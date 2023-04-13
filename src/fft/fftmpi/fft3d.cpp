/* ----------------------------------------------------------------------
   fftMPI - library for computing 3d/2d FFTs in parallel
   http://fftmpi.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright 2018 National Technology & Engineering Solutions of
   Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
   NTESS, the U.S. Government retains certain rights in this software.
   This software is distributed under the modified Berkeley Software
   Distribution (BSD) License.

   See the README file in the top-level fftMPI directory.
------------------------------------------------------------------------- */

// FFT3d class

/* ----------------------------------------------------------------------
   Contributing authors: Steve Plimpton (Sandia)
                         Jim Shepherd (GA Tech) added SGI SCSL support
                         Axel Kohlmeyer (Temple U) added support for
                           FFTW3, KISSFFT, Dfti/MKL, and ACML
                         Phil Blood (PSC) added single precision FFT
                         Paul Coffman (IBM) added MPI collective remap
------------------------------------------------------------------------- */

#include <mpi.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fft3d.h"
#include "remap3d.h"
#include "version.h"
#include "memory.h"
#include "error.h"

using namespace FFTMPI_NS;

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

#define NFACTOR 50
#define BIG 1.0e20

#ifdef FFT_LONGLONG_TO_LONG
#define MPI_FFT_BIGINT MPI_LONG
#else
#define MPI_FFT_BIGINT MPI_LONG_LONG
#endif

typedef int64_t bigint;

/* ----------------------------------------------------------------------
   data layout for 3d FFTs:

   data set of Nfast x Nmid x Nslow elements is owned by P procs
   on input, each proc owns a subsectioffn of the elements
   on output, each proc will own a (possibly different) subsection
   my subsection must not overlap with any other proc's subsection,
     i.e. the union of all proc's input (or output) subsections must
     exactly tile the global Nfast x Nmid x Nslow data set
   when called from C, all subsection indices are
     C-style from 0 to N-1 where N = Nfast or Nmid or Nslow
   when called from F77, all subsection indices are
     F77-style from 1 to N where N = Nfast or Nmid or Nslow
   a proc can own 0 elements on input or output
     by specifying hi index < lo index
   on both input and output, data is stored contiguously on a processor
     with a fast-varying, mid-varying, and slow-varying index

   flags caller can set before setup()
     collective = 0/1/2 (default = 2)
       collective MPI operations for remapping data
       0 = point to point comm for all remaps
       1 = MPI all2all for all remaps
       2 = point for pencil2brick remaps, all2all for pencil2pencil remaps
     exchange = 0/1 (default = 0)
       style of data exchanges
       0 = pencil to pencil
       1 = brick to pencil and pencil to brick
     packflag = 0/1/2 (default = 2)
       style of pack/unpack methods for remapping data
       0 = array
       1 = pointer
       2 = memcpy
     memoryflag = 0/1 (default = 1)
       caller provides remap memory or system does
       0 = user provides memory
       1 = system provides memory internally

   flags caller can set before compute()
     scaling = 0/1 (default = 1)
       0 = no scaling
       1 = scaling of forward FFT
     remaponly = 0/1 (default = 0)
       0 = regular FFT
       1 = only the remap operations, no 1d FFTs
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   instantiate a 3d FFT
   user_comm = MPI communicator for the P procs which own the data
   user_precision = caller precision for its FFT data, 1=single, 2=double
------------------------------------------------------------------------- */

FFT3d::FFT3d(MPI_Comm user_comm, int user_precision)
{
  world = user_comm;
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // default settings
  // user must change them before setup()

  collective = 2;
  exchange = 0;
  packflag = 2;
  memoryflag = 1;

  // default settings
  // user can change them before compute()

  scaled = 1;
  remaponly = 0;

  // tuning results

  ntrial = npertrial = 0;
  cbest = ebest = pbest = -1;
  besttime = 0.0;

  // Memory and Error classes

  memory = new Memory();
  error = new Error(world);

  // error if caller and lib precision do not match

  if (user_precision != 1 && user_precision != 2)
    error->all("Invalid precision setting for FFT lib");

  if ((user_precision == 1 && sizeof(FFT_SCALAR) != 4) ||
      (user_precision == 2 && sizeof(FFT_SCALAR) != 8))
    error->all("Precision mis-match with FFT lib");

  if (sizeof(FFT_SCALAR) == 4) precision = "single";
  if (sizeof(FFT_SCALAR) == 8) precision = "double";

  // allowed prime factors for each FFT grid dimension

  nprime = 3;
  primes = new int[nprime];
  primes[0] = 2;
  primes[1] = 3;
  primes[2] = 5;

  // initialize memory allocations

  factors = new int[NFACTOR];
  nfactor = 0;

  remap_prefast = remap_fastmid = remap_midslow = remap_postslow = NULL;
  remap_preslow = remap_slowmid = remap_midfast = remap_postfast = NULL;
  remap_premid = remap_fastslow = remap_slowfast = remap_postmid = NULL;
  fft_fast = fft_mid = fft_slow = NULL;

  memusage = 0;
  sendbuf = recvbuf = NULL;

  setupflag = 0;
  setup_memory_flag = 0;
}

/* ----------------------------------------------------------------------
   delete a 3d FFT
------------------------------------------------------------------------- */

FFT3d::~FFT3d()
{
  delete memory;
  delete error;

  delete [] primes;
  delete [] factors;

  if (setupflag) deallocate_setup();
  if (memoryflag) deallocate_setup_memory();
}

/* ----------------------------------------------------------------------
   create plan for performing a 3d FFT

   inputs:
   nfast,nmid,nslow     size of global 3d matrix
   in_ilo,in_ihi        input bounds of data I own in fast index
   in_jlo,in_jhi        input bounds of data I own in mid index
   in_klo,in_khi        input bounds of data I own in slow index
   out_ilo,out_ihi      output bounds of data I own in fast index
   out_jlo,out_jhi      output bounds of data I own in mid index
   out_klo,out_khi      output bounds of data I own in slow index
   permute              permutation in storage order of indices on output
                          0 = no permutation
                          1 = permute once = mid->fast, slow->mid, fast->slow
                          2 = permute twice = slow->fast, fast->mid, mid->slow

   outputs:
   fftsize = size of in/out FFT arrays required from caller
     can be greater than insize, due to remaps
   sendsize = size of send buffer, caller may choose to provide it
   recvsize = size of recv buffer, caller may choose to provide it
------------------------------------------------------------------------- */

void FFT3d::setup(int user_nfast, int user_nmid, int user_nslow,
                  int user_in_ilo, int user_in_ihi, int user_in_jlo,
                  int user_in_jhi, int user_in_klo, int user_in_khi,
                  int user_out_ilo, int user_out_ihi, int user_out_jlo,
                  int user_out_jhi, int user_out_klo, int user_out_khi,
                  int user_permute,
                  int &user_fftsize, int &user_sendsize, int &user_recvsize)
{
  int flag,allflag;

  if (setupflag) error->all("FFT is already setup");
  setupflag = 1;

  // internal copies of input params

  nfast = user_nfast;
  nmid = user_nmid;
  nslow = user_nslow;

  in_ilo = user_in_ilo; in_ihi = user_in_ihi;
  in_jlo = user_in_jlo; in_jhi = user_in_jhi;
  in_klo = user_in_klo; in_khi = user_in_khi;
  out_ilo = user_out_ilo; out_ihi = user_out_ihi;
  out_jlo = user_out_jlo; out_jhi = user_out_jhi;
  out_klo = user_out_klo; out_khi = user_out_khi;

  permute = user_permute;

  // all dimensions must be >= 2
  // all dimensions must be factorable

  if (nfast < 2 || nmid < 2 || nslow < 2)
    error->all("Each FFT dimension must be >= 2");

  if (!prime_factorable(nfast)) error->all("Invalid nfast");
  if (!prime_factorable(nmid)) error->all("Invalid nmid");
  if (!prime_factorable(nslow)) error->all("Invalid nslow");

  // error checks in indices and tiling

  flag = 0;
  if (in_ilo > in_ihi+1 || in_jlo > in_jhi+1 || in_klo > in_khi+1) flag = 1;
  if (in_ilo < 0 || in_jlo < 0 || in_klo < 0) flag = 1;
  if (in_ihi >= nfast || in_jhi >= nmid || in_khi >= nslow) flag = 1;

  MPI_Allreduce(&flag,&allflag,1,MPI_INT,MPI_MAX,world);

  if (allflag) error->all("FFT setup in/out indices are invalid");

  bigint n = (bigint) (in_ihi-in_ilo+1) * (in_jhi-in_jlo+1) * (in_khi-in_klo+1);
  bigint nall;
  MPI_Allreduce(&n,&nall,1,MPI_FFT_BIGINT,MPI_SUM,world);

  if (nall != ((bigint) nfast * nmid*nslow))
    error->all("FFT setup in/out indices do not tile grid");

  // set collective flags for different remap operations
  // bp = brick2pencil or pencil2brick, pp = pencel2pencil

  if (collective == 0) collective_bp = collective_pp = 0;
  else if (collective == 1) collective_bp = collective_pp = 1;
  else {
    collective_bp = 0;
    collective_pp = 1;
  }

  // inout_layout_same is set only if:
  // in/out indices are same on every proc and permute = 0

  flag = 0;
  if (in_ilo != out_ilo || in_ihi != out_ihi ||
      in_jlo != out_jlo || in_jhi != out_jhi ||
      in_klo != out_klo || in_khi != out_khi) flag = 1;
  if (permute) flag = 1;
  MPI_Allreduce(&flag,&allflag,1,MPI_INT,MPI_MAX,world);
  if (allflag) inout_layout_same = 0;
  else inout_layout_same = 1;

  // compute partitioning of FFT grid across procs for each pencil layout
  // if exchange set, also partition in 3d for brick layout
  // np = # of procs in each dimension
  // ip = my location in each dimension

  factor(nprocs);
  if (nfactor > NFACTOR) error->all("Nprocs is too large to factor");

  procfactors(1,nmid,nslow,
              npfast1,npfast2,npfast3,ipfast1,ipfast2,ipfast3);
  procfactors(nfast,1,nslow,
              npmid1,npmid2,npmid3,ipmid1,ipmid2,ipmid3);
  procfactors(nfast,nmid,1,
              npslow1,npslow2,npslow3,ipslow1,ipslow2,ipslow3);

  if (exchange)
    procfactors(nfast,nmid,nslow,
                npbrick1,npbrick2,npbrick3,ipbrick1,ipbrick2,ipbrick3);
  else npbrick1 = npbrick2 = npbrick3 = 0;

  // remap from initial layout to fast pencil layout
  // remap_preflag = 1 if remap is needed, else 0
  // not needed if all procs own entire fast dimension initially
  // fast indices = data layout before/after 1st set of FFTs

  if (in_ilo == 0 && in_ihi == nfast-1) flag = 0;
  else flag = 1;
  MPI_Allreduce(&flag,&allflag,1,MPI_INT,MPI_MAX,world);

  if (allflag == 0) {
    remap_preflag = 0;
    fast_ilo = in_ilo;
    fast_ihi = in_ihi;
    fast_jlo = in_jlo;
    fast_jhi = in_jhi;
    fast_klo = in_klo;
    fast_khi = in_khi;
  } else {
    remap_preflag = 1;
    fast_ilo = 0;
    fast_ihi = nfast - 1;
    fast_jlo = ipfast2*nmid/npfast2;
    fast_jhi = (ipfast2+1)*nmid/npfast2 - 1;
    fast_klo = ipfast3*nslow/npfast3;
    fast_khi = (ipfast3+1)*nslow/npfast3 - 1;
  }

  // remap from fast pencil layout to mid pencil layout
  // always needed, b/c permutation changes
  // mid indices = data layout before/after 2nd set of FFTs

  mid_ilo = ipmid1*nfast/npmid1;
  mid_ihi = (ipmid1+1)*nfast/npmid1 - 1;
  mid_jlo = 0;
  mid_jhi = nmid - 1;
  mid_klo = ipmid3*nslow/npmid3;
  mid_khi = (ipmid3+1)*nslow/npmid3 - 1;

  // remap from mid pencil layout to slow pencil layout
  // always needed, b/c permutation changes
  // slow indices = data layout before/after 3rd set of FFTs
  // if final layout is slow pencil with permute=2, set slow = out

  if (permute == 2 && out_klo == 0 && out_khi == nslow-1) flag = 0;
  else flag = 1;
  MPI_Allreduce(&flag,&allflag,1,MPI_INT,MPI_MAX,world);

  if (allflag == 0) {
    slow_ilo = out_ilo;
    slow_ihi = out_ihi;
    slow_jlo = out_jlo;
    slow_jhi = out_jhi;
    slow_klo = out_klo;
    slow_khi = out_khi;
  } else {
    slow_ilo = ipslow1*nfast/npslow1;
    slow_ihi = (ipslow1+1)*nfast/npslow1 - 1;
    slow_jlo = ipslow2*nmid/npslow2;
    slow_jhi = (ipslow2+1)*nmid/npslow2 - 1;
    slow_klo = 0;
    slow_khi = nslow - 1;
  }

  // remap from slow pencil layout to final layout
  // remap_postflag = 1 if remap is needed, else 0
  // not needed if permute=2 and slow = out already

  if (permute == 2 &&
      out_ilo == slow_ilo && out_ihi == slow_ihi &&
      out_jlo == slow_jlo && out_jhi == slow_jhi &&
      out_klo == slow_klo && out_khi == slow_khi) flag = 0;
  else flag = 1;
  MPI_Allreduce(&flag,&allflag,1,MPI_INT,MPI_MAX,world);

  if (allflag == 0) remap_postflag = 0;
  else remap_postflag = 1;

  // if exchange is set, then remap for fast/mid and mid/slow
  // remap will be two stages, with brick layout and brick indices inbetween

  if (exchange) {
    brick_ilo = ipbrick1*nfast/npbrick1;
    brick_ihi = (ipbrick1+1)*nfast/npbrick1 - 1;
    brick_jlo = ipbrick2*nmid/npbrick2;
    brick_jhi = (ipbrick2+1)*nmid/npbrick2 - 1;
    brick_klo = ipbrick3*nslow/npbrick3;
    brick_khi = (ipbrick3+1)*nslow/npbrick3 - 1;
  }

  // create Remap instances for 4 forward remaps
  // likewise for inverse remaps if in/out layout is not the same
  // create calls return max size of send/recv buffers needed by remaps

  sendsize = recvsize = 0;
  remap_forward_create(sendsize,recvsize);
  if (!inout_layout_same) remap_inverse_create(sendsize,recvsize);

  // insize/outsize = # of FFT data points in initial/final layout
  // fastsize/midsize/slowsize = # of data points in fast/mid/slow layout
  // maxsize = max of all these sizes, returned to caller

  insize = (in_ihi-in_ilo+1) * (in_jhi-in_jlo+1) *
    (in_khi-in_klo+1);
  outsize = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1) *
    (out_khi-out_klo+1);

  fastsize = (fast_ihi-fast_ilo+1) * (fast_jhi-fast_jlo+1) *
    (fast_khi-fast_klo+1);
  midsize = (mid_ihi-mid_ilo+1) * (mid_jhi-mid_jlo+1) *
    (mid_khi-mid_klo+1);
  slowsize = (slow_ihi-slow_ilo+1) * (slow_jhi-slow_jlo+1) *
    (slow_khi-slow_klo+1);
  if (exchange)
    bricksize = (brick_ihi-brick_ilo+1) * (brick_jhi-brick_jlo+1) *
      (brick_khi-brick_klo+1);

  fftsize = MAX(insize,outsize);
  fftsize = MAX(fftsize,fastsize);
  fftsize = MAX(fftsize,midsize);
  fftsize = MAX(fftsize,slowsize);
  if (exchange) fftsize = MAX(fftsize,bricksize);

  // setup for 3 sets of 1d FFTs, also scaling normalization
  // outsize must be already set for setup_ffts() to use to setup scaling
  // norm must allow for nfast*nmid*nslow to exceed a 4-byte int (2B)

  fft_fast = new FFT1d;
  fft_mid = new FFT1d;
  fft_slow = new FFT1d;

  fft_fast->length = nfast;
  fft_fast->n = (fast_jhi-fast_jlo+1) * (fast_khi-fast_klo+1);
  fft_fast->total = fft_fast->n * fft_fast->length;

  fft_mid->length = nmid;
  fft_mid->n = (mid_ihi-mid_ilo+1) * (mid_khi-mid_klo+1);
  fft_mid->total = fft_mid->n * fft_mid->length;

  fft_slow->length = nslow;
  fft_slow->n = (slow_ihi-slow_ilo+1) * (slow_jhi-slow_jlo+1);
  fft_slow->total = fft_slow->n * fft_slow->length;

  setup_ffts();

  norm = 1.0/((bigint) nfast * nmid*nslow);
  normnum = outsize;

  // allocate sendbuf, recvbuf arrays to max sizes needed by any remap

  if (memoryflag) {
    setup_memory_flag = 1;
    if (sendsize) {
      sendbuf = (FFT_SCALAR *) memory->smalloc(sendsize*sizeof(FFT_SCALAR));
      if (!sendbuf) error->one("Could not allocate sendbuf array");
    }
    if (recvsize) {
      recvbuf = (FFT_SCALAR *) memory->smalloc(recvsize*sizeof(FFT_SCALAR));
      if (!recvbuf) error->one("Could not allocate recvbuf array");
    }
  }

  // return buffer sizes to caller

  user_fftsize = fftsize;
  user_sendsize = sendsize;
  user_recvsize = recvsize;

  // set memusage for FFT and Remap memory

  memusage = 0;

  if (memoryflag) {
    memusage += (int64_t) sendsize * sizeof(FFT_SCALAR);
    memusage += (int64_t) recvsize * sizeof(FFT_SCALAR);
  }

  memusage += remap_memory();
}

/* ----------------------------------------------------------------------
   deallocate memory allocated by setup()
------------------------------------------------------------------------- */

void FFT3d::deallocate_setup()
{
  setupflag = 0;

  deallocate_remap(remap_prefast);
  deallocate_remap(remap_fastmid);
  deallocate_remap(remap_midslow);
  deallocate_remap(remap_postslow);

  deallocate_remap(remap_preslow);
  deallocate_remap(remap_slowmid);
  deallocate_remap(remap_midfast);
  deallocate_remap(remap_postfast);

  deallocate_remap(remap_premid);
  deallocate_remap(remap_fastslow);
  deallocate_remap(remap_slowfast);
  deallocate_remap(remap_postmid);

  deallocate_ffts();
  delete fft_fast;
  delete fft_mid;
  delete fft_slow;

  remap_prefast = remap_fastmid = remap_midslow = remap_postslow = NULL;
  remap_preslow = remap_slowmid = remap_midfast = remap_postfast = NULL;
  remap_premid = remap_fastslow = remap_slowfast = remap_postmid = NULL;

  fft_fast = fft_mid = fft_slow = NULL;
}

/* ----------------------------------------------------------------------
   pass in user memory for Remap send/recv operations
   user_sendbuf = send buffer of length user_sendsize
   user_recvbuf = send buffer of length user_recvsize
------------------------------------------------------------------------- */

void FFT3d::setup_memory(FFT_SCALAR *user_sendbuf, FFT_SCALAR *user_recvbuf)
{
  if (!setupflag) error->all("Cannot setup FFT memory before setup");
  if (memoryflag) error-> all("Cannot setup FFT memory with memoryflag set");

  setup_memory_flag = 1;
  sendbuf = user_sendbuf;
  recvbuf = user_recvbuf;
}

/* ----------------------------------------------------------------------
   deallocate memory allocated internally for send/recv
   only called if allocated internally
------------------------------------------------------------------------- */

void FFT3d::deallocate_setup_memory()
{
  setup_memory_flag = 0;
  memory->sfree(sendbuf);
  memory->sfree(recvbuf);
  sendbuf = recvbuf = NULL;
}

/* ----------------------------------------------------------------------
   perform a 3d FFT
   in           address of input data on this proc
   out          address of output data for this proc, can be same as in
   flag         1 for forward FFT, -1 for inverse FFT
------------------------------------------------------------------------- */

void FFT3d::compute(FFT_SCALAR *in, FFT_SCALAR *out, int flag)
{
  if (!setupflag) error->all("Cannot compute FFT before setup");
  if (!setup_memory_flag) error->all("Cannot compute FFT before setup_memory");

  FFT_SCALAR *data = out;

  if (flag == 1 || inout_layout_same) {

    if (remap_prefast) remap(in,out,remap_prefast);
    else if (in != out) memcpy(out,in,insize*sizeof(FFT_SCALAR));

    if (remaponly) {
      if (remap_fastmid) remap(data,data,remap_fastmid);
      if (remap_midslow) remap(data,data,remap_midslow);
    } else {
      perform_ffts((FFT_DATA *) data,flag,fft_fast);
      if (remap_fastmid) remap(data,data,remap_fastmid);
      perform_ffts((FFT_DATA *) data,flag,fft_mid);
      if (remap_midslow) remap(data,data,remap_midslow);
      perform_ffts((FFT_DATA *) data,flag,fft_slow);
    }

    if (remap_postslow) remap(data,data,remap_postslow);

    if (flag == 1 && scaled && !remaponly) scale_ffts((FFT_DATA *) data);

  } else {

    if (remap_preslow) remap(in,out,remap_preslow);
    else if (in != out) memcpy(out,in,outsize*sizeof(FFT_SCALAR));

    if (remaponly) {
      if (remap_slowmid) remap(data,data,remap_slowmid);
      if (remap_midfast) remap(data,data,remap_midfast);
    } else {
      perform_ffts((FFT_DATA *) data,flag,fft_slow);
      if (remap_slowmid) remap(data,data,remap_slowmid);
      perform_ffts((FFT_DATA *) data,flag,fft_mid);
      if (remap_midfast) remap(data,data,remap_midfast);
      perform_ffts((FFT_DATA *) data,flag,fft_fast);
    }

    if (remap_postfast) remap(in,out,remap_postfast);
  }
}

/* ----------------------------------------------------------------------
   perform just the 1d FFTs needed by a 3d FFT, no data movement
     useful for timing purposes
   in           starting address of input data on this proc, all set to 0.0
   flag         1 for forward FFT, -1 for inverse FFT
------------------------------------------------------------------------- */

void FFT3d::only_1d_ffts(FFT_SCALAR *in, int flag)
{
  if (!setupflag) error->all("Cannot compute 1d FFTs before setup");

  FFT_DATA *data = (FFT_DATA *) in;
  perform_ffts((FFT_DATA *) data,flag,fft_fast);
  perform_ffts((FFT_DATA *) data,flag,fft_mid);
  perform_ffts((FFT_DATA *) data,flag,fft_slow);
}

/* ----------------------------------------------------------------------
   perform all the remaps in a 3d FFT, but no 1d FFTs
     useful for debugging and timing purposes
   in           starting address of input data on this proc
   out          starting address of where output data for this proc
                  will be placed (can be same as in)
   flag         1 for forward FFT, -1 for inverse FFT
------------------------------------------------------------------------- */

void FFT3d::only_remaps(FFT_SCALAR *in, FFT_SCALAR *out, int flag)
{
  if (!setupflag) error->all("Cannot perform FFT remap before setup");
  if (!setup_memory_flag)
    error->all("Cannot perform FFT remap before setup_memory");

  FFT_SCALAR *data = out;

  if (flag == 1 || inout_layout_same) {

    if (remap_prefast) remap(in,out,remap_prefast);
    else if (in != out) memcpy(out,in,insize*sizeof(FFT_SCALAR));

    if (remap_fastmid) remap(data,data,remap_fastmid);
    if (remap_midslow) remap(data,data,remap_midslow);

    if (remap_postslow) remap(data,data,remap_postslow);

  } else {

    if (remap_preslow) remap(in,out,remap_preslow);
    else if (in != out) memcpy(out,in,outsize*sizeof(FFT_SCALAR));

    if (remap_slowmid) remap(data,data,remap_slowmid);
    if (remap_midfast) remap(data,data,remap_midfast);

    if (remap_postfast) remap(data,data,remap_postfast);
  }
}

/* ----------------------------------------------------------------------
   perform just a single remap operation
     useful for debugging and timing purposes
   in           starting address of input data on this proc, all set to 0.0
   out          starting address of where output data for this proc
                  will be placed (can be same as in)
   flag         1 for forward FFT, -1 for inverse FFT
   which        which remap to perform = 1,2,3,4
------------------------------------------------------------------------- */

void FFT3d::only_one_remap(FFT_SCALAR *in, FFT_SCALAR *out, int flag, int which)
{
  if (!setupflag) error->all("Cannot perform an FFT remap before setup");
  if (!setup_memory_flag)
    error->all("Cannot perform an FFT remap before setup_memory");

  if (flag == 1 || inout_layout_same) {
    if (which == 1) {
      if (remap_prefast) remap(in,out,remap_prefast);
      else if (in != out) memcpy(out,in,insize*sizeof(FFT_SCALAR));
    } else if (which == 2) {
      if (remap_fastmid) remap(in,out,remap_fastmid);
    } else if (which == 3) {
      if (remap_midslow) remap(in,out,remap_midslow);
    } else if (which == 4) {
      if (remap_postslow) remap(in,out,remap_postslow);
    }

  } else {
    if (which == 4) {
      if (remap_preslow) remap(in,out,remap_preslow);
      else if (in != out) memcpy(out,in,outsize*sizeof(FFT_SCALAR));
    } else if (which == 3) {
      if (remap_slowmid) remap(in,out,remap_slowmid);
    } else if (which == 2) {
      if (remap_midfast) remap(in,out,remap_midfast);
    } else if (which == 1) {
      if (remap_postfast) remap(in,out,remap_postfast);
    }
  }
}

/* ----------------------------------------------------------------------
   tune settings for fastest FFT: collective, exchange, pack flags
   this SETS collective/exchange/packflag to optimal values
   initial args are the same as for setup()
   additional args:
   flag = 1 or -1 to tune just forward or inverse, 0 to tune both together
   niter = perform niter FFTs in each tuning trial
   tmax = do not exceed tmax CPU secs for all tuning trials
     tmax = 0.0 to allow unlimited time
   tflag = 1 to also time 1d FFTs and remap operations, 0 if not
   output:
     collective/exchange/packflag are set to optimal values
     these internal variables can be accessed by caller:
       besttime = fastest CPU time for a single FFT
       ntrial = # of time trials performed
       npertrial = # of FFTs performed per trial (2x this for flag = 0)
       cbest,ebest,pbest = collective/exchange/pack settings for fastest
       cflags = collective setting for each trial
       eflags = exchange setting for each trial
       pflags = packflag setting for each trial
       t3d = CPU time for single 3d FFT for each trial
       t1d = CPU time for 1d FFTs for each trial
       tremap = CPU time for all remaps for each trial
       tremap1234 = CPU time for one remap (1-4) for each trial
------------------------------------------------------------------------- */

void FFT3d::tune(int user_nfast, int user_nmid, int user_nslow,
     int user_in_ilo, int user_in_ihi, int user_in_jlo,
     int user_in_jhi, int user_in_klo, int user_in_khi,
     int user_out_ilo, int user_out_ihi, int user_out_jlo,
     int user_out_jhi, int user_out_klo, int user_out_khi,
     int user_permute,
     int &user_fftsize, int &user_sendsize, int &user_recvsize,
     int flag, int niter, double tmax, int tflag)
{
  if (setupflag) error->all("FFT is already setup");
  if (flag < -1 || flag > 1) error->all("Invalid flag arg for FFT tune");
  if (niter <= 0 || tmax < 0.0)
    error->all("Invalid niter/tmax args for FFT tune");
  if (tflag < 0 || tflag > 1) error->all("Invalid tflag arg for FFT tune");

  // preserve user memoryflag setting since will override below

  int user_memoryflag = memoryflag;

  // local data memory for tune operation

  int maxfftsize = 0;
  FFT_SCALAR *data = NULL;

  // perform a setup and single FFT iteration using current settings
  // do not use tflag for this tune_trial()

  memoryflag = 1;

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  setup(user_nfast,user_nmid,user_nslow,
        user_in_ilo,user_in_ihi,user_in_jlo,
        user_in_jhi,user_in_klo,user_in_khi,
        user_out_ilo,user_out_ihi,user_out_jlo,
        user_out_jhi,user_out_klo,user_out_khi,
        user_permute,user_fftsize,user_sendsize,user_recvsize);

  MPI_Barrier(world);
  double time2 = MPI_Wtime();
  double timesetup = time2-time1;

  if (user_fftsize > maxfftsize) {
    data = (FFT_SCALAR *)
      memory->srealloc(data,user_fftsize*2*sizeof(FFT_SCALAR));
    for (int i = 2*maxfftsize; i < 2*user_fftsize; i++) data[i] = 0.0;
    maxfftsize = user_fftsize;
  }

  ntrial = 0;
  tune_trial(data,1,flag,0,tfft[ntrial],t1d[ntrial],tremap[ntrial],
             tremap1[ntrial],tremap2[ntrial],tremap3[ntrial],tremap4[ntrial]);
  double timefft = tfft[ntrial];
  cflags[ntrial] = collective;
  eflags[ntrial] = exchange;
  pflags[ntrial] = packflag;

  deallocate_setup();
  if (memoryflag) deallocate_setup_memory();

  // insure all procs use exact same timesetup & timefft for setup

  MPI_Bcast(&timesetup,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&timefft,1,MPI_DOUBLE,0,world);

  // for tmax = time limit > 0.0
  // use single iteration timing to setup sequence of tuning runs
  // reset niter if necessary to limit tuning time to tmax
  // ntrial = # of test runs, limit if even niter = 1 exceeds tmax
  // 10 = initial-run + (3-collective * 2-exchange) + 3-packflag

  double timetotal;
  int nruns = 10;
  if (tflag) timetotal = nruns * (timesetup + 3*niter*timefft);
  else timetotal = nruns * (timesetup + niter*timefft);

  if (tmax > 0.0 && timetotal > tmax) {
    if (tflag) niter = static_cast<int> ((tmax/nruns-timesetup) / timefft/3);
    else niter = static_cast<int> ((tmax/nruns-timesetup) / timefft);
    if (niter <= 0) {
      tflag = 0;
      niter = static_cast<int> ((tmax/nruns-timesetup) / timefft);
      if (niter <= 0) {
        niter = 1;
        nruns = static_cast<int> (tmax/(timesetup+timefft));
      }
    }
  }

  // use nruns to set bounds for collective/exchange/packflag tests

  int cdefault = 2;
  int edefault = 0;
  int pdefault = 2;
  ntrial = 0;

  int cstart,cstop,estart,estop,pstart,pstop;

  if (nruns-1 >= 9) {
    cstart = 0;
    cstop = 3;
    estart = 0;
    estop = 2;
    pstart = 0;
    pstop = 3;
  } else if (nruns-1 >= 6) {
    cstart = 0;
    cstop = 3;
    estart = 0;
    estop = 2;
    pstart = -1;
    pbest = pdefault;
  } else if (nruns-1 >= 3) {
    cstart = 0;
    cstop = 3;
    estart = edefault;
    estop = edefault+1;
    pstart = -1;
    pbest = pdefault;
  } else if (nruns-1 == 2) {
    cstart = 1;
    cstop = 3;
    estart = edefault;
    estop = edefault+1;
    pstart = -1;
    pbest = pdefault;

  // if only time for one more run
  //   and initial run was not with defaults:
  // perform additional default run

  } else if (nruns-1 == 1 && (collective != cdefault ||
                              exchange != edefault || packflag != pdefault)) {
    ntrial = 1;
    cstart = cdefault;
    cstop = cdefault+1;
    estart = edefault;
    estop = edefault+1;
    pstart = -1;
    pbest = pdefault;

  // initial run is all there is time to perform

  } else {
    ntrial = 1;
    cstart = estart = pstart = -1;
    besttime = timefft;
    cbest = cflags[0] = collective;
    ebest = eflags[0] = exchange;
    pbest = pflags[0] = packflag;
  }

  // find best time for collective/exchange settings
  // perform a new setup() for each pair of settings
  // use packflag = pdefault, may optimize it below
  // only do this operation if cstart >= 0

  if (cstart >= 0) {
    besttime = BIG;

    for (int cflag = cstart; cflag < cstop; cflag++) {
      for (int eflag = estart; eflag < estop; eflag++) {
        collective = cflags[ntrial] = cflag;
        exchange = eflags[ntrial] = eflag;
        packflag = pflags[ntrial] = pdefault;
        memoryflag = 1;
        setup(user_nfast,user_nmid,user_nslow,
              user_in_ilo,user_in_ihi,user_in_jlo,
              user_in_jhi,user_in_klo,user_in_khi,
              user_out_ilo,user_out_ihi,user_out_jlo,
              user_out_jhi,user_out_klo,user_out_khi,
              user_permute,user_fftsize,user_sendsize,user_recvsize);
        if (user_fftsize > maxfftsize) {
          data = (FFT_SCALAR *)
            memory->srealloc(data,user_fftsize*2*sizeof(FFT_SCALAR));
          for (int i = 2*maxfftsize; i < 2*user_fftsize; i++) data[i] = 0.0;
          maxfftsize = user_fftsize;
        }

        tune_trial(data,niter,flag,tflag,tfft[ntrial],t1d[ntrial],tremap[ntrial],
                   tremap1[ntrial],tremap2[ntrial],
                   tremap3[ntrial],tremap4[ntrial]);

        deallocate_setup();
        if (memoryflag) deallocate_setup_memory();

        if (tfft[ntrial] < besttime) {
          besttime = tfft[ntrial];
          cbest = cflag;
          ebest = eflag;
        }

        ntrial++;
      }
    }
  }

  // find best time for packflag setting, using best collective/exchange
  // perform a new setup() for each setting
  // only do this operation if pstart >= 0

  if (pstart >= 0) {
    besttime = BIG;

    for (int pflag = pstart; pflag < pstop; pflag++) {
      collective = cflags[ntrial] = cbest;
      exchange = eflags[ntrial] = ebest;
      packflag = pflags[ntrial] = pflag;
      memoryflag = 1;
      setup(user_nfast,user_nmid,user_nslow,
            user_in_ilo,user_in_ihi,user_in_jlo,
            user_in_jhi,user_in_klo,user_in_khi,
            user_out_ilo,user_out_ihi,user_out_jlo,
            user_out_jhi,user_out_klo,user_out_khi,
            user_permute,user_fftsize,user_sendsize,user_recvsize);
      if (user_fftsize > maxfftsize) {
        data = (FFT_SCALAR *)
          memory->srealloc(data,user_fftsize*2*sizeof(FFT_SCALAR));
        for (int i = 2*maxfftsize; i < 2*user_fftsize; i++) data[i] = 0.0;
        maxfftsize = user_fftsize;
      }

      tune_trial(data,niter,flag,tflag,tfft[ntrial],t1d[ntrial],tremap[ntrial],
                 tremap1[ntrial],tremap2[ntrial],
                 tremap3[ntrial],tremap4[ntrial]);

      deallocate_setup();
      if (memoryflag) deallocate_setup_memory();

      if (tfft[ntrial] < besttime) {
        besttime = tfft[ntrial];
        pbest = pflag;
      }

      ntrial++;
    }
  }

  // scale all times for a single FFT

  npertrial = niter;
  double scale = 1.0/npertrial;
  if (flag == 0) scale *= 0.5;

  besttime *= scale;
  for (int i = 0; i < ntrial; i++) {
    tfft[i] *= scale;
    t1d[i] *= scale;
    tremap[i] *= scale;
    tremap1[i] *= scale;
    tremap2[i] *= scale;
    tremap3[i] *= scale;
    tremap4[i] *= scale;
  }

  // free local data

  memory->sfree(data);

  // final setup with optimal settings
  // restore original memoryflag
  // setuptime = CPU time to perform final setup

  collective = cbest;
  exchange = ebest;
  packflag = pbest;
  memoryflag = user_memoryflag;

  MPI_Barrier(world);
  time1 = MPI_Wtime();

  setup(user_nfast, user_nmid, user_nslow,
  user_in_ilo, user_in_ihi, user_in_jlo,
  user_in_jhi, user_in_klo, user_in_khi,
  user_out_ilo, user_out_ihi, user_out_jlo,
  user_out_jhi, user_out_klo, user_out_khi,
  user_permute, user_fftsize, user_sendsize, user_recvsize);

  MPI_Barrier(world);
  time2 = MPI_Wtime();
  setuptime = time2-time1;
}

/* ----------------------------------------------------------------------
   perform a timing trial for tune() method
   inputs:
   data = input data, all set to zeroes
   nper = # of iterations in trial
   flag = 1,-1 for forward/reverse FFT, 0 for both
   tflag = 1 if also test 1d FFTs and remaps
   return:
   time for 3d FFT, 1d FFTs, all remaps, each of 4 remaps
------------------------------------------------------------------------- */

void FFT3d::tune_trial(FFT_SCALAR *data, int nper, int flag, int tflag,
                       double &time3d, double &time1d, double &timeremap,
                       double &timeremap1, double &timeremap2,
                       double &timeremap3, double &timeremap4)
{
  double time1,time2;

  MPI_Barrier(world);
  time1 = MPI_Wtime();

  if (flag)
    for (int i = 0; i < nper; i++)
      compute(data,data,flag);
  else {
    for (int i = 0; i < nper; i++) {
      compute(data,data,1);
      compute(data,data,-1);
    }
  }

  MPI_Barrier(world);
  time2 = MPI_Wtime();
  time3d = time2-time1;

  if (tflag) {
    if (flag)
      for (int i = 0; i < nper; i++)
        only_1d_ffts(data,flag);
    else {
      for (int i = 0; i < nper; i++) {
        only_1d_ffts(data,1);
        only_1d_ffts(data,-1);
      }
    }

    MPI_Barrier(world);
    time1 = MPI_Wtime();
    time1d = time1-time2;

    if (flag)
      for (int i = 0; i < nper; i++)
        only_remaps(data,data,flag);
    else {
      for (int i = 0; i < nper; i++) {
        only_remaps(data,data,1);
        only_remaps(data,data,-1);
      }
    }

    MPI_Barrier(world);
    time2 = MPI_Wtime();
    timeremap = time2-time1;

    for (int which = 1; which <= 4; which++) {
      if (flag)
        for (int i = 0; i < nper; i++)
          only_one_remap(data,data,flag,which);
      else {
        for (int i = 0; i < nper; i++) {
          only_one_remap(data,data,1,which);
          only_one_remap(data,data,-1,which);
        }
      }

      MPI_Barrier(world);
      time1 = MPI_Wtime();
      if (which == 1) timeremap1 = time1-time2;
      else if (which == 2) timeremap2 = time1-time2;
      else if (which == 3) timeremap3 = time1-time2;
      else if (which == 4) timeremap4 = time1-time2;
      time2 = time1;
    }

  } else time1d = timeremap =
           timeremap1 = timeremap2 = timeremap3 = timeremap4 = 0.0;
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------
// private methods
// -------------------------------------------------------------------
// -------------------------------------------------------------------

/* ----------------------------------------------------------------------
   perform a 3d remap
   in           address of input data on this proc
   out          address of output data for this proc, can be same as in
   plan         plan for remap
------------------------------------------------------------------------- */

void FFT3d::remap(FFT_SCALAR *in, FFT_SCALAR *out, Remap *remap)
{
  remap->remap3d->remap(in,out,sendbuf,recvbuf);
  if (remap->remap3d_extra)
    remap->remap3d_extra->remap(in,out,sendbuf,recvbuf);
}

/* ----------------------------------------------------------------------
   dellocate a Remap and its contents
------------------------------------------------------------------------- */

void FFT3d::deallocate_remap(Remap *remap)
{
  if (remap == NULL) return;
  delete remap->remap3d;
  delete remap->remap3d_extra;
  delete remap;
}

/* ----------------------------------------------------------------------
   create Remap3d instances for forward FFT
------------------------------------------------------------------------- */

void FFT3d::remap_forward_create(int &sendsize, int &recvsize)
{
  int ssize,rsize;

  // remap uses I=fast, J=mid, K=slow, b/c current permute=0

  if (remap_preflag) {
    remap_prefast = new Remap;
    remap_prefast->remap3d = new Remap3d(world);
    remap_prefast->remap3d->collective = collective_bp;
    remap_prefast->remap3d->packflag = packflag;
    remap_prefast->remap3d->
      setup(in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
            fast_ilo,fast_ihi,fast_jlo,fast_jhi,fast_klo,fast_khi,
            2,0,0,ssize,rsize);
    sendsize = MAX(sendsize,ssize);
    recvsize = MAX(recvsize,rsize);
    remap_prefast->remap3d_extra = NULL;
  }

  // if exchange = 0, remap direct from pencil to pencil
  // if exchange = 1, two remaps from pencil to brick, then brick to pencil
  // remap uses I=fast, J=mid, K=slow, b/c current permute=0

  remap_fastmid = new Remap;
  if (exchange == 0) {
    remap_fastmid->remap3d = new Remap3d(world);
    remap_fastmid->remap3d->collective = collective_pp;
    remap_fastmid->remap3d->packflag = packflag;
    remap_fastmid->remap3d->
      setup(fast_ilo,fast_ihi,fast_jlo,fast_jhi,fast_klo,fast_khi,
            mid_ilo,mid_ihi,mid_jlo,mid_jhi,mid_klo,mid_khi,
            2,1,0,ssize,rsize);
    sendsize = MAX(sendsize,ssize);
    recvsize = MAX(recvsize,rsize);
    remap_fastmid->remap3d_extra = NULL;

  } else {
    remap_fastmid->remap3d = new Remap3d(world);
    remap_fastmid->remap3d->collective = collective_bp;
    remap_fastmid->remap3d->packflag = packflag;
    remap_fastmid->remap3d->
      setup(fast_ilo,fast_ihi,fast_jlo,fast_jhi,fast_klo,fast_khi,
            brick_ilo,brick_ihi,brick_jlo,brick_jhi,brick_klo,brick_khi,
            2,0,0,ssize,rsize);
    sendsize = MAX(sendsize,ssize);
    recvsize = MAX(recvsize,rsize);
    remap_fastmid->remap3d_extra = new Remap3d(world);
    remap_fastmid->remap3d_extra->collective = collective_bp;
    remap_fastmid->remap3d_extra->packflag = packflag;
    remap_fastmid->remap3d_extra->
      setup(brick_ilo,brick_ihi,brick_jlo,brick_jhi,brick_klo,brick_khi,
            mid_ilo,mid_ihi,mid_jlo,mid_jhi,mid_klo,mid_khi,
            2,1,0,ssize,rsize);
    sendsize = MAX(sendsize,ssize);
    recvsize = MAX(recvsize,rsize);
  }

  // if exchange = 0, remap direct from pencil to pencil
  // if exchange = 1, two remaps from pencil to brick, then brick to pencil
  // remap uses J=fast, K=mid, I=slow, b/c current permute=1

  remap_midslow = new Remap;
  if (exchange == 0) {
    remap_midslow->remap3d = new Remap3d(world);
    remap_midslow->remap3d->collective = collective_pp;
    remap_midslow->remap3d->packflag = packflag;
    remap_midslow->remap3d->
      setup(mid_jlo,mid_jhi,mid_klo,mid_khi,mid_ilo,mid_ihi,
            slow_jlo,slow_jhi,slow_klo,slow_khi,slow_ilo,slow_ihi,
            2,1,0,ssize,rsize);
    sendsize = MAX(sendsize,ssize);
    recvsize = MAX(recvsize,rsize);
    remap_midslow->remap3d_extra = NULL;

  } else {
    remap_midslow->remap3d = new Remap3d(world);
    remap_midslow->remap3d->collective = collective_bp;
    remap_midslow->remap3d->packflag = packflag;
    remap_midslow->remap3d->
      setup(mid_jlo,mid_jhi,mid_klo,mid_khi,mid_ilo,mid_ihi,
            brick_jlo,brick_jhi,brick_klo,brick_khi,brick_ilo,brick_ihi,
            2,0,0,ssize,rsize);
    sendsize = MAX(sendsize,ssize);
    recvsize = MAX(recvsize,rsize);
    remap_midslow->remap3d_extra = new Remap3d(world);
    remap_midslow->remap3d_extra->collective = collective_bp;
    remap_midslow->remap3d_extra->packflag = packflag;
    remap_midslow->remap3d_extra->
      setup(brick_jlo,brick_jhi,brick_klo,brick_khi,brick_ilo,brick_ihi,
            slow_jlo,slow_jhi,slow_klo,slow_khi,slow_ilo,slow_ihi,
            2,1,0,ssize,rsize);
    sendsize = MAX(sendsize,ssize);
    recvsize = MAX(recvsize,rsize);
  }

  // remap uses K=fast, I=mid, J=slow, b/c current permute=2
  // newpermute is from current permute=2 to desired permute=user_permute

  if (remap_postflag) {
    remap_postslow = new Remap;
    int newpermute;
    if (permute == 0) newpermute = 1;
    if (permute == 1) newpermute = 2;
    if (permute == 2) newpermute = 0;
    remap_postslow->remap3d = new Remap3d(world);
    remap_postslow->remap3d->collective = collective_bp;
    remap_postslow->remap3d->packflag = packflag;
    remap_postslow->remap3d->
      setup(slow_klo,slow_khi,slow_ilo,slow_ihi,slow_jlo,slow_jhi,
            out_klo,out_khi,out_ilo,out_ihi,out_jlo,out_jhi,
            2,newpermute,0,ssize,rsize);
    sendsize = MAX(sendsize,ssize);
    recvsize = MAX(recvsize,rsize);
    remap_postslow->remap3d_extra = NULL;
  }

  remap_premid = new Remap;
  remap_premid->remap3d = new Remap3d(world);
  remap_premid->remap3d->collective = collective_bp;
  remap_premid->remap3d->packflag = packflag;
  remap_premid->remap3d->
    setup(in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
          mid_ilo,mid_ihi,mid_jlo,mid_jhi,mid_klo,mid_khi,
          2,1,0,ssize,rsize);
  sendsize = MAX(sendsize,ssize);
  recvsize = MAX(recvsize,rsize);
  remap_premid->remap3d_extra = NULL;

  remap_fastslow = new Remap;
  remap_fastslow->remap3d = new Remap3d(world);
  remap_fastslow->remap3d->collective = collective_pp;
  remap_fastslow->remap3d->packflag = packflag;
  remap_fastslow->remap3d->
    setup(fast_ilo,fast_ihi,fast_jlo,fast_jhi,fast_klo,fast_khi,
          slow_ilo,slow_ihi,slow_jlo,slow_jhi,slow_klo,slow_khi,
          2,2,0,ssize,rsize);
  sendsize = MAX(sendsize,ssize);
  recvsize = MAX(recvsize,rsize);
  remap_fastslow->remap3d_extra = NULL;
}

/* ----------------------------------------------------------------------
   create inverted Remap3d instances for inverse FFT if needed
------------------------------------------------------------------------- */

void FFT3d::remap_inverse_create(int &sendsize, int &recvsize)
{
  int ssize,rsize;

  // if current permute=0. remap uses I=fast, J=mid, K=slow
  // if current permute=1, remap uses J=fast, K=mid, I=slow
  // if current permute=2, remap uses K=fast, I=mid, J=slow

  if (remap_postflag) {
    remap_preslow = new Remap();
    if (permute == 0) {
      remap_preslow->remap3d = new Remap3d(world);
      remap_preslow->remap3d->collective = collective_bp;
      remap_preslow->remap3d->packflag = packflag;
      remap_preslow->remap3d->
        setup(out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,
              slow_ilo,slow_ihi,slow_jlo,slow_jhi,slow_klo,slow_khi,
              2,2,0,ssize,rsize);
      sendsize = MAX(sendsize,ssize);
      recvsize = MAX(recvsize,rsize);
    } else if (permute == 1) {
      remap_preslow->remap3d = new Remap3d(world);
      remap_preslow->remap3d->collective = collective_bp;
      remap_preslow->remap3d->packflag = packflag;
      remap_preslow->remap3d->
        setup(out_jlo,out_jhi,out_klo,out_khi,out_ilo,out_ihi,
              slow_jlo,slow_jhi,slow_klo,slow_khi,slow_ilo,slow_ihi,
              2,1,0,ssize,rsize);
      sendsize = MAX(sendsize,ssize);
      recvsize = MAX(recvsize,rsize);
    } else if (permute == 2) {
      remap_preslow->remap3d = new Remap3d(world);
      remap_preslow->remap3d->collective = collective_bp;
      remap_preslow->remap3d->packflag = packflag;
      remap_preslow->remap3d->
        setup(out_klo,out_khi,out_ilo,out_ihi,out_jlo,out_jhi,
              slow_klo,slow_khi,slow_ilo,slow_ihi,slow_jlo,slow_jhi,
              2,0,0,ssize,rsize);
      sendsize = MAX(sendsize,ssize);
      recvsize = MAX(recvsize,rsize);
    }
    remap_preslow->remap3d_extra = NULL;
  }

  // if exchange = 0, remap direct from pencil to pencil
  // if exchange = 1, two remaps from pencil to brick, then brick to pencil
  // remap uses K=fast, I=mid, J=slow, b/c current permute=2

  remap_slowmid = new Remap;
  if (exchange == 0) {
    remap_slowmid->remap3d = new Remap3d(world);
    remap_slowmid->remap3d->collective = collective_pp;
    remap_slowmid->remap3d->packflag = packflag;
    remap_slowmid->remap3d->
      setup(slow_klo,slow_khi,slow_ilo,slow_ihi,slow_jlo,slow_jhi,
            mid_klo,mid_khi,mid_ilo,mid_ihi,mid_jlo,mid_jhi,
            2,2,0,ssize,rsize);
    sendsize = MAX(sendsize,ssize);
    recvsize = MAX(recvsize,rsize);
    remap_slowmid->remap3d_extra = NULL;
  } else {
    remap_slowmid->remap3d = new Remap3d(world);
    remap_slowmid->remap3d->collective = collective_bp;
    remap_slowmid->remap3d->packflag = packflag;
    remap_slowmid->remap3d->
      setup(slow_klo,slow_khi,slow_ilo,slow_ihi,slow_jlo,slow_jhi,
            brick_klo,brick_khi,brick_ilo,brick_ihi,brick_jlo,brick_jhi,
            2,0,0,ssize,rsize);
    sendsize = MAX(sendsize,ssize);
    recvsize = MAX(recvsize,rsize);
    remap_slowmid->remap3d_extra = new Remap3d(world);
    remap_slowmid->remap3d_extra->collective = collective_bp;
    remap_slowmid->remap3d_extra->packflag = packflag;
    remap_slowmid->remap3d_extra->
      setup(brick_klo,brick_khi,brick_ilo,brick_ihi,brick_jlo,brick_jhi,
            mid_klo,mid_khi,mid_ilo,mid_ihi,mid_jlo,mid_jhi,
            2,2,0,ssize,rsize);
    sendsize = MAX(sendsize,ssize);
    recvsize = MAX(recvsize,rsize);
  }

  // if exchange = 0, remap direct from pencil to pencil
  // if exchange = 1, two remaps from pencil to brick, then brick to pencil
  // remap uses J=fast, K=mid, I=slow, b/c current permute=1

  remap_midfast = new Remap;
  if (exchange == 0) {
    remap_midfast->remap3d = new Remap3d(world);
    remap_midfast->remap3d->collective = collective_pp;
    remap_midfast->remap3d->packflag = packflag;
    remap_midfast->remap3d->
      setup(mid_jlo,mid_jhi,mid_klo,mid_khi,mid_ilo,mid_ihi,
            fast_jlo,fast_jhi,fast_klo,fast_khi,fast_ilo,fast_ihi,
            2,2,0,ssize,rsize);
    sendsize = MAX(sendsize,ssize);
    recvsize = MAX(recvsize,rsize);
    remap_midfast->remap3d_extra = NULL;
  } else {
    remap_midfast->remap3d = new Remap3d(world);
    remap_midfast->remap3d->collective = collective_bp;
    remap_midfast->remap3d->packflag = packflag;
    remap_midfast->remap3d->
      setup(mid_jlo,mid_jhi,mid_klo,mid_khi,mid_ilo,mid_ihi,
            brick_jlo,brick_jhi,brick_klo,brick_khi,brick_ilo,brick_ihi,
            2,0,0,ssize,rsize);
    sendsize = MAX(sendsize,ssize);
    recvsize = MAX(recvsize,rsize);
    remap_midfast->remap3d_extra = new Remap3d(world);
    remap_midfast->remap3d_extra->collective = collective_bp;
    remap_midfast->remap3d_extra->packflag = packflag;
    remap_midfast->remap3d_extra->
      setup(brick_jlo,brick_jhi,brick_klo,brick_khi,brick_ilo,brick_ihi,
            fast_jlo,fast_jhi,fast_klo,fast_khi,fast_ilo,fast_ihi,
            2,2,0,ssize,rsize);
  }

  // remap uses I=fast, J=mid, K=slow, b/c current permute=0

  if (remap_preflag) {
    remap_postfast = new Remap;
    remap_postfast->remap3d = new Remap3d(world);
    remap_postfast->remap3d->collective = collective_bp;
    remap_postfast->remap3d->packflag = packflag;
    remap_postfast->remap3d->
      setup(fast_ilo,fast_ihi,fast_jlo,fast_jhi,fast_klo,fast_khi,
            in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
            2,0,0,ssize,rsize);
    sendsize = MAX(sendsize,ssize);
    recvsize = MAX(recvsize,rsize);
    remap_postfast->remap3d_extra = NULL;
  }

  remap_slowfast = new Remap;
  remap_slowfast->remap3d = new Remap3d(world);
  remap_slowfast->remap3d->collective = collective_pp;
  remap_slowfast->remap3d->packflag = packflag;
  remap_slowfast->remap3d->
    setup(slow_klo,slow_khi,slow_ilo,slow_ihi,slow_jlo,slow_jhi,
          fast_klo,fast_khi,fast_ilo,fast_ihi,fast_jlo,fast_jhi,
          2,1,0,ssize,rsize);
  sendsize = MAX(sendsize,ssize);
  recvsize = MAX(recvsize,rsize);
  remap_slowfast->remap3d_extra = NULL;

  remap_postmid = new Remap;
  remap_postmid->remap3d = new Remap3d(world);
  remap_postmid->remap3d->collective = collective_bp;
  remap_postmid->remap3d->packflag = packflag;
  remap_postmid->remap3d->
    setup(mid_jlo,mid_jhi,mid_klo,mid_khi,mid_ilo,mid_ihi,
          in_jlo,in_jhi,in_klo,in_khi,in_ilo,in_ihi,
          2,2,0,ssize,rsize);
  sendsize = MAX(sendsize,ssize);
  recvsize = MAX(recvsize,rsize);
  remap_postmid->remap3d_extra = NULL;
}

/* ----------------------------------------------------------------------
   tally memory used by all Remap3d instances
------------------------------------------------------------------------- */

int64_t FFT3d::remap_memory()
{
  int64_t memusage = 0;

  if (remap_prefast) {
    memusage += remap_prefast->remap3d->memusage;
    if (remap_prefast->remap3d_extra)
      memusage += remap_prefast->remap3d_extra->memusage;
  }
  if (remap_fastmid) {
    memusage += remap_fastmid->remap3d->memusage;
    if (remap_fastmid->remap3d_extra)
      memusage += remap_fastmid->remap3d_extra->memusage;
  }
  if (remap_midslow) {
    memusage += remap_midslow->remap3d->memusage;
    if (remap_midslow->remap3d_extra)
      memusage += remap_midslow->remap3d_extra->memusage;
  }
  if (remap_postslow) {
    memusage += remap_postslow->remap3d->memusage;
    if (remap_postslow->remap3d_extra)
      memusage += remap_postslow->remap3d_extra->memusage;
  }

  if (remap_preslow) {
    memusage += remap_preslow->remap3d->memusage;
    if (remap_preslow->remap3d_extra)
      memusage += remap_preslow->remap3d_extra->memusage;
  }
  if (remap_slowmid) {
    memusage += remap_slowmid->remap3d->memusage;
    if (remap_slowmid->remap3d_extra)
      memusage += remap_slowmid->remap3d_extra->memusage;
  }
  if (remap_midfast) {
    memusage += remap_midfast->remap3d->memusage;
    if (remap_midfast->remap3d_extra)
      memusage += remap_midfast->remap3d_extra->memusage;
  }
  if (remap_postfast) {
    memusage += remap_postfast->remap3d->memusage;
    if (remap_postfast->remap3d_extra)
      memusage += remap_postfast->remap3d_extra->memusage;
  }

  return memusage;
}

// -------------------------------------------------------------------
// FFTW3 FFTs
// -------------------------------------------------------------------

void FFT3d::setup_ffts()
{
  fft1d = "FFTW3";

  int n = fft_fast->n;
  fft_fast->plan_forward =
    FFTW_API(plan_many_dft)(1,&nfast,n,NULL,&nfast,1,nfast,NULL,&nfast,1,nfast,
                            FFTW_FORWARD,FFTW_ESTIMATE);
  fft_fast->plan_backward =
    FFTW_API(plan_many_dft)(1,&nfast,n,NULL,&nfast,1,nfast,NULL,&nfast,1,nfast,
                            FFTW_BACKWARD,FFTW_ESTIMATE);
  n = fft_mid->n;
  fft_mid->plan_forward =
    FFTW_API(plan_many_dft)(1,&nmid,n,NULL,&nmid,1,nmid,NULL,&nmid,1,nmid,
                            FFTW_FORWARD,FFTW_ESTIMATE);
  fft_mid->plan_backward =
    FFTW_API(plan_many_dft)(1,&nmid,n,NULL,&nmid,1,nmid,NULL,&nmid,1,nmid,
                            FFTW_BACKWARD,FFTW_ESTIMATE);
  n = fft_slow->n;
  fft_slow->plan_forward =
    FFTW_API(plan_many_dft)(1,&nslow,n,NULL,&nslow,1,nslow,NULL,&nslow,1,nslow,
                            FFTW_FORWARD,FFTW_ESTIMATE);
  fft_slow->plan_backward =
    FFTW_API(plan_many_dft)(1,&nslow,n,NULL,&nslow,1,nslow,NULL,&nslow,1,nslow,
                            FFTW_BACKWARD,FFTW_ESTIMATE);
}

void FFT3d::perform_ffts(FFT_DATA *data, int flag, FFT1d *plan)
{
  if (flag == -1) FFTW_API(execute_dft)(plan->plan_forward,data,data);
  else FFTW_API(execute_dft)(plan->plan_backward,data,data);
}

void FFT3d::scale_ffts(FFT_DATA *data)
{
  FFT_SCALAR fnorm = norm;
  FFT_SCALAR *data_ptr = (FFT_SCALAR *) data;
  for (int i = 0; i < normnum; i++) {
    *(data_ptr++) *= fnorm;
    *(data_ptr++) *= fnorm;
  }
}

void FFT3d::deallocate_ffts()
{
  FFTW_API(destroy_plan)(fft_fast->plan_forward);
  FFTW_API(destroy_plan)(fft_fast->plan_backward);
  FFTW_API(destroy_plan)(fft_mid->plan_forward);
  FFTW_API(destroy_plan)(fft_mid->plan_backward);
  FFTW_API(destroy_plan)(fft_slow->plan_forward);
  FFTW_API(destroy_plan)(fft_slow->plan_backward);
}

/* ----------------------------------------------------------------------
   check if all prime factors of N are in list of prime factors
   return 1 if yes, 0 if no
------------------------------------------------------------------------- */

int FFT3d::prime_factorable(int n)
{
  int i;

  while (n > 1) {
    for (i = 0; i < nprime; i++) {
      if (n % primes[i] == 0) {
        n /= primes[i];
        break;
      }
    }
    if (i == nprime) return 0;
  }

  return 1;
}

/* ----------------------------------------------------------------------
   computes factors of N up to std::sqrt(N)
   store ascending list in pre-allocated factors
   return nfactor
------------------------------------------------------------------------- */

void FFT3d::factor(int n)
{
  int sqroot = (int) std::sqrt(n) + 1;
  if (sqroot*sqroot > n) sqroot--;

  nfactor = 0;
  for (int i = 1; i <= sqroot; i++) {
    if (n % i) continue;
    if (nfactor < NFACTOR) factors[nfactor++] = i;
  }
}

/* ----------------------------------------------------------------------
   compute proc grid that is best match to FFT grid: Nx by Ny by Nz
   best = minimum surface area
   caller sets Nx or Ny or Nz = 1 if a 2d proc grid is desired
   else 3d if returned
   return npx,npy,npz = proc grid
   return ipx,ipy,ipz = my location in proc grid
------------------------------------------------------------------------- */

void FFT3d::procfactors(int nx, int ny, int nz,
                       int &npx, int &npy, int &npz,
                       int &ipx, int &ipy, int &ipz)
{
  int i,j,jk,ifac,jfac,kfac;
  double newarea;

  int sqroot = (int) std::sqrt(nprocs) + 1;
  if (sqroot*sqroot > nprocs) sqroot--;

  double minarea = 2.0*nx*ny + 2.0*ny*nz + 2.0*nx*nz;

  // find 3d factorization of nprocs with min surface area for (Nx,Ny,Nz) grid
  // loop over all combinations of (ifac,jfac,kfac)
  // where ifac <= jfac and jfac <= kfac
  // then do surface-area test of all 6 permutations of (ifac,jfac,kfac)

  for (i = 0; i < nfactor; i++) {
    ifac = factors[i];
    jk = nprocs/ifac;
    for (j = i; j < nfactor; j++) {
      jfac = factors[j];
      kfac = jk/jfac;
      if (ifac*jfac*kfac != nprocs) continue;
      if (ifac > jfac || jfac > kfac) continue;

      newarea = surfarea(ifac,jfac,kfac,nx,ny,nz);
      if (newarea < minarea) {
        minarea = newarea;
        npx = ifac;
        npy = jfac;
        npz = kfac;
      }

      newarea = surfarea(ifac,kfac,jfac,nx,ny,nz);
      if (newarea < minarea) {
        minarea = newarea;
        npx = ifac;
        npy = kfac;
        npz = jfac;
      }

      newarea = surfarea(jfac,ifac,kfac,nx,ny,nz);
      if (newarea < minarea) {
        minarea = newarea;
        npx = jfac;
        npy = ifac;
        npz = kfac;
      }

      newarea = surfarea(jfac,kfac,ifac,nx,ny,nz);
      if (newarea < minarea) {
        minarea = newarea;
        npx = jfac;
        npy = kfac;
        npz = ifac;
      }

      newarea = surfarea(kfac,ifac,jfac,nx,ny,nz);
      if (newarea < minarea) {
        minarea = newarea;
        npx = kfac;
        npy = ifac;
        npz = jfac;
      }

      newarea = surfarea(kfac,jfac,ifac,nx,ny,nz);
      if (newarea < minarea) {
        minarea = newarea;
        npx = kfac;
        npy = jfac;
        npz = ifac;
      }
    }
  }

  // my location in 3d proc grid

  ipx = me % npx;
  ipy = (me/npx) % npy;
  ipz = me / (npx*npy);
}

/* ----------------------------------------------------------------------
   compute per-proc surface area for I,J,K proc grid and a Nx,Ny,Nz FFT grid
   if Nx or Ny or Nz = 1, force corresponding I,J,K to be 1, else return BIG
------------------------------------------------------------------------- */

double FFT3d::surfarea(int i, int j, int k, int nx, int ny, int nz)
{
  if (nx == 1 && i != 1) return BIG;
  if (ny == 1 && j != 1) return BIG;
  if (nz == 1 && k != 1) return BIG;

  double dx = 1.0*nx/i;
  double dy = 1.0*ny/j;
  double dz = 1.0*nz/k;
  return dx*dy + dy*dz + dx*dz;
}
