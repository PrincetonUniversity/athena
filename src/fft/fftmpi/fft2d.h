// FFT2d class

#ifndef FFT_FFT2D_H
#define FFT_FFT2D_H

#include <mpi.h>
#include "ffttype.h"

namespace FFTMPI_NS {/* ----------------------------------------------------------------------
   fftMPI - library for computing 3d/2d FFTs in parallel
   http://fftmpi.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright 2018 National Technology & Engineering Solutions of
   Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
   NTESS, the U.S. Government retains certain rights in this software.
   This software is distributed under the modified Berkeley Software
   Distribution (BSD) License.

   See the README file in the top-level fftMPI directory.
---------------------
---------------------------------------------------- */


class FFT2d {
 public:
  int scaled,remaponly;
  int permute,collective,exchange,packflag,memoryflag;
  int collective_bp,collective_pp;
  int64_t memusage;                 // memory usage in bytes

  int npfast1,npfast2;             // size of pencil decomp in fast dim
  int npslow1,npslow2;             // ditto for slow dim
  int npbrick1,npbrick2;           // size of brick decomp in 2 dims

  int ntrial;                            // # of tuning trial runs
  int npertrial;                         // # of FFTs per trial
  int cbest,ebest,pbest;                 // fastest setting for coll,exch,pack
  int cflags[10],eflags[10],pflags[10];  // same 3 settings for each trial
  double besttime;                       // fastest single 2d FFT time
  double setuptime;                      // setup() time after tuning
  double tfft[10];                       // single 2d FFT time for each trial
  double t1d[10];                        // 1d FFT time for each trial
  double tremap[10];                     // total remap time for each trial
  double tremap1[10],tremap2[10],tremap3[10];    // per-remap time for each trial

  const char *fft1d;                // name of 1d FFT lib
  const char *precision;            // precision of FFTs, "single" or "double"

  FFT2d(MPI_Comm,int);
  ~FFT2d();
  void setup(int, int, int, int, int, int, int, int, int, int,
             int, int &, int &, int &);
  void setup_memory(FFT_SCALAR *, FFT_SCALAR *);
  void compute(FFT_SCALAR *, FFT_SCALAR *, int);
  void only_1d_ffts(FFT_SCALAR *, int);
  void only_remaps(FFT_SCALAR *, FFT_SCALAR *, int);
  void only_one_remap(FFT_SCALAR *, FFT_SCALAR *, int, int);
  void tune(int, int, int, int, int, int, int, int, int, int, int,
	    int &, int &, int &, int, int, double, int);

 private:
  MPI_Comm world;
  int me,nprocs;
  int setupflag,setup_memory_flag;

  class Memory *memory;
  class Error *error;
  
  int normnum;                      // # of values to rescale
  double norm;                      // normalization factor for rescaling

  int nprime,nfactor;
  int *primes,*factors;

  int nfast,nslow;

  int in_ilo,in_ihi,in_jlo,in_jhi;
  int out_ilo,out_ihi,out_jlo,out_jhi;
  int fast_ilo,fast_ihi,fast_jlo,fast_jhi;
  int slow_ilo,slow_ihi,slow_jlo,slow_jhi;
  int brick_ilo,brick_ihi,brick_jlo,brick_jhi;

  int ipfast1,ipfast2;              // my loc in pencil decomp in fast dim
  int ipslow1,ipslow2;              // diito for slow dim
  int ipbrick1,ipbrick2;            // my loc in brick decomp in 3 dims

  int insize,outsize;
  int fastsize,slowsize,bricksize;
  int fftsize,sendsize,recvsize;

  int inout_layout_same;            // 1 if initial layout = final layout

  // Remap data structs

  struct Remap {
    class Remap2d *remap2d;
    class Remap2d *remap2d_extra;
  };

  Remap *remap_prefast,*remap_fastslow,*remap_postslow;
  Remap *remap_preslow,*remap_slowfast,*remap_postfast;
  int remap_preflag,remap_postflag;

  FFT_SCALAR *sendbuf;              // buffer for remap sends
  FFT_SCALAR *recvbuf;              // buffer for remap recvs

  struct FFT1d {
    int n,length,total;
    FFTW_API(plan) plan_forward;
    FFTW_API(plan) plan_backward;
  };

  FFT1d *fft_fast,*fft_slow;

  // private methods

  void deallocate_setup();
  void deallocate_setup_memory();

  void tune_trial(FFT_SCALAR *, int, int, int,
                  double &, double &, double &, double &, double &, double &);

  void remap(FFT_SCALAR *, FFT_SCALAR *, Remap *);
  void deallocate_remap(Remap *);
  void remap_forward_create(int &, int &);
  void remap_inverse_create(int &, int &);
  int64_t remap_memory();

  void setup_ffts();
  void perform_ffts(FFT_DATA *, int, FFT1d *);
  void scale_ffts(FFT_DATA *);
  void deallocate_ffts();

  int prime_factorable(int);
  void factor(int);
  void procfactors(int, int, int &, int &, int &, int &);
  double surfarea(int, int, int, int);
};

}

#endif

