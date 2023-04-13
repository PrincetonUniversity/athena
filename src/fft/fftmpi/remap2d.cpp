// Remap2d class

#include <stdio.h>
#include <stdlib.h>

#include "remap2d.h"
#include "pack2d.h"
#include "version.h"
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

#include "memory.h"
#include "error.h"

using namespace FFTMPI_NS;

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

/* ----------------------------------------------------------------------
   data layout for 2d remaps:

   data set of Nfast x Nslow elements is owned by P procs
   each element = nqty contiguous datums
   on input, each proc owns a subsection of the elements
   on output, each proc will own a (presumably different) subsection
   my subsection must not overlap with any other proc's subsection,
     i.e. the union of all proc's input (or output) subsections must
     exactly tile the global Nfast x Nslow data set
   when called from C, all subsection indices are
     C-style from 0 to N-1 where N = Nfast or Nmid or Nslow
   when called from F77, all subsection indices are
     F77-style from 1 to N where N = Nfast or Nmid or Nslow
   a proc can own 0 elements on input or output
     by specifying hi index < lo index
   on both input and output, data is stored contiguously on a processor
     with a fast-varying and slow-varying index

   flags caller can set before setup()
     collective = 0/1 (default = 1)
       collective MPI operations for remapping data
       0 = point to point comm
       1 = MPI all2all
     packflag = 0/1/2 (default = 0)
       style of pack/unpack methods
       0 = array
       1 = pointer
       2 = memcpy
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   instantiate a 2d Remap
   user_comm = MPI communicator for the P procs which own the data
------------------------------------------------------------------------- */

Remap2d::Remap2d(MPI_Comm user_comm)
{
  world = user_comm;
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // default settings
  // user can change them before setup()

  collective = 1;
  packflag = 0;

  // Memory and Error classes

  memory = new Memory();
  error = new Error(world);

  // initialize memory allocations

  send_offset = send_size = send_proc = NULL;
  packplan = NULL;

  recv_offset = recv_size = recv_proc = recv_bufloc = NULL;
  request = NULL;
  unpackplan = NULL;

  memusage = 0;
  sendbuf = recvbuf = NULL;

  setupflag = 0;
}

/* ----------------------------------------------------------------------
   delete a 2d remap plan
------------------------------------------------------------------------- */

Remap2d::~Remap2d()
{
  delete memory;
  delete error;

  // free new MPI communicator for collective comm

  if (collective) {
    if (newcomm != MPI_COMM_NULL) MPI_Comm_free(&newcomm);
    memory->sfree(pgroup);
  }

  // free internal arrays for point-to-point comm
  // also allocated for collective comm

  memory->sfree(send_offset);
  memory->sfree(send_size);
  memory->sfree(send_proc);
  memory->sfree(packplan);

  memory->sfree(recv_offset);
  memory->sfree(recv_size);
  memory->sfree(recv_proc);
  memory->sfree(recv_bufloc);
  memory->sfree(request);
  memory->sfree(unpackplan);

  // free internal arrays for collective commm

  if (collective) {
    memory->sfree(sendcnts);
    memory->sfree(recvcnts);
    memory->sfree(senddispls);
    memory->sfree(recvdispls);
    memory->sfree(recvmap);
    memory->sfree(sendmap);
  }

  // free buffers if internal

  if (memoryflag) {
    memory->sfree(sendbuf);
    memory->sfree(recvbuf);
  }
}

/* ----------------------------------------------------------------------
   create plan for performing a 2d remap

   inputs:
   in_ilo,in_ihi        input bounds of data I own in fast index
   in_jlo,in_jhi        input bounds of data I own in slow index
   out_ilo,out_ihi      output bounds of data I own in fast index
   out_jlo,out_jhi      output bounds of data I own in slow index
   nqty                 # of datums per element
   permute              permutation in storage order of indices on output
                          0 = no permutation
                          1 = permute = slow->fast, fast->slow
   memoryflag           user provides buffer memory or system does
                          0 = caller will provide memory
                          1 = system provides memory internally

   outputs:
   sendsize = size of send buffer, caller may choose to provide it
   recvsize = size of recv buffer, caller may choose to provide it
------------------------------------------------------------------------- */

void Remap2d::setup(int in_ilo, int in_ihi, int in_jlo, int in_jhi,
                    int out_ilo, int out_ihi, int out_jlo, int out_jhi,
                    int nqty, int user_permute, int user_memoryflag,
                    int &user_sendsize, int &user_recvsize)
{
  int i,iproc,ibuf,sendsize,recvsize;
  struct extent_2d in,out,overlap;
  struct extent_2d *inarray,*outarray;

  setupflag = 1;

  permute = user_permute;
  memoryflag = user_memoryflag;

  // store parameters in local data structs

  in.ilo = in_ilo;
  in.ihi = in_ihi;
  in.isize = in.ihi - in.ilo + 1;

  in.jlo = in_jlo;
  in.jhi = in_jhi;
  in.jsize = in.jhi - in.jlo + 1;

  out.ilo = out_ilo;
  out.ihi = out_ihi;
  out.isize = out.ihi - out.ilo + 1;

  out.jlo = out_jlo;
  out.jhi = out_jhi;
  out.jsize = out.jhi - out.jlo + 1;

  // combine output extents across all procs

  inarray = (struct extent_2d *)
    memory->smalloc(nprocs*sizeof(struct extent_2d));
  if (!inarray) error->one("Could not allocate inarray");

  outarray = (struct extent_2d *)
    memory->smalloc(nprocs*sizeof(struct extent_2d));
  if (!outarray) error->one("Could not allocate outarray");

  MPI_Allgather(&out,sizeof(struct extent_2d),MPI_BYTE,
                outarray,sizeof(struct extent_2d),MPI_BYTE,world);

  // count send collides, including self

  nsend = 0;
  iproc = me;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    nsend += collide(&in,&outarray[iproc],&overlap);
  }

  // malloc space for send info

  if (nsend) {
    if (packflag == 0) pack = pack_2d_array;
    else if (packflag == 1) pack = pack_2d_pointer;
    else if (packflag == 2) pack = pack_2d_memcpy;
    send_offset = (int *) memory->smalloc(nsend*sizeof(int));
    send_size = (int *) memory->smalloc(nsend*sizeof(int));
    send_proc = (int *) memory->smalloc(nsend*sizeof(int));
    packplan = (struct pack_plan_2d *)
      memory->smalloc(nsend*sizeof(struct pack_plan_2d));
    if (!send_offset || !send_size || !send_proc || !packplan)
      error->one("Could not allocate remap send info");
  }

  // store send info, with self as last entry

  nsend = 0;
  iproc = me;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    if (collide(&in,&outarray[iproc],&overlap)) {
      send_proc[nsend] = iproc;
      send_offset[nsend] = nqty *
        ((overlap.jlo-in.jlo)*in.isize + overlap.ilo-in.ilo);
      packplan[nsend].nfast = nqty*overlap.isize;
      packplan[nsend].nslow = overlap.jsize;
      packplan[nsend].nstride = nqty*in.isize;
      packplan[nsend].nqty = nqty;
      send_size[nsend] = nqty*overlap.isize*overlap.jsize;
      nsend++;
    }
  }

  // nsend = # of sends not including self
  // for collective mode include self in nsend list

  if (nsend && send_proc[nsend-1] == me && !collective) nsend--;

  // combine input extents across all procs

  MPI_Allgather(&in,sizeof(struct extent_2d),MPI_BYTE,
                inarray,sizeof(struct extent_2d),MPI_BYTE,world);

  // count recv collides, including self

  nrecv = 0;
  iproc = me;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    nrecv += collide(&out,&inarray[iproc],&overlap);
  }

  // malloc space for recv info

  if (nrecv) {
    if (permute == 0) {
      if (packflag == 0) unpack = unpack_2d_array;
      else if (packflag == 1) unpack = unpack_2d_pointer;
      else if (packflag == 2) unpack = unpack_2d_memcpy;
    } else if (permute == 1) {
      if (nqty == 1) {
        if (packflag == 0) unpack = unpack_2d_permute_1_array;
        else if (packflag == 1) unpack = unpack_2d_permute_1_pointer;
        else if (packflag == 2) unpack = unpack_2d_permute_1_memcpy;
      } else if (nqty == 2) {
        if (packflag == 0) unpack = unpack_2d_permute_2_array;
        else if (packflag == 1) unpack = unpack_2d_permute_2_pointer;
        else if (packflag == 2) unpack = unpack_2d_permute_2_memcpy;
      } else {
        if (packflag == 0) unpack = unpack_2d_permute_n_array;
        else if (packflag == 1) unpack = unpack_2d_permute_n_pointer;
        else if (packflag == 2) unpack = unpack_2d_permute_n_memcpy;
      }
    }

    recv_offset = (int *) memory->smalloc(nrecv*sizeof(int));
    recv_size = (int *) memory->smalloc(nrecv*sizeof(int));
    recv_proc = (int *) memory->smalloc(nrecv*sizeof(int));
    recv_bufloc = (int *) memory->smalloc(nrecv*sizeof(int));
    request = (MPI_Request *) memory->smalloc(nrecv*sizeof(MPI_Request));
    unpackplan = (struct pack_plan_2d *)
      memory->smalloc(nrecv*sizeof(struct pack_plan_2d));
    if (!recv_offset || !recv_size || !recv_proc || !recv_bufloc ||
        !request || !unpackplan)
      error->one("Could not allocate remap recv info");
  }

  // store recv info, with self as last entry

  ibuf = 0;
  nrecv = 0;
  iproc = me;

  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    if (collide(&out,&inarray[iproc],&overlap)) {
      recv_proc[nrecv] = iproc;
      recv_bufloc[nrecv] = ibuf;

      if (permute == 0) {
        recv_offset[nrecv] = nqty *
          ((overlap.jlo-out.jlo)*out.isize + (overlap.ilo-out.ilo));
        unpackplan[nrecv].nfast = nqty*overlap.isize;
        unpackplan[nrecv].nslow = overlap.jsize;
        unpackplan[nrecv].nstride = nqty*out.isize;
        unpackplan[nrecv].nqty = nqty;
      }
      else if (permute == 1) {
        recv_offset[nrecv] = nqty *
          ((overlap.ilo-out.ilo)*out.jsize + (overlap.jlo-out.jlo));
        unpackplan[nrecv].nfast = overlap.isize;
        unpackplan[nrecv].nslow = overlap.jsize;
        unpackplan[nrecv].nstride = nqty*out.jsize;
        unpackplan[nrecv].nqty = nqty;
      }

      recv_size[nrecv] = nqty*overlap.isize*overlap.jsize;
      ibuf += recv_size[nrecv];
      nrecv++;
    }
  }

  // nrecv = # of recvs not including self
  // for collectives include self in nrecv list

  int nrecv_original = nrecv;
  if (nrecv && recv_proc[nrecv-1] == me && !collective) nrecv--;

  // self = 1 if send/recv data to self

  if (nrecv == nrecv_original) self = 0;
  else self = 1;

  // for point-to-point comm
  // find biggest send message (not including self) and malloc space for it
  // if requested, allocate internal scratch space for recvs,
  // only need it if I will receive any data (including self)

  if (!collective) {
    sendsize = 0;
    for (i = 0; i < nsend; i++) sendsize = MAX(sendsize,send_size[i]);
    recvsize = nqty * out.isize*out.jsize;

    if (memoryflag && sendsize) {
      sendbuf = (FFT_SCALAR *) memory->smalloc(sendsize*sizeof(FFT_SCALAR));
      if (!sendbuf) error->one("Could not allocate sendbuf array");
    }
    if (memoryflag && recvsize) {
      recvbuf = (FFT_SCALAR *) memory->smalloc(recvsize*sizeof(FFT_SCALAR));
      if (!recvbuf) error->one("Could not allocate recvbuf array");
    }
  }

  // setup for collective communication
  // pgroup = list of procs I communicate with during remap
  // ngroup = # of procs in pgroup

  if (collective) {

    // pflag = 1 if proc is in group
    // allocate pgroup as large as all procs

    int *pflag = (int *) memory->smalloc(nprocs*sizeof(int));
    for (i = 0; i < nprocs; i++) pflag[i] = 0;

    pgroup = (int *) memory->smalloc(nprocs*sizeof(int));
    ngroup = 0;

    // add procs to pgroup that I send to and recv from, including self

    for (i = 0; i < nsend; i++) {
      if (pflag[send_proc[i]]) continue;
      pflag[send_proc[i]] = 1;
      pgroup[ngroup++] = send_proc[i];
    }

    for (i = 0; i < nrecv; i++) {
      if (pflag[recv_proc[i]]) continue;
      pflag[recv_proc[i]] = 1;
      pgroup[ngroup++] = recv_proc[i];
    }

    // loop over procs in pgroup
    // collide each inarray extent with all Nprocs output extents
    // collide each outarray extent with all Nprocs input extents
    // add any new collision to pgroup
    // keep iterating until nothing is added to pgroup

    int ngroup_extra;

    int active = 1;
    while (active) {
      active = 0;
      ngroup_extra = ngroup;
      for (int i = 0; i < ngroup; i++) {
        iproc = pgroup[i];
        for (int jproc = 0; jproc < nprocs; jproc++) {
          if (pflag[jproc]) continue;
          if (collide(&inarray[iproc],&outarray[jproc],&overlap)) {
            pflag[jproc] = 1;
            pgroup[ngroup_extra++] = jproc;
            active = 1;
          }
          if (pflag[jproc]) continue;
          if (collide(&outarray[iproc],&inarray[jproc],&overlap)) {
            pflag[jproc] = 1;
            pgroup[ngroup_extra++] = jproc;
            active = 1;
          }
        }
      }
      ngroup = ngroup_extra;
    }

    // resize pgroup to final size
    // recreate sorted pgroup from pflag

    pgroup = (int *) memory->srealloc(pgroup,ngroup*sizeof(int));

    ngroup = 0;
    for (i = 0; i < nprocs; i++)
      if (pflag[i]) pgroup[ngroup++] = i;

    memory->sfree(pflag);

    // create all2all communicators for the remap
    // based on the group each proc belongs to

    MPI_Group orig_group,new_group;
    MPI_Comm_group(world,&orig_group);
    MPI_Group_incl(orig_group,ngroup,pgroup,&new_group);
    MPI_Comm_create(world,new_group,&newcomm);
    MPI_Group_free(&orig_group);
    MPI_Group_free(&new_group);

    // create send and recv buffers for AlltoAllv collective

    sendsize = 0;
    for (int i = 0; i < nsend; i++) sendsize += send_size[i];
    recvsize = 0;
    for (int i = 0; i < nrecv; i++) recvsize += recv_size[i];

    if (memoryflag && sendsize) {
      sendbuf = (FFT_SCALAR *) memory->smalloc(sendsize*sizeof(FFT_SCALAR));
      if (!sendbuf) error->one("Could not allocate sendbuf array");
    }
    if (memoryflag && recvsize) {
      recvbuf = (FFT_SCALAR *) memory->smalloc(recvsize*sizeof(FFT_SCALAR));
      if (!recvbuf) error->one("Could not allocate recvbuf array");
    }

    sendcnts = (int *) memory->smalloc(sizeof(int)*ngroup);
    senddispls = (int *) memory->smalloc(sizeof(int)*ngroup);
    sendmap = (int *) memory->smalloc(sizeof(int)*ngroup);
    recvcnts = (int *) memory->smalloc(sizeof(int)*ngroup);
    recvdispls = (int *) memory->smalloc(sizeof(int)*ngroup);
    recvmap = (int *) memory->smalloc(sizeof(int)*ngroup);

    if (!sendcnts || !senddispls || !sendmap ||
        !recvcnts || !recvdispls || !recvmap)
      if (ngroup) error->one("Could not allocate all2all args");

    // populate sendcnts and recvdispls vectors
    // order and size of proc group is different than send_proc
    // sendmap[i] = index into send info for Ith proc in pgroup

    int offset = 0;
    for (int isend = 0; isend < ngroup; isend++) {
      sendcnts[isend] = 0;
      senddispls[isend] = 0;
      sendmap[isend] = -1;
      for (int i = 0; i < nsend; i++) {
        if (send_proc[i] != pgroup[isend]) continue;
        sendcnts[isend] = send_size[i];
        senddispls[isend] = offset;
        offset += send_size[i];
        sendmap[isend] = i;
        break;
      }
    }

    // populate recvcnts and recvdispls vectors
    // order and size of proc group is different than recv_proc
    // recvmap[i] = index into recv info for Ith proc in pgroup

    offset = 0;
    for (int irecv = 0; irecv < ngroup; irecv++) {
      recvcnts[irecv] = 0;
      recvdispls[irecv] = 0;
      recvmap[irecv] = -1;
      for (int i = 0; i < nrecv; i++) {
        if (recv_proc[i] != pgroup[irecv]) continue;
        recvcnts[irecv] = recv_size[i];
        recvdispls[irecv] = offset;
        offset += recv_size[i];
        recvmap[irecv] = i;
        break;
      }
    }
  }

  // free allocated extents

  memory->sfree(inarray);
  memory->sfree(outarray);

  // return sizes for send and recv buffers

  user_sendsize = sendsize;
  user_recvsize = recvsize;

  // set memusage
  // note there was also temporary allocation of
  //   inarray,outarray = Nprocs * sizeof(struc extent_2d)

  memusage = 0;

  // allocated for both point-to-point and collective comm
  // 3 send vectors and packplan
  // 4 recv vectors, request, and unpackplan
  // send and recv bufs if caller doesn't allocate them

  memusage += 3*nsend * sizeof(int);
  memusage += nsend * sizeof(struct pack_plan_2d);

  memusage += 4*nrecv * sizeof(int);
  memusage += nrecv * sizeof(MPI_Request *);
  memusage += nrecv * sizeof(struct pack_plan_2d);

  if (memoryflag) {
    memusage += (int64_t) sendsize * sizeof(FFT_SCALAR);
    memusage += (int64_t) recvsize * sizeof(FFT_SCALAR);
  }

  // allocated only for collective commm

  if (collective) memusage += 7*ngroup * sizeof(int);
}

/* ----------------------------------------------------------------------
   perform a 2d remap

   in           starting address of input data on this proc
   out          starting address of where output data for this proc
                  will be placed (can be same as in)
   buf          extra memory required for remap
                if memoryflag=0 was used in call to setup()
                  user_sendbuf and user_recvbuf are used
                  size was returned to caller by setup()
                if memoryflag=1 was used in call to setup()
                  user_sendbuf and user_recvbuf are not used, can be NULL
------------------------------------------------------------------------- */

void Remap2d::remap(FFT_SCALAR *in, FFT_SCALAR *out,
                    FFT_SCALAR *user_sendbuf, FFT_SCALAR *user_recvbuf)
{
  int isend,irecv;

  if (!setupflag) error->all("Cannot perform remap before setup");

  if (!memoryflag) {
    sendbuf = user_sendbuf;
    recvbuf = user_recvbuf;
  }

  // point-to-point remap communication

  if (!collective) {

    // post all recvs into scratch space

    for (irecv = 0; irecv < nrecv; irecv++)
      MPI_Irecv(&recvbuf[recv_bufloc[irecv]],recv_size[irecv],
                MPI_FFT_SCALAR,recv_proc[irecv],0,
                world,&request[irecv]);

    // send all messages to other procs

    for (isend = 0; isend < nsend; isend++) {
      pack(&in[send_offset[isend]],sendbuf,&packplan[isend]);
      MPI_Send(sendbuf,send_size[isend],MPI_FFT_SCALAR,
               send_proc[isend],0,world);
    }

    // copy in -> recvbuf -> out for self data

    if (self) {
      isend = nsend;
      pack(&in[send_offset[isend]],&recvbuf[recv_bufloc[nrecv]],
           &packplan[isend]);
      unpack(&recvbuf[recv_bufloc[nrecv]],&out[recv_offset[nrecv]],
             &unpackplan[nrecv]);
    }

    // unpack all messages from mybuf -> out

    for (int i = 0; i < nrecv; i++) {
      MPI_Waitany(nrecv,request,&irecv,MPI_STATUS_IGNORE);
      unpack(&recvbuf[recv_bufloc[irecv]],&out[recv_offset[irecv]],
             &unpackplan[irecv]);
    }

  // All2Allv collective for remap communication

  } else {

    // pack the data into SendBuffer from in

    int offset = 0;
    for (int igroup = 0; igroup < ngroup; igroup++) {
      if (sendmap[igroup] >= 0) {
        isend = sendmap[igroup];
        pack(&in[send_offset[isend]],&sendbuf[offset],&packplan[isend]);
        offset += send_size[isend];
      }
    }

    // perform All2All

    if (newcomm != MPI_COMM_NULL)
      MPI_Alltoallv(sendbuf,sendcnts,senddispls,MPI_FFT_SCALAR,
                    recvbuf,recvcnts,recvdispls,MPI_FFT_SCALAR,
                    newcomm);

    // unpack the data from recvbuf into out

    offset = 0;
    for (int igroup = 0; igroup < ngroup; igroup++) {
      if (recvmap[igroup] >= 0) {
        irecv = recvmap[igroup];
        unpack(&recvbuf[offset],&out[recv_offset[irecv]],&unpackplan[irecv]);
        offset += recv_size[irecv];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   collide 2 sets of indices to determine overlap
   compare bounds of block1 with block2 to see if they overlap
   return 1 if they do and put bounds of overlapping section in overlap
   return 0 if they do not overlap
------------------------------------------------------------------------- */

int Remap2d::collide(struct extent_2d *block1, struct extent_2d *block2,
                     struct extent_2d *overlap)

{
  overlap->ilo = MAX(block1->ilo,block2->ilo);
  overlap->ihi = MIN(block1->ihi,block2->ihi);
  overlap->jlo = MAX(block1->jlo,block2->jlo);
  overlap->jhi = MIN(block1->jhi,block2->jhi);

  if (overlap->ilo > overlap->ihi || overlap->jlo > overlap->jhi) return 0;

  overlap->isize = overlap->ihi - overlap->ilo + 1;
  overlap->jsize = overlap->jhi - overlap->jlo + 1;

  return 1;
}
