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

// Remap2d class

#ifndef FFT_REMAP2D_H
#define FFT_REMAP2D_H

#include <mpi.h>
#include "fftdata.h"

namespace FFTMPI_NS {

class Remap2d {
 public:
  int collective;         // 0 = point-to-point MPI, 1 = collective all2all
  int packflag;           // 0 = array, 1 = pointer, 2 = memcpy
  int64_t memusage;       // memory usage in bytes

  Remap2d(MPI_Comm);
  ~Remap2d();
  void setup(int, int, int, int, int, int, int, int,
             int, int, int, int &, int &);
  void remap(FFT_SCALAR *, FFT_SCALAR *, FFT_SCALAR *, FFT_SCALAR *);

 private:
  MPI_Comm world;
  int me,nprocs;
  int setupflag;

  class Memory *memory;
  class Error *error;

  // details of how to perform a 2d remap

  int permute;                      // permutation setting = 0,1
  int memoryflag;                   // 0 = user-provided bufs, 1 = internal

  // point to point communication

  int nrecv;                        // # of recvs from other procs
  int nsend;                        // # of sends to other procs
  int self;                         // whether I send/recv with myself

  int *send_offset;                 // extraction loc for each send
  int *send_size;                   // size of each send message
  int *send_proc;                   // proc to send each message to
  struct pack_plan_2d *packplan;    // pack plan for each send message
  int *recv_offset;                 // insertion loc for each recv
  int *recv_size;                   // size of each recv message
  int *recv_proc;                   // proc to recv each message from
  int *recv_bufloc;                 // offset in scratch buf for each recv
  MPI_Request *request;             // MPI request for each posted recv
  struct pack_plan_2d *unpackplan;  // unpack plan for each recv message

  // collective communication

  int ngroup;                       // # of procs in my collective comm group
  int *pgroup;                      // list of ranks in my comm group
  MPI_Comm newcomm;                 // communicator for my comm group

  int *sendcnts;                    // args for MPI_All2all()
  int *senddispls;
  int *sendmap;
  int *recvcnts;
  int *recvdispls;
  int *recvmap;

  // memory for remap sends and recvs and All2all
  // either provided by caller or allocated internally

  FFT_SCALAR *sendbuf;              // send buffer
  FFT_SCALAR *recvbuf;              // recv buffer

  // which pack & unpack functions to use

  void (*pack)(FFT_SCALAR *, FFT_SCALAR *, struct pack_plan_2d *);
  void (*unpack)(FFT_SCALAR *, FFT_SCALAR *, struct pack_plan_2d *);

  // collision between 2 regions
  
  struct extent_2d {
    int ilo,ihi,isize;
    int jlo,jhi,jsize;
  };

  int collide(struct extent_2d *, struct extent_2d *, struct extent_2d *);
};

}

#endif
