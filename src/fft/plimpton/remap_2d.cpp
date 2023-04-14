/* parallel remap functions - 1998, 1999

   Steve Plimpton, MS 1111, Dept 9221, Sandia National Labs
   (505) 845-7873
   sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory of the distribution.
*/

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include "pack_2d.h"
#include "remap_2d.h"
#include "assert.h"

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ------------------------------------------------------------------- */
/* Data layout for 2d remaps:

   data set of Nfast x Nslow elements is owned by P procs
   each element = nqty contiguous datums
   on input, each proc owns a subsection of the elements
   on output, each proc will own a (presumably different) subsection
   my subsection must not overlap with any other proc's subsection,
     i.e. the union of all proc's input (or output) subsections must
     exactly tile the global Nfast x Nslow data set
   when called from C, all subsection indices are
     C-style from 0 to N-1 where N = Nfast or Nslow
   when called from F77, all subsection indices are
     F77-style from 1 to N where N = Nfast or Nslow
   a proc can own 0 elements on input or output
     by specifying hi index < lo index
   on both input and output, data is stored contiguously on a processor
     with a fast-varying and slow-varying index
*/
/* ------------------------------------------------------------------- */

/* ------------------------------------------------------------------- */
/* Perform 2d remap */

/* Arguments:

   in           starting address of input data on this proc
   out          starting address of where output data for this proc
                  will be placed (can be same as in)
   buf          extra memory required for remap
                if memory=0 was used in call to remap_2d_create_plan
		  then buf must be big enough to hold output result
		  i.e. nqty * (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1)
		if memory=1 was used in call to remap_2d_create_plan
		  then buf is not used, can just be a dummy pointer
   plan         plan returned by previous call to remap_2d_create_plan
*/

void remap_2d(double *in, double *out, double *buf,
	      struct remap_plan_2d *plan)

{
  MPI_Status status;
  int i,isend,irecv;
  double *scratch;

  if (plan->memory == 0)
    scratch = buf;
  else
    scratch = plan->scratch;

/* post all recvs into scratch space */

  for (irecv = 0; irecv < plan->nrecv; irecv++)
    MPI_Irecv(&scratch[plan->recv_bufloc[irecv]],plan->recv_size[irecv],
	      MPI_DOUBLE,plan->recv_proc[irecv],0,
	      plan->comm,&plan->request[irecv]);

/* send all messages to other procs */

  for (isend = 0; isend < plan->nsend; isend++) {
    plan->pack(&in[plan->send_offset[isend]],
	       plan->sendbuf,&plan->packplan[isend]);
    MPI_Send(plan->sendbuf,plan->send_size[isend],MPI_DOUBLE,
	     plan->send_proc[isend],0,plan->comm);
  }

/* copy in -> scratch -> out for self data */

  if (plan->self) {
    isend = plan->nsend;
    irecv = plan->nrecv;
    plan->pack(&in[plan->send_offset[isend]],
	       &scratch[plan->recv_bufloc[irecv]],
	       &plan->packplan[isend]);
    plan->unpack(&scratch[plan->recv_bufloc[irecv]],
		 &out[plan->recv_offset[irecv]],&plan->unpackplan[irecv]);
  }

/* unpack all messages from scratch -> out */

  for (i = 0; i < plan->nrecv; i++) {
    MPI_Waitany(plan->nrecv,plan->request,&irecv,&status);
    plan->unpack(&scratch[plan->recv_bufloc[irecv]],
		 &out[plan->recv_offset[irecv]],&plan->unpackplan[irecv]);
  }
}

/* ------------------------------------------------------------------- */
/* Create plan for performing a 2d remap */

/* Arguments:

   comm                 MPI communicator for the P procs which own the data
   in_ilo,in_ihi        input bounds of data I own in fast index
   in_jlo,in_jhi        input bounds of data I own in slow index
   out_ilo,out_ihi      output bounds of data I own in fast index
   out_jlo,out_jhi      output bounds of data I own in slow index
   nqty                 # of datums per element
   permute              permutation in storage order of indices on output
                          0 = no permutation
			  1 = permute = slow->fast, fast->slow
   memory               user provides buffer memory for remap or system does
                          0 = user provides memory
			  1 = system provides memory
   precision            precision of data
                          1 = single precision (4 bytes per datum)
			  2 = double precision (8 bytes per datum)
*/

struct remap_plan_2d *remap_2d_create_plan(
       MPI_Comm comm,
       int in_ilo, int in_ihi, int in_jlo, int in_jhi,
       int out_ilo, int out_ihi, int out_jlo, int out_jhi,
       int nqty, int permute, int memory, int precision)

{
  struct remap_plan_2d *plan;
  MPI_Comm newcomm;
  struct extent_2d *array;
  struct extent_2d in,out,overlap;
  int i,iproc,nsend,nrecv,ibuf,size,me,nprocs;

/* query MPI info */

  MPI_Comm_rank(comm,&me);
  MPI_Comm_size(comm,&nprocs);

/* single precision not yet supported */

  if (precision == 1) {
    if (me == 0) printf("Single precision not supported\n");
    return NULL;
  }

/* allocate memory for plan data struct */

  plan = (struct remap_plan_2d *) malloc(sizeof(struct remap_plan_2d));
  if (plan == NULL) return NULL;

/* store parameters in local data structs */

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

/* combine output extents across all procs */

  array = (struct extent_2d *) malloc(nprocs*sizeof(struct extent_2d));
  if (array == NULL) return NULL;

  MPI_Allgather(&out,sizeof(struct extent_2d),MPI_BYTE,
		array,sizeof(struct extent_2d),MPI_BYTE,comm);

/* count send collides, including self */

  nsend = 0;
  iproc = me;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    nsend += remap_2d_collide(&in,&array[iproc],&overlap);
  }

/* malloc space for send info */

  if (nsend) {
    if (precision == 1)
      plan->pack = NULL;
    else
      plan->pack = pack_2d;

    plan->send_offset = (int *) malloc(nsend*sizeof(int));
    plan->send_size = (int *) malloc(nsend*sizeof(int));
    plan->send_proc = (int *) malloc(nsend*sizeof(int));
    plan->packplan = (struct pack_plan_2d *)
      malloc(nsend*sizeof(struct pack_plan_2d));

    if (plan->send_offset == NULL || plan->send_size == NULL ||
	plan->send_proc == NULL || plan->packplan == NULL) return NULL;
  }

/* store send info, with self as last entry */

  nsend = 0;
  iproc = me;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    if (remap_2d_collide(&in,&array[iproc],&overlap)) {
      plan->send_proc[nsend] = iproc;
      plan->send_offset[nsend] = nqty * ((overlap.jlo-in.jlo)*in.isize +
					(overlap.ilo-in.ilo));
      plan->packplan[nsend].nfast = nqty*overlap.isize;
      plan->packplan[nsend].nslow = overlap.jsize;
      plan->packplan[nsend].nstride = nqty*in.isize;
      plan->packplan[nsend].nqty = nqty;
      plan->send_size[nsend] = nqty*overlap.isize*overlap.jsize;
      nsend++;
    }
  }

/* plan->nsend = # of sends not including self */

  if (nsend && plan->send_proc[nsend-1] == me)
    plan->nsend = nsend - 1;
  else
    plan->nsend = nsend;

/* combine input extents across all procs */

  MPI_Allgather(&in,sizeof(struct extent_2d),MPI_BYTE,
		array,sizeof(struct extent_2d),MPI_BYTE,comm);

/* count recv collides, including self */

  nrecv = 0;
  iproc = me;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    nrecv += remap_2d_collide(&out,&array[iproc],&overlap);
  }

/* malloc space for recv info */

  if (nrecv) {
    if (precision == 1) {
      if (permute == 0)
	plan->unpack = NULL;
      else if (nqty == 1)
	plan->unpack = NULL;
      else if (nqty == 2)
	plan->unpack = NULL;
      else
	plan->unpack = NULL;
    } else if (precision == 2) {
      if (permute == 0)
	plan->unpack = unpack_2d;
      else if (nqty == 1)
	plan->unpack = unpack_2d_permute_1;
      else if (nqty == 2)
	plan->unpack = unpack_2d_permute_2;
      else
	plan->unpack = unpack_2d_permute_n;
    }

    plan->recv_offset = (int *) malloc(nrecv*sizeof(int));
    plan->recv_size = (int *) malloc(nrecv*sizeof(int));
    plan->recv_proc = (int *) malloc(nrecv*sizeof(int));
    plan->recv_bufloc = (int *) malloc(nrecv*sizeof(int));
    plan->request = (MPI_Request *) malloc(nrecv*sizeof(MPI_Request));
    plan->unpackplan = (struct pack_plan_2d *)
      malloc(nrecv*sizeof(struct pack_plan_2d));

    if (plan->recv_offset == NULL || plan->recv_size == NULL ||
	plan->recv_proc == NULL || plan->recv_bufloc == NULL ||
	plan->request == NULL || plan->unpackplan == NULL) return NULL;
  }

/* store recv info, with self as last entry */

  ibuf = 0;
  nrecv = 0;
  iproc = me;

  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    if (remap_2d_collide(&out,&array[iproc],&overlap)) {
      plan->recv_proc[nrecv] = iproc;
      plan->recv_bufloc[nrecv] = ibuf;

      if (permute == 0) {
	plan->recv_offset[nrecv] = nqty * ((overlap.jlo-out.jlo)*out.isize +
					  (overlap.ilo-out.ilo));
	plan->unpackplan[nrecv].nfast = nqty*overlap.isize;
	plan->unpackplan[nrecv].nslow = overlap.jsize;
	plan->unpackplan[nrecv].nstride = nqty*out.isize;
	plan->unpackplan[nrecv].nqty = nqty;
      } else {
	plan->recv_offset[nrecv] = nqty * ((overlap.ilo-out.ilo)*out.jsize +
					  (overlap.jlo-out.jlo));
	plan->unpackplan[nrecv].nfast = overlap.isize;
	plan->unpackplan[nrecv].nslow = overlap.jsize;
	plan->unpackplan[nrecv].nstride = nqty*out.jsize;
	plan->unpackplan[nrecv].nqty = nqty;
      }

      plan->recv_size[nrecv] = nqty*overlap.isize*overlap.jsize;
      ibuf += plan->recv_size[nrecv];
      nrecv++;
    }
  }

/* plan->nrecv = # of recvs not including self */

  if (nrecv && plan->recv_proc[nrecv-1] == me)
    plan->nrecv = nrecv - 1;
  else
    plan->nrecv = nrecv;

/* init remaining fields in remap plan */

  plan->memory = memory;

  if (nrecv == plan->nrecv)
    plan->self = 0;
  else
    plan->self = 1;

/* free locally malloced space */

  free(array);

/* find biggest send message (not including self) and malloc space for it */

  plan->sendbuf = NULL;

  size = 0;
  for (nsend = 0; nsend < plan->nsend; nsend++)
    size = MAX(size,plan->send_size[nsend]);

  if (size) {
    if (precision == 1)
      plan->sendbuf = NULL;
    else
      plan->sendbuf = (double *) malloc(size*sizeof(double));
    if (plan->sendbuf == NULL) return NULL;
  }

/* if requested, allocate internal scratch space for recvs,
   only need it if I will receive any data (including self) */

  plan->scratch = NULL;

  if (memory == 1) {
    if (nrecv > 0) {
      if (precision == 1)
	plan->scratch = NULL;
      else
	plan->scratch =
	  (double *) malloc(nqty*out.isize*out.jsize*sizeof(double));
      if (plan->scratch == NULL) return NULL;
    }
  }

/* create new MPI communicator for remap */

  MPI_Comm_dup(comm,&plan->comm);

/* return pointer to plan */

  return plan;
}

/* ------------------------------------------------------------------- */
/* Destroy a 2d remap plan */

/* Arguments:

   plan         plan returned by previous call to remap_2d_create_plan
*/

void remap_2d_destroy_plan(struct remap_plan_2d *plan)

{
  /* free MPI communicator */

  MPI_Comm_free(&plan->comm);

  /* free internal arrays */

  if (plan->nsend || plan->self) {
    free(plan->send_offset);
    free(plan->send_size);
    free(plan->send_proc);
    free(plan->packplan);
    if (plan->sendbuf) free(plan->sendbuf);
  }

  if (plan->nrecv || plan->self) {
    free(plan->recv_offset);
    free(plan->recv_size);
    free(plan->recv_proc);
    free(plan->recv_bufloc);
    free(plan->request);
    free(plan->unpackplan);
    if (plan->scratch) free(plan->scratch);
  }

  /* free plan itself */

  free(plan);
}

/* ------------------------------------------------------------------- */
/* collide 2 sets of indices to determine overlap */

/* compare bounds of block1 with block2 to see if they overlap
   return 1 if they do and put bounds of overlapping section in overlap
   return 0 if they do not overlap */

int remap_2d_collide(struct extent_2d *block1, struct extent_2d *block2,
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
