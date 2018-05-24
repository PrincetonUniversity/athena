/* parallel FFT functions - 1998, 1999

   Steve Plimpton, MS 1111, Dept 9221, Sandia National Labs
   (505) 845-7873
   sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory of the distribution.

   Modified March 12, 2007 by Nicole Lemaster
     Rewritten to work only with FFTW 3.x
*/

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include "pack_2d.h"
#include "remap_2d.h"
#include "fft_2d.h"

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ------------------------------------------------------------------- */
/* Data layout for 2d FFTs:

   data set of Nfast x Nslow elements is owned by P procs
   on input, each proc owns a subsection of the elements
   on output, each proc will own a (possibly different) subsection
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
/* Perform 2d FFT */

/* Arguments:

   in           starting address of input data on this proc
   out          starting address of where output data for this proc
                  will be placed (can be same as in)
   flag         1 for forward FFT, -1 for inverse FFT
   plan         plan returned by previous call to fft_2d_create_plan
*/

void fft_2d(FFT_DATA *in, FFT_DATA *out, int flag, struct fft_plan_2d *plan)

{
  int i,offset,num;
  double norm;
  FFT_DATA *data,*copy;

/* pre-remap to prepare for 1st FFTs if needed
   copy = loc for remap result */

  if (plan->pre_plan) {
    if (plan->pre_target == 0)
      copy = out;
    else
      copy = plan->copy;
    remap_2d((double *) in, (double *) copy, (double *) plan->scratch,
	     plan->pre_plan);
    data = copy;
  } else
    data = in;

/* 1d FFTs along fast axis */

  if (flag == -1)
    fftw_execute_dft(plan->plan_fast_forward,data,data);
  else
    fftw_execute_dft(plan->plan_fast_backward,data,data);

/* mid-remap to prepare for 2nd FFTs
   copy = loc for remap result */

  if (plan->mid_target == 0)
    copy = out;
  else
    copy = plan->copy;
  remap_2d((double *) data, (double *) copy, (double *) plan->scratch,
	   plan->mid_plan);
  data = copy;

/* 1d FFTs along slow axis */

  if (flag == -1)
    fftw_execute_dft(plan->plan_slow_forward,data,data);
  else
    fftw_execute_dft(plan->plan_slow_backward,data,data);

/* post-remap to put data in output format if needed
   destination is always out */

  if (plan->post_plan)
    remap_2d((double *) data, (double *) out, (double *) plan->scratch,
	     plan->post_plan);

/* scaling if required */

  if (flag == 1 && plan->scaled) {
    norm = plan->norm;
    num = plan->normnum;
    for (i = 0; i < num; i++) {
      out[i][0] *= norm;
      out[i][1] *= norm;
    }
  }

}

/* ------------------------------------------------------------------- */
/* Create plan for performing a 2d FFT */

/* Arguments:

   comm                 MPI communicator for the P procs which own the data
   nfast,nslow          size of global 2d matrix
   in_ilo,in_ihi        input bounds of data I own in fast index
   in_jlo,in_jhi        input bounds of data I own in slow index
   out_ilo,out_ihi      output bounds of data I own in fast index
   out_jlo,out_jhi      output bounds of data I own in slow index
   scaled               0 = no scaling of result, 1 = scaling
   permute              permutation in storage order of indices on output
                          0 = no permutation
			  1 = permute = slow->fast, fast->slow
   nbuf                 returns size of internal storage buffers used by FFT
*/

struct fft_plan_2d *fft_2d_create_plan(
       MPI_Comm comm, int nfast, int nslow,
       int in_ilo, int in_ihi, int in_jlo, int in_jhi,
       int out_ilo, int out_ihi, int out_jlo, int out_jhi,
       int scaled, int permute, int *nbuf)

{
  struct fft_plan_2d *plan;
  int me,nprocs;
  int i,num,flag,remapflag,fftflag;
  int first_ilo,first_ihi,first_jlo,first_jhi;
  int second_ilo,second_ihi,second_jlo,second_jhi;
  int out_size,first_size,second_size,copy_size,scratch_size;
  int list[50];

/* query MPI info */

  MPI_Comm_rank(comm,&me);
  MPI_Comm_size(comm,&nprocs);

/* allocate memory for plan data struct */

  plan = (struct fft_plan_2d *) malloc(sizeof(struct fft_plan_2d));
  if (plan == NULL) return NULL;

/* remap from initial distribution to layout needed for 1st set of 1d FFTs
   not needed if all procs own entire fast axis initially
   first indices = distribution after 1st set of FFTs */

  if (in_ilo == 0 && in_ihi == nfast-1)
    flag = 0;
  else
    flag = 1;

  MPI_Allreduce(&flag,&remapflag,1,MPI_INT,MPI_MAX,comm);

  if (remapflag == 0) {
    first_ilo = in_ilo;
    first_ihi = in_ihi;
    first_jlo = in_jlo;
    first_jhi = in_jhi;
    plan->pre_plan = NULL;
  } else {
    first_ilo = 0;
    first_ihi = nfast - 1;
    first_jlo = me*nslow/nprocs;
    first_jhi = (me+1)*nslow/nprocs - 1;
    plan->pre_plan =
      remap_2d_create_plan(comm,in_ilo,in_ihi,in_jlo,in_jhi,
			   first_ilo,first_ihi,first_jlo,first_jhi,
			   FFT_PRECISION,0,0,2);
    if (plan->pre_plan == NULL) return NULL;
  }

/* 1d FFTs along fast axis */

  plan->length1 = nfast;
  plan->total1 = nfast * (first_jhi-first_jlo+1);

/* remap from 1st to 2nd FFT
   if final distribution is permute=1 with all procs owning entire slow axis
     then this remapping goes directly to final distribution
   second indices = distribution after 2nd set of FFTs */

  if (permute == 1 && out_jlo == 0 && out_jhi == nslow-1)
    flag = 0;
  else
    flag = 1;

  MPI_Allreduce(&flag,&remapflag,1,MPI_INT,MPI_MAX,comm);

  if (remapflag == 0) {
    second_ilo = out_ilo;
    second_ihi = out_ihi;
    second_jlo = out_jlo;
    second_jhi = out_jhi;
  } else {
    second_ilo = me*nfast/nprocs;
    second_ihi = (me+1)*nfast/nprocs - 1;
    second_jlo = 0;
    second_jhi = nslow - 1;
  }

  plan->mid_plan =
    remap_2d_create_plan(comm,first_ilo,first_ihi,first_jlo,first_jhi,
			 second_ilo,second_ihi,second_jlo,second_jhi,
			 FFT_PRECISION,1,0,2);
  if (plan->mid_plan == NULL) return NULL;

/* 1d FFTs along slow axis */

  plan->length2 = nslow;
  plan->total2 = (second_ihi-second_ilo+1) * nslow;

/* remap from 2nd FFT to final distribution
   not needed if permute = 1 and second indices = out indices on all procs */

  if (permute == 1 &&
      out_ilo == second_ilo && out_ihi == second_ihi &&
      out_jlo == second_jlo && out_jhi == second_jhi)
    flag = 0;
  else
    flag = 1;

  MPI_Allreduce(&flag,&remapflag,1,MPI_INT,MPI_MAX,comm);

  if (remapflag == 0) {
    plan->post_plan = NULL;
  } else {
    plan->post_plan =
      remap_2d_create_plan(comm,second_jlo,second_jhi,second_ilo,second_ihi,
			   out_jlo,out_jhi,out_ilo,out_ihi,
			   FFT_PRECISION,(permute+1)%2,0,2);
    if (plan->post_plan == NULL) return NULL;
  }

/* configure plan memory pointers and allocate work space
   out_size = amount of memory given to FFT by user
   first/second_size = amount of memory needed after pre, mid remaps
   copy_size = amount needed internally for extra copy of data
   scratch_size = amount needed internally for remap scratch space
   for each remap:
     use out space for result if big enough, else require copy buffer
     accumulate largest required remap scratch space */

  out_size = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1);
  first_size = (first_ihi-first_ilo+1) * (first_jhi-first_jlo+1);
  second_size = (second_ihi-second_ilo+1) * (second_jhi-second_jlo+1);

  copy_size = 0;
  scratch_size = 0;

  if (plan->pre_plan) {
    if (first_size <= out_size)
      plan->pre_target = 0;
    else {
      plan->pre_target = 1;
      copy_size = MAX(copy_size,first_size);
    }
    scratch_size = MAX(scratch_size,first_size);
  }

  if (plan->mid_plan) {
    if (second_size <= out_size)
      plan->mid_target = 0;
    else {
      plan->mid_target = 1;
      copy_size = MAX(copy_size,second_size);
    }
    scratch_size = MAX(scratch_size,second_size);
  }

  if (plan->post_plan)
    scratch_size = MAX(scratch_size,out_size);

  *nbuf = copy_size + scratch_size;

  if (copy_size) {
    plan->copy = (FFT_DATA *) malloc(copy_size*sizeof(FFT_DATA));
    if (plan->copy == NULL) return NULL;
  } else plan->copy = NULL;

  if (scratch_size) {
    plan->scratch = (FFT_DATA *) malloc(scratch_size*sizeof(FFT_DATA));
    if (plan->scratch == NULL) return NULL;
  } else plan->scratch = NULL;

/* system specific pre-computation of 1d FFT coeffs
   and scaling normalization */

  plan->plan_fast_forward =
    fftw_plan_many_dft(1,&(plan->length1),plan->total1/plan->length1,
                       plan->scratch,NULL,1,plan->length1,plan->scratch,
                       NULL,1,plan->length1,FFTW_FORWARD,FFTW_ESTIMATE);
  plan->plan_fast_backward =
    fftw_plan_many_dft(1,&(plan->length1),plan->total1/plan->length1,
                       plan->scratch,NULL,1,plan->length1,plan->scratch,
                       NULL,1,plan->length1,FFTW_BACKWARD,FFTW_ESTIMATE);

  if (plan->length2 == plan->length1) {
    plan->plan_slow_forward = plan->plan_fast_forward;
    plan->plan_slow_backward = plan->plan_fast_backward;
  } else {
    plan->plan_slow_forward =
      fftw_plan_many_dft(1,&(plan->length2),plan->total2/plan->length2,
                         plan->scratch,NULL,1,plan->length2,plan->scratch,
                         NULL,1,plan->length2,FFTW_FORWARD,FFTW_ESTIMATE);
    plan->plan_slow_backward =
      fftw_plan_many_dft(1,&(plan->length2),plan->total2/plan->length2,
                         plan->scratch,NULL,1,plan->length2,plan->scratch,
                         NULL,1,plan->length2,FFTW_BACKWARD,FFTW_ESTIMATE);
  }

  if (scaled == 0)
    plan->scaled = 0;
  else {
    plan->scaled = 1;
    plan->norm = 1.0/(nfast*nslow);
    plan->normnum = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1);
  }

  return plan;
}

/* ------------------------------------------------------------------- */
/* Destroy a 2d fft plan */

void fft_2d_destroy_plan(struct fft_plan_2d *plan)

{
  if (plan->pre_plan) remap_2d_destroy_plan(plan->pre_plan);
  if (plan->mid_plan) remap_2d_destroy_plan(plan->mid_plan);
  if (plan->post_plan) remap_2d_destroy_plan(plan->post_plan);

  if (plan->copy) free(plan->copy);
  if (plan->scratch) free(plan->scratch);

  if (plan->plan_slow_forward != plan->plan_fast_forward) {
    fftw_destroy_plan(plan->plan_slow_forward);
    fftw_destroy_plan(plan->plan_slow_backward);
  }
  fftw_destroy_plan(plan->plan_fast_forward);
  fftw_destroy_plan(plan->plan_fast_backward);

  free(plan);
}
