/* parallel pack functions - 1998, 1999

   Steve Plimpton, MS 1111, Dept 9221, Sandia National Labs
   (505) 845-7873
   sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory of the distribution.
*/

#include "pack_2d.h"

#if !defined(PACK_POINTER) && !defined(PACK_MEMCPY)
#define PACK_ARRAY
#endif

/* ------------------------------------------------------------------- */
/* Pack and unpack functions:

   pack routines copy strided values from data into contiguous locs in buf
   unpack routines copy contiguous values from buf into strided locs in data
   different versions of unpack depending on permutation
     and # of values/element
   PACK_ARRAY routines work via array indices (default)
   PACK_POINTER routines work via pointers
   PACK_MEMCPY routines work via pointers and memcpy function
*/
/* ------------------------------------------------------------------- */

/* ------------------------------------------------------------------- */
/* pack/unpack with array indices */
/* ------------------------------------------------------------------- */

#ifdef PACK_ARRAY

/* ------------------------------------------------------------------- */
/* pack from data -> buf */

void pack_2d(double *data, double *buf, struct pack_plan_2d *plan)

{
  register int in,out,fast,slow;
  register int nfast,nslow,nstride;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;

  in = 0;
  for (slow = 0; slow < nslow; slow++) {
    out = slow*nstride;
    for (fast = 0; fast < nfast; fast++)
      buf[in++] = data[out++];
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data */

void unpack_2d(double *buf, double *data, struct pack_plan_2d *plan)

{
  register int in,out,fast,slow;
  register int nfast,nslow,nstride;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;

  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    in = slow*nstride;
    for (fast = 0; fast < nfast; fast++)
      data[in++] = buf[out++];
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data, axis permutation, 1 value/element */

void unpack_2d_permute_1(double *buf, double *data, struct pack_plan_2d *plan)

{
  register int in,out,fast,slow;
  register int nfast,nslow,nstride;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;

  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    in = slow;
    for (fast = 0; fast < nfast; fast++, in += nstride)
      data[in] = buf[out++];
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data, axis permutation, 2 values/element */

void unpack_2d_permute_2(double *buf, double *data, struct pack_plan_2d *plan)

{
  register int in,out,fast,slow;
  register int nfast,nslow,nstride;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;

  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    in = 2*slow;
    for (fast = 0; fast < nfast; fast++, in += nstride) {
      data[in] = buf[out++];
      data[in+1] = buf[out++];
    }
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data, axis permutation, nqty values/element */

void unpack_2d_permute_n(double *buf, double *data, struct pack_plan_2d *plan)

{
  register int in,out,iqty,instart,fast,slow;
  register int nfast,nslow,nstride,nqty;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;
  nqty = plan->nqty;

  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    instart = nqty*slow;
    for (fast = 0; fast < nfast; fast++, instart += nstride) {
      in = instart;
      for (iqty = 0; iqty < nqty; iqty++) data[in++] = buf[out++];
    }
  }
}

#endif

/* ------------------------------------------------------------------- */
/* pack/unpack with pointers */
/* ------------------------------------------------------------------- */

#ifdef PACK_POINTER

/* ------------------------------------------------------------------- */
/* pack from data -> buf */

void pack_2d(double *data, double *buf, struct pack_plan_2d *plan)

{
  register double *in,*out,*begin,*end;
  register int slow;
  register int nfast,nslow,nstride;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;

  in = buf;
  for (slow = 0; slow < nslow; slow++) {
    begin = &(data[slow*nstride]);
    end = begin + nfast;
    for (out = begin; out < end; out++)
      *(in++) = *out;
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data */

void unpack_2d(double *buf, double *data, struct pack_plan_2d *plan)

{
  register double *in,*out,*begin,*end;
  register int slow;
  register int nfast,nslow,nstride;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    begin = &(data[slow*nstride]);
    end = begin + nfast;
    for (in = begin; in < end; in++)
      *in = *(out++);
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data, axis permutation, 1 value/element */

void unpack_2d_permute_1(double *buf, double *data, struct pack_plan_2d *plan)

{
  register double *in,*out,*begin,*end;
  register int slow;
  register int nfast,nslow,nstride;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    begin = &(data[slow]);
    end = begin + nfast*nstride;
    for (in = begin; in < end; in += nstride)
      *in = *(out++);
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data, axis permutation, 2 values/element */

void unpack_2d_permute_2(double *buf, double *data, struct pack_plan_2d *plan)

{
  register double *in,*out,*begin,*end;
  register int slow;
  register int nfast,nslow,nstride;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    begin = &(data[2*slow]);
    end = begin + nfast*nstride;
    for (in = begin; in < end; in += nstride) {
      *in = *(out++);
      *(in+1) = *(out++);
    }
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data, axis permutation, nqty values/element */

void unpack_2d_permute_n(double *buf, double *data, struct pack_plan_2d *plan)

{
  register double *in,*out,*instart,*begin,*end;
  register int iqty,slow;
  register int nfast,nslow,nstride,nqty;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;
  nqty = plan->nqty;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    begin = &(data[nqty*slow]);
    end = begin + nfast*nstride;
    for (instart = begin; instart < end; instart += nstride) {
      in = instart;
      for (iqty = 0; iqty < nqty; iqty++) *(in++) = *(out++);
    }
  }
}

#endif

/* ------------------------------------------------------------------- */
/* pack/unpack with pointers and memcpy function
   no memcpy version of unpack_permute routines,
     just use PACK_POINTER versions */
/* ------------------------------------------------------------------- */

#ifdef PACK_MEMCPY

/* ------------------------------------------------------------------- */
/* pack from data -> buf */

void pack_2d(double *data, double *buf, struct pack_plan_2d *plan)

{
  register double *in,*out;
  register int slow,size;
  register int nfast,nslow,nstride;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;

  size = nfast*sizeof(double);
  for (slow = 0; slow < nslow; slow++) {
    in = &(buf[slow*nfast]);
    out = &(data[slow*nstride]);
    memcpy(in,out,size);
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data */

void unpack_2d(double *buf, double *data, struct pack_plan_2d *plan)

{
  register double *in,*out;
  register int slow,size;
  register int nfast,nslow,nstride;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;

  size = nfast*sizeof(double);
  for (slow = 0; slow < nslow; slow++) {
    in = &(data[slow*nstride]);
    out = &(buf[slow*nfast]);
    memcpy(in,out,size);
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data, axis permutation, 1 value/element */

void unpack_2d_permute_1(double *buf, double *data, struct pack_plan_2d *plan)

{
  register double *in,*out,*begin,*end;
  register int slow;
  register int nfast,nslow,nstride;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    begin = &(data[slow]);
    end = begin + nfast*nstride;
    for (in = begin; in < end; in += nstride)
      *in = *(out++);
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data, axis permutation, 2 values/element */

void unpack_2d_permute_2(double *buf, double *data, struct pack_plan_2d *plan)

{
  register double *in,*out,*begin,*end;
  register int slow;
  register int nfast,nslow,nstride;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    begin = &(data[2*slow]);
    end = begin + nfast*nstride;
    for (in = begin; in < end; in += nstride) {
      *in = *(out++);
      *(in+1) = *(out++);
    }
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data, axis permutation, nqty values/element */

void unpack_2d_permute_n(double *buf, double *data, struct pack_plan_2d *plan)

{
  register double *in,*out,*instart,*begin,*end;
  register int iqty,slow;
  register int nfast,nslow,nstride,nqty;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;
  nqty = plan->nqty;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    begin = &(data[nqty*slow]);
    end = begin + nfast*nstride;
    for (instart = begin; instart < end; instart += nstride) {
      in = instart;
      for (iqty = 0; iqty < nqty; iqty++) *(in++) = *(out++);
    }
  }
}

#endif
