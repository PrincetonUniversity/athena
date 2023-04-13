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

// 2d pack/unpack library

#ifndef FFT_PACK2D_H
#define FFT_PACK2D_H

#include <string.h>
#include "fftdata.h"

namespace FFTMPI_NS {

// loop counters for doing a pack/unpack

struct pack_plan_2d {
  int nfast;                 // # of elements in fast index
  int nslow;                 // # of elements in slow index
  int nstride;               // stride between successive slow indices
  int nqty;                  // # of values/element
};

#ifndef PACK_DATA
#define PACK_DATA FFT_SCALAR
#endif

/* ----------------------------------------------------------------------
   Pack and unpack functions:

   pack routines copy strided values from data into contiguous locs in buf
   unpack routines copy contiguous values from buf into strided locs in data
   different versions of unpack depending on permutation
     and # of values/element
   ARRAY methods work via array indices
   POINTER methods work via pointers
   MEMCPY methods work via pointers and memcpy function
------------------------------------------------------------------------- */

// ----------------------------------------------------------------------
// pack/unpack with array indices
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   pack from data -> buf
------------------------------------------------------------------------- */

static void pack_2d_array(PACK_DATA *data, PACK_DATA *buf, 
                          struct pack_plan_2d *plan)
{
  int in,out,fast,slow;
  int nfast,nslow,nstride;

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

/* ----------------------------------------------------------------------
   unpack from buf -> data
------------------------------------------------------------------------- */

static void unpack_2d_array(PACK_DATA *buf, PACK_DATA *data, 
                            struct pack_plan_2d *plan)
{
  int in,out,fast,slow;
  int nfast,nslow,nstride;

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

/* ----------------------------------------------------------------------
   unpack from buf -> data, axis permutation, 1 value/element
------------------------------------------------------------------------- */

static void unpack_2d_permute_1_array(PACK_DATA *buf, PACK_DATA *data, 
                                      struct pack_plan_2d *plan)
{
  int in,out,fast,slow;
  int nfast,nslow,nstride;

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

/* ----------------------------------------------------------------------
   unpack from buf -> data, axis permutation, 2 values/element
------------------------------------------------------------------------- */

static void unpack_2d_permute_2_array(PACK_DATA *buf, PACK_DATA *data, 
                                      struct pack_plan_2d *plan)
{
  int in,out,fast,slow;
  int nfast,nslow,nstride;

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

/* ----------------------------------------------------------------------
   unpack from buf -> data, axis permutation, nqty values/element
------------------------------------------------------------------------- */

static void unpack_2d_permute_n_array(PACK_DATA *buf, PACK_DATA *data, 
                                      struct pack_plan_2d *plan)

{
  int in,out,iqty,instart,fast,slow;
  int nfast,nslow,nstride,nqty;

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

// ----------------------------------------------------------------------
// pack/unpack with pointers
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   pack from data -> buf
------------------------------------------------------------------------- */

static void pack_2d_pointer(PACK_DATA *data, PACK_DATA *buf, 
                            struct pack_plan_2d *plan)

{
  PACK_DATA *in,*out,*begin,*end;
  int slow;
  int nfast,nslow,nstride;

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

/* ----------------------------------------------------------------------
   unpack from buf -> data
------------------------------------------------------------------------- */

static void unpack_2d_pointer(PACK_DATA *buf, PACK_DATA *data, 
                              struct pack_plan_2d *plan)

{
  PACK_DATA *in,*out,*begin,*end;
  int slow;
  int nfast,nslow,nstride;

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

/* ----------------------------------------------------------------------
   unpack from buf -> data, axis permutation, 1 value/element
------------------------------------------------------------------------- */

static void unpack_2d_permute_1_pointer(PACK_DATA *buf, PACK_DATA *data, 
                                        struct pack_plan_2d *plan)

{
  PACK_DATA *in,*out,*begin,*end;
  int slow;
  int nfast,nslow,nstride;

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

/* ----------------------------------------------------------------------
   unpack from buf -> data, axis permutation, 2 values/element
------------------------------------------------------------------------- */

static void unpack_2d_permute_2_pointer(PACK_DATA *buf, PACK_DATA *data, 
                                        struct pack_plan_2d *plan)

{
  PACK_DATA *in,*out,*begin,*end;
  int slow;
  int nfast,nslow,nstride;

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

/* ----------------------------------------------------------------------
   unpack from buf -> data, axis permutation, nqty values/element
------------------------------------------------------------------------- */

static void unpack_2d_permute_n_pointer(PACK_DATA *buf, PACK_DATA *data, 
                                        struct pack_plan_2d *plan)

{
  PACK_DATA *in,*out,*instart,*begin,*end;
  int iqty,slow;
  int nfast,nslow,nstride,nqty;

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

// ----------------------------------------------------------------------
// pack/unpack with pointers and memcpy function
// no memcpy version of unpack_permute methods, just use POINTER version
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   pack from data -> buf
------------------------------------------------------------------------- */

static void pack_2d_memcpy(PACK_DATA *data, PACK_DATA *buf, 
                           struct pack_plan_2d *plan)

{
  PACK_DATA *in,*out;
  int slow,size;
  int nfast,nslow,nstride;

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

/* ----------------------------------------------------------------------
   unpack from buf -> data
------------------------------------------------------------------------- */

static void unpack_2d_memcpy(PACK_DATA *buf, PACK_DATA *data, 
                             struct pack_plan_2d *plan)

{
  PACK_DATA *in,*out;
  int slow,size;
  int nfast,nslow,nstride;

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

/* ----------------------------------------------------------------------
   unpack from buf -> data, one axis permutation, 1 value/element
------------------------------------------------------------------------- */

static void unpack_2d_permute_1_memcpy(PACK_DATA *buf, PACK_DATA *data, 
                                       struct pack_plan_2d *plan)

{
  PACK_DATA *in,*out,*begin,*end;
  int slow;
  int nfast,nslow,nstride;

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

/* ----------------------------------------------------------------------
   unpack from buf -> data, axis permutation, 2 values/element
------------------------------------------------------------------------- */

static void unpack_2d_permute_2_memcpy(PACK_DATA *buf, PACK_DATA *data, 
                                       struct pack_plan_2d *plan)

{
  PACK_DATA *in,*out,*begin,*end;
  int slow;
  int nfast,nslow,nstride;

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

/* ----------------------------------------------------------------------
   unpack from buf -> data, axis permutation, nqty values/element
------------------------------------------------------------------------- */

static void unpack_2d_permute_n_memcpy(PACK_DATA *buf, PACK_DATA *data, 
                                       struct pack_plan_2d *plan)

{
  PACK_DATA *in,*out,*instart,*begin,*end;
  int iqty,slow;
  int nfast,nslow,nstride,nqty;

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

/* ---------------------------------------------------------------------- */

}

#endif
