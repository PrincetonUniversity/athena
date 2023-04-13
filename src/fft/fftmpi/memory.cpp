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

// Memory class

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "memory.h"

// memory alignment settings

#define FFT_INTEL_NO_TBB
#if defined(__INTEL_COMPILER)
#ifndef FFT_INTEL_NO_TBB
#define FFT_USE_TBB_ALLOCATOR
#include "tbb/scalable_allocator.h"
#else
#include <malloc.h>
size_t smalloc_size(const void *p) {
    return malloc_usable_size((void*)p);
}
#endif
#else
// platform dependent malloced size return function
// https://stackoverflow.com/questions/1281686/determine-size-of-dynamically-allocated-memory-in-c
#if defined(__linux__)
#include <malloc.h>
size_t smalloc_size(const void *p) {
    return malloc_usable_size((void*)p);
}
#elif defined(__APPLE__)
#include <malloc/malloc.h>
size_t smalloc_size(const void *p) {
    return malloc_size(p);
}
#endif
#endif

#if !defined(FFT_MEMALIGN)
#define FFT_MEMALIGN 64
#endif

using namespace FFTMPI_NS;

/* ----------------------------------------------------------------------
   safe malloc
------------------------------------------------------------------------- */

void *Memory::smalloc(size_t nbytes)
{
  if (nbytes == 0) return NULL;

#if defined(FFT_MEMALIGN)
  void *ptr;

#if defined(FFT_USE_TBB_ALLOCATOR)
  ptr = scalable_aligned_malloc(nbytes, FFT_MEMALIGN);
#else
  int retval = posix_memalign(&ptr,FFT_MEMALIGN,nbytes);
  if (retval) ptr = NULL;
#endif

#else
  void *ptr = malloc(nbytes);
#endif

  return ptr;
}

/* ----------------------------------------------------------------------
   safe realloc
------------------------------------------------------------------------- */

void *Memory::srealloc(void *ptr, size_t nbytes)
{
  if (nbytes == 0) {
    sfree(ptr);
    return NULL;
  }

#if defined(FFT_USE_TBB_ALLOCATOR)
  ptr = scalable_aligned_realloc(ptr, nbytes, FFT_MEMALIGN);
#elif defined(FFT_INTEL_NO_TBB) && defined(FFT_MEMALIGN)
  ptr = realloc(ptr, nbytes);
  uintptr_t offset = ((uintptr_t)(const void *)(ptr)) % FFT_MEMALIGN;
  if (offset) {
    void *optr = ptr;
    ptr = smalloc(nbytes);
    if (nbytes < smalloc_size(optr)) memcpy(ptr,optr,nbytes);
    else memcpy(ptr,optr,smalloc_size(optr));
    free(optr);
  }
#else
  ptr = realloc(ptr,nbytes);
#endif

  return ptr;
}

/* ----------------------------------------------------------------------
   safe free
------------------------------------------------------------------------- */

void Memory::sfree(void *ptr)
{
  if (ptr == NULL) return;

  #if defined(FFT_USE_TBB_ALLOCATOR)
  scalable_aligned_free(ptr);
  #else
  free(ptr);
  #endif

  ptr = NULL;
}
