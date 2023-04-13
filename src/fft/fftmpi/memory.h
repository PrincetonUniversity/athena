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

#ifndef FFT_MEMORY_H
#define FFT_MEMORY_H

#include <stdint.h>

namespace FFTMPI_NS {

class Memory {
 public:
  Memory() {}
  void *smalloc(size_t n);
  void *srealloc(void *, size_t n);
  void sfree(void *);
};

}

#endif
