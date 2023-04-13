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

// Error class

#ifndef FFT_ERROR_H
#define FFT_ERROR_H

#include <mpi.h>

namespace FFTMPI_NS {

class Error {
 public:
  Error(MPI_Comm);
  void all(const char *);
  void one(const char *);
  void warning(const char *);

 private:
  MPI_Comm world;
};

}

#endif
