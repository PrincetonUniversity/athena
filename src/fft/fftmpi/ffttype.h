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

#ifndef FFT_FFTTYPE_H
#define FFT_FFTTYPE_H

// FFT_PRECISION = 1 for single-precision complex (4-byte real, 4-byte imag)
// FFT_PRECISION = 2 for double-precision complex (8-byte real, 8-byte imag)

#include "fftdata.h"

#define FFT_PRECISION 2

// -------------------------------------------------------------------------
// data types for double-precision complex
#include "fftw3.h"
typedef fftw_complex FFT_DATA;
#define FFTW_API(function)  fftw_ ## function

#endif
