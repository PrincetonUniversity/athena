//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gravity.cpp
//  \brief implementation of functions in class Field

// C/C++ headers
#include <iostream>
#include <sstream>

// Athena++ headers
#include "athena_fft.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"

void AthenaFFT::MpiCleanup()
{
  if(FFT_ENABLED) {
#ifdef MPI_PARALLEL
    fftw_mpi_cleanup();
#endif
  }
}

// plan 2D fft
AthenaFFTPlan *AthenaFFT::CreatePlan(AthenaFFTInt nx1, AthenaFFTInt nx2, 
                                     AthenaFFTComplex *data, enum AthenaFFTDirection dir)
{
  AthenaFFTPlan *plan;

  if(FFT_ENABLED){
    plan = new AthenaFFTPlan;
    plan->dir = dir;
    plan->dim = 2;
#ifdef MPI_PARALLEL
    if(dir == AthenaFFTForward)
      plan->plan = fftw_mpi_plan_dft_2d(nx1,nx2,data, data, 
                                  MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE);
    else
      plan->plan = fftw_mpi_plan_dft_2d(nx1,nx2,data, data, 
                                  MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE);
#else // MPI_PARALLEL
    if(dir == AthenaFFTForward)
      plan->plan = fftw_plan_dft_2d(nx1,nx2,data, data, FFTW_FORWARD, FFTW_MEASURE);
    else
      plan->plan = fftw_plan_dft_2d(nx1,nx2,data, data, FFTW_BACKWARD, FFTW_MEASURE);
#endif
  }

  return plan;
}

// plan 3D fft
AthenaFFTPlan *AthenaFFT::CreatePlan(AthenaFFTInt nx1, AthenaFFTInt nx2, AthenaFFTInt nx3,
                                     AthenaFFTComplex *data, enum AthenaFFTDirection dir)
{
  AthenaFFTPlan *plan;

  if(FFT_ENABLED){
    plan = new AthenaFFTPlan;
    plan->dir = dir;
    plan->dim = 3;
#ifdef MPI_PARALLEL
    if(dir == AthenaFFTForward)
      plan->plan = fftw_mpi_plan_dft_3d(nx1,nx2,nx3,data, data, 
                                  MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE);
    else
      plan->plan = fftw_mpi_plan_dft_3d(nx1,nx2,nx3,data, data, 
                                  MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE);
#else // MPI_PARALLEL
    if(dir == AthenaFFTForward){
      plan->plan = fftw_plan_dft_3d(nx1,nx2,nx3,data, data, FFTW_FORWARD, FFTW_MEASURE);
    } else {
      plan->plan = fftw_plan_dft_3d(nx1,nx2,nx3,data, data, FFTW_BACKWARD, FFTW_MEASURE);
    }
#endif
  }

  return plan;
}

void AthenaFFT::Execute(AthenaFFTPlan *plan)
{
  if(FFT_ENABLED){
    fftw_execute(plan->plan);
  }
}

void AthenaFFT::Execute(AthenaFFTPlan *plan, AthenaFFTComplex *data)
{
  if(FFT_ENABLED){
    fftw_execute_dft(plan->plan, data, data);
  }
}

void AthenaFFT::Initialize()
{
  if(FFT_ENABLED){
#ifdef MPI_PARALLEL
    fftw_mpi_init();
#endif

#ifdef OPENMP_PARALLEL
    int threads_ok;
    threads_ok = fftw_init_threads();
    if(threads_ok) fftw_plan_with_nthreads(nthreads_);
#endif
  }
}
