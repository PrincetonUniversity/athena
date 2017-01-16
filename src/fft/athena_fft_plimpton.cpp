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

// plimpton's
#ifdef MPI_PARALLEL
#include "mpi.h"
#include "plimpton/fft_2d.h"
#include "plimpton/fft_3d.h"
#endif

// destructor

void AthenaFFT::MpiCleanup()
{
}

// plan 1D fft
AthenaFFTPlan *AthenaFFT::CreatePlan(AthenaFFTInt nx1, AthenaFFTComplex *data, 
                                     enum AthenaFFTDirection dir)
{
  AthenaFFTPlan *plan;
  if(FFT_ENABLED){
    plan = new AthenaFFTPlan;
    plan->dim = 1;
    if(dir == AthenaFFTForward){
      plan->dir = FFTW_FORWARD;
      plan->plan = fftw_plan_dft_1d(nx1, data, data, FFTW_FORWARD, FFTW_ESTIMATE);
    } else {
      plan->dir = FFTW_BACKWARD;
      plan->plan = fftw_plan_dft_1d(nx1, data, data, FFTW_BACKWARD, FFTW_ESTIMATE);
    }
  }
  return plan;
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
    int nbuf;
    plan->plan2d = fft_2d_create_plan(MPI_COMM_WORLD,nx2,nx1,
                                      gjs_,gje_,gis_,gie_,
                                      gjs_,gje_,gis_,gie_,
                                      0, 0, &nbuf);
    if(dir == AthenaFFTForward) plan->dir = FFTW_FORWARD;
    else plan->dir = FFTW_BACKWARD;
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
    int nbuf;
    plan->plan3d = fft_3d_create_plan(MPI_COMM_WORLD,nx3,nx2,nx1,
                                     gks_,gke_,gjs_,gje_,gis_,gie_,
                                     gks_,gke_,gjs_,gje_,gis_,gie_,
                                     0, 0, &nbuf);
    if(dir == AthenaFFTForward) plan->dir = FFTW_FORWARD;
    else plan->dir = FFTW_BACKWARD;
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
    AthenaFFTComplex *data = work;
#ifdef MPI_PARALLEL
    if(plan->dim == 3) fft_3d(data, data, plan->dir, plan->plan3d);
    if(plan->dim == 2) fft_2d(data, data, plan->dir, plan->plan2d);
#else
    fftw_execute(plan->plan);
#endif
  }
}

void AthenaFFT::Execute(AthenaFFTPlan *plan, AthenaFFTComplex *data)
{
  if(FFT_ENABLED){
#ifdef MPI_PARALLEL
    if(plan->dim == 3) fft_3d(data, data, plan->dir, plan->plan3d);
    if(plan->dim == 2) fft_2d(data, data, plan->dir, plan->plan2d);
#else
    fftw_execute_dft(plan->plan, data, data);
#endif
  }
}

void AthenaFFT::Initialize()
{
  if(FFT_ENABLED){
#ifdef OPENMP_PARALLEL
    int threads_ok;
    threads_ok = fftw_init_threads();
    if(threads_ok) fftw_plan_with_nthreads(nthreads_);
#endif
  }
}

