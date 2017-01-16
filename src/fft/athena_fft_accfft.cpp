//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gravity.cpp
//  \brief implementation of functions in class Field

// C/C++ headers
#include <iostream>
#include <sstream> //sstream
#include <stdexcept>  // runtime_error


// Athena++ headers
#include "athena_fft.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"

// destructor

void AthenaFFT::MpiCleanup()
{
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
    std::cout << "ACCFFT does not support 2D MPI FFT!" << std::endl;
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
    int gnx[3] = {nx1, nx2, nx3};
    if(dir == AthenaFFTForward) plan->dir = ACCFFT_FORWARD;
    else plan->dir = ACCFFT_BACKWARD;
    //std::cout << "gnx: " << gnx[0] << "x" << gnx[1] << "x" << gnx[2] << std::endl;
    plan->plan3d = accfft_plan_dft_3d_c2c(gnx,data,data,comm_,ACCFFT_MEASURE);
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
    if(plan->dim == 3) accfft_execute_c2c(plan->plan3d, plan->dir, data, data);
#else
    fftw_execute(plan->plan);
#endif
  }
}

void AthenaFFT::Execute(AthenaFFTPlan *plan, AthenaFFTComplex *data)
{
  if(FFT_ENABLED){
#ifdef MPI_PARALLEL
    if(plan->dim == 3) accfft_execute_c2c(plan->plan3d, plan->dir, data, data);
#else
    fftw_execute_dft(plan->plan, data, data);
#endif
  }
}

void AthenaFFT::Initialize()
{
  if(FFT_ENABLED){
#ifdef MPI_PARALLEL
    std::stringstream msg;
    int np[2];
    if(np1_ == 1) {
      np[0] = np2_; np[1] = np3_;
    } else if (np2_ == 1) {
      np[0] = np1_; np[1] = np3_;
    } else if (np3_ == 1) {
      np[0] = np1_; np[1] = np2_;
    } else {
      msg << "ACCFFT requires pencil decomposition!" << std:: endl
          << "Current MPI Configureation is " 
          << np1_ << " x " << np2_ << " x " << np3_ << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    std::cout << np1_ << " x " << np2_ << " x " << np3_ << std::endl;
    std::cout << "np: " << np[0] << "x" << np[1] << std::endl;
    accfft_create_comm(MPI_COMM_WORLD, np, &comm_);
#endif

#ifdef OPENMP_PARALLEL
    int threads_ok;
    threads_ok = fftw_init_threads();
    if(threads_ok) fftw_plan_with_nthreads(nthreads_);
#endif
  }
}

