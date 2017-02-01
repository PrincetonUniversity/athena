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
// plan 1D fft
AthenaFFTPlan *AthenaFFT::CreatePlan(AthenaFFTInt nfast, AthenaFFTComplex *data, 
                                     enum AthenaFFTDirection dir)
{
  AthenaFFTPlan *plan;
  if(FFT_ENABLED){
    plan = new AthenaFFTPlan;
    plan->dir = dir;
    plan->dim = dim;
    if(dir == AthenaFFTForward)
      plan->plan = fftw_plan_dft_1d(nfast, data, data, FFTW_FORWARD, FFTW_ESTIMATE);
    else
      plan->plan = fftw_plan_dft_1d(nfast, data, data, FFTW_BACKWARD, FFTW_ESTIMATE);
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
    plan->dim = dim;
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
    plan->dim = dim;
#ifdef MPI_PARALLEL
    int gnx[3] = {nx1, nx2, nx3};
    if(dir == AthenaFFTForward) plan->dir = ACCFFT_FORWARD;
    else plan->dir = ACCFFT_BACKWARD;
    //std::cout << "gnx: " << gnx[0] << "x" << gnx[1] << "x" << gnx[2] << std::endl;
    plan->plan3d = accfft_plan_dft_3d_c2c(gnx,data,data,cart_2d_comm_,ACCFFT_MEASURE);
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

void AthenaFFT::Execute(AthenaFFTPlan *plan, AthenaFFTComplex *in, AthenaFFTComplex *out)
{
  if(FFT_ENABLED){
#ifdef MPI_PARALLEL
    if(plan->dim == 3) accfft_execute_c2c(plan->plan3d, plan->dir, in, out);
#else
    fftw_execute_dft(plan->plan, in, out);
#endif
  }
}

void AthenaFFT::Initialize()
{
  if(FFT_ENABLED){
#ifdef MPI_PARALLEL
    std::stringstream msg;
    int np[2];
    if(np1_ == 1) { // decomposed y and z axes
      np[0] = np2_; np[1] = np3_; decomp_ = 1;
    } else if (np2_ == 1) { // decomposed x and z axes
      np[0] = np1_; np[1] = np3_; decomp_ = 2;
    } else if (np3_ == 1) { // decomposed x and y axes
      np[0] = np1_; np[1] = np2_; decomp_ = 3;
    } else {
      msg << "ACCFFT requires pencil decomposition!" << std:: endl
          << "Current MPI Configuration is " 
          << np1_ << " x " << np2_ << " x " << np3_ << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }

    std::cout << "MPI configuration : " << np[0] << "x" << np[1] << std::endl;
    accfft_create_comm(MPI_COMM_WORLD, np, &cart_2d_comm_);
    int alloc_local,Ni[3],No[3],iis[3],ios[3];
    AthenaFFTInt Nx[3]={gnx3_,gnx2_,gnx1_};
    if (decomp_ == 1){ // decomposed in x-dir; swap x<->y in k-space
    } else if (decomp_ == 2){ // decomposed in y-dir; swap y<->x in k-space
      Nx[0]=gnx3_; Nx[1]=gnx1_; Nx[2] = gnx2_;
    } else if (decomp_ == 3){ // decomposed in z-dir; swap x<->z in k-space
      Nx[0]=gnx2_; Nx[1]=gnx1_; Nx[2] = gnx3_;
    }
    alloc_local = accfft_local_size_dft_c2c(Nx, Ni, iis, No, ios, cart_2d_comm_);
    int procid;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);

    std::cout << procid << " Athena's decomposition: " << std::endl
              << " MeshBlock: " << nx1 << " " << nx2 << " " << nx3 
              << " starts at " << idisp << " " << jdisp << " " << kdisp << std::endl
              << " ACCFFT's decomposition: " << std::endl
              << "  " << Ni[0] << "x" << Ni[1] << "x" << Ni[2] 
              << " starts at " << iis[0] << " " << iis[1] << " " << iis[2] << std::endl
              << " ACCFFT decomposition in k-space: " << std::endl
              << "  " << No[0] << "x" << No[1] << "x" << No[2] << std::endl
              << " starts at " << ios[0] << " " << ios[1] << " " << ios[2] << std::endl;

    if(decomp_ == 1){
      knx3 = No[0]; knx2 = No[1]; knx1 = No[2];
      kdisp_k = ios[0]; jdisp_k = ios[1]; idisp_k = ios[2];
    } else if (decomp_ == 2) {
      knx3 = No[0]; knx1 = No[1]; knx2 = No[2];
      kdisp_k = ios[0]; idisp_k = ios[1]; jdisp_k = ios[2];
    } else {
      knx2 = No[0]; knx1 = No[1]; knx3 = No[2];
      jdisp_k = ios[0]; idisp_k = ios[1]; kdisp_k = ios[2];
    }
#endif

#ifdef OPENMP_PARALLEL
    int threads_ok;
    threads_ok = fftw_init_threads();
    if(threads_ok) fftw_plan_with_nthreads(nthreads_);
#endif
  }
}

long int AthenaFFT::GetIndex(const int i, const int j, const int k)
{
  if (decomp_ == 1) return i + nx1 * ( j + nx2 * k );
  else if (decomp_ == 2) return j + nx2 * ( i + nx1 * k );
  else  return k + nx3 * ( i + nx1 * j );
}

long int AthenaFFT::GetFreq(const int i, const int j, const int k)
{
  if (decomp_ == 1) return i + knx1 * ( j + knx2 * k);
  else if (decomp_ == 2) return j + knx2 * ( i + knx1 * k );
  else  return k + knx3 * ( i + knx1 * j );
}

long int AthenaFFT::GetIndex(const int i, const int j, const int k, bool swap)
{
  if (decomp_ == 1) return i + nx1 * ( j + nx2 * k );
  else if (decomp_ == 2) return j + nx2 * ( i + nx1 * k );
  else  return k + nx3 * ( i + nx1 * j );
}

long int AthenaFFT::GetFreq(const int i, const int j, const int k, bool swap)
{
  if (decomp_ == 1) return i + knx1 * ( j + knx2 * k);
  else if (decomp_ == 2) return j + knx2 * ( i + knx1 * k );
  else  return k + knx3 * ( i + knx1 * j );
}
AthenaFFTPlan *AthenaFFT::QuickCreatePlan(enum AthenaFFTDirection dir)
{
  AthenaFFTInt Nx[3]={gnx3_,gnx2_,gnx1_};
  if (decomp_ == 1){ // decomposed in x-dir; swap x<->y in k-space
  } else if (decomp_ == 2){ // decomposed in y-dir; swap y<->x in k-space
    Nx[0]=gnx3_; Nx[1]=gnx1_; Nx[2] = gnx2_;
  } else if (decomp_ == 3){ // decomposed in z-dir; swap x<->z in k-space
    Nx[0]=gnx2_; Nx[1]=gnx1_; Nx[2] = gnx3_;
  }
  if(gnx3_ > 1) return CreatePlan(Nx[0],Nx[1],Nx[2],work,dir);
  else if(gnx2_ > 1) return CreatePlan(Nx[0],Nx[1],work,dir);
  else  return CreatePlan(gnx1_,work,dir);
}


