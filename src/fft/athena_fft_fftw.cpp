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
#include <stdexcept>  // runtime_error

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
    plan->dim = dim;
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
    plan->dim = dim;
#ifdef MPI_PARALLEL
    if(dir == AthenaFFTForward)
      plan->plan = fftw_mpi_plan_dft_3d(nx1,nx2,nx3,data, data, 
                                  MPI_COMM_WORLD, FFTW_FORWARD,
                                  FFTW_MEASURE | FFTW_MPI_TRANSPOSED_OUT);
    else
      plan->plan = fftw_mpi_plan_dft_3d(nx1,nx2,nx3,data, data, 
                                  MPI_COMM_WORLD, FFTW_BACKWARD, 
                                  FFTW_MEASURE | FFTW_MPI_TRANSPOSED_IN);
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

void AthenaFFT::CompatabilityCheck(int verbose)
{
  if(FFT_ENABLED){
#ifdef MPI_PARALLEL
    std::stringstream msg;
    if(np1_ == 1 && np2_ == 1) decomp_=3;
    else if(np1_ == 1 && np3_ == 1) decomp_=2;
    else if(np2_ == 1 && np3_ == 1) decomp_=1;
    else {
      msg << "## FATAL ERROR in FFT domain decomposition" << std::endl
          << "MPI configuration: " << np1_ << " x " << np2_ << " x " << np3_ 
          << " is not allowed." << std::endl
          << "FFTW only supports slab decomposition along z-axis!" << std:: endl;
      throw std::runtime_error(msg.str().c_str());
    }

    AthenaFFTInt Nx[3]={gnx1_,gnx2_,gnx3_};
    if (decomp_ == 1){ // decomposed in x-dir; swap x<->y in k-space
      knx1 = gnx1_; knx2 = gnx2_/np1_; knx3 = gnx3_;
      idisp_k = jdisp; jdisp_k = idisp/nx1*knx2; kdisp_k = kdisp;
    } else if (decomp_ == 2){ // decomposed in y-dir; swap y<->x in k-space
      Nx[0]=gnx2_; Nx[1]=gnx1_;
      knx1 = gnx1_/np2_; knx2 = gnx2_; knx3 = gnx3_;
      idisp_k = jdisp/nx2*knx1; jdisp_k = idisp; kdisp_k = kdisp;
    } else if (decomp_ == 3){ // decomposed in z-dir; swap z<->y in k-space
      Nx[0]=gnx3_; Nx[2]=gnx1_;
      knx1 = gnx1_; knx2 = gnx2_/np3_; knx3 = gnx3_;
      idisp_k = idisp; jdisp_k = kdisp/nx3*knx2; kdisp_k = jdisp;
    }

    if(vergose) {
      AthenaFFTInt local_n0,local_0_start;
      AthenaFFTInt local_n1,local_1_start;
      AthenaFFTInt alloc_local=
                fftw_mpi_local_size_3d_transposed(Nx[0],Nx[1],Nx[2],MPI_COMM_WORLD,
                                                  &local_n0,&local_0_start,
                                                  &local_n1,&local_1_start);

 
      int procid;
      MPI_Comm_rank(MPI_COMM_WORLD, &procid);
 
      std::cout << "MPI rank: " << procid << std::endl
                << "Athena Mesh " << gnx1_ << "x" << gnx2_ << "x" << gnx3_
                << " Athena MeshBlock " << nx1 << "x" << nx2 << "x" << nx3 
                << " starts at " << gis_ << " " << gjs_ << " " << gks_ << std::endl
                << "FFTW-MPI decompsed axis " << decomp_ << " with " << Nx[0] 
                << " into " << local_n0 << " starts at " << local_0_start << std::endl
                << "In k-space, MeshBlock becomes " << knx1  << "x" << knx2 << "x" << knx3 
                << " starts at " << idisp_k << " " << jdisp_k << " " << kdisp_k << std::endl
                << "FFTW-MPI transposed out gives " << local_n1 << " starts at " 
                << local_1_start << std::endl;
    }
#endif
  }
}

#ifdef MPI_PARALLEL
AthenaFFTPlan *AthenaFFT::QuickCreatePlan(enum AthenaFFTDirection dir)
{
  AthenaFFTInt Nx[3]={gnx1_,gnx2_,gnx3_};
  if (decomp_ == 1){ // decomposed in x-dir; swap x<->y in k-space
  } else if (decomp_ == 2){ // decomposed in y-dir; swap y<->x in k-space
    Nx[0]=gnx2_; Nx[1]=gnx1_;
  } else if (decomp_ == 3){ // decomposed in z-dir; swap x<->z in k-space
    Nx[0]=gnx3_; Nx[2]=gnx1_;
  }
  if(gnx3_ > 1) return CreatePlan(Nx[0],Nx[1],Nx[2],work,dir);
  else if(gnx2_ > 1) return CreatePlan(Nx[0],Nx[1],work,dir);
  else  return CreatePlan(gnx1_,work,dir);
}

long int AthenaFFT::GetIndex(const int i, const int j, const int k)
{
  if (decomp_ == 1) return k + nx3 * ( j + nx2 * i );
  else if (decomp_ == 2) return k + nx3 * ( i + nx1 * j);
  else return i + nx1 * ( j + nx2 * k );
}

long int AthenaFFT::GetFreq(const int i, const int j, const int k)
{
  if (decomp_ == 1) return k + knx3 * ( i + knx1 * j);
  else if (decomp_ == 2) return k + knx3 * ( j + knx2 * i);
  else return i + knx1 * ( k + knx3 * j);
}
#endif
