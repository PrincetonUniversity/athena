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
#include <string.h>

// Athena++ headers
#include "athena_fft.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../globals.hpp"


// destructor

void AthenaFFT::MpiCleanup()
{
  if(FFT_ENABLED){
  }
}

// plan 1D fft
AthenaFFTPlan *AthenaFFT::CreatePlan(int nfast, AthenaFFTComplex *data, 
                                     enum AthenaFFTDirection dir)
{
  AthenaFFTPlan *plan;
  if(FFT_ENABLED){
    plan = new AthenaFFTPlan;
    plan->dir = dir;
    plan->dim = dim_;
    if(dir == AthenaFFTForward)
      plan->plan = fftw_plan_dft_1d(nfast, data, data, FFTW_FORWARD, FFTW_ESTIMATE);
    else
      plan->plan = fftw_plan_dft_1d(nfast, data, data, FFTW_BACKWARD, FFTW_ESTIMATE);
  }
  return plan;
}
// plan 2D fft
AthenaFFTPlan *AthenaFFT::CreatePlan(int nfast, int nslow, 
                                     AthenaFFTComplex *data, enum AthenaFFTDirection dir)
{
  AthenaFFTPlan *plan;

  if(FFT_ENABLED){
    plan = new AthenaFFTPlan;
    plan->dir = dir;
    plan->dim = dim_;
#ifdef MPI_PARALLEL
    int nbuf;
    if(dir == AthenaFFTForward){
      plan->dir = FFTW_FORWARD;
      plan->plan2d = fft_2d_create_plan(MPI_COMM_WORLD,nfast,nslow,
                                        f_in->is[0],f_in->ie[0],
                                        f_in->is[1],f_in->ie[1],
                                        f_out->is[f_in->iloc[0]],f_out->ie[f_in->iloc[0]],
                                        f_out->is[f_in->iloc[1]],f_out->ie[f_in->iloc[1]],
                                        0, permute1_, &nbuf);
    } else {
      plan->dir = FFTW_BACKWARD;
      plan->plan2d = fft_2d_create_plan(MPI_COMM_WORLD,nfast,nslow,
                                        b_in->is[0],b_in->ie[0],
                                        b_in->is[1],b_in->ie[1],
                                        b_out->is[b_in->iloc[0]],b_out->ie[b_in->iloc[0]],
                                        b_out->is[b_in->iloc[1]],b_out->ie[b_in->iloc[1]],
                                        0, permute2_, &nbuf);
    }
#else // MPI_PARALLEL
    if(dir == AthenaFFTForward)
      plan->plan = fftw_plan_dft_2d(nfast,nslow,data, data, FFTW_FORWARD, FFTW_MEASURE);
    else
      plan->plan = fftw_plan_dft_2d(nfast,nslow,data, data, FFTW_BACKWARD, FFTW_MEASURE);
#endif
  }

  return plan;
}

// plan 3D fft
AthenaFFTPlan *AthenaFFT::CreatePlan(int nfast, int nmid, int nslow,
                                     AthenaFFTComplex *data, enum AthenaFFTDirection dir)
{
  AthenaFFTPlan *plan;

  if(FFT_ENABLED){
    plan = new AthenaFFTPlan;
    plan->dir = dir;
    plan->dim = dim_;
#ifdef MPI_PARALLEL
    int nbuf;
    int ois[3], oie[3];
    if(dir == AthenaFFTForward){
      for(int l=0; l<dim_; l++){
        ois[l]=f_out->is[(l+(dim_-permute1_)) % dim_];
        oie[l]=f_out->ie[(l+(dim_-permute1_)) % dim_];
      }
      plan->dir = FFTW_FORWARD;
      plan->plan3d = fft_3d_create_plan(MPI_COMM_WORLD,nfast,nmid,nslow,
                                        f_in->is[0],f_in->ie[0],
                                        f_in->is[1],f_in->ie[1],
                                        f_in->is[2],f_in->ie[2],
                                        ois[0],oie[0],
                                        ois[1],oie[1],
                                        ois[2],oie[2],
                                        0, permute1_, &nbuf);
    } else {
      for(int l=0; l<dim_; l++){
        ois[l]=b_out->is[(l+(dim_-permute2_)) % dim_];
        oie[l]=b_out->ie[(l+(dim_-permute2_)) % dim_];
      }
      plan->dir = FFTW_BACKWARD;
      plan->plan3d = fft_3d_create_plan(MPI_COMM_WORLD,nfast,nmid,nslow,
                                        b_in->is[0],b_in->ie[0],
                                        b_in->is[1],b_in->ie[1],
                                        b_in->is[2],b_in->ie[2],
                                        ois[0],oie[0],
                                        ois[1],oie[1],
                                        ois[2],oie[2],
                                        0, permute2_, &nbuf);
    }
#else // MPI_PARALLEL
    if(dir == AthenaFFTForward){
      plan->plan = fftw_plan_dft_3d(nfast,nmid,nslow,data, data, FFTW_FORWARD, FFTW_MEASURE);
    } else {
      plan->plan = fftw_plan_dft_3d(nfast,nmid,nslow,data, data, FFTW_BACKWARD, FFTW_MEASURE);
    }
#endif
  }

  return plan;
}

void AthenaFFT::Execute(AthenaFFTPlan *plan)
{
  if(FFT_ENABLED){
#ifdef MPI_PARALLEL
    if(plan->dim == 3) fft_3d(in, out, plan->dir, plan->plan3d);
    if(plan->dim == 2) fft_2d(in, out, plan->dir, plan->plan2d);
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

void AthenaFFT::Execute(AthenaFFTPlan *plan, AthenaFFTComplex *in_data, AthenaFFTComplex *out_data)
{
  if(FFT_ENABLED){
#ifdef MPI_PARALLEL
    if(plan->dim == 3) fft_3d(in_data, out_data, plan->dir, plan->plan3d);
    if(plan->dim == 2) fft_2d(in_data, out_data, plan->dir, plan->plan2d);
#else
    fftw_execute_dft(plan->plan, in_data, out_data);
#endif
  }
}

void AthenaFFT::MpiInitialize()
{
  std::stringstream msg;
  if(FFT_ENABLED){
#ifdef OPENMP_PARALLEL
    int threads_ok;
    threads_ok = fftw_init_threads();
    if(threads_ok) fftw_plan_with_nthreads(nthreads_);
#endif
    if(pdim_ < dim_){
// To achieve best performance with 2D-pencil decomposition,
// (1) if the "long"-axis (undecomposed-axis) is not the "fast"-axis (x-axis),
//     one needs to permute the axes to make it fast by setting "permute0":
//       yz_decomp (long-axis = x) --> permute0=0 (i,j,k)
//       xz_decomp (long-axis = y) --> permute0=1 (j,k,i)
//       xy_decomp (long-axis = z) --> permute0=2 (k,i,j)
// (2) swap mid<->slow axes for input array and permute twice for forward FFT. 
//       permute1=2 & swap1=true 
//       yz_decomp (i,k,j) --> (j,i,k)
//       xz_decomp (j,i,k) --> (k,j,i)
//       xy_decomp (k,j,i) --> (i,k,j)
// (3) swap mid<->slow axes from the output of forward FFT to prepare backward FFT.
//       permute2=2 & swap2=true
//       yz_decomp (j,k,i) --> (i,j,k)
//       xz_decomp (k,i,j) --> (j,k,i)
//       xy_decomp (i,j,k) --> (k,i,j)
// (4) final outcome is the same with original input before swapping. 
//     assign it back to original Athena array with permutation
//
// swap1=swap2=true; permute1=permute2=2; permute0 depends on the decomposition

      swap1_ = true; swap2_ = true;
      permute1_ = 2; permute2_ = 2;
      {using namespace DecompositionNames;
      if(decomp_ == x_decomp){
        permute0_ = 1;
      } else if(decomp_ == y_decomp){
        permute0_ = 2;
      } else if(decomp_ == z_decomp){
        permute0_ = 0;
      } else if(decomp_ == xy_decomp){
        permute0_ = 2;
      } else if(decomp_ == yz_decomp){
        permute0_ = 0;
      } else if(decomp_ == xz_decomp){
        permute0_ = 1;
      } else {
        msg << "Something wrong with " << pdim_ << "D decomposition!" << std::endl
        << "Current MPI Configuration is " 
        << orig_idx->np[0] << " x " << orig_idx->np[1] 
        << " x " << orig_idx->np[2] << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }}
    } else {
// For 3D block decompsition, simply set indices as in original Athena Array.
// two additional remapping will be performed to prepare and recover indices.
        swap1_ = false; swap2_ = false;
        permute0_ = 0; permute1_ = 0; permute2_ = 0;
    }

// permute axes and procs & swap mid <-> slow indices to prepare forward FFT
    f_in = new AthenaFFTIndex(orig_idx);
    f_in->PermuteAxis(permute0_);
    f_in->PermuteProc(permute0_);
    if(swap1_){
      f_in->SwapAxis(0);
      f_in->SwapProc(0);
    }
    f_in->SetLocalIndex();

// set output indices of forward FFT; 
// keep global mesh size as input, 
// reverse permutation for MPI configurations to get correct indices
    f_out = new AthenaFFTIndex(f_in);
    f_out->PermuteAxis(permute1_);
    f_out->SetLocalIndex();

// prepare backward FFT; 
// now permute fast, mid, and slow axes twice
    b_in = new AthenaFFTIndex(f_out);
    if(swap2_){
      b_in->SwapAxis(0);
      b_in->SwapProc(0);
    }
    b_in->SetLocalIndex();

// set output indices of backward FFT; 
// keep global mesh size as input, 
// reverse permutation for MPI configurations to get correct indices
    b_out = new AthenaFFTIndex(b_in);
    b_out->PermuteAxis(permute2_);
    b_out->SetLocalIndex();

    if(Globals::my_rank == 0){
      std::cout << "iloc: " << orig_idx->iloc[0] << " " << orig_idx->iloc[1] << " " << orig_idx->iloc[2] << std::endl
                << "iloc: " << f_in->iloc[0] << " " << f_in->iloc[1] << " " << f_in->iloc[2] << std::endl
                << "iloc: " << f_out->iloc[0] << " " << f_out->iloc[1] << " " << f_out->iloc[2] << std::endl
                << "iloc: " << b_in->iloc[0] << " " << b_in->iloc[1] << " " << b_in->iloc[2] << std::endl
                << "iloc: " << b_out->iloc[0] << " " << b_out->iloc[1] << " " << b_out->iloc[2] << std::endl;
    }
    if(Globals::my_rank == Globals::nranks-1){
      std::cout << "f_in  ";
      f_in->PrintIndex();
      std::cout << "f_out ";
      f_out->PrintIndex();
      std::cout << "b_in  ";
      b_in->PrintIndex();
      std::cout << "b_out ";
      b_out->PrintIndex();
    }

  }
}

long int AthenaFFT::GetIndex(const int i, const int j, const int k, AthenaFFTIndex *pidx)
{
  int old_idx[3]={i,j,k};
  int new_idx[3],nx[3];
  new_idx[0]=old_idx[pidx->iloc[0]];
  new_idx[1]=old_idx[pidx->iloc[1]];
  new_idx[2]=old_idx[pidx->iloc[2]];

  return new_idx[0] + pidx->nx[0] * (new_idx[1] + pidx->nx[1] * (new_idx[2]));
}
