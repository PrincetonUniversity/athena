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

// local functions
inline void SetIndex(AthenaFFTIndex *pidx);
inline void SwapIndex(int in[], int dim, int axis);
inline void PermuteIndex(int in[], int out[], int dim, int permute);

// destructor

void AthenaFFT::MpiCleanup()
{
}

AthenaFFTPlan *AthenaFFT::QuickCreatePlan(enum AthenaFFTDirection dir)
{
  AthenaFFTInt nfast,nmid,nslow;
  if(dir == AthenaFFTForward){
    nfast = f_in.N[0]; nmid = f_in.N[1]; nslow = f_in.N[2];
  } else {
    nfast = b_in.N[0]; nmid = b_in.N[1]; nslow = b_in.N[2];
  }
  if(gnx3_ > 1) return CreatePlan(nfast,nmid,nslow,work,dir);
  else if(gnx2_ > 1) return CreatePlan(nfast,nslow,work,dir);
  else  return CreatePlan(nfast,work,dir);
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
AthenaFFTPlan *AthenaFFT::CreatePlan(AthenaFFTInt nfast, AthenaFFTInt nslow, 
                                     AthenaFFTComplex *data, enum AthenaFFTDirection dir)
{
  AthenaFFTPlan *plan;

  if(FFT_ENABLED){
    plan = new AthenaFFTPlan;
    plan->dir = dir;
    plan->dim = dim;
#ifdef MPI_PARALLEL
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
AthenaFFTPlan *AthenaFFT::CreatePlan(AthenaFFTInt nfast, 
                                     AthenaFFTInt nmid, 
                                     AthenaFFTInt nslow,
                                     AthenaFFTComplex *data, enum AthenaFFTDirection dir)
{
  AthenaFFTPlan *plan;

  if(FFT_ENABLED){
    plan = new AthenaFFTPlan;
    plan->dir = dir;
    plan->dim = dim;
#ifdef MPI_PARALLEL
    int nbuf;
    if(dir == AthenaFFTForward){
      plan->dir = FFTW_FORWARD;
      plan->plan3d = fft_3d_create_plan(MPI_COMM_WORLD,nfast,nmid,nslow,
                                        f_in.is[0],f_in.ie[0],f_in.is[1],
                                        f_in.ie[1],f_in.is[2],f_in.ie[2],
                                        f_out.is[0],f_out.ie[0],f_out.is[1],
                                        f_out.ie[1],f_out.is[2],f_out.ie[2],
                                        0, permute1, &nbuf);
    } else {
      plan->dir = FFTW_BACKWARD;
      plan->plan3d = fft_3d_create_plan(MPI_COMM_WORLD,nfast,nmid,nslow,
                                        b_in.is[0],b_in.ie[0],b_in.is[1],
                                        b_in.ie[1],b_in.is[2],b_in.ie[2],
                                        b_out.is[0],b_out.ie[0],b_out.is[1],
                                        b_out.ie[1],b_out.is[2],b_out.ie[2],
                                        0, permute2, &nbuf);
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

void AthenaFFT::Execute(AthenaFFTPlan *plan, AthenaFFTComplex *in, AthenaFFTComplex *out)
{
  if(FFT_ENABLED){
#ifdef MPI_PARALLEL
    if(plan->dim == 3) fft_3d(in, out, plan->dir, plan->plan3d);
    if(plan->dim == 2) fft_2d(in, out, plan->dir, plan->plan2d);
#else
    fftw_execute_dft(plan->plan, in, out);
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
    std::stringstream msg;
// Check decomposition
    {using namespace DecompositionNames;

    decomp = 0;
    int pdim = 0;
    if(np1_ > 1){
      decomp = decomp | x_decomp;
      pdim++;
    }
    if(np2_ > 1){
      decomp = decomp | y_decomp;
      pdim++;
    }
    if(np3_ > 1){
      decomp = decomp | z_decomp;
      pdim++;
    }

    LogicalLocation& loc = pmy_block->loc;
// set default Mesh and MPI information
    int N[]={gnx1_,gnx2_,gnx3_};
    int np[]={np1_,np2_,np3_};
    int ip[]={loc.lx1,loc.lx2,loc.lx3};

    switch(pdim){
      case 1:
// For 1D-slab decomposition, it performs best if fast and slow axes are in the local grid
//   x_decomp permute twice (i,j,k) --> (k,i,j)
//   y_decomp no permutation(i,j,k) --> (i,j,k)
//   z_decomp permute once  (i,j,k) --> (j,k,i)
        if(decomp == x_decomp) std::cout << "1D slab decompostion in x-dir." << std::endl;
        if(decomp == y_decomp) std::cout << "1D slab decompostion in y-dir." << std::endl;
        if(decomp == z_decomp) std::cout << "1D slab decompostion in z-dir." << std::endl;
        break;
      case 2:
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

        swap1 = true; swap2 = true;
        permute1 = 2; permute2 = 2;
        if(decomp == xy_decomp){
          permute0 = 2;
        } else if(decomp == yz_decomp){
          permute0 = 0;
        } else if(decomp == xz_decomp){
          permute0 = 1;
        } else {
          msg << "Something wrong with 2D pencil decomposition!" << std:: endl
          << "Current MPI Configuration is " 
          << np1_ << " x " << np2_ << " x " << np3_ << std::endl;
          throw std::runtime_error(msg.str().c_str());
        }

// permute axes and procs & swap mid <-> slow indices to prepare forward FFT
        PermuteIndex(N,f_in.N,dim,permute0);
        PermuteIndex(np,f_in.np,dim,permute0);
        PermuteIndex(ip,f_in.ip,dim,permute0);
        if(swap1){
          SwapIndex(f_in.N,dim,0);
          SwapIndex(f_in.np,dim,0);
          SwapIndex(f_in.ip,dim,0);
        }
        SetIndex(&f_in);

// set output indices of forward FFT; 
// keep global mesh size as input, 
// reverse permutation for MPI configurations to get correct indices
        PermuteIndex(f_in.N,f_out.N,dim,0);
        PermuteIndex(f_in.np,f_out.np,dim,dim-permute1);
        PermuteIndex(f_in.ip,f_out.ip,dim,dim-permute1);
        SetIndex(&f_out);

// prepare backward FFT; 
// now permute fast, mid, and slow axes twice
        PermuteIndex(f_out.N,b_in.N,dim,permute1);
        PermuteIndex(f_out.np,b_in.np,dim,permute1);
        PermuteIndex(f_out.ip,b_in.ip,dim,permute1);
        if(swap2){
          SwapIndex(b_in.N,dim,0);
          SwapIndex(b_in.np,dim,0);
          SwapIndex(b_in.ip,dim,0);
        }
        SetIndex(&b_in);

// set output indices of backward FFT; 
// keep global mesh size as input, 
// reverse permutation for MPI configurations to get correct indices
        PermuteIndex(b_in.N,b_out.N,dim,0);
        PermuteIndex(b_in.np,b_out.np,dim,dim-permute2);
        PermuteIndex(b_in.ip,b_out.ip,dim,dim-permute2);
        SetIndex(&b_out);

        break;
      case 3:
// For 3D block decompsition, simply set indices as in original Athena Array.
// two additional remapping will be performed to prepare and recover indices.
  
        swap1 = false; swap2 = false;
        permute0 = 0; permute1 = 0; permute2 = 0;
        PermuteIndex(N,f_in.N,dim,permute0);
        PermuteIndex(np,f_in.np,dim,permute0);
        PermuteIndex(ip,f_in.ip,dim,permute0);

        SetIndex(&f_in); 
        f_out = f_in;
        b_in = f_in;
        b_out = f_in;

        break;
      default:
        msg << "no MPI decomposition!" << std:: endl
            << "Current MPI Configuration is " 
            << np1_ << " x " << np2_ << " x " << np3_ << std::endl;
        throw std::runtime_error(msg.str().c_str());
    }

    int idx_in[]={0,1,2}; 
    int idx1[3],idx2[3],idx3[3],idx4[3];
    PermuteIndex(idx_in,idx1,dim,permute0);
    if(swap1) SwapIndex(idx1,dim,0);
    PermuteIndex(idx1,idx2,dim,permute1);
    PermuteIndex(idx2,idx3,dim,0);
    if(swap1) SwapIndex(idx3,dim,0);
    PermuteIndex(idx3,idx4,dim,permute2);

    std::cout << idx3[0] << idx3[1] << idx3[2] << std::endl;
    int knx[3], disp[3];
    Real dk[3];
    knx[idx3[0]] = b_in.nx[0];
    knx[idx3[1]] = b_in.nx[1];
    knx[idx3[2]] = b_in.nx[2];
    disp[idx3[0]] = b_in.is[0];
    disp[idx3[1]] = b_in.is[1];
    disp[idx3[2]] = b_in.is[2];
    dk[idx3[0]] = 2*PI/(Real)b_in.N[0];
    dk[idx3[1]] = 2*PI/(Real)b_in.N[1];
    dk[idx3[2]] = 2*PI/(Real)b_in.N[2];
    knx1 = knx[0];
    knx2 = knx[1];
    knx3 = knx[2];
    idisp_k = disp[0];
    jdisp_k = disp[1];
    kdisp_k = disp[2];
    dkx = dk[0];
    dky = dk[1];
    dkz = dk[2];
    } // end of using namespace block
  }
}

long int AthenaFFT::GetIndex(const int i, const int j, const int k, bool swap)
{
  {using namespace DecompositionNames;

  if(swap){
    if(decomp == xz_decomp) // fast = j, mid = i, slow = k
      return j + nx2 * ( i + nx1 * k);
    else if (decomp == yz_decomp) // fast = i, mid = k, slow = j
      return i + nx1 * ( k + nx3 * j);
    else if (decomp == xy_decomp) // fast = k, mid = j, slow = i
      return k + nx3 * ( j + nx2 * i);
  } else {
    if(decomp == xz_decomp) // fast = j, mid = k, slow = i
      return j + nx2 * ( k + nx3 * i);
    else if (decomp == yz_decomp) // fast = i, mid = j, slow = k
      return i + nx1 * ( j + nx2 * k);
    else if (decomp == xy_decomp) // fast = k, mid = i, slow = j
      return k + nx3 * ( i + nx1 * j);
  }

  return i + nx1 * ( j + nx2 * k);

  }
}

long int AthenaFFT::GetFreq(const int i, const int j, const int k, bool swap)
{
  {using namespace DecompositionNames;

  if(swap){
    if(decomp == xz_decomp) // fast = k, mid = j, slow = i
      return k + knx3 * ( j + knx2 * i);
    else if (decomp == yz_decomp) // fast = j, mid = i, slow = k
      return j + knx2 * ( i + knx1 * k);
    else if (decomp == xy_decomp) // fast = i, mid = k, slow = j
      return i + knx1 * ( k + knx3 * j);
  } else {
    if(decomp == xz_decomp) // fast = k, mid = i, slow = j
      return k + knx3 * ( i + knx1 * j);
    else if (decomp == yz_decomp) // fast = j, mid = k, slow = i
      return j + knx2 * ( k + knx3 * i);
    else if (decomp == xy_decomp) // fast = i, mid = j, slow = k
      return i + knx1 * ( j + knx2 * k);
  }

  return i + knx1 * ( j + knx2 * k);
  }
}


inline void SetIndex(AthenaFFTIndex *pidx){
  for(int i=0; i<3; i++){
    pidx->nx[i] = pidx->N[i]/pidx->np[i]; 
    pidx->is[i] = pidx->ip[i]*pidx->nx[i]; 
    pidx->ie[i] = pidx->is[i] + pidx->nx[i] - 1; 
  }
}

inline void SwapIndex(int in[], int dim, int axis){
  int tmp;
  int axis1=(axis+1) % dim, axis2=axis+2 % dim;
  tmp=in[axis1];
  in[axis1]=in[axis2];
  in[axis2]=tmp;
}

inline void PermuteIndex(int in[], int out[], int dim, int permute){
  for(int i=0; i<dim; i++) out[i] = in[(i+permute) % dim];
}
