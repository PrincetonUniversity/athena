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
#include <string.h>

// Athena++ headers
#include "athena_fft.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"

// local functions
void SetIndex(AthenaFFTIndex *pidx);
void SwapIndex(AthenaFFTIndex *pidx,int axis);
void PermuteIndex(int in[], int out[], int dim, int permute);

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
                                    f_in.fs,f_in.fe,f_in.ms,f_in.me,f_in.ss,f_in.se,
                                    f_out.fs,f_out.fe,f_out.ms,f_out.me,f_out.ss,f_out.se,
                                    0, permute1, &nbuf);
    } else {
      plan->dir = FFTW_BACKWARD;
      plan->plan3d = fft_3d_create_plan(MPI_COMM_WORLD,nfast,nmid,nslow,
                                    b_in.fs,b_in.fe,b_in.ms,b_in.me,b_in.ss,b_in.se,
                                    b_out.fs,b_out.fe,b_out.ms,b_out.me,b_out.ss,b_out.se,
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
// Check decomposition
    {using namespace DecompositionNames;

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

    std::cout << pdim << " " << np1_ << "x" << np2_ << "x" << np3_ << "->" << decomp << std::endl;
   
    LogicalLocation& loc = pmy_block->loc;
    int N[]={gnx1_,gnx2_,gnx3_};
    int np[]={np1_,np2_,np3_};
    int ip[]={loc.lx1,loc.lx2,loc.lx3};

    switch(pdim){
      case 1:
        if(decomp == x_decomp) std::cout << "1D slab decompostion in x-dir." << std::endl;
        if(decomp == y_decomp) std::cout << "1D slab decompostion in y-dir." << std::endl;
        if(decomp == z_decomp) std::cout << "1D slab decompostion in z-dir." << std::endl;
        break;
      case 2:
        if(decomp == xy_decomp){
          std::cout << "2D pencil decompostion in xy-dir." << std::endl;

          PermuteIndex(N,f_in.N,dim,1);
          PermuteIndex(np,f_in.np,dim,1);
          PermuteIndex(ip,f_in.ip,dim,1);

          swap1 = false;
          swap2 = false;
          permute1 = 2;
          permute2 = 1;
        }
        if(decomp == yz_decomp){
          std::cout << "2D pencil decompostion in yz-dir. " << std::endl;
          if(loc.lx1 != 0 ) std::cout << "Something wrong..." << std::endl;

          PermuteIndex(N,f_in.N,dim,0);
          PermuteIndex(np,f_in.np,dim,0);
          PermuteIndex(ip,f_in.ip,dim,0);

          swap1 = true;
          swap2 = true;
          permute1 = 2;
          permute2 = 2;

/*
          f_in.gnfast = gnx1_; 
          f_in.gnmid = gnx2_; 
          f_in.gnslow = gnx3_;
          f_in.np1 = np2_; f_in.np2 = np3_;
          f_in.ip1 = loc.lx2; f_in.ip2 = loc.lx3;

          swap1 = false;
          swap2 = false;
          permute1 = 2;
          permute2 = 1;
*/
        }
        if(decomp == xz_decomp){
          std::cout << "2D pencil decompostion in xz-dir." << std::endl;

          PermuteIndex(N,f_in.N,dim,2);
          PermuteIndex(np,f_in.np,dim,2);
          PermuteIndex(ip,f_in.ip,dim,2);

          swap1 = false;
          swap2 = false;
          permute1 = 2;
          permute2 = 1;
        }
        if(swap1) SwapIndex(&f_in,0);
        SetIndex(&f_in);

        PermuteIndex(f_in.N,f_out.N,dim,0);
        PermuteIndex(f_in.np,f_out.np,dim,permute1);
        PermuteIndex(f_in.ip,f_out.ip,dim,permute1);
        SetIndex(&f_out);

        PermuteIndex(f_in.N,b_in.N,dim,permute1);
        PermuteIndex(f_in.np,b_in.np,dim,0);
        PermuteIndex(f_in.ip,b_in.ip,dim,0);
        if(swap2) SwapIndex(&b_in,0);
        SetIndex(&b_in);

        PermuteIndex(b_in.N,b_out.N,dim,0);
        PermuteIndex(b_in.np,b_out.np,dim,permute2);
        PermuteIndex(b_in.ip,b_out.ip,dim,permute2);
        SetIndex(&b_out);

        break;
      case 3:
        std::cout << "3D block mpi decompisiton at " 
                  << loc.lx1 << loc.lx2 << loc.lx3 << std::endl;

        PermuteIndex(N,f_in.N,dim,0);
        PermuteIndex(np,f_in.np,dim,0);
        PermuteIndex(ip,f_in.ip,dim,0);

        SetIndex(&f_in); 
        f_out = f_in;
        b_in = f_in;
        b_out = f_in;
  
        swap1 = false;
        swap2 = false;
        permute1 = 0;
        permute2 = 0;

        break;
      default:
        std::cout << "no mpi decompisiton" << std::endl;
    }

    knx1 = b_in.nfast;
    knx2 = b_in.nmid;
    knx3 = b_in.nslow;
    idisp_k = b_in.fs;
    jdisp_k = b_in.ms;
    kdisp_k = b_in.ss;
    dkx = 2*PI/(Real)b_in.N[0];
    dky = 2*PI/(Real)b_in.N[1];
    dkz = 2*PI/(Real)b_in.N[2];
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
  if(swap) return i + knx1 * ( k + knx3 * j);
  else return i + knx1 * ( j + knx2 * k);
}


void SetIndex(AthenaFFTIndex *pidx){
  pidx->nfast = pidx->N[0]/pidx->np[0]; 
  pidx->nmid  = pidx->N[1]/pidx->np[1]; 
  pidx->nslow = pidx->N[2]/pidx->np[2];
  pidx->fs = pidx->ip[0]*pidx->nfast ;
  pidx->ms = pidx->ip[1]*pidx->nmid ;
  pidx->ss = pidx->ip[2]*pidx->nslow;
  pidx->fe = pidx->fs + pidx->nfast - 1;
  pidx->me = pidx->ms + pidx->nmid  - 1;
  pidx->se = pidx->ss + pidx->nslow - 1;
}

void SwapIndex(AthenaFFTIndex *pidx, int axis){
  int tmp;
  int dim=3;
  int axis1=(axis+1) % dim,axis2=axis+2 % dim;
  tmp=pidx->N[axis1];
  pidx->N[axis1] = pidx->N[axis2];
  pidx->N[axis2] = tmp;

  tmp=pidx->np[axis1];
  pidx->np[axis1] = pidx->np[axis2];
  pidx->np[axis2] = tmp;

  tmp=pidx->ip[axis1];
  pidx->ip[axis1] = pidx->ip[axis2];
  pidx->ip[axis2] = tmp;
}

void PermuteIndex(int in[], int out[], int dim, int permute){
  for(int i=0; i<dim; i++) out[(i+permute) % dim]=in[i];
}
