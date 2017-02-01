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

// destructor

void AthenaFFT::MpiCleanup()
{
}

AthenaFFTPlan *AthenaFFT::QuickCreatePlan(enum AthenaFFTDirection dir)
{
  AthenaFFTInt nfast,nmid,nslow;
  if(dir == AthenaFFTForward){
    nfast = f_in.gnfast; nmid = f_in.gnmid; nslow = f_in.gnslow;
  } else {
    nfast = b_in.gnfast; nmid = b_in.gnmid; nslow = b_in.gnslow;
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

    switch(pdim){
      case 1:
        if(decomp == x_decomp) std::cout << "1D slab decompostion in x-dir." << std::endl;
        if(decomp == y_decomp) std::cout << "1D slab decompostion in y-dir." << std::endl;
        if(decomp == z_decomp) std::cout << "1D slab decompostion in z-dir." << std::endl;
        break;
      case 2:
        if(decomp == xy_decomp){
          std::cout << "2D pencil decompostion in xy-dir." << std::endl;

          f_in.gnfast = gnx3_; 
          f_in.gnmid = gnx1_; 
          f_in.gnslow = gnx2_;
          f_in.np1 = np1_; f_in.np2 = np2_;
          f_in.ip1 = loc.lx1; f_in.ip2 = loc.lx2;

          swap1 = true;
          swap2 = false;
        }
        if(decomp == yz_decomp){
          std::cout << "2D pencil decompostion in yz-dir. " << std::endl;
          if(loc.lx1 != 0 ) std::cout << "Something wrong..." << std::endl;
/*
          f_in.gnfast = gnx1_; 
          f_in.gnmid = gnx3_; 
          f_in.gnslow = gnx2_;
          f_in.np1 = np3_; f_in.np2 = np2_;
          f_in.ip1 = loc.lx3; f_in.ip2 = loc.lx2;
          f_in.nfast = f_in.gnfast; 
          f_in.nmid  = f_in.gnmid/f_in.np1; 
          f_in.nslow = f_in.gnslow/f_in.np2;
          f_in.fs = 0;
          f_in.ms = f_in.ip1*f_in.nmid ;
          f_in.ss = f_in.ip2*f_in.nslow;
          f_in.fe = f_in.fs + f_in.nfast - 1;
          f_in.me = f_in.ms + f_in.nmid  - 1;
          f_in.se = f_in.ss + f_in.nslow - 1;

          f_out = f_in;
          f_out.nfast = f_out.gnfast/f_out.np1;
          f_out.nmid = f_out.gnmid/f_out.np2;
          f_out.nslow = f_out.gnslow;
          f_out.fs = f_out.ip1*f_out.nfast;
          f_out.ms = f_out.ip2*f_out.nmid;
          f_out.ss = 0;
          f_out.fe = f_out.fs + f_out.nfast - 1;
          f_out.me = f_out.ms + f_out.nmid  - 1;
          f_out.se = f_out.ss + f_out.nslow - 1;

          b_in.gnfast = f_out.gnslow;
          b_in.gnmid  = f_out.gnmid;
          b_in.gnslow = f_out.gnfast;
          b_in.np1 = f_in.np2; b_in.np2 = f_in.np1;
          b_in.ip1 = f_in.ip2; b_in.ip2 = f_in.ip1;
          b_in.nfast = b_in.gnfast;
          b_in.nmid = b_in.gnmid/b_in.np1;
          b_in.nslow = b_in.gnslow/b_in.np2;
          b_in.fs = 0;
          b_in.ms = b_in.ip1*b_in.nmid ;
          b_in.ss = b_in.ip2*b_in.nslow;
          b_in.fe = b_in.fs + b_in.nfast - 1;
          b_in.me = b_in.ms + b_in.nmid  - 1;
          b_in.se = b_in.ss + b_in.nslow - 1;

          b_out = b_in;
          b_out.nfast = b_out.gnfast/b_out.np1;
          b_out.nmid = b_out.gnmid/b_out.np2;
          b_out.nslow = b_out.gnslow;
          b_out.fs = b_out.ip1*b_out.nfast;
          b_out.ms = b_out.ip2*b_out.nmid;
          b_out.ss = 0;
          b_out.fe = b_out.fs + b_out.nfast - 1;
          b_out.me = b_out.ms + b_out.nmid  - 1;
          b_out.se = b_out.ss + b_out.nslow - 1;

          swap1 = true;
          swap2 = true;
          permute1 = 2;
          permute2 = 2;
*/

          f_in.gnfast = gnx1_; 
          f_in.gnmid = gnx2_; 
          f_in.gnslow = gnx3_;
          f_in.np1 = np2_; f_in.np2 = np3_;
          f_in.ip1 = loc.lx2; f_in.ip2 = loc.lx3;

          swap1 = false;
          swap2 = false;
        }
        if(decomp == xz_decomp){
          std::cout << "2D pencil decompostion in xz-dir." << std::endl;

          f_in.gnfast = gnx2_; 
          f_in.gnmid = gnx3_; 
          f_in.gnslow = gnx1_;
          f_in.np1 = np3_; f_in.np2 = np1_;
          f_in.ip1 = loc.lx3; f_in.ip2 = loc.lx1;
          swap1 = true;
          swap2 = false;
        }

        f_in.nfast = f_in.gnfast; 
        f_in.nmid  = f_in.gnmid/f_in.np1; 
        f_in.nslow = f_in.gnslow/f_in.np2;
        f_in.fs = 0;
        f_in.ms = f_in.ip1*f_in.nmid ;
        f_in.ss = f_in.ip2*f_in.nslow;
        f_in.fe = f_in.fs + f_in.nfast - 1;
        f_in.me = f_in.ms + f_in.nmid  - 1;
        f_in.se = f_in.ss + f_in.nslow - 1;

        f_out = f_in;
        f_out.nfast = f_out.gnfast/f_out.np1;
        f_out.nmid = f_out.gnmid/f_out.np2;
        f_out.nslow = f_out.gnslow;
        f_out.fs = f_out.ip1*f_out.nfast;
        f_out.ms = f_out.ip2*f_out.nmid;
        f_out.ss = 0;
        f_out.fe = f_out.fs + f_out.nfast - 1;
        f_out.me = f_out.ms + f_out.nmid  - 1;
        f_out.se = f_out.ss + f_out.nslow - 1;

        b_in = f_out;
        b_in.gnfast = f_out.gnslow;
        b_in.gnmid  = f_out.gnfast;
        b_in.gnslow = f_out.gnmid;
        b_in.nfast = b_in.gnfast;
        b_in.nmid = b_in.gnmid/b_in.np1;
        b_in.nslow = b_in.gnslow/b_in.np2;
        b_in.fs = 0;
        b_in.ms = b_in.ip1*b_in.nmid ;
        b_in.ss = b_in.ip2*b_in.nslow;
        b_in.fe = b_in.fs + b_in.nfast - 1;
        b_in.me = b_in.ms + b_in.nmid  - 1;
        b_in.se = b_in.ss + b_in.nslow - 1;

        b_out = b_in;
        b_out.nfast = b_out.gnfast/b_out.np1;
        b_out.nmid = b_out.gnmid;
        b_out.nslow = b_out.gnslow/b_out.np2;
        b_out.fs = b_out.ip1*b_out.nfast;
        b_out.ms = 0;
        b_out.ss = b_out.ip2*b_out.nslow;
        b_out.fe = b_out.fs + b_out.nfast - 1;
        b_out.me = b_out.ms + b_out.nmid  - 1;
        b_out.se = b_out.ss + b_out.nslow - 1;

        permute1 = 2;
        permute2 = 1;

        break;
      case 3:
        std::cout << "3D block mpi decompisiton at " 
                  << loc.lx1 << loc.lx2 << loc.lx3 << std::endl;

        f_in.gnfast = gnx1_; f_in.gnmid = gnx2_; f_in.gnslow = gnx3_;
        f_in.nfast = nx1; 
        f_in.nmid  = nx2;
        f_in.nslow = nx3;
        f_in.fs = idisp;
        f_in.ms = jdisp;
        f_in.ss = kdisp;
        f_in.fe = f_in.fs + f_in.nfast - 1;
        f_in.me = f_in.ms + f_in.nmid  - 1;
        f_in.se = f_in.ss + f_in.nslow - 1;
  
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
    dkx = 2*PI/(Real)b_in.gnfast;
    dky = 2*PI/(Real)b_in.gnmid;
    dkz = 2*PI/(Real)b_in.gnslow;
    } // end of using namespace block
  }
}

long int AthenaFFT::GetIndex(const int i, const int j, const int k)
{
  return k + nx3 * ( j + nx2 * i);
}


long int AthenaFFT::GetIndex(const int i, const int j, const int k, bool swap)
{
  {using namespace DecompositionNames;

  if(swap){
    if(decomp == xz_decomp)
      return j + nx2 * ( k + nx3 * i);
    else if (decomp == yz_decomp)
      return i + nx1 * ( j + nx2 * k);
    else if (decomp == xy_decomp)
      return k + nx3 * ( i + nx1 * j);
  } else return i + nx1 * ( j + nx2 * k);

  }
}

long int AthenaFFT::GetFreq(const int i, const int j, const int k, bool swap)
{
  if(swap) return j + knx2 * ( k + knx3 * i);
  else return i + knx1 * ( j + knx2 * k);
}


