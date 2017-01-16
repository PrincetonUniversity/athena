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

// constructor, initializes data structures and parameters

AthenaFFT::AthenaFFT(MeshBlock *pmb)
{
  pmy_block = pmb;  
  if (FFT_ENABLED) {
    Mesh *pm=pmy_block->pmy_mesh;
    RegionSize& mesh_size  = pmy_block->pmy_mesh->mesh_size;
    RegionSize& block_size = pmy_block->block_size;
    LogicalLocation& loc = pmy_block->loc;
    int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
    int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
    nx1_ = block_size.nx1;
    nx2_ = block_size.nx2;
    nx3_ = block_size.nx3;
 
    np1_ = pm->nrbx1;
    np2_ = pm->nrbx2;
    np3_ = pm->nrbx3;
 
    int dim = 1;
    if(mesh_size.nx2 > 1) dim=2;
    if(mesh_size.nx3 > 1) dim=3;
 
    int idisp = loc.lx1*nx1_;
    int jdisp = loc.lx2*nx2_;
    int kdisp = loc.lx3*nx3_;
 
    gis_ = idisp; gie_= idisp+nx1_; 
    gjs_ = jdisp; gje_= jdisp+nx2_; 
    gks_ = kdisp; gke_= kdisp+nx3_; 
 
    gnx1_ = mesh_size.nx1;
    gnx2_ = mesh_size.nx2;
    gnx3_ = mesh_size.nx3;
 
    cnt = nx1_*nx2_*nx3_;
    gcnt = gnx1_*gnx2_*gnx3_;
 
//    work = (AthenaFFTComplex *) fftw_malloc(sizeof(AthenaFFTComplex *) * cnt);
    work = new AthenaFFTComplex[cnt];
//    work->NewAthenaArray(nx3,nx2,nx1);
    fplan = new AthenaFFTPlan;
    bplan = new AthenaFFTPlan;
 
    nthreads_ = pmy_block->pmy_mesh->GetNumMeshThreads();
  }
}

// destructor

AthenaFFT::~AthenaFFT()
{
  if(FFT_ENABLED) {
    MpiCleanup();
    delete work;
    delete fplan;
    delete bplan;
#ifdef OPENMP_PARALLEL
    fftw_cleanup_threads();
#endif
    fftw_cleanup();
  }
}

AthenaFFTPlan *AthenaFFT::QuickCreatePlan(enum AthenaFFTDirection dir)
{
  if(FFT_ENABLED){
    if(gnx3_ > 1) return CreatePlan(gnx1_,gnx2_,gnx3_,work,dir);
    else if(gnx2_ > 1) return CreatePlan(gnx1_,gnx2_,work,dir);
    else  return CreatePlan(gnx1_,work,dir);
  }
}


int AthenaFFT::GetIndex(const int i)
{
  return i;
}
int AthenaFFT::GetIndex(const int j, const int i)
{
  return j + nx2_*i;
//  return i + nx1_*j;
}
int AthenaFFT::GetIndex(const int k, const int j, const int i)
{
  return k + nx3_*(j + nx2_*i);
//  return i + nx1_*(j + nx2_*k);
}

// plan 1D fft
AthenaFFTPlan *AthenaFFT::CreatePlan(AthenaFFTInt nx1, AthenaFFTComplex *data, 
                                     enum AthenaFFTDirection dir)
{
  AthenaFFTPlan *plan;
  if(FFT_ENABLED){
    plan = new AthenaFFTPlan;
    plan->dir = dir;
    plan->dim = 1;
    if(dir == AthenaFFTForward)
      plan->plan = fftw_plan_dft_1d(nx1, data, data, FFTW_FORWARD, FFTW_ESTIMATE);
    else
      plan->plan = fftw_plan_dft_1d(nx1, data, data, FFTW_BACKWARD, FFTW_ESTIMATE);
  }
  return plan;
}


