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
    std::stringstream msg;
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
 
    gis_ = idisp+is; gie_= idisp+ie; 
    gjs_ = jdisp+js; gje_= jdisp+je; 
    gks_ = kdisp+ks; gke_= kdisp+ke; 
 
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
    fftw_free(work);
#ifdef MPI_PARALLEL
    if(MPIFFT_MODE == "fftw-mpi"){
      fftw_destroy_plan(fplan->plan);
      fftw_destroy_plan(bplan->plan);
      fftw_mpi_cleanup();
    } else if (MPIFFT_MODE == "accfft"){
      accfft_destroy_plan(fplan->plan);
      accfft_destroy_plan(bplan->plan);
      MPI_Comm_free(&comm_);
    } else if (MPIFFT_MODE == "plimpton"){
      fft_3d_destroy_plan(plan->plan);
    }
#else
    fftw_destroy_plan(fplan->plan);
    fftw_destroy_plan(bplan->plan);
#endif
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

// plan 1D fft
AthenaFFTPlan *AthenaFFT::CreatePlan(int nx1, AthenaFFTComplex *data, enum AthenaFFTDirection dir)
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

// plan 2D fft
AthenaFFTPlan *AthenaFFT::CreatePlan(int nx1, int nx2, AthenaFFTComplex *data, enum AthenaFFTDirection dir)
{
  AthenaFFTPlan *plan;

  if(FFT_ENABLED){
    plan = new AthenaFFTPlan;
    plan->dir = dir;
    plan->dim = 2;
#ifdef MPI_PARALLEL
    if(MPIFFT_MODE == "plimpton"){
      int nbuf;
      plan->plan = fft_2d_create_plan(MPI_COMM_WORLD,nx2,nx1,
                                gjs_,gje_,gis_,gie_,
                                gjs_,gje_,gis_,gie_,
                                0, 0, &nbuf);
    } else if (MPIFFT_MODE == "fftw-mpi"){
      if(dir == AthenaFFTForward)
        plan->plan = fftw_mpi_plan_dft_2d(nx1,nx2,data, data, 
                                    MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE);
      else
        plan->plan = fftw_mpi_plan_dft_2d(nx1,nx2,data, data, 
                                    MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE);
    } else if (MPIFFT_MODE == "accfft"){
      int gnx[2] = {nx1, nx2};
      if(dir == AthenaFFTForward) plan->dir = ACCFFT_FORWARD;
      else plan->dir = ACCFFT_BACKWARD;
      plan->plan = accfft_plan_dft_2d_c2c(gnx,data,data,comm_,ACCFFT_MEASURE);
    }
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
AthenaFFTPlan *AthenaFFT::CreatePlan(int nx1, int nx2, int nx3, AthenaFFTComplex *data, enum AthenaFFTDirection dir)
{
  AthenaFFTPlan *plan;

  if(FFT_ENABLED){
    plan = new AthenaFFTPlan;
    plan->dir = dir;
    plan->dim = 3;
#ifdef MPI_PARALLEL
    if(MPIFFT_MODE == "plimpton"){
      int nbuf;
      plan->plan = fft_3d_create_plan(MPI_COMM_WORLD,nx3,nx2,nx1,
                                gks_,gke_,gjs_,gje_,gis_,gie_,
                                gks_,gke_,gjs_,gje_,gis_,gie_,
                                0, 0, &nbuf);
    } else if (MPIFFT_MODE == "fftw-mpi"){
      if(dir == AthenaFFTForward)
        plan->plan = fftw_mpi_plan_dft_3d(nx1,nx2,nx3,data, data, 
                                    MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE);
      else
        plan->plan = fftw_mpi_plan_dft_3d(nx1,nx2,nx3,data, data, 
                                    MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE);
    } else if (MPIFFT_MODE == "accfft"){
      int gnx[3] = {nx1, nx2, nx3};
      if(dir == AthenaFFTForward) plan->dir = ACCFFT_FORWARD;
      else plan->dir = ACCFFT_BACKWARD;
      plan->plan = accfft_plan_dft_3d_c2c(gnx,data,data,comm_,ACCFFT_MEASURE);
    }
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
    if(MPIFFT_MODE == "plimpton"){
      if(plan->dim == 3) fft_3d(data, data, plan->dir, plan->plan);
      if(plan->dim == 2) fft_2d(data, data, plan->dir, plan->plan);
    } else if (MPIFFT_MODE == "fftw-mpi"){
      fftw_execute_dft(plan->plan, data, data);
    } else if (MPIFFT_MODE == "accfft"){
      accfft_execute_c2c(plan->plan, plan->dir, data, data);
    }
#else
    fftw_execute_dft(plan->plan, data, data);
#endif
  }
}

void AthenaFFT::Execute(AthenaFFTPlan *plan, AthenaFFTComplex *data)
{
  if(FFT_ENABLED){
#ifdef MPI_PARALLEL
    if(MPIFFT_MODE == "plimpton"){
      if(plan->dim == 3) fft_3d(data, data, plan->dir, plan->plan);
      if(plan->dim == 2) fft_2d(data, data, plan->dir, plan->plan);
    } else if (MPIFFT_MODE == "fftw-mpi"){
      fftw_execute_dft(plan->plan, data, data);
    } else if (MPIFFT_MODE == "accfft"){
      accfft_execute_c2c(plan->plan, plan->dir, data, data);
    }
#else
    fftw_execute_dft(plan->plan, data, data);
#endif
  }
}

void AthenaFFT::Initialize()
{
  if(FFT_ENABLED){
#ifdef MPI_PARALLEL
    if(MPIFFT_MODE == "fftw-mpi"){
      fftw_mpi_init();
    } else if (MPIFFT_MODE == "accfft"){
      int np[2];
      if(np1_ == 1) 
        np = { np2_, np3_ };
      else if (np2_ == 1)
        np = { np1_, np3_ };
      else if (np3_ == 1)
        np = { np1_, np2_ };
      else {
        msg << "ACCFFT requires pencil decomposition!" << std:: endl
            << "Current MPI Configureation is " 
            << np1_ << " x " << np2_ << " x " << np3_ << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
 
      AthenaFFTInt c_dims[2] = { np[0], np[1] };
      accfft_create_comm(MPI_COMM_WORLD, c_dims, &comm_);
    }
#endif

#ifdef OPENMP_PARALLEL
    int threads_ok;
    threads_ok = fftw_init_threads();
    if(threads_ok) fftw_plan_with_nthreads(nthreads_);
#endif
  }
}

int AthenaFFT::GetIndex(const int i)
{
  return i;
}
int AthenaFFT::GetIndex(const int j, const int i)
{
  return j + nx2_*i;
}
int AthenaFFT::GetIndex(const int k, const int j, const int i)
{
  return k + nx3_*(j + nx2_*i);
}
