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

    dim_ = 1;
    if(mesh_size.nx2 > 1) dim_=2;
    if(mesh_size.nx3 > 1) dim_=3;
 
    orig_idx = new AthenaFFTIndex(dim_,pmb);

    cnt = block_size.nx1*block_size.nx2*block_size.nx3;
    gcnt = mesh_size.nx1*mesh_size.nx2*mesh_size.nx3;
 
    fplan = new AthenaFFTPlan;
    bplan = new AthenaFFTPlan;
 
    in = new AthenaFFTComplex[cnt];
    out = new AthenaFFTComplex[cnt];

    nthreads_ = pmy_block->pmy_mesh->GetNumMeshThreads();

#ifdef MPI_PARALLEL
    {using namespace DecompositionNames;
    decomp_ = 0; pdim_ = 0;
    if(pm->nrbx1 > 1){
      decomp_ = decomp_ | x_decomp;
      pdim_++;
    }
    if(pm->nrbx2 > 1){
      decomp_ = decomp_ | y_decomp;
      pdim_++;
    }
    if(pm->nrbx3 > 1){
      decomp_ = decomp_ | z_decomp;
      pdim_++;
    }

    MpiInitialize();
    }
#else
    f_in = new AthenaFFTIndex(orig_idx);
    f_out = new AthenaFFTIndex(orig_idx);
    b_in = new AthenaFFTIndex(orig_idx);
    b_out = new AthenaFFTIndex(orig_idx);
#endif
    for(int i=0;i<dim_;i++){
      Nx[f_in->iloc[i]]=f_in->Nx[i];
      nx[f_in->iloc[i]]=f_in->nx[i];
      disp[f_in->iloc[i]]=f_in->is[i];
      knx[b_in->iloc[i]]=b_in->nx[i];
      kdisp[b_in->iloc[i]]=b_in->is[i];
      dkx[b_in->iloc[i]]=2*PI/(Real)b_in->Nx[i];
    }
  }
}

// destructor

AthenaFFT::~AthenaFFT()
{
  if(FFT_ENABLED){
    MpiCleanup();
    delete fplan;
    delete bplan;
    delete[] in;
    delete[] out;
  }
}

AthenaFFTPlan *AthenaFFT::QuickCreatePlan(enum AthenaFFTDirection dir,
                                          AthenaFFTComplex *work)
{
  if(FFT_ENABLED){
    int nfast,nmid,nslow;
    if(dir == AthenaFFTForward){
      nfast = f_in->Nx[0]; nmid = f_in->Nx[1]; nslow = f_in->Nx[2];
    } else {
      nfast = b_in->Nx[0]; nmid = b_in->Nx[1]; nslow = b_in->Nx[2];
    }
    if(dim_==3) return CreatePlan(nfast,nmid,nslow,work,dir);
    else if(dim_==2) return CreatePlan(nfast,nslow,work,dir);
    else  return CreatePlan(nfast,work,dir);
  }
}

long int AthenaFFT::GetIndex(const int i, const int j, const int k)
{
  return i + nx[0] * ( j + nx[1] * k);
}

long int AthenaFFT::GetGlobalIndex(const int i, const int j, const int k)
{
  return i + disp[0] + Nx[0]* ( j + disp[1] + Nx[1]* ( k + disp[2]) );
}

AthenaFFTIndex::AthenaFFTIndex(int dim, MeshBlock *pmb)
{
  pmy_block = pmb;
  Mesh *pm=pmy_block->pmy_mesh;
  RegionSize& mesh_size  = pmy_block->pmy_mesh->mesh_size;
  RegionSize& block_size = pmy_block->block_size;
  LogicalLocation& loc = pmy_block->loc;
  dim_=dim;

  Nx = new int[dim_];
  np = new int[dim_];
  ip = new int[dim_];

  nx = new int[dim_];
  is = new int[dim_];
  ie = new int[dim_];

  Nx[0] = mesh_size.nx1;
  np[0] = mesh_size.nx1/block_size.nx1;
  ip[0] = loc.lx1;
  if(dim_ > 1){
    Nx[1] = mesh_size.nx2;
    np[1] = mesh_size.nx2/block_size.nx2;
    ip[1] = loc.lx2; 
  }
  if(dim_ > 2){
    Nx[2] = mesh_size.nx3;
    np[2] = mesh_size.nx3/block_size.nx3;
    ip[2] = loc.lx3;
  }

  iloc[0]=0;
  iloc[1]=1;
  iloc[2]=2;

  ploc[0]=0;
  ploc[1]=1;
  ploc[2]=2;

  SetLocalIndex();
}

// copy constructor
AthenaFFTIndex::AthenaFFTIndex(const AthenaFFTIndex *psrc){
  dim_ = psrc->dim_;
  Nx = new int[dim_];
  np = new int[dim_];
  ip = new int[dim_];

  nx = new int[dim_];
  is = new int[dim_];
  ie = new int[dim_];

  for(int i=0; i<dim_; i++){
    Nx[i]=psrc->Nx[i];
    np[i]=psrc->np[i];
    ip[i]=psrc->ip[i];
    iloc[i] = psrc->iloc[i];
    ploc[i] = psrc->ploc[i];
  }

  SetLocalIndex();
}

AthenaFFTIndex::~AthenaFFTIndex()
{
  delete[] Nx;
  delete[] np;
  delete[] ip;
  delete[] nx;
  delete[] is;
  delete[] ie;
  delete[] iloc;
  delete[] ploc;
}

void AthenaFFTIndex::SetLocalIndex(){
  for(int i=0; i<dim_; i++){
    nx[i] = Nx[i]/np[i];
    is[i] = ip[i]*nx[i];
    ie[i] = is[i]+nx[i]-1;
  }
//  PrintIndex();
}

void AthenaFFTIndex::Swap_(int loc[],int ref_axis){
  int tmp;
  int axis1=(ref_axis+1) % dim_, axis2=ref_axis+2 % dim_;
  tmp=loc[axis1]; 
  loc[axis1]=loc[axis2]; 
  loc[axis2]=tmp;
}

void AthenaFFTIndex::SwapAxis(int ref_axis){
  Swap_(iloc,ref_axis);
  Swap_(Nx,ref_axis);
}

void AthenaFFTIndex::SwapProc(int ref_axis){
  Swap_(ploc,ref_axis);
  Swap_(ip,ref_axis);
  Swap_(np,ref_axis);
}

void AthenaFFTIndex::Permute_(int loc[], int npermute){
  int tmp;
  for(int i=0; i<npermute; i++){
    tmp=loc[0];
    loc[0]=loc[1];
    loc[1]=loc[2];
    loc[2]=tmp;
  } 
}

void AthenaFFTIndex::PermuteAxis(int npermute){
  Permute_(iloc,npermute);
  Permute_(Nx,npermute);
}

void AthenaFFTIndex::PermuteProc(int npermute){
  Permute_(ploc,npermute);
  Permute_(np,npermute);
  Permute_(ip,npermute);
}

void AthenaFFTIndex::RemapArray_(int arr[], int loc[], int dir){
  int tmp[dim_];
  for(int i=0; i<dim_; i++) tmp[i]=arr[i];
  for(int i=0; i<dim_; i++){
    if(dir>0) arr[loc[i]]=tmp[i];
    else arr[i]=tmp[loc[i]];
  }
}

void AthenaFFTIndex::RemapAxis(int dir){
  RemapArray_(Nx,iloc,dir);
}

void AthenaFFTIndex::RemapProc(int dir){
  RemapArray_(np,ploc,dir);
  RemapArray_(ip,ploc,dir);
}

void AthenaFFTIndex::PrintIndex(void){
  std::cout << "Nx:" << Nx[0] << " "  << Nx[1] << " " << Nx[2] << std::endl
            << "np:" << np[0] << " "  << np[1] << " " << np[2] << std::endl
            << "ip:" << ip[0] << " "  << ip[1] << " " << ip[2] << std::endl
            << "nx:" << nx[0] << " "  << nx[1] << " " << nx[2] << std::endl
            << "is:" << is[0] << " "  << is[1] << " " << is[2] << std::endl
            << "ie:" << ie[0] << " "  << ie[1] << " " << ie[2] << std::endl
            << "iloc:" << iloc[0] << " "  << iloc[1] << " " << iloc[2] << std::endl
            << "ploc:" << ploc[0] << " "  << ploc[1] << " " << ploc[2] << std::endl;
}
