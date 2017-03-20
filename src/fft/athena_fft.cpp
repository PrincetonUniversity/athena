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
#include "../globals.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"

// constructor, initializes data structures and parameters

AthenaFFT::AthenaFFT(MeshBlock *pmb)
{
  pmy_block = pmb;  
#ifdef FFT
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
#ifdef OPENMP_PARALLEL
    int threads_ok;
    threads_ok = fftw_init_threads();
    if(threads_ok) fftw_plan_with_nthreads(nthreads_);
#endif

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
    f_in  = new AthenaFFTIndex(dim_,pmb);
    f_out = new AthenaFFTIndex(dim_,pmb);
    b_in  = new AthenaFFTIndex(dim_,pmb);
    b_out = new AthenaFFTIndex(dim_,pmb);
#endif
    for(int i=0;i<3;i++){
      Nx[f_in->iloc[i]]=f_in->Nx[i];
      nx[f_in->iloc[i]]=f_in->nx[i];
      disp[f_in->iloc[i]]=f_in->is[i];
      knx[b_in->iloc[i]]=b_in->nx[i];
      kdisp[b_in->iloc[i]]=b_in->is[i];
      dkx[b_in->iloc[i]]=2*PI/(Real)b_in->Nx[i];
    }
#endif
}

// destructor

AthenaFFT::~AthenaFFT()
{
#ifdef FFT
    delete fplan;
    delete bplan;
    delete[] in;
    delete[] out;
    delete orig_idx;
    delete f_in;
    delete f_out;
    delete b_in;
    delete b_out;
#endif
}

AthenaFFTPlan *AthenaFFT::QuickCreatePlan(enum AthenaFFTDirection dir,
                                          AthenaFFTComplex *work)
{
  int nfast,nmid,nslow;
  if(dir == AthenaFFTForward){
    nfast = f_in->Nx[0]; nmid = f_in->Nx[1]; nslow = f_in->Nx[2];
  } else {
    nfast = b_in->Nx[0]; nmid = b_in->Nx[1]; nslow = b_in->Nx[2];
  }
  if(dim_==3) return CreatePlan(nfast,nmid,nslow,work,dir);
  else if(dim_==2) return CreatePlan(nfast,nmid,work,dir);
  else  return CreatePlan(nfast,work,dir);
}

long int AthenaFFT::GetIndex(const int i, const int j, const int k)
{
  return i + nx[0] * ( j + nx[1] * k);
}

long int AthenaFFT::GetGlobalIndex(const int i, const int j, const int k)
{
  return i + disp[0] + Nx[0]* ( j + disp[1] + Nx[1]* ( k + disp[2]) );
}

long int AthenaFFT::GetIndex(const int i, const int j, const int k, AthenaFFTIndex *pidx)
{
  int old_idx[3]={i,j,k};
  int new_idx[3];
  new_idx[0]=old_idx[pidx->iloc[0]];
  new_idx[1]=old_idx[pidx->iloc[1]];
  new_idx[2]=old_idx[pidx->iloc[2]];

  return new_idx[0] + pidx->nx[0] * (new_idx[1] + pidx->nx[1] * (new_idx[2]));
}

// plan 1D fft
AthenaFFTPlan *AthenaFFT::CreatePlan(int nfast, AthenaFFTComplex *data, 
                                     enum AthenaFFTDirection dir)
{
  AthenaFFTPlan *plan;
#ifdef FFT
  plan = new AthenaFFTPlan;
  plan->dir = dir;
  plan->dim = dim_;
  if(dir == AthenaFFTForward)
    plan->plan = fftw_plan_dft_1d(nfast, data, data, FFTW_FORWARD, FFTW_ESTIMATE);
  else
    plan->plan = fftw_plan_dft_1d(nfast, data, data, FFTW_BACKWARD, FFTW_ESTIMATE);
#endif
  return plan;
}

// plan 2D fft
AthenaFFTPlan *AthenaFFT::CreatePlan(int nfast, int nslow, 
                                     AthenaFFTComplex *data, enum AthenaFFTDirection dir)
{
  AthenaFFTPlan *plan;

#ifdef FFT
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
    plan->plan = fftw_plan_dft_2d(nslow,nfast,data, data, FFTW_FORWARD, FFTW_MEASURE);
  else
    plan->plan = fftw_plan_dft_2d(nslow,nfast,data, data, FFTW_BACKWARD, FFTW_MEASURE);
#endif
#endif // FFT

  return plan;
}

// plan 3D fft
AthenaFFTPlan *AthenaFFT::CreatePlan(int nfast, int nmid, int nslow,
                                     AthenaFFTComplex *data, enum AthenaFFTDirection dir)
{
  AthenaFFTPlan *plan;

#ifdef FFT
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
    plan->plan = fftw_plan_dft_3d(nslow,nmid,nfast,data, data, FFTW_FORWARD, FFTW_MEASURE);
  } else {
    plan->plan = fftw_plan_dft_3d(nslow,nmid,nfast,data, data, FFTW_BACKWARD, FFTW_MEASURE);
  }
#endif
#endif // FFT

  return plan;
}

void AthenaFFT::Execute(AthenaFFTPlan *plan)
{
#ifdef FFT
#ifdef MPI_PARALLEL
    if(plan->dim == 3) fft_3d(in, out, plan->dir, plan->plan3d);
    if(plan->dim == 2) fft_2d(in, out, plan->dir, plan->plan2d);
#else
    fftw_execute_dft(plan->plan, in, out);
#endif
#endif //FFT
}

void AthenaFFT::Execute(AthenaFFTPlan *plan, AthenaFFTComplex *data)
{
#ifdef FFT
#ifdef MPI_PARALLEL
    if(plan->dim == 3) fft_3d(data, data, plan->dir, plan->plan3d);
    if(plan->dim == 2) fft_2d(data, data, plan->dir, plan->plan2d);
#else
    fftw_execute_dft(plan->plan, data, data);
#endif
#endif // FFT
}

void AthenaFFT::Execute(AthenaFFTPlan *plan, AthenaFFTComplex *in_data, AthenaFFTComplex *out_data)
{
#ifdef FFT
#ifdef MPI_PARALLEL
  if(plan->dim == 3) fft_3d(in_data, out_data, plan->dir, plan->plan3d);
  if(plan->dim == 2) fft_2d(in_data, out_data, plan->dir, plan->plan2d);
#else
  fftw_execute_dft(plan->plan, in_data, out_data);
#endif
#endif
}

void AthenaFFT::MpiInitialize()
{
#ifdef MPI_PARALLEL
  std::stringstream msg;
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
#endif
}
//---------------------------------------------------------------------------------------
// AthenaFFTIndex class:
AthenaFFTIndex::AthenaFFTIndex(int dim, MeshBlock *pmb)
{
  pmy_block = pmb;
  Mesh *pm=pmy_block->pmy_mesh;
  RegionSize& mesh_size  = pmy_block->pmy_mesh->mesh_size;
  RegionSize& block_size = pmy_block->block_size;
  LogicalLocation& loc = pmy_block->loc;
  dim_=3;

  Nx = new int[dim_];
  np = new int[dim_];
  ip = new int[dim_];

  nx = new int[dim_];
  is = new int[dim_];
  ie = new int[dim_];

  Nx[0] = mesh_size.nx1;
  np[0] = mesh_size.nx1/block_size.nx1;
  ip[0] = loc.lx1;
  iloc[0]=0;
  ploc[0]=0;
  if(dim_ > 1){
    Nx[1] = mesh_size.nx2;
    np[1] = mesh_size.nx2/block_size.nx2;
    ip[1] = loc.lx2; 
    iloc[1]=1;
    ploc[1]=1;
  }
  if(dim_ > 2){
    Nx[2] = mesh_size.nx3;
    np[2] = mesh_size.nx3/block_size.nx3;
    ip[2] = loc.lx3;
    iloc[2]=2;
    ploc[2]=2;
  }

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
