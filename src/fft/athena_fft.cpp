//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file athena_fft.cpp
//  \brief

// C/C++ headers
#include <iostream>
#include <sstream>
#include <stdexcept>  // runtime_error

// Athena++ headers
#include "../globals.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../mesh/meshblock_tree.hpp"
#include "athena_fft.hpp"

// constructor, initializes data structures and parameters

FFTBlock::FFTBlock(FFTDriver *pfd, LogicalLocation iloc, int igid,
            RegionSize msize, RegionSize bsize) {
  pmy_driver_ = pfd;
  loc_=iloc;
  gid_=igid;
  msize_=msize;
  bsize_=bsize;
  dx1=(msize_.x1max-msize_.x1min)/static_cast<Real>(msize_.nx1);
  dx2=(msize_.x2max-msize_.x2min)/static_cast<Real>(msize_.nx2);
  dx3=(msize_.x3max-msize_.x3min)/static_cast<Real>(msize_.nx3);
  fplan_=NULL;
  bplan_=NULL;

  cnt_ = bsize_.nx1*bsize_.nx2*bsize_.nx3;
  gcnt_ = pmy_driver_->gcnt_;

  dim_=pmy_driver_->dim_;

  norm_factor_ = 1.0;

  in_ = new AthenaFFTComplex[cnt_];
  out_ = new AthenaFFTComplex[cnt_];

  orig_idx_ = new AthenaFFTIndex(dim_,loc_,msize_,bsize_);

#ifdef MPI_PARALLEL
  decomp_=pmy_driver_->decomp_;
  pdim_=pmy_driver_->pdim_;
  MpiInitialize();
#else
  f_in_  = new AthenaFFTIndex(orig_idx_);
  f_out_ = new AthenaFFTIndex(orig_idx_);
  b_in_  = new AthenaFFTIndex(orig_idx_);
  b_out_ = new AthenaFFTIndex(orig_idx_);
#endif

//  f_in_->PrintIndex();
#ifdef FFT
  for (int i=0;i<3;i++) {
    Nx[f_in_->iloc[i]]=f_in_->Nx[i];
    nx[f_in_->iloc[i]]=f_in_->nx[i];
    disp[f_in_->iloc[i]]=f_in_->is[i];
    kNx[b_in_->iloc[i]]=b_in_->Nx[i];
    knx[b_in_->iloc[i]]=b_in_->nx[i];
    kdisp[b_in_->iloc[i]]=b_in_->is[i];
    dkx[b_in_->iloc[i]]=2*PI/b_in_->Lx[i];
  }
#endif
}

// destructor

FFTBlock::~FFTBlock() {
    delete[] in_;
    delete[] out_;
    delete orig_idx_;
    delete f_in_;
    delete f_out_;
    delete b_in_;
    delete b_out_;
    if (fplan_!=NULL) DestroyPlan(fplan_);
    if (bplan_!=NULL) DestroyPlan(bplan_);
}

void FFTBlock::DestroyPlan(AthenaFFTPlan *plan) {
#ifdef FFT
#ifdef MPI_PARALLEL
  if (plan->plan3d != NULL) fft_3d_destroy_plan(plan->plan3d);
  if (plan->plan2d != NULL) fft_2d_destroy_plan(plan->plan2d);
#endif
  if (plan->plan != NULL) fftw_destroy_plan(plan->plan);
  delete plan;
#endif
}

void FFTBlock::PrintSource(int in) {
  std::cout << Nx[0] << "x" << Nx[1] << "x" << Nx[2] << std::endl;

  for (int k=0; k<Nx[2]; ++k) {
    for (int j=0; j<Nx[1]; ++j) {
      for (int i=0; i<Nx[0]; ++i) {
        int64_t idx=GetIndex(i,j,k,f_in_);
        if (in == 1) std::cout << in_[idx][0] << " ";
        if (in == -1) std::cout << out_[idx][0] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}
int64_t FFTBlock::GetIndex(const int i, const int j, const int k) {
  return i + nx[0] * ( j + nx[1] * k);
}

int64_t FFTBlock::GetGlobalIndex(const int i, const int j, const int k) {
  return i + disp[0] + Nx[0]* ( j + disp[1] + Nx[1]* ( k + disp[2]) );
}

int64_t FFTBlock::GetIndex(const int i, const int j, const int k, AthenaFFTIndex *pidx) {
  int old_idx[3]={i,j,k};
  int new_idx[3];
  new_idx[0]=old_idx[pidx->iloc[0]];
  new_idx[1]=old_idx[pidx->iloc[1]];
  new_idx[2]=old_idx[pidx->iloc[2]];

  return new_idx[0] + pidx->nx[0] * (new_idx[1] + pidx->nx[1] * (new_idx[2]));
}

//----------------------------------------------------------------------------------------
//! \fn void FFTBlock::RetrieveResult(const AthenaArray<Real> &src, int ns)
//  \brief Fill the result in the active zone

void FFTBlock::RetrieveResult(AthenaArray<Real> &dst, int ns, int ngh,
                              LogicalLocation loc, RegionSize bsize) {
  const AthenaFFTComplex *src=out_;
  int is, ie, js, je, ks, ke;
  // possible loss of accuracy from int64_t loc.lx1 to int is, e.g.
  is = static_cast<int>(loc.lx1*bsize.nx1-loc_.lx1*bsize_.nx1);
  js = static_cast<int>(loc.lx2*bsize.nx2-loc_.lx2*bsize_.nx2);
  ks = static_cast<int>(loc.lx3*bsize.nx3-loc_.lx3*bsize_.nx3);
  ie = is+bsize.nx1-1;
  je = bsize.nx2>1 ? js+bsize.nx2-1:js;
  ke = bsize.nx3>1 ? ks+bsize.nx3-1:ks;
  int jl = bsize.nx2>1 ? ngh:0;
  int kl = bsize.nx3>1 ? ngh:0;

  for (int n=0; n<ns; n++) {
    for (int k=kl, mk=ks; mk<=ke; k++, mk++) {
      for (int j=jl, mj=js; mj<=je; j++, mj++) {
        for (int i=ngh, mi=is; mi<=ie; i++, mi++) {
          int64_t idx=GetIndex(mi,mj,mk,b_out_);
          if (ns == 1) {
            dst(k,j,i)=src[idx][0]*norm_factor_;
          } else {
            dst(n,k,j,i)=src[idx][n]*norm_factor_;
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FFTBlock::LoadSource(const AthenaArray<Real> &src, int ns)
//  \brief Fill the source in the active zone

void FFTBlock::LoadSource(const AthenaArray<Real> &src, int ns, int ngh,
                          LogicalLocation loc, RegionSize bsize) {
  AthenaFFTComplex *dst=in_;
  int is, ie, js, je, ks, ke;
  // possible loss of accuracy from int64_t loc.lx1 to int is, e.g.
  is = static_cast<int>(loc.lx1*bsize.nx1-loc_.lx1*bsize_.nx1);
  js = static_cast<int>(loc.lx2*bsize.nx2-loc_.lx2*bsize_.nx2);
  ks = static_cast<int>(loc.lx3*bsize.nx3-loc_.lx3*bsize_.nx3);
  ie = is+bsize.nx1-1;
  je = bsize.nx2>1 ? js+bsize.nx2-1:js;
  ke = bsize.nx3>1 ? ks+bsize.nx3-1:ks;
  int jl = bsize.nx2>1 ? ngh:0;
  int kl = bsize.nx3>1 ? ngh:0;

  for (int n=0; n<ns; n++) {
    for (int k=kl, mk=ks; mk<=ke; k++, mk++) {
      for (int j=jl, mj=js; mj<=je; j++, mj++) {
        for (int i=ngh, mi=is; mi<=ie; i++, mi++) {
          int64_t idx=GetIndex(mi,mj,mk,f_in_);
          if (ns == 1) {
            dst[idx][0]=src(n,k,j,i);
            dst[idx][1]=0.0;
          } else {
            dst[idx][n]=src(n,k,j,i);
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FFTBlock::ApplyKernel(int mode)
//  \brief Apply kernel

void FFTBlock::ApplyKernel(int mode) {
  for (int k=0; k<knx[2]; k++) {
    for (int j=0; j<knx[1]; j++) {
      for (int i=0; i<knx[0]; i++) {
        int64_t idx_in=GetIndex(i,j,k,b_in_);
        int64_t idx_out=GetIndex(i,j,k,f_out_);
        in_[idx_in][0] = out_[idx_out][0];
        in_[idx_in][1] = out_[idx_out][1];
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn AthenaFFTPlan *FFTBlock::QuickCreatePlan(AthenaFFTComplex *data,
//                                                enum AthenaFFTDirection dir)
//  \brief initialize FFT plan using mesh information

AthenaFFTPlan *FFTBlock::QuickCreatePlan(AthenaFFTComplex *data,
                                          enum AthenaFFTDirection dir) {
  int nfast,nmid,nslow;
  if (dir == AthenaFFTForward) {
    nfast = f_in_->Nx[0]; nmid = f_in_->Nx[1]; nslow = f_in_->Nx[2];
  } else {
    nfast = b_in_->Nx[0]; nmid = b_in_->Nx[1]; nslow = b_in_->Nx[2];
  }
  if (dim_==3) return CreatePlan(nfast,nmid,nslow,data,dir);
  else if (dim_==2) return CreatePlan(nfast,nmid,data,dir);
  else  return CreatePlan(nfast,data,dir);
}

//----------------------------------------------------------------------------------------
//! \fn AthenaFFTPlan *FFTBlock::CreatePlan(int nfast, AthenaFFTComplex *data,
//                                           enum AthenaFFTDirection dir)
//  \brief initialize FFT plan for 1D FFT
AthenaFFTPlan *FFTBlock::CreatePlan(int nfast, AthenaFFTComplex *data,
                                     enum AthenaFFTDirection dir) {
  AthenaFFTPlan *plan = NULL;
#ifdef FFT
  plan = new AthenaFFTPlan;
  plan->dir = dir;
  plan->dim = dim_;
  if (dir == AthenaFFTForward)
    plan->plan = fftw_plan_dft_1d(nfast, data, data, FFTW_FORWARD, FFTW_ESTIMATE);
  else
    plan->plan = fftw_plan_dft_1d(nfast, data, data, FFTW_BACKWARD, FFTW_ESTIMATE);
#endif
  return plan;
}

//----------------------------------------------------------------------------------------
//! \fn AthenaFFTPlan *FFTBlock::CreatePlan(int nfast, int nslow,
//                                           AthenaFFTComplex *data,
//                                           enum AthenaFFTDirection dir)
//  \brief initialize FFT plan for 2D FFT
AthenaFFTPlan *FFTBlock::CreatePlan(int nfast, int nslow,
                                     AthenaFFTComplex *data,
                                     enum AthenaFFTDirection dir) {
  AthenaFFTPlan *plan = NULL;

#ifdef FFT
  plan = new AthenaFFTPlan;
  plan->dir = dir;
  plan->dim = dim_;
#ifdef MPI_PARALLEL
  int nbuf;
  if (dir == AthenaFFTForward) {
    plan->dir = FFTW_FORWARD;
    plan->plan2d = fft_2d_create_plan(MPI_COMM_WORLD,nfast,nslow,
                                      f_in_->is[0],f_in_->ie[0],
                                      f_in_->is[1],f_in_->ie[1],
                                      f_out_->is[f_in_->iloc[0]],
                                      f_out_->ie[f_in_->iloc[0]],
                                      f_out_->is[f_in_->iloc[1]],
                                      f_out_->ie[f_in_->iloc[1]],
                                      0, permute1_, &nbuf);
  } else {
    plan->dir = FFTW_BACKWARD;
    plan->plan2d = fft_2d_create_plan(MPI_COMM_WORLD,nfast,nslow,
                                      b_in_->is[0],b_in_->ie[0],
                                      b_in_->is[1],b_in_->ie[1],
                                      b_out_->is[b_in_->iloc[0]],
                                      b_out_->ie[b_in_->iloc[0]],
                                      b_out_->is[b_in_->iloc[1]],
                                      b_out_->ie[b_in_->iloc[1]],
                                      0, permute2_, &nbuf);
  }
  plan->plan3d=NULL;
  plan->plan=NULL;
#else // MPI_PARALLEL
  if (dir == AthenaFFTForward)
    plan->plan = fftw_plan_dft_2d(nslow,nfast,data,data,FFTW_FORWARD,FFTW_MEASURE);
  else
    plan->plan = fftw_plan_dft_2d(nslow,nfast,data,data,FFTW_BACKWARD,FFTW_MEASURE);
#endif
#endif // FFT

  return plan;
}

//----------------------------------------------------------------------------------------
//! \fn AthenaFFTPlan *FFTBlock::CreatePlan(int nfast, int nmid, int nslow,
//                                           AthenaFFTComplex *data,
//                                           enum AthenaFFTDirection dir)
//  \brief initialize FFT plan for 3D FFT

AthenaFFTPlan *FFTBlock::CreatePlan(int nfast, int nmid, int nslow,
                                     AthenaFFTComplex *data,
                                     enum AthenaFFTDirection dir) {
  AthenaFFTPlan *plan = NULL;

#ifdef FFT
  plan = new AthenaFFTPlan;
  plan->dir = dir;
  plan->dim = dim_;
#ifdef MPI_PARALLEL
  int nbuf;
  int ois[3], oie[3];
  if (dir == AthenaFFTForward) {
    for (int l=0; l<dim_; l++) {
      ois[l]=f_out_->is[(l+(dim_-permute1_)) % dim_];
      oie[l]=f_out_->ie[(l+(dim_-permute1_)) % dim_];
    }
    plan->dir = FFTW_FORWARD;
    plan->plan3d = fft_3d_create_plan(MPI_COMM_WORLD,nfast,nmid,nslow,
                                      f_in_->is[0],f_in_->ie[0],
                                      f_in_->is[1],f_in_->ie[1],
                                      f_in_->is[2],f_in_->ie[2],
                                      ois[0],oie[0],
                                      ois[1],oie[1],
                                      ois[2],oie[2],
                                      0, permute1_, &nbuf);
  } else {
    for (int l=0; l<dim_; l++) {
      ois[l]=b_out_->is[(l+(dim_-permute2_)) % dim_];
      oie[l]=b_out_->ie[(l+(dim_-permute2_)) % dim_];
    }
    plan->dir = FFTW_BACKWARD;
    plan->plan3d = fft_3d_create_plan(MPI_COMM_WORLD,nfast,nmid,nslow,
                                      b_in_->is[0],b_in_->ie[0],
                                      b_in_->is[1],b_in_->ie[1],
                                      b_in_->is[2],b_in_->ie[2],
                                      ois[0],oie[0],
                                      ois[1],oie[1],
                                      ois[2],oie[2],
                                      0, permute2_, &nbuf);
  }
  plan->plan2d=NULL;
  plan->plan=NULL;
#else // MPI_PARALLEL
  if (dir == AthenaFFTForward) {
    plan->plan = fftw_plan_dft_3d(nslow,nmid,nfast,data,data,FFTW_FORWARD,FFTW_MEASURE);
  } else {
    plan->plan = fftw_plan_dft_3d(nslow,nmid,nfast,data,data,FFTW_BACKWARD,FFTW_MEASURE);
  }
#endif
#endif // FFT

  return plan;
}

//----------------------------------------------------------------------------------------
//! \fn void FFTBlock::Execute(AthenaFFTPlan *plan)
//  \brief excute FFT using private vars

void FFTBlock::Execute(AthenaFFTPlan *plan) {
#ifdef FFT
#ifdef MPI_PARALLEL
    if (plan->dim == 3) fft_3d(in_, out_, plan->dir, plan->plan3d);
    if (plan->dim == 2) fft_2d(in_, out_, plan->dir, plan->plan2d);
#else
    fftw_execute_dft(plan->plan, in_, out_);
#endif
#endif // FFT
}

//----------------------------------------------------------------------------------------
//! \fn void FFTBlock::Execute(AthenaFFTPlan *plan, AthenaFFTComplex *data)
//  \brief excute in-place FFT

void FFTBlock::Execute(AthenaFFTPlan *plan, AthenaFFTComplex *data) {
#ifdef FFT
#ifdef MPI_PARALLEL
    if (plan->dim == 3) fft_3d(data, data, plan->dir, plan->plan3d);
    if (plan->dim == 2) fft_2d(data, data, plan->dir, plan->plan2d);
#else
    fftw_execute_dft(plan->plan, data, data);
#endif
#endif // FFT
}

//----------------------------------------------------------------------------------------
//! \fn void FFTBlock::Execute(AthenaFFTPlan *plan,
//                              AthenaFFTComplex *in_data,AthenaFFTComplex *out_data)
//  \brief excute out-place FFT

void FFTBlock::Execute(AthenaFFTPlan *plan, AthenaFFTComplex *in_data,
                        AthenaFFTComplex *out_data) {
#ifdef FFT
#ifdef MPI_PARALLEL
  if (plan->dim == 3) fft_3d(in_data, out_data, plan->dir, plan->plan3d);
  if (plan->dim == 2) fft_2d(in_data, out_data, plan->dir, plan->plan2d);
#else
  fftw_execute_dft(plan->plan, in_data, out_data);
#endif
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void FFTBlock::MpiInitialize(AthenaFFTPlan *plan,
//                              AthenaFFTComplex *in_data,AthenaFFTComplex *out_data)
//  \brief excute out-place FFT

void FFTBlock::MpiInitialize() {
#ifdef MPI_PARALLEL
  std::stringstream msg;
  if ((pdim_ == 2 || pdim_ ==1) && dim_ == 3) {
// To achieve best performance with 2D-pencil decomposition,
// (1) if the "long"-axis (undecomposed-axis) is not the "slow"-axis (x-axis),
//     one needs to permute the axes to make it fast by setting "permute0":
//       yz_decomp (long-axis = x) (i,j,k) --> (i,j,k) permute0=0
//       xz_decomp (long-axis = y) (i,j,k) --> (j,k,i) permute0=1
//       xy_decomp (long-axis = z) (i,j,k) --> (k,i,j) permute0=2
// (2) swap axes for input array (swap1=true)
//     forward FFT with permute1=2 option.
//       yz_decomp (i,j,k) --> (i,k,j) --> (j,i,k)
//       xz_decomp (j,k,i) --> (j,i,k) --> (k,j,i)
//       xy_decomp (k,i,j) --> (k,j,i) --> (i,k,j)
// (3) swap axes from the output of forward FFT to prepare backward FFT (swap2==treu).
//     excute backward FFT with permute2=2 option
//       yz_decomp (j,i,k) --> (j,k,i) --> (i,j,k)
//       xz_decomp (k,j,i) --> (k,i,j) --> (j,k,i)
//       xy_decomp (i,k,j) --> (i,j,k) --> (k,i,j)
// (4) final outcome is the same with original input before swapping.
//     assign it back to original Athena array with permutation
//
// swap1=swap2=true; permute1=permute2=2; permute0 depends on the decomposition

    swap1_ = true; swap2_ = true;
    permute1_ = 2; permute2_ = 2;
    if (decomp_ == DecompositionNames::x_decomp) {
      permute0_ = 1;
    } else if (decomp_ == DecompositionNames::y_decomp) {
      permute0_ = 2;
    } else if (decomp_ == DecompositionNames::z_decomp) {
      permute0_ = 0;
    } else if (decomp_ == DecompositionNames::xy_decomp) {
      permute0_ = 2;
    } else if (decomp_ == DecompositionNames::yz_decomp) {
      permute0_ = 0;
    } else if (decomp_ == DecompositionNames::xz_decomp) {
      permute0_ = 1;
    } else {
      msg << "Something wrong with " << pdim_ << "D decomposition!" << std::endl
          << "Current MPI Configuration is "
          << orig_idx_->np[0] << " x " << orig_idx_->np[1]
          << " x " << orig_idx_->np[2] << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
  } else {
    // For 3D block decompsition, simply set indices as in original Athena Array.
    // two additional remapping will be performed to prepare and recover indices.
    swap1_ = false; swap2_ = false;
    permute0_ = 0; permute1_ = 0; permute2_ = 0;
  }

  // permute axes and procs & swap mid <-> slow indices to prepare forward FFT
  f_in_ = new AthenaFFTIndex(orig_idx_);
  f_in_->PermuteAxis(permute0_);
  f_in_->PermuteProc(permute0_);
  if (swap1_) {
    f_in_->SwapAxis(0);
    f_in_->SwapProc(0);
  }
  f_in_->SetLocalIndex();

  // set output indices of forward FFT;
  // keep global mesh size as input,
  // reverse permutation for MPI configurations to get correct indices
  f_out_ = new AthenaFFTIndex(f_in_);
  f_out_->PermuteAxis(permute1_);
  f_out_->SetLocalIndex();

  // prepare backward FFT;
  // now permute fast, mid, and slow axes twice
  b_in_ = new AthenaFFTIndex(f_out_);
  if (swap2_) {
    b_in_->SwapAxis(0);
    b_in_->SwapProc(0);
  }
  b_in_->SetLocalIndex();

  // set output indices of backward FFT;
  // keep global mesh size as input,
  // reverse permutation for MPI configurations to get correct indices
  b_out_ = new AthenaFFTIndex(b_in_);
  b_out_->PermuteAxis(permute2_);
  b_out_->SetLocalIndex();

#endif
}


//---------------------------------------------------------------------------------------
// AthenaFFTIndex class:

AthenaFFTIndex::AthenaFFTIndex(int dim, LogicalLocation loc, RegionSize msize,
                               RegionSize bsize) {
  dim_=dim;
  // loc.lxi are int64_t in general, but w/o AMR, they are unilikely to overflow int32_t
  Lx[0] = msize.x1max-msize.x1min;
  Nx[0] = msize.nx1;
  np[0] = msize.nx1/bsize.nx1;
  ip[0] = static_cast<int>(loc.lx1);
  iloc[0]=0;
  ploc[0]=0;
  Lx[1] = msize.x2max-msize.x2min;
  Nx[1] = msize.nx2;
  np[1] = msize.nx2/bsize.nx2;
  ip[1] = static_cast<int>(loc.lx2);
  iloc[1]=1;
  ploc[1]=1;
  Lx[2] = msize.x3max-msize.x3min;
  Nx[2] = msize.nx3;
  np[2] = msize.nx3/bsize.nx3;
  ip[2] = static_cast<int>(loc.lx3);
  iloc[2]=2;
  ploc[2]=2;

  SetLocalIndex();
}

// copy constructor
AthenaFFTIndex::AthenaFFTIndex(const AthenaFFTIndex *psrc) {
  dim_ = psrc->dim_;

  for (int i=0; i<3; i++) {
    Lx[i]=psrc->Lx[i];
    Nx[i]=psrc->Nx[i];
    np[i]=psrc->np[i];
    ip[i]=psrc->ip[i];
    iloc[i] = psrc->iloc[i];
    ploc[i] = psrc->ploc[i];
  }

  SetLocalIndex();
}

AthenaFFTIndex::~AthenaFFTIndex() {

}

void AthenaFFTIndex::SetLocalIndex() {
  for (int i=0; i<3; i++) {
    nx[i] = Nx[i]/np[i];
    is[i] = ip[i]*nx[i];
    ie[i] = is[i]+nx[i]-1;
  }
  //  PrintIndex();
}

void AthenaFFTIndex::Swap_(int loc[],int ref_axis) {
  int tmp;
  int axis1=(ref_axis+1) % dim_, axis2=ref_axis+2 % dim_;
  tmp=loc[axis1];
  loc[axis1]=loc[axis2];
  loc[axis2]=tmp;
}

void AthenaFFTIndex::SwapAxis(int ref_axis) {
  Swap_(iloc,ref_axis);
  Swap_(Nx,ref_axis);
}

void AthenaFFTIndex::SwapProc(int ref_axis) {
  Swap_(ploc,ref_axis);
  Swap_(ip,ref_axis);
  Swap_(np,ref_axis);
}

void AthenaFFTIndex::Permute_(int loc[], int npermute) {
  int tmp;
  for (int i=0; i<npermute; i++) {
    tmp=loc[0];
    loc[0]=loc[1];
    loc[1]=loc[2];
    loc[2]=tmp;
  }
}

void AthenaFFTIndex::PermuteAxis(int npermute) {
  Permute_(iloc,npermute);
  Permute_(Nx,npermute);
}

void AthenaFFTIndex::PermuteProc(int npermute) {
  Permute_(ploc,npermute);
  Permute_(np,npermute);
  Permute_(ip,npermute);
}

void AthenaFFTIndex::RemapArray_(int arr[], int loc[], int dir) {
  int tmp[3];
  for (int i=0; i<dim_; i++) tmp[i]=arr[i];
  for (int i=0; i<dim_; i++) {
    if (dir>0) arr[loc[i]]=tmp[i];
    else arr[i]=tmp[loc[i]];
  }
}

void AthenaFFTIndex::RemapAxis(int dir) {
  RemapArray_(Nx,iloc,dir);
}

void AthenaFFTIndex::RemapProc(int dir) {
  RemapArray_(np,ploc,dir);
  RemapArray_(ip,ploc,dir);
}

void AthenaFFTIndex::PrintIndex(void) {
  std::cout << "Lx:" << Lx[0] << " "  << Lx[1] << " " << Lx[2] << std::endl
            << "Nx:" << Nx[0] << " "  << Nx[1] << " " << Nx[2] << std::endl
            << "np:" << np[0] << " "  << np[1] << " " << np[2] << std::endl
            << "ip:" << ip[0] << " "  << ip[1] << " " << ip[2] << std::endl
            << "nx:" << nx[0] << " "  << nx[1] << " " << nx[2] << std::endl
            << "is:" << is[0] << " "  << is[1] << " " << is[2] << std::endl
            << "ie:" << ie[0] << " "  << ie[1] << " " << ie[2] << std::endl
            << "iloc:" << iloc[0] << " "  << iloc[1] << " " << iloc[2] << std::endl
            << "ploc:" << ploc[0] << " "  << ploc[1] << " " << ploc[2] << std::endl;
}
