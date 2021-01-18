//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file athena_fft.cpp
//  \brief

// C headers

// C++ headers
#include <complex>
#include <iostream>
#include <sstream>
#include <stdexcept>  // runtime_error

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../globals.hpp"
#include "../mesh/mesh.hpp"
#include "../mesh/meshblock_tree.hpp"
#include "athena_fft.hpp"

// constructor, initializes data structures and parameters

FFTBlock::FFTBlock(FFTDriver *pfd, LogicalLocation iloc, int igid,
                   RegionSize msize, RegionSize bsize) :
    // public:
    dx1((msize.x1max - msize.x1min)/static_cast<Real>(msize.nx1)),
    dx2((msize.x2max - msize.x2min)/static_cast<Real>(msize.nx2)),
    dx3((msize.x3max - msize.x3min)/static_cast<Real>(msize.nx3)),
    // protected:
    pmy_driver_(pfd),
    cnt_(bsize.nx1*bsize.nx2*bsize.nx3), gcnt_(pmy_driver_->gcnt_),
    gid_(igid), fplan_{}, bplan_{},
    norm_factor_(1.0), dim_(pmy_driver_->dim_),
    loc_(iloc), msize_(msize), bsize_(bsize),
    orig_idx_{dim_, loc_, msize_, bsize_} {
  in_ = new std::complex<Real>[cnt_];
  out_ = new std::complex<Real>[cnt_];
#ifdef MPI_PARALLEL
  decomp_ = pmy_driver_->decomp_;
  pdim_ = pmy_driver_->pdim_;
  InitializeMPI();
#else
  f_in_  = new AthenaFFTIndex(&orig_idx_);
  f_out_ = new AthenaFFTIndex(&orig_idx_);
  b_in_  = new AthenaFFTIndex(&orig_idx_);
  b_out_ = new AthenaFFTIndex(&orig_idx_);
#endif

  //  f_in_->PrintIndex();
#ifdef FFT
  for (int i=0; i<3; i++) {
    Nx[f_in_->iloc[i]] = f_in_->Nx[i];
    nx[f_in_->iloc[i]] = f_in_->nx[i];
    disp[f_in_->iloc[i]] = f_in_->is[i];
    kNx[b_in_->iloc[i]] = b_in_->Nx[i];
    knx[b_in_->iloc[i]] = b_in_->nx[i];
    kdisp[b_in_->iloc[i]] = b_in_->is[i];
    dkx[b_in_->iloc[i]] = TWO_PI/b_in_->Lx[i];
  }
#endif
}

// destructor

FFTBlock::~FFTBlock() {
  delete[] in_;
  delete[] out_;
  delete f_in_;
  delete f_out_;
  delete b_in_;
  delete b_out_;
  if (fplan_ != nullptr) DestroyPlan(fplan_);
  if (bplan_ != nullptr) DestroyPlan(bplan_);
#ifdef FFT
  fftw_cleanup();
#endif
}

void FFTBlock::DestroyPlan(AthenaFFTPlan *plan) {
#ifdef FFT
#ifdef MPI_PARALLEL
  if (plan->plan3d != nullptr) fft_3d_destroy_plan(plan->plan3d);
  if (plan->plan2d != nullptr) fft_2d_destroy_plan(plan->plan2d);
#endif
  if (plan->plan != nullptr) fftw_destroy_plan(plan->plan);
  delete plan;
#endif
}

void FFTBlock::PrintSource(int in) {
  std::cout << Nx[0] << "x" << Nx[1] << "x" << Nx[2] << std::endl;

  for (int k=0; k<Nx[2]; ++k) {
    for (int j=0; j<Nx[1]; ++j) {
      for (int i=0; i<Nx[0]; ++i) {
        std::int64_t idx = GetIndex(i, j, k, f_in_);
        if (in == 1) std::cout << std::real(in_[idx]) << " ";
        if (in == -1) std::cout << std::real(out_[idx]) << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}
std::int64_t FFTBlock::GetIndex(const int i, const int j, const int k) {
  return i + nx[0] * ( j + nx[1] * k);
}

std::int64_t FFTBlock::GetGlobalIndex(const int i, const int j, const int k) {
  return i + disp[0] + Nx[0]* ( j + disp[1] + Nx[1]* ( k + disp[2]) );
}

std::int64_t FFTBlock::GetIndex(const int i, const int j, const int k,
                                AthenaFFTIndex *pidx) {
  int old_idx[3] = {i, j, k};
  int new_idx[3];
  new_idx[0] = old_idx[pidx->iloc[0]];
  new_idx[1] = old_idx[pidx->iloc[1]];
  new_idx[2] = old_idx[pidx->iloc[2]];

  return new_idx[0] + pidx->nx[0] * (new_idx[1] + pidx->nx[1] * (new_idx[2]));
}

//----------------------------------------------------------------------------------------
//! \fn void FFTBlock::RetrieveResult(const AthenaArray<Real> &src, bool nu,
//!                                   int ngh, LogicalLocation loc, RegionSize bsize)
//! \brief Fill the result in the active zone

void FFTBlock::RetrieveResult(AthenaArray<Real> &dst, bool nu, int ngh,
                              LogicalLocation loc, RegionSize bsize) {
  const std::complex<Real> *src = out_;
  int is, ie, js, je, ks, ke;
  // KGF: possible overflow from std::int64_t loc.lx1 to int is, e.g.
  is = static_cast<int>(loc.lx1*bsize.nx1 - loc_.lx1*bsize_.nx1);
  js = static_cast<int>(loc.lx2*bsize.nx2 - loc_.lx2*bsize_.nx2);
  ks = static_cast<int>(loc.lx3*bsize.nx3 - loc_.lx3*bsize_.nx3);
  ie = is + bsize.nx1 - 1;
  je = bsize.nx2 > 1 ? js + bsize.nx2 - 1:js;
  ke = bsize.nx3 > 1 ? ks + bsize.nx3 - 1:ks;
  int jl = bsize.nx2 > 1 ? ngh:0;
  int kl = bsize.nx3 > 1 ? ngh:0;

  for (int n=0; n<=nu; n++) {
    for (int k=kl, mk=ks; mk<=ke; k++, mk++) {
      for (int j=jl, mj=js; mj<=je; j++, mj++) {
        for (int i=ngh, mi=is; mi<=ie; i++, mi++) {
          std::int64_t idx = GetIndex(mi, mj, mk, b_out_);
          if (n == 0) {
            dst(k,j,i) = std::real(src[idx])*norm_factor_;
          } else {
            dst(n,k,j,i) = std::imag(src[idx])*norm_factor_;
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FFTBlock::LoadSource(const AthenaArray<Real> &src, bool nu, int ngh,
//!                               LogicalLocation loc, RegionSize bsize)
//! \brief Fill the source in the active zone

void FFTBlock::LoadSource(const AthenaArray<Real> &src, bool nu, int ngh,
                          LogicalLocation loc, RegionSize bsize) {
  std::complex<Real> *dst = in_;
  int is, ie, js, je, ks, ke;
  // KGF: possible overflow from std::int64_t loc.lx1 to int is, e.g.
  is = static_cast<int>(loc.lx1*bsize.nx1 - loc_.lx1*bsize_.nx1);
  js = static_cast<int>(loc.lx2*bsize.nx2 - loc_.lx2*bsize_.nx2);
  ks = static_cast<int>(loc.lx3*bsize.nx3 - loc_.lx3*bsize_.nx3);
  ie = is + bsize.nx1 - 1;
  je = bsize.nx2 > 1 ? js + bsize.nx2 - 1:js;
  ke = bsize.nx3 > 1 ? ks + bsize.nx3 - 1:ks;
  int jl = bsize.nx2 > 1 ? ngh:0;
  int kl = bsize.nx3 > 1 ? ngh:0;

  for (int n=0; n<=nu; n++) {
    for (int k=kl, mk=ks; mk<=ke; k++, mk++) {
      for (int j=jl, mj=js; mj<=je; j++, mj++) {
        for (int i=ngh, mi=is; mi<=ie; i++, mi++) {
          std::int64_t idx = GetIndex(mi, mj, mk, f_in_);
          if (n == 0) {
            // copy-list initialization (since C++11)
            dst[idx] = {src(n,k,j,i), 0.0};
          } else {
            dst[idx].imag(src(n,k,j,i));
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FFTBlock::ApplyKernel(int mode)
//! \brief Apply kernel

void FFTBlock::ApplyKernel(int mode) {
  for (int k=0; k<knx[2]; k++) {
    for (int j=0; j<knx[1]; j++) {
      for (int i=0; i<knx[0]; i++) {
        std::int64_t idx_in = GetIndex(i,j,k,b_in_);
        std::int64_t idx_out = GetIndex(i,j,k,f_out_);
        in_[idx_in] = out_[idx_out];
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn AthenaFFTPlan *FFTBlock::QuickCreatePlan(std::complex<Real> *data,
//!                                               AthenaFFTDirection dir)
//! \brief initialize FFT plan using mesh information

AthenaFFTPlan *FFTBlock::QuickCreatePlan(std::complex<Real> *data,
                                         AthenaFFTDirection dir) {
  int nfast, nmid, nslow;
  if (dir == AthenaFFTDirection::forward) {
    nfast = f_in_->Nx[0]; nmid = f_in_->Nx[1]; nslow = f_in_->Nx[2];
  } else {
    nfast = b_in_->Nx[0]; nmid = b_in_->Nx[1]; nslow = b_in_->Nx[2];
  }
  if (dim_ == 3) return CreatePlan(nfast, nmid, nslow, data, dir);
  else if (dim_ == 2) return CreatePlan(nfast, nmid, data, dir);
  else  return CreatePlan(nfast, data, dir);
}

//----------------------------------------------------------------------------------------
//! \fn AthenaFFTPlan *FFTBlock::CreatePlan(int nfast, std::complex<Real> *data,
//!                                          AthenaFFTDirection dir)
//! \brief initialize FFT plan for 1D FFT
AthenaFFTPlan *FFTBlock::CreatePlan(int nfast, std::complex<Real> *data,
                                    AthenaFFTDirection dir) {
  AthenaFFTPlan *plan = nullptr;
#ifdef FFT
  plan = new AthenaFFTPlan;
  plan->dir = static_cast<int>(dir);
  plan->dim = dim_;
  if (dir == AthenaFFTDirection::forward)
    plan->plan = fftw_plan_dft_1d(nfast,
                                  reinterpret_cast<fftw_complex *>(data),
                                  reinterpret_cast<fftw_complex *>(data),
                                  FFTW_FORWARD, FFTW_ESTIMATE);
  else
    plan->plan = fftw_plan_dft_1d(nfast,
                                  reinterpret_cast<fftw_complex *>(data),
                                  reinterpret_cast<fftw_complex *>(data),
                                  FFTW_BACKWARD, FFTW_ESTIMATE);
#endif
  return plan;
}

//----------------------------------------------------------------------------------------
//! \fn AthenaFFTPlan *FFTBlock::CreatePlan(int nfast, int nslow,
//!                                          std::complex<Real> *data,
//!                                          AthenaFFTDirection dir)
//! \brief initialize FFT plan for 2D FFT
AthenaFFTPlan *FFTBlock::CreatePlan(int nfast, int nslow,
                                    std::complex<Real> *data,
                                    AthenaFFTDirection dir) {
  AthenaFFTPlan *plan = nullptr;
#ifdef FFT
  plan = new AthenaFFTPlan;
  plan->dir = static_cast<int>(dir);
  plan->dim = dim_;
#ifdef MPI_PARALLEL
  int nbuf;
  if (dir == AthenaFFTDirection::forward) {
    plan->dir = FFTW_FORWARD;
    plan->plan2d = fft_2d_create_plan(MPI_COMM_WORLD, nfast, nslow,
                                      f_in_->is[0], f_in_->ie[0],
                                      f_in_->is[1], f_in_->ie[1],
                                      f_out_->is[f_in_->iloc[0]],
                                      f_out_->ie[f_in_->iloc[0]],
                                      f_out_->is[f_in_->iloc[1]],
                                      f_out_->ie[f_in_->iloc[1]],
                                      0, permute1_, &nbuf);
  } else {
    plan->dir = FFTW_BACKWARD;
    plan->plan2d = fft_2d_create_plan(MPI_COMM_WORLD, nfast, nslow,
                                      b_in_->is[0], b_in_->ie[0],
                                      b_in_->is[1], b_in_->ie[1],
                                      b_out_->is[b_in_->iloc[0]],
                                      b_out_->ie[b_in_->iloc[0]],
                                      b_out_->is[b_in_->iloc[1]],
                                      b_out_->ie[b_in_->iloc[1]],
                                      0, permute2_, &nbuf);
  }
  plan->plan3d = nullptr;
  plan->plan = nullptr;
#else // MPI_PARALLEL
  if (dir == AthenaFFTDirection::forward)
    plan->plan = fftw_plan_dft_2d(nslow, nfast,
                                  reinterpret_cast<fftw_complex *>(data),
                                  reinterpret_cast<fftw_complex *>(data),
                                  FFTW_FORWARD, FFTW_MEASURE);
  else
    plan->plan = fftw_plan_dft_2d(nslow, nfast,
                                  reinterpret_cast<fftw_complex *>(data),
                                  reinterpret_cast<fftw_complex *>(data),
                                  FFTW_BACKWARD, FFTW_MEASURE);
#endif
#endif // FFT

  return plan;
}

//----------------------------------------------------------------------------------------
//! \fn AthenaFFTPlan *FFTBlock::CreatePlan(int nfast, int nmid, int nslow,
//!                                          std::complex<Real> *data,
//!                                          AthenaFFTDirection dir)
//! \brief initialize FFT plan for 3D FFT

AthenaFFTPlan *FFTBlock::CreatePlan(int nfast, int nmid, int nslow,
                                    std::complex<Real> *data,
                                    AthenaFFTDirection dir) {
  AthenaFFTPlan *plan = nullptr;
#ifdef FFT
  plan = new AthenaFFTPlan;
  plan->dir = static_cast<int>(dir);
  plan->dim = dim_;
#ifdef MPI_PARALLEL
  int nbuf;
  int ois[3], oie[3];
  if (dir == AthenaFFTDirection::forward) {
    for (int l=0; l<dim_; l++) {
      ois[l] = f_out_->is[(l+(dim_-permute1_)) % dim_];
      oie[l] = f_out_->ie[(l+(dim_-permute1_)) % dim_];
    }
    plan->dir = FFTW_FORWARD;
    plan->plan3d = fft_3d_create_plan(MPI_COMM_WORLD, nfast, nmid, nslow,
                                      f_in_->is[0], f_in_->ie[0],
                                      f_in_->is[1], f_in_->ie[1],
                                      f_in_->is[2], f_in_->ie[2],
                                      ois[0], oie[0],
                                      ois[1], oie[1],
                                      ois[2], oie[2],
                                      0, permute1_, &nbuf);
  } else {
    for (int l=0; l<dim_; l++) {
      ois[l] = b_out_->is[(l+(dim_-permute2_)) % dim_];
      oie[l] = b_out_->ie[(l+(dim_-permute2_)) % dim_];
    }
    plan->dir = FFTW_BACKWARD;
    plan->plan3d = fft_3d_create_plan(MPI_COMM_WORLD, nfast, nmid, nslow,
                                      b_in_->is[0], b_in_->ie[0],
                                      b_in_->is[1], b_in_->ie[1],
                                      b_in_->is[2], b_in_->ie[2],
                                      ois[0], oie[0],
                                      ois[1], oie[1],
                                      ois[2], oie[2],
                                      0, permute2_, &nbuf);
  }
  plan->plan2d = nullptr;
  plan->plan = nullptr;
#else // MPI_PARALLEL
  if (dir == AthenaFFTDirection::forward) {
    plan->plan = fftw_plan_dft_3d(nslow, nmid, nfast,
                                  reinterpret_cast<fftw_complex *>(data),
                                  reinterpret_cast<fftw_complex *>(data), FFTW_FORWARD,
                                  FFTW_MEASURE);
  } else {
    plan->plan = fftw_plan_dft_3d(nslow, nmid, nfast,
                                  reinterpret_cast<fftw_complex *>(data),
                                  reinterpret_cast<fftw_complex *>(data), FFTW_BACKWARD,
                                  FFTW_MEASURE);
  }
#endif
#endif // FFT

  return plan;
}

//----------------------------------------------------------------------------------------
//! \fn void FFTBlock::Execute(AthenaFFTPlan *plan)
//! \brief execute FFT using private vars

void FFTBlock::Execute(AthenaFFTPlan *plan) {
#ifdef FFT
#ifdef MPI_PARALLEL
  if (plan->dim == 3) fft_3d(reinterpret_cast<fftw_complex *>(in_),
                             reinterpret_cast<fftw_complex *>(out_),
                             plan->dir, plan->plan3d);
  if (plan->dim == 2) fft_2d(reinterpret_cast<fftw_complex *>(in_),
                             reinterpret_cast<fftw_complex *>(out_),
                             plan->dir, plan->plan2d);
#else
  fftw_execute_dft(plan->plan, reinterpret_cast<fftw_complex *>(in_),
                   reinterpret_cast<fftw_complex *>(out_));
#endif
#endif // FFT
}

//----------------------------------------------------------------------------------------
//! \fn void FFTBlock::Execute(AthenaFFTPlan *plan, std::complex<Real> *data)
//! \brief execute in-place FFT

void FFTBlock::Execute(AthenaFFTPlan *plan, std::complex<Real> *data) {
#ifdef FFT
#ifdef MPI_PARALLEL
  if (plan->dim == 3) fft_3d(reinterpret_cast<fftw_complex *>(data),
                             reinterpret_cast<fftw_complex *>(data),
                             plan->dir, plan->plan3d);
  if (plan->dim == 2) fft_2d(reinterpret_cast<fftw_complex *>(data),
                             reinterpret_cast<fftw_complex *>(data),
                             plan->dir, plan->plan2d);
#else
  fftw_execute_dft(plan->plan, reinterpret_cast<fftw_complex *>(data),
                   reinterpret_cast<fftw_complex *>(data));
#endif
#endif // FFT
}

//----------------------------------------------------------------------------------------
//! \fn void FFTBlock::Execute(AthenaFFTPlan *plan,
//!                             std::complex<Real> *in_data,std::complex<Real> *out_data)
//! \brief execute out-of-place FFT

void FFTBlock::Execute(AthenaFFTPlan *plan, std::complex<Real> *in_data,
                       std::complex<Real> *out_data) {
#ifdef FFT
#ifdef MPI_PARALLEL
  if (plan->dim == 3) fft_3d(reinterpret_cast<fftw_complex *>(in_data),
                             reinterpret_cast<fftw_complex *>(out_data),
                             plan->dir, plan->plan3d);
  if (plan->dim == 2) fft_2d(reinterpret_cast<fftw_complex *>(in_data),
                             reinterpret_cast<fftw_complex *>(out_data),
                             plan->dir, plan->plan2d);
#else
  fftw_execute_dft(plan->plan, reinterpret_cast<fftw_complex *>(in_data),
                   reinterpret_cast<fftw_complex *>(out_data));
#endif
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void FFTBlock::InitializeMPI(AthenaFFTPlan *plan,
//!                             std::complex<Real> *in_data,std::complex<Real> *out_data)
//! \brief initialize for paralle FFT
//!
//! \note
//! To achieve best performance with 2D-pencil decomposition,
//! (1) if the "long"-axis (undecomposed-axis) is not the "slow"-axis (x-axis),
//!     one needs to permute the axes to make it fast by setting "permute0":
//!       yz_decomp (long-axis = x) (i,j,k) --> (i,j,k) permute0=0
//!       xz_decomp (long-axis = y) (i,j,k) --> (j,k,i) permute0=1
//!       xy_decomp (long-axis = z) (i,j,k) --> (k,i,j) permute0=2
//! (2) swap axes for input array (swap1=true)
//!     forward FFT with permute1=2 option.
//!       yz_decomp (i,j,k) --> (i,k,j) --> (j,i,k)
//!       xz_decomp (j,k,i) --> (j,i,k) --> (k,j,i)
//!       xy_decomp (k,i,j) --> (k,j,i) --> (i,k,j)
//! (3) swap axes from the output of forward FFT to prepare backward FFT (swap2==treu).
//!     excute backward FFT with permute2=2 option
//!       yz_decomp (j,i,k) --> (j,k,i) --> (i,j,k)
//!       xz_decomp (k,j,i) --> (k,i,j) --> (j,k,i)
//!       xy_decomp (i,k,j) --> (i,j,k) --> (k,i,j)
//! (4) final outcome is the same with original input before swapping.
//!     assign it back to original Athena array with permutation
//!
//! swap1=swap2=true; permute1=permute2=2; permute0 depends on the decomposition

void FFTBlock::InitializeMPI() {
#ifdef MPI_PARALLEL
  std::stringstream msg;
  if ((pdim_ == 2 || pdim_ ==1) && dim_ == 3) {
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
          << orig_idx_.np[0] << " x " << orig_idx_.np[1]
          << " x " << orig_idx_.np[2] << std::endl;
      ATHENA_ERROR(msg);
    }
  } else {
    // For 3D block decompsition, simply set indices as in original Athena Array.
    // two additional remapping will be performed to prepare and recover indices.
    swap1_ = false; swap2_ = false;
    permute0_ = 0; permute1_ = 0; permute2_ = 0;
  }

  // permute axes and procs & swap mid <-> slow indices to prepare forward FFT
  f_in_ = new AthenaFFTIndex(&orig_idx_);
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
  dim_ = dim;
  // KGF: loc.lxi are std::int64_t in general, but w/o AMR, they are unilikely to overflow
  // std::int32_t type limits
  Lx[0] = msize.x1max-msize.x1min;
  Nx[0] = msize.nx1;
  np[0] = msize.nx1/bsize.nx1;
  ip[0] = static_cast<int>(loc.lx1);
  iloc[0] = 0;
  ploc[0] = 0;
  Lx[1] = msize.x2max-msize.x2min;
  Nx[1] = msize.nx2;
  np[1] = msize.nx2/bsize.nx2;
  ip[1] = static_cast<int>(loc.lx2);
  iloc[1] = 1;
  ploc[1] = 1;
  Lx[2] = msize.x3max-msize.x3min;
  Nx[2] = msize.nx3;
  np[2] = msize.nx3/bsize.nx3;
  ip[2] = static_cast<int>(loc.lx3);
  iloc[2] = 2;
  ploc[2] = 2;

  SetLocalIndex();
}

// copy constructor
AthenaFFTIndex::AthenaFFTIndex(const AthenaFFTIndex *psrc) {
  dim_ = psrc->dim_;

  for (int i=0; i<3; i++) {
    Lx[i] = psrc->Lx[i];
    Nx[i] = psrc->Nx[i];
    np[i] = psrc->np[i];
    ip[i] = psrc->ip[i];
    iloc[i] = psrc->iloc[i];
    ploc[i] = psrc->ploc[i];
  }

  SetLocalIndex();
}

void AthenaFFTIndex::SetLocalIndex() {
  for (int i=0; i<3; i++) {
    nx[i] = Nx[i]/np[i];
    is[i] = ip[i]*nx[i];
    ie[i] = is[i]+nx[i]-1;
  }
  //  PrintIndex();
}

template <typename T> void AthenaFFTIndex::Swap_(T loc[], int ref_axis) {
  T tmp;
  int axis1 = (ref_axis+1) % dim_, axis2 = ref_axis+2 % dim_;
  tmp = loc[axis1];
  loc[axis1] = loc[axis2];
  loc[axis2] = tmp;
}

void AthenaFFTIndex::SwapAxis(int ref_axis) {
  Swap_(iloc,ref_axis);
  Swap_(Nx,ref_axis);
  Swap_(Lx,ref_axis);
}

void AthenaFFTIndex::SwapProc(int ref_axis) {
  Swap_(ploc,ref_axis);
  Swap_(ip,ref_axis);
  Swap_(np,ref_axis);
}

template <typename T> void AthenaFFTIndex::Permute_(T loc[], int npermute) {
  T tmp;
  for (int i=0; i<npermute; i++) {
    tmp = loc[0];
    loc[0] = loc[1];
    loc[1] = loc[2];
    loc[2] = tmp;
  }
}

// For safety when linking to other TUs which might use these function templates in the
// future, provide explicit instantiations (not using template argument deduction)
template void AthenaFFTIndex::Swap_<Real>(Real loc[], int ref_axis);
template void AthenaFFTIndex::Swap_<int>(int loc[], int ref_axis);
template void AthenaFFTIndex::Permute_<Real>(Real loc[], int npermute);
template void AthenaFFTIndex::Permute_<int>(int loc[], int npermute);

void AthenaFFTIndex::PermuteAxis(int npermute) {
  Permute_(iloc,npermute);
  Permute_(Nx,npermute);
  Permute_(Lx,npermute);
}

void AthenaFFTIndex::PermuteProc(int npermute) {
  Permute_(ploc,npermute);
  Permute_(np,npermute);
  Permute_(ip,npermute);
}

void AthenaFFTIndex::PrintIndex() {
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
