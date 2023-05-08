//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file dc.cpp
//! \brief piecewise constant (donor cell) reconstruction

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "reconstruction.hpp"

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::DonorCellX1(const int k, const int j, const int il, const int iu,
//!                              const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
//!                              AthenaArray<Real> &wl, AthenaArray<Real> &wr)
//! \brief reconstruct L/R surfaces of the i-th cells

void Reconstruction::DonorCellX1(const int k, const int j, const int il, const int iu,
                                 const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
                                 AthenaArray<Real> &wl, AthenaArray<Real> &wr) {
  // compute L/R states for each variable
  for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      wl(n,i+1) =  wr(n,i) = w(n,k,j,i);
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      wl(IBY,i+1) = wr(IBY,i) = bcc(IB2,k,j,i);
    }
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      wl(IBZ,i+1) = wr(IBZ,i) = bcc(IB3,k,j,i);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::DonorCellX2(const int k, const int j, const int il, const int iu,
//!                              const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
//!                              AthenaArray<Real> &wl, AthenaArray<Real> &wr)
//! \brief


void Reconstruction::DonorCellX2(const int k, const int j, const int il, const int iu,
                                 const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
                                 AthenaArray<Real> &wl, AthenaArray<Real> &wr) {
  // compute L/R states for each variable
  for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      wl(n,i) = wr(n,i) = w(n,k,j,i);
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      wl(IBY,i) = wr(IBY,i) = bcc(IB3,k,j,i);
    }
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      wl(IBZ,i) = wr(IBZ,i) = bcc(IB1,k,j,i);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::DonorCellX3(const int k, const int j, const int il, const int iu,
//!                              const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
//!                              AthenaArray<Real> &wl, AthenaArray<Real> &wr)
//! \brief

void Reconstruction::DonorCellX3(const int k, const int j, const int il, const int iu,
                                 const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
                                 AthenaArray<Real> &wl, AthenaArray<Real> &wr) {
  // compute L/R states for each variable
  for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      wl(n,i) = wr(n,i) = w(n,k,j,i);
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      wl(IBY,i) = wr(IBY,i) = bcc(IB1,k,j,i);
    }
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      wl(IBZ,i) = wr(IBZ,i) = bcc(IB2,k,j,i);
    }
  }
  return;
}
