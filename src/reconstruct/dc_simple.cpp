//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file dc_simple.cpp
//! \brief piecewise constant (donor cell) reconstruction
//! Operates on the entire nx4 range of a single AthenaArray<Real> input (no MHD).
//! No assumptions of hydrodynamic fluid variable input; no characteristic projection.

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "reconstruction.hpp"

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::DonorCellX1(const int k, const int j,
//!                              const int il, const int iu,
//!                              const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
//!                              AthenaArray<Real> &wl, AthenaArray<Real> &wr)
//! \brief reconstruct L/R surfaces of the i-th cells

void Reconstruction::DonorCellX1(const int k, const int j, const int il, const int iu,
                                 const AthenaArray<Real> &q,
                                 AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  const int nu = q.GetDim4() - 1;

  // compute L/R states for each variable
  for (int n=0; n<=nu; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      ql(n,i+1) =  qr(n,i) = q(n,k,j,i);
    }
  }
  return;
}

// function for arrays with different order
void Reconstruction::DonorCellX1(const int k, const int j, const int il, const int iu,
                                 AthenaArray<Real> &q, const int array_order,
                                 AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  if (array_order < 0) {
    const int nu = q.GetDim1() - 1;
    // compute L/R states for each variable
    for (int i=il; i<=iu; ++i) {
      Real *qn = &(q(k,j,i,0));
      Real *qln = &(ql(i+1,0));
      Real *qrn = &(qr(i,0));
      for (int n=0; n<=nu; ++n) {
        qln[n] =  qrn[n] = qn[n];
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::DonorCellX2(const int k, const int j,
//!                              const int il, const int iu,
//!                              const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
//!                              AthenaArray<Real> &wl, AthenaArray<Real> &wr)
//! \brief

void Reconstruction::DonorCellX2(const int k, const int j, const int il, const int iu,
                                 const AthenaArray<Real> &q,
                                 AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  const int nu = q.GetDim4() - 1;
  // compute L/R states for each variable
  for (int n=0; n<=nu; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      ql(n,i) = qr(n,i) = q(n,k,j,i);
    }
  }
  return;
}

// function for arrys with different order
void Reconstruction::DonorCellX2(const int k, const int j, const int il, const int iu,
                                 AthenaArray<Real> &q, const int array_order,
                                 AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  if (array_order < 0) {
    const int nu = q.GetDim1() - 1;
    // compute L/R states for each variable
    for (int i=il; i<=iu; ++i) {
      Real *qln = &(ql(i,0));
      Real *qrn = &(qr(i,0));
      Real *qn = &(q(k,j,i,0));
      for (int n=0; n<=nu; ++n) {
      qln[n] = qrn[n] = qn[n];
    }
  }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::DonorCellX3(const int k, const int j,
//!                              const int il, const int iu,
//!                              const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
//!                              AthenaArray<Real> &wl, AthenaArray<Real> &wr)
//! \brief

void Reconstruction::DonorCellX3(const int k, const int j, const int il, const int iu,
                                 const AthenaArray<Real> &q,
                                 AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  const int nu = q.GetDim4() - 1;
  // compute L/R states for each variable
  for (int n=0; n<=nu; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      ql(n,i) = qr(n,i) = q(n,k,j,i);
    }
  }
  return;
}



void Reconstruction::DonorCellX3(const int k, const int j, const int il, const int iu,
                                 AthenaArray<Real> &q, const int array_order,
                                 AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  const int nu = q.GetDim1() - 1;
  // compute L/R states for each variable
  for (int i=il; i<=iu; ++i) {
    Real *qln = &(ql(i,0));
    Real *qrn = &(qr(i,0));
    Real *qn = &(q(k,j,i,0));
    for (int n=0; n<=nu; ++n) {
      qln[n] = qrn[n] = qn[n];
    }
  }
  return;
}

void Reconstruction::DonorCellZeta(
    NRRadiation *prad, const int zs, const int ze,
    AthenaArray<Real> &q,
    AthenaArray<Real> &ql, AthenaArray<Real> &qr) {

    Real *qln = &(ql(zs+1));
    Real *qrn = &(qr(zs));
    Real *qn = &(q(zs));
    for (int n=0; n<=ze-zs; ++n) {
      qln[n] =  qrn[n] = qn[n];
    }

  return;
}


void Reconstruction::DonorCellPsi(
    NRRadiation *prad, const int ps, const int pe,
    AthenaArray<Real> &q,
    AthenaArray<Real> &ql, AthenaArray<Real> &qr) {

    Real *qln = &(ql(ps+1));
    Real *qrn = &(qr(ps));
    Real *qn = &(q(ps));
    for (int m=0; m<=pe-ps; ++m) {
      // renamed dw* -> dq* from plm.cpp
      qln[m] =  qrn[m] = qn[m];
    }
  return;
}
