//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file plm_simple.cpp
//  \brief  piecewise linear reconstruction for both uniform and non-uniform meshes
//  Operates on the entire nx4 range of a single AthenaArray<Real> input (no MHD).
//  No assumptions of hydrodynamic fluid variable input; no characteristic projection.

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "reconstruction.hpp"

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::PiecewiseLinearX1()
//  \brief

void Reconstruction::PiecewiseLinearX1(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &q,
    AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  Coordinates *pco = pmy_block_->pcoord;
  // set work arrays to shallow copies of scratch arrays
  AthenaArray<Real> &qc = scr1_ni_, &dql = scr2_ni_, &dqr = scr3_ni_,
                   &dqm = scr4_ni_;
  const int nu = q.GetDim4() - 1;

  // compute L/R slopes for each variable
  for (int n=0; n<=nu; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      dql(n,i) = (q(n,k,j,i  ) - q(n,k,j,i-1));
      dqr(n,i) = (q(n,k,j,i+1) - q(n,k,j,i  ));
      qc(n,i) = q(n,k,j,i);
    }
  }

  // Apply van Leer limiter for uniform grid
  if (uniform_limiter[X1DIR]) {
    for (int n=0; n<=nu; ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        Real dq2 = dql(n,i)*dqr(n,i);
        dqm(n,i) = 2.0*dq2/(dql(n,i) + dqr(n,i));
        if (dq2 <= 0.0) dqm(n,i) = 0.0;
      }
    }

    // Apply Mignone limiter for non-uniform grid
  } else {
    for (int n=0; n<=nu; ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        Real dq2 = dql(n,i)*dqr(n,i);
        Real cf = pco->dx1v(i  )/(pco->x1f(i+1) - pco->x1v(i));
        Real cb = pco->dx1v(i-1)/(pco->x1v(i  ) - pco->x1f(i));
        dqm(n,i) = (dq2*(cf*dql(n,i) + cb*dqr(n,i))/
                    (SQR(dql(n,i)) + SQR(dqr(n,i)) + dq2*(cf + cb - 2.0)));
        if (dq2 <= 0.0) dqm(n,i) = 0.0;
      }
    }
  }

  // compute ql_(i+1/2) and qr_(i-1/2) using monotonized slopes
  for (int n=0; n<=nu; ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
    for (int i=il; i<=iu; ++i) {
      ql(n,i+1) = qc(n,i) + ((pco->x1f(i+1)-pco->x1v(i))/pco->dx1f(i))*dqm(n,i);
      qr(n,i  ) = qc(n,i) - ((pco->x1v(i  )-pco->x1f(i))/pco->dx1f(i))*dqm(n,i);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::PiecewiseLinearX2()
//  \brief

void Reconstruction::PiecewiseLinearX2(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &q,
    AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  Coordinates *pco = pmy_block_->pcoord;
  // set work arrays to shallow copies of scratch arrays
  AthenaArray<Real> &qc = scr1_ni_, &dql = scr2_ni_,
                   &dqr = scr3_ni_, &dqm = scr4_ni_;
  const int nu = q.GetDim4() - 1;

  // compute L/R slopes for each variable
  for (int n=0; n<=nu; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      dql(n,i) = (q(n,k,j  ,i) - q(n,k,j-1,i));
      dqr(n,i) = (q(n,k,j+1,i) - q(n,k,j  ,i));
      qc(n,i) = q(n,k,j,i);
    }
  }

  // Apply van Leer limiter for uniform grid
  if (uniform_limiter[X2DIR]) {
    for (int n=0; n<=nu; ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        Real dq2 = dql(n,i)*dqr(n,i);
        dqm(n,i) = 2.0*dq2/(dql(n,i) + dqr(n,i));
        if (dq2 <= 0.0) dqm(n,i) = 0.0;
      }
    }

    // Apply Mignone limiter for non-uniform grid
  } else {
    Real cf = pco->dx2v(j  )/(pco->x2f(j+1) - pco->x2v(j));
    Real cb = pco->dx2v(j-1)/(pco->x2v(j  ) - pco->x2f(j));
    for (int n=0; n<=nu; ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        Real dq2 = dql(n,i)*dqr(n,i);
        dqm(n,i) = (dq2*(cf*dql(n,i) + cb*dqr(n,i))/
                    (SQR(dql(n,i)) + SQR(dqr(n,i)) + dq2*(cf + cb - 2.0)));
        if (dq2 <= 0.0) dqm(n,i) = 0.0;
      }
    }
  }

  // compute ql_(j+1/2) and qr_(j-1/2) using monotonized slopes
  Real dxp = (pco->x2f(j+1)-pco->x2v(j))/pco->dx2f(j);
  Real dxm = (pco->x2v(j  )-pco->x2f(j))/pco->dx2f(j);
  for (int n=0; n<=nu; ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
    for (int i=il; i<=iu; ++i) {
      ql(n,i) = qc(n,i) + dxp*dqm(n,i);
      qr(n,i) = qc(n,i) - dxm*dqm(n,i);
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::PiecewiseLinearX3()
//  \brief

void Reconstruction::PiecewiseLinearX3(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &q,
    AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  Coordinates *pco = pmy_block_->pcoord;
  // set work arrays to shallow copies of scratch arrays
  AthenaArray<Real> &qc = scr1_ni_, &dql = scr2_ni_, &dqr = scr3_ni_,
                   &dqm = scr4_ni_;
  const int nu = q.GetDim4() - 1;

  // compute L/R slopes for each variable
  for (int n=0; n<=nu; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      dql(n,i) = (q(n,k  ,j,i) - q(n,k-1,j,i));
      dqr(n,i) = (q(n,k+1,j,i) - q(n,k  ,j,i));
      qc(n,i) = q(n,k,j,i);
    }
  }

  // Apply van Leer limiter for uniform grid
  if (uniform_limiter[X3DIR]) {
    for (int n=0; n<=nu; ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        Real dq2 = dql(n,i)*dqr(n,i);
        dqm(n,i) = 2.0*dq2/(dql(n,i) + dqr(n,i));
        if (dq2 <= 0.0) dqm(n,i) = 0.0;
      }
    }

    // Apply Mignone limiter for non-uniform grid
  } else {
    for (int n=0; n<=nu; ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        Real dq2 = dql(n,i)*dqr(n,i);
        Real cf = pco->dx3v(k  )/(pco->x3f(k+1) - pco->x3v(k));
        Real cb = pco->dx3v(k-1)/(pco->x3v(k  ) - pco->x3f(k));
        dqm(n,i) = (dq2*(cf*dql(n,i) + cb*dqr(n,i))/
                    (SQR(dql(n,i)) + SQR(dqr(n,i)) + dq2*(cf + cb - 2.0)));
        if (dq2 <= 0.0) dqm(n,i) = 0.0;
      }
    }
  }

  // compute ql_(k+1/2) and qr_(k-1/2) using monotonized slopes
  Real dxp = (pco->x3f(k+1)-pco->x3v(k))/pco->dx3f(k);
  Real dxm = (pco->x3v(k  )-pco->x3f(k))/pco->dx3f(k);
  for (int n=0; n<=nu; ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
    for (int i=il; i<=iu; ++i) {
      ql(n,i) = qc(n,i) + dxp*dqm(n,i);
      qr(n,i) = qc(n,i) - dxm*dqm(n,i);
    }
  }
  return;
}
