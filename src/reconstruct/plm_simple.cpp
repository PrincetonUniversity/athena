//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file plm_simple.cpp
//! \brief  piecewise linear reconstruction for both uniform and non-uniform meshes
//! Operates on the entire nx4 range of a single AthenaArray<Real> input (no MHD).
//! No assumptions of hydrodynamic fluid variable input; no characteristic projection.
//!
//! REFERENCES:
//! - (Mignone) A. Mignone, "High-order conservative reconstruction schemes for finite
//!   volume methods in cylindrical and spherical coordinates", JCP, 270, 784 (2014)
//========================================================================================

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "reconstruction.hpp"

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::PiecewiseLinearX1(const int k, const int j,
//!                              const int il, const int iu,
//!                              const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
//!                              AthenaArray<Real> &wl, AthenaArray<Real> &wr)
//! \brief

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
      // renamed dw* -> dq* from plm.cpp
      dql(n,i) = (q(n,k,j,i  ) - q(n,k,j,i-1));
      dqr(n,i) = (q(n,k,j,i+1) - q(n,k,j,i  ));
      qc(n,i) = q(n,k,j,i);
    }
  }

  // Apply simplified van Leer (VL) limiter expression for a Cartesian-like coordinate
  // with uniform mesh spacing
  if (uniform[X1DIR] && !curvilinear[X1DIR]) {
    for (int n=0; n<=nu; ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        Real dq2 = dql(n,i)*dqr(n,i);
        dqm(n,i) = 2.0*dq2/(dql(n,i) + dqr(n,i));
        if (dq2 <= 0.0) dqm(n,i) = 0.0;
      }
    }

    // Apply general VL limiter expression w/ the Mignone correction for a Cartesian-like
    // coordinate with nonuniform mesh spacing or for any curvilinear coordinate spacing
  } else {
    for (int n=0; n<=nu; ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        Real dqF =  dqr(n,i)*pco->dx1f(i)/pco->dx1v(i);
        Real dqB =  dql(n,i)*pco->dx1f(i)/pco->dx1v(i-1);
        Real dq2 = dqF*dqB;
        // cf, cb -> 2 (uniform Cartesian mesh / original VL value) w/ vanishing curvature
        // (may not exactly hold for nonuniform meshes, but converges w/ smooth
        // nonuniformity)
        Real cf = pco->dx1v(i  )/(pco->x1f(i+1) - pco->x1v(i)); // (Mignone eq 33)
        Real cb = pco->dx1v(i-1)/(pco->x1v(i  ) - pco->x1f(i));
        // (modified) VL limiter (Mignone eq 37)
        // (dQ^F term from eq 31 pulled into eq 37, then multiply by (dQ^F/dQ^F)^2)
        dqm(n,i) = (dq2*(cf*dqB + cb*dqF)/
                    (SQR(dqB) + SQR(dqF) + dq2*(cf + cb - 2.0)));
        if (dq2 <= 0.0) dqm(n,i) = 0.0; // ---> no concern for divide-by-0 in above line

        // Real v = dqB/dqF;
        // monotoniced central (MC) limiter (Mignone eq 38)
        // (std::min calls should avoid issue if divide-by-zero causes v=Inf)
        //dqm(n,i) = dqF*std::max(0.0, std::min(0.5*(1.0 + v), std::min(cf, cb*v)));
      }
    }
  }

  // compute ql_(i+1/2) and qr_(i-1/2) using limited slopes
  for (int n=0; n<=nu; ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
    for (int i=il; i<=iu; ++i) {
      // Mignone equation 30
      ql(n,i+1) = qc(n,i) + ((pco->x1f(i+1) - pco->x1v(i))/pco->dx1f(i))*dqm(n,i);
      qr(n,i  ) = qc(n,i) - ((pco->x1v(i  ) - pco->x1f(i))/pco->dx1f(i))*dqm(n,i);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::PiecewiseLinearX2(const int k, const int j,
//!                              const int il, const int iu,
//!                              const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
//!                              AthenaArray<Real> &wl, AthenaArray<Real> &wr)
//! \brief

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
      // renamed dw* -> dq* from plm.cpp
      dql(n,i) = (q(n,k,j  ,i) - q(n,k,j-1,i));
      dqr(n,i) = (q(n,k,j+1,i) - q(n,k,j  ,i));
      qc(n,i) = q(n,k,j,i);
    }
  }

  // Apply simplified van Leer (VL) limiter expression for a Cartesian-like coordinate
  // with uniform mesh spacing
  if (uniform[X2DIR] && !curvilinear[X2DIR]) {
    for (int n=0; n<=nu; ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        Real dq2 = dql(n,i)*dqr(n,i);
        dqm(n,i) = 2.0*dq2/(dql(n,i) + dqr(n,i));
        if (dq2 <= 0.0) dqm(n,i) = 0.0;
      }
    }

    // Apply general VL limiter expression w/ the Mignone correction for a Cartesian-like
    // coordinate with nonuniform mesh spacing or for any curvilinear coordinate spacing
  } else {
    Real cf = pco->dx2v(j  )/(pco->x2f(j+1) - pco->x2v(j));
    Real cb = pco->dx2v(j-1)/(pco->x2v(j  ) - pco->x2f(j));
    Real dxF = pco->dx2f(j)/pco->dx2v(j); // dimensionless, not technically a dx quantity
    Real dxB = pco->dx2f(j)/pco->dx2v(j-1);
    for (int n=0; n<=nu; ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        Real dqF =  dqr(n,i)*dxF;
        Real dqB =  dql(n,i)*dxB;
        Real dq2 = dqF*dqB;
        // (modified) VL limiter (Mignone eq 37)
        dqm(n,i) = (dq2*(cf*dqB + cb*dqF)/
                    (SQR(dqB) + SQR(dqF) + dq2*(cf + cb - 2.0)));
        if (dq2 <= 0.0) dqm(n,i) = 0.0; // ---> no concern for divide-by-0 in above line

        // Real v = dqB/dqF;
        // // monotoniced central (MC) limiter (Mignone eq 38)
        // // (std::min calls should avoid issue if divide-by-zero causes v=Inf)
        // dqm(n,i) = dqF*std::max(0.0, std::min(0.5*(1.0 + v), std::min(cf, cb*v)));
      }
    }
  }

  // compute ql_(j+1/2) and qr_(j-1/2) using limited slopes
  // dimensionless, not technically a "dx" quantity
  Real dxp = (pco->x2f(j+1) - pco->x2v(j))/pco->dx2f(j);
  Real dxm = (pco->x2v(j  ) - pco->x2f(j))/pco->dx2f(j);
  for (int n=0; n<=nu; ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
    for (int i=il; i<=iu; ++i) {
      ql(n,i) = qc(n,i) + dxp*dqm(n,i);
      qr(n,i) = qc(n,i) - dxm*dqm(n,i);
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::PiecewiseLinearX3(const int k, const int j,
//!                              const int il, const int iu,
//!                              const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
//!                              AthenaArray<Real> &wl, AthenaArray<Real> &wr)
//! \brief

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
      // renamed dw* -> dq* from plm.cpp
      dql(n,i) = (q(n,k  ,j,i) - q(n,k-1,j,i));
      dqr(n,i) = (q(n,k+1,j,i) - q(n,k  ,j,i));
      qc(n,i) = q(n,k,j,i);
    }
  }

  // Apply simplified van Leer (VL) limiter expression for a Cartesian-like coordinate
  // with uniform mesh spacing
  if (uniform[X3DIR]) {
    for (int n=0; n<=nu; ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        Real dq2 = dql(n,i)*dqr(n,i);
        dqm(n,i) = 2.0*dq2/(dql(n,i) + dqr(n,i));
        if (dq2 <= 0.0) dqm(n,i) = 0.0;
      }
    }

    // Apply original VL limiter's general expression for a Cartesian-like coordinate with
    // nonuniform mesh spacing
  } else {
    Real dxF = pco->dx3f(k)/pco->dx3v(k);
    Real dxB = pco->dx3f(k)/pco->dx3v(k-1);
    for (int n=0; n<=nu; ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        Real dqF =  dqr(n,i)*dxF;
        Real dqB =  dql(n,i)*dxB;
        Real dq2 = dqF*dqB;
        // original VL limiter (Mignone eq 36)
        dqm(n,i) = 2.0*dq2/(dqF + dqB);
        // dq2 > 0 ---> dqF, dqB are nonzero and have the same sign ----> no risk for
        // (dqF + dqB) = 0 cancellation causing a divide-by-0 in the above line
        if (dq2 <= 0.0) dqm(n,i) = 0.0;
      }
    }
  }

  // compute ql_(k+1/2) and qr_(k-1/2) using limited slopes
  Real dxp = (pco->x3f(k+1) - pco->x3v(k))/pco->dx3f(k);
  Real dxm = (pco->x3v(k  ) - pco->x3f(k))/pco->dx3f(k);
  for (int n=0; n<=nu; ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
    for (int i=il; i<=iu; ++i) {
      ql(n,i) = qc(n,i) + dxp*dqm(n,i);
      qr(n,i) = qc(n,i) - dxm*dqm(n,i);
    }
  }
  return;
}
