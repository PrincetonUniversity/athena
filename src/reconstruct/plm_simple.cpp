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
#include "../nr_radiation/radiation.hpp"
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



void Reconstruction::PiecewiseLinearX1(
    const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &q, const int array_order,
    AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  Coordinates *pco = pmy_block_->pcoord;
  if (array_order < 0) {
    // set work arrays to shallow copies of scratch arrays
    AthenaArray<Real> &qc = scr1_in2_, &dql = scr2_in2_, &dqr = scr3_in2_,
                     &dqm = scr4_in2_;
    const int nu = q.GetDim1() - 1;

    // compute L/R slopes for each variable
    for (int i=il; i<=iu; ++i) {
      Real *dqln = &(dql(i,0));
      Real *dqrn = &(dqr(i,0));
      Real *qn = &(q(k,j,i,0));
      Real *q1n = &(q(k,j,i+1,0));
      Real *q2n = &(q(k,j,i-1,0));
      Real *qcn = &(qc(i,0));
      for (int n=0; n<=nu; ++n) {
        dqln[n] = (qn[n] - q2n[n]);
        dqrn[n] = (q1n[n] - qn[n]);
        qcn[n] = qn[n];
        // renamed dw* -> dq* from plm.cpp
      }
    }

    // Apply simplified van Leer (VL) limiter expression for a Cartesian-like coordinate
    // with uniform mesh spacing
    if (uniform[X1DIR] && !curvilinear[X1DIR]) {
      for (int i=il; i<=iu; ++i) {
#pragma omp simd simdlen(SIMD_WIDTH)
        for (int n=0; n<=nu; ++n) {
          Real dq2 = dql(i,n)*dqr(i,n);
          dqm(i,n) = 2.0*dq2/(dql(i,n) + dqr(i,n));
          if (dq2 <= 0.0) dqm(i,n) = 0.0;
        }
      }
      // coordinate with nonuniform mesh spacing or for any curvilinear coordinate spacing
    } else {
      for (int i=il; i<=iu; ++i) {
        // variables independent n
        // cf, cb -> 2 (uniform Cartesian mesh / original VL value) w/ vanishing curvature
        // (may not exactly hold for nonuniform meshes, but converges w/ smooth
        // nonuniformity)
        Real cf = pco->dx1v(i  )/(pco->x1f(i+1) - pco->x1v(i)); // (Mignone eq 33)
        Real cb = pco->dx1v(i-1)/(pco->x1v(i  ) - pco->x1f(i));
        Real dxF = pco->dx1f(i)/pco->dx1v(i);
        Real dxB = pco->dx1f(i)/pco->dx1v(i-1);
        Real *dqrn = &(dqr(i,0));
        Real *dqln = &(dql(i,0));
        Real *dqmn = &(dqm(i,0));
        for (int n=0; n<=nu; ++n) {
          Real dqF =  dqrn[n]*dxF;
          Real dqB =  dqln[n]*dxB;
          Real dq2 = dqF*dqB;
          // (modified) VL limiter (Mignone eq 37)
          // (dQ^F term from eq 31 pulled into eq 37, then multiply by (dQ^F/dQ^F)^2)
          dqmn[n] = (dq2*(cf*dqB + cb*dqF)/
                     (SQR(dqB) + SQR(dqF) + dq2*(cf + cb - 2.0)));
          if (dq2 <= 0.0) dqmn[n] = 0.0; // ---> no concern for divide-by-0 in above line
          // Real v = dqB/dqF;
          // monotoniced central (MC) limiter (Mignone eq 38)
          // (std::min calls should avoid issue if divide-by-zero causes v=Inf)
          //dqm(n,i) = dqF*std::max(0.0, std::min(0.5*(1.0 + v), std::min(cf, cb*v)));
        }
      }
    }

    // compute ql_(i+1/2) and qr_(i-1/2) using limited slopes
    for (int i=il; i<=iu; ++i) {
      Real ratio_l=(pco->x1f(i+1) - pco->x1v(i))/pco->dx1f(i);
      Real ratio_r=(pco->x1v(i  ) - pco->x1f(i))/pco->dx1f(i);
      Real *qln = &(ql(i+1,0));
      Real *qrn = &(qr(i,0));
      Real *qcn = &(qc(i,0));
      Real *dqmn = &(dqm(i,0));
      for (int n=0; n<=nu; ++n) {
        // Mignone equation 30
        qln[n] = qcn[n] + ratio_l*dqmn[n];
        qrn[n] = qcn[n] - ratio_r*dqmn[n];
      }
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


void Reconstruction::PiecewiseLinearX2(
    const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &q, const int array_order,
    AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  Coordinates *pco = pmy_block_->pcoord;
  if (array_order < 0) {
    // set work arrays to shallow copies of scratch arrays
    AthenaArray<Real> &qc = scr1_in2_, &dql = scr2_in2_,
                     &dqr = scr3_in2_, &dqm = scr4_in2_;
    const int nu = q.GetDim1() - 1;

    // compute L/R slopes for each variable
    for (int i=il; i<=iu; ++i) {
      Real *dqln = &(dql(i,0));
      Real *dqrn = &(dqr(i,0));
      Real *qcn = &(qc(i,0));
      Real *qn  = &(q(k,j  ,i,0));
      Real *q1n = &(q(k,j+1,i,0));
      Real *q2n = &(q(k,j-1,i,0));
      for (int n=0; n<=nu; ++n) {
        // renamed dw* -> dq* from plm.cpp
        dqln[n] = (qn[n] - q2n[n]);
        dqrn[n] = (q1n[n] - qn[n]);
        qcn[n] = qn[n];
      }
    }

    // Apply simplified van Leer (VL) limiter expression for a Cartesian-like coordinate
    // with uniform mesh spacing
    if (uniform[X2DIR] && !curvilinear[X2DIR]) {
      for (int i=il; i<=iu; ++i) {
        Real *dqln = &(dql(i,0));
        Real *dqrn = &(dqr(i,0));
        Real *dqmn = &(dqm(i,0));
        for (int n=0; n<=nu; ++n) {
          Real dq2 = dqln[n]*dqrn[n];
          dqmn[n] = 2.0*dq2/(dqln[n] + dqrn[n]);
          if (dq2 <= 0.0) dqmn[n] = 0.0;
        }
      }

      // Apply general VL limiter expression w/ the Mignone correction for a
      // Cartesian-like coordinate sys with nonuniform mesh spacing or for
      // any curvilinear coordinate spacing
    } else {
      Real cf = pco->dx2v(j  )/(pco->x2f(j+1) - pco->x2v(j));
      Real cb = pco->dx2v(j-1)/(pco->x2v(j  ) - pco->x2f(j));
      Real dxF = pco->dx2f(j)/pco->dx2v(j); // dimensionless, not technically a dx qty
      Real dxB = pco->dx2f(j)/pco->dx2v(j-1);

      for (int i=il; i<=iu; ++i) {
        Real *dqrn = &(dqr(i,0));
        Real *dqln = &(dql(i,0));
        Real *dqmn = &(dqm(i,0));
        for (int n=0; n<=nu; ++n) {
          Real dqF =  dqrn[n]*dxF;
          Real dqB =  dqln[n]*dxB;
          Real dq2 = dqF*dqB;
          // (modified) VL limiter (Mignone eq 37)
          dqmn[n] = (dq2*(cf*dqB + cb*dqF)/
                     (SQR(dqB) + SQR(dqF) + dq2*(cf + cb - 2.0)));
          if (dq2 <= 0.0) dqmn[n] = 0.0; // ---> no concern for divide-by-0 in above line

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
    for (int i=il; i<=iu; ++i) {
      Real *qln = &(ql(i,0));
      Real *qrn = &(qr(i,0));
      Real *dqmn = &(dqm(i,0));
      Real *qcn = &(qc(i,0));
      for (int n=0; n<=nu; ++n) {
        qln[n] = qcn[n] + dxp*dqmn[n];
        qrn[n] = qcn[n] - dxm*dqmn[n];
      }
    }
  }// End array_order
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



void Reconstruction::PiecewiseLinearX3(
    const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &q, const int array_order,
    AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  Coordinates *pco = pmy_block_->pcoord;
  if (array_order < 0) {
    // set work arrays to shallow copies of scratch arrays
    AthenaArray<Real> &qc = scr1_in2_, &dql = scr2_in2_, &dqr = scr3_in2_,
                     &dqm = scr4_in2_;
    const int nu = q.GetDim1() - 1;

    // compute L/R slopes for each variable
    for (int i=il; i<=iu; ++i) {
      Real *dqln = &(dql(i,0));
      Real *dqrn = &(dqr(i,0));
      Real *qcn = &(qc(i,0));
      Real *qn  = &(q(k,j  ,i,0));
      Real *q1n = &(q(k+1,j,i,0));
      Real *q2n = &(q(k-1,j,i,0));
      for (int n=0; n<=nu; ++n) {
        // renamed dw* -> dq* from plm.cpp
        dqln[n] = (qn[n] - q2n[n]);
        dqrn[n] = (q1n[n] - qn[n]);
        qcn[n] = qn[n];
      }
    }

    // Apply simplified van Leer (VL) limiter expression for a Cartesian-like coordinate
    // with uniform mesh spacing
    if (uniform[X3DIR]) {
      for (int i=il; i<=iu; ++i) {
        Real *dqln = &(dql(i,0));
        Real *dqrn = &(dqr(i,0));
        Real *dqmn = &(dqm(i,0));
        for (int n=0; n<=nu; ++n) {
          Real dq2 = dqln[n]*dqrn[n];
          dqmn[n] = 2.0*dq2/(dqln[n] + dqrn[n]);
          if (dq2 <= 0.0) dqmn[n] = 0.0;
        }
      }

      // Apply original VL limiter's general expression for a Cartesian-like
      // coordinate system with nonuniform mesh spacing
    } else {
      Real dxF = pco->dx3f(k)/pco->dx3v(k);
      Real dxB = pco->dx3f(k)/pco->dx3v(k-1);

      for (int i=il; i<=iu; ++i) {
        Real *dqmn = &(dqm(i,0));
        Real *dqrn = &(dqr(i,0));
        Real *dqln = &(dql(i,0));
        for (int n=0; n<=nu; ++n) {
          Real dqF =  dqrn[n]*dxF;
          Real dqB =  dqln[n]*dxB;
          Real dq2 = dqF*dqB;
          // original VL limiter (Mignone eq 36)
          dqmn[n] = 2.0*dq2/(dqF + dqB);
          // dq2 > 0 ---> dqF, dqB are nonzero and have the same sign ----> no risk for
          // (dqF + dqB) = 0 cancellation causing a divide-by-0 in the above line
          if (dq2 <= 0.0) dqmn[n] = 0.0;
        }
      }
    }

    // compute ql_(k+1/2) and qr_(k-1/2) using limited slopes
    Real dxp = (pco->x3f(k+1) - pco->x3v(k))/pco->dx3f(k);
    Real dxm = (pco->x3v(k  ) - pco->x3f(k))/pco->dx3f(k);

    for (int i=il; i<=iu; ++i) {
      Real *qln = &(ql(i,0));
      Real *qrn = &(qr(i,0));
      Real *qcn = &(qc(i,0));
      Real *dqmn = &(dqm(i,0));
      for (int n=0; n<=nu; ++n) {
        qln[n] = qcn[n] + dxp*dqmn[n];
        qrn[n] = qcn[n] - dxm*dqmn[n];
      }
    }
  }
  return;
}

// reconstruction in the angular space
// q, ql, and qr are 2D array in zeta, psi
// zs, ze are starting/end index for zeta
// ps, pe are starting/end index for psi
// q is always in the order (zeta, psi)
// q, ql and qr can change
void Reconstruction::PiecewiseLinearZeta(
    NRRadiation *prad, const int zs, const int ze,
    AthenaArray<Real> &q,
    AthenaArray<Real> &ql, AthenaArray<Real> &qr) {

  // set work arrays to shallow copies of scratch arrays
  AthenaArray<Real> &qc = scr1_nn_, &dql = scr2_nn_, &dqr = scr3_nn_,
                   &dqm = scr4_nn_;
#pragma omp simd
  for (int n=zs; n<=ze; ++n) {
    // renamed dw* -> dq* from plm.cpp
    dql(n) = (q(n) - q(n-1));
    dqr(n) = (q(n+1) - q(n));
    qc(n) = q(n);
  }
  for (int n=zs; n<=ze; ++n) {
    Real cf = prad->dzeta_v(n)/(prad->zeta_f_full(n+1)
                                - prad->zeta_v_full(n)); // (Mignone eq 33)
    Real cb = prad->dzeta_v(n-1)/(prad->zeta_v_full(n) - prad->zeta_f_full(n));
    Real dqF =  dqr(n)*prad->dzeta_f(n)/prad->dzeta_v(n);
    Real dqB =  dql(n)*prad->dzeta_f(n)/prad->dzeta_v(n-1);
    Real dq2 = dqF*dqB;
    // (modified) VL limiter (Mignone eq 37)
    // (dQ^F term from eq 31 pulled into eq 37, then multiply by (dQ^F/dQ^F)^2)
    dqm(n) = (dq2*(cf*dqB + cb*dqF)/
              (SQR(dqB) + SQR(dqF) + dq2*(cf + cb - 2.0)));
    if (dq2 <= 0.0) dqm(n) = 0.0; // ---> no concern for divide-by-0 in above line
    // Real v = dqB/dqF;
    // monotoniced central (MC) limiter (Mignone eq 38)
    // (std::min calls should avoid issue if divide-by-zero causes v=Inf)
    //dqm(n,i) = dqF*std::max(0.0, std::min(0.5*(1.0 + v), std::min(cf, cb*v)));
  }
  for (int n=zs; n<=ze; ++n) {
    Real ratio_l=(prad->zeta_f_full(n+1) - prad->zeta_v_full(n))/prad->dzeta_f(n);
    Real ratio_r=(prad->zeta_v_full(n) - prad->zeta_f_full(n))/prad->dzeta_f(n);
    // Mignone equation 30
    ql(n+1) = qc(n) + ratio_l*dqm(n);
    qr(n  ) = qc(n) - ratio_r*dqm(n);
  }
  return;
}

void Reconstruction::PiecewiseLinearPsi(
    NRRadiation *prad, const int ps, const int pe,
    AthenaArray<Real> &q,
    AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  // set work arrays to shallow copies of scratch arrays
  AthenaArray<Real> &qc = scr1_nn_, &dql = scr2_nn_, &dqr = scr3_nn_,
                   &dqm = scr4_nn_;
#pragma omp simd
  for (int m=ps; m<=pe; ++m) {
    // renamed dw* -> dq* from plm.cpp
    dql(m) = (q(m) - q(m-1));
    dqr(m) = (q(m+1) - q(m));
    qc(m) = q(m);
  }

#pragma omp simd
  for (int m=ps; m<=pe; ++m) {
    Real cf = prad->dpsi_v(m)/(prad->psi_f_full(m+1)
                               - prad->psi_v_full(m)); // (Mignone eq 33)
    Real cb = prad->dpsi_v(m-1)/(prad->psi_v_full(m) - prad->psi_f_full(m));
    Real dqF =  dqr(m)*prad->dpsi_f(m)/prad->dpsi_v(m);
    Real dqB =  dql(m)*prad->dpsi_f(m)/prad->dpsi_v(m-1);
    Real dq2 = dqF*dqB;
    // (modified) VL limiter (Mignone eq 37)
    // (dQ^F term from eq 31 pulled into eq 37, then multiply by (dQ^F/dQ^F)^2)
    dqm(m) = (dq2*(cf*dqB + cb*dqF)/
              (SQR(dqB) + SQR(dqF) + dq2*(cf + cb - 2.0)));
    if (dq2 <= 0.0) dqm(m) = 0.0; // ---> no concern for divide-by-0 in above line
    // Real v = dqB/dqF;
    // monotoniced central (MC) limiter (Mignone eq 38)
    // (std::min calls should avoid issue if divide-by-zero causes v=Inf)
    //dqm(n,i) = dqF*std::max(0.0, std::min(0.5*(1.0 + v), std::min(cf, cb*v)));
  }

#pragma omp simd
  for (int m=ps; m<=pe; ++m) {
    Real ratio_l=(prad->psi_f_full(m+1) - prad->psi_v_full(m))/prad->dpsi_f(m);
    Real ratio_r=(prad->psi_v_full(m) - prad->psi_f_full(m))/prad->dpsi_f(m);
    // Mignone equation 30
    ql(m+1) = qc(m) + ratio_l*dqm(m);
    qr(m  ) = qc(m) - ratio_r*dqm(m);
  }
  return;
}
