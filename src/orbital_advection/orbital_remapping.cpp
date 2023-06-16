//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file orbital_remapping.cpp
//! \brief functions for remapping in the orbital direction

// C/C++ headers
#include <algorithm>  // min()
#include <cfloat>     // FLT_MAX
#include <cmath>      // abs
#include <cstring>    // memcpy
#include <iostream>   // cout, endl
#include <sstream>    //
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"

// this class header
#include "./orbital_advection.hpp"

//----------------------------------------------------------------------------------------
//! \fn void OrbitalAdvection::RemapFluxPlm(AthenaArray<Real> &pflux_,
//!                                         const AthenaArray<Real> &pbuf_,
//!                                         const Real eps_, const int osgn_,
//!                                         const int k, const int j, const int il,
//!                                         const int iu, const int shift_)
//! \brief Remap & Calculate flux with plm

void OrbitalAdvection::RemapFluxPlm(AthenaArray<Real> &pflux_,
                                    const AthenaArray<Real> &pbuf_,
                                    const Real eps_, const int osgn_, const int k,
                                    const int j, const int il, const int iu,
                                    const int shift_) {
  const Real odi     = osgn_-0.5;
  const Real coeff   = odi*(1.0-2.0*odi*eps_);
#pragma omp simd simdlen(SIMD_WIDTH)
  for(int i = il; i <= iu; i++) {
    Real dul  = pbuf_(k,j,i+shift_)-pbuf_(k,j,i+shift_-1);
    Real dur  = pbuf_(k,j,i+shift_+1)-pbuf_(k,j,i+shift_);
    Real du2  = dul*dur;
    Real dum  = 2.0*du2/(dul+dur);
    if (du2 <= 0.0) dum = 0.0;
    pflux_(i) = eps_*(pbuf_(k,j,i+shift_)+coeff*dum);
  }
  return;
}

void OrbitalAdvection::RemapFluxPlm(AthenaArray<Real> &pflux_,
                                    const AthenaArray<Real> &pbuf_,
                                    const Real eps_, const int osgn_, const int k,
                                    const int j, const int il, const int iu,
                                    const int nl, const int nu, const int shift_) {
  const Real odi     = osgn_-0.5;
  const Real coeff   = odi*(1.0-2.0*odi*eps_);
  for (int i = il; i <= iu; i++) {
    for (int n=nl; n<=nu; ++n) {
      Real dul  = pbuf_(k,j,i+shift_,n)-pbuf_(k,j,i+shift_-1,n);
      Real dur  = pbuf_(k,j,i+shift_+1,n)-pbuf_(k,j,i+shift_,n);
      Real du2  = dul*dur;
      Real dum  = 2.0*du2/(dul+dur);
      if (du2 <= 0.0) dum = 0.0;
      pflux_(i,n) = eps_*(pbuf_(k,j,i+shift_,n)+coeff*dum);
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void OrbitalAdvection::RemapFluxPpm(AthenaArray<Real> &pflux_,
//!                                         AthenaArray<Real> &pbuf_,
//!                                         const Real eps_, const int osgn_,
//!                                         const int k, const int j, const int il,
//!                                         const int iu, const int shift_)
//! \brief Remap & Calculate flux with ppm

void OrbitalAdvection::RemapFluxPpm(AthenaArray<Real> &pflux_,
                                    AthenaArray<Real> &pbuf_,
                                    const Real eps_, const int osgn_,
                                    const int k, const int j, const int il,
                                    const int iu, const int shift_) {
  const Real C2      = 1.25;
  const Real odi     = 2.0*(osgn_-0.5);
  AthenaArray<Real> &q_m2  = s_src[0], &q_m1  = s_src[1], &q     = s_src[2],
                    &q_p1  = s_src[3], &q_p2  = s_src[4];
  AthenaArray<Real> &dd    = d_src[0], &dd_m1 = d_src[1], &dd_p1 = d_src[2],
                    &dph     = d_src[3] , &dph_p1  = d_src[4],
                    &d2qc_m1 = d_src[5] , &d2qc    = d_src[6],
                    &d2qc_p1 = d_src[7] , &d2qf    = d_src[8],
                    &qplus   = d_src[9] , &qminus = d_src[10],
                    &dqf_plus = d_src[11], &dqf_minus = d_src[12];
  q_m2.ShallowSlice3DToPencil(pbuf_, k, j, shift_-2, iu+1);
  q_m1.ShallowSlice3DToPencil(pbuf_, k, j, shift_-1, iu+1);
     q.ShallowSlice3DToPencil(pbuf_, k, j, shift_  , iu+1);
  q_p1.ShallowSlice3DToPencil(pbuf_, k, j, shift_+1, iu+1);
  q_p2.ShallowSlice3DToPencil(pbuf_, k, j, shift_+2, iu+1);

//--- Step 1. -------------------------------------------------------------------------
#pragma omp simd
  for (int i = il; i <= iu; i++) {
    Real qa = q   (i) -q_m1(i);
    Real qb = q_p1(i) -q   (i);
    dd_m1(i) = 0.5*(qa+q_m1(i)-q_m2(i));
    dd   (i) = 0.5*(qb+qa);
    dd_p1(i) = 0.5*(q_p2(i)-q_p1(i)+qb);

    dph  (i) = 0.5*(q_m1(i)+q(i))+(dd_m1(i)-dd(i))/6.0;
    dph_p1(i)= 0.5*(q(i)+q_p1(i))+(dd(i)-dd_p1(i))/6.0;
  }
//--- Step 2. -----------------------------------------------------------------------
#pragma omp simd simdlen(SIMD_WIDTH)
  for (int i = il; i <= iu; i++) {
    d2qc_m1(i) = q_m2(i)+q   (i)-2.0*q_m1(i);
    d2qc   (i) = q_m1(i)+q_p1(i)-2.0*q   (i);
    d2qc_p1(i) = q   (i)+q_p2(i)-2.0*q_p1(i);
  }
  // i-1/2
#pragma omp simd simdlen(SIMD_WIDTH)
  for (int i = il; i <= iu; i++) {
    Real qa_tmp = dph(i) - q_m1(i); // (CD eq 84a)
    Real qb_tmp = q  (i) - dph (i);     // (CD eq 84b)
    Real qa = 3.0*(q_m1(i) + q(i) - 2.0*dph(i));  // (CD eq 85b)
    Real qb = d2qc_m1(i);    // (CD eq 85a) (no 1/2)
    Real qc = d2qc(i);   // (CD eq 85c) (no 1/2)
    Real qd = 0.0;
    if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)) {
      qd = SIGN(qa)* std::min(C2*std::abs(qb),
                              std::min(C2*std::abs(qc), std::abs(qa)));
    }
    Real dph_tmp = 0.5*(q_m1(i) + q(i)) - qd/6.0;
    if (qa_tmp*qb_tmp < 0.0) {  // Local extrema detected at i-1/2 face
      dph(i) = dph_tmp;
    }
  }
  // i+1/2
#pragma omp simd simdlen(SIMD_WIDTH)
  for (int i=il; i<=iu; ++i) {
    Real qa_tmp = dph_p1(i) - q     (i);       // (CD eq 84a)
    Real qb_tmp = q_p1  (i) - dph_p1(i);   // (CD eq 84b)
    Real qa = 3.0*(q(i) + q_p1(i) - 2.0*dph_p1(i));  // (CD eq 85b)
    Real qb = d2qc(i);            // (CD eq 85a) (no 1/2)
    Real qc = d2qc_p1(i);   // (CD eq 85c) (no 1/2)
    Real qd = 0.0;
    if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)) {
      qd = SIGN(qa)* std::min(C2*std::abs(qb),
                              std::min(C2*std::abs(qc), std::abs(qa)));
    }
    Real dphp1_tmp = 0.5*(q(i) + q_p1(i)) - qd/6.0;
    if (qa_tmp*qb_tmp < 0.0) { // Local extrema detected at k+1/2 face
      dph_p1(i) = dphp1_tmp;
    }
  }
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    d2qf(i) = 6.0*(dph(i) + dph_p1(i) - 2.0*q(i)); // a6 coefficient * -2
  }

  // Cache Riemann states for both non-/uniform limiters
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    qminus(i) = dph(i);
    qplus(i)  = dph_p1(i);
  }

//--- Step 3. --------------------------------------------------------------------
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    dqf_minus(i) = q(i) - qminus(i); // (CS eq 25)
    dqf_plus(i)  = qplus(i) - q(i);
  }
//--- Step 4. -----------------------------------------------------------------------
#pragma omp simd simdlen(SIMD_WIDTH)
  for (int i=il; i<=iu; ++i) {
    Real qa_tmp = dqf_minus(i)*dqf_plus(i);
    Real qb_tmp = (q_p1(i) - q(i))*(q(i) - q_m1(i));

    // Check if extrema is smooth
    Real qa = d2qc_m1(i);
    Real qb = d2qc(i);
    Real qc = d2qc_p1(i);
    Real qd = d2qf(i);
    Real qe = 0.0;
    if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc) && SIGN(qa) == SIGN(qd)) {
      // Extrema is smooth
      qe = SIGN(qd)* std::min(std::min(C2*std::abs(qa),C2*std::abs(qb)),
                              std::min(C2*std::abs(qc),std::abs(qd))); // (CS eq 22)
    }

    // Check if 2nd derivative is close to roundoff error
    qa = std::max(std::abs(q_m1(i)),std::abs(q_m2(i)));
    qb = std::max(std::max(std::abs(q(i)),
                           std::abs(q_p1(i))), std::abs(q_p2(i)));

    Real rho = 0.0;
    if (std::abs(qd) > (1.0e-12)*std::max(qa,qb)) {
      // Limiter is not sensitive to roundoff. Use limited ratio (MC eq 27)
      rho = qe/qd;
    }

    Real tmp_m = q(i) - rho*dqf_minus(i);
    Real tmp_p = q(i) + rho*dqf_plus(i);
    Real tmp2_m = q(i) - 2.0*dqf_plus(i);
    Real tmp2_p = q(i) + 2.0*dqf_minus(i);

    // Check for local extrema
    if (qa_tmp <= 0.0 || qb_tmp <= 0.0 ) {
      // Check if relative change in limited 2nd deriv is > roundoff
      if (rho <= (1.0 - (1.0e-12))) {
        // Limit smooth extrema
        qminus(i) = tmp_m; // (CS eq 23)
        qplus(i) = tmp_p;
      }

      // No extrema detected
    } else {
      // Overshoot k-1/2,R / k,(-) state
      if (std::abs(dqf_minus(i)) >= 2.0*std::abs(dqf_plus(i))) {
        qminus(i) = tmp2_m;
      }
      // Overshoot k+1/2,L / k,(+) state
      if (std::abs(dqf_plus(i)) >= 2.0*std::abs(dqf_minus(i))) {
        qplus(i) = tmp2_p;
      }
    }
  }

//--- Step 5. -----------------------------------------------------------------------
  Real qx = odi*TWO_3RD*eps_;
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    Real dU = qplus(i)-qminus(i);
    Real U6 = 6.0*(q(i)-0.5*(qplus(i)+qminus(i)));
    pflux_(i) = eps_*(osgn_*qplus(i)+(1-osgn_)*qminus(i)
                      -0.5*eps_*(dU-odi*(1.0-qx)*U6));
  }
  return;
}

void OrbitalAdvection::RemapFluxPpm(AthenaArray<Real> &pflux_,
                                    AthenaArray<Real> &pbuf_,
                                    const Real eps_, const int osgn_,
                                    const int k, const int j, const int il,
                                    const int iu, const int nl, const int nu,
                                    const int shift_) {
  // will implement later
  return;
}
