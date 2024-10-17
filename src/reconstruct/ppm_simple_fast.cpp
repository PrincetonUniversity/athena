//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file ppm_simple.cpp
//! \brief piecewise parabolic reconstruction with modified McCorquodale/Colella limiter
//!        for a Cartesian-like coordinate with uniform spacing, Mignone modified original
//!        PPM limiter for nonuniform and/or curvilinear coordinate.
//! Operates on the entire nx4 range of a single AthenaArray<Real> input (no MHD).
//! No assumptions of hydrodynamic fluid variable input; no characteristic projection.
//!
//! REFERENCES:
//! - (CW) P. Colella & P. Woodward, "The Piecewise Parabolic Method (PPM) for Gas-
//!   Dynamical Simulations", JCP, 54, 174 (1984)
//! - (CS) P. Colella & M. Sekora, "A limiter for PPM that preserves accuracy at smooth
//!   extrema", JCP, 227, 7069 (2008)
//! - (MC) P. McCorquodale & P. Colella,  "A high-order finite-volume method for
//!   conservation laws on locally refined grids", CAMCoS, 6, 1 (2011)
//! - (CD) P. Colella, M.R. Dorr, J. Hittinger, D. Martin, "High-order, finite-volume
//!   methods in mapped coordinates", JCP, 230, 2952 (2011)
//! - (Mignone) A. Mignone, "High-order conservative reconstruction schemes for finite
//!   volume methods in cylindrical and spherical coordinates", JCP, 270, 784 (2014)
//========================================================================================

// C headers

// C++ headers
#include <algorithm>    // max()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "reconstruction.hpp"

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::PiecewiseParabolicFastX1(const int k, const int j,
//!                              const int il, const int iu,
//!                              const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
//!                              AthenaArray<Real> &wl, AthenaArray<Real> &wr)
//! \brief Returns L/R interface values in X1-dir constructed using fourth-order PPM and
//!        Colella-Sekora or Mignone limiting over [kl,ku][jl,ju][il,iu]

void Reconstruction::PiecewiseParabolicFastX1(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &q,
    AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  const int nu = q.GetDim4() - 1;

  // CS08 constant used in second derivative limiter, >1 , independent of h
  const Real C2 = 1.25;

  // TODO(felker): renumber scratch array references; not using 2x from ppm.cpp
  // bx (MHD) and wc (characteristic projection)

  // set work arrays used for primitive/characterstic cell-averages to scratch
  AthenaArray<Real> &q_i = scr4_ni_;

  // cache the x1-sliced primitive states for eigensystem calculation
  for (int n=0; n<=nu; ++n) {
#pragma omp simd
    for (int i=il-2; i<=iu+2; ++i) {
      q_i  (n,i) = q(n,k,j,i  );
    }
  }

  //--- Step 1. --------------------------------------------------------------------------
  // Reconstruct interface averages <a>_{i-1/2} and <a>_{i+1/2}
  for (int n=0; n<=nu; ++n) {
    // Compute average slope in i-1, i, i+1 zones
#pragma omp simd simdlen(SIMD_WIDTH)
    for (int i=il; i<=iu; ++i) {
      // nonuniform or uniform Cartesian-like coord reconstruction from volume averages:
      Real dph_i, dph_ip1_i;
      Real q_i0_i  = q_i(n,i);
      Real q_im1_i = q_i(n,i-1);
      Real q_im2_i = q_i(n,i-2);
      Real q_ip1_i = q_i(n,i+1);
      Real q_ip2_i = q_i(n,i+2);

      if (!curvilinear_[X1DIR]) {
        Real qa = (q_i0_i - q_im1_i);
        Real qb = (q_ip1_i - q_i0_i);
        Real dd_im1_i = c1i(i-1)*qa + c2i(i-1)*(q_im1_i - q_im2_i);
        Real dd_i     = c1i(i  )*qb + c2i(i  )*qa;
        Real dd_ip1_i = c1i(i+1)*(q_ip2_i - q_ip1_i) + c2i(i+1)*qb;

        // Approximate interface average at i-1/2 and i+1/2 using PPM (CW eq 1.6)
        // KGF: group the biased stencil quantities to preserve FP symmetry
        dph_i = (c3i(i)*q_im1_i + c4i(i)*q_i0_i) +
                 (c5i(i)*dd_im1_i + c6i(i)*dd_i);
        dph_ip1_i = (c3i(i+1)*q_i0_i + c4i(i+1)*q_ip1_i) +
                     (c5i(i+1)*dd_i + c6i(i+1)*dd_ip1_i );
      } else { // radial coordinate
        dph_i = c1i(i)*q_im2_i + c2i(i)*q_im1_i + c3i(i)*q_i0_i
                 + c4i(i)*q_ip1_i;
        dph_ip1_i = c1i(i+1)*q_im1_i + c2i(i+1)*q_i0_i + c3i(i+1)*q_ip1_i
                     + c4i(i+1)*q_ip2_i;
      }


      //--- Step 2a. --------------------------------------------------------------------
      // Uniform Cartesian-like coordinate:
      // limit interpolated interface states (CD 4.3.1)
      if (uniform_[X1DIR] && !curvilinear_[X1DIR]) {
        // approximate second derivative at interfaces for smooth extrema preservation

        // KGF: add the off-centered quantities first to preserve FP symmetry
        Real d2qc_im1_i = q_im2_i + q_i0_i  - 2.0*q_im1_i;
        // (CD eq 85a) (no 1/2)
        Real d2qc_i     = q_im1_i + q_ip1_i - 2.0*q_i0_i;
        Real d2qc_ip1_i = q_i0_i  + q_ip2_i - 2.0*q_ip1_i;

        // i-1/2
        Real qa_tmp = dph_i - q_im1_i;   // (CD eq 84a)
        Real qb_tmp = q_i0_i - dph_i;     // (CD eq 84b)
        // KGF: add the off-centered quantities first to preserve FP symmetry
        Real qa = 3.0*(q_im1_i + q_i0_i  - 2.0*dph_i);  // (CD eq 85b)
        Real qb = d2qc_im1_i;    // (CD eq 85a) (no 1/2)
        Real qc = d2qc_i;   // (CD eq 85c) (no 1/2)
        Real qd = 0.0;
        if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)) {
          qd = SIGN(qa)* std::min(C2*std::abs(qb),
                                  std::min(C2*std::abs(qc), std::abs(qa)));
        }
        Real dph_tmp = 0.5*(q_im1_i + q_i0_i) - qd/6.0;
        if (qa_tmp*qb_tmp < 0.0) { // Local extrema detected at i-1/2 face
          dph_i = dph_tmp;
        }


        // i+1/2
        qa_tmp = dph_ip1_i - q_i0_i;     // (CD eq 84a)
        qb_tmp = q_ip1_i - dph_ip1_i;   // (CD eq 84b)
        // KGF: add the off-centered quantities first to preserve FP symmetry
        qa = 3.0*(q_i0_i + q_ip1_i - 2.0*dph_ip1_i);  // (CD eq 85b)
        qb = d2qc_i;            // (CD eq 85a) (no 1/2)
        qc = d2qc_ip1_i;   // (CD eq 85c) (no 1/2)
        qd = 0.0;
        if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)) {
          qd = SIGN(qa)* std::min(C2*std::abs(qb),
                                  std::min(C2*std::abs(qc), std::abs(qa)));
        }
        Real dphip1_tmp = 0.5*(q_i0_i + q_ip1_i) - qd/6.0;
        if (qa_tmp*qb_tmp < 0.0) { // Local extrema detected at i+1/2 face
          dph_ip1_i = dphip1_tmp;
        }

      //--- Step 2b. ---------------------------------------------------------------------
      // Nonuniform and/or curvilinear coordinate: apply strict monotonicity constraints
      // (Mignone eq 45)
      } else {
        dph_i     = std::min(dph_i    , std::max(q_i0_i,q_im1_i));
        dph_ip1_i = std::min(dph_ip1_i, std::max(q_i0_i,q_ip1_i));

        dph_i     = std::max(dph_i    , std::min(q_i0_i,q_im1_i));
        dph_ip1_i = std::max(dph_ip1_i, std::min(q_i0_i,q_ip1_i));
      }

      // Cache Riemann states for both non-/uniform limiters
      Real qminus_i = dph_i;
      Real qplus_i  = dph_ip1_i;

      //--- Step 3. ---------------------------------------------------------------------
      // Compute cell-centered difference stencils (MC section 2.4.1)
      Real dqf_minus_i = q_i0_i - qminus_i; // (CS eq 25) = -dQ^- in Mignone's notation
      Real dqf_plus_i  = qplus_i - q_i0_i;

      //--- Step 4b. --------------------------------------------------------------------
      // Nonuniform and/or curvilinear coordinate: apply Mignone limiters to parabolic
      // interpolant. Note, the Mignone limiter does not check for cell-averaged extrema:
      Real qa = dqf_minus_i*dqf_plus_i;
      if (qa <= 0.0) { // Local extrema detected
        qminus_i = q_i0_i;
        qplus_i = q_i0_i;
      } else { // No extrema detected
        // Overshoot i-1/2,R / i,(-) state
        if (std::abs(dqf_minus_i) >= hplus_ratio_i(i)*std::abs(dqf_plus_i)) {
          qminus_i = q_i0_i - hplus_ratio_i(i)*dqf_plus_i;
        }
        // Overshoot i+1/2,L / i,(+) state
        if (std::abs(dqf_plus_i) >= hminus_ratio_i(i)*std::abs(dqf_minus_i)) {
          qplus_i = q_i0_i + hminus_ratio_i(i)*dqf_minus_i;
        }
      }

      // compute ql_(i+1/2) and qr_(i-1/2)
      ql(n,i+1) = qplus_i;
      qr(n,i  ) = qminus_i;
    }
  } // end char PPM loop over =nu

  return;
}

//-------------------------------------------------------------------------------------
//! \fn Reconstruction::PiecewiseParabolicFastX2(const int k, const int j,
//!                              const int il, const int iu,
//!                              const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
//!                              AthenaArray<Real> &wl, AthenaArray<Real> &wr)
//! \brief Returns L/R interface values in X2-dir constructed using fourth-order PPM and
//!         Colella-Sekora or Mignone limiting over [kl,ku][jl,ju][il,iu]

void Reconstruction::PiecewiseParabolicFastX2(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &q,
    AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  const int nu = q.GetDim4() - 1;
  // CS08 constant used in second derivative limiter, >1 , independent of h
  const Real C2 = 1.25;

  // set work arrays used for primitive/characterstic cell-averages to scratch
  AthenaArray<Real> &q_jm2 = scr2_ni_, &q_jm1 = scr3_ni_,
                    &q_j = scr4_ni_, &q_jp1 = scr5_ni_, &q_jp2 = scr6_ni_;

  // cache the x1-sliced primitive states for eigensystem calculation
  if (j==(pmy_block_->js-1)) {
    for (int n=0; n<=nu; ++n) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        q_j  (n,i) = q(n,k,j  ,i);
        q_jm2(n,i) = q(n,k,j-2,i);
        q_jm1(n,i) = q(n,k,j-1,i);
        q_jp1(n,i) = q(n,k,j+1,i);
        q_jp2(n,i) = q(n,k,j+2,i);
      }
    }
  } else { // reuse data
    q_jm2.SwapAthenaArray(q_jm1);
    q_jm1.SwapAthenaArray(q_j);
    q_j.SwapAthenaArray(q_jp1);
    q_jp1.SwapAthenaArray(q_jp2);
    for (int n=0; n<=nu; ++n) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        q_jp2(n,i) = q(n,k,j+2,i);
      }
    }
  }

  //--- Step 1. ------------------------------------------------------------------------
  // Reconstruct interface averages <a>_{j-1/2} and <a>_{j+1/2}
  for (int n=0; n<=nu; ++n) {
    // Compute average slope in j-1, j, j+1 zones
#pragma omp simd simdlen(SIMD_WIDTH)
    for (int i=il; i<=iu; ++i) {
      // nonuniform or uniform Cartesian-like coord reconstruction from volume averages:
      Real dph_i, dph_jp1_i;
      Real q_j0_i  = q_j(n,i);
      Real q_jm1_i = q_jm1(n,i);
      Real q_jm2_i = q_jm2(n,i);
      Real q_jp1_i = q_jp1(n,i);
      Real q_jp2_i = q_jp2(n,i);

      if (!curvilinear_[X2DIR]) {
        Real qa = (q_j0_i - q_jm1_i);
        Real qb = (q_jp1_i - q_j0_i);
        Real dd_jm1_i = c1j(j-1)*qa + c2j(j-1)*(q_jm1_i - q_jm2_i);
        Real dd_i     = c1j(j  )*qb + c2j(j  )*qa;
        Real dd_jp1_i = c1j(j+1)*(q_jp2_i - q_jp1_i) + c2j(j+1)*qb;

        // Approximate interface average at j-1/2 and j+1/2 using PPM (CW eq 1.6)
        // KGF: group the biased stencil quantities to preserve FP symmetry
        dph_i = (c3j(j)*q_jm1_i + c4j(j)*q_j0_i) +
                 (c5j(j)*dd_jm1_i + c6j(j)*dd_i);
        dph_jp1_i = (c3j(j+1)*q_j0_i + c4j(j+1)*q_jp1_i) +
                     (c5j(j+1)*dd_i + c6j(j+1)*dd_jp1_i);
      } else { // spherical-polar coordinates meridional x2 direction:
        // TODO(felker): implement weight calculation in reconstruction.cpp
        dph_i = c1j(j)*q_jm2_i + c2j(j)*q_jm1_i + c3j(j)*q_j0_i
                 + c4j(j)*q_jp1_i;
        dph_jp1_i = c1j(j+1)*q_jm1_i + c2j(j+1)*q_j0_i + c3j(j+1)*q_jp1_i
                     + c4j(j+1)*q_jp2_i;
      }


      //--- Step 2a. ---------------------------------------------------------------------
      // Uniform Cartesian-like coordinate:
      // limit interpolated interface states (CD 4.3.1)
      if (uniform_[X2DIR] && !curvilinear_[X2DIR]) {
        // approximate second derivative at interfaces for smooth extrema preservation

        // KGF: add the off-centered quantities first to preserve FP symmetry
        Real d2qc_jm1_i = q_jm2_i + q_j0_i  - 2.0*q_jm1_i;
        //(CD eq 85a) (no 1/2)
        Real d2qc_i     = q_jm1_i + q_jp1_i - 2.0*q_j0_i;
        Real d2qc_jp1_i = q_j0_i  + q_jp2_i - 2.0*q_jp1_i;

        // j-1/2
        Real qa_tmp = dph_i - q_jm1_i;    // (CD eq 84a)
        Real qb_tmp = q_j0_i - dph_i;     // (CD eq 84b)
        // KGF: add the off-centered quantities first to preserve FP symmetry
        Real qa = 3.0*(q_jm1_i + q_j0_i - 2.0*dph_i);  // (CD eq 85b)
        Real qb = d2qc_jm1_i;    // (CD eq 85a) (no 1/2)
        Real qc = d2qc_i;   // (CD eq 85c) (no 1/2)
        Real qd = 0.0;
        if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)) {
          qd = SIGN(qa)* std::min(C2*std::abs(qb),
                                  std::min(C2*std::abs(qc), std::abs(qa)));
        }
        Real dph_tmp = 0.5*(q_jm1_i + q_j0_i) - qd/6.0;
        if (qa_tmp*qb_tmp < 0.0) { // Local extrema detected at j-1/2 face
          dph_i = dph_tmp;
        }

        // j+1/2
        qa_tmp = dph_jp1_i - q_j0_i;    // (CD eq 84a)
        qb_tmp = q_jp1_i - dph_jp1_i;   // (CD eq 84b)
        // KGF: add the off-centered quantities first to preserve FP symmetry
        qa = 3.0*(q_j0_i + q_jp1_i  - 2.0*dph_jp1_i);  // (CD eq 85b)
        qb = d2qc_i;            // (CD eq 85a) (no 1/2)
        qc = d2qc_jp1_i;   // (CD eq 85c) (no 1/2)
        qd = 0.0;
        if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)) {
          qd = SIGN(qa)* std::min(C2*std::abs(qb),
                                  std::min(C2*std::abs(qc), std::abs(qa)));
        }
        Real dphjp1_tmp = 0.5*(q_j0_i + q_jp1_i) - qd/6.0;
        if (qa_tmp*qb_tmp < 0.0) { // Local extrema detected at j+1/2 face
          dph_jp1_i = dphjp1_tmp;
        }

      //--- Step 2b. -------------------------------------------------------------------
      // Nonuniform and/or curvilinear coordinate: apply strict monotonicity constraints
      // (Mignone eq 45)
      } else {
        dph_i     = std::min(dph_i    , std::max(q_j0_i,q_jm1_i));
        dph_jp1_i = std::min(dph_jp1_i, std::max(q_j0_i,q_jp1_i));

        dph_i     = std::max(dph_i    , std::min(q_j0_i,q_jm1_i));
        dph_jp1_i = std::max(dph_jp1_i, std::min(q_j0_i,q_jp1_i));
      }

      // Cache Riemann states for both non-/uniform limiters
      Real qminus_i = dph_i;
      Real qplus_i =  dph_jp1_i;

      //--- Step 3. ----------------------------------------------------------------------
      // Compute cell-centered difference stencils (MC section 2.4.1)
      Real dqf_minus_i = q_j0_i - qminus_i; // (CS eq 25)
      Real dqf_plus_i  = qplus_i - q_j0_i;

      //--- Step 4b. -------------------------------------------------------------------
      // Nonuniform and/or curvilinear coordinate: apply Mignone limiters to parabolic
      // interpolant. Note, the Mignone limiter does not check for cell-averaged extrema:

      Real qa = dqf_minus_i*dqf_plus_i;
      if (qa <= 0.0) { // Local extrema detected
        qminus_i = q_j0_i;
        qplus_i = q_j0_i;
      } else { // No extrema detected
        // Overshoot j-1/2,R / j,(-) state
        if (std::abs(dqf_minus_i) >= hplus_ratio_j(j)*std::abs(dqf_plus_i)) {
          qminus_i = q_j0_i - hplus_ratio_j(j)*dqf_plus_i;
        }
        // Overshoot j+1/2,L / j,(+) state
        if (std::abs(dqf_plus_i) >= hminus_ratio_j(j)*std::abs(dqf_minus_i)) {
          qplus_i = q_j0_i + hminus_ratio_j(j)*dqf_minus_i;
        }
      }

      // compute ql_(j+1/2) and qr_(j-1/2)
      ql(n,i) = qplus_i;
      qr(n,i) = qminus_i;
    }
  } // end char PPM loop over =nu
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::PiecewiseParabolicFastX3(const int k, const int j,
//!                              const int il, const int iu,
//!                              const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
//!                              AthenaArray<Real> &wl, AthenaArray<Real> &wr)
//! \brief Returns L/R interface values in X3-dir constructed using fourth-order PPM and
//!         Colella-Sekora or Mignone limiting over [kl,ku][jl,ju][il,iu]

void Reconstruction::PiecewiseParabolicFastX3(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &q,
    AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  const int nu = q.GetDim4() - 1;
  // CS08 constant used in second derivative limiter, >1 , independent of h
  const Real C2 = 1.25;

  // set work arrays used for primitive/characterstic cell-averages to scratch
  AthenaArray<Real> &q_km2 = scr2_ni_, &q_km1 = scr3_ni_,
                     &q_k = scr4_ni_, &q_kp1 = scr5_ni_, &q_kp2 = scr6_ni_;

  // cache the x1-sliced primitive states for eigensystem calculation
  if (k==(pmy_block_->ks-1)) {
    for (int n=0; n<=nu; ++n) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        q_k  (n,i) = q(n,k  ,j,i);
        q_km2(n,i) = q(n,k-2,j,i);
        q_km1(n,i) = q(n,k-1,j,i);
        q_kp1(n,i) = q(n,k+1,j,i);
        q_kp2(n,i) = q(n,k+2,j,i);
      }
    }
  } else { // reuse data
    q_km2.SwapAthenaArray(q_km1);
    q_km1.SwapAthenaArray(q_k);
    q_k.SwapAthenaArray(q_kp1);
    q_kp1.SwapAthenaArray(q_kp2);
    for (int n=0; n<=nu; ++n) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        q_kp2(n,i) = q(n,k+2,j,i);
      }
    }
  }

  //--- Step 1. -------------------------------------------------------------------------
  // Reconstruct interface averages <a>_{k-1/2} and <a>_{k+1/2}
  for (int n=0; n<=nu; ++n) {
    // Compute average slope in k-1, k, k+1 zones
#pragma omp simd simdlen(SIMD_WIDTH)
    for (int i=il; i<=iu; ++i) {
      // nonuniform or uniform Cartesian-like coord reconstruction from volume averages:
      Real dph_i, dph_kp1_i;
      Real q_k0_i  = q_k(n,i);
      Real q_km1_i = q_km1(n,i);
      Real q_km2_i = q_km2(n,i);
      Real q_kp1_i = q_kp1(n,i);
      Real q_kp2_i = q_kp2(n,i);

      Real qa = (q_k0_i - q_km1_i);
      Real qb = (q_kp1_i - q_k0_i);
      Real dd_km1_i = c1k(k-1)*qa + c2k(k-1)*(q_km1_i - q_km2_i);
      Real dd_i     = c1k(k  )*qb + c2k(k  )*qa;
      Real dd_kp1_i = c1k(k+1)*(q_kp2_i - q_kp1_i) + c2k(k+1)*qb;

      // Approximate interface average at k-1/2 and k+1/2 using PPM (CW eq 1.6)
      // KGF: group the biased stencil quantities to preserve FP symmetry
      dph_i = (c3k(k)*q_km1_i + c4k(k)*q_k0_i) +
               (c5k(k)*dd_km1_i + c6k(k)*dd_i);
      dph_kp1_i = (c3k(k+1)*q_k0_i + c4k(k+1)*q_kp1_i) +
                   (c5k(k+1)*dd_i + c6k(k+1)*dd_kp1_i);

      //--- Step 2a. --------------------------------------------------------------------
      // Uniform Cartesian-like coordinate:
      // limit interpolated interface states (CD 4.3.1)
      if (uniform_[X3DIR]) {
        // approximate second derivative at interfaces for smooth extrema preservation

        // KGF: add the off-centered quantities first to preserve FP symmetry
        Real d2qc_km1_i = q_km2_i + q_k0_i  - 2.0*q_km1_i;
        // (CD eq 85a) (no 1/2)
        Real d2qc_i     = q_km1_i + q_kp1_i - 2.0*q_k0_i;
        Real d2qc_kp1_i = q_k0_i  + q_kp2_i - 2.0*q_kp1_i;

        // k-1/2
        Real qa_tmp = dph_i - q_km1_i; // (CD eq 84a)
        Real qb_tmp = q_k0_i - dph_i;     // (CD eq 84b)
        // KGF: add the off-centered quantities first to preserve FP symmetry
        Real qa = 3.0*(q_km1_i + q_k0_i - 2.0*dph_i);  // (CD eq 85b)
        Real qb = d2qc_km1_i;    // (CD eq 85a) (no 1/2)
        Real qc = d2qc_i;   // (CD eq 85c) (no 1/2)
        Real qd = 0.0;
        if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)) {
          qd = SIGN(qa)* std::min(C2*std::abs(qb),
                                  std::min(C2*std::abs(qc), std::abs(qa)));
        }
        Real dph_tmp = 0.5*(q_km1_i + q_k0_i) - qd/6.0;
        if (qa_tmp*qb_tmp < 0.0) {  // Local extrema detected at k-1/2 face
          dph_i = dph_tmp;
        }

        // k+1/2
        qa_tmp = dph_kp1_i - q_k0_i;       // (CD eq 84a)
        qb_tmp = q_kp1_i - dph_kp1_i;   // (CD eq 84b)
        // KGF: add the off-centered quantities first to preserve FP symmetry
        qa = 3.0*(q_k0_i + q_kp1_i - 2.0*dph_kp1_i);  // (CD eq 85b)
        qb = d2qc_i;            // (CD eq 85a) (no 1/2)
        qc = d2qc_kp1_i;   // (CD eq 85c) (no 1/2)
        qd = 0.0;
        if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)) {
          qd = SIGN(qa)* std::min(C2*std::abs(qb),
                                  std::min(C2*std::abs(qc), std::abs(qa)));
        }
        Real dphkp1_tmp = 0.5*(q_k0_i + q_kp1_i) - qd/6.0;
        if (qa_tmp*qb_tmp < 0.0) { // Local extrema detected at k+1/2 face
          dph_kp1_i = dphkp1_tmp;
        }

      //--- Step 2b. ---------------------------------------------------------------------
      // Nonuniform coordinate spacing: apply strict monotonicity constraints
      // (Mignone eq 45)
      } else {
        dph_i     = std::min(dph_i    , std::max(q_k0_i,q_km1_i));
        dph_kp1_i = std::min(dph_kp1_i, std::max(q_k0_i,q_kp1_i));

        dph_i     = std::max(dph_i    , std::min(q_k0_i,q_km1_i));
        dph_kp1_i = std::max(dph_kp1_i, std::min(q_k0_i,q_kp1_i));
      }

      // Cache Riemann states for both non-/uniform limiters
      Real qminus_i = dph_i;
      Real qplus_i =  dph_kp1_i;

      //--- Step 3. --------------------------------------------------------------------
      // Compute cell-centered difference stencils (MC section 2.4.1)
      Real dqf_minus_i = q_k0_i - qminus_i; // (CS eq 25)
      Real dqf_plus_i  = qplus_i - q_k0_i;

      //--- Step 4b. ---------------------------------------------------------------------
      // Nonuniform coordinate spacing: apply Mignone limiters to parabolic interpolant
      // Note, the Mignone limiter does not check for cell-averaged extrema:
      qa = dqf_minus_i*dqf_plus_i;
      if (qa <= 0.0) { // Local extrema detected
        qminus_i = q_k0_i;
        qplus_i = q_k0_i;
      } else { // No extrema detected
        // TODO(felker): could delete hplus_ratio_k() arrays for curvilinear PPMx3
        // Overshoot k-1/2,R / k,(-) state
        if (std::abs(dqf_minus_i) >= hplus_ratio_k(k)*std::abs(dqf_plus_i)) {
          qminus_i = q_k0_i - hplus_ratio_k(k)*dqf_plus_i;
        }
        // Overshoot k+1/2,L / k,(+) state
        if (std::abs(dqf_plus_i) >= hminus_ratio_k(k)*std::abs(dqf_minus_i)) {
          qplus_i = q_k0_i + hminus_ratio_k(k)*dqf_minus_i;
        }
      }

      // compute ql_(k+1/2) and qr_(k-1/2)
      ql(n,i) = qplus_i;
      qr(n,i) = qminus_i;
    }
  }
  return;
}
