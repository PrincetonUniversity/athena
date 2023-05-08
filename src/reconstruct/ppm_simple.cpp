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
//! \fn Reconstruction::PiecewiseParabolicX1(const int k, const int j,
//!                              const int il, const int iu,
//!                              const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
//!                              AthenaArray<Real> &wl, AthenaArray<Real> &wr)
//! \brief Returns L/R interface values in X1-dir constructed using fourth-order PPM and
//!        Colella-Sekora or Mignone limiting over [kl,ku][jl,ju][il,iu]

void Reconstruction::PiecewiseParabolicX1(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &q,
    AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  const int nu = q.GetDim4() - 1;

  // CS08 constant used in second derivative limiter, >1 , independent of h
  const Real C2 = 1.25;

  // TODO(felker): renumber scratch array references; not using 2x from ppm.cpp
  // bx (MHD) and wc (characteristic projection)

  // set work arrays used for primitive/characterstic cell-averages to scratch
  AthenaArray<Real> &q_im2 = scr2_ni_, &q_im1 = scr3_ni_,
                     &q_i = scr4_ni_, &q_ip1 = scr5_ni_, &q_ip2 = scr6_ni_,
                &qr_imh = scr7_ni_, &ql_iph = scr8_ni_;

  // set work PPM work arrays to shallow copies of scratch arrays:
  AthenaArray<Real> &dd = scr02_i_, &dd_im1 = scr03_i_, &dd_ip1 = scr04_i_,
                   &dph = scr05_i_, &dph_ip1 = scr06_i_;

  AthenaArray<Real> &d2qc_im1 = scr07_i_, &d2qc = scr08_i_, &d2qc_ip1 = scr09_i_,
                        &d2qf = scr10_i_;

  AthenaArray<Real> &qplus = scr11_i_, &qminus = scr12_i_, &dqf_plus = scr13_i_,
                &dqf_minus = scr14_i_;

  // cache the x1-sliced primitive states for eigensystem calculation
  for (int n=0; n<=nu; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      q_i  (n,i) = q(n,k,j,i  );
      q_im2(n,i) = q(n,k,j,i-2);
      q_im1(n,i) = q(n,k,j,i-1);
      q_ip1(n,i) = q(n,k,j,i+1);
      q_ip2(n,i) = q(n,k,j,i+2);
    }
  }

  //--- Step 1. --------------------------------------------------------------------------
  // Reconstruct interface averages <a>_{i-1/2} and <a>_{i+1/2}
  for (int n=0; n<=nu; ++n) {
    // Compute average slope in i-1, i, i+1 zones
#pragma omp simd simdlen(SIMD_WIDTH)
    for (int i=il; i<=iu; ++i) {
      // nonuniform or uniform Cartesian-like coord reconstruction from volume averages:
      if (!curvilinear[X1DIR]) {
        Real qa = (q_i(n,i) - q_im1(n,i));
        Real qb = (q_ip1(n,i) - q_i(n,i));
        dd_im1(i) = c1i(i-1)*qa + c2i(i-1)*(q_im1(n,i) - q_im2(n,i));
        dd    (i) = c1i(i  )*qb + c2i(i  )*qa;
        dd_ip1(i) = c1i(i+1)*(q_ip2(n,i) - q_ip1(n,i)) + c2i(i+1)*qb;

        // Approximate interface average at i-1/2 and i+1/2 using PPM (CW eq 1.6)
        // KGF: group the biased stencil quantities to preserve FP symmetry
        dph(i) = (c3i(i)*q_im1(n,i) + c4i(i)*q_i(n,i)) +
                 (c5i(i)*dd_im1(i) + c6i(i)*dd(i));
        dph_ip1(i) = (c3i(i+1)*q_i(n,i) + c4i(i+1)*q_ip1(n,i)) +
                     (c5i(i+1)*dd(i) + c6i(i+1)*dd_ip1(i) );
      } else { // radial coordinate
        dph(i) = c1i(i)*q_im2(n,i) + c2i(i)*q_im1(n,i) + c3i(i)*q_i(n,i)
                 + c4i(i)*q_ip1(n,i);
        dph_ip1(i) = c1i(i+1)*q_im1(n,i) + c2i(i+1)*q_i(n,i) + c3i(i+1)*q_ip1(n,i)
                     + c4i(i+1)*q_ip2(n,i);
      }
    }

    //--- Step 2a. -----------------------------------------------------------------------
    // Uniform Cartesian-like coordinate: limit interpolated interface states (CD 4.3.1)
    if (uniform[X1DIR] && !curvilinear[X1DIR]) {
      // approximate second derivative at interfaces for smooth extrema preservation
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu+1; ++i) {
        // KGF: add the off-centered quantities first to preserve FP symmetry
        d2qc_im1(i) = q_im2(n,i) + q_i  (n,i) - 2.0*q_im1(n,i);
        d2qc    (i) = q_im1(n,i) + q_ip1(n,i) - 2.0*q_i  (n,i); // (CD eq 85a) (no 1/2)
        d2qc_ip1(i) = q_i  (n,i) + q_ip2(n,i) - 2.0*q_ip1(n,i);
      }

      // i-1/2
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=(iu+1); ++i) {
        Real qa_tmp = dph(i) - q_im1(n,i); // (CD eq 84a)
        Real qb_tmp = q_i(n,i) - dph(i);     // (CD eq 84b)
        // KGF: add the off-centered quantities first to preserve FP symmetry
        Real qa = 3.0*(q_im1(n,i) + q_i(n,i)  - 2.0*dph(i));  // (CD eq 85b)
        Real qb = d2qc_im1(i);    // (CD eq 85a) (no 1/2)
        Real qc = d2qc(i);   // (CD eq 85c) (no 1/2)
        Real qd = 0.0;
        if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)) {
          qd = SIGN(qa)* std::min(C2*std::abs(qb),
                                  std::min(C2*std::abs(qc), std::abs(qa)));
        }
        Real dph_tmp = 0.5*(q_im1(n,i) + q_i(n,i)) - qd/6.0;
        if (qa_tmp*qb_tmp < 0.0) { // Local extrema detected at i-1/2 face
          dph(i) = dph_tmp;
        }
      }

      // i+1/2
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=(iu+1); ++i) {
        Real qa_tmp = dph_ip1(i) - q_i(n,i);       // (CD eq 84a)
        Real qb_tmp = q_ip1(n,i) - dph_ip1(i);   // (CD eq 84b)
        // KGF: add the off-centered quantities first to preserve FP symmetry
        Real qa = 3.0*(q_i(n,i) + q_ip1(n,i) - 2.0*dph_ip1(i));  // (CD eq 85b)
        Real qb = d2qc(i);            // (CD eq 85a) (no 1/2)
        Real qc = d2qc_ip1(i);   // (CD eq 85c) (no 1/2)
        Real qd = 0.0;
        if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)) {
          qd = SIGN(qa)* std::min(C2*std::abs(qb),
                                  std::min(C2*std::abs(qc), std::abs(qa)));
        }
        Real dphip1_tmp = 0.5*(q_i(n,i) + q_ip1(n,i)) - qd/6.0;
        if (qa_tmp*qb_tmp < 0.0) { // Local extrema detected at i+1/2 face
          dph_ip1(i) = dphip1_tmp;
        }
      }

#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        // KGF: add the off-centered quantities first to preserve FP symmetry
        d2qf(i) = 6.0*(dph(i) + dph_ip1(i) - 2.0*q_i(n,i)); // a6 coefficient * -2
      }

      //--- Step 2b. ---------------------------------------------------------------------
      // Nonuniform and/or curvilinear coordinate: apply strict monotonicity constraints
      // (Mignone eq 45)
    } else {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        dph    (i) = std::min(dph    (i), std::max(q_i(n,i),q_im1(n,i)));
        dph_ip1(i) = std::min(dph_ip1(i), std::max(q_i(n,i),q_ip1(n,i)));

        dph    (i) = std::max(dph    (i), std::min(q_i(n,i),q_im1(n,i)));
        dph_ip1(i) = std::max(dph_ip1(i), std::min(q_i(n,i),q_ip1(n,i)));
      }
    }

    // Cache Riemann states for both non-/uniform limiters
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      qminus(i) = dph(i  );
      qplus(i) =  dph_ip1(i );
    }

    //--- Step 3. ------------------------------------------------------------------------
    // Compute cell-centered difference stencils (MC section 2.4.1)
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      dqf_minus(i) = q_i(n,i) - qminus(i); // (CS eq 25) = -dQ^- in Mignone's notation
      dqf_plus(i)  = qplus(i) - q_i(n,i);
    }

    //--- Step 4a. -----------------------------------------------------------------------
    // For uniform Cartesian-like coordinate: apply CS limiters to parabolic interpolant
    if (uniform[X1DIR] && !curvilinear[X1DIR]) {
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        Real qa_tmp = dqf_minus(i)*dqf_plus(i);
        Real qb_tmp = (q_ip1(n,i) - q_i(n,i))*(q_i(n,i) - q_im1(n,i));

        Real qa = d2qc_im1(i);
        Real qb = d2qc(i);
        Real qc = d2qc_ip1(i);
        Real qd = d2qf(i);
        Real qe = 0.0;
        if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc) && SIGN(qa) == SIGN(qd)) {
          // Extrema is smooth
          qe = SIGN(qd)* std::min(std::min(C2*std::abs(qa), C2*std::abs(qb)),
                                  std::min(C2*std::abs(qc),
                                           std::abs(qd))); // (CS eq 22)
        }

        // Check if 2nd derivative is close to roundoff error
        qa = std::max(std::abs(q_im1(n,i)), std::abs(q_im2(n,i)));
        qb = std::max(std::max(std::abs(q_i(n,i)), std::abs(q_ip1(n,i))),
                      std::abs(q_ip2(n,i)));

        Real rho = 0.0;
        if (std::abs(qd) > (1.0e-12)*std::max(qa,qb)) {
          // Limiter is not sensitive to roundoff. Use limited ratio (MC eq 27)
          rho = qe/qd;
        }

        Real tmp_m = q_i(n,i) - rho*dqf_minus(i);
        Real tmp_p = q_i(n,i) + rho*dqf_plus(i);
        Real tmp2_m = q_i(n,i) - 2.0*dqf_plus(i);
        Real tmp2_p = q_i(n,i) + 2.0*dqf_minus(i);

        // Check for local extrema
        if ((qa_tmp <= 0.0 || qb_tmp <= 0.0)) {
          // Check if relative change in limited 2nd deriv is > roundoff
          if (rho <= (1.0 - (1.0e-12))) {
            // Limit smooth extrema
            qminus(i) = tmp_m;// (CS eq 23)
            qplus(i) = tmp_p;
          }
          // No extrema detected
        } else {
          // Overshoot i-1/2,R / i,(-) state
          if (std::abs(dqf_minus(i)) >= 2.0*std::abs(dqf_plus(i))) {
            qminus(i) = tmp2_m;
          }
          // Overshoot i+1/2,L / i,(+) state
          if (std::abs(dqf_plus(i)) >= 2.0*std::abs(dqf_minus(i))) {
            qplus(i) = tmp2_p;
          }
        }
      }

      //--- Step 4b. ---------------------------------------------------------------------
      // Nonuniform and/or curvilinear coordinate: apply Mignone limiters to parabolic
      // interpolant. Note, the Mignone limiter does not check for cell-averaged extrema:
    } else {
      for (int i=il; i<=iu; ++i) {
        Real qa = dqf_minus(i)*dqf_plus(i);
        if (qa <= 0.0) { // Local extrema detected
          qminus(i) = q_i(n,i);
          qplus(i) = q_i(n,i);
        } else { // No extrema detected
          // Overshoot i-1/2,R / i,(-) state
          if (std::abs(dqf_minus(i)) >= hplus_ratio_i(i)*std::abs(dqf_plus(i))) {
            qminus(i) = q_i(n,i) - hplus_ratio_i(i)*dqf_plus(i);
          }
          // Overshoot i+1/2,L / i,(+) state
          if (std::abs(dqf_plus(i)) >= hminus_ratio_i(i)*std::abs(dqf_minus(i))) {
            qplus(i) = q_i(n,i) + hminus_ratio_i(i)*dqf_minus(i);
          }
        }
      }
    }
    //--- Step 5. ------------------------------------------------------------------------
    // Convert limited cell-centered values to interface-centered L/R Riemann states
    // both L/R values defined over [il,iu]
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      ql_iph(n,i ) = qplus(i);
      qr_imh(n,i ) = qminus(i);
    }
  } // end char PPM loop over =nu

  // compute ql_(i+1/2) and qr_(i-1/2)
  for (int n=0; n<=nu; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      ql(n,i+1) = ql_iph(n,i);
      qr(n,i  ) = qr_imh(n,i);
    }
  }
  return;
}

//-------------------------------------------------------------------------------------
//! \fn Reconstruction::PiecewiseParabolicX2(const int k, const int j,
//!                              const int il, const int iu,
//!                              const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
//!                              AthenaArray<Real> &wl, AthenaArray<Real> &wr)
//! \brief Returns L/R interface values in X2-dir constructed using fourth-order PPM and
//!         Colella-Sekora or Mignone limiting over [kl,ku][jl,ju][il,iu]

void Reconstruction::PiecewiseParabolicX2(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &q,
    AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  const int nu = q.GetDim4() - 1;
  // CS08 constant used in second derivative limiter, >1 , independent of h
  const Real C2 = 1.25;

  // set work arrays used for primitive/characterstic cell-averages to scratch
  AthenaArray<Real> &q_jm2 = scr2_ni_, &q_jm1 = scr3_ni_,
                     &q_j = scr4_ni_, &q_jp1 = scr5_ni_, &q_jp2 = scr6_ni_,
                &qr_jmh = scr7_ni_, &ql_jph = scr8_ni_;

  // set work PPM work arrays to shallow copies of scratch arrays:
  AthenaArray<Real> &dd = scr02_i_, &dd_jm1 = scr03_i_, &dd_jp1 = scr04_i_,
                   &dph = scr05_i_, &dph_jp1 = scr06_i_;

  AthenaArray<Real> &d2qc_jm1 = scr07_i_, &d2qc = scr08_i_, &d2qc_jp1 = scr09_i_,
                        &d2qf = scr10_i_;

  AthenaArray<Real> &qplus = scr11_i_, &qminus = scr12_i_, &dqf_plus = scr13_i_,
                &dqf_minus = scr14_i_;

  // cache the x1-sliced primitive states for eigensystem calculation
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

  //--- Step 1. ------------------------------------------------------------------------
  // Reconstruct interface averages <a>_{j-1/2} and <a>_{j+1/2}
  for (int n=0; n<=nu; ++n) {
    // Compute average slope in j-1, j, j+1 zones
#pragma omp simd simdlen(SIMD_WIDTH)
    for (int i=il; i<=iu; ++i) {
      // nonuniform or uniform Cartesian-like coord reconstruction from volume averages:
      if (!curvilinear[X2DIR]) {
        Real qa = (q_j(n,i) - q_jm1(n,i));
        Real qb = (q_jp1(n,i) - q_j(n,i));
        dd_jm1(i) = c1j(j-1)*qa + c2j(j-1)*(q_jm1(n,i) - q_jm2(n,i));
        dd    (i) = c1j(j  )*qb + c2j(j  )*qa;
        dd_jp1(i) = c1j(j+1)*(q_jp2(n,i) - q_jp1(n,i)) + c2j(j+1)*qb;

        // Approximate interface average at j-1/2 and j+1/2 using PPM (CW eq 1.6)
        // KGF: group the biased stencil quantities to preserve FP symmetry
        dph(i) = (c3j(j)*q_jm1(n,i) + c4j(j)*q_j(n,i)) +
                 (c5j(j)*dd_jm1(i) + c6j(j)*dd(i));
        dph_jp1(i) = (c3j(j+1)*q_j(n,i) + c4j(j+1)*q_jp1(n,i)) +
                     (c5j(j+1)*dd(i) + c6j(j+1)*dd_jp1(i));
      } else { // spherical-polar coordinates meridional x2 direction:
        // TODO(felker): implement weight calculation in reconstruction.cpp
        dph(i) = c1j(j)*q_jm2(n,i) + c2j(j)*q_jm1(n,i) + c3j(j)*q_j(n,i)
                 + c4j(j)*q_jp1(n,i);
        dph_jp1(i) = c1j(j+1)*q_jm1(n,i) + c2j(j+1)*q_j(n,i) + c3j(j+1)*q_jp1(n,i)
                     + c4j(j+1)*q_jp2(n,i);
      }
    }

    //--- Step 2a. ---------------------------------------------------------------------
    // Uniform Cartesian-like coordinate: limit interpolated interface states (CD 4.3.1)
    if (uniform[X2DIR] && !curvilinear[X2DIR]) {
      // approximate second derivative at interfaces for smooth extrema preservation
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        // KGF: add the off-centered quantities first to preserve FP symmetry
        d2qc_jm1(i) = q_jm2(n,i) + q_j  (n,i) - 2.0*q_jm1(n,i);
        d2qc    (i) = q_jm1(n,i) + q_jp1(n,i) - 2.0*q_j  (n,i); //(CD eq 85a) (no 1/2)
        d2qc_jp1(i) = q_j  (n,i) + q_jp2(n,i) - 2.0*q_jp1(n,i);
      }

      // j-1/2
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        Real qa_tmp = dph(i) - q_jm1(n,i); // (CD eq 84a)
        Real qb_tmp = q_j(n,i) - dph(i);     // (CD eq 84b)
        // KGF: add the off-centered quantities first to preserve FP symmetry
        Real qa = 3.0*(q_jm1(n,i) + q_j(n,i) - 2.0*dph(i));  // (CD eq 85b)
        Real qb = d2qc_jm1(i);    // (CD eq 85a) (no 1/2)
        Real qc = d2qc(i);   // (CD eq 85c) (no 1/2)
        Real qd = 0.0;
        if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)) {
          qd = SIGN(qa)* std::min(C2*std::abs(qb),
                                  std::min(C2*std::abs(qc), std::abs(qa)));
        }
        Real dph_tmp = 0.5*(q_jm1(n,i) + q_j(n,i)) - qd/6.0;
        if (qa_tmp*qb_tmp < 0.0) { // Local extrema detected at j-1/2 face
          dph(i) = dph_tmp;
        }
      }
      // j+1/2
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        Real qa_tmp = dph_jp1(i) - q_j(n,i);       // (CD eq 84a)
        Real qb_tmp = q_jp1(n,i) - dph_jp1(i);   // (CD eq 84b)
        // KGF: add the off-centered quantities first to preserve FP symmetry
        Real qa = 3.0*(q_j(n,i) + q_jp1(n,i)  - 2.0*dph_jp1(i));  // (CD eq 85b)
        Real qb = d2qc(i);            // (CD eq 85a) (no 1/2)
        Real qc = d2qc_jp1(i);   // (CD eq 85c) (no 1/2)
        Real qd = 0.0;
        if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)) {
          qd = SIGN(qa)* std::min(C2*std::abs(qb),
                                  std::min(C2*std::abs(qc), std::abs(qa)));
        }
        Real dphjp1_tmp = 0.5*(q_j(n,i) + q_jp1(n,i)) - qd/6.0;
        if (qa_tmp*qb_tmp < 0.0) { // Local extrema detected at j+1/2 face
          dph_jp1(i) = dphjp1_tmp;
        }
      }

#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        // KGF: add the off-centered quantities first to preserve FP symmetry
        d2qf(i) = 6.0*(dph(i) + dph_jp1(i) - 2.0*q_j(n,i)); // a6 coefficient * -2
      }

      //--- Step 2b. -------------------------------------------------------------------
      // Nonuniform and/or curvilinear coordinate: apply strict monotonicity constraints
      // (Mignone eq 45)
    } else {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        dph    (i) = std::min(dph    (i), std::max(q_j(n,i),q_jm1(n,i)));
        dph_jp1(i) = std::min(dph_jp1(i), std::max(q_j(n,i),q_jp1(n,i)));

        dph    (i) = std::max(dph    (i), std::min(q_j(n,i),q_jm1(n,i)));
        dph_jp1(i) = std::max(dph_jp1(i), std::min(q_j(n,i),q_jp1(n,i)));
      }
    }

    // Cache Riemann states for both non-/uniform limiters
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      qminus(i) = dph(i  );
      qplus(i) =  dph_jp1(i );
    }

    //--- Step 3. ----------------------------------------------------------------------
    // Compute cell-centered difference stencils (MC section 2.4.1)
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      dqf_minus(i) = q_j(n,i) - qminus(i); // (CS eq 25)
      dqf_plus(i)  = qplus(i) - q_j(n,i);
    }

    //--- Step 4a. ---------------------------------------------------------------------
    // For uniform Cartesian-like coordinate: apply CS limiters to parabolic interpolant
    if (uniform[X2DIR] && !curvilinear[X2DIR]) {
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        Real qa_tmp = dqf_minus(i)*dqf_plus(i);
        Real qb_tmp = (q_jp1(n,i) - q_j(n,i))*(q_j(n,i) - q_jm1(n,i));

        Real qa = d2qc_jm1(i);
        Real qb = d2qc(i);
        Real qc = d2qc_jp1(i);
        Real qd = d2qf(i);
        Real qe = 0.0;
        if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc) && SIGN(qa) == SIGN(qd)) {
          // Extrema is smooth
          qe = SIGN(qd)* std::min(std::min(C2*std::abs(qa), C2*std::abs(qb)),
                                  std::min(C2*std::abs(qc), std::abs(qd)));// (CS eq 22)
        }

        // Check if 2nd derivative is close to roundoff error
        qa = std::max(std::abs(q_jm1(n,i)), std::abs(q_jm2(n,i)));
        qb = std::max(std::max(std::abs(q_j(n,i)),
                               std::abs(q_jp1(n,i))), std::abs(q_jp2(n,i)));

        Real rho = 0.0;
        if (std::abs(qd) > (1.0e-12)*std::max(qa, qb)) {
          // Limiter is not sensitive to roundoff. Use limited ratio (MC eq 27)
          rho = qe/qd;
        }

        Real tmp_m = q_j(n,i) - rho*dqf_minus(i);
        Real tmp_p = q_j(n,i) + rho*dqf_plus(i);
        Real tmp2_m = q_j(n,i) - 2.0*dqf_plus(i);
        Real tmp2_p = q_j(n,i) + 2.0*dqf_minus(i);

        // Check if relative change in limited 2nd deriv is > roundoff
        // Check for local extrema
        if ((qa_tmp <= 0.0 || qb_tmp <= 0.0)) {
          if (rho <= (1.0 - (1.0e-12))) {
            // Limit smooth extrema
            qminus(i) = tmp_m; // (CS eq 23)
            qplus(i) = tmp_p;
          }
          // No extrema detected
        } else {
          // Overshoot j-1/2,R / j,(-) state
          if (std::abs(dqf_minus(i)) >= 2.0*std::abs(dqf_plus(i))) {
            qminus(i) = tmp2_m;
          }
          // Overshoot j+1/2,L / j,(+) state
          if (std::abs(dqf_plus(i)) >= 2.0*std::abs(dqf_minus(i))) {
            qplus(i) = tmp2_p;
          }
        }
      }

      //--- Step 4b. -------------------------------------------------------------------
      // Nonuniform and/or curvilinear coordinate: apply Mignone limiters to parabolic
      // interpolant. Note, the Mignone limiter does not check for cell-averaged extrema:
    } else {
      for (int i=il; i<=iu; ++i) {
        Real qa = dqf_minus(i)*dqf_plus(i);
        if (qa <= 0.0) { // Local extrema detected
          qminus(i) = q_j(n,i);
          qplus(i) = q_j(n,i);
        } else { // No extrema detected
          // Overshoot j-1/2,R / j,(-) state
          if (std::abs(dqf_minus(i)) >= hplus_ratio_j(j)*std::abs(dqf_plus(i))) {
            qminus(i) = q_j(n,i) - hplus_ratio_j(j)*dqf_plus(i);
          }
          // Overshoot j+1/2,L / j,(+) state
          if (std::abs(dqf_plus(i)) >= hminus_ratio_j(j)*std::abs(dqf_minus(i))) {
            qplus(i) = q_j(n,i) + hminus_ratio_j(j)*dqf_minus(i);
          }
        }
      }
    }
    //--- Step 5. ----------------------------------------------------------------------
    // Convert limited cell-centered values to interface-centered L/R Riemann states
    // both L/R values defined over [il,iu]
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      ql_jph(n,i ) = qplus(i);
      qr_jmh(n,i ) = qminus(i);
    }
  } // end char PPM loop over =nu


  // compute ql_(j+1/2) and qr_(j-1/2)
  for (int n=0; n<=nu; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      ql(n,i) = ql_jph(n,i);
      qr(n,i) = qr_jmh(n,i);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::PiecewiseParabolicX3(const int k, const int j,
//!                              const int il, const int iu,
//!                              const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
//!                              AthenaArray<Real> &wl, AthenaArray<Real> &wr)
//! \brief Returns L/R interface values in X3-dir constructed using fourth-order PPM and
//!         Colella-Sekora or Mignone limiting over [kl,ku][jl,ju][il,iu]

void Reconstruction::PiecewiseParabolicX3(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &q,
    AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  const int nu = q.GetDim4() - 1;
  // CS08 constant used in second derivative limiter, >1 , independent of h
  const Real C2 = 1.25;

  // set work arrays used for primitive/characterstic cell-averages to scratch
  AthenaArray<Real> &q_km2 = scr2_ni_, &q_km1 = scr3_ni_,
                     &q_k = scr4_ni_, &q_kp1 = scr5_ni_, &q_kp2 = scr6_ni_,
                &qr_kmh = scr7_ni_, &ql_kph = scr8_ni_;

  // set work PPM work arrays to shallow copies of scratch arrays:
  AthenaArray<Real> &dd = scr02_i_, &dd_km1 = scr03_i_, &dd_kp1 = scr04_i_,
                   &dph = scr05_i_, &dph_kp1 = scr06_i_;

  AthenaArray<Real> &d2qc_km1 = scr07_i_, &d2qc = scr08_i_,
                    &d2qc_kp1 = scr09_i_, &d2qf = scr10_i_;

  AthenaArray<Real> &qplus = scr11_i_, &qminus = scr12_i_, &dqf_plus = scr13_i_,
                &dqf_minus = scr14_i_;

  // cache the x1-sliced primitive states for eigensystem calculation
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

  //--- Step 1. -------------------------------------------------------------------------
  // Reconstruct interface averages <a>_{k-1/2} and <a>_{k+1/2}
  for (int n=0; n<=nu; ++n) {
    // Compute average slope in k-1, k, k+1 zones
#pragma omp simd simdlen(SIMD_WIDTH)
    for (int i=il; i<=iu; ++i) {
      // nonuniform or uniform Cartesian-like coord reconstruction from volume averages:
      Real qa = (q_k(n,i) - q_km1(n,i));
      Real qb = (q_kp1(n,i) - q_k(n,i));
      dd_km1(i) = c1k(k-1)*qa + c2k(k-1)*(q_km1(n,i) - q_km2(n,i));
      dd    (i) = c1k(k  )*qb + c2k(k  )*qa;
      dd_kp1(i) = c1k(k+1)*(q_kp2(n,i) - q_kp1(n,i)) + c2k(k+1)*qb;

      // Approximate interface average at k-1/2 and k+1/2 using PPM (CW eq 1.6)
      // KGF: group the biased stencil quantities to preserve FP symmetry
      dph(i) = (c3k(k)*q_km1(n,i) + c4k(k)*q_k(n,i)) +
               (c5k(k)*dd_km1(i) + c6k(k)*dd(i));
      dph_kp1(i) = (c3k(k+1)*q_k(n,i) + c4k(k+1)*q_kp1(n,i)) +
                   (c5k(k+1)*dd(i) + c6k(k+1)*dd_kp1(i));
    }

    //--- Step 2a. -----------------------------------------------------------------------
    // Uniform Cartesian-like coordinate: limit interpolated interface states (CD 4.3.1)
    if (uniform[X3DIR]) {
      // approximate second derivative at interfaces for smooth extrema preservation
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        // KGF: add the off-centered quantities first to preserve FP symmetry
        d2qc_km1(i) = q_km2(n,i) + q_k  (n,i) - 2.0*q_km1(n,i);
        d2qc    (i) = q_km1(n,i) + q_kp1(n,i) - 2.0*q_k  (n,i); // (CD eq 85a) (no 1/2)
        d2qc_kp1(i) = q_k  (n,i) + q_kp2(n,i) - 2.0*q_kp1(n,i);
      }

      // k-1/2
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        Real qa_tmp = dph(i) - q_km1(n,i); // (CD eq 84a)
        Real qb_tmp = q_k(n,i) - dph(i);     // (CD eq 84b)
        // KGF: add the off-centered quantities first to preserve FP symmetry
        Real qa = 3.0*(q_km1(n,i) + q_k(n,i) - 2.0*dph(i));  // (CD eq 85b)
        Real qb = d2qc_km1(i);    // (CD eq 85a) (no 1/2)
        Real qc = d2qc(i);   // (CD eq 85c) (no 1/2)
        Real qd = 0.0;
        if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)) {
          qd = SIGN(qa)* std::min(C2*std::abs(qb),
                                  std::min(C2*std::abs(qc), std::abs(qa)));
        }
        Real dph_tmp = 0.5*(q_km1(n,i) + q_k(n,i)) - qd/6.0;
        if (qa_tmp*qb_tmp < 0.0) {  // Local extrema detected at k-1/2 face
          dph(i) = dph_tmp;
        }
      }
      // k+1/2
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        Real qa_tmp = dph_kp1(i) - q_k(n,i);       // (CD eq 84a)
        Real qb_tmp = q_kp1(n,i) - dph_kp1(i);   // (CD eq 84b)
        // KGF: add the off-centered quantities first to preserve FP symmetry
        Real qa = 3.0*(q_k(n,i) + q_kp1(n,i) - 2.0*dph_kp1(i));  // (CD eq 85b)
        Real qb = d2qc(i);            // (CD eq 85a) (no 1/2)
        Real qc = d2qc_kp1(i);   // (CD eq 85c) (no 1/2)
        Real qd = 0.0;
        if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)) {
          qd = SIGN(qa)* std::min(C2*std::abs(qb),
                                  std::min(C2*std::abs(qc), std::abs(qa)));
        }
        Real dphkp1_tmp = 0.5*(q_k(n,i) + q_kp1(n,i)) - qd/6.0;
        if (qa_tmp*qb_tmp < 0.0) { // Local extrema detected at k+1/2 face
          dph_kp1(i) = dphkp1_tmp;
        }
      }

#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        // KGF: add the off-centered quantities first to preserve FP symmetry
        d2qf(i) = 6.0*(dph(i) + dph_kp1(i) - 2.0*q_k(n,i)); // a6 coefficient * -2
      }
      //--- Step 2b. ---------------------------------------------------------------------
      // Nonuniform coordinate spacing: apply strict monotonicity constraints
      // (Mignone eq 45)
    } else {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        dph    (i) = std::min(dph    (i), std::max(q_k(n,i),q_km1(n,i)));
        dph_kp1(i) = std::min(dph_kp1(i), std::max(q_k(n,i),q_kp1(n,i)));

        dph    (i) = std::max(dph    (i), std::min(q_k(n,i),q_km1(n,i)));
        dph_kp1(i) = std::max(dph_kp1(i), std::min(q_k(n,i),q_kp1(n,i)));
      }
    }

    // Cache Riemann states for both non-/uniform limiters
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      qminus(i) = dph(i  );
      qplus(i) =  dph_kp1(i );
    }

    //--- Step 3. --------------------------------------------------------------------
    // Compute cell-centered difference stencils (MC section 2.4.1)
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      dqf_minus(i) = q_k(n,i) - qminus(i); // (CS eq 25)
      dqf_plus(i)  = qplus(i) - q_k(n,i);
    }

    //--- Step 4a. -----------------------------------------------------------------------
    // For uniform Cartesian-like coordinate: apply CS limiters to parabolic interpolant
    if (uniform[X3DIR]) {
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        Real qa_tmp = dqf_minus(i)*dqf_plus(i);
        Real qb_tmp = (q_kp1(n,i) - q_k(n,i))*(q_k(n,i) - q_km1(n,i));

        // Check if extrema is smooth
        Real qa = d2qc_km1(i);
        Real qb = d2qc(i);
        Real qc = d2qc_kp1(i);
        Real qd = d2qf(i);
        Real qe = 0.0;
        if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc) && SIGN(qa) == SIGN(qd)) {
          // Extrema is smooth
          qe = SIGN(qd)* std::min(std::min(C2*std::abs(qa),C2*std::abs(qb)),
                                  std::min(C2*std::abs(qc),std::abs(qd))); // (CS eq 22)
        }

        // Check if 2nd derivative is close to roundoff error
        qa = std::max(std::abs(q_km1(n,i)),std::abs(q_km2(n,i)));
        qb = std::max(std::max(std::abs(q_k(n,i)),
                               std::abs(q_kp1(n,i))), std::abs(q_kp2(n,i)));

        Real rho = 0.0;
        if (std::abs(qd) > (1.0e-12)*std::max(qa,qb)) {
          // Limiter is not sensitive to roundoff. Use limited ratio (MC eq 27)
          rho = qe/qd;
        }

        Real tmp_m = q_k(n,i) - rho*dqf_minus(i);
        Real tmp_p = q_k(n,i) + rho*dqf_plus(i);
        Real tmp2_m = q_k(n,i) - 2.0*dqf_plus(i);
        Real tmp2_p = q_k(n,i) + 2.0*dqf_minus(i);

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
      //--- Step 4b. ---------------------------------------------------------------------
      // Nonuniform coordinate spacing: apply Mignone limiters to parabolic interpolant
      // Note, the Mignone limiter does not check for cell-averaged extrema:
    } else {
      for (int i=il; i<=iu; ++i) {
        Real qa = dqf_minus(i)*dqf_plus(i);
        if (qa <= 0.0) { // Local extrema detected
          qminus(i) = q_k(n,i);
          qplus(i) = q_k(n,i);
        } else { // No extrema detected
          // TODO(felker): could delete hplus_ratio_k() arrays for curvilinear PPMx3
          // Overshoot k-1/2,R / k,(-) state
          if (std::abs(dqf_minus(i)) >= hplus_ratio_k(k)*std::abs(dqf_plus(i))) {
            qminus(i) = q_k(n,i) - hplus_ratio_k(k)*dqf_plus(i);
          }
          // Overshoot k+1/2,L / k,(+) state
          if (std::abs(dqf_plus(i)) >= hminus_ratio_k(k)*std::abs(dqf_minus(i))) {
            qplus(i) = q_k(n,i) + hminus_ratio_k(k)*dqf_minus(i);
          }
        }
      }
    }
    //--- Step 5. --------------------------------------------------------------------
    // Convert limited cell-centered values to interface-centered L/R Riemann states
    // both L/R values defined over [il,iu]
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      ql_kph(n,i ) = qplus(i);
      qr_kmh(n,i ) = qminus(i);
    }
  } // end char PPM loop over =nu

  // compute ql_(k+1/2) and qr_(k-1/2)
  for (int n=0; n<=nu; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      ql(n,i) = ql_kph(n,i);
      qr(n,i) = qr_kmh(n,i);
    }
  }
  return;
}
