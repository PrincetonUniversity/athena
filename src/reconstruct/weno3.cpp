//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file weno3.cpp
//  \brief third-order accurate weighted essentially non-oscillatory (WENO_3) scheme
//
// REFERENCES:
// (Mignone) A. Mignone, "High-order conservative reconstruction schemes for finite volume
// methods in cylindrical and spherical coordinates", JCP, 270, 784 (2014)
//
// (Classical WENO3) G.C. Jiang, C.W. Shu, Efficient implementation of weighted ENO
// schemes, J. Comput. Phys. 126 (1996) 202–228
//
// (ESWENO) N.K. Yamaleev, M.H. Carpenter, Third-order energy stable WENO scheme,
// J. Comput. Phys. 228 (2009) 3025–3047
//
// A. Mignone, P. Tzeferacos, G. Bodo, High-order conservative finite difference GLM-MHD
// schemes for cell-centered MHD, J. Comput. Phys. 229 (2010) 5896–5920
//======================================================================================

// C headers

// C++ headers
#include <algorithm>    // max()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "reconstruction.hpp"

// number of stencils per cell; WENO3 consists of 2x two-point substencils per cell
#define NSUBSTENCIL 2

namespace {
// Optimal linear weights \gamma, for convex combination of stencils
constexpr Real dminus_k[NSUBSTENCIL] = {TWO_3RD, ONE_3RD};
constexpr Real dplus_k[NSUBSTENCIL] = {ONE_3RD, TWO_3RD};
// coefficient of smooth reference solution
constexpr int c_ref = 20.0;
} // namespace

//--------------------------------------------------------------------------------------
//! \fn Reconstruction::WenoThirdX1()
//  \brief Uniform Cartesian mesh WENO3 implementation. Using nonlinear weights from
// Mignone (2014)

void Reconstruction::WenoThirdX1(const int k, const int j, const int il, const int iu,
                                 const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
                                 AthenaArray<Real> &wl, AthenaArray<Real> &wr) {
  // set work arrays used for primitive/characterstic cell-averages to scratch
  AthenaArray<Real> &bx = scr01_i_, // 1D
                          // 2D:
                    &wc = scr1_ni_, // linearization pts for characteristic projection
                 &q_im1 = scr2_ni_,
                     &q = scr3_ni_, &q_ip1 = scr4_ni_;
  AthenaArray<Real> &qr_imh = scr5_ni_, &ql_iph = scr6_ni_;

  // cache the x1-sliced primitive states for eigensystem calculation
  for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      wc(n,i) = w(n,k,j,i);
      q    (n,i) = w(n,k,j,i  );
      //q_im2(n,i) = w(n,k,j,i-2);
      q_im1(n,i) = w(n,k,j,i-1);
      q_ip1(n,i) = w(n,k,j,i+1);
      //q_ip2(n,i) = w(n,k,j,i+2);
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      bx(i) = bcc(IB1,k,j,i);

      wc(IBY,i) = bcc(IB2,k,j,i);
      q    (IBY,i) = bcc(IB2,k,j,i  );
      //q_im2(IBY,i) = bcc(IB2,k,j,i-2);
      q_im1(IBY,i) = bcc(IB2,k,j,i-1);
      q_ip1(IBY,i) = bcc(IB2,k,j,i+1);
      //q_ip2(IBY,i) = bcc(IB2,k,j,i+2);

      wc(IBZ,i) = bcc(IB3,k,j,i);
      q    (IBZ,i) = bcc(IB3,k,j,i  );
      //q_im2(IBZ,i) = bcc(IB3,k,j,i-2);
      q_im1(IBZ,i) = bcc(IB3,k,j,i-1);
      q_ip1(IBZ,i) = bcc(IB3,k,j,i+1);
      //q_ip2(IBZ,i) = bcc(IB3,k,j,i+2);
    }
  }
  // Project cell-averages to characteristic variables, if necessary
  // Note order of characteristic fields in output vect corresponds to (IVX,IVY,IVZ)
  if (characteristic_projection) {
    //LeftEigenmatrixDotVector(IVX, il, iu, bx, wc, q_im2);
    LeftEigenmatrixDotVector(IVX, il, iu, bx, wc, q_im1);
    LeftEigenmatrixDotVector(IVX, il, iu, bx, wc, q);
    LeftEigenmatrixDotVector(IVX, il, iu, bx, wc, q_ip1);
    //LeftEigenmatrixDotVector(IVX, il, iu, bx, wc, q_ip2);
  }


  // WENO3 PARAMETERS
  // -----------------------SUBSTENCILS FOR INTERPOLATING INTERFACES-------------
  // Approximations to upper and lower interfaces per cell
  Real qminus_k[NSUBSTENCIL];
  Real qplus_k[NSUBSTENCIL]; //  = (number of stencils)

  // not indenxing another dimension of array by number of cell interfaces since it is
  // fixed =2 for all methods

  // (number of stencils) x (number of interpolating cells)
  // WENO3: 2 interfaces, 2x 2-cell stencils
  Real sminus_weight_k[NSUBSTENCIL][NSUBSTENCIL];
  Real splus_weight_k[NSUBSTENCIL][NSUBSTENCIL];
  // r=2 case of Table 1 in Jiang & Shu (1996)
  // see also Mignone (2014) eqs 28 and 29 (assume uniform mesh) and Yamaleev eq 15
  // omega_k^+ =
  // a_(k=0, l=0,1)^+ =
  splus_weight_k[0][0] = -0.5, splus_weight_k[0][1] = 1.5;
  // a_(k=1, l=0,1)^+ =
  splus_weight_k[1][0] = 0.5, splus_weight_k[1][1] = 0.5;
  //
  // omega_k^- =
  // (inferred from a^+_{k,l} by symmetry w.r.t. relative position of interface)
  sminus_weight_k[0][0] = 0.5, sminus_weight_k[0][1] = 0.5;
  sminus_weight_k[1][0] = 1.5, sminus_weight_k[1][1] = -0.5;

  // ----------------- NONLINEAR WEIGHTING OF SUBSTENCILS------------------------
  // Optimal linear weights \gamma, for convex combination of stencils
  // Real dminus_k[NSUBSTENCIL], dplus_k[NSUBSTENCIL]; // = gamma in Shu notation
  // dminus_k[0] = TWO_3RD, dminus_k[1] = ONE_3RD;
  // dplus_k[0] = ONE_3RD, dplus_k[1] = TWO_3RD;

  // Mesh spacing information
  Coordinates *pco = pmy_block_->pcoord;

  // Smoothness indicator per stencil (per cell per cycle)
  Real beta_k[NSUBSTENCIL];
  // Alpha nonlinear weight functions
  Real alpha_minus_k[NSUBSTENCIL], alpha_plus_k[NSUBSTENCIL];
  // Final omega nonlinear weight functions
  Real w_minus_k[NSUBSTENCIL], w_plus_k[NSUBSTENCIL];

  for (int n=0; n<NWAVE; ++n) {
    // Start outermost spatial cell loop at the first ghost cell to reconstruct w^L at
    // lowest real interface (extend to iu=ie+1 for upper-most real interaface)
    for (int i=il; i<=iu; ++i) {
      // Compute smoothness reference value for this cell and cycle
      // (does not depend on the specific variable stencil)
      Real h = pco->dx1v(i); // TODO(felker): consider replacing with direct div. by Nx1
      Real q_ref = c_ref*h*std::max(std::abs(q_im1(n,i)),
                                    std::max(abs(q(n,i)), std::abs(q_ip1(n,i))));

      // Sum of nonlinear alpha weights for normalization of nonlinear stencil weights
      Real alpha_minus_sum = 0.0;
      Real alpha_plus_sum = 0.0;

      // Loop over all stencils: compute both interface approximations and smoothness
      for (int sk=0; sk<NSUBSTENCIL; ++sk) {
        // First stencil starts at index i-(k-1) = i-k+1 cell index, final stencil starts
        // at i cell index
        int stencil_il = i-NSUBSTENCIL+1+sk;
        Real q_0 = 0.0;
        Real q_1 = 0.0;
        if (sk == 0) {
          q_0 = q_im1(n,i);
          q_1 = q(n,i);
        } else {
          q_0 = q(n,i);
          q_1 = q_ip1(n,i);
        }
        // ------- Compute reconstruction at interfaces for all substencils  --------
        // TODO(felker): use += operator w/ loop to generalize to any ncell size stencil
        // qminus_k[sk] = sminus_weight_k[sk][0]*q(n,stencil_il)
        //                + sminus_weight_k[sk][1]*q(n,stencil_il+1);
        // // Stencil shape and reconstructed polynomial are identical for omega_k^+/-,
        // // different coefficient weights are due to different positions
        // // of polynomial evaluation
        // qplus_k[sk] = splus_weight_k[sk][0]*q(n,stencil_il)
        //               + splus_weight_k[sk][1]*q(n,stencil_il+1);

        qminus_k[sk] = sminus_weight_k[sk][0]*q_0 + sminus_weight_k[sk][1]*q_1;
        qplus_k[sk] = splus_weight_k[sk][0]*q_0 + splus_weight_k[sk][1]*q_1;

        // ----- Compute substencil's smoothness indicator at this cell and cycle ----
        // (smoothness indicator is substencil-based, but equivalent for both interfaces)
        // WENO3 indicators: backward difference squared for lower stencil S_0,
        //                   forward difference squared for upper stencil S_1
        //beta_k[sk] = SQR(q(n,stencil_il) - q(n,stencil_il+1));
        beta_k[sk] = SQR(q_0 - q_1);
        // note: the flipped sign is squared

        //--- Compute nonlinear alpha weights for every stencil, ESWENO variation ----
        // AVOID DIVIDING BY ZERO!
        if ((beta_k[sk] + q_ref*q_ref) != 0.0) {
          alpha_minus_k[sk] = dminus_k[sk]*(
              1.0 + SQR(q_ip1(n,i) - 2.0*q(n,i) + q_im1(n,i)
                        )/(beta_k[sk] + SQR(q_ref)));
          alpha_plus_k[sk] = dplus_k[sk]*(
              1.0 + SQR(q_ip1(n,i) - 2.0*q(n,i) + q_im1(n,i)
                        )/(beta_k[sk] + SQR(q_ref)));
        } else {
          alpha_minus_k[sk] = 0.0;
          alpha_plus_k[sk] = 0.0;
        }
        // Initialize final omega nonlinear weights for every stencil
        w_minus_k[sk] = alpha_minus_k[sk];
        w_plus_k[sk] = alpha_plus_k[sk];
        // Build sums of alpha weights for normalization
        alpha_minus_sum += alpha_minus_k[sk];
        alpha_plus_sum += alpha_plus_k[sk];
      }
      // Initialize two Riemann states (on opposite interfaces)
      qr_imh(n,i) = 0.0;
      ql_iph(n,i) = 0.0;
      for (int sk=0; sk<NSUBSTENCIL; ++sk) {
        // Normalize omega nonlinear substencil weights at both interfaces
        // AVOID DIVIDING BY ZERO, specifically in the case of a trivial u=0 profile
        // (need to have alpha weights for all stencils)
        if (alpha_minus_sum != 0.0) {
          w_minus_k[sk] /= alpha_minus_sum;
        } else {
          w_minus_k[sk] = 0.0;
        }
        if (alpha_plus_sum != 0.0) {
          w_plus_k[sk] /= alpha_plus_sum;
        } else {
          w_plus_k[sk] = 0.0;
        }

        // Build nonlinear convex combination of approximations
        // equivalent to Q^-_i, Q^+_i, respectively, in Mignone notation
        qr_imh(n,i) += w_minus_k[sk]*qminus_k[sk];
        ql_iph(n,i) += w_plus_k[sk]*qplus_k[sk];
        // Interface states are automatically discontinuous due to cell-centering of
        // stencils. This also causes a mismatch of indices relative to Athena++
        // face-centering of ql,qr
      }
    }
  } // end loop over primitive variables

  // Project limited slope back to primitive variables, if necessary
  if (characteristic_projection) {
    RightEigenmatrixDotVector(IVX, il, iu, bx, wc, ql_iph);
    RightEigenmatrixDotVector(IVX, il, iu, bx, wc, qr_imh);
  }

  // compute ql_(i+1/2) and qr_(i-1/2)
  for (int n=0; n<NWAVE; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      wl(n,i+1) = ql_iph(n,i);
      wr(n,i  ) = qr_imh(n,i);
    }
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn Reconstruction::WenoThirdX2()
//  \brief

void Reconstruction::WenoThirdX2(const int k, const int j, const int il, const int iu,
                                 const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
                                 AthenaArray<Real> &wl, AthenaArray<Real> &wr) {
  return;
}

//--------------------------------------------------------------------------------------
//! \fn Reconstruction::WenoThirdX3()
//  \brief

void Reconstruction::WenoThirdX3(const int k, const int j, const int il, const int iu,
                                 const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
                                 AthenaArray<Real> &wl, AthenaArray<Real> &wr) {
  return;
}
