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

//--------------------------------------------------------------------------------------
//! \fn Reconstruction::WenoThirdX1()
//  \brief Uniform Cartesian mesh WENO3 implementation. Using nonlinear weights from
// Mignone (2014)

void Reconstruction::WenoThirdX1(const int k, const int j, const int il, const int iu,
                                 const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
                                 AthenaArray<Real> &wl, AthenaArray<Real> &wr) {
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
  // u_k^- =
  sminus_weight_k[0][0] = 0.5, sminus_weight_k[0][1] = 0.5;
  sminus_weight_k[1][0] = 1.5, sminus_weight_k[1][1] = -0.5;
  // u_k^+ =
  splus_weight_k[0][0] = -0.5, splus_weight_k[0][1] = 1.5;
  splus_weight_k[1][0] = 0.5, splus_weight_k[1][1] = 0.5;

  // ----------------- NONLINEAR WEIGHTING OF SUBSTENCILS------------------------
  // Optimal linear weights \gamma, for convex combination of stencils
  Real dminus_k[NSUBSTENCIL], dplus_k[NSUBSTENCIL]; // = gamma in Shu notation
  dminus_k[0] = TWO_3RD, dminus_k[1] = ONE_3RD;
  dplus_k[0] = ONE_3RD, dplus_k[1] = TWO_3RD;
  // coefficient of smooth reference solution
  constexpr int c_ref = 20.0;
  // Mesh spacing information
  Coordinates *pco = pmy_block_->pcoord;

  // Smoothness indicator per stencil (per cell per cycle)
  Real beta_k[NSUBSTENCIL];
  // Alpha nonlinear weight functions
  Real alpha_minus_k[NSUBSTENCIL], alpha_plus_k[NSUBSTENCIL];
  // Final omega nonlinear weight functions
  Real w_minus_k[NSUBSTENCIL], w_plus_k[NSUBSTENCIL];

  for (int n=0; n<NHYDRO; ++n) {
    // Start outermost spatial cell loop at the first ghost cell to reconstruct w^L at
    // lowest real interface
    for (int i=il-1; i<=iu; ++i) { // recall, iu=ie+1 for upper-most real interaface
      // Compute smoothness reference value for this cell and cycle
      // (does not depend on the specific stencil)
      Real h = pco->dx1v(i); // h = 1/N always?
      Real q_ref = c_ref*h*std::max(abs(q(n,k,j,i-1)),
                                    std::max(abs(q(n,k,j,i)), abs(q(n,k,j,i+1))));

      // Sum of nonlinear alpha weights for normalization of nonlinear stencil weights
      Real alpha_minus_sum = 0.0;
      Real alpha_plus_sum = 0.0;

      // Loop over all stencils: compute both interface approximations and smoothness
      for (int sk=0; sk<NSUBSTENCIL; ++sk) {
        // First stencil starts at index i-(k-1) = i-k+1 cell index, final stencil starts
        // at i cell index
        int stencil_il = i-NSUBSTENCIL+1+sk;
        // ------- Compute reconstruction at interfaces for all substencils  --------
        // TODO(felker): use += operator w/ loop to generalize to any ncell size stencil
        qminus_k[sk] = sminus_weight_k[sk][0]*q(n,k,j,stencil_il)
                       + sminus_weight_k[sk][1]*q(n,k,j,stencil_il+1);
        // Stencil shape and reconstructed polynomial are identical for u_k^+/-,
        // different coefficient weights are due to different positions
        // of polynomial evaluation
        qplus_k[sk] = splus_weight_k[sk][0]*q(n,k,j,stencil_il)
                      + splus_weight_k[sk][1]*q(n,k,j,stencil_il+1);
        // ----- Compute substencil's smoothness indicator at this cell and cycle ----
        // (smoothness indicator is substencil-based = same for both interfaces)
        // WENO3 indicators: backward difference squared for lower stencil S_0,
        //                   forward difference squared for upper stencil S_1
        beta_k[sk] = SQR(q(n,k,j,stencil_il) - q(n,k,j,stencil_il+1));
        // note: the flipped sign is squared

        //--- Compute nonlinear alpha weights for every stencil, ESWENO variation ----
        // AVOID DIVIDING BY ZERO!
        if ((beta_k[sk] + q_ref*q_ref) != 0.0) {
          alpha_minus_k[sk] = dminus_k[sk]*(
              1.0 + SQR(q(n,k,j,i+1) - 2.0*q(n,k,j,i) + q(n,k,j,i-1)
                        )/(beta_k[sk] + SQR(q_ref)));
          alpha_plus_k[sk] = dplus_k[sk]*(
              1.0 + SQR(q(n,k,j,i+1) - 2.0*q(n,k,j,i) + q(n,k,j,i-1)
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
      qr(n,i) = 0.0;
      ql(n,i+1) = 0.0;
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
        qr(n,i) += w_minus_k[sk]*qminus_k[sk]; // equivalent to Q^-_i in Mignone notation
        ql(n,i+1) += w_plus_k[sk]*qplus_k[sk]; // equivalent to Q^+_i in Mignone notation
        // Interface states are automatically discontinuous due to cell-centering of
        // stencils. This also causes a mismatch of indices relative to Athena++
        // face-centering of ql,qr
      }
    }
  } // end loop over primitive variables
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
