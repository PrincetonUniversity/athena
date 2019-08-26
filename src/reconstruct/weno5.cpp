//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file weno5.cpp
//  \brief fifth-order accurate weighted essentially non-oscillatory (WENO_5) scheme
//
// REFERENCES:
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

// number of stencils per cell; WENO5 consists of 3x three-point substencils per cell
#define NSUBSTENCIL 3

//--------------------------------------------------------------------------------------
//! \fn Reconstruction::WenoFifthX1()
//  \brief Uniform Cartesian mesh WENO5 implementation

void Reconstruction::WenoFifthX1(const int k, const int j, const int il, const int iu,
                                 const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
                                 AthenaArray<Real> &wl, AthenaArray<Real> &wr) {
  // WENO5 PARAMETERS
  //  -----------------------SUBSTENCILS FOR INTERPOLATING INTERFACES-------------
  // Approximations to upper and lower interfaces per cell
  Real qminus_k[NSUBSTENCIL];
  Real qplus_k[NSUBSTENCIL];

  // (number of stencils) x (number of interpolating cells)
  // WENO5: for both x1 interfaces, 3x 3-cell stencils
  Real sminus_weight_k[NSUBSTENCIL][NSUBSTENCIL];
  Real splus_weight_k[NSUBSTENCIL][NSUBSTENCIL];

  // Mirror symmetry relative to x_i cell center of entire reconstruction process:
  // Given reconstruction weights w+_k for all k stencils (cell-averages ordered from
  // lowest to highest u_j), to get w-_k, simply:
  // 1) Swap weight vectors/sets relative to central stencil (cyclic shift of weight
  // vectors). If k=odd, w+/-_{k/2+1} has same set of weights
  // 2) Flip weight vector (cyclic shift of u_j weights)

  // Never any sign flipping!
  // Or, just use my mignone_reconstructions Mathematica script to see the w+/- for all
  // il,ir combinations such that ir+il+1=k

  // u_k^- =
  sminus_weight_k[0][0] = -1.0/6.0, sminus_weight_k[0][1] = 5.0/6.0,
  sminus_weight_k[0][2] = ONE_3RD;
  sminus_weight_k[1][0] = ONE_3RD, sminus_weight_k[1][1] = 5.0/6.0,
  sminus_weight_k[1][2] = -1.0/6.0;
  sminus_weight_k[2][0] = 11.0/6.0, sminus_weight_k[2][1] = -7.0/6.0,
  sminus_weight_k[2][2] = ONE_3RD;
  // u_k^+ =
  splus_weight_k[0][0] = ONE_3RD, splus_weight_k[0][1] = -7.0/6.0,
  splus_weight_k[0][2] = 11.0/6.0;
  splus_weight_k[1][0] = -1.0/6.0,  splus_weight_k[1][1] = 5.0/6.0,
  splus_weight_k[1][2] = ONE_3RD;
  splus_weight_k[2][0] = ONE_3RD,  splus_weight_k[2][1] = 5.0/6.0,
  splus_weight_k[2][2] = -1.0/6.0;

  // ----------------- NONLINEAR WEIGHTING OF SUBSTENCILS------------------------
  // Optimal linear weights \gamma, for convex combination of stencils
  Real dminus_k[NSUBSTENCIL], dplus_k[NSUBSTENCIL]; // gamma in Shu notation
  dminus_k[0] = 3.0/10.0, dminus_k[1] = 6.0/10.0, dminus_k[2] = 1.0/10.0;
  dplus_k[0] = 1.0/10.0, dplus_k[1] = 6.0/10.0, dplus_k[2] = 3.0/10;

  // Smoothness indicator per substencil (per cell per cycle)
  Real beta_k[NSUBSTENCIL];
  Real beta_w_k[NSUBSTENCIL][NSUBSTENCIL];
  // per cell weights for second squared term in each k stencil's smoothness indicator
  beta_w_k[0][0] = 1.0, beta_w_k[0][1] = -4.0, beta_w_k[0][2] = 3.0;
  beta_w_k[1][0] = 1.0, beta_w_k[1][1] = 0.0, beta_w_k[1][2] = -1.0;
  beta_w_k[2][0] = 3.0, beta_w_k[2][1] = -4.0, beta_w_k[2][2] = 1.0;

  // Final omega nonlinear weight functions
  Real w_minus_k[NSUBSTENCIL], w_plus_k[NSUBSTENCIL];
  Real epsilon = 1e-6;
  // Should look for improved WENO5 nonlinear weight functions, epsilon() values that
  // depend on relative size of variable, etc.

  for (int n=0; n<NHYDRO; ++n) {
    // Start outermost spatial cell loop at the first ghost cell to reconstruct w^L at
    // lowest real interface
    for (int i=il-1; i<=iu; ++i) { // recall, iu=ie+1 for upper-most real interaface
      // Sum of nonlinear alpha weights for normalization of nonlinear stencil weights
      Real w_minus_sum =0.0;
      Real w_plus_sum =0.0;

      // Loop over all stencils: compute both interface approximations and smoothness
      for (int sk=0; sk<NSUBSTENCIL; ++sk) {
        // First stencil starts at index i-(k-1) = i-k+1 cell index, final stencil starts
        // at i cell index

        int stencil_il = i-NSUBSTENCIL+1+sk;
        // ------- Compute reconstruction at interfaces for all substencils  --------//
        qminus_k[sk] = sminus_weight_k[sk][0]*q(n,k,j,stencil_il)
                       + sminus_weight_k[sk][1]*q(n,k,j,stencil_il+1)
                       + sminus_weight_k[sk][2]*q(n,k,j,stencil_il+2);
        // Stencil shape and reconstructed polynomial are identical for u_k^+/-, different
        // coefficient weights are due to different positions of evaluation of polynomial
        qplus_k[sk] = splus_weight_k[sk][0]*q(n,k,j,stencil_il)
                      + splus_weight_k[sk][1]*q(n,k,j,stencil_il+1)
                      + splus_weight_k[sk][2]*q(n,k,j,stencil_il+2);
        // ---- Compute substencil's smoothness indicator at this cell and cycle ----//
        // (smoothness indicator is stencil-based = same for both interfaces)
        beta_k[sk] = 13.0/12.0*SQR(
            q(n,k,j,stencil_il) - 2.0*q(n,k,j,stencil_il+1) + q(n,k,j,stencil_il+2)
                                   )
                     + 0.25*SQR(beta_w_k[sk][0]*q(n,k,j,stencil_il)
                                + beta_w_k[sk][1]*q(n,k,j,stencil_il+1)
                                + beta_w_k[sk][2]*q(n,k,j,stencil_il+2)
                                ); // no real need to store NSUBSTENCIL indicator values
        // ---- Compute traditional WENO5/WENO3 nonlinear weight functions ---------//
        w_minus_k[sk] = dminus_k[sk]/SQR(epsilon + beta_k[sk]);
        w_plus_k[sk] = dplus_k[sk]/SQR(epsilon + beta_k[sk]);
        // Build sums of weights for normalization
        w_minus_sum += w_minus_k[sk];
        w_plus_sum += w_plus_k[sk];
      }
      // Initialize two Riemann states (on opposite interfaces)
      qr(n,i) = 0.0;
      ql(n,i+1) = 0.0;
      for (int sk=0; sk<NSUBSTENCIL; ++sk) {
        // Normalize omega nonlinear substencil weights at both interfaces
        // AVOID DIVIDING BY ZERO
        // Should be impossible in traditional WENO nonlinear weighting due to epsilon
        if (w_minus_sum != 0.0) {
          w_minus_k[sk] /= w_minus_sum;
        } else {
          w_minus_k[sk] = 0.0;
        }
        if (w_plus_sum != 0.0) {
          w_plus_k[sk] /= w_plus_sum;
        } else {
          w_plus_k[sk] = 0.0;
        }

        // Build nonlinear convex combination of approximations
        qr(n,i) += w_minus_k[sk]*qminus_k[sk]; // equivalent to Q^-_i in Mignone notation
        ql(n,i+1) += w_plus_k[sk]*qplus_k[sk]; // equivalent to Q^+_i in Mignone notation
      }
    }
  } // end loop over primitive variables
  return;
}

//--------------------------------------------------------------------------------------
//! \fn Reconstruction::WenoFifthX2()
//  \brief Uniform Cartesian mesh WENO5 implementation

void Reconstruction::WenoFifthX2(const int k, const int j, const int il, const int iu,
                                 const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
                                 AthenaArray<Real> &wl, AthenaArray<Real> &wr) {
  return;
}

//--------------------------------------------------------------------------------------
//! \fn Reconstruction::WenoFifthX3()
//  \brief Uniform Cartesian mesh WENO5 implementation

void Reconstruction::WenoFifthX3(const int k, const int j, const int il, const int iu,
                                 const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
                                 AthenaArray<Real> &wl, AthenaArray<Real> &wr) {
  return;
}
