//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file radiation.cpp
//  \brief implementation of flux-related functions in class Radiation

// C++ headers
#include <algorithm>  // max, min
#include <cmath>      // expm1, sqrt

// Athena++ headers
#include "radiation.hpp"
#include "../athena.hpp"                   // Real, indices
#include "../athena_arrays.hpp"            // AthenaArray
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../eos/eos.hpp"                  // EquationOfState
#include "../hydro/hydro.hpp"              // Hydro
#include "../mesh/mesh.hpp"                // MeshBlock

//----------------------------------------------------------------------------------------
// Function for calculating spatial and angular fluxes
// Inputs:
//   prim_rad: primitive intensity
//   prim_hydro: primitive hydrodynamic variables
//   order: reconstruction order
//   dt: timestep
// Outputs:
//   this->flux_x, this->flux_a: fluxes set

void Radiation::CalculateFluxes(AthenaArray<Real> &prim_rad,
    const AthenaArray<Real> &prim_hydro, int order, Real dt) {

  // Calculate x1-fluxes
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {

      // Reconstruct radiation
      if (order == 1) {
        RadiationDonorCellX1(prim_rad, k, j);
      } else {
        RadiationPiecewiseLinearX1(prim_rad, k, j);
      }

      // Calculate pure radiation fluxes
      if (not (coupled_to_matter and tau_multiplier > 0.0)) {
        for (int l = zs; l <= ze; ++l) {
          for (int m = ps; m <= pe; ++m) {
            int lm = AngleInd(l, m);
            for (int i = is; i <= ie+1; ++i) {
              Real n1 = n1_1_(l,m,k,j,i);
              Real n1_n_0 = n1_n_mu_(0,l,m,k,j,i);
              flux_x[X1DIR](lm,k,j,i) =
                  n1_n_0 * (n1 > 0.0 ? ii_l_(lm,i) : ii_r_(lm,i));
            }
          }
        }
      }

      // Calculate radiation-hydrodynamic fluxes
      else {

        // Calculate negative fluid velocity in tetrad frame
        pmy_block->pcoord->Face1Metric(k, j, is, ie + 1, g_, gi_);
        for (int i = is; i <= ie+1; ++i) {
          Real dx_l = pmy_block->pcoord->x1f(i) - pmy_block->pcoord->x1v(i-1);
          Real dx_r = pmy_block->pcoord->x1v(i) - pmy_block->pcoord->x1f(i);
          Real uu1 = (dx_r * prim_hydro(IVX,k,j,i-1) + dx_l * prim_hydro(IVX,k,j,i))
              / (dx_l + dx_r);
          Real uu2 = (dx_r * prim_hydro(IVY,k,j,i-1) + dx_l * prim_hydro(IVY,k,j,i))
              / (dx_l + dx_r);
          Real uu3 = (dx_r * prim_hydro(IVZ,k,j,i-1) + dx_l * prim_hydro(IVZ,k,j,i))
              / (dx_l + dx_r);
          Real temp_var = g_(I11,i) * SQR(uu1) + 2.0 * g_(I12,i) * uu1 * uu2
              + 2.0 * g_(I13,i) * uu1 * uu3 + g_(I22,i) * SQR(uu2)
              + 2.0 * g_(I23,i) * uu2 * uu3 + g_(I33,i) * SQR(uu3);
          Real uu0 = std::sqrt(1.0 + temp_var);
          u_tet_(0,i) = norm_to_tet_1_(0,0,k,j,i) * uu0 + norm_to_tet_1_(0,1,k,j,i) * uu1
              + norm_to_tet_1_(0,2,k,j,i) * uu2 + norm_to_tet_1_(0,3,k,j,i) * uu3;
        }

        // Calculate fluid velocity in tetrad frame to left of interface
        pmy_block->pcoord->CellMetric(k, j, is - 1, ie + 1, g_, gi_);
        for (int i = is; i <= ie+1; ++i) {
          Real uu1 = prim_hydro(IVX,k,j,i-1);
          Real uu2 = prim_hydro(IVY,k,j,i-1);
          Real uu3 = prim_hydro(IVZ,k,j,i-1);
          Real temp_var = g_(I11,i) * SQR(uu1) + 2.0 * g_(I12,i) * uu1 * uu2
              + 2.0 * g_(I13,i) * uu1 * uu3 + g_(I22,i) * SQR(uu2)
              + 2.0 * g_(I23,i) * uu2 * uu3 + g_(I33,i) * SQR(uu3);
          Real uu0 = std::sqrt(1.0 + temp_var);
          u_tet_l_(0,i) = norm_to_tet_(0,0,k,j,i-1) * uu0
              + norm_to_tet_(0,1,k,j,i-1) * uu1 + norm_to_tet_(0,2,k,j,i-1) * uu2
              + norm_to_tet_(0,3,k,j,i-1) * uu3;
          u_tet_l_(1,i) = -norm_to_tet_(1,0,k,j,i-1) * uu0
              - norm_to_tet_(1,1,k,j,i-1) * uu1 - norm_to_tet_(1,2,k,j,i-1) * uu2
              - norm_to_tet_(1,3,k,j,i-1) * uu3;
          u_tet_l_(2,i) = -norm_to_tet_(2,0,k,j,i-1) * uu0
              - norm_to_tet_(2,1,k,j,i-1) * uu1 - norm_to_tet_(2,2,k,j,i-1) * uu2
              - norm_to_tet_(2,3,k,j,i-1) * uu3;
          u_tet_l_(3,i) = -norm_to_tet_(3,0,k,j,i-1) * uu0
              - norm_to_tet_(3,1,k,j,i-1) * uu1 - norm_to_tet_(3,2,k,j,i-1) * uu2
              - norm_to_tet_(3,3,k,j,i-1) * uu3;
        }

        // Calculate fluid velocity in tetrad frame to right of interface
        for (int i = is; i <= ie+1; ++i) {
          Real uu1 = prim_hydro(IVX,k,j,i);
          Real uu2 = prim_hydro(IVY,k,j,i);
          Real uu3 = prim_hydro(IVZ,k,j,i);
          Real temp_var = g_(I11,i) * SQR(uu1) + 2.0 * g_(I12,i) * uu1 * uu2
              + 2.0 * g_(I13,i) * uu1 * uu3 + g_(I22,i) * SQR(uu2)
              + 2.0 * g_(I23,i) * uu2 * uu3 + g_(I33,i) * SQR(uu3);
          Real uu0 = std::sqrt(1.0 + temp_var);
          u_tet_r_(0,i) = norm_to_tet_(0,0,k,j,i) * uu0 + norm_to_tet_(0,1,k,j,i) * uu1
              + norm_to_tet_(0,2,k,j,i) * uu2 + norm_to_tet_(0,3,k,j,i) * uu3;
          u_tet_r_(1,i) = -norm_to_tet_(1,0,k,j,i) * uu0 - norm_to_tet_(1,1,k,j,i) * uu1
              - norm_to_tet_(1,2,k,j,i) * uu2 - norm_to_tet_(1,3,k,j,i) * uu3;
          u_tet_r_(2,i) = -norm_to_tet_(2,0,k,j,i) * uu0 - norm_to_tet_(2,1,k,j,i) * uu1
              - norm_to_tet_(2,2,k,j,i) * uu2 - norm_to_tet_(2,3,k,j,i) * uu3;
          u_tet_r_(3,i) = -norm_to_tet_(3,0,k,j,i) * uu0 - norm_to_tet_(3,1,k,j,i) * uu1
              - norm_to_tet_(3,2,k,j,i) * uu2 - norm_to_tet_(3,3,k,j,i) * uu3;
        }

        // Calculate optical depth
        for (int i = is; i <= ie+1; ++i) {
          Real dx_l = pmy_block->pcoord->x1f(i) - pmy_block->pcoord->x1v(i-1);
          Real dx_r = pmy_block->pcoord->x1v(i) - pmy_block->pcoord->x1f(i);
          Real k_a = (dx_r * opacity(OPAA,k,j,i-1) * prim_hydro(IDN,k,j,i-1)
              + dx_l * opacity(OPAA,k,j,i) * prim_hydro(IDN,k,j,i)) / (dx_l + dx_r);
          Real k_s = (dx_r * opacity(OPAS,k,j,i-1) * prim_hydro(IDN,k,j,i-1)
              + dx_l * opacity(OPAS,k,j,i) * prim_hydro(IDN,k,j,i)) / (dx_l + dx_r);
          Real k_tot = k_a + k_s;
          Real tau = tau_multiplier * k_tot * dt / u_tet_(0,i);
          tau_factor_(i) = -std::expm1(-tau * tau);
        }

        // Calculate modified fluxes
        for (int l = zs; l <= ze; ++l) {
          for (int m = ps; m <= pe; ++m) {
            int lm = AngleInd(l, m);
            for (int i = is; i <= ie+1; ++i) {
              Real u_n_left = -u_tet_l_(0,i) * nh_cc_(0,l,m)
                  + u_tet_l_(1,i) * nh_cc_(1,l,m) + u_tet_l_(2,i) * nh_cc_(2,l,m)
                  + u_tet_l_(3,i) * nh_cc_(3,l,m);
              Real u_n_right = -u_tet_r_(0,i) * nh_cc_(0,l,m)
                  + u_tet_r_(1,i) * nh_cc_(1,l,m) + u_tet_r_(2,i) * nh_cc_(2,l,m)
                  + u_tet_r_(3,i) * nh_cc_(3,l,m);
              Real advection_factor_left =
                  -tau_factor_(i) * (1.0 + 1.0 / (u_n_left * u_n_left * u_n_left));
              Real advection_factor_right =
                  -tau_factor_(i) * (1.0 + 1.0 / (u_n_right * u_n_right * u_n_right));
              Real ii_1_left = ii_l_(lm,i) * (1.0 - advection_factor_left);
              Real ii_1_right = ii_r_(lm,i) * (1.0 - advection_factor_right);
              Real ii_2_left = advection_factor_left;
              Real ii_2_right = advection_factor_right;
              Real v_1 = n1_n_mu_(0,l,m,k,j,i);
              Real dx_l = pmy_block->pcoord->x1f(i) - pmy_block->pcoord->x1v(i-1);
              Real dx_r = pmy_block->pcoord->x1v(i) - pmy_block->pcoord->x1f(i);
              Real v_2 = v_1 * (dx_r * ii_2_left + dx_l * ii_2_right) / (dx_l + dx_r);
              Real n1 = n1_1_(l,m,k,j,i);
              flux_x[X1DIR](lm,k,j,i) = v_1 * (n1 > 0.0 ? ii_1_left : ii_1_right);
              flux_x[X1DIR](lm,k,j,i) += v_2 * (n1 > 0.0 ? ii_l_(lm,i) : ii_r_(lm,i));
            }
          }
        }
      }
    }
  }

  // Calculate x2-fluxes
  if (js != je) {
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je+1; ++j) {

        // Reconstruct radiation
        if (order == 1) {
          RadiationDonorCellX2(prim_rad, k, j);
        } else {
          RadiationPiecewiseLinearX2(prim_rad, k, j);
        }

        // Calculate pure radiation fluxes
        if (not (coupled_to_matter and tau_multiplier > 0.0)) {
          for (int l = zs; l <= ze; ++l) {
            for (int m = ps; m <= pe; ++m) {
              int lm = AngleInd(l, m);
              for (int i = is; i <= ie; ++i) {
                Real n2 = n2_2_(l,m,k,j,i);
                Real n2_n_0 = n2_n_mu_(0,l,m,k,j,i);
                flux_x[X2DIR](lm,k,j,i) =
                    n2_n_0 * (n2 > 0.0 ? ii_l_(lm,i) : ii_r_(lm,i));
              }
            }
          }
        }

        // Calculate radiation-hydrodynamic fluxes
        else {

          // Calculate negative fluid velocity in tetrad frame
          pmy_block->pcoord->Face2Metric(k, j, is, ie, g_, gi_);
          Real dx_l = pmy_block->pcoord->x2f(j) - pmy_block->pcoord->x2v(j-1);
          Real dx_r = pmy_block->pcoord->x2v(j) - pmy_block->pcoord->x2f(j);
          for (int i = is; i <= ie; ++i) {
            Real uu1 = (dx_r * prim_hydro(IVX,k,j-1,i) + dx_l * prim_hydro(IVX,k,j,i))
                / (dx_l + dx_r);
            Real uu2 = (dx_r * prim_hydro(IVY,k,j-1,i) + dx_l * prim_hydro(IVY,k,j,i))
                / (dx_l + dx_r);
            Real uu3 = (dx_r * prim_hydro(IVZ,k,j-1,i) + dx_l * prim_hydro(IVZ,k,j,i))
                / (dx_l + dx_r);
            Real temp_var = g_(I11,i) * SQR(uu1) + 2.0 * g_(I12,i) * uu1 * uu2
                + 2.0 * g_(I13,i) * uu1 * uu3 + g_(I22,i) * SQR(uu2)
                + 2.0 * g_(I23,i) * uu2 * uu3 + g_(I33,i) * SQR(uu3);
            Real uu0 = std::sqrt(1.0 + temp_var);
            u_tet_(0,i) = norm_to_tet_2_(0,0,k,j,i) * uu0
                + norm_to_tet_2_(0,1,k,j,i) * uu1 + norm_to_tet_2_(0,2,k,j,i) * uu2
                + norm_to_tet_2_(0,3,k,j,i) * uu3;
          }

          // Calculate fluid velocity in tetrad frame to left of interface
          pmy_block->pcoord->CellMetric(k, j - 1, is, ie, g_, gi_);
          for (int i = is; i <= ie; ++i) {
            Real uu1 = prim_hydro(IVX,k,j-1,i);
            Real uu2 = prim_hydro(IVY,k,j-1,i);
            Real uu3 = prim_hydro(IVZ,k,j-1,i);
            Real temp_var = g_(I11,i) * SQR(uu1) + 2.0 * g_(I12,i) * uu1 * uu2
                + 2.0 * g_(I13,i) * uu1 * uu3 + g_(I22,i) * SQR(uu2)
                + 2.0 * g_(I23,i) * uu2 * uu3 + g_(I33,i) * SQR(uu3);
            Real uu0 = std::sqrt(1.0 + temp_var);
            u_tet_l_(0,i) = norm_to_tet_(0,0,k,j-1,i) * uu0
                + norm_to_tet_(0,1,k,j-1,i) * uu1 + norm_to_tet_(0,2,k,j-1,i) * uu2
                + norm_to_tet_(0,3,k,j-1,i) * uu3;
            u_tet_l_(1,i) = -norm_to_tet_(1,0,k,j-1,i) * uu0
                - norm_to_tet_(1,1,k,j-1,i) * uu1 - norm_to_tet_(1,2,k,j-1,i) * uu2
                - norm_to_tet_(1,3,k,j-1,i) * uu3;
            u_tet_l_(2,i) = -norm_to_tet_(2,0,k,j-1,i) * uu0
                - norm_to_tet_(2,1,k,j-1,i) * uu1 - norm_to_tet_(2,2,k,j-1,i) * uu2
                - norm_to_tet_(2,3,k,j-1,i) * uu3;
            u_tet_l_(3,i) = -norm_to_tet_(3,0,k,j-1,i) * uu0
                - norm_to_tet_(3,1,k,j-1,i) * uu1 - norm_to_tet_(3,2,k,j-1,i) * uu2
                - norm_to_tet_(3,3,k,j-1,i) * uu3;
          }

          // Calculate fluid velocity in tetrad frame to right of interface
          pmy_block->pcoord->CellMetric(k, j, is, ie, g_, gi_);
          for (int i = is; i <= ie; ++i) {
            Real uu1 = prim_hydro(IVX,k,j,i);
            Real uu2 = prim_hydro(IVY,k,j,i);
            Real uu3 = prim_hydro(IVZ,k,j,i);
            Real temp_var = g_(I11,i) * SQR(uu1) + 2.0 * g_(I12,i) * uu1 * uu2
                + 2.0 * g_(I13,i) * uu1 * uu3 + g_(I22,i) * SQR(uu2)
                + 2.0 * g_(I23,i) * uu2 * uu3 + g_(I33,i) * SQR(uu3);
            Real uu0 = std::sqrt(1.0 + temp_var);
            u_tet_r_(0,i) = norm_to_tet_(0,0,k,j,i) * uu0 + norm_to_tet_(0,1,k,j,i) * uu1
                + norm_to_tet_(0,2,k,j,i) * uu2 + norm_to_tet_(0,3,k,j,i) * uu3;
            u_tet_r_(1,i) = -norm_to_tet_(1,0,k,j,i) * uu0 - norm_to_tet_(1,1,k,j,i) * uu1
                - norm_to_tet_(1,2,k,j,i) * uu2 - norm_to_tet_(1,3,k,j,i) * uu3;
            u_tet_r_(2,i) = -norm_to_tet_(2,0,k,j,i) * uu0 - norm_to_tet_(2,1,k,j,i) * uu1
                - norm_to_tet_(2,2,k,j,i) * uu2 - norm_to_tet_(2,3,k,j,i) * uu3;
            u_tet_r_(3,i) = -norm_to_tet_(3,0,k,j,i) * uu0 - norm_to_tet_(3,1,k,j,i) * uu1
                - norm_to_tet_(3,2,k,j,i) * uu2 - norm_to_tet_(3,3,k,j,i) * uu3;
          }

          // Calculate optical depth
          for (int i = is; i <= ie; ++i) {
            Real k_a = (dx_r * opacity(OPAA,k,j-1,i) * prim_hydro(IDN,k,j-1,i)
                + dx_l * opacity(OPAA,k,j,i) * prim_hydro(IDN,k,j,i)) / (dx_l + dx_r);
            Real k_s = (dx_r * opacity(OPAS,k,j-1,i) * prim_hydro(IDN,k,j-1,i)
                + dx_l * opacity(OPAS,k,j,i) * prim_hydro(IDN,k,j,i)) / (dx_l + dx_r);
            Real k_tot = k_a + k_s;
            Real tau = tau_multiplier * k_tot * dt / u_tet_(0,i);
            tau_factor_(i) = -std::expm1(-tau * tau);
          }

          // Calculate modified fluxes
          for (int l = zs; l <= ze; ++l) {
            for (int m = ps; m <= pe; ++m) {
              int lm = AngleInd(l, m);
              for (int i = is; i <= ie; ++i) {
                Real u_n_left = -u_tet_l_(0,i) * nh_cc_(0,l,m)
                    + u_tet_l_(1,i) * nh_cc_(1,l,m) + u_tet_l_(2,i) * nh_cc_(2,l,m)
                    + u_tet_l_(3,i) * nh_cc_(3,l,m);
                Real u_n_right = -u_tet_r_(0,i) * nh_cc_(0,l,m)
                    + u_tet_r_(1,i) * nh_cc_(1,l,m) + u_tet_r_(2,i) * nh_cc_(2,l,m)
                    + u_tet_r_(3,i) * nh_cc_(3,l,m);
                Real advection_factor_left =
                    -tau_factor_(i) * (1.0 + 1.0 / (u_n_left * u_n_left * u_n_left));
                Real advection_factor_right =
                    -tau_factor_(i) * (1.0 + 1.0 / (u_n_right * u_n_right * u_n_right));
                Real ii_1_left = ii_l_(lm,i) * (1.0 - advection_factor_left);
                Real ii_1_right = ii_r_(lm,i) * (1.0 - advection_factor_right);
                Real ii_2_left = advection_factor_left;
                Real ii_2_right = advection_factor_right;
                Real v_1 = n2_n_mu_(0,l,m,k,j,i);
                Real v_2 = v_1 * (dx_r * ii_2_left + dx_l * ii_2_right) / (dx_l + dx_r);
                Real n2 = n2_2_(l,m,k,j,i);
                flux_x[X2DIR](lm,k,j,i) = v_1 * (n2 > 0.0 ? ii_1_left : ii_1_right);
                flux_x[X2DIR](lm,k,j,i) += v_2 * (n2 > 0.0 ? ii_l_(lm,i) : ii_r_(lm,i));
              }
            }
          }
        }
      }
    }
  }

  // Calculate x3-fluxes
  if (ks != ke) {
    for (int k = ks; k <= ke+1; ++k) {
      for (int j = js; j <= je; ++j) {

        // Reconstruct radiation
        if (order == 1) {
          RadiationDonorCellX3(prim_rad, k, j);
        } else {
          RadiationPiecewiseLinearX3(prim_rad, k, j);
        }

        // Calculate pure radiation fluxes
        if (not (coupled_to_matter and tau_multiplier > 0.0)) {
          for (int l = zs; l <= ze; ++l) {
            for (int m = ps; m <= pe; ++m) {
              int lm = AngleInd(l, m);
              for (int i = is; i <= ie; ++i) {
                Real n3 = n3_3_(l,m,k,j,i);
                Real n3_n_0 = n3_n_mu_(0,l,m,k,j,i);
                flux_x[X3DIR](lm,k,j,i) =
                    n3_n_0 * (n3 > 0.0 ? ii_l_(lm,i) : ii_r_(lm,i));
              }
            }
          }
        }

        // Calculate radiation-hydrodynamic fluxes
        else {

          // Calculate negative fluid velocity in tetrad frame
          pmy_block->pcoord->Face3Metric(k, j, is, ie, g_, gi_);
          Real dx_l = pmy_block->pcoord->x3f(k) - pmy_block->pcoord->x3v(k-1);
          Real dx_r = pmy_block->pcoord->x3v(k) - pmy_block->pcoord->x3f(k);
          for (int i = is; i <= ie; ++i) {
            Real uu1 = (dx_r * prim_hydro(IVX,k-1,j,i) + dx_l * prim_hydro(IVX,k,j,i))
                / (dx_l + dx_r);
            Real uu2 = (dx_r * prim_hydro(IVY,k-1,j,i) + dx_l * prim_hydro(IVY,k,j,i))
                / (dx_l + dx_r);
            Real uu3 = (dx_r * prim_hydro(IVZ,k-1,j,i) + dx_l * prim_hydro(IVZ,k,j,i))
                / (dx_l + dx_r);
            Real temp_var = g_(I11,i) * SQR(uu1) + 2.0 * g_(I12,i) * uu1 * uu2
                + 2.0 * g_(I13,i) * uu1 * uu3 + g_(I22,i) * SQR(uu2)
                + 2.0 * g_(I23,i) * uu2 * uu3 + g_(I33,i) * SQR(uu3);
            Real uu0 = std::sqrt(1.0 + temp_var);
            u_tet_(0,i) = norm_to_tet_3_(0,0,k,j,i) * uu0
                + norm_to_tet_3_(0,1,k,j,i) * uu1 + norm_to_tet_3_(0,2,k,j,i) * uu2
                + norm_to_tet_3_(0,3,k,j,i) * uu3;
          }

          // Calculate fluid velocity in tetrad frame to left of interface
          pmy_block->pcoord->CellMetric(k - 1, j, is, ie, g_, gi_);
          for (int i = is; i <= ie; ++i) {
            Real uu1 = prim_hydro(IVX,k-1,j,i);
            Real uu2 = prim_hydro(IVY,k-1,j,i);
            Real uu3 = prim_hydro(IVZ,k-1,j,i);
            Real temp_var = g_(I11,i) * SQR(uu1) + 2.0 * g_(I12,i) * uu1 * uu2
                + 2.0 * g_(I13,i) * uu1 * uu3 + g_(I22,i) * SQR(uu2)
                + 2.0 * g_(I23,i) * uu2 * uu3 + g_(I33,i) * SQR(uu3);
            Real uu0 = std::sqrt(1.0 + temp_var);
            u_tet_l_(0,i) = norm_to_tet_(0,0,k-1,j,i) * uu0
                + norm_to_tet_(0,1,k-1,j,i) * uu1 + norm_to_tet_(0,2,k-1,j,i) * uu2
                + norm_to_tet_(0,3,k-1,j,i) * uu3;
            u_tet_l_(1,i) = -norm_to_tet_(1,0,k-1,j,i) * uu0
                - norm_to_tet_(1,1,k-1,j,i) * uu1 - norm_to_tet_(1,2,k-1,j,i) * uu2
                - norm_to_tet_(1,3,k-1,j,i) * uu3;
            u_tet_l_(2,i) = -norm_to_tet_(2,0,k-1,j,i) * uu0
                - norm_to_tet_(2,1,k-1,j,i) * uu1 - norm_to_tet_(2,2,k-1,j,i) * uu2
                - norm_to_tet_(2,3,k-1,j,i) * uu3;
            u_tet_l_(3,i) = -norm_to_tet_(3,0,k-1,j,i) * uu0
                - norm_to_tet_(3,1,k-1,j,i) * uu1 - norm_to_tet_(3,2,k-1,j,i) * uu2
                - norm_to_tet_(3,3,k-1,j,i) * uu3;
          }

          // Calculate fluid velocity in tetrad frame to right of interface
          pmy_block->pcoord->CellMetric(k, j, is, ie, g_, gi_);
          for (int i = is; i <= ie; ++i) {
            Real uu1 = prim_hydro(IVX,k,j,i);
            Real uu2 = prim_hydro(IVY,k,j,i);
            Real uu3 = prim_hydro(IVZ,k,j,i);
            Real temp_var = g_(I11,i) * SQR(uu1) + 2.0 * g_(I12,i) * uu1 * uu2
                + 2.0 * g_(I13,i) * uu1 * uu3 + g_(I22,i) * SQR(uu2)
                + 2.0 * g_(I23,i) * uu2 * uu3 + g_(I33,i) * SQR(uu3);
            Real uu0 = std::sqrt(1.0 + temp_var);
            u_tet_r_(0,i) = norm_to_tet_(0,0,k,j,i) * uu0 + norm_to_tet_(0,1,k,j,i) * uu1
                + norm_to_tet_(0,2,k,j,i) * uu2 + norm_to_tet_(0,3,k,j,i) * uu3;
            u_tet_r_(1,i) = -norm_to_tet_(1,0,k,j,i) * uu0 - norm_to_tet_(1,1,k,j,i) * uu1
                - norm_to_tet_(1,2,k,j,i) * uu2 - norm_to_tet_(1,3,k,j,i) * uu3;
            u_tet_r_(2,i) = -norm_to_tet_(2,0,k,j,i) * uu0 - norm_to_tet_(2,1,k,j,i) * uu1
                - norm_to_tet_(2,2,k,j,i) * uu2 - norm_to_tet_(2,3,k,j,i) * uu3;
            u_tet_r_(3,i) = -norm_to_tet_(3,0,k,j,i) * uu0 - norm_to_tet_(3,1,k,j,i) * uu1
                - norm_to_tet_(3,2,k,j,i) * uu2 - norm_to_tet_(3,3,k,j,i) * uu3;
          }

          // Calculate optical depth
          for (int i = is; i <= ie; ++i) {
            Real k_a = (dx_r * opacity(OPAA,k-1,j,i) * prim_hydro(IDN,k-1,j,i)
                + dx_l * opacity(OPAA,k,j,i) * prim_hydro(IDN,k,j,i)) / (dx_l + dx_r);
            Real k_s = (dx_r * opacity(OPAS,k-1,j,i) * prim_hydro(IDN,k-1,j,i)
                + dx_l * opacity(OPAS,k,j,i) * prim_hydro(IDN,k,j,i)) / (dx_l + dx_r);
            Real k_tot = k_a + k_s;
            Real tau = tau_multiplier * k_tot * dt / u_tet_(0,i);
            tau_factor_(i) = -std::expm1(-tau * tau);
          }

          // Calculate modified fluxes
          for (int l = zs; l <= ze; ++l) {
            for (int m = ps; m <= pe; ++m) {
              int lm = AngleInd(l, m);
              for (int i = is; i <= ie; ++i) {
                Real u_n_left = -u_tet_l_(0,i) * nh_cc_(0,l,m)
                    + u_tet_l_(1,i) * nh_cc_(1,l,m) + u_tet_l_(2,i) * nh_cc_(2,l,m)
                    + u_tet_l_(3,i) * nh_cc_(3,l,m);
                Real u_n_right = -u_tet_r_(0,i) * nh_cc_(0,l,m)
                    + u_tet_r_(1,i) * nh_cc_(1,l,m) + u_tet_r_(2,i) * nh_cc_(2,l,m)
                    + u_tet_r_(3,i) * nh_cc_(3,l,m);
                Real advection_factor_left =
                    -tau_factor_(i) * (1.0 + 1.0 / (u_n_left * u_n_left * u_n_left));
                Real advection_factor_right =
                    -tau_factor_(i) * (1.0 + 1.0 / (u_n_right * u_n_right * u_n_right));
                Real ii_1_left = ii_l_(lm,i) * (1.0 - advection_factor_left);
                Real ii_1_right = ii_r_(lm,i) * (1.0 - advection_factor_right);
                Real ii_2_left = advection_factor_left;
                Real ii_2_right = advection_factor_right;
                Real v_1 = n3_n_mu_(0,l,m,k,j,i);
                Real v_2 = v_1 * (dx_r * ii_2_left + dx_l * ii_2_right) / (dx_l + dx_r);
                Real n3 = n3_3_(l,m,k,j,i);
                flux_x[X3DIR](lm,k,j,i) = v_1 * (n3 > 0.0 ? ii_1_left : ii_1_right);
                flux_x[X3DIR](lm,k,j,i) += v_2 * (n3 > 0.0 ? ii_l_(lm,i) : ii_r_(lm,i));
              }
            }
          }
        }
      }
    }
  }

  // Calculate zeta-fluxes
  if (use_angular_fluxes) {
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {

        // Reconstruct radiation
        if (order == 1) {
          RadiationDonorCellA1(prim_rad, k, j);
        } else {
          RadiationPiecewiseLinearA1(prim_rad, k, j);
        }

        // Calculate radiation fluxes
        for (int l = zs; l <= ze+1; ++l) {
          for (int m = ps; m <= pe; ++m) {
            int lm = AngleInd(l, m, true, false);
            for (int i = is; i <= ie; ++i) {
              Real na1 = na1_(l,m,k,j,i);
              Real na1_n_0 = na1_n_0_(l,m,k,j,i);
              flux_a[ZETADIR](lm,k,j,i) =
                  na1_n_0 * (na1 > 0.0 ? ii_l_(lm,i) : ii_r_(lm,i));
            }
          }
        }
      }
    }
  }

  // Calculate psi-fluxes
  if (use_angular_fluxes) {
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {

        // Reconstruct radiation
        if (order == 1) {
          RadiationDonorCellA2(prim_rad, k, j);
        } else {
          RadiationPiecewiseLinearA2(prim_rad, k, j);
        }

        // Calculate radiation fluxes
        for (int l = zs; l <= ze; ++l) {
          for (int m = ps; m <= pe+1; ++m) {
            int lm = AngleInd(l, m, false, true);
            for (int i = is; i <= ie; ++i) {
              Real na2 = na2_(l,m,k,j,i);
              Real na2_n_0 = na2_n_0_(l,m,k,j,i);
              flux_a[PSIDIR](lm,k,j,i) =
                  na2_n_0 * (na2 > 0.0 ? ii_l_(lm,i) : ii_r_(lm,i));
            }
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for updating conserved quantities
// Inputs:
//   prim_in: primitive intensity
//   dt: timestep for stage of integration
// Outputs:
//   cons_out: conserved values updated

void Radiation::AddFluxDivergenceToAverage(AthenaArray<Real> &prim_in, const Real dt,
    AthenaArray<Real> &cons_out) {

  // Extract Coordinates and timestep
  Coordinates *pcoord = pmy_block->pcoord;

  // Calculate angle index range (including some but not all unnecessary ghost zones)
  int lms = AngleInd(zs, ps);
  int lme = AngleInd(ze, pe);

  // Go through all cells
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {

      // Determine poles
      bool left_pole = pcoord->IsPole(j);
      bool right_pole = pcoord->IsPole(j+1);

      // Calculate x1-divergence
      pcoord->Face1Area(k, j, is, ie+1, area_l_);
      for (int lm = lms; lm <= lme; ++lm) {
        for (int i = is; i <= ie; ++i) {
          flux_div_(lm,i) = area_l_(i+1) * flux_x[X1DIR](lm,k,j,i+1)
              - area_l_(i) * flux_x[X1DIR](lm,k,j,i);
        }
      }

      // Add x2-divergence
      if (js != je) {
        pcoord->Face2Area(k, j, is, ie, area_l_);
        pcoord->Face2Area(k, j+1, is, ie, area_r_);
        for (int lm = lms; lm <= lme; ++lm) {
          for (int i = is; i <= ie; ++i) {
            Real left_flux = left_pole ? 0.0 : -area_l_(i) * flux_x[X2DIR](lm,k,j,i);
            Real right_flux = right_pole ? 0.0 : area_r_(i) * flux_x[X2DIR](lm,k,j+1,i);
            flux_div_(lm,i) += left_flux + right_flux;
          }
        }
      }

      // Add x3-divergence
      if (ks != ke) {
        pcoord->Face3Area(k, j, is, ie, area_l_);
        pcoord->Face3Area(k+1, j, is, ie, area_r_);
        for (int lm = lms; lm <= lme; ++lm) {
          for (int i = is; i <= ie; ++i) {
            flux_div_(lm,i) += area_r_(i) * flux_x[X3DIR](lm,k+1,j,i)
                - area_l_(i) * flux_x[X3DIR](lm,k,j,i);
          }
        }
      }

      // Update conserved variables
      pcoord->CellVolume(k, j, is, ie, vol_);
      for (int lm = lms; lm <= lme; ++lm) {
        for (int i = is; i <= ie; ++i) {
          cons_out(lm,k,j,i) -= dt * flux_div_(lm,i) / vol_(i);
        }
      }
    }
  }

  // Go through all angles
  if (use_angular_fluxes) {
    for (int l = zs; l <= ze; ++l) {
      for (int m = ps; m <= pe; ++m) {

        // Determine poles
        bool left_pole = l == zs;
        bool right_pole = l == ze;

        // Calculate angle lengths and solid angles
        int lm = AngleInd(l, m);
        int lm_lc = AngleInd(l, m, true, false);
        int lm_rc = AngleInd(l + 1, m, true, false);
        int lm_cl = AngleInd(l, m, false, true);
        int lm_cr = AngleInd(l, m + 1, false, true);
        Real zeta_length_m = zeta_length(l,m);
        Real zeta_length_p = zeta_length(l,m+1);
        Real psi_length_m = psi_length(l,m);
        Real psi_length_p = psi_length(l+1,m);
        Real omega = solid_angle(l,m);

        // Go through all cells
        for (int k = ks; k <= ke; ++k) {
          for (int j = js; j <= je; ++j) {
            for (int i = is; i <= ie; ++i) {

              // Calculate zeta-divergence
              Real left_flux =
                  left_pole ? 0.0 : -psi_length_m * flux_a[ZETADIR](lm_lc,k,j,i);
              Real right_flux =
                  right_pole ? 0.0 : psi_length_p * flux_a[ZETADIR](lm_rc,k,j,i);
              Real flux_div = left_flux + right_flux;

              // Add psi-divergence
              flux_div += zeta_length_p * flux_a[PSIDIR](lm_cr,k,j,i)
                  - zeta_length_m * flux_a[PSIDIR](lm_cl,k,j,i);

              // Update conserved variables
              cons_out(lm,k,j,i) -= dt * flux_div / omega;
            }
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for averaging intensities according to integrator weights
// Inputs:
//   cons_out, cons_in_1, cons_in_2: conserved intensity arrays, possibly uninitialized
//   weights: integrator weights
// Outputs:
//   cons_out: weighted intensity
// Notes:
//   Same procedure as in MeshBlock::WeightedAve().

void Radiation::WeightedAve(AthenaArray<Real> &cons_out, AthenaArray<Real> &cons_in_1,
    AthenaArray<Real> &cons_in_2, const Real weights[3]) {

  // Apply averaging based on which weights are 0
  if (weights[2] != 0.0) {
    for (int n = 0; n < nang; ++n) {
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
          #pragma omp simd
          for (int i = is; i <= ie; ++i) {
            cons_out(n,k,j,i) = weights[0] * cons_out(n,k,j,i)
                + weights[1] * cons_in_1(n,k,j,i) + weights[2] * cons_in_2(n,k,j,i);
          }
        }
      }
    }
  } else if (weights[1] != 0.0) {
    for (int n = 0; n < nang; ++n) {
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
          #pragma omp simd
          for (int i = is; i <= ie; ++i) {
            cons_out(n,k,j,i) =
                weights[0] * cons_out(n,k,j,i) + weights[1] * cons_in_1(n,k,j,i);
          }
        }
      }
    }
  } else if (weights[0] != 0.0) {
    for (int n = 0; n < nang; ++n) {
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
          #pragma omp simd
          for (int i = is; i <= ie; ++i) {
            cons_out(n,k,j,i) = weights[0] * cons_out(n,k,j,i);
          }
        }
      }
    }
  } else {
    for (int n = 0; n < nang; ++n) {
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
          #pragma omp simd
          for (int i = is; i <= ie; ++i) {
            cons_out(n,k,j,i) = 0.0;
          }
        }
      }
    }
  }
  return;
}
