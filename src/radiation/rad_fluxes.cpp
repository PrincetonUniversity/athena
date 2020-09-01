//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file radiation.cpp
//  \brief implementation of flux-related functions in class Radiation

// Athena++ headers
#include "radiation.hpp"
#include "../athena_arrays.hpp"            // AthenaArray
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../mesh/mesh.hpp"                // MeshBlock

//----------------------------------------------------------------------------------------
// Function for calculating spatial and angular fluxes
// Inputs:
//   prim_rad: primitive intensity
//   prim_hydro: primitive hydrodynamic variables
//   order: reconstruction order
// Outputs:
//   this->flux_x, this->flux_a: fluxes set

void Radiation::CalculateFluxes(AthenaArray<Real> &prim_rad,
    const AthenaArray<Real> &prim_hydro, int order) {

  // Calculate x1-fluxes
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {

      // Reconstruct radiation
      if (order == 1) {
        RadiationDonorCellX1(prim_rad, k, j);
      } else {
        RadiationPiecewiseLinearX1(prim_rad, k, j);
      }

      // Calculate radiation fluxes
      for (int l = zs; l <= ze; ++l) {
        for (int m = ps; m <= pe; ++m) {
          int lm = AngleInd(l, m);
          for (int i = is; i <= ie+1; ++i) {
            Real n1_n_0 = n1_n_0_(l,m,k,j,i);
            if (n1_n_0 < 0.0) {
              flux_x[X1DIR](lm,k,j,i) = n1_n_0 * rad_l_(lm,i);
            } else {
              flux_x[X1DIR](lm,k,j,i) = n1_n_0 * rad_r_(lm,i);
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

        // Calculate radiation fluxes
        for (int l = zs; l <= ze; ++l) {
          for (int m = ps; m <= pe; ++m) {
            int lm = AngleInd(l, m);
            for (int i = is; i <= ie; ++i) {
              Real n2_n_0 = n2_n_0_(l,m,k,j,i);
              if (n2_n_0 < 0.0) {
                flux_x[X2DIR](lm,k,j,i) = n2_n_0 * rad_l_(lm,i);
              } else {
                flux_x[X2DIR](lm,k,j,i) = n2_n_0 * rad_r_(lm,i);
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

        // Calculate radiation fluxes
        for (int l = zs; l <= ze; ++l) {
          for (int m = ps; m <= pe; ++m) {
            int lm = AngleInd(l, m);
            for (int i = is; i <= ie; ++i) {
              Real n3_n_0 = n3_n_0_(l,m,k,j,i);
              if (n3_n_0 < 0.0) {
                flux_x[X3DIR](lm,k,j,i) = n3_n_0 * rad_l_(lm,i);
              } else {
                flux_x[X3DIR](lm,k,j,i) = n3_n_0 * rad_r_(lm,i);
              }
            }
          }
        }
      }
    }
  }

  // Calculate zeta-fluxes
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
          int lm_c = AngleInd(l, m, true, false);
          for (int i = is; i <= ie; ++i) {
            Real na1_n_0 = na1_n_0_(l,m,k,j,i);
            if (na1_n_0 < 0.0) {
              flux_a[ZETADIR](lm_c,k,j,i) = na1_n_0 * rad_l_(lm_c,i);
            } else {
              flux_a[ZETADIR](lm_c,k,j,i) = na1_n_0 * rad_r_(lm_c,i);
            }
          }
        }
      }
    }
  }

  // Calculate psi-fluxes
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
          int lm_c = AngleInd(l, m, false, true);
          for (int i = is; i <= ie; ++i) {
            Real na2_n_0 = na2_n_0_(l,m,k,j,i);
            if (na2_n_0 < 0.0) {
              flux_a[PSIDIR](lm_c,k,j,i) = na2_n_0 * rad_l_(lm_c,i);
            } else {
              flux_a[PSIDIR](lm_c,k,j,i) = na2_n_0 * rad_r_(lm_c,i);
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
  return;
}
