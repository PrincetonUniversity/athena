//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file radiation.cpp
//  \brief implementation of flux-related functions in class Radiation

// Athena++ headers
#include "radiation.hpp"
#include "../athena.hpp"                   // Real, indices
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

      // Calculate pure radiation fluxes
      if (not coupled_to_matter) {
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

      // Calculate radiation-hydrodynamic fluxes
      if (coupled_to_matter) {
        pmy_block->pcoord->CellMetric(k, j, is - 1, ie + 1, g_, gi_);
        pmy_block->pcoord->CenterWidth1(k, j, is - 1, ie + 1, widths_l_);
        for (int i = is-1; i <= ie+1; ++i) {
          ee_l_(i) = 0.0;
        }
        for (int l = zs; l <= ze; ++l) {
          for (int m = ps; m <= pe; ++m) {
            int lm = AngleInd(l, m);
            for (int i = is-1; i <= ie+1; ++i) {
              ee_l_(i) += n0_n_mu_(0,l,m,k,j,i) * prim_rad(lm,k,j,i) * solid_angle(l,m);
            }
          }
        }
        for (int i = is-1; i <= ie+1; ++i) {
          Real uu1 = prim_hydro(IVX,k,j,i);
          Real uu2 = prim_hydro(IVY,k,j,i);
          Real uu3 = prim_hydro(IVZ,k,j,i);
          Real gamma_sq = 1.0 + g_(I11,i) * SQR(uu1) + 2.0 * g_(I12,i) * uu1 * uu2
              + 2.0 * g_(I13,i) * uu1 * uu3 + g_(I22,i) * SQR(uu2)
              + 2.0 * g_(I23,i) * uu2 * uu3 + g_(I33,i) * SQR(uu3);
          Real gamma = std::sqrt(gamma_sq);
          Real alpha = std::sqrt(-1.0 / gi_(I00,i));
          Real u0 = gamma / alpha;
          Real u1 = uu1 - alpha * gamma * gi_(I01,i);
          v_l_(i) = u1 / u0;
        }
        for (int i = is; i <= ie+1; ++i) {
          Real dx_l = pmy_block->pcoord->x1f(i) - pmy_block->pcoord->x1v(i-1);
          Real dx_r = pmy_block->pcoord->x1v(i) - pmy_block->pcoord->x1f(i);
          Real ee = 0.5 * (ee_l_(i-1) + ee_l_(i));
          Real grad_ee =
              2.0 * std::abs(ee_l_(i-1) - ee_l_(i)) / (widths_l_(i-1) + widths_l_(i));
          Real kappa_l = opacity(OPAA,k,j,i-1) + opacity(OPAS,k,j,i-1);
          Real kappa_r = opacity(OPAA,k,j,i) + opacity(OPAS,k,j,i);
          Real tau_factor_l = kappa_l > 0.0 ? 1.0 : 0.0;
          Real tau_factor_r = kappa_r > 0.0 ? 1.0 : 0.0;
          if (grad_ee > 0.0) {
            Real tau_l = ee / grad_ee * prim_hydro(IDN,k,j,i-1) * kappa_l;
            Real tau_r = ee / grad_ee * prim_hydro(IDN,k,j,i) * kappa_r;
            tau_factor_l = -std::expm1(-SQR(tau_l));
            tau_factor_r = -std::expm1(-SQR(tau_r));
          }
          for (int l = zs; l <= ze; ++l) {
            for (int m = ps; m <= pe; ++m) {
              int lm = AngleInd(l, m);
              Real ii_diffusive_l = rad_l_(lm,i) * (1.0 - tau_factor_l);
              Real ii_diffusive_r = rad_r_(lm,i) * (1.0 - tau_factor_r);
              Real v_diffusive = n1_n_0_(l,m,k,j,i);
              if (v_diffusive < 0.0) {
                flux_x[X1DIR](lm,k,j,i) = v_diffusive * ii_diffusive_l;
              } else {
                flux_x[X1DIR](lm,k,j,i) = v_diffusive * ii_diffusive_r;
              }
              Real ii_advective_l = rad_l_(lm,i) * tau_factor_l;
              Real ii_advective_r = rad_r_(lm,i) * tau_factor_r;
              Real v_advective = (dx_r * v_l_(i-1) * n0_n_mu_(0,l,m,k,j,i-1)
                  / nmu_(0,l,m,k,j,i-1) + dx_l * v_l_(i) * n0_n_mu_(0,l,m,k,j,i)
                  / nmu_(0,l,m,k,j,i)) / (dx_l + dx_r);
              if (v_advective < 0.0) {
                flux_x[X1DIR](lm,k,j,i) += v_advective * ii_advective_l;
              } else {
                flux_x[X1DIR](lm,k,j,i) += v_advective * ii_advective_r;
              }
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
