//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file radiation.cpp
//  \brief implementation of source-related functions in class Radiation

// C++ headers
#include <algorithm>  // max, min
#include <cmath>      // abs, isnan, pow, sqrt

// Athena++ headers
#include "radiation.hpp"
#include "../athena.hpp"                   // Real, indices, function prototypes
#include "../athena_arrays.hpp"            // AthenaArray
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../eos/eos.hpp"                  // EquationOfState
#include "../mesh/mesh.hpp"                // MeshBlock

// Declarations
bool FourthPolyRoot(const Real coef4, const Real tconst, Real *root);

//----------------------------------------------------------------------------------------
// Function for adding all source terms beyond those induced by coordinates
// Inputs:
//   time: time of simulation
//   dt: simulation timestep
//   prim_rad: primitive intensity at beginning of stage
//   prim_hydro: primitive hydro variables at beginning of step
//   prim_hydro_alt: primitive hydro variables at beginning of stage
//   cons_rad: conserved intensity after stage integration
//   cons_hydro: conserved hydro variables after stage integration
// Outputs:
//   cons: conserved intensity updated
//   cons_hydro: conserved hydro variables updated

void Radiation::AddSourceTerms(const Real time, const Real dt,
    const AthenaArray<Real> &prim_rad, const AthenaArray<Real> &prim_hydro,
    const AthenaArray<Real> &prim_hydro_alt, AthenaArray<Real> &cons_rad,
    AthenaArray<Real> &cons_hydro) {

  // Get adiabatic index
  Real gamma_adi = pmy_block->peos->GetGamma();

  // Go through outer loops of cells
  if (coupled_to_matter) {
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        pmy_block->pcoord->CellMetric(k, j, is, ie, g_, gi_);

        // Calculate zeroth and first moments of radiation before coupling
        if (affect_fluid) {
          for (int n = 0; n < 4; ++n) {
            for (int i = is; i <= ie; ++i) {
              moments_old_(n,i) = 0.0;
            }
          }
          for (int l = zs; l <= ze; ++l) {
            for (int m = ps; m <= pe; ++m) {
              int lm = AngleInd(l, m);
              for (int n = 0; n < 4; ++n) {
                for (int i = is; i <= ie; ++i) {
                  moments_old_(n,i) += n0_n_mu_(n,l,m,k,j,i) * cons_rad(lm,k,j,i)
                      / n0_n_mu_(0,l,m,k,j,i) * solid_angle(l,m);
                }
              }
            }
          }
        }

        // Calculate fluid velocity in tetrad frame
        for (int i = is; i <= ie; ++i) {
          Real uu1 = prim_hydro_alt(IVX,k,j,i);
          Real uu2 = prim_hydro_alt(IVY,k,j,i);
          Real uu3 = prim_hydro_alt(IVZ,k,j,i);
          Real temp_var = g_(I11,i) * SQR(uu1) + 2.0 * g_(I12,i) * uu1 * uu2
              + 2.0 * g_(I13,i) * uu1 * uu3 + g_(I22) * SQR(uu2)
              + 2.0 * g_(I23,i) * uu2 * uu3 + g_(I33,i) * SQR(uu3);
          Real uu0 = std::sqrt(1.0 + temp_var);
          u_tet_(0,i) = norm_to_tet_(0,0,k,j,i) * uu0 + norm_to_tet_(0,1,k,j,i) * uu1
              + norm_to_tet_(0,2,k,j,i) * uu2 + norm_to_tet_(0,3,k,j,i) * uu3;
          u_tet_(1,i) = norm_to_tet_(1,0,k,j,i) * uu0 + norm_to_tet_(1,1,k,j,i) * uu1
              + norm_to_tet_(1,2,k,j,i) * uu2 + norm_to_tet_(1,3,k,j,i) * uu3;
          u_tet_(2,i) = norm_to_tet_(2,0,k,j,i) * uu0 + norm_to_tet_(2,1,k,j,i) * uu1
              + norm_to_tet_(2,2,k,j,i) * uu2 + norm_to_tet_(2,3,k,j,i) * uu3;
          u_tet_(3,i) = norm_to_tet_(3,0,k,j,i) * uu0 + norm_to_tet_(3,1,k,j,i) * uu1
              + norm_to_tet_(3,2,k,j,i) * uu2 + norm_to_tet_(3,3,k,j,i) * uu3;
        }

        // Calculate radiation energy density in fluid frame
        for (int i = is; i <= ie; ++i) {
          ee_f_minus_(i) = 0.0;
          for (int l = zs; l <= ze; ++l) {
            for (int m = ps; m <= pe; ++m) {
              int lm = AngleInd(l, m);
              Real u_n = -u_tet_(0,i) * nh_cc_(0,l,m) + u_tet_(1,i) * nh_cc_(1,l,m)
                  + u_tet_(2,i) * nh_cc_(2,l,m) + u_tet_(3,i) * nh_cc_(3,l,m);
              ee_f_minus_(i) += prim_rad(lm,k,j,i) * SQR(u_n) * solid_angle(l,m);
            }
          }
        }

        // Calculate velocity modification effective timestep
        // TODO(CJW): below only works for VL2
        Real dt_velocity = dt;
        if (prim_hydro.data() != prim_hydro_alt.data()) {
          dt_velocity /= 2.0;
        }

        // Modify velocity in radiation-dominated regime
        for (int i = is; i <= ie; ++i) {

          // Check that radiation dominates fluid
          Real rho = prim_hydro_alt(IDN,k,j,i);
          Real pgas = prim_hydro_alt(IPR,k,j,i);
          Real egas = rho + pgas / (gamma_adi - 1.);
          if (ee_f_minus_(i) <= egas) {
            continue;
          }

          // Calculate radiation moments
          Real rr_tet00 = 0.0;
          Real rr_tet01 = 0.0;
          Real rr_tet02 = 0.0;
          Real rr_tet03 = 0.0;
          Real rr_tet11 = 0.0;
          Real rr_tet12 = 0.0;
          Real rr_tet13 = 0.0;
          Real rr_tet22 = 0.0;
          Real rr_tet23 = 0.0;
          Real rr_tet33 = 0.0;
          for (int l = zs; l <= ze; ++l) {
            for (int m = ps; m <= pe; ++m) {
              int lm = AngleInd(l, m);
              rr_tet00 += prim_rad(lm,k,j,i) * solid_angle(l,m);
              rr_tet01 += prim_rad(lm,k,j,i) * nh_cc_(1,l,m) * solid_angle(l,m);
              rr_tet02 += prim_rad(lm,k,j,i) * nh_cc_(2,l,m) * solid_angle(l,m);
              rr_tet03 += prim_rad(lm,k,j,i) * nh_cc_(3,l,m) * solid_angle(l,m);
              rr_tet11 += prim_rad(lm,k,j,i) * SQR(nh_cc_(1,l,m)) * solid_angle(l,m);
              rr_tet12 +=
                  prim_rad(lm,k,j,i) * nh_cc_(1,l,m) * nh_cc_(2,l,m) * solid_angle(l,m);
              rr_tet13 +=
                  prim_rad(lm,k,j,i) * nh_cc_(1,l,m) * nh_cc_(3,l,m) * solid_angle(l,m);
              rr_tet22 += prim_rad(lm,k,j,i) * SQR(nh_cc_(2,l,m)) * solid_angle(l,m);
              rr_tet23 +=
                  prim_rad(lm,k,j,i) * nh_cc_(2,l,m) * nh_cc_(3,l,m) * solid_angle(l,m);
              rr_tet33 += prim_rad(lm,k,j,i) * SQR(nh_cc_(3,l,m)) * solid_angle(l,m);
            }
          }

          // Calculate radiation-frame velocity
          Real vrad_tet1 = rr_tet01 / rr_tet00;
          Real vrad_tet2 = rr_tet02 / rr_tet00;
          Real vrad_tet3 = rr_tet03 / rr_tet00;
          Real vrad_sq = SQR(vrad_tet1) + SQR(vrad_tet2) + SQR(vrad_tet3);
          if (vrad_sq > v_sq_max) {
            Real ratio = std::sqrt(v_sq_max / vrad_sq);
            vrad_tet1 *= ratio;
            vrad_tet2 *= ratio;
            vrad_tet3 *= ratio;
            vrad_sq = v_sq_max;
          }
          Real urad_tet0 = 1.0 / std::sqrt(1.0 - vrad_sq);
          Real urad_tet1 = urad_tet0 * vrad_tet1;
          Real urad_tet2 = urad_tet0 * vrad_tet2;
          Real urad_tet3 = urad_tet0 * vrad_tet3;

          // Calculate current fluid momentum
          Real wgas = egas + pgas;
          Real mgas_tet1 = wgas * u_tet_(0,i) * u_tet_(1,i);
          Real mgas_tet2 = wgas * u_tet_(0,i) * u_tet_(2,i);
          Real mgas_tet3 = wgas * u_tet_(0,i) * u_tet_(3,i);

          // Calculate fluid momentum if accelerated to radiation frame
          Real mgas_rad_tet1 = wgas * urad_tet0 * urad_tet1;
          Real mgas_rad_tet2 = wgas * urad_tet0 * urad_tet2;
          Real mgas_rad_tet3 = wgas * urad_tet0 * urad_tet3;

          // Estimate change in fluid momentum from source terms
          Real tt = pgas / rho;
          Real k_a = rho * opacity(OPAA,k,j,i);
          Real k_s = rho * opacity(OPAS,k,j,i);
          Real gg_tet1 = -(k_a * arad * SQR(SQR(tt)) + k_s * ee_f_minus_(i))
              * u_tet_(1,i) - (k_a + k_s) * (-u_tet_(0,i) * rr_tet01
              + u_tet_(1,i) * rr_tet11 + u_tet_(2,i) * rr_tet12
              + u_tet_(3,i) * rr_tet13);
          Real gg_tet2 = -(k_a * arad * SQR(SQR(tt)) + k_s * ee_f_minus_(i))
              * u_tet_(2,i) - (k_a + k_s) * (-u_tet_(0,i) * rr_tet02
              + u_tet_(1,i) * rr_tet12 + u_tet_(2,i) * rr_tet22
              + u_tet_(3,i) * rr_tet23);
          Real gg_tet3 = -(k_a * arad * SQR(SQR(tt)) + k_s * ee_f_minus_(i))
              * u_tet_(3,i) - (k_a + k_s) * (-u_tet_(0,i) * rr_tet03
              + u_tet_(1,i) * rr_tet13 + u_tet_(2,i) * rr_tet23
              + u_tet_(3,i) * rr_tet33);
          Real dmgas_tet1 = dt_velocity * gg_tet1;
          Real dmgas_tet2 = dt_velocity * gg_tet2;
          Real dmgas_tet3 = dt_velocity * gg_tet3;

          // Estimate new fluid velocity
          Real frac1 = mgas_rad_tet1 == mgas_tet1 ? 0.0
              : dmgas_tet1 / (mgas_rad_tet1 - mgas_tet1);
          Real frac2 = mgas_rad_tet2 == mgas_tet2 ? 0.0
              : dmgas_tet2 / (mgas_rad_tet2 - mgas_tet2);
          Real frac3 = mgas_rad_tet3 == mgas_tet3 ? 0.0
              : dmgas_tet3 / (mgas_rad_tet3 - mgas_tet3);
          frac1 = std::min(std::max(frac1, 0.0), 1.0);
          frac2 = std::min(std::max(frac2, 0.0), 1.0);
          frac3 = std::min(std::max(frac3, 0.0), 1.0);
          u_tet_(1,i) = (1.0 - frac1) * u_tet_(1,i) + frac1 * urad_tet1;
          u_tet_(2,i) = (1.0 - frac2) * u_tet_(2,i) + frac2 * urad_tet2;
          u_tet_(3,i) = (1.0 - frac3) * u_tet_(3,i) + frac3 * urad_tet3;
          u_tet_(0,i) =
              std::sqrt(1.0 + SQR(u_tet_(1,i)) + SQR(u_tet_(2,i)) + SQR(u_tet_(3,i)));
        }

        // Calculate quartic coefficients
        for (int i = is; i <= ie; ++i) {
          Real rho = prim_hydro(IDN,k,j,i);
          Real tt_minus = prim_hydro(IPR,k,j,i) / rho;
          Real k_a = rho * opacity(OPAA,k,j,i);
          Real k_s = rho * opacity(OPAS,k,j,i);
          Real k_tot = k_a + k_s;
          Real var_a = 0.0;
          Real var_b = 0.0;
          ee_f_minus_(i) = 0.0;
          for (int l = zs; l <= ze; ++l) {
            for (int m = ps; m <= pe; ++m) {
              int lm = AngleInd(l, m);
              Real ii_minus = prim_rad(lm,k,j,i);
              Real u_n = -u_tet_(0,i) * nh_cc_(0,l,m) + u_tet_(1,i) * nh_cc_(1,l,m)
                  + u_tet_(2,i) * nh_cc_(2,l,m) + u_tet_(3,i) * nh_cc_(3,l,m);
              Real denominator = 1.0 - dt * k_tot * u_n;
              var_a += ii_minus * SQR(u_n) / denominator * solid_angle(l,m);
              var_b += 1.0 / u_n / denominator * solid_angle(l,m);
              ee_f_minus_(i) += ii_minus * SQR(u_n) * solid_angle(l,m);
            }
          }
          var_b *= dt / (4.0 * PI);
          coefficients_(0,i) =
              -(gamma_adi - 1.0) / rho * var_b * k_a * arad / (1.0 + var_b * k_s);
          coefficients_(1,i) = -tt_minus - (gamma_adi - 1.0) / rho * ee_f_minus_(i)
              + (gamma_adi - 1.0) / rho * var_a / (1.0 + var_b * k_s);
        }

        // Calculate new gas temperature
        for (int i = is; i <= ie; ++i) {
          bad_cell_(i) = false;
          if (std::abs(coefficients_(0,i)) > TINY_NUMBER) {
            bool quartic_flag =
                FourthPolyRoot(coefficients_(0,i), coefficients_(1,i), &tt_plus_(i));
            if (not quartic_flag or std::isnan(tt_plus_(i))) {
              bad_cell_(i) = true;
              tt_plus_(i) = prim_hydro(IPR,k,j,i) / prim_hydro(IDN,k,j,i);
            }
          } else {
            tt_plus_(i) = -coefficients_(1,i);
          }
        }

        // Calculate new radiation energy density
        for (int i = is; i <= ie; ++i) {
          if (not bad_cell_(i)) {
            Real rho = prim_hydro(IDN,k,j,i);
            Real tt_minus = prim_hydro(IPR,k,j,i) / rho;
            ee_f_plus_(i) =
                ee_f_minus_(i) + rho / (gamma_adi - 1.0) * (tt_minus - tt_plus_(i));
            ee_f_plus_(i) = std::max(ee_f_plus_(i), 0.0);
          }
        }

        // Calculate new intensity
        for (int i = is; i <= ie; ++i) {
          if (not bad_cell_(i)) {
            Real rho = prim_hydro(IDN,k,j,i);
            Real k_a = rho * opacity(OPAA,k,j,i);
            Real k_s = rho * opacity(OPAS,k,j,i);
            Real k_tot = k_a + k_s;
            for (int l = zs; l <= ze; ++l) {
              for (int m = ps; m <= pe; ++m) {
                int lm = AngleInd(l, m);
                Real ii_minus = prim_rad(lm,k,j,i);
                Real u_n = -u_tet_(0,i) * nh_cc_(0,l,m) + u_tet_(1,i) * nh_cc_(1,l,m)
                    + u_tet_(2,i) * nh_cc_(2,l,m) + u_tet_(3,i) * nh_cc_(3,l,m);
                Real ii_plus = (ii_minus - dt / (4.0 * PI) / u_n / SQR(u_n)
                    * (k_a * arad * SQR(SQR(tt_plus_(i))) + k_s * ee_f_plus_(i)))
                    / (1.0 - dt * k_tot * u_n);
                cons_rad(lm,k,j,i) += (ii_plus - ii_minus) * n0_n_mu_(0,l,m,k,j,i);
                cons_rad(lm,k,j,i) = std::min(cons_rad(lm,k,j,i), 0.0);
              }
            }
          }
        }

        // Calculate zeroth and first moments of radiation after coupling
        if (affect_fluid) {
          for (int n = 0; n < 4; ++n) {
            for (int i = is; i <= ie; ++i) {
              moments_new_(n,i) = 0.0;
            }
          }
          for (int l = zs; l <= ze; ++l) {
            for (int m = ps; m <= pe; ++m) {
              int lm = AngleInd(l, m);
              for (int n = 0; n < 4; ++n) {
                for (int i = is; i <= ie; ++i) {
                  moments_new_(n,i) += n0_n_mu_(n,l,m,k,j,i) * cons_rad(lm,k,j,i)
                      / n0_n_mu_(0,l,m,k,j,i) * solid_angle(l,m);
                }
              }
            }
          }
        }

        // Apply radiation-fluid coupling to fluid
        if (affect_fluid) {
          for (int i = is; i <= ie; ++i) {
            cons_hydro(IEN,k,j,i) += moments_old_(0,i) - moments_new_(0,i);
            cons_hydro(IM1,k,j,i) += moments_old_(1,i) - moments_new_(1,i);
            cons_hydro(IM2,k,j,i) += moments_old_(2,i) - moments_new_(2,i);
            cons_hydro(IM3,k,j,i) += moments_old_(3,i) - moments_new_(3,i);
          }
        }
      }
    }
  }

  // Apply user source terms
  if (UserSourceTerm != nullptr) {
    UserSourceTerm(pmy_block, time, dt, prim_rad, cons_rad);
  }
  return;
}

//----------------------------------------------------------------------------------------
// Opacity enrollment
// Inputs:
//   MyOpacityFunction: user-defined function from problem generator
// Outputs: (none)
// Notes:
//   If nothing else enrolled, default function keeps opacities (not absorption
//     coefficients) at their initial values.

void Radiation::EnrollOpacityFunction(OpacityFunc MyOpacityFunction)
{
  UpdateOpacity = MyOpacityFunction;
  return;
}

//----------------------------------------------------------------------------------------
// Exact solution for fourth order polynomial
// Inputs:
//   coef4: quartic coefficient
//   tconst: constant coefficient
// Outputs:
//   root: solution to equation
//   returned value: flag indicating success
// Notes:
//   Polynomial has the form coef4 * x^4 + x + tconst = 0.

bool FourthPolyRoot(const Real coef4, const Real tconst, Real *root) {

  // Calculate real root of z^3 - 4*tconst/coef4 * z - 1/coef4^2 = 0
  Real asquar = coef4 * coef4;
  Real acubic = coef4 * asquar;
  Real ccubic = tconst * tconst * tconst;
  Real delta1 = 0.25 - 64.0 * ccubic * coef4 / 27.0;
  if (delta1 < 0.0) {
    return false;
  }
  delta1 = std::sqrt(delta1);
  if (delta1 < 0.5) {
    return false;
  }
  Real zroot;
  if (delta1 > 1.0e11) {  // to avoid small number cancellation
    zroot = std::pow(delta1, -2.0/3.0) / 3.0;
  } else {
    zroot = std::pow(0.5 + delta1, 1.0/3.0) - std::pow(-0.5 + delta1, 1.0/3.0);
  }
  if (zroot < 0.0) {
    return false;
  }
  zroot *= std::pow(coef4, -2.0/3.0);

  // Calculate quartic root using cubic root
  Real rcoef = std::sqrt(zroot);
  Real delta2 = -zroot + 2.0 / (coef4 * rcoef);
  if (delta2 < 0.0) {
    return false;
  }
  delta2 = std::sqrt(delta2);
  *root = 0.5 * (delta2 - rcoef);
  if (*root < 0.0) {
    return false;
  }
  return true;
}
