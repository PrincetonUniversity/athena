//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file adiabatic_mhd_sr.cpp
//! \brief Implements functions for going between primitive and conserved variables in
//! special-relativistic MHD, as well as for computing wavespeeds.

// C headers

// C++ headers
#include <algorithm>  // max, min
#include <cmath>      // abs, acos, cbrt, cos, isfinite, sqrt

// Athena++ headers
#include "../athena.hpp"                   // enums, macros
#include "../athena_arrays.hpp"            // AthenaArray
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../field/field.hpp"              // FaceField
#include "../mesh/mesh.hpp"                // MeshBlock
#include "../parameter_input.hpp"          // ParameterInput
#include "eos.hpp"

namespace {
// Declarations
bool ConservedToPrimitiveNormal(
    const AthenaArray<Real> &dd_vals, const AthenaArray<Real> &ee_vals,
    const AthenaArray<Real> &mm_vals, const AthenaArray<Real> &bb_vals,
    const AthenaArray<Real> &tt_vals, Real gamma_adi, Real pgas_old, int k, int j, int i,
    AthenaArray<Real> &prim, Real *p_gamma_lor, Real *p_pmag);
} // namespace

//----------------------------------------------------------------------------------------
// Constructor
// Inputs:
//   pmb: pointer to MeshBlock
//   pin: pointer to runtime inputs

EquationOfState::EquationOfState(MeshBlock *pmb, ParameterInput *pin) :
    pmy_block_(pmb),
    gamma_{pin->GetReal("hydro", "gamma")},
    density_floor_{pin->GetOrAddReal("hydro", "dfloor", std::sqrt(1024*float_min))},
    pressure_floor_{pin->GetOrAddReal("hydro", "pfloor", std::sqrt(1024*float_min))},
    scalar_floor_{pin->GetOrAddReal("hydro", "sfloor", std::sqrt(1024*float_min))},
    sigma_max_{pin->GetOrAddReal("hydro", "sigma_max",  0.0)},
    beta_min_{pin->GetOrAddReal("hydro", "beta_min", 0.0)},
    gamma_max_{pin->GetOrAddReal("hydro", "gamma_max", 1000.0)} {
  int nc1 = pmb->ncells1;
  normal_dd_.NewAthenaArray(nc1);
  normal_ee_.NewAthenaArray(nc1);
  normal_mm_.NewAthenaArray(4,nc1);
  normal_bb_.NewAthenaArray(4,nc1);
  normal_tt_.NewAthenaArray(nc1);
}

//----------------------------------------------------------------------------------------
// Variable inverter
// Inputs:
//   cons: conserved quantities
//   prim_old: primitive quantities from previous half timestep
//   bb: face-centered magnetic field
//   pco: pointer to Coordinates
//   il, iu, jl, ju, kl, ku: index bounds of region to be updated
// Outputs:
//   prim: primitives
//   bb_cc: cell-centered magnetic field
// Notes:
//   Simpler version without magnetic fields found in adiabatic_hydro_sr.cpp.
//   More complex version for GR found in adiabatic_mhd_gr.cpp.

void EquationOfState::ConservedToPrimitive(
    AthenaArray<Real> &cons, const AthenaArray<Real> &prim_old, const FaceField &bb,
    AthenaArray<Real> &prim, AthenaArray<Real> &bb_cc, Coordinates *pco, int il, int iu,
    int jl, int ju, int kl, int ku) {
  // Parameters
  const Real mm_sq_ee_sq_max = 1.0 - 1.0e-12;  // max. of squared momentum over energy

  // Extract ratio of specific heats
  const Real &gamma_adi = gamma_;

  // Interpolate magnetic field from faces to cell centers
  pmy_block_->pfield->CalculateCellCenteredField(bb, bb_cc, pco, il, iu, jl, ju, kl, ku);

  // Go through all rows
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      // Extract conserved variables
      for (int i=il; i<=iu; ++i) {
        normal_dd_(i) = cons(IDN,k,j,i);
        normal_ee_(i) = cons(IEN,k,j,i);
        normal_mm_(0,i) =
            SQR(cons(IM1,k,j,i)) + SQR(cons(IM2,k,j,i)) + SQR(cons(IM3,k,j,i));
        normal_mm_(1,i) = cons(IM1,k,j,i);
        normal_mm_(2,i) = cons(IM2,k,j,i);
        normal_mm_(3,i) = cons(IM3,k,j,i);
        normal_bb_(0,i) =
            SQR(bb_cc(IB1,k,j,i)) + SQR(bb_cc(IB2,k,j,i)) + SQR(bb_cc(IB3,k,j,i));
        normal_bb_(1,i) = bb_cc(IB1,k,j,i);
        normal_bb_(2,i) = bb_cc(IB2,k,j,i);
        normal_bb_(3,i) = bb_cc(IB3,k,j,i);
        normal_tt_(i) = cons(IM1,k,j,i) * bb_cc(IB1,k,j,i)
            + cons(IM2,k,j,i) * bb_cc(IB2,k,j,i) + cons(IM3,k,j,i) * bb_cc(IB3,k,j,i);
      }

      // Go through cells
      for (int i=il; i<=iu; ++i) {
        // Set flag indicating conserved values need adjusting at end
        bool fixed = false;

        // Calculate floors for density and pressure
        Real density_floor_local = density_floor_;
        Real pressure_floor_local = pressure_floor_;

        // Ensure conserved density is large enough
        Real dd_min = density_floor_local;
        if (normal_dd_(i) < dd_min) {
          normal_dd_(i) = dd_min;
          fixed = true;
        }

        // Ensure conserved energy is large enough
        Real ee_min = density_floor_local + pressure_floor_local/(gamma_adi-1.0)
                      + 0.5*normal_bb_(0,i);
        if (normal_ee_(i) < ee_min) {
          normal_ee_(i) = ee_min;
          fixed = true;
        }

        // Ensure conserved momentum is not too large given energy
        Real mm_sq_max = mm_sq_ee_sq_max * SQR(normal_ee_(i));
        if (normal_mm_(0,i) > mm_sq_max) {
          Real factor = std::sqrt(mm_sq_max/normal_mm_(0,i));
          normal_mm_(0,i) = mm_sq_max;
          normal_mm_(1,i) *= factor;
          normal_mm_(2,i) *= factor;
          normal_mm_(3,i) *= factor;
          normal_tt_(i) *= factor;
          fixed = true;
        }

        // Define fallback state
        Real rho_old = std::max(prim(IDN,k,j,i), density_floor_local);
        Real pgas_old = std::max(prim(IPR,k,j,i), pressure_floor_local);
        Real u1_old = prim(IVX,k,j,i);
        Real u2_old = prim(IVY,k,j,i);
        Real u3_old = prim(IVZ,k,j,i);

        // Set primitives
        Real gamma, pmag;
        bool success = ConservedToPrimitiveNormal(normal_dd_, normal_ee_, normal_mm_,
                                                  normal_bb_, normal_tt_, gamma_adi,
                                                  pgas_old, k, j, i, prim,
                                                  &gamma, &pmag);

        // Handle failures
        if (!success) {
          prim(IDN,k,j,i) = rho_old;
          prim(IPR,k,j,i) = pgas_old;
          prim(IVX,k,j,i) = u1_old;
          prim(IVY,k,j,i) = u2_old;
          prim(IVZ,k,j,i) = u3_old;
          fixed = true;
        }

        // Apply density and gas pressure floors in normal frame
        if (sigma_max_ > 0.0) {
          density_floor_local = std::max(density_floor_local, 2.0*pmag/sigma_max_);
        }
        if (beta_min_ > 0.0) {
          pressure_floor_local = std::max(pressure_floor_local, beta_min_*pmag);
        }
        Real rho_add = std::max(density_floor_local-prim(IDN,k,j,i),
                                                static_cast<Real>(0.0));
        Real pgas_add = std::max(pressure_floor_local-prim(IPR,k,j,i),
                                                static_cast<Real>(0.0));
        if (success && (rho_add > 0.0 || pgas_add > 0.0)) {
          // Adjust conserved density and energy
          Real wgas_add = rho_add + gamma_adi/(gamma_adi-1.0) * pgas_add;
          normal_dd_(i) += rho_add * gamma;
          normal_ee_(i) += wgas_add * SQR(gamma) + pgas_add;

          // Recalculate primitives
          success = ConservedToPrimitiveNormal(normal_dd_, normal_ee_, normal_mm_,
                                               normal_bb_, normal_tt_, gamma_adi,
                                               pgas_old, k, j, i, prim, &gamma,
                                               &pmag);

          // Handle failures
          if (!success) {
            prim(IDN,k,j,i) = rho_old;
            prim(IPR,k,j,i) = pgas_old;
            prim(IVX,k,j,i) = u1_old;
            prim(IVY,k,j,i) = u2_old;
            prim(IVZ,k,j,i) = u3_old;
          }
          fixed = true;
        }

        // Apply velocity ceiling
        Real &u1 = prim(IVX,k,j,i);
        Real &u2 = prim(IVY,k,j,i);
        Real &u3 = prim(IVZ,k,j,i);
        if (!success) {
          gamma = std::sqrt(1.0 + SQR(u1) + SQR(u2) + SQR(u3));
        }
        bool velocity_ceiling = false;
        if (gamma > gamma_max_) {
          Real factor = std::sqrt((SQR(gamma_max_)-1.0) / (SQR(gamma)-1.0));
          u1 *= factor;
          u2 *= factor;
          u3 *= factor;
          fixed = true;
          velocity_ceiling = true;
        }

        // Recalculate density and pressure floors given new velocity
        if (velocity_ceiling) {
          Real b0 = u1 * bb_cc(IB1,k,j,i) + u2 * bb_cc(IB2,k,j,i) + u3 * bb_cc(IB3,k,j,i);
          pmag = 0.5 * (normal_bb_(0,i)/SQR(gamma) + SQR(b0/gamma));
        }
        density_floor_local = density_floor_;
        if (sigma_max_ > 0.0) {
          density_floor_local = std::max(density_floor_local, 2.0*pmag/sigma_max_);
        }
        pressure_floor_local = pressure_floor_;
        if (beta_min_ > 0.0) {
          pressure_floor_local = std::max(pressure_floor_local, beta_min_*pmag);
        }

        // Apply density and gas pressure floors in fluid frame
        Real &rho = prim(IDN,k,j,i);
        Real &pgas = prim(IPR,k,j,i);
        if (rho < density_floor_local) {
          rho = density_floor_local;
          fixed = true;
        }
        if (pgas < pressure_floor_local) {
          pgas = pressure_floor_local;
          fixed = true;
        }
        if (!success) {
          rho = density_floor_local;
          pgas = pressure_floor_local;
          u1 = u2 = u3 = 0.0;
        }

        // Ensure conserved variables match primitives
        if (fixed) {
          gamma = std::sqrt(1.0 + SQR(u1) + SQR(u2) + SQR(u3));
          Real bb1 = bb_cc(IB1,k,j,i);
          Real bb2 = bb_cc(IB2,k,j,i);
          Real bb3 = bb_cc(IB3,k,j,i);
          Real b0 = bb1 * u1 + bb2 * u2 + bb3 * u3;
          Real b1 = (bb1 + b0 * u1) / gamma;
          Real b2 = (bb2 + b0 * u2) / gamma;
          Real b3 = (bb3 + b0 * u3) / gamma;
          Real b_sq = -SQR(b0) + SQR(b1) + SQR(b2) + SQR(b3);
          Real wgas = rho + gamma_adi / (gamma_adi - 1.0) * pgas;
          Real wtot = wgas + b_sq;
          pmag = 0.5 * b_sq;
          cons(IDN,k,j,i) = rho * gamma;
          cons(IEN,k,j,i) = wgas * SQR(gamma) - SQR(b0) - (pgas + pmag);
          cons(IM1,k,j,i) = wgas * gamma * u1 - b0 * b1;
          cons(IM2,k,j,i) = wgas * gamma * u2 - b0 * b2;
          cons(IM3,k,j,i) = wgas * gamma * u3 - b0 * b3;
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for converting all primitives to conserved variables
// Inputs:
//   prim: 3D array of primitives
//   bc: 3D array of cell-centered magnetic fields
//   pco: pointer to Coordinates
//   il,iu,jl,ju,kl,ku: index bounds of region to be updated
// Outputs:
//   cons: 3D array of conserved variables

void EquationOfState::PrimitiveToConserved(
    const AthenaArray<Real> &prim, const AthenaArray<Real> &bc,
    AthenaArray<Real> &cons, Coordinates *pco,
    int il, int iu, int jl, int ju, int kl, int ku) {
  // Calculate reduced ratio of specific heats
  Real gamma_adi_red = gamma_/(gamma_-1.0);

  // Go through all cells
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        // Extract primitives and magnetic fields
        Real rho = prim(IDN,k,j,i);
        Real pgas = prim(IPR,k,j,i);
        Real u1 = prim(IVX,k,j,i);
        Real u2 = prim(IVY,k,j,i);
        Real u3 = prim(IVZ,k,j,i);
        Real bb1 = bc(IB1,k,j,i);
        Real bb2 = bc(IB2,k,j,i);
        Real bb3 = bc(IB3,k,j,i);

        // Calculate Lorentz factor
        Real u0 = std::sqrt(1.0 + SQR(u1) + SQR(u2) + SQR(u3));

        // Calculate 4-magnetic field
        Real b0 = bb1*u1 + bb2*u2 + bb3*u3;
        Real b1 = (bb1 + b0 * u1) / u0;
        Real b2 = (bb2 + b0 * u2) / u0;
        Real b3 = (bb3 + b0 * u3) / u0;
        Real b_sq = -SQR(b0) + SQR(b1) + SQR(b2) + SQR(b3);

        // Set conserved quantities
        Real wtot_u02 = (rho + gamma_adi_red * pgas + b_sq) * u0 * u0;
        Real dd = rho * u0;
        Real ee = wtot_u02 - b0 * b0 - (pgas + 0.5*b_sq);
        Real m1 = wtot_u02 * u1 / u0 - b0 * b1;
        Real m2 = wtot_u02 * u2 / u0 - b0 * b2;
        Real m3 = wtot_u02 * u3 / u0 - b0 * b3;
        cons(IDN,k,j,i) = dd;
        cons(IEN,k,j,i) = ee;
        cons(IM1,k,j,i) = m1;
        cons(IM2,k,j,i) = m2;
        cons(IM3,k,j,i) = m3;
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for calculating relativistic fast wavespeeds in arbitrary coordinates
// Inputs:
//   wgas: gas enthalpy
//   pgas: gas pressure
//   u0, u1: contravariant components of 4-velocity
//   b_sq: b_\mu b^\mu
//   g00, g01, g11: contravariant components of metric (-1, 0, 1 in SR)
// Outputs:
//   p_lambda_plus: value set to most positive wavespeed
//   p_lambda_minus: value set to most negative wavespeed
// Notes:
//   Follows same general procedure as vchar() in phys.c in Harm.
//   Variables are named as though 1 is normal direction.
//   Same function as in adiabatic_mhd_gr.cpp.

void EquationOfState::FastMagnetosonicSpeedsGR(Real wgas, Real pgas, Real u0, Real u1,
                                               Real b_sq, Real g00, Real g01, Real g11,
                                               Real *p_lambda_plus,
                                               Real *p_lambda_minus) {
  // Calculate comoving fast magnetosonic speed
  const Real gamma_adi = gamma_;
  Real cs_sq = gamma_adi * pgas / wgas;
  Real va_sq = b_sq / (b_sq + wgas);
  Real cms_sq = cs_sq + va_sq - cs_sq * va_sq;

  // Set fast magnetosonic speeds in appropriate coordinates
  Real a = SQR(u0) - (g00 + SQR(u0)) * cms_sq;
  Real b = -2.0 * (u0 * u1 - (g01 + u0 * u1) * cms_sq);
  Real c = SQR(u1) - (g11 + SQR(u1)) * cms_sq;
  Real a1 = b / a;
  Real a0 = c / a;
  Real s = std::sqrt(std::max(SQR(a1) - 4.0 * a0, static_cast<Real>(0.0)));
  *p_lambda_plus = (a1 >= 0.0) ? -2.0 * a0 / (a1 + s) : (-a1 + s) / 2.0;
  *p_lambda_minus = (a1 >= 0.0) ? (-a1 - s) / 2.0 : -2.0 * a0 / (a1 - s);
  return;
}

namespace {

//----------------------------------------------------------------------------------------
// Function for calculating primitives in normal observer frame
// Inputs:
//   dd_vals: array of conserved densities
//   ee_vals: array of conserved energies
//   mm_vals: array of conserved momenta \mathcal{M}^2, M^i
//   bb_vals: array of magnetic fields \mathcal{B{^2, B^i
//   tt_vals: array of M_i B^i values
//   gamma_adi: ratio of specific heats
//   pgas_old: previous value of p_{gas} used to initialize iteration
//   k, j, i: indices of cell
// Outputs:
//   returned value: true for successful convergence, false otherwise
//   prim: all values set in given cell
//   p_gamma_lor: normal-frame Lorentz factor
//   p_pmag: magnetic pressure
// Notes:
//   Implements Newman & Hamlin 2014, SIAM J. Sci. Comput. 36(4) B661 (NH).
//   Notation here largely follows (NH).
//   Symbols:
//     tt: \mathcal{T}
//     ee: E (NH: e)
//     mm: M (NH: m)
//     eee: \mathcal{E}
//     ll: \mathcal{L}
//     wgas: w_{gas} (NH: w)
//     rr: \mathcal{R}
//   Simpler version without magnetic fields found in adiabatic_hydro_sr.cpp.
//   Also found in adiabatic_mhd_gr.cpp.

bool ConservedToPrimitiveNormal(
    const AthenaArray<Real> &dd_vals, const AthenaArray<Real> &ee_vals,
    const AthenaArray<Real> &mm_vals, const AthenaArray<Real> &bb_vals,
    const AthenaArray<Real> &tt_vals, Real gamma_adi, Real pgas_old, int k, int j, int i,
    AthenaArray<Real> &prim, Real *p_gamma_lor, Real *p_pmag) {
  // Parameters
  const int max_iterations = 15;
  const Real tol = 1.0e-12;
  const Real pgas_uniform_min = 1.0e-12;
  const Real a_min = 1.0e-12;
  const Real v_sq_max = 1.0 - 1.0e-12;
  const Real rr_max = 1.0 - 1.0e-12;

  // Extract conserved values
  const Real &dd = dd_vals(i);
  const Real &ee = ee_vals(i);
  const Real &mm_sq = mm_vals(0,i);
  const Real &mm1 = mm_vals(1,i);
  const Real &mm2 = mm_vals(2,i);
  const Real &mm3 = mm_vals(3,i);
  const Real &bb_sq = bb_vals(0,i);
  const Real &bb1 = bb_vals(1,i);
  const Real &bb2 = bb_vals(2,i);
  const Real &bb3 = bb_vals(3,i);
  const Real &tt = tt_vals(i);

  // Calculate functions of conserved quantities
  Real d = 0.5 * (mm_sq * bb_sq - SQR(tt));                  // (NH 5.7)
  d = std::max(d, 0.0);
  Real pgas_min = std::cbrt(27.0/4.0 * d) - ee - 0.5*bb_sq;
  pgas_min = std::max(pgas_min, pgas_uniform_min);

  // Iterate until convergence
  Real pgas[3];
  pgas[0] = std::max(pgas_old, pgas_min);
  int n;
  for (n = 0; n < max_iterations; ++n) {
    // Step 1: Calculate cubic coefficients
    Real a;
    if (n%3 != 2) {
      a = ee + pgas[n%3] + 0.5*bb_sq;  // (NH 5.7)
      a = std::max(a, a_min);
    }

    // Step 2: Calculate correct root of cubic equation
    Real phi, eee, ll, v_sq;
    if (n%3 != 2) {
      phi = std::acos(1.0/a * std::sqrt(27.0*d/(4.0*a)));                     // (NH 5.10)
      eee = a/3.0 - 2.0/3.0 * a * std::cos(2.0/3.0 * (phi+PI));               // (NH 5.11)
      ll = eee - bb_sq;                                                       // (NH 5.5)
      v_sq = (mm_sq*SQR(ll) + SQR(tt)*(bb_sq+2.0*ll)) / SQR(ll * (bb_sq+ll)); // (NH 5.2)
      v_sq = std::min(std::max(v_sq, static_cast<Real>(0.0)), v_sq_max);
      Real gamma_sq = 1.0/(1.0-v_sq);                                         // (NH 3.1)
      Real gamma = std::sqrt(gamma_sq);                                       // (NH 3.1)
      Real wgas = ll/gamma_sq;                                                // (NH 5.1)
      Real rho = dd/gamma;                                                    // (NH 4.5)
      pgas[(n+1)%3] = (gamma_adi-1.0)/gamma_adi * (wgas - rho);               // (NH 4.1)
      pgas[(n+1)%3] = std::max(pgas[(n+1)%3], pgas_min);
    }

    // Step 3: Check for convergence
    if (n%3 != 2) {
      if (pgas[(n+1)%3] > pgas_min && std::abs(pgas[(n+1)%3]-pgas[n%3]) < tol) {
        break;
      }
    }

    // Step 4: Calculate Aitken accelerant and check for convergence
    if (n%3 == 2) {
      Real rr = (pgas[2] - pgas[1]) / (pgas[1] - pgas[0]);  // (NH 7.1)
      if (!std::isfinite(rr) || std::abs(rr) > rr_max) {
        continue;
      }
      pgas[0] = pgas[1] + (pgas[2] - pgas[1]) / (1.0 - rr);  // (NH 7.2)
      pgas[0] = std::max(pgas[0], pgas_min);
      if (pgas[0] > pgas_min && std::abs(pgas[0]-pgas[2]) < tol) {
        break;
      }
    }
  }

  // Step 5: Set primitives
  if (n == max_iterations) {
    return false;
  }
  prim(IPR,k,j,i) = pgas[(n+1)%3];
  if (!std::isfinite(prim(IPR,k,j,i))) {
    return false;
  }
  Real a = ee + prim(IPR,k,j,i) + 0.5*bb_sq;                      // (NH 5.7)
  a = std::max(a, a_min);
  Real phi = std::acos(1.0/a * std::sqrt(27.0*d/(4.0*a)));        // (NH 5.10)
  Real eee = a/3.0 - 2.0/3.0 * a * std::cos(2.0/3.0 * (phi+PI));  // (NH 5.11)
  Real ll = eee - bb_sq;                                          // (NH 5.5)
  Real v_sq = (mm_sq*SQR(ll) + SQR(tt)*(bb_sq+2.0*ll))
              / SQR(ll * (bb_sq+ll));                             // (NH 5.2)
  v_sq = std::min(std::max(v_sq, static_cast<Real>(0.0)), v_sq_max);
  Real gamma_sq = 1.0/(1.0-v_sq);                                 // (NH 3.1)
  Real gamma = std::sqrt(gamma_sq);                               // (NH 3.1)
  prim(IDN,k,j,i) = dd/gamma;                                     // (NH 4.5)
  if (!std::isfinite(prim(IDN,k,j,i))) {
    return false;
  }
  Real ss = tt/ll;                          // (NH 4.8)
  Real v1 = (mm1 + ss*bb1) / (ll + bb_sq);  // (NH 4.6)
  Real v2 = (mm2 + ss*bb2) / (ll + bb_sq);  // (NH 4.6)
  Real v3 = (mm3 + ss*bb3) / (ll + bb_sq);  // (NH 4.6)
  prim(IVX,k,j,i) = gamma*v1;               // (NH 3.3)
  prim(IVY,k,j,i) = gamma*v2;               // (NH 3.3)
  prim(IVZ,k,j,i) = gamma*v3;               // (NH 3.3)
  if (!std::isfinite(prim(IVX,k,j,i))
      || !std::isfinite(prim(IVY,k,j,i))
      || !std::isfinite(prim(IVZ,k,j,i))) {
    return false;
  }
  *p_gamma_lor = gamma;
  *p_pmag = 0.5 * (bb_sq/gamma_sq + SQR(ss));  // (NH 3.7, 3.11)
  return true;
}

} // namespace

//---------------------------------------------------------------------------------------
// \!fn void EquationOfState::ApplyPrimitiveFloors(AthenaArray<Real> &prim, int k, int j,
//                                                 int i)
// \brief Apply density and pressure floors to reconstructed L/R cell interface states

void EquationOfState::ApplyPrimitiveFloors(AthenaArray<Real> &prim, int k, int j, int i) {
  Real& w_d  = prim(IDN,i);
  Real& w_p  = prim(IPR,i);
  // Eventually, may want to check that small field errors don't overwhelm gas floor
  // Real density_floor_local = density_floor_;
  // if (sigma_max_ > 0.0) {
  //   density_floor_local = std::max(density_floor_local, 2.0*pmag/sigma_max_);
  // }
  // Real pressure_floor_local = pressure_floor_;
  // if (beta_min_ > 0.0) {
  //   pressure_floor_local = std::max(pressure_floor_local, beta_min_*pmag);
  // }
  // w_d = (w_d > density_floor_local) ?  w_d : density_floor_local;
  // w_p = (w_p > pressure_floor_local) ?  w_p : pressure_floor_local;

  // apply density floor
  w_d = (w_d > density_floor_) ?  w_d : density_floor_;
  // apply pressure floor
  w_p = (w_p > pressure_floor_) ?  w_p : pressure_floor_;

  return;
}
