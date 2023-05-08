//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file adiabatic_hydro_sr.cpp
//! \brief Implements functions for going between primitive and conserved variables in
//!  special-relativistic hydrodynamics, as well as for computing wavespeeds.

// C headers

// C++ headers
#include <algorithm>  // max, min
#include <cmath>      // abs, isfinite, sqrt

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
    const AthenaArray<Real> &mm_vals, Real gamma_adi, Real pgas_old, int k, int j, int i,
    AthenaArray<Real> &prim, Real *p_gamma_lor);
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
    gamma_max_{pin->GetOrAddReal("hydro", "gamma_max", 1000.0)} {
  int nc1 = pmb->ncells1;
  normal_dd_.NewAthenaArray(nc1);
  normal_ee_.NewAthenaArray(nc1);
  normal_mm_.NewAthenaArray(4,nc1);
}

//----------------------------------------------------------------------------------------
// Variable inverter
// Inputs:
//   cons: conserved quantities
//   prim_old: primitive quantities from previous half timestep (not used)
//   bb: face-centered magnetic field (not used)
//   pco: pointer to Coordinates
//   il, iu, jl, ju, kl, ku: index bounds of region to be updated
// Outputs:
//   prim: primitives
//   bb_cc: cell-centered magnetic field (not set)
// Notes:
//   More complex version with magnetic fields found in adiabatic_mhd_sr.cpp.
//   More complex version for GR found in adiabatic_hydro_gr.cpp.

void EquationOfState::ConservedToPrimitive(
    AthenaArray<Real> &cons, const AthenaArray<Real> &prim_old, const FaceField &bb,
    AthenaArray<Real> &prim, AthenaArray<Real> &bb_cc, Coordinates *pco, int il, int iu,
    int jl, int ju, int kl, int ku) {
  // Parameters
  const Real mm_sq_ee_sq_max = 1.0 - 1.0e-12;  // max. of squared momentum over energy

  // Extract ratio of specific heats
  const Real &gamma_adi = gamma_;

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
      }

      // Go through cells
      for (int i=il; i<=iu; ++i) {
        // Set flag indicating conserved values need adjusting at end
        bool fixed = false;

        // Ensure conserved density is large enough
        Real dd_min = density_floor_;
        if (normal_dd_(i) < dd_min) {
          normal_dd_(i) = dd_min;
          fixed = true;
        }

        // Ensure conserved energy is large enough
        Real ee_min = density_floor_ + pressure_floor_/(gamma_adi-1.0);
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
          fixed = true;
        }

        // Define fallback state
        Real rho_old = std::max(prim(IDN,k,j,i), density_floor_);
        Real pgas_old = std::max(prim(IPR,k,j,i), pressure_floor_);
        Real u1_old = prim(IVX,k,j,i);
        Real u2_old = prim(IVY,k,j,i);
        Real u3_old = prim(IVZ,k,j,i);

        // Set primitives
        Real gamma;
        bool success = ConservedToPrimitiveNormal(normal_dd_, normal_ee_, normal_mm_,
                                                  gamma_adi, pgas_old, k, j, i, prim,
                                                  &gamma);

        // Handle failures
        if (!success) {
          prim(IDN,k,j,i) = rho_old;
          prim(IPR,k,j,i) = pgas_old;
          prim(IVX,k,j,i) = u1_old;
          prim(IVY,k,j,i) = u2_old;
          prim(IVZ,k,j,i) = u3_old;
          fixed = true;
        }

        // Apply density and gas pressure floors in coordinate frame
        Real rho_add = std::max(density_floor_-prim(IDN,k,j,i),
                                                static_cast<Real>(0.0));
        Real pgas_add = std::max(pressure_floor_-prim(IPR,k,j,i),
                                                static_cast<Real>(0.0));
        if (success && (rho_add > 0.0 || pgas_add > 0.0)) {
          // Adjust conserved density and energy
          Real wgas_add = rho_add + gamma_adi/(gamma_adi-1.0) * pgas_add;
          normal_dd_(i) += rho_add * gamma;
          normal_ee_(i) += wgas_add * SQR(gamma) + pgas_add;

          // Recalculate primitives
          success = ConservedToPrimitiveNormal(normal_dd_, normal_ee_, normal_mm_,
                                               gamma_adi, pgas_old, k, j, i, prim,
                                               &gamma);

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

        // Apply density and gas pressure floors in fluid frame
        Real &rho = prim(IDN,k,j,i);
        Real &pgas = prim(IPR,k,j,i);
        if (rho < density_floor_) {
          rho = density_floor_;
          fixed = true;
        }
        if (pgas < pressure_floor_) {
          pgas = pressure_floor_;
          fixed = true;
        }
        if (!success) {
          rho = density_floor_;
          pgas = pressure_floor_;
          u1 = u2 = u3 = 0.0;
        }

        // Ensure conserved variables match primitives
        if (fixed) {
          Real wgas = rho + gamma_adi / (gamma_adi - 1.0) * pgas;
          gamma = std::sqrt(1.0 + SQR(u1) + SQR(u2) + SQR(u3));
          cons(IDN,k,j,i) = rho * gamma;
          cons(IEN,k,j,i) = wgas * SQR(gamma) - pgas;
          cons(IM1,k,j,i) = wgas * gamma * u1;
          cons(IM2,k,j,i) = wgas * gamma * u2;
          cons(IM3,k,j,i) = wgas * gamma * u3;
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for converting all primitives to conserved variables
// Inputs:
//   prim: primitives
//   bb_cc: cell-centered magnetic field (unused)
//   pco: pointer to Coordinates
//   il,iu,jl,ju,kl,ku: index bounds of region to be updated
// Outputs:
//   cons: conserved variables

void EquationOfState::PrimitiveToConserved(
    const AthenaArray<Real> &prim, const AthenaArray<Real> &bb_cc,
    AthenaArray<Real> &cons, Coordinates *pco,
    int il, int iu, int jl, int ju, int kl, int ku) {
  // Calculate reduced ratio of specific heats
  Real gamma_prime = gamma_/(gamma_-1.0);

  // Go through all cells
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        // Extract primitives
        Real rho = prim(IDN,k,j,i);
        Real pgas = prim(IPR,k,j,i);
        Real u1 = prim(IVX,k,j,i);
        Real u2 = prim(IVY,k,j,i);
        Real u3 = prim(IVZ,k,j,i);

        // Calculate Lorentz factor
        Real u0 = std::sqrt(1.0 + SQR(u1) + SQR(u2) + SQR(u3));

        // Set conserved quantities
        Real wgas_u0 = (rho + gamma_prime * pgas) * u0;
        Real d = rho * u0;
        Real e = wgas_u0 * u0 - pgas;
        Real m1 = wgas_u0 * u1;
        Real m2 = wgas_u0 * u2;
        Real m3 = wgas_u0 * u3;
        cons(IDN,k,j,i) = d;
        cons(IEN,k,j,i) = e;
        cons(IM1,k,j,i) = m1;
        cons(IM2,k,j,i) = m2;
        cons(IM3,k,j,i) = m3;
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for calculating relativistic sound speeds
// Inputs:
//   rho_h: enthalpy per unit volume
//   pgas: gas pressure
//   vx: 3-velocity component v^x
//   gamma_lorentz_sq: Lorentz factor \gamma^2
// Outputs:
//   plambda_plus: value set to most positive wavespeed
//   plambda_minus: value set to most negative wavespeed
// Notes:
//   Same function as in adiabatic_hydro_gr.cpp.
//   References Mignone & Bodo 2005, MNRAS 364 126 (MB).

void EquationOfState::SoundSpeedsSR(Real rho_h, Real pgas, Real vx, Real gamma_lorentz_sq,
                                    Real *plambda_plus, Real *plambda_minus) {
  const Real gamma_adi = gamma_;
  Real cs_sq = gamma_adi * pgas / rho_h;  // (MB 4)
  Real sigma_s = cs_sq / (gamma_lorentz_sq * (1.0-cs_sq));
  Real relative_speed = std::sqrt(sigma_s * (1.0 + sigma_s - SQR(vx)));
  Real sigma_s_tmp = 1.0/(1.0+sigma_s);
  *plambda_plus = sigma_s_tmp * (vx + relative_speed);  // (MB 23)
  *plambda_minus = sigma_s_tmp * (vx - relative_speed);  // (MB 23)
  return;
}

namespace {

//----------------------------------------------------------------------------------------
// Function for calculating primitives
// Inputs:
//   dd_vals: array of conserved densities
//   ee_vals: array of conserved energies
//   mm_vals: array of conserved momenta \mathcal{M}^2, M^i
//   gamma_adi: ratio of specific heats
//   pgas_old: previous value of p_{gas} used to initialize iteration
//   k, j, i: indices of cell
// Outputs:
//   returned value: true for successful convergence, false otherwise
//   prim: all values set in given cell
//   p_gamma_lor: normal-frame Lorentz factor
// Notes:
//   Specializes Newman & Hamlin 2014, SIAM J. Sci. Comput. 36(4) B661 (NH).
//     Omits magnetic fields.
//   Notation here largely follows (NH).
//   Symbols:
//     ee: E (NH: e)
//     mm: M (NH: m)
//     wgas: w_{gas} (NH: w)
//     rr: \mathcal{R}
//   More complex version with magnetic fields found in adiabatic_mhd_sr.cpp.
//   Also found in adiabatic_hydro_gr.cpp.

bool ConservedToPrimitiveNormal(
    const AthenaArray<Real> &dd_vals, const AthenaArray<Real> &ee_vals,
    const AthenaArray<Real> &mm_vals, Real gamma_adi, Real pgas_old, int k, int j, int i,
    AthenaArray<Real> &prim, Real *p_gamma_lor) {
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

  // Calculate functions of conserved quantities
  Real pgas_min = -ee;
  pgas_min = std::max(pgas_min, pgas_uniform_min);

  // Iterate until convergence
  Real pgas[3];
  pgas[0] = std::max(pgas_old, pgas_min);
  int n;
  for (n = 0; n < max_iterations; ++n) {
    // Step 1: Calculate cubic coefficients
    Real a;
    if (n%3 != 2) {
      a = ee + pgas[n%3];      // (NH 5.7)
      a = std::max(a, a_min);
    }

    // Step 2: Calculate correct root of cubic equation
    Real v_sq;
    if (n%3 != 2) {
      v_sq = mm_sq / SQR(a);                                     // (NH 5.2)
      v_sq = std::min(std::max(v_sq, static_cast<Real>(0.0)), v_sq_max);
      Real gamma_sq = 1.0/(1.0-v_sq);                            // (NH 3.1)
      Real gamma = std::sqrt(gamma_sq);                          // (NH 3.1)
      Real wgas = a/gamma_sq;                                    // (NH 5.1)
      Real rho = dd/gamma;                                       // (NH 4.5)
      pgas[(n+1)%3] = (gamma_adi-1.0)/gamma_adi * (wgas - rho);  // (NH 4.1)
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
  Real a = ee + prim(IPR,k,j,i);                   // (NH 5.7)
  a = std::max(a, a_min);
  Real v_sq = mm_sq / SQR(a);                      // (NH 5.2)
  v_sq = std::min(std::max(v_sq, static_cast<Real>(0.0)), v_sq_max);
  Real gamma_sq = 1.0/(1.0-v_sq);                  // (NH 3.1)
  Real gamma = std::sqrt(gamma_sq);                // (NH 3.1)
  prim(IDN,k,j,i) = dd/gamma;                      // (NH 4.5)
  if (!std::isfinite(prim(IDN,k,j,i))) {
    return false;
  }
  Real v1 = mm1 / a;           // (NH 4.6)
  Real v2 = mm2 / a;           // (NH 4.6)
  Real v3 = mm3 / a;           // (NH 4.6)
  prim(IVX,k,j,i) = gamma*v1;  // (NH 3.3)
  prim(IVY,k,j,i) = gamma*v2;  // (NH 3.3)
  prim(IVZ,k,j,i) = gamma*v3;  // (NH 3.3)
  if (!std::isfinite(prim(IVX,k,j,i))
      || !std::isfinite(prim(IVY,k,j,i))
      || !std::isfinite(prim(IVZ,k,j,i))) {
    return false;
  }
  *p_gamma_lor = gamma;
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

  // apply density floor
  w_d = (w_d > density_floor_) ?  w_d : density_floor_;
  // apply pressure floor
  w_p = (w_p > pressure_floor_) ?  w_p : pressure_floor_;

  return;
}
