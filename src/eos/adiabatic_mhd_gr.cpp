//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file adiabatic_mhd_gr.cpp
//  \brief Implements functions for going between primitive and conserved variables in
//  general-relativistic MHD, as well as for computing wavespeeds.

// C headers

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // abs(), cbrt(), isfinite(), isnan(), NAN, pow(), sqrt()
#include <limits>

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
void CalculateNormalConserved(
    const AthenaArray<Real> &cons,
    const AthenaArray<Real> &bb,
    const AthenaArray<Real> &g, const AthenaArray<Real> &gi,
    int k, int j, int il, int iu,
    AthenaArray<Real> &dd, AthenaArray<Real> &ee,
    AthenaArray<Real> &mm, AthenaArray<Real> &bbb,
    AthenaArray<Real> &tt);
#pragma omp declare simd simdlen(SIMD_WIDTH)                            \
  uniform(num_iterations,dd_vals,ee_vals,mm_vals,                       \
          bb_vals,tt_vals,gamma_adi,k,j,prim,gamma_vals,pmag_vals) linear(i)
bool ConservedToPrimitiveNormal(
    const int num_iterations,
    const AthenaArray<Real> &dd_vals,
    const AthenaArray<Real> &ee_vals,
    const AthenaArray<Real> &mm_vals,
    const AthenaArray<Real> &bb_vals,
    const AthenaArray<Real> &tt_vals,
    Real gamma_adi, Real pgas_old,
    int k, int j, int i, AthenaArray<Real> &prim,
    AthenaArray<Real> &gamma_vals,
    AthenaArray<Real> &pmag_vals);
#pragma omp declare simd simdlen(SIMD_WIDTH)                    \
  uniform(prim,gamma_adi,bb_cc,g, gi,k,j,cons,pco) linear(i)
void PrimitiveToConservedSingle(
    const AthenaArray<Real> &prim, Real gamma_adi,
    const AthenaArray<Real> &bb_cc,
    const AthenaArray<Real> &g, const AthenaArray<Real> &gi,
    int k, int j, int i,
    AthenaArray<Real> &cons, Coordinates *pco);
} // namespace

//----------------------------------------------------------------------------------------
// Constructor
// Inputs:
//   pmb: pointer to MeshBlock
//   pin: pointer to runtime inputs

EquationOfState::EquationOfState(MeshBlock *pmb, ParameterInput *pin) {
  pmy_block_ = pmb;
  gamma_ = pin->GetReal("hydro", "gamma");
  Real float_min = std::numeric_limits<float>::min();
  density_floor_ = pin->GetOrAddReal("hydro", "dfloor", std::sqrt(1024*(float_min)) );
  pressure_floor_ = pin->GetOrAddReal("hydro", "pfloor", std::sqrt(1024*(float_min)) );
  rho_min_ = pin->GetOrAddReal("hydro", "rho_min", density_floor_);
  rho_pow_ = pin->GetOrAddReal("hydro", "rho_pow", 0.0);
  pgas_min_ = pin->GetOrAddReal("hydro", "pgas_min", pressure_floor_);
  pgas_pow_ = pin->GetOrAddReal("hydro", "pgas_pow", 0.0);
  sigma_max_ = pin->GetOrAddReal("hydro", "sigma_max",  0.0);
  beta_min_ = pin->GetOrAddReal("hydro", "beta_min", 0.0);
  gamma_max_ = pin->GetOrAddReal("hydro", "gamma_max", 1000.0);
  int nc1 = pmb->ncells1;
  g_.NewAthenaArray(NMETRIC, nc1);
  g_inv_.NewAthenaArray(NMETRIC, nc1);
  fixed_.NewAthenaArray(nc1);
  success_.NewAthenaArray(nc1);
  normal_dd_.NewAthenaArray(nc1);
  normal_ee_.NewAthenaArray(nc1);
  normal_mm_.NewAthenaArray(4, nc1);
  normal_bb_.NewAthenaArray(4, nc1);
  normal_tt_.NewAthenaArray(nc1);
  dens_floor_local_.NewAthenaArray(nc1);
  press_floor_local_.NewAthenaArray(nc1);
  normal_gamma_.NewAthenaArray(nc1);
  pmag_.NewAthenaArray(nc1);
}

//----------------------------------------------------------------------------------------
// Variable inverter
// Inputs:
//   cons: conserved quantities
//   prim_old: primitive quantities from previous half timestep
//   bb: face-centered magnetic field
//   pco: pointer to Coordinates
//   il,iu,jl,ju,kl,ku: index bounds of region to be updated
// Outputs:
//   prim: primitives
//   bb_cc: cell-centered magnetic field

void EquationOfState::ConservedToPrimitive(
    AthenaArray<Real> &cons, const AthenaArray<Real> &prim_old, const FaceField &bb,
    AthenaArray<Real> &prim, AthenaArray<Real> &bb_cc, Coordinates *pco,
    int il, int iu, int jl, int ju, int kl, int ku) {

  // Parameters
  const Real mm_sq_ee_sq_max = 1.0 - 1.0e-12;  // max. of squared momentum over energy

  // Extract ratio of specific heats
  const Real gamma_adi = gamma_;
  const Real gamma_prime = gamma_adi/(gamma_adi-1.0);

  // Interpolate magnetic field from faces to cell centers
  pmy_block_->pfield->CalculateCellCenteredField(bb, bb_cc, pco, il, iu, jl, ju, kl, ku);

  // Go through all rows
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      // Calculate metric
      pco->CellMetric(k, j, il, iu, g_, g_inv_);

      // Cast problem into normal frame
      CalculateNormalConserved(cons, bb_cc, g_, g_inv_, k, j, il, iu, normal_dd_,
                               normal_ee_, normal_mm_, normal_bb_, normal_tt_);

      // Pass 1: Initial attempt at inversion
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        // Set flag indicating conserved values need adjusting at end
        fixed_(i) = false;

        // Calculate floors for density and pressure
        dens_floor_local_(i) = density_floor_;
        if (rho_pow_ != 0.0) {
          dens_floor_local_(i) =
              std::max(dens_floor_local_(i), rho_min_ * std::pow(pco->x1v(i), rho_pow_));
        }
        press_floor_local_(i) = pressure_floor_;
        if (pgas_pow_ != 0.0) {
          press_floor_local_(i) = std::max(press_floor_local_(i),
                                           pgas_min_ * std::pow(pco->x1v(i), pgas_pow_));
        }

        // Ensure conserved density is large enough
        Real dd_min = dens_floor_local_(i);
        if (normal_dd_(i) < dd_min) {
          normal_dd_(i) = dd_min;
          fixed_(i) = true;
        }

        // Ensure conserved energy is large enough
        Real ee_min = dens_floor_local_(i) + press_floor_local_(i)/(gamma_adi-1.0)
                      + 0.5*normal_bb_(0,i);
        if (normal_ee_(i) < ee_min) {
          normal_ee_(i) = ee_min;
          fixed_(i) = true;
        }

        // Ensure conserved momentum is not too large given energy
        Real mm_sq_max = mm_sq_ee_sq_max * SQR(normal_ee_(i));
        Real factor = std::sqrt(mm_sq_max/normal_mm_(0,i));
        if (normal_mm_(0,i) > mm_sq_max) {
          normal_mm_(0,i) = mm_sq_max;
          normal_mm_(1,i) *= factor;
          normal_mm_(2,i) *= factor;
          normal_mm_(3,i) *= factor;
          normal_tt_(i) *= factor;
          fixed_(i) = true;
        }

        // Set primitives
        success_(i) = ConservedToPrimitiveNormal(3, normal_dd_, normal_ee_, normal_mm_,
                                                 normal_bb_, normal_tt_, gamma_adi,
                                                 prim_old(IPR,k,j,i),
                                                 k, j, i, prim, normal_gamma_, pmag_);
      }

      // Pass 2: Cleanup (most cells should be skipped)
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        // Reapply iteration procedure in case convergence not attained (possibly slow)
        if (!success_(i)) {
          success_(i) = ConservedToPrimitiveNormal(7, normal_dd_, normal_ee_, normal_mm_,
                                                   normal_bb_, normal_tt_, gamma_adi,
                                                   prim(IPR,k,j,i),
                                                   k, j, i, prim, normal_gamma_, pmag_);
        }
      }

      // Pass 3: Density and pressure floors, and possible re-inversion
//#pragma omp simd simdlen(SIMD_WIDTH) // this leads to memory corruption
      for (int i=il; i<=iu; ++i) {
        // Handle failures
        if (!success_(i)) {
          for (int n = 0; n < NHYDRO; ++n) {
            prim(n,k,j,i) = prim_old(n,k,j,i);
          }
          fixed_(i) = true;
        }

        // Apply primitive density and gas pressure floors
        if (sigma_max_ > 0.0) {
          dens_floor_local_(i) = std::max(dens_floor_local_(i), 2.0*pmag_(i)/sigma_max_);
        }
        if (beta_min_ > 0.0) {
          press_floor_local_(i) = std::max(press_floor_local_(i), beta_min_*pmag_(i));
        }
        Real rho_add = std::max(dens_floor_local_(i) - prim(IDN,k,j,i), 0.0);
        Real pgas_add = std::max(press_floor_local_(i) - prim(IPR,k,j,i), 0.0);
        Real wgas_add = rho_add + gamma_prime * pgas_add;

        // Adjust normal frame conserved quantities, and recalculate primitives
        if (success_(i) && (rho_add > 0.0 || pgas_add > 0.0)) {
          // Adjust conserved density and energy
          normal_dd_(i) += rho_add * normal_gamma_(i);
          normal_ee_(i) += wgas_add * SQR(normal_gamma_(i)) + pgas_add;

          // Recalculate primitives
          success_(i) = ConservedToPrimitiveNormal(10, normal_dd_, normal_ee_, normal_mm_,
                                                   normal_bb_, normal_tt_, gamma_adi,
                                                   prim_old(IPR,k,j,i), k, j, i, prim,
                                                   normal_gamma_, pmag_);

          // Handle failures
          if (!success_(i)) {
            for (int n = 0; n < NHYDRO; ++n) {
              prim(n,k,j,i) = prim_old(n,k,j,i);
            }
          }
          fixed_(i) = true;
        }
      }

      // Pass 4: Velocity ceiling
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        // Apply velocity ceiling
        Real &uu1 = prim(IVX,k,j,i);
        Real &uu2 = prim(IVY,k,j,i);
        Real &uu3 = prim(IVZ,k,j,i);
        Real tmp = g_(I11,i)*SQR(uu1) + 2.0*g_(I12,i)*uu1*uu2 + 2.0*g_(I13,i)*uu1*uu3
                   + g_(I22,i)*SQR(uu2) + 2.0*g_(I23,i)*uu2*uu3
                   + g_(I33,i)*SQR(uu3);
        if (!success_(i)) {
          normal_gamma_(i) = std::sqrt(1.0 + tmp);
        }
        bool velocity_ceiling = false;
        Real factor = std::sqrt((SQR(gamma_max_)-1.0) / (SQR(normal_gamma_(i))-1.0));
        if (normal_gamma_(i) > gamma_max_) {
          uu1 *= factor;
          uu2 *= factor;
          uu3 *= factor;
          fixed_(i) = true;
          velocity_ceiling = true;
        }

        // Calculate b^0 (only needed if velocity_ceiling == true)
        Real alpha = std::sqrt(-1.0/g_inv_(I00,i));
        Real u0 = normal_gamma_(i)/alpha;
        Real u1 = uu1 - alpha * normal_gamma_(i) * g_inv_(I01,i);
        Real u2 = uu2 - alpha * normal_gamma_(i) * g_inv_(I02,i);
        Real u3 = uu3 - alpha * normal_gamma_(i) * g_inv_(I03,i);
        const Real &bb1 = bb_cc(IB1,k,j,i);
        const Real &bb2 = bb_cc(IB2,k,j,i);
        const Real &bb3 = bb_cc(IB3,k,j,i);
        Real b0 = g_(I01,i)*u0*bb1 + g_(I02,i)*u0*bb2 + g_(I03,i)*u0*bb3
                  + g_(I11,i)*u1*bb1 + g_(I12,i)*u1*bb2 + g_(I13,i)*u1*bb3
                  + g_(I12,i)*u2*bb1 + g_(I22,i)*u2*bb2 + g_(I23,i)*u2*bb3
                  + g_(I13,i)*u3*bb1 + g_(I23,i)*u3*bb2 + g_(I33,i)*u3*bb3;

        // Recalculate density and pressure floors given new velocity
        if (velocity_ceiling) {
          pmag_(i) = 0.5 * (normal_bb_(0,i)/SQR(normal_gamma_(i)) + SQR(b0/u0));
        }
        dens_floor_local_(i) = density_floor_;
        if (rho_pow_ != 0.0) {
          dens_floor_local_(i) =
              std::max(dens_floor_local_(i), rho_min_ * std::pow(pco->x1v(i), rho_pow_));
        }
        if (sigma_max_ > 0.0) {
          dens_floor_local_(i) = std::max(dens_floor_local_(i), 2.0*pmag_(i)/sigma_max_);
        }
        press_floor_local_(i) = pressure_floor_;
        if (pgas_pow_ != 0.0) {
          press_floor_local_(i) = std::max(press_floor_local_(i),
                                           pgas_min_ * std::pow(pco->x1v(i), pgas_pow_));
        }
        if (beta_min_ > 0.0) {
          press_floor_local_(i) = std::max(press_floor_local_(i), beta_min_*pmag_(i));
        }

        // Apply density and gas pressure floors in fluid frame
        Real &rho = prim(IDN,k,j,i);
        Real &pgas = prim(IPR,k,j,i);
        if (rho < dens_floor_local_(i)) {
          rho = dens_floor_local_(i);
          fixed_(i) = true;
        }
        if (pgas < press_floor_local_(i)) {
          pgas = press_floor_local_(i);
          fixed_(i) = true;
        }
        if (!success_(i)) {
          rho = dens_floor_local_(i);
          pgas = press_floor_local_(i);
          uu1 = uu2 = uu3 = 0.0;
        }
      }

      // Pass 5: Consistency check
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        if (fixed_(i)) {
          PrimitiveToConservedSingle(prim, gamma_adi, bb_cc, g_, g_inv_, k, j, i, cons,
                                     pco);
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
//   bb_cc: cell-centered magnetic field
//   pco: pointer to Coordinates
//   il,iu,jl,ju,kl,ku: index bounds of region to be updated
// Outputs:
//   cons: conserved variables
// Notes:
//   single-cell function exists for other purposes; call made to that function rather
//       than having duplicate code

void EquationOfState::PrimitiveToConserved(
    const AthenaArray<Real> &prim,
    const AthenaArray<Real> &bb_cc,
    AthenaArray<Real> &cons, Coordinates *pco,
    int il, int iu, int jl,
    int ju, int kl, int ku) {
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      pco->CellMetric(k, j, il, iu, g_, g_inv_);
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        PrimitiveToConservedSingle(prim, gamma_, bb_cc, g_, g_inv_, k, j, i, cons, pco);
      }
    }
  }
  return;
}


namespace {

//----------------------------------------------------------------------------------------
// Function for casting quantities into normal observer frame
// Inputs:
//   cons: conserved quantities rho u^0, T^0_\mu
//   bb: cell-centered magnetic field B^i
//   g,gi: 1D arrays of metric covariant and contravariant coefficients
//   k,j,il,iu: indices and index bounds of 1D array to use
// Outputs:
//   dd: normal density D
//   ee: normal energy E
//   mm: squared normal momentum M^2 = g_{ij} M^i M^j and its components M^i
//   bbb: squared normal field \mathcal{B}^2 = g_{ij} \mathcal{B}^i \mathcal{B}^j and
//       its components \mathcal{B}^i
//   tt: projection of momentum onto field \mathcal{T} = g_{ij} M^i \mathcal{B}^j
// Notes:
//   references Noble et al. 2006, ApJ 641 626 (N)
//   symbols:
//     t: T
//     qq_n: Q \cdot n
//     bbb: \mathcal{B}
//     tt: \mathcal{T}

void CalculateNormalConserved(
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bb,
    const AthenaArray<Real> &g, const AthenaArray<Real> &gi,
    int k, int j, int il, int iu, AthenaArray<Real> &dd, AthenaArray<Real> &ee,
    AthenaArray<Real> &mm, AthenaArray<Real> &bbb, AthenaArray<Real> &tt) {
  // Go through row
#pragma omp simd simdlen(SIMD_WIDTH)
  for (int i=il; i<=iu; ++i) {
    // Extract metric
    Real g_11 = g(I11,i), g_12 = g(I12,i), g_13 = g(I13,i),
         g_21 = g(I12,i), g_22 = g(I22,i), g_23 = g(I23,i),
         g_31 = g(I13,i), g_32 = g(I23,i), g_33 = g(I33,i);
    Real g00 = gi(I00,i), g01 = gi(I01,i), g02 = gi(I02,i), g03 = gi(I03,i),
         g10 = gi(I01,i), g11 = gi(I11,i), g12 = gi(I12,i), g13 = gi(I13,i),
         g20 = gi(I02,i), g21 = gi(I12,i), g22 = gi(I22,i), g23 = gi(I23,i),
         g30 = gi(I03,i), g31 = gi(I13,i), g32 = gi(I23,i), g33 = gi(I33,i);

    // Calculate unit timelike normal
    const Real alpha = std::sqrt(-1.0/g00);
    const Real alpha_2 = alpha * alpha;

    // Calculate projection operator
    const Real j10 = g10 + alpha_2*g01*g00, j20 = g20 + alpha_2*g02*g00,
               j30 = g30 + alpha_2*g03*g00, j11 = g11 + alpha_2*g01*g01,
               j21 = g21 + alpha_2*g02*g01, j31 = g31 + alpha_2*g03*g01,
               j12 = g12 + alpha_2*g01*g02, j22 = g22 + alpha_2*g02*g02,
               j32 = g32 + alpha_2*g03*g02, j13 = g13 + alpha_2*g01*g03,
               j23 = g23 + alpha_2*g02*g03, j33 = g33 + alpha_2*g03*g03;

    // Extract conserved quantities
    Real rho_u0 = cons(IDN,k,j,i);
    Real t0_0 = cons(IEN,k,j,i);
    Real t0_1 = cons(IVX,k,j,i);
    Real t0_2 = cons(IVY,k,j,i);
    Real t0_3 = cons(IVZ,k,j,i);
    Real bb1 = bb(IB1,k,j,i);
    Real bb2 = bb(IB2,k,j,i);
    Real bb3 = bb(IB3,k,j,i);

    // Calculate projected momentum densities Q_\mu = -n_\nu T^\nu_\mu (N 17)
    const Real qq_n = -alpha_2*(t0_0*g00 + t0_1*g01 + t0_2*g02 + t0_3*g03);

    // Calculate projected momentum M^i = j^{i\mu} Q_\mu
    const Real mm1 = alpha*(j10*t0_0 + j11*t0_1 + j12*t0_2 + j13*t0_3);
    const Real mm2 = alpha*(j20*t0_0 + j21*t0_1 + j22*t0_2 + j23*t0_3);
    const Real mm3 = alpha*(j30*t0_0 + j31*t0_1 + j32*t0_2 + j33*t0_3);

    // Calculate projected field \mathcal{B}^i = alpha B^i (N 5)
    const Real bbb1 = alpha * bb1;
    const Real bbb2 = alpha * bb2;
    const Real bbb3 = alpha * bb3;

    // Set normal conserved quantities
    dd(i) = alpha * rho_u0;  // (N 21)
    ee(i) = -qq_n;
    mm(0,i) = g_11*SQR(mm1) + 2.0*g_12*mm1*mm2 + 2.0*g_13*mm1*mm3
              + g_22*SQR(mm2) + 2.0*g_23*mm2*mm3
              + g_33*SQR(mm3);
    mm(1,i) = mm1;
    mm(2,i) = mm2;
    mm(3,i) = mm3;
    bbb(0,i) = g_11*SQR(bbb1) + 2.0*g_12*bbb1*bbb2 + 2.0*g_13*bbb1*bbb3
               + g_22*SQR(bbb2) + 2.0*g_23*bbb2*bbb3
               + g_33*SQR(bbb3);
    bbb(1,i) = bbb1;
    bbb(2,i) = bbb2;
    bbb(3,i) = bbb3;
    tt(i) = g_11*mm1*bbb1 + g_12*mm1*bbb2 + g_13*mm1*bbb3
            + g_21*mm2*bbb1 + g_22*mm2*bbb2 + g_23*mm2*bbb3
            + g_31*mm3*bbb1 + g_32*mm3*bbb2 + g_33*mm3*bbb3;
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for calculating primitives in normal observer frame: cleanup passes
// Inputs:
//   dd_vals: array of conserved densities
//   ee_vals: array of conserved energies
//   mm_vals: array of conserved momenta \mathcal{M}^2, M^i
//   bb_vals: array of magnetic fields \mathcal{B}^2, B^i
//   tt_vals: array of M_i B^i values
//   gamma_adi: ratio of specific heats
//   pgas_old: previous value of p_{gas} used to initialize iteration
//   k,j,i: indices of cell
// Outputs:
//   returned value: true for successful convergence, false otherwise
//   prim: all values set in given cell
//   gamma_vals: normal-frame Lorentz factor set in given cell
//   pmag_vals: magnetic pressure set in given cell
// Notes:
//   generalizes Newman & Hamlin 2014, SIAM J. Sci. Comput. 36(4) B661 (NH)
//     like SR, but all 3-vector operations done with respect to g_{ij} rather than
//         \eta_{ij}
//   notation here largely follows (NH), so for example writing B^i for what is really
//       \mathcal{B}^i
//   symbols:
//     tt: \mathcal{T}
//     ee: E (NH: e)
//     mm: M (NH: m)
//     eee: \mathcal{E}
//     ll: \mathcal{L}
//     wgas: w_{gas} (NH: w)
//     rr: \mathcal{R}
//   same exact function as ConservedToPrimitiveNormalInitial(), except this does 7
//       iterations

bool ConservedToPrimitiveNormal(
    const int num_iterations,
    const AthenaArray<Real> &dd_vals,
    const AthenaArray<Real> &ee_vals,
    const AthenaArray<Real> &mm_vals,
    const AthenaArray<Real> &bb_vals,
    const AthenaArray<Real> &tt_vals,
    Real gamma_adi, Real pgas_old, int k, int j, int i,
    AthenaArray<Real> &prim, AthenaArray<Real> &gamma_vals,
    AthenaArray<Real> &pmag_vals) {
  // Parameters
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

  // Precompute function of adiabatic index
  Real gamma_prime_inv = (gamma_adi-1.0) / gamma_adi;

  // Calculate functions of conserved quantities
  Real d = 0.5 * (mm_sq * bb_sq - SQR(tt));                  // (NH 5.7)
  d = std::max(d, 0.0);
  Real pgas_min = std::cbrt(27.0/4.0 * d) - ee - 0.5*bb_sq;
  pgas_min = std::max(pgas_min, pgas_uniform_min);

  // Iterate until convergence
  Real pgas[3];
  pgas[0] = std::max(pgas_old, pgas_min);
  for (int n = 0; n < num_iterations; ++n) {
    // Iterate normally for 2 out of every 3 times
    if (n%3 != 2) {
      // Step 1: Calculate cubic coefficients
      Real a = ee + pgas[n%3] + 0.5*bb_sq;  // (NH 5.7)
      a = std::max(a, a_min);

      // Step 2: Calculate correct root of cubic equation
      Real phi = std::acos(1.0/a * std::sqrt(27.0*d/(4.0*a)));        // (NH 5.10)
      Real eee = a/3.0 - 2.0/3.0 * a * std::cos(2.0/3.0 * (phi+PI));  // (NH 5.11)
      Real ll = eee - bb_sq;                                          // (NH 5.5)
      Real v_sq = (mm_sq*SQR(ll) + SQR(tt)*(bb_sq+2.0*ll))
                  / SQR(ll * (bb_sq+ll));                             // (NH 5.2)
      v_sq = std::min(std::max(v_sq, 0.0), v_sq_max);
      Real gamma_sq = 1.0/(1.0-v_sq);                                 // (NH 3.1)
      Real gamma = std::sqrt(gamma_sq);                               // (NH 3.1)
      pgas[(n+1)%3] = gamma_prime_inv * (ll/gamma_sq - dd/gamma);     // (NH 4.1)
      pgas[(n+1)%3] = std::max(pgas[(n+1)%3], pgas_min);

      // Every third iteration, apply Aitken acceleration instead
    } else {
      // Step 4: Calculate Aitken accelerant and check for convergence
      Real rr = (pgas[2] - pgas[1]) / (pgas[1] - pgas[0]);  // (NH 7.1)
      if (std::abs(rr) < rr_max + tol) {
        pgas[0] = pgas[1] + (pgas[2] - pgas[1]) / (1.0 - rr);  // (NH 7.2)
        pgas[0] = std::max(pgas[0], pgas_min);
      }
    }
  }

  // Step 3: Check for convergence
  bool success=true;
  if (pgas[num_iterations%3] < pgas_min ||
      std::abs(pgas[num_iterations%3]-pgas[(num_iterations-1)%3]) > tol) {
    success = false;
  }

  // Step 5: Set primitives
  prim(IPR,k,j,i) = pgas[num_iterations%3];
  Real a = ee + prim(IPR,k,j,i) + 0.5*bb_sq;                      // (NH 5.7)
  a = std::max(a, a_min);
  Real phi = std::acos(1.0/a * std::sqrt(27.0*d/(4.0*a)));        // (NH 5.10)
  Real eee = a/3.0 - 2.0/3.0 * a * std::cos(2.0/3.0 * (phi+PI));  // (NH 5.11)
  Real ll = eee - bb_sq;                                          // (NH 5.5)
  Real v_sq = (mm_sq*SQR(ll) + SQR(tt)*(bb_sq+2.0*ll))
              / SQR(ll * (bb_sq+ll));                                     // (NH 5.2)
  v_sq = std::min(std::max(v_sq, 0.0), v_sq_max);
  Real gamma_sq = 1.0/(1.0-v_sq);                                 // (NH 3.1)
  Real gamma = std::sqrt(gamma_sq);                               // (NH 3.1)
  prim(IDN,k,j,i) = dd/gamma;
  Real ss = tt/ll;                                                // (NH 4.8)
  Real eee_inv = 1.0/eee;
  Real gamma_v1 = gamma * (mm1 + ss*bb1) * eee_inv;               // (NH 4.6, 5.5)
  Real gamma_v2 = gamma * (mm2 + ss*bb2) * eee_inv;               // (NH 4.6, 5.5)
  Real gamma_v3 = gamma * (mm3 + ss*bb3) * eee_inv;               // (NH 4.6, 5.5)
  prim(IVX,k,j,i) = gamma_v1;                                     // (NH 3.3)
  prim(IVY,k,j,i) = gamma_v2;                                     // (NH 3.3)
  prim(IVZ,k,j,i) = gamma_v3;                                     // (NH 3.3)
  if (!std::isfinite(prim(IPR,k,j,i))
      || !std::isfinite(prim(IDN,k,j,i))
      || !std::isfinite(prim(IVX,k,j,i))
      || !std::isfinite(prim(IVY,k,j,i))
      || !std::isfinite(prim(IVZ,k,j,i))) {
    success = false;
  }
  gamma_vals(i) = gamma;
  pmag_vals(i) = 0.5 * (bb_sq/gamma_sq + SQR(ss));                // (NH 3.7, 3.11)
  return success;
}

//----------------------------------------------------------------------------------------
// Function for converting primitives to conserved variables in a single cell
// Inputs:
//   prim: 1D array of primitives
//   gamma_adi: ratio of specific heats
//   bb_cc: 3D array of cell-centered magnetic field
//   g,gi: 1D arrays of metric covariant and contravariant coefficients
//   k,j,i: indices of cell
//   pco: pointer to Coordinates
// Outputs:
//   cons: conserved variables set in desired cell

void PrimitiveToConservedSingle(
    const AthenaArray<Real> &prim, Real gamma_adi, const AthenaArray<Real> &bb_cc,
    const AthenaArray<Real> &g, const AthenaArray<Real> &gi,
    int k, int j, int i, AthenaArray<Real> &cons, Coordinates *pco) {
  // Extract primitives and magnetic fields
  const Real &rho = prim(IDN,k,j,i);
  const Real &pgas = prim(IPR,k,j,i);
  const Real &uu1 = prim(IVX,k,j,i);
  const Real &uu2 = prim(IVY,k,j,i);
  const Real &uu3 = prim(IVZ,k,j,i);
  const Real &bb1 = bb_cc(IB1,k,j,i);
  const Real &bb2 = bb_cc(IB2,k,j,i);
  const Real &bb3 = bb_cc(IB3,k,j,i);

  Real gamma_prime =  gamma_adi/(gamma_adi-1.0);
  // Calculate 4-velocity
  Real alpha = std::sqrt(-1.0/gi(I00,i));
  Real tmp = g(I11,i)*uu1*uu1 + 2.0*g(I12,i)*uu1*uu2 + 2.0*g(I13,i)*uu1*uu3
             + g(I22,i)*uu2*uu2 + 2.0*g(I23,i)*uu2*uu3
             + g(I33,i)*uu3*uu3;
  Real gamma = std::sqrt(1.0 + tmp);
  Real u0 = gamma/alpha;
  Real u1 = uu1 - alpha * gamma * gi(I01,i);
  Real u2 = uu2 - alpha * gamma * gi(I02,i);
  Real u3 = uu3 - alpha * gamma * gi(I03,i);
  Real u_0, u_1, u_2, u_3;
  pco->LowerVectorCell(u0, u1, u2, u3, k, j, i, &u_0, &u_1, &u_2, &u_3);

  // Calculate 4-magnetic field
  Real b0 = g(I01,i)*u0*bb1 + g(I02,i)*u0*bb2 + g(I03,i)*u0*bb3
            + g(I11,i)*u1*bb1 + g(I12,i)*u1*bb2 + g(I13,i)*u1*bb3
            + g(I12,i)*u2*bb1 + g(I22,i)*u2*bb2 + g(I23,i)*u2*bb3
            + g(I13,i)*u3*bb1 + g(I23,i)*u3*bb2 + g(I33,i)*u3*bb3;
  Real b1 = (bb1 + b0 * u1) / u0;
  Real b2 = (bb2 + b0 * u2) / u0;
  Real b3 = (bb3 + b0 * u3) / u0;
  Real b_0, b_1, b_2, b_3;
  pco->LowerVectorCell(b0, b1, b2, b3, k, j, i, &b_0, &b_1, &b_2, &b_3);
  Real b_sq = b0*b_0 + b1*b_1 + b2*b_2 + b3*b_3;

  // Extract conserved quantities
  Real &rho_u0 = cons(IDN,k,j,i);
  Real &t0_0 = cons(IEN,k,j,i);
  Real &t0_1 = cons(IM1,k,j,i);
  Real &t0_2 = cons(IM2,k,j,i);
  Real &t0_3 = cons(IM3,k,j,i);

  // Set conserved quantities
  Real wtot_u0 = (rho + gamma_prime * pgas + b_sq) * u0;
  Real ptot = pgas + 0.5 * b_sq;
  rho_u0 = rho * u0;
  t0_0 = wtot_u0 * u_0 - b0 * b_0 + ptot;
  t0_1 = wtot_u0 * u_1 - b0 * b_1;
  t0_2 = wtot_u0 * u_2 - b0 * b_2;
  t0_3 = wtot_u0 * u_3 - b0 * b_3;
  return;
}
} // namespace

//----------------------------------------------------------------------------------------
// Function for calculating relativistic fast wavespeeds
// Inputs:
//   prim: 1D array of primitive states
//   bbx_vals: 1D array of B^x
//   k,j: x3- and x2-indices
//   il,iu: lower and upper x1-indices
//   ivx: type of interface (IVX for x1, IVY for x2, IVZ for x3)
// Outputs:
//   lambdas_p,lambdas_m: 1D arrays set to +/- wavespeeds
// Notes:
//   follows advice in NR for avoiding large cancellations in solving quadratics
//   uses same approximation as FastMagnetosonicSpeedsGR()
//   almost same function as in adiabatic_mhd_sr.cpp

void EquationOfState::FastMagnetosonicSpeedsSR(
    const AthenaArray<Real> &prim, const AthenaArray<Real> &bbx_vals,
    int k, int j, int il, int iu, int ivx,
    AthenaArray<Real> &lambdas_p, AthenaArray<Real> &lambdas_m) {

  // Calculate cyclic permutations of indices
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  // Calculate ratio of specific heats
  const Real gamma_adi = gamma_;
  const Real gamma_adi_red = gamma_adi/(gamma_adi-1.0);

  // Go through states
#pragma omp simd simdlen(SIMD_WIDTH)
  for (int i = il; i <= iu; ++i) {
    // Extract primitives
    Real rho = prim(IDN,i);
    Real pgas = prim(IPR,i);
    Real u1 = prim(ivx,i);
    Real u2 = prim(ivy,i);
    Real u3 = prim(ivz,i);
    Real bb1 = bbx_vals(i);
    Real bb2 = prim(IBY,i);
    Real bb3 = prim(IBZ,i);

    // Calculate Lorentz factor and (twice) magnetic pressure
    Real u0 = std::sqrt(1.0 + SQR(u1) + SQR(u2) + SQR(u3));
    Real b0 = bb1 * u1 + bb2 * u2 + bb3 * u3;
    Real b1 = (bb1 + b0 * u1) / u0;
    Real b2 = (bb2 + b0 * u2) / u0;
    Real b3 = (bb3 + b0 * u3) / u0;
    Real b_sq = -SQR(b0) + SQR(b1) + SQR(b2) + SQR(b3);

    // Calculate comoving fast magnetosonic speed
    Real w_gas = rho + gamma_adi_red * pgas;
    Real cs_sq = gamma_adi * pgas / w_gas;
    Real va_sq = b_sq / (w_gas + b_sq);
    Real cms_sq = cs_sq + va_sq - cs_sq * va_sq;

    // Set fast magnetosonic speeds in appropriate coordinates
    Real a = SQR(u0) - (SQR(u0) - 1.0) * cms_sq;
    Real b = -2.0 * u0 * u1 * (1.0 - cms_sq);
    Real c = SQR(u1) - (SQR(u1) + 1.0) * cms_sq;
    Real d = std::max(SQR(b) - 4.0 * a * c, 0.0);
    Real d_sqrt = std::sqrt(d);
    Real root_1 = (-b + d_sqrt) / (2.0*a);
    Real root_2 = (-b - d_sqrt) / (2.0*a);
    if (root_1 > root_2) {
      lambdas_p(i) = root_1;
      lambdas_m(i) = root_2;
    } else {
      lambdas_p(i) = root_2;
      lambdas_m(i) = root_1;
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for calculating relativistic fast wavespeeds in arbitrary coordinates
// Inputs:
//   rho_h: gas enthalpy
//   pgas: gas pressure
//   u0,u1: contravariant components of 4-velocity
//   b_sq: b_\mu b^\mu
//   g00,g01,g11: contravariant components of metric
// Outputs:
//   plambda_plus: value set to most positive wavespeed
//   plambda_minus: value set to most negative wavespeed
// Notes:
//   follows same general procedure as vchar() in phys.c in Harm
//   variables are named as though 1 is normal direction

void EquationOfState::FastMagnetosonicSpeedsGR(Real rho_h, Real pgas, Real u0, Real u1,
                                               Real b_sq, Real g00, Real g01, Real g11,
                                               Real *plambda_plus, Real *plambda_minus) {
  // Parameters and constants
  const Real gamma_adi = gamma_;

  // Calculate comoving fast magnetosonic speed
  Real cs_sq = gamma_adi * pgas / rho_h;
  Real va_sq = b_sq / (b_sq + rho_h);
  Real cms_sq = cs_sq + va_sq - cs_sq * va_sq;

  // Set fast magnetosonic speeds in appropriate coordinates
  Real a = SQR(u0) - (g00 + SQR(u0)) * cms_sq;
  Real b = -2.0 * (u0*u1 - (g01 + u0*u1) * cms_sq);
  Real c = SQR(u1) - (g11 + SQR(u1)) * cms_sq;
  Real d = std::max(SQR(b) - 4.0*a*c, 0.0);
  Real d_sqrt = std::sqrt(d);
  Real root_1 = (-b + d_sqrt) / (2.0*a);
  Real root_2 = (-b - d_sqrt) / (2.0*a);
  if (root_1 > root_2) {
    *plambda_plus = root_1;
    *plambda_minus = root_2;
  } else {
    *plambda_plus = root_2;
    *plambda_minus = root_1;
  }
  return;
}

//---------------------------------------------------------------------------------------
// \!fn void EquationOfState::ApplyPrimitiveFloors(AthenaArray<Real> &prim,
//           int k, int j, int i)
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

  // Not applying position-dependent floors here in GR, nor using rho_min
  // apply density floor
  w_d = (w_d > density_floor_) ?  w_d : density_floor_;
  // apply pressure floor
  w_p = (w_p > pressure_floor_) ?  w_p : pressure_floor_;

  return;
}
