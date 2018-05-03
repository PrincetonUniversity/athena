//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file adiabatic_mhd_gr.cpp
//  \brief Implements functions for going between primitive and conserved variables in
//  general-relativistic MHD, as well as for computing wavespeeds.

// C++ headers
#include <algorithm>  // max(), min()
#include <cfloat>     // FLT_MIN
#include <cmath>      // abs(), cbrt(), isfinite(), isnan(), NAN, pow(), sqrt()

// Athena++ headers
#include "eos.hpp"
#include "../athena.hpp"                   // enums, macros
#include "../athena_arrays.hpp"            // AthenaArray
#include "../parameter_input.hpp"          // ParameterInput
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../field/field.hpp"              // FaceField
#include "../mesh/mesh.hpp"                // MeshBlock

// Declarations
static void CalculateNormalConserved(const AthenaArray<Real> &cons,
    const AthenaArray<Real> &bb, const AthenaArray<Real> &g, const AthenaArray<Real> &gi,
    int k, int j, int il, int iu, AthenaArray<Real> &dd, AthenaArray<Real> &ee,
    AthenaArray<Real> &mm, AthenaArray<Real> &bbb, AthenaArray<Real> &tt);
static bool ConservedToPrimitiveNormal(const AthenaArray<Real> &dd_vals,
    const AthenaArray<Real> &ee_vals, const AthenaArray<Real> &mm_vals,
    const AthenaArray<Real> &bb_vals, const AthenaArray<Real> &tt_vals, Real gamma_adi,
    Real pgas_old, int k, int j, int i, AthenaArray<Real> &prim, Real *p_gamma_lor,
    Real *p_pmag);
static void PrimitiveToConservedSingle(const AthenaArray<Real> &prim, Real gamma_adi,
    const AthenaArray<Real> &bb_cc, const AthenaArray<Real> &g,
    const AthenaArray<Real> &gi, int k, int j, int i, AthenaArray<Real> &cons,
    Coordinates *pco);

//----------------------------------------------------------------------------------------
// Constructor
// Inputs:
//   pmb: pointer to MeshBlock
//   pin: pointer to runtime inputs

EquationOfState::EquationOfState(MeshBlock *pmb, ParameterInput *pin) {
  pmy_block_ = pmb;
  gamma_ = pin->GetReal("hydro", "gamma");
  density_floor_ = pin->GetOrAddReal("hydro", "dfloor", std::sqrt(1024*(FLT_MIN)) );
  pressure_floor_ = pin->GetOrAddReal("hydro", "pfloor", std::sqrt(1024*(FLT_MIN)) );
  rho_min_ = pin->GetOrAddReal("hydro", "rho_min", density_floor_);
  rho_pow_ = pin->GetOrAddReal("hydro", "rho_pow", 0.0);
  pgas_min_ = pin->GetOrAddReal("hydro", "pgas_min", pressure_floor_);
  pgas_pow_ = pin->GetOrAddReal("hydro", "pgas_pow", 0.0);
  sigma_max_ = pin->GetOrAddReal("hydro", "sigma_max",  0.0);
  beta_min_ = pin->GetOrAddReal("hydro", "beta_min", 0.0);
  gamma_max_ = pin->GetOrAddReal("hydro", "gamma_max", 1000.0);
  int ncells1 = pmb->block_size.nx1 + 2*NGHOST;
  g_.NewAthenaArray(NMETRIC, ncells1);
  g_inv_.NewAthenaArray(NMETRIC, ncells1);
  normal_dd_.NewAthenaArray(ncells1);
  normal_ee_.NewAthenaArray(ncells1);
  normal_mm_.NewAthenaArray(4,ncells1);
  normal_bb_.NewAthenaArray(4,ncells1);
  normal_tt_.NewAthenaArray(ncells1);
}

//----------------------------------------------------------------------------------------
// Destructor

EquationOfState::~EquationOfState() {
  g_.DeleteAthenaArray();
  g_inv_.DeleteAthenaArray();
  normal_dd_.DeleteAthenaArray();
  normal_ee_.DeleteAthenaArray();
  normal_mm_.DeleteAthenaArray();
  normal_bb_.DeleteAthenaArray();
  normal_tt_.DeleteAthenaArray();
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

void EquationOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
    const AthenaArray<Real> &prim_old, const FaceField &bb, AthenaArray<Real> &prim,
    AthenaArray<Real> &bb_cc, Coordinates *pco, int il, int iu, int jl, int ju, int kl,
    int ku) {
  // Parameters
  const Real mm_sq_ee_sq_max = 1.0 - 1.0e-12;  // max. of squared momentum over energy

  // Extract ratio of specific heats
  const Real &gamma_adi = gamma_;

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

      // Go through cells
      for (int i=il; i<=iu; ++i) {

        // Set flag indicating conserved values need adjusting at end
        bool fixed = false;

        // Calculate floors for density and pressure
        Real density_floor_local = density_floor_;
        if (rho_pow_ != 0.0) {
          density_floor_local =
              std::max(density_floor_local, rho_min_ * std::pow(pco->x1v(i), rho_pow_));
        }
        Real pressure_floor_local = pressure_floor_;
        if (pgas_pow_ != 0.0) {
          pressure_floor_local = std::max(pressure_floor_local,
              pgas_min_ * std::pow(pco->x1v(i), pgas_pow_));
        }

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

        // Set primitives
        Real gamma, pmag;
        bool success = ConservedToPrimitiveNormal(normal_dd_, normal_ee_, normal_mm_,
            normal_bb_, normal_tt_, gamma_adi, prim_old(IPR,k,j,i), k, j, i, prim,
            &gamma, &pmag);

        // Handle failures
        if (not success) {
          for (int n = 0; n < NHYDRO; ++n) {
            prim(n,k,j,i) = prim_old(n,k,j,i);
          }
          fixed = true;
        }

        // Apply density and gas pressure floors in normal frame
        if (sigma_max_ > 0.0) {
          density_floor_local = std::max(density_floor_local, 2.0*pmag/sigma_max_);
        }
        if (beta_min_ > 0.0) {
          pressure_floor_local = std::max(pressure_floor_local, beta_min_*pmag);
        }
        Real rho_add = std::max(density_floor_local-prim(IDN,k,j,i), 0.0);
        Real pgas_add = std::max(pressure_floor_local-prim(IPR,k,j,i), 0.0);
        if (success and (rho_add > 0.0 or pgas_add > 0.0)) {

          // Adjust conserved density and energy
          Real wgas_add = rho_add + gamma_adi/(gamma_adi-1.0) * pgas_add;
          normal_dd_(i) += rho_add * gamma;
          normal_ee_(i) += wgas_add * SQR(gamma) + pgas_add;

          // Recalculate primitives
          success = ConservedToPrimitiveNormal(normal_dd_, normal_ee_, normal_mm_,
              normal_bb_, normal_tt_, gamma_adi, prim_old(IPR,k,j,i), k, j, i, prim,
              &gamma, &pmag);

          // Handle failures
          if (not success) {
            for (int n = 0; n < NHYDRO; ++n) {
              prim(n,k,j,i) = prim_old(n,k,j,i);
            }
          }
          fixed = true;
        }

        // Apply velocity ceiling
        Real &uu1 = prim(IVX,k,j,i);
        Real &uu2 = prim(IVY,k,j,i);
        Real &uu3 = prim(IVZ,k,j,i);
        if (not success) {
          Real tmp = g_(I11,i)*SQR(uu1) + 2.0*g_(I12,i)*uu1*uu2 + 2.0*g_(I13,i)*uu1*uu3
                   + g_(I22,i)*SQR(uu2) + 2.0*g_(I23,i)*uu2*uu3
                   + g_(I33,i)*SQR(uu3);
          gamma = std::sqrt(1.0 + tmp);
        }
        bool velocity_ceiling = false;
        if (gamma > gamma_max_) {
          Real factor = std::sqrt((SQR(gamma_max_)-1.0) / (SQR(gamma)-1.0));
          uu1 *= factor;
          uu2 *= factor;
          uu3 *= factor;
          fixed = true;
          velocity_ceiling = true;
        }

        // Recalculate density and pressure floors given new velocity
        if (velocity_ceiling) {
          Real alpha = std::sqrt(-1.0/g_inv_(I00,i));
          Real u0 = gamma/alpha;
          Real u1 = uu1 - alpha * gamma * g_inv_(I01,i);
          Real u2 = uu2 - alpha * gamma * g_inv_(I02,i);
          Real u3 = uu3 - alpha * gamma * g_inv_(I03,i);
          const Real &bb1 = bb_cc(IB1,k,j,i);
          const Real &bb2 = bb_cc(IB2,k,j,i);
          const Real &bb3 = bb_cc(IB3,k,j,i);
          Real b0 = g_(I01,i)*u0*bb1 + g_(I02,i)*u0*bb2 + g_(I03,i)*u0*bb3
                  + g_(I11,i)*u1*bb1 + g_(I12,i)*u1*bb2 + g_(I13,i)*u1*bb3
                  + g_(I12,i)*u2*bb1 + g_(I22,i)*u2*bb2 + g_(I23,i)*u2*bb3
                  + g_(I13,i)*u3*bb1 + g_(I23,i)*u3*bb2 + g_(I33,i)*u3*bb3;
          pmag = 0.5 * (normal_bb_(0,i)/SQR(gamma) + SQR(b0/u0));
        }
        density_floor_local = density_floor_;
        if (rho_pow_ != 0.0) {
          density_floor_local =
              std::max(density_floor_local, rho_min_ * std::pow(pco->x1v(i), rho_pow_));
        }
        if (sigma_max_ > 0.0) {
          density_floor_local = std::max(density_floor_local, 2.0*pmag/sigma_max_);
        }
        pressure_floor_local = pressure_floor_;
        if (pgas_pow_ != 0.0) {
          pressure_floor_local = std::max(pressure_floor_local,
              pgas_min_ * std::pow(pco->x1v(i), pgas_pow_));
        }
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
        if (not success) {
          rho = density_floor_local;
          pgas = pressure_floor_local;
          uu1 = uu2 = uu3 = 0.0;
        }

        // Ensure conserved variables match primitives
        if (fixed) {
          PrimitiveToConservedSingle(prim, gamma_adi, bb_cc, g_, g_inv_, k, j, i, cons,
              pco);
        }
      }
    }
  }
  return;
}

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

static void CalculateNormalConserved(const AthenaArray<Real> &cons,
    const AthenaArray<Real> &bb, const AthenaArray<Real> &g, const AthenaArray<Real> &gi,
    int k, int j, int il, int iu, AthenaArray<Real> &dd, AthenaArray<Real> &ee,
    AthenaArray<Real> &mm, AthenaArray<Real> &bbb, AthenaArray<Real> &tt) {
  // Go through row
  for (int i=il; i<=iu; ++i) {

    // Extract metric
    const Real &g_11 = g(I11,i), &g_12 = g(I12,i), &g_13 = g(I13,i),
               &g_21 = g(I12,i), &g_22 = g(I22,i), &g_23 = g(I23,i),
               &g_31 = g(I13,i), &g_32 = g(I23,i), &g_33 = g(I33,i);
    const Real &g00 = gi(I00,i), &g01 = gi(I01,i), &g02 = gi(I02,i), &g03 = gi(I03,i),
               &g10 = gi(I01,i), &g11 = gi(I11,i), &g12 = gi(I12,i), &g13 = gi(I13,i),
               &g20 = gi(I02,i), &g21 = gi(I12,i), &g22 = gi(I22,i), &g23 = gi(I23,i),
               &g30 = gi(I03,i), &g31 = gi(I13,i), &g32 = gi(I23,i), &g33 = gi(I33,i);

    // Calculate unit timelike normal
    const Real alpha = std::sqrt(-1.0/g00);
    const Real n0 = -alpha * g00;
    const Real n1 = -alpha * g01;
    const Real n2 = -alpha * g02;
    const Real n3 = -alpha * g03;

    // Calculate projection operator
    const Real j10 = g10 + n1*n0, j20 = g20 + n2*n0, j30 = g30 + n3*n0;
    const Real j11 = g11 + n1*n1, j21 = g21 + n2*n1, j31 = g31 + n3*n1;
    const Real j12 = g12 + n1*n2, j22 = g22 + n2*n2, j32 = g32 + n3*n2;
    const Real j13 = g13 + n1*n3, j23 = g23 + n2*n3, j33 = g33 + n3*n3;

    // Extract conserved quantities
    const Real &rho_u0 = cons(IDN,k,j,i);
    const Real &t0_0 = cons(IEN,k,j,i);
    const Real &t0_1 = cons(IVX,k,j,i);
    const Real &t0_2 = cons(IVY,k,j,i);
    const Real &t0_3 = cons(IVZ,k,j,i);
    const Real &bb1 = bb(IB1,k,j,i);
    const Real &bb2 = bb(IB2,k,j,i);
    const Real &bb3 = bb(IB3,k,j,i);

    // Calculate projected momentum densities Q_\mu = -n_\nu T^\nu_\mu (N 17)
    const Real qq_0 = alpha * t0_0;
    const Real qq_1 = alpha * t0_1;
    const Real qq_2 = alpha * t0_2;
    const Real qq_3 = alpha * t0_3;
    const Real qq_n = qq_0*n0 + qq_1*n1 + qq_2*n2 + qq_3*n3;

    // Calculate projected momentum M^i = j^{i\mu} Q_\mu
    const Real mm1 = j10*qq_0 + j11*qq_1 + j12*qq_2 + j13*qq_3;
    const Real mm2 = j20*qq_0 + j21*qq_1 + j22*qq_2 + j23*qq_3;
    const Real mm3 = j30*qq_0 + j31*qq_1 + j32*qq_2 + j33*qq_3;

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
// Function for calculating primitives in normal observer frame
// Inputs:
//   dd_vals: array of conserved densities
//   ee_vals: array of conserved energies
//   mm_vals: array of conserved momenta \mathcal{M}^2, M^i
//   bb_vals: array of magnetic fields \mathcal{B{^2, B^i
//   tt_vals: array of M_i B^i values
//   gamma_adi: ratio of specific heats
//   pgas_old: previous value of p_{gas} used to initialize iteration
//   k,j,i: indices of cell
// Outputs:
//   returned value: true for successful convergence, false otherwise
//   prim: all values set in given cell
//   p_gamma_lor: normal-frame Lorentz factor
//   p_pmag: magnetic pressure
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

static bool ConservedToPrimitiveNormal(const AthenaArray<Real> &dd_vals,
    const AthenaArray<Real> &ee_vals, const AthenaArray<Real> &mm_vals,
    const AthenaArray<Real> &bb_vals, const AthenaArray<Real> &tt_vals, Real gamma_adi,
    Real pgas_old, int k, int j, int i, AthenaArray<Real> &prim, Real *p_gamma_lor,
    Real *p_pmag) {
  // Parameters
  const int max_iterations = 10;
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
  Real d = 0.5 * (mm_sq * bb_sq - SQR(tt));             // (NH 5.7)
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
      v_sq = std::min(std::max(v_sq, 0.0), v_sq_max);
      Real gamma_sq = 1.0/(1.0-v_sq);                                         // (NH 3.1)
      Real gamma = std::sqrt(gamma_sq);                                       // (NH 3.1)
      Real wgas = ll/gamma_sq;                                                // (NH 5.1)
      Real rho = dd/gamma;                                                    // (NH 4.5)
      pgas[(n+1)%3] = (gamma_adi-1.0)/gamma_adi * (wgas - rho);               // (NH 4.1)
      pgas[(n+1)%3] = std::max(pgas[(n+1)%3], pgas_min);
    }

    // Step 3: Check for convergence
    if (n%3 != 2) {
      if (pgas[(n+1)%3] > pgas_min and std::abs(pgas[(n+1)%3]-pgas[n%3]) < tol) {
        break;
      }
    }

    // Step 4: Calculate Aitken accelerant and check for convergence
    if (n%3 == 2) {
      Real rr = (pgas[2] - pgas[1]) / (pgas[1] - pgas[0]);  // (NH 7.1)
      if (not std::isfinite(rr) or std::abs(rr) > rr_max) {
        continue;
      }
      pgas[0] = pgas[1] + (pgas[2] - pgas[1]) / (1.0 - rr);  // (NH 7.2)
      pgas[0] = std::max(pgas[0], pgas_min);
      if (pgas[0] > pgas_min and std::abs(pgas[0]-pgas[2]) < tol) {
        break;
      }
    }
  }

  // Step 5: Set primitives
  if (n == max_iterations) {
    return false;
  }
  prim(IPR,k,j,i) = pgas[(n+1)%3];
  if (not std::isfinite(prim(IPR,k,j,i))) {
    return false;
  }
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
  Real wgas = ll/gamma_sq;                                        // (NH 5.1)
  prim(IDN,k,j,i) = dd/gamma;                                     // (NH 4.5)
  if (not std::isfinite(prim(IDN,k,j,i))) {
    return false;
  }
  Real ss = tt/ll;                          // (NH 4.8)
  Real v1 = (mm1 + ss*bb1) / (ll + bb_sq);  // (NH 4.6)
  Real v2 = (mm2 + ss*bb2) / (ll + bb_sq);  // (NH 4.6)
  Real v3 = (mm3 + ss*bb3) / (ll + bb_sq);  // (NH 4.6)
  prim(IVX,k,j,i) = gamma*v1;               // (NH 3.3)
  prim(IVY,k,j,i) = gamma*v2;               // (NH 3.3)
  prim(IVZ,k,j,i) = gamma*v3;               // (NH 3.3)
  if (not std::isfinite(prim(IVX,k,j,i))
      or not std::isfinite(prim(IVY,k,j,i))
      or not std::isfinite(prim(IVZ,k,j,i))) {
    return false;
  }
  *p_gamma_lor = gamma;
  *p_pmag = 0.5 * (bb_sq/gamma_sq + SQR(ss));  // (NH 3.7, 3.11)
  return true;
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

void EquationOfState::PrimitiveToConserved(const AthenaArray<Real> &prim,
     const AthenaArray<Real> &bb_cc, AthenaArray<Real> &cons, Coordinates *pco,
     int il, int iu, int jl, int ju, int kl, int ku) {
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      pco->CellMetric(k, j, il, iu, g_, g_inv_);
      //#pragma omp simd // fn is too long to inline
      for (int i=il; i<=iu; ++i) {
        PrimitiveToConservedSingle(prim, gamma_, bb_cc, g_, g_inv_, k, j, i, cons, pco);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for converting primitives to conserved variables in a single cell
// Inputs:
//   prim: 3D array of primitives
//   gamma_adi: ratio of specific heats
//   bb_cc: 3D array of cell-centered magnetic field
//   g,gi: 1D arrays of metric covariant and contravariant coefficients
//   k,j,i: indices of cell
//   pco: pointer to Coordinates
// Outputs:
//   cons: conserved variables set in desired cell

static void PrimitiveToConservedSingle(const AthenaArray<Real> &prim, Real gamma_adi,
    const AthenaArray<Real> &bb_cc, const AthenaArray<Real> &g,
    const AthenaArray<Real> &gi, int k, int j, int i, AthenaArray<Real> &cons,
    Coordinates *pco) {
  // Extract primitives and magnetic fields
  const Real &rho = prim(IDN,k,j,i);
  const Real &pgas = prim(IPR,k,j,i);
  const Real &uu1 = prim(IVX,k,j,i);
  const Real &uu2 = prim(IVY,k,j,i);
  const Real &uu3 = prim(IVZ,k,j,i);
  const Real &bb1 = bb_cc(IB1,k,j,i);
  const Real &bb2 = bb_cc(IB2,k,j,i);
  const Real &bb3 = bb_cc(IB3,k,j,i);

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
  Real wtot = rho + gamma_adi/(gamma_adi-1.0) * pgas + b_sq;
  Real ptot = pgas + 0.5 * b_sq;
  rho_u0 = rho * u0;
  t0_0 = wtot * u0 * u_0 - b0 * b_0 + ptot;
  t0_1 = wtot * u0 * u_1 - b0 * b_1;
  t0_2 = wtot * u0 * u_2 - b0 * b_2;
  t0_3 = wtot * u0 * u_3 - b0 * b_3;
  return;
}

//----------------------------------------------------------------------------------------
// Function for calculating relativistic fast wavespeeds
// Inputs:
//   prim: 3D array of primitive states
//   bbx_vals: 1D array of B^x
//   k,j: x3- and x2-indices
//   il,iu: lower and upper x1-indices
//   ivx: type of interface (IVX for x1, IVY for x2, IVZ for x3)
// Outputs:
//   lambdas_p,lambdas_m: 1D arrays set to +/- wavespeeds
// Notes:
//   references Mignone & Bodo 2005, MNRAS 364 126 (MB2005)
//   references Mignone & Bodo 2006, MNRAS 368 1040 (MB2006)
//   references Numerical Recipes, 3rd ed. (NR)
//   follows advice in NR for avoiding large cancellations in solving quadratics
//   almost same function as in adiabatic_mhd_sr.cpp

void EquationOfState::FastMagnetosonicSpeedsSR(const AthenaArray<Real> &prim,
    const AthenaArray<Real> &bbx_vals, int k, int j, int il, int iu, int ivx,
    AthenaArray<Real> &lambdas_p, AthenaArray<Real> &lambdas_m) {
  // Parameters
  const double v_limit = 1.0e-12;  // squared velocities less than this are considered 0
  const double b_limit = 1.0e-14;  // squared B^x less than this is considered 0

  // Calculate cyclic permutations of indices
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  // Calculate ratio of specific heats
  const Real gamma_adi = gamma_;
  const Real gamma_adi_red = gamma_adi/(gamma_adi-1.0);

  // Go through states
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {

    // Extract primitives
    const Real &rho = prim(IDN,k,j,i);
    const Real &pgas = prim(IPR,k,j,i);
    Real u[4];
    u[1] = prim(ivx,k,j,i);
    u[2] = prim(ivy,k,j,i);
    u[3] = prim(ivz,k,j,i);
    u[0] = std::sqrt(1.0 + SQR(u[1]) + SQR(u[2]) + SQR(u[3]));
    const Real &bbx = bbx_vals(i);
    const Real &bby = prim(IBY,k,j,i);
    const Real &bbz = prim(IBZ,k,j,i);

    // Calculate 3-velocity
    Real vx = u[1]/u[0];
    Real vy = u[2]/u[0];
    Real vz = u[3]/u[0];

    // Calculate contravariant magnetic field
    Real b[4];
    b[0] = bbx*u[1] + bby*u[2] + bbz*u[3];
    b[1] = (bbx + b[0] * u[1]) / u[0];
    b[2] = (bby + b[0] * u[2]) / u[0];
    b[3] = (bbz + b[0] * u[3]) / u[0];

    // Calculate intermediate quantities
    Real v_sq = SQR(vx) + SQR(vy) + SQR(vz);
    Real gamma_rel_sq = 1.0/(1.0-v_sq);
    Real w_gas = rho + gamma_adi_red * pgas;
    Real cs_sq = gamma_adi * pgas / w_gas;                       // (MB2005 4)
    Real b_sq = -SQR(b[0]) + SQR(b[1]) + SQR(b[2]) + SQR(b[3]);
    Real bbx_sq = SQR(bbx);

    // Calculate wavespeeds in vanishing velocity case (MB2006 57)
    Real lambda_plus_no_v, lambda_minus_no_v;
    {
      Real w_tot = w_gas + b_sq;
      Real a1 = -(b_sq + cs_sq * (w_gas + bbx_sq)) / w_tot;
      Real a0 = cs_sq * bbx_sq / w_tot;
      Real s2 = SQR(a1) - 4.0*a0;
      Real s = (s2 > 0.0) ? std::sqrt(s2) : 0.0;
      Real lambda_sq = 0.5 * (-a1 + s);
      lambda_plus_no_v = std::sqrt(lambda_sq);
      lambda_minus_no_v = -lambda_plus_no_v;
    }

    // Calculate wavespeeds in vanishing normal field case (MB2006 58)
    Real lambda_plus_no_bbx, lambda_minus_no_bbx;
    {
      Real vx_sq = SQR(vx);
      Real v_dot_bb_perp = vy*bby + vz*bbz;
      Real q = b_sq - cs_sq*SQR(v_dot_bb_perp);
      Real denominator = w_gas * (cs_sq + gamma_rel_sq*(1.0-cs_sq)) + q;
      Real a1 = -2.0 * w_gas * gamma_rel_sq * vx * (1.0-cs_sq) / denominator;
      Real a0 = (w_gas * (-cs_sq + gamma_rel_sq*vx_sq*(1.0-cs_sq)) - q) / denominator;
      Real s2 = SQR(a1) - 4.0*a0;
      Real s = (s2 > 0.0) ? std::sqrt(s2) : 0.0;
      lambda_plus_no_bbx = (a1 >= 0.0) ? -2.0*a0/(a1+s) : (-a1+s)/2.0;
      lambda_minus_no_bbx = (a1 >= 0.0) ? (-a1-s)/2.0 : -2.0*a0/(a1-s);
    }

    // Calculate wavespeeds in general case (MB2006 56)
    Real lambda_plus, lambda_minus;
    {
      // Calculate quartic coefficients
      Real vx2 = SQR(vx);
      Real vx3 = vx2 * vx;
      Real vx4 = SQR(vx2);
      Real bt_sq = SQR(b[0]);
      Real bx_sq = SQR(b[1]);
      Real tmp1 = SQR(gamma_rel_sq) * w_gas * (1.0-cs_sq);
      Real tmp2 = gamma_rel_sq * (b_sq + w_gas * cs_sq);
      Real denominator = tmp1 + tmp2 - cs_sq * bt_sq;
      Real a3 = (-(4.0*tmp1+2.0*tmp2)*vx + 2.0*cs_sq*b[0]*b[1]) / denominator;
      Real a2 = (6.0*tmp1*vx2 + tmp2*(vx2-1.0) + cs_sq*(bt_sq-bx_sq)) / denominator;
      Real a1 = (-4.0*tmp1*vx3 + 2.0*tmp2*vx - 2.0*cs_sq*b[0]*b[1]) / denominator;
      Real a0 = (tmp1*vx4 - tmp2*vx2 + cs_sq*bx_sq) / denominator;

      // Calculate reduced quartic coefficients
      Real b2 = a2 - 3.0/8.0*SQR(a3);
      Real b1 = a1 - 1.0/2.0*a2*a3 + 1.0/8.0*a3*SQR(a3);
      Real b0 = a0 - 1.0/4.0*a1*a3 + 1.0/16.0*a2*SQR(a3) - 3.0/256.0*SQR(SQR(a3));

      // Solve reduced quartic equation
      Real y1, y2, y3, y4;
      {
        // Calculate resolvent cubic coefficients
        Real c2 = -b2;
        Real c1 = -4.0*b0;
        Real c0 = 4.0*b0*b2 - SQR(b1);

        // Solve resolvent cubic equation
        Real q = (c2*c2 - 3.0*c1) / 9.0;                       // (NR 5.6.10)
        Real r = (2.0*c2*c2*c2 - 9.0*c1*c2 + 27.0*c0) / 54.0;  // (NR 5.6.10)
        Real q3 = q*q*q;
        Real r2 = SQR(r);
        Real s2 = r2 - q3;
        Real z0;
        if (s2 < 0.0) {
          Real theta = std::acos(r/std::sqrt(q3));             // (NR 5.6.11)
          z0 = -2.0 * std::sqrt(q) * cos(theta/3.0) - c2/3.0;  // (NR 5.6.12)
        } else {
          Real s = std::sqrt(s2);
          Real aa = -copysign(1.0, r) * std::cbrt(std::abs(r) + s);  // (NR 5.6.15)
          Real bb = (aa != 0.0) ? q/aa : 0.0;                   // (NR 5.6.16)
          z0 = aa + bb - c2/3.0;
        }

        // Calculate quadratic coefficients
        Real d1 = (z0-b2 > 0.0) ? std::sqrt(z0-b2) : 0.0;
        Real e1 = -d1;
        s2 = SQR(z0)/4.0 - b0;
        Real s = (s2 < 0.0) ? 0.0 : std::sqrt(s2);
        Real d0 = (b1 < 0) ? 0.5*z0+s : 0.5*z0-s;
        Real e0 = (b1 < 0) ? 0.5*z0-s : 0.5*z0+s;

        // Solve quadratic equations
        s2 = SQR(d1) - 4.0*d0;
        s = (s2 < 0.0) ? 0.0 : std::sqrt(s2);
        y1 = (d1 >= 0.0) ? (-d1-s)/2.0 : -2.0*d0/(d1-s);
        y2 = (d1 >= 0.0) ? -2.0*d0/(d1+s) : (-d1+s)/2.0;
        s2 = SQR(e1) - 4.0*e0;
        s = (s2 < 0.0) ? 0.0 : std::sqrt(s2);
        y3 = (e1 >= 0.0) ? (-e1-s)/2.0 : -2.0*e0/(e1-s);
        y4 = (e1 >= 0.0) ? -2.0*e0/(e1+s) : (-e1+s)/2.0;
      }

      // Calculate extremal original quartic roots
      lambda_minus = std::min(y1, y3) - a3/4.0;
      lambda_plus = std::max(y2, y4) - a3/4.0;

      // Ensure wavespeeds are not superluminal
      if (not std::isfinite(lambda_minus) or lambda_minus < -1.0) {
        lambda_minus = -1.0;
      }
      if (not std::isfinite(lambda_plus) or lambda_plus > 1.0) {
        lambda_plus = 1.0;
      }
    }

    // Set wavespeeds based on velocity and magnetic field
    if (v_sq < v_limit) {
      lambdas_p(i) = lambda_plus_no_v;
      lambdas_m(i) = lambda_minus_no_v;
    } else if (bbx_sq < b_limit) {
      lambdas_p(i) = lambda_plus_no_bbx;
      lambdas_m(i) = lambda_minus_no_bbx;
    } else {
      lambdas_p(i) = lambda_plus;
      lambdas_m(i) = lambda_minus;
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
    Real b_sq, Real g00, Real g01, Real g11, Real *plambda_plus, Real *plambda_minus) {
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
  Real& w_d  = prim(IDN,k,j,i);
  Real& w_p  = prim(IPR,k,j,i);
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
