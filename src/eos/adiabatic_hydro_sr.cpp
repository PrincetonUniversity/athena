//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file adiabatic_hydro_sr.cpp
//  \brief Implements functions for going between primitive and conserved variables in
//  special-relativistic hydrodynamics, as well as for computing wavespeeds.

// C++ headers
#include <cmath>   // atan2(), cbrt(), cos(), sqrt()
#include <cfloat>  // FLT_MIN

// Athena++ headers
#include "eos.hpp"
#include "../athena.hpp"                   // enums, macros
#include "../athena_arrays.hpp"            // AthenaArray
#include "../parameter_input.hpp"          // ParameterInput
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../field/field.hpp"              // FaceField
#include "../mesh/mesh.hpp"                // MeshBlock


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
  gamma_max_ = pin->GetOrAddReal("hydro", "gamma_max", 1000.0);
}

//----------------------------------------------------------------------------------------
// Destructor

EquationOfState::~EquationOfState() {}

//----------------------------------------------------------------------------------------
// Variable inverter
// Inputs:
//   cons: conserved quantities
//   prim_old: primitive quantities from previous half timestep (not used)
//   bb: face-centered magnetic field (not used)
//   pco: pointer to Coordinates
//   il,iu,jl,ju,kl,ku: index bounds of region to be updated
// Outputs:
//   prim: primitives
//   bb_cc: cell-centered magnetic field (not used)
// Notes:
//   solves quartic equation for velocity explicitly
//   follows relativistic hydro routine in old Athena's convert_var.c:
//     1) given quartic |v|^4 + a3 |v|^3 + a2 |v|^2 + a1 |v| + a0
//     2) construct resolvent cubic x^3 + b2 x^2 + b1 x + b0:
//          b2 = -a2
//          b1 = -4 a0 + a1 a3
//          b0 = 4 a0 a2 - a1^2 - a0 a3^2
//     3) eliminate quadratic term from cubic y^3 + (3 c1) y - 2 c2 = 0:
//          c1 = b1/3 - b2^2/9
//          c2 = -(1/2) b0 + (1/6) b1 b2 - (1/27) b2^3
//          c3 = c1^3 + c2^2
//     4) find real root of new cubic:
//          if c3 >= 0, y0 = (c2 + sqrt(c3))^(1/3) + (c2 - sqrt(c3))^(1/3)
//          otherwise use y0 = 2 (c2^2 + c3)^(1/6) cos((1/3) atan2(sqrt(-c2), c3))
//          formulas are equivalent except for implicit assumptions about branch cuts
//     5) find real root of original (resolvent) cubic:
//          x0 = y0 - b2/3
//     6) solve for (correct) root of original quartic:
//          d1 = 1/2 * (a3 + sqrt(4 x0 - 4 a2 + a3^2))
//          d0 = 1/2 * (x0 - sqrt(x0^2 - 4 a0))
//          then |v|^2 + d1 |v| + d0 = 0
//          |v| = 1/2 * (-d1 + sqrt(d1^2 - 4 d0))

void EquationOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
  const AthenaArray<Real> &prim_old, const FaceField &bb, AthenaArray<Real> &prim,
  AthenaArray<Real> &bb_cc, Coordinates *pco, int il, int iu, int jl, int ju, int kl,
  int ku) {
  // Parameters
  const Real max_velocity = std::sqrt(1.0 - 1.0/SQR(gamma_max_));

  // Extract ratio of specific heats
  const Real gamma_adi = GetGamma();
  const Real gamma_adi_minus_1 = gamma_adi - 1.0;

  const Real TINY_NUMBER_p2 = SQR(TINY_NUMBER);

  // Make array copies for performance reasons
  AthenaArray<Real> cons_copy,prim_copy;
  cons_copy.InitWithShallowCopy(cons);
  prim_copy.InitWithShallowCopy(prim);

  // Go through cells
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {

        // Extract conserved quantities
        Real d = cons_copy(IDN,k,j,i);
        Real e = cons_copy(IEN,k,j,i);

        Real mx = cons_copy(IVX,k,j,i);
        Real my = cons_copy(IVY,k,j,i);
        Real mz = cons_copy(IVZ,k,j,i);

        // Extract primitives quantities
        Real rho = prim_copy(IDN,k,j,i);
        Real pgas = prim_copy(IPR,k,j,i);
        Real vx = prim_copy(IVX,k,j,i);
        Real vy = prim_copy(IVY,k,j,i);
        Real vz = prim_copy(IVZ,k,j,i);

        // Calculate total momentum
        Real m_sq = SQR(mx) + SQR(my) + SQR(mz);

        // Step 1: Prepare quartic coefficients
        Real m_abs = std::sqrt(m_sq);
        Real d_p2 = SQR(d);
        Real gamma_adi_minus_1_p2 = SQR(gamma_adi_minus_1);
        Real denom_inv = 1.0 / (gamma_adi_minus_1_p2 * (d_p2 + m_sq));
        Real a3 = -2.0 * gamma_adi * gamma_adi_minus_1 * m_abs * e * denom_inv;
        Real a2 = (SQR(gamma_adi) * SQR(e) + 2.0 * gamma_adi_minus_1 * m_sq -
                   gamma_adi_minus_1_p2 * d_p2) * denom_inv;
        Real a1 = -2.0 * gamma_adi * m_abs * e * denom_inv;
        Real a0 = m_sq * denom_inv;

        // Step 2: Find resolvent cubic coefficients
        Real b2 = -a2;
        Real b1 = -4.0*a0 + a1*a3;
        Real b0 = 4.0*a0*a2 - SQR(a1) - a0*SQR(a3);

        // Step 3: Eliminate quadratic term from cubic
        Real b2_p2 = SQR(b2);
        Real c1 = b1/3.0 - b2_p2/9.0;
        Real c2 = -b0/2.0 + b1*b2/6.0 - b2_p2*b2/27.0;
        Real c3 = SQR(c1)*c1 + SQR(c2);

        // Step 4: Find real root of new cubic
        Real y0;

        if (c3 >= 0.0) {
          y0 = std::cbrt(c2 + std::sqrt(c3)) + std::cbrt(c2 - std::sqrt(c3));
        } else {
          y0 = 2.0 * std::cbrt(SQR(c2) + c3)
            * std::cos(std::atan2(std::sqrt(-c3), c2) / 3.0);
        }

        // Step 5: Find real root of original (resolvent) cubic:
        Real x0 = y0 - b2/3.0;

        // Step 6: Solve for (correct) root of original quartic
        Real d1 = 0.5 * (a3 + std::sqrt(4.0*x0 - 4.0*a2 + SQR(a3)));
        Real d0 = 0.5 * (x0 - std::sqrt(SQR(x0) - 4.0*a0));
        Real v_abs = 0.5 * (-d1 + std::sqrt(SQR(d1) - 4.0*d0));

        // Ensure velocity is physical
        v_abs = (v_abs > 0.0) ? v_abs : 0.0;                    // sets NaN to 0
        v_abs = (v_abs < max_velocity) ? v_abs : max_velocity;

        // Set density, correcting only conserved density if floor applied
        Real gamma_rel = 1.0 / std::sqrt(1.0 - SQR(v_abs));
        Real v_over_m = v_abs / m_abs;
        Real e_tmp = m_sq * v_over_m;

        // Case out based on whether momentum vanishes
        if (m_sq < TINY_NUMBER_p2 ) { // zero velocity case
          gamma_rel = 1.0;
          v_over_m = 0.0;
          e_tmp = 0.0;
        }

        rho = d / gamma_rel;
        if (rho < density_floor_) {
          rho = density_floor_;
          d = gamma_rel * rho;
        }

        // Set velocity
        vx = mx * v_over_m;
        vy = my * v_over_m;
        vz = mz * v_over_m;

        // Set pressure, correcting only energy if floor applied
        pgas = gamma_adi_minus_1 * (e - e_tmp - rho);
        if (pgas < pressure_floor_) {
          pgas = pressure_floor_;
          e = pgas/gamma_adi_minus_1 + e_tmp + rho;
        }

        cons_copy(IDN,k,j,i) = d;
        cons_copy(IEN,k,j,i) = e;
        prim_copy(IDN,k,j,i) = rho;
        prim_copy(IPR,k,j,i) = pgas;
        prim_copy(IVX,k,j,i) = vx;
        prim_copy(IVY,k,j,i) = vy;
        prim_copy(IVZ,k,j,i) = vz;

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
// Notes:
//   single-cell function exists for other purposes; call made to that function rather
//       than having duplicate code

void EquationOfState::PrimitiveToConserved(const AthenaArray<Real> &prim,
     const AthenaArray<Real> &bb_cc, AthenaArray<Real> &cons, Coordinates *pco, int il,
     int iu, int jl, int ju, int kl, int ku) {
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
        Real v1 = prim(IVX,k,j,i);
        Real v2 = prim(IVY,k,j,i);
        Real v3 = prim(IVZ,k,j,i);

        // Calculate 4-velocity
        Real u0 = 1.0 / std::sqrt(1.0 - SQR(v1) - SQR(v2) - SQR(v3));
        Real u1 = u0 * v1;
        Real u2 = u0 * v2;
        Real u3 = u0 * v3;

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
//   same function as in adiabatic_hydro_gr.cpp
//   references Mignone & Bodo 2005, MNRAS 364 126 (MB)

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

//---------------------------------------------------------------------------------------
// \!fn void EquationOfState::ApplyPrimitiveFloors(AthenaArray<Real> &prim,
//           int k, int j, int i)
// \brief Apply density and pressure floors to reconstructed L/R cell interface states

void EquationOfState::ApplyPrimitiveFloors(AthenaArray<Real> &prim, int k, int j, int i) {
  Real& w_d  = prim(IDN,k,j,i);
  Real& w_p  = prim(IPR,k,j,i);

  // apply density floor
  w_d = (w_d > density_floor_) ?  w_d : density_floor_;
  // apply pressure floor
  w_p = (w_p > pressure_floor_) ?  w_p : pressure_floor_;

  return;
}
