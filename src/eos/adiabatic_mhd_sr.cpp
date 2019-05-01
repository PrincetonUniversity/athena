//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file adiabatic_mhd_sr.cpp
//  \brief Implements functions for going between primitive and conserved variables in
//  special-relativistic MHD, as well as for computing wavespeeds.

// C headers

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // abs(), cbrt(), isfinite(), NAN, sqrt()
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
Real EResidual(Real w_guess, Real dd, Real ee, Real m_sq, Real bb_sq, Real ss_sq,
               Real gamma_prime);
Real EResidualPrime(Real w_guess, Real dd, Real m_sq, Real bb_sq, Real ss_sq,
                    Real gamma_prime);
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
  sigma_max_ = pin->GetOrAddReal("hydro", "sigma_max",  0.0);
  beta_min_ = pin->GetOrAddReal("hydro", "beta_min", 0.0);
  gamma_max_ = pin->GetOrAddReal("hydro", "gamma_max", 1000.0);
}

//----------------------------------------------------------------------------------------
// Variable inverter
// Inputs:
//   cons: conserved quantities
//   prim_old: primitive quantities from previous half timestep (not used)
//   bb: face-centered magnetic field
//   pco: pointer to Coordinates
//   il,iu,jl,ju,kl,ku: index bounds of region to be updated
// Outputs:
//   prim: primitives
//   bb_cc: cell-centered magnetic field
// Notes:
//   follows Mignone & McKinney 2007, MNRAS 378 1118 (MM)
//   follows hlld_sr.c in Athena 4.2 in using W and E rather than W' and E'

void EquationOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
                                           const AthenaArray<Real> &prim_old,
                                           const FaceField &bb, AthenaArray<Real> &prim,
                                           AthenaArray<Real> &bb_cc, Coordinates *pco,
                                           int il, int iu, int jl,
                                           int ju, int kl, int ku) {
  // Parameters
  const Real gamma_prime = gamma_/(gamma_-1.0);
  const Real v_sq_max = 1.0 - 1.0/SQR(gamma_max_);
  const Real max_velocity = std::sqrt(v_sq_max);

  // Interpolate magnetic field from faces to cell centers
  pmy_block_->pfield->CalculateCellCenteredField(bb, bb_cc, pco, il, iu, jl, ju, kl, ku);

  // Go through cells
  for (int k=kl; k<=ku; k++) {
    for (int j=jl; j<=ju; j++) {
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        // Extract conserved quantities
        Real &dd = cons(IDN,k,j,i);
        Real &ee = cons(IEN,k,j,i);

        Real mx = cons(IM1,k,j,i);
        Real my = cons(IM2,k,j,i);
        Real mz = cons(IM3,k,j,i);

        // Extract primitives
        Real rho = prim(IDN,k,j,i);
        Real pgas = prim(IPR,k,j,i);
        Real vx = prim(IVX,k,j,i);
        Real vy = prim(IVY,k,j,i);
        Real vz = prim(IVZ,k,j,i);

        // Extract cell-centered magnetic field
        Real bbx = bb_cc(IB1,k,j,i);
        Real bby = bb_cc(IB2,k,j,i);
        Real bbz = bb_cc(IB3,k,j,i);

        // Calculate variations on conserved quantities
        Real m_sq = SQR(mx) + SQR(my) + SQR(mz);
        Real bb_sq = SQR(bbx) + SQR(bby) + SQR(bbz);
        Real m_dot_bb = mx*bbx + my*bby + mz*bbz;
        Real ss_sq = SQR(m_dot_bb);

        // Construct initial guess for enthalpy W (cf. MM A26-A27)
        Real a1 = 4.0/3.0 * (bb_sq - ee);
        Real a0 = ONE_3RD * (m_sq + bb_sq * (bb_sq - 2.0*ee));
        Real s2 = SQR(a1) - 4.0*a0;
        Real s = (s2 < 0.0) ? 0.0 : std::sqrt(s2);
        Real w_init = (s2 >= 0.0 && a1 >= 0.0) ? -2.0*a0/(a1+s) : 0.5*(-a1+s);

        // Apply Newton-Raphson method to find new W
        const int num_iterations = 5;
        Real w_new = w_init;
        Real res_new = EResidual(w_new, dd, ee, m_sq, bb_sq, ss_sq, gamma_prime);
        for (int n = 0; n < num_iterations; ++n) {
          Real w_old = w_new;
          Real res_old = res_new;
          Real derivative = EResidualPrime(w_old, dd, m_sq, bb_sq, ss_sq, gamma_prime);
          Real delta = -res_old / derivative;
          w_new = w_old + delta;
          res_new = EResidual(w_new, dd, ee, m_sq, bb_sq, ss_sq, gamma_prime);
        }
        Real w_true = w_new;

        // Set velocity
        vx = (mx + m_dot_bb/w_true * bbx) / (w_true + bb_sq);  // (MM A10)
        vy = (my + m_dot_bb/w_true * bby) / (w_true + bb_sq);  // (MM A10)
        vz = (mz + m_dot_bb/w_true * bbz) / (w_true + bb_sq);  // (MM A10)
        Real v_sq = SQR(vx) + SQR(vy) + SQR(vz);
        Real vel_ratio = max_velocity / std::sqrt(v_sq);
        if (v_sq > v_sq_max) {
          vx *= vel_ratio;
          vy *= vel_ratio;
          vz *= vel_ratio;
          v_sq = v_sq_max;
        }
        Real gamma_sq = 1.0/(1.0-v_sq);
        Real gamma_lorentz = std::sqrt(gamma_sq);

        // Calculate magnetic pressure
        Real v_dot_bb = vx*bbx + vy*bby + vz*bbz;
        Real b_sq = bb_sq * (1.0-v_sq) + SQR(v_dot_bb);
        Real pmag = 0.5*b_sq;

        // Calculate floors for density and pressure
        Real density_floor_local = density_floor_;
        if (sigma_max_ > 0.0) {
          density_floor_local = std::max(density_floor_local, 2.0*pmag/sigma_max_);
        }
        Real pressure_floor_local = pressure_floor_;
        if (beta_min_ > 0.0) {
          pressure_floor_local = std::max(pressure_floor_local, beta_min_*pmag);
        }

        // Set density, correcting only conserved density if floor applied
        rho = dd/gamma_lorentz;  // (MM A12)
        if (rho < density_floor_local) {
          rho = density_floor_local;
          dd = gamma_lorentz * rho;
        }

        // Set pressure, correcting only energy if floor applied
        Real chi = (1.0 - v_sq) * (w_true - gamma_lorentz * dd);  // (cf. MM A11)
        pgas = chi/gamma_prime;  // (MM A17
        Real bt = gamma_lorentz * v_dot_bb;
        if (pgas < pressure_floor_local) {
          pgas = pressure_floor_local;
          Real w = rho + gamma_prime * pgas + b_sq;
          ee = gamma_sq * w - SQR(bt) - (pgas + pmag);
        }

        // cons(IDN,k,j,i) = dd;
        // cons(IEN,k,j,i) = ee;

        prim(IDN,k,j,i) = rho;
        prim(IPR,k,j,i) = pgas;
        prim(IVX,k,j,i) = vx;
        prim(IVY,k,j,i) = vy;
        prim(IVZ,k,j,i) = vz;
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

void EquationOfState::PrimitiveToConserved(const AthenaArray<Real> &prim,
                                           const AthenaArray<Real> &bc,
                                           AthenaArray<Real> &cons, Coordinates *pco,
                                           int il, int iu, int jl,
                                           int ju, int kl, int ku) {
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
        Real v1 = prim(IVX,k,j,i);
        Real v2 = prim(IVY,k,j,i);
        Real v3 = prim(IVZ,k,j,i);

        Real bb1 = bc(IB1,k,j,i);
        Real bb2 = bc(IB2,k,j,i);
        Real bb3 = bc(IB3,k,j,i);

        // Calculate velocity
        Real u0 = 1.0 / std::sqrt(1.0 - SQR(v1) - SQR(v2) - SQR(v3));

        // Calculate 4-magnetic field
        Real b0 = (bb1*v1 + bb2*v2 + bb3*v3) * u0;
        Real b1 = bb1 / u0 + b0 * v1;
        Real b2 = bb2 / u0 + b0 * v2;
        Real b3 = bb3 / u0 + b0 * v3;
        Real b_sq = -SQR(b0) + SQR(b1) + SQR(b2) + SQR(b3);

        // Set conserved quantities
        Real wtot_u02 = (rho + gamma_adi_red * pgas + b_sq) * u0 * u0;
        Real dd = rho * u0;
        Real ee = wtot_u02 - b0 * b0 - (pgas + 0.5*b_sq);
        Real m1 = wtot_u02 * v1 - b0 * b1;
        Real m2 = wtot_u02 * v2 - b0 * b2;
        Real m3 = wtot_u02 * v3 - b0 * b3;

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
//   applies approximation that may slightly overestimate wavespeeds
//   almost same function as in adiabatic_mhd_gr.cpp

void EquationOfState::FastMagnetosonicSpeedsSR(const AthenaArray<Real> &prim,
                                               const AthenaArray<Real> &bbx_vals,
                                               int k, int j, int il, int iu, int ivx,
                                               AthenaArray<Real> &lambdas_p,
                                               AthenaArray<Real> &lambdas_m) {
  // Calculate cyclic permutations of indices
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  // Calculate ratio of specific heats
  const Real gamma_adi = gamma_;
  const Real gamma_adi_red = gamma_adi/(gamma_adi-1.0);

  // Go through states
#pragma omp simd simdlen(SIMD_WIDTH)
  for (int i=il; i<=iu; ++i) {
    // Extract primitives
    Real rho = prim(IDN,i);
    Real pgas = prim(IPR,i);
    Real v1 = prim(ivx,i);
    Real v2 = prim(ivy,i);
    Real v3 = prim(ivz,i);
    Real bb1 = bbx_vals(i);
    Real bb2 = prim(IBY,i);
    Real bb3 = prim(IBZ,i);

    // Calculate Lorentz factor and (twice) magnetic pressure
    Real v_sq = SQR(v1) + SQR(v2) + SQR(v3);
    Real u0 = std::sqrt(1.0 / (1.0 - v_sq));
    Real b0 = (bb1 * v1 + bb2 * v2 + bb3 * v3) * u0;
    Real b1 = bb1 / u0 + b0 * v1;
    Real b2 = bb2 / u0 + b0 * v2;
    Real b3 = bb3 / u0 + b0 * v3;
    Real b_sq = -SQR(b0) + SQR(b1) + SQR(b2) + SQR(b3);

    // Calculate comoving fast magnetosonic speed
    Real w_gas = rho + gamma_adi_red * pgas;
    Real cs_sq = gamma_adi * pgas / w_gas;
    Real va_sq = b_sq / (w_gas + b_sq);
    Real cms_sq = cs_sq + va_sq - cs_sq * va_sq;

    // Set fast magnetosonic speeds
    Real a = SQR(u0) - (SQR(u0) - 1.0) * cms_sq;
    Real b = -2.0 * SQR(u0) * v1 * (1.0 - cms_sq);
    Real c = SQR(u0 * v1) - (SQR(u0 * v1) + 1.0) * cms_sq;
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

namespace {
//----------------------------------------------------------------------------------------
// Function whose value vanishes for correct enthalpy
// Inputs:
//   w_guess: guess for total enthalpy W
//   dd: relativistic density D
//   ee: total energy E
//   m_sq: square magnitude of momentum \vec{m}
//   bb_sq: square magnitude of magnetic field \vec{B}
//   ss_sq: (\vec{m} \cdot \vec{B})^2
//   gamma_prime: reduced adiabatic gas constant Gamma' = Gamma/(Gamma-1)
// Outputs:
//   returned value: calculated minus given value of E
// Notes:
//   follows Mignone & McKinney 2007, MNRAS 378 1118 (MM)
//   implementation follows that of hlld_sr.c in Athena 4.2
//   same function as in hlld_rel.cpp

Real EResidual(Real w_guess, Real dd, Real ee, Real m_sq, Real bb_sq, Real ss_sq,
               Real gamma_prime) {
  Real v_sq = (m_sq + ss_sq/SQR(w_guess) * (2.0*w_guess + bb_sq))
              / SQR(w_guess + bb_sq);                                      // (cf. MM A3)
  Real gamma_sq = 1.0/(1.0-v_sq);
  Real gamma_lorentz = std::sqrt(gamma_sq);
  Real chi = (1.0 - v_sq) * (w_guess - gamma_lorentz * dd);        // (cf. MM A11)
  Real pgas = chi/gamma_prime;                                     // (MM A17)
  Real ee_calc = w_guess - pgas + 0.5*bb_sq * (1.0+v_sq)
                 - ss_sq*0.5/SQR(w_guess);                                    // (MM A1)
  return ee_calc - ee;
}

//----------------------------------------------------------------------------------------
// Derivative of EResidual()
// Inputs:
//   w_guess: guess for total enthalpy W
//   dd: relativistic density D
//   m_sq: square magnitude of momentum \vec{m}
//   bb_sq: square magnitude of magnetic field \vec{B}
//   ss_sq: (\vec{m} \cdot \vec{B})^2
//   gamma_prime: reduced adiabatic gas constant Gamma' = Gamma/(Gamma-1)
// Outputs:
//   returned value: derivative of calculated value of E
// Notes:
//   follows Mignone & McKinney 2007, MNRAS 378 1118 (MM)
//   implementation follows that of hlld_sr.c in Athena 4.2
//   same function as in hlld_mhd_rel.cpp

Real EResidualPrime(Real w_guess, Real dd, Real m_sq, Real bb_sq, Real ss_sq,
                    Real gamma_prime) {
  Real v_sq = (m_sq + ss_sq/SQR(w_guess) * (2.0*w_guess + bb_sq))
              / SQR(w_guess + bb_sq);                                 // (cf. MM A3)
  Real gamma_sq = 1.0/(1.0-v_sq);
  Real gamma_lorentz = std::sqrt(gamma_sq);
  Real chi = (1.0 - v_sq) * (w_guess - gamma_lorentz * dd);           // (cf. MM A11)
  Real w_cu = SQR(w_guess) * w_guess;
  Real w_b_cu = SQR(w_guess + bb_sq) * (w_guess + bb_sq);
  Real dv_sq_dw = -2.0 / (w_cu*w_b_cu)                                // (MM A16)
                  * (ss_sq * (3.0*w_guess*(w_guess+bb_sq) + SQR(bb_sq)) + m_sq*w_cu);
  Real dchi_dw = 1.0 - v_sq                                           // (cf. MM A14)
                 - 0.5*gamma_lorentz * (dd + 2.0*gamma_lorentz*chi) * dv_sq_dw;
  Real drho_dw = -0.5*gamma_lorentz*dd*dv_sq_dw;                      // (MM A15)
  Real dpgas_dchi = 1.0/gamma_prime;                                  // (MM A18)
  Real dpgas_drho = 0.0;                                              // (MM A18)
  Real dpgas_dw = dpgas_dchi * dchi_dw + dpgas_drho * drho_dw;
  return 1.0 - dpgas_dw + 0.5*bb_sq * dv_sq_dw + ss_sq/w_cu;
}
} // namespace

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

  // apply density floor
  w_d = (w_d > density_floor_) ?  w_d : density_floor_;
  // apply pressure floor
  w_p = (w_p > pressure_floor_) ?  w_p : pressure_floor_;

  return;
}
