//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file adiabatic_mhd_sr.cpp
//  \brief Implements functions for going between primitive and conserved variables in
//  special-relativistic MHD, as well as for computing wavespeeds.

// C++ headers
#include <algorithm>  // max(), min()
#include <cfloat>     // FLT_MIN
#include <cmath>      // abs(), cbrt(), isfinite(), NAN, sqrt()

// Athena++ headers
#include "eos.hpp"
#include "../athena.hpp"                   // enums, macros
#include "../athena_arrays.hpp"            // AthenaArray
#include "../parameter_input.hpp"          // ParameterInput
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../field/field.hpp"              // FaceField
#include "../mesh/mesh.hpp"                // MeshBlock

// Declarations
static Real EResidual(Real w_guess, Real dd, Real ee, Real m_sq, Real bb_sq, Real ss_sq,
    Real gamma_prime);
static Real EResidualPrime(Real w_guess, Real dd, Real m_sq, Real bb_sq, Real ss_sq,
    Real gamma_prime);

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
  sigma_max_ = pin->GetOrAddReal("hydro", "sigma_max",  0.0);
  beta_min_ = pin->GetOrAddReal("hydro", "beta_min", 0.0);
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
    const AthenaArray<Real> &prim_old, const FaceField &bb, AthenaArray<Real> &prim,
    AthenaArray<Real> &bb_cc, Coordinates *pco, int il, int iu, int jl, int ju, int kl,
    int ku) {
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
        Real dd = cons(IDN,k,j,i);
        Real ee = cons(IEN,k,j,i);

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
        Real w_init = (s2 >= 0.0 and a1 >= 0.0) ? -2.0*a0/(a1+s) : 0.5*(-a1+s);

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

        cons(IDN,k,j,i) = dd;
        cons(IEN,k,j,i) = ee;

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
     const AthenaArray<Real> &bc, AthenaArray<Real> &cons, Coordinates *pco, int il,
     int iu, int jl, int ju, int kl, int ku) {
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
//   almost same function as in adiabatic_mhd_gr.cpp

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
#pragma omp simd simdlen(SIMD_WIDTH)
  for (int i=il; i<=iu; ++i) {

    // Extract primitives
    Real rho = prim(IDN,k,j,i);
    Real pgas = prim(IPR,k,j,i);
    Real vx = prim(ivx,k,j,i);
    Real vy = prim(ivy,k,j,i);
    Real vz = prim(ivz,k,j,i);
    Real bbx = bbx_vals(i);
    Real bby = prim(IBY,k,j,i);
    Real bbz = prim(IBZ,k,j,i);

    // Calculate velocity
    Real v_sq = SQR(vx) + SQR(vy) + SQR(vz);
    Real u0 = std::sqrt(1.0 / (1.0 - v_sq));

    // Calculate contravariant magnetic field
    Real b[4];
    b[0] = (bbx*vx + bby*vy + bbz*vz) * u0;
    b[1] = bbx / u0 + b[0] * vx;
    b[2] = bby / u0 + b[0] * vy;
    b[3] = bbz / u0 + b[0] * vz;

    // Calculate intermediate quantities
    Real gamma_rel_sq = u0 * u0;
    Real w_gas = rho + gamma_adi_red * pgas;
    Real cs_sq = gamma_adi * pgas / w_gas;                       // (MB2005 4)
    Real b_sq = -SQR(b[0]) + SQR(b[1]) + SQR(b[2]) + SQR(b[3]);
    Real bbx_sq = SQR(bbx);
    Real vx_sq = SQR(vx);

    // Calculate wavespeeds in vanishing velocity case (MB2006 57)
    Real lambda_plus_no_v, lambda_minus_no_v;
    {
      Real w_tot_inv = 1.0 / (w_gas + b_sq);
      Real a1 = -(b_sq + cs_sq * (w_gas + bbx_sq)) * w_tot_inv;
      Real a0 = cs_sq * bbx_sq * w_tot_inv;
      Real s2 = SQR(a1) - 4.0*a0;
      Real s = (s2 < 0.0) ? 0.0 : std::sqrt(s2);
      Real lambda_sq = 0.5 * (-a1 + s);
      lambda_plus_no_v = std::sqrt(lambda_sq);
      lambda_minus_no_v = -lambda_plus_no_v;
    }

    // Calculate wavespeeds in vanishing normal field case (MB2006 58)
    Real lambda_plus_no_bbx, lambda_minus_no_bbx;
    {
      Real v_dot_bb_perp = vy*bby + vz*bbz;
      Real q = b_sq - cs_sq*SQR(v_dot_bb_perp);
      Real denominator_inv = 1.0 / (w_gas * (cs_sq + gamma_rel_sq*(1.0-cs_sq)) + q);
      Real a1 = -2.0 * w_gas * gamma_rel_sq * vx * (1.0-cs_sq) * denominator_inv;
      Real a0 = (w_gas * (-cs_sq + gamma_rel_sq*vx_sq*(1.0-cs_sq)) - q) * denominator_inv;
      Real s2 = SQR(a1) - 4.0*a0;
      Real s = (s2 < 0.0) ? 0.0 : std::sqrt(s2);
      lambda_plus_no_bbx = (s2 >= 0.0 and a1 >= 0.0) ? -2.0*a0/(a1+s) : 0.5*(-a1+s);
      lambda_minus_no_bbx = (s2 >= 0.0 and a1 < 0.0) ? -2.0*a0/(a1-s) : 0.5*(-a1-s);
    }

    // Calculate wavespeeds in general case (MB2006 56)
    Real lambda_plus, lambda_minus;
    {
      // Calculate quartic coefficients
      Real bt_sq = SQR(b[0]);
      Real bx_sq = SQR(b[1]);
      Real tmp1 = SQR(gamma_rel_sq) * w_gas * (1.0-cs_sq);
      Real tmp2 = gamma_rel_sq * (b_sq + w_gas * cs_sq);
      Real denominator_inv = 1.0 / (tmp1 + tmp2 - cs_sq * bt_sq);
      Real a3 = (-(4.0*tmp1+2.0*tmp2)*vx + 2.0*cs_sq*b[0]*b[1]) * denominator_inv;
      Real a2 = (6.0*tmp1*vx_sq + tmp2*(vx_sq-1.0) + cs_sq*(bt_sq-bx_sq))
        * denominator_inv;
      Real a1 = (-4.0*tmp1*vx_sq*vx + 2.0*tmp2*vx - 2.0*cs_sq*b[0]*b[1])
        * denominator_inv;
      Real a0 = (tmp1*SQR(vx_sq) - tmp2*vx_sq + cs_sq*bx_sq) * denominator_inv;

      // Calculate reduced quartic coefficients
      Real b2 = a2 - 0.375*SQR(a3);
      Real b1 = a1 - 0.5*a2*a3 + 0.125*a3*SQR(a3);
      Real b0 = a0 - 0.25*a1*a3 + 0.0625*a2*SQR(a3) - 3.0/256.0*SQR(SQR(a3));

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
        s2 = 0.25*SQR(z0) - b0;
        Real s = (s2 < 0.0) ? 0.0 : std::sqrt(s2);
        Real d0 = (b1 < 0) ? 0.5*z0+s : 0.5*z0-s;
        Real e0 = (b1 < 0) ? 0.5*z0-s : 0.5*z0+s;

        // Solve quadratic equations
        s2 = SQR(d1) - 4.0*d0;
        s = (s2 < 0.0) ? 0.0 : std::sqrt(s2);
        y1 = (s2 >= 0.0 and d1 < 0.0) ? -2.0*d0/(d1-s) : 0.5*(-d1-s);
        y2 = (s2 >= 0.0 and d1 >= 0.0) ? -2.0*d0/(d1+s) : 0.5*(-d1+s);
        s2 = SQR(e1) - 4.0*e0;
        s = (s2 < 0.0) ? 0.0 : std::sqrt(s2);
        y3 = (s2 >= 0.0 and e1 < 0.0) ? -2.0*e0/(e1-s) : 0.5*(-e1-s);
        y4 = (s2 >= 0.0 and e1 >= 0.0) ? -2.0*e0/(e1+s) : 0.5*(-e1+s);
      }

      // Calculate extremal original quartic roots
      lambda_minus = std::min(y1, y3) - 0.25*a3;
      lambda_plus = std::max(y2, y4) - 0.25*a3;

      // Ensure wavespeeds are not superluminal
      lambda_minus = std::max(lambda_minus, -1.0);
      lambda_plus = std::min(lambda_plus, 1.0);
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

static Real EResidual(Real w_guess, Real dd, Real ee, Real m_sq, Real bb_sq, Real ss_sq,
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

static Real EResidualPrime(Real w_guess, Real dd, Real m_sq, Real bb_sq, Real ss_sq,
    Real gamma_prime) {
  Real v_sq = (m_sq + ss_sq/SQR(w_guess) * (2.0*w_guess + bb_sq))
      / SQR(w_guess + bb_sq);                                         // (cf. MM A3)
  Real gamma_sq = 1.0/(1.0-v_sq);
  Real gamma_lorentz = std::sqrt(gamma_sq);
  Real chi = (1.0 - v_sq) * (w_guess - gamma_lorentz * dd);           // (cf. MM A11)
  Real w_cu = SQR(w_guess) * w_guess;
  Real w_b_cu = SQR(w_guess + bb_sq) * (w_guess + bb_sq);
  Real dv_sq_dw = -2.0 / (w_cu*w_b_cu) * (ss_sq
      * (3.0*w_guess*(w_guess+bb_sq) + SQR(bb_sq)) + m_sq*w_cu);      // (MM A16)
  Real dchi_dw = 1.0 - v_sq
      - 0.5*gamma_lorentz * (dd + 2.0*gamma_lorentz*chi) * dv_sq_dw;  // (cf. MM A14)
  Real drho_dw = -0.5*gamma_lorentz*dd*dv_sq_dw;                      // (MM A15)
  Real dpgas_dchi = 1.0/gamma_prime;                                  // (MM A18)
  Real dpgas_drho = 0.0;                                              // (MM A18)
  Real dpgas_dw = dpgas_dchi * dchi_dw + dpgas_drho * drho_dw;
  return 1.0 - dpgas_dw + 0.5*bb_sq * dv_sq_dw + ss_sq/w_cu;
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

  // apply density floor
  w_d = (w_d > density_floor_) ?  w_d : density_floor_;
  // apply pressure floor
  w_p = (w_p > pressure_floor_) ?  w_p : pressure_floor_;

  return;
}
