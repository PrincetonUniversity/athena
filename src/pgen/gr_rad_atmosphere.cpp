//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gr_rad_beam.cpp
//  \brief Problem generator for GR radiation isothermal Schwarzschild atmosphere

// C++ headers
#include <cmath>      // exp, pow, sqrt
#include <cstdlib>    // exit (needed for defs.hpp)
#include <iostream>   // cout (needed for defs.hpp), endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error (needed for defs.hpp)
#include <string>     // c_str, strcmp, string (some needed for defs.hpp)

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "../athena.hpp"                   // Real, enums
#include "../athena_arrays.hpp"            // AthenaArray
#include "../parameter_input.hpp"          // ParameterInput
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../eos/eos.hpp"                  // EquationOfState
#include "../hydro/hydro.hpp"              // Hydro
#include "../radiation/radiation.hpp"      // Radiation

// Configuration checking
#if not RADIATION_ENABLED
#error "This problem generator must be used with radiation"
#endif
#if not GENERAL_RELATIVITY
#error "This problem generator must be used with general relativity"
#endif
#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator does not support magnetic fields"
#endif

// Declarations
void FixedLeft(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, AthenaArray<Real> &prim_rad, Real time, Real dt, int il, int iu,
    int jl, int ju, int kl, int ku, int ngh);
void FixedRight(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, AthenaArray<Real> &prim_rad, Real time, Real dt, int il, int iu,
    int jl, int ju, int kl, int ku, int ngh);

// Global variables
namespace {
Real kb_cgs, mp_cgs, arad_cgs;  // physical constants
Real rho_unit, p_unit;          // conversion factors
Real mu;                        // molecular weight
Real gamma_adi;                 // adiabatic index
Real alpha_in;                  // lapse at inner edge
Real tt_inf_cgs;                // temperature at infinity in K
Real kappa;                     // opacity in code units
Real ptot_in_cgs;               // total pressure at inner edge in dyne/cm^2
Real var_a, var_b;              // solution parameters
}  // namespace

//----------------------------------------------------------------------------------------
// Function for preparing Mesh
// Inputs:
//   pin: input parameters
// Outputs: (none)

void Mesh::InitUserMeshData(ParameterInput *pin) {

  // Define constants
  Real c_cgs = 2.99792458e10;
  kb_cgs = 1.380649e-16;
  mp_cgs = 1.67262192369e-24;
  Real gg_msun_cgs = 1.32712440018e26;
  arad_cgs = 4.0 * 5.670374419e-5 / c_cgs;

  // Check coordinates
  if (std::strcmp(COORDINATE_SYSTEM, "schwarzschild") != 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Problem Generator\n"
        << "gr_rad_atmosphere only supports Schwarzschild coordinates" << std::endl;
    ATHENA_ERROR(msg);
    return;
  }

  // Read parameters from profile file
  Real r_in = pin->GetReal("mesh", "x1min");
  Real r_out = pin->GetReal("mesh", "x1max");
  rho_unit = pin->GetReal("radiation", "density_cgs");
  mu = pin->GetReal("radiation", "mol_weight");
  gamma_adi = pin->GetReal("hydro", "gamma");
  Real m_msun = pin->GetReal("problem", "m_msun");
  Real h_r = pin->GetReal("problem", "h_r");
  Real beta_r = pin->GetReal("problem", "beta_r");
  Real tau = pin->GetReal("problem", "tau");

  // Calculate solution parameters
  p_unit = rho_unit * SQR(c_cgs);
  Real rg_cgs = gg_msun_cgs * m_msun / SQR(c_cgs);
  Real r_in_cgs = r_in * rg_cgs;
  Real r_out_cgs = r_out * rg_cgs;
  alpha_in = std::sqrt(1.0 - 2.0 / r_in);
  tt_inf_cgs = alpha_in * mu * mp_cgs / kb_cgs * SQR(c_cgs) / r_in * h_r;
  Real tt_in_cgs = tt_inf_cgs / alpha_in;
  Real rho_in_cgs =
      mu * mp_cgs / kb_cgs * arad_cgs / (3.0 * beta_r) * tt_in_cgs * SQR(tt_in_cgs);
  Real rho_in = rho_in_cgs / rho_unit;
  kappa = tau / (rho_in * (r_out - r_in));
  Real pgas_in_cgs = kb_cgs / (mu * mp_cgs) * tt_in_cgs * rho_in_cgs;
  Real prad_in_cgs = 1.0 / 3.0 * arad_cgs * SQR(SQR(tt_in_cgs));
  ptot_in_cgs = pgas_in_cgs + prad_in_cgs;
  var_a = SQR(c_cgs) * mu * mp_cgs / (kb_cgs * tt_inf_cgs);
  var_b = 1.0 / 3.0 * arad_cgs * SQR(SQR(tt_inf_cgs));

  // Enroll boundary functions
  if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, FixedLeft);
  }
  if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, FixedRight);
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for setting initial conditions
// Inputs:
//   pin: parameters (unused)
// Outputs: (none)
// Notes:
//   initializes constant state with no radiation

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  // Calculate index bounds
  int il = is - NGHOST;
  int iu = ie + NGHOST;
  int jl = js - (ncells2 > 1 ? NGHOST : 0);
  int ju = je + (ncells2 > 1 ? NGHOST : 0);
  int kl = ks - (ncells3 > 1 ? NGHOST : 0);
  int ku = ke + (ncells3 > 1 ? NGHOST : 0);

  // Prepare scratch arrays
  AthenaArray<Real> gcov, gcon;
  gcov.NewAthenaArray(NMETRIC, iu + 1);
  gcon.NewAthenaArray(NMETRIC, iu + 1);

  // Go through cells
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      pcoord->CellMetric(k, j, il, iu, gcov, gcon);
      for (int i = il; i <= iu; ++i) {

        // Calculate exact solution
        Real r = pcoord->x1v(i);
        Real alpha = std::sqrt(1.0 - 2.0 / r);
        Real ptot_cgs = var_b / SQR(SQR(alpha)) + std::exp(var_a * (alpha_in - alpha))
            * std::pow(alpha_in / alpha, gamma_adi / (gamma_adi - 1.0))
            * (ptot_in_cgs - var_b / SQR(SQR(alpha_in)));
        Real tt_cgs = tt_inf_cgs / alpha;
        Real prad_cgs = 1.0 / 3.0 * arad_cgs * SQR(SQR(tt_cgs));
        Real pgas_cgs = ptot_cgs - prad_cgs;
        Real rho_cgs = mu * mp_cgs / kb_cgs * pgas_cgs / tt_cgs;
        Real rho = rho_cgs / rho_unit;
        Real pgas = pgas_cgs / p_unit;

        // Set hydrodynamical primitives
        phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) = rho;
        phydro->w(IPR,k,j,i) = phydro->w1(IPR,k,j,i) = pgas;
        phydro->w(IVX,k,j,i) = phydro->w1(IVX,k,j,i) = 0.0;
        phydro->w(IVY,k,j,i) = phydro->w1(IVY,k,j,i) = 0.0;
        phydro->w(IVZ,k,j,i) = phydro->w1(IVZ,k,j,i) = 0.0;

        // Initialize radiation
        Real e_rad = prad->arad * SQR(SQR(pgas / rho));
        prad->CalculateRadiationInCellLinear(e_rad, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, k, j, i,
            gcov, prad->cons);
      }
    }
  }

  // Calculate conserved values
  AthenaArray<Real> bb;
  peos->PrimitiveToConserved(phydro->w, bb, phydro->u, pcoord, il, iu, jl, ju, kl, ku);

  // Initialize opacity (never changed)
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      for (int i = il; i <= iu; ++i) {
        prad->opacity(OPAS,k,j,i) = 0.0;
        prad->opacity(OPAA,k,j,i) = kappa;
        prad->opacity(OPAP,k,j,i) = 0.0;
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Fixed boundary condition for left side
// Inputs:
//   pmb: pointer to MeshBlock
//   pcoord: pointer to Coordinates
//   il, iu, jl, ju, kl, ku: indices demarkating active region
//   time, dt: current time and time step (unused)
//   ngh: effective number of ghost zones to populate
// Outputs:
//   prim: primitives set in ghost zones
//   bb: face-centered magnetic field set in ghost zones (unused)
//   prim_rad: radiation primitives set in ghost zones

void FixedLeft(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, AthenaArray<Real> &prim_rad, Real time, Real dt, int il, int iu,
    int jl, int ju, int kl, int ku, int ngh) {

  // Allocate scratch arrays
  AthenaArray<Real> gcov, gcon;
  gcov.NewAthenaArray(NMETRIC, il);
  gcon.NewAthenaArray(NMETRIC, il);

  // Go through cells
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      pcoord->CellMetric(k, j, il - ngh, il - 1, gcov, gcon);
      for (int i = il - ngh; i <= il - 1; ++i) {

        // Calculate exact solution
        Real r = pcoord->x1v(i);
        Real alpha = std::sqrt(1.0 - 2.0 / r);
        Real ptot_cgs = var_b / SQR(SQR(alpha)) + std::exp(var_a * (alpha_in - alpha))
            * std::pow(alpha_in / alpha, gamma_adi / (gamma_adi - 1.0))
            * (ptot_in_cgs - var_b / SQR(SQR(alpha_in)));
        Real tt_cgs = tt_inf_cgs / alpha;
        Real prad_cgs = 1.0 / 3.0 * arad_cgs * SQR(SQR(tt_cgs));
        Real pgas_cgs = ptot_cgs - prad_cgs;
        Real rho_cgs = mu * mp_cgs / kb_cgs * pgas_cgs / tt_cgs;
        Real rho = rho_cgs / rho_unit;
        Real pgas = pgas_cgs / p_unit;

        // Set fluid
        prim(IDN,k,j,i) = rho;
        prim(IPR,k,j,i) = pgas;
        prim(IVX,k,j,i) = 0.0;
        prim(IVY,k,j,i) = 0.0;
        prim(IVZ,k,j,i) = 0.0;

        // Set radiation
        Real e_rad = pmb->prad->arad * SQR(SQR(pgas / rho));
        pmb->prad->CalculateRadiationInCellLinear(e_rad, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, k,
            j, i, gcov, pmb->prad->cons);
      }
    }
  }
  pmb->prad->ConservedToPrimitive(pmb->prad->cons, prim_rad, il - ngh, il - 1, jl, ju, kl,
      ku);
  return;
}

//----------------------------------------------------------------------------------------
// Fixed boundary condition for right side
// Inputs:
//   pmb: pointer to MeshBlock
//   pcoord: pointer to Coordinates
//   il, iu, jl, ju, kl, ku: indices demarkating active region
//   time, dt: current time and time step (unused)
//   ngh: effective number of ghost zones to populate
// Outputs:
//   prim: primitives set in ghost zones
//   bb: face-centered magnetic field set in ghost zones (unused)
//   prim_rad: radiation primitives set in ghost zones

void FixedRight(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, AthenaArray<Real> &prim_rad, Real time, Real dt, int il, int iu,
    int jl, int ju, int kl, int ku, int ngh) {

  // Allocate scratch arrays
  AthenaArray<Real> gcov, gcon;
  gcov.NewAthenaArray(NMETRIC, iu + ngh + 1);
  gcon.NewAthenaArray(NMETRIC, iu + ngh + 1);

  // Go through cells
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      pcoord->CellMetric(k, j, iu + 1, iu + ngh, gcov, gcon);
      for (int i = iu + 1; i <= iu + ngh; ++i) {

        // Calculate exact solution
        Real r = pcoord->x1v(i);
        Real alpha = std::sqrt(1.0 - 2.0 / r);
        Real ptot_cgs = var_b / SQR(SQR(alpha)) + std::exp(var_a * (alpha_in - alpha))
            * std::pow(alpha_in / alpha, gamma_adi / (gamma_adi - 1.0))
            * (ptot_in_cgs - var_b / SQR(SQR(alpha_in)));
        Real tt_cgs = tt_inf_cgs / alpha;
        Real prad_cgs = 1.0 / 3.0 * arad_cgs * SQR(SQR(tt_cgs));
        Real pgas_cgs = ptot_cgs - prad_cgs;
        Real rho_cgs = mu * mp_cgs / kb_cgs * pgas_cgs / tt_cgs;
        Real rho = rho_cgs / rho_unit;
        Real pgas = pgas_cgs / p_unit;

        // Set fluid
        prim(IDN,k,j,i) = rho;
        prim(IPR,k,j,i) = pgas;
        prim(IVX,k,j,i) = 0.0;
        prim(IVY,k,j,i) = 0.0;
        prim(IVZ,k,j,i) = 0.0;

        // Set radiation
        Real e_rad = pmb->prad->arad * SQR(SQR(pgas / rho));
        pmb->prad->CalculateRadiationInCellLinear(e_rad, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, k,
            j, i, gcov, pmb->prad->cons);
      }
    }
  }
  pmb->prad->ConservedToPrimitive(pmb->prad->cons, prim_rad, iu + 1, iu + ngh, jl, ju, kl,
      ku);
  return;
}
