//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gr_rad_beam.cpp
//  \brief Problem generator for 1D GR radiation shocks

// C++ headers
#include <cstdlib>    // exit (needed for defs.hpp)
#include <iostream>   // cout (needed for defs.hpp), endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error (needed for defs.hpp)
#include <string>     // c_str, strcmp, string (needed for defs.hpp)

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "../athena.hpp"                   // Real, enums
#include "../athena_arrays.hpp"            // AthenaArray
#include "../parameter_input.hpp"          // ParameterInput
#include "../bvals/bvals_interfaces.hpp"   // enums
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
Real x_shock;                // position of initial interface
Real rho_left, rho_right;    // initial gas density
Real pgas_left, pgas_right;  // initial gas pressure
Real ux_left, ux_right;      // initial gas x-velocity
Real uy_left, uy_right;      // initial gas y-velocity
Real uz_left, uz_right;      // initial gas z-velocity
Real erad_left, erad_right;  // initial fluid-frame radiation energy density
Real fx_left, fx_right;      // initial fluid-frame radiation x-flux
Real fy_left, fy_right;      // initial fluid-frame radiation y-flux
Real fz_left, fz_right;      // initial fluid-frame radiation z-flux
Real kappa_a, kappa_s;       // absorption and scattering opacities
}  // namespace

//----------------------------------------------------------------------------------------
// Function for preparing Mesh
// Inputs:
//   pin: input parameters
// Outputs: (none)

void Mesh::InitUserMeshData(ParameterInput *pin) {

  // Check coordinates
  if (std::strcmp(COORDINATE_SYSTEM, "minkowski") != 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Problem Generator\n"
        << "gr_rad_constant only supports Minkowski coordinates" << std::endl;
    ATHENA_ERROR(msg);
    return;
  }

  // Read parameters from input file
  x_shock = pin->GetReal("problem", "x_shock");
  rho_left = pin->GetReal("problem", "rho_left");
  rho_right = pin->GetReal("problem", "rho_right");
  pgas_left = pin->GetReal("problem", "pgas_left");
  pgas_right = pin->GetReal("problem", "pgas_right");
  ux_left = pin->GetReal("problem", "ux_left");
  ux_right = pin->GetReal("problem", "ux_right");
  uy_left = pin->GetReal("problem", "uy_left");
  uy_right = pin->GetReal("problem", "uy_right");
  uz_left = pin->GetReal("problem", "uz_left");
  uz_right = pin->GetReal("problem", "uz_right");
  erad_left = pin->GetReal("problem", "erad_left");
  erad_right = pin->GetReal("problem", "erad_right");
  fx_left = pin->GetReal("problem", "fx_left");
  fx_right = pin->GetReal("problem", "fx_right");
  fy_left = pin->GetReal("problem", "fy_left");
  fy_right = pin->GetReal("problem", "fy_right");
  fz_left = pin->GetReal("problem", "fz_left");
  fz_right = pin->GetReal("problem", "fz_right");
  kappa_a = pin->GetReal("problem", "kappa_a");
  kappa_s = pin->GetReal("problem", "kappa_s");

  // Enroll boundary functions
  if (mesh_bcs[BoundaryFace::inner_x1] == BoundaryFlag::user) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, FixedLeft);
  }
  if (mesh_bcs[BoundaryFace::outer_x1] == BoundaryFlag::user) {
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

  // Initialize fluid
  int il = is - NGHOST;
  int iu = ie + NGHOST;
  int jl = js - (ncells2 > 1 ? NGHOST : 0);
  int ju = je + (ncells2 > 1 ? NGHOST : 0);
  int kl = ks - (ncells3 > 1 ? NGHOST : 0);
  int ku = ke + (ncells3 > 1 ? NGHOST : 0);
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      for (int i = il; i <= iu; ++i) {
        Real x = pcoord->x1v(i);
        if (x <= x_shock) {
          phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) = rho_left;
          phydro->w(IPR,k,j,i) = phydro->w1(IPR,k,j,i) = pgas_left;
          phydro->w(IVX,k,j,i) = phydro->w1(IVX,k,j,i) = ux_left;
          phydro->w(IVY,k,j,i) = phydro->w1(IVY,k,j,i) = uy_left;
          phydro->w(IVZ,k,j,i) = phydro->w1(IVZ,k,j,i) = uz_left;
        } else {
          phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) = rho_right;
          phydro->w(IPR,k,j,i) = phydro->w1(IPR,k,j,i) = pgas_right;
          phydro->w(IVX,k,j,i) = phydro->w1(IVX,k,j,i) = ux_right;
          phydro->w(IVY,k,j,i) = phydro->w1(IVY,k,j,i) = uy_right;
          phydro->w(IVZ,k,j,i) = phydro->w1(IVZ,k,j,i) = uz_right;
        }
      }
    }
  }
  AthenaArray<Real> bb;
  peos->PrimitiveToConserved(phydro->w, bb, phydro->u, pcoord, il, iu, jl, ju, kl, ku);

  // Initialize radiation
  AthenaArray<Real> gcov, gcon;
  gcov.NewAthenaArray(NMETRIC, iu + 1);
  gcon.NewAthenaArray(NMETRIC, iu + 1);
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      pcoord->CellMetric(k, j, il, iu, gcov, gcon);
      for (int i = il; i <= iu; ++i) {
        Real x = pcoord->x1v(i);
        Real erad = erad_left;
        Real fx = fx_left;
        Real fy = fy_left;
        Real fz = fz_left;
        if (x > x_shock) {
          erad = erad_right;
          fx = fx_right;
          fy = fy_right;
          fz = fz_right;
        }
        prad->CalculateRadiationInCellLinear(erad, fx, fy, fz, phydro->w(IVX,k,j,i),
            phydro->w(IVY,k,j,i), phydro->w(IVZ,k,j,i), k, j, i, gcov, prad->cons);
      }
    }
  }

  // Initialize opacity (never changed)
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      for (int i = il; i <= iu; ++i) {
        prad->opacity(OPAA,k,j,i) = kappa_a;
        prad->opacity(OPAS,k,j,i) = kappa_s;
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
//   pcoord: pointer to Coordinates (unused)
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
  AthenaArray<Real> gcov, gcon;
  gcov.NewAthenaArray(NMETRIC, il);
  gcon.NewAthenaArray(NMETRIC, il);
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      pcoord->CellMetric(k, j, il - ngh, il - 1, gcov, gcon);
      for (int i = il - ngh; i <= il - 1; ++i) {
        prim(IDN,k,j,i) = rho_left;
        prim(IPR,k,j,i) = pgas_left;
        prim(IVX,k,j,i) = ux_left;
        prim(IVY,k,j,i) = uy_left;
        prim(IVZ,k,j,i) = uz_left;
        pmb->prad->CalculateRadiationInCellLinear(erad_left, fx_left, fy_left, fz_left,
            ux_left, uy_left, uz_left, k, j, i, gcov, pmb->prad->cons);
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
//   pcoord: pointer to Coordinates (unused)
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
  AthenaArray<Real> gcov, gcon;
  gcov.NewAthenaArray(NMETRIC, iu + ngh + 1);
  gcon.NewAthenaArray(NMETRIC, iu + ngh + 1);
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      pcoord->CellMetric(k, j, iu + 1, iu + ngh, gcov, gcon);
      for (int i = iu + 1; i <= iu + ngh; ++i) {
        prim(IDN,k,j,i) = rho_right;
        prim(IPR,k,j,i) = pgas_right;
        prim(IVX,k,j,i) = ux_right;
        prim(IVY,k,j,i) = uy_right;
        prim(IVZ,k,j,i) = uz_right;
        pmb->prad->CalculateRadiationInCellLinear(erad_right, fx_right, fy_right,
            fz_right, ux_right, uy_right, uz_right, k, j, i, gcov, pmb->prad->cons);
      }
    }
  }
  pmb->prad->ConservedToPrimitive(pmb->prad->cons, prim_rad, iu + 1, iu + ngh, jl, ju, kl,
      ku);
  return;
}
