//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gr_rad_beam.cpp
//  \brief Problem generator for GR radiation with spatially constant initial conditions

// C++ headers
#include <cstdlib>    // exit (needed for defs.hpp)
#include <cmath>      // abs, fmin
#include <iostream>   // cout (needed for defs.hpp), endl
#include <limits>     // numeric_limits
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error (needed for defs.hpp)
#include <string>     // c_str, strcmp, string (needed for defs.hpp)

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

// Global variables
namespace {
Real rho, pgas;               // initial thermodynamic variables for fluid
Real ux, uy, uz;              // initial spatial components of fluid 4-velocity
Real e_rad;                   // initial coordinate-frame radiation energy density
Real ux_rad, uy_rad, uz_rad;  // initial spatial components of isotropic radiation frame
Real step_limit;              // factor to use in limiting timestep
}  // namespace

// Declarations
Real CouplingTimestep(MeshBlock *pmb);
void TestOpacity(MeshBlock *pmb, const AthenaArray<Real> &prim);

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
  if (pin->GetString("coord", "rad_tetrad") != "cartesian") {
    std::stringstream msg;
    msg << "### FATAL ERROR in Problem Generator\n"
        << "gr_rad_constant only supports Cartesian tetrad" << std::endl;
    ATHENA_ERROR(msg);
    return;
  }

  // Read parameters from input file
  rho = pin->GetReal("problem", "rho");
  pgas = pin->GetReal("problem", "pgas");
  ux = pin->GetReal("problem", "ux");
  uy = pin->GetReal("problem", "uy");
  uz = pin->GetReal("problem", "uz");
  e_rad = pin->GetReal("problem", "e_rad");
  ux_rad = pin->GetReal("problem", "ux_rad");
  uy_rad = pin->GetReal("problem", "uy_rad");
  uz_rad = pin->GetReal("problem", "uz_rad");
  step_limit = pin->GetReal("problem", "step_limit");

  // Enroll timestep limiter
  if (step_limit > 0.0) {
    EnrollUserTimeStepFunction(CouplingTimestep);
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for preparing MeshBlock
// Inputs:
//   pin: input parameters (unused)
// Outputs: (none)

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  prad->EnrollOpacityFunction(TestOpacity);
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
  int kl = ks - (ncells3 > 1 ? NGHOST : 0);
  int ku = ke + (ncells3 > 1 ? NGHOST : 0);
  int jl = js - (ncells2 > 1 ? NGHOST : 0);
  int ju = je + (ncells2 > 1 ? NGHOST : 0);
  int il = is - NGHOST;
  int iu = ie + NGHOST;
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      for (int i = il; i <= iu; ++i) {
        phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) = rho;
        phydro->w(IPR,k,j,i) = phydro->w1(IPR,k,j,i) = pgas;
        phydro->w(IVX,k,j,i) = phydro->w1(IVX,k,j,i) = ux;
        phydro->w(IVY,k,j,i) = phydro->w1(IVY,k,j,i) = uy;
        phydro->w(IVZ,k,j,i) = phydro->w1(IVZ,k,j,i) = uz;
      }
    }
  }
  AthenaArray<Real> bb;
  peos->PrimitiveToConserved(phydro->w, bb, phydro->u, pcoord, il, iu, jl, ju, kl, ku);

  // Initialize radiation
  prad->CalculateConstantRadiation(e_rad, ux_rad, uy_rad, uz_rad, prad->cons);
  return;
}

//----------------------------------------------------------------------------------------
// Timestep limiter
// Inputs:
//   pmb: pointer to MeshBlock
// Outputs:
//   returned value: minimum timescale for gas internal energy to change
// Notes:
//   Returns min(u_gas / |d(u_gas)/dt|) multiplied by step_limit.
//   Uses 0th moment: d(u_gas)/dt = kappa_a * rho * c * a_rad * (T_gas^4 - T_rad^4).

Real CouplingTimestep(MeshBlock *pmb)
{
  Real gamma_adi = pmb->peos->GetGamma();
  Real arad = pmb->prad->arad;
  Real dt_min = std::numeric_limits<Real>::max();
  if (pmb->prad->coupled_to_matter) {
    for (int k = pmb->ks; k <= pmb->ke; ++k) {
      for (int j = pmb->js; j <= pmb->je; ++j) {
        for (int i = pmb->is; i <= pmb->ie; ++i) {

          // Extract values
          Real rho = pmb->phydro->w(IDN,k,j,i);
          Real pgas = pmb->phydro->w(IPR,k,j,i);
          Real erad = pmb->prad->moments_fluid(0,k,j,i);
          Real k_a = pmb->prad->opacity(OPAA,k,j,i);

          // Calculate timescale
          Real ugas = pgas / (gamma_adi - 1.0);
          Real ttgas = pgas / rho;
          Real ttgas4 = SQR(SQR(ttgas));
          Real ttrad4 = erad / arad;
          Real dugas_dt = k_a * arad * (ttrad4 - ttgas4);
          Real dt = ugas / std::abs(dugas_dt);

          // Minimize timescale
          if (dt > 0.0) {
            dt_min = std::fmin(dt_min, dt);
          }
        }
      }
    }
  }
  return dt_min * step_limit;
}

//----------------------------------------------------------------------------------------
// Opacity Function
// Inputs:
//   pmb: pointer to MeshBlock
//   prim: primitive variables
// Outputs: (none)
// Notes:
//   Sets prad->opacity.

void TestOpacity(MeshBlock *pmb, const AthenaArray<Real> &prim)
{
  // Prepare index bounds
  int il = pmb->is - NGHOST;
  int iu = pmb->ie + NGHOST;
  int jl = pmb->js;
  int ju = pmb->je;
  int kl = pmb->ks;
  int ku = pmb->ke;
  if (ju > jl){
    jl -= NGHOST;
    ju += NGHOST;
  }
  if (ku > kl){
    kl -= NGHOST;
    ku += NGHOST;
  }

  // Set coefficients for electron scattering
  Real kappas = 0.0;
  Real kappaa = pmb->prad->kappa;

  // Calculate opacity
  Radiation *prad = pmb->prad;
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      for (int i = il; i <= iu; ++i) {
        Real rho = prim(IDN,k,j,i);
        prad->opacity(OPAS,k,j,i) = rho * kappas;
        prad->opacity(OPAA,k,j,i) = rho * kappaa;
        prad->opacity(OPAP,k,j,i) = rho * kappaa;
      }
    }
  }
  return;
}
