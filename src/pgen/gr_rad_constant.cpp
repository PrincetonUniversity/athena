//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gr_rad_beam.cpp
//  \brief Problem generator for GR radiation with spatially constant initial conditions

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
Real zs, ze;                  // index bounds on zeta
Real ps, pe;                  // index bounds on psi
}  // namespace

// Declarations
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
  zs = NGHOST;
  ze = zs + pin->GetInteger("radiation", "n_polar");
  ps = NGHOST;
  pe = ps + pin->GetInteger("radiation", "n_azimuthal");
  return;
}

//----------------------------------------------------------------------------------------
// Function for preparing MeshBlock
// Inputs:
//   pin: input parameters (unused)
// Outputs: (none)

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {

  // Enroll opacity
  prad->EnrollOpacityFunction(TestOpacity);

  // Prepare moment outputs
  AllocateUserOutputVariables(4);
  SetUserOutputVariableName(0, "E");
  SetUserOutputVariableName(1, "M1");
  SetUserOutputVariableName(2, "M2");
  SetUserOutputVariableName(3, "M3");
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
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) = rho;
        phydro->w(IPR,k,j,i) = phydro->w1(IPR,k,j,i) = pgas;
        phydro->w(IVX,k,j,i) = phydro->w1(IVX,k,j,i) = ux;
        phydro->w(IVY,k,j,i) = phydro->w1(IVY,k,j,i) = uy;
        phydro->w(IVZ,k,j,i) = phydro->w1(IVZ,k,j,i) = uz;
      }
    }
  }
  AthenaArray<Real> bb;
  peos->PrimitiveToConserved(phydro->w, bb, phydro->u, pcoord, is, ie, js, je, ks, ke);

  // Initialize radiation
  prad->CalculateConstantRadiation(e_rad, ux_rad, uy_rad, uz_rad, prad->cons);
  return;
}

//----------------------------------------------------------------------------------------
// Function for preparing output
// Inputs:
//   pin: parameters (unused)
// Outputs: (none)
// Notes:
//   sets user_out_var array

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  prad->SetMoments(user_out_var);
  return;
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
  Real kappaa = 1.0;

  // Calculate opacity
  Radiation *prad = pmb->prad;
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      for (int i = il; i <= iu; ++i) {
        Real rho = prim(IDN,k,j,i);
        prad->opacity(OPAS,k,j,i) = rho * kappas;
        prad->opacity(OPAA,k,j,i) = rho * kappaa;
        prad->opacity(OPAP,k,j,i) = 0.0;
      }
    }
  }
  return;
}

