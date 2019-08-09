//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gr_rad_beam.cpp
//  \brief Problem generator for GR radiation with spatially constant initial conditions

// C++ headers
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str, string

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

// Global variables
static Real rho, pgas;   // initial thermodynamic variables for fluid
static Real ux, uy, uz;  // initial spatial components of fluid 4-velocity
static Real zs, ze;      // index bounds on zeta
static Real ps, pe;      // index bounds on psi

//----------------------------------------------------------------------------------------
// Function for preparing Mesh
// Inputs:
//   pin: input parameters
// Outputs: (none)

void Mesh::InitUserMeshData(ParameterInput *pin) {

  // Read parameters from input file
  rho = pin->GetReal("problem", "rho");
  pgas = pin->GetReal("problem", "pgas");
  ux = pin->GetReal("problem", "ux");
  uy = pin->GetReal("problem", "uy");
  uz = pin->GetReal("problem", "uz");
  zs = NGHOST;
  ze = zs + pin->GetInteger("radiation", "n_polar");
  ps = NGHOST;
  pe = ps + pin->GetInteger("radiation", "n_azimuthal");

  // Check coordinates
  if (COORDINATE_SYSTEM != std::string("minkowski")) {
    std::stringstream msg;
    msg << "### FATAL ERROR in problem generator\n";
    msg << "unsupported coordinate system\n";
    throw std::runtime_error(msg.str().c_str());
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for preparing MeshBlock
// Inputs:
//   pin: input parameters (unused)
// Outputs: (none)

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
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
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = prad->AngleInd(l, m);
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
          for (int i = is; i <= ie; ++i) {
            prad->prim(lm,k,j,i) = 0.0;
            prad->cons(lm,k,j,i) = 0.0;
          }
        }
      }
    }
  }
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
