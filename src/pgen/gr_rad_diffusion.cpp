//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gr_rad_beam.cpp
//  \brief Problem generator for GR radiation and diffusion of pulse

// C++ headers
#include <cstdlib>    // exit (needed for defs.hpp)
#include <cmath>      // exp
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
Real x0, t0;                  // initial location of peak and time since delta function
Real e_rad_back, e_rad_peak;  // initial radiation energy density parameters
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
  x0 = pin->GetReal("problem", "x0");
  t0 = pin->GetReal("problem", "t0");
  e_rad_back = pin->GetReal("problem", "e_rad_back");
  e_rad_peak = pin->GetReal("problem", "e_rad_peak");
  return;
}

//----------------------------------------------------------------------------------------
// Function for setting initial conditions
// Inputs:
//   pin: parameters (unused)
// Outputs: (none)

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
  Real kappa = prad->kappa;
  Real dd = 1.0 / (3.0 * kappa * rho);
  for (int i = is; i <= ie; ++i) {
    Real x = pcoord->x1v(i);
    Real e_rad = e_rad_peak * (e_rad_back + std::exp(-SQR(x) / (4.0 * dd * t0)));
    Real ii = e_rad / (4.0*PI);
    for (int n = 0; n < prad->nang; ++n) {
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
          prad->prim(n,k,j,i) = ii;
        }
      }
    }
  }
  prad->PrimitiveToConserved(prad->prim, prad->cons, pcoord, is, ie, js, je, ks, ke);

  // Initialize opacity (never changed)
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        prad->opacity(OPAS,k,j,i) = rho * kappa;
        prad->opacity(OPAA,k,j,i) = 0.0;
        prad->opacity(OPAP,k,j,i) = 0.0;
      }
    }
  }
  return;
}
