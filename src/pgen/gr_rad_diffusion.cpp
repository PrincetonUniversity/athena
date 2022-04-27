//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gr_rad_beam.cpp
//  \brief Problem generator for GR radiation and diffusion of pulse

// C++ headers
#include <algorithm>  // max
#include <cstdlib>    // exit (needed for defs.hpp)
#include <cmath>      // exp, sqrt
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
Real rho, pgas;   // initial thermodynamic variables for fluid
Real ut, ux;      // initial components of fluid 4-velocity
Real x0;          // initial location of peak
Real sigma;       // initial width of Gaussian
Real e_rad_peak;  // initial peak radiation energy density above background
Real e_rad_back;  // initial background radiation energy density
Real k_s;         // scattering coefficient
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
  x0 = pin->GetReal("problem", "x0");
  sigma = pin->GetReal("problem", "sigma");
  e_rad_peak = pin->GetReal("problem", "e_rad_peak");
  e_rad_back = pin->GetReal("problem", "e_rad_back");
  k_s = pin->GetReal("problem", "k_s");

  // Calculate Lorentz factor
  ut = std::sqrt(1.0 + SQR(ux));
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
        phydro->w(IVY,k,j,i) = phydro->w1(IVY,k,j,i) = 0.0;
        phydro->w(IVZ,k,j,i) = phydro->w1(IVZ,k,j,i) = 0.0;
      }
    }
  }
  AthenaArray<Real> bb;
  peos->PrimitiveToConserved(phydro->w, bb, phydro->u, pcoord, is, ie, js, je, ks, ke);

  // Initialize radiation
  AthenaArray<Real> g, gi;
  g.NewAthenaArray(NMETRIC, ie + 1);
  gi.NewAthenaArray(NMETRIC, ie + 1);
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      pcoord->CellMetric(k, j, is, ie, g, gi);
      for (int i = is; i <= ie; ++i) {

        // Locate cell in coordinate frame
        Real t = 0.0;
        Real x = pcoord->x1v(i);

        // Locate cell in fluid frame
        Real tp = ut * t - ux * x;
        Real xp = ut * x - ux * t;

        // Calculate fluid-frame moments
        Real aa = 0.0;
        if (tp > -3.0 * k_s * SQR(sigma) / 2.0) {
          aa = 1.0 / std::sqrt(1.0 + 2.0 * tp / (3.0 * k_s * SQR(sigma)));
        }
        Real ef_rad =
            e_rad_peak * aa * std::exp(-0.5 * SQR(aa * (xp - x0) / sigma)) + e_rad_back;
        Real ff_rad = SQR(aa / sigma) * (xp - x0) / (3.0 * k_s) * (ef_rad - e_rad_back);
        prad->CalculateRadiationInCellLinear(ef_rad, ff_rad, 0.0, 0.0, ux, 0.0, 0.0, k, j,
            i, g, prad->cons);
      }
    }
  }

  // Initialize opacity
  int il = is - NGHOST;
  int iu = ie + NGHOST;
  int jl = js;
  int ju = je;
  if (jl != ju) {
    jl -= NGHOST;
    ju += NGHOST;
  }
  int kl = ks;
  int ku = ke;
  if (kl != ku) {
    kl -= NGHOST;
    ku += NGHOST;
  }
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      for (int i = il; i <= iu; ++i) {
        prad->opacity(OPAS,k,j,i) = k_s / rho;
        prad->opacity(OPAA,k,j,i) = 0.0;
        prad->opacity(OPAP,k,j,i) = 0.0;
      }
    }
  }
  return;
}
