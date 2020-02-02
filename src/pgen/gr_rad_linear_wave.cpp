//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gr_rad_beam.cpp
//  \brief Problem generator for GR radiation with sinusoidal initial conditions

// C++ headers
#include <algorithm>  // max, min
#include <cstdlib>    // exit (needed for defs.hpp)
#include <cmath>      // cos, sin, sqrt
#include <iostream>   // cout (needed for defs.hpp), endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error (needed for defs.hpp)
#include <string>     // c_str, strcmp, string (all needed for defs.hpp)

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "../athena.hpp"                   // Real, enums
#include "../athena_arrays.hpp"            // AthenaArray
#include "../parameter_input.hpp"          // ParameterInput
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../eos/eos.hpp"                  // EquationOfState
#include "../field/field.hpp"              // Field
#include "../hydro/hydro.hpp"              // Hydro
#include "../radiation/radiation.hpp"      // Radiation

// Configuration checking
#if not RADIATION_ENABLED
#error "This problem generator must be used with radiation"
#endif
#if not GENERAL_RELATIVITY
#error "This problem generator must be used with general relativity"
#endif
#if not MAGNETIC_FIELDS_ENABLED
#error "This problem generator must be used with magnetic fields"
#endif

// Global variables
namespace {
Real rho, pgas;                 // background thermodynamic variables
Real ux, uy, uz;                // background velocity
Real bbx, bby, bbz;             // background magnetic field
Real erad;                      // background radiation energy density
Real fxrad, fyrad, fzrad;       // background radiation flux
Real drho_real, drho_imag;      // density perturbation
Real dpgas_real, dpgas_imag;    // pressure perturbation
Real dux_real, dux_imag;        // x-velocity perturbation
Real duy_real, duy_imag;        // y-velocity perturbation
Real duz_real, duz_imag;        // z-velocity perturbation
Real dbbx_real, dbbx_imag;      // x-field perturbation
Real dbby_real, dbby_imag;      // y-field perturbation
Real dbbz_real, dbbz_imag;      // z-field perturbation
Real derad_real, derad_imag;    // radiation energy density perturbation
Real dfxrad_real, dfxrad_imag;  // x-flux perturbation
Real dfyrad_real, dfyrad_imag;  // y-flux perturbation
Real dfzrad_real, dfzrad_imag;  // z-flux perturbation
Real delta;                     // amplitude of perturbation
Real lambda;                    // wavelength
}  // namespace

// Declarations
void Opacity(MeshBlock *pmb, const AthenaArray<Real> &prim);

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
  bbx = pin->GetReal("problem", "bbx");
  bby = pin->GetReal("problem", "bby");
  bbz = pin->GetReal("problem", "bbz");
  erad = pin->GetReal("problem", "erad");
  fxrad = pin->GetReal("problem", "fxrad");
  fyrad = pin->GetReal("problem", "fyrad");
  fzrad = pin->GetReal("problem", "fzrad");
  drho_real = pin->GetReal("problem", "drho_real");
  drho_imag = pin->GetReal("problem", "drho_imag");
  dpgas_real = pin->GetReal("problem", "dpgas_real");
  dpgas_imag = pin->GetReal("problem", "dpgas_imag");
  dux_real = pin->GetReal("problem", "dux_real");
  dux_imag = pin->GetReal("problem", "dux_imag");
  duy_real = pin->GetReal("problem", "duy_real");
  duy_imag = pin->GetReal("problem", "duy_imag");
  duz_real = pin->GetReal("problem", "duz_real");
  duz_imag = pin->GetReal("problem", "duz_imag");
  dbbx_real = pin->GetReal("problem", "dbbx_real");
  dbbx_imag = pin->GetReal("problem", "dbbx_imag");
  dbby_real = pin->GetReal("problem", "dbby_real");
  dbby_imag = pin->GetReal("problem", "dbby_imag");
  dbbz_real = pin->GetReal("problem", "dbbz_real");
  dbbz_imag = pin->GetReal("problem", "dbbz_imag");
  derad_real = pin->GetReal("problem", "derad_real");
  derad_imag = pin->GetReal("problem", "derad_imag");
  dfxrad_real = pin->GetReal("problem", "dfxrad_real");
  dfxrad_imag = pin->GetReal("problem", "dfxrad_imag");
  dfyrad_real = pin->GetReal("problem", "dfyrad_real");
  dfyrad_imag = pin->GetReal("problem", "dfyrad_imag");
  dfzrad_real = pin->GetReal("problem", "dfzrad_real");
  dfzrad_imag = pin->GetReal("problem", "dfzrad_imag");
  delta = pin->GetReal("problem", "delta");
  lambda = pin->GetReal("mesh", "x1max") - pin->GetReal("mesh", "x1min");
  return;
}

//----------------------------------------------------------------------------------------
// Function for preparing MeshBlock
// Inputs:
//   pin: input parameters (unused)
// Outputs: (none)

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  prad->EnrollOpacityFunction(Opacity);
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

  // Initialize primitive hydro variables
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        Real x = pcoord->x1v(i);
        Real s = std::sin(2.0*PI * x / lambda);
        Real c = std::cos(2.0*PI * x / lambda);
        phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) =
            rho + delta * (drho_real * c - drho_imag * s);
        phydro->w(IPR,k,j,i) = phydro->w1(IPR,k,j,i) =
            pgas + delta * (dpgas_real * c - dpgas_imag * s);
        phydro->w(IVX,k,j,i) = phydro->w1(IVX,k,j,i) =
            ux + delta * (dux_real * c - dux_imag * s);
        phydro->w(IVY,k,j,i) = phydro->w1(IVY,k,j,i) =
            uy + delta * (duy_real * c - duy_imag * s);
        phydro->w(IVZ,k,j,i) = phydro->w1(IVZ,k,j,i) =
            uz + delta * (duz_real * c - duz_imag * s);
      }
    }
  }

  // Initialize magnetic fields
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie+1; ++i) {
        Real x = pcoord->x1f(i);
        Real s = std::sin(2.0*PI * x / lambda);
        Real c = std::cos(2.0*PI * x / lambda);
        pfield->b.x1f(k,j,i) = bbx + delta * (dbbx_real * c - dbbx_imag * s);
      }
    }
  }
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je+1; ++j) {
      for (int i = is; i <= ie; ++i) {
        Real x = pcoord->x1v(i);
        Real s = std::sin(2.0*PI * x / lambda);
        Real c = std::cos(2.0*PI * x / lambda);
        pfield->b.x2f(k,j,i) = bby + delta * (dbby_real * c - dbby_imag * s);
      }
    }
  }
  for (int k = ks; k <= ke+1; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        Real x = pcoord->x1v(i);
        Real s = std::sin(2.0*PI * x / lambda);
        Real c = std::cos(2.0*PI * x / lambda);
        pfield->b.x3f(k,j,i) = bbz + delta * (dbbz_real * c - dbbz_imag * s);
      }
    }
  }
  pfield->CalculateCellCenteredField(pfield->b, pfield->bcc, pcoord, is, ie, js, je, ks,
      ke);

  // Initialize conserved hydro variables
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je,
      ks, ke);

  // Initialize radiation
  AthenaArray<Real> g, gi;
  g.NewAthenaArray(NMETRIC, ie + 1);
  gi.NewAthenaArray(NMETRIC, ie + 1);
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      pcoord->CellMetric(k, j, is, ie, g, gi);
      for (int i = is; i <= ie; ++i) {
        Real x = pcoord->x1v(i);
        Real s = std::sin(2.0*PI * x / lambda);
        Real c = std::cos(2.0*PI * x / lambda);
        Real e = erad + delta * (derad_real * c - derad_imag * s);
        Real fx = fxrad + delta * (dfxrad_real * c - dfxrad_imag * s);
        Real fy = fyrad + delta * (dfyrad_real * c - dfyrad_imag * s);
        Real fz = fzrad + delta * (dfzrad_real * c - dfzrad_imag * s);
        Real f = std::sqrt(SQR(fx) + SQR(fy) + SQR(fz)) / e;
        f = std::min(std::max(f, 0.0), 1.0);
        Real ux = 0.0, uy = 0.0, uz = 0.0;
        if (f > 0.0) {
          Real v = (2.0 - std::sqrt(4.0 - 3.0 * SQR(f))) / f;
          v = std::min(std::max(v, 0.0), 1.0);
          Real gamma = 1.0 / std::sqrt(1.0 - SQR(v));
          ux = gamma * v * fx / (f * e);
          uy = gamma * v * fy / (f * e);
          uz = gamma * v * fz / (f * e);
        }
        prad->CalculateRadiationInCell(e, ux, uy, uz, k, j, i, g, prad->cons);
      }
    }
  }
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

void Opacity(MeshBlock *pmb, const AthenaArray<Real> &prim)
{
  // Prepare index bounds
  int il = pmb->is;
  int iu = pmb->ie;
  int jl = pmb->js;
  int ju = pmb->je;
  int kl = pmb->ks;
  int ku = pmb->ke;

  // Set coefficients
  Radiation *prad = pmb->prad;
  Real kappaa = prad->kappa;
  Real kappas = 0.0;

  // Calculate opacity
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      for (int i = il; i <= iu; ++i) {
        Real rho = prim(IDN,k,j,i);
        prad->opacity(OPAA,k,j,i) = rho * kappaa;
        prad->opacity(OPAS,k,j,i) = rho * kappas;
      }
    }
  }
  return;
}
