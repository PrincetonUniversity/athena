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
Real rho, pgas;                        // background thermodynamic variables
Real ux, uy, uz;                       // background velocity
Real bbx, bby, bbz;                    // background magnetic field
Real erad;                             // background radiation energy density
Real fxrad, fyrad, fzrad;              // background radiation flux
Real kappa_a, kappa_s;                 // absorption and scattering opacities
Real drho_real, drho_imag;             // density perturbation
Real dpgas_real, dpgas_imag;           // pressure perturbation
Real dux_real, dux_imag;               // x-velocity perturbation
Real duy_real, duy_imag;               // y-velocity perturbation
Real duz_real, duz_imag;               // z-velocity perturbation
Real dbbx_real, dbbx_imag;             // x-field perturbation
Real dbby_real, dbby_imag;             // y-field perturbation
Real dbbz_real, dbbz_imag;             // z-field perturbation
AthenaArray<Real> dii_real, dii_imag;  // radiation intensity perturbation
Real delta;                            // amplitude of perturbation
Real lambda;                           // wavelength
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
  kappa_a = pin->GetReal("problem", "kappa_a");
  kappa_s = pin->GetReal("problem", "kappa_s");
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
  int n_polar = pin->GetInteger("radiation", "n_polar");
  int n_azimuthal = pin->GetInteger("radiation", "n_azimuthal");
  int n_ang = n_polar * n_azimuthal;
  dii_real.NewAthenaArray(n_ang);
  dii_imag.NewAthenaArray(n_ang);
  for (int n = 0; n < n_ang; ++n) {
    std::stringstream real_name, imag_name;
    real_name << "dii" << n + 1 << "_real";
    imag_name << "dii" << n + 1 << "_imag";
    dii_real(n) = pin->GetReal("problem", real_name.str().c_str());
    dii_imag(n) = pin->GetReal("problem", imag_name.str().c_str());
  }
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
  Real gamma_sq = 1.0 + SQR(ux) + SQR(uy) + SQR(uz);
  Real erad = gamma_sq * prad->arad * SQR(SQR(pgas / rho));
  prad->CalculateConstantRadiation(erad, ux, uy, uz, prad->cons);
  prad->ConservedToPrimitive(prad->cons, prad->prim, is, ie, js, je, ks, ke);
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        Real x = pcoord->x1v(i);
        Real s = std::sin(2.0*PI * x / lambda);
        Real c = std::cos(2.0*PI * x / lambda);
        for (int l = prad->zs, n = 0; l <= prad->ze; ++l) {
          for (int m = prad->ps; m <= prad->pe; ++m, ++n) {
            int lm = prad->AngleInd(l, m);
            prad->prim(lm,k,j,i) += delta * (dii_real(n) * c - dii_imag(n) * s);
          }
        }
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
        prad->opacity(OPAA,k,j,i) = kappa_a;
        prad->opacity(OPAS,k,j,i) = kappa_s;
        prad->opacity(OPAP,k,j,i) = 0.0;
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
//   Keeps absorption coefficient, rather than opacity, constant in time.

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
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      for (int i = il; i <= iu; ++i) {
        Real rho_new = prim(IDN,k,j,i);
        prad->opacity(OPAA,k,j,i) = kappa_a * rho / rho_new;
        prad->opacity(OPAS,k,j,i) = kappa_s * rho / rho_new;
      }
    }
  }
  return;
}
