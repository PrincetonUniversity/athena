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
#include <cmath>      // cos, hypot, sin
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
Real kappa_a, kappa_s;          // absorption and scattering opacities
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
Real lx, ly;                    // domain side lengths
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
        << "gr_rad_linear_wave only supports Minkowski coordinates" << std::endl;
    ATHENA_ERROR(msg);
    return;
  }
  if (pin->GetString("coord", "rad_tetrad") != "cartesian") {
    std::stringstream msg;
    msg << "### FATAL ERROR in Problem Generator\n"
        << "gr_rad_linear_wave only supports Cartesian tetrad" << std::endl;
    ATHENA_ERROR(msg);
    return;
  }
  if (pin->GetInteger("mesh", "nx3") > 1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Problem Generator\n"
        << "gr_rad_linear_wave only supports 1D and 2D" << std::endl;
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
  kappa_a = pin->GetReal("problem", "kappa_a");
  kappa_s = pin->GetReal("problem", "kappa_s");
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
  lx = pin->GetReal("mesh", "x1max") - pin->GetReal("mesh", "x1min");
  ly = 0.0;
  if (pin->GetInteger("mesh", "nx2") > 1) {
    ly = pin->GetReal("mesh", "x2max") - pin->GetReal("mesh", "x2min");
  }
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
//   Makes wavevector normal to box diagonal.
//   Assumes given x-direction is along wave.

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  // Calculate wavenumber
  Real kx = 2.0 * PI / lx;
  Real ky = 0.0;
  if (ly > 0.0) {
    kx = 2.0 * PI / ly;
    ky = 2.0 * PI / lx;
  }

  // Calculate rotation matrix
  Real lxy = std::hypot(lx, ly);
  Real sa = ly / lxy;
  Real ca = lx / lxy;

  // Initialize primitive hydro variables
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);
        Real s = std::sin(kx * x + ky * y);
        Real c = std::cos(kx * x + ky * y);
        Real rho_val = rho + delta * (drho_real * c - drho_imag * s);
        Real pgas_val = pgas + delta * (dpgas_real * c - dpgas_imag * s);
        Real ux_val_wave = ux + delta * (dux_real * c - dux_imag * s);
        Real uy_val_wave = uy + delta * (duy_real * c - duy_imag * s);
        Real uz_val_wave = uz + delta * (duz_real * c - duz_imag * s);
        Real ux_val = ca * ux_val_wave - sa * uy_val_wave;
        Real uy_val = sa * ux_val_wave + ca * uy_val_wave;
        Real uz_val = uz_val_wave;
        phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) = rho_val;
        phydro->w(IPR,k,j,i) = phydro->w1(IPR,k,j,i) = pgas_val;
        phydro->w(IVX,k,j,i) = phydro->w1(IVX,k,j,i) = ux_val;
        phydro->w(IVY,k,j,i) = phydro->w1(IVY,k,j,i) = uy_val;
        phydro->w(IVZ,k,j,i) = phydro->w1(IVZ,k,j,i) = uz_val;
      }
    }
  }

  // Initialize magnetic fields
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie+1; ++i) {
        Real x = pcoord->x1f(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);
        Real s = std::sin(kx * x + ky * y);
        Real c = std::cos(kx * x + ky * y);
        Real bbx_val_wave = bbx + delta * (dbbx_real * c - dbbx_imag * s);
        Real bby_val_wave = bby + delta * (dbby_real * c - dbby_imag * s);
        Real bbx_val = ca * bbx_val_wave - sa * bby_val_wave;
        pfield->b.x1f(k,j,i) = bbx_val;
      }
    }
  }
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je+1; ++j) {
      for (int i = is; i <= ie; ++i) {
        Real x = pcoord->x1v(i);
        Real y = pcoord->x2f(j);
        Real z = pcoord->x3v(k);
        Real s = std::sin(kx * x + ky * y);
        Real c = std::cos(kx * x + ky * y);
        Real bbx_val_wave = bbx + delta * (dbbx_real * c - dbbx_imag * s);
        Real bby_val_wave = bby + delta * (dbby_real * c - dbby_imag * s);
        Real bby_val = sa * bbx_val_wave + ca * bby_val_wave;
        pfield->b.x2f(k,j,i) = bby_val;
      }
    }
  }
  for (int k = ks; k <= ke+1; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3f(k);
        Real s = std::sin(kx * x + ky * y);
        Real c = std::cos(kx * x + ky * y);
        Real bbz_val_wave = bbz + delta * (dbbz_real * c - dbbz_imag * s);
        Real bbz_val = bbz_val_wave;
        pfield->b.x3f(k,j,i) = bbz_val;
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

        // Calculate location
        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);
        Real s = std::sin(kx * x + ky * y);
        Real c = std::cos(kx * x + ky * y);

        // Calculate wave-aligned coordinate-frame fluid velocity
        Real u_wave[4];
        u_wave[1] = ux + delta * (dux_real * c - dux_imag * s);
        u_wave[2] = uy + delta * (duy_real * c - duy_imag * s);
        u_wave[3] = uz + delta * (duz_real * c - duz_imag * s);
        u_wave[0] =
            std::hypot(1.0, std::hypot(u_wave[1], std::hypot(u_wave[2], u_wave[3])));

        // Calculate coordinate-frame fluid velocity
        Real u[4];
        u[0] = u_wave[0];
        u[1] = ca * u_wave[1] - sa * u_wave[2];
        u[2] = sa * u_wave[1] + ca * u_wave[2];
        u[3] = u_wave[3];

        // Calculate wave-aligned fluid-frame radiation moments
        Real rf_wave[4][4];
        rf_wave[0][0] = erad + delta * (derad_real * c - derad_imag * s);
        rf_wave[0][1] = rf_wave[1][0] = fxrad + delta * (dfxrad_real * c - dfxrad_imag * s);
        rf_wave[0][2] = rf_wave[2][0] = fyrad + delta * (dfyrad_real * c - dfyrad_imag * s);
        rf_wave[0][3] = rf_wave[3][0] = fzrad + delta * (dfzrad_real * c - dfzrad_imag * s);
        rf_wave[1][1] = 1.0 / 3.0 * rf_wave[0][0];
        rf_wave[1][2] = rf_wave[2][1] = 0.0;
        rf_wave[1][3] = rf_wave[3][1] = 0.0;
        rf_wave[2][2] = 1.0 / 3.0 * rf_wave[0][0];
        rf_wave[2][3] = rf_wave[3][2] = 0.0;
        rf_wave[3][3] = 1.0 / 3.0 * rf_wave[0][0];

        // Calculate wave-aligned coordinate-frame radiation moments
        Real lambda_c_f_wave[4][4];
        lambda_c_f_wave[0][0] = u_wave[0];
        lambda_c_f_wave[0][1] = lambda_c_f_wave[1][0] = u_wave[1];
        lambda_c_f_wave[0][2] = lambda_c_f_wave[2][0] = u_wave[2];
        lambda_c_f_wave[0][3] = lambda_c_f_wave[3][0] = u_wave[3];
        lambda_c_f_wave[1][1] = 1.0 + 1.0 / (1.0 + u_wave[0]) * SQR(u_wave[1]);
        lambda_c_f_wave[1][2] = lambda_c_f_wave[2][1] =
            1.0 / (1.0 + u_wave[0]) * u_wave[1] * u_wave[2];
        lambda_c_f_wave[1][3] = lambda_c_f_wave[3][1] =
            1.0 / (1.0 + u_wave[0]) * u_wave[1] * u_wave[3];
        lambda_c_f_wave[2][2] = 1.0 + 1.0 / (1.0 + u_wave[0]) * SQR(u_wave[2]);
        lambda_c_f_wave[2][3] = lambda_c_f_wave[3][2] =
            1.0 / (1.0 + u_wave[0]) * u_wave[2] * u_wave[3];
        lambda_c_f_wave[3][3] = 1.0 + 1.0 / (1.0 + u_wave[0]) * SQR(u_wave[3]);
        Real r_wave[4][4];
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            r_wave[a][b] = 0.0;
            for (int c = 0; c < 4; ++c) {
              for (int d = 0; d < 4; ++d) {
                r_wave[a][b] +=
                    lambda_c_f_wave[a][c] * lambda_c_f_wave[b][d] * rf_wave[c][d];
              }
            }
          }
        }

        // Calculate coordinate-frame radiation moments
        Real r[4][4];
        r[0][0] = r_wave[0][0];
        r[0][1] = r[1][0] = ca * r_wave[0][1] - sa * r_wave[0][2];
        r[0][2] = r[2][0] = sa * r_wave[0][1] + ca * r_wave[0][2];
        r[0][3] = r[3][0] = r_wave[0][3];
        r[1][1] = SQR(ca) * r_wave[1][1] - 2.0 * sa * ca * r_wave[1][2]
            + SQR(sa) * r_wave[2][2];
        r[1][2] = r[2][1] = sa * ca * r_wave[1][1] + (SQR(ca) - SQR(sa)) * r_wave[1][2]
            - sa * ca * r_wave[2][2];
        r[1][3] = r[3][1] = ca * r_wave[1][3] - sa * r_wave[2][3];
        r[2][2] = SQR(sa) * r_wave[1][1] + 2.0 * sa * ca * r_wave[1][2]
            + SQR(ca) * r_wave[2][2];
        r[2][3] = r[3][2] = sa * r_wave[1][3] + ca * r_wave[2][3];
        r[3][3] = r_wave[3][3];

        // Calculate fluid-frame radiation moments
        Real lambda_f_c[4][4];
        lambda_f_c[0][0] = u[0];
        lambda_f_c[0][1] = lambda_f_c[1][0] = -u[1];
        lambda_f_c[0][2] = lambda_f_c[2][0] = -u[2];
        lambda_f_c[0][3] = lambda_f_c[3][0] = -u[3];
        lambda_f_c[1][1] = 1.0 + 1.0 / (1.0 + u[0]) * SQR(u[1]);
        lambda_f_c[1][2] = lambda_f_c[2][1] = 1.0 / (1.0 + u[0]) * u[1] * u[2];
        lambda_f_c[1][3] = lambda_f_c[3][1] = 1.0 / (1.0 + u[0]) * u[1] * u[3];
        lambda_f_c[2][2] = 1.0 + 1.0 / (1.0 + u[0]) * SQR(u[2]);
        lambda_f_c[2][3] = lambda_f_c[3][2] = 1.0 / (1.0 + u[0]) * u[2] * u[3];
        lambda_f_c[3][3] = 1.0 + 1.0 / (1.0 + u[0]) * SQR(u[3]);
        Real rf[4][4];
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            rf[a][b] = 0.0;
            for (int c = 0; c < 4; ++c) {
              for (int d = 0; d < 4; ++d) {
                rf[a][b] += lambda_f_c[a][c] * lambda_f_c[b][d] * r[c][d];
              }
            }
          }
        }

        // Initialize radiation based on fluid velocity and fluid-frame quantities
        prad->CalculateRadiationInCellLinear(rf[0][0], rf[0][1], rf[0][2], rf[0][3], u[1],
            u[2], u[3], k, j, i, g, prad->cons);
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
  int il = pmb->is - NGHOST;
  int iu = pmb->ie + NGHOST;
  int jl = pmb->js;
  int ju = pmb->je;
  if (jl != ju) {
    jl -= NGHOST;
    ju += NGHOST;
  }
  int kl = pmb->ks;
  int ku = pmb->ke;
  if (kl != ku) {
    kl -= NGHOST;
    ku += NGHOST;
  }

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
