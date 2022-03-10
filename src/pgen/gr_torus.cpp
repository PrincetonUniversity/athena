//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gr_torus.cpp
//! \brief Problem generator for magnetized Fishbone-Moncrief torus.

// C headers

// C++ headers
#include <algorithm>  // max, min
#include <cmath>      // abs, acos, atan2, cos, exp, hypot, log, NAN, pow, sin, sqrt
#include <cstdlib>    // exit (needed for defs.hpp)
#include <cstring>    // strcmp
#include <iostream>   // cout (needed for defs.hpp), endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error (needed for defs.hpp)
#include <string>     // string (needed for defs.hpp)

// Athena++ headers
#include "../athena.hpp"                   // macros, enums, FaceField
#include "../athena_arrays.hpp"            // AthenaArray
#include "../bvals/bvals_interfaces.hpp"   // BoundaryFace
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../eos/eos.hpp"                  // EquationOfState
#include "../field/field.hpp"              // Field
#include "../hydro/hydro.hpp"              // Hydro
#include "../mesh/mesh.hpp"                // Mesh, MeshBlock
#include "../parameter_input.hpp"          // ParameterInput

// Configuration checking
#if !GENERAL_RELATIVITY
#error "This problem generator must be used with general relativity"
#endif

// Declarations
enum class MagneticFieldConfigs {density, loops};
void InflowBoundary(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ku,
    int ngh);
void OutflowBoundary(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ku,
    int ngh);
Real ThetaGrid(Real x, RegionSize rs);
Real HistorySum(MeshBlock *pmb, int iout);

// File declarations
namespace {
void GetKerrSchildCoordinates(Real x1, Real x2, Real x3, Real *p_r, Real *p_th,
    Real *p_ph);
void GetBoyerLindquistCoordinates(Real x1, Real x2, Real x3, Real *p_r, Real *p_th,
    Real *p_ph);
void TransformContravariantFromBoyerLindquist(Real at_bl, Real ar_bl, Real ath_bl,
    Real aph_bl, Real x1, Real x2, Real x3, Real *p_a0, Real *p_a1, Real *p_a2,
    Real *p_a3);
void TransformCovariantFromKerrSchild(Real a_r, Real a_th, Real a_ph, Real x1, Real x2,
    Real x3, Real *p_a_1, Real *p_a_2, Real *p_a_3);
Real CalculateLFromRPeak(Real r);
Real CalculateRPeakFromL(Real l_target);
Real LogHAux(Real r, Real sth);
void CalculateVelocityInTorus(Real r, Real sth, Real *p_ut, Real *p_uph);
void CalculateVelocityInTiltedTorus(Real r, Real th, Real ph, Real *p_ut, Real *p_ur,
    Real *p_uth, Real *p_uph);
Real IntegratedA1(Real x1_m, Real x1_p, Real x2, Real x3);
Real IntegratedA2(Real x1, Real x2_m, Real x2_p, Real x3);
Real IntegratedA3(Real x1, Real x2, Real x3_m, Real x3_p);
void VectorPotential(Real x1, Real x2, Real x3, Real *p_a_1, Real *p_a_2, Real *p_a_3);
} // namespace

// File variables
namespace {
Real m, a;                                    // black hole parameters
Real h_grid;                                  // grid compression parameter
Real gamma_adi;                               // adiabatic index
Real rho_min, rho_pow, pgas_min, pgas_pow;    // background parameters
bool prograde;                                // flag indicating disk is prograde
Real r_edge, r_peak, l, r_peak_max, rho_max;  // torus parameters
Real tilt;                                    // tilt angle
Real pert_amp, pert_kr, pert_kz;              // initial perturbations parameters
MagneticFieldConfigs field_config;            // type of magnetic field
Real pot_r_pow;                               // density vector potential parameters
Real pot_rho_pow, pot_rho_cutoff;             // density vector potential parameters
Real pot_r_min, pot_r_max, pot_r_num;         // loops vector potential parameters
Real pot_theta_min, pot_theta_num;            // loops vector potential parameters
Real pot_pgas_pow, pot_pgas_cutoff;           // loops vector potential parameters
Real pot_samples;                             // number of sample points for integrating
Real pot_amp;                                 // vector potential amplitude
Real log_h_edge, log_h_peak;                  // calculated torus parameters
Real pgas_over_rho_peak, rho_amp;             // calculated torus parameters
Real sin_tilt, cos_tilt;                      // calculated tilt parameters
int num_flux_radii;                           // number of spheres to use
Real *flux_radii;                             // locations to calculate fluxes
} // namespace

//----------------------------------------------------------------------------------------
// Function for preparing Mesh
// Inputs:
//   pin: input parameters (unused)
// Outputs: (none)

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Check for Kerr-Schild coordinates
  if (std::strcmp(COORDINATE_SYSTEM, "kerr-schild") != 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Problem Generator\n"
        << "GR torus only supports Kerr-Schild coordinates" << std::endl;
    ATHENA_ERROR(msg);
    return;
  }

  // Read coordinate parameters from input file
  m = pin->GetReal("coord", "m");
  a = pin->GetReal("coord", "a");
  Real x2rat = pin->GetReal("mesh", "x2rat");
  if (x2rat < 0.0) {
    h_grid = pin->GetOrAddReal("coord", "h", 1.0);
  }

  // Read fluid parameters from input file
  gamma_adi = pin->GetReal("hydro", "gamma");

  // Read torus parameters from input file
  prograde = pin->GetBoolean("problem", "prograde");
  r_edge = pin->GetReal("problem", "r_edge");
  r_peak = pin->GetReal("problem", "r_peak");
  if (r_peak < 0.0) {
    r_peak_max = pin->GetReal("problem", "r_peak_max");
    l = pin->GetReal("problem", "l");
  }
  rho_max = pin->GetReal("problem", "rho_max");

  // Read tilt parameter from input file
  tilt = pin->GetOrAddReal("problem", "tilt_angle", 0.0) * PI/180.0;

  // Read background parameters from input file
  rho_min = pin->GetReal("hydro", "rho_min");
  rho_pow = pin->GetReal("hydro", "rho_pow");
  pgas_min = pin->GetReal("hydro", "pgas_min");
  pgas_pow = pin->GetReal("hydro", "pgas_pow");

  // Read magnetic field parameters from input file
  if (MAGNETIC_FIELDS_ENABLED) {
    std::string field_config_str = pin->GetString("problem", "field_config");
    if (field_config_str == "density") {
      field_config = MagneticFieldConfigs::density;
      pot_r_pow = pin->GetReal("problem", "pot_r_pow");
      pot_rho_pow = pin->GetReal("problem", "pot_rho_pow");
      pot_rho_cutoff = pin->GetReal("problem", "pot_rho_cutoff");
      if (tilt != 0.0) {
        pot_samples = pin->GetReal("problem", "pot_samples");
      }
    } else if (field_config_str == "loops") {
      field_config = MagneticFieldConfigs::loops;
      pot_r_min = pin->GetReal("problem", "pot_r_min");
      pot_r_max = pin->GetReal("problem", "pot_r_max");
      pot_r_num = pin->GetReal("problem", "pot_r_num");
      pot_theta_min = pin->GetReal("problem", "pot_theta_min");
      pot_theta_num = pin->GetReal("problem", "pot_theta_num");
      pot_pgas_pow = pin->GetReal("problem", "pot_pgas_pow");
      pot_pgas_cutoff = pin->GetReal("problem", "pot_pgas_cutoff");
      if (tilt != 0.0) {
        pot_samples = pin->GetReal("problem", "pot_samples");
      }
    } else {
      std::stringstream msg;
      msg << "### FATAL ERROR in Problem Generator\n"
          << "unrecognized field_config input" << std::endl;
      ATHENA_ERROR(msg);
    }
    pot_amp = pin->GetReal("problem", "pot_amp");
  }

  // Read perturbation parameters from input file
  pert_amp = pin->GetOrAddReal("problem", "pert_amp", 0.0);
  pert_kr = pin->GetOrAddReal("problem", "pert_kr", 0.0);
  pert_kz = pin->GetOrAddReal("problem", "pert_kz", 0.0);

  // Read flux parameters from input file
  num_flux_radii = pin->GetOrAddInteger("problem", "num_flux_radii", 0);
  if (num_flux_radii > 0) {
    flux_radii = new Real[num_flux_radii];
    for (int n = 0; n < num_flux_radii; ++n) {
      std::stringstream flux_name;
      flux_name << "flux_radius_" << n + 1;
      flux_radii[n] = pin->GetReal("problem", flux_name.str().c_str());
    }
  }

  // Read boundary conditions from file
  std::string inner_boundary = pin->GetString("mesh", "ix1_bc");
  std::string outer_boundary = pin->GetString("mesh", "ox1_bc");

  // Enroll user-defined functions
  if (inner_boundary == "user") {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, InflowBoundary);
  }
  if (outer_boundary == "user") {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, OutflowBoundary);
  }
  if (x2rat < 0.0) {
    EnrollUserMeshGenerator(X2DIR, ThetaGrid);
  }
  if (num_flux_radii > 0) {
    AllocateUserHistoryOutput(num_flux_radii * 4);
    for (int n = 0; n < num_flux_radii; ++n) {
      std::stringstream mdot_name, edot_name, jdot_name, phi_name;
      mdot_name << "mdot_" << n + 1;
      edot_name << "edot_" << n + 1;
      jdot_name << "jdot_" << n + 1;
      phi_name << "phi_" << n + 1;
      EnrollUserHistoryOutput(n * 4, HistorySum, mdot_name.str().c_str());
      EnrollUserHistoryOutput(n * 4 + 1, HistorySum, edot_name.str().c_str());
      EnrollUserHistoryOutput(n * 4 + 2, HistorySum, jdot_name.str().c_str());
      EnrollUserHistoryOutput(n * 4 + 3, HistorySum, phi_name.str().c_str());
    }
  }

  // Calculate global constants describing torus and tilt
  if (r_peak >= 0.0) {
    l = CalculateLFromRPeak(r_peak);
  } else {
    r_peak = CalculateRPeakFromL(l);
  }
  log_h_edge = LogHAux(r_edge, 1.0);
  log_h_peak = LogHAux(r_peak, 1.0) - log_h_edge;
  pgas_over_rho_peak = (gamma_adi - 1.0) / gamma_adi * std::expm1(log_h_peak);
  rho_amp = rho_max / std::pow(pgas_over_rho_peak, 1.0 / (gamma_adi - 1.0));
  sin_tilt = std::sin(tilt);
  cos_tilt = std::cos(tilt);
  return;
}

//----------------------------------------------------------------------------------------
// Function for cleaning up Mesh
// Inputs:
//   pin: input parameters (unused)
// Outputs: (none)

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  if (num_flux_radii > 0) {
    delete[] flux_radii;
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for preparing MeshBlock
// Inputs:
//   pin: input parameters (unused)
// Outputs: (none)
// Notes:
//   User arrays are metric and its inverse.

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  // Allocate space for user output variables
  if (MAGNETIC_FIELDS_ENABLED) {
    AllocateUserOutputVariables(2);
    SetUserOutputVariableName(0, "gamma");
    SetUserOutputVariableName(1, "pmag");
  } else {
    AllocateUserOutputVariables(1);
    SetUserOutputVariableName(0, "gamma");
  }

  // Allocate space for scratch arrays
  AllocateRealUserMeshBlockDataField(num_flux_radii > 0 ? 5 : 2);
  ruser_meshblock_data[0].NewAthenaArray(NMETRIC, ie + NGHOST + 1);
  ruser_meshblock_data[1].NewAthenaArray(NMETRIC, ie + NGHOST + 1);

  // Allocate space for history variable computation
  if (num_flux_radii > 0) {
    ruser_meshblock_data[2].NewAthenaArray(num_flux_radii);
    ruser_meshblock_data[3].NewAthenaArray(num_flux_radii, je + 1);
    ruser_meshblock_data[4].NewAthenaArray(num_flux_radii, 4);
    AllocateIntUserMeshBlockDataField(1);
    iuser_meshblock_data[0].NewAthenaArray(num_flux_radii);
    for (int n = 0; n < num_flux_radii; ++n) {
      if (flux_radii[n] >= pcoord->x1f(is) && flux_radii[n] < pcoord->x1f(ie+1)) {
        int i;
        for (i = is; i <= ie + 1; ++i) {
          if (pcoord->x1v(i) > flux_radii[n]) {
            break;
          }
        }
        iuser_meshblock_data[0](n) = i - 1;
        ruser_meshblock_data[2](n) =
            (flux_radii[n] - pcoord->x1v(i-1)) / pcoord->dx1v(i-1);
        for (int j = js; j <= je; ++j) {
          Real cos_theta_m = std::cos(pcoord->x2f(j));
          Real cos_theta_p = std::cos(pcoord->x2f(j+1));
          ruser_meshblock_data[3](n,j) =
              1.0/3.0 * std::abs(cos_theta_m - cos_theta_p) * pcoord->dx3f(ks)
              * (3.0 * SQR(flux_radii[n]) + SQR(a) * (SQR(cos_theta_m)
              + cos_theta_m * cos_theta_p + SQR(cos_theta_p)));
        }
      } else {
        iuser_meshblock_data[0](n) = -1;
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for setting initial conditions
// Inputs:
//   pin: parameters
// Outputs: (none)
// Notes:
//   Initializes Fishbone-Moncrief torus.
//     Sets both primitive and conserved variables.
//   Defines and enrolls fixed r- and theta-direction boundary conditions.
//   References:
//     Fishbone & Moncrief 1976, ApJ 207 962 (FM)
//     Fishbone 1977, ApJ 215 323 (F)
//   Assumes x3 is axisymmetric direction.

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Prepare index bounds
  int il = is - NGHOST;
  int iu = ie + NGHOST;
  int jl = js;
  int ju = je;
  if (block_size.nx2 > 1) {
    jl -= NGHOST;
    ju += NGHOST;
  }
  int kl = ks;
  int ku = ke;
  if (block_size.nx3 > 1) {
    kl -= NGHOST;
    ku += NGHOST;
  }

  // Prepare scratch arrays
  AthenaArray<Real> &g = ruser_meshblock_data[0];
  AthenaArray<Real> &gi = ruser_meshblock_data[1];

  // Initialize primitive values
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      pcoord->CellMetric(k, j, il, iu, g, gi);
      for (int i = il; i <= iu; ++i) {
        // Extract preferred coordinates of cell
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real x3 = pcoord->x3v(k);

        // Calculate spherical Kerr-Schild coordinates of cell
        Real r, th, ph;
        GetKerrSchildCoordinates(x1, x2, x3, &r, &th, &ph);

        // Calculate Boyer-Lindquist coordinates of cell
        Real r_bl, th_bl, ph_bl;
        GetBoyerLindquistCoordinates(x1, x2, x3, &r_bl, &th_bl, &ph_bl);
        Real sth_bl = std::sin(th_bl);
        Real cth_bl = std::cos(th_bl);
        Real sph_bl = std::sin(ph_bl);
        Real cph_bl = std::cos(ph_bl);

        // Account for tilt
        Real sth_bl_t, cth_bl_t; // ph_bl_t;
        if (tilt == 0.0) {
          sth_bl_t = std::abs(sth_bl);
          cth_bl_t = cth_bl;
          //ph_bl_t = (sth_bl < 0.0) ? ph_bl - PI : ph_bl;
        } else {
          Real x_bl = sth_bl * cph_bl;
          Real y_bl = sth_bl * sph_bl;
          Real z_bl = cth_bl;
          Real x_bl_t = cos_tilt * x_bl - sin_tilt * z_bl;
          Real y_bl_t = y_bl;
          Real z_bl_t = sin_tilt * x_bl + cos_tilt * z_bl;
          sth_bl_t = std::hypot(x_bl_t, y_bl_t);
          cth_bl_t = z_bl_t;
          //ph_bl_t = std::atan2(y_bl_t, x_bl_t);
        }

        // Determine if we are in the torus
        Real log_h;
        bool in_torus = false;
        if (r_bl >= r_edge) {
          log_h = LogHAux(r_bl, sth_bl_t) - log_h_edge;  // (FM 3.6)
          if (log_h >= 0.0) {
            in_torus = true;
          }
        }

        // Calculate background primitives
        Real rho = rho_min * std::pow(r, rho_pow);
        Real pgas = pgas_min * std::pow(r, pgas_pow);
        Real uu1 = 0.0;
        Real uu2 = 0.0;
        Real uu3 = 0.0;

        // Overwrite primitives inside torus
        if (in_torus) {
          // Calculate thermodynamic variables
          Real pgas_over_rho = (gamma_adi - 1.0) / gamma_adi * std::expm1(log_h);
          rho = rho_amp * std::pow(pgas_over_rho, 1.0 / (gamma_adi - 1.0));
          pgas = pgas_over_rho * rho;

          // Calculate velocities in Boyer-Lindquist coordinates
          Real u0_bl, u1_bl, u2_bl, u3_bl;
          CalculateVelocityInTiltedTorus(r_bl, th_bl, ph_bl, &u0_bl, &u1_bl, &u2_bl,
              &u3_bl);

          // Transform to preferred coordinates
          Real u0, u1, u2, u3;
          TransformContravariantFromBoyerLindquist(u0_bl, 0.0, u2_bl, u3_bl, x1, x2, x3,
              &u0, &u1, &u2, &u3);
          uu1 = u1 - gi(I01,i) / gi(I00,i) * u0;
          uu2 = u2 - gi(I02,i) / gi(I00,i) * u0;
          uu3 = u3 - gi(I03,i) / gi(I00,i) * u0;
        }

        // Set primitive values, including cylindrically symmetric radial velocity
        // perturbations
        Real rr_bl = r_bl * sth_bl_t;
        Real z_bl = r_bl * cth_bl_t;
        Real amp_rel = 0.0;
        if (in_torus) {
          amp_rel = pert_amp * std::sin(pert_kr * rr_bl) * std::cos(pert_kz * z_bl);
        }
        Real amp_abs = amp_rel * uu3;
        Real pert_uur = rr_bl / r_bl * amp_abs;
        Real pert_uutheta = cth_bl / r_bl * amp_abs;
        phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) = rho;
        phydro->w(IPR,k,j,i) = phydro->w1(IPR,k,j,i) = pgas;
        phydro->w(IVX,k,j,i) = phydro->w1(IVX,k,j,i) = uu1 + pert_uur;
        phydro->w(IVY,k,j,i) = phydro->w1(IVY,k,j,i) = uu2 + pert_uutheta;
        phydro->w(IVZ,k,j,i) = phydro->w1(IVZ,k,j,i) = uu3;
      }
    }
  }

  // Initialize magnetic fields
  if (MAGNETIC_FIELDS_ENABLED) {
    // Initialize B^1
    for (int k = kl; k <= ku; ++k) {
      for (int j = jl; j <= ju; ++j) {
        for (int i = il; i <= iu + 1; ++i) {
          Real x1 = pcoord->x1f(i);
          Real x2_m = pcoord->x2f(j);
          Real x2_p = pcoord->x2f(j+1);
          Real x3_m = pcoord->x3f(k);
          Real x3_p = pcoord->x3f(k+1);
          Real a_2_3m = IntegratedA2(x1, x2_m, x2_p, x3_m);
          Real a_2_3p = IntegratedA2(x1, x2_m, x2_p, x3_p);
          Real a_3_2m = IntegratedA3(x1, x2_m, x3_m, x3_p);
          Real a_3_2p = IntegratedA3(x1, x2_p, x3_m, x3_p);
          Real area = pcoord->GetFace1Area(k,j,i);
          if (area > 0.0 && !(std::isnan(a_2_3m) || std::isnan(a_2_3p)
                || std::isnan(a_3_2m) || std::isnan(a_3_2p))) {
            pfield->b.x1f(k,j,i) = (a_3_2p - a_3_2m - a_2_3p + a_2_3m) / area;
          } else {
            pfield->b.x1f(k,j,i) = 0.0;
          }
        }
      }
    }

    // Initialize B^2
    for (int k = kl; k <= ku; ++k) {
      for (int j = jl; j <= ju + 1; ++j) {
        for (int i = il; i <= iu; ++i) {
          Real x1_m = pcoord->x1f(i);
          Real x1_p = pcoord->x1f(i+1);
          Real x2 = pcoord->x2f(j);
          Real x3_m = pcoord->x3f(k);
          Real x3_p = pcoord->x3f(k+1);
          Real a_1_3m = IntegratedA1(x1_m, x1_p, x2, x3_m);
          Real a_1_3p = IntegratedA1(x1_m, x1_p, x2, x3_p);
          Real a_3_1m = IntegratedA3(x1_m, x2, x3_m, x3_p);
          Real a_3_1p = IntegratedA3(x1_p, x2, x3_m, x3_p);
          Real area = pcoord->GetFace2Area(k,j,i);
          if (area > 0.0 && !(std::isnan(a_1_3m) || std::isnan(a_1_3p)
                || std::isnan(a_3_1m) || std::isnan(a_3_1p))) {
            pfield->b.x2f(k,j,i) = (a_1_3p - a_1_3m - a_3_1p + a_3_1m) / area;
          } else {
            pfield->b.x2f(k,j,i) = 0.0;
          }
        }
      }
    }

    // Initialize B^3
    for (int k = kl; k <= ku + 1; ++k) {
      for (int j = jl; j <= ju; ++j) {
        for (int i = il; i <= iu; ++i) {
          Real x1_m = pcoord->x1f(i);
          Real x1_p = pcoord->x1f(i+1);
          Real x2_m = pcoord->x2f(j);
          Real x2_p = pcoord->x2f(j+1);
          Real x3 = pcoord->x3f(k);
          Real a_1_2m = IntegratedA1(x1_m, x1_p, x2_m, x3);
          Real a_1_2p = IntegratedA1(x1_m, x1_p, x2_p, x3);
          Real a_2_1m = IntegratedA2(x1_m, x2_m, x2_p, x3);
          Real a_2_1p = IntegratedA2(x1_p, x2_m, x2_p, x3);
          Real area = pcoord->GetFace3Area(k,j,i);
          if (area > 0.0 && !(std::isnan(a_1_2m) || std::isnan(a_1_2p)
                || std::isnan(a_2_1m) || std::isnan(a_2_1p))) {
            pfield->b.x3f(k,j,i) = (a_2_1p - a_2_1m - a_1_2p + a_1_2m) / area;
          } else {
            pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
    }
    // Calculate cell-centered magnetic field
    pfield->CalculateCellCenteredField(pfield->b, pfield->bcc, pcoord, il, iu, jl, ju, kl,
                                       ku);
  }

  // Initialize conserved values
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, il, iu, jl, ju,
      kl, ku);

  // Call user work function to set output variables
  UserWorkInLoop();
  return;
}

//----------------------------------------------------------------------------------------
// Function responsible for storing history quantities for output
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Writes to ruser_meshblock_data[4] array the following quantities:
//     n,0: mdot (mass flux across specified radius)
//     n,1: edot (energy flux across specified radius)
//     n,2: jdot (angular momentum flux across specified radius)
//     n,3: phi (magnetic flux through specified radius)

void MeshBlock::UserWorkInLoop() {
  // Prepare scratch arrays
  AthenaArray<Real> &g = ruser_meshblock_data[0];
  AthenaArray<Real> &gi = ruser_meshblock_data[1];

  if (num_flux_radii > 0) {
    // Clear storage for history output accumulation
    for (int n = 0; n < num_flux_radii; ++n) {
      for (int m = 0; m < 4; ++m) {
        ruser_meshblock_data[4](n,m) = 0.0;
      }
    }
    // Go through all cells
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        pcoord->CellMetric(k, j, is, ie, g, gi);
        for (int i = is; i <= ie; ++i) {
          // Calculate normal-frame Lorentz factor
          Real uu1 = phydro->w(IVX,k,j,i);
          Real uu2 = phydro->w(IVY,k,j,i);
          Real uu3 = phydro->w(IVZ,k,j,i);
          Real tmp = g(I11,i) * SQR(uu1) + 2.0 * g(I12,i) * uu1 * uu2
              + 2.0 * g(I13,i) * uu1 * uu3 + g(I22,i) * SQR(uu2)
              + 2.0 * g(I23,i) * uu2 * uu3 + g(I33,i) * SQR(uu3);
          Real gamma = std::sqrt(1.0 + tmp);

          // Calculate 4-velocity
          Real alpha = std::sqrt(-1.0 / gi(I00,i));
          Real u0 = gamma / alpha;
          Real u1 = uu1 - alpha * gamma * gi(I01,i);
          Real u2 = uu2 - alpha * gamma * gi(I02,i);
          Real u3 = uu3 - alpha * gamma * gi(I03,i);
          Real u_0, u_1, u_2, u_3;
          pcoord->LowerVectorCell(u0, u1, u2, u3, k, j, i, &u_0, &u_1, &u_2, &u_3);

          // Calculate 4-magnetic field
          Real bb1 = 0.0, bb2 = 0.0, bb3 = 0.0;
          Real b0 = 0.0, b1 = 0.0, b2 = 0.0, b3 = 0.0;
          Real b_0 = 0.0, b_1 = 0.0, b_2 = 0.0, b_3 = 0.0;
          if (MAGNETIC_FIELDS_ENABLED) {
            bb1 = pfield->bcc(IB1,k,j,i);
            bb2 = pfield->bcc(IB2,k,j,i);
            bb3 = pfield->bcc(IB3,k,j,i);
            b0 = u_1 * bb1 + u_2 * bb2 + u_3 * bb3;
            b1 = (bb1 + b0 * u1) / u0;
            b2 = (bb2 + b0 * u2) / u0;
            b3 = (bb3 + b0 * u3) / u0;
            pcoord->LowerVectorCell(b0, b1, b2, b3, k, j, i, &b_0, &b_1, &b_2, &b_3);
          }

          // Calculate magnetic pressure
          Real b_sq = b0 * b_0 + b1 * b_1 + b2 * b_2 + b3 * b_3;

          // Check for contribution to history output
          for (int n = 0; n < num_flux_radii; ++n) {
            if (iuser_meshblock_data[0](n) == i || iuser_meshblock_data[0](n) + 1 == i) {
              Real factor = ruser_meshblock_data[2](n);
              if (iuser_meshblock_data[0](n) == i) {
                factor = 1.0 - ruser_meshblock_data[2](n);
              }
              Real rho = phydro->w(IDN,k,j,i);
              Real pgas = phydro->w(IPR,k,j,i);
              Real t1_0 = (rho + gamma_adi / (gamma_adi - 1.0) * pgas + b_sq) * u1 * u_0
                  - b1 * b_0;
              Real t1_3 = (rho + gamma_adi / (gamma_adi - 1.0) * pgas + b_sq) * u1 * u_3
                  - b1 * b_3;
              Real area = ruser_meshblock_data[3](n,j);
              ruser_meshblock_data[4](n,0) -= factor * rho * u1 * area;
              ruser_meshblock_data[4](n,1) += factor * t1_0 * area;
              ruser_meshblock_data[4](n,2) += factor * t1_3 * area;
              ruser_meshblock_data[4](n,3) +=
                  factor * 0.5 * std::sqrt(4.0*PI) * std::abs(bb1) * area;
            }
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function responsible for storing cell quantities for output
// Inputs:
//   pin: input parameters (unused)
// Outputs: (none)
// Notes:
//   Writes to user_out_var array the following quantities:
//     0,k,j,i: gamma (normal-frame Lorentz factor)
//     1,k,j,i: pmag (fluid-frame magnetic pressure, if magnetic fields enabled)

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  // Prepare scratch arrays
  AthenaArray<Real> &g = ruser_meshblock_data[0];
  AthenaArray<Real> &gi = ruser_meshblock_data[1];

  // Go through all cells
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      pcoord->CellMetric(k, j, is, ie, g, gi);
      for (int i = is; i <= ie; ++i) {
        // Calculate normal-frame Lorentz factor
        Real uu1 = phydro->w(IVX,k,j,i);
        Real uu2 = phydro->w(IVY,k,j,i);
        Real uu3 = phydro->w(IVZ,k,j,i);
        Real tmp = g(I11,i) * SQR(uu1) + 2.0 * g(I12,i) * uu1 * uu2
            + 2.0 * g(I13,i) * uu1 * uu3 + g(I22,i) * SQR(uu2)
            + 2.0 * g(I23,i) * uu2 * uu3 + g(I33,i) * SQR(uu3);
        Real gamma = std::sqrt(1.0 + tmp);
        user_out_var(0,k,j,i) = gamma;

        // Calculate 4-velocity
        Real alpha = std::sqrt(-1.0 / gi(I00,i));
        Real u0 = gamma / alpha;
        Real u1 = uu1 - alpha * gamma * gi(I01,i);
        Real u2 = uu2 - alpha * gamma * gi(I02,i);
        Real u3 = uu3 - alpha * gamma * gi(I03,i);
        Real u_0, u_1, u_2, u_3;
        pcoord->LowerVectorCell(u0, u1, u2, u3, k, j, i, &u_0, &u_1, &u_2, &u_3);

        // Calculate 4-magnetic field
        Real bb1 = 0.0, bb2 = 0.0, bb3 = 0.0;
        Real b0 = 0.0, b1 = 0.0, b2 = 0.0, b3 = 0.0;
        Real b_0 = 0.0, b_1 = 0.0, b_2 = 0.0, b_3 = 0.0;
        if (MAGNETIC_FIELDS_ENABLED) {
          bb1 = pfield->bcc(IB1,k,j,i);
          bb2 = pfield->bcc(IB2,k,j,i);
          bb3 = pfield->bcc(IB3,k,j,i);
          b0 = u_1 * bb1 + u_2 * bb2 + u_3 * bb3;
          b1 = (bb1 + b0 * u1) / u0;
          b2 = (bb2 + b0 * u2) / u0;
          b3 = (bb3 + b0 * u3) / u0;
          pcoord->LowerVectorCell(b0, b1, b2, b3, k, j, i, &b_0, &b_1, &b_2, &b_3);
        }

        // Calculate magnetic pressure
        Real b_sq = b0 * b_0 + b1 * b_1 + b2 * b_2 + b3 * b_3;
        if (MAGNETIC_FIELDS_ENABLED) {
          user_out_var(1,k,j,i) = b_sq / 2.0;
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Inflow boundary condition
// Inputs:
//   pmb: pointer to MeshBlock
//   pcoord: pointer to Coordinates
//   is,ie,js,je,ks,ke: indices demarkating active region
// Outputs:
//   prim: primitives set in ghost zones
//   bb: face-centered magnetic field set in ghost zones
// Notes:
//   Density and pressure are copied from last active cell into ghost cells.
//   Velocity is copied in the same way, except any primitive radial velocity (u^1 in the
//       normal observer frame) is made 0 if it is positive. This limiter is more
//       aggressive than acting on u^1 in the coordinate frame, but it is much simpler and
//       it guarantees a valid extrapolated state.
//   Magnetic field in the 2- and 3-directions is copied so as to preserve (exactly using
//       pointwise formulas, approximately in a finite-area sense) flux per unit area. For
//       example, 1-flux per unit area is
//           B^1 \sqrt{-det(g)} / \sqrt{g_{22} g_{33} - g_{23} g_{23}},
//       which is equal to
//           B^1 (g^{01} g^{01} - g^{00} g^{11})^{-1/2}
//       by Cramer's rule.
//   Magnetic field in the 1-direction is set so as to make the finite-volume
//       representation of divergence exactly 0 in each ghost cell.

void InflowBoundary(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ku,
    int ngh) {
  // Set hydro variables
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      for (int i = il - ngh; i <= il - 1; ++i) {
        prim(IDN,k,j,i) = prim(IDN,k,j,il);
        prim(IPR,k,j,i) = prim(IPR,k,j,il);
        prim(IVX,k,j,i) = std::min(prim(IVX,k,j,il), static_cast<Real>(0.0));
        prim(IVY,k,j,i) = prim(IVY,k,j,il);
        prim(IVZ,k,j,i) = prim(IVZ,k,j,il);
      }
    }
  }
  if (!MAGNETIC_FIELDS_ENABLED) {
    return;
  }

  // Prepare scratch arrays
  AthenaArray<Real> &g = pmb->ruser_meshblock_data[0];
  AthenaArray<Real> &gi = pmb->ruser_meshblock_data[1];

  // Set B^2
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju + 1; ++j) {
      pcoord->Face2Metric(k, j, il - ngh, il, g, gi);
      Real factor_active = 1.0 / std::sqrt(SQR(gi(I02,il)) - gi(I00,il) * gi(I22,il));
      for (int i = il - ngh; i <= il - 1; ++i) {
        Real factor_ghost = 1.0 / std::sqrt(SQR(gi(I02,i)) - gi(I00,i) * gi(I22,i));
        bb.x2f(k,j,i) = bb.x2f(k,j,il) * factor_active / factor_ghost;
      }
    }
  }

  // Set B^3
  for (int k = kl; k <= ku + 1; ++k) {
    for (int j = jl; j <= ju; ++j) {
      pcoord->Face3Metric(k, j, il - ngh, il, g, gi);
      Real factor_active = 1.0 / std::sqrt(SQR(gi(I03,il)) - gi(I00,il) * gi(I33,il));
      for (int i = il - ngh; i <= il - 1; ++i) {
        Real factor_ghost = 1.0 / std::sqrt(SQR(gi(I03,i)) - gi(I00,i) * gi(I33,i));
        bb.x3f(k,j,i) = bb.x3f(k,j,il) * factor_active / factor_ghost;
      }
    }
  }

  // Set B^1
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      for (int i = il - 1; i >= il - ngh; --i) {
        bb.x1f(k,j,i) = (pcoord->GetFace1Area(k,j,i+1) * bb.x1f(k,j,i+1)
            - pcoord->GetFace2Area(k,j,i) * bb.x2f(k,j,i)
            + pcoord->GetFace2Area(k,j+1,i) * bb.x2f(k,j+1,i)
            - pcoord->GetFace3Area(k,j,i) * bb.x3f(k,j,i)
            + pcoord->GetFace3Area(k+1,j,i) * bb.x3f(k+1,j,i))
            / pcoord->GetFace1Area(k,j,i);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Outflow boundary condition
// Inputs:
//   pmb: pointer to MeshBlock
//   pcoord: pointer to Coordinates
//   is,ie,js,je,ks,ke: indices demarkating active region
// Outputs:
//   prim: primitives set in ghost zones
//   bb: face-centered magnetic field set in ghost zones
// Notes:
//   Density and pressure are copied from last active cell into ghost cells.
//   Velocity is copied in the same way, except any primitive radial velocity (u^1 in the
//       normal observer frame) is made 0 if it is negative. This limiter is more
//       aggressive than acting on u^1 in the coordinate frame, but it is much simpler and
//       it guarantees a valid extrapolated state.
//   Magnetic field in the 2- and 3-directions is copied so as to preserve (exactly using
//       pointwise formulas, approximately in a finite-area sense) flux per unit area. For
//       example, 1-flux per unit area is
//           B^1 \sqrt{-det(g)} / \sqrt{g_{22} g_{33} - g_{23} g_{23}},
//       which is equal to
//           B^1 (g^{01} g^{01} - g^{00} g^{11})^{-1/2}
//       by Cramer's rule.
//   Magnetic field in the 1-direction is set so as to make the finite-volume
//       representation of divergence exactly 0 in each ghost cell.

void OutflowBoundary(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ku,
    int ngh) {
  // Set hydro variables
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      for (int i = iu + 1; i <= iu + ngh; ++i) {
        prim(IDN,k,j,i) = prim(IDN,k,j,iu);
        prim(IPR,k,j,i) = prim(IPR,k,j,iu);
        prim(IVX,k,j,i) = std::max(prim(IVX,k,j,iu), static_cast<Real>(0.0));
        prim(IVY,k,j,i) = prim(IVY,k,j,iu);
        prim(IVZ,k,j,i) = prim(IVZ,k,j,iu);
      }
    }
  }
  if (!MAGNETIC_FIELDS_ENABLED) {
    return;
  }

  // Prepare scratch arrays
  AthenaArray<Real> &g = pmb->ruser_meshblock_data[0];
  AthenaArray<Real> &gi = pmb->ruser_meshblock_data[1];

  // Set B^2
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju + 1; ++j) {
      pcoord->Face2Metric(k, j, iu, iu + ngh, g, gi);
      Real factor_active = 1.0 / std::sqrt(SQR(gi(I02,iu)) - gi(I00,iu) * gi(I22,iu));
      for (int i = iu + 1; i <= iu + ngh; ++i) {
        Real factor_ghost = 1.0 / std::sqrt(SQR(gi(I02,i)) - gi(I00,i) * gi(I22,i));
        bb.x2f(k,j,i) = bb.x2f(k,j,iu) * factor_active / factor_ghost;
      }
    }
  }

  // Set B^3
  for (int k = kl; k <= ku + 1; ++k) {
    for (int j = jl; j <= ju; ++j) {
      pcoord->Face3Metric(k, j, iu, iu + ngh, g, gi);
      Real factor_active = 1.0 / std::sqrt(SQR(gi(I03,iu)) - gi(I00,iu) * gi(I33,iu));
      for (int i = iu + 1; i <= iu + ngh; ++i) {
        Real factor_ghost = 1.0 / std::sqrt(SQR(gi(I03,i)) - gi(I00,i) * gi(I33,i));
        bb.x3f(k,j,i) = bb.x3f(k,j,iu) * factor_active / factor_ghost;
      }
    }
  }

  // Set B^1
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      for (int i = iu + 2; i <= iu + ngh + 1; ++i) {
        bb.x1f(k,j,i) = (pcoord->GetFace1Area(k,j,i-1) * bb.x1f(k,j,i-1)
            + pcoord->GetFace2Area(k,j,i-1) * bb.x2f(k,j,i-1)
            - pcoord->GetFace2Area(k,j+1,i-1) * bb.x2f(k,j+1,i-1)
            + pcoord->GetFace3Area(k,j,i-1) * bb.x3f(k,j,i-1)
            - pcoord->GetFace3Area(k+1,j,i-1) * bb.x3f(k+1,j,i-1))
            / pcoord->GetFace1Area(k,j,i);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Theta-spacing function
// Inputs:
//   x2: internal theta-like coordinate, scaled from 0 to 1
//   rs: region size struct
// Outputs:
//   returned value: corresponding theta-value
// Notes:
//   Implements remapping from Gammie, McKinney, & Toth 2003, ApJ 589 444.

Real ThetaGrid(Real x2, RegionSize rs) {
  return PI * x2 + (1.0 - h_grid) / 2.0 * std::sin(2.0*PI * x2);
}

//----------------------------------------------------------------------------------------
// History extraction
// Inputs:
//   pmb: pointer to MeshBlock
//   iout: index of history output
// Outputs:
//   returned value: block sum of corresponding variable

Real HistorySum(MeshBlock *pmb, int iout) {
  return pmb->ruser_meshblock_data[4](iout/4,iout%4);
}

//----------------------------------------------------------------------------------------
// Function for returning corresponding spherical Kerr-Schild coordinates of point
// Inputs:
//   x1,x2,x3: global coordinates to be converted
// Outputs:
//   p_r,p_th,p_ph: variables pointed to set to Kerr-Schild coordinates
// Notes:
//   Conversion is trivial in all currently implemented coordinate systems.

namespace {
void GetKerrSchildCoordinates(Real x1, Real x2, Real x3, Real *p_r, Real *p_th,
    Real *p_ph) {
  if (std::strcmp(COORDINATE_SYSTEM, "kerr-schild") == 0) {
    *p_r = x1;
    *p_th = x2;
    *p_ph = x3;
  } else {
    *p_r = NAN;
    *p_th = NAN;
    *p_ph = NAN;
  }
  return;
}
} // namespace

//----------------------------------------------------------------------------------------
// Function for returning corresponding Boyer-Lindquist coordinates of point
// Inputs:
//   x1,x2,x3: global coordinates to be converted
// Outputs:
//   p_r,p_th,p_ph: variables pointed to set to Boyer-Lindquist coordinates
// Notes:
//   Conversion is trivial in all currently implemented coordinate systems.

namespace {
void GetBoyerLindquistCoordinates(Real x1, Real x2, Real x3, Real *p_r, Real *p_th,
    Real *p_ph) {
  if (std::strcmp(COORDINATE_SYSTEM, "kerr-schild") == 0) {
    *p_r = x1;
    *p_th = x2;
    *p_ph = x3;
  } else {
    *p_r = NAN;
    *p_th = NAN;
    *p_ph = NAN;
  }
  return;
}
} // namespace

//----------------------------------------------------------------------------------------
// Function for transforming vector from Boyer-Lindquist to desired coordinates
// Inputs:
//   ar_bl,ar_bl,ath_bl,aph_bl: contravariant components in Boyer-Lindquist coordinates
//   x1,x2,x3: preferred coordinates of point
// Outputs:
//   p_a0,p_a1,p_a2,p_a3: contravariant components pointed to set in desired coordinates

namespace {
void TransformContravariantFromBoyerLindquist(Real at_bl, Real ar_bl, Real ath_bl,
    Real aph_bl, Real x1, Real x2, Real x3, Real *p_a0, Real *p_a1, Real *p_a2,
    Real *p_a3) {
  if (std::strcmp(COORDINATE_SYSTEM, "kerr-schild") == 0) {
    Real delta = SQR(x1) - 2.0 * m * x1 + SQR(a);
    *p_a0 = at_bl + 2.0 * m * x1 / delta * ar_bl;
    *p_a1 = ar_bl;
    *p_a2 = ath_bl;
    *p_a3 = aph_bl + a / delta * ar_bl;
  } else {
    *p_a0 = NAN;
    *p_a1 = NAN;
    *p_a2 = NAN;
    *p_a3 = NAN;
  }
  return;
}
} // namespace

//----------------------------------------------------------------------------------------
// Function for transforming covector from Kerr-Schild to desired coordinates
// Inputs:
//   a_r,a_th,a_ph: covariant components in Kerr-Schild coordinates
//   x1,x2,x3: preferred coordinates of point
// Outputs:
//   p_a_1,p_a_2,p_a_3: covariant components pointed to set in desired coordinates
// Notes:
//   Ignores time components. This is only valid if dt/dx^i = 0.

namespace {
void TransformCovariantFromKerrSchild(Real a_r, Real a_th, Real a_ph, Real x1, Real x2,
    Real x3, Real *p_a_1, Real *p_a_2, Real *p_a_3) {
  if (std::strcmp(COORDINATE_SYSTEM, "kerr-schild") == 0) {
    // Real delta = SQR(x1) - 2.0 * m * x1 + SQR(a);
    *p_a_1 = a_r;
    *p_a_2 = a_th;
    *p_a_3 = a_ph;
  } else {
    *p_a_1 = NAN;
    *p_a_2 = NAN;
    *p_a_3 = NAN;
  }
  return;
}
} // namespace

//----------------------------------------------------------------------------------------
// Function for calculating angular momentum variable l
// Inputs:
//   r: desired radius of pressure maximum
// Outputs:
//   returned value: l = u^t u_\phi such that pressure maximum occurs at r_peak
// Notes:
//   Beware many different definitions of l abound; this is *not* -u_phi/u_t.
//   Harm has a similar function: lfish_calc() in init.c
//     Harm's function assumes M = 1 and that corotation is desired.
//     It is equivalent to this, though seeing this requires much manipulation.
//   Implements (3.8) from Fishbone & Moncrief 1976, ApJ 207 962.
//   Allows for retrograde torus.
//   See CalculateRPeakFromL().

namespace {
Real CalculateLFromRPeak(Real r) {
  Real num, denom;
  if (prograde) {
    num = SQR(SQR(r)) + SQR(a * r) - 2.0 * m * SQR(a) * r
        - a * (SQR(r) - SQR(a)) * std::sqrt(m * r);
    denom = SQR(r) - 3.0 * m * r + 2.0 * a * std::sqrt(m * r);
  } else {
    num = -SQR(SQR(r)) - SQR(a * r) + 2.0 * m * SQR(a) * r
        - a * (SQR(r) - SQR(a)) * std::sqrt(m * r);
    denom = SQR(r) - 3.0 * m * r - 2.0 * a * std::sqrt(m * r);
  }
  return 1.0 / r * std::sqrt(m / r) * num / denom;
}
} // namespace

//----------------------------------------------------------------------------------------
// Function for calculating pressure maximum radius r_peak
// Inputs:
//   l_target: desired u^t u_\phi
// Outputs:
//   returned value: location of pressure maximum given l_target
// Notes:
//   Beware many different definitions of l abound; this is *not* -u_phi/u_t.
//   Uses (3.8) from Fishbone & Moncrief 1976, ApJ 207 962.
//   Allows for retrograde torus.
//   Uses bisection to find r such that formula for l agrees with given value.
//   Proceeds until either absolute tolerance is met.
//   Returns best value after max_iterations reached if tolerances not met.
//   Returns NAN in case of failure (e.g. root not bracketed).
//   See CalculateLFromRPeak().

namespace {
Real CalculateRPeakFromL(Real l_target) {
  // Parameters
  const Real tol_r = 1.0e-10;      // absolute tolerance on abscissa r_peak
  const Real tol_l = 1.0e-10;      // absolute tolerance on ordinate l
  const int max_iterations = 100;  // maximum number of iterations before best res

  // Prepare initial values
  Real r_a = r_edge;
  Real r_b = r_peak_max;
  Real r_c = 0.5 * (r_a + r_b);
  Real l_a = CalculateLFromRPeak(r_a);
  Real l_b = CalculateLFromRPeak(r_b);
  Real l_c = CalculateLFromRPeak(r_c);
  if (!((l_a < l_target && l_b > l_target) || (l_a > l_target && l_b < l_target))) {
    return NAN;
  }

  // Find root
  for (int n = 0; n < max_iterations; ++n) {
    if (std::abs(r_b - r_a) <= 2.0 * tol_r || std::abs(l_c - l_target) <= tol_l) {
      break;
    }
    if ((l_a < l_target && l_c < l_target) || (l_a > l_target && l_c > l_target)) {
      r_a = r_c;
      l_a = l_c;
    } else {
      r_b = r_c;
      //l_b = l_c;
    }
    r_c = 0.5 * (r_a + r_b);
    l_c = CalculateLFromRPeak(r_c);
  }
  return r_c;
}
} // namespace

//----------------------------------------------------------------------------------------
// Function for helping to calculate enthalpy
// Inputs:
//   r: radial Boyer-Lindquist coordinate
//   sth: sine of polar Boyer-Lindquist coordinate
// Outputs:
//   returned value: log(h)
// Notes:
//   Enthalpy defined here as h = p_gas/rho.
//   References Fishbone & Moncrief 1976, ApJ 207 962 (FM).
//   Implements first half of (FM 3.6).

namespace {
Real LogHAux(Real r, Real sth) {
  Real sth2 = SQR(sth);
  Real cth2 = 1.0 - sth2;
  Real delta = SQR(r) - 2.0 * m * r + SQR(a);                // \Delta
  Real sigma = SQR(r) + SQR(a) * cth2;                       // \Sigma
  Real aa = SQR(SQR(r) + SQR(a)) - delta * SQR(a) * sth2;    // A
  Real exp_2nu = sigma * delta / aa;                         // exp(2\nu) (FM 3.5)
  Real exp_2psi = aa / sigma * sth2;                         // exp(2\psi) (FM 3.5)
  Real exp_neg2chi = exp_2nu / exp_2psi;                     // exp(-2\chi) (cf. FM 2.15)
  Real omega = 2.0 * m * a * r / aa;                         // \omega (FM 3.5)
  Real var_a = std::sqrt(1.0 + 4.0 * SQR(l) * exp_neg2chi);
  Real var_b = 0.5 * std::log((1.0 + var_a)
                              / (sigma * delta / aa));
  Real var_c = -0.5 * var_a;
  Real var_d = -l * omega;
  return var_b + var_c + var_d;                              // (FM 3.4)
}
} // namespace

//----------------------------------------------------------------------------------------
// Function for computing 4-velocity components at a given position inside untilted torus
// Inputs:
//   r: radial Boyer-Lindquist coordinate
//   sth: sine of polar Boyer-Lindquist coordinate
// Outputs:
//   p_ut: Boyer-Lindquist u^t set
//   p_uph: Boyer-Lindquist u^\phi set
// Notes:
//   The formula for u^3 as a function of u_{(\phi)} is tedious to derive, but this
//       matches the formula used in Harm (init.c).

namespace {
void CalculateVelocityInTorus(Real r, Real sth, Real *p_ut, Real *p_uph) {
  Real sth2 = SQR(sth);
  Real cth2 = 1.0 - sth2;
  Real delta = SQR(r) - 2.0 * m * r + SQR(a);                // \Delta
  Real sigma = SQR(r) + SQR(a) * cth2;                       // \Sigma
  Real aa = SQR(SQR(r) + SQR(a)) - delta * SQR(a) * sth2;    // A
  Real exp_2nu = sigma * delta / aa;                         // exp(2\nu) (FM 3.5)
  Real exp_2psi = aa / sigma * sth2;                         // exp(2\psi) (FM 3.5)
  Real exp_neg2chi = exp_2nu / exp_2psi;                     // exp(-2\chi) (cf. FM 2.15)
  Real u_phi_proj_a = 1.0 + 4.0 * SQR(l) * exp_neg2chi;
  Real u_phi_proj_b = -1.0 + std::sqrt(u_phi_proj_a);
  Real u_phi_proj = std::sqrt(0.5 * u_phi_proj_b);           // (FM 3.3)
  u_phi_proj *= prograde ? 1.0 : -1.0;
  Real u3_a =
      (1.0 + SQR(u_phi_proj)) / (aa * sigma * delta);
  Real u3_b = 2.0 * m * a * r * std::sqrt(u3_a);
  Real u3_c = std::sqrt(sigma / aa) / sth;
  Real u3 = u3_b + u3_c * u_phi_proj;
  Real g_00 = -(1.0 - 2.0 * m * r / sigma);
  Real g_03 = -2.0 * m * a * r / sigma * sth2;
  Real g_33 = (sigma + (1.0 + 2.0 * m * r / sigma) * SQR(a)
               * sth2) * sth2;
  Real u0_a = (SQR(g_03) - g_00 * g_33) * SQR(u3);
  Real u0_b = std::sqrt(u0_a - g_00);
  Real u0 = -1.0 / g_00 * (g_03 * u3 + u0_b);
  *p_ut = u0;
  *p_uph = u3;
  return;
}
} // namespace

//----------------------------------------------------------------------------------------
// Function for computing 4-velocity components at a given position inside tilted torus
// Inputs:
//   r: radial Boyer-Lindquist coordinate
//   th,ph: polar and azimuthal Boyer-Lindquist coordinates aligned with black hole
// Outputs:
//   p_ut,p_ur,p_uth,p_uph: u^\mu set (Boyer-Lindquist coordinates)
// Notes:
//   First finds corresponding location in untilted torus.
//   Next calculates velocity at that point in untilted case.
//   Finally transforms that velocity into coordinates in which torus is tilted.

namespace {
void CalculateVelocityInTiltedTorus(Real r, Real th, Real ph, Real *p_ut, Real *p_ur,
    Real *p_uth, Real *p_uph) {
  // Calculate corresponding location
  Real sth = std::sin(th);
  Real cth = std::cos(th);
  Real sph = std::sin(ph);
  Real cph = std::cos(ph);
  Real sth_t, cth_t, ph_t;
  if (tilt == 0.0) {
    sth_t = std::abs(sth);
    cth_t = cth;
    ph_t = (sth < 0.0) ? ph - PI : ph;
  } else {
    Real x = sth * cph;
    Real y = sth * sph;
    Real z = cth;
    Real x_t = cos_tilt * x - sin_tilt * z;
    Real y_t = y;
    Real z_t = sin_tilt * x + cos_tilt * z;
    sth_t = std::hypot(x_t, y_t);
    cth_t = z_t;
    ph_t = std::atan2(y_t, x_t);
  }
  Real sph_t = std::sin(ph_t);
  Real cph_t = std::cos(ph_t);

  // Calculate untilted velocity
  Real u0_t, u3_t;
  CalculateVelocityInTorus(r, sth_t, &u0_t, &u3_t);
  Real u1_t = 0.0;
  Real u2_t = 0.0;

  // Account for tilt
  *p_ut = u0_t;
  *p_ur = u1_t;
  if (tilt == 0.0) {
    *p_uth = u2_t;
    *p_uph = u3_t;
  } else {
    Real dth_dth_t = (cos_tilt * sth_t + sin_tilt * cth_t * cph_t) / sth;
    Real dth_dph_t = -sin_tilt * sth_t * sph_t / sth;
    Real dph_dth_t = sin_tilt * sph_t / SQR(sth);
    Real dph_dph_t = sth_t / SQR(sth) * (cos_tilt * sth_t + sin_tilt * cth_t * cph_t);
    *p_uth = dth_dth_t * u2_t + dth_dph_t * u3_t;
    *p_uph = dph_dth_t * u2_t + dph_dph_t * u3_t;
  }
  if (sth < 0.0) {
    *p_uth *= -1.0;
    *p_uph *= -1.0;
  }
  return;
}
} // namespace

//----------------------------------------------------------------------------------------
// Function for calculating integrated 1-component of vector potential
// Inputs:
//   x1_m,x1_p: x^1 limits of integration
//   x2: x^2
//   x3: x^3
// Outputs:
//   returned value: integrated value
// Notes:
//   Evaluates \int A_1 dx^1.

namespace {
Real IntegratedA1(Real x1_m, Real x1_p, Real x2, Real x3) {
  // Multiple loops and density isocontour configurations
  if (field_config == MagneticFieldConfigs::density
      || field_config == MagneticFieldConfigs::loops) {
    // Closed form in spherical Kerr-Schild coordinates
    if (std::strcmp(COORDINATE_SYSTEM, "kerr-schild") == 0) {
      return 0.0;

    // General form for other coordinate systems
    } else {
      Real a_1_sum = 0.0;
      for (int n = 0; n < pot_samples; ++n) {
        Real x1 = x1_m
            + static_cast<Real>(n) / static_cast<Real>(pot_samples - 1) * (x1_p - x1_m);
        Real a_1, a_2, a_3;
        VectorPotential(x1, x2, x3, &a_1, &a_2, &a_3);
        a_1_sum += (n == 0 || n == pot_samples - 1) ? 0.5 * a_1 : a_1;
      }
      return a_1_sum * (x1_p - x1_m) / (pot_samples - 1);
    }

  // Unknown configuration
  } else {
    return NAN;
  }
}
} // namespace

//----------------------------------------------------------------------------------------
// Function for calculating integrated 2-component of vector potential
// Inputs:
//   x1: x^1
//   x2_m,x2_p: x^2 limits of integration
//   x3: x^3
// Outputs:
//   returned value: integrated value
// Notes:
//   Evaluates \int A_2 dx^2.

namespace {
Real IntegratedA2(Real x1, Real x2_m, Real x2_p, Real x3) {
  // Multiple loops and density isocontour configurations
  if (field_config == MagneticFieldConfigs::density
      || field_config == MagneticFieldConfigs::loops) {
    // Closed form for untilted disk in spherical Kerr-Schild coordinates
    if (std::strcmp(COORDINATE_SYSTEM, "kerr-schild") == 0 && tilt == 0.0) {
      return 0.0;

    // General form for tilted disk or other coordinate systems
    } else {
      Real a_2_sum = 0.0;
      for (int n = 0; n < pot_samples; ++n) {
        Real x2 = x2_m
            + static_cast<Real>(n) / static_cast<Real>(pot_samples - 1) * (x2_p - x2_m);
        Real a_1, a_2, a_3;
        VectorPotential(x1, x2, x3, &a_1, &a_2, &a_3);
        a_2_sum += (n == 0 || n == pot_samples - 1) ? 0.5 * a_2 : a_2;
      }
      return a_2_sum * (x2_p - x2_m) / (pot_samples - 1);
    }

  // Unknown configuration
  } else {
    return NAN;
  }
}
} // namespace

//----------------------------------------------------------------------------------------
// Function for calculating integrated 3-component of vector potential
// Inputs:
//   x1: x^1
//   x2: x^2
//   x3_m,x3_p: x^3 limits of integration
// Outputs:
//   returned value: integrated value
// Notes:
//   Evaluates \int A_3 dx^3.

namespace {
Real IntegratedA3(Real x1, Real x2, Real x3_m, Real x3_p) {
  // Multiple loops and density isocontour configurations
  if (field_config == MagneticFieldConfigs::density
      || field_config == MagneticFieldConfigs::loops) {
    // Closed form for untilted disk in spherical Kerr-Schild coordinates
    if (std::strcmp(COORDINATE_SYSTEM, "kerr-schild") == 0 && tilt == 0.0) {
      Real a_r, a_th, a_ph;
      VectorPotential(x1, x2, 0.5 * (x3_m + x3_p), &a_r, &a_th, &a_ph);
      return a_ph * (x3_p - x3_m);

    // General form for tilted disk or other coordinate systems
    } else {
      Real a_3_sum = 0.0;
      for (int n = 0; n < pot_samples; ++n) {
        Real x3 = x3_m
            + static_cast<Real>(n) / static_cast<Real>(pot_samples - 1) * (x3_p - x3_m);
        Real a_1, a_2, a_3;
        VectorPotential(x1, x2, x3, &a_1, &a_2, &a_3);
        a_3_sum += (n == 0 || n == pot_samples - 1) ? 0.5 * a_3 : a_3;
      }
      return a_3_sum * (x3_p - x3_m) / (pot_samples - 1);
    }

  // Unknown configuration
  } else {
    return NAN;
  }
}
} // namespace

//----------------------------------------------------------------------------------------
// Function for calculating covariant spatial components of vector potential
// Inputs:
//   x1,x2,x3: x^1, x^2, x^3
// Outputs:
//   p_a_1,p_a_2,p_a_3: values of A_1, A_2, and A_3 set
// Notes:
//   References Fishbone & Moncrief 1976, ApJ 207 962 (FM).
//   field_config == density: In Kerr-Schild coordinates aligned with the torus:
//     A_r = 0,
//     A_{\theta'} = 0,
//     A_{\phi'} = pot_amp * \max(\rho - pot_rho_cutoff, 0)^pot_rho_pow * r^pot_r_pow.
//   field_config == loops: In Kerr-Schild coordinates aligned with the torus:
//     A_r = 0,
//     A_{\theta'} = 0,
//     A_{\phi'} = pot_amp * \max(p_{gas} - pot_pgas_cutoff, 0)^pot_pgas_pow * r^2
//         * \sin(\theta') * \sin(\pi * pot_r_num * L(r; pot_r_min, pot_r_max))
//         * \sin(\pi * pot_theta_num * L(\theta'; pot_theta_min, \pi - pot_theta_min)),
//     L(x; x_{min}, x_{max}) = \min(\max((x - x_{\min}) / (x_{\max} - x_{min}), 0), 1).

namespace {
void VectorPotential(Real x1, Real x2, Real x3, Real *p_a_1, Real *p_a_2, Real *p_a_3) {
  // Calculate spherical Kerr-Schild coordinates of cell
  Real r, th, ph;
  GetKerrSchildCoordinates(x1, x2, x3, &r, &th, &ph);
  Real sth = std::sin(th);
  Real cth = std::cos(th);
  Real sph = std::sin(ph);
  Real cph = std::cos(ph);

  // Calculate Boyer-Lindquist coordinates of cell
  Real r_bl, th_bl, ph_bl;
  GetBoyerLindquistCoordinates(x1, x2, x3, &r_bl, &th_bl, &ph_bl);
  Real sth_bl = std::sin(th_bl);
  Real cth_bl = std::cos(th_bl);
  Real sph_bl = std::sin(ph_bl);
  Real cph_bl = std::cos(ph_bl);

  // Account for tilt
  Real th_t, sth_t, sth_bl_t;
  if (tilt == 0.0) {
    th_t = th;
    sth_t = std::abs(sth);
    sth_bl_t = std::abs(sth_bl);
  } else {
    Real x = sth * cph;
    Real y = sth * sph;
    Real z = cth;
    Real x_t = cos_tilt * x - sin_tilt * z;
    Real y_t = y;
    Real z_t = sin_tilt * x + cos_tilt * z;
    th_t = std::acos(z_t);
    sth_t = std::hypot(x_t, y_t);
    Real x_bl = sth_bl * cph_bl;
    Real y_bl = sth_bl * sph_bl;
    Real z_bl = cth_bl;
    Real x_bl_t = cos_tilt * x_bl - sin_tilt * z_bl;
    Real y_bl_t = y_bl;
    sth_bl_t = std::hypot(x_bl_t, y_bl_t);
  }

  // Determine if we are in the magnetized torus
  if (r_bl < r_edge) {
    *p_a_1 = 0.0;
    *p_a_2 = 0.0;
    *p_a_3 = 0.0;
    return;
  }
  Real log_h = LogHAux(r_bl, sth_bl_t) - log_h_edge;  // (FM 3.6)
  if (log_h < 0.0) {
    *p_a_1 = 0.0;
    *p_a_2 = 0.0;
    *p_a_3 = 0.0;
    return;
  }
  if (field_config == MagneticFieldConfigs::loops && (r <= pot_r_min || r >= pot_r_max
        || th_t <= pot_theta_min || th_t >= PI - pot_theta_min)) {
    *p_a_1 = 0.0;
    *p_a_2 = 0.0;
    *p_a_3 = 0.0;
    return;
  }

  // Calculate thermodynamic variables
  Real pgas_over_rho = (gamma_adi - 1.0) / gamma_adi * std::expm1(log_h);
  Real rho = rho_amp * std::pow(pgas_over_rho, 1.0 / (gamma_adi - 1.0));
  Real pgas = pgas_over_rho * rho;

  // Density isocontour configuration
  Real a_r, a_th_t, a_ph_t;
  if (field_config == MagneticFieldConfigs::density) {
    Real rho_cut = std::max(rho - pot_rho_cutoff, static_cast<Real>(0.0));
    a_r = 0.0;
    a_th_t = 0.0;
    a_ph_t = std::pow(rho_cut, pot_rho_pow) * std::pow(r, pot_r_pow);

  // Multiple loops configuration
  } else if (field_config == MagneticFieldConfigs::loops) {
    Real pgas_cut = std::max(pgas - pot_pgas_cutoff, static_cast<Real>(0.0));
    Real arg_r = PI * pot_r_num * (r - pot_r_min) / (pot_r_max - pot_r_min);
    Real arg_th =
        PI * pot_theta_num * (th_t - pot_theta_min) / (PI - 2.0 * pot_theta_min);
    a_r = 0.0;
    a_th_t = 0.0;
    a_ph_t = std::pow(pgas_cut, pot_pgas_pow) * SQR(r) * sth_t * std::sin(arg_r)
        * std::sin(arg_th);

  // Unknown configuration
  } else {
    *p_a_1 = NAN;
    *p_a_2 = NAN;
    *p_a_3 = NAN;
    return;
  }

  // Account for tilt
  Real a_th, a_ph;
  if (tilt == 0.0) {
    a_th = a_th_t;
    a_ph = a_ph_t;
  } else {
    Real x = sth * cph;
    Real y = sth * sph;
    Real z = cth;
    Real x_t = cos_tilt * x - sin_tilt * z;
    Real y_t = y;
    // Real z_t = sin_tilt * x + cos_tilt * z;
    Real sth_tilt = std::hypot(x_t, y_t);
    Real dth_t_dth = (-sin_tilt * cth * cph + cos_tilt * sth) / sth_tilt;
    Real dth_t_dph = sin_tilt * sth * sph / sth_tilt;
    Real dph_t_dth = -sin_tilt * sph / SQR(sth_tilt);
    Real dph_t_dph = (-sin_tilt * cth * sth * cph + cos_tilt * SQR(sth)) / SQR(sth_tilt);
    a_th = dth_t_dth * a_th_t + dph_t_dth * a_ph_t;
    a_ph = dth_t_dph * a_th_t + dph_t_dph * a_ph_t;
  }

  // Account for normalization
  a_r *= pot_amp;
  a_th *= pot_amp;
  a_ph *= pot_amp;

  // Transform to preferred coordinates
  TransformCovariantFromKerrSchild(a_r, a_th, a_ph, x1, x2, x3, p_a_1, p_a_2, p_a_3);
  return;
}
} // namespace
