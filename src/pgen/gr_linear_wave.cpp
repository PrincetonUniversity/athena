//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gr_linear_wave.cpp
//  \brief Problem generator for linear waves in special and general relativity.

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // abs(), cbrt(), sin(), sqrt()
#include <cstdio>     // fopen(), fprintf(), freopen()
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // string

// Athena++ headers
#include "../athena.hpp"                   // macros, enums, FaceField
#include "../athena_arrays.hpp"            // AthenaArray
#include "../bvals/bvals.hpp"              // BoundaryValues
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../eos/eos.hpp"                  // EquationOfState
#include "../field/field.hpp"              // Field
#include "../globals.hpp"                  // Globals
#include "../hydro/hydro.hpp"              // Hydro
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"          // ParameterInput

// Configuration checking
#if not RELATIVISTIC_DYNAMICS
#error "This problem generator must be used with relativity"
#endif

// Declarations
static void GetMinkowskiCoordinates(Real x0, Real x1, Real x2, Real x3, Real *pt,
    Real *px, Real *py, Real *pz);
static void TransformVector(Real at, Real ax, Real ay, Real az, Real x, Real y, Real z,
    Real *pa0, Real *pa1, Real *pa2, Real *pa3);
static Real QuadraticRoot(Real a1, Real a0, bool greater_root);
static Real CubicRootReal(Real a2, Real a1, Real a0);
static void QuarticRoots(Real a3, Real a2, Real a1, Real a0, Real *px1, Real *px2,
    Real *px3, Real *px4);

// Global variables
Real amp;                     // amplitude of wave
bool compute_error;           // flag indicating L1 errors should be computed and saved
Real gamma_adi_red;           // reduced adiabatic index \Gamma/(\Gamma-1)
Real rho, pgas;               // thermodynamic quantities
Real vx, vy, vz;              // 3-velocity components
Real bx;                      // longitudinal magnetic field
Real u[4], b[4];              // contravariant quantities
Real delta_rho, delta_pgas;   // perturbations to thermodynamic quantities
Real delta_u[4], delta_b[4];  // perturbations to contravariant quantities
Real delta_v[4];              // perturbations to 3-velocity
Real lambda;                  // wavespeed
Real wavenumber;              // wavenumber
AthenaArray<Real> g, gi;      // metric and inverse
AthenaArray<Real> bcc;        // cell-centered initial magnetic fields
AthenaArray<Real> initial;    // initial conditions
AthenaArray<Real> volume;     // 1D array of volumes

//----------------------------------------------------------------------------------------
// Function for initializing global variables and allocating arrays
// Inputs:
//   pin: input parameters (unused)
// Outputs: (none)

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Read information regarding desired wave and check input
  int wave_flag = pin->GetInteger("problem", "wave_flag");
  if (wave_flag < 0 or wave_flag > NWAVE-1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Problem Generator\n"
        << "wave_flag=" << wave_flag << " must be between 0 and " << NWAVE-1 << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  amp = pin->GetReal("problem", "amp");
  compute_error = pin->GetOrAddBoolean("problem", "compute_error", false);

  // Get ratio of specific heats
  Real gamma_adi = pin->GetReal("hydro", "gamma");
  gamma_adi_red = gamma_adi / (gamma_adi - 1.0);

  // Read background state
  rho = pin->GetReal("problem", "rho");
  pgas = pin->GetReal("problem", "pgas");
  vx = pin->GetReal("problem", "vx");
  vy = pin->GetReal("problem", "vy");
  vz = pin->GetReal("problem", "vz");
  bx = 0.0;
  Real by = 0.0, bz = 0.0;
  if (MAGNETIC_FIELDS_ENABLED) {
    bx = pin->GetReal("problem", "Bx");
    by = pin->GetReal("problem", "By");
    bz = pin->GetReal("problem", "Bz");
  }

  // Calculate background 4-vectors
  Real v_sq = SQR(vx) + SQR(vy) + SQR(vz);
  u[0] = 1.0 / std::sqrt(1.0 - v_sq);
  u[1] = u[0]*vx;
  u[2] = u[0]*vy;
  u[3] = u[0]*vz;
  b[0] = bx*u[1] + by*u[2] + bz*u[3];
  b[1] = 1.0/u[0] * (bx + b[0]*u[1]);
  b[2] = 1.0/u[0] * (by + b[0]*u[2]);
  b[3] = 1.0/u[0] * (bz + b[0]*u[3]);

  // Calculate useful background scalars
  Real b_sq = -SQR(b[0]) + SQR(b[1]) + SQR(b[2]) + SQR(b[3]);
  Real wgas = rho + gamma_adi_red * pgas;
  Real wtot = wgas + b_sq;
  Real cs_sq = gamma_adi * pgas / wgas;
  Real cs = std::sqrt(cs_sq);

  // Calculate desired perturbation
  if (MAGNETIC_FIELDS_ENABLED) {
    switch (wave_flag) {
      case 3: {  // entropy (A 46)
        lambda = vx;
        delta_rho = 1.0;
        delta_pgas = 0.0;
        for (int mu = 0; mu < 4; ++mu) {
          delta_u[mu] = 0.0;
          delta_b[mu] = 0.0;
        }
        break;
      }
      case 1: case 5: {  // Alfven (A 65)

        // Calculate wavespeed
        Real lambda_ap = (b[1] + std::sqrt(wtot) * u[1])
            / (b[0] + std::sqrt(wtot) * u[0]);            // (A 38)
        Real lambda_am = (b[1] - std::sqrt(wtot) * u[1])
            / (b[0] - std::sqrt(wtot) * u[0]);            // (A 38)
        Real sign = 1.0;
        if (lambda_ap > lambda_am) {  // \lambda_{a,\pm} = \lambda_a^\pm
          if (wave_flag == 1) {  // leftgoing
            sign = -1.0;
          }
        } else {  // lambda_{a,\pm} = \lambda_a^\mp
          if (wave_flag == 5) {  // rightgoing
            sign = -1.0;
          }
        }
        if (sign > 0) {  // want \lambda_{a,+}
          lambda = lambda_ap;
        } else {  // want \lambda_{a,-} instead
          lambda = lambda_am;
        }

        // Prepare auxiliary quantities
        Real alpha_1[4], alpha_2[4];
        alpha_1[0] = u[3];                                              // (A 58)
        alpha_1[1] = lambda * u[3];                                     // (A 58)
        alpha_1[2] = 0.0;                                               // (A 58)
        alpha_1[3] = u[0] - lambda * u[1];                              // (A 58)
        alpha_2[0] = -u[2];                                             // (A 59)
        alpha_2[1] = -lambda * u[2];                                    // (A 59)
        alpha_2[2] = lambda * u[1] - u[0];                              // (A 59)
        alpha_2[3] = 0.0;                                               // (A 59)
        Real g_1 = 1.0/u[0] * (by + lambda*vy / (1.0-lambda*vx) * bx);  // (A 60)
        Real g_2 = 1.0/u[0] * (bz + lambda*vz / (1.0-lambda*vx) * bx);  // (A 61)
        Real f_1, f_2;
        if (g_1 == 0.0 and g_2 == 0.0) {
          f_1 = f_2 = ONE_OVER_SQRT2;  // (A 67)
        } else {
          f_1 = g_1 / std::sqrt(SQR(g_1) + SQR(g_2));  // (A 66)
          f_2 = g_2 / std::sqrt(SQR(g_1) + SQR(g_2));  // (A 66)
        }

        // Set perturbation
        delta_rho = 0.0;
        delta_pgas = 0.0;
        for (int mu = 0; mu < 4; ++mu) {
          delta_u[mu] = f_1 * alpha_1[mu] + f_2 * alpha_2[mu];
          delta_b[mu] = -sign * std::sqrt(wtot) * delta_u[mu];
        }
        break;
      }
      default: {  // magnetosonic (A 71)

        // Calculate wavespeed
        Real factor_a = wgas * (1.0/cs_sq - 1.0);
        Real factor_b = -(wgas + b_sq/cs_sq);
        Real gamma_2 = SQR(u[0]);
        Real gamma_4 = SQR(gamma_2);
        Real coeff_4 = factor_a * gamma_4
                     - factor_b * gamma_2
                     - SQR(b[0]);
        Real coeff_3 = -factor_a * 4.0 * gamma_4 * vx
                     + factor_b * 2.0 * gamma_2 * vx
                     + 2.0 * b[0] * b[1];
        Real coeff_2 = factor_a * 6.0 * gamma_4 * SQR(vx)
                     + factor_b * gamma_2 * (1.0-SQR(vx))
                     + SQR(b[0]) - SQR(b[1]);
        Real coeff_1 = -factor_a * 4.0 * gamma_4 * vx*SQR(vx)
                     - factor_b * 2.0 * gamma_2 * vx
                     - 2.0 * b[0] * b[1];
        Real coeff_0 = factor_a * gamma_4 * SQR(SQR(vx))
                     + factor_b * gamma_2 * SQR(vx)
                     + SQR(b[1]);
        Real lambda_fl, lambda_sl, lambda_sr, lambda_fr;
        QuarticRoots(coeff_3/coeff_4, coeff_2/coeff_4, coeff_1/coeff_4, coeff_0/coeff_4,
            &lambda_fl, &lambda_sl, &lambda_sr, &lambda_fr);
        Real lambda_other_ms;
        if (wave_flag == 0) {
          lambda = lambda_fl;
          lambda_other_ms = lambda_sl;
        }
        if (wave_flag == 2) {
          lambda = lambda_sl;
          lambda_other_ms = lambda_fl;
        }
        if (wave_flag == 4) {
          lambda = lambda_sr;
          lambda_other_ms = lambda_fr;
        }
        if (wave_flag == 6) {
          lambda = lambda_fr;
          lambda_other_ms = lambda_sr;
        }

        // Determine which sign to use
        Real lambda_ap = (b[1] + std::sqrt(wtot) * u[1])
            / (b[0] + std::sqrt(wtot) * u[0]);            // (A 38)
        Real lambda_am = (b[1] - std::sqrt(wtot) * u[1])
            / (b[0] - std::sqrt(wtot) * u[0]);            // (A 38)
        Real lambda_a = lambda_ap;
        Real sign = 1.0;
        if (lambda_ap > lambda_am) {  // \lambda_{a,\pm} = \lambda_a^\pm
          if (wave_flag < 3) {  // leftgoing
            lambda_a = lambda_am;
            sign = -1.0;
          }
        } else {  // lambda_{a,\pm} = \lambda_a^\mp
          if (wave_flag > 3) {  // rightgoing
            lambda_a = lambda_am;
            sign = -1.0;
          }
        }

        // Prepare auxiliary quantities
        Real a = u[0] * (vx - lambda);                                       // (A 39)
        Real g = 1.0 - SQR(lambda);                                          // (A 41)
        Real b_over_a = -sign * std::sqrt(-factor_b - factor_a * SQR(a)/g);  // (A 68)
        Real alpha_1[4], alpha_2[4];
        alpha_1[0] = u[3];                                                   // (A 58)
        alpha_1[1] = lambda * u[3];                                          // (A 58)
        alpha_1[2] = 0.0;                                                    // (A 58)
        alpha_1[3] = u[0] - lambda * u[1];                                   // (A 58)
        alpha_2[0] = -u[2];                                                  // (A 59)
        alpha_2[1] = -lambda * u[2];                                         // (A 59)
        alpha_2[2] = lambda * u[1] - u[0];                                   // (A 59)
        alpha_2[3] = 0.0;                                                    // (A 59)
        Real alpha_11 = -SQR(alpha_1[0]);
        Real alpha_12 = -alpha_1[0] * alpha_2[0];
        Real alpha_22 = -SQR(alpha_2[0]);
        for (int i = 1; i < 4; ++i) {
          alpha_11 += SQR(alpha_1[i]);
          alpha_12 += alpha_1[i] * alpha_2[i];
          alpha_22 += SQR(alpha_2[i]);
        }
        Real g_1 = 1.0/u[0] * (by + lambda*vy / (1.0-lambda*vx) * bx);       // (A 60)
        Real g_2 = 1.0/u[0] * (bz + lambda*vz / (1.0-lambda*vx) * bx);       // (A 61)
        Real c_1 = (g_1*alpha_12 + g_2*alpha_22)
            / (alpha_11*alpha_22 - SQR(alpha_12)) * u[0] * (1.0-lambda*vx);  // (A 63)
        Real c_2 = -(g_1*alpha_11 + g_2*alpha_12)
            / (alpha_11*alpha_22 - SQR(alpha_12)) * u[0] * (1.0-lambda*vx);  // (A 63)
        Real b_t[4];
        for (int mu = 0; mu < 4; ++mu) {
          b_t[mu] = c_1 * alpha_1[mu] + c_2 * alpha_2[mu];  // (A 62)
        }
        Real f_1, f_2;
        if (g_1 == 0.0 and g_2 == 0.0) {
          f_1 = f_2 = ONE_OVER_SQRT2;  // (A 67)
        } else {
          f_1 = g_1 / std::sqrt(SQR(g_1) + SQR(g_2));  // (A 66)
          f_2 = g_2 / std::sqrt(SQR(g_1) + SQR(g_2));  // (A 66)
        }
        Real phi_plus_a_u[4];
        for (int mu = 0; mu < 4; ++mu) {
          phi_plus_a_u[mu] = a * u[mu];
        }
        phi_plus_a_u[0] += lambda;
        phi_plus_a_u[1] += 1.0;

        // Set perturbation
        if (std::abs(lambda-lambda_a)                 // using closer magnetosonic wave...
            <= std::abs(lambda_other_ms-lambda_a)) {  // ...to the associated Alfven wave
          Real b_t_normalized[4];
          Real denom = std::sqrt((alpha_11*alpha_22 - SQR(alpha_12))
              * (SQR(f_1)*alpha_11 + 2.0*f_1*f_2*alpha_12 + SQR(f_2)*alpha_22));
          for (int mu = 0; mu < 4; ++mu) {
            b_t_normalized[mu] = ((f_1*alpha_12+f_2*alpha_22) * alpha_1[mu]
                - (f_1*alpha_11+f_2*alpha_12) * alpha_2[mu]) / denom;        // (A 75)
          }
          Real b_t_norm = -SQR(b_t[0]);
          for (int i = 1; i < 4; ++i) {
            b_t_norm += SQR(b_t[i]);
          }
          b_t_norm = std::sqrt(b_t_norm);
          denom = SQR(a) - (g+SQR(a)) * cs_sq;
          if (denom == 0.0) {
            delta_pgas = 0.0;
          } else {
            delta_pgas = -(g+SQR(a)) * cs_sq / denom * b_t_norm;  // (A 74)
          }
          delta_rho = rho / (gamma_adi*pgas) * delta_pgas;
          for (int mu = 0; mu < 4; ++mu) {
            delta_u[mu] =
                -a*delta_pgas / (wgas*cs_sq*(g+SQR(a))) * phi_plus_a_u[mu]
                - b_over_a / wgas * b_t_normalized[mu];                     // (A 72)
            delta_b[mu] = -b_over_a * delta_pgas/wgas * u[mu]
                - (1.0+SQR(a)/g) * b_t_normalized[mu];                      // (A 73)
          }
        } else {  // using more distant magnetosonic wave
          delta_pgas = -1.0;                                // (A 78)
          delta_rho = rho / (gamma_adi*pgas) * delta_pgas;
          Real b_t_reduced[4] = {0.0};                      // (A 79)
          Real denom = wgas * SQR(a) - b_sq * g;
          if (denom != 0.0) {
            for (int mu = 0; mu < 4; ++mu) {
              b_t_reduced[mu] = b_t[mu] / denom;
            }
          }
          for (int mu = 0; mu < 4; ++mu) {
            delta_u[mu] = a / (wgas*cs_sq*(g+SQR(a))) * phi_plus_a_u[mu]
                - b_over_a * g/wgas * b_t_reduced[mu];                    // (A 76)
            delta_b[mu] = b_over_a / wgas * u[mu]
                - (1.0+SQR(a)/g) * g * b_t_reduced[mu];                   // (A 77)
          }
        }
      }
    }
  } else {  // hydro

    // Calculate perturbation in 4-velocity components (Q of FK)
    switch (wave_flag) {
      case 1:  // entropy 1/3
        lambda = vx;
        delta_rho = 1.0;
        delta_pgas = 0.0;
        delta_u[1] = delta_u[2] = delta_u[3] = 0.0;
        break;
      case 2:  // entropy 2/3
        lambda = vx;
        delta_rho = 0.0;
        delta_pgas = 0.0;
        delta_u[1] = vx * vy / (1.0 - SQR(vx));
        delta_u[2] = 1.0;
        delta_u[3] = 0.0;
        break;
      case 3:  // entropy 3/3
        lambda = vx;
        delta_rho = 0.0;
        delta_pgas = 0.0;
        delta_u[1] = vx * vz / (1.0 - SQR(vx));
        delta_u[2] = 0.0;
        delta_u[3] = 1.0;
        break;
      default:  // sound
        Real delta = SQR(u[0]) * (1.0-cs_sq) + cs_sq;
        Real v_minus_lambda_a = vx * cs_sq;
        Real v_minus_lambda_b =
            cs * std::sqrt(SQR(u[0]) * (1.0-cs_sq) * (1.0-SQR(vx)) + cs_sq);
        Real v_minus_lambda;
        if (wave_flag == 0) {  // leftgoing
          v_minus_lambda = (v_minus_lambda_a + v_minus_lambda_b) / delta;  // (FK A1)
        } else {  // rightgoing
          v_minus_lambda = (v_minus_lambda_a - v_minus_lambda_b) / delta;  // (FK A1)
        }
        lambda = vx - v_minus_lambda;
        delta_rho = rho;
        delta_pgas = wgas * cs_sq;
        delta_u[1] = -cs_sq * u[1] - cs_sq / u[0] / v_minus_lambda;
        delta_u[2] = -cs_sq * u[2];
        delta_u[3] = -cs_sq * u[3];
    }

    // Calculate perturbation in 3-velocity components (P of FK)
    delta_v[1] = (1.0-SQR(vx)) * delta_u[1] - vx*vy * delta_u[2] - vx*vz * delta_u[3];
    delta_v[2] = -vx*vy * delta_u[1] + (1.0-SQR(vy)) * delta_u[2] - vy*vz * delta_u[3];
    delta_v[3] = -vx*vz * delta_u[1] - vy*vz * delta_u[2] + (1.0-SQR(vz)) * delta_u[3];
    for (int i = 1; i < 4; ++i) {
      delta_v[i] /= u[0];
    }
  }

  // Renormalize perturbation to unit L^2 norm
  if (MAGNETIC_FIELDS_ENABLED) {
    Real perturbation_size = SQR(delta_rho) + SQR(delta_pgas);
    for (int mu = 0; mu < 4; ++mu) {
      perturbation_size += SQR(delta_u[mu]) + SQR(delta_b[mu]);
    }
    perturbation_size = std::sqrt(perturbation_size);
    delta_rho /= perturbation_size;
    delta_pgas /= perturbation_size;
    for (int mu = 0; mu < 4; ++mu) {
      delta_u[mu] /= perturbation_size;
      delta_b[mu] /= perturbation_size;
    }
  } else {  // hydro
    Real perturbation_size = SQR(delta_rho) + SQR(delta_pgas);
    for (int i = 1; i < 4; ++i) {
      perturbation_size += SQR(delta_v[i]);
    }
    perturbation_size = std::sqrt(perturbation_size);
    delta_rho /= perturbation_size;
    delta_pgas /= perturbation_size;
    for (int i = 1; i < 4; ++i) {
      delta_v[i] /= perturbation_size;
    }
  }

  // Prepare arrays to hold metric
  if (GENERAL_RELATIVITY) {
    int ncells1 = mesh_size.nx1/nrbx1 + 2*NGHOST;
    g.NewAthenaArray(NMETRIC, ncells1);
    gi.NewAthenaArray(NMETRIC, ncells1);
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for evaluating errors and freeing arrays
// Inputs:
//   pin: parameters
// Outputs: (none)

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  // Free metric
  if (GENERAL_RELATIVITY) {
    g.DeleteAthenaArray();
    gi.DeleteAthenaArray();
  }

  // Calculate L1 error against initial conditions
  if (compute_error) {

    // Prepare error calculation variables
    Real errors[(NHYDRO+NFIELD)+1];
    for (int n = 0; n < (NHYDRO+NFIELD)+1; ++n) {
      errors[n] = 0.0;
    }

    // Go through blocks to calculate errors
    MeshBlock *pmb = pblock;
    while (pmb != NULL) {
      for (int k = pmb->ks; k <= pmb->ke; ++k) {
        for (int j = pmb->js; j <= pmb->je; ++j) {
          pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, volume);
          for (int i = pmb->is; i <= pmb->ie; ++i) {
            for (int n = 0; n < NHYDRO; ++n) {
              errors[n] += std::abs(pmb->phydro->u(n,k,j,i) - initial(pmb->lid,n,k,j,i))
                  * volume(i);
            }
            if (MAGNETIC_FIELDS_ENABLED) {
              for (int n = IB1; n <= IB3; ++n) {
                errors[NHYDRO+n] += std::abs(pmb->pfield->bcc(n,k,j,i)
                    - initial(pmb->lid,NHYDRO+n,k,j,i)) * volume(i);
              }
            }
            errors[(NHYDRO+NFIELD)] += volume(i);
          }
        }
      }
      pmb = pmb->next;
    }

    // Reduce errors across ranks
    #ifdef MPI_PARALLEL
    {
      if (Globals::my_rank == 0) {
        MPI_Reduce(MPI_IN_PLACE, errors, (NHYDRO+NFIELD)+1, MPI_ATHENA_REAL, MPI_SUM, 0,
            MPI_COMM_WORLD);
      } else {
        MPI_Reduce(errors, 0, (NHYDRO+NFIELD)+1, MPI_ATHENA_REAL, MPI_SUM, 0,
                   MPI_COMM_WORLD);
      }
    }
    #endif

    // Write errors to file if root
    if (Globals::my_rank == 0) {

      // Divide volume-weighted errors by total volume
      for (int n = 0; n < (NHYDRO+NFIELD); ++n) {
        errors[n] /= errors[(NHYDRO+NFIELD)];
      }

      // Calculate RMS of volume-averaged errors
      Real total_error = 0.0;
      for (int n = 0; n < (NHYDRO+NFIELD); ++n) {
        total_error += SQR(errors[n]);
      }
      total_error = std::sqrt(total_error/(NHYDRO+NFIELD));

      // Prepare output file
      std::string filename;
      filename.assign("linearwave-errors.dat");
      std::stringstream msg;
      FILE *pfile;

      // Open file
      pfile = fopen(filename.c_str(), "r");
      if (pfile != NULL) {  // file exists
        pfile = freopen(filename.c_str(), "a", pfile);
        if (pfile == NULL) {
          msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]\n"
              << "Error output file could not be opened" << std::endl;
          throw std::runtime_error(msg.str().c_str());
        }
      } else {  // file does not exist
        pfile = fopen(filename.c_str(), "w");
        if (pfile == NULL) {
          msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]\n"
              << "Error output file could not be opened" << std::endl;
          throw std::runtime_error(msg.str().c_str());
        }
        fprintf(pfile, "# Nx1  Nx2  Nx3  Ncycle  RMS-Error  D  E  M1  M2  M3");
        if (MAGNETIC_FIELDS_ENABLED) {
          fprintf(pfile, "  B1c  B2c  B3c");
        }
        fprintf(pfile, "\n");
      }

      // Write errors
      fprintf(pfile, "%d  %d  %d  %d  %e", mesh_size.nx1, mesh_size.nx2, mesh_size.nx3,
          ncycle, total_error);
      fprintf(pfile, "  %e  %e  %e  %e  %e", errors[IDN], errors[IEN], errors[IM1],
          errors[IM2], errors[IM3]);
      if (MAGNETIC_FIELDS_ENABLED) {
        fprintf(pfile,"  %e  %e  %e", errors[NHYDRO+IB1], errors[NHYDRO+IB2],
            errors[NHYDRO+IB3]);
      }
      fprintf(pfile, "\n");

      // Close file
      fclose(pfile);
    }

    // Free initial conditions arrays
    if (MAGNETIC_FIELDS_ENABLED) {
      bcc.DeleteAthenaArray();
    }
    initial.DeleteAthenaArray();
    volume.DeleteAthenaArray();
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for setting initial conditions
// Inputs:
//   pin: parameters
// Outputs: (none)
// Notes:
//   initializes linear wave with sinusoidal variation
//     sets both primitive and conserved variables
//   references Anton et al. 2010, ApJS 188 1 (A, MHD)
//              Falle & Komissarov 1996, MNRAS 278 586 (FK, hydro)

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

  // Calculate wavenumber such that wave has single period over domain
  Real x1_min = pmy_mesh->mesh_size.x1min;
  Real x1_max = pmy_mesh->mesh_size.x1max;
  Real x2_min = pmy_mesh->mesh_size.x2min;
  Real x3_min = pmy_mesh->mesh_size.x3min;
  Real arg_min, arg_max;
  #if GENERAL_RELATIVITY
  {
    Real t_left, x_left, y_left, z_left;
    Real t_right, x_right, y_right, z_right;
    GetMinkowskiCoordinates(0.0, x1_min, x2_min, x3_min, &t_left, &x_left, &y_left,
        &z_left);
    GetMinkowskiCoordinates(0.0, x1_max, x2_min, x3_min, &t_right, &x_right, &y_right,
        &z_right);
    arg_min = x_left - lambda * t_left;
    arg_max = x_right - lambda * t_right;
  }
  #else  // SR
  {
    arg_min = x1_min;
    arg_max = x1_max;
  }
  #endif  // GENERAL_RELATIVITY
  wavenumber = 2.0*PI / (arg_max - arg_min);

  // Initialize hydro variables
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      #if GENERAL_RELATIVITY
      {
        pcoord->CellMetric(k, j, il, iu, g, gi);
      }
      #endif  // GENERAL_RELATIVITY
      for (int i = il; i <= iu; ++i) {

        // Find location of cell in spacetime
        Real t, x, y, z;
        #if GENERAL_RELATIVITY
        {
          GetMinkowskiCoordinates(0.0, pcoord->x1v(i), pcoord->x2v(j), pcoord->x3v(k), &t,
              &x, &y, &z);
        }
        #else  // SR
        {
          t = 0.0;
          x = pcoord->x1v(i);
        }
        #endif  // GENERAL_RELATIVITY

        // Calculate scalar perturbations
        Real local_amp = amp * std::sin(wavenumber * (x - lambda * t));
        Real rho_local = rho + local_amp * delta_rho;
        Real pgas_local = pgas + local_amp * delta_pgas;

        // Calculate vector perturbations
        Real u_mink[4];
        Real b_mink[4] = {0.0};
        if (MAGNETIC_FIELDS_ENABLED) {
          for (int mu = 0; mu < 4; ++mu) {
            u_mink[mu] = u[mu] + local_amp * delta_u[mu];
            b_mink[mu] = b[mu] + local_amp * delta_b[mu];
          }
        } else {  // hydro
          Real vx_mink = vx + local_amp * delta_v[1];
          Real vy_mink = vy + local_amp * delta_v[2];
          Real vz_mink = vz + local_amp * delta_v[3];
          u_mink[0] = 1.0 / std::sqrt(1.0 - SQR(vx_mink) - SQR(vy_mink) - SQR(vz_mink));
          u_mink[1] = u_mink[0] * vx_mink;
          u_mink[2] = u_mink[0] * vy_mink;
          u_mink[3] = u_mink[0] * vz_mink;
        }

        // Transform vector perturbations
        Real u_local[4], b_local[4], u_local_low[4], b_local_low[4];
        #if GENERAL_RELATIVITY
        {
          TransformVector(u_mink[0], u_mink[1], u_mink[2], u_mink[3], x, y, z,
              &u_local[0], &u_local[1], &u_local[2], &u_local[3]);
          TransformVector(b_mink[0], b_mink[1], b_mink[2], b_mink[3], x, y, z,
              &b_local[0], &b_local[1], &b_local[2], &b_local[3]);
          pcoord->LowerVectorCell(u_local[0], u_local[1], u_local[2], u_local[3], k, j, i,
              &u_local_low[0], &u_local_low[1], &u_local_low[2], &u_local_low[3]);
          pcoord->LowerVectorCell(b_local[0], b_local[1], b_local[2], b_local[3], k, j, i,
              &b_local_low[0], &b_local_low[1], &b_local_low[2], &b_local_low[3]);
        }
        #else  // SR
        {
          for (int mu = 0; mu < 4; ++mu) {
            u_local[mu] = u_mink[mu];
            b_local[mu] = b_mink[mu];
            u_local_low[mu] = (mu == 0 ? -1.0 : 1.0) * u_local[mu];
            b_local_low[mu] = (mu == 0 ? -1.0 : 1.0) * b_local[mu];
          }
        }
        #endif  // GENERAL_RELATIVITY

        // Calculate useful local scalars
        Real b_sq_local = 0.0;
        for (int mu = 0; mu < 4; ++mu) {
          b_sq_local += b_local[mu] * b_local_low[mu];
        }
        Real wtot_local = rho_local + gamma_adi_red * pgas_local + b_sq_local;
        Real ptot_local = pgas_local + 0.5*b_sq_local;

        // Set primitive hydro variables
        phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) = rho_local;
        phydro->w(IPR,k,j,i) = phydro->w1(IPR,k,j,i) = pgas_local;
        if (GENERAL_RELATIVITY) {
          Real uu1 = u_local[1] - gi(I01,i)/gi(I00,i) * u_local[0];
          Real uu2 = u_local[2] - gi(I02,i)/gi(I00,i) * u_local[0];
          Real uu3 = u_local[3] - gi(I03,i)/gi(I00,i) * u_local[0];
          phydro->w(IVX,k,j,i) = phydro->w1(IVX,k,j,i) = uu1;
          phydro->w(IVY,k,j,i) = phydro->w1(IVY,k,j,i) = uu2;
          phydro->w(IVZ,k,j,i) = phydro->w1(IVZ,k,j,i) = uu3;
        } else {  // SR
          phydro->w(IVX,k,j,i) = phydro->w1(IVX,k,j,i) = u_local[1] / u_local[0];
          phydro->w(IVY,k,j,i) = phydro->w1(IVY,k,j,i) = u_local[2] / u_local[0];
          phydro->w(IVZ,k,j,i) = phydro->w1(IVZ,k,j,i) = u_local[3] / u_local[0];
        }

        // Set conserved hydro variables
        phydro->u(IDN,k,j,i) = u_local[0] * rho_local;
        if (GENERAL_RELATIVITY) {
          phydro->u(IEN,k,j,i) = wtot_local*u_local[0]*u_local_low[0]
              - b_local[0]*b_local_low[0] + ptot_local;
          phydro->u(IM1,k,j,i) = wtot_local*u_local[0]*u_local_low[1]
              - b_local[0]*b_local_low[1];
          phydro->u(IM2,k,j,i) = wtot_local*u_local[0]*u_local_low[2]
              - b_local[0]*b_local_low[2];
          phydro->u(IM3,k,j,i) = wtot_local*u_local[0]*u_local_low[3]
              - b_local[0]*b_local_low[3];
        } else {  // SR
          phydro->u(IEN,k,j,i) = wtot_local*u_local[0]*u_local[0] - b_local[0]*b_local[0]
            - ptot_local;
          phydro->u(IM1,k,j,i) = wtot_local*u_local[0]*u_local[1] - b_local[0]*b_local[1];
          phydro->u(IM2,k,j,i) = wtot_local*u_local[0]*u_local[2] - b_local[0]*b_local[2];
          phydro->u(IM3,k,j,i) = wtot_local*u_local[0]*u_local[3] - b_local[0]*b_local[3];
        }
      }
    }
  }

  // Initialize magnetic fields
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k = kl; k <= ku+1; ++k) {
      for (int j = jl; j <= ju+1; ++j) {
        for (int i = il; i <= iu+1; ++i) {
          #if GENERAL_RELATIVITY
          {
            // Set B^1 if needed
            if (j != ju+1 and k != ku+1) {
              Real t, x, y, z;
              GetMinkowskiCoordinates(0.0, pcoord->x1f(i), pcoord->x2v(j), pcoord->x3v(k),
                  &t, &x, &y, &z);
              Real local_amp = amp * std::sin(wavenumber * (x - lambda * t));
              Real u_mink[4], b_mink[4];
              for (int mu = 0; mu < 4; ++mu) {
                u_mink[mu] = u[mu] + local_amp * delta_u[mu];
                b_mink[mu] = b[mu] + local_amp * delta_b[mu];
              }
              Real u_local[4], b_local[4];
              TransformVector(u_mink[0], u_mink[1], u_mink[2], u_mink[3], x, y, z,
                  &u_local[0], &u_local[1], &u_local[2], &u_local[3]);
              TransformVector(b_mink[0], b_mink[1], b_mink[2], b_mink[3], x, y, z,
                  &b_local[0], &b_local[1], &b_local[2], &b_local[3]);
              pfield->b.x1f(k,j,i) = b_local[1] * u_local[0] - b_local[0] * u_local[1];
            }

            // Set B^2 if needed
            if (i != iu+1 and k != ku+1) {
              Real t, x, y, z;
              GetMinkowskiCoordinates(0.0, pcoord->x1v(i), pcoord->x2f(j), pcoord->x3v(k),
                  &t, &x, &y, &z);
              Real local_amp = amp * std::sin(wavenumber * (x - lambda * t));
              Real u_mink[4], b_mink[4];
              for (int mu = 0; mu < 4; ++mu) {
                u_mink[mu] = u[mu] + local_amp * delta_u[mu];
                b_mink[mu] = b[mu] + local_amp * delta_b[mu];
              }
              Real u_local[4], b_local[4];
              TransformVector(u_mink[0], u_mink[1], u_mink[2], u_mink[3], x, y, z,
                  &u_local[0], &u_local[1], &u_local[2], &u_local[3]);
              TransformVector(b_mink[0], b_mink[1], b_mink[2], b_mink[3], x, y, z,
                  &b_local[0], &b_local[1], &b_local[2], &b_local[3]);
              pfield->b.x2f(k,j,i) = b_local[2] * u_local[0] - b_local[0] * u_local[2];
            }

            // Set B^3 if needed
            if (i != iu+1 and j != ju+1) {
              Real t, x, y, z;
              GetMinkowskiCoordinates(0.0, pcoord->x1v(i), pcoord->x2v(j), pcoord->x3f(k),
                  &t, &x, &y, &z);
              Real local_amp = amp * std::sin(wavenumber * (x - lambda * t));
              Real u_mink[4], b_mink[4];
              for (int mu = 0; mu < 4; ++mu) {
                u_mink[mu] = u[mu] + local_amp * delta_u[mu];
                b_mink[mu] = b[mu] + local_amp * delta_b[mu];
              }
              Real u_local[4], b_local[4];
              TransformVector(u_mink[0], u_mink[1], u_mink[2], u_mink[3], x, y, z,
                  &u_local[0], &u_local[1], &u_local[2], &u_local[3]);
              TransformVector(b_mink[0], b_mink[1], b_mink[2], b_mink[3], x, y, z,
                  &b_local[0], &b_local[1], &b_local[2], &b_local[3]);
              pfield->b.x3f(k,j,i) = b_local[3] * u_local[0] - b_local[0] * u_local[3];
            }
          }
          #else  // SR
          {
            Real local_amp = amp * std::sin(wavenumber * pcoord->x1v(i));
            Real u_local[4], b_local[4];
            for (int mu = 0; mu < 4; ++mu) {
              u_local[mu] = u[mu] + local_amp * delta_u[mu];
              b_local[mu] = b[mu] + local_amp * delta_b[mu];
            }
            Real by_local = b_local[2]*u_local[0] - b_local[0]*u_local[2];
            Real bz_local = b_local[3]*u_local[0] - b_local[0]*u_local[3];
            if (j != ju+1 and k != ku+1) {
              pfield->b.x1f(k,j,i) = bx;
            }
            if (i != iu+1 and k != ku+1) {
              pfield->b.x2f(k,j,i) = by_local;
            }
            if (i != iu+1 and j != ju+1) {
              pfield->b.x3f(k,j,i) = bz_local;
            }
          }
          #endif  // GENERAL_RELATIVITY
        }
      }
    }
  }

  // Prepare arrays for comparing to initial conditions (only once per Mesh)
  if (compute_error and lid == 0) {
    int num_blocks = pmy_mesh->GetNumMeshBlocksThisRank(Globals::my_rank);
    int nx1 = block_size.nx1;
    int nx2 = block_size.nx2;
    int nx3 = block_size.nx3;
    if (MAGNETIC_FIELDS_ENABLED) {
      bcc.NewAthenaArray(NFIELD, nx3+NGHOST, nx2+NGHOST, nx1+NGHOST);
    }
    initial.NewAthenaArray(num_blocks, (NHYDRO+NFIELD), nx3+NGHOST, nx2+NGHOST,
                           nx1+NGHOST);
    volume.NewAthenaArray(nx1+NGHOST);
  }

  // Record initial conditions
  if (compute_error) {
    for (int n = 0; n < NHYDRO; ++n) {
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
          for (int i = is; i <= ie; ++i) {
            initial(lid,n,k,j,i) = phydro->u(n,k,j,i);
          }
        }
      }
    }
    if (MAGNETIC_FIELDS_ENABLED) {
      pfield->CalculateCellCenteredField(pfield->b, bcc, pcoord, is, ie, js, je, ks, ke);
      for (int n = IB1; n <= IB3; ++n) {
        for (int k = ks; k <= ke; ++k) {
          for (int j = js; j <= je; ++j) {
            for (int i = is; i <= ie; ++i) {
              initial(lid,NHYDRO+n,k,j,i) = bcc(n,k,j,i);
            }
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for returning corresponding Minkowski coordinates of point
// Inputs:
//   x0,x1,x2,x3: global coordinates to be converted
// Outputs:
//   pt,px,py,pz: variables pointed to set to Minkowski coordinates
// Notes:
//   conversion is trivial
//   useful to have if other coordinate systems for Minkowski space are developed

static void GetMinkowskiCoordinates(Real x0, Real x1, Real x2, Real x3, Real *pt,
    Real *px, Real *py, Real *pz) {
  if (COORDINATE_SYSTEM == "minkowski") {
    *pt = x0;
    *px = x1;
    *py = x2;
    *pz = x3;
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for transforming 4-vector from Minkowski to desired coordinates
// Inputs:
//   at,ax,ay,az: upper 4-vector components in Minkowski coordinates
//   x,y,z: Minkowski coordinates of point
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to upper 4-vector components in desired coordinates
// Notes:
//   conversion is trivial
//   useful to have if other coordinate systems for Minkowski space are developed

static void TransformVector(Real at, Real ax, Real ay, Real az, Real x, Real y, Real z,
    Real *pa0, Real *pa1, Real *pa2, Real *pa3) {
  if (COORDINATE_SYSTEM == "minkowski") {
    *pa0 = at;
    *pa1 = ax;
    *pa2 = ay;
    *pa3 = az;
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for finding root of monic quadratic equation
// Inputs:
//   a1: linear coefficient
//   a0: constant coefficient
//   greater_root: flag indicating that larger root is to be returned
//     "larger" does not mean absolute value
// Outputs:
//   returned value: desired root
// Notes:
//   same function as in adiabatic_mhd_sr.cpp, adiabatic_mhd_gr.cpp, and hlld_rel.cpp
//   solves x^2 + a_1 x + a_0 = 0 for x
//   returns abscissa of vertex if there are no real roots
//   follows advice in Numerical Recipes, 3rd ed. (5.6) for avoiding large cancellations

static Real QuadraticRoot(Real a1, Real a0, bool greater_root) {
  if (a1*a1 < 4.0*a0) {  // no real roots
    return -a1/2.0;
  }
  if (greater_root) {
    if (a1 >= 0.0) {
      return -2.0*a0 / (a1 + std::sqrt(a1*a1 - 4.0*a0));
    } else {
      return (-a1 + std::sqrt(a1*a1 - 4.0*a0)) / 2.0;
    }
  } else {
    if (a1 >= 0.0) {
      return (-a1 - std::sqrt(a1*a1 - 4.0*a0)) / 2.0;
    } else {
      return -2.0*a0 / (a1 - std::sqrt(a1*a1 - 4.0*a0));
    }
  }
}

//----------------------------------------------------------------------------------------
// Function for finding real root of monic cubic equation
// Inputs:
//   a2: quadratic coefficient
//   a1: linear coefficient
//   a0: constant coefficient
// Outputs:
//   returned value: a real root
// Notes:
//   solves x^3 + a_2 x^2 + a_1 x + a_0 = 0 for x
//   same function as in adiabatic_mhd_sr.cpp and adiabatic_mhd_gr.cpp
//   references Numerical Recipes, 3rd ed. (NR)

static Real CubicRootReal(Real a2, Real a1, Real a0) {
  Real q = (a2*a2 - 3.0*a1) / 9.0;                       // (NR 5.6.10)
  Real r = (2.0*a2*a2*a2 - 9.0*a1*a2 + 27.0*a0) / 54.0;  // (NR 5.6.10)
  if (r*r - q*q*q < 0.0) {
    Real theta = acos(r/std::sqrt(q*q*q));                 // (NR 5.6.11)
    return -2.0 * std::sqrt(q) * cos(theta/3.0) - a2/3.0;  // (NR 5.6.12)
  } else {
    Real a = -copysign(1.0, r)
        * std::cbrt(std::abs(r) + std::sqrt(r*r - q*q*q));  // (NR 5.6.15)
    Real b = (a != 0.0) ? q/a : 0.0;                   // (NR 5.6.16)
    return a + b - a2/3.0;
  }
}

//----------------------------------------------------------------------------------------
// Function for finding extremal real roots of monic quartic equation
// Inputs:
//   a3: cubic coefficient
//   a2: quadratic coefficient
//   a1: linear coefficient
//   a0: constant coefficient
// Outputs:
//   px1: value set to least real root
//   px2: value set to second least real root
//   px3: value set to second greatest real root
//   px4: value set to greatest real root
// Notes:
//   solves x^4 + a3 x^3 + a2 x^2 + a1 x + a0 = 0 for x
//   uses following procedure:
//     1) eliminate cubic term y^4 + b2 y^2 + b1 y + b0
//     2) construct resolvent cubic z^3 + c2 z^2 + c1 z + c0
//     3) find real root z0 of cubic
//     4) construct quadratics:
//          y^2 + d1 y + d0
//          y^2 + e1 y + e0
//     5) find roots of quadratics
//   similar function to those in adiabatic_mhd_sr.cpp and adiabatic_mhd_gr.cpp

static void QuarticRoots(Real a3, Real a2, Real a1, Real a0, Real *px1, Real *px2,
    Real *px3, Real *px4) {
  // Step 1: Find reduced quartic coefficients
  Real b2 = a2 - 3.0/8.0*SQR(a3);
  Real b1 = a1 - 1.0/2.0*a2*a3 + 1.0/8.0*a3*SQR(a3);
  Real b0 = a0 - 1.0/4.0*a1*a3 + 1.0/16.0*a2*SQR(a3) - 3.0/256.0*SQR(SQR(a3));

  // Step 2: Find resolvent cubic coefficients
  Real c2 = -b2;
  Real c1 = -4.0*b0;
  Real c0 = 4.0*b0*b2 - SQR(b1);

  // Step 3: Solve cubic
  Real z0 = CubicRootReal(c2, c1, c0);

  // Step 4: Find quadratic coefficients
  Real d1 = (z0 - b2 > 0.0) ? std::sqrt(z0 - b2) : 0.0;
  Real e1 = -d1;
  Real d0, e0;
  if (b1 < 0) {
    d0 = z0/2.0 + std::sqrt(SQR(z0)/4.0 - b0);
    e0 = z0/2.0 - std::sqrt(SQR(z0)/4.0 - b0);
  } else {
    d0 = z0/2.0 - std::sqrt(SQR(z0)/4.0 - b0);
    e0 = z0/2.0 + std::sqrt(SQR(z0)/4.0 - b0);
  }

  // Step 5: Solve quadratics
  Real y1 = QuadraticRoot(d1, d0, false);
  Real y2 = QuadraticRoot(d1, d0, true);
  Real y3 = QuadraticRoot(e1, e0, false);
  Real y4 = QuadraticRoot(e1, e0, true);

  // Step 6: Set original quartic roots
  *px1 = std::min(y1, y3) - a3/4.0;
  Real mid_1 = std::max(y1, y3) - a3/4.0;
  *px4 = std::max(y2, y4) - a3/4.0;
  Real mid_2 = std::min(y2, y4) - a3/4.0;
  *px2 = std::min(mid_1, mid_2);
  *px3 = std::max(mid_1, mid_2);
  return;
}
