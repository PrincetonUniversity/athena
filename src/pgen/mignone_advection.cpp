//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mignone_advection.cpp
//! \brief 1D passive scalar advection tests in curvilinear coordinates from:
//! Mignone, A. (2014). High-order conservative reconstruction schemes for finite volume
//! methods in cylindrical and spherical coordinates. Journal of Computational Physics,
//! 270, 784-814.
//========================================================================================

// C++ headers
#include <algorithm>  // min, max
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"
#include "../utils/gl_quadrature.hpp"

// Parameters which define initial solution -- made global so that they can be shared
namespace {
constexpr int N_gl = 12;
constexpr Real d0 = 1.0;
constexpr bool use_gl_quadrature = true;
Real (*IntegrandInitial)(Real x1);
Real (*IntegrandFinal)(Real x1);
Real a_width, b_center, alpha;
Real iso_cs;
int m_coord;
Real t_final;
int iprob;

// pointwise analytic initial condition for Gaussian bell curve
Real InitialGaussianProfile(Real x1);
// pointwise analytic exact solution at any t >=0 for a linear x1 velocity profile
Real FinalGaussianProfile(Real x1);
Real InitialGaussianCylindricalIntegrand(Real x1);
Real FinalGaussianCylindricalIntegrand(Real x1);
Real InitialGaussianSphericalIntegrand(Real x1);
Real FinalGaussianSphericalIntegrand(Real x1);

Real InitialCosineProfile(Real x2);
Real FinalCosineProfile(Real x2);
Real InitialCosineSphericalIntegrand(Real x2);
Real FinalCosineSphericalIntegrand(Real x2);
} // namespace

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // read and initialize global parameters
  iso_cs = pin->GetReal("hydro", "iso_sound_speed");
  // initial condition; see Mignone (2014) section 5.1.1
  alpha = pin->GetOrAddReal("problem", "alpha", 1.0);
  a_width = pin->GetOrAddReal("problem", "a_width", 10.0);
  b_center = pin->GetOrAddReal("problem", "b_center", 0.0);
  t_final = pin->GetOrAddReal("time", "tlim", 1.0);
  iprob = pin->GetInteger("problem", "iprob");

  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    if (iprob == 1) {
      m_coord = 0;
      IntegrandInitial = &InitialGaussianProfile;
      IntegrandFinal = &FinalGaussianProfile;
    } else if (iprob == 2) {
      IntegrandInitial = &InitialCosineProfile;
      IntegrandFinal = &FinalCosineProfile;
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    if (iprob == 1) {
      m_coord = 1;
      IntegrandInitial = &InitialGaussianCylindricalIntegrand;
      IntegrandFinal = &FinalGaussianCylindricalIntegrand;
    } else if (iprob == 2) {
      // iprob=2 is for spherical-polar coordinates only
    }
  } else { // if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0)
    if (iprob == 1) {
      m_coord = 2;
      IntegrandInitial = &InitialGaussianSphericalIntegrand;
      IntegrandFinal = &FinalGaussianSphericalIntegrand;
    } else if (iprob == 2) {
      IntegrandInitial = &InitialCosineSphericalIntegrand;
      IntegrandFinal = &FinalCosineSphericalIntegrand;
    }
  }
  // Restrict to: 1) no-MHD 2) isothermal EOS only 3) hydro/active=background
  // TODO(felker): add explicit user-safety checks for these conditions
  return;
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  if (!pin->GetOrAddBoolean("problem", "compute_error", false)) return;

  if (NSCALARS > 0) { // TODO(felker): error-out at compile time if NSCALARS < 1
    // Initialize errors to zero
    // (temp. workaround for zero-length array prohibition from ISO C++; replace with
    // dynamically-sized arrays)
    Real l1_err[(NSCALARS > 0 ? NSCALARS : 1)]{},
        max_err[(NSCALARS > 0 ? NSCALARS : 1)]{};
    Real total_vol = 0.0;
    MeshBlock *pmb = my_blocks(0);
    AthenaArray<Real> vol(pmb->ncells1);
    // recalculate initial condition from ProblemGenerator on final Mesh configuration:
    // (may have changed due to AMR)
    for (int b=0; b<nblocal; ++b) {
      pmb = my_blocks(b);
      int il = pmb->is, iu = pmb->ie, jl = pmb->js, ju = pmb->je,
          kl = pmb->ks, ku = pmb->ke;
      // only interested in error of the evolved passive scalar profiles
      constexpr int scalar_norm = NSCALARS > 0 ? NSCALARS : 1.0;

      for (int n=0; n<NSCALARS; ++n) {
        for (int k=kl; k<=ku; k++) {
          for (int j=jl; j<=ju; j++) {
            pmb->pcoord->CellVolume(k, j, il, iu, vol);
            for (int i=il; i<=iu; i++) {
              Real cell_ave;
              if (iprob == 1) {
                if (use_gl_quadrature) {
                  Real xl, xu;
                  xl = pmb->pcoord->x1f(i);
                  xu = pmb->pcoord->x1f(i+1);

                  // GL implementation returns total integral, not ave. Divide by delta_r
                  Real cell_quad = GaussLegendre::integrate(N_gl, IntegrandFinal, xl, xu);
                  cell_ave = cell_quad/vol(i);
                  // assuming that the Gaussian profile is 1D in radial coordinate to pull
                  // out 2x integrals from the triple volume integral
                  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0 ||    // dy*dz
                      std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {  // dz*dphi
                    cell_ave *= pmb->pcoord->dx2f(j);
                    cell_ave *= pmb->pcoord->dx3f(k);
                  } else { // if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0)
                    // sin(theta)*dtheta*dphi
                    cell_ave *= std::cos(pmb->pcoord->x2f(j)) -
                                std::cos(pmb->pcoord->x2f(j+1));
                    cell_ave *= pmb->pcoord->dx3f(k);
                  }
                } else {
                  // Use standard midpoint approximation with cell-centered coords:
                  cell_ave = FinalGaussianProfile(pmb->pcoord->x1v(i));
                }
              } else if (iprob == 2) {
                if (use_gl_quadrature) {
                  Real xl, xu;
                  xl = pmb->pcoord->x2f(j);
                  xu = pmb->pcoord->x2f(j+1);
                  Real cell_quad = GaussLegendre::integrate(N_gl, IntegrandFinal, xl, xu);
                  cell_ave = cell_quad*pmb->pcoord->dx3f(k)/vol(i);
                  cell_ave *= ONE_3RD*(SQR(pmb->pcoord->x1f(i+1))*pmb->pcoord->x1f(i+1)
                                       - SQR(pmb->pcoord->x1f(i))*pmb->pcoord->x1f(i));
                } else {
                  cell_ave = FinalCosineProfile(pmb->pcoord->x2v(j));
                }
              } // end if iprob == 2
              Real sol = 1.0/scalar_norm*cell_ave;
              Real abs_diff = std::abs(sol - pmb->pscalars->s(n,k,j,i));
              l1_err[n] += abs_diff*vol(i);
              max_err[n] = std::max(abs_diff, max_err[n]);
              total_vol += vol(i);
            }
          }
        }
      }
    }
#ifdef MPI_PARALLEL
    if (Globals::my_rank == 0) {
      MPI_Reduce(MPI_IN_PLACE, &l1_err, NSCALARS, MPI_ATHENA_REAL, MPI_SUM, 0,
                 MPI_COMM_WORLD);
      MPI_Reduce(MPI_IN_PLACE, &total_vol, 1, MPI_ATHENA_REAL, MPI_SUM, 0,
                 MPI_COMM_WORLD);
      MPI_Reduce(MPI_IN_PLACE, &max_err, NSCALARS, MPI_ATHENA_REAL, MPI_MAX, 0,
                 MPI_COMM_WORLD);
    } else {
      MPI_Reduce(&l1_err, &l1_err, NSCALARS, MPI_ATHENA_REAL, MPI_SUM, 0,
                 MPI_COMM_WORLD);
      MPI_Reduce(&total_vol, &total_vol, 1, MPI_ATHENA_REAL, MPI_SUM, 0,
                 MPI_COMM_WORLD);
      MPI_Reduce(&max_err, &max_err, NSCALARS, MPI_ATHENA_REAL, MPI_MAX, 0,
                 MPI_COMM_WORLD);
    }
#endif

    // only the root process outputs the data
    if (Globals::my_rank == 0) {
      // normalize errors by total domain volume
      for (int i=0; i<NSCALARS; ++i) l1_err[i] = l1_err[i]/total_vol;
      // open output file and write out errors
      std::stringstream msg;
      std::string fname;
      if (iprob == 1) {
        fname.assign("mignone_radial-errors.dat");
      } else {
        fname.assign("mignone_meridional-errors.dat");
      }
      FILE *pfile;

      // The file exists -- reopen the file in append mode
      if ((pfile = std::fopen(fname.c_str(), "r")) != nullptr) {
        if ((pfile = std::freopen(fname.c_str(), "a", pfile)) == nullptr) {
          msg << "### FATAL ERROR in function Mesh::UserWorkAfterLoop"
              << std::endl << "Error output file could not be opened" <<std::endl;
          ATHENA_ERROR(msg);
        }

        // The file does not exist -- open the file in write mode and add headers
      } else {
        if ((pfile = std::fopen(fname.c_str(), "w")) == nullptr) {
          msg << "### FATAL ERROR in function Mesh::UserWorkAfterLoop"
              << std::endl << "Error output file could not be opened" <<std::endl;
          ATHENA_ERROR(msg);
        }
        std::fprintf(pfile, "# Nx1  Nx2  Nx3  Ncycle  ");
        for (int n=0; n<NSCALARS; ++n)
          std::fprintf(pfile, "s%d_L1  ", n);
        for (int n=0; n<NSCALARS; ++n)
          std::fprintf(pfile, "s%d_max  ", n);
        std::fprintf(pfile, "\n");
      }

      // write errors
      std::fprintf(pfile, "%d  %d", mesh_size.nx1, mesh_size.nx2);
      std::fprintf(pfile, "  %d  %d", mesh_size.nx3, ncycle);
      for (int n=0; n<NSCALARS; ++n)
        std::fprintf(pfile, "  %e", l1_err[n]);
      for (int n=0; n<NSCALARS; ++n)
        std::fprintf(pfile, "  %e", max_err[n]);
      std::fprintf(pfile, "\n");
      std::fclose(pfile);
    }
  } // if NSCALARS > 0
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)

//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  AthenaArray<Real> vol(ncells1);
  // initialize conserved variables
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      pcoord->CellVolume(k, j, is, ie, vol);
      for (int i=is; i<=ie; i++) {
        // background fluid:
        phydro->u(IDN,k,j,i) = d0;
        phydro->u(IM3,k,j,i) = 0.0;
        // assuming isothermal EOS:
        //  phydro->u(IEN,k,j,i) =
        Real cell_ave;

        //--- iprob=1
        if (iprob == 1) {
          // Note: GL quadrature is only ever used for the passive scalar profile (even
          // though it would be straightforward to compute the analytic integral for the
          // specific Gaussian profile) for generality with any initial condition
          Real xl, xu;
          xl = pcoord->x1f(i);
          xu = pcoord->x1f(i+1);
          phydro->u(IM2,k,j,i) = 0.0;
          // The cell-averaged linear velocity profile is always computed from the
          // analytic integral expression:
          if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0 ||    // dy*dz
              std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {  // dz*dphi
            phydro->u(IM1,k,j,i) =
                d0*alpha*ONE_3RD*(pcoord->x1f(i+1)*SQR(pcoord->x1f(i+1)) -
                                  pcoord->x1f(i)*SQR(pcoord->x1f(i)));
            phydro->u(IM1,k,j,i) *= pcoord->dx2f(j)*pcoord->dx3f(k)/vol(i);
          } else { // if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0)
            // sin(theta)*dtheta*dphi
            phydro->u(IM1,k,j,i) =
                d0*alpha*0.25*(SQR(SQR(pcoord->x1f(i+1))) - SQR(SQR(pcoord->x1f(i))));
            phydro->u(IM1,k,j,i) *= pcoord->dx3f(k)/vol(i);
            phydro->u(IM1,k,j,i) *= std::cos(pcoord->x2f(j)) - std::cos(pcoord->x2f(j+1));
          }
          // vs. midpoint approximation for intiialization of linear velocity profile:
          // phydro->u(IM1,k,j,i) = d0*alpha*pcoord->x1v(i);

          if (use_gl_quadrature) {
            // Use Gauss-Legendre quadrature rules to compute cell-averaged passive scalar
            // initial condition based on the pointwise analytic formula
            // GL implementation returns total integral, not average. Divide by delta_r
            Real cell_quad = GaussLegendre::integrate(N_gl, IntegrandInitial, xl, xu);
            cell_ave = cell_quad/vol(i);
            // assume that the Gaussian profile is 1D in radial coordinate to pull
            // out 2x integrals from the triple volume integral
            if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0 ||    // dy*dz
                std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {  // dz*dphi
              cell_ave *= pcoord->dx2f(j);
              cell_ave *= pcoord->dx3f(k);
            } else { // if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0)
              // sin(theta)*dtheta*dphi
              cell_ave *= std::cos(pcoord->x2f(j)) - std::cos(pcoord->x2f(j+1));
              cell_ave *= pcoord->dx3f(k);
            }
          } else {
            // Use standard midpoint approximation with cell centered coords:
            cell_ave = InitialGaussianProfile(pcoord->x1v(i));
          }
        } else if (iprob == 2) {
          Real xl, xu;
          xl = pcoord->x2f(j);
          xu = pcoord->x2f(j+1);
          phydro->u(IM1,k,j,i) = 0.0;
          // exact integral of cell-averaged linear velocity profile:
          phydro->u(IM2,k,j,i) =
              d0*alpha*(std::sin(xu) - xu*std::cos(xu)
                        - (std::sin(xl) - xl*std::cos(xl)));
          phydro->u(IM2,k,j,i) *= pcoord->dx3f(k)/vol(i);
          phydro->u(IM2,k,j,i) *= 0.25*(SQR(SQR(pcoord->x1f(i+1)))
                                        - SQR(SQR(pcoord->x1f(i))));
          // vs. midpoint approximation for intiialization of linear velocity profile:
          // (note, Mignone assumes r=1 for the 1D variant of this problem; need to scale
          // v_theta with radius in this 2D variant)
          // phydro->u(IM2,k,j,i) = d0*alpha*pcoord->x2v(j)*pcoord->x1f(i);
          if (use_gl_quadrature) {
            Real cell_quad = GaussLegendre::integrate(N_gl, IntegrandInitial, xl, xu);
            cell_ave = cell_quad*pcoord->dx3f(k)/vol(i);
            cell_ave *= ONE_3RD*(SQR(pcoord->x1f(i+1))*pcoord->x1f(i+1)
                                 - SQR(pcoord->x1f(i))*pcoord->x1f(i));
          } else {
            cell_ave = InitialCosineProfile(pcoord->x2v(j));
          }
        } // end if iprob == 2

        // uniformly fill all scalars to have equal concentration
        constexpr int scalar_norm = NSCALARS > 0 ? NSCALARS : 1.0;
        if (NSCALARS > 0) {
          for (int n=0; n<NSCALARS; ++n) {
            pscalars->s(n,k,j,i) = 1.0/scalar_norm*cell_ave*d0;
          }
        }
      }
    }
  }
  return;
}

namespace {
Real InitialGaussianProfile(Real x1) {
  return std::exp(-SQR(a_width)*SQR(x1 - b_center));  // Mignone eq 73
}


Real FinalGaussianProfile(Real x1) {
  // hardcoding t=t_final to maintain compatiblity with GL quadrature routines
  Real x_initial = x1*std::exp(-alpha*t_final);
  Real q_initial = InitialGaussianProfile(x_initial);
  Real amp = std::exp(-(m_coord + 1)*alpha*t_final);
  return amp*q_initial;  // Mignone eq 72
}


Real InitialCosineProfile(Real x2) {
  Real shift_x2 = x2 - b_center;
  if (std::abs(shift_x2) < PI/a_width) {
    return SQR(0.5 + 0.5*std::cos(a_width*shift_x2));  // Mignone eq 77
  } else {
    return 0.0;
  }
}


Real FinalCosineProfile(Real x2) {
  Real x2_initial = x2*std::exp(-alpha*t_final);
  Real q_initial = InitialCosineProfile(x2_initial);
  Real amp = std::exp(-alpha*t_final)*std::sin(x2_initial)/sin(x2);
  return amp*q_initial;  // Mignone eq 76
}


Real InitialGaussianCylindricalIntegrand(Real x1) {
  return x1*InitialGaussianProfile(x1);
}


Real FinalGaussianCylindricalIntegrand(Real x1) {
  return x1*FinalGaussianProfile(x1);
}


Real InitialGaussianSphericalIntegrand(Real x1) {
  return SQR(x1)*InitialGaussianProfile(x1);
}


Real FinalGaussianSphericalIntegrand(Real x1) {
  return SQR(x1)*FinalGaussianProfile(x1);
}


Real InitialCosineSphericalIntegrand(Real x2) {
  return std::sin(x2)*InitialCosineProfile(x2);
}


Real FinalCosineSphericalIntegrand(Real x2) {
  return std::sin(x2)*FinalCosineProfile(x2);
}
} // namespace
