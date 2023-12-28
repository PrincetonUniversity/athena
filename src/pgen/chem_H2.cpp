//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file chem_H2.cpp
//! \brief problem generator, H2 chemistry
//======================================================================================

// C headers

// C++ headers
#include <algorithm>  // std::find()
#include <cstdio>
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // std::runtime_error()
#include <string>     // c_str()
#include <vector>     // vector container

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../chem_rad/chem_rad.hpp"
#include "../chem_rad/integrators/rad_integrators.hpp"
#include "../chemistry/utils/thermo.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"
#include "../units/units.hpp"

Real threshold;
int RefinementCondition(MeshBlock *pmb);
namespace {
Real HistoryT(MeshBlock *pmb, int iout); // average temperature output
} // namespace

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//========================================================================================
void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (adaptive) {
    EnrollUserRefinementCondition(RefinementCondition);
    threshold = pin->GetReal("problem", "thr");
  }
  AllocateUserHistoryOutput(1);
  EnrollUserHistoryOutput(0, HistoryT, "T");
  return;
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief initialize problem by reading in vtk file.
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // dimensions of meshblock
  //const int Nx = ie - is + 1;
  //const int Ny = je - js + 1;
  //const int Nz = ke - ks + 1;

  // read input parameters
  const Real nH = pin->GetReal("problem", "nH"); // density
  const Real vx = pin->GetOrAddReal("problem", "vx_kms", 0); // velocity x
  const Real fH1 = pin->GetOrAddReal("problem", "fH1", 0.); // H abundance at x>1
  // mean and std of the initial gaussian profile
  const Real gaussian_mean = pin->GetOrAddReal("problem", "gaussian_mean", 0.5);
  const Real gaussian_std = pin->GetOrAddReal("problem", "gaussian_std", 0.1);
  const Real iso_cs = pin->GetReal("hydro", "iso_sound_speed");
  const Real pres = nH*SQR(iso_cs);
  const Real gm1  = peos->GetGamma() - 1.0;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        // density
        phydro->u(IDN, k, j, i) = nH;
        // velocity, x direction
        phydro->u(IM1, k, j, i) = nH*vx;
        // energy
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN, k, j, i) = pres/gm1 + 0.5*nH*SQR(vx);
        }
      }
    }
  }

  // intialize chemical species
  if (NSPECIES > 0) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          for (int ispec=0; ispec < NSPECIES; ++ispec) {
            if (CHEMISTRY_ENABLED) {
              Real x1 = pcoord->x1v(i);
              // gaussian initial H abundance in [0, 1), and no H in [1, 2]
              if (x1 <= 1) {
                pscalars->s(0, k, j, i) = std::exp( -SQR(x1-gaussian_mean)
                                                 /(2.*SQR(gaussian_std)) )*nH; // H
                pscalars->s(1, k, j, i) = 0.5*(nH - pscalars->s(0, k, j, i)); // H2
              } else {
                pscalars->s(0, k, j, i) = fH1*nH; // H
                pscalars->s(1, k, j, i) = (1.-fH1)*0.5*nH; // H2
              }
            }
          }
        }
      }
    }
  }

  return;
}

//======================================================================================
//! \fn int RefinementCondition(MeshBlock *pmb)
//! \brief refinement condition: maximum gradient of each passive scalar profile
//======================================================================================

int RefinementCondition(MeshBlock *pmb) {
  int f2 = pmb->pmy_mesh->f2, f3 = pmb->pmy_mesh->f3;
  AthenaArray<Real> &r = pmb->pscalars->r;
  Real maxeps = 0.0;
  if (f3) {
    for (int n=0; n<NSPECIES; ++n) {
      for (int k=pmb->ks-1; k<=pmb->ke+1; k++) {
        for (int j=pmb->js-1; j<=pmb->je+1; j++) {
          for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
            Real eps = std::sqrt(SQR(0.5*(r(n,k,j,i+1) - r(n,k,j,i-1)))
                                 + SQR(0.5*(r(n,k,j+1,i) - r(n,k,j-1,i)))
                                 + SQR(0.5*(r(n,k+1,j,i) - r(n,k-1,j,i))));
            // /r(n,k,j,i); Do not normalize by scalar, since (unlike IDN and IPR) there
            // are are no physical floors / r=0 might be allowed. Compare w/ blast.cpp.
            maxeps = std::max(maxeps, eps);
          }
        }
      }
    }
  } else if (f2) {
    int k = pmb->ks;
    for (int n=0; n<NSPECIES; ++n) {
      for (int j=pmb->js-1; j<=pmb->je+1; j++) {
        for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
          Real eps = std::sqrt(SQR(0.5*(r(n,k,j,i+1) - r(n,k,j,i-1)))
                               + SQR(0.5*(r(n,k,j+1,i) - r(n,k,j-1,i)))); // /r(n,k,j,i);
          maxeps = std::max(maxeps, eps);
        }
      }
    }
  } else {
    return 0;
  }

  if (maxeps > threshold) return 1;
  if (maxeps < 0.25*threshold) return -1;
  return 0;
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  if (!pin->GetOrAddBoolean("problem", "compute_error", false)) return;

  // read input parameters
  const Real nH = pin->GetReal("problem", "nH"); // density
  const Real vx = pin->GetOrAddReal("problem", "vx_kms", 0); // velocity x
  const Real fH1 = pin->GetOrAddReal("problem", "fH1", 0.); // H abundance at x>1
  const Real gaussian_mean = pin->GetOrAddReal("problem", "gaussian_mean", 0.5);
  const Real gaussian_std = pin->GetOrAddReal("problem", "gaussian_std", 0.1);
  // chemistry parameters
  Units *punit = my_blocks(0)->pmy_mesh->punit;
  const Real xi_cr = pin->GetOrAddReal("chemistry", "xi_cr", 2e-16);
  const Real kcr = xi_cr * 3.;
  const Real kgr = 3e-17;
  const Real a1 = kcr + 2.*nH*kgr;
  const Real a2 = kcr;
  // cooling parameters
  const Real iso_cs = pin->GetReal("hydro", "iso_sound_speed");
  const Real gm = pin->GetReal("hydro", "gamma");
  const Real ED0  = SQR(iso_cs) / (gm - 1.0);
  const Real CvHI = Thermo::CvCold(0., 0.1, 0.);
  const Real T0 =  (ED0 * punit->code_energydensity_cgs)  / CvHI;
  const Real tdust = 2 * CvHI / (3.2e-34 * nH);

  // end of the simulation time
  const Real tchem = time*punit->code_time_cgs;
  const Real mu = gaussian_mean + vx*time;
  const Real xg_min = vx*time;
  const Real xg_max = xg_min + 1.;
  // only compute error if the Gaussian profile did not travel outside of the
  // simulation domain at the end of the simulation
  if (xg_max > mesh_size.x1max) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function Mesh::UserWorkAfterLoop"
      << std::endl << "Gaussian profile outside of the simulation domain" <<std::endl;
    ATHENA_ERROR(msg);
  }

  // Initialize errors to zero
  Real l1_err[NSPECIES]{}, max_err[NSPECIES]{}, cons_err[1]{},
       l1_err_T[1]{}, max_err_T[1]{};
  Real T1_a, T1_s;

  for (int b=0; b<nblocal; ++b) {
    MeshBlock *pmb = my_blocks(b);
    int il = pmb->is, iu = pmb->ie, jl = pmb->js, ju = pmb->je,
        kl = pmb->ks, ku = pmb->ke;
    //  Compute errors at cell centers
    for (int k=kl; k<=ku; k++) {
      for (int j=jl; j<=ju; j++) {
        for (int i=il; i<=iu; i++) {
          Real x = pmb->pcoord->x1v(i);
          Real fH0 = 0;
          if ( (x < xg_min) || (x > xg_max) ) {
            fH0 = fH1;
          } else {
            fH0 = std::exp( -SQR(x-mu)/(2.*SQR(gaussian_std)) );
          }
          Real fH = (fH0 - a2/a1)*std::exp(-a1*tchem) + a2/a1;
          Real fH2 = 0.5*(1. - fH);
          T1_a = 1. / SQR( tchem/tdust + 1./std::sqrt(T0) ); // analytic T
          T1_s = pmb->phydro->w(IPR,k,j,i)/pmb->phydro->w(IDN,k,j,i)/(gm-1)
                        * punit->code_energydensity_cgs / CvHI; // simulation T
          // Weight l1 error by cell volume
          Real vol = pmb->pcoord->GetCellVolume(k, j, i);
          l1_err[0] += std::abs(fH - pmb->pscalars->r(0,k,j,i))*vol;
          max_err[0] = std::max(
              static_cast<Real>(std::abs(fH - pmb->pscalars->r(0,k,j,i))),
              max_err[0]);
          l1_err[1] += std::abs(fH2 - pmb->pscalars->r(1,k,j,i))*vol;
          max_err[1] = std::max(
              static_cast<Real>(std::abs(fH2 - pmb->pscalars->r(1,k,j,i))),
              max_err[1]);
          if (NON_BAROTROPIC_EOS) {
            l1_err_T[0] += std::abs(T1_a - T1_s)*vol;
            max_err_T[0] = std::max(
                static_cast<Real>(std::abs(T1_a - T1_s)), max_err_T[0]);
          }
          cons_err[0] += std::abs(pmb->pscalars->r(0,k,j,i) +
                                  2*pmb->pscalars->r(1,k,j,i) - 1.)*vol;
        }
      }
    }
  }

#ifdef MPI_PARALLEL
  if (Globals::my_rank == 0) {
    MPI_Reduce(MPI_IN_PLACE, &l1_err, NSPECIES, MPI_ATHENA_REAL, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &max_err, NSPECIES, MPI_ATHENA_REAL, MPI_MAX, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &cons_err, 1, MPI_ATHENA_REAL, MPI_SUM, 0,
               MPI_COMM_WORLD);
    if (NON_BAROTROPIC_EOS) {
      MPI_Reduce(MPI_IN_PLACE, &l1_err_T, 1, MPI_ATHENA_REAL, MPI_SUM, 0,
                 MPI_COMM_WORLD);
      MPI_Reduce(MPI_IN_PLACE, &max_err_T, 1, MPI_ATHENA_REAL, MPI_MAX, 0,
                 MPI_COMM_WORLD);
    }
  } else {
    MPI_Reduce(&l1_err, &l1_err, NSPECIES, MPI_ATHENA_REAL, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&max_err, &max_err, NSPECIES, MPI_ATHENA_REAL, MPI_MAX, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&cons_err, &cons_err, 1, MPI_ATHENA_REAL, MPI_SUM, 0,
               MPI_COMM_WORLD);
    if (NON_BAROTROPIC_EOS) {
      MPI_Reduce(&l1_err, &l1_err_T, 1, MPI_ATHENA_REAL, MPI_SUM, 0,
                 MPI_COMM_WORLD);
      MPI_Reduce(&max_err, &max_err_T, 1, MPI_ATHENA_REAL, MPI_MAX, 0,
                 MPI_COMM_WORLD);
    }
  }
#endif

  // only the root process outputs the data
  if (Globals::my_rank == 0) {
    // normalize errors by number of cells
    Real vol= (mesh_size.x1max - mesh_size.x1min)*(mesh_size.x2max - mesh_size.x2min)
              *(mesh_size.x3max - mesh_size.x3min);
    for (int i=0; i<NSPECIES; ++i) {
      l1_err[i] = l1_err[i]/vol;
    }
    cons_err[0] = cons_err[0]/vol;
    if (NON_BAROTROPIC_EOS) {
      l1_err_T[0] = l1_err_T[0]/vol;
    }

    // open output file and write out errors
    std::stringstream msg;
    std::string fname;
    fname.assign("chem_H2-errors.dat");
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
      for (int n=0; n<NSPECIES; ++n) {
        std::fprintf(pfile, "r%d_L1  ", n);
        std::fprintf(pfile, "r%d_max  ", n);
      }
      if (NON_BAROTROPIC_EOS) {
        std::fprintf(pfile, "T_L1  T_max  ");
      }
      std::fprintf(pfile, "cons_L1  \n");
    }

    // write errors
    std::fprintf(pfile, "%d  %d", mesh_size.nx1, mesh_size.nx2);
    std::fprintf(pfile, "  %d  %d", mesh_size.nx3, ncycle);
    for (int n=0; n<NSPECIES; ++n) {
      std::fprintf(pfile, "  %e", l1_err[n]);
      std::fprintf(pfile, "  %e", max_err[n]);
    }
    if (NON_BAROTROPIC_EOS) {
      std::fprintf(pfile, "  %e  %e", l1_err_T[0], max_err_T[0]);
    }
    std::fprintf(pfile, "  %e", cons_err[0]);
    std::fprintf(pfile, "\n");
    std::fclose(pfile);
  }
  return;
}

namespace {

Real HistoryT(MeshBlock *pmb, int iout) {
  const Real gm1  = pmb->peos->GetGamma() - 1.0;
  const Real CvHI = Thermo::CvCold(0., 0.1, 0.);
  int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je, ks = pmb->ks, ke = pmb->ke;
  Real T = 0;
  AthenaArray<Real> volume; // 1D array of volumes
  // allocate 1D array for cell volume used in usr def history
  volume.NewAthenaArray(pmb->ncells1);
  // total volume
  Real vol_tot= (pmb->pmy_mesh->mesh_size.x1max - pmb->pmy_mesh->mesh_size.x1min)
    *(pmb->pmy_mesh->mesh_size.x2max - pmb->pmy_mesh->mesh_size.x2min)
    *(pmb->pmy_mesh->mesh_size.x3max - pmb->pmy_mesh->mesh_size.x3min);

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      pmb->pcoord->CellVolume(k,j,pmb->is,pmb->ie,volume);
      for (int i=is; i<=ie; i++) {
        T += volume(i) * pmb->phydro->w(IPR,k,j,i)/pmb->phydro->w(IDN,k,j,i)/gm1
                        * pmb->pmy_mesh->punit->code_energydensity_cgs / CvHI;
      }
    }
  }
  T /= vol_tot;
  return T;
}

} // namespace
