//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file msa.cpp
//! \brief Modified swing amplification problem generator.
//! REFERENCE: Kim & Ostriker, ApJ, 2001
//======================================================================================

// C headers

// C++ headers
#include <cmath>      // sqrt()
#include <fstream>    // ofstream
#include <iomanip>    // setprecision
#include <iostream>   // cout, endl
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
#include "../orbital_advection/orbital_advection.hpp"
#include "../outputs/outputs.hpp"
#include "../parameter_input.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

namespace {
Real cs, gm1, d0, p0, gconst;
Real Q, nJ, scaleH, beta, amp;
int nwx, nwy; // wavenumbers
Real x1size,x2size,x3size;
Real qshear, Omega0; // shear parameters
bool strat;
} // namespace

//======================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Init the Mesh properties
//======================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (!shear_periodic) {
    std::stringstream msg;
    msg << "### FATAL ERROR in msa.cpp ProblemGenerator" << std::endl
        << "This problem generator requires shearing box." << std::endl;
    ATHENA_ERROR(msg);
  }

  if (MAGNETIC_FIELDS_ENABLED) {
    std::stringstream msg;
    msg << "### FATAL ERROR in msa.cpp ProblemGenerator" << std::endl
        << "Magnetic field is not yet included." << std::endl;
    ATHENA_ERROR(msg);
  }

  if (mesh_size.nx2 == 1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in msa.cpp ProblemGenerator" << std::endl
        << "This problem does NOT work on a 1D grid." << std::endl;
    ATHENA_ERROR(msg);
  }

  x1size = mesh_size.x1max - mesh_size.x1min;
  x2size = mesh_size.x2max - mesh_size.x2min;
  x3size = mesh_size.x3max - mesh_size.x3min;

  // shearing box parameters
  qshear = pin->GetReal("orbital_advection","qshear");
  Omega0 = pin->GetReal("orbital_advection","Omega0");

  // hydro parameters
  if (NON_BAROTROPIC_EOS) {
    gm1 = (pin->GetReal("hydro","gamma") - 1.0);
  }

  // MSA parameters
  Q = pin->GetReal("problem","Q");
  nJ = pin->GetReal("problem","nJ");
  beta = pin->GetReal("problem","beta");
  amp = pin->GetReal("problem","amp");
  nwx = pin->GetInteger("problem","nwx");
  nwy = pin->GetInteger("problem","nwy");
  strat = pin->GetBoolean("problem","strat");
  cs = std::sqrt(4.0-2.0*qshear)/PI/nJ/Q;
  scaleH = 1.0/std::sqrt(TWO_PI*nJ); // scale height for hyperbolic secant^2 profile
  d0 = 1.0; // midplane density
  if (NON_BAROTROPIC_EOS) {
    p0 = SQR(cs)*d0/(gm1+1.0); // midplane pressure
  }

  if (SELF_GRAVITY_ENABLED) {
    gconst = nJ*SQR(cs);
    SetGravitationalConstant(gconst);
    Real eps = pin->GetOrAddReal("self_gravity","grav_eps", 0.0);
    SetGravityThreshold(eps);
  }
  return;
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  if (gid == 0) {
    std::cout << "cs = " << cs << std::endl;
    std::cout << "G = " << gconst << std::endl;
    std::cout << "[msa.cpp]: [Lx,Ly,Lz] = [" <<x1size <<","<<x2size
              <<","<<x3size<<"]"<<std::endl;
  }

  // set wavenumbers
  Real kx = (TWO_PI/x1size)*(static_cast<Real>(nwx));
  Real ky = (TWO_PI/x2size)*(static_cast<Real>(nwy));

  Real x1, x2, x3, rd, rp, rvx, rvy;
  Real den, prs;
  // update the physical variables as initial conditions
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        x1 = pcoord->x1v(i);
        x2 = pcoord->x2v(j);
        x3 = pcoord->x3v(k);

        if (strat) {
          den = d0*SQR(1.0/std::cosh(x3/scaleH));
          rd = amp*SQR(1.0/std::cosh(x3/scaleH))*std::cos(kx*x1 + ky*x2);
        } else {
          den = d0;
          rd = amp*std::cos(kx*x1 + ky*x2);
        }
        rvx = amp*kx/ky*std::sin(kx*x1 + ky*x2);
        rvy = amp*std::sin(kx*x1 + ky*x2);
        if (NON_BAROTROPIC_EOS) {
          prs = p0*std::pow(den/d0, gm1+1.0);
          rp = SQR(cs)*rd;
        }

        phydro->u(IDN,k,j,i) = (den+rd);
        phydro->u(IM1,k,j,i) = (den+rd)*rvx;
        phydro->u(IM2,k,j,i) = (den+rd)*rvy;
        if (pmy_mesh->shear_periodic) {
          phydro->u(IM2,k,j,i) -= (den+rd)*(qshear*Omega0*x1);
        }
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = (prs+rp)/gm1 + 0.5*(SQR(phydro->u(IM1,k,j,i)) +
                                                SQR(phydro->u(IM2,k,j,i)) +
                                                SQR(phydro->u(IM3,k,j,i))
                                                ) / phydro->u(IDN,k,j,i);
        }
      }
    }
  }
  return;
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief Compute L1 error in shearing waves and output to file
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  if (!pin->GetOrAddBoolean("problem", "compute_error", false)) return;

  if (!SELF_GRAVITY_ENABLED) {
    std::stringstream msg;
    msg << "### FATAL ERROR in msa.cpp ProblemGenerator" << std::endl
        << "Errors are valid only when SELF_GRAVITY_ENABLED" << std::endl;
    ATHENA_ERROR(msg);
  }

  // linear perturbation amplitudes at t=2
  Real tf = pin->GetReal("time","tlim");
  if ((tf!=2.0)||(Q!=2.0)||(nJ!=2.5)||(nwx!=-3.0)||(nwy!=1)||(amp!=1e-6)) {
    std::stringstream msg;
    msg << "### FATAL ERROR in msa.cpp ProblemGenerator" << std::endl
        << "The errors can be computed only with the following input parameters:"
        << std::endl
        << "time/tlim = 2.0" << std::endl
        << "problem/Q = 2.0" << std::endl
        << "problem/nJ = 2.5" << std::endl
        << "problem/nwx = -3" << std::endl
        << "problem/nxy = 1" << std::endl
        << "problem/amp = 1e-6" << std::endl;
    ATHENA_ERROR(msg);
  }
  Real d1a = 1.56657750393522897e-06;
  Real v1a = 6.99809836090060095e-07;
  Real v2a = -8.69560918031967031e-07;

  // set wavenumbers
  Real kx = (TWO_PI/x1size)*(static_cast<Real>(nwx));
  Real ky = (TWO_PI/x2size)*(static_cast<Real>(nwy));
  Real kxt = kx + qshear*ky*tf;

  // Initialize errors to zero
  Real l1_err[NHYDRO+NFIELD]{}, max_err[NHYDRO+NFIELD]{};

  for (int b=0; b<nblocal; ++b) {
    MeshBlock *pmb = my_blocks(b);
    int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je,
        ks = pmb->ks, ke = pmb->ke;

    //  Compute errors at cell centers
    Real d1, m1, m2, m3, e0, b1, b2, b3;
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          Real x1 = pmb->pcoord->x1v(i);
          Real x2 = pmb->pcoord->x2v(j);
          d1 = 1.0 + 2.0*d1a*std::cos(kxt*x1 + ky*x2);
          m1 = d1*2.0*v1a*std::sin(kxt*x1 + ky*x2);
          m2 = d1*(-qshear*Omega0*x1 + 2.0*v2a*std::sin(kxt*x1 + ky*x2));
          m3 = 0.0;
          if (NON_BAROTROPIC_EOS) {
            // (background pressure) = (sound speed)^2*(background density)/(Gamma)
            // (perturbed pressure) = (sound speed)^2*(perturbed density)
            Real p0 = SQR(cs)/(gm1+1.0) + SQR(cs)*(d1-1.0);
            e0 = p0/gm1 + 0.5*(SQR(m1)+SQR(m2)+SQR(m3))/d1;
            if (MAGNETIC_FIELDS_ENABLED) {
              // TODO(SMOON) set magnetic field
              b1 = 0.0;
              b2 = 0.0;
              b3 = 0.0;
              e0 += 0.5*(b1*b1 + b2*b2 + b3*b3);
            }
          }
          // Weight l1 error by cell volume
          Real vol = pmb->pcoord->GetCellVolume(k, j, i);

          l1_err[IDN] += std::abs(d1 - pmb->phydro->u(IDN,k,j,i))*vol;
          max_err[IDN] = std::max(
              static_cast<Real>(std::abs(d1 - pmb->phydro->u(IDN,k,j,i))),
              max_err[IDN]);
          l1_err[IM1] += std::abs(m1 - pmb->phydro->u(IM1,k,j,i))*vol;
          l1_err[IM2] += std::abs(m2 - pmb->phydro->u(IM2,k,j,i))*vol;
          l1_err[IM3] += std::abs(m3 - pmb->phydro->u(IM3,k,j,i))*vol;
          max_err[IM1] = std::max(
              static_cast<Real>(std::abs(m1 - pmb->phydro->u(IM1,k,j,i))),
              max_err[IM1]);
          max_err[IM2] = std::max(
              static_cast<Real>(std::abs(m2 - pmb->phydro->u(IM2,k,j,i))),
              max_err[IM2]);
          max_err[IM3] = std::max(
              static_cast<Real>(std::abs(m3 - pmb->phydro->u(IM3,k,j,i))),
              max_err[IM3]);

          if (NON_BAROTROPIC_EOS) {
            l1_err[IEN] += std::abs(e0 - pmb->phydro->u(IEN,k,j,i))*vol;
            max_err[IEN] = std::max(
                static_cast<Real>(std::abs(e0-pmb->phydro->u(IEN,k,j,i))),
                max_err[IEN]);
          }

          if (MAGNETIC_FIELDS_ENABLED) {
            Real db1 = std::abs(b1 - pmb->pfield->bcc(IB1,k,j,i));
            Real db2 = std::abs(b2 - pmb->pfield->bcc(IB2,k,j,i));
            Real db3 = std::abs(b3 - pmb->pfield->bcc(IB3,k,j,i));

            l1_err[NHYDRO + IB1] += db1*vol;
            l1_err[NHYDRO + IB2] += db2*vol;
            l1_err[NHYDRO + IB3] += db3*vol;
            max_err[NHYDRO + IB1] = std::max(db1, max_err[NHYDRO+IB1]);
            max_err[NHYDRO + IB2] = std::max(db2, max_err[NHYDRO+IB2]);
            max_err[NHYDRO + IB3] = std::max(db3, max_err[NHYDRO+IB3]);
          }
        }
      }
    }
  }
  Real rms_err = 0.0, max_max_over_l1 = 0.0;

#ifdef MPI_PARALLEL
  if (Globals::my_rank == 0) {
    MPI_Reduce(MPI_IN_PLACE, &l1_err, (NHYDRO+NFIELD), MPI_ATHENA_REAL, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &max_err, (NHYDRO+NFIELD), MPI_ATHENA_REAL, MPI_MAX, 0,
               MPI_COMM_WORLD);
  } else {
    MPI_Reduce(&l1_err, &l1_err, (NHYDRO+NFIELD), MPI_ATHENA_REAL, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&max_err, &max_err, (NHYDRO+NFIELD), MPI_ATHENA_REAL, MPI_MAX, 0,
               MPI_COMM_WORLD);
  }
#endif

  // only the root process outputs the data
  if (Globals::my_rank == 0) {
    // normalize errors by number of cells
    Real vol= (mesh_size.x1max - mesh_size.x1min)*(mesh_size.x2max - mesh_size.x2min)
              *(mesh_size.x3max - mesh_size.x3min);
    for (int i=0; i<(NHYDRO+NFIELD); ++i) l1_err[i] = l1_err[i]/vol;
    // compute rms error
    for (int i=0; i<(NHYDRO+NFIELD); ++i) {
      rms_err += SQR(l1_err[i]);
      max_max_over_l1 = std::max(max_max_over_l1, (max_err[i]/l1_err[i]));
    }
    rms_err = std::sqrt(rms_err);

    // open output file and write out errors
    std::string fname;
    fname.assign("msa-errors.dat");
    std::stringstream msg;
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
      std::fprintf(pfile, "RMS-L1-Error  d_L1  M1_L1  M2_L1  M3_L1  E_L1 ");
      if (MAGNETIC_FIELDS_ENABLED) std::fprintf(pfile, "  B1c_L1  B2c_L1  B3c_L1");
      std::fprintf(pfile, "  Largest-Max/L1  d_max  M1_max  M2_max  M3_max  E_max ");
      if (MAGNETIC_FIELDS_ENABLED) std::fprintf(pfile, "  B1c_max  B2c_max  B3c_max");
      std::fprintf(pfile, "\n");
    }

    // write errors
    std::fprintf(pfile, "%d  %d", mesh_size.nx1, mesh_size.nx2);
    std::fprintf(pfile, "  %d  %d", mesh_size.nx3, ncycle);
    std::fprintf(pfile, "  %e  %e", rms_err, l1_err[IDN]);
    std::fprintf(pfile, "  %e  %e  %e", l1_err[IM1], l1_err[IM2], l1_err[IM3]);
    if (NON_BAROTROPIC_EOS)
      std::fprintf(pfile, "  %e", l1_err[IEN]);
    if (MAGNETIC_FIELDS_ENABLED) {
      std::fprintf(pfile, "  %e", l1_err[NHYDRO+IB1]);
      std::fprintf(pfile, "  %e", l1_err[NHYDRO+IB2]);
      std::fprintf(pfile, "  %e", l1_err[NHYDRO+IB3]);
    }
    std::fprintf(pfile, "  %e  %e  ", max_max_over_l1, max_err[IDN]);
    std::fprintf(pfile, "%e  %e  %e", max_err[IM1], max_err[IM2], max_err[IM3]);
    if (NON_BAROTROPIC_EOS)
      std::fprintf(pfile, "  %e", max_err[IEN]);
    if (MAGNETIC_FIELDS_ENABLED) {
      std::fprintf(pfile, "  %e", max_err[NHYDRO+IB1]);
      std::fprintf(pfile, "  %e", max_err[NHYDRO+IB2]);
      std::fprintf(pfile, "  %e", max_err[NHYDRO+IB3]);
    }
    std::fprintf(pfile, "\n");
    std::fclose(pfile);
  }

  return;
}
