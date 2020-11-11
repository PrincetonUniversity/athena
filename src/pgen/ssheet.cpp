//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file ssheet.cpp
//  \brief Shearing wave problem generator.
//  Several different initial conditions:
//  - ipert = 1  pure shearing background flow
//  - ipert = 2  epicycle motion (0.1 c_s initial kick in radial)
//  - ipert = 3  shwave test for compressible flow
//               JG: Johnson & Gammie 2005, ApJ, 626, 978
//
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
#include "../bvals/bvals.hpp"
#include "../bvals/bvals_interfaces.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../orbital_advection/orbital_advection.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"

#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator requires NOT MHD"
#endif

namespace {
Real amp, nwx, nwy; // amplitude, Wavenumbers
int ipert; // initial pattern
Real gm1,iso_cs;
Real x1size,x2size,x3size;
Real Omega_0,qshear;
int shboxcoord;

Real Historydvyc(MeshBlock *pmb, int iout);
} // namespace

//======================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Init the Mesh properties
//======================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // read ipert parameter
  ipert = pin->GetInteger("problem","ipert");

  if (ipert == 3) {
    amp = pin->GetReal("problem","amp");
    nwx = pin->GetInteger("problem","nwx");
    nwy = pin->GetInteger("problem","nwy");
    if (nwx == 0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in ssheet.cpp ProblemGenerator"   << std::endl
          << "Parameter nwx must be non-zero." << std::endl;
      ATHENA_ERROR(msg);
    }
    AllocateUserHistoryOutput(1);
    EnrollUserHistoryOutput(0, Historydvyc, "dvyc", UserHistoryOperation::sum);
  } else if (ipert == 1 || ipert == 2) {
    amp = 0.0;
    nwx = 0;
    nwy = 0;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in ssheet.cpp ProblemGenerator"   << std::endl
        << "This problem requires that ipert is from 1 to 3." << std::endl;
    ATHENA_ERROR(msg);
  }

  if (!shear_periodic) {
    std::stringstream msg;
    msg << "### FATAL ERROR in ssheet.cpp ProblemGenerator" << std::endl
        << "This problem generator requires shearing box"   << std::endl;
    ATHENA_ERROR(msg);
  }

  return;
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // shearing sheet parameter
  Omega_0 = porb->Omega0;
  qshear  = porb->qshear;
  shboxcoord = pin->GetOrAddInteger("problem","shboxcoord",1);

  // conditional error branch
  if (pmy_mesh->mesh_size.nx2 == 1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in ssheet.cpp ProblemGenerator" << std::endl
        << "Shearing wave sheet does NOT work on a 1D grid" << std::endl;
    ATHENA_ERROR(msg);
  }

  int il = is - NGHOST; int iu = ie + NGHOST;
  int jl = js - NGHOST; int ju = je + NGHOST;
  int kl = ks;          int ku = ke;
  if (block_size.nx3 > 1) {
    kl = ks - NGHOST;
    ku = ke + NGHOST;
  }

  Real d0 = 1.0;
  Real p0 = 1e-6;

  if (NON_BAROTROPIC_EOS) {
    gm1 = (peos->GetGamma() - 1.0);
    iso_cs = std::sqrt((gm1+1.0)*p0/d0);
  } else {
    iso_cs = peos->GetIsoSoundSpeed();
    p0 = d0*SQR(iso_cs);
  }
  x1size = pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min;
  x2size = pmy_mesh->mesh_size.x2max - pmy_mesh->mesh_size.x2min;
  x3size = pmy_mesh->mesh_size.x3max - pmy_mesh->mesh_size.x3min;

  if (gid == 0) {
    std::cout << "iso_cs = " << iso_cs << std::endl;
    std::cout << "d0 = " << d0 << std::endl;
    std::cout << "p0 = " << p0 << std::endl;
    std::cout << "ipert  = " << ipert  << std::endl;
    std::cout << "[ssheet.cpp]: [Lx,Ly,Lz] = [" <<x1size <<","<<x2size
              <<","<<x3size<<"]"<<std::endl;
  }

  // calculate wave number just for ipert = 3
  Real kx = (TWO_PI/x1size)*(static_cast<Real>(nwx));
  Real ky = (TWO_PI/x2size)*(static_cast<Real>(nwy));

  Real x1, x2, rd, rp, rvx, rvy;
  // update the physical variables as initial conditions
  for (int k=kl; k<=ku; k++) {
    for (int j=jl; j<=ju; j++) {
      for (int i=il; i<=iu; i++) {
        x1 = pcoord->x1v(i);
        x2 = pcoord->x2v(j);
        rd = d0;
        rp = p0;
        if (ipert == 1) {
          // 1) pure shear bg flow:
          phydro->u(IDN,k,j,i) = rd;
          phydro->u(IM1,k,j,i) = 0.0;
          phydro->u(IM2,k,j,i) = 0.0;
          if(!porb->orbital_advection_defined)
            phydro->u(IM2,k,j,i) -= rd*qshear*Omega_0*x1;
          phydro->u(IM3,k,j,i) = 0.0;
        } else if (ipert == 2) {
          // 2) epicyclic oscillation
          if (shboxcoord == 1) { // x-y shear
            rvx = 0.1*iso_cs;
            rvy = 0.0;
            phydro->u(IDN,k,j,i) = rd;
            phydro->u(IM1,k,j,i) = rd*rvx;
            phydro->u(IM2,k,j,i) = -rd*rvy;
            if(!porb->orbital_advection_defined)
            phydro->u(IM2,k,j,i) -= rd*qshear*Omega_0*x1;
            phydro->u(IM3,k,j,i) = 0.0;
          } else { // x-z plane
            rvx = 0.1*iso_cs;
            rvy = 0.0;
            phydro->u(IDN,k,j,i) = rd;
            phydro->u(IM1,k,j,i) = rd*rvx;
            phydro->u(IM2,k,j,i) = 0.0;
            phydro->u(IM3,k,j,i) = -rd*(rvy+qshear*Omega_0*x1);
          }
        } else if (ipert == 3) {
          // 3) JG HD shwave test
          rvx = amp*iso_cs*std::cos(kx*x1 + ky*x2);
          rvy = amp*iso_cs*(ky/kx)*std::cos(kx*x1 + ky*x2);
          phydro->u(IDN,k,j,i) = rd;
          phydro->u(IM1,k,j,i) = -rd*rvx;
          phydro->u(IM2,k,j,i) = -rd*rvy;
          if(!porb->orbital_advection_defined)
            phydro->u(IM2,k,j,i) -= rd*qshear*Omega_0*x1;
          phydro->u(IM3,k,j,i) = 0.0;
        } else {
          std::stringstream msg;
          msg << "### FATAL ERROR in ssheet.cpp ProblemGenerator" << std::endl
              << "Shearing wave sheet ipert=" << ipert << " is unrecognized" << std::endl;
          ATHENA_ERROR(msg);
        }
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = rp/gm1 + 0.5*(SQR(phydro->u(IM1,k,j,i)) +
                                               SQR(phydro->u(IM2,k,j,i)) +
                                               SQR(phydro->u(IM3,k,j,i))
                                              ) / phydro->u(IDN,k,j,i);
        }
      }
    }
  }

  return;
}

namespace {

Real Historydvyc(MeshBlock *pmb, int iout) {
  Real kx = (TWO_PI/x1size)*(static_cast<Real>(nwx));
  Real ky = (TWO_PI/x2size)*(static_cast<Real>(nwy));
  kx += qshear*Omega_0*pmb->pmy_mesh->time*ky;
  Real dvyc = 0.0;
  AthenaArray<Real> volume; // 1D array of volumes
  volume.NewAthenaArray(pmb->ncells1);
  Real tvol = x1size*x2size*x3size;
  for (int k=pmb->ks; k<=pmb->ke  ; k++) {
    for (int j=pmb->js; j<=pmb->je  ; j++) {
      pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, volume);
      for (int i=pmb->is; i<=pmb->ie  ; i++) {
        Real x1 = pmb->pcoord->x1v(i);
        Real x2 = pmb->pcoord->x2v(j);
        Real CS = std::cos(kx*x1+ky*x2);
        Real dvy = pmb->phydro->w(IVY,k,j,i);
        if(!pmb->porb->orbital_advection_defined)
          dvy += qshear*Omega_0*x1;
        dvyc += volume(i)*2.0*dvy*CS/tvol;
      }
    }
  }
  return dvyc;
}
} // namespace
