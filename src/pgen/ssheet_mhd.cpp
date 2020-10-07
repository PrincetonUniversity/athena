//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file ssheet.cpp
//  \brief MHD shwave problem generator (Johnson 2007, ApJ, 660, 1375).
//
// Code must be configured using -b
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
#include "../parameter_input.hpp"

#if !MAGNETIC_FIELDS_ENABLED
#error "This problem generator requires magnetic fields"
#endif

namespace {
Real iso_cs, gm1, p0;
Real Omega_0, qshear; // shear
Real nwx, nwy, nwz; // Wavenumbers
Real x1size,x2size,x3size;
Real d0, epsilon, beta; // parameters
int shboxcoord;

Real HistoryDby(MeshBlock *pmb, int iout);
} // namespace

//======================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Init the Mesh properties
//======================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // initialize global variables
  // wave number
  nwx     = pin->GetInteger("problem","nwx");
  nwy     = pin->GetInteger("problem","nwy");
  nwz     = pin->GetInteger("problem","nwz");
  // parameters
  d0    = pin->GetReal("problem","d0");
  if (NON_BAROTROPIC_EOS)
    p0    = pin->GetReal("problem","p0");
  epsilon = pin->GetReal("problem","epsilon");
  beta    = pin->GetReal("problem","beta");
  shboxcoord = pin->GetOrAddInteger("problem","shboxcoord",1);

  AllocateUserHistoryOutput(1);
  EnrollUserHistoryOutput(0,HistoryDby, "dby", UserHistoryOperation::sum);

  if (!shear_periodic) {
    std::stringstream msg;
    msg << "### FATAL ERROR in hb3.cpp ProblemGenerator" << std::endl
        << "This problem generator requires shearing box" << std::endl;
    ATHENA_ERROR(msg);
  }

  return;
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  if (pmy_mesh->mesh_size.nx2 == 1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in ssheet_mhd.cpp ProblemGenerator" << std::endl
        << "Shwave for MHD works on 2D or 3D grid." << std::endl;
    ATHENA_ERROR(msg);
  }

  if (shboxcoord != 1 && pmy_mesh->mesh_size.nx3==1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in ssheet_mhd.cpp ProblemGenerator" << std::endl
        << "Shwave for MHD is incompatible with shearing box in x-z plane." << std::endl;
    ATHENA_ERROR(msg);
  }

  if (NON_BAROTROPIC_EOS) {
    gm1 = (peos->GetGamma() - 1.0);
    iso_cs = std::sqrt((gm1+1.0)*p0/d0);
  } else {
    iso_cs = peos->GetIsoSoundSpeed();
    p0 = d0*SQR(iso_cs);
  }

  // shearing sheet parameter
  Omega_0 = porb->Omega0;
  qshear  = porb->qshear;

  x1size = pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min;
  x2size = pmy_mesh->mesh_size.x2max - pmy_mesh->mesh_size.x2min;
  x3size = pmy_mesh->mesh_size.x3max - pmy_mesh->mesh_size.x3min;

  // calculate k = [kx, ky, kz]
  Real kx = (TWO_PI/x1size)*(static_cast<Real>(nwx));
  Real ky = (TWO_PI/x2size)*(static_cast<Real>(nwy));
  Real kz = (TWO_PI/x3size)*(static_cast<Real>(nwz));
  if (pmy_mesh->mesh_size.nx3 == 1) {
    kz = 0.0;
  }
  Real k2  =  SQR(kx)+SQR(ky)+SQR(kz);

  Real B2 =  d0*SQR(iso_cs)/beta;
  Real alpha = std::sqrt(B2/(SQR(kx)+SQR(ky)));
  Real Bx =  alpha*ky;
  Real By = -alpha*kx;
  Real Bz =  0.0;

  Real sch = iso_cs/Omega_0;
  Real cf1 = std::sqrt(B2*(1.0+beta));
  Real cf2 = epsilon*std::sqrt(sch*std::sqrt(k2*beta/(1+beta)));
  Real vd = cf1/std::sqrt(k2)*cf2;

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real x3 = pcoord->x3v(k);
        Real CS = std::cos(kx*x1+ky*x2+kz*x3);
        Real rd = d0*(1.0+cf2*CS);
        phydro->u(IDN,k,j,i) = rd;
        phydro->u(IM1,k,j,i) = rd*vd*kx*CS;
        phydro->u(IM2,k,j,i) = rd*vd*ky*CS;
        if(!porb->orbital_advection_defined)
          phydro->u(IM2,k,j,i) -= rd*qshear*Omega_0*x1;
        phydro->u(IM3,k,j,i) = rd*vd*kz*CS;
      }
    }
  }

  // nxN != ncellsN, in general. Allocate to extend through 2*ghost, regardless # dim
  int nx1 = block_size.nx1 + 2*NGHOST;
  int nx2 = block_size.nx2 + 2*NGHOST;
  int nx3 = block_size.nx3 + 2*NGHOST;
  AthenaArray<Real> dax, day, daz;
  dax.NewAthenaArray(nx3, nx2, nx1);
  day.NewAthenaArray(nx3, nx2, nx1);
  daz.NewAthenaArray(nx3, nx2, nx1);

  for (int k=ks; k<=ke+1; k++) {
    for (int j=js; j<=je+1; j++) {
      for (int i=is; i<=ie+1; i++) {
        Real x1 = pcoord->x1f(i);
        Real x2 = pcoord->x2f(j);
        Real x3 = pcoord->x3f(k);
        Real SN = std::sin(kx*x1+ky*x2+kz*x3);
        Real temp = cf2*SN/k2;
        dax(k,j,i) = temp*(By*kz-Bz*ky);
        day(k,j,i) = temp*(Bz*kx-Bx*kz);
        daz(k,j,i) = temp*(Bx*ky-By*kx);
      }
    }
  }

  // initialize interface B
  for (int k=ks; k<=ke  ; k++) {
    for (int j=js; j<=je  ; j++) {
      for (int i=is; i<=ie+1; i++) {
        pfield->b.x1f(k,j,i) = Bx +
                               0.5*(daz(k,j+1,i)+daz(k+1,j+1,i)
                                    -daz(k,j,i)-daz(k+1,j,i))/pcoord->dx2f(j) -
                               0.5*(day(k+1,j,i)+day(k+1,j+1,i)
                                    -day(k,j,i)-day(k,j+1,i))/pcoord->dx3f(k);
      }
    }
  }
  for (int k=ks; k<=ke  ; k++) {
    for (int j=js; j<=je+1; j++) {
      for (int i=is; i<=ie  ; i++) {
        pfield->b.x2f(k,j,i) = By +
                               0.5*(dax(k+1,j,i)+dax(k+1,j,i+1)
                                    -dax(k,j,i)-dax(k,j,i+1))/pcoord->dx3f(k) -
                               0.5*(daz(k,j,i+1)+daz(k+1,j,i+1)
                                    -daz(k,j,i)-daz(k+1,j,i))/pcoord->dx1f(i);
      }
    }
  }

  for (int k=ks; k<=ke+1; k++) {
    for (int j=js; j<=je  ; j++) {
      for (int i=is; i<=ie  ; i++) {
        pfield->b.x3f(k,j,i) = 0.5*(day(k,j,i+1)+day(k,j+1,i+1)
                                    -day(k,j,i)-day(k,j+1,i))/pcoord->dx1f(i) -
                               0.5*(dax(k,j+1,i)+dax(k,j+1,i+1)
                                    -dax(k,j,i)-dax(k,j,i+1))/pcoord->dx2f(j);
      }
    }
  }
  dax.DeleteAthenaArray();
  day.DeleteAthenaArray();
  daz.DeleteAthenaArray();

  // initialize total energy
  if (NON_BAROTROPIC_EOS) {
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          phydro->u(IEN,k,j,i) =
            p0/gm1 +
            0.5*(SQR(0.5*(pfield->b.x1f(k,j,i) + pfield->b.x1f(k,j,i+1)))  +
                 SQR(0.5*(pfield->b.x2f(k,j,i) + pfield->b.x2f(k,j+1,i)))  +
                 SQR(0.5*(pfield->b.x3f(k,j,i) + pfield->b.x3f(k+1,j,i)))) +
            0.5*(SQR(phydro->u(IM1,k,j,i)) +
                 SQR(phydro->u(IM2,k,j,i)) +
                 SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
        }
      }
    }
  }

  return;
}

namespace {

Real HistoryDby(MeshBlock *pmb, int iout) {
  Real kx = (TWO_PI/x1size)*(static_cast<Real>(nwx));
  Real ky = (TWO_PI/x2size)*(static_cast<Real>(nwy));
  Real kz = (TWO_PI/x3size)*(static_cast<Real>(nwz));
  if (pmb->pmy_mesh->mesh_size.nx3 == 1)
    kz = 0.0;
  Real alpha = std::sqrt(d0*SQR(iso_cs)/(beta*(SQR(kx)+SQR(ky))));
  kx += qshear*Omega_0*pmb->pmy_mesh->time*(TWO_PI/static_cast<Real>(x1size));
  Real By  = -alpha*kx;
  Real dby = 0.0;
  AthenaArray<Real> volume; // 1D array of volumes
  volume.NewAthenaArray(pmb->ncells1);
  Real tvol = x1size*x2size*x3size;
  for (int k=pmb->ks; k<=pmb->ke  ; k++) {
    for (int j=pmb->js; j<=pmb->je  ; j++) {
      pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, volume);
      for (int i=pmb->is; i<=pmb->ie  ; i++) {
        Real x1 = pmb->pcoord->x1v(i);
        Real x2 = pmb->pcoord->x2v(j);
        Real x3 = pmb->pcoord->x3v(k);
        Real CS = std::cos(kx*x1+ky*x2+kz*x3);
        dby += volume(i)*2.0*(pmb->pfield->bcc(IB2,k,j,i)-By)*CS/tvol;
      }
    }
  }
  return dby;
}
} // namespace
