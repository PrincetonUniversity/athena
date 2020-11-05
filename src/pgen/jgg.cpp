//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file jgg.cpp
//
//  \brief Problem generator for linear MHD in shearing sheet.
//
// PURPOSE:  Problem generator for linear MHD in shearing sheet. Based on the initial
//   conditions described in "ORBITAL ADVECTION BY INTERPOLATION: A FAST AND ACCURATE
//   NUMERICAL SCHEME FOR SUPER-FAST MHD FLOWS" by Johnson, Guan, & Gammie, or JGG.
//
// Two kinds of perturbations are possible:
//
//- ipert = 1 - MHD simple shwave test of JGG -- their figures 5 - 7
//- ipert = 2 - MHD compressive shwave test of JGG -- their figure 11
//
// REFERENCE: Johnson, Guan, & Gammie, ApJS, 177, 373 (2008)
//============================================================================

// C headers

// C++ headers
#include <cmath>      // sqrt()
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
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../orbital_advection/orbital_advection.hpp"
#include "../parameter_input.hpp"

#if !MAGNETIC_FIELDS_ENABLED
#error "This problem generator requires magnetic fields"
#endif

namespace {
Real iso_cs, gm1, d0, p0;
Real nwx, nwy, nwz; // Wavenumbers
Real Lx, Ly, Lz; // root grid size, global to share with output functions
Real Omega_0, qshear;
Real amp, beta;
int ipert, shboxcoord;

Real HistoryBxAmp(MeshBlock *pmb, int iout);
Real HistoryByAmp(MeshBlock *pmb, int iout);
Real HistoryBzAmp(MeshBlock *pmb, int iout);
Real HistorydBy(MeshBlock *pmb, int iout);
} // namespace

// ===================================================================================
void Mesh::InitUserMeshData(ParameterInput *pin) {
  ipert = pin->GetOrAddInteger("problem","ipert", 1);
  if (ipert == 1) {
    if (mesh_size.nx3 > 1) { // 3D
      AllocateUserHistoryOutput(3);
      EnrollUserHistoryOutput(0, HistoryBxAmp, "BxAmp");
      EnrollUserHistoryOutput(1, HistoryByAmp, "ByAmp");
      EnrollUserHistoryOutput(2, HistoryBzAmp, "BzAmp");
    } else { // 2D
      AllocateUserHistoryOutput(2);
      EnrollUserHistoryOutput(0, HistoryBxAmp, "BxAmp");
      EnrollUserHistoryOutput(1, HistoryByAmp, "ByAmp");
    }
  } else if (ipert == 2) {
    AllocateUserHistoryOutput(1);
    EnrollUserHistoryOutput(0, HistorydBy, "dBy");
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in jgg.cpp ProblemGenerator" << std::endl
        << "Parameter problem/ipert should be chosen from 1  or 2." << std::endl;
    ATHENA_ERROR(msg);
  }

  if (!shear_periodic) {
    std::stringstream msg;
    msg << "### FATAL ERROR in jgg.cpp ProblemGenerator" << std::endl
        << "This problem generator requires shearing box" << std::endl;
    ATHENA_ERROR(msg);
  }

  if (mesh_size.nx2 == 1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in jgg.cpp ProblemGenerator" << std::endl
        << "This problem generator works only in 2D or 3D." << std::endl;
    ATHENA_ERROR(msg);
  }
  return;
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief linear MHD waves in shearing box
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // gamma, press, sound speed
  Real gamma = 1.0;
  d0    = pin->GetOrAddReal("problem","d0", 1.0);
  if (NON_BAROTROPIC_EOS) {
    p0    = pin->GetReal("problem","p0");
    gamma  = peos->GetGamma();
    iso_cs = std::sqrt(gamma*p0/d0);
  } else {
    iso_cs = peos->GetIsoSoundSpeed();
    p0 = d0*SQR(iso_cs);
  }

  // shearing box parameter
  shboxcoord = pin->GetOrAddInteger("problem","shboxcoord",1);
  if (shboxcoord != 1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in jgg.cpp ProblemGenerator" << std::endl
        << "This problem generator requires shearing box in x-y plane." << std::endl;
    ATHENA_ERROR(msg);
  }
  Omega_0 = porb->Omega0;
  qshear  = porb->qshear;

  // Read problem parameters for initial conditions
  amp     = pin->GetReal("problem","amp");

  // Initialize boxsize
  Lx = pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min;
  Ly = pmy_mesh->mesh_size.x2max - pmy_mesh->mesh_size.x2min;
  Lz = pmy_mesh->mesh_size.x3max - pmy_mesh->mesh_size.x3min;

  // initialize wavenumbers
  nwx = pin->GetOrAddInteger("problem","nwx",1);
  nwy = pin->GetOrAddInteger("problem","nwy",1);
  if (block_size.nx3 == 1)
    nwz = pin->GetOrAddInteger("problem","nwz",1);
  else
    nwz = pin->GetOrAddInteger("problem","nwz",0);
  Real kx = (TWO_PI/Lx)*(static_cast<Real>(nwx));
  Real ky = (TWO_PI/Ly)*(static_cast<Real>(nwy));
  Real kz = (TWO_PI/Lz)*(static_cast<Real>(nwz));
  if (nwx == 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in jgg.cpp ProblemGenerator" << std::endl
        << "Parameterproblem/nwx must be non-zero." << std::endl;
    ATHENA_ERROR(msg);
  }
  if (block_size.nx3 == 1 && nwz != 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in jgg.cpp ProblemGenerator" << std::endl
        << "In 2D, parameer problem/nwz must vanish." << std::endl;
    ATHENA_ERROR(msg);
  }

  // Initialize perturbations
  // parameters
  Real x1, x2, x3;
  Real rp(0.0);
  Real rbx(0.0), rby(0.0), rbz(0.0);
  // Calculate magnetic fields using the vector potential.
  AthenaArray<Real> rax, ray, raz;
  int nz = (ncells3>1)?ncells3:2;
  rax.NewAthenaArray(nz, ncells2, ncells1);
  ray.NewAthenaArray(nz, ncells2, ncells1);
  raz.NewAthenaArray(nz, ncells2, ncells1);

  // ipert = 1 - MHD linear shwave test of JGG -- their figures 5 - 7
  // ipert = 2 - MHD compressive shwave test of JGG -- their figure 11
  if (ipert == 1) {
    // hydro
    Real rd = d0;
    rp = p0;
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          phydro->u(IDN,k,j,i) = rd;
          phydro->u(IM1,k,j,i) = 0.0;
          phydro->u(IM2,k,j,i) = 0.0;
          if(!porb->orbital_advection_defined)
            phydro->u(IM2,k,j,i) -= rd*qshear*Omega_0*pcoord->x1v(i);
          phydro->u(IM3,k,j,i) = 0.0;
        }
      }
    }
    // vector potential
    rby = amp*ky/std::fabs(ky);
    rbz = amp*kz/std::fabs(ky);
    rbx = -(rby*ky+rbz*kz)/kx;
    // set rax
    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie; i++) {
          rax(k,j,i) = 0.0;
        }
      }
    }
    // set ray
    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie+1; i++) {
          x1 = pcoord->x1f(i);
          x2 = pcoord->x2v(j);
          x3 = pcoord->x3f(k);
          ray(k,j,i) = rbz/kx*std::sin(kx*x1+ky*x2+kz*x3);
        }
      }
    }
    // set raz
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie+1; i++) {
          x1 = pcoord->x1f(i);
          x2 = pcoord->x2f(j);
          x3 = pcoord->x3v(k);
          raz(k,j,i) = -rby/kx*std::sin(kx*x1+ky*x2+kz*x3);
        }
      }
    }
    rbx = 0.0;
    rby = 0.0;
    rbz = 0.0;
  } else { // ipert == 2
    Real rd(0.0);
    rp = p0;
    // In JGG
    // amp (epsilon in JGG)  = 1.0e-6
    // Omega_0 = d0 = iso_cs = H = 1.0
    // Lx = Ly = Lz = 0.5H
    // nwx = -2, nwy = nwz = 1
    beta    = pin->GetReal("problem","beta");
    Real B02  = static_cast<Real>(p0/beta);
    Real k2    =  SQR(kx)+SQR(ky)+SQR(kz);
    Real alpha =  std::sqrt(B02/(kx*kx+ky*ky));
    rbx   =  alpha*ky;
    rby   = -alpha*kx;
    rbz   =  0.0;

    Real sch   =  iso_cs/Omega_0; // H = 1
    Real cf1   =  std::sqrt(B02*(1.0+beta)); // va*sqrt(1+beta)
    Real cf2   =  amp*std::sqrt(sch*std::sqrt(k2*beta/(1.0+beta)));
    Real vd    =  cf1/std::sqrt(k2)*cf2;

    // hydro
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          x1 = pcoord->x1v(i);
          x2 = pcoord->x2v(j);
          x3 = pcoord->x3v(k);
          Real CS = std::cos(kx*x1+ky*x2+kz*x3);
          rd = d0*(1.0+cf2*CS);
          phydro->u(IDN,k,j,i) = rd;
          phydro->u(IM1,k,j,i) = rd*vd*kx*CS;
          phydro->u(IM2,k,j,i) = rd*vd*ky*CS;
          if(!porb->orbital_advection_defined)
            phydro->u(IM2,k,j,i) -= rd*qshear*Omega_0*x1;
          phydro->u(IM3,k,j,i) = rd*vd*kz*CS;
        }
      }
    }

    // vector potential
    // set rax
    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie; i++) {
          x1 = pcoord->x1v(i);
          x2 = pcoord->x2f(j);
          x3 = pcoord->x3f(k);
          Real temp = cf2/k2*std::sin(kx*x1+ky*x2+kz*x3);
          rax(k,j,i) = temp*(rby*kz-rbz*ky);
        }
      }
    }

    // set ray
    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie+1; i++) {
          x1 = pcoord->x1f(i);
          x2 = pcoord->x2v(j);
          x3 = pcoord->x3f(k);
          Real temp = cf2/k2*std::sin(kx*x1+ky*x2+kz*x3);
          ray(k,j,i) = temp*(rbz*kx-rbx*kz);
        }
      }
    }

    // set raz
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie+1; i++) {
          x1 = pcoord->x1f(i);
          x2 = pcoord->x2f(j);
          x3 = pcoord->x3v(k);
          Real temp = cf2/k2*std::sin(kx*x1+ky*x2+kz*x3);
          raz(k,j,i) = temp*(rbx*ky-rby*kx);
        }
      }
    }
  }

  // initialize interface B
  // set bx
  for (int k=ks; k<=ke  ; k++) {
    for (int j=js; j<=je  ; j++) {
      for (int i=is; i<=ie+1; i++) {
        pfield->b.x1f(k,j,i) = rbx
                               +(raz(k,j+1,i)-raz(k,j,i))/pcoord->dx2f(j)
                               -(ray(k+1,j,i)-ray(k,j,i))/pcoord->dx3f(k);
      }
    }
  }

  // set by
  for (int k=ks; k<=ke  ; k++) {
    for (int j=js; j<=je+1; j++) {
      for (int i=is; i<=ie  ; i++) {
        pfield->b.x2f(k,j,i) = rby
                               +(rax(k+1,j,i)-rax(k,j,i))/pcoord->dx3f(k)
                               -(raz(k,j,i+1)-raz(k,j,i))/pcoord->dx1f(i);
      }
    }
  }

  // set bz
  for (int k=ks; k<=ke+1; k++) {
    for (int j=js; j<=je  ; j++) {
      for (int i=is; i<=ie  ; i++) {
        pfield->b.x3f(k,j,i) = rbz
                               +(ray(k,j,i+1)-ray(k,j,i))/pcoord->dx1f(i)
                               -(rax(k,j+1,i)-rax(k,j,i))/pcoord->dx2f(j);
      }
    }
  }
  rax.DeleteAthenaArray();
  ray.DeleteAthenaArray();
  raz.DeleteAthenaArray();

  // initialize total energy
  if (NON_BAROTROPIC_EOS) {
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          phydro->u(IEN,k,j,i) =
              rp/(gamma-1.0) +
              0.5*(SQR(0.5*(pfield->b.x1f(k,j,i) + pfield->b.x1f(k,j,i+1))) +
                   SQR(0.5*(pfield->b.x2f(k,j,i) + pfield->b.x2f(k,j+1,i))) +
                   SQR(0.5*(pfield->b.x3f(k,j,i) + pfield->b.x3f(k+1,j,i)))) + (0.5)*
              (SQR(phydro->u(IM1,k,j,i)) + SQR(phydro->u(IM2,k,j,i))
               + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
        }
      }
    }
  }

  return;
}

namespace {

Real HistoryBxAmp(MeshBlock *pmb, int iout) {
  Real bx_amp = 0.0;
  Real x1, x2, x3;
  AthenaArray<Real> volume; // 1D array of volumes
  volume.NewAthenaArray(pmb->ncells1);
  int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je, ks = pmb->ks, ke = pmb->ke;
  AthenaArray<Real> &b = pmb->pfield->bcc;

  Real kx = (TWO_PI/Lx)
            *(static_cast<Real>(nwx)+qshear*Omega_0*pmb->pmy_mesh->time);
  Real ky = (TWO_PI/Ly)*(static_cast<Real>(nwy));
  Real kz = (TWO_PI/Lz)*(static_cast<Real>(nwz));
  Real total_volume = Lx*Ly*Lz;

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      pmb->pcoord->CellVolume(k, j, is, ie, volume);
      for (int i=is; i<=ie; i++) {
        x1 = pmb->pcoord->x1v(i);
        x2 = pmb->pcoord->x2v(j);
        x3 = pmb->pcoord->x3v(k);
        Real CS = std::cos(kx*x1+ky*x2+kz*x3);
        bx_amp += volume(i)*2.0*b(IB1,k,j,i)*CS/total_volume;
      }
    }
  }
  return bx_amp;
}

Real HistoryByAmp(MeshBlock *pmb, int iout) {
  Real by_amp = 0.0;
  Real x1, x2, x3;
  AthenaArray<Real> volume; // 1D array of volumes
  volume.NewAthenaArray(pmb->ncells1);
  int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je, ks = pmb->ks, ke = pmb->ke;
  AthenaArray<Real> &b = pmb->pfield->bcc;

  Real kx = (TWO_PI/Lx)
            *(static_cast<Real>(nwx)+qshear*Omega_0*pmb->pmy_mesh->time);
  Real ky = (TWO_PI/Ly)*(static_cast<Real>(nwy));
  Real kz = (TWO_PI/Lz)*(static_cast<Real>(nwz));
  Real total_volume = Lx*Ly*Lz;

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      pmb->pcoord->CellVolume(k, j, is, ie, volume);
      for (int i=is; i<=ie; i++) {
        x1 = pmb->pcoord->x1v(i);
        x2 = pmb->pcoord->x2v(j);
        x3 = pmb->pcoord->x3v(k);
        Real CS = std::cos(kx*x1+ky*x2+kz*x3);
        by_amp += volume(i)*2.0*b(IB2,k,j,i)*CS/total_volume;
      }
    }
  }
  return by_amp;
}

Real HistoryBzAmp(MeshBlock *pmb, int iout) {
  Real bz_amp = 0.0;
  Real x1, x2, x3;
  AthenaArray<Real> volume; // 1D array of volumes
  volume.NewAthenaArray(pmb->ncells1);
  int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je, ks = pmb->ks, ke = pmb->ke;
  AthenaArray<Real> &b = pmb->pfield->bcc;

  Real kx = (TWO_PI/Lx)
            *(static_cast<Real>(nwx)+qshear*Omega_0*pmb->pmy_mesh->time);
  Real ky = (TWO_PI/Ly)*(static_cast<Real>(nwy));
  Real kz = (TWO_PI/Lz)*(static_cast<Real>(nwz));
  Real total_volume = Lx*Ly*Lz;

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      pmb->pcoord->CellVolume(k, j, is, ie, volume);
      for (int i=is; i<=ie; i++) {
        x1 = pmb->pcoord->x1v(i);
        x2 = pmb->pcoord->x2v(j);
        x3 = pmb->pcoord->x3v(k);
        Real CS = std::cos(kx*x1+ky*x2+kz*x3);
        bz_amp += volume(i)*2.0*b(IB3,k,j,i)*CS/total_volume;
      }
    }
  }
  return bz_amp;
}

Real HistorydBy(MeshBlock *pmb, int iout) {
  Real dby = 0.0;
  Real x1, x2, x3;
  AthenaArray<Real> volume; // 1D array of volumes
  volume.NewAthenaArray(pmb->ncells1);
  int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je, ks = pmb->ks, ke = pmb->ke;
  AthenaArray<Real> &b = pmb->pfield->bcc;

  Real kx = (TWO_PI/Lx)
            *(static_cast<Real>(nwx)+qshear*Omega_0*pmb->pmy_mesh->time);
  Real ky = (TWO_PI/Ly)*(static_cast<Real>(nwy));
  Real kz = (TWO_PI/Lz)*(static_cast<Real>(nwz));
  Real by_analytic = -kx*std::sqrt(d0*SQR(iso_cs)/(beta*(SQR(kx)+SQR(ky))));
  Real total_volume = Lx*Ly*Lz;

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      pmb->pcoord->CellVolume(k, j, is, ie, volume);
      for (int i=is; i<=ie; i++) {
        x1 = pmb->pcoord->x1v(i);
        x2 = pmb->pcoord->x2v(j);
        x3 = pmb->pcoord->x3v(k);
        Real CS = std::cos(kx*x1+ky*x2+kz*x3);
        dby += volume(i)*2.0*(b(IB2,k,j,i)-by_analytic)*CS/total_volume;
      }
    }
  }
  return dby;
}
} // namespace
