//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hb3.cpp
//! \brief Problem generator for 2D MRI simulations using the shearing sheet
//!  based on "A powerful local shear instability in weakly magnetized disks"
//!
//! PURPOSE: Problem generator for 2D MRI simulations using the shearing sheet
//!   based on "A powerful local shear instability in weakly magnetized disks.
//!   III - Long-term evolution in a shearing sheet" by Hawley & Balbus.  This
//!   is the third of the HB papers on the MRI, thus hb3.
//!
//! Several different perturbations and field configurations are possible:
//! - ipert = 1 - isentropic perturbations to P & d [default]
//! - ipert = 2 - uniform Vx=amp, sinusoidal density
//!
//! - ifield = 1 - Bz=B0 sin(x1) field with zero-net-flux [default]
//! - ifield = 2 - uniform Bz
//!
//! PRIVATE FUNCTION PROTOTYPES:
//! - ran2() - random number generator from NR
//!
//! REFERENCE: Hawley, J. F. & Balbus, S. A., ApJ 400, 595-609 (1992).*/
//======================================================================================

// C headers

// C++ headers
#include <cmath>      // sqrt()
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
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../orbital_advection/orbital_advection.hpp"
#include "../parameter_input.hpp"
#include "../utils/utils.hpp" // ran2()


#if !MAGNETIC_FIELDS_ENABLED
#error "This problem generator requires magnetic fields"
#endif

namespace {
Real amp, nwx, nwy; // amplitude, Wavenumbers
int ShBoxCoord, ipert,ifield; // initial pattern
Real beta, B0;
Real gm1, iso_cs;
Real x1size, x2size, x3size;
Real Omega_0, qshear;
AthenaArray<Real> volume; // 1D array of volumes
Real HistoryBxBy(MeshBlock *pmb, int iout);
} // namespace

//======================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Init the Mesh properties
//======================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // initialize global variables
  amp    = pin->GetReal("problem","amp");
  beta   = pin->GetReal("problem","beta");
  nwx = pin->GetOrAddInteger("problem","nwx",1);
  nwy = pin->GetOrAddInteger("problem","nwy",1);
  ipert  = pin->GetOrAddInteger("problem","ipert",1);
  ifield = pin->GetOrAddInteger("problem","ifield",1);

  if (!shear_periodic) {
    std::stringstream msg;
    msg << "### FATAL ERROR in hb3.cpp ProblemGenerator" << std::endl
        << "This problem generator requires shearing box" << std::endl;
    ATHENA_ERROR(msg);
  }

  // enroll new history variables
  AllocateUserHistoryOutput(1);
  EnrollUserHistoryOutput(0, HistoryBxBy, "<-BxBy>");
  return;
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  if (pmy_mesh->mesh_size.nx2 == 1 || pmy_mesh->mesh_size.nx3 > 1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in hb3.cpp ProblemGenerator" << std::endl
        << "Shearing sheet only works on a 2D grid" << std::endl;
    ATHENA_ERROR(msg);
  }

  if (porb->orbital_advection_defined) {
    std::stringstream msg;
    msg << "### FATAL ERROR in hb3.cpp ProblemGenerator" << std::endl
        << "This problem does NOT work with orbital advection." << std::endl;
    ATHENA_ERROR(msg);
  }

  // shearing sheet parameter
  Omega_0 = porb->Omega0;
  qshear  = porb->qshear;
  ShBoxCoord = porb->shboxcoord;

  if (ShBoxCoord != 2) {
    std::stringstream msg;
    msg << "### FATAL ERROR in hb3.cpp ProblemGenerator" << std::endl
        << "Shearing sheet only works for x-z plane with shboxcoord=2" << std::endl
        << "Check <orbital_advection> shboxcoord parameter." <<std::endl;
    ATHENA_ERROR(msg);
  }

  // allocate 1D array for cell volume used in usr def history
  volume.NewAthenaArray(ncells1);

  Real d0 = 1.0;
  Real p0 = 1e-5;

  if (NON_BAROTROPIC_EOS) {
    gm1 = (peos->GetGamma() - 1.0);
    iso_cs = std::sqrt((gm1+1.0)*p0/d0);
    std::cout << "gamma  = " << peos->GetGamma() << std::endl;
  } else {
    iso_cs = peos->GetIsoSoundSpeed();
    p0 = d0*SQR(iso_cs);
    std::cout << "iso_cs = " << iso_cs << std::endl;
  }

  B0 = std::sqrt(static_cast<Real>(2.0*p0/beta));
  std::cout << "d0     = " << d0     << std::endl;
  std::cout << "p0     = " << p0     << std::endl;
  std::cout << "B0     = " << B0     << std::endl;
  std::cout << "ipert  = " << ipert  << std::endl;
  std::cout << "ifield = " << ifield << std::endl;
  std::cout << "beta   = " << beta   << std::endl;


  x1size = pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min;
  x2size = pmy_mesh->mesh_size.x2max - pmy_mesh->mesh_size.x2min;
  x3size = pmy_mesh->mesh_size.x3max - pmy_mesh->mesh_size.x3min;
  std::cout << "[hb3.cpp]: [Lx,Lz,Ly] = [" <<x1size <<","<<x2size<<","<<x3size<<"]"
            << std::endl;

  Real kx = (TWO_PI/x1size)*(static_cast<Real>(nwx));
  // Real kz = (TWO_PI/x2size)*(static_cast<Real>(nwy));

  Real x1, x2; //, x3;
  Real rd, rp, rval;
  Real rvx, rvy, rvz;
  std::int64_t iseed = -1-gid; // Initialize on the first call to ran2
  // Initialize perturbations
  //   ipert = 1 - isentropic perturbations to P & d [default]
  //   ipert = 2 - uniform Vx=amp, sinusoidal density
  for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      x1 = pcoord->x1v(i);
      x2 = pcoord->x2v(j);
      rd = d0;
      rp = p0;
      rvx = 0.0;
      rvy = 0.0;
      rvz = 0.0;
      if (ipert == 1) {
        rval = 1.0 + amp*(ran2(&iseed) - 0.5);
        if (NON_BAROTROPIC_EOS) {
          rp = rval*p0;
          rd = d0;
        } else {
          rd = rval*d0;
        }
        rvx = 0.0;
      } else if (ipert == 2) {
        rp = p0;
        rd = d0*(1.0+0.1*std::sin(static_cast<Real>(kx)*x1));
        if (NON_BAROTROPIC_EOS) {
          rvx = amp*std::sqrt((gm1+1.0)*p0/d0);
        } else {
          rvx = amp*std::sqrt(p0/d0);
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in hb3.cpp ProblemGenerator" << std::endl
            << "Shearing sheet ipert=" << ipert << " is invalid." << std::endl;
        ATHENA_ERROR(msg);
      }
      phydro->u(IDN,ks,j,i) = rd;
      phydro->u(IM1,ks,j,i) = rd*rvx;
      phydro->u(IM2,ks,j,i) = rd*rvy;
      phydro->u(IM3,ks,j,i) = rd*rvz;
      phydro->u(IM3,ks,j,i) -= rd*qshear*Omega_0*x1;
      if (NON_BAROTROPIC_EOS) {
        phydro->u(IEN,ks,j,i) = rp/gm1 +
                                0.5*(SQR(phydro->u(IM1,ks,j,i)) +
                                     SQR(phydro->u(IM2,ks,j,i)) +
                                     SQR(phydro->u(IM3,ks,j,i)))/rd;
      }

      // Initialize magnetic field.  For 2D shearing box
      // B1=Bx, B2=Bz, B3=By
      // ifield = 1 - Bz=B0 std::sin(x1) field with zero-net-flux [default]
      // ifield = 2 - uniform Bz
      // ifield = 3 - sinusiodal modes (Nordita workshop test)
      if (MAGNETIC_FIELDS_ENABLED) {
        if (ifield == 1) {
          pfield->b.x1f(ks,j,i) = 0.0;
          pfield->b.x2f(ks,j,i) = B0*(std::sin(static_cast<Real>(kx)*x1));
          pfield->b.x3f(ks,j,i) = 0.0;
          if (i==ie) pfield->b.x1f(ks,j,ie+1) = 0.0;
          if (j==je) pfield->b.x2f(ks,je+1,i) = B0*(std::sin(static_cast<Real>(kx)*x1));
        } else if (ifield == 2) {
          pfield->b.x1f(ks,j,i) = 0.0;
          pfield->b.x2f(ks,j,i) = B0;
          pfield->b.x3f(ks,j,i) = 0.0;
          if (i==ie) pfield->b.x1f(ks,j,ie+1) = 0.0;
          if (j==je) pfield->b.x2f(ks,je+1,i) = B0;
        } else {
          std::stringstream msg;
          msg << "### FATAL ERROR in hb3.cpp ProblemGenerator" << std::endl
              << "Shearing sheet ifield=" << ifield << " is unrecognized." << std::endl;
          ATHENA_ERROR(msg);
        }
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,ks,j,i) += 0.5*(
              SQR(0.5*(pfield->b.x1f(ks,j,i) + pfield->b.x1f(ks,j,i+1))) +
              SQR(0.5*(pfield->b.x2f(ks,j,i) + pfield->b.x2f(ks,j+1,i))) +
              SQR(pfield->b.x3f(ks,j,i)));
        }
      }
    }
  }
  return;
}


//======================================================================================
//! \fn void MeshBlock::UserWorkInLoop()
//  \brief User-defined work function for every time step
//======================================================================================
void MeshBlock::UserWorkInLoop() {
  // nothing to do
  return;
}

namespace {
Real HistoryBxBy(MeshBlock *pmb, int iout) {
  Real bxby=0;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  AthenaArray<Real> &b = pmb->pfield->bcc;

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      pmb->pcoord->CellVolume(k,j,pmb->is,pmb->ie,volume);
      for (int i=is; i<=ie; i++) {
        bxby-=volume(i)*b(IB1,k,j,i)*b(IB3,k,j,i);
      }
    }
  }

  return bxby;
}
} // namespace
