//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file ssheet.cpp
//  \brief Shearing wave problem generator for 2D/3D problems.
//  Several different initial conditions:
//  - ipert = 0  pure shearing background flow
//  - ipert = 1  shearing wave perturbed in velocity
//  - ipert = 2  epicycle motion (0.1 c_s initial kick in radial)
//
// Code must be configured using -shear
//======================================================================================

// C headers
#include <stdlib.h>   // exit

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
#include "../parameter_input.hpp"

#if !SHEARING_BOX
#error "This problem generator requires shearing box"
#endif

static Real amp, nwx, nwy; // amplitude, Wavenumbers
static int ipert; // initial pattern
static Real gm1,iso_cs;
static Real x1size,x2size,x3size;
static Real Omega_0,qshear;
static int shboxcoord;
AthenaArray<Real> volume; // 1D array of volumes

//======================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Init the Mesh properties
//======================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // initialize global variables
  amp = pin->GetReal("problem","amp");
  nwx = pin->GetInteger("problem","nwx");
  nwy = pin->GetInteger("problem","nwy");
  ipert = pin->GetInteger("problem","ipert");
  Omega_0 = pin->GetOrAddReal("problem","Omega0",0.001);
  qshear  = pin->GetOrAddReal("problem","qshear",1.5);
  shboxcoord = pin->GetOrAddInteger("problem","shboxcoord",1);

  return;
}



//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Linear wave problem generator for 1D/2D/3D problems.
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  if (pmy_mesh->mesh_size.nx2 == 1 || pmy_mesh->mesh_size.nx3 > 1) {
      std::cout << "[ssheet.cpp]: only works on 2D grid" << std::endl;
      exit(0);
  }

  if (MAGNETIC_FIELDS_ENABLED) {
      std::cout << "[ssheet.cpp]: only works for hydro alone" << std::endl;
      exit(0);
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
  std::cout << "iso_cs = " << iso_cs << std::endl;
  std::cout << "d0 = " << d0 << std::endl;
  std::cout << "p0 = " << p0 << std::endl;
  std::cout << "ipert  = " << ipert  << std::endl;


  x1size = pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min;
  x2size = pmy_mesh->mesh_size.x2max - pmy_mesh->mesh_size.x2min;
  x3size = pmy_mesh->mesh_size.x3max - pmy_mesh->mesh_size.x3min;
  std::cout << "[ssheet.cpp]: [Lx,Ly,Lz] = [" <<x1size <<","<<x2size
            <<","<<x3size<<"]"<<std::endl;

  Real kx = (2.0*PI/x1size)*(static_cast<Real>(nwx));
  Real ky = (2.0*PI/x2size)*(static_cast<Real>(nwy));

  Real x1,x2,rd,rp,rvx,rvy;
// update the physical variables as initial conditions
  //int nx1 = (ie-is)+1 + 2*(NGHOST);
  //int nx2 = (je-js)+1 + 2*(NGHOST);
  //int nx3 = (ke-ks)+1 + 2*(NGHOST);

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      x1 = pcoord->x1v(i);
      x2 = pcoord->x2v(j);
      rd = d0;
      rp = p0;
      if (ipert == 0) {
        // 1) pure shear bg flow:
        phydro->u(IDN,k,j,i) = rd;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) -= rd*(qshear*Omega_0*x1);
        phydro->u(IM3,k,j,i) = 0.0;
      } else if (ipert == 1) {
        // 2) initialize with shwave in velocity
        rvx = amp*iso_cs*sin(kx*x1 + ky*x2);
        rvy = amp*iso_cs*(kx/ky)*sin(kx*x1 + ky*x2);
        phydro->u(IDN,k,j,i) = rd;
        phydro->u(IM1,k,j,i) = rd*rvx;
        phydro->u(IM2,k,j,i) -= rd*(rvy + qshear*Omega_0*x1);
        phydro->u(IM3,k,j,i) = 0.0;
      } else if (ipert == 2) {
        // 3) epicyclic oscillation
        if (shboxcoord == 1) { // x-y shear
          rvx = 0.1*iso_cs;
          rvy = 0.0;
          phydro->u(IDN,k,j,i) = rd;
          phydro->u(IM1,k,j,i) = rd*rvx;
          phydro->u(IM2,k,j,i) -= rd*(rvy + qshear*Omega_0*x1);
          phydro->u(IM3,k,j,i) = 0.0;
        } else { // x-z plane
          rvx = 0.1*iso_cs;
          rvy = 0.0;
          phydro->u(IDN,k,j,i) = rd;
          phydro->u(IM1,k,j,i) = rd*rvx;
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = -rd*(rvy + qshear*Omega_0*x1);
        }
      } else {
          std::cout << "[ssheet.cpp] ipert = " << ipert
                    << " is unrecognized " <<std::endl;
          exit(0);
      }
      if (NON_BAROTROPIC_EOS) {
        phydro->u(IEN,k,j,i) = rp/gm1 + 0.5*(SQR(phydro->u(IM1,k,j,i)) +
                                             SQR(phydro->u(IM2,k,j,i)) +
                                             SQR(phydro->u(IM3,k,j,i))
                                             ) / phydro->u(IDN,k,j,i);
      }
    }
  }}


  return;
}


//======================================================================================
//! \fn void MeshBlock::UserWorkInLoop(void)
//  \brief User-defined work function for every time step
//======================================================================================
void MeshBlock::UserWorkInLoop(void) {
  // nothing to do
  return;
}
