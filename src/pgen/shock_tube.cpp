//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 *
 * This program is free software: you can redistribute and/or modify it under the terms
 * of the GNU General Public License (GPL) as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of GNU GPL in the file LICENSE included in the code
 * distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

// Primary header
#include "../fluid/fluid.hpp"

// C++ headers
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena headers
#include "../athena.hpp"           // enums, Real
#include "../athena_arrays.hpp"    // AthenaArray
#include "../mesh.hpp"             // MeshBlock
#include "../parameter_input.hpp"  // ParameterInput
#include "../fluid/eos/eos.hpp"    // ParameterInput

//======================================================================================
/*! \file shock_tube.cpp
 *  \brief Problem generator for shock tube problems.  
 *
 * Problem generator for shock tube (1-D Riemann) problems. Initializes plane-parallel
 * shock along x1 (in 1D, 2D, 3D), along x2 (in 2D, 3D), and along x3 (in 3D).
 *====================================================================================*/

void Fluid::InitFluid(ParameterInput *pin)
{
  MeshBlock *pb = pmy_block;
  std::stringstream msg;

  int is = pb->is; int js = pb->js; int ks = pb->ks;
  int ie = pb->ie; int je = pb->je; int ke = pb->ke;

// parse shock direction: {1,2,3} -> {x1,x2,x3}

  int shk_dir = pin->GetInteger("problem","shock_dir"); 

// parse shock location (must be inside grid)

  Real xshock = pin->GetReal("problem","xshock"); 
  if (shk_dir == 1 && (xshock < pb->pmy_domain->pmy_mesh->mesh_size.x1min ||
                       xshock > pb->pmy_domain->pmy_mesh->mesh_size.x1max)) {
    msg << "### FATAL ERROR in Problem Generator" << std::endl << "xshock="
        << xshock << " lies outside x1 domain for shkdir=" << shk_dir << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (shk_dir == 2 && (xshock < pb->pmy_domain->pmy_mesh->mesh_size.x2min ||
                       xshock > pb->pmy_domain->pmy_mesh->mesh_size.x2max)) {
    msg << "### FATAL ERROR in Problem Generator" << std::endl << "xshock="
        << xshock << " lies outside x2 domain for shkdir=" << shk_dir << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (shk_dir == 3 && (xshock < pb->pmy_domain->pmy_mesh->mesh_size.x3min ||
                       xshock > pb->pmy_domain->pmy_mesh->mesh_size.x3max)) {
    msg << "### FATAL ERROR in Problem Generator" << std::endl << "xshock="
        << xshock << " lies outside x3 domain for shkdir=" << shk_dir << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

// Parse left state read from input file: dl,ul,vl,wl,[pl]

  Real wl[NVAR];
  wl[IDN] = pin->GetReal("problem","dl");
  wl[IVX] = pin->GetReal("problem","ul");
  wl[IVY] = pin->GetReal("problem","vl");
  wl[IVZ] = pin->GetReal("problem","wl");
  if (NON_BAROTROPIC_EOS) wl[IEN] = pin->GetReal("problem","pl");

// Parse right state read from input file: dr,ur,vr,wr,[pr]

  Real wr[NVAR];
  wr[IDN] = pin->GetReal("problem","dr");
  wr[IVX] = pin->GetReal("problem","ur");
  wr[IVY] = pin->GetReal("problem","vr");
  wr[IVZ] = pin->GetReal("problem","wr");
  if (NON_BAROTROPIC_EOS) wr[IEN] = pin->GetReal("problem","pr");

// Initialize the discontinuity

  switch(shk_dir) {

//------ shock in 1-direction ----------------------------------------------------------
  case 1:
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        if (pb->x1v(i) < xshock) {
          u(IDN,k,j,i) = wl[IDN];
          u(IM1,k,j,i) = wl[IVX]*wl[IDN];
          u(IM2,k,j,i) = wl[IVY]*wl[IDN];
          u(IM3,k,j,i) = wl[IVZ]*wl[IDN];
          if (NON_BAROTROPIC_EOS) u(IEN,k,j,i)= wl[IEN]/(pf_eos->GetGamma() - 1.0) + 
            0.5*wl[IDN]*(wl[IVX]*wl[IVX] + wl[IVY]*wl[IVY] + wl[IVZ]*wl[IVZ]);
        } else {
          u(IDN,k,j,i) = wr[IDN];
          u(IM1,k,j,i) = wr[IVX]*wr[IDN];
          u(IM2,k,j,i) = wr[IVY]*wr[IDN];
          u(IM3,k,j,i) = wr[IVZ]*wr[IDN];
          if (NON_BAROTROPIC_EOS) u(IEN,k,j,i)= wr[IEN]/(pf_eos->GetGamma() - 1.0) + 
            0.5*wr[IDN]*(wr[IVX]*wr[IVX] + wr[IVY]*wr[IVY] + wr[IVZ]*wr[IVZ]);
        }
      }
    }}
    break;

//------ shock in 2-direction ----------------------------------------------------------
  case 2:
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      if (pb->x2v(j) < xshock) {
        for (int i=is; i<=ie; ++i) {
          u(IDN,k,j,i) = wl[IDN];
          u(IM2,k,j,i) = wl[IVX]*wl[IDN];
          u(IM3,k,j,i) = wl[IVY]*wl[IDN];
          u(IM1,k,j,i) = wl[IVZ]*wl[IDN];
          if (NON_BAROTROPIC_EOS) u(IEN,k,j,i)= wl[IEN]/(pf_eos->GetGamma() - 1.0) + 
            0.5*wl[IDN]*(wl[IVX]*wl[IVX] + wl[IVY]*wl[IVY] + wl[IVZ]*wl[IVZ]);
        }
      } else {
        for (int i=is; i<=ie; ++i) {
          u(IDN,k,j,i) = wr[IDN];
          u(IM2,k,j,i) = wr[IVX]*wr[IDN];
          u(IM3,k,j,i) = wr[IVY]*wr[IDN];
          u(IM1,k,j,i) = wr[IVZ]*wr[IDN];
          if (NON_BAROTROPIC_EOS) u(IEN,k,j,i)= wr[IEN]/(pf_eos->GetGamma() - 1.0) + 
            0.5*wr[IDN]*(wr[IVX]*wr[IVX] + wr[IVY]*wr[IVY] + wr[IVZ]*wr[IVZ]);
        }
      }
    }}
    break;

//------ shock in 3-direction ----------------------------------------------------------

  case 3:
    for (int k=ks; k<=ke; ++k) {
      if (pb->x3v(k) < xshock) {
        for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          u(IDN,k,j,i) = wl[IDN];
          u(IM3,k,j,i) = wl[IVX]*wl[IDN];
          u(IM1,k,j,i) = wl[IVY]*wl[IDN];
          u(IM2,k,j,i) = wl[IVZ]*wl[IDN];
          if (NON_BAROTROPIC_EOS) u(IEN,k,j,i)= wl[IEN]/(pf_eos->GetGamma() - 1.0) + 
            0.5*wl[IDN]*(wl[IVX]*wl[IVX] + wl[IVY]*wl[IVY] + wl[IVZ]*wl[IVZ]);
        }}
      } else {
        for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          u(IDN,k,j,i) = wr[IDN];
          u(IM3,k,j,i) = wr[IVX]*wr[IDN];
          u(IM1,k,j,i) = wr[IVY]*wr[IDN];
          u(IM2,k,j,i) = wr[IVZ]*wr[IDN];
          if (NON_BAROTROPIC_EOS) u(IEN,k,j,i)= wr[IEN]/(pf_eos->GetGamma() - 1.0) + 
            0.5*wr[IDN]*(wr[IVX]*wr[IVX] + wr[IVY]*wr[IVY] + wr[IVZ]*wr[IVZ]);
        }}
      }
    }
    break;

  default:
    msg << "### FATAL ERROR in Problem Generator" << std::endl
        << "shock_dir=" << shk_dir << " must be either 1,2, or 3" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  return;
}
