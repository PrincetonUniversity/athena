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
 * You should have received a copy of GNU GPL in the file LICENSE included in
 * the code distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

#include <sstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <math.h>

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh.hpp"
#include "../fluid.hpp"

//======================================================================================
/*! \file shock_tube.cpp
 *  \brief Problem generator for shock tube problems.  
 *
 * Problem generator for shock tube (1-D Riemann) problems. Initializes plane-parallel
 * shock along x1 (in 1D, 2D, 3D), along x2 (in 2D, 3D), and along x3 (in 3D).
 *====================================================================================*/

void Fluid::InitProblem(ParameterInput *pin)
{
  Block *pb = pparent_block;
  std::stringstream msg;

  int is = pb->is; int js = pb->js; int ks = pb->ks;
  int ie = pb->ie; int je = pb->je; int ke = pb->ke;

// Read parameters from input file

  gamma_ = pin->GetReal("fluid","gamma");
  Real gm1 = (GetGamma()) - 1.0;

// shock direction: {1,2,3} -> {x1,x2,x3}

  int shk_dir = pin->GetInteger("problem","shock_dir"); 

// shock location (must be inside grid)

  Real xshock = pin->GetReal("problem","xshock"); 
  if (shk_dir == 1 && (xshock < pb->pparent_domain->pparent_mesh->mesh_size.x1min ||
                       xshock > pb->pparent_domain->pparent_mesh->mesh_size.x1max)) {
    msg << "### FATAL ERROR in Problem Generator" << std::endl << "xshock="
        << xshock << " lies outside x1 domain for shkdir=" << shk_dir << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (shk_dir == 2 && (xshock < pb->pparent_domain->pparent_mesh->mesh_size.x2min ||
                       xshock > pb->pparent_domain->pparent_mesh->mesh_size.x2max)) {
    msg << "### FATAL ERROR in Problem Generator" << std::endl << "xshock="
        << xshock << " lies outside x2 domain for shkdir=" << shk_dir << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (shk_dir == 3 && (xshock < pb->pparent_domain->pparent_mesh->mesh_size.x3min ||
                       xshock > pb->pparent_domain->pparent_mesh->mesh_size.x3max)) {
    msg << "### FATAL ERROR in Problem Generator" << std::endl << "xshock="
        << xshock << " lies outside x3 domain for shkdir=" << shk_dir << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

// Parse left state read from input file: dl,pl,ul,vl,wl

  Real dl = pin->GetReal("problem","dl");
  Real pl = pin->GetReal("problem","pl");
  Real ul = pin->GetReal("problem","ul");
  Real vl = pin->GetReal("problem","vl");
  Real wl = pin->GetReal("problem","wl");

// Parse right state read from input file: dr,pr,ur,vr,wr

  Real dr = pin->GetReal("problem","dr");
  Real pr = pin->GetReal("problem","pr");
  Real ur = pin->GetReal("problem","ur");
  Real vr = pin->GetReal("problem","vr");
  Real wr = pin->GetReal("problem","wr");

// Initialize the discontinuity

  switch(shk_dir) {

//------ shock in 1-direction ----------------------------------------------------------
  case 1:
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        if (pb->x1v(i) < xshock) {
          u(IDN,k,j,i) = dl;
          u(IM1,k,j,i) = ul*dl;
          u(IM2,k,j,i) = vl*dl;
          u(IM3,k,j,i) = wl*dl;
          u(IEN,k,j,i) = pl/gm1 + 0.5*dl*(ul*ul + vl*vl + wl*wl);
        } else {
          u(IDN,k,j,i) = dr;
          u(IM1,k,j,i) = ur*dr;
          u(IM2,k,j,i) = vr*dr;
          u(IM3,k,j,i) = wr*dr;
          u(IEN,k,j,i) = pr/gm1 + 0.5*dr*(ur*ur + vr*vr + wr*wr);
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
          u(IDN,k,j,i) = dl;
          u(IM1,k,j,i) = wl*dl;
          u(IM2,k,j,i) = ul*dl;
          u(IM3,k,j,i) = vl*dl;
          u(IEN,k,j,i) = pl/gm1 + 0.5*dl*(ul*ul + vl*vl + wl*wl);
        }
      } else {
        for (int i=is; i<=ie; ++i) {
          u(IDN,k,j,i) = dr;
          u(IM1,k,j,i) = wr*dr;
          u(IM2,k,j,i) = ur*dr;
          u(IM3,k,j,i) = vr*dr;
          u(IEN,k,j,i) = pr/gm1 + 0.5*dr*(ur*ur + vr*vr + wr*wr);
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
          u(IDN,k,j,i) = dl;
          u(IM1,k,j,i) = vl*dl;
          u(IM2,k,j,i) = wl*dl;
          u(IM3,k,j,i) = ul*dl;
          u(IEN,k,j,i) = pl/gm1 + 0.5*dl*(ul*ul + vl*vl + wl*wl);
        }}
      } else {
        for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          u(IDN,k,j,i) = dr;
          u(IM1,k,j,i) = vr*dr;
          u(IM2,k,j,i) = wr*dr;
          u(IM3,k,j,i) = ur*dr;
          u(IEN,k,j,i) = pr/gm1 + 0.5*dr*(ur*ur + vr*vr + wr*wr);
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
