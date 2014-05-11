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

#include <iostream>
#include <string>
#include <math.h>
#include <float.h>

#include "athena.hpp"
#include "athena_arrays.hpp"
#include "parameter_input.hpp"
#include "mesh.hpp"
#include "fluid.hpp"
#include "convert_var/convert_var.hpp"

//======================================================================================
//! \file fluid.cpp
//  \brief implementation of functions in class Fluid
//======================================================================================

// constructor, initializes data structures and parameters, calls problem generator

Fluid::Fluid(ParameterInput *pin, Block *pb)
{
  pmy_block = pb;

// Read some parameters from input file

  gamma_ = pin->GetReal("fluid","gamma");
  cfl_number_ = pin->GetReal("time","cfl_number");

// Allocate memory for primitive/conserved variables

  int ncells1 = pmy_block->block_size.nx1 + 2*(NGHOST);
  int ncells2 = pmy_block->block_size.nx2 + 2*(NGHOST);
  int ncells3 = pmy_block->block_size.nx3 + 2*(NGHOST);

  u.NewAthenaArray(NVAR,ncells3,ncells2,ncells1);
  w.NewAthenaArray(NVAR,ncells3,ncells2,ncells1);

// Allocate memory for primitive/conserved variables at half-time step, and scratch

  u1.NewAthenaArray(NVAR,ncells3,ncells2,ncells1);
  w1.NewAthenaArray(NVAR,ncells3,ncells2,ncells1);

  wl_.NewAthenaArray(NVAR,ncells1);
  wr_.NewAthenaArray(NVAR,ncells1);
  flx_.NewAthenaArray(NVAR,ncells1);

//  Problem(pin);

}

// destructor

Fluid::~Fluid()
{
  u.DeleteAthenaArray();
  w.DeleteAthenaArray();
  u1.DeleteAthenaArray();
  w1.DeleteAthenaArray();
  
  wl_.DeleteAthenaArray();
  wr_.DeleteAthenaArray();
  flx_.DeleteAthenaArray();
}

//--------------------------------------------------------------------------------------
// \!fn 
// \brief

void Fluid::NewTimeStep(Block *pb)
{
  Real min_dt1=(FLT_MAX), min_dt2=(FLT_MAX), min_dt3=(FLT_MAX);
  int is = pb->is; int js = pb->js; int ks = pb->ks;
  int ie = pb->ie; int je = pb->je; int ke = pb->ke;
  Real gam = GetGamma();

  AthenaArray<Real> w = pb->pfluid->w.ShallowCopy();

  for (int k=ks; k<=ke; ++k){
  for (int j=js; j<=je; ++j){
    Real& dx2 = pb->dx2f(j);
    Real& dx3 = pb->dx3f(k);
#pragma simd
    for (int i=is; i<=ie; ++i){
      Real& w_d  = w(IDN,k,j,i);
      Real& w_v1 = w(IVX,k,j,i);
      Real& w_v2 = w(IVY,k,j,i);
      Real& w_v3 = w(IVZ,k,j,i);
      Real& w_p  = w(IEN,k,j,i);
      Real& dx1  = pb->dx1f(i);

      Real cs = sqrt(gam*w_p/((gam-1.0)*w_d));
      min_dt1 = std::max(min_dt1,(dx1/(fabs(w_v1) + cs)));
      if (pb->block_size.nx2 > 1) min_dt2 = std::min(min_dt2,(dx2/(fabs(w_v2) + cs)));
      if (pb->block_size.nx3 > 1) min_dt3 = std::min(min_dt3,(dx3/(fabs(w_v3) + cs)));
    }
  }}

  Real min_dt = std::min(min_dt1,min_dt2);
  min_dt = std::min(min_dt ,min_dt3);

  Real old_dt = pb->pmy_domain->pmy_mesh->dt;
  pb->pmy_domain->pmy_mesh->dt = std::min((cfl_number_*min_dt), (2.0*old_dt));

  return;

}
