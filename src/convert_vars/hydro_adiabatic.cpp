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
#include <algorithm>
#include <stdio.h>

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh.hpp"
#include "../fluid.hpp"


//======================================================================================
/*! \file vl2.cpp
 *  \brief van-Leer (MUSCL-Hancock) second-order integrator
 *====================================================================================*/

void Fluid::ConservedToPrimitive(Domain *pd, AthenaArray<Real> &u, AthenaArray<Real> &w)
{
  Block *pb = pd->pblock;
  int is = pb->is; int js = pb->js; int ks = pb->ks;
  int ie = pb->ie; int je = pb->je; int ke = pb->ke;

/*
  AthenaArray<Real> ut = u;
  AthenaArray<Real> wt = w;
*/

//--------------------------------------------------------------------------------------
// Convert to Primitives

  for (int k=ks-2; k<=ke+2; ++k){
  for (int j=js-2; j<=je+2; ++j){
#pragma simd
    for (int i=is-2; i<=ie+2; ++i){
      Real& u_d  = u(IDN,k,j,i);
      Real& u_m1 = u(IVX,k,j,i);
      Real& u_m2 = u(IVY,k,j,i);
      Real& u_m3 = u(IVZ,k,j,i);
      Real& u_e  = u(IEN,k,j,i);

      Real& w_d  = w(IDN,k,j,i);
      Real& w_m1 = w(IVX,k,j,i);
      Real& w_m2 = w(IVY,k,j,i);
      Real& w_m3 = w(IVZ,k,j,i);
      Real& w_e  = w(IEN,k,j,i);

      Real di = 1.0/u_d;
      w_d  = u_d;
      w_m1 = u_m1*di;
      w_m2 = u_m2*di;
      w_m3 = u_m3*di;

      Real qa = u_m1*u_m1;
      qa += u_m2*u_m2;
      qa += u_m3*u_m3;
      w_e = u_e - 0.5*qa*di;
    }
  }}

  return;
}
 
