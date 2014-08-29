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

// Primary header
#include "../fluid.hpp"

// Athena headers
#include "../athena.hpp"         // enums, macros, Real
#include "../athena_arrays.hpp"  // AthenaArray
#include "../mesh.hpp"           // MeshBlock

//======================================================================================
/*! \file adiabatic_hydro.cpp
 *  \brief converts conserved to primitive variables in adiabatic hydrodynamics`
 *====================================================================================*/

void Fluid::ConservedToPrimitive(AthenaArray<Real> &c, AthenaArray<Real> &p_old,
    AthenaArray<Real> &p)
{
  MeshBlock *pb = pmy_block;
  int is = pb->is; int ie = pb->ie;
  int jl = pb->js; int ju = pb->je;
  int kl = pb->ks; int ku = pb->ke;
  if (pb->block_size.nx2 > 1) {
    jl -= (NGHOST);
    ju += (NGHOST);
  }
  if (pb->block_size.nx3 > 1) {
    kl -= (NGHOST);
    ku += (NGHOST);
  }
  Real gamma = pb->pfluid->GetGamma();

  AthenaArray<Real> lc = c.ShallowCopy();
  AthenaArray<Real> lp = p.ShallowCopy();

//--------------------------------------------------------------------------------------
// Convert to Primitives

  for (int k=kl; k<=ku; ++k){
  for (int j=jl; j<=ju; ++j){
#pragma simd
    for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i){
      Real& u_d  = lc(IDN,k,j,i);
      Real& u_m1 = lc(IVX,k,j,i);
      Real& u_m2 = lc(IVY,k,j,i);
      Real& u_m3 = lc(IVZ,k,j,i);
      Real& u_e  = lc(IEN,k,j,i);

      Real& w_d  = lp(IDN,k,j,i);
      Real& w_m1 = lp(IVX,k,j,i);
      Real& w_m2 = lp(IVY,k,j,i);
      Real& w_m3 = lp(IVZ,k,j,i);
      Real& w_p  = lp(IEN,k,j,i);

      Real di = 1.0/u_d;
      w_d  = u_d;
      w_m1 = u_m1*di;
      w_m2 = u_m2*di;
      w_m3 = u_m3*di;

      w_p = u_e - 0.5*di*(u_m1*u_m1 + u_m2*u_m2 + u_m3*u_m3);
      w_p *= (gamma - 1.0);
    }
  }}

  return;
}
 
