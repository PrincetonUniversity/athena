//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file adiabatic_hydro.cpp
//  \brief implements functions in class HydroEqnOfState for adiabatic hydrodynamics`
//======================================================================================

// C/C++ headers
#include <cmath>   // sqrt()
#include <cfloat>  // FLT_MIN

// Athena++ headers
#include "../fluid.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh.hpp"
#include "../../parameter_input.hpp"
#include "../../field/field.hpp"

// this class header
#include "eos.hpp"

// HydroEqnOfState constructor

HydroEqnOfState::HydroEqnOfState(Hydro *pf, ParameterInput *pin)
{
  pmy_fluid_ = pf;
  gamma_ = pin->GetReal("fluid","gamma");
  density_floor_  = pin->GetOrAddReal("fluid","dfloor",(1024*(FLT_MIN)));
  pressure_floor_ = pin->GetOrAddReal("fluid","pfloor",(1024*(FLT_MIN)));
}

// destructor

HydroEqnOfState::~HydroEqnOfState()
{
}

//--------------------------------------------------------------------------------------
// \!fn void HydroEqnOfState::ConservedToPrimitive(const AthenaArray<Real> &cons,
//  const AthenaArray<Real> &prim_old, const InterfaceField &b, AthenaArray<Real> &prim,
//  AthenaArray<Real> &bcc)
// \brief Converts conserved into primitive variables in adiabatic hydro.

void HydroEqnOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
  const AthenaArray<Real> &prim_old, const InterfaceField &b, AthenaArray<Real> &prim,
  AthenaArray<Real> &bcc)
{
  MeshBlock *pmb = pmy_fluid_->pmy_block;
  int jl = pmb->js; int ju = pmb->je;
  int kl = pmb->ks; int ku = pmb->ke;
  if (pmb->block_size.nx2 > 1) {
    jl -= (NGHOST);
    ju += (NGHOST);
  }
  if (pmb->block_size.nx3 > 1) {
    kl -= (NGHOST);
    ku += (NGHOST);
  }
  Real gm1 = GetGamma() - 1.0;

  int nthreads = pmb->pmy_mesh->GetNumMeshThreads();
#pragma omp parallel default(shared) num_threads(nthreads)
{
  for (int k=kl; k<=ku; ++k){
#pragma omp for schedule(dynamic)
  for (int j=jl; j<=ju; ++j){
#pragma simd
    for (int i=pmb->is-(NGHOST); i<=pmb->ie+(NGHOST); ++i){
      Real& u_d  = cons(IDN,k,j,i);
      Real& u_m1 = cons(IM1,k,j,i);
      Real& u_m2 = cons(IM2,k,j,i);
      Real& u_m3 = cons(IM3,k,j,i);
      Real& u_e  = cons(IEN,k,j,i);

      Real& w_d  = prim(IDN,k,j,i);
      Real& w_vx = prim(IVX,k,j,i);
      Real& w_vy = prim(IVY,k,j,i);
      Real& w_vz = prim(IVZ,k,j,i);
      Real& w_p  = prim(IEN,k,j,i);

      // apply density floor, without changing momentum or energy
      u_d = (u_d > density_floor_) ?  u_d : density_floor_;
      w_d = u_d;

      Real di = 1.0/u_d;
      w_vx = u_m1*di;
      w_vy = u_m2*di;
      w_vz = u_m3*di;

      Real ke = 0.5*di*(SQR(u_m1) + SQR(u_m2) + SQR(u_m3));
      w_p = gm1*(u_e - ke);

      // apply pressure floor, correct total energy
      u_e = (w_p > pressure_floor_) ?  u_e : ((pressure_floor_/gm1) + ke);
      w_p = (w_p > pressure_floor_) ?  w_p : pressure_floor_;
    }
  }}
}

  return;
}

//--------------------------------------------------------------------------------------
// \!fn Real HydroEqnOfState::SoundSpeed(Real prim[NFLUID])
// \brief returns adiabatic sound speed given vector of primitive variables

Real HydroEqnOfState::SoundSpeed(const Real prim[NFLUID])
{
  return sqrt(gamma_*prim[IEN]/prim[IDN]);
}
