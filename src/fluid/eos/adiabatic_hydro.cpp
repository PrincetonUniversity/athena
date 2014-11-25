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

// Primary header
#include "eos.hpp"

// C++ headers
#include <cmath>   // sqrt()
#include <cfloat>  // FLT_MIN

// Athena headers
#include "../fluid.hpp"               // Fluid
#include "../../athena.hpp"           // enums, macros, Real
#include "../../athena_arrays.hpp"    // AthenaArray
#include "../../mesh.hpp"             // MeshBlock
#include "../../parameter_input.hpp"  // GetReal()
#include "../../field/field.hpp"      // InterfaceField

//======================================================================================
//! \file adiabatic_hydro.cpp
//  \brief implements functions in class FluidEqnOfState for adiabatic hydrodynamics`
//======================================================================================

// FluidEqnOfState constructor

FluidEqnOfState::FluidEqnOfState(Fluid *pf, ParameterInput *pin)
{
  pmy_fluid_ = pf;
  gamma_ = pin->GetReal("fluid","gamma");
  density_floor_  = pin->GetOrAddReal("fluid","dfloor",(1024*(FLT_MIN)));
  pressure_floor_ = pin->GetOrAddReal("fluid","pfloor",(1024*(FLT_MIN)));
}

// destructor

FluidEqnOfState::~FluidEqnOfState()
{
}

//--------------------------------------------------------------------------------------
// \!fn void FluidEqnOfState::ConservedToPrimitive(const AthenaArray<Real> &cons,
//  const AthenaArray<Real> &prim_old, const InterfaceField &b, AthenaArray<Real> &prim,
//  AthenaArray<Real> &bcc)
// \brief Converts conserved into primitive variables in adiabatic hydro.

void FluidEqnOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
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

// Convert to Primitives

  int max_nthreads = pmb->pmy_mesh->nthreads_mesh;

// #pragma omp parallel default(shared) num_threads(max_nthreads)
{
  for (int k=kl; k<=ku; ++k){
// #pragma omp for schedule(dynamic)
  for (int j=jl; j<=ju; ++j){
#pragma simd
    for (int i=pmb->is-(NGHOST); i<=pmb->ie+(NGHOST); ++i){
      const Real& u_d  = cons(IDN,k,j,i);
      const Real& u_m1 = cons(IVX,k,j,i);
      const Real& u_m2 = cons(IVY,k,j,i);
      const Real& u_m3 = cons(IVZ,k,j,i);
      const Real& u_e  = cons(IEN,k,j,i);

// apply density floor, without changing momentum or energy
      cons(IDN,k,j,i) = std::max(cons(IDN,k,j,i), density_floor_);
      prim(IDN,k,j,i) = u_d;

      Real di = 1.0/u_d;
      prim(IVX,k,j,i) = u_m1*di;
      prim(IVY,k,j,i) = u_m2*di;
      prim(IVZ,k,j,i) = u_m3*di;

      prim(IEN,k,j,i) = gm1*( u_e - 0.5*di*(SQR(u_m1) + SQR(u_m2) + SQR(u_m3)) );

// apply pressure floor, correct total energy
      if (prim(IEN,k,j,i) < pressure_floor_) {
        prim(IEN,k,j,i) = pressure_floor_;
        cons(IEN,k,j,i) = (pressure_floor_/gm1) +0.5*di*(SQR(u_m1)+SQR(u_m2)+SQR(u_m3));
      }
    }
  }}
}

  return;
}

//--------------------------------------------------------------------------------------
// \!fn Real FluidEqnOfState::SoundSpeed(Real prim[NFLUID])
// \brief returns adiabatic sound speed given vector of primitive variables

Real FluidEqnOfState::SoundSpeed(const Real prim[NFLUID])
{
  return sqrt(GetGamma()*prim[IEN]/prim[IDN]);
}
