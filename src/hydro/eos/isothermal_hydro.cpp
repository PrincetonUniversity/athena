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
//! \file isothermal_hydro.cpp
//  \brief implements functions in class HydroEqnOfState for isothermal hydrodynamics`
//======================================================================================

// C/C++ headers
#include <cfloat>  // FLT_MIN

// Athena++ headers
#include "../hydro.hpp"
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
  pmy_hydro_ = pf;
  iso_sound_speed_ = pin->GetReal("hydro","iso_sound_speed"); // error if missing!
  density_floor_  = pin->GetOrAddReal("hydro","dfloor",(1024*(FLT_MIN)));
}

// destructor

HydroEqnOfState::~HydroEqnOfState()
{
}

//--------------------------------------------------------------------------------------
// \!fn void HydroEqnOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
//           const AthenaArray<Real> &prim_old, const FaceField &b,
//           AthenaArray<Real> &prim, AthenaArray<Real> &bcc, Coordinates *pco,
//           int is, int ie, int js, int je, int ks, int ke)
// \brief Converts conserved into primitive variables in adiabatic hydro.

void HydroEqnOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
  const AthenaArray<Real> &prim_old, const FaceField &b, AthenaArray<Real> &prim,
  AthenaArray<Real> &bcc, Coordinates *pco, int is, int ie, int js, int je, int ks, int ke)
{
  MeshBlock *pmb = pmy_hydro_->pmy_block;

  int nthreads = pmb->pmy_mesh->GetNumMeshThreads();
#pragma omp parallel default(shared) num_threads(nthreads)
{
  for (int k=ks; k<=ke; ++k){
#pragma omp for schedule(dynamic)
  for (int j=js; j<=je; ++j){
#pragma simd
    for (int i=is; i<=ie; ++i){
      Real& u_d  = cons(IDN,k,j,i);
      Real& u_m1 = cons(IM1,k,j,i);
      Real& u_m2 = cons(IM2,k,j,i);
      Real& u_m3 = cons(IM3,k,j,i);

      Real& w_d  = prim(IDN,k,j,i);
      Real& w_vx = prim(IVX,k,j,i);
      Real& w_vy = prim(IVY,k,j,i);
      Real& w_vz = prim(IVZ,k,j,i);

      // apply density floor, without changing momentum or energy
      u_d = (u_d > density_floor_) ?  u_d : density_floor_;
      w_d = u_d;

      Real di = 1.0/u_d;
      w_vx = u_m1*di;
      w_vy = u_m2*di;
      w_vz = u_m3*di;
    }
  }}
}

  return;
}

//--------------------------------------------------------------------------------------
// \!fn void HydroEqnOfState::PrimitiveToConserved(const AthenaArray<Real> &prim,
//           const AthenaArray<Real> &bc, AthenaArray<Real> &cons, Coordinates *pco,
//           int is, int ie, int js, int je, int ks, int ke);
// \brief Converts primitive variables into conservative variables

void HydroEqnOfState::PrimitiveToConserved(const AthenaArray<Real> &prim,
     const AthenaArray<Real> &bc, AthenaArray<Real> &cons, Coordinates *pco,
     int is, int ie, int js, int je, int ks, int ke)
{
  MeshBlock *pmb = pmy_hydro_->pmy_block;
  Real igm1 = 1.0/(GetGamma() - 1.0);

  int nthreads = pmb->pmy_mesh->GetNumMeshThreads();
#pragma omp parallel default(shared) num_threads(nthreads)
{
  for (int k=ks; k<=ke; ++k){
#pragma omp for schedule(dynamic)
  for (int j=js; j<=je; ++j){
#pragma simd
    for (int i=is; i<=ie; ++i){
      Real& u_d  = cons(IDN,k,j,i);
      Real& u_m1 = cons(IM1,k,j,i);
      Real& u_m2 = cons(IM2,k,j,i);
      Real& u_m3 = cons(IM3,k,j,i);

      const Real& w_d  = prim(IDN,k,j,i);
      const Real& w_vx = prim(IVX,k,j,i);
      const Real& w_vy = prim(IVY,k,j,i);
      const Real& w_vz = prim(IVZ,k,j,i);

      u_d = w_d;
      u_m1 = w_vx*w_d;
      u_m2 = w_vy*w_d;
      u_m3 = w_vz*w_d;
    }
  }}
}
  return;
}

//--------------------------------------------------------------------------------------
// \!fn Real HydroEqnOfState::SoundSpeed(Real dummy_arg[NHYDRO])
// \brief returns isothermal sound speed

Real HydroEqnOfState::SoundSpeed(const Real dummy_arg[NHYDRO])
{
  return iso_sound_speed_;
}
