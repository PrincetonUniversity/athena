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
#include "../../field/field.hpp"      // BFields

//======================================================================================
//! \file isothermal_mhd.cpp
//  \brief implements functions in class FluidEqnOfState for isothermal MHD
//======================================================================================

// FluidEqnOfState constructor

FluidEqnOfState::FluidEqnOfState(Fluid *pf, ParameterInput *pin)
{
  pmy_fluid_ = pf;
  iso_sound_speed_ = pin->GetReal("fluid","iso_sound_speed"); // error if missing!
  density_floor_  = pin->GetOrAddReal("fluid","dfloor",(1024*(FLT_MIN)));
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

  // Convert to Primitives
  for (int k=kl; k<=ku; ++k){
  for (int j=jl; j<=ju; ++j){
#pragma simd
    for (int i=pmb->is-(NGHOST); i<=pmb->ie+(NGHOST); ++i){
      Real& u_d  = cons(IDN,k,j,i);
      Real& u_m1 = cons(IVX,k,j,i);
      Real& u_m2 = cons(IVY,k,j,i);
      Real& u_m3 = cons(IVZ,k,j,i);

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

      const Real& b1_i   = b.x1f(k,j,i  );
      const Real& b1_ip1 = b.x1f(k,j,i+1);
      const Real& b2_j   = b.x2f(k,j  ,i);
      const Real& b2_jp1 = b.x2f(k,j+1,i);
      const Real& b3_k   = b.x3f(k  ,j,i);
      const Real& b3_kp1 = b.x3f(k+1,j,i);

      Real& bcc1 = bcc(IB1,k,j,i);
      Real& bcc2 = bcc(IB2,k,j,i);
      Real& bcc3 = bcc(IB3,k,j,i);

      // cell center B-fields are defined as spatial interpolation at the volume center
      const Real& x1f_i  = pmb->x1f(i);
      const Real& x1f_ip = pmb->x1f(i+1);
      const Real& x1v_i  = pmb->x1v(i);
      const Real& dx1_i  = pmb->dx1f(i);
      Real lw=(x1f_ip-x1v_i)/dx1_i;
      Real rw=(x1v_i -x1f_i)/dx1_i;
      bcc1 = lw*b1_i + rw*b1_ip1;
      const Real& x2f_j  = pmb->x2f(j);
      const Real& x2f_jp = pmb->x2f(j+1);
      const Real& x2v_j  = pmb->x2v(j);
      const Real& dx2_j  = pmb->dx2f(j);
      lw=(x2f_jp-x2v_j)/dx2_j;
      rw=(x2v_j -x2f_j)/dx2_j;
      bcc2 = lw*b2_j + rw*b2_jp1;
      const Real& x3f_k  = pmb->x3f(k);
      const Real& x3f_kp = pmb->x3f(k+1);
      const Real& x3v_k  = pmb->x3v(k);
      const Real& dx3_k  = pmb->dx3f(k);
      lw=(x3f_kp-x3v_k)/dx3_k;
      rw=(x3v_k -x3f_k)/dx3_k;
      bcc3 = lw*b3_k + rw*b3_kp1;
    }
  }}

  return;
}

//--------------------------------------------------------------------------------------
// \!fn Real FluidEqnOfState::SoundSpeed(Real prim[NFLUID])
// \brief returns adiabatic sound speed given vector of primitive variables

Real FluidEqnOfState::SoundSpeed(const Real prim[NFLUID])
{
  return iso_sound_speed_;
}

//--------------------------------------------------------------------------------------
// \!fn Real FluidEqnOfState::FastMagnetosonicSpeed()
// \brief returns fast magnetosonic speed given vector of primitive variables
// Note the formula for (C_f)^2 is positive definite, so this func never returns a NaN

Real FluidEqnOfState::FastMagnetosonicSpeed(const Real prim[(NWAVE)], const Real bx)
{
  Real asq = (iso_sound_speed_*iso_sound_speed_);
  Real vaxsq = bx*bx/prim[IDN];
  Real ct2 = (prim[IBY]*prim[IBY] + prim[IBZ]*prim[IBZ])/prim[IDN];
  Real qsq = vaxsq + ct2 + asq;
  Real tmp = vaxsq + ct2 - asq;
  return sqrt(0.5*(qsq + sqrt(tmp*tmp + 4.0*asq*ct2)));
}
