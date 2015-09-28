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
//! \file adiabatic_mhd.cpp
//  \brief implements functions in class HydroEqnOfState for adiabatic MHD
//======================================================================================

// C++ headers
#include <cmath>   // sqrt()
#include <cfloat>  // FLT_MIN

// Athena++ headers
#include "../hydro.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh.hpp"
#include "../../parameter_input.hpp"
#include "../../field/field.hpp"
#include "../../coordinates/coordinates.hpp"

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
// \brief For the Hydro, converts conserved into primitive variables in adiabatic MHD.
//  For the Field, computes cell-centered from face-centered magnetic field.

void HydroEqnOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
  const AthenaArray<Real> &prim_old, const InterfaceField &b, AthenaArray<Real> &prim,
  AthenaArray<Real> &bcc)
{
  MeshBlock *pmb = pmy_fluid_->pmy_block;
  Coordinates *pco = pmb->pcoord;
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

  // calc cell centered fields first
#pragma simd
    for (int i=pmb->is-(NGHOST); i<=pmb->ie+(NGHOST); ++i){
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
      const Real& x1f_i  = pco->x1f(i);
      const Real& x1f_ip = pco->x1f(i+1);
      const Real& x1v_i  = pco->x1v(i);
      const Real& dx1_i  = pco->dx1f(i);
      Real lw=(x1f_ip-x1v_i)/dx1_i;
      Real rw=(x1v_i -x1f_i)/dx1_i;
      bcc1 = lw*b1_i + rw*b1_ip1;
      const Real& x2f_j  = pco->x2f(j);
      const Real& x2f_jp = pco->x2f(j+1);
      const Real& x2v_j  = pco->x2v(j);
      const Real& dx2_j  = pco->dx2f(j);
      lw=(x2f_jp-x2v_j)/dx2_j;
      rw=(x2v_j -x2f_j)/dx2_j;
      bcc2 = lw*b2_j + rw*b2_jp1;
      const Real& x3f_k  = pco->x3f(k);
      const Real& x3f_kp = pco->x3f(k+1);
      const Real& x3v_k  = pco->x3v(k);
      const Real& dx3_k  = pco->dx3f(k);
      lw=(x3f_kp-x3v_k)/dx3_k;
      rw=(x3v_k -x3f_k)/dx3_k;
      bcc3 = lw*b3_k + rw*b3_kp1;
    }

#pragma simd
    for (int i=pmb->is-(NGHOST); i<=pmb->ie+(NGHOST); ++i){
      Real& u_d  = cons(IDN,k,j,i);
      Real& u_m1 = cons(IVX,k,j,i);
      Real& u_m2 = cons(IVY,k,j,i);
      Real& u_m3 = cons(IVZ,k,j,i);
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

      const Real& bcc1 = bcc(IB1,k,j,i);
      const Real& bcc2 = bcc(IB2,k,j,i);
      const Real& bcc3 = bcc(IB3,k,j,i);

      Real pb = 0.5*(SQR(bcc1) + SQR(bcc2) + SQR(bcc3));
      Real ke = 0.5*di*(SQR(u_m1) + SQR(u_m2) + SQR(u_m3));
      w_p = gm1*(u_e - ke - pb);

      // apply pressure floor, correct total energy
      u_e = (w_p > pressure_floor_) ?  u_e : ((pressure_floor_/gm1) + ke + pb);
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
  return sqrt(GetGamma()*prim[IEN]/prim[IDN]);
}

//--------------------------------------------------------------------------------------
// \!fn Real HydroEqnOfState::FastMagnetosonicSpeed(const Real prim[], const Real bx)
// \brief returns fast magnetosonic speed given vector of primitive variables
// Note the formula for (C_f)^2 is positive definite, so this func never returns a NaN 

Real HydroEqnOfState::FastMagnetosonicSpeed(const Real prim[(NWAVE)], const Real bx)
{
  Real asq = GetGamma()*prim[IEN]/prim[IDN];
  Real vaxsq = bx*bx/prim[IDN];
  Real ct2 = (prim[IBY]*prim[IBY] + prim[IBZ]*prim[IBZ])/prim[IDN];
  Real qsq = vaxsq + ct2 + asq;
  Real tmp = vaxsq + ct2 - asq;
  return sqrt(0.5*(qsq + sqrt(tmp*tmp + 4.0*asq*ct2)));
}
