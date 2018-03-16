//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file adiabatic_mhd.cpp
//  \brief implements functions in class EquationOfState for adiabatic MHD

// C++ headers
#include <cmath>   // sqrt()
#include <cfloat>  // FLT_MIN

// Athena++ headers
#include "eos.hpp"
#include "../hydro/hydro.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../field/field.hpp"
#include "../coordinates/coordinates.hpp"

// EquationOfState constructor

EquationOfState::EquationOfState(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block_ = pmb;
  gamma_ = pin->GetReal("hydro","gamma");
  density_floor_  = pin->GetOrAddReal("hydro","dfloor",(1024*(FLT_MIN)));
  pressure_floor_ = pin->GetOrAddReal("hydro","pfloor",(1024*(FLT_MIN)));
}

// destructor

EquationOfState::~EquationOfState()
{
}

//----------------------------------------------------------------------------------------
// \!fn void EquationOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
//    const AthenaArray<Real> &prim_old, const FaceField &b,
//    AthenaArray<Real> &prim, AthenaArray<Real> &bcc, Coordinates *pco,
//    int is, int ie, int js, int je, int ks, int ke);
// \brief For the Hydro, converts conserved into primitive variables in adiabatic MHD.
//  For the Field, computes cell-centered from face-centered magnetic field.

void EquationOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
    const AthenaArray<Real> &prim_old, const FaceField &b, AthenaArray<Real> &prim,
    AthenaArray<Real> &bcc, Coordinates *pco,
    int is, int ie, int js, int je, int ks, int ke)
{
  Real gm1 = GetGamma() - 1.0;

  pmy_block_->pfield->CalculateCellCenteredField(b,bcc,pco,is,ie,js,je,ks,ke);

  for (int k=ks; k<=ke; ++k){
  for (int j=js; j<=je; ++j){
#pragma omp simd
    for (int i=is; i<=ie; ++i){
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

  return;
}

//----------------------------------------------------------------------------------------
// \!fn void EquationOfState::PrimitiveToConserved(const AthenaArray<Real> &prim,
//           const AthenaArray<Real> &bc, AthenaArray<Real> &cons, Coordinates *pco,
//           int is, int ie, int js, int je, int ks, int ke);
// \brief Converts primitive variables into conservative variables
//        Note that this function assumes cell-centered fields are already calculated

void EquationOfState::PrimitiveToConserved(const AthenaArray<Real> &prim,
     const AthenaArray<Real> &bc, AthenaArray<Real> &cons, Coordinates *pco,
     int is, int ie, int js, int je, int ks, int ke)
{
  Real igm1 = 1.0/(GetGamma() - 1.0);

  for (int k=ks; k<=ke; ++k){
  for (int j=js; j<=je; ++j){
#pragma omp simd
    for (int i=is; i<=ie; ++i){
      Real& u_d  = cons(IDN,k,j,i);
      Real& u_m1 = cons(IM1,k,j,i);
      Real& u_m2 = cons(IM2,k,j,i);
      Real& u_m3 = cons(IM3,k,j,i);
      Real& u_e  = cons(IEN,k,j,i);

      const Real& w_d  = prim(IDN,k,j,i);
      const Real& w_vx = prim(IVX,k,j,i);
      const Real& w_vy = prim(IVY,k,j,i);
      const Real& w_vz = prim(IVZ,k,j,i);
      const Real& w_p  = prim(IEN,k,j,i);

      const Real& bcc1 = bc(IB1,k,j,i);
      const Real& bcc2 = bc(IB2,k,j,i);
      const Real& bcc3 = bc(IB3,k,j,i);

      u_d = w_d;
      u_m1 = w_vx*w_d;
      u_m2 = w_vy*w_d;
      u_m3 = w_vz*w_d;
      u_e = w_p*igm1 + 0.5*(w_d*(SQR(w_vx) + SQR(w_vy) + SQR(w_vz))
            + (SQR(bcc1) + SQR(bcc2) + SQR(bcc3)));
    }
  }}

  return;
}

//----------------------------------------------------------------------------------------
// \!fn Real EquationOfState::SoundSpeed(Real prim[NHYDRO])
// \brief returns adiabatic sound speed given vector of primitive variables

Real EquationOfState::SoundSpeed(const Real prim[NHYDRO])
{
  return sqrt(GetGamma()*prim[IEN]/prim[IDN]);
}

//----------------------------------------------------------------------------------------
// \!fn Real EquationOfState::FastMagnetosonicSpeed(const Real prim[], const Real bx)
// \brief returns fast magnetosonic speed given vector of primitive variables
// Note the formula for (C_f)^2 is positive definite, so this func never returns a NaN 

Real EquationOfState::FastMagnetosonicSpeed(const Real prim[(NWAVE)], const Real bx)
{
  Real asq = GetGamma()*prim[IEN]/prim[IDN];
  Real vaxsq = bx*bx/prim[IDN];
  Real ct2 = (prim[IBY]*prim[IBY] + prim[IBZ]*prim[IBZ])/prim[IDN];
  Real qsq = vaxsq + ct2 + asq;
  Real tmp = vaxsq + ct2 - asq;
  return sqrt(0.5*(qsq + sqrt(tmp*tmp + 4.0*asq*ct2)));
}
