//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file adiabatic_hydro.cpp
//  \brief implements functions in class EquationOfState for adiabatic hydrodynamics`

// C/C++ headers
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

// EquationOfState constructor

EquationOfState::EquationOfState(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block_ = pmb;
  density_floor_  = pin->GetOrAddReal("hydro","dfloor",(1024*(FLT_MIN)));
#if EOS_TABLE_ENABLED
  if (pin->DoesParameterExist("hydro","efloor")){
    energy_floor_ = pin->GetReal("hydro","efloor");
  }
  else{
    energy_floor_ = pin->GetOrAddReal("hydro","pfloor",(1024*(FLT_MIN)));
    energy_floor_ /= pin->GetOrAddReal("hydro","gamma", 2.) - 1.;
  }
  GetEosFn = NULL;
  PrepEOS(pin);
  pressure_floor_ = sqrt(-1);
  gamma_ = sqrt(-1);
#else
  pressure_floor_ = pin->GetOrAddReal("hydro","pfloor",(1024*(FLT_MIN)));
  gamma_ = pin->GetReal("hydro","gamma");
#endif
}

// destructor

EquationOfState::~EquationOfState()
{
  #if EOS_TABLE_ENABLED
  CleanEOS();
  #endif
}

//----------------------------------------------------------------------------------------
// \!fn void EquationOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
//           const AthenaArray<Real> &prim_old, const FaceField &b,
//           AthenaArray<Real> &prim, AthenaArray<Real> &bcc, Coordinates *pco,
//           int is, int ie, int js, int je, int ks, int ke)
// \brief Converts conserved into primitive variables in adiabatic hydro.

void EquationOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
  const AthenaArray<Real> &prim_old, const FaceField &b, AthenaArray<Real> &prim,
  AthenaArray<Real> &bcc, Coordinates *pco, int is,int ie, int js,int je, int ks,int ke)
{
#if !EOS_TABLE_ENABLED
  Real gm1 = GetGamma() - 1.0;
#endif

  for (int k=ks; k<=ke; ++k){
  for (int j=js; j<=je; ++j){
#pragma omp simd
    for (int i=is; i<=ie; ++i){
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

      // apply pressure/energy floor, correct total energy
#if EOS_TABLE_ENABLED
      u_e = (u_e - ke > energy_floor_) ?  u_e : energy_floor_ + ke;
      w_p = GetEosData(u_d, u_e - ke, axisEgas, iPresEOS) * (u_e - ke);
#else
      w_p = gm1*(u_e - ke);
      u_e = (w_p > pressure_floor_) ?  u_e : ((pressure_floor_/gm1) + ke);
      w_p = (w_p > pressure_floor_) ?  w_p : pressure_floor_;
#endif
    }
  }}

  return;
}


//----------------------------------------------------------------------------------------
// \!fn void EquationOfState::PrimitiveToConserved(const AthenaArray<Real> &prim,
//           const AthenaArray<Real> &bc, AthenaArray<Real> &cons, Coordinates *pco,
//           int is, int ie, int js, int je, int ks, int ke);
// \brief Converts primitive variables into conservative variables

void EquationOfState::PrimitiveToConserved(const AthenaArray<Real> &prim,
     const AthenaArray<Real> &bc, AthenaArray<Real> &cons, Coordinates *pco,
     int is, int ie, int js, int je, int ks, int ke)
{
#if !EOS_TABLE_ENABLED
  Real igm1 = 1.0/(GetGamma() - 1.0);
#endif


  #pragma omp simd
  for (int k=ks; k<=ke; ++k){
  for (int j=js; j<=je; ++j){
    //#pragma omp simd
    #pragma novector
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

      u_d = w_d;
      u_m1 = w_vx*w_d;
      u_m2 = w_vy*w_d;
      u_m3 = w_vz*w_d;
#if EOS_TABLE_ENABLED
      u_e = GetEosData(u_d, w_p, axisPres, iPresEOS) * w_p + 0.5*w_d*(SQR(w_vx) + SQR(w_vy) + SQR(w_vz));
#else
      u_e = w_p*igm1 + 0.5*w_d*(SQR(w_vx) + SQR(w_vy) + SQR(w_vz));
#endif
    }
  }}

  return;
}

//----------------------------------------------------------------------------------------
// \!fn Real EquationOfState::SoundSpeed(Real prim[NHYDRO])
// \brief returns adiabatic sound speed given vector of primitive variables

Real EquationOfState::SoundSpeed(const Real prim[NHYDRO])
{
#if EOS_TABLE_ENABLED
  return sqrt(GetASqFromRhoPres(prim[IDN], prim[IEN]));
#else
  return sqrt(gamma_*prim[IEN]/prim[IDN]);
#endif
}
