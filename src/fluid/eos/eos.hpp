#ifndef EOS_HPP
#define EOS_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file eos.hpp
//  \brief defines class FluidEqnOfState
//  Contains data and functions that implement the equation of state for the fluid
//======================================================================================

// Athena headers
#include "../../athena.hpp"         // Real
#include "../../athena_arrays.hpp"  // AthenaArray

// Declarations
class Fluid;
class ParameterInput;
struct InterfaceField;

//! \class FluidEqnOfState
//  \brief data and functions that implement EoS for fluid

class FluidEqnOfState {
public:
  FluidEqnOfState(Fluid *pf, ParameterInput *pin);
  ~FluidEqnOfState();

  void ConservedToPrimitive(AthenaArray<Real> &cons, const AthenaArray<Real> &prim_old,
    const InterfaceField &b, AthenaArray<Real> &prim, AthenaArray<Real> &bcc);

  Real SoundSpeed(const Real prim[NFLUID]); 
  Real FastMagnetosonicSpeed(const Real prim[((NFLUID)+(NFIELD)-1)], const Real bx); 
  void FastMagnetosonicSpeedsRelativistic(Real rho, Real pgas,
      Real vx, Real vy, Real vz,
      Real ut, Real ux, Real uy, Real uz,
      Real bx, Real by, Real bz,
      Real bcovt, Real bcovx, Real bcovy, Real bcovz,
      Real &lambda_plus, Real &lambda_minus);
  Real GetGamma() const {return gamma_;}

private:
  Fluid *pmy_fluid_;             // ptr to Fluid containing this EqnOfState
  Real iso_sound_speed_, gamma_; // isothermal Cs, ratio of specific heats
  AthenaArray<Real> g_, g_inv_;  // metric and its inverse, used for cons->prim in GR
  Real density_floor_, pressure_floor_; // density and pressure floors
};
#endif
