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

  Real SoundSpeed(const Real prim[(NFLUID)]); 
  Real FastMagnetosonicSpeed(const Real prim[(NWAVE)], const Real bx); 
  void SoundSpeedsSR(Real rho_h, Real pgas, Real vx, Real gamma_lorentz_sq,
      Real *plambda_plus, Real *plambda_minus);
  void FastMagnetosonicSpeedsSR(Real rho, Real pgas, const Real u[4], const Real b[4],
      Real *plambda_plus, Real *plambda_minus);
  void SoundSpeedsGR(Real rho_h, Real pgas, Real u0, Real u1,
      Real g00, Real g01, Real g11,
      Real *plambda_plus, Real *plambda_minus);
  void FastMagnetosonicSpeedsGR(Real rho_h, Real pgas, Real u0, Real u1, Real b_sq,
      Real g00, Real g01, Real g11,
      Real *plambda_plus, Real *plambda_minus);
  Real GetGamma() const {return gamma_;}
  Real GetIsoSoundSpeed() const {return iso_sound_speed_;}
  Real GetDensityFloor() const {return density_floor_;}
  Real GetPressureFloor() const {return pressure_floor_;}

private:
  Fluid *pmy_fluid_;             // ptr to Fluid containing this EqnOfState
  Real iso_sound_speed_, gamma_; // isothermal Cs, ratio of specific heats
  AthenaArray<Real> g_, g_inv_;  // metric and its inverse, used for cons->prim in GR
  Real density_floor_, pressure_floor_; // density and pressure floors
};
#endif
