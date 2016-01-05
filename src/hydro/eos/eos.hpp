#ifndef EOS_HPP
#define EOS_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file eos.hpp
//  \brief defines class HydroEqnOfState
//  Contains data and functions that implement the equation of state for hydro
//======================================================================================

// Athena headers
#include "../../athena.hpp"         // Real
#include "../../athena_arrays.hpp"  // AthenaArray
#include "../../coordinates/coordinates.hpp" // Coordinates

// Declarations
class Hydro;
class ParameterInput;
struct FaceField;

//! \class HydroEqnOfState
//  \brief data and functions that implement EoS for hydro

class HydroEqnOfState {
public:
  HydroEqnOfState(Hydro *pf, ParameterInput *pin);
  ~HydroEqnOfState();

  void ConservedToPrimitive(AthenaArray<Real> &cons, const AthenaArray<Real> &prim_old,
    const FaceField &b, AthenaArray<Real> &prim, AthenaArray<Real> &bcc,
    Coordinates *pco, int is, int ie, int js, int je, int ks, int ke);
  void PrimitiveToConserved(const AthenaArray<Real> &prim, const AthenaArray<Real> &bc,
       AthenaArray<Real> &cons, Coordinates *pco,
       int is, int ie, int js, int je, int ks, int ke);

  // Sound speed functions in different regimes
  #if !RELATIVISTIC_DYNAMICS  // Newtonian: SR, GR defined as no-op
    Real SoundSpeed(const Real prim[(NHYDRO)]);
    #if !MAGNETIC_FIELDS_ENABLED  // hydro: MHD defined as no-op
      Real FastMagnetosonicSpeed(const Real [], const Real) {return 0.0;}
    #else  // MHD
      Real FastMagnetosonicSpeed(const Real prim[(NWAVE)], const Real bx);
    #endif  // !MAGNETIC_FIELDS_ENABLED
    void SoundSpeedsSR(Real, Real, Real, Real, Real *, Real *) {return;}
    void FastMagnetosonicSpeedsSR(const AthenaArray<Real> &,
        const AthenaArray<Real> &, int, int, int, AthenaArray<Real> &,
        AthenaArray<Real> &) {return;}
    void SoundSpeedsGR(Real, Real, Real, Real, Real, Real, Real, Real *, Real *)
        {return;}
    void FastMagnetosonicSpeedsGR(Real, Real, Real, Real, Real, Real, Real, Real, Real *,
        Real *) {return;}
  #elif !GENERAL_RELATIVITY  // SR: Newtonian, GR defined as no-op
    Real SoundSpeed(const Real []) {return 0.0;}
    Real FastMagnetosonicSpeed(const Real [], const Real) {return 0.0;}
    #if !MAGNETIC_FIELDS_ENABLED  // hydro: MHD defined as no-op
      void SoundSpeedsSR(Real rho_h, Real pgas, Real vx, Real gamma_lorentz_sq,
          Real *plambda_plus, Real *plambda_minus);
      void FastMagnetosonicSpeedsSR(const AthenaArray<Real> &,
          const AthenaArray<Real> &, int, int, int, AthenaArray<Real> &,
          AthenaArray<Real> &) {return;}
    #else  // MHD: hydro defined as no-op
      void SoundSpeedsSR(Real, Real, Real, Real, Real *, Real *) {return;}
      void FastMagnetosonicSpeedsSR(const AthenaArray<Real> &prim,
          const AthenaArray<Real> &bbx_vals, int il, int iu, int ivx,
          AthenaArray<Real> &lambdas_p, AthenaArray<Real> &lambdas_m);
    #endif  // !MAGNETIC_FIELDS_ENABLED
    void SoundSpeedsGR(Real, Real, Real, Real, Real, Real, Real, Real *, Real *)
        {return;}
    void FastMagnetosonicSpeedsGR(Real, Real, Real, Real, Real, Real, Real, Real,
        Real *, Real *) {return;}
  #else  // GR: Newtonian defined as no-op
    Real SoundSpeed(const Real []) {return 0.0;}
    Real FastMagnetosonicSpeed(const Real [], const Real) {return 0.0;}
    #if !MAGNETIC_FIELDS_ENABLED  // hydro: MHD defined as no-op
      void SoundSpeedsSR(Real rho_h, Real pgas, Real vx, Real gamma_lorentz_sq,
          Real *plambda_plus, Real *plambda_minus);
      void FastMagnetosonicSpeedsSR(const AthenaArray<Real> &,
          const AthenaArray<Real> &, int, int, int, AthenaArray<Real> &,
          AthenaArray<Real> &) {return;}
      void SoundSpeedsGR(Real rho_h, Real pgas, Real u0, Real u1,
          Real g00, Real g01, Real g11,
          Real *plambda_plus, Real *plambda_minus);
      void FastMagnetosonicSpeedsGR(Real, Real, Real, Real, Real, Real, Real, Real,
          Real *, Real *) {return;}
    #else  // MHD: hydro defined as no-op
      void SoundSpeedsSR(Real, Real, Real, Real, Real *, Real *) {return;}
      void FastMagnetosonicSpeedsSR(const AthenaArray<Real> &prim,
          const AthenaArray<Real> &bbx_vals, int il, int iu, int ivx,
          AthenaArray<Real> &lambdas_p, AthenaArray<Real> &lambdas_m);
      void SoundSpeedsGR(Real, Real, Real, Real, Real, Real, Real, Real *, Real *)
          {return;}
      void FastMagnetosonicSpeedsGR(Real rho_h, Real pgas, Real u0, Real u1, Real b_sq,
          Real g00, Real g01, Real g11,
          Real *plambda_plus, Real *plambda_minus);
    #endif  // !MAGNETIC_FIELDS_ENABLED
  #endif  // !RELATIVISTIC_DYNAMICS

  Real GetGamma() const {return gamma_;}
  Real GetIsoSoundSpeed() const {return iso_sound_speed_;}
  Real GetDensityFloor() const {return density_floor_;}
  Real GetPressureFloor() const {return pressure_floor_;}

private:
  Hydro *pmy_hydro_;                     // ptr to Hydro containing this EqnOfState
  Real iso_sound_speed_, gamma_;         // isothermal Cs, ratio of specific heats
  Real density_floor_, pressure_floor_;  // density and pressure floors
  Real rho_pmag_min_, u_pmag_min_;       // floors on ratios of gas quantities to pmag
  Real gamma_max_;                       // maximum Lorentz factor
  Real rho_min_, rho_pow_;               // variables to control power-law denity floor
  Real u_min_, u_pow_;                   // variables to control power-law energy floor
  AthenaArray<Real> g_, g_inv_;          // metric and its inverse, used in GR
  AthenaArray<bool> fixed_;              // array for flagging fixed cells, used in GR
};

#endif
