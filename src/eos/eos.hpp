#ifndef EOS_EOS_HPP_
#define EOS_EOS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file eos.hpp
//  \brief defines class EquationOfState
//  Contains data and functions that implement the equation of state

// Athena headers
#include "../athena.hpp"         // Real
#include "../athena_arrays.hpp"  // AthenaArray
#include "../coordinates/coordinates.hpp" // Coordinates

// Declarations
class Hydro;
class ParameterInput;
struct FaceField;

//! \class EquationOfState
//  \brief data and functions that implement EoS

class EquationOfState {
public:
  EquationOfState(MeshBlock *pmb, ParameterInput *pin);
  ~EquationOfState();

  void ConservedToPrimitive(AthenaArray<Real> &cons, const AthenaArray<Real> &prim_old,
    const FaceField &b, AthenaArray<Real> &prim, AthenaArray<Real> &bcc,
    Coordinates *pco, int il, int iu, int jl, int ju, int kl, int ku);
  void PrimitiveToConserved(const AthenaArray<Real> &prim, const AthenaArray<Real> &bc,
       AthenaArray<Real> &cons, Coordinates *pco,
       int il, int iu, int jl, int ju, int kl, int ku);
#pragma omp declare simd simdlen(SIMD_WIDTH) uniform(this,prim,k,j) linear(i)
  void ApplyPrimitiveFloors(AthenaArray<Real> &prim, int k, int j, int i);

  // Sound speed functions in different regimes
  #if !RELATIVISTIC_DYNAMICS  // Newtonian: SR, GR defined as no-op
#pragma omp declare simd simdlen(SIMD_WIDTH) uniform(this)
    Real SoundSpeed(const Real prim[(NHYDRO)]);
    #if !MAGNETIC_FIELDS_ENABLED  // hydro: MHD defined as no-op
      Real FastMagnetosonicSpeed(const Real[], const Real) {return 0.0;}
    #else  // MHD
#pragma omp declare simd simdlen(SIMD_WIDTH) uniform(this)
      Real FastMagnetosonicSpeed(const Real prim[(NWAVE)], const Real bx);
    #endif  // !MAGNETIC_FIELDS_ENABLED
    void SoundSpeedsSR(Real, Real, Real, Real, Real *, Real *) {return;}
    void FastMagnetosonicSpeedsSR(const AthenaArray<Real> &,
        const AthenaArray<Real> &, int, int, int, int, int, AthenaArray<Real> &,
        AthenaArray<Real> &) {return;}
    void SoundSpeedsGR(Real, Real, Real, Real, Real, Real, Real, Real *, Real *)
        {return;}
    void FastMagnetosonicSpeedsGR(Real, Real, Real, Real, Real, Real, Real, Real, Real *,
        Real *) {return;}
  #elif !GENERAL_RELATIVITY  // SR: Newtonian, GR defined as no-op
    Real SoundSpeed(const Real[]) {return 0.0;}
    Real FastMagnetosonicSpeed(const Real[], const Real) {return 0.0;}
    #if !MAGNETIC_FIELDS_ENABLED  // hydro: MHD defined as no-op
      void SoundSpeedsSR(Real rho_h, Real pgas, Real vx, Real gamma_lorentz_sq,
          Real *plambda_plus, Real *plambda_minus);
      void FastMagnetosonicSpeedsSR(const AthenaArray<Real> &,
          const AthenaArray<Real> &, int, int, int, int, int, AthenaArray<Real> &,
          AthenaArray<Real> &) {return;}
    #else  // MHD: hydro defined as no-op
      void SoundSpeedsSR(Real, Real, Real, Real, Real *, Real *) {return;}
      void FastMagnetosonicSpeedsSR(const AthenaArray<Real> &prim,
          const AthenaArray<Real> &bbx_vals, int k, int j, int il, int iu, int ivx,
          AthenaArray<Real> &lambdas_p, AthenaArray<Real> &lambdas_m);
    #endif  // !MAGNETIC_FIELDS_ENABLED
    void SoundSpeedsGR(Real, Real, Real, Real, Real, Real, Real, Real *, Real *)
        {return;}
    void FastMagnetosonicSpeedsGR(Real, Real, Real, Real, Real, Real, Real, Real,
        Real *, Real *) {return;}
  #else  // GR: Newtonian defined as no-op
    Real SoundSpeed(const Real[]) {return 0.0;}
    Real FastMagnetosonicSpeed(const Real[], const Real) {return 0.0;}
    #if !MAGNETIC_FIELDS_ENABLED  // hydro: MHD defined as no-op
      void SoundSpeedsSR(Real rho_h, Real pgas, Real vx, Real gamma_lorentz_sq,
          Real *plambda_plus, Real *plambda_minus);
      void FastMagnetosonicSpeedsSR(const AthenaArray<Real> &,
          const AthenaArray<Real> &, int, int, int, int, int, AthenaArray<Real> &,
          AthenaArray<Real> &) {return;}
      void SoundSpeedsGR(Real rho_h, Real pgas, Real u0, Real u1,
          Real g00, Real g01, Real g11,
          Real *plambda_plus, Real *plambda_minus);
      void FastMagnetosonicSpeedsGR(Real, Real, Real, Real, Real, Real, Real, Real,
          Real *, Real *) {return;}
    #else  // MHD: hydro defined as no-op
      void SoundSpeedsSR(Real, Real, Real, Real, Real *, Real *) {return;}
      void FastMagnetosonicSpeedsSR(const AthenaArray<Real> &prim,
          const AthenaArray<Real> &bbx_vals, int k, int j, int il, int iu, int ivx,
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
  MeshBlock *pmy_block_;                 // ptr to MeshBlock containing this EOS
  Real iso_sound_speed_, gamma_;         // isothermal Cs, ratio of specific heats
  Real density_floor_, pressure_floor_;  // density and pressure floors
  Real sigma_max_, beta_min_;            // limits on ratios of gas quantities to pmag
  Real gamma_max_;                       // maximum Lorentz factor
  Real rho_min_, rho_pow_;               // variables to control power-law denity floor
  Real pgas_min_, pgas_pow_;             // variables to control power-law pressure floor
  AthenaArray<Real> g_, g_inv_;          // metric and its inverse, used in GR
  AthenaArray<Real> fixed_;              // cells with problems, used in GR hydro
  AthenaArray<Real> normal_dd_;          // normal-frame densities, used in GR MHD
  AthenaArray<Real> normal_ee_;          // normal-frame energies, used in GR MHD
  AthenaArray<Real> normal_mm_;          // normal-frame momenta, used in GR MHD
  AthenaArray<Real> normal_bb_;          // normal-frame fields, used in GR MHD
  AthenaArray<Real> normal_tt_;          // normal-frame M.B, used in GR MHD
};

#endif // EOS_EOS_HPP_
