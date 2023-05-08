#ifndef EOS_EOS_HPP_
#define EOS_EOS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file eos.hpp
//! \brief defines class EquationOfState
//!
//!  Contains data and functions that implement the equation of state

// C headers

// C++ headers
#include <limits>     // std::numeric_limits<float>

// Athena++ headers
#include "../athena.hpp"         // Real
#include "../athena_arrays.hpp"  // AthenaArray
#include "../coordinates/coordinates.hpp" // Coordinates
#include "../utils/interp_table.hpp"

// Declarations
class Hydro;
class ParameterInput;
struct FaceField;

//! \class EquationOfState
//! \brief data and functions that implement EoS

class EquationOfState {
  friend class Hydro;
 public:
  EquationOfState(MeshBlock *pmb, ParameterInput *pin);

  void ConservedToPrimitive(
      AthenaArray<Real> &cons, const AthenaArray<Real> &prim_old, const FaceField &b,
      AthenaArray<Real> &prim, AthenaArray<Real> &bcc,
      Coordinates *pco, int il, int iu, int jl, int ju, int kl, int ku);
  void PrimitiveToConserved(const AthenaArray<Real> &prim, const AthenaArray<Real> &bc,
                            AthenaArray<Real> &cons, Coordinates *pco,
                            int il, int iu, int jl, int ju, int kl, int ku);
  void ConservedToPrimitiveCellAverage(
      AthenaArray<Real> &cons, const AthenaArray<Real> &prim_old, const FaceField &b,
      AthenaArray<Real> &prim, AthenaArray<Real> &bcc,
      Coordinates *pco, int il, int iu, int jl, int ju, int kl, int ku);

  // void PrimitiveToConservedCellAverage(const AthenaArray<Real> &prim,
  //   const AthenaArray<Real> &bc, AthenaArray<Real> &cons, Coordinates *pco, int il,
  //   int iu, int jl, int ju, int kl, int ku);

  void PassiveScalarConservedToPrimitive(
      AthenaArray<Real> &s, const AthenaArray<Real> &u, const AthenaArray<Real> &r_old,
      AthenaArray<Real> &r,
      Coordinates *pco, int il, int iu, int jl, int ju, int kl, int ku);
  void PassiveScalarPrimitiveToConserved(
    const AthenaArray<Real> &r, const AthenaArray<Real> &u,
    AthenaArray<Real> &s, Coordinates *pco,
    int il, int iu, int jl, int ju, int kl, int ku);
  void PassiveScalarConservedToPrimitiveCellAverage(
    AthenaArray<Real> &s, const AthenaArray<Real> &r_old, AthenaArray<Real> &r,
    Coordinates *pco, int il, int iu, int jl, int ju, int kl, int ku);

  // pass k, j, i to following 2x functions even though x1-sliced input array is expected
  // in order to accomodate position-dependent floors
#pragma omp declare simd simdlen(SIMD_WIDTH) uniform(this,prim,k,j) linear(i)
  void ApplyPrimitiveFloors(AthenaArray<Real> &prim, int k, int j, int i);

#pragma omp declare simd simdlen(SIMD_WIDTH) uniform(this,s,n,k,j) linear(i)
  void ApplyPassiveScalarFloors(AthenaArray<Real> &s, int n, int k, int j, int i);

#pragma omp declare simd simdlen(SIMD_WIDTH) uniform(this,s,w,r,n,k,j) linear(i)
  void ApplyPassiveScalarPrimitiveConservedFloors(
    AthenaArray<Real> &s, const AthenaArray<Real> &w, AthenaArray<Real> &r,
    int n, int k, int j, int i);

  // Sound speed functions in different regimes
#if !RELATIVISTIC_DYNAMICS  // Newtonian: SR, GR defined as no-op
#pragma omp declare simd simdlen(SIMD_WIDTH) uniform(this)
  Real SoundSpeed(const Real prim[(NHYDRO)]);
  // Define flooring function for fourth-order EOS as no-op for SR, GR regimes
#pragma omp declare simd simdlen(SIMD_WIDTH) uniform(this,prim,cons,bcc,k,j) linear(i)
  void ApplyPrimitiveConservedFloors(
      AthenaArray<Real> &prim, AthenaArray<Real> &cons, AthenaArray<Real> &bcc,
      int k, int j, int i);
#if !MAGNETIC_FIELDS_ENABLED  // Newtonian hydro: Newtonian MHD defined as no-op
  Real FastMagnetosonicSpeed(const Real[], const Real) {return 0.0;}
#else  // Newtonian MHD
#pragma omp declare simd simdlen(SIMD_WIDTH) uniform(this)
  Real FastMagnetosonicSpeed(const Real prim[(NWAVE)], const Real bx);
#endif  // !MAGNETIC_FIELDS_ENABLED
  void SoundSpeedsSR(Real, Real, Real, Real, Real *, Real *) {return;}
  void SoundSpeedsGR(Real, Real, Real, Real, Real, Real, Real, Real *, Real *)
  {return;}
  void FastMagnetosonicSpeedsGR(Real, Real, Real, Real, Real, Real, Real, Real, Real *,
                                Real *) {return;}
#elif !GENERAL_RELATIVITY  // SR: Newtonian, GR defined as no-op
  Real SoundSpeed(const Real[]) {return 0.0;}
  Real FastMagnetosonicSpeed(const Real[], const Real) {return 0.0;}
  void ApplyPrimitiveConservedFloors(
      AthenaArray<Real> &, AthenaArray<Real> &, AthenaArray<Real> &,
      int, int, int) {return;}
#if !MAGNETIC_FIELDS_ENABLED  // SR hydro: SR MHD defined as no-op
#pragma omp declare simd simdlen(SIMD_WIDTH) uniform(this)
  void SoundSpeedsSR(Real rho_h, Real pgas, Real vx, Real gamma_lorentz_sq,
                     Real *plambda_plus, Real *plambda_minus);
  void FastMagnetosonicSpeedsGR(Real, Real, Real, Real, Real, Real, Real, Real,
                                Real *, Real *) {return;}
#else  // SR MHD: SR hydro defined as no-op
  void SoundSpeedsSR(Real, Real, Real, Real, Real *, Real *) {return;}
  void FastMagnetosonicSpeedsGR(Real wgas, Real pgas, Real u0, Real u1, Real b_sq,
                                Real g00, Real g01, Real g11,
                                Real *p_lambda_plus, Real *p_lambda_minus);
#endif  // !MAGNETIC_FIELDS_ENABLED
  void SoundSpeedsGR(Real, Real, Real, Real, Real, Real, Real, Real *, Real *)
  {return;}
#else  // GR: Newtonian defined as no-op
  Real SoundSpeed(const Real[]) {return 0.0;}
  Real FastMagnetosonicSpeed(const Real[], const Real) {return 0.0;}
  void ApplyPrimitiveConservedFloors(
      AthenaArray<Real> &, AthenaArray<Real> &, AthenaArray<Real> &,
      int, int, int) {return;}
#if !MAGNETIC_FIELDS_ENABLED  // GR hydro: GR+SR MHD defined as no-op
#pragma omp declare simd simdlen(SIMD_WIDTH) uniform(this)
  void SoundSpeedsSR(Real rho_h, Real pgas, Real vx, Real gamma_lorentz_sq,
                     Real *plambda_plus, Real *plambda_minus);
#pragma omp declare simd simdlen(SIMD_WIDTH) uniform(this)
  void SoundSpeedsGR(Real rho_h, Real pgas, Real u0, Real u1,
                     Real g00, Real g01, Real g11,
                     Real *plambda_plus, Real *plambda_minus);
  void FastMagnetosonicSpeedsGR(Real, Real, Real, Real, Real, Real, Real, Real,
                                Real *, Real *) {return;}
#else  // GR MHD: GR+SR hydro defined as no-op
  void SoundSpeedsSR(Real, Real, Real, Real, Real *, Real *) {return;}
  void SoundSpeedsGR(Real, Real, Real, Real, Real, Real, Real, Real *, Real *)
  {return;}
#pragma omp declare simd simdlen(SIMD_WIDTH) uniform(this)
  void FastMagnetosonicSpeedsGR(Real wgas, Real pgas, Real u0, Real u1, Real b_sq,
                                Real g00, Real g01, Real g11,
                                Real *p_lambda_plus, Real *p_lambda_minus);
#endif  // !MAGNETIC_FIELDS_ENABLED (GR)
#endif  // #else (#if !RELATIVISTIC_DYNAMICS, #elif !GENERAL_RELATIVITY)

  Real PresFromRhoEg(Real rho, Real egas);
  Real EgasFromRhoP(Real rho, Real pres);
  Real AsqFromRhoP(Real rho, Real pres);
  Real GetIsoSoundSpeed() const {return iso_sound_speed_;}
  Real GetDensityFloor() const {return density_floor_;}
  Real GetPressureFloor() const {return pressure_floor_;}
  EosTable* ptable; // pointer to EOS table data
#if GENERAL_EOS
  Real GetGamma();
#else // not GENERAL_EOS
  Real GetGamma() const {return gamma_;}
#endif

 private:
  // (C++11) in-class Default Member Initializer (fallback option):
  const Real float_min{std::numeric_limits<float>::min()};
  MeshBlock *pmy_block_;                 // ptr to MeshBlock containing this EOS
  Real iso_sound_speed_, gamma_;         // isothermal Cs, ratio of specific heats
  Real density_floor_, pressure_floor_;  // density and pressure floors
  Real energy_floor_;                    // energy floor
  Real scalar_floor_; // dimensionless concentration floor
  Real sigma_max_, beta_min_;            // limits on ratios of gas quantities to pmag
  Real gamma_max_;                       // maximum Lorentz factor
  Real rho_min_, rho_pow_;               // variables to control power-law denity floor
  Real pgas_min_, pgas_pow_;             // variables to control power-law pressure floor
  Real rho_unit_, inv_rho_unit_;         // physical unit/sim unit for mass density
  Real egas_unit_, inv_egas_unit_;       // physical unit/sim unit for energy density
  Real vsqr_unit_, inv_vsqr_unit_;       // physical unit/sim unit for speed^2
  AthenaArray<Real> g_, g_inv_;          // metric and its inverse, used in GR
  AthenaArray<Real> normal_dd_;          // normal-frame densities, used in relativity
  AthenaArray<Real> normal_ee_;          // normal-frame energies, used in relativity
  AthenaArray<Real> normal_mm_;          // normal-frame momenta, used in relativity
  AthenaArray<Real> normal_bb_;          // normal-frame fields, used in relativistic MHD
  AthenaArray<Real> normal_tt_;          // normal-frame M.B, used in relativistic MHD
  void InitEosConstants(ParameterInput *pin);
};

#endif // EOS_EOS_HPP_
