#ifndef NR_RADIATION_RADIATION_HPP_
#define NR_RADIATION_RADIATION_HPP_
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file radiation.hpp
//  \brief definitions for Radiation class
//======================================================================================

// C headers

// C++ headers
#include <cstdint>     // int64_t
#include <functional>  // reference_wrapper
#include <string>
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/cc/nr_radiation/bvals_rad.hpp"
#include "../parameter_input.hpp"

class MeshBlock;
class ParameterInput;
class RadIntegrator;

//! \class Radiation
//  \brief radiation data and functions

// prototype for user-defined opacity function for radiative transfer

// Array indices for radiation moments
enum {IER=0, IFR1=1, IFR2=2, IFR3=3, IPR11=4, IPR22=5, IPR33=6, IPR12=7,
      IPR13=8, IPR23=9, IPR21=10, IPR31=11, IPR32=12};

enum {OPAS=0, OPAA=1, OPAP=2}; // scattering, absorption, Planck, opacity

class NRRadiation {
  friend class RadIntegrator;
  friend class Mesh;
 public:
  NRRadiation(MeshBlock *pmb, ParameterInput *pin);
  ~NRRadiation();

  MeshBlock* pmy_block;    // ptr to MeshBlock containing this Fluid
  AthenaArray<Real> ir, ir1, ir2, ir_old; // radiation specific intensity
  AthenaArray<Real> ir_gray;
  AthenaArray<Real> rad_mom; // frequency integrated radiation moments
  AthenaArray<Real> rad_mom_cm; // co-moving frame Er, Frx, Fry, Frz
  AthenaArray<Real> rad_mom_nu, rad_mom_cm_nu; // multi_group radiation moments
  AthenaArray<Real> sigma_s, sigma_a; //   opacity
                  //scattering and fluxes weighted Rosseland mean
  AthenaArray<Real> sigma_p, sigma_pe;
                  // Planck mean and radiation energy weighted mean
  AthenaArray<Real> output_sigma; // frequency integrated opacity
  AthenaArray<Real> mu, wmu; // angles and weight

  AthenaArray<Real> flux[3]; // store transport flux, also need for refinement

  AthenaArray<Real> coarse_ir_;
  int refinement_idx{-1};

  Real prat, crat; // prat=aT^4/P_0, crat=c/c_s
  Real vmax;
  Real reduced_c; // reduced speed of light
  Real tunit, telectron; // gas temperature cgs unit,
                         // effective electron scattering temperature
  Real rhounit; // density unit
  Real mol_weight; // mean molecular weight
  Real lunit;  // length unit

  Real sum_diff;
  Real sum_full; // store

  int nang, noct, n_fre_ang; // n_fre_ang=nang*nfreq
  int angle_flag;
  int polar_angle;
  // variables related to the angular space transport
  int nzeta, npsi;
  AthenaArray<Real> coszeta_v, zeta_v_full, zeta_f_full, dzeta_v, dzeta_f,
                    coszeta_f, len_zeta;
  AthenaArray<Real> psi_v, psi_f, len_psi, psi_v_full, psi_f_full,
                    dpsi_v, dpsi_f, sin_psi_f, cot_theta;

  // The frequency grid
  Real nu_min, nu_max;
  // mininum and maximum frequencies, and number of frequency bins
  int nfreq; // number of frequency bins
  int restart_from_gray; //
  Real fre_ratio; // ratio between neighboring frequency bins
  // frequency grid, center of each frequency bin
  AthenaArray<Real> nu_grid, nu_cen, delta_nu;
  // gas emission term in each frequency bin relative to a_rT^4
  AthenaArray<Real> emission_spec;
  FrequencyFunc UserFrequency; // user defined frequency grid
  void EnrollFrequencyFunction(FrequencyFunc MyFrequencyFunction);
  EmissionFunc UserEmissionSpec;
  void EnrollEmissionFunction(EmissionFunc MyEmissionSpec);

  //  int ir_output; // the number of specific intensity to dump
  //  AthenaArray<int> ir_index; // the array
  //  AthenaArray<Real> dump_ir;

  RadBoundaryVariable rad_bvar;

  RadIntegrator *pradintegrator;

  // The function pointer for the opacity
  OpacityFunc UpdateOpacity;

  Real kappa_es; // the frequency independent electron scattering opacity

  int rotate_theta; // flag to rotate the boundary
  int rotate_phi;
  int set_source_flag; // flag to add radiation source term or not

  // Functions
  // Function in problem generators to update opacity
  void EnrollOpacityFunction(OpacityFunc MyOpacityFunction);

  void CalculateMoment(AthenaArray<Real> &ir_in);
  void CalculateComMoment();

  void AngularGrid(int angle_flag, int nmu);
  void AngularGrid(int angle_flag, int nzeta, int npsi);

  void FrequencyGrid();

  Real FitIntPlanckFunc(Real nu_t);

  Real IntPlanckFunc(Real nu_min, Real nu_max);
  Real EffectiveBlackBody(Real intensity, Real nu);
  Real EffectiveBlackBodyNNu2(Real n_nu2, Real nu);
  Real IntegrateBBNuJ(Real nu_t); // integral of
  Real IntegrateBBJONuSq(Real nu_t); //\integral of (j/\nu)^2d\nu
  Real IntegrateBBNNu2(Real nu_t); // ingral of n\nu^2d\nu
  Real ConvertBBJNNu2(Real &bb_j, Real &nu_f);
  // Convert from n\nu^2 dnu to J
  Real InverseConvertBBJNNu2(Real &nnu2, Real &nu_f);
  // Convert J to \int (J/nu)^2
  Real BBJToJONuSq(Real &bb_j, Real &nu_f);
  // Convert J to n(nu_f)
  Real BBJtoNnu(Real &bb_j, Real &nu_f);
  // Convert J to \int J\nu
  Real BBJtoJnu(Real &bb_j, Real &nu_f);
  Real DBBjDNNu2(Real &bb_j, Real &nu_f);

  // convert j to different quantities assuming Wien spectrum
  void ConvertBBJWien(Real &bb_j, Real &nu_f, Real &tgas,
                                 Real &nuj, Real &jonusq);
  void ConvertBBJWien2(Real &bb_j, Real &nu_f, Real &tgas,
                                 Real &nnu2, Real &n_nuf);
  Real InverseConvertBBJNNu2Wien(Real &nnu2, Real &nu_f, Real &tgas);


  AthenaArray<Real> t_floor_, t_ceiling_; // temperature floor

 private:
  int user_unit_;
  // temporary arrays for co-moving moments
  AthenaArray<Real> cosx_cm_, cosy_cm_, cosz_cm_;

  friend class BoundaryValues;

  // to use log frequency spaceing or not
  int log_fre_;
};

#endif // NR_RADIATION_RADIATION_HPP_
