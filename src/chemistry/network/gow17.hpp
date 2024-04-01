#ifndef CHEMISTRY_NETWORK_GOW17_HPP_
#define CHEMISTRY_NETWORK_GOW17_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gow17.hpp
//! \brief definitions for chemical network in Gong, Osriker and Wolfire 2017

// C headers

// C++ headers
#include <array>
#include <string>

// Athena++ classes headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "network.hpp"

class Units;

//! \class ChemNetwork
//! \brief Chemical Network that defines the reaction rates between species.
//!
//!  Note: This is a template for chemistry network.
//!  When implementing a new chemistry network, all public functions should be
//!  in the same form.
//!  The internal calculations are in cgs units. The input and
//!  return of RHS and Edot must be in code units.
class ChemNetwork : public NetworkWrapper {
  friend class ChemRadiation; // number of frequencies n_freq_
  friend class ChemRadIntegrator; // shielding
  // OutputProperties in problem generator called by Mesh::UserWorkAfterLoop.
  friend class Mesh;
  // It would be convenient to know the species names in
  // initialization of chemical species in problem
  friend class MeshBlock;
 public:
  ChemNetwork(MeshBlock *pmb, ParameterInput *pin);
  ~ChemNetwork();

  // a list of species name, used in output
  static const std::array<std::string, NSPECIES> species_names;

  void InitializeNextStep(const int k, const int j, const int i);
  void RHS(const Real t, const Real *y, const Real ED, Real *ydot);
  Real Edot(const Real t, const Real *y, const Real ED);

 private:
  PassiveScalars *pmy_spec_;
  MeshBlock *pmy_mb_;

  // constants
  static constexpr int ngs_ = 6;
  static constexpr int n_cr_ = 7;
  static constexpr int n_2body_ = 31;
  static constexpr int n_ph_ = 6;
  static constexpr int n_gr_ = 5;
  static constexpr int nE_ = 15;
  static constexpr int n_freq_ = n_ph_ + 2;
  static constexpr int index_gpe_ = n_ph_;
  static constexpr int index_cr_ = n_ph_ + 1;
  // other variables
  static const std::array<std::string, ngs_> ghost_species_names_;
  std::array<std::string, NSCALARS+ngs_> species_names_all_;
  Real nH_; // density, updated at InitializeNextStep
  Real temperature_;
  Real temp_max_heat_;
  Real temp_min_cool_;
  Real temp_max_cool_nm_;
  Real temp_min_rates_;
  Real temp_max_rates_;
  bool is_H2_rovib_cooling_; // whether to include H2 rovibrational cooling
  // CR shielding
  bool is_cr_shielding_;
  // parameters of the netowork
  Real zdg_;
  Real xHe_;
  Real xC_std_;
  Real xO_std_;
  Real xSi_std_;
  Real xC_;
  Real xO_;
  Real xSi_;
  // index of species
  static const int iHeplus_;
  static const int iOHx_;
  static const int iCHx_;
  static const int iCO_;
  static const int iCplus_;
  static const int iHCOplus_;
  static const int iH2_;
  static const int iHplus_;
  static const int iH3plus_;
  static const int iH2plus_;
  static const int iOplus_;
  static const int iSiplus_;
  static const int iE_;
  // index of ghost species
  static const int igSi_;
  static const int igC_;
  static const int igO_;
  static const int igHe_;
  static const int ige_;
  static const int igH_;
  // species needed for calculating shielding
  static constexpr int n_cols_ = 4;
  static constexpr int iNHtot_ = 0;
  static constexpr int iNH2_ = 1;
  static constexpr int iNCO_ = 2;
  static constexpr int iNC_ = 3;
  //-------------------chemical network---------------------
  // number of different reactions
  // cosmic ray reactions
  static const int icr_H2_;
  static const int icr_He_;
  static const int icr_H_;
  static const int incr_[n_cr_]; // reactant
  static const int outcr_[n_cr_]; // product
  static const Real kcr_base_[n_cr_]; // coefficents of rates relative to H
  Real kcr_[n_cr_]; // rates for cosmic-ray reactions  // NOLINT (runtime/arrays)
  // 2body reactions
  static const int i2body_H2_H;
  static const int i2body_H2_H2;
  static const int i2body_H_e;
  static const int in2body1_[n_2body_];
  static const int in2body2_[n_2body_];
  static const int out2body1_[n_2body_];
  static const int out2body2_[n_2body_];
  static const Real k2Texp_[n_2body_]; // exponent of temp dependance
  static const Real k2body_base_[n_2body_]; // base rate coefficents
  Real k2body_[n_2body_]; // rates for 2 body reactions  // NOLINT (runtime/arrays)
  static const Real A_kCHx_;
  static const Real n_kCHx_;
  static const Real c_kCHx_[4];
  static const Real Ti_kCHx_[4];
  // (15) H2 + *H -> 3 *H  and (16) H2 + H2 -> H2 + 2 *H
  // temperature above which collisional dissociation (15), (16) and (17)
  // will be importatant k>~1e-30
  static const Real temp_coll_;
  // photo reactions
  // radiation related variables, used in ChemRadIntegrator
  // radiation field in unit of Draine 1987 field (G0), updated at InitializeNextStep
  // vector: [Gph, GPE]
  Real rad_[n_freq_];   // NOLINT (runtime/arrays)
  Real kph_[n_ph_]; // rates for photo-reactions  // NOLINT (runtime/arrays)
  static const int iph_C_;
  static const int iph_CHx_;
  static const int iph_CO_;
  static const int iph_OHx_;
  static const int iph_H2_;
  static const int iph_Si_;
  static const int inph_[n_ph_];
  static const int outph1_[n_ph_];
  static const Real kph_base_[n_ph_];   // base rate of photo reaction
  static const Real kph_avfac_[n_ph_];  // exponent factor in front of AV
  // grain assisted reactions
  static const int igr_H_;
  static const int ingr_[n_gr_];
  static const int outgr_[n_gr_];
  Real kgr_[n_gr_]; // rates for grain assisted reactions  // NOLINT (runtime/arrays)
  // constants for H+, C+, He+ grain recombination.
  // See Draine ISM book table 14.9 page 158, or Weingartner+Draine2001.
  static const Real cHp_[7];
  static const Real cCp_[7];
  static const Real cHep_[7];
  static const Real cSip_[7];
  // factor to calculate psi in H+ recombination on grain
  Real psi_gr_fac_;
  // parameters related to CO cooling
  // these are needed for LVG approximation
  Real gradv_; // absolute value of velocity gradient in cgs, >0
  Real Leff_CO_max_; // maximum effective length for CO cooling
  // a small number to avoid divide by zero.
  static const Real small_;

  // private functions
  void GetGhostSpecies(const Real *y, Real *yall);
  Real CII_rec_rate_(const Real temp);
  void UpdateRates(const Real *y, const Real E);
  void OutputRates(FILE *pf) const;
  Real GetStddev(Real arr[], const int len);
  void SetGrad_v(const int k, const int j, const int i);
};

#endif // CHEMISTRY_NETWORK_GOW17_HPP_
