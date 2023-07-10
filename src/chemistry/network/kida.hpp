#ifndef CHEMISTRY_NETWORK_KIDA_HPP_
#define CHEMISTRY_NETWORK_KIDA_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file kida.hpp
//! \brief definitions for a general chemical network in KIDA format

// C headers

// C++ headers
#include <array>
#include <map>
#include <string>
#include <vector>

// Athena++ classes headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../utils/kida_reaction.hpp"
#include "../utils/kida_species.hpp"
#include "network.hpp"

class Units;

//reaction types
enum class ReactionType {none, cr, crp, photo, twobody, twobodytr,
                         grain_implicit, special, grain_collision};

//! \class ChemNetwork
//! \brief General chemical network in KIDA format
class ChemNetwork : public NetworkWrapper {
  // It would be convenient to know the species names in
  // initialization of chemical species in problem
  friend class MeshBlock;
 public:
  ChemNetwork(MeshBlock *pmb, ParameterInput *pin);
  ~ChemNetwork();

  // a list of species name, used in output
  std::array<std::string, NSPECIES> species_names;

  void InitializeNextStep(const int k, const int j, const int i);
  void RHS(const Real t, const Real *y, const Real ED, Real *ydot);
  Real Edot(const Real t, const Real *y, const Real ED);
  void Jacobian_isothermal(const Real t, const Real *y, const Real *ydot,
                           AthenaArray<Real> &jac);

 private:
  PassiveScalars *pmy_spec_;
  MeshBlock *pmy_mb_;

  std::map<std::string, int> ispec_map_;
  std::string network_dir_;
  std::vector<KidaSpecies> species_;
  std::vector<KidaReaction> reactions_;
  AthenaArray<int> id_ices_;// species index for ices
  std::vector<int> id_2bodytr_;
  int nr_; // number of reactions
  int nices_;// number of ice species

  // physical quantities
  static constexpr Real xHe_ = 0.1; // Helium abundance
  // whether to cap temperature if the reaction is outside of the temperature range
  // only for 2body reactions. Default is false, which means extrapolation
  bool is_Tcap_2body_;
  // whether to fix dust metallicity for implicit grain assisted
  // reactions (opposed to reading from the initial grain abundance from the
  // input file)
  bool is_fixed_Zd_;
  // whether to fix PAH abundance for PE heating and Rec cooling
  // if so, use PE heating and Rec cooling from Draine & Weingartner
  // if not, use PE heating and Rec cooling from Wolfire+ 2003
  bool is_fixed_PAH_;
  // flag for in isothermal EOS only temperature related rates are
  // calculated only once in the beginning
  bool flag_T_rates_;
  Real Z_g_; // gas metallicity relative to solar, default 1.
  Real Z_d_; // larger dust grain metallicity relative to solar, default 1.
  Real Z_PAH_; // PAH abundance, for PE heating and Rec cooling
  Real phi_PAH_; // PAH recombination efficiency, default 0.4
  Real a_d_; // size of the dust grain in cm, default 1e-5 (0.1 micron)
  Real rho_d_; // density of grain in cgs, default 2 g/cm3
  Real m_d_; // mass of the dust grain in g
  Real x_d_; // relative abundance of all dust
  Real nH_; // density, updated at InitializeNextStep from hydro variable
  Real o2pH2_;// ortho to para H2 ratio, default 3:1
  Real Yi_;// Yield for crp desorption, default 1e-3
  Real temperature_; // temperature of the gas if isothermal
  Real temp_min_rates_; // temperature floor for reaction rates
  Real temp_min_cool_; // temperature minimum for cooling
  Real temp_max_cool_nm_; // temperature maximum for neutral medium cooling
  Real temp_dust_thermo_; // dust temperature for dust thermo cooling and desorption

  // reaction constants
  // special rates index map
  int id7max_;
  AthenaArray<int> id7map_;
  AthenaArray<ReactionType> id7type_;
  // direct cosmic-ray ionization
  int n_cr_;
  AthenaArray<int> incr_;
  AthenaArray<int> outcr1_;
  AthenaArray<int> outcr2_;
  AthenaArray<int> outcr3_;
  AthenaArray<Real> kcr_base_;
  AthenaArray<Real> kcr_;
  // index for cr ionization rates for CR heating
  int icr_H_;
  int icr_H2_;
  int icr_He_;
  // cosmic-ray induced photo ionization
  int n_crp_;
  AthenaArray<int> incrp_;
  AthenaArray<int> outcrp1_;
  AthenaArray<int> outcrp2_;
  AthenaArray<Real> kcrp_base_;
  AthenaArray<Real> kcrp_;
  // photo reactions
  int n_ph_;
  AthenaArray<int> inph_;
  AthenaArray<int> outph1_;
  AthenaArray<int> outph2_;
  // reactant species name map, for radiation calculation
  std::map<std::string, int> smap_ph_;
  AthenaArray<Real> kph_base_;
  AthenaArray<Real> kph_avfac_;
  AthenaArray<Real> kph_;
  // index for H2 photodissociation for H2 UV pumping and dissociation heating
  int iph_H2_;
  // 2body reactions
  int n_2body_;
  static constexpr int n_in2body_ = 2;
  static constexpr int n_out2body_ = 4;
  AthenaArray<int> in2body_;
  AthenaArray<int> out2body_;
  AthenaArray<int> frml_2body_;
  AthenaArray<Real> a2body_; // alpha
  AthenaArray<Real> b2body_; // beta
  AthenaArray<Real> c2body_; // gamma
  AthenaArray<Real> Tmin_2body_; // minimum temperature for reaction rates
  AthenaArray<Real> Tmax_2body_; // maximum temperature for reaction rates
  AthenaArray<Real> k2body_;
  int i2body_H2_H_; // index for H2+H collisional dissociation, for cooling
  int i2body_H2_H2_; // index for H2+H2 collisional dissociation, for cooling
  int i2body_H_e_; // index for H+e collisional ionization, for cooling
  // 2body reactions with temperature dependent rates
  // they have to be the same reaction, same ID, and arranged next to each other
  // with ascending temperature ranges
  int n_2bodytr_;
  static constexpr int n_range_ = 3; // maximum number of temperature ranges
  AthenaArray<int> in2bodytr1_;
  AthenaArray<int> in2bodytr2_;
  AthenaArray<int> out2bodytr1_;
  AthenaArray<int> out2bodytr2_;
  AthenaArray<int> out2bodytr3_;
  AthenaArray<int> nr_2bodytr_;
  AthenaArray<int> frml_2bodytr_;
  AthenaArray<Real> a2bodytr_; // alpha
  AthenaArray<Real> b2bodytr_; // beta
  AthenaArray<Real> c2bodytr_; // gamma
  AthenaArray<Real> Tmin_2bodytr_; // minimum temperature for reaction rates
  AthenaArray<Real> Tmax_2bodytr_; // maximum temperature for reaction rates
  AthenaArray<Real> k2bodytr_;
  // grain assisted reactions
  int n_gr_;
  AthenaArray<int> ingr1_;
  AthenaArray<int> ingr2_;
  AthenaArray<int> outgr_;
  AthenaArray<int> frml_gr_; // 7: special, 10: desorption
  AthenaArray<Real> kgr_;
  AthenaArray<Real> TDgr_;   // desorption/binding energy/temperature in K
  AthenaArray<Real> nu0gr_;  // vibrational frequency for desorption
  int igr_H_;  // index for gr fromation of H2 for its heating
  // special reactions
  int n_sr_;
  static constexpr int n_insr_ = 3;
  static constexpr int n_outsr_ = 5;
  AthenaArray<int> insr_;
  AthenaArray<int> outsr_;
  AthenaArray<Real> ksr_;
  // grain collision: grain - electron/ion reactions
  int n_gc_;
  static constexpr int n_ingc_ = 2;
  static constexpr int n_outgc_ = 3;
  AthenaArray<int> ingc_;
  AthenaArray<int> outgc_;
  AthenaArray<int> nu_gc_; // nu={0, -1} for rate formula
  // effective rate at 1K: pi a_g^2 s_i sqrt( 8k_B/(pi m_i) )*branch_ratio
  AthenaArray<Real> r1_gc_;
  AthenaArray<Real> t1_gc_; // tau at 1K: a_g k_B/qi^2
  AthenaArray<Real> kgc_;

  // radiation related reactions and variables
  int n_freq_;
  int index_gpe_;
  int index_cr_;
  AthenaArray<Real> rad_;

  // parameters related to CO cooling
  // these are needed for LVG approximation
  Real Leff_CO_max_;  // maximum effective length in cm for CO cooling
  Real gradv_;        // absolute value of velocity gradient in cgs, >0

  // private functions
  void InitializeReactions();
  void UpdateRates(const Real *y, const Real E);
  ReactionType SortReaction(KidaReaction* pr) const;
  void CheckReaction(KidaReaction reaction);
  void PrintProperties() const;
  void OutputRates(FILE *pf) const;
  void OutputJacobian(FILE *pf, const AthenaArray<Real> &jac) const;
  void Jacobian_isothermal_numerical(const Real t, const Real *y,
                                     const Real *ydot,
                                     AthenaArray<Real> &jac);
  void SetGrad_v(const int k, const int j, const int i);
  void UpdateRatesSpecial(const Real *y, const Real E);
  void UpdateJacobianSpecial(const Real *y, const Real E, AthenaArray<Real> &jac);
};

#endif // CHEMISTRY_NETWORK_KIDA_HPP_
