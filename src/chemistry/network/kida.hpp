#ifndef KIDA_HPP
#define KIDA_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file H2.hpp
//  \brief definitions for a very simple chemical network with H2 formation on grains,
//  and H2 distruction by CR. This has an analytic solution.
//======================================================================================
//c++ headers
#include <string> //std::string
#include <vector> // vector container
#include <map>    //map

// Athena++ classes headers
#include "network.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../utils/kida_species.hpp"
#include "../utils/kida_reaction.hpp"

//reaction types
enum class ReactionType {none, cr, crp, photo, twobody, twobodytr, grain, special};

//! \class ChemNetwork
//  \brief Chemical Network that defines the reaction rates between species.
//  Note: This is a template for chemistry network.
//  When implementing a new chemistry network, all public functions should be
//  in the same form.
//  The internal calculations are in cgs units. The input and 
//  return of RHS and Edot must be in code units.
class ChemNetwork : public NetworkWrapper {
  //It would be convenient to know the species names in
  //initialization of chemical species in problem
  friend class MeshBlock; 
public:
  ChemNetwork(MeshBlock *pmb, ParameterInput *pin);
  ~ChemNetwork();

	//a list of species name, used in output
	std::string species_names[NSCALARS];

	//Set the rates of chemical reactions, eg. through density and radiation field.
  //k, j, i are the corresponding index of the grid
  void InitializeNextStep(const int k, const int j, const int i);

  //RHS: right-hand-side of ODE. dy/dt = ydot(t, y). Here y are the abundance
  //of species. details see CVODE package documentation.
  //all input/output variables are in code units
  void RHS(const Real t, const Real y[NSCALARS], const Real ED,
           Real ydot[NSCALARS]);
  
  //energy equation dED/dt, all input/output variables are in code units
  //(ED is the energy density)
  Real Edot(const Real t, const Real y[NSCALARS], const Real ED);

private:
  PassiveScalars *pmy_spec_;
	MeshBlock *pmy_mb_;

  std::map<std::string, int> ispec_map_;
  std::string network_dir_;
  std::vector<KidaSpecies> species_;
  std::vector<KidaReaction> reactions_;
  std::vector<int> id_2bodytr_;
  int nr_; //number of reactions

  //physical quantities
  const Real xHe_ = 0.1;//Helium abundance
  //whether to cap temperature if the reaction is outside of the temperature range
  //only for 2body reactions. Default is false, which means extrapolation
  bool is_Tcap_2body_; 
	Real zdg_; //dust and gas metallicity relative to solar, default 1.
	Real nH_; //density, updated at InitializeNextStep from hydro variable
  Real o2pH2_;//ortho to para H2 ratio, default 3:1
  Real temperature_; //temperature of the gas if isothermal 
  Real temp_min_rates_; //temperature floor for reaction rates
  Real temp_min_cool_; //temperature minimum for cooling
  Real temp_dust_thermo_; //dust temperature for dust thermo cooling 

	//units 
	Real unit_density_in_nH_; //read from input
	Real unit_length_in_cm_; //read from input
	Real unit_vel_in_cms_; //read from input
	Real unit_time_in_s_; //from length and velocity units
  Real unit_radiation_in_draine1987_; //FUV radiation unit
  //unit of energy density, in erg cm-3, from density and velocity units
  Real unit_E_in_cgs_; 

  //reaction constants
  //special rates index map
  int id7max_;
  AthenaArray<int> id7map_;
  AthenaArray<ReactionType> id7type_;
  //direct cosmic-ray ionization
  int n_cr_;
  AthenaArray<int> incr_;
  AthenaArray<int> outcr1_;
  AthenaArray<int> outcr2_;
  AthenaArray<int> outcr3_;
  AthenaArray<Real> kcr_base_;
  AthenaArray<Real> kcr_;
  //index for cr ionization rates for CR heating
  int icr_H_;  
  int icr_H2_;
  int icr_He_;
  //cosmic-ray induced photo ionization
  int n_crp_;
  AthenaArray<int> incrp_;
  AthenaArray<int> outcrp1_;
  AthenaArray<int> outcrp2_;
  AthenaArray<Real> kcrp_base_;
  AthenaArray<Real> kcrp_;
  //photo reactions
  int n_ph_;
  AthenaArray<int> inph_;
  AthenaArray<int> outph1_;
  AthenaArray<int> outph2_;
  //reactant species name map, for radiation calculation
  std::map<std::string, int> smap_ph_; 
  AthenaArray<Real> kph_base_;
  AthenaArray<Real> kph_avfac_;
  AthenaArray<Real> kph_;
  //index for H2 photodissociation for H2 UV pumping and dissociation heating
  int iph_H2_; 
  //2body reactions
  int n_2body_;
  AthenaArray<int> in2body1_;
  AthenaArray<int> in2body2_;
  AthenaArray<int> out2body1_;
  AthenaArray<int> out2body2_;
  AthenaArray<int> out2body3_;
  AthenaArray<int> out2body4_;
  AthenaArray<int> frml_2body_;
  AthenaArray<Real> a2body_; //alpha
  AthenaArray<Real> b2body_; //beta
  AthenaArray<Real> c2body_; //gamma
  AthenaArray<Real> Tmin_2body_; //minimum temperature for reaction rates
  AthenaArray<Real> Tmax_2body_; //maximum temperature for reaction rates
  AthenaArray<Real> k2body_;
  int i2body_H2_H_; //index for H2+H collisional dissociation, for cooling
  int i2body_H2_H2_; //index for H2+H2 collisional dissociation, for cooling
  int i2body_H_e_; //index for H+e collisional ionization, for cooling
  //2body reactions with temperature dependent rates
  //they have to be the same reaction, same ID, and arranged next to each other
  //with accending temperature ranges
  int n_2bodytr_;
  const int n_range_ = 3; //maximum number of temperature ranges
  AthenaArray<int> in2bodytr1_;
  AthenaArray<int> in2bodytr2_;
  AthenaArray<int> out2bodytr1_;
  AthenaArray<int> out2bodytr2_;
  AthenaArray<int> out2bodytr3_;
  AthenaArray<int> nr_2bodytr_;
  AthenaArray<int> frml_2bodytr_;
  AthenaArray<Real> a2bodytr_; //alpha
  AthenaArray<Real> b2bodytr_; //beta
  AthenaArray<Real> c2bodytr_; //gamma
  AthenaArray<Real> Tmin_2bodytr_; //minimum temperature for reaction rates
  AthenaArray<Real> Tmax_2bodytr_; //maximum temperature for reaction rates
  AthenaArray<Real> k2bodytr_;
  //grain assisted reactions
  int n_gr_;
  AthenaArray<int> ingr1_;
  AthenaArray<int> ingr2_;
  AthenaArray<int> outgr_;
  AthenaArray<Real> kgr_;
  int igr_H_;//index for gr fromation of H2 for its heating
  //special reactions
  int n_sr_;
  const int n_insr_ = 3;
  const int n_outsr_ = 5;
  AthenaArray<int> insr_;
  AthenaArray<int> outsr_;
  AthenaArray<Real> ksr_;

  //radiation related reactions and variables
  int n_freq_;
  int index_gpe_;
  int index_cr_;
	AthenaArray<Real> rad_;

	//parameters related to CO cooling
	//these are needed for LVG approximation
	Real Leff_CO_max_; //maximum effective length in cm for CO cooling
	Real gradv_; //abosolute value of velocity gradient in cgs, >0

  //functions
  void InitializeReactions(); //set up coefficients of reactions
  void UpdateRates(const Real y[NSCALARS], const Real E); //reaction rates
  void UpdateRatesSpecial(const Real y[NSCALARS], const Real E); //formula = 7
  ReactionType SortReaction(KidaReaction* pr) const;
  void CheckReaction(KidaReaction reaction);
  void PrintProperties() const; //print out reactions and rates, for debug
  void OutputRates(FILE *pf) const;//output reaction rates
  //set gradients of v and nH for CO cooling
  void SetGrad_v(const int k, const int j, const int i); 

};


#endif // KIDA_HPP
