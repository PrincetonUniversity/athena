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
enum class ReactionType {cr, crp, photo, twobody, grain, special};

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
  int nr_; //number of reactions

  //physical quantities
  const Real xHe_ = 0.1;//Helium abundance
	Real zdg_; //dust and gas metallicity relative to solar, default 1.
	Real nH_; //density, updated at InitializeNextStep from hydro variable
  Real temperature_; //temperature of the gas if isothermal 
  Real temp_min_rates_; //temperature floor for reaction rates

	//units 
	Real unit_density_in_nH_; //read from input
	Real unit_length_in_cm_; //read from input
	Real unit_vel_in_cms_; //read from input
	Real unit_time_in_s_; //from length and velocity units
  Real unit_radiation_in_draine1987_; //FUV radiation unit
  //unit of energy density, in erg cm-3, from density and velocity units
  Real unit_E_in_cgs_; 

  //reaction constants
  //direct cosmic-ray ionization
  int n_cr_;
  std::vector<int> incr_;
  std::vector<int> outcr1_;
  std::vector<int> outcr2_;
  std::map<int, int> idmap_cr_;
  std::vector<Real> kcr_base_;
  std::vector<Real> kcr_;
  //cosmic-ray induced photo ionization
  int n_crp_;
  std::vector<int> incrp_;
  std::vector<int> outcrp1_;
  std::vector<int> outcrp2_;
  std::map<int, int> idmap_crp_;
  std::vector<Real> kcrp_base_;
  std::vector<Real> kcrp_;
  //photo reactions
  int n_ph_;
  std::vector<int> inph_;
  std::vector<int> outph1_;
  std::vector<int> outph2_;
  //reactant species name map, for radiation calculation
  std::map<std::string, int> smap_ph_; 
  std::vector<Real> kph_base_;
  std::vector<Real> kph_avfac_;
  std::vector<Real> kph_;
  //2body reactions
  int n_2body_;
  std::vector<int> in2body1_;
  std::vector<int> in2body2_;
  std::vector<int> out2body1_;
  std::vector<int> out2body2_;
  std::vector<int> out2body3_;
  std::vector<int> frml_2body_;
  std::map<int, int> idmap_2body_;
  std::vector<Real> a2body_; //alpha
  std::vector<Real> b2body_; //beta
  std::vector<Real> c2body_; //gamma
  std::vector<Real> k2body_;
  //grain assisted reactions
  int n_gr_;
  std::vector<int> ingr1_;
  std::vector<int> ingr2_;
  std::vector<int> outgr_;
  std::map<int, int> idmap_gr_;
  std::vector<Real> kgr_;
  //special reactions
  int n_sr_;
  const int n_insr_ = 3;
  const int n_outsr_ = 5;
  std::vector<KidaReaction*> pr_sr_;
  AthenaArray<int> insr_;
  AthenaArray<int> outsr_;
  std::map<int, int> idmap_sr_;
  std::vector<Real> ksr_;

  //radiation related reactions and variables
  int n_freq_;
  int index_gpe_;
  int index_cr_;
	AthenaArray<Real> rad_;

  //functions
  void InitializeReactions(); //set up coefficients of reactions
  void UpdateRates(const Real y[NSCALARS], const Real E); //reaction rates
  void UpdateRatesSpecial(const Real y[NSCALARS], const Real E); //formula = 7
  ReactionType SortReaction(KidaReaction* pr) const;
  void CheckReaction(KidaReaction reaction);
  void PrintProperties() const; //print out reactions and rates, for debug
  void OutputRates(FILE *pf) const;//output reaction rates

};


#endif // KIDA_HPP
