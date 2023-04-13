#ifndef G14Sod_HPP
#define G14Sod_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file G14Sod.hpp
//  \brief implementation of functions in Grassi 2014 Fig. 21
//======================================================================================

//c++ headers
#include <string> //std::string

// Athena++ classes headers
#include "network.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
class ChemNetwork : public NetworkWrapper {
  //OutputProperties in problem generator called by Mesh::UserWorkAfterLoop.
  friend class Mesh; 
  //It would be convenient to know the species names in
  //initialization of chemical species in problem
  friend class MeshBlock; 
public:
	static const int ngs_     = 1;

  ChemNetwork(MeshBlock *pmb, ParameterInput *pin);
  ~ChemNetwork();

	//a list of species name, used in output
	static const std::string species_names[NSCALARS];

  //Set the rates of chemical reactions, eg. through density.
  //k, j, i are the corresponding index of the grid
  void InitializeNextStep(const int k, const int j, const int i);

  //RHS: right-hand-side of ODE. dy/dt = ydot(t, y). Here y are the abundance
  //of species. details see CVODE package documentation.
  //all input/output variables are in code units
  void RHS(const Real t, const Real y[NSCALARS], const Real ED,
           Real ydot[NSCALARS]);
  
  //energy equation dE/dt, all input/output variables are in code units
  //(ED is the energy density)
  Real Edot(const Real t, const Real y[NSCALARS], const Real ED);


  void GetGhostSpecies(const Real *y, Real yghost[NSCALARS+ngs_]);
  void UpdateRates(const Real y[NSCALARS+ngs_], const Real E);
  void OutputRates(FILE *pf) const;
private:
  PassiveScalars *pmy_spec_;
	MeshBlock *pmy_mb_;

	//constants
	static const int n_2body_ = 20;
	static const std::string ghost_species_names_[ngs_];
	const Real gamma = 5./3.;
	const Real mu = 1.25;
	std::string species_names_all_[NSCALARS+ngs_];//all species
	
	Real nH_;   //hydrogen number, updated at InitializeNextStep
	Real rho;   //density, updated at InitializeNextStep

  //	Real indi,indj,indk;
  //	static const Real gm1 = 2/3; // factor gamma - 1
	//units of density and radiation
	Real unit_density_in_nH_;
	Real unit_length_in_cm_;
	Real unit_vel_in_cms_;
	Real unit_time_in_s_;
	Real unit_E_in_cgs_;//unit of energy density, in erg cm-3
	Real rate_in_cgs_ ; // convert the rate unit from SI to cgs
	Real temperature_;
	Real temp_max_heat_; 
	Real temp_min_cool_; 
	Real temp_min_rates_; 
	Real temp_max_rates_; 
	//parameters of the netowork
	Real xHplus_  ;
	Real xHe_     ;
	Real xHeplus_ ;
	Real xHe2plus_;
	Real xHmin_   ;
	Real xH2_     ;
	Real xH2P_    ;
	//index of species
	static const int iH_;
	static const int iHplus_;
	static const int iHe_;
	static const int iHeplus_;
	static const int iHe2plus_;
	static const int iHmin_;
	static const int iH2_;
	static const int iH2plus_;
	static const int ige_;
	static const int igr_;
	//-------------------chemical network---------------------
	//2body reactions
	static const int in2body1_[n_2body_];
	static const int in2body2_[n_2body_];
	static const int out2body1_[n_2body_];
	static const int out2body2_[n_2body_];
	static const Real stoich_in2body1[n_2body_];
	static const Real stoich_in2body2[n_2body_];
	static const Real stoich_out2body1[n_2body_];
	static const Real stoich_out2body2[n_2body_];
	Real k2body_[n_2body_]; //rates for 2 body reacrtions.
};

#endif // G14SodHPP
