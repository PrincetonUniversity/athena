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

// Athena++ classes headers
#include "network.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"

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
  std::string network_dir_;

  //physical quantities
	Real xi_cr_; //primary CRIR in s-1 H-1, read from input file, default 2e-16.
	Real nH_; //density, updated at InitializeNextStep from hydro variable

	//units 
	Real unit_density_in_nH_; //read from input
	Real unit_length_in_cm_; //read from input
	Real unit_vel_in_cms_; //read from input
	Real unit_time_in_s_; //from length and velocity units
  //unit of energy density, in erg cm-3, from density and velocity units
  Real unit_E_in_cgs_; 
};

#endif // KIDA_HPP
