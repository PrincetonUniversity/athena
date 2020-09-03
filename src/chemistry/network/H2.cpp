//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file H2.cpp
//  \brief implementation of functions in class ChemNetwork, using the simple
//  network for H2 formation and destruction.
//======================================================================================

// this class header
#include "H2.hpp"

//athena++ header
#include "network.hpp"
#include "../../scalars/scalars.hpp"
#include "../../parameter_input.hpp"       //ParameterInput
#include "../../mesh/mesh.hpp"
#include "../../hydro/hydro.hpp"
#include "../utils/chemistry_utils.hpp"
#include "../../defs.hpp"
#include "../../eos/eos.hpp"
#include "../utils/thermo.hpp"

//c++ header
#include <sstream>    // stringstream
#include <iostream>   // endl
#include <limits>    //inf

//constants
const Real ChemNetwork::kgr_ = 3e-17;

//species names
const std::string ChemNetwork::species_names[NSCALARS] = 
{"H", "H2"};

const int ChemNetwork::iH_ =
  ChemistryUtility::FindStrIndex(species_names, NSCALARS, "H");
const int ChemNetwork::iH2_ =
  ChemistryUtility::FindStrIndex(species_names, NSCALARS, "H2");

//flag for Cv
static bool is_const_Cv;


ChemNetwork::ChemNetwork(MeshBlock *pmb, ParameterInput *pin) {
	//number of species and a list of name of species
  pmy_spec_ = pmb->pscalars;
	pmy_mb_ = pmb;

	//set the parameters from input file
	xi_cr_ = pin->GetOrAddReal("chemistry", "xi_cr", 2e-16);
  kcr_ = xi_cr_ * 3.;
  //units
	unit_density_in_nH_ = pin->GetReal("chemistry", "unit_density_in_nH");
	unit_length_in_cm_ = pin->GetReal("chemistry", "unit_length_in_cm");
	unit_vel_in_cms_ = pin->GetReal("chemistry", "unit_vel_in_cms");
  unit_time_in_s_ = unit_length_in_cm_/unit_vel_in_cms_;
  unit_E_in_cgs_ = 1.67e-24 * 1.4 * unit_density_in_nH_
                           * unit_vel_in_cms_ * unit_vel_in_cms_;
  //set Cv: constant or H2 abundance dependent
  is_const_Cv = pin->GetOrAddBoolean("problem", "is_const_Cv", true);
}

ChemNetwork::~ChemNetwork() {}

void ChemNetwork::InitializeNextStep(const int k, const int j, const int i) {
  Real rho, rho_floor;
  //density
  rho = pmy_mb_->phydro->w(IDN, k, j, i);
  //apply density floor
  rho_floor = pmy_mb_->peos->GetDensityFloor();
  rho = (rho > rho_floor) ?  rho : rho_floor;
  //hydrogen atom number density
  nH_ =  rho * unit_density_in_nH_;
  return;
}

void ChemNetwork::RHS(const Real t, const Real y[NSCALARS], const Real ED,
                      Real ydot[NSCALARS]){
  const Real rate_cr = kcr_ * y[iH2_];
  const Real rate_gr = kgr_ * nH_ * y[iH_];
  ydot[iH2_] = rate_gr - rate_cr;
  ydot[iH_] = -2*rate_gr + 2*rate_cr;
	for (int i=0; i<NSCALARS; i++) {
    //return in code units
		ydot[i] *= unit_time_in_s_;
	}
  return;
}

Real ChemNetwork::Edot(const Real t, const Real y[NSCALARS], const Real ED){
  //isothermal
  if (!NON_BAROTROPIC_EOS) {
    return 0;
  }
  const Real x_He = 0.1;
  const Real x_e = 0.;
  Real x_H2;
  if (is_const_Cv) {
    x_H2 = 0;
  } else {
    x_H2 = y[iH2_];
  }
  const Real T_floor = 1.;//temperature floor for cooling
  //ernergy per hydrogen atom
  const Real E_ergs = ED * unit_E_in_cgs_ / nH_; 
  Real T = E_ergs / Thermo::CvCold(x_H2, x_He, x_e);
	if (T < T_floor) {
    return 0;
  }
  Real dEdt = - Thermo::alpha_GD_ * nH_ * sqrt(T) * T;
  //return in code units
  Real dEDdt = (dEdt * nH_ / unit_E_in_cgs_) * unit_time_in_s_;
  return dEDdt;
}
