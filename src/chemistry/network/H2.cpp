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
#include "../../defs.hpp"
#include "../../eos/eos.hpp"

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

ChemNetwork::ChemNetwork(MeshBlock *pmb, ParameterInput *pin) {
	//number of species and a list of name of species
  pmy_spec_ = pmb->pscalars;
	pmy_mb_ = pmb;

	//set the parameters from input file
	xi_CR_ = pin->GetOrAddReal("chemistry", "xi_CR", 2e-16);
  //units
	unit_density_in_nH_ = pin->GetReal("chemistry", "unit_density_in_nH");
	unit_length_in_cm_ = pin->GetReal("chemistry", "unit_length_in_cm");
	unit_vel_in_cms_ = pin->GetReal("chemistry", "unit_vel_in_cms");
  unit_time_in_s_ = unit_length_in_cm_/unit_vel_in_cms_;
  unit_E_in_cgs_ = 1.67e-24 * 1.4 * unit_density_in_nH_
                           * unit_vel_in_cms_ * unit_vel_in_cms_;
}

ChemNetwork::~ChemNetwork() {}
