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
//! \file kida.cpp
//  \brief implementation of functions in class ChemNetwork, using the simple
//  network for kida style network files.
//======================================================================================

// this class header
#include "kida.hpp"

//athena++ header
#include "network.hpp"
#include "../../scalars/scalars.hpp"
#include "../../parameter_input.hpp"       //ParameterInput
#include "../../mesh/mesh.hpp"
#include "../../hydro/hydro.hpp"
#include "../utils/chemistry_utils.hpp"
#include "../utils/kida_species.hpp"
#include "../../utils/string_utils.hpp"
#include "../../defs.hpp"
#include "../../eos/eos.hpp"

//c++ header
#include <sstream>    // stringstream
#include <iostream>   // endl
#include <limits>    //inf
#include <fstream>   //file()
#include <stdio.h>    // c style file

ChemNetwork::ChemNetwork(MeshBlock *pmb, ParameterInput *pin) {
	//number of species and a list of name of species
  pmy_spec_ = pmb->pscalars;
	pmy_mb_ = pmb;

	//set the parameters from input file
	xi_cr_ = pin->GetOrAddReal("chemistry", "xi_cr", 2e-16);
  //units
	unit_density_in_nH_ = pin->GetReal("chemistry", "unit_density_in_nH");
	unit_length_in_cm_ = pin->GetReal("chemistry", "unit_length_in_cm");
	unit_vel_in_cms_ = pin->GetReal("chemistry", "unit_vel_in_cms");
  unit_time_in_s_ = unit_length_in_cm_/unit_vel_in_cms_;
  unit_E_in_cgs_ = 1.67e-24 * 1.4 * unit_density_in_nH_
                           * unit_vel_in_cms_ * unit_vel_in_cms_;
  //folder of the network
  network_dir_ = pin->GetString("chemistry", "network_dir");

  //read in the species
  std::string species_file_name = network_dir_ + "/species.dat";
  std::ifstream species_file(species_file_name);
  if (!species_file) {
    std::stringstream msg; //error message
    msg << "### FATAL ERROR in ChemNetwork constructor [ChemNetwork]" << std::endl
      << "Cannot open file" << species_file_name << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  std::string line;
  int nline = 0;
  while (getline(species_file, line)) {
    //trim white spaces
    StringUtils::trim(line);
    //skip blank lines and comments
    if(line.empty() || (line.find("!") == 0)) {
        continue;
    }
    KidaSpecies si(line, nline);
    species_names[nline] = si.name;
    nline++;
    //test: print each species.
    std::cout << "name=" << si.name << ", index=" << si.index << ", charge="
      << si.charge_ << std::endl;
    std::cout << "atom_count_ = ";
    for (int i=0; i < si.natom_; i++) {
      std::cout << si.atom_count_[i] << " ";
    }
    std::cout << std::endl;
  }
  if (nline != NSCALARS) {
    std::stringstream msg; //error message
    msg << "### FATAL ERROR in ChemNetwork constructor [ChemNetwork]" << std::endl
      << "number of species in species.dat does not match the number of scalars (" 
      << NSCALARS << ")" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
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
  for (int i=0; i<NSCALARS; i++) {
    ydot[i] = 0;
  }
  return;
}

Real ChemNetwork::Edot(const Real t, const Real y[NSCALARS], const Real ED){
  const Real dEDdt = 0;
  return dEDdt;
}
