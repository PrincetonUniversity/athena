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
//! \file kida_species.cpp
//  \brief implementation of functions in class KidaSpecies
//======================================================================================
#include "kida_species.hpp"
#include "../../utils/string_utils.hpp"
#include <iostream>   // endl
#include <stdexcept>  // std::runtime_error()

//atom masses
const Real KidaSpecies::ma_atom_[natom_] = 
{1., 4., 12., 14., 16., 28., 32., 55.8, 23., 24.3, 35.5, 31., 19. };
//H  He   C    N    O   Si   S    Fe    Na   Mg    Cl    P    F

KidaSpecies::KidaSpecies(std::string line, int index) :
  index(index) {
  std::stringstream msg; //error message
  std::vector<std::string> fields = StringUtils::split(line, ' ');
  if ( fields.size() != (natom_ + 2) ) {
    msg << "### FATAL ERROR in KidaSpecies constructor [KidaSpecies]" << std::endl
      << "number of fields (" << fields.size() << ")does not match " 
      << natom_ + 2 << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  name = fields[0];
  charge_ = std::stoi(fields[1]);
  mass_ = 0.;
  for (int i=0; i < natom_; i++) {
    atom_count_[i] = std::stoi(fields[i+2]);
    mass_ += float(atom_count_[i]) * ma_atom_[i] * ChemistryUtility::mH;
  }
  //electron mass
  mass_ += float(-charge_) * ChemistryUtility::me;
}

