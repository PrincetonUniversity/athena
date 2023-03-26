//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file kida_species.cpp
//! \brief implementation of functions in class KidaSpecies

//this class header
#include "kida_species.hpp"

//c++ header
#include <iostream>   // endl
#include <stdexcept>  // std::runtime_error()
#include <vector>

//athena++ header
#include "../../units/units.hpp"
#include "../../utils/string_utils.hpp"

//atom masses
const Real KidaSpecies::ma_atom_[natom_] =
{1., 4., 12., 14., 16., 28., 32., 55.8, 23., 24.3, 35.5, 31., 19. };
//H  He   C    N    O   Si   S    Fe    Na   Mg    Cl    P    F

//----------------------------------------------------------------------------------------
//! constructor
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
    mass_ += static_cast<float>(atom_count_[i]) * ma_atom_[i]
      * Constants::hydrogen_mass_cgs;
  }
  //electron mass
  mass_ += static_cast<float>(-charge_) * ChemistryUtility::me;
}

//----------------------------------------------------------------------------------------
//! \fn void KidaSpecies::SetMass(Real mass)
//! \brief set the mass of species. used for grains.
void KidaSpecies::SetMass(Real mass) {
  mass_ = mass;
  return;
}

