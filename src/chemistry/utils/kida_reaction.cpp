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
//! \file kida_reaction.cpp
//  \brief implementation of functions in class KidaReaction
//======================================================================================
#include "kida_reaction.hpp"
#include "../../utils/string_utils.hpp"
#include <iostream>   // endl
#include <stdexcept>  // std::runtime_error()
#include <string>
#include <stdio.h> //sscanf()

KidaReaction::KidaReaction(std::string line) {
  const int nfields = 13;
  Real Tmin, Tmax; 
  int r;

  //read the reaction
  reactants_ = StringUtils::split(line.substr(0, 33), ' ');
  products_ = StringUtils::split(line.substr(34, 55), ' ');
  std::vector<std::string> fields = StringUtils::split(line.substr(91, 85), ' ');

  if ( fields.size() != nfields ) {
    std::stringstream msg; 
    msg << "### FATAL ERROR in KidaReaction constructor [KidaReaction]" << std::endl
      << "number of fields (" << fields.size() << ")does not match " 
      << nfields << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  alpha_ = std::stof(fields[0]);
  beta_ = std::stof(fields[1]);
  gamma_ = std::stof(fields[2]);
  itype_ = std::stoi(fields[6]);
  Tmin = std::stof(fields[7]);
  Tmax = std::stof(fields[8]);
  formula_ = std::stoi(fields[9]);
  id_ = std::stoi(fields[10]);
  r = std::stoi(fields[12]);
  
  //Temperature range warning
  if ( (Tmin > 0.) || (Tmax < 9998.) ) {
    std::cout << "### WARNING KidaReaction constructor [KidaReaction]" << std::endl
      << "Temperature range (" << Tmin << "," << Tmax << ") for reaction (ID="  
      << id_ << "). Extrapolation outside of range will be used" << std::endl;
  }

  //recommendation error
  if (r == 0) {
    std::stringstream msg; 
    msg << "### FATAL ERROR in KidaReaction constructor [KidaReaction]" << std::endl
      << "reaction (ID=" << id_ << ") not recommended (r=0)." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

}

void KidaReaction::Print() const {
  std::cout << "reaction ID=" << id_ << ": ";
  for (int i=0; i<reactants_.size()-1; i++) {
    std::cout << reactants_[i] << " + ";
  }
  std::cout << reactants_[reactants_.size()-1]
    << " -> ";

  for (int i=0; i<products_.size()-1; i++) {
    std::cout << products_[i] << " + ";
  }
  std::cout << products_[products_.size()-1]
    << std::endl;

  return;
}
