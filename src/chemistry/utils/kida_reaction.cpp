//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file kida_reaction.cpp
//! \brief implementation of functions in class KidaReaction

// this class header
#include "kida_reaction.hpp"

// C headers

// C++ headers
#include <iostream>   // endl
#include <stdexcept>  // std::runtime_error()
#include <string>

// Athena++ headers
#include "../../utils/string_utils.hpp"

//----------------------------------------------------------------------------------------
//! constructor
KidaReaction::KidaReaction(std::string line) {
  const int nfields = 13;
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
  Tmin_ = std::stof(fields[7]);
  Tmax_ = std::stof(fields[8]);
  formula_ = std::stoi(fields[9]);
  id_ = std::stoi(fields[10]);
  r = std::stoi(fields[12]);

  //Temperature range warning
  //if ( (Tmin_ > 0.) || (Tmax_ < 9998.) ) {
  //  std::cout << "### WARNING KidaReaction constructor [KidaReaction]" << std::endl
  //    << "Temperature range (" << Tmin_ << "," << Tmax_ << ") for reaction (ID="
  //    << id_ << "). Extrapolation or temperature cap will be used" << std::endl;
  //}

  //recommendation error
  if (r == 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in KidaReaction constructor [KidaReaction]" << std::endl
      << "reaction (ID=" << id_ << ") not recommended (r=0)." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
}

//----------------------------------------------------------------------------------------
//! \fn void KidaReaction::Print() const
//! \brief Print reaction to standard output
void KidaReaction::Print() const {
  std::cout << "reaction ID=" << id_ << ": ";
  for (std::size_t i=0; i<reactants_.size()-1; i++) {
    std::cout << reactants_[i] << " + ";
  }
  std::cout << reactants_[reactants_.size()-1]
    << " -> ";

  for (std::size_t i=0; i<products_.size()-1; i++) {
    std::cout << products_[i] << " + ";
  }
  std::cout << products_[products_.size()-1]
    << std::endl;

  return;
}
