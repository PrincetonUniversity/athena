#ifndef CHEMISTRY_UTILS_KIDA_REACTION_HPP_
#define CHEMISTRY_UTILS_KIDA_REACTION_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file kida_reaction.hpp
//! \brief definitions for reactions in kida network format

//c++ header
#include <sstream>    // stringstream
#include <string>     // string
#include <vector>     // vector container

// Athena++ classes headers
#include "../../athena.hpp"

//! \class KidaReaction
//! \brief Chemical reactions in kida network format
class KidaReaction{
  friend class ChemNetwork;
 public:
  explicit KidaReaction(std::string line);
  void Print() const;

 private:
  std::vector<std::string> reactants_;
  std::vector<std::string> products_;

  int id_; //id number, has to be unique but doesn't have to be in order

  int itype_; //type of reaction
  int formula_; //type of formula

  //reaction rates coefficients;
  Real alpha_;
  Real beta_;
  Real gamma_;
  Real Tmin_;
  Real Tmax_;
};

#endif //CHEMISTRY_UTILS_KIDA_REACTION_HPP_
