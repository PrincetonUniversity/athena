#ifndef CHEMISTRY_UTILS_CHEMISTRY_UTILS_HPP_
#define CHEMISTRY_UTILS_CHEMISTRY_UTILS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file chemistry_utils.hpp
//! \brief utility functions and constants for chemistry

// C headders

// C++ header
#include <string>

// Athena++ headers
#include "../../athena.hpp"

//! \namespace ChemistryUtility
//! \brief Utility functions for chemistry
namespace ChemistryUtility {
  // physical constants
  const Real me = 9.1093897e-28; // electron mass
  const Real mCO = 4.68e-23;
  int FindStrIndex(const std::string *str_arr, const int len,
                   const std::string name);
} // namespace ChemistryUtility
#endif // CHEMISTRY_UTILS_CHEMISTRY_UTILS_HPP_
