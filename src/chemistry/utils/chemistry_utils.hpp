#ifndef CHEMISTRY_UTILS_CHEMISTRY_UTILS_HPP_
#define CHEMISTRY_UTILS_CHEMISTRY_UTILS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file chemistry_utils.hpp
//! \brief utility functions for TIGRESS simulations by C.G. Kim

//c++ header
#include <string>

//athena++ header
#include "../../athena.hpp"

//! \namespace ChemistryUtility
//! \brief Utility functions for chemistry
//TODO (Munan Gong) Perhaps change this to a class
namespace ChemistryUtility {
  const Real mumin=0.6182, mumax=1.295;
  const Real muH = 1.4271;
  //physical constants
  const Real kB = 1.380658e-16;
  const Real mH = 1.6733e-24;
  const Real me = 9.1093897e-28;//electron mass
  const Real mCO = 4.68e-23;
  const Real pc = 3.085678e18; //parsec in cm
  const Real qe = 4.803206e-10; //elementary charge
  //units
  const Real unitL = pc;
  const Real unitD = muH * mH;
  const Real unitV = 1.0e5;
  const Real unitT = unitD * unitV * unitV / (kB * muH);
  const Real unitE = unitD * unitV * unitV;
  Real get_temp_from_t1(const Real t1);
  Real get_temp(Real pressure, Real dens);
  int FindStrIndex(const std::string *str_arr, const int len,
                   const std::string name);
  int GetOppositeDirection(const int direction);
} // namespace ChemistryUtility
#endif // CHEMISTRY_UTILS_CHEMISTRY_UTILS_HPP_
