#ifndef CHEMISTRY_UTILS_SHIELDING_HPP_
#define CHEMISTRY_UTILS_SHIELDING_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file shielding.hpp
//! \brief definition of shielding class.

//c header
#include <math.h> //a^x = pow(a,x)

//athena++ header
#include "../../athena.hpp"

//! \class Shielding
//! \brief Shielding functions.
class Shielding {
 public:
  Shielding();
  static Real fShield_CO_vDB(const Real NCO, const Real NH2);
  static Real fShield_CO_V09(const Real NCO, const Real NH2);
  static Real fShield_H2(const Real NH2, const Real bH2);
  static Real fShield_C(const Real NC, const Real NH2);
  static Real fShield_CO_C(const Real NC);

 private:
  //CO column density for DB table
  static const int len_NCO_DB_ = 8;
  static const Real logNCOvDB_[len_NCO_DB_];
  // H2 column densities for DB table
  static const int len_NH2_DB_ = 6;
  static const Real logNH2vDB_[len_NH2_DB_];
  // Tabulated shielding factors
  static const Real ThetavDB_[len_NH2_DB_][len_NCO_DB_];

  //Visser+ 2009 Table 5, b(CO)=0.3, Tex(CO)=5K
  static const int len_NCO_V09_ = 47;
  static const Real logNCOV09_[len_NCO_V09_];
  static const int len_NH2_V09_ = 42;
  static const Real logNH2V09_[len_NH2_V09_];
  static const Real ThetaV09_[len_NH2_V09_][len_NCO_V09_];
};

#endif //CHEMISTRY_UTILS_SHIELDING_HPP_
