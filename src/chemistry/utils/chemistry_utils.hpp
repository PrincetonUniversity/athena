#ifndef CHEMISTRY_UTILS_HPP
#define CHEMISTRY_UTILS_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file chemistry_utils.hpp
//  \brief prototypes of utility functions for Chang-Goo Kim's galactic disk
//  simulations.
//======================================================================================

#include <string>
#include "../../athena.hpp"

namespace ChemistryUtility
{
  const Real mumin=0.6182, mumax=1.295;
	const Real muH = 1.4271;
	//physical constants
	const Real kB = 1.380658e-16;
	const Real mH = 1.6733e-24; 
	const Real mCO = 4.68e-23;
  const Real pc = 3.085678e18; //parsec in cm
	//units
	const Real unitL = pc;
	const Real unitD = muH * mH;
	const Real unitV = 1.0e5;
	const Real unitT = unitD * unitV * unitV / (kB * muH);
	const Real unitE = unitD * unitV * unitV;
  //calculate tempereature in Kelvin from t1=m_H P / (rho * k_B)
  Real get_temp_from_t1(const Real t1);
  //calculate temperature in Kelvin
  //pressure: in units of [dens][v]^2, can be found in phydro.w(IEN, k, j, i). 
  //dens: density in units of [dens] = mu_H m_H cm^-3,
  //in phydro.w(IDN, k, j, i) and phydro.u(IDN, k, j, i)
  //veolcity units: [v] = km/s
  Real get_temp(Real pressure, Real dens);
	//find index of string
  int FindStrIndex(const std::string *str_arr, const int len,
		               const std::string name);
  int GetOppositeDirection(const int direction);
}
#endif // CHEMISTRY_UTILS_HPP
