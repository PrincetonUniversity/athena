#ifndef UTILS_UNITS_HPP_
#define UTILS_UNITS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file units.hpp
//! \brief prototypes of unit and constant classes

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

//! \brief Physical constants defined in c.g.s.
namespace Constants {
static const Real G     = 6.67259e-8;
static const Real Msun  = 1.9891e+33;
static const Real Lsun  = 3.8268e+33;
static const Real Myr   = 3.155815e+13;
static const Real pc    = 3.085678e+18;
static const Real kpc   = 3.085678e+21;
static const Real kms   = 1.0e+5;
static const Real mH    = 1.6733e-24;
static const Real aR    = 7.5646e-15;
static const Real kB    = 1.380658e-16;
static const Real c     = 2.99792458e+10;
static const Real e     = 4.80320427e-10;
} // namespace Constants

//! \brief Class for units
class Units {
 public:
  Units(Real dunit, Real lunit, Real vunit);
  Units(Real dunit, Real lunit, Real vunit, Real mu0);

  void SetUnitsConstants();
  void PrintCodeUnits();
  void PrintConstantsInCodeUnits();

  bool fixed_mu;

  Real Mass, Length, Time;
  Real Volume, Density, Velocity;
  Real EnergyDensity, Pressure;
  Real MagneticField;
  Real Temperature, mu;

  Real gram, cm, second, dyne, erg, Kelvin, Gauss;

  Real G_in_code;
  Real Msun_in_code;
  Real Lsun_in_code;
  Real Myr_in_code;
  Real pc_in_code;
  Real kpc_in_code;
  Real kms_in_code;
  Real mH_in_code;
  Real aR_in_code;
  Real kB_in_code;
  Real c_in_code;
  Real e_in_code;
  Real Bethe_in_code;
};
#endif // UTILS_UNITS_HPP_
