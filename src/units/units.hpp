#ifndef UNITS_UNITS_HPP_
#define UNITS_UNITS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file units.hpp
//! \brief prototypes of unit and constant classes

// C headers

// C++ headers
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

//! \brief Physical constants defined in c.g.s.
namespace Constants {
static const Real grav_const_cgs       = 6.67259e-8;
static const Real solar_mass_cgs       = 1.9891e+33;
static const Real solar_lum_cgs        = 3.8268e+33;
static const Real yr_cgs               = 3.155815e+7;
static const Real million_yr_cgs       = 3.155815e+13;
static const Real pc_cgs               = 3.085678e+18;
static const Real kpc_cgs              = 3.085678e+21;
static const Real km_s_cgs             = 1.0e+5;
static const Real hydrogen_mass_cgs    = 1.6733e-24;
static const Real radiation_aconst_cgs = 7.5646e-15;
static const Real k_boltzmann_cgs      = 1.380658e-16;
static const Real speed_of_light_cgs   = 2.99792458e+10;
static const Real echarge_cgs          = 4.80320427e-10;
static const Real kelvin_cgs           = 1;
} // namespace Constants

//! \brief Class for units
class Units {
 public:
  explicit Units(ParameterInput *pin);

  void SetUnitsConstants();
  void PrintCodeUnits();
  void PrintConstantsInCodeUnits();

  // unit system
  std::string unit_system;

  // code units in c.g.s.
  // i.e. multiply this to convert quantities in c.g.s.
  // public accessor for default MLT units
  // ideally this should be accessor methods, but
  // no on-the-fly change in unit system is assumed
  Real code_mass_cgs, code_length_cgs, code_time_cgs;

  Real code_volume_cgs, code_density_cgs, code_velocity_cgs;
  Real code_energydensity_cgs, code_pressure_cgs;
  Real code_magneticfield_cgs;
  Real code_temperature_mu_cgs; // T/mu

  // c.g.s. units in code units
  Real gram_code, cm_code, second_code, dyne_code, erg_code, kelvin_code;

  // physical constants in code units
  Real grav_const_code;
  Real solar_mass_code;
  Real solar_lum_code;
  Real yr_code;
  Real million_yr_code;
  Real pc_code;
  Real kpc_code;
  Real km_s_code;
  Real hydrogen_mass_code;
  Real radiation_aconst_code; // aR
  Real k_boltzmann_code; // k_B
  Real speed_of_light_code;
  Real echarge_code;
  Real bethe_code; // 1.e51 erg

 private:
  // code MLT units in c.g.s.
  Real code_mass_cgs_, code_length_cgs_, code_time_cgs_;
};
#endif // UNITS_UNITS_HPP_
