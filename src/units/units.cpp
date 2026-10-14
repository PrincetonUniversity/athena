//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file units.cpp
//! \brief define unit class and physical constants

// C headers

// C++ headers
#include <iomanip>    // setprecision
#include <iostream>
#include <limits>     // numeric_limits
#include <sstream>    // stringstream
#include <stdexcept> // throw exceptions
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "units.hpp"

//========================================================================================
//! \fn void Units::Units(ParameterInput *pin)
//! \brief default unit constructor from the parameter input (default is for ISM problems)
//!        temperature units are not set (due to the mu dependence)
//!
//!        Additional basis units should be added to ToCgs().
//========================================================================================
Units::Units(ParameterInput *pin) :
  unit_system(pin->GetOrAddString("units", "unit_system", "ism")),
  mass_per_hydrogen(pin->GetOrAddReal("units", "mass_per_hydrogen", 1.4)) {
  // define units for given unit system
  // (note the slight change in the default mass_per_hydrogen 1.4271->1.4)
  // 1) ism unit system adopted in TIGRESS
  //     basis_length   = 1.0 pc
  //     basis_velocity = 1.0 km/s
  //     basis_ndensity = 1.0 n/cm^3
  //   which gives MLT units
  //     [mass]   = mass_per_hydrogen*m_h*(pc/cm)^3
  //     [length] = pc
  //     [time]   = (pc/km)*s
  // 2) galaxy unit system
  //     basis_length   = 1.0 kpc
  //     basis_time     = 1.0 Myr
  //     basis_ndensity = 1.0 n/cm^3
  //   which gives MLT units
  //     [mass]   = mass_per_hydrogen*m_h*(kpc/cm)^3
  //     [length] = kpc
  //     [time]   = Myr
  // 3) galaxypc unit system
  //     basis_length   = 1.0 pc
  //     basis_time     = 1.0 Myr
  //     basis_ndensity = 1.0 n/cm^3
  //   which gives MLT units
  //     [mass]   = mass_per_hydrogen*m_h*(pc/cm)^3
  //     [length] = pc
  //     [time]   = Myr
  // 4) ism_SR unit system (for special relativity)
  //     basis_length   = 1.0 pc
  //     basis_velocity = c
  //     basis_ndensity = 1.0 n/cm^3
  // 5) cgs unit system
  //     basis_length   = 1.0 cm
  //     basis_time     = 1.0 s
  //     basis_mass     = 1.0 g
  // 6) SI unit system
  //     basis_length   = 1.0 m
  //     basis_time     = 1.0 s
  //     basis_mass     = 1.0 kg
  // 7) custom_basis
  //    User chooses a 'length' basis (with possible units):
  //      length (pc, kpc, au, cm, m, km)
  //    Then either 'time' or 'velocity' (with possible units):
  //      time (yr, Myr, s)
  //      velocity (km/s, cm/s, m/s)
  //    Then either 'ndensity' or 'mass' (with possible units):
  //      ndensity (n/cm^3, n/m^3)
  //      mass (Msun, g, kg)
  // 8) custom
  //    User sets MLT units directly in c.g.s. with mass_cgs, length_cgs, time_cgs
  //
  // add other default units system here
  // The basis values that are not part of the chosen basis are placeholders with
  // default units; CompleteBasis() overwrites their values.
  if (unit_system.compare("ism") == 0) {
    basis_length   = std::make_tuple(1.0, "pc");
    basis_velocity = std::make_tuple(1.0, "km/s");
    basis_ndensity = std::make_tuple(1.0, "n/cm^3");
    basis_time     = std::make_tuple(1.0, "Myr");
    basis_mass     = std::make_tuple(1.0, "Msun");
    velocity_basis_ = true;
  } else if (unit_system.compare("galaxy") == 0) {
    basis_length   = std::make_tuple(1.0, "kpc");
    basis_time     = std::make_tuple(1.0, "Myr");
    basis_ndensity = std::make_tuple(1.0, "n/cm^3");
    basis_velocity = std::make_tuple(1.0, "km/s");
    basis_mass     = std::make_tuple(1.0, "Msun");
  } else if (unit_system.compare("galaxypc") == 0) {
    basis_length   = std::make_tuple(1.0, "pc");
    basis_time     = std::make_tuple(1.0, "Myr");
    basis_ndensity = std::make_tuple(1.0, "n/cm^3");
    basis_velocity = std::make_tuple(1.0, "km/s");
    basis_mass     = std::make_tuple(1.0, "Msun");
  } else if (unit_system.compare("ism_SR") == 0) {
    basis_length   = std::make_tuple(1.0, "pc");
    basis_velocity = std::make_tuple(Constants::speed_of_light_cgs/Constants::km_s_cgs,
                                     "km/s");
    basis_ndensity = std::make_tuple(1.0, "n/cm^3");
    basis_time     = std::make_tuple(1.0, "Myr");
    basis_mass     = std::make_tuple(1.0, "Msun");
    velocity_basis_ = true;
  } else if (unit_system.compare("cgs") == 0) {
    basis_length   = std::make_tuple(1.0, "cm");
    basis_time     = std::make_tuple(1.0, "s");
    basis_mass     = std::make_tuple(1.0, "g");
    basis_velocity = std::make_tuple(1.0, "cm/s");
    basis_ndensity = std::make_tuple(1.0, "n/cm^3");
    mass_basis_ = true;
  } else if (unit_system.compare("SI") == 0) {
    basis_length   = std::make_tuple(1.0, "m");
    basis_time     = std::make_tuple(1.0, "s");
    basis_mass     = std::make_tuple(1.0, "kg");
    basis_velocity = std::make_tuple(1.0, "m/s");
    basis_ndensity = std::make_tuple(1.0, "n/m^3");
    mass_basis_ = true;
  } else if (unit_system.compare("custom_basis") == 0) {
    basis_length   = std::make_tuple(1.0, pin->GetOrAddString("units", "length_unit",
                                                              "pc"));
    basis_time     = std::make_tuple(1.0, pin->GetOrAddString("units", "time_unit",
                                                              "Myr"));
    basis_velocity = std::make_tuple(1.0, pin->GetOrAddString("units", "velocity_unit",
                                                              "km/s"));
    basis_ndensity = std::make_tuple(1.0, pin->GetOrAddString("units", "ndensity_unit",
                                                              "n/cm^3"));
    basis_mass     = std::make_tuple(1.0, pin->GetOrAddString("units", "mass_unit",
                                                              "Msun"));

    bool has_time = pin->DoesParameterExist("units", "time");
    bool has_velocity = pin->DoesParameterExist("units", "velocity");
    bool has_ndensity = pin->DoesParameterExist("units", "ndensity");
    bool has_mass = pin->DoesParameterExist("units", "mass");
    if (!pin->DoesParameterExist("units", "length") || has_time == has_velocity
        || has_ndensity == has_mass) {
      std::stringstream msg;
      msg << "### FATAL ERROR in Units constructor" << std::endl
          << "  unit_system=custom_basis requires length, exactly one of" << std::endl
          << "  time or velocity, and exactly one of ndensity or mass" << std::endl;
      ATHENA_ERROR(msg);
    }

    std::get<0>(basis_length) = pin->GetReal("units", "length");
    if (has_time) {
      std::get<0>(basis_time) = pin->GetReal("units", "time");
    } else {
      std::get<0>(basis_velocity) = pin->GetReal("units", "velocity");
      velocity_basis_ = true;
    }
    if (has_ndensity) {
      std::get<0>(basis_ndensity) = pin->GetReal("units", "ndensity");
    } else {
      std::get<0>(basis_mass) = pin->GetReal("units", "mass");
      mass_basis_ = true;
    }
  } else if (unit_system.compare("custom") == 0) {
    // this must raise error if MLT units are not given in the input file
    code_mass_cgs_ = pin->GetReal("units", "mass_cgs");
    code_length_cgs_ = pin->GetReal("units", "length_cgs");
    code_time_cgs_ = pin->GetReal("units", "time_cgs");

    // express the MLT units in the default basis units
    basis_length   = std::make_tuple(code_length_cgs_/Constants::pc_cgs, "pc");
    basis_time     = std::make_tuple(code_time_cgs_/Constants::million_yr_cgs, "Myr");
    basis_velocity = std::make_tuple(
        code_length_cgs_/code_time_cgs_/Constants::km_s_cgs, "km/s");
    basis_ndensity = std::make_tuple(
        code_mass_cgs_/(mass_per_hydrogen*Constants::hydrogen_mass_cgs
                        *CUBE(code_length_cgs_)), "n/cm^3");
    basis_mass     = std::make_tuple(code_mass_cgs_/Constants::solar_mass_cgs, "Msun");
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in Units constructor" << std::endl
        << "  unit_system=" << unit_system << " is not valid unit system " << std::endl
        << "  choose one of the default unit systems in" << std::endl
        << "  [ism, galaxy, galaxypc, ism_SR, cgs, SI], or set" << std::endl
        << "  unit_system=custom_basis or custom and define the units manually"
        << std::endl;
    ATHENA_ERROR(msg);
  }

  // compute MLT units from the basis and write them back to the input file
  if (unit_system.compare("custom") != 0) {
    CompleteBasis();
    pin->SetReal("units","mass_cgs",code_mass_cgs_);
    pin->SetReal("units","length_cgs",code_length_cgs_);
    pin->SetReal("units","time_cgs",code_time_cgs_);
  }

  // calculate default unit conversion factors
  SetUnitsConstants();
}

//========================================================================================
//! \fn void Units::SetUnitsConstants()
//! \brief calculate default unit conversion factors, constants in code units
//========================================================================================
void Units::SetUnitsConstants() {
  // Multiply a code value by code_X_cgs to get the value in c.g.s.,
  // e.g. mass (code units) * code_mass_cgs = mass (g).
  // Multiply a code value by std::get<0>(basis_X) to get the value in the basis unit,
  // e.g. mass (code units) * std::get<0>(basis_mass) = mass (Msun).
  // Divide a code value by X_code to get the value in units of X,
  // e.g. mass (code units) / solar_mass_code = mass (Msun).

  // set public MLT unit variable
  code_mass_cgs = code_mass_cgs_;
  code_length_cgs = code_length_cgs_;
  code_time_cgs = code_time_cgs_;

  // variable units in cgs
  code_volume_cgs = CUBE(code_length_cgs_);
  code_density_cgs = code_mass_cgs_/code_volume_cgs;
  code_velocity_cgs = code_length_cgs_/code_time_cgs_;

  code_energy_cgs = code_mass_cgs_*SQR(code_velocity_cgs);
  code_energydensity_cgs = code_pressure_cgs = code_density_cgs*SQR(code_velocity_cgs);

  code_magneticfield_cgs = std::sqrt(4.*PI*code_pressure_cgs);

  code_temperature_mu_cgs = code_pressure_cgs/code_density_cgs
                           *Constants::hydrogen_mass_cgs/Constants::k_boltzmann_cgs;

  // constans in code units
  cm_code = 1.0/code_length_cgs_;
  gram_code = 1.0/code_mass_cgs_;
  second_code = 1.0/code_time_cgs_;
  dyne_code = gram_code*cm_code/(second_code*second_code);
  erg_code = dyne_code*cm_code;
  kelvin_code = 1.0; // (changgoo) in principle, this should be 1/[code temperature]
                     // but this is what has been adopted in Athena (not sure why)

  grav_const_code = Constants::grav_const_cgs
                     *cm_code*cm_code*cm_code/(gram_code*second_code*second_code);
  solar_mass_code = Constants::solar_mass_cgs*gram_code;
  solar_lum_code = Constants::solar_lum_cgs*erg_code/second_code;
  yr_code = Constants::yr_cgs*second_code;
  million_yr_code = Constants::million_yr_cgs*second_code;
  pc_code = Constants::pc_cgs*cm_code;
  kpc_code = Constants::kpc_cgs*cm_code;
  km_s_code = Constants::km_s_cgs*cm_code/second_code;
  hydrogen_mass_code = Constants::hydrogen_mass_cgs*gram_code;
  radiation_aconst_code = Constants::radiation_aconst_cgs*erg_code
                         /(cm_code*cm_code*cm_code
                          *kelvin_code*kelvin_code*kelvin_code*kelvin_code);
  k_boltzmann_code = Constants::k_boltzmann_cgs*erg_code/kelvin_code;
  speed_of_light_code = Constants::speed_of_light_cgs*cm_code/second_code;
  echarge_code = Constants::echarge_cgs*std::sqrt(dyne_code*4*PI)*cm_code;
  bethe_code = 1.e51 * erg_code;
}

//========================================================================================
//! \fn void Units::PrintBasisUnits()
//! \brief print basis parameters in chosen units
//========================================================================================
void Units::PrintBasisUnits() {
  std::stringstream ss;
  ss << std::scientific << "============ Unit Basis ============" << std::endl;
  ss << "basis length = " << std::get<0>(basis_length) << " "
     << std::get<1>(basis_length) << std::endl;
  ss << "basis time = " << std::get<0>(basis_time) << " "
     << std::get<1>(basis_time) << std::endl;
  ss << "basis mass = " << std::get<0>(basis_mass) << " "
     << std::get<1>(basis_mass) << std::endl;
  ss << "basis velocity = " << std::get<0>(basis_velocity) << " "
     << std::get<1>(basis_velocity) << std::endl;
  ss << "basis ndensity = " << std::get<0>(basis_ndensity) << " "
     << std::get<1>(basis_ndensity) << std::endl;
  ss << "====================================" << std::endl;
  std::cout << ss.str();
}

//========================================================================================
//! \fn void Units::PrintCodeUnits()
//! \brief print code units in c.g.s.
//========================================================================================
void Units::PrintCodeUnits() {
  std::stringstream ss;
  ss << std::scientific << "============ Code Units ============" << std::endl;
  ss << "code Mass = " << code_mass_cgs << " g" << std::endl;
  ss << "code Length = " << code_length_cgs << " cm" << std::endl;
  ss << "code Time = " << code_time_cgs << " s" << std::endl;

  ss << "code density = " << code_density_cgs << " g/cm^3" << std::endl;
  ss << "code velocity = " << code_velocity_cgs << " cm/s" << std::endl;
  ss << "code energy = " << code_energy_cgs << " erg" << std::endl;
  ss << "code pressure = " << code_pressure_cgs << " erg/cm^3" << std::endl;
  ss << "code temperature/mu = " << code_temperature_mu_cgs << " K" << std::endl;
  ss << "====================================" << std::endl;
  std::cout << ss.str();
}

//========================================================================================
//! \fn void Units::PrintConstantsInCodeUnits()
//! \brief print physical constatns in code units
//========================================================================================
void Units::PrintConstantsInCodeUnits() {
  std::stringstream ss;
  ss << std::scientific << "===== Constants  in Code Units =====" << std::endl;
  ss << "dyne in code = " << dyne_code << std::endl;
  ss << "erg in code = " << erg_code << std::endl;

  ss << "Gconst in code = " << grav_const_code << std::endl;
  ss << "Msun in code = " << solar_mass_code << std::endl;
  ss << "Lsun in code = " << solar_lum_code << std::endl;
  ss << "Myr in code = " << million_yr_code << std::endl;
  ss << "pc in code = " << pc_code << std::endl;
  ss << "km/s in code = " << km_s_code << std::endl;
  ss << "m_H in code = " << hydrogen_mass_code << std::endl;
  ss << "kB in code = " << k_boltzmann_code << std::endl;
  ss << "c in code = " << speed_of_light_code << std::endl;
  ss << "e in code = " << echarge_code << std::endl;
  ss << "====================================" << std::endl;
  std::cout << ss.str();
}

//========================================================================================
//! \fn void Units::PrintYtUnitsOverride()
//! \brief print the units_override dictionary for loading HDF5 outputs into yt
//========================================================================================
void Units::PrintYtUnitsOverride() {
  std::stringstream ss;
  ss << std::scientific << std::setprecision(std::numeric_limits<Real>::max_digits10)
     << "{\"length_unit\":(" << std::get<0>(basis_length) << ",\""
     << std::get<1>(basis_length)
     << "\"),\"time_unit\":(" << std::get<0>(basis_time) << ",\""
     << std::get<1>(basis_time)
     << "\"),\"mass_unit\":(" << std::get<0>(basis_mass) << ",\""
     << std::get<1>(basis_mass) << "\")}";
  std::cout << "yt units_override for loading HDF5 outputs:" << std::endl
            << ss.str() << std::endl;
}

//========================================================================================
//! \fn Real Units::ToCgs(const std::string &parameter, Real value,
//!                       const std::string &unit)
//! \brief Converts a basis value given in unit into c.g.s.
//!
//!        Additional units should be added here.
//========================================================================================
Real Units::ToCgs(const std::string &parameter, Real value, const std::string &unit) {
  std::string allowed;
  if (parameter == "length") {
    if (unit == "pc") return Constants::pc_cgs*value;
    if (unit == "kpc") return Constants::kpc_cgs*value;
    if (unit == "au") return Constants::au_cgs*value;
    if (unit == "cm") return value;
    if (unit == "m") return 1.0e2*value;
    if (unit == "km") return 1.0e5*value;
    allowed = "pc, kpc, au, cm, m, km";
  } else if (parameter == "time") {
    if (unit == "yr") return Constants::yr_cgs*value;
    if (unit == "Myr") return Constants::million_yr_cgs*value;
    if (unit == "s") return value;
    allowed = "yr, Myr, s";
  } else if (parameter == "velocity") {
    if (unit == "km/s") return Constants::km_s_cgs*value;
    if (unit == "m/s") return 1.0e2*value;
    if (unit == "cm/s") return value;
    allowed = "km/s, m/s, cm/s";
  } else if (parameter == "ndensity") {
    if (unit == "n/cm^3") return value;
    if (unit == "n/m^3") return 1.0e-6*value;
    allowed = "n/cm^3, n/m^3";
  } else if (parameter == "mass") {
    if (unit == "Msun") return Constants::solar_mass_cgs*value;
    if (unit == "g") return value;
    if (unit == "kg") return 1.0e3*value;
    allowed = "Msun, g, kg";
  }
  std::stringstream msg;
  msg << "### FATAL ERROR in Units::ToCgs" << std::endl
      << "  " << parameter << " unit " << unit << " is not a valid unit." << std::endl
      << "  Allowed units are: " << allowed << std::endl;
  ATHENA_ERROR(msg);
  return 0.0;
}

//========================================================================================
//! \fn void Units::CompleteBasis()
//! \brief Sets the code MLT units in c.g.s. from the chosen basis, and fills in the
//!        values of the basis quantities that were not part of the chosen basis
//========================================================================================
void Units::CompleteBasis() {
  code_length_cgs_ = ToCgs("length", std::get<0>(basis_length),
                           std::get<1>(basis_length));

  // time from either velocity or time basis
  if (velocity_basis_) {
    Real velocity_cgs = ToCgs("velocity", std::get<0>(basis_velocity),
                                   std::get<1>(basis_velocity));
    code_time_cgs_ = code_length_cgs_/velocity_cgs;
    std::get<0>(basis_time) = code_time_cgs_/ToCgs("time", 1.0, std::get<1>(basis_time));
  } else {
    code_time_cgs_ = ToCgs("time", std::get<0>(basis_time), std::get<1>(basis_time));
    Real velocity_cgs = code_length_cgs_/code_time_cgs_;
    std::get<0>(basis_velocity) =
        velocity_cgs/ToCgs("velocity", 1.0, std::get<1>(basis_velocity));
  }

  // mass from either mass or ndensity basis
  Real mass_per_particle = mass_per_hydrogen*Constants::hydrogen_mass_cgs;
  if (mass_basis_) {
    code_mass_cgs_ = ToCgs("mass", std::get<0>(basis_mass), std::get<1>(basis_mass));
    Real ndensity_cgs = code_mass_cgs_/(mass_per_particle*CUBE(code_length_cgs_));
    std::get<0>(basis_ndensity) =
        ndensity_cgs/ToCgs("ndensity", 1.0, std::get<1>(basis_ndensity));
  } else {
    Real ndensity_cgs = ToCgs("ndensity", std::get<0>(basis_ndensity),
                                   std::get<1>(basis_ndensity));
    code_mass_cgs_ = mass_per_particle*CUBE(code_length_cgs_)*ndensity_cgs;
    std::get<0>(basis_mass) = code_mass_cgs_/ToCgs("mass", 1.0, std::get<1>(basis_mass));
  }
}
