//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file units.cpp
//! \brief define unit class and physical constants

// C headers

// C++ headers
#include <iostream>
#include <sstream>    // stringstream
#include <stdexcept> //throw exceptions

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "units.hpp"

//========================================================================================
//! \fn void Units::Units(ParameterInput *pin)
//! \brief default unit constructor from the parameter input (default is for ISM problems)
//!        temperature units are not set (due to the mu dependence)
//========================================================================================
Units::Units(ParameterInput *pin) :
  unit_system(pin->GetOrAddString("units", "unit_system", "ism")) {
  // define units for given unit system
  // 1) ism unit system adopted in TIGRESS
  //   (note the slight change in the density units 1.4271->1.4):
  //     [rho] = 1.4 m_h/cm^3
  //     [length] = pc
  //     [velocity] = km/s
  //   which gives MLT units
  //     [mass] = 1.4 m_h*(pc/cm)^3
  //     [length] = pc
  //     [time] = (pc/km)*s
  // 2) galaxy unit system (example)
  //     [mass] = 1.4 m_h*(pc/cm)^3
  //     [length] = kpc
  //     [time] = Myr
  // add other default units system here
  // set private MLT unit variables here
  if (unit_system.compare("ism") == 0) {
    code_length_cgs_ = Constants::pc_cgs;
    code_mass_cgs_ = 1.4*Constants::hydrogen_mass_cgs*CUBE(code_length_cgs_);
    code_time_cgs_ = Constants::pc_cgs/Constants::km_s_cgs;
  } else if (unit_system.compare("galaxy") == 0) {
    code_length_cgs_ = Constants::kpc_cgs;
    code_mass_cgs_ = 1.4*Constants::hydrogen_mass_cgs*CUBE(code_length_cgs_);
    code_time_cgs_ = Constants::million_yr_cgs;
  } else if (unit_system.compare("custom") == 0) {
    // this must raise error if MLT units are not given in the input file
    code_mass_cgs_ = pin->GetReal("units", "mass_cgs");
    code_length_cgs_ = pin->GetReal("units", "length_cgs");
    code_time_cgs_ = pin->GetReal("units", "time_cgs");
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in Units constructor" << std::endl
        << "  unit_system=" << unit_system << " is not valid unit system " << std::endl
        << "  choose one of the default unit systems in [ism] or " << std::endl
        << "  set unit_system=custom and define MLT units manually" << std::endl;
    ATHENA_ERROR(msg);
  }

  // write units back to the input file
  if (unit_system.compare("custom") != 0) {
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
  // set public MLT unit variable
  code_mass_cgs = code_mass_cgs_;
  code_length_cgs = code_length_cgs_;
  code_time_cgs = code_time_cgs_;

  // variable units in cgs
  code_volume_cgs = CUBE(code_length_cgs_);
  code_density_cgs = code_mass_cgs_/code_volume_cgs;
  code_velocity_cgs = code_length_cgs_/code_time_cgs_;

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
//! \fn void Units::PrintCodeUnits()
//! \brief print code units in c.g.s.
//========================================================================================
void Units::PrintCodeUnits() {
  std::cout << std::scientific << "============ Code Units ============" << std::endl;
  std::cout << "code Mass = " << code_mass_cgs << " g" << std::endl;
  std::cout << "code Length = " << code_length_cgs << " cm" << std::endl;
  std::cout << "code Time = " << code_time_cgs << " s" << std::endl;

  std::cout << "code density = " << code_density_cgs << " g/cm^3" << std::endl;
  std::cout << "code velocity = " << code_velocity_cgs << " cm/s" << std::endl;
  std::cout << "code pressure = " << code_pressure_cgs << " erg/cm^3" << std::endl;
  std::cout << "====================================" << std::endl;
}

//========================================================================================
//! \fn void Units::PrintConstantsInCodeUnits()
//! \brief print physical constatns in code units
//========================================================================================
void Units::PrintConstantsInCodeUnits() {
  std::cout << std::scientific << "===== Constants  in Code Units =====" << std::endl;
  std::cout << "dyne in code = " << dyne_code << std::endl;
  std::cout << "erg in code = " << erg_code << std::endl;

  std::cout << "Gconst in code = " << grav_const_code << std::endl;
  std::cout << "Msun in code = " << solar_mass_code << std::endl;
  std::cout << "Lsun in code = " << solar_lum_code << std::endl;
  std::cout << "Myr in code = " << million_yr_code << std::endl;
  std::cout << "kB in code = " << k_boltzmann_code << std::endl;
  std::cout << "c in code = " << speed_of_light_code << std::endl;
  std::cout << "e in code = " << echarge_code << std::endl;
  std::cout << "====================================" << std::endl;
}
