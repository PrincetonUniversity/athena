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
#include <stdexcept> // throw exceptions

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "units.hpp"

//========================================================================================
//! \fn void Units::Units(ParameterInput *pin)
//! \brief default unit constructor from the parameter input (default is for ISM problems)
//!        temperature units are not set (due to the mu (mean_weight) dependence)
//!
//!        Additional units should be added to Returncgs().
//========================================================================================
Units::Units(ParameterInput *pin) :
  unit_system(pin->GetOrAddString("units", "unit_system", "ism")) {
  //   (note the slight change in the density units 1.4271->1.4):
  mean_weight = pin->GetOrAddReal("units", "mean_weight", 1.4);
  // define units for given unit system
  // 1) ism unit system adopted in TIGRESS
  //     basis_length   = 1.0 pc
  //     basis_velocity = 1.0 km/s
  //     basis_ndensity = 1.0 n/cm^3
  //   which gives MLT units
  //     [mass]   = mean_weight*m_h*(pc/cm)^3
  //     [length] = pc
  //     [time]   = (pc/km)*s
  // 2) galaxy unit system
  //     basis_length   = 1.0 kpc
  //     basis_time     = 1.0 Myr
  //     basis_ndensity = 1.0 n/cm^3
  //   which gives MLT units
  //     [mass]   = n*mean_weight*m_h*(pc/cm)^3
  //     [length] = kpc
  //     [time]   = Myr
  // 3) galaxypc unit system
  //     basis_length   = 1.0 pc
  //     basis_time     = 1.0 Myr
  //     basis_ndensity = 1.0 n/cm^3
  //   which gives MLT units
  //     [mass]   = n*mean_weight*m_h*(pc/cm)^3
  //     [length] = pc
  //     [time]   = Myr
  // 4) ism_SR unit system (for special relativity)
  //     basis_length   = 1.0 pc
  //     basis_velocity = 2.99792458e5 km/s (speed of light)
  //     basis_ndensity = 1.0 n/cm^3
  // 5) cgs unit system
  //     basis_length   = 1.0 cm
  //     basis_mass     = 1.0 g
  //     basis_time     = 1.0 s
  //     basis_velocity = 1.0 cm/s
  //     basis_ndensity = 1.0 n/cm^3
  // 6) SI unit system
  //     basis_length   = 1.0 m
  //     basis_mass     = 1.0 kg
  //     basis_time     = 1.0 s
  //     basis_velocity = 1.0 m/s
  //     basis_ndensity = 1.0 n/m^3
  // 7) custom_basis
  //    User chooses a 'length' basis (with possible units):
  //      length (pc, kpc, au, cm, m, km)
  //    Then either 'time' or 'velocity' (with possible units):
  //      time (yr, Myr, s)
  //      velocity (km/s, cm/s, m/s)
  //    Then either 'ndensity' or 'mass' (with possible units):
  //      ndensity (n/cm^3, n/m^3)
  //      mass (Msun, g, kg)
  //
  // add other default units system here
  // set private MLT unit variables here
  if (unit_system.compare("ism") == 0) {
    basis_length   = std::make_tuple(1.0,"pc");
    basis_velocity = std::make_tuple(1.0,"km/s");
    basis_ndensity = std::make_tuple(1.0,"n/cm^3");

    // Complete the basis set with dummy values and default units
    basis_time = std::make_tuple(1.0,"Myr");
    basis_mass = std::make_tuple(1.0,"Msun");

    velocity_basis = true;

  } else if (unit_system.compare("galaxy") == 0) {
    basis_length   = std::make_tuple(1.0,"kpc");
    basis_time     = std::make_tuple(1.0,"Myr");
    basis_ndensity = std::make_tuple(1.0,"n/cm^3");

    // Complete the basis set with dummy values and default units
    basis_velocity = std::make_tuple(1.0,"km/s");
    basis_mass = std::make_tuple(1.0,"Msun");

  } else if (unit_system.compare("galaxypc") == 0) {
    basis_length   = std::make_tuple(1.0,"pc");
    basis_time     = std::make_tuple(1.0,"Myr");
    basis_ndensity = std::make_tuple(1.0,"n/cm^3");

    // Complete the basis set with dummy values and default units
    basis_velocity = std::make_tuple(1.0,"km/s");
    basis_mass = std::make_tuple(1.0,"Msun");

  } else if (unit_system.compare("ism_SR") == 0) {
    // ism unit system, but for special relativity
    basis_length   = std::make_tuple(1.0,"pc");
    basis_velocity = std::make_tuple(2.99792458e5,"km/s");
    basis_ndensity = std::make_tuple(1.0,"n/cm^3");

    // Complete the basis set with dummy values and default units
    basis_time = std::make_tuple(1.0,"Myr");
    basis_mass = std::make_tuple(1.0,"Msun");

    velocity_basis = true;
  } else if (unit_system.compare("cgs") == 0) {
    // ism unit system, but for special relativity
    basis_length = std::make_tuple(1.0,"cm");
    basis_time   = std::make_tuple(1.0,"s");
    basis_mass   = std::make_tuple(1.0,"g");

    // Complete the basis set with dummy values and default units
    basis_velocity = std::make_tuple(1.0,"cm/s");
    basis_ndensity = std::make_tuple(1.0,"n/cm^3");

    mass_basis = true;
  } else if (unit_system.compare("SI") == 0) {
    // ism unit system, but for special relativity
    basis_length = std::make_tuple(1.0,"m");
    basis_time   = std::make_tuple(1.0,"s");
    basis_mass   = std::make_tuple(1.0,"kg");

    // Complete the basis set with dummy values and default units
    basis_velocity = std::make_tuple(1.0,"m/s");
    basis_ndensity = std::make_tuple(1.0,"n/m^3");

    mass_basis = true;
  } else if (unit_system.compare("custom_basis") == 0) {
    // New unit basis systems should be added here
    int basis_count = 0;

    // Create the tuples with dummy values and default units if not given
    basis_length   = std::make_tuple(1.0,pin->GetOrAddString("units",
                                                             "length_unit","pc"));
    basis_time     = std::make_tuple(1.0,pin->GetOrAddString("units",
                                                             "time_unit","Myr"));
    basis_velocity = std::make_tuple(1.0,pin->GetOrAddString("units",
                                                             "velocity_unit","km/s"));
    basis_ndensity = std::make_tuple(1.0,pin->GetOrAddString("units",
                                                             "ndensity_unit","n/cm^3"));
    basis_mass     = std::make_tuple(1.0,pin->GetOrAddString("units",
                                                             "mass_unit","Msun"));

    // Process values from input file
    if (pin->DoesParameterExist("units", "length")) {
      std::get<0>(basis_length) = pin->GetReal("units", "length");
      basis_count++;
    }

    // Either 'time' or 'velocity' but not both
    if (pin->DoesParameterExist("units", "time")) {
      std::get<0>(basis_time) = pin->GetReal("units", "time");
      basis_count++;
    } else if (pin->DoesParameterExist("units", "velocity")) {
      std::get<0>(basis_velocity) = pin->GetReal("units", "velocity");
      velocity_basis = true;
      basis_count++;
    }

    // Either 'ndensity' or 'mass' but not both
    if (pin->DoesParameterExist("units", "ndensity")) {
      std::get<0>(basis_ndensity) = pin->GetReal("units", "ndensity");
      basis_count++;
    } else if (pin->DoesParameterExist("units", "mass")) {
      std::get<0>(basis_mass) = pin->GetReal("units", "mass");
      mass_basis = true;
      basis_count++;
    }
    if (basis_count != 3) {
      std::stringstream msg;
      msg << "### FATAL ERROR in Units constructor" << std::endl
          << "  Counted " << basis_count << " unit basis inputs." << std::endl
          << "  3 unit basis inputs needed." << std::endl
          << "  Choose length, and either time or velocity," << std::endl
          << "  and either ndensity or mass." << std::endl;
      ATHENA_ERROR(msg);
    }
  } else if (unit_system.compare("custom") == 0) {
    // This is the original "custom" unit constructor
    // This must raise error if MLT units are not given in the input file
    code_mass_cgs_   = pin->GetReal("units", "mass_cgs");
    code_length_cgs_ = pin->GetReal("units", "length_cgs");
    code_time_cgs_   = pin->GetReal("units", "time_cgs");

    // Set the basis using the input values
    basis_length   = std::make_tuple(code_length_cgs_/Constants::pc_cgs,"pc");
    basis_time     = std::make_tuple(code_time_cgs_/Constants::Myr_cgs,"Myr");
    basis_velocity = std::make_tuple(code_length_cgs_/code_time_cgs_/Constants::km_s_cgs,
                                     "km/s");
    basis_ndensity = std::make_tuple((code_mass_cgs_/(
        mean_weight*Constants::H_mass_cgs))/CUBE(Constants::pc_cgs),"n/cm^3");
    basis_mass     = std::make_tuple(code_mass_cgs_/Constants::Msun_cgs,"Msun");
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in Units constructor" << std::endl
        << "  unit_system=" << unit_system << " is not valid unit system " << std::endl
        << "  choose one of the default unit systems in [ism] or " << std::endl
        << "  set unit_system=custom and define MLT units manually" << std::endl;
    ATHENA_ERROR(msg);
  }

  // If a basis unit system was used calculate code_mass_cgs_,
  // code_length_cgs_, and code_time_cgs_.
  if (unit_system.compare("custom") != 0) {
    CompleteBasis();

    // Write cgs basis back to the input file
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
  // For constants named code_X_cgs, multiply the code value
  // by the constant to get real (physical) world value.
  // i.e. mass (code units) * code_mass_cgs = mass (g)
  //
  // The same is true for the basis values (basis_X).
  // i.e. mass (code units) * basis_mass = mass (Msun)
  //
  // For constants named X_code, divide the code value by
  // the constant to get the real (pysical) world value.
  // i.e. mass (code units) / Msun_code = mass (Msun)
  // or energy (code units) / erg_code = energy (ergs)
  // set public MLT unit variable
  code_length_cgs = code_length_cgs_;
  code_time_cgs   = code_time_cgs_;
  code_mass_cgs   = code_mass_cgs_;

  // variable units in cgs
  code_volume_cgs   = CUBE(code_length_cgs);
  code_density_cgs  = code_mass_cgs/code_volume_cgs;
  code_velocity_cgs = code_length_cgs/code_time_cgs;

  code_energy_cgs = code_mass_cgs*SQR(code_velocity_cgs);
  code_energydensity_cgs = code_pressure_cgs = code_density_cgs*SQR(code_velocity_cgs);

  code_magneticfield_cgs = std::sqrt(4.*PI*code_pressure_cgs);

  code_temp_cgs = code_pressure_cgs/code_density_cgs
                    *Constants::H_mass_cgs/Constants::k_B_cgs;

  // constants in code units
  cm_code     = 1.0/code_length_cgs;
  gram_code   = 1.0/code_mass_cgs;
  second_code = 1.0/code_time_cgs;
  dyne_code   = gram_code*cm_code/(second_code*second_code);
  erg_code    = dyne_code*cm_code;
  kelvin_code = 1.0; // (changgoo) in principle, this should be 1/[code temperature]
                     // but this is what has been adopted in Athena (not sure why)

  Gconst_code = Constants::grav_const_cgs
                     *cm_code*cm_code*cm_code/(gram_code*second_code*second_code);
  Msun_code = Constants::Msun_cgs*gram_code;
  solar_lum_code  = Constants::solar_lum_cgs*erg_code/second_code;

  yr_code   = Constants::yr_cgs*second_code;
  Myr_code  = Constants::Myr_cgs*second_code;
  pc_code   = Constants::pc_cgs*cm_code;
  kpc_code  = Constants::kpc_cgs*cm_code;
  km_s_code = Constants::km_s_cgs*cm_code/second_code;

  H_mass_code = Constants::H_mass_cgs*gram_code;
  rad_aconst_code = Constants::rad_aconst_cgs*erg_code
                         /(cm_code*cm_code*cm_code
                          *kelvin_code*kelvin_code*kelvin_code*kelvin_code);
  k_B_code = Constants::k_B_cgs*erg_code/kelvin_code;
  c_code = Constants::c_cgs*cm_code/second_code;
  echarge_code = Constants::echarge_cgs*std::sqrt(dyne_code*4*PI)*cm_code;
  bethe_code = 1.e51 * erg_code;
}

//========================================================================================
//! \fn void Units::PrintBasisUnits()
//! \brief print basis parameters in chosen units
//========================================================================================
void Units::PrintBasisUnits() {
  Real lval = std::get<0>(basis_length);
  Real tval = std::get<0>(basis_time);
  Real mval = std::get<0>(basis_mass);
  Real vval = std::get<0>(basis_velocity);
  Real nval = std::get<0>(basis_ndensity);
  std::string lunit = std::get<1>(basis_length);
  std::string tunit = std::get<1>(basis_time);
  std::string munit = std::get<1>(basis_mass);
  std::string vunit = std::get<1>(basis_velocity);
  std::string nunit = std::get<1>(basis_ndensity);
  std::cout << "=========== Unit Basis ===========" << std::endl;
  std::cout << "basis_length   = " << lval << " " << lunit << std::endl;
  std::cout << "basis_time     = " << tval << " " << tunit << std::endl;
  std::cout << "basis_mass     = " << mval << " " << munit << std::endl;
  std::cout << "basis_velocity = " << vval << " " << vunit << std::endl;
  std::cout << "basis_ndensity = " << nval << " " << nunit << std::endl;
  std::cout << "==================================" << std::endl;
}

//========================================================================================
//! \fn void Units::PrintCodeUnits()
//! \brief print code units in c.g.s.
//========================================================================================
void Units::PrintCodeUnits() {
  std::cout << "For constants named code_X_cgs, multiply the code value" << std::endl;
  std::cout << "by the constant to get real (physical) world value." << std::endl;
  std::cout << "i.e. mass (code units) * code_mass_cgs = mass (g)" << std::endl;
  std::cout << std::scientific << "=========== Code Units (cgs) ===========" << std::endl;
  std::cout << "code_length_cgs   = " << code_length_cgs   << " cm" << std::endl;
  std::cout << "code_time_cgs     = " << code_time_cgs     << " s" << std::endl;
  std::cout << "code_mass_cgs     = " << code_mass_cgs     << " g" << std::endl;
  std::cout << "code_density_cgs  = " << code_density_cgs  << " g/cm^3" << std::endl;
  std::cout << "code_velocity_cgs = " << code_velocity_cgs << " cm/s" << std::endl;
  std::cout << "code_energy_cgs   = " << code_energy_cgs   << " erg" << std::endl;
  std::cout << "code_pressure_cgs = " << code_pressure_cgs << " erg/cm^2" << std::endl;
  std::cout << "code_temp_cgs     = " << code_temp_cgs     << " K" << std::endl;
  std::cout << "========================================" << std::endl;
}

//========================================================================================
//! \fn void Units::PrintConstantsInCodeUnits()
//! \brief print physical constatns in code units
//========================================================================================
void Units::PrintConstantsInCodeUnits() {
  std::cout << "For constants named X_code, divide the code value by" << std::endl;
  std::cout << "the constant to get the real (pysical) world value." << std::endl;
  std::cout << "i.e. mass (code units) / Msun_code = mass (Msun)" << std::endl;
  std::cout << "or energy (code units) / erg_code = energy (ergs)" << std::endl;
  std::cout << std::scientific << "==== Constants in Code Units ====" << std::endl;
  std::cout << "Myr_code       = " << Myr_code << std::endl;
  std::cout << "pc_code        = " << pc_code << std::endl;
  std::cout << "km_s_code      = " << km_s_code << std::endl;
  std::cout << "Msun_code      = " << Msun_code << std::endl;
  std::cout << "dyne_code      = " << dyne_code << std::endl;
  std::cout << "erg_code       = " << erg_code << std::endl;
  std::cout << "Gconst_code    = " << Gconst_code << std::endl;
  std::cout << "solar_lum_code = " << solar_lum_code << std::endl;
  std::cout << "H_mass_code    = " << H_mass_code << std::endl;
  std::cout << "k_B_code       = " << k_B_code << std::endl;
  std::cout << "c_code         = " << c_code << std::endl;
  std::cout << "echarge_code   = " << echarge_code << std::endl;
  std::cout << "=================================" << std::endl;
}

//========================================================================================
//! \fn void Units::PrintytUnitsOverride()
//! \brief print units override string for yt input
//========================================================================================
void Units::PrintytUnitsOverride() {
  Real lval = std::get<0>(basis_length);
  Real tval = std::get<0>(basis_time);
  Real mval = std::get<0>(basis_mass);
  std::string lunit = std::get<1>(basis_length);
  std::string tunit = std::get<1>(basis_time);
  std::string munit = std::get<1>(basis_mass);
  std::cout << "======= Units Override for yt ======" << std::endl;
  std::cout << "Units override for importing HDF5 files into yt:" << std::endl;
  std::cout << std::fixed << std::setprecision(6)
            << "{\"length_unit\":("  << lval << ",\"" << lunit
            << "\"),\"time_unit\":(" << tval << ",\"" << tunit
            << "\"),\"mass_unit\":(" << mval << ",\"" << munit << "\")}" << std::endl;
  std::cout << "====================================" << std::endl;
}

//========================================================================================
//! \fn void Units::ConvertTimeFromInputFile()
//! \brief Converts time values from input file using code
//!        basis unit conversions
//========================================================================================
void Units::ConvertTimeFromInputFile(ParameterInput *pin) {
  InputBlock *pib = pin->pfirst_block;

  if (pin->DoesParameterExist("time", "start_time")) {
    Real start_time = pin->GetReal("time", "start_time");
    pin->SetReal("time", "start_time", start_time/std::get<0>(basis_time));
  }
  Real tlim = pin->GetReal("time", "tlim");
  pin->SetReal("time", "tlim", tlim/std::get<0>(basis_time));

  // Convert dt from all output blocks
  while (pib != nullptr) {
    if (pib->block_name.compare(0, 6, "output") == 0) {
      std::string block_name;
      block_name.assign(pib->block_name);
      Real dt = pin->GetReal(block_name, "dt");
      pin->SetReal(block_name, "dt", dt/std::get<0>(basis_time));
    }
    pib = pib->pnext;  // move to next input block name
  }
}

//========================================================================================
//! \fn void Units::ConvertMeshFromInputFile()
//! \brief Converts mesh box size values from input file
//!        using code basis unit conversions
//========================================================================================
void Units::ConvertMeshFromInputFile(ParameterInput *pin) {
  InputBlock *pib = pin->pfirst_block;

  Real x1min = pin->GetReal("mesh", "x1min");
  Real x2min = pin->GetReal("mesh", "x2min");
  Real x3min = pin->GetReal("mesh", "x3min");
  Real x1max = pin->GetReal("mesh", "x1max");
  Real x2max = pin->GetReal("mesh", "x2max");
  Real x3max = pin->GetReal("mesh", "x3max");
  pin->SetReal("mesh", "x1min", x1min/std::get<0>(basis_length));
  pin->SetReal("mesh", "x2min", x2min/std::get<0>(basis_length));
  pin->SetReal("mesh", "x3min", x3min/std::get<0>(basis_length));
  pin->SetReal("mesh", "x1max", x1max/std::get<0>(basis_length));
  pin->SetReal("mesh", "x2max", x2max/std::get<0>(basis_length));
  pin->SetReal("mesh", "x3max", x3max/std::get<0>(basis_length));

  // Needed for SMR
  while (pib != nullptr) {
    if (pib->block_name.compare(0, 10, "refinement") == 0) {
      std::string block_name;
      block_name.assign(pib->block_name);
      x1min = pin->GetReal(block_name, "x1min");
      x2min = pin->GetReal(block_name, "x2min");
      x3min = pin->GetReal(block_name, "x3min");
      x1max = pin->GetReal(block_name, "x1max");
      x2max = pin->GetReal(block_name, "x2max");
      x3max = pin->GetReal(block_name, "x3max");
      pin->SetReal(block_name, "x1min", x1min/std::get<0>(basis_length));
      pin->SetReal(block_name, "x2min", x2min/std::get<0>(basis_length));
      pin->SetReal(block_name, "x3min", x3min/std::get<0>(basis_length));
      pin->SetReal(block_name, "x1max", x1max/std::get<0>(basis_length));
      pin->SetReal(block_name, "x2max", x2max/std::get<0>(basis_length));
      pin->SetReal(block_name, "x3max", x3max/std::get<0>(basis_length));
    }
    pib = pib->pnext;  // move to next input block name
  }
}

//========================================================================================
//! \fn void Units::ConvertInputFile()
//! \brief Converts values from input file using code unit conversions
//========================================================================================
void Units::ConvertInputFile(ParameterInput *pin) {
  ConvertTimeFromInputFile(pin);
  ConvertMeshFromInputFile(pin);
}

//========================================================================================
//! \fn void Units::Returncgs()
//! \brief Converts basis value into cgs value
//
//    Additional units should be added here.
//========================================================================================
Real Units::Returncgs(std::string parameter, Real value, std::string unit) {
  Real code_cgs_;
  // Only needed if unit not recognized
  std::stringstream msg;
  msg << "### FATAL ERROR in Units constructor" << std::endl
      << "  " << parameter << " unit: " << unit << " is not valid unit." << std::endl
      << "  Additional units should be added to units.cpp." << std::endl;

  if (parameter == "basis_length") {
    if (unit == "pc") {
      code_cgs_ = Constants::pc_cgs*value;
    } else if (unit == "kpc") {
      code_cgs_ = Constants::kpc_cgs*value;
    } else if (unit == "au") {
      code_cgs_ = Constants::au_cgs*value;
    } else if (unit == "cm") {
      code_cgs_ = value;
    } else if (unit == "m") {
      code_cgs_ = 100*value;
    } else if (unit == "km") {
      code_cgs_ = 1000*value;
    } else { // If more length units are added, they should be added here
      msg << "  Allowed units are: pc, kpc, au, cm, m, or km" << std::endl;
      ATHENA_ERROR(msg);
    }
  } else if (parameter == "basis_time") {
    if (unit == "yr") {
      code_cgs_ = Constants::yr_cgs*value;
    } else if (unit == "Myr") {
      code_cgs_ = Constants::Myr_cgs*value;
    } else if (unit == "s") {
      code_cgs_ = value;
    } else { // If more time units are added, they should be added here
      msg << "  Allowed units are: yr, Myr, or s" << std::endl;
      ATHENA_ERROR(msg);
    }
  } else if (parameter == "basis_velocity") {
    if (unit == "km/s") {
      code_cgs_ = Constants::km_s_cgs*value;
    } else if (unit == "m/s") {
      code_cgs_ = 100*value;
    } else if (unit == "cm/s") {
      code_cgs_ = value;
    } else { // If more velocity units are added, they should be added here
      msg << "  Allowed units are: km/s, m/s, or cm/s" << std::endl;
      ATHENA_ERROR(msg);
    }
  } else if (parameter == "basis_ndensity") {
    if (unit == "n/cm^3") {
      code_cgs_ = value;
    } else if (unit == "n/m^3") {
      code_cgs_ = value/CUBE(100.0);
    } else { // If more ndensity units are added, they should be added here
      msg << "  Allowed units are: n/cm^3, or n/m^3" << std::endl;
      ATHENA_ERROR(msg);
    }
  } else if (parameter == "basis_mass") {
    if (unit == "Msun") {
      code_cgs_ = Constants::Msun_cgs*value;
    } else if (unit == "g") {
      code_cgs_ = value;
    } else if (unit == "kg") {
      code_cgs_ = 1000*value;
    } else { // If more mass units are added, they should be added here
      msg << "  Allowed units are: Msun, g, or kg" << std::endl;
      ATHENA_ERROR(msg);
    }
  } else {
    code_cgs_ = 0.0;
    std::stringstream msg2;
    msg2 << "### FATAL ERROR in Units constructor" << std::endl
         << "   Parameter: " << parameter << " is not a valid choice" << std::endl;
      ATHENA_ERROR(msg2);
  }
  return code_cgs_;
}

//========================================================================================
//! \fn void Units::CompleteBasis()
//! \brief Completes the basis values and sets code_cgs_
//========================================================================================
void Units::CompleteBasis() {
  // Set code_length_cgs_
  code_length_cgs_ = Returncgs("basis_length",
                      std::get<0>(basis_length),std::get<1>(basis_length));

  // Set code_time_cgs_
  // If velocity was given as the basis, else use time basis
  if (velocity_basis) {
    code_velocity_cgs_ = Returncgs("basis_velocity",
                          std::get<0>(basis_velocity),std::get<1>(basis_velocity));
    code_time_cgs_ = code_length_cgs_/code_velocity_cgs_;
    // Set corresponding time basis
    Real time_conv = Returncgs("basis_time",1.0,std::get<1>(basis_time));
    std::get<0>(basis_time) = code_time_cgs_/time_conv;
  } else {
    code_time_cgs_ = Returncgs("basis_time",
                      std::get<0>(basis_time),std::get<1>(basis_time));
    code_velocity_cgs_ = code_length_cgs_/code_time_cgs_;
    // Set corresponding velocity basis
    Real vel_conv = Returncgs("basis_velocity",1.0,std::get<1>(basis_velocity));
    std::get<0>(basis_velocity) = code_velocity_cgs_/vel_conv;
  }

  // Set code_mass_cgs_
  // If mass was given as the basis, else use ndensity basis
  if (mass_basis) {
    code_mass_cgs_ = Returncgs("basis_mass",
                      std::get<0>(basis_mass),std::get<1>(basis_mass));
    code_ndensity_cgs_ = code_mass_cgs_/(
        mean_weight*Constants::H_mass_cgs*CUBE(code_length_cgs_));
    // Set corresponding ndensity basis
    Real nden_conv = Returncgs("basis_ndensity",1.0,std::get<1>(basis_ndensity));
    std::get<0>(basis_ndensity) = code_ndensity_cgs_/nden_conv;
  } else {
    code_ndensity_cgs_ = Returncgs("basis_ndensity",
                          std::get<0>(basis_ndensity),std::get<1>(basis_ndensity));
    // Using length and ndensity to get mass basis
    code_mass_cgs_ = mean_weight*Constants::H_mass_cgs*
                      CUBE(code_length_cgs_)*code_ndensity_cgs_;
    Real mass_conv = Returncgs("basis_mass",1.0,std::get<1>(basis_mass));
    std::get<0>(basis_mass) = code_mass_cgs_/mass_conv;
  }
}
