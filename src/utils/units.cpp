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

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "units.hpp"

//! \brief ctor without mu
Units::Units(Real dunit, Real lunit, Real vunit) {
  Density = dunit;
  Length = lunit;
  Velocity = vunit;
  mu = 1.0;
  fixed_mu = false;

  SetUnitsConstants();
}

//! \brief ctor with fixed mu
Units::Units(Real dunit, Real lunit, Real vunit, Real mu0) {
  Density = dunit;
  Length = lunit;
  Velocity = vunit;
  mu = mu0;
  fixed_mu = true;

  SetUnitsConstants();
}

void Units::SetUnitsConstants() {
  Volume = Length*Length*Length;
  Mass = Density*Volume;
  Time = Length/Velocity;

  EnergyDensity = Pressure = Density*Velocity*Velocity;

  MagneticField = std::sqrt(4.*PI*Pressure);

  Temperature = Pressure/Density*mu*Constants::mH/Constants::kB;

  cm = 1.0/Length;
  gram = 1.0/Mass;
  second = 1.0/Time;
  dyne = gram*cm/(second*second);
  erg = gram*cm*cm/(second*second);
  Kelvin = 1.0;

  G_in_code = Constants::G * cm*cm*cm/(gram*second*second);
  Msun_in_code = Constants::Msun * gram;
  Lsun_in_code = Constants::Lsun * erg/second;
  Myr_in_code = Constants::Myr * second;
  pc_in_code = Constants::pc * cm;
  kpc_in_code = Constants::kpc * cm;
  kms_in_code = Constants::kms * cm/second;
  mH_in_code = Constants::mH * gram;
  aR_in_code = Constants::aR * erg/(cm*cm*cm*Kelvin*Kelvin*Kelvin*Kelvin);
  kB_in_code = Constants::kB * erg/Kelvin;
  c_in_code = Constants::c * cm/second;
  e_in_code = Constants::e * std::sqrt(dyne*4*PI)*cm;
  Bethe_in_code = 1.e51 * erg;
}

void Units::PrintCodeUnits() {
  std::cout << std::scientific << "============ Code Units ============" << std::endl;

  std::cout << "code Density = " << Density << " c.g.s." << std::endl;
  std::cout << "code Length = " << Length << " c.g.s." << std::endl;
  std::cout << "code Velocity = " << Velocity << " c.g.s." << std::endl;

  std::cout << "code Mass = " << Mass << " c.g.s." << std::endl;
  std::cout << "code Time = " << Time << " c.g.s." << std::endl;
  std::cout << "code Pressure = " << Pressure << " c.g.s." << std::endl;

  std::cout << "code Temperature = " << Temperature << " c.g.s." << std::endl;
  if (fixed_mu)
    std::cout << "    with fixed mu0 = " << mu << std::endl;
  else
    std::cout << "    for mu = 1.0 (temperature conversion is non-trivial for general mu)"
              << std::endl;
  std::cout << "====================================" << std::endl;
}

void Units::PrintConstantsInCodeUnits() {
  std::cout << std::scientific << "===== Constants  in Code Units =====" << std::endl;

  std::cout << "dyne in code = " << dyne << std::endl;
  std::cout << "erg in code = " << erg << std::endl;

  std::cout << "Gconst in code = " << G_in_code << std::endl;
  std::cout << "Msun in code = " << Msun_in_code << std::endl;
  std::cout << "Lsun in code = " << Lsun_in_code << std::endl;
  std::cout << "Myr in code = " << Myr_in_code << std::endl;
  std::cout << "kB in code = " << kB_in_code << std::endl;
  std::cout << "c in code = " << c_in_code << std::endl;
  std::cout << "e in code = " << e_in_code << std::endl;

  std::cout << "P/kB conversion = " << Pressure/Constants::kB << std::endl;
  std::cout << "====================================" << std::endl;
}
