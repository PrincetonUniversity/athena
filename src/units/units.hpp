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
#include <iomanip>
#include <string>
#include <tuple>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

//! \brief Physical constants defined in c.g.s.
namespace Constants {
static const Real grav_const_cgs = 6.67259e-8;
// Changed from 1.9891e+33 to agree with Astropy:
static const Real Msun_cgs       = 1.9884099e33;
static const Real solar_lum_cgs  = 3.8268e33;
// Changed from 3.155815e7 to agree with Astropy and the FITS standard:
static const Real yr_cgs         = 3.15576e7;
static const Real Myr_cgs        = 3.15576e13;
// Changed from 3.085678e+18 to agree with Astropy:
static const Real pc_cgs         = 3.08567758e18;
static const Real kpc_cgs        = 3.08567758e21;
static const Real au_cgs         = 1.495978707e13;
static const Real km_s_cgs       = 1.0e5;
static const Real H_mass_cgs     = 1.6733e-24;
static const Real rad_aconst_cgs = 7.5646e-15;
static const Real k_B_cgs        = 1.380658e-16;
static const Real c_cgs          = 2.99792458e10;
static const Real echarge_cgs    = 4.80320427e-10;
static const Real kelvin_cgs     = 1;
static const Real yrinmyr        = 1e6;
static const Real cminkm         = 100000;
static const Real thomson        = 6.6524e-25;       // Thompson cross section (cm^2)
static const Real mp             = 1.6726231e-24;    // Mass of a proton (g)
} // namespace Constants

//! \brief Class for units
class Units {
 public:
  explicit Units(ParameterInput *pin);

  void SetUnitsConstants();
  void PrintCodeUnits();
  void PrintBasisUnits();
  void PrintConstantsInCodeUnits();
  void PrintytUnitsOverride();
  Real Returncgs(std::string parameter,Real value,std::string unit);
  void CompleteBasis();
  void ConvertInputFile(ParameterInput *pin);
  void ConvertTimeFromInputFile(ParameterInput *pin);
  void ConvertMeshFromInputFile(ParameterInput *pin);

  // unit system
  std::string unit_system;

  // Code to physical units conversion basis
  // Value and unit stored in a tuple (val,unit)
  // To access basis value use: std::get<0>(basis_X)
  // To access basis unit use : std::get<1>(basis_X)
  //
  // User chooses a 'length' basis (with possible units):
  //   length (pc, kpc, au, cm, m, km)
  // Then either 'time' or 'velocity' (with possible units):
  //   time (yr, Myr, s)
  //   velocity (km/s, cm/s, m/s)
  // Then either 'ndensity' or 'mass' (with possible units):
  //   ndensity (n/cm^3, n/m^3)
  //   mass (Msun, g, kg)

  std::tuple<Real, std::string> basis_length;
  std::tuple<Real, std::string> basis_time;
  std::tuple<Real, std::string> basis_velocity;
  std::tuple<Real, std::string> basis_ndensity;
  std::tuple<Real, std::string> basis_mass;

  // code units in c.g.s.
  // i.e. multiply this to convert quantities in c.g.s.
  // public accessor for default MLT units
  // ideally this should be accessor methods, but
  // no on-the-fly change in unit system is assumed

  Real code_mass_cgs, code_length_cgs, code_time_cgs;

  Real mean_weight;

  Real code_volume_cgs, code_density_cgs, code_velocity_cgs;
  Real code_ndensity;
  Real code_energy_cgs, code_energydensity_cgs, code_pressure_cgs;
  Real code_magneticfield_cgs;
  Real code_temp_cgs; // T/mu

  // c.g.s. units in code units
  Real gram_code, cm_code, second_code, dyne_code, erg_code, kelvin_code;

  // physical constants in code units
  Real Gconst_code;
  Real Msun_code;
  Real solar_lum_code;
  Real yr_code;
  Real Myr_code; // Changed from million_yr to Myr to agree with the FITS standard
  Real pc_code;
  Real kpc_code;
  Real km_s_code;
  Real H_mass_code;
  Real rad_aconst_code; // aR
  Real k_B_code; // k_B
  Real c_code;
  Real echarge_code;
  Real bethe_code; // 1.e51 erg

 private:
  // code MLT units in c.g.s.
  Real code_mass_cgs_, code_length_cgs_, code_time_cgs_;
  Real code_velocity_cgs_, code_ndensity_cgs_;

  bool velocity_basis = false;
  bool mass_basis = false;
};
#endif // UNITS_UNITS_HPP_
