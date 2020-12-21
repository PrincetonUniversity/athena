//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file eos_table.cpp
//! \brief implements functions in class EquationOfState for an EOS lookup table
//======================================================================================

// C headers

// C++ headers
#include <cmath>   // sqrt()
#include <fstream>
#include <iostream> // ifstream
#include <sstream>
#include <stdexcept> // std::invalid_argument
#include <string>

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../field/field.hpp"
#include "../../parameter_input.hpp"
#include "../../utils/interp_table.hpp"
#include "../eos.hpp"

namespace {
Real dens_pow = -1.0;

//----------------------------------------------------------------------------------------
//! \fn Real GetEosData(EosTable *ptable, int kOut, Real var, Real rho)
//! \brief Gets interpolated data from EOS table assuming 'var' has dimensions
//!        of energy per volume.
inline Real GetEosData(EosTable *ptable, int kOut, Real var, Real rho) {
  Real x1 = std::log10(rho * ptable->rhoUnit);
  Real x2 = std::log10(var * ptable->EosRatios(kOut) * ptable->eUnit) + dens_pow * x1;
  return std::pow((Real)10, ptable->table.interpolate(kOut, x2, x1));
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::PresFromRhoEg(Real rho, Real egas)
//! \brief Return interpolated gas pressure
Real EquationOfState::PresFromRhoEg(Real rho, Real egas) {
  return GetEosData(ptable, 0, egas, rho) * egas;
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::EgasFromRhoP(Real rho, Real pres)
//! \brief Return interpolated internal energy density
Real EquationOfState::EgasFromRhoP(Real rho, Real pres) {
  return GetEosData(ptable, 1, pres, rho) * pres;
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::AsqFromRhoP(Real rho, Real pres)
//! \brief Return interpolated adiabatic sound speed squared
Real EquationOfState::AsqFromRhoP(Real rho, Real pres) {
  return GetEosData(ptable, 2, pres, rho) * pres / rho;
}

//----------------------------------------------------------------------------------------
//! void EquationOfState::InitEosConstants(ParameterInput* pin)
//! \brief Initialize constants for EOS
void EquationOfState::InitEosConstants(ParameterInput* pin) {
  dens_pow = pin->GetOrAddReal("hydro", "dens_pow", dens_pow);
  return;
}
