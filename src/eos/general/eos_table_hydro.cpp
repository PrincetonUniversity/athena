//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file eos_table_hydro.cpp
//  \brief implements functions in class EquationOfState for an EOS lookup table
//======================================================================================

// C headers

// C++ headers
#include <cfloat>  // FLT_MIN
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

void EosTestLoop(EquationOfState *peos);

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::SimplePres(Real rho, Real egas)
//  \brief Return interpolated gas pressure
Real EquationOfState::SimplePres(Real rho, Real egas) {
  return ptable->GetEosData(0, egas, rho) * egas;
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::SimpleEgas(Real rho, Real pres)
//  \brief Return interpolated internal energy density
Real EquationOfState::SimpleEgas(Real rho, Real pres) {
  return ptable->GetEosData(1, pres, rho) * pres;
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::SimpleAsq(Real rho, Real pres)
//  \brief Return interpolated adiabatic sound speed squared
Real EquationOfState::SimpleAsq(Real rho, Real pres) {
  return ptable->GetEosData(2, pres, rho) * pres / rho;
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::RiemannAsq(Real rho, Real hint)
//  \brief Return interpolated adiabatic sound speed squared for use in
//         Riemann solver.
Real EquationOfState::RiemannAsq(Real rho, Real hint) {
  return std::pow(static_cast<Real>(10.0),
                  ptable->table.interpolate(3, std::log10(hint*ptable->hUnit),
                                            std::log10(rho*ptable->rhoUnit))) * hint;
}
//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::PrepEOS(ParameterInput *pin)
//  \brief Read data and initialize interpolated table.
void EquationOfState::PrepEOS(ParameterInput *pin) {
}

//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::CleanEOS()
//  \brief Clear memory/objects used to store EOS data
void EquationOfState::CleanEOS() {
}
