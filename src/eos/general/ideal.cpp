//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file ideal.cpp
//! \brief implements ideal EOS in general EOS framework, mostly for debuging
//======================================================================================

// C headers

// C++ headers

// Athena++ headers
#include "../eos.hpp"

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::PresFromRhoEg(Real rho, Real egas, Real* s)
//! \brief Return gas pressure
Real EquationOfState::PresFromRhoEg(Real rho, Real egas, Real* s) {
  return (gamma_ - 1.) * egas;
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::EgasFromRhoP(Real rho, Real pres, Real* r)
//! \brief Return internal energy density
Real EquationOfState::EgasFromRhoP(Real rho, Real pres, Real* r) {
  return pres / (gamma_ - 1.);
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::AsqFromRhoP(Real rho, Real pres, const Real* r)
//! \brief Return adiabatic sound speed squared
Real EquationOfState::AsqFromRhoP(Real rho, Real pres, const Real* r) {
  return gamma_ * pres / rho;
}

//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::InitEosConstants(ParameterInput* pin)
//! \brief Initialize constants for EOS
void EquationOfState::InitEosConstants(ParameterInput *pin) {
  return;
}
