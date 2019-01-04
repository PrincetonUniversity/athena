//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file noop.cpp
//  \brief Implements no-op versions of the general eos functions

// C/C++ headers
#include <cmath>   // sqrt()
#include <cfloat>  // FLT_MIN

// Athena++ headers
#include "../eos.hpp"

Real EquationOfState::RiemannAsq(Real rho, Real hint) {
  return -1.0;
}
Real EquationOfState::SimplePres(Real rho, Real egas) {
  return -1.0;
}
Real EquationOfState::SimpleEgas(Real rho, Real pres) {
  return -1.0;
}
Real EquationOfState::SimpleAsq(Real rho, Real pres) {
  return -1.0;
}
