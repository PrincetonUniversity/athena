//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file hydrogen.cpp
//  \brief implements functions in class EquationOfState for an EOS lookup table
//======================================================================================

// C++ headers
#include <cmath>   // sqrt()
#include <cfloat>  // FLT_MIN
#include <iostream> // ifstream
#include <fstream>
#include <sstream>
#include <stdexcept> // std::invalid_argument

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../parameter_input.hpp"
#include "../../field/field.hpp"
#include "../../coordinates/coordinates.hpp"

// this class header
#include "../eos.hpp"

const Real my_1pe = 1. + FLT_EPSILON;

//Real phi_(Real T) {
//  return std::exp(1. / T) * std::pow(T, -1.5);
//}

Real x_(Real rho, Real T) {
  return 2. /(1. + std::sqrt(1. + 4. * std::exp(1. / T - 1.5 * std::log(T) + std::log(rho))));
}

//Real x_T_(Real rho, Real T){
//  Real x = x_(rho, T);
//  return std::pow(T, 3) / (2. + x) * std::exp(1. / T) * std::pow(T, -3.5) * (1. + 1.5 * T) * rho;
//}

Real P_of_rho_T(Real rho, Real T) {
  return rho * T * (1. + x_(rho, T));
}

Real e_of_rho_T(Real rho, Real T) {
  Real x = x_(rho, T);
  return x * rho + 1.5 * rho * T * (1. + x);
}

Real h_of_rho_T(Real rho, Real T){
  Real x = x_(rho, T);
  return x + 2.5 * T * (1. + x);
}

Real asq_(Real rho, Real T) {
  Real x = x_(rho, T);
  Real lt = std::log(T);
  Real b = 8. * rho * std::exp(-1.25 * lt - .5 / T);
  Real c = std::exp(1.5 * lt - .5 / T);
  c = (std::sqrt(c) + std::sqrt(c + 4. * rho));
  b /= c*c*c;
  c = 2. + x - SQR(x);
  return (1.+x)*T*(b * (4. + 20. * T + 15. * SQR(T)) + 10. * c)/(b * SQR(2. + 3. * T) + 6. * c);
}

Real invert(Real(*f) (Real, Real), Real rho, Real sol, Real T0, Real T1) {
  Real f0, f1, fa, fb;
  Real Ta=-1;
  Real Tb=-1;
  const Real prec=1e-12;
  f0 = f(rho, T0) / sol - 1.;
  f1 = f(rho, T1) / sol - 1.;

  if ( f0 < 0 ) {
    fa = f0;
    Ta = T0;
  } else if ( f0 > 0 ) {
    fb = f0;
    Tb = T0;
  }
  else return T0;

  if ( f1 < 0 ) {
    fa = f1;
    Ta = T1;
  } else if ( f1 > 0 ) {
    fb = f1;
    Tb = T1;
  }
  else return T1;

  if ( (Ta < 0) or (Tb < 0) ) {
    std::stringstream msg;
    msg << "### FATAL ERROR in EquationOfState inversion"
        << std::endl << "Root not bracketed" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  while ( (std::fabs(Ta - Tb) >= prec) and (fb - fa >= prec) ){
    T0 = .5 * Ta + .5 * Tb;
    f0 = f(rho, T0) / sol - 1.;
    if ( f0 < 0 ) {
      fa = f0;
      Ta = T0;
    } else if ( f0 > 0 ) {
      fb = f0;
      Tb = T0;
    }
    else return T0;
  }

  return T0;
}

Real EquationOfState::RiemannAsq(Real rho, Real hint) {
  Real T = invert(*h_of_rho_T, rho, hint, .2 * std::max(hint - 1., .1 * hint), my_1pe * .4 * hint);
  return asq_(rho, T);
}

Real EquationOfState::SimplePres(Real rho, Real egas) {
  Real es = egas / rho;
  Real T = invert(*e_of_rho_T, rho, egas, std::max(es - 1., .1 * es)/3., my_1pe * 2. * es / 3.);
  return P_of_rho_T(rho, T);
}

Real EquationOfState::SimpleEgas(Real rho, Real pres) {
  Real ps = pres / rho;
  Real T = invert(*P_of_rho_T, rho, pres, .5 * ps, my_1pe * ps);
  return e_of_rho_T(rho, T);
}

Real EquationOfState::SimpleAsq(Real rho, Real pres) {
  Real ps = pres / rho;
  Real T = invert(*P_of_rho_T, rho, pres, .5 * ps, my_1pe * ps);
  return asq_(rho, T);
}
