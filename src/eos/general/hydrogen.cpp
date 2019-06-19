//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file hydrogen_hydro.cpp
//  \brief implements functions in class EquationOfState for simple hydrogen EOS
//======================================================================================

// C headers

// C++ headers
#include <algorithm>
#include <cmath>   // sqrt()
#include <fstream>
#include <iostream> // ifstream
#include <limits>   // std::numeric_limits<float>::epsilon()
#include <sstream>
#include <stdexcept> // std::invalid_argument

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../field/field.hpp"
#include "../../parameter_input.hpp"
#include "../eos.hpp"

namespace {
const Real float_eps = std::numeric_limits<float>::epsilon();
const Real my_1pe = 1.0 + float_eps;

//----------------------------------------------------------------------------------------
//! \fn Real x_(Real rho, Real T) {
//  \brief compute ionization fraction
Real x_(Real rho, Real T) {
  return 2./(1.+std::sqrt(1.+4.*std::exp(1./T-1.5*std::log(T)+std::log(rho))));
}

//----------------------------------------------------------------------------------------
//! \fn Real P_of_rho_T(Real rho, Real T) {
//  \brief compute gas pressure
Real P_of_rho_T(Real rho, Real T) {
  return rho * T * (1. + x_(rho, T));
}

//----------------------------------------------------------------------------------------
//! \fn Real e_of_rho_T(Real rho, Real T) {
//  \brief compute internal energy density
Real e_of_rho_T(Real rho, Real T) {
  Real x = x_(rho, T);
  return x * rho + 1.5 * rho * T * (1. + x);
}

//----------------------------------------------------------------------------------------
//! \fn Real h_of_rho_T(Real rho, Real T) {
//  \brief compute specific enthalpy
Real h_of_rho_T(Real rho, Real T) {
  Real x = x_(rho, T);
  return x + 2.5 * T * (1. + x);
}

//----------------------------------------------------------------------------------------
//! \fn Real asq_(Real rho, Real T) {
//  \brief compute adiabatic sound speed squared
Real asq_(Real rho, Real T) {
  Real x = x_(rho, T);
  Real lt = std::log(T);
  Real b = 8. * rho * std::exp(-1.25 * lt - .5 / T);
  Real c = std::exp(1.5 * lt - .5 / T);
  c = (std::sqrt(c) + std::sqrt(c + 4. * rho));
  b /= c*c*c;
  c = 2. + x - SQR(x);
  return (1.+x)*T*(b*(4.+20.*T+15.*SQR(T))+10.*c)/(b*SQR(2.+3.*T)+6.*c);
}

//----------------------------------------------------------------------------------------
//! \fn Real invert(Real(*f) (Real, Real), Real rho, Real sol, Real T0, Real T1)
//  \brief Invert an EOS function at a given density and return temperature.
//         Uses simple bisection for inversion.
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
  } else {
    return T0;
  }

  if ( f1 < 0 ) {
    fa = f1;
    Ta = T1;
  } else if ( f1 > 0 ) {
    fb = f1;
    Tb = T1;
  } else {
    return T1;
  }

  if ( (Ta < 0) || (Tb < 0) ) {
    std::stringstream msg;
    msg << "### FATAL ERROR in EquationOfState inversion"
        << std::endl << "Root not bracketed" << std::endl;
    ATHENA_ERROR(msg);
  }

  while ( (std::fabs(Ta - Tb) >= prec) && (fb - fa >= prec) ) {
    T0 = 0.5 * Ta + 0.5 * Tb;
    f0 = f(rho, T0) / sol - 1.0;
    if ( f0 < 0 ) {
      fa = f0;
      Ta = T0;
    } else if ( f0 > 0 ) {
      fb = f0;
      Tb = T0;
    } else {
      return T0;
    }
  }
  return T0;
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::RiemannAsq(Real rho, Real hint)
//  \brief Return adiabatic sound speed squared for use in Riemann solver.
Real EquationOfState::RiemannAsq(Real rho, Real hint) {
  Real T = invert(*h_of_rho_T, rho, hint, 0.2*std::max(hint - 1.0, 0.1*hint),
                  my_1pe*0.4*hint);
  return asq_(rho, T);
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::PresFromRhoEg(Real rho, Real egas)
//  \brief Return gas pressure
Real EquationOfState::PresFromRhoEg(Real rho, Real egas) {
  Real es = egas / rho;
  Real T = invert(*e_of_rho_T, rho, egas, std::max(es - 1.0, 0.1*es)/3.0,
                  my_1pe*2.0*es/3.0);
  return P_of_rho_T(rho, T);
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::EgasFromRhoP(Real rho, Real pres)
//  \brief Return internal energy density
Real EquationOfState::EgasFromRhoP(Real rho, Real pres) {
  Real ps = pres / rho;
  Real T = invert(*P_of_rho_T, rho, pres, 0.5*ps, my_1pe*ps);
  return e_of_rho_T(rho, T);
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::AsqFromRhoP(Real rho, Real pres)
//  \brief Return adiabatic sound speed squared
Real EquationOfState::AsqFromRhoP(Real rho, Real pres) {
  Real ps = pres / rho;
  Real T = invert(*P_of_rho_T, rho, pres, 0.5*ps, my_1pe*ps);
  return asq_(rho, T);
}
