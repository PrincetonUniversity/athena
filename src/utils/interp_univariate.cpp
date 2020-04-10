//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file interp_univariate.cpp
//  \brief Collection of univariate interpolators

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "interp_univariate.hpp"

template<>
Real const InterpolateLagrangeUniform<1>::coeff[3][2] = {
  1, 0,  // for injection
  0.5, 0.5,
  0, 1,  // for injection
};
template<>
Real const InterpolateLagrangeUniform<2>::coeff[3][4] = {
  0, 1, 0, 0,
  -1./16., 9./16., 9./16., -1./16.,
  0, 0, 1, 0,
};
template<>
Real const InterpolateLagrangeUniform<3>::coeff[3][6] = {
  0, 0, 1, 0, 0, 0,
  3./256., -25./256., 75./128., 75./128., -25./256., 3./256.,
  0, 0, 0, 1, 0, 0,
};
template<>
Real const InterpolateLagrangeUniform<4>::coeff[3][8] = {
  0, 0, 0, 1, 0, 0, 0, 0,
  -5./2048., 49./2048., -245./2048., 1225./2048., 1225./2048., -245./2048., 49./2048., -5./2048.,
  0, 0, 0, 0, 1, 0, 0, 0,
};
template<>
Real const InterpolateLagrangeUniform<5>::coeff[3][10] = {
  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
  35./65536., -405./65536., 567./16384., -2205./16384., 19845./32768., 19845./32768., -2205./16384., 567./16384., -405./65536., 35./65536.,
  0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
};