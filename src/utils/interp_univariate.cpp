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

// centered stencils
template<>
Real const InterpolateLagrangeUniform<1>::coeff[3][2] = {
  1, 0,         // for injection [interp. @ poly root]
  1./2., 1./2., // interp. to midpoint of stencil
  0, 1,         // injected
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

// right biased stencils
template<>
Real const InterpolateLagrangeUniformBiasR<2>::coeff[] = {
  2., -1.
};
template<>
Real const InterpolateLagrangeUniformBiasR<3>::coeff[] = {
  3., -3., 1.
};
template<>
Real const InterpolateLagrangeUniformBiasR<4>::coeff[] = {
  4., -6., 4., -1.
};
template<>
Real const InterpolateLagrangeUniformBiasR<5>::coeff[] = {
  5., -10., 10., -5., 1.
};
template<>
Real const InterpolateLagrangeUniformBiasR<6>::coeff[] = {
  6., -15., 20., -15., 6., -1.
};
template<>
Real const InterpolateLagrangeUniformBiasR<7>::coeff[] = {
  7., -21., 35., -35., 21., -7., 1.
};
template<>
Real const InterpolateLagrangeUniformBiasR<8>::coeff[] = {
  8., -28., 56., -70., 56., -28., 8., -1.
};
template<>
Real const InterpolateLagrangeUniformBiasR<9>::coeff[] = {
  9., -36., 84., -126., 126., -84., 36., -9., 1.
};
template<>
Real const InterpolateLagrangeUniformBiasR<10>::coeff[] = {
  10., -45., 120., -210., 252., -210., 120., -45., 10., -1.
};

// left biased stencils
template<>
Real const InterpolateLagrangeUniformBiasL<2>::coeff[] = {
  -1., 2.
};
template<>
Real const InterpolateLagrangeUniformBiasL<3>::coeff[] = {
  1., -3., 3.
};
template<>
Real const InterpolateLagrangeUniformBiasL<4>::coeff[] = {
  -1., 4., -6., 4.
};
template<>
Real const InterpolateLagrangeUniformBiasL<5>::coeff[] = {
  1., -5., 10., -10., 5.
};
template<>
Real const InterpolateLagrangeUniformBiasL<6>::coeff[] = {
  -1., 6., -15., 20., -15., 6.
};
template<>
Real const InterpolateLagrangeUniformBiasL<7>::coeff[] = {
  1., -7., 21., -35., 35., -21., 7.
};
template<>
Real const InterpolateLagrangeUniformBiasL<8>::coeff[] = {
  -1., 8., -28., 56., -70., 56., -28., 8.
};
template<>
Real const InterpolateLagrangeUniformBiasL<9>::coeff[] = {
  1., -9., 36., -84., 126., -126., 84., -36., 9.
};
template<>
Real const InterpolateLagrangeUniformBiasL<10>::coeff[] = {
  -1., 10., -45., 120., -210., 252., -210., 120., -45., 10.
};


// centered stencils
template<>
Real const InterpolateLagrangeUniform_opt<1>::coeff[1] = {
  1./2., // interp. to midpoint of stencil
};
template<>
Real const InterpolateLagrangeUniform_opt<2>::coeff[2] = {
  -1./16., 9./16.,
};
template<>
Real const InterpolateLagrangeUniform_opt<3>::coeff[3] = {
  3./256., -25./256., 75./128.,
};
template<>
Real const InterpolateLagrangeUniform_opt<4>::coeff[4] = {
  -5./2048., 49./2048., -245./2048., 1225./2048.,
};
template<>
Real const InterpolateLagrangeUniform_opt<5>::coeff[5] = {
  35./65536., -405./65536., 567./16384., -2205./16384., 19845./32768.,
};