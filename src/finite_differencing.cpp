//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file finite_differencing.cpp
//  \brief High-performance finite differencing kernel

#include "finite_differencing.hpp"

// Centered finite differencing 1st derivative
template<>
Real const FDCenteredStencil<1, 1>::coeff[] = {
  -1./2., 0., 1./2.,
};

template<>
Real const FDCenteredStencil<1, 2>::coeff[] = {
  1./12., -2./3., 0., 2./3., -1./12.,
};

template<>
Real const FDCenteredStencil<1, 3>::coeff[] = {
  -1./60., 3./20., -3./4., 0., 3./4., -3./20., 1/60
};

template<>
Real const FDCenteredStencil<1, 4>::coeff[] = {
  1./280., -4./105., 1./5., -4./5., 0., 4./5., -1./5., 4./105., -1./280.
};

// Centered finite differencing 2nd derivative
template<>
Real const FDCenteredStencil<2, 1>::coeff[] = {
  1., -2., 1.,
};

template<>
Real const FDCenteredStencil<2, 2>::coeff[] = {
  -1./12., 4./3., -5./2., 4./3., -1./12.
};

template<>
Real const FDCenteredStencil<2, 3>::coeff[] = {
  1./90., -3./20., 3./2., -49./18., 3./2., -3./20., 1./90.
};

template<>
Real const FDCenteredStencil<2, 4>::coeff[] = {
  -1./560., 8./315., -1./5., 8./5., -205./72., 8./5., -1./5., 8./315., -1./560.
};

// High order derivative operators for Kreiss-Oliger dissipation
template<>
Real const FDCenteredStencil<4, 2>::coeff[] = {
  1., -4., 6., -4., 1.,
};

template<>
Real const FDCenteredStencil<6, 3>::coeff[] = {
  1., -6., 15., -20., 15., -6., 1.,
};

template<>
Real const FDCenteredStencil<8, 4>::coeff[] = {
  1., -8., 28., -56., 70., -56., 28., -8., 1.,
};

template<>
Real const FDCenteredStencil<10, 5>::coeff[] = {
  1., -10., 45., -120., 210., -252., 210., -120., 45., -10., 1.,
};

// Left-biased finite differencing 1st derivative
template<>
Real const FDLeftBiasedStencil<1, 1, 0>::coeff[] = {
  -1., 1.,
};

template<>
Real const FDLeftBiasedStencil<1, 2, 0>::coeff[] = {
  0.5, -2.0, 1.5,
};

template<>
Real const FDLeftBiasedStencil<1, 2, 1>::coeff[] = {
  1./6., -1., 1./2., 1./3.,
};

template<>
Real const FDLeftBiasedStencil<1, 3, 1>::coeff[] = {
  -1./12, 6./12, -18./12, +10./12, +3./12,
};

template<>
Real const FDLeftBiasedStencil<1, 3, 2>::coeff[] = {
  -1./30., 1./4., -1., 1./3., 1./2., -1./20.,
};

template<>
Real const FDLeftBiasedStencil<1, 4, 2>::coeff[] = {
  1./60., -2./15., 1./2., -4./3., 7./12., 2./5., -1./30.,
};

template<>
Real const FDLeftBiasedStencil<1, 4, 3>::coeff[] = {
  1./140., -1./15., 3./10., -1., 1./4., 3./5., -1./10., 1./105.,
};

template<>
Real const FDLeftBiasedStencil<1, 5, 3>::coeff[] = {
  -1./280., 1./28., -1./6., 1./2., -5./4., 9./20., 1./2., -1./14., 1./168.,
};

template<>
Real const FDLeftBiasedStencil<1, 5, 4>::coeff[] ={
  -1./630., 1./56., -2./21., 1./3., -1., 1./5., 2./3., -1./7., 1./42., -1./504.,
};

// Right-biased finite differencing 1st derivative
template<>
Real const FDRightBiasedStencil<1, 1, 0>::coeff[] = {
  -1., 1.,
};

template<>
Real const FDRightBiasedStencil<1, 2, 0>::coeff[] = {
  -1.5, 2.0, -0.5,
};

template<>
Real const FDRightBiasedStencil<1, 2, 1>::coeff[] = {
  -1./3., -1./2., 1., -1./6.,
};

template<>
Real const FDRightBiasedStencil<1, 3, 1>::coeff[] = {
   -3./12., -10./12., 18./12., -6./12., 1./12.,
};

template<>
Real const FDRightBiasedStencil<1, 3, 2>::coeff[] = {
  1./20., -1./2., -1./3., 1., -1./4., 1./30.,
};

template<>
Real const FDRightBiasedStencil<1, 4, 2>::coeff[] = {
  1./30., -2./5., -7./12., 4./3., -1./2., 2./15., -1./60.,
};

template<>
Real const FDRightBiasedStencil<1, 4, 3>::coeff[] = {
  -1./105., 1./10., -3./5., -1./4., 1., -3./10., 1./15., -1./140.,
};

template<>
Real const FDRightBiasedStencil<1, 5, 3>::coeff[] = {
  -1./168, 1./14, -1./2, -9./20, 5./4, -1./2, 1./6, -1./28, 1./280
};

template<>
Real const FDRightBiasedStencil<1, 5, 4>::coeff[] = {
  1./504., -1./42., 1./7., -2./3., -1./5., 1., -1./3., 2./21., -1./56., 1./630.,
};
