#ifndef INTERP_UNIVARIATE_HPP_
#define INTERP_UNIVARIATE_HPP_

//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file interp_univariate.hpp
//  \brief Collection of univariate interpolators

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

// uniform grid assumed; bias_ controls location of interp. point
// bias_ = 0 corresponds to a point exactly centered
template<int half_stencil_size_>
class InterpolateLagrangeUniform {
  public:
    // order of convergence (in spacing) for "derivative-dominated" functions
    enum {interpolation_order = 2 * half_stencil_size_ - 1};
    enum {npoints = 2 * half_stencil_size_};
    static Real const coeff[3][npoints];
};

// uniform grid; stencil fully biased towards right (target left)
// use for extrapolation
template<int stencil_size_>
class InterpolateLagrangeUniformBiasR {
  public:
    enum {interpolation_order = stencil_size_ - 1};
    enum {npoints = stencil_size_};
    static Real const coeff[npoints];
};
// fully bias left (target right)
template<int stencil_size_>
class InterpolateLagrangeUniformBiasL {
  public:
    enum {interpolation_order = stencil_size_ - 1};
    enum {npoints = stencil_size_};
    static Real const coeff[npoints];
};

template<int half_stencil_size_>
class InterpolateLagrangeUniform_opt {
  public:
    enum {interpolation_order = 2 * half_stencil_size_ - 1};
    enum {npoints = half_stencil_size_};
    static Real const coeff[npoints];
};

#endif