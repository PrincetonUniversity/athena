#ifndef UTILS_GL_QUADRATURE_HPP_
#define UTILS_GL_QUADRATURE_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gl_quadrature.hpp
//! \brief functions for computing the Gauss-Legendre (GL) quadrature of a given 1D, 2D,
//!  or 3D function. Provided for convenience / intended for use in pgen/*.cpp files.
//!
//! \todo (felker):
//! - add other Gaussian quadratures, or alternative approaches for computing
//!   the initial condition that outperform GL quadrature for a discontinuous function

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

namespace GaussLegendre {
// TODO(felker): for more complicate f(), use functors/lambdas/pass AthenaArray to both fn
// 1D f(x1)
Real integrate(const int n, Real (*f)(Real), Real x1l, Real x1u);
// 2D f(x1, x2)
Real integrate(const int n, Real (*f)(Real, Real),
               Real x1l, Real x1u, Real x2l, Real x2u);
// 3D f(x1, x2, x3)
Real integrate(const int n, Real (*f)(Real, Real, Real),
               Real x1l, Real x1u, Real x2l, Real x2u, Real x3l, Real x3u);
} // namespace GaussLegendre
#endif // UTILS_GL_QUADRATURE_HPP_
