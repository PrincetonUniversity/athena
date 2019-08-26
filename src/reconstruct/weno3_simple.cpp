//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file weno3_simple.cpp
//  \brief third-order accurate weighted essentially non-oscillatory (WENO_3) scheme
//
//  Operates on the entire nx4 range of a single AthenaArray<Real> input (no MHD).
//  No assumptions of hydrodynamic fluid variable input; no characteristic projection.
//
// REFERENCES:
// (Mignone) A. Mignone, "High-order conservative reconstruction schemes for finite volume
// methods in cylindrical and spherical coordinates", JCP, 270, 784 (2014)
//
// (Classical WENO3) G.C. Jiang, C.W. Shu, Efficient implementation of weighted ENO
// schemes, J. Comput. Phys. 126 (1996) 202–228
//
// (ESWENO) N.K. Yamaleev, M.H. Carpenter, Third-order energy stable WENO scheme,
// J. Comput. Phys. 228 (2009) 3025–3047
//
// A. Mignone, P. Tzeferacos, G. Bodo, High-order conservative finite difference GLM-MHD
// schemes for cell-centered MHD, J. Comput. Phys. 229 (2010) 5896–5920
//======================================================================================

// C headers

// C++ headers
#include <algorithm>    // max()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "reconstruction.hpp"

//--------------------------------------------------------------------------------------
//! \fn Reconstruction::WenoThirdX1()
//  \brief Uniform Cartesian mesh WENO3 implementation. Using nonlinear weights from
// Mignone (2014)

void Reconstruction::WenoThirdX1(const int k, const int j, const int il, const int iu,
                                 const AthenaArray<Real> &q,
                                 AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  return;
}

//--------------------------------------------------------------------------------------
//! \fn Reconstruction::WenoThirdX2()
//  \brief

void Reconstruction::WenoThirdX2(const int k, const int j, const int il, const int iu,
                                 const AthenaArray<Real> &q,
                                 AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  return;
}

//--------------------------------------------------------------------------------------
//! \fn Reconstruction::WenoThirdX3()
//  \brief

void Reconstruction::WenoThirdX3(const int k, const int j, const int il, const int iu,
                                 const AthenaArray<Real> &q,
                                 AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  return;
}
