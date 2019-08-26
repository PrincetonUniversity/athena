//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file weno5.cpp
//  \brief fifth-order accurate weighted essentially non-oscillatory (WENO_5) scheme
//
//  Operates on the entire nx4 range of a single AthenaArray<Real> input (no MHD).
//  No assumptions of hydrodynamic fluid variable input; no characteristic projection.
//
// REFERENCES:
//
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

// number of stencils per cell; WENO5 consists of 3x three-point substencils per cell
#define NSUBSTENCIL 3

//--------------------------------------------------------------------------------------
//! \fn Reconstruction::WenoFifthX1()
//  \brief Uniform Cartesian mesh WENO5 implementation

void Reconstruction::WenoFifthX1(const int k, const int j, const int il, const int iu,
                                 const AthenaArray<Real> &q,
                                 AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  return;
}

//--------------------------------------------------------------------------------------
//! \fn Reconstruction::WenoFifthX2()
//  \brief Uniform Cartesian mesh WENO5 implementation

void Reconstruction::WenoFifthX2(const int k, const int j, const int il, const int iu,
                                 const AthenaArray<Real> &q,
                                 AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  return;
}

//--------------------------------------------------------------------------------------
//! \fn Reconstruction::WenoFifthX3()
//  \brief Uniform Cartesian mesh WENO5 implementation

void Reconstruction::WenoFifthX3(const int k, const int j, const int il, const int iu,
                                 const AthenaArray<Real> &q,
                                 AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  return;
}
