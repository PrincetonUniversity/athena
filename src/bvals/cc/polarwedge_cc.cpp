//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file polarwedge_cc.cpp
//! \brief implementation of polar wedge BCs in x2 direction

// C headers

// C++ headers

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "bvals_cc.hpp"

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::PolarWedgeInnerX2(
//!         Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh)
//! \brief polar wedge boundary conditions, inner x2 boundary

void CellCenteredBoundaryVariable::PolarWedgeInnerX2(
    Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh) {
  for (int n=0; n<=nu_; ++n) {
    Real sign = 1.0;
    if (flip_across_pole_ != nullptr) sign = flip_across_pole_[n] ? -1.0 : 1.0;
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          (*var_cc)(n,k,jl-j,i) = sign * (*var_cc)(n,k,jl+j-1,i);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::PolarWedgeOuterX2(
//!         Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh)
//! \brief polar wedge boundary conditions, outer x2 boundary

void CellCenteredBoundaryVariable::PolarWedgeOuterX2(
    Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh) {
  for (int n=0; n<=nu_; ++n) {
    Real sign = 1.0;
    if (flip_across_pole_ != nullptr) sign = flip_across_pole_[n] ? -1.0 : 1.0;
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          (*var_cc)(n,k,ju+j,i) = sign * (*var_cc)(n,k,ju-j+1,i);
        }
      }
    }
  }
  return;
}
