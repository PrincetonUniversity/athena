//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file polarwedge_rad.cpp
//  \brief implementation of polar wedge radiation BCs in x2 direction

// C headers

// C++ headers

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../mesh/mesh.hpp"
#include "../../../nr_radiation/radiation.hpp"
#include "bvals_rad.hpp"

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::PolarWedgeInnerX2(
//          Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh)
//  \brief polar wedge boundary conditions, inner x2 boundary

// do not flip the sign of specific intensity
void RadBoundaryVariable::PolarWedgeInnerX2(
    Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh) {
  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {
#pragma omp simd
        for (int n=0; n<=nu_; ++n) {
          (*var_cc)(k,jl-j,i,n) = (*var_cc)(k,jl+j-1,i,n);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::PolarWedgeOuterX2(
//          Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh)
//  \brief polar wedge boundary conditions, outer x2 boundary

void RadBoundaryVariable::PolarWedgeOuterX2(
    Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh) {
  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {
#pragma omp simd
        for (int n=0; n<=nu_; ++n) {
          (*var_cc)(k,ju+j,i,n) = (*var_cc)(k,ju-j+1,i,n);
        }
      }
    }
  }
  return;
}
