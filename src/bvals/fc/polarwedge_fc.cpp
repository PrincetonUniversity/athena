//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file polarwedge_fc.cpp
//! \brief implementation of polar wedge BCs in x2 direction

// C headers

// C++ headers

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "bvals_fc.hpp"

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::PolarWedgeInnerX2(
//!              Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh)
//! \brief polar wedge boundary conditions, inner x2 boundary

void FaceCenteredBoundaryVariable::PolarWedgeInnerX2(
    Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh) {
  // copy face-centered magnetic fields into ghost zones, reflecting b2
  Real sign = flip_across_pole_field[IB1] ? -1.0 : 1.0;
  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu+1; ++i) {
        (*var_fc).x1f(k,(jl-j),i) = sign * (*var_fc).x1f(k,(jl+j-1),i);
      }
    }
  }
  sign = flip_across_pole_field[IB2] ? -1.0 : 1.0;
  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        (*var_fc).x2f(k,(jl-j),i) = sign * (*var_fc).x2f(k,(jl+j  ),i);
      }
    }
  }
  for (int k=kl; k<=ku; ++k) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      (*var_fc).x2f(k,jl,i) = 0.0;
    }
  }
  sign = flip_across_pole_field[IB3] ? -1.0 : 1.0;
  for (int k=kl; k<=ku+1; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        (*var_fc).x3f(k,(jl-j),i) = sign * (*var_fc).x3f(k,(jl+j-1),i);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::PolarWedgeOuterX2(
//!              Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh)
//! \brief polar wedge boundary conditions, outer x2 boundary

void FaceCenteredBoundaryVariable::PolarWedgeOuterX2(
    Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh) {
  // copy face-centered magnetic fields into ghost zones, reflecting b2
  Real sign = flip_across_pole_field[IB1] ? -1.0 : 1.0;
  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu+1; ++i) {
        (*var_fc).x1f(k,(ju+j  ),i) = sign * (*var_fc).x1f(k,(ju-j+1),i);
      }
    }
  }
  sign = flip_across_pole_field[IB2] ? -1.0 : 1.0;
  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        (*var_fc).x2f(k,(ju+j+1),i) = sign * (*var_fc).x2f(k,(ju-j+1),i);
      }
    }
  }
  for (int k=kl; k<=ku; ++k) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      (*var_fc).x2f(k,(ju+1),i) = 0.0;
    }
  }
  sign = flip_across_pole_field[IB3] ? -1.0 : 1.0;
  for (int k=kl; k<=ku+1; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        (*var_fc).x3f(k,(ju+j  ),i) =  sign * (*var_fc).x3f(k,(ju-j+1),i);
      }
    }
  }
  return;
}
