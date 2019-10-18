//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file outflow.cpp
//  \brief implementation of extrapolated outflow BCs in each dimension for cell-centered
//  AthenaArray

// C, C++ headers
#include <iostream>

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "bvals_cc.hpp"

// BD: TODO - extrapolation order should probably be modified based on ghosts

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::ExtrapolateOutflowInnerX1(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief OUTFLOW boundary conditions with extrapolation, inner x1 boundary

void CellCenteredBoundaryVariable::ExtrapolateOutflowInnerX1(
    Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh) {
  for (int n=0; n<=nu_; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i = il-1; i >= il-ngh; --i) {
          // extrapolate variables at 4th order
          (*var_cc)(n,k,j,i) = 4.*(*var_cc)(n,k,j,i+1) - 6.*(*var_cc)(n,k,j,i+2) +
                               4.*(*var_cc)(n,k,j,i+3) - 1.*(*var_cc)(n,k,j,i+4);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::ExtrapolateOutflowOuterX1(
//          Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh)
//  \brief OUTFLOW boundary conditions with extrapolation, outer x1 boundary

void CellCenteredBoundaryVariable::ExtrapolateOutflowOuterX1(
    Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int n=0; n<=nu_; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=iu+1; i<=iu+ngh; ++i) {
          // extrapolate variables at 4th order
          (*var_cc)(n,k,j,i) = 4.*(*var_cc)(n,k,j,i-1) - 6.*(*var_cc)(n,k,j,i-2) +
                               4.*(*var_cc)(n,k,j,i-3) - 1.*(*var_cc)(n,k,j,i-4);
        }

      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::ExtrapolateOutflowInnerX2(
//          Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh)
//  \brief OUTFLOW boundary conditions with extrapolation, inner x2 boundary

void CellCenteredBoundaryVariable::ExtrapolateOutflowInnerX2(
    Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh) {
  for (int n=0; n<=nu_; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl-1; j>=jl-ngh; --j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          // extrapolate variables at 4th order
          (*var_cc)(n,k,j,i) = 4.*(*var_cc)(n,k,j+1,i) - 6.*(*var_cc)(n,k,j+2,i) +
                               4.*(*var_cc)(n,k,j+3,i) - 1.*(*var_cc)(n,k,j+4,i);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::ExtrapolateOutflowOuterX2(
//          Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh)
//  \brief OUTFLOW boundary conditions with extrapolation, outer x2 boundary

void CellCenteredBoundaryVariable::ExtrapolateOutflowOuterX2(
    Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh) {
  for (int n=0; n<=nu_; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=ju+1; j<=ju+ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          // extrapolate variables at 4th order
          (*var_cc)(n,k,j,i) = 4.*(*var_cc)(n,k,j-1,i) - 6.*(*var_cc)(n,k,j-2,i) +
                               4.*(*var_cc)(n,k,j-3,i) - 1.*(*var_cc)(n,k,j-4,i);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::ExtrapolateOutflowInnerX3(
//          Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ngh)
//  \brief OUTFLOW boundary conditions with extrapolation, inner x3 boundary

void CellCenteredBoundaryVariable::ExtrapolateOutflowInnerX3(
    Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ngh) {
  for (int n=0; n<=nu_; ++n) {
    for (int k=kl-1; k>=kl-ngh; --k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          // extrapolate variables at 4th order
          (*var_cc)(n,k,j,i) = 4.*(*var_cc)(n,k+1,j,i) - 6.*(*var_cc)(n,k+2,j,i) +
                               4.*(*var_cc)(n,k+3,j,i) - 1.*(*var_cc)(n,k+4,j,i);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::ExtrapolateOutflowOuterX3(
//          Real time, Real dt, int il, int iu, int jl, int ju, int ku, int ngh)
//  \brief OUTFLOW boundary conditions with extrapolation, outer x3 boundary

void CellCenteredBoundaryVariable::ExtrapolateOutflowOuterX3(
    Real time, Real dt, int il, int iu, int jl, int ju, int ku, int ngh) {
  for (int n=0; n<=nu_; ++n) {
    for (int k=ku+1; k<=ku+ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          // extrapolate variables at 4th order
          (*var_cc)(n,k,j,i) = 4.*(*var_cc)(n,k-1,j,i) - 6.*(*var_cc)(n,k-2,j,i) +
                               4.*(*var_cc)(n,k-3,j,i) - 1.*(*var_cc)(n,k-4,j,i);
        }
      }
    }
  }
  return;
}
