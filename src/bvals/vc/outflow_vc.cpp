//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file outflow_vc.cpp
//  \brief implementation of outflow BCs in each dimension for vertex-centered AthenaArray
//
//  Notes:
//  - data is copied from extremal [boundary] vertices out to ghosts

// C headers

// C++ headers

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "bvals_vc.hpp"

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::OutflowInnerX1(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief OUTFLOW boundary conditions, inner x1 boundary

void VertexCenteredBoundaryVariable::OutflowInnerX1(
    Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh) {

  for (int n=0; n<=nu_; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          (*var_vc)(n,k,j,il-i) = (*var_vc)(n,k,j,il);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::OutflowOuterX1(
//          Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh)
//  \brief OUTFLOW boundary conditions, outer x1 boundary

void VertexCenteredBoundaryVariable::OutflowOuterX1(
    Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh) {

  for (int n=0; n<=nu_; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          (*var_vc)(n,k,j,iu+i) = (*var_vc)(n,k,j,iu);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::OutflowInnerX2(
//          Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh)
//  \brief OUTFLOW boundary conditions, inner x2 boundary

void VertexCenteredBoundaryVariable::OutflowInnerX2(
    Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh) {

  for (int n=0; n<=nu_; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          (*var_vc)(n,k,jl-j,i) = (*var_vc)(n,k,jl,i);
        }
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::OutflowOuterX2(
//          Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh)
//  \brief OUTFLOW boundary conditions, outer x2 boundary

void VertexCenteredBoundaryVariable::OutflowOuterX2(
    Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh) {

  for (int n=0; n<=nu_; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          (*var_vc)(n,k,ju+j,i) = (*var_vc)(n,k,ju,i);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::OutflowInnerX3(
//          Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ngh)
//  \brief OUTFLOW boundary conditions, inner x3 boundary

void VertexCenteredBoundaryVariable::OutflowInnerX3(
    Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ngh) {

  for (int n=0; n<=nu_; ++n) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          (*var_vc)(n,kl-k,j,i) = (*var_vc)(n,kl,j,i);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::OutflowOuterX3(
//          Real time, Real dt, int il, int iu, int jl, int ju, int ku, int ngh)
//  \brief OUTFLOW boundary conditions, outer x3 boundary

void VertexCenteredBoundaryVariable::OutflowOuterX3(
    Real time, Real dt, int il, int iu, int jl, int ju, int ku, int ngh) {

  for (int n=0; n<=nu_; ++n) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          (*var_vc)(n,ku+k,j,i) = (*var_vc)(n,ku,j,i);
        }
      }
    }
  }
  return;
}
