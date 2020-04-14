//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file extrapolate_outflow_vc.cpp
//  \brief implementation of extrapolated outflow BCs in each dimension for
//         vertex-centered AthenaArray
//
//  Notes:
//  - data is copied from extremal [boundary] vertices out to ghosts
//  - OuterXn conditions expect the same input idx arguments as outflow_cc but
//    are adjusted to work with VC

// C, C++ headers
#include <iostream>

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "bvals_vc.hpp"

// BD: TODO - extrapolation order should probably be modified based on ghosts

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::ExtrapolateOutflowInnerX1(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief OUTFLOW boundary conditions with extrapolation, inner x1 boundary

void VertexCenteredBoundaryVariable::ExtrapolateOutflowInnerX1(
    Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh) {

  ju = IncrementIfNonzero(ju);
  ku = IncrementIfNonzero(ku);

  for (int n=0; n<=nu_; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i = il-1; i >= il-ngh; --i) {
          // extrapolate variables at 4th order
          (*var_vc)(n,k,j,i) = 4.*(*var_vc)(n,k,j,i+1) - 6.*(*var_vc)(n,k,j,i+2) +
                               4.*(*var_vc)(n,k,j,i+3) - 1.*(*var_vc)(n,k,j,i+4);

          // extrapolate variables at 6th order
          // (*var_vc)(n,k,j,i) = 6.*(*var_vc)(n,k,j,i+1) - 15.*(*var_vc)(n,k,j,i+2) +
          //                     20.*(*var_vc)(n,k,j,i+3) - 15.*(*var_vc)(n,k,j,i+4) +
          //                      6.*(*var_vc)(n,k,j,i+5) - 1.*(*var_vc)(n,k,j,i+6);

        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::ExtrapolateOutflowOuterX1(
//          Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh)
//  \brief OUTFLOW boundary conditions with extrapolation, outer x1 boundary

void VertexCenteredBoundaryVariable::ExtrapolateOutflowOuterX1(
    Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh) {

  iu = IncrementIfNonzero(iu);
  ju = IncrementIfNonzero(ju);
  ku = IncrementIfNonzero(ku);

  for (int n=0; n<=nu_; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=iu+1; i<=iu+ngh; ++i) {
          // extrapolate variables at 4th order
          (*var_vc)(n,k,j,i) = 4.*(*var_vc)(n,k,j,i-1) - 6.*(*var_vc)(n,k,j,i-2) +
                               4.*(*var_vc)(n,k,j,i-3) - 1.*(*var_vc)(n,k,j,i-4);

          // extrapolate variables at 6th order
          // (*var_vc)(n,k,j,i) = 6.*(*var_vc)(n,k,j,i-1) - 15.*(*var_vc)(n,k,j,i-2) +
          //                     20.*(*var_vc)(n,k,j,i-3) - 15.*(*var_vc)(n,k,j,i-4) +
          //                      6.*(*var_vc)(n,k,j,i-5) - 1.*(*var_vc)(n,k,j,i-6);
        }

      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::ExtrapolateOutflowInnerX2(
//          Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh)
//  \brief OUTFLOW boundary conditions with extrapolation, inner x2 boundary

void VertexCenteredBoundaryVariable::ExtrapolateOutflowInnerX2(
    Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh) {

  iu = IncrementIfNonzero(iu);
  ku = IncrementIfNonzero(ku);

  for (int n=0; n<=nu_; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl-1; j>=jl-ngh; --j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          // extrapolate variables at 4th order
          (*var_vc)(n,k,j,i) = 4.*(*var_vc)(n,k,j+1,i) - 6.*(*var_vc)(n,k,j+2,i) +
                               4.*(*var_vc)(n,k,j+3,i) - 1.*(*var_vc)(n,k,j+4,i);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::ExtrapolateOutflowOuterX2(
//          Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh)
//  \brief OUTFLOW boundary conditions with extrapolation, outer x2 boundary

void VertexCenteredBoundaryVariable::ExtrapolateOutflowOuterX2(
    Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh) {

  iu = IncrementIfNonzero(iu);
  ju = IncrementIfNonzero(ju);
  ku = IncrementIfNonzero(ku);

  for (int n=0; n<=nu_; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=ju+1; j<=ju+ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          // extrapolate variables at 4th order
          (*var_vc)(n,k,j,i) = 4.*(*var_vc)(n,k,j-1,i) - 6.*(*var_vc)(n,k,j-2,i) +
                               4.*(*var_vc)(n,k,j-3,i) - 1.*(*var_vc)(n,k,j-4,i);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::ExtrapolateOutflowInnerX3(
//          Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ngh)
//  \brief OUTFLOW boundary conditions with extrapolation, inner x3 boundary

void VertexCenteredBoundaryVariable::ExtrapolateOutflowInnerX3(
    Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ngh) {

  iu = IncrementIfNonzero(iu);
  ju = IncrementIfNonzero(ju);

  for (int n=0; n<=nu_; ++n) {
    for (int k=kl-1; k>=kl-ngh; --k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          // extrapolate variables at 4th order
          (*var_vc)(n,k,j,i) = 4.*(*var_vc)(n,k+1,j,i) - 6.*(*var_vc)(n,k+2,j,i) +
                               4.*(*var_vc)(n,k+3,j,i) - 1.*(*var_vc)(n,k+4,j,i);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::ExtrapolateOutflowOuterX3(
//          Real time, Real dt, int il, int iu, int jl, int ju, int ku, int ngh)
//  \brief OUTFLOW boundary conditions with extrapolation, outer x3 boundary

void VertexCenteredBoundaryVariable::ExtrapolateOutflowOuterX3(
    Real time, Real dt, int il, int iu, int jl, int ju, int ku, int ngh) {

  iu = IncrementIfNonzero(iu);
  ju = IncrementIfNonzero(ju);
  ku = IncrementIfNonzero(ku);

  for (int n=0; n<=nu_; ++n) {
    for (int k=ku+1; k<=ku+ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          // extrapolate variables at 4th order
          (*var_vc)(n,k,j,i) = 4.*(*var_vc)(n,k-1,j,i) - 6.*(*var_vc)(n,k-2,j,i) +
                               4.*(*var_vc)(n,k-3,j,i) - 1.*(*var_vc)(n,k-4,j,i);
        }
      }
    }
  }
  return;
}
