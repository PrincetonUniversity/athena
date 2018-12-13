//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file outflow.cpp
//  \brief implementation of outflow BCs in each dimension for cell-centered AthenaArray

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "bvals.hpp"

//----------------------------------------------------------------------------------------
//! \fn void OutflowInnerX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
//                          AthenaArray<Real> &arr_cc, int il, int iu, int jl, int ju,
//                          int kl, int ku, int nl, int nu) {
//  \brief OUTFLOW boundary conditions, inner x1 boundary

void OutflowInnerX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                    AthenaArray<Real> &arr_cc, int il, int iu, int jl, int ju,
                    int kl, int ku, int nl, int nu) {
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=nl; i<=nu; ++i) {
        arr_cc(n,k,j,il-i) = arr_cc(n,k,j,il);
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
//                          AthenaArray<Real> &arr_cc, int il, int iu, int jl, int ju,
//                          int kl, int ku, int nl, int nu) {
//  \brief OUTFLOW boundary conditions, outer x1 boundary

void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                    AthenaArray<Real> &arr_cc, int il, int iu, int jl, int ju,
                    int kl, int ku, int nl, int nu) {
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=nl; i<=nu; ++i) {
        arr_cc(n,k,j,iu+i) = arr_cc(n,k,j,iu);
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OutflowInnerX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
//                          AthenaArray<Real> &arr_cc, int il, int iu, int jl, int ju,
//                          int kl, int ku, int nl, int nu) {
//  \brief OUTFLOW boundary conditions, inner x2 boundary

void OutflowInnerX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                    AthenaArray<Real> &arr_cc, int il, int iu, int jl, int ju,
                    int kl, int ku, int nl, int nu) {
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=kl; k<=ku; ++k) {
    for (int j=nl; j<=nu; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        arr_cc(n,k,jl-j,i) = arr_cc(n,k,jl,i);
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OutflowOuterX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
//                          AthenaArray<Real> &arr_cc, int il, int iu, int jl, int ju,
//                          int kl, int ku, int nl, int nu) {
//  \brief OUTFLOW boundary conditions, outer x2 boundary

void OutflowOuterX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                    AthenaArray<Real> &arr_cc, int il, int iu, int jl, int ju,
                    int kl, int ku, int nl, int nu) {
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=kl; k<=ku; ++k) {
    for (int j=nl; j<=nu; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        arr_cc(n,k,ju+j,i) = arr_cc(n,k,ju,i);
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OutflowInnerX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
//                          AthenaArray<Real> &arr_cc, int il, int iu, int jl, int ju,
//                          int kl, int ku, int nl, int nu) {
//  \brief OUTFLOW boundary conditions, inner x3 boundary

void OutflowInnerX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                    AthenaArray<Real> &arr_cc, int il, int iu, int jl, int ju,
                    int kl, int ku, int nl, int nu) {
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=nl; k<=nu; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        arr_cc(n,kl-k,j,i) = arr_cc(n,kl,j,i);
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OutflowOuterX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
//                          AthenaArray<Real> &arr_cc, int il, int iu, int jl, int ju,
//                          int kl, int ku, int nl, int nu) {
//  \brief OUTFLOW boundary conditions, outer x3 boundary

void OutflowOuterX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                    AthenaArray<Real> &arr_cc, int il, int iu, int jl, int ju,
                    int kl, int ku, int nl, int nu) {
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=nl; k<=nu; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        arr_cc(n,ku+k,j,i) = arr_cc(n,ku,j,i);
      }
    }}
  }

  return;
}
