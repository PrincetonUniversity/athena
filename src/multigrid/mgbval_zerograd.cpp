//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mgbval_zerograd.cpp
//! \brief 6x zero gradient (outflow) boundary functions for Multigrid

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../defs.hpp"
#include "multigrid.hpp"


//----------------------------------------------------------------------------------------
//! \fn void::MGZeroGradientInnerX1(AthenaArray<Real> &dst, Real time,
//!                 int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//!                 const MGCoordinates &coord)
//! \brief Zero gradient boundary condition in the inner-X1 direction

void MGZeroGradientInnerX1(AthenaArray<Real> &dst, Real time, int nvar,
                           int is, int ie, int js, int je, int ks, int ke, int ngh,
                           const MGCoordinates &coord) {
  for (int n=0; n<nvar; n++) {
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=0; i<ngh; i++)
          dst(n,k,j,is-i-1) = dst(n,k,j,is);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void::MGZeroGradientOuterX1(AthenaArray<Real> &dst, Real time,
//!                 int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//!                 const MGCoordinates &coord)
//! \brief Zero gradient boundary condition in the outer-X1 direction

void MGZeroGradientOuterX1(AthenaArray<Real> &dst, Real time, int nvar,
                           int is, int ie, int js, int je, int ks, int ke, int ngh,
                           const MGCoordinates &coord) {
  for (int n=0; n<nvar; n++) {
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=0; i<ngh; i++)
          dst(n,k,j,ie+i+1) = dst(n,k,j,ie);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void::MGZeroGradientInnerX2(AthenaArray<Real> &dst, Real time,
//!                 int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//!                 const MGCoordinates &coord)
//! \brief Zero gradient boundary condition in the inner-X2 direction

void MGZeroGradientInnerX2(AthenaArray<Real> &dst, Real time, int nvar,
                           int is, int ie, int js, int je, int ks, int ke, int ngh,
                           const MGCoordinates &coord) {
  for (int n=0; n<nvar; n++) {
    for (int k=ks; k<=ke; k++) {
      for (int j=0; j<ngh; j++) {
        for (int i=is; i<=ie; i++)
          dst(n,k,js-j-1,i) = dst(n,k,js,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void::MGZeroGradientOuterX2(AthenaArray<Real> &dst, Real time,
//!                 int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//!                 const MGCoordinates &coord)
//! \brief Zero gradient boundary condition in the outer-X2 direction

void MGZeroGradientOuterX2(AthenaArray<Real> &dst, Real time, int nvar,
                           int is, int ie, int js, int je, int ks, int ke, int ngh,
                           const MGCoordinates &coord) {
  for (int n=0; n<nvar; n++) {
    for (int k=ks; k<=ke; k++) {
      for (int j=0; j<ngh; j++) {
        for (int i=is; i<=ie; i++)
          dst(n,k,je+j+1,i) = dst(n,k,je,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void::MGZeroGradientInnerX3(AthenaArray<Real> &dst, Real time,
//!                 int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//!                 const MGCoordinates &coord)
//! \brief Zero gradient boundary condition in the inner-X3 direction

void MGZeroGradientInnerX3(AthenaArray<Real> &dst, Real time, int nvar,
                           int is, int ie, int js, int je, int ks, int ke, int ngh,
                           const MGCoordinates &coord) {
  for (int n=0; n<nvar; n++) {
    for (int k=0; k<ngh; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++)
          dst(n,ks-k-1,j,i) = dst(n,ks,j,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void::MGZeroGradientOuterX3(AthenaArray<Real> &dst, Real time,
//!                 int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//!                 const MGCoordinates &coord)
//! \brief Zero gradient boundary condition in the outer-X3 direction

void MGZeroGradientOuterX3(AthenaArray<Real> &dst, Real time, int nvar,
                           int is, int ie, int js, int je, int ks, int ke, int ngh,
                           const MGCoordinates &coord) {
  for (int n=0; n<nvar; n++) {
    for (int k=0; k<ngh; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++)
          dst(n,ke+k+1,j,i) = dst(n,ke,j,i);
      }
    }
  }
  return;
}


