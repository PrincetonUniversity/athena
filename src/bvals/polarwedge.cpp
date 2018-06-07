//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file polarwedge.cpp
//  \brief implementation of polar wedge BCs in x2 direction

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "bvals.hpp"

//----------------------------------------------------------------------------------------
//! \fn void PolarWedgeInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                         FaceField &b, const Real time, const Real dt,
//                         int is, int ie, int js, int je, int ks, int ke, int ngh)
//  \brief polar wedge boundary conditions, inner x2 boundary

void PolarWedgeInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
                    int is, int ie, int js, int je, int ks, int ke, int ngh) {
  // copy hydro variables into ghost zones, reflecting v2
  for (int n=0; n<(NHYDRO); ++n) {
    Real sign = flip_across_pole_hydro[n] ? -1.0 : 1.0;
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          prim(n,k,js-j,i) = sign * prim(n,k,js+j-1,i);
        }
      }
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b2
  if (MAGNETIC_FIELDS_ENABLED) {
    Real sign = flip_across_pole_field[IB1] ? -1.0 : 1.0;
    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        b.x1f(k,(js-j),i) = sign * b.x1f(k,(js+j-1),i);
      }
    }}

    sign = flip_across_pole_field[IB2] ? -1.0 : 1.0;
    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x2f(k,(js-j),i) = sign * b.x2f(k,(js+j  ),i);
      }
    }}
    for (int k=ks; k<=ke; ++k) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x2f(k,js,i) = 0.0;
      }
    }

    sign = flip_across_pole_field[IB3] ? -1.0 : 1.0;
    for (int k=ks; k<=ke+1; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x3f(k,(js-j),i) = sign * b.x3f(k,(js+j-1),i);
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void PolarWedgeOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                          FaceField &b, const Real time, const Real dt,
//                          int is, int ie, int js, int je, int ks, int ke, int ngh)
//  \brief polar wedge boundary conditions, outer x2 boundary

void PolarWedgeOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
                    int is, int ie, int js, int je, int ks, int ke, int ngh) {
  // copy hydro variables into ghost zones, reflecting v2
  for (int n=0; n<(NHYDRO); ++n) {
    Real sign = flip_across_pole_hydro[n] ? -1.0 : 1.0;
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          prim(n,k,je+j,i) = sign * prim(n,k,je-j+1,i);
        }
      }
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b2
  if (MAGNETIC_FIELDS_ENABLED) {
    Real sign = flip_across_pole_field[IB1] ? -1.0 : 1.0;
    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        b.x1f(k,(je+j  ),i) = sign * b.x1f(k,(je-j+1),i);
      }
    }}

    sign = flip_across_pole_field[IB2] ? -1.0 : 1.0;
    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x2f(k,(je+j+1),i) = sign * b.x2f(k,(je-j+1),i);
      }
    }}
    for (int k=ks; k<=ke; ++k) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x2f(k,(je+1),i) = 0.0;
      }
    }


    sign = flip_across_pole_field[IB3] ? -1.0 : 1.0;
    for (int k=ks; k<=ke+1; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x3f(k,(je+j  ),i) =  sign * b.x3f(k,(je-j+1),i);
      }
    }}
  }

  return;
}
