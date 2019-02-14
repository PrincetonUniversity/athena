//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file polarwedge.cpp
//  \brief implementation of polar wedge BCs in x2 direction

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "bvals.hpp"

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::PolarWedgeInnerX2(
//                         FaceField &b, const Real time, const Real dt,
//                         int il, int iu, int jl, int ju, int kl, int ku, int ngh)
//  \brief polar wedge boundary conditions, inner x2 boundary

void CellCenteredBoundaryVariable::PolarWedgeInnerX2(
                       FaceField &b, Real time, Real dt,
                       int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // copy hydro variables into ghost zones, reflecting v2
  for (int n=0; n<(NHYDRO); ++n) {
    Real sign = flip_across_pole_hydro[n] ? -1.0 : 1.0;
    for (int k=kl; k<=ku; ++k) {
      for (int j=nl; j<=nu; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          prim(n,k,jl-j,i) = sign * prim(n,k,jl+j-1,i);
        }
      }
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b2
  if (MAGNETIC_FIELDS_ENABLED) {
    Real sign = flip_across_pole_field[IB1] ? -1.0 : 1.0;
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu+1; ++i) {
          b.x1f(k,(jl-j),i) = sign * b.x1f(k,(jl+j-1),i);
        }
      }
    }

    sign = flip_across_pole_field[IB2] ? -1.0 : 1.0;
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x2f(k,(jl-j),i) = sign * b.x2f(k,(jl+j  ),i);
        }
      }
    }
    for (int k=kl; k<=ku; ++k) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        b.x2f(k,jl,i) = 0.0;
      }
    }

    sign = flip_across_pole_field[IB3] ? -1.0 : 1.0;
    for (int k=kl; k<=ku+1; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x3f(k,(jl-j),i) = sign * b.x3f(k,(jl+j-1),i);
        }
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::PolarWedgeOuterX2(
//                        MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                        FaceField &b, const Real time, const Real dt,
//                        int il, int iu, int jl, int ju, int kl, int ku, int ngh)
//  \brief polar wedge boundary conditions, outer x2 boundary

void CellCenteredBoundaryVariable::PolarWedgeOuterX2(
    MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                       FaceField &b, Real time, Real dt,
                       int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // copy hydro variables into ghost zones, reflecting v2
  for (int n=0; n<(NHYDRO); ++n) {
    Real sign = flip_across_pole_hydro[n] ? -1.0 : 1.0;
    for (int k=kl; k<=ku; ++k) {
      for (int j=nl; j<=nu; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          prim(n,k,ju+j,i) = sign * prim(n,k,ju-j+1,i);
        }
      }
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b2
  if (MAGNETIC_FIELDS_ENABLED) {
    Real sign = flip_across_pole_field[IB1] ? -1.0 : 1.0;
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu+1; ++i) {
          b.x1f(k,(ju+j  ),i) = sign * b.x1f(k,(ju-j+1),i);
        }
      }
    }

    sign = flip_across_pole_field[IB2] ? -1.0 : 1.0;
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x2f(k,(ju+j+1),i) = sign * b.x2f(k,(ju-j+1),i);
        }
      }
    }
    for (int k=kl; k<=ku; ++k) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        b.x2f(k,(ju+1),i) = 0.0;
      }
    }


    sign = flip_across_pole_field[IB3] ? -1.0 : 1.0;
    for (int k=kl; k<=ku+1; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x3f(k,(ju+j  ),i) =  sign * b.x3f(k,(ju-j+1),i);
        }
      }
    }
  }

  return;
}
