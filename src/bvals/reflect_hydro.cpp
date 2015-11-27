//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file reflect_hydro.cpp
//  \brief implements reflecting BCs in each dimension for primitive hydro variables
//======================================================================================

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh.hpp"

// this class header
#include "bvals.hpp"

//--------------------------------------------------------------------------------------
//! \fn void ReflectInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief  REFLECTING boundary conditions primitive vars, inner x1 boundary (ix1_bc=1)

void ReflectInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                    int is, int ie, int js, int je, int ks, int ke)
{
  for (int n=0; n<(NHYDRO); ++n) {
    if (n==(IVX)) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma simd
          for (int i=1; i<=(NGHOST); ++i)
            a(IVX,k,j,is-i) = -a(IVX,k,j,(is+i-1));  // reflect 1-velocity
        }
      }
    } else {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma simd
          for (int i=1; i<=(NGHOST); ++i)
            a(n,k,j,is-i) = a(n,k,j,(is+i-1));
        }
      }
    }
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void ReflectOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief  REFLECTING boundary conditions primitive vars, outer x1 boundary (ox1_bc=1)

void ReflectOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                    int is, int ie, int js, int je, int ks, int ke)
{
  for (int n=0; n<(NHYDRO); ++n) {
    if (n==(IVX)) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma simd
          for (int i=1; i<=(NGHOST); ++i)
            a(IVX,k,j,ie+i) = -a(IVX,k,j,(ie-i+1));  // reflect 1-velocity
        }
      }
    }
    else {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma simd
          for (int i=1; i<=(NGHOST); ++i)
          a(n,k,j,ie+i) = a(n,k,j,(ie-i+1));
        }
      }
    }
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void ReflecInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief  REFLECTING boundary conditions primitive vars, inner x2 boundary (ix2_bc=1)

void ReflectInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                    int is, int ie, int js, int je, int ks, int ke)
{
  for (int n=0; n<(NHYDRO); ++n) {
    if (n==(IVY)) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
          for (int i=is; i<=ie; ++i)
            a(IVY,k,js-j,i) = -a(IVY,k,js+j-1,i);  // reflect 2-velocity
        }
      }
    }
    else {
#pragma simd
      for (int k=ks; k<=ke; ++k) {
        for (int j=1; j<=(NGHOST); ++j) {
          for (int i=is; i<=ie; ++i)
          a(n,k,js-j,i) = a(n,k,js+j-1,i);
        }
      }

    }
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void ReflectOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief  REFLECTING boundary conditions primitive vars, outer x2 boundary (ox2_bc=1)

void ReflectOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                    int is, int ie, int js, int je, int ks, int ke)
{
  for (int n=0; n<(NHYDRO); ++n) {
    if (n==(IVY)) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
          for (int i=is; i<=ie; ++i)
            a(IVY,k,je+j,i) = -a(IVY,k,je-j+1,i);  // reflect 2-velocity
        }
      }
    }
    else {
      for (int k=ks; k<=ke; ++k) {
        for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
          for (int i=is; i<=ie; ++i)
            a(n,k,je+j,i) = a(n,k,je-j+1,i);
        }
      }
    }
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void ReflectInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief  REFLECTING boundary conditions primitive vars, inner x3 boundary (ix3_bc=1)

void ReflectInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                    int is, int ie, int js, int je, int ks, int ke)
{
  for (int n=0; n<(NHYDRO); ++n) {
    if (n==(IVZ)) {
      for (int k=1; k<=(NGHOST); ++k) {
        for (int j=js; j<=je; ++j) {
#pragma simd
          for (int i=is; i<=ie; ++i)
            a(IVZ,ks-k,j,i) = -a(IVZ,ks+k-1,j,i);  // reflect 3-velocity
        }
      }
    }
    else {
      for (int k=1; k<=(NGHOST); ++k) {
        for (int j=js; j<=je; ++j) {
#pragma simd
          for (int i=is; i<=ie; ++i)
            a(n,ks-k,j,i) = a(n,ks+k-1,j,i);
        }
      }
    }
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void ReflectOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief  REFLECTING boundary conditions primitive vars, outer x3 boundary (ox3_bc=1)

void ReflectOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                    int is, int ie, int js, int je, int ks, int ke)
{
  for (int n=0; n<(NHYDRO); ++n) {
    if (n==(IVZ)) {
      for (int k=1; k<=(NGHOST); ++k) {
        for (int j=js; j<=je; ++j) {
#pragma simd
          for (int i=is; i<=ie; ++i)
            a(IVZ,ke+k,j,i) = -a(IVZ,ke-k+1,j,i);  // reflect 3-velocity
        }
      }
    }
    else {
#pragma simd
      for (int k=1; k<=(NGHOST); ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i)
            a(n,ke+k,j,i) = a(n,ke-k+1,j,i);
        }
      }
    }
  }

  return;
}
