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

// Primary header
#include "bvals.hpp"

// Athena headers
#include "../athena.hpp"         // macros, Real
#include "../athena_arrays.hpp"  // AthenaArray
#include "../mesh.hpp"           // MeshBlock

//======================================================================================
//! \file reflect_fluid.cpp
//  \brief implements reflecting BCs in each dimension for conserved fluid variables
//======================================================================================
//--------------------------------------------------------------------------------------
//! \fn void ReflectInnerX1(MeshBlock *pmb, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief  REFLECTING boundary conditions conserved vars, inner x1 boundary (ix1_bc=1)

void ReflectInnerX1(MeshBlock *pmb, AthenaArray<Real> &a,
                    int is, int ie, int js, int je, int ks, int ke)
{
  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for (int n=0; n<(NFLUID); ++n) {

      if (n==(IM1)) {
#pragma simd
        for (int i=1; i<=(NGHOST); ++i) {
          a(IM1,k,j,is-i) = -a(IM1,k,j,(is+i-1));  // reflect 1-mom
        }

      } else {
#pragma simd
        for (int i=1; i<=(NGHOST); ++i) {
          a(n,k,j,is-i) = a(n,k,j,(is+i-1));
        }
      }

    }
  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void ReflectOuterX1(MeshBlock *pmb, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief  REFLECTING boundary conditions conserved vars, outer x1 boundary (ox1_bc=1)

void ReflectOuterX1(MeshBlock *pmb, AthenaArray<Real> &a,
                    int is, int ie, int js, int je, int ks, int ke)
{
  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for (int n=0; n<(NFLUID); ++n) {

      if (n==(IM1)) {
#pragma simd
        for (int i=1; i<=(NGHOST); ++i) {
          a(IM1,k,j,ie+i) = -a(IM1,k,j,(ie-i+1));  // reflect 1-mom
        }

      } else {
#pragma simd
        for (int i=1; i<=(NGHOST); ++i) {
          a(n,k,j,ie+i) = a(n,k,j,(ie-i+1));
        }
      }

    }
  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void ReflecInnerX2(MeshBlock *pmb, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief  REFLECTING boundary conditions conserved vars, inner x2 boundary (ix2_bc=1)

void ReflectInnerX2(MeshBlock *pmb, AthenaArray<Real> &a,
                    int is, int ie, int js, int je, int ks, int ke)
{
  for (int k=ks; k<=ke; ++k) {
  for (int j=1; j<=(NGHOST); ++j) {
    for (int n=0; n<(NFLUID); ++n) {

      if (n==(IM2)) {
#pragma simd
        for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
          a(IM2,k,js-j,i) = -a(IM2,k,js+j-1,i);  // reflect 2-mom
        }

      } else {
#pragma simd
        for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
          a(n,k,js-j,i) = a(n,k,js+j-1,i);
        }
      }

    }
  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void ReflectOuterX2(MeshBlock *pmb, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief  REFLECTING boundary conditions conserved vars, outer x2 boundary (ox2_bc=1)

void ReflectOuterX2(MeshBlock *pmb, AthenaArray<Real> &a,
                    int is, int ie, int js, int je, int ks, int ke)
{
  for (int k=ks; k<=ke; ++k) {
  for (int j=1; j<=(NGHOST); ++j) {
    for (int n=0; n<(NFLUID); ++n) {

      if (n==(IM2)) {
#pragma simd
        for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
          a(IM2,k,je+j,i) = -a(IM2,k,je-j+1,i);  // reflect 2-mom
        }

      } else {
#pragma simd
        for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
          a(n,k,je+j,i) = a(n,k,je-j+1,i);
        }
      }

    }
  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void ReflectInnerX3(MeshBlock *pmb)
//  \brief  REFLECTING boundary conditions conserved vars, inner x3 boundary (ix3_bc=1)

void ReflectInnerX3(MeshBlock *pmb, AthenaArray<Real> &a,
                    int is, int ie, int js, int je, int ks, int ke)
{
  for (int k=1; k<=(NGHOST); ++k) {
  for (int j=js-(NGHOST); j<=je+(NGHOST); ++j) {
    for (int n=0; n<(NFLUID); ++n) {

      if (n==(IM3)) {
#pragma simd
        for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
          a(IM3,ks-k,j,i) = -a(IM3,ks+k-1,j,i);  // reflect 3-mom
        }

      } else {
#pragma simd
        for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
          a(n,ks-k,j,i) = a(n,ks+k-1,j,i);
        }
      }

    }
  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void ReflectOuterX3(MeshBlock *pmb, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief  REFLECTING boundary conditions conserved vars, outer x3 boundary (ox3_bc=1)

void ReflectOuterX3(MeshBlock *pmb, AthenaArray<Real> &a,
                    int is, int ie, int js, int je, int ks, int ke)
{
  for (int k=1; k<=(NGHOST); ++k) {
  for (int j=js-(NGHOST); j<=je+(NGHOST); ++j) {
    for (int n=0; n<(NFLUID); ++n) {

      if (n==(IM3)) {
#pragma simd
        for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
          a(IM3,ke+k,j,i) = -a(IM3,ke-k+1,j,i);  // reflect 3-mom
        }

      } else {
#pragma simd
        for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
          a(n,ke+k,j,i) = a(n,ke-k+1,j,i);
        }
      }

    }
  }}

  return;
}
