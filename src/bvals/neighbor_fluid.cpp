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

// *** temporary for test, replace this adequately for MPI ***
int myrank=0;

//======================================================================================
//! \file Neighbor_fluid.cpp
//  \brief implements Neighbor BCs in each dimension for conserved fluid variables
//======================================================================================
//--------------------------------------------------------------------------------------
//! \fn void NeighborInnerX1(MeshBlock *pmb)
//  \brief  Neighbor boundary conditions conserved vars, inner x1 boundary (ix1_bc=4)

void NeighborInnerX1(MeshBlock *pmb, AthenaArray<Real> &a)
{
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;
  
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int n=0; n<(NFLUID); ++n) {
#pragma simd
        for (int i=1; i<=(NGHOST); ++i) {
          a(n,k,j,is-i) = a(n,k,j,(ie-i+1));
        }
      }
    }
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void NeighborOuterX1(MeshBlock *pmb)
//  \brief  Neighbor boundary conditions conserved vars, outer x1 boundary (ox1_bc=4)

void NeighborOuterX1(MeshBlock *pmb, AthenaArray<Real> &a)
{
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for (int n=0; n<(NFLUID); ++n) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        a(n,k,j,ie+i) = a(n,k,j,(is+i-1));
      }
    }
  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void NeighborInnerX2(MeshBlock *pmb)
//  \brief  Neighbor boundary conditions conserved vars, inner x2 boundary (ix2_bc=4)

void NeighborInnerX2(MeshBlock *pmb, AthenaArray<Real> &a)
{
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  for (int k=ks; k<=ke; ++k) {
  for (int j=1; j<=(NGHOST); ++j) {
    for (int n=0; n<(NFLUID); ++n) {
#pragma simd
      for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
        a(n,k,js-j,i) = a(n,k,je-j+1,i);
      }
    }
  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void NeighborOuterX2(MeshBlock *pmb)
//  \brief  Neighbor boundary conditions conserved vars, outer x2 boundary (ox2_bc=4)

void NeighborOuterX2(MeshBlock *pmb, AthenaArray<Real> &a)
{
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  for (int k=ks; k<=ke; ++k) {
  for (int j=1; j<=(NGHOST); ++j) {
    for (int n=0; n<(NFLUID); ++n) {
#pragma simd
      for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
        a(n,k,je+j,i) = a(n,k,js+j-1,i);
      }
    }
  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void NeighborInnerX3(MeshBlock *pmb)
//  \brief  Neighbor boundary conditions conserved vars, inner x3 boundary (ix3_bc=4)

void NeighborInnerX3(MeshBlock *pmb, AthenaArray<Real> &a)
{
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  for (int k=1; k<=(NGHOST); ++k) {
  for (int j=js-(NGHOST); j<=je+(NGHOST); ++j) {
    for (int n=0; n<(NFLUID); ++n) {
#pragma simd
      for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
        a(n,ks-k,j,i) = a(n,ke-k+1,j,i);
      }
    }
  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void NeighborOuterX3(MeshBlock *pmb)
//  \brief  Neighbor boundary conditions conserved vars, outer x3 boundary (ox3_bc=4)

void NeighborOuterX3(MeshBlock *pmb, AthenaArray<Real> &a)
{
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  for (int k=1; k<=(NGHOST); ++k) {
  for (int j=js-(NGHOST); j<=je+(NGHOST); ++j) {
    for (int n=0; n<(NFLUID); ++n) {
#pragma simd
      for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
        a(n,ke+k,j,i) = a(n,ks+k-1,j,i);
      }
    }
  }}

  return;
}
