//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 *
 * This program is free software: you can redistribute and/or modify it under the terms
 * of the GNU General Public License (GPL) as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of GNU GPL in the file LICENSE included in the code
 * distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

// Primary header
#include "bvals.hpp"

// Athena headers
#include "../../athena.hpp"         // macros, Real
#include "../../athena_arrays.hpp"  // AthenaArray
#include "../../mesh.hpp"           // MeshBlock

//======================================================================================
/*! \file outflow_fluid.cpp
 *  \brief implements outflow BCs in each dimension for conserved fluid variables
 *====================================================================================*/
//--------------------------------------------------------------------------------------
//! \fn void OutflowInnerX1(MeshBlock *pmb)
//  \brief  OUTFLOW  boundary conditions conserved vars, inner x1 boundary (ix1_bc=2)

void OutflowInnerX1(MeshBlock *pmb, AthenaArray<Real> &a)
{
  int is = pmb->is;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;
  AthenaArray<Real> la = a.ShallowCopy();

  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for (int n=0; n<(NVAR); ++n) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        la(n,k,j,is-i) = la(n,k,j,is);
      }
    }
  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void OutflowOuterX1(MeshBlock *pmb)
//  \brief  OUTFLOW  boundary conditions conserved vars, outer x1 boundary (ox1_bc=2)

void OutflowOuterX1(MeshBlock *pmb, AthenaArray<Real> &a)
{
  int ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;
  AthenaArray<Real> la = a.ShallowCopy();

  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for (int n=0; n<(NVAR); ++n) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        la(n,k,j,ie+i) = la(n,k,j,ie);
      }
    }
  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void OutflowInnerX2(MeshBlock *pmb)
//  \brief  OUTFLOW  boundary conditions conserved vars, inner x2 boundary (ix2_bc=2)

void OutflowInnerX2(MeshBlock *pmb, AthenaArray<Real> &a)
{
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js;
  int ks = pmb->ks, ke = pmb->ke;
  AthenaArray<Real> la = a.ShallowCopy();

  for (int k=ks; k<=ke; ++k) {
  for (int j=1; j<=(NGHOST); ++j) {
    for (int n=0; n<(NVAR); ++n) {
#pragma simd
      for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
        la(n,k,js-j,i) = la(n,k,js,i);
      }
    }
  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void OutflowOuterX2(MeshBlock *pmb)
//  \brief  OUTFLOW  boundary conditions conserved vars, outer x2 boundary (ox2_bc=2)

void OutflowOuterX2(MeshBlock *pmb, AthenaArray<Real> &a)
{
  int is = pmb->is, ie = pmb->ie;
  int je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;
  AthenaArray<Real> la = a.ShallowCopy();

  for (int k=ks; k<=ke; ++k) {
  for (int j=1; j<=(NGHOST); ++j) {
    for (int n=0; n<(NVAR); ++n) {
#pragma simd
      for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
        la(n,k,je+j,i) = la(n,k,je,i);
      }
    }
  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void OutflowInnerX3(MeshBlock *pmb)
//  \brief  OUTFLOW  boundary conditions conserved vars, inner x3 boundary (ix3_bc=2)

void OutflowInnerX3(MeshBlock *pmb, AthenaArray<Real> &a)
{
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks;
  AthenaArray<Real> la = a.ShallowCopy();

  for (int k=1; k<=(NGHOST); ++k) {
  for (int j=js-(NGHOST); j<=je+(NGHOST); ++j) {
    for (int n=0; n<(NVAR); ++n) {
#pragma simd
      for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
        la(n,ks-k,j,i) = la(n,ks,j,i);
      }
    }
  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void OutflowOuterX3(MeshBlock *pmb)
//  \brief  OUTFLOW  boundary conditions conserved vars, outer x3 boundary (ox3_bc=2)

void OutflowOuterX3(MeshBlock *pmb, AthenaArray<Real> &a)
{
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ke = pmb->ke;
  AthenaArray<Real> la = a.ShallowCopy();

  for (int k=1; k<=(NGHOST); ++k) {
  for (int j=js-(NGHOST); j<=je+(NGHOST); ++j) {
    for (int n=0; n<(NVAR); ++n) {
#pragma simd
      for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
        la(n,ke+k,j,i) = la(n,ke,j,i);
      }
    }
  }}

  return;
}
