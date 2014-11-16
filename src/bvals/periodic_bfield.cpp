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
#include "../field/field.hpp"          // InterfaceField

//======================================================================================
//! \file periodic_bfield.cpp
//  \brief implements periodic BCs in each dimension for interface B-field
//======================================================================================
//--------------------------------------------------------------------------------------
//! \fn void PeriodicInnerX1()
//  \brief  PERIODIC boundary conditions interface B, inner x1 boundary (ix1_bc=1)

void PeriodicInnerX1(MeshBlock *pmb, InterfaceField &a)
{
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {

#pragma simd
    for (int i=1; i<=(NGHOST); ++i) {
      a.x1f(k,j,(is-i)) = a.x1f(k,j,(ie-i+1));
      a.x2f(k,j,(is-i)) = a.x2f(k,j,(ie-i+1));
      a.x3f(k,j,(is-i)) = a.x3f(k,j,(ie-i+1));
    }

  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void PeriodicOuterX1()
//  \brief  PERIODIC boundary conditions interface B, outer x1 boundary (ox1_bc=1)

void PeriodicOuterX1(MeshBlock *pmb, InterfaceField &a)
{
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {

#pragma simd
    for (int i=1; i<=(NGHOST); ++i) {
      a.x1f(k,j,(ie+i+1)) = a.x1f(k,j,(is+i  ));
      a.x2f(k,j,(ie+i  )) = a.x2f(k,j,(is+i-1));
      a.x3f(k,j,(ie+i  )) = a.x3f(k,j,(is+i-1));
    }

  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void PeriodicInnerX2()
//  \brief  PERIODIC boundary conditions interface B, inner x2 boundary (ix2_bc=1)

void PeriodicInnerX2(MeshBlock *pmb, InterfaceField &a)
{
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  for (int k=ks; k<=ke; ++k) {
  for (int j=1; j<=(NGHOST); ++j) {

#pragma simd
    for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
      a.x1f(k,(js-j),i) = a.x1f(k,(je-j+1),i);
      a.x2f(k,(js-j),i) = a.x2f(k,(je-j+1),i);
      a.x3f(k,(js-j),i) = a.x3f(k,(je-j+1),i);
    }

  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void PeriodicOuterX2()
//  \brief  PERIODIC boundary conditions interface B, outer x2 boundary (ox2_bc=1)

void PeriodicOuterX2(MeshBlock *pmb, InterfaceField &a)
{
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  for (int k=ks; k<=ke; ++k) {
  for (int j=1; j<=(NGHOST); ++j) {

#pragma simd
    for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
      a.x1f(k,(je+j  ),i) = a.x1f(k,(js+j-1),i);
      a.x2f(k,(je+j+1),i) = a.x2f(k,(js+j  ),i);
      a.x3f(k,(je+j  ),i) = a.x3f(k,(js+j-1),i);
    }

  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void PeriodicInnerX3()
//  \brief  PERIODIC boundary conditions interface B, inner x3 boundary (ix3_bc=1)

void PeriodicInnerX3(MeshBlock *pmb, InterfaceField &a)
{
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  for (int k=1; k<=(NGHOST); ++k) {
  for (int j=js-(NGHOST); j<=je+(NGHOST); ++j) {

#pragma simd
    for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
      a.x1f((ks-k),j,i) = a.x1f((ke-k+1),j,i);
      a.x2f((ks-k),j,i) = a.x2f((ke-k+1),j,i);
      a.x3f((ks-k),j,i) = a.x3f((ke-k+1),j,i);
    }

  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void PeriodicOuterX3()
//  \brief  PERIODIC boundary conditions interface B, outer x3 boundary (ox3_bc=1)

void PeriodicOuterX3(MeshBlock *pmb, InterfaceField &a)
{
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  for (int k=1; k<=(NGHOST); ++k) {
  for (int j=js-(NGHOST); j<=je+(NGHOST); ++j) {

#pragma simd
    for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
      a.x1f((ke+k  ),j,i) = a.x1f((ks+k-1),j,i);
      a.x2f((ke+k  ),j,i) = a.x2f((ks+k-1),j,i);
      a.x3f((ke+k+1),j,i) = a.x3f((ks+k  ),j,i);
    }

  }}

  return;
}
