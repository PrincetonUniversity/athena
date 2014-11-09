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
//! \file outflow_bfield.cpp
//  \brief implements outflow BCs in each dimension for interface B-field
//======================================================================================
//--------------------------------------------------------------------------------------
//! \fn void OutflowInnerX1()
//  \brief  OUTFLOW boundary conditions interface B, inner x1 boundary (ix1_bc=1)

void OutflowInnerX1(MeshBlock *pmb, InterfaceField &a)
{
  int is = pmb->is;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {

#pragma simd
    for (int i=1; i<=(NGHOST); ++i) {
      a.x1(k,j,(is-i)) = a.x1(k,j,is);
      a.x2(k,j,(is-i)) = a.x2(k,j,is);
      a.x3(k,j,(is-i)) = a.x3(k,j,is);
    }

  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void OutflowOuterX1()
//  \brief  OUTFLOW boundary conditions interface B, outer x1 boundary (ox1_bc=1)

void OutflowOuterX1(MeshBlock *pmb, InterfaceField &a)
{
  int ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {

#pragma simd
    for (int i=1; i<=(NGHOST); ++i) {
      a.x1(k,j,(ie+i+1)) = a.x1(k,j,(ie+1));
      a.x2(k,j,(ie+i  )) = a.x2(k,j,(ie  ));
      a.x3(k,j,(ie+i  )) = a.x3(k,j,(ie  ));
    }

  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void OutflowInnerX2()
//  \brief  OUTFLOW boundary conditions interface B, inner x2 boundary (ix2_bc=1)

void OutflowInnerX2(MeshBlock *pmb, InterfaceField &a)
{
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js;
  int ks = pmb->ks, ke = pmb->ke;

  for (int k=ks; k<=ke; ++k) {
  for (int j=1; j<=(NGHOST); ++j) {

#pragma simd
    for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
      a.x1(k,(js-j),i) = a.x1(k,js,i);
      a.x2(k,(js-j),i) = a.x2(k,js,i);
      a.x3(k,(js-j),i) = a.x3(k,js,i);
    }

  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void OutflowOuterX2()
//  \brief  OUTFLOW boundary conditions interface B, outer x2 boundary (ox2_bc=1)

void OutflowOuterX2(MeshBlock *pmb, InterfaceField &a)
{
  int is = pmb->is, ie = pmb->ie;
  int je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  for (int k=ks; k<=ke; ++k) {
  for (int j=1; j<=(NGHOST); ++j) {

#pragma simd
    for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
      a.x1(k,(je+j  ),i) = a.x1(k,(je  ),i);
      a.x2(k,(je+j+1),i) = a.x2(k,(je+1),i);
      a.x3(k,(je+j  ),i) = a.x3(k,(je  ),i);
    }

  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void OutflowInnerX3()
//  \brief  OUTFLOW boundary conditions interface B, inner x3 boundary (ix3_bc=1)

void OutflowInnerX3(MeshBlock *pmb, InterfaceField &a)
{
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks;

  for (int k=1; k<=(NGHOST); ++k) {
  for (int j=js-(NGHOST); j<=je+(NGHOST); ++j) {

#pragma simd
    for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
      a.x1((ks-k),j,i) = a.x1(ks,j,i);
      a.x2((ks-k),j,i) = a.x2(ks,j,i);
      a.x3((ks-k),j,i) = a.x3(ks,j,i);
    }

  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void OutflowOuterX3()
//  \brief  OUTFLOW boundary conditions interface B, outer x3 boundary (ox3_bc=1)

void OutflowOuterX3(MeshBlock *pmb, InterfaceField &a)
{
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ke = pmb->ke;

  for (int k=1; k<=(NGHOST); ++k) {
  for (int j=js-(NGHOST); j<=je+(NGHOST); ++j) {

#pragma simd
    for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
      a.x1((ke+k  ),j,i) = a.x1((ke  ),j,i);
      a.x2((ke+k  ),j,i) = a.x2((ke  ),j,i);
      a.x3((ke+k+1),j,i) = a.x3((ke+1),j,i);
    }

  }}

  return;
}
