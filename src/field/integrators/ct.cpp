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
#include "../field.hpp"             // Field
#include "field_integrator.hpp"

// C++ headers
#include <algorithm>  // max(), min()

// Athena headers
#include "../../athena.hpp"         // enums, macros, Real
#include "../../athena_arrays.hpp"  // AthenaArray
#include "../../mesh.hpp"           // MeshBlock
#include "../../coordinates/coordinates.hpp"    // Areas/Lengths

//======================================================================================
//! \file ct.cpp
//  \brief
//======================================================================================

//--------------------------------------------------------------------------------------
//! \fn  void FieldIntegrator::CT
//  \brief

void FieldIntegrator::CT(MeshBlock *pmb, InterfaceField &bin, InterfaceField &bout,
  Real dt)
{
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  ComputeCornerEMFs(pmb);

//  AthenaArray<Real> b_x1f = b.x1f.ShallowCopy();
//  AthenaArray<Real> b_x2f = b.x2f.ShallowCopy();
//  AthenaArray<Real> b_x3f = b.x3f.ShallowCopy();
  AthenaArray<Real> emf1 = pmb->pfield->emf1.ShallowCopy();
  AthenaArray<Real> emf2 = pmb->pfield->emf2.ShallowCopy();
  AthenaArray<Real> emf3 = pmb->pfield->emf3.ShallowCopy();

  AthenaArray<Real> area = face_area_.ShallowCopy();
  AthenaArray<Real> len  = edge_length_.ShallowCopy();

//---- 1-D update

  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je+1; ++j) {
    pmb->pcoord->Face2Area(k,j,is,ie,area);
    pmb->pcoord->Edge3Length(k,j,is,ie+1,len);

    for (int i=is; i<=ie; ++i) {
      bout.x2f(k,j,i) = bin.x2f(k,j,i) +
        (dt/area(i))*(len(i+1)*emf3(k,j,i+1) - len(i)*emf3(k,j,i));
    }

  }}

  for (int k=ks; k<=ke+1; ++k) {
  for (int j=js; j<=je; ++j) {
    pmb->pcoord->Face3Area(k,j,is,ie,area);
    pmb->pcoord->Edge2Length(k,j,is,ie+1,len);

    for (int i=is; i<=ie; ++i) {
      bout.x3f(k,j,i) = bin.x3f(k,j,i) -
        (dt/area(i))*(len(i+1)*emf2(k,j,i+1) - len(i)*emf2(k,j,i));
    }

  }}

//---- 2-D update

  if (pmb->block_size.nx2 > 1) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      pmb->pcoord->Face1Area(k,j,is,ie,area);
      pmb->pcoord->Edge3Length(k,j,is,ie,len);

      for (int i=is; i<=ie+1; ++i) {
        bout.x1f(k,j,i) -= (dt/area(i))*(emf3(k,j+1,i) - emf3(k,j,i));
      }
    }}
    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
      pmb->pcoord->Face3Area(k,j,is,ie,area);
      pmb->pcoord->Edge1Length(k,j,is,ie,len);

      for (int i=is; i<=ie; ++i) {
        bout.x3f(k,j,i) += (dt/area(i))*(emf1(k,j+1,i) - emf1(k,j,i));
      }
    }}
  }

//---- 3-D update

  if (pmb->block_size.nx3 > 1) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      pmb->pcoord->Face1Area(k,j,is,ie,area);
      pmb->pcoord->Edge2Length(k,j,is,ie,len);

      for (int i=is; i<=ie+1; ++i) {
        bout.x1f(k,j,i) += (dt/area(i))*(emf2(k+1,j,i) - emf2(k,j,i));
      }
    }}
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
      pmb->pcoord->Face2Area(k,j,is,ie,area);
      pmb->pcoord->Edge1Length(k,j,is,ie,len);

      for (int i=is; i<=ie; ++i) {
        bout.x2f(k,j,i) -= (dt/area(i))*(emf1(k+1,j,i) - emf1(k,j,i));
      }
    }}
  }

  return;
}
