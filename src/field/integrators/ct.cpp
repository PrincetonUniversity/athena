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
//  \brief Constrained Transport implementation of dB/dt = -Curl(E), where E=-(v X B)

void FieldIntegrator::CT(MeshBlock *pmb, InterfaceField &b, AthenaArray<Real> &w,
  AthenaArray<Real> &bcc, Real dt)
{
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  ComputeCornerE(pmb, w, bcc);

  AthenaArray<Real> e1 = pmb->pfield->e1.ShallowCopy();
  AthenaArray<Real> e2 = pmb->pfield->e2.ShallowCopy();
  AthenaArray<Real> e3 = pmb->pfield->e3.ShallowCopy();

  AthenaArray<Real> area = face_area_.ShallowCopy();
  AthenaArray<Real> len   = edge_length_.ShallowCopy();
  AthenaArray<Real> lenp1 = edge_lengthp1_.ShallowCopy();

// First copy input into output

/**
  bout.x1f = bin.x1f;
  bout.x2f = bin.x2f;
  bout.x3f = bin.x3f;
**/

//---- 1-D update

  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je+1; ++j) {
    pmb->pcoord->Face2Area(k,j,is,ie,area);
    pmb->pcoord->Edge3Length(k,j,is,ie+1,len);

    for (int i=is; i<=ie; ++i) {
      b.x2f(k,j,i) += (dt/area(i))*(len(i+1)*e3(k,j,i+1) - len(i)*e3(k,j,i));
    }

  }}

  for (int k=ks; k<=ke+1; ++k) {
  for (int j=js; j<=je; ++j) {
    pmb->pcoord->Face3Area(k,j,is,ie,area);
    pmb->pcoord->Edge2Length(k,j,is,ie+1,len);

    for (int i=is; i<=ie; ++i) {
      b.x3f(k,j,i) -= (dt/area(i))*(len(i+1)*e2(k,j,i+1) - len(i)*e2(k,j,i));
    }

  }}

//---- 2-D update

  if (pmb->block_size.nx2 > 1) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      pmb->pcoord->Face1Area(k,j,is,ie+1,area);
      pmb->pcoord->Edge3Length(k,j  ,is,ie+1,len);
      pmb->pcoord->Edge3Length(k,j+1,is,ie+1,lenp1);

      for (int i=is; i<=ie+1; ++i) {
        b.x1f(k,j,i) -= (dt/area(i))*(lenp1(i)*e3(k,j+1,i) - len(i)*e3(k,j,i));
      }
    }}
    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
      pmb->pcoord->Face3Area(k,j,is,ie,area);
      pmb->pcoord->Edge1Length(k,j  ,is,ie,len);
      pmb->pcoord->Edge1Length(k,j+1,is,ie,lenp1);

      for (int i=is; i<=ie; ++i) {
        b.x3f(k,j,i) += (dt/area(i))*(lenp1(i)*e1(k,j+1,i) - len(i)*e1(k,j,i));
      }
    }}
  }

//---- 3-D update

  if (pmb->block_size.nx3 > 1) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      pmb->pcoord->Face1Area(k,j,is,ie+1,area);
      pmb->pcoord->Edge2Length(k  ,j,is,ie+1,len);
      pmb->pcoord->Edge2Length(k+1,j,is,ie+1,lenp1);

      for (int i=is; i<=ie+1; ++i) {
        b.x1f(k,j,i) += (dt/area(i))*(lenp1(i)*e2(k+1,j,i) - len(i)*e2(k,j,i));
      }
    }}
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
      pmb->pcoord->Face2Area(k,j,is,ie,area);
      pmb->pcoord->Edge2Length(k  ,j,is,ie,len);
      pmb->pcoord->Edge2Length(k+1,j,is,ie,lenp1);

      for (int i=is; i<=ie; ++i) {
        b.x2f(k,j,i) -= (dt/area(i))*(lenp1(i)*e1(k+1,j,i) - len(i)*e1(k,j,i));
      }
    }}
  }

  return;
}
