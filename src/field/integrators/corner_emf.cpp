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
#include "../field.hpp"    // Field
#include "field_integrator.hpp"

// C++ headers
#include <algorithm>  // max(), min()

// Athena headers
#include "../../athena.hpp"         // enums, macros, Real
#include "../../athena_arrays.hpp"  // AthenaArray
#include "../../mesh.hpp"           // MeshBlock
#include "../../fluid/fluid.hpp"    // Fluid

//======================================================================================
//! \file corner_emf.cpp
//  \brief
//======================================================================================

//--------------------------------------------------------------------------------------
//! \fn  void FieldIntegrator::ComputeCornerEMFs
//  \brief

void FieldIntegrator::ComputeCornerE(MeshBlock *pmb, AthenaArray<Real> &w, 
  AthenaArray<Real> &bcc)
{
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  AthenaArray<Real> e1 = pmb->pfield->e1.ShallowCopy();
  AthenaArray<Real> e2 = pmb->pfield->e2.ShallowCopy();
  AthenaArray<Real> e3 = pmb->pfield->e3.ShallowCopy();
  AthenaArray<Real> e_x1f = pmb->pfield->e.x1f.ShallowCopy();
  AthenaArray<Real> e_x2f = pmb->pfield->e.x2f.ShallowCopy();
  AthenaArray<Real> e_x3f = pmb->pfield->e.x3f.ShallowCopy();
  AthenaArray<Real> w_x1f = pmb->pfield->wght.x1f.ShallowCopy();
  AthenaArray<Real> w_x2f = pmb->pfield->wght.x2f.ShallowCopy();
  AthenaArray<Real> w_x3f = pmb->pfield->wght.x3f.ShallowCopy();
 
//---- 1-D update:
//  copy face-centered E-fields to edges and return.

  if (pmb->block_size.nx2 == 1) {
    for (int i=is; i<=ie+1; ++i) {
      e2(ks  ,js  ,i) = e_x1f(X1E2,ks,js,i);
      e2(ke+1,js  ,i) = e_x1f(X1E2,ks,js,i);
      e3(ks  ,js  ,i) = e_x1f(X1E3,ks,js,i);
      e3(ks  ,je+1,i) = e_x1f(X1E3,ks,js,i);
    }
    return;
  }

//---- 2-D/3-D update:
// Set BCs for face-centered E3, compute cell-centered E3=-(v X B)=VyBx-VxBy

  BoundaryValuesFaceCenteredE3(pmb);

  for (int k=ks; k<=ke; ++k) {
  for (int j=js-1; j<=je+1; ++j) {
    for (int i=is-1; i<=ie+1; ++i) {
      cc_e3_(k,j,i) = w(IVY,k,j,i)*bcc(IB1,k,j,i) - w(IVX,k,j,i)*bcc(IB2,k,j,i);
    }
  }}

// integrate E3 to corner using SG07

  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je+1; ++j) {
    for (int i=is; i<=ie+1; ++i) {
      Real de3_l2 = (1.0-w_x1f(k,j-1,i))*(e_x2f(X2E3,k,j,i  ) - cc_e3_(k,j-1,i  )) +
                    (    w_x1f(k,j-1,i))*(e_x2f(X2E3,k,j,i-1) - cc_e3_(k,j-1,i-1));

      Real de3_r2 = (1.0-w_x1f(k,j  ,i))*(e_x2f(X2E3,k,j,i  ) - cc_e3_(k,j,i  )) +
                    (    w_x1f(k,j  ,i))*(e_x2f(X2E3,k,j,i-1) - cc_e3_(k,j,i-1));

      Real de3_l1 = (1.0-w_x2f(k,j,i-1))*(e_x1f(X1E3,k,j  ,i) - cc_e3_(k,j  ,i-1)) +
                    (    w_x2f(k,j,i-1))*(e_x1f(X1E3,k,j-1,i) - cc_e3_(k,j-1,i-1));

      Real de3_r1 = (1.0-w_x2f(k,j,i  ))*(e_x1f(X1E3,k,j  ,i) - cc_e3_(k,j  ,i)) +
                    (    w_x2f(k,j,i  ))*(e_x1f(X1E3,k,j-1,i) - cc_e3_(k,j-1,i));

      e3(k,j,i) = 0.25*(de3_l1 + de3_r1 + de3_l2 + de3_r2 + e_x2f(X2E3,k,j,i-1) +
        e_x2f(X2E3,k,j,i) + e_x1f(X1E3,k,j-1,i) + e_x1f(X1E3,k,j,i));
    }
  }}

// for 2D: copy E1 and E2 to edges and return

  if (pmb->block_size.nx3 == 1) {
    for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie+1; ++i) {
      e2(ks  ,j,i) = e_x1f(X1E2,ks,j,i);
      e2(ke+1,j,i) = e_x1f(X1E2,ks,j,i);
    }}
    for (int j=js; j<=je+1; ++j) {
    for (int i=is; i<=ie; ++i) {
      e1(ks  ,j,i) = e_x2f(X2E1,ks,j,i);
      e1(ke+1,j,i) = e_x2f(X2E1,ks,j,i);
    }}
    return;
  }

//---- 3-D update:
// Set BCs for face-centered E1,E2, compute cell-centered E1,E2

  BoundaryValuesFaceCenteredE1(pmb);
  BoundaryValuesFaceCenteredE2(pmb);

// E1=-(v X B)=VzBy-VyBz
// E2=-(v X B)=VxBz-VzBx

  for (int k=ks-1; k<=ke+1; ++k) {
  for (int j=js-1; j<=je+1; ++j) {
    for (int i=is; i<=ie; ++i) {
      cc_e1_(k,j,i) = w(IVZ,k,j,i)*bcc(IB2,k,j,i) - w(IVY,k,j,i)*bcc(IB3,k,j,i);
    }
  }}

  for (int k=ks-1; k<=ke+1; ++k) {
  for (int j=js; j<=je; ++j) {
    for (int i=is-1; i<=ie+1; ++i) {
      cc_e2_(k,j,i) = w(IVX,k,j,i)*bcc(IB3,k,j,i) - w(IVZ,k,j,i)*bcc(IB1,k,j,i);
    }
  }}

// integrate E1 and E2 to corners using GS07 (E3 already done above)

  for (int k=ks; k<=ke+1; ++k) {
  for (int j=js; j<=je+1; ++j) {
    for (int i=is; i<=ie; ++i) {
      Real de1_l3 = (1.0-w_x2f(k-1,j,i))*(e_x3f(X3E1,k,j  ,i) - cc_e1_(k-1,j  ,i)) +
                    (    w_x2f(k-1,j,i))*(e_x3f(X3E1,k,j-1,i) - cc_e1_(k-1,j-1,i));

      Real de1_r3 = (1.0-w_x2f(k  ,j,i))*(e_x3f(X3E1,k,j  ,i) - cc_e1_(k,j  ,i)) +
                    (    w_x2f(k  ,j,i))*(e_x3f(X3E1,k,j-1,i) - cc_e1_(k,j-1,i));

      Real de1_l2 = (1.0-w_x3f(k,j-1,i))*(e_x2f(X2E1,k  ,j,i) - cc_e1_(k  ,j-1,i)) +
                    (    w_x3f(k,j-1,i))*(e_x2f(X2E1,k-1,j,i) - cc_e1_(k-1,j-1,i));

      Real de1_r2 = (1.0-w_x3f(k,j  ,i))*(e_x2f(X2E1,k  ,j,i) - cc_e1_(k  ,j,i)) +
                    (    w_x3f(k,j  ,i))*(e_x2f(X2E1,k-1,j,i) - cc_e1_(k-1,j,i));

      e1(k,j,i) = 0.25*(de1_l3 + de1_r3 + de1_l2 + de1_r2 + e_x2f(X2E1,k-1,j,i) +
        e_x2f(X2E1,k,j,i) + e_x3f(X3E1,k,j-1,i) + e_x3f(X3E1,k,j,i));
    }
  }}

  for (int k=ks; k<=ke+1; ++k) {
  for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie+1; ++i) {
      Real de2_l3 = (1.0-w_x1f(k-1,j,i))*(e_x3f(X3E2,k,j,i  ) - cc_e2_(k-1,j,i  )) +
                    (    w_x1f(k-1,j,i))*(e_x3f(X3E2,k,j,i-1) - cc_e2_(k-1,j,i-1));

      Real de2_r3 = (1.0-w_x1f(k,j  ,i))*(e_x3f(X3E2,k,j,i  ) - cc_e2_(k,j,i  )) +
                    (    w_x1f(k,j  ,i))*(e_x3f(X3E2,k,j,i-1) - cc_e2_(k,j,i-1));

      Real de2_l1 = (1.0-w_x3f(k,j,i-1))*(e_x1f(X1E2,k  ,j,i) - cc_e2_(k  ,j,i-1)) +
                    (    w_x3f(k,j,i-1))*(e_x1f(X1E2,k-1,j,i) - cc_e2_(k-1,j,i-1));

      Real de2_r1 = (1.0-w_x3f(k,j,i  ))*(e_x1f(X1E2,k  ,j,i) - cc_e2_(k  ,j,i)) +
                    (    w_x3f(k,j,i  ))*(e_x1f(X1E2,k-1,j,i) - cc_e2_(k-1,j,i));

      e2(k,j,i) = 0.25*(de2_l3 + de2_r3 + de2_l1 + de2_r1 + e_x3f(X3E2,k,j,i-1) +
        e_x3f(X3E2,k,j,i) + e_x1f(X1E2,k-1,j,i) + e_x1f(X1E2,k,j,i));
    }
  }}

  return;
}
