//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file corner_emf.cpp
//  \brief

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()

// Athena++ headers
#include "field.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro.hpp"
#include "field_diffusion/field_diffusion.hpp"
//----------------------------------------------------------------------------------------
//! \fn  void Field::ComputeCornerEMFs
//  \brief

void Field::ComputeCornerE(AthenaArray<Real> &w, AthenaArray<Real> &bcc) {
  MeshBlock *pmb = pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  AthenaArray<Real> e1,e2,e3,w_x1f,w_x2f,w_x3f;
  e1.InitWithShallowCopy(pmb->pfield->e.x1e);
  e2.InitWithShallowCopy(pmb->pfield->e.x2e);
  e3.InitWithShallowCopy(pmb->pfield->e.x3e);
  w_x1f.InitWithShallowCopy(pmb->pfield->wght.x1f);
  w_x2f.InitWithShallowCopy(pmb->pfield->wght.x2f);
  w_x3f.InitWithShallowCopy(pmb->pfield->wght.x3f);

//---- 1-D update:
//  copy face-centered E-fields to edges and return.

  if (pmb->block_size.nx2 == 1) {
    for (int i=is; i<=ie+1; ++i) {
      e2(ks  ,js  ,i) = e2_x1f(ks,js,i);
      e2(ke+1,js  ,i) = e2_x1f(ks,js,i);
      e3(ks  ,js  ,i) = e3_x1f(ks,js,i);
      e3(ks  ,je+1,i) = e3_x1f(ks,js,i);
    }
    return;
  }

//---- 2-D/3-D update:
  // E3=-(v X B)=VyBx-VxBy
  for (int k=ks; k<=ke; ++k) {
  for (int j=js-1; j<=je+1; ++j) {

#if GENERAL_RELATIVITY==1
    pmb->pcoord->CellMetric(k, j, is-1, ie+1, g_, gi_);
#pragma omp simd
    for (int i=is-1; i<=ie+1; ++i) {
      const Real &uu1 = w(IVX,k,j,i);
      const Real &uu2 = w(IVY,k,j,i);
      const Real &uu3 = w(IVZ,k,j,i);
      const Real &bb1 = bcc(IB1,k,j,i);
      const Real &bb2 = bcc(IB2,k,j,i);
      const Real &bb3 = bcc(IB3,k,j,i);
      Real alpha = std::sqrt(-1.0/gi_(I00,i));
      Real tmp = g_(I11,i)*SQR(uu1) + 2.0*g_(I12,i)*uu1*uu2 + 2.0*g_(I13,i)*uu1*uu3
               + g_(I22,i)*SQR(uu2) + 2.0*g_(I23,i)*uu2*uu3
               + g_(I33,i)*SQR(uu3);
      Real gamma = std::sqrt(1.0 + tmp);
      Real u0 = gamma / alpha;
      Real u1 = uu1 - alpha * gamma * gi_(I01,i);
      Real u2 = uu2 - alpha * gamma * gi_(I02,i);
      Real u3 = uu3 - alpha * gamma * gi_(I03,i);
      Real b0 = bb1 * (g_(I01,i)*u0 + g_(I11,i)*u1 + g_(I12,i)*u2 + g_(I13,i)*u3)
              + bb2 * (g_(I02,i)*u0 + g_(I12,i)*u1 + g_(I22,i)*u2 + g_(I23,i)*u3)
              + bb3 * (g_(I03,i)*u0 + g_(I13,i)*u1 + g_(I23,i)*u2 + g_(I33,i)*u3);
      Real b1 = (bb1 + b0 * u1) / u0;
      Real b2 = (bb2 + b0 * u2) / u0;
      Real b3 = (bb3 + b0 * u3) / u0;
      cc_e_(k,j,i) = b1 * u2 - b2 * u1;
    }
#else
#pragma omp simd
    for (int i=is-1; i<=ie+1; ++i) {
      cc_e_(k,j,i) = w(IVY,k,j,i)*bcc(IB1,k,j,i) - w(IVX,k,j,i)*bcc(IB2,k,j,i);
    }
#endif // GENERAL_RELATIVITY

  }}

  // integrate E3 to corner using SG07
  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je+1; ++j) {
#pragma omp simd
    for (int i=is; i<=ie+1; ++i) {
      Real de3_l2 = (1.0-w_x1f(k,j-1,i))*(e3_x2f(k,j,i  ) - cc_e_(k,j-1,i  )) +
                    (    w_x1f(k,j-1,i))*(e3_x2f(k,j,i-1) - cc_e_(k,j-1,i-1));

      Real de3_r2 = (1.0-w_x1f(k,j  ,i))*(e3_x2f(k,j,i  ) - cc_e_(k,j  ,i  )) +
                    (    w_x1f(k,j  ,i))*(e3_x2f(k,j,i-1) - cc_e_(k,j  ,i-1));

      Real de3_l1 = (1.0-w_x2f(k,j,i-1))*(e3_x1f(k,j  ,i) - cc_e_(k,j  ,i-1)) +
                    (    w_x2f(k,j,i-1))*(e3_x1f(k,j-1,i) - cc_e_(k,j-1,i-1));

      Real de3_r1 = (1.0-w_x2f(k,j,i  ))*(e3_x1f(k,j  ,i) - cc_e_(k,j  ,i  )) +
                    (    w_x2f(k,j,i  ))*(e3_x1f(k,j-1,i) - cc_e_(k,j-1,i  ));

      e3(k,j,i) = 0.25*(de3_l1 + de3_r1 + de3_l2 + de3_r2 + e3_x2f(k,j,i-1) +
        e3_x2f(k,j,i) + e3_x1f(k,j-1,i) + e3_x1f(k,j,i));
    }
  }}

  // for 2D: copy E1 and E2 to edges and return
  if (pmb->block_size.nx3 == 1) {
    for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie+1; ++i) {
      e2(ks  ,j,i) = e2_x1f(ks,j,i);
      e2(ke+1,j,i) = e2_x1f(ks,j,i);
    }}
    for (int j=js; j<=je+1; ++j) {
    for (int i=is; i<=ie; ++i) {
      e1(ks  ,j,i) = e1_x2f(ks,j,i);
      e1(ke+1,j,i) = e1_x2f(ks,j,i);
    }}
  } else {

//---- 3-D update:
    // integrate E1 to corners using GS07 (E3 already done above)
    // E1=-(v X B)=VzBy-VyBz
    for (int k=ks-1; k<=ke+1; ++k) {
    for (int j=js-1; j<=je+1; ++j) {

#if GENERAL_RELATIVITY==1
      pmb->pcoord->CellMetric(k, j, is, ie, g_, gi_);
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        const Real &uu1 = w(IVX,k,j,i);
        const Real &uu2 = w(IVY,k,j,i);
        const Real &uu3 = w(IVZ,k,j,i);
        const Real &bb1 = bcc(IB1,k,j,i);
        const Real &bb2 = bcc(IB2,k,j,i);
        const Real &bb3 = bcc(IB3,k,j,i);
        Real alpha = std::sqrt(-1.0/gi_(I00,i));
        Real tmp = g_(I11,i)*SQR(uu1) + 2.0*g_(I12,i)*uu1*uu2 + 2.0*g_(I13,i)*uu1*uu3
                 + g_(I22,i)*SQR(uu2) + 2.0*g_(I23,i)*uu2*uu3
                 + g_(I33,i)*SQR(uu3);
        Real gamma = std::sqrt(1.0 + tmp);
        Real u0 = gamma / alpha;
        Real u1 = uu1 - alpha * gamma * gi_(I01,i);
        Real u2 = uu2 - alpha * gamma * gi_(I02,i);
        Real u3 = uu3 - alpha * gamma * gi_(I03,i);
        Real b0 = bb1 * (g_(I01,i)*u0 + g_(I11,i)*u1 + g_(I12,i)*u2 + g_(I13,i)*u3)
                + bb2 * (g_(I02,i)*u0 + g_(I12,i)*u1 + g_(I22,i)*u2 + g_(I23,i)*u3)
                + bb3 * (g_(I03,i)*u0 + g_(I13,i)*u1 + g_(I23,i)*u2 + g_(I33,i)*u3);
        Real b1 = (bb1 + b0 * u1) / u0;
        Real b2 = (bb2 + b0 * u2) / u0;
        Real b3 = (bb3 + b0 * u3) / u0;
        cc_e_(k,j,i) = b2 * u3 - b3 * u2;
      }
#else
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        cc_e_(k,j,i) = w(IVZ,k,j,i)*bcc(IB2,k,j,i) - w(IVY,k,j,i)*bcc(IB3,k,j,i);
      }
#endif // GENERAL_RELATIVITY

    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        Real de1_l3 = (1.0-w_x2f(k-1,j,i))*(e1_x3f(k,j  ,i) - cc_e_(k-1,j  ,i)) +
                      (    w_x2f(k-1,j,i))*(e1_x3f(k,j-1,i) - cc_e_(k-1,j-1,i));

        Real de1_r3 = (1.0-w_x2f(k  ,j,i))*(e1_x3f(k,j  ,i) - cc_e_(k  ,j  ,i)) +
                      (    w_x2f(k  ,j,i))*(e1_x3f(k,j-1,i) - cc_e_(k  ,j-1,i));

        Real de1_l2 = (1.0-w_x3f(k,j-1,i))*(e1_x2f(k  ,j,i) - cc_e_(k  ,j-1,i)) +
                      (    w_x3f(k,j-1,i))*(e1_x2f(k-1,j,i) - cc_e_(k-1,j-1,i));

        Real de1_r2 = (1.0-w_x3f(k,j  ,i))*(e1_x2f(k  ,j,i) - cc_e_(k  ,j  ,i)) +
                      (    w_x3f(k,j  ,i))*(e1_x2f(k-1,j,i) - cc_e_(k-1,j  ,i));

        e1(k,j,i) = 0.25*(de1_l3 + de1_r3 + de1_l2 + de1_r2 + e1_x2f(k-1,j,i) +
          e1_x2f(k,j,i) + e1_x3f(k,j-1,i) + e1_x3f(k,j,i));
      }
    }}

    // integrate E2 to corners using GS07 (E3 already done above)
    // E2=-(v X B)=VxBz-VzBx
    for (int k=ks-1; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {

#if GENERAL_RELATIVITY==1
      pmb->pcoord->CellMetric(k, j, is-1, ie+1, g_, gi_);
#pragma omp simd
      for (int i=is-1; i<=ie+1; ++i) {
        const Real &uu1 = w(IVX,k,j,i);
        const Real &uu2 = w(IVY,k,j,i);
        const Real &uu3 = w(IVZ,k,j,i);
        const Real &bb1 = bcc(IB1,k,j,i);
        const Real &bb2 = bcc(IB2,k,j,i);
        const Real &bb3 = bcc(IB3,k,j,i);
        Real alpha = std::sqrt(-1.0/gi_(I00,i));
        Real tmp = g_(I11,i)*SQR(uu1) + 2.0*g_(I12,i)*uu1*uu2 + 2.0*g_(I13,i)*uu1*uu3
                 + g_(I22,i)*SQR(uu2) + 2.0*g_(I23,i)*uu2*uu3
                 + g_(I33,i)*SQR(uu3);
        Real gamma = std::sqrt(1.0 + tmp);
        Real u0 = gamma / alpha;
        Real u1 = uu1 - alpha * gamma * gi_(I01,i);
        Real u2 = uu2 - alpha * gamma * gi_(I02,i);
        Real u3 = uu3 - alpha * gamma * gi_(I03,i);
        Real b0 = bb1 * (g_(I01,i)*u0 + g_(I11,i)*u1 + g_(I12,i)*u2 + g_(I13,i)*u3)
                + bb2 * (g_(I02,i)*u0 + g_(I12,i)*u1 + g_(I22,i)*u2 + g_(I23,i)*u3)
                + bb3 * (g_(I03,i)*u0 + g_(I13,i)*u1 + g_(I23,i)*u2 + g_(I33,i)*u3);
        Real b1 = (bb1 + b0 * u1) / u0;
        Real b2 = (bb2 + b0 * u2) / u0;
        Real b3 = (bb3 + b0 * u3) / u0;
        cc_e_(k,j,i) = b3 * u1 - b1 * u3;
      }
#else
#pragma omp simd
      for (int i=is-1; i<=ie+1; ++i) {
        cc_e_(k,j,i) = w(IVX,k,j,i)*bcc(IB3,k,j,i) - w(IVZ,k,j,i)*bcc(IB1,k,j,i);
      }
#endif // GENERAL_RELATIVITY

    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        Real de2_l3 = (1.0-w_x1f(k-1,j,i))*(e2_x3f(k,j,i  ) - cc_e_(k-1,j,i  )) +
                      (    w_x1f(k-1,j,i))*(e2_x3f(k,j,i-1) - cc_e_(k-1,j,i-1));

        Real de2_r3 = (1.0-w_x1f(k,j  ,i))*(e2_x3f(k,j,i  ) - cc_e_(k  ,j,i  )) +
                      (    w_x1f(k,j  ,i))*(e2_x3f(k,j,i-1) - cc_e_(k  ,j,i-1));

        Real de2_l1 = (1.0-w_x3f(k,j,i-1))*(e2_x1f(k  ,j,i) - cc_e_(k  ,j,i-1)) +
                      (    w_x3f(k,j,i-1))*(e2_x1f(k-1,j,i) - cc_e_(k-1,j,i-1));

        Real de2_r1 = (1.0-w_x3f(k,j,i  ))*(e2_x1f(k  ,j,i) - cc_e_(k  ,j,i  )) +
                      (    w_x3f(k,j,i  ))*(e2_x1f(k-1,j,i) - cc_e_(k-1,j,i  ));

        e2(k,j,i) = 0.25*(de2_l3 + de2_r3 + de2_l1 + de2_r1 + e2_x3f(k,j,i-1) +
          e2_x3f(k,j,i) + e2_x1f(k-1,j,i) + e2_x1f(k,j,i));
      }
    }}
  }

  // add diffusion flux
  if (pfdif->field_diffusion_defined) pfdif->AddEMF(pfdif->e_oa, e);

  return;
}
