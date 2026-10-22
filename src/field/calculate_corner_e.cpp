//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file calculate_corner_e.cpp
//! \brief

// C headers

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "field.hpp"
#include "field_diffusion/field_diffusion.hpp"

//----------------------------------------------------------------------------------------
//! \fn  void Field::ComputeCornerE
//! \brief calculate the corner EMFs

void Field::ComputeCornerE(AthenaArray<Real> &w, AthenaArray<Real> &bcc) {
  MeshBlock *pmb = pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  AthenaArray<Real> &e1 = e.x1e, &e2 = e.x2e, &e3 = e.x3e,
                 &w_x1f = wght.x1f, &w_x2f = wght.x2f, &w_x3f = wght.x3f;
  //---- 1-D update:
  //  copy face-centered E-fields to edges and return.

  if (pmb->block_size.nx2 == 1) {
    for (int i=is; i<=ie+1; ++i) {
      e2(ks  ,js  ,i) = e2_x1f(ks,js,i);
      e2(ke+1,js  ,i) = e2_x1f(ks,js,i);
      e3(ks  ,js  ,i) = e3_x1f(ks,js,i);
      e3(ks  ,je+1,i) = e3_x1f(ks,js,i);
    }
    if (!STS_ENABLED) // add diffusion flux
      if (fdif.field_diffusion_defined) fdif.AddEMF(fdif.e_oa, e);
    return;
  }

  if (pmb->block_size.nx3 == 1) {
    //---- 2-D update - cc_e_ is 3D array
    for (int k=ks; k<=ke; ++k) {
      for (int j=js-1; j<=je+1; ++j) {
        // E3=-(v X B)=VyBx-VxBy
#if GENERAL_RELATIVITY==1  // GR
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
#elif RELATIVISTIC_DYNAMICS==1  // SR
#pragma omp simd
        for (int i=is-1; i<=ie+1; ++i) {
          const Real &u1 = w(IVX,k,j,i);
          const Real &u2 = w(IVY,k,j,i);
          const Real &u3 = w(IVZ,k,j,i);
          Real u0 = std::sqrt(1.0 + SQR(u1) + SQR(u2) + SQR(u3));
          cc_e_(k,j,i) = (u2 * bcc(IB1,k,j,i) - u1 * bcc(IB2,k,j,i)) / u0;
        }
#else  // Newtonian
#pragma omp simd
        for (int i=is-1; i<=ie+1; ++i) {
          cc_e_(k,j,i) = w(IVY,k,j,i)*bcc(IB1,k,j,i) - w(IVX,k,j,i)*bcc(IB2,k,j,i);
        }
#endif // GENERAL_RELATIVITY
      }
    }
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie+1; ++i) {
        e2(ke+1,j,i) = e2(ks  ,j,i) = e2_x1f(ks,j,i);
      }
    }
    for (int j=js; j<=je+1; ++j) {
      for (int i=is; i<=ie; ++i) {
        e1(ke+1,j,i) = e1(ks  ,j,i) = e1_x2f(ks,j,i);
      }
    }

    // integrate E3 to corner using SG07
    // The upwind weights are stored as a signed velocity v in [-0.5, 0.5], and the
    // sums are grouped so that each pair of terms exchanged by a reflection is added
    // first; this keeps the corner EMF bit-exactly mirror symmetric (#165).
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          Real de3_l2 = (0.5-w_x1f(k,j-1,i))*(e3_x2f(k,j,i  ) - cc_e_(k,j-1,i  )) +
                        (0.5+w_x1f(k,j-1,i))*(e3_x2f(k,j,i-1) - cc_e_(k,j-1,i-1));
          Real de3_r2 = (0.5-w_x1f(k,j  ,i))*(e3_x2f(k,j,i  ) - cc_e_(k,j  ,i  )) +
                        (0.5+w_x1f(k,j  ,i))*(e3_x2f(k,j,i-1) - cc_e_(k,j  ,i-1));
          Real de3_l1 = (0.5-w_x2f(k,j,i-1))*(e3_x1f(k,j  ,i) - cc_e_(k,j  ,i-1)) +
                        (0.5+w_x2f(k,j,i-1))*(e3_x1f(k,j-1,i) - cc_e_(k,j-1,i-1));
          Real de3_r1 = (0.5-w_x2f(k,j,i  ))*(e3_x1f(k,j  ,i) - cc_e_(k,j  ,i  )) +
                        (0.5+w_x2f(k,j,i  ))*(e3_x1f(k,j-1,i) - cc_e_(k,j-1,i  ));

          e3(k,j,i) = 0.25*(((de3_l1 + de3_r1) + (de3_l2 + de3_r2)) +
                            ((e3_x2f(k,j,i-1) + e3_x2f(k,j,i)) +
                             (e3_x1f(k,j-1,i) + e3_x1f(k,j,i))));
        }
      }
    }
  } else {
    // 3-D updates - cc_e_ is 4D array
    for (int k=ks-1; k<=ke+1; ++k) {
      for (int j=js-1; j<=je+1; ++j) {
        // E1=-(v X B)=VzBy-VyBz
        // E2=-(v X B)=VxBz-VzBx
        // E3=-(v X B)=VyBx-VxBy
#if GENERAL_RELATIVITY==1  // GR
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
          cc_e_(IB1,k,j,i) = b2 * u3 - b3 * u2;
          cc_e_(IB2,k,j,i) = b3 * u1 - b1 * u3;
          cc_e_(IB3,k,j,i) = b1 * u2 - b2 * u1;
        }
#elif RELATIVISTIC_DYNAMICS==1  // SR
#pragma omp simd
        for (int i=is-1; i<=ie+1; ++i) {
          const Real &u1 = w(IVX,k,j,i);
          const Real &u2 = w(IVY,k,j,i);
          const Real &u3 = w(IVZ,k,j,i);
          Real u0 = std::sqrt(1.0 + SQR(u1) + SQR(u2) + SQR(u3));
          cc_e_(IB1,k,j,i) = (u3 * bcc(IB2,k,j,i) - u2 * bcc(IB3,k,j,i)) / u0;
          cc_e_(IB2,k,j,i) = (u1 * bcc(IB3,k,j,i) - u3 * bcc(IB1,k,j,i)) / u0;
          cc_e_(IB3,k,j,i) = (u2 * bcc(IB1,k,j,i) - u1 * bcc(IB2,k,j,i)) / u0;
        }
#else  // Newtonian
#pragma omp simd
        for (int i=is-1; i<=ie+1; ++i) {
          cc_e_(IB1,k,j,i) = w(IVZ,k,j,i)*bcc(IB2,k,j,i) - w(IVY,k,j,i)*bcc(IB3,k,j,i);
          cc_e_(IB2,k,j,i) = w(IVX,k,j,i)*bcc(IB3,k,j,i) - w(IVZ,k,j,i)*bcc(IB1,k,j,i);
          cc_e_(IB3,k,j,i) = w(IVY,k,j,i)*bcc(IB1,k,j,i) - w(IVX,k,j,i)*bcc(IB2,k,j,i);
        }
#endif // GENERAL_RELATIVITY
      }
    }

    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          // integrate E1,E2,E3 to corner using SG07
          Real de1_l3 = (0.5-w_x2f(k-1,j,i))*(e1_x3f(k,j  ,i) - cc_e_(IB1,k-1,j  ,i)) +
                        (0.5+w_x2f(k-1,j,i))*(e1_x3f(k,j-1,i) - cc_e_(IB1,k-1,j-1,i));
          Real de1_r3 = (0.5-w_x2f(k  ,j,i))*(e1_x3f(k,j  ,i) - cc_e_(IB1,k  ,j  ,i)) +
                        (0.5+w_x2f(k  ,j,i))*(e1_x3f(k,j-1,i) - cc_e_(IB1,k  ,j-1,i));
          Real de1_l2 = (0.5-w_x3f(k,j-1,i))*(e1_x2f(k  ,j,i) - cc_e_(IB1,k  ,j-1,i)) +
                        (0.5+w_x3f(k,j-1,i))*(e1_x2f(k-1,j,i) - cc_e_(IB1,k-1,j-1,i));
          Real de1_r2 = (0.5-w_x3f(k,j  ,i))*(e1_x2f(k  ,j,i) - cc_e_(IB1,k  ,j  ,i)) +
                        (0.5+w_x3f(k,j  ,i))*(e1_x2f(k-1,j,i) - cc_e_(IB1,k-1,j  ,i));

          e1(k,j,i) = 0.25*(((de1_l3 + de1_r3) + (de1_l2 + de1_r2)) +
                            ((e1_x2f(k-1,j,i) + e1_x2f(k,j,i)) +
                             (e1_x3f(k,j-1,i) + e1_x3f(k,j,i))));

          Real de2_l3 = (0.5-w_x1f(k-1,j,i))*(e2_x3f(k,j,i  ) - cc_e_(IB2,k-1,j,i  )) +
                        (0.5+w_x1f(k-1,j,i))*(e2_x3f(k,j,i-1) - cc_e_(IB2,k-1,j,i-1));
          Real de2_r3 = (0.5-w_x1f(k,j  ,i))*(e2_x3f(k,j,i  ) - cc_e_(IB2,k  ,j,i  )) +
                        (0.5+w_x1f(k,j  ,i))*(e2_x3f(k,j,i-1) - cc_e_(IB2,k  ,j,i-1));
          Real de2_l1 = (0.5-w_x3f(k,j,i-1))*(e2_x1f(k  ,j,i) - cc_e_(IB2,k  ,j,i-1)) +
                        (0.5+w_x3f(k,j,i-1))*(e2_x1f(k-1,j,i) - cc_e_(IB2,k-1,j,i-1));
          Real de2_r1 = (0.5-w_x3f(k,j,i  ))*(e2_x1f(k  ,j,i) - cc_e_(IB2,k  ,j,i  )) +
                        (0.5+w_x3f(k,j,i  ))*(e2_x1f(k-1,j,i) - cc_e_(IB2,k-1,j,i  ));

          e2(k,j,i) = 0.25*(((de2_l3 + de2_r3) + (de2_l1 + de2_r1)) +
                            ((e2_x3f(k,j,i-1) + e2_x3f(k,j,i)) +
                             (e2_x1f(k-1,j,i) + e2_x1f(k,j,i))));

          Real de3_l2 = (0.5-w_x1f(k,j-1,i))*(e3_x2f(k,j,i  ) - cc_e_(IB3,k,j-1,i  )) +
                        (0.5+w_x1f(k,j-1,i))*(e3_x2f(k,j,i-1) - cc_e_(IB3,k,j-1,i-1));
          Real de3_r2 = (0.5-w_x1f(k,j  ,i))*(e3_x2f(k,j,i  ) - cc_e_(IB3,k,j  ,i  )) +
                        (0.5+w_x1f(k,j  ,i))*(e3_x2f(k,j,i-1) - cc_e_(IB3,k,j  ,i-1));
          Real de3_l1 = (0.5-w_x2f(k,j,i-1))*(e3_x1f(k,j  ,i) - cc_e_(IB3,k,j  ,i-1)) +
                        (0.5+w_x2f(k,j,i-1))*(e3_x1f(k,j-1,i) - cc_e_(IB3,k,j-1,i-1));
          Real de3_r1 = (0.5-w_x2f(k,j,i  ))*(e3_x1f(k,j  ,i) - cc_e_(IB3,k,j  ,i  )) +
                        (0.5+w_x2f(k,j,i  ))*(e3_x1f(k,j-1,i) - cc_e_(IB3,k,j-1,i  ));

          e3(k,j,i) = 0.25*(((de3_l1 + de3_r1) + (de3_l2 + de3_r2)) +
                            ((e3_x2f(k,j,i-1) + e3_x2f(k,j,i)) +
                             (e3_x1f(k,j-1,i) + e3_x1f(k,j,i))));
        }
      }
    }
  }

  if (!STS_ENABLED) // add diffusion flux
    if (fdif.field_diffusion_defined) fdif.AddEMF(fdif.e_oa, e);

  return;
}

//----------------------------------------------------------------------------------------
//! \fn  void Field::ComputeCornerE_STS
//! \brief Compute corner E for STS

void Field::ComputeCornerE_STS() {
  // add diffusion flux
  if (fdif.field_diffusion_defined) fdif.AddEMF(fdif.e_oa, e);
  return;
}

namespace {
// Point value at the edge midpoint from the line averages on an edge and its two
// neighbors along the edge direction (uniform spacing): <q> - (h^2/24) d2q/ds2 with a
// central second difference, written as (q_m + q_p) - 2q so that it is bit-exactly
// mirror symmetric. EdgeAverage is the same correction applied to point values.
inline Real EdgePointValue(const Real qm, const Real q, const Real qp) {
  return q - ((qm + qp) - 2.0*q)*(1.0/24.0);
}
inline Real EdgeAverage(const Real qm, const Real q, const Real qp) {
  return q + ((qm + qp) - 2.0*q)*(1.0/24.0);
}
// ...for the 4th-dim slice n of a 4D array (n = 0 for a 3D array), along x1, x2, x3
inline Real PointX1(const AthenaArray<Real> &q, int n, int k, int j, int i) {
  return EdgePointValue(q(n,k,j,i-1), q(n,k,j,i), q(n,k,j,i+1));
}
inline Real PointX2(const AthenaArray<Real> &q, int n, int k, int j, int i) {
  return EdgePointValue(q(n,k,j-1,i), q(n,k,j,i), q(n,k,j+1,i));
}
inline Real PointX3(const AthenaArray<Real> &q, int n, int k, int j, int i) {
  return EdgePointValue(q(n,k-1,j,i), q(n,k,j,i), q(n,k+1,j,i));
}
}  // namespace

//----------------------------------------------------------------------------------------
//! \fn  void Field::ComputeCornerE_UCT4
//! \brief calculate the corner (edge-averaged) EMFs for the fourth-order UCT scheme
//!
//! Computes the multidimensional HLL average of the four one-sided corner states plus
//! two field-jump diffusion terms, following Londrillo & Del Zanna (2004) eq. 56 and
//! Felker & Stone (2018) section 4.2. The corner-reconstructed velocity/field states
//! and the alpha^{+,-} fast magnetosonic wavespeeds are computed during
//! Hydro::CalculateFluxes when xorder=4.
//!
//! In 2D the corner states are point values and the corner EMF is the edge average.
//! In 3D the corner states are line averages along the edge (transverse PPM of the
//! face-averaged states), so they are first deconvolved to point values at the edge
//! midpoint, q_pt = <q> - (h^2/24) d2<q>/ds2, the corner EMF is evaluated from the
//! point states, and the result is re-averaged along the edge,
//! <E> = E_pt + (h^2/24) d2E_pt/ds2, with central second differences along the edge
//! direction s (EdgePointValue/EdgeAverage). This is the edge analog of the face-flux
//! Laplacian correction in Hydro::CalculateFluxes and keeps the edge-averaged EMF
//! fourth-order accurate when the EMF is a nonlinear function of the varying
//! quantities; the product of the line averages alone is O(h^2) (mhd4/rotating_column).
//! The point EMF is needed on one extra row at each end of the edge and the corner
//! states on two (kl_e3, il_e1, jl_e2 in Hydro::CalculateFluxes).
//!
//! Sign convention of the diffusion terms: for each transverse field component B_t
//! differenced across direction d at the edge, the term enters with the sign such that
//! the corner EMF reduces to the 1D HLL flux of B_t along d in the plane-symmetric
//! limit. Since (dB/dt)_t = -(curl E)_t, the EMF component is (+/-) the flux of B_t
//! along d, which fixes the sign of -alpha^+ alpha^- (B_t^R - B_t^L)/(alpha^+ + alpha^-)
//! per component (see the E3 comments below; E1, E2 follow by cyclic permutation).

void Field::ComputeCornerE_UCT4() {
  MeshBlock *pmb = pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  AthenaArray<Real> &e1 = e.x1e, &e2 = e.x2e, &e3 = e.x3e;

  //---- 1-D update:
  //  copy face-centered E-fields to edges and return.

  if (pmb->block_size.nx2 == 1) {
    for (int i=is; i<=ie+1; ++i) {
      e2(ks  ,js  ,i) = e2_x1f(ks,js,i);
      e2(ke+1,js  ,i) = e2_x1f(ks,js,i);
      e3(ks  ,js  ,i) = e3_x1f(ks,js,i);
      e3(ks  ,je+1,i) = e3_x1f(ks,js,i);
    }
    if (!STS_ENABLED) // add diffusion flux
      if (fdif.field_diffusion_defined) fdif.AddEMF(fdif.e_oa, e);
    return;
  }

  if (pmb->block_size.nx3 == 1) {
    //---- 2-D update: the corner states are point values
    // E3=-(v X B)=VyBx-VxBy
    // integrate E3 to corner using UCT
    // corner state naming: E/W = L1/R1, N/S = L2/R2 one-sided limits w.r.t. the corner
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          // corner wavespeeds: take the max of the face-centered estimates on the two
          // faces adjacent to the corner (symmetric w.r.t. the corner, following LDZ04)
          Real alpha_plus_x = std::max(std::abs(alpha_plus_x1_(k,j,i)),
                                       std::abs(alpha_plus_x1_(k,j-1,i)));
          Real alpha_minus_x = std::max(std::abs(alpha_minus_x1_(k,j,i)),
                                        std::abs(alpha_minus_x1_(k,j-1,i)));

          Real alpha_plus_y = std::max(std::abs(alpha_plus_x2_(k,j,i)),
                                       std::abs(alpha_plus_x2_(k,j,i-1)));
          Real alpha_minus_y = std::max(std::abs(alpha_minus_x2_(k,j,i)),
                                        std::abs(alpha_minus_x2_(k,j,i-1)));

          // Londrillo & Del Zanna (2004) eq. 56; L states are weighted by alpha^+
          Real e3_NE, e3_SE, e3_NW, e3_SW;
          Real bx_S_minus_bx_N, by_W_minus_by_E;
          e3_NE = v_NE(1,k,j,i)*bx_N(k,j,i) - v_NE(0,k,j,i)*by_E(k,j,i);
          e3_SE = v_SE(1,k,j,i)*bx_S(k,j,i) - v_SE(0,k,j,i)*by_E(k,j,i);
          e3_NW = v_NW(1,k,j,i)*bx_N(k,j,i) - v_NW(0,k,j,i)*by_W(k,j,i);
          e3_SW = v_SW(1,k,j,i)*bx_S(k,j,i) - v_SW(0,k,j,i)*by_W(k,j,i);

          // diffusion terms: +E3 is the x2-flux of B1 (hence "-" on the B1 jump across
          // x2), and -E3 is the x1-flux of B2 (hence "+" on the B2 jump across x1)
          bx_S_minus_bx_N = bx_S(k,j,i) - bx_N(k,j,i);
          by_W_minus_by_E = by_W(k,j,i) - by_E(k,j,i);
          // (the two corner states exchanged by an x1 reflection are added first, so that
          // the sum is mirror symmetric, #165)
          e3(k,j,i) = (
              (alpha_plus_x*alpha_plus_y*e3_NE + alpha_minus_x*alpha_plus_y*e3_NW) +
              (alpha_plus_x*alpha_minus_y*e3_SE + alpha_minus_x*alpha_minus_y*e3_SW)
                       ) / ((alpha_plus_x + alpha_minus_x)*(alpha_plus_y + alpha_minus_y))
                      - alpha_plus_y*alpha_minus_y*bx_S_minus_bx_N / (
                          alpha_plus_y + alpha_minus_y)
                      + alpha_plus_x*alpha_minus_x*by_W_minus_by_E/(
                          alpha_plus_x + alpha_minus_x);
        }
      }
    }


    // copy E1 and E2 to edges
    // (e1_x2f and e2_x1f are the 4th-order face-averaged fluxes of B3 in 2D)
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie+1; ++i) {
        e2(ks  ,j,i) = e2_x1f(ks,j,i);
        e2(ke+1,j,i) = e2_x1f(ks,j,i);
      }
    }
    for (int j=js; j<=je+1; ++j) {
      for (int i=is; i<=ie; ++i) {
        e1(ks  ,j,i) = e1_x2f(ks,j,i);
        e1(ke+1,j,i) = e1_x2f(ks,j,i);
      }
    }
  } else {
    //---- 3-D update:
    // The corner states are line averages along the edge (transverse PPM of the
    // face-averaged Riemann states and of the face fields). Deconvolve them to point
    // values at the edge midpoint, evaluate the UCT corner EMF from the point states on
    // one extra row at each end of the edge, and re-average the EMF along the edge.
    // The point EMFs use the Field conversion scratch, which is free during this task.
    AthenaArray<Real> &e3_pt = scr1_kji_x1fc_, &e1_pt = scr2_kji_x2fc_,
                      &e2_pt = scr3_kji_x3fc_;

    // E3=-(v X B)=VyBx-VxBy: point EMF on x3 rows ks-1..ke+1
    // corner state naming: E/W = L1/R1, N/S = L2/R2 one-sided limits w.r.t. the corner
    for (int k=ks-1; k<=ke+1; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          // corner wavespeeds: take the max of the face-centered estimates on the two
          // faces adjacent to the corner (symmetric w.r.t. the corner, following LDZ04)
          Real alpha_plus_x = std::max(std::abs(alpha_plus_x1_(k,j,i)),
                                       std::abs(alpha_plus_x1_(k,j-1,i)));
          Real alpha_minus_x = std::max(std::abs(alpha_minus_x1_(k,j,i)),
                                        std::abs(alpha_minus_x1_(k,j-1,i)));

          Real alpha_plus_y = std::max(std::abs(alpha_plus_x2_(k,j,i)),
                                       std::abs(alpha_plus_x2_(k,j,i-1)));
          Real alpha_minus_y = std::max(std::abs(alpha_minus_x2_(k,j,i)),
                                        std::abs(alpha_minus_x2_(k,j,i-1)));

          // point values of the corner states at the edge midpoint (deconvolved in x3)
          Real vx_NE = PointX3(v_NE, 0, k, j, i), vy_NE = PointX3(v_NE, 1, k, j, i);
          Real vx_SE = PointX3(v_SE, 0, k, j, i), vy_SE = PointX3(v_SE, 1, k, j, i);
          Real vx_NW = PointX3(v_NW, 0, k, j, i), vy_NW = PointX3(v_NW, 1, k, j, i);
          Real vx_SW = PointX3(v_SW, 0, k, j, i), vy_SW = PointX3(v_SW, 1, k, j, i);
          Real bxN = PointX3(bx_N, 0, k, j, i), bxS = PointX3(bx_S, 0, k, j, i);
          Real byE = PointX3(by_E, 0, k, j, i), byW = PointX3(by_W, 0, k, j, i);

          // Londrillo & Del Zanna (2004) eq. 56; L states are weighted by alpha^+
          Real e3_NE = vy_NE*bxN - vx_NE*byE;
          Real e3_SE = vy_SE*bxS - vx_SE*byE;
          Real e3_NW = vy_NW*bxN - vx_NW*byW;
          Real e3_SW = vy_SW*bxS - vx_SW*byW;

          // diffusion terms: +E3 is the x2-flux of B1 (hence "-" on the B1 jump across
          // x2), and -E3 is the x1-flux of B2 (hence "+" on the B2 jump across x1)
          Real bx_S_minus_bx_N = bxS - bxN;
          Real by_W_minus_by_E = byW - byE;
          // (the two corner states exchanged by an x1 reflection are added first, so
          // that the sum is mirror symmetric, #165)
          e3_pt(k,j,i) = (
              (alpha_plus_x*alpha_plus_y*e3_NE + alpha_minus_x*alpha_plus_y*e3_NW) +
              (alpha_plus_x*alpha_minus_y*e3_SE + alpha_minus_x*alpha_minus_y*e3_SW)
                         ) / ((alpha_plus_x + alpha_minus_x)
                              *(alpha_plus_y + alpha_minus_y))
                        - alpha_plus_y*alpha_minus_y*bx_S_minus_bx_N / (
                            alpha_plus_y + alpha_minus_y)
                        + alpha_plus_x*alpha_minus_x*by_W_minus_by_E/(
                            alpha_plus_x + alpha_minus_x);
        }
      }
    }
    // edge average of E3 along x3
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          e3(k,j,i) = EdgeAverage(e3_pt(k-1,j,i), e3_pt(k,j,i), e3_pt(k+1,j,i));
        }
      }
    }

    // E1=-(v X B)=VzBy-VyBz: point EMF on x1 rows is-1..ie+1
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=is-1; i<=ie+1; ++i) {
          // edge wavespeeds: max of the face-centered estimates on the two faces
          // adjacent to the edge (symmetric w.r.t. the edge, as for E3 above)
          Real alpha_plus_z = std::max(std::abs(alpha_plus_x3_(k,j,i)),
                                       std::abs(alpha_plus_x3_(k,j-1,i)));
          Real alpha_minus_z = std::max(std::abs(alpha_minus_x3_(k,j,i)),
                                        std::abs(alpha_minus_x3_(k,j-1,i)));

          Real alpha_plus_y = std::max(std::abs(alpha_plus_x2_(k,j,i)),
                                       std::abs(alpha_plus_x2_(k-1,j,i)));
          Real alpha_minus_y = std::max(std::abs(alpha_minus_x2_(k,j,i)),
                                        std::abs(alpha_minus_x2_(k-1,j,i)));

          // point values of the corner states at the edge midpoint (deconvolved in x1)
          Real vy_L3L2 = PointX1(v_L3L2, 1, k, j, i);
          Real vz_L3L2 = PointX1(v_L3L2, 2, k, j, i);
          Real vy_L3R2 = PointX1(v_L3R2, 1, k, j, i);
          Real vz_L3R2 = PointX1(v_L3R2, 2, k, j, i);
          Real vy_R3L2 = PointX1(v_R3L2, 1, k, j, i);
          Real vz_R3L2 = PointX1(v_R3L2, 2, k, j, i);
          Real vy_R3R2 = PointX1(v_R3R2, 1, k, j, i);
          Real vz_R3R2 = PointX1(v_R3R2, 2, k, j, i);
          Real byL3 = PointX1(by_L3, 0, k, j, i), byR3 = PointX1(by_R3, 0, k, j, i);
          Real bzL2 = PointX1(bz_L2, 0, k, j, i), bzR2 = PointX1(bz_R2, 0, k, j, i);

          Real e1_L3L2 = vz_L3L2*byL3 - vy_L3L2*bzL2;
          Real e1_L3R2 = vz_L3R2*byL3 - vy_L3R2*bzR2;
          Real e1_R3L2 = vz_R3L2*byR3 - vy_R3L2*bzL2;
          Real e1_R3R2 = vz_R3R2*byR3 - vy_R3R2*bzR2;

          // diffusion terms: -E1 is the x2-flux of B3 ("+" on the B3 jump across x2),
          // and +E1 is the x3-flux of B2 ("-" on the B2 jump across x3)
          Real by_R3_minus_by_L3 = byR3 - byL3;
          Real bz_R2_minus_bz_L2 = bzR2 - bzL2;
          e1_pt(k,j,i) = ((alpha_plus_z*alpha_plus_y*e1_L3L2 +
                           alpha_minus_z*alpha_plus_y*e1_R3L2) +
                          (alpha_plus_z*alpha_minus_y*e1_L3R2 +
                           alpha_minus_z*alpha_minus_y*e1_R3R2))
              / ((alpha_plus_z + alpha_minus_z)*(alpha_plus_y + alpha_minus_y))
              + alpha_plus_y*alpha_minus_y*bz_R2_minus_bz_L2 /
              (alpha_plus_y + alpha_minus_y)
              - alpha_plus_z*alpha_minus_z*by_R3_minus_by_L3 /
              (alpha_plus_z + alpha_minus_z);
        }
      }
    }
    // edge average of E1 along x1
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          e1(k,j,i) = EdgeAverage(e1_pt(k,j,i-1), e1_pt(k,j,i), e1_pt(k,j,i+1));
        }
      }
    }

    // E2=-(v X B)=VxBz-VzBx: point EMF on x2 rows js-1..je+1
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js-1; j<=je+1; ++j) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          Real alpha_plus_z = std::max(std::abs(alpha_plus_x3_(k,j,i)),
                                       std::abs(alpha_plus_x3_(k,j,i-1)));
          Real alpha_minus_z = std::max(std::abs(alpha_minus_x3_(k,j,i)),
                                        std::abs(alpha_minus_x3_(k,j,i-1)));

          Real alpha_plus_x = std::max(std::abs(alpha_plus_x1_(k,j,i)),
                                       std::abs(alpha_plus_x1_(k-1,j,i)));
          Real alpha_minus_x = std::max(std::abs(alpha_minus_x1_(k,j,i)),
                                        std::abs(alpha_minus_x1_(k-1,j,i)));

          // point values of the corner states at the edge midpoint (deconvolved in x2)
          Real vx_L3L1 = PointX2(v_L3L1, 0, k, j, i);
          Real vz_L3L1 = PointX2(v_L3L1, 2, k, j, i);
          Real vx_L3R1 = PointX2(v_L3R1, 0, k, j, i);
          Real vz_L3R1 = PointX2(v_L3R1, 2, k, j, i);
          Real vx_R3L1 = PointX2(v_R3L1, 0, k, j, i);
          Real vz_R3L1 = PointX2(v_R3L1, 2, k, j, i);
          Real vx_R3R1 = PointX2(v_R3R1, 0, k, j, i);
          Real vz_R3R1 = PointX2(v_R3R1, 2, k, j, i);
          Real bxL3 = PointX2(bx_L3, 0, k, j, i), bxR3 = PointX2(bx_R3, 0, k, j, i);
          Real bzL1 = PointX2(bz_L1, 0, k, j, i), bzR1 = PointX2(bz_R1, 0, k, j, i);

          Real e2_L3L1 = - vz_L3L1*bxL3 + vx_L3L1*bzL1;
          Real e2_L3R1 = - vz_L3R1*bxL3 + vx_L3R1*bzR1;
          Real e2_R3L1 = - vz_R3L1*bxR3 + vx_R3L1*bzL1;
          Real e2_R3R1 = - vz_R3R1*bxR3 + vx_R3R1*bzR1;

          // diffusion terms: +E2 is the x1-flux of B3 ("-" on the B3 jump across x1),
          // and -E2 is the x3-flux of B1 ("+" on the B1 jump across x3)
          Real bx_R3_minus_bx_L3 = bxR3 - bxL3;
          Real bz_R1_minus_bz_L1 = bzR1 - bzL1;
          e2_pt(k,j,i) = ((alpha_plus_z*alpha_plus_x*e2_L3L1 +
                           alpha_minus_z*alpha_plus_x*e2_R3L1) +
                          (alpha_plus_z*alpha_minus_x*e2_L3R1 +
                           alpha_minus_z*alpha_minus_x*e2_R3R1))
              / ((alpha_plus_z + alpha_minus_z)*(alpha_plus_x + alpha_minus_x))
              - alpha_plus_x*alpha_minus_x*bz_R1_minus_bz_L1 /
              (alpha_plus_x + alpha_minus_x)
              + alpha_plus_z*alpha_minus_z*bx_R3_minus_bx_L3 /
              (alpha_plus_z + alpha_minus_z);
        }
      }
    }
    // edge average of E2 along x2
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          e2(k,j,i) = EdgeAverage(e2_pt(k,j-1,i), e2_pt(k,j,i), e2_pt(k,j+1,i));
        }
      }
    }
  } // end if 3D

  // add diffusion EMFs (e.g. Ohmic resistivity), computed at second-order accuracy;
  // consistent with the O(dx^2) treatment of the explicit diffusive fluxes
  if (!STS_ENABLED)
    if (fdif.field_diffusion_defined) fdif.AddEMF(fdif.e_oa, e);

  return;
}
