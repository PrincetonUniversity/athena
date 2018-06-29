//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hlle_mhd_rel_no_transform.cpp
//  \brief Implements HLLE Riemann solver for relativistic MHD in pure GR.

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()

// Athena++ headers
#include "../../hydro.hpp"
#include "../../../athena.hpp"                   // enums, macros
#include "../../../athena_arrays.hpp"            // AthenaArray
#include "../../../coordinates/coordinates.hpp"  // Coordinates
#include "../../../eos/eos.hpp"                  // EquationOfState
#include "../../../mesh/mesh.hpp"                // MeshBlock

//----------------------------------------------------------------------------------------
// Riemann solver
// Inputs:
//   kl,ku,jl,ju,il,iu: lower and upper x1-, x2-, and x3-indices
//   bb: 3D array of normal magnetic fields
//   prim_l,prim_r: 3D arrays of left and right primitive states
// Outputs:
//   flux: 3D array of hydrodynamical fluxes across interfaces
//   ey,ez: 3D arrays of magnetic fluxes (electric fields) across interfaces
// Notes:
//   implements HLLE algorithm similar to that of fluxcalc() in step_ch.c in Harm
//   cf. HLLENonTransforming() in hlle_mhd_rel.cpp and hlld_rel.cpp

void Hydro::RiemannSolver(const int kl, const int ku, const int jl, const int ju,
    const int il, const int iu, const int ivx, const AthenaArray<Real> &bb,
    AthenaArray<Real> &prim_l, AthenaArray<Real> &prim_r, AthenaArray<Real> &flux,
    AthenaArray<Real> &ey, AthenaArray<Real> &ez) {
  // Calculate cyclic permutations of indices
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  // Extract ratio of specific heats
  const Real gamma_adi = pmy_block->peos->GetGamma();
  const Real gamma_prime = gamma_adi/(gamma_adi-1.0);

  Real cons_l[NWAVE][SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  Real flux_l[NWAVE][SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  Real cons_r[NWAVE][SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  Real flux_r[NWAVE][SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  Real flux_hll[NWAVE][SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  Real flux_interface[NWAVE][SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));

  // Go through 1D arrays of interfaces
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {

      // Get metric components
      switch (ivx) {
        case IVX:
          pmy_block->pcoord->Face1Metric(k, j, il, iu, g_, gi_);
          break;
        case IVY:
          pmy_block->pcoord->Face2Metric(k, j, il, iu, g_, gi_);
          break;
        case IVZ:
          pmy_block->pcoord->Face3Metric(k, j, il, iu, g_, gi_);
          break;
      }

      // Go through each interface
      for (int i = il; i <= iu; i+=SIMD_WIDTH) {
#pragma omp simd simdlen(SIMD_WIDTH)
        for (int m=0; m<std::min(SIMD_WIDTH, iu-i+1); m++) {
          int ipm = i+m;
          // Extract metric
          Real g_00 = g_(I00,ipm), g_01 = g_(I01,ipm), g_02 = g_(I02,ipm), g_03 = g_(I03,ipm),
            g_10 = g_(I01,ipm), g_11 = g_(I11,ipm), g_12 = g_(I12,ipm), g_13 = g_(I13,ipm),
            g_20 = g_(I02,ipm), g_21 = g_(I12,ipm), g_22 = g_(I22,ipm), g_23 = g_(I23,ipm),
            g_30 = g_(I03,ipm), g_31 = g_(I13,ipm), g_32 = g_(I23,ipm), g_33 = g_(I33,ipm);
          Real g00 = gi_(I00,ipm), g01 = gi_(I01,ipm), g02 = gi_(I02,ipm), g03 = gi_(I03,ipm),
            g10 = gi_(I01,ipm), g11 = gi_(I11,ipm), g12 = gi_(I12,ipm), g13 = gi_(I13,ipm),
            g20 = gi_(I02,ipm), g21 = gi_(I12,ipm), g22 = gi_(I22,ipm), g23 = gi_(I23,ipm),
            g30 = gi_(I03,ipm), g31 = gi_(I13,ipm), g32 = gi_(I23,ipm), g33 = gi_(I33,ipm);
          Real alpha = std::sqrt(-1.0/g00);

          Real gii, g0i;
          if(ivx==IVX) {
            gii = g11;
            g0i = g01;
          } else if (ivx==IVY) {
            gii = g22;
            g0i = g02;
          } else if (ivx==IVZ) {
            gii = g33;
            g0i = g03;
          }

          // Extract left primitives
          Real rho_l = prim_l(IDN,k,j,ipm);
          Real pgas_l = prim_l(IPR,k,j,ipm);
          Real uu1_l = prim_l(IVX,k,j,ipm);
          Real uu2_l = prim_l(IVY,k,j,ipm);
          Real uu3_l = prim_l(IVZ,k,j,ipm);
          Real bb1_l, bb2_l, bb3_l;
          if (ivx==IVX) {
            bb1_l = bb(k,j,ipm);
            bb2_l = prim_l(IBY,k,j,ipm);
            bb3_l = prim_l(IBZ,k,j,ipm);
          } else if (ivx==IVY) {
            bb2_l = bb(k,j,ipm);
            bb3_l = prim_l(IBY,k,j,ipm);
            bb1_l = prim_l(IBZ,k,j,ipm);
          } else if (ivx==IVZ) {
            bb3_l = bb(k,j,ipm);
            bb1_l = prim_l(IBY,k,j,ipm);
            bb2_l = prim_l(IBZ,k,j,ipm);
          }

          // Extract right primitives
          Real rho_r = prim_r(IDN,k,j,ipm);
          Real pgas_r = prim_r(IPR,k,j,ipm);
          Real uu1_r = prim_r(IVX,k,j,ipm);
          Real uu2_r = prim_r(IVY,k,j,ipm);
          Real uu3_r = prim_r(IVZ,k,j,ipm);
          Real bb1_r, bb2_r, bb3_r;
          if (ivx==IVZ) {
            bb1_r = bb(k,j,ipm);
            bb2_r = prim_r(IBY,k,j,ipm);
            bb3_r = prim_r(IBZ,k,j,ipm);
          } else if (ivx==IVY) {
            bb2_r = bb(k,j,ipm);
            bb3_r = prim_r(IBY,k,j,ipm);
            bb1_r = prim_r(IBZ,k,j,ipm);
          } else if (ivx==IVZ) {
            bb3_r = bb(k,j,ipm);
            bb1_r = prim_r(IBY,k,j,ipm);
            bb2_r = prim_r(IBZ,k,j,ipm);
          }

          // Calculate 4-velocity in left state
          Real ucon_l[4], ucov_l[4];
          Real tmp = g_11*SQR(uu1_l) + 2.0*g_12*uu1_l*uu2_l + 2.0*g_13*uu1_l*uu3_l
            + g_22*SQR(uu2_l) + 2.0*g_23*uu2_l*uu3_l
            + g_33*SQR(uu3_l);
          Real gamma_l = std::sqrt(1.0 + tmp);
          ucon_l[0] = gamma_l / alpha;
          ucon_l[1] = uu1_l - alpha * gamma_l * g01;
          ucon_l[2] = uu2_l - alpha * gamma_l * g02;
          ucon_l[3] = uu3_l - alpha * gamma_l * g03;
          ucov_l[0] = g_00*ucon_l[0] + g_01*ucon_l[1] + g_02*ucon_l[2] + g_03*ucon_l[3];
          ucov_l[1] = g_10*ucon_l[0] + g_11*ucon_l[1] + g_12*ucon_l[2] + g_13*ucon_l[3];
          ucov_l[2] = g_20*ucon_l[0] + g_21*ucon_l[1] + g_22*ucon_l[2] + g_23*ucon_l[3];
          ucov_l[3] = g_30*ucon_l[0] + g_31*ucon_l[1] + g_32*ucon_l[2] + g_33*ucon_l[3];

          // Calculate 4-velocity in right state
          Real ucon_r[4], ucov_r[4];
          tmp = g_11*SQR(uu1_r) + 2.0*g_12*uu1_r*uu2_r + 2.0*g_13*uu1_r*uu3_r
            + g_22*SQR(uu2_r) + 2.0*g_23*uu2_r*uu3_r
            + g_33*SQR(uu3_r);
          Real gamma_r = std::sqrt(1.0 + tmp);
          ucon_r[0] = gamma_r / alpha;
          ucon_r[1] = uu1_r - alpha * gamma_r * g01;
          ucon_r[2] = uu2_r - alpha * gamma_r * g02;
          ucon_r[3] = uu3_r - alpha * gamma_r * g03;
          ucov_r[0] = g_00*ucon_r[0] + g_01*ucon_r[1] + g_02*ucon_r[2] + g_03*ucon_r[3];
          ucov_r[1] = g_10*ucon_r[0] + g_11*ucon_r[1] + g_12*ucon_r[2] + g_13*ucon_r[3];
          ucov_r[2] = g_20*ucon_r[0] + g_21*ucon_r[1] + g_22*ucon_r[2] + g_23*ucon_r[3];
          ucov_r[3] = g_30*ucon_r[0] + g_31*ucon_r[1] + g_32*ucon_r[2] + g_33*ucon_r[3];

          // Calculate 4-magnetic field in left state
          Real bcon_l[4], bcov_l[4];
          bcon_l[0] = ucon_l[0] * (g_01*bb1_l + g_02*bb2_l + g_03*bb3_l)
            + ucon_l[1] * (g_11*bb1_l + g_12*bb2_l + g_13*bb3_l)
            + ucon_l[2] * (g_21*bb1_l + g_22*bb2_l + g_23*bb3_l)
            + ucon_l[3] * (g_31*bb1_l + g_32*bb2_l + g_33*bb3_l);
          bcon_l[1] = (bb1_l + bcon_l[0] * ucon_l[1]) / ucon_l[0];
          bcon_l[2] = (bb2_l + bcon_l[0] * ucon_l[2]) / ucon_l[0];
          bcon_l[3] = (bb3_l + bcon_l[0] * ucon_l[3]) / ucon_l[0];
          bcov_l[0] = g_00*bcon_l[0] + g_01*bcon_l[1] + g_02*bcon_l[2] + g_03*bcon_l[3];
          bcov_l[1] = g_10*bcon_l[0] + g_11*bcon_l[1] + g_12*bcon_l[2] + g_13*bcon_l[3];
          bcov_l[2] = g_20*bcon_l[0] + g_21*bcon_l[1] + g_22*bcon_l[2] + g_23*bcon_l[3];
          bcov_l[3] = g_30*bcon_l[0] + g_31*bcon_l[1] + g_32*bcon_l[2] + g_33*bcon_l[3];
          Real b_sq_l = bcon_l[0]*bcov_l[0] + bcon_l[1]*bcov_l[1] + bcon_l[2]*bcov_l[2]
            + bcon_l[3]*bcov_l[3];

          // Calculate 4-magnetic field in right state
          Real bcon_r[4], bcov_r[4];
          bcon_r[0] = ucon_r[0] * (g_01*bb1_r + g_02*bb2_r + g_03*bb3_r)
            + ucon_r[1] * (g_11*bb1_r + g_12*bb2_r + g_13*bb3_r)
            + ucon_r[2] * (g_21*bb1_r + g_22*bb2_r + g_23*bb3_r)
            + ucon_r[3] * (g_31*bb1_r + g_32*bb2_r + g_33*bb3_r);
          bcon_r[1] = (bb1_r + bcon_r[0] * ucon_r[1]) / ucon_r[0];
          bcon_r[2] = (bb2_r + bcon_r[0] * ucon_r[2]) / ucon_r[0];
          bcon_r[3] = (bb3_r + bcon_r[0] * ucon_r[3]) / ucon_r[0];
          bcov_r[0] = g_00*bcon_r[0] + g_01*bcon_r[1] + g_02*bcon_r[2] + g_03*bcon_r[3];
          bcov_r[1] = g_10*bcon_r[0] + g_11*bcon_r[1] + g_12*bcon_r[2] + g_13*bcon_r[3];
          bcov_r[2] = g_20*bcon_r[0] + g_21*bcon_r[1] + g_22*bcon_r[2] + g_23*bcon_r[3];
          bcov_r[3] = g_30*bcon_r[0] + g_31*bcon_r[1] + g_32*bcon_r[2] + g_33*bcon_r[3];
          Real b_sq_r = bcon_r[0]*bcov_r[0] + bcon_r[1]*bcov_r[1] + bcon_r[2]*bcov_r[2]
            + bcon_r[3]*bcov_r[3];

          // Calculate wavespeeds in left state
          Real lambda_p_l, lambda_m_l;
          Real wgas_l = rho_l + gamma_prime * pgas_l;
          pmy_block->peos->FastMagnetosonicSpeedsGR(wgas_l, pgas_l, ucon_l[0], ucon_l[ivx],
                                                    b_sq_l, g00, g0i, gii, &lambda_p_l, &lambda_m_l);

          // Calculate wavespeeds in right state
          Real lambda_p_r, lambda_m_r;
          Real wgas_r = rho_r + gamma_prime * pgas_r;
          pmy_block->peos->FastMagnetosonicSpeedsGR(wgas_r, pgas_r, ucon_r[0], ucon_r[ivx],
                                                    b_sq_r, g00, g0i, gii, &lambda_p_r, &lambda_m_r);

          // Calculate extremal wavespeeds
          Real lambda_l = std::min(lambda_m_l, lambda_m_r);
          Real lambda_r = std::max(lambda_p_l, lambda_p_r);

          // Calculate conserved quantities in L region
          // (rho u^0, T^0_\mu, and B^j = *F^{j0}, where j != ivx)
          Real wtot_l = wgas_l + b_sq_l;
          Real ptot_l = pgas_l + 0.5*b_sq_l;
          cons_l[IDN][m] = rho_l * ucon_l[0];
          cons_l[IEN][m] = wtot_l * ucon_l[0] * ucov_l[0] - bcon_l[0] * bcov_l[0] + ptot_l;
          cons_l[IVX][m] = wtot_l * ucon_l[0] * ucov_l[1] - bcon_l[0] * bcov_l[1];
          cons_l[IVY][m] = wtot_l * ucon_l[0] * ucov_l[2] - bcon_l[0] * bcov_l[2];
          cons_l[IVZ][m] = wtot_l * ucon_l[0] * ucov_l[3] - bcon_l[0] * bcov_l[3];
          cons_l[IBY][m] = bcon_l[ivy] * ucon_l[0] - bcon_l[0] * ucon_l[ivy];
          cons_l[IBZ][m] = bcon_l[ivz] * ucon_l[0] - bcon_l[0] * ucon_l[ivz];

          // Calculate fluxes in L region
          // (rho u^i, T^i_\mu, and *F^{ji}, where i = ivx and j != ivx)
          flux_l[IDN][m] = rho_l * ucon_l[ivx];
          flux_l[IEN][m] = wtot_l * ucon_l[ivx] * ucov_l[0] - bcon_l[ivx] * bcov_l[0];
          flux_l[IVX][m] = wtot_l * ucon_l[ivx] * ucov_l[1] - bcon_l[ivx] * bcov_l[1];
          flux_l[IVY][m] = wtot_l * ucon_l[ivx] * ucov_l[2] - bcon_l[ivx] * bcov_l[2];
          flux_l[IVZ][m] = wtot_l * ucon_l[ivx] * ucov_l[3] - bcon_l[ivx] * bcov_l[3];
          flux_l[ivx][m] += ptot_l;
          flux_l[IBY][m] = bcon_l[ivy] * ucon_l[ivx] - bcon_l[ivx] * ucon_l[ivy];
          flux_l[IBZ][m] = bcon_l[ivz] * ucon_l[ivx] - bcon_l[ivx] * ucon_l[ivz];

        // Calculate conserved quantities in R region
        // (rho u^0, T^0_\mu, and B^j = *F^{j0}, where j != ivx)
          Real wtot_r = wgas_r + b_sq_r;
          Real ptot_r = pgas_r + 0.5*b_sq_r;
          cons_r[IDN][m] = rho_r * ucon_r[0];
          cons_r[IEN][m] = wtot_r * ucon_r[0] * ucov_r[0] - bcon_r[0] * bcov_r[0] + ptot_r;
          cons_r[IVX][m] = wtot_r * ucon_r[0] * ucov_r[1] - bcon_r[0] * bcov_r[1];
          cons_r[IVY][m] = wtot_r * ucon_r[0] * ucov_r[2] - bcon_r[0] * bcov_r[2];
          cons_r[IVZ][m] = wtot_r * ucon_r[0] * ucov_r[3] - bcon_r[0] * bcov_r[3];
          cons_r[IBY][m] = bcon_r[ivy] * ucon_r[0] - bcon_r[0] * ucon_r[ivy];
          cons_r[IBZ][m] = bcon_r[ivz] * ucon_r[0] - bcon_r[0] * ucon_r[ivz];

          // Calculate fluxes in R region
          // (rho u^i, T^i_\mu, and *F^{ji}, where i = ivx and j != ivx)
          flux_r[IDN][m] = rho_r * ucon_r[ivx];
          flux_r[IEN][m] = wtot_r * ucon_r[ivx] * ucov_r[0] - bcon_r[ivx] * bcov_r[0];
          flux_r[IVX][m] = wtot_r * ucon_r[ivx] * ucov_r[1] - bcon_r[ivx] * bcov_r[1];
          flux_r[IVY][m] = wtot_r * ucon_r[ivx] * ucov_r[2] - bcon_r[ivx] * bcov_r[2];
          flux_r[IVZ][m] = wtot_r * ucon_r[ivx] * ucov_r[3] - bcon_r[ivx] * bcov_r[3];
          flux_r[ivx][m] += ptot_r;
          flux_r[IBY][m] = bcon_r[ivy] * ucon_r[ivx] - bcon_r[ivx] * ucon_r[ivy];
          flux_r[IBZ][m] = bcon_r[ivz] * ucon_r[ivx] - bcon_r[ivx] * ucon_r[ivz];

          Real lambda_diff_inv = 1.0 / (lambda_r-lambda_l);
          // Calculate fluxes in HLL region
          for (int n = 0; n < NWAVE; ++n) {
            flux_hll[n][m] = (lambda_r*flux_l[n][m] - lambda_l*flux_r[n][m]
                           + lambda_r*lambda_l * (cons_r[n][m] - cons_l[n][m]))
              * lambda_diff_inv;
          }

          // Determine region of waveface
          if (lambda_l >= 0.0) {  // L region
            for (int n = 0; n < NWAVE; ++n) {
              flux_interface[n][m] = flux_l[n][m];
            }
          } else if (lambda_r <= 0.0) { // R region
            for (int n = 0; n < NWAVE; ++n) {
              flux_interface[n][m] = flux_r[n][m];
            }
          } else {  // HLL region
            for (int n = 0; n < NWAVE; ++n) {
              flux_interface[n][m] = flux_hll[n][m];
            }
          }

          // Set fluxes
          for (int n = 0; n < NHYDRO; ++n) {
            flux(n,k,j,ipm) = flux_interface[n][m];
          }
          ey(k,j,ipm) = -flux_interface[IBY][m];
          ez(k,j,ipm) = flux_interface[IBZ][m];

        }
      }
    }
  }
  return;
}
