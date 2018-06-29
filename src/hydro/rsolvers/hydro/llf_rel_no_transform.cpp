//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file llf_rel_no_transform.cpp
//  \brief Implements local Lax-Friedrichs Riemann solver for relativistic hydrodynamics
//  in pure GR.

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
//   ivx: type of interface (IVX for x1, IVY for x2, IVZ for x3)
//   bb: 3D array of normal magnetic fields (not used)
//   prim_l,prim_r: 3D arrays of left and right primitive states
// Outputs:
//   flux: 3D array of hydrodynamical fluxes across interfaces
//   ey,ez: 3D arrays of magnetic fluxes (electric fields) across interfaces (not used)
// Notes:
//   implements LLF algorithm similar to that of fluxcalc() in step_ch.c in Harm
//   cf. LLFNonTransforming() in llf_rel.cpp

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
          Real g_00 = g_(I00,ipm), g_01 = g_(I01,ipm), g_02 = g_(I02,ipm),
            g_03 = g_(I03,ipm), g_10 = g_(I01,ipm), g_11 = g_(I11,ipm),
            g_12 = g_(I12,ipm), g_13 = g_(I13,ipm), g_20 = g_(I02,ipm),
            g_21 = g_(I12,ipm), g_22 = g_(I22,ipm), g_23 = g_(I23,ipm),
            g_30 = g_(I03,ipm), g_31 = g_(I13,ipm), g_32 = g_(I23,ipm),
            g_33 = g_(I33,ipm);
          Real g00 = gi_(I00,ipm), g01 = gi_(I01,ipm), g02 = gi_(I02,ipm),
            g03 = gi_(I03,ipm), g10 = gi_(I01,ipm), g11 = gi_(I11,ipm),
            g12 = gi_(I12,ipm), g13 = gi_(I13,ipm), g20 = gi_(I02,ipm),
            g21 = gi_(I12,ipm), g22 = gi_(I22,ipm), g23 = gi_(I23,ipm),
            g30 = gi_(I03,ipm), g31 = gi_(I13,ipm), g32 = gi_(I23,ipm),
            g33 = gi_(I33,ipm);
          Real alpha = std::sqrt(-1.0/g00);
          Real gii, g0i;
          if (ivx==IVX) {
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

          // Extract right primitives
          Real rho_r = prim_r(IDN,k,j,ipm);
          Real pgas_r = prim_r(IPR,k,j,ipm);
          Real uu1_r = prim_r(IVX,k,j,ipm);
          Real uu2_r = prim_r(IVY,k,j,ipm);
          Real uu3_r = prim_r(IVZ,k,j,ipm);

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

          // Calculate wavespeeds in left state
          Real lambda_p_l, lambda_m_l;
          Real wgas_l = rho_l + gamma_prime * pgas_l;
          pmy_block->peos->SoundSpeedsGR(wgas_l, pgas_l, ucon_l[0], ucon_l[ivx], g00, g0i,
                                         gii, &lambda_p_l, &lambda_m_l);

          // Calculate wavespeeds in right state
          Real lambda_p_r, lambda_m_r;
          Real wgas_r = rho_r + gamma_prime * pgas_r;
          pmy_block->peos->SoundSpeedsGR(wgas_r, pgas_r, ucon_r[0], ucon_r[ivx], g00, g0i,
                                         gii, &lambda_p_r, &lambda_m_r);

          // Calculate extremal wavespeed
          Real lambda_l = std::min(lambda_m_l, lambda_m_r);
          Real lambda_r = std::max(lambda_p_l, lambda_p_r);
          Real lambda = std::max(lambda_r, -lambda_l);

          // Calculate conserved quantities in L region (rho u^0 and T^0_\mu)
          cons_l[IDN][m] = rho_l * ucon_l[0];
          cons_l[IEN][m] = wgas_l * ucon_l[0] * ucov_l[0] + pgas_l;
          cons_l[IVX][m] = wgas_l * ucon_l[0] * ucov_l[1];
          cons_l[IVY][m] = wgas_l * ucon_l[0] * ucov_l[2];
          cons_l[IVZ][m] = wgas_l * ucon_l[0] * ucov_l[3];

          // Calculate fluxes in L region (rho u^i and T^i_\mu, where i = ivx)
          flux_l[IDN][m] = rho_l * ucon_l[ivx];
          flux_l[IEN][m] = wgas_l * ucon_l[ivx] * ucov_l[0];
          flux_l[IVX][m] = wgas_l * ucon_l[ivx] * ucov_l[1];
          flux_l[IVY][m] = wgas_l * ucon_l[ivx] * ucov_l[2];
          flux_l[IVZ][m] = wgas_l * ucon_l[ivx] * ucov_l[3];
          flux_l[ivx][m] += pgas_l;

          // Calculate conserved quantities in R region (rho u^0 and T^0_\mu)
          cons_r[IDN][m] = rho_r * ucon_r[0];
          cons_r[IEN][m] = wgas_r * ucon_r[0] * ucov_r[0] + pgas_r;
          cons_r[IVX][m] = wgas_r * ucon_r[0] * ucov_r[1];
          cons_r[IVY][m] = wgas_r * ucon_r[0] * ucov_r[2];
          cons_r[IVZ][m] = wgas_r * ucon_r[0] * ucov_r[3];

          // Calculate fluxes in R region (rho u^i and T^i_\mu, where i = ivx)
          flux_r[IDN][m] = rho_r * ucon_r[ivx];
          flux_r[IEN][m] = wgas_r * ucon_r[ivx] * ucov_r[0];
          flux_r[IVX][m] = wgas_r * ucon_r[ivx] * ucov_r[1];
          flux_r[IVY][m] = wgas_r * ucon_r[ivx] * ucov_r[2];
          flux_r[IVZ][m] = wgas_r * ucon_r[ivx] * ucov_r[3];
          flux_r[ivx][m] += pgas_r;

          // Set fluxes
          for (int n = 0; n < NHYDRO; ++n) {
            flux(n,k,j,ipm) =
              0.5 * (flux_l[n][m] + flux_r[n][m] -
                     lambda * (cons_r[n][m] - cons_l[n][m]));
          }
        }
      }
    }
  }
  return;
}
