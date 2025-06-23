//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file llf_mhd_rel_no_transform.cpp
//! \brief Implements local Lax-Friedrichs Riemann solver for relativistic MHD in pure GR.

// C headers

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()

// Athena++ headers
#include "../../../athena.hpp"                   // enums, macros
#include "../../../athena_arrays.hpp"            // AthenaArray
#include "../../../coordinates/coordinates.hpp"  // Coordinates
#include "../../../eos/eos.hpp"                  // EquationOfState
#include "../../../mesh/mesh.hpp"                // MeshBlock
#include "../../hydro.hpp"

//----------------------------------------------------------------------------------------
// Riemann solver
// Inputs:
//   k,j: x3- and x2-indices
//   il,iu: lower and upper x1-indices
//   ivx: type of interface (IVX for x1, IVY for x2, IVZ for x3)
//   bb: 3D array of normal magnetic fields
//   prim_l,prim_r: 1D arrays of left and right primitive states
//   dxw: 1D array of mesh spacing in the x-direction
// Outputs:
//   flux: 3D array of hydrodynamical fluxes across interfaces
//   ey,ez: 3D arrays of magnetic fluxes (electric fields) across interfaces
//   wct: 3D array of weighting factors for CT

void Hydro::RiemannSolver(const int k, const int j, const int il, const int iu,
                          const int ivx, const AthenaArray<Real> &bb,
                          AthenaArray<Real> &prim_l, AthenaArray<Real> &prim_r,
                          AthenaArray<Real> &flux,
                          AthenaArray<Real> &ey, AthenaArray<Real> &ez,
                          AthenaArray<Real> &wct, const AthenaArray<Real> &dxw) {
  // Calculate cyclic permutations of indices
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  // Extract ratio of specific heats
  const Real gamma_adi = pmy_block->peos->GetGamma();
  Real dt = pmy_block->pmy_mesh->dt;

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
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    // Extract metric
    const Real
        &g_00 = g_(I00,i), &g_01 = g_(I01,i), &g_02 = g_(I02,i), &g_03 = g_(I03,i),
        &g_10 = g_(I01,i), &g_11 = g_(I11,i), &g_12 = g_(I12,i), &g_13 = g_(I13,i),
        &g_20 = g_(I02,i), &g_21 = g_(I12,i), &g_22 = g_(I22,i), &g_23 = g_(I23,i),
        &g_30 = g_(I03,i), &g_31 = g_(I13,i), &g_32 = g_(I23,i), &g_33 = g_(I33,i);
    const Real
        &g00 = gi_(I00,i), &g01 = gi_(I01,i), &g02 = gi_(I02,i), &g03 = gi_(I03,i),
        &g10 = gi_(I01,i), &g11 = gi_(I11,i), &g12 = gi_(I12,i), &g13 = gi_(I13,i),
        &g20 = gi_(I02,i), &g21 = gi_(I12,i), &g22 = gi_(I22,i), &g23 = gi_(I23,i),
        &g30 = gi_(I03,i), &g31 = gi_(I13,i), &g32 = gi_(I23,i), &g33 = gi_(I33,i);
    Real alpha = std::sqrt(-1.0/g00);
    Real gii, g0i;
    switch (ivx) {
      case IVX:
        gii = g11;
        g0i = g01;
        break;
      case IVY:
        gii = g22;
        g0i = g02;
        break;
      case IVZ:
        gii = g33;
        g0i = g03;
        break;
    }

    // Extract left primitives
    const Real &rho_l = prim_l(IDN,i);
    const Real &pgas_l = prim_l(IPR,i);
    const Real &uu1_l = prim_l(IVX,i);
    const Real &uu2_l = prim_l(IVY,i);
    const Real &uu3_l = prim_l(IVZ,i);
    Real bb1_l, bb2_l, bb3_l;
    switch (ivx) {
      case IVX:
        bb1_l = bb(k,j,i);
        bb2_l = prim_l(IBY,i);
        bb3_l = prim_l(IBZ,i);
        break;
      case IVY:
        bb2_l = bb(k,j,i);
        bb3_l = prim_l(IBY,i);
        bb1_l = prim_l(IBZ,i);
        break;
      case IVZ:
        bb3_l = bb(k,j,i);
        bb1_l = prim_l(IBY,i);
        bb2_l = prim_l(IBZ,i);
        break;
    }

    // Extract right primitives
    const Real &rho_r = prim_r(IDN,i);
    const Real &pgas_r = prim_r(IPR,i);
    const Real &uu1_r = prim_r(IVX,i);
    const Real &uu2_r = prim_r(IVY,i);
    const Real &uu3_r = prim_r(IVZ,i);
    Real bb1_r, bb2_r, bb3_r;
    switch (ivx) {
      case IVX:
        bb1_r = bb(k,j,i);
        bb2_r = prim_r(IBY,i);
        bb3_r = prim_r(IBZ,i);
        break;
      case IVY:
        bb2_r = bb(k,j,i);
        bb3_r = prim_r(IBY,i);
        bb1_r = prim_r(IBZ,i);
        break;
      case IVZ:
        bb3_r = bb(k,j,i);
        bb1_r = prim_r(IBY,i);
        bb2_r = prim_r(IBZ,i);
        break;
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
    Real wgas_l = rho_l + gamma_adi/(gamma_adi-1.0) * pgas_l;
    pmy_block->peos->FastMagnetosonicSpeedsGR(wgas_l, pgas_l, ucon_l[0], ucon_l[ivx],
                                              b_sq_l, g00, g0i, gii,
                                              &lambda_p_l, &lambda_m_l);

    // Calculate wavespeeds in right state
    Real lambda_p_r, lambda_m_r;
    Real wgas_r = rho_r + gamma_adi/(gamma_adi-1.0) * pgas_r;
    pmy_block->peos->FastMagnetosonicSpeedsGR(wgas_r, pgas_r, ucon_r[0], ucon_r[ivx],
                                              b_sq_r, g00, g0i, gii,
                                              &lambda_p_r, &lambda_m_r);

    // Calculate extremal wavespeed
    Real lambda_l = std::min(lambda_m_l, lambda_m_r);
    Real lambda_r = std::max(lambda_p_l, lambda_p_r);
    Real lambda = std::max(lambda_r, -lambda_l);

    // Calculate conserved quantities in L region
    // (rho u^0, T^0_\mu, and B^j = *F^{j0}, where j != ivx)
    Real cons_l[NWAVE];
    Real wtot_l = wgas_l + b_sq_l;
    Real ptot_l = pgas_l + 0.5*b_sq_l;
    cons_l[IDN] = rho_l * ucon_l[0];
    cons_l[IEN] = wtot_l * ucon_l[0] * ucov_l[0] - bcon_l[0] * bcov_l[0] + ptot_l;
    cons_l[IVX] = wtot_l * ucon_l[0] * ucov_l[1] - bcon_l[0] * bcov_l[1];
    cons_l[IVY] = wtot_l * ucon_l[0] * ucov_l[2] - bcon_l[0] * bcov_l[2];
    cons_l[IVZ] = wtot_l * ucon_l[0] * ucov_l[3] - bcon_l[0] * bcov_l[3];
    cons_l[IBY] = bcon_l[ivy] * ucon_l[0] - bcon_l[0] * ucon_l[ivy];
    cons_l[IBZ] = bcon_l[ivz] * ucon_l[0] - bcon_l[0] * ucon_l[ivz];

    // Calculate fluxes in L region
    // (rho u^i, T^i_\mu, and *F^{ji}, where i = ivx and j != ivx)
    Real flux_l[NWAVE];
    flux_l[IDN] = rho_l * ucon_l[ivx];
    flux_l[IEN] = wtot_l * ucon_l[ivx] * ucov_l[0] - bcon_l[ivx] * bcov_l[0];
    flux_l[IVX] = wtot_l * ucon_l[ivx] * ucov_l[1] - bcon_l[ivx] * bcov_l[1];
    flux_l[IVY] = wtot_l * ucon_l[ivx] * ucov_l[2] - bcon_l[ivx] * bcov_l[2];
    flux_l[IVZ] = wtot_l * ucon_l[ivx] * ucov_l[3] - bcon_l[ivx] * bcov_l[3];
    flux_l[ivx] += ptot_l;
    flux_l[IBY] = bcon_l[ivy] * ucon_l[ivx] - bcon_l[ivx] * ucon_l[ivy];
    flux_l[IBZ] = bcon_l[ivz] * ucon_l[ivx] - bcon_l[ivx] * ucon_l[ivz];

    // Calculate conserved quantities in R region
    // (rho u^0, T^0_\mu, and B^j = *F^{j0}, where j != ivx)
    Real cons_r[NWAVE];
    Real wtot_r = wgas_r + b_sq_r;
    Real ptot_r = pgas_r + 0.5*b_sq_r;
    cons_r[IDN] = rho_r * ucon_r[0];
    cons_r[IEN] = wtot_r * ucon_r[0] * ucov_r[0] - bcon_r[0] * bcov_r[0] + ptot_r;
    cons_r[IVX] = wtot_r * ucon_r[0] * ucov_r[1] - bcon_r[0] * bcov_r[1];
    cons_r[IVY] = wtot_r * ucon_r[0] * ucov_r[2] - bcon_r[0] * bcov_r[2];
    cons_r[IVZ] = wtot_r * ucon_r[0] * ucov_r[3] - bcon_r[0] * bcov_r[3];
    cons_r[IBY] = bcon_r[ivy] * ucon_r[0] - bcon_r[0] * ucon_r[ivy];
    cons_r[IBZ] = bcon_r[ivz] * ucon_r[0] - bcon_r[0] * ucon_r[ivz];

    // Calculate fluxes in R region
    // (rho u^i, T^i_\mu, and *F^{ji}, where i = ivx and j != ivx)
    Real flux_r[NWAVE];
    flux_r[IDN] = rho_r * ucon_r[ivx];
    flux_r[IEN] = wtot_r * ucon_r[ivx] * ucov_r[0] - bcon_r[ivx] * bcov_r[0];
    flux_r[IVX] = wtot_r * ucon_r[ivx] * ucov_r[1] - bcon_r[ivx] * bcov_r[1];
    flux_r[IVY] = wtot_r * ucon_r[ivx] * ucov_r[2] - bcon_r[ivx] * bcov_r[2];
    flux_r[IVZ] = wtot_r * ucon_r[ivx] * ucov_r[3] - bcon_r[ivx] * bcov_r[3];
    flux_r[ivx] += ptot_r;
    flux_r[IBY] = bcon_r[ivy] * ucon_r[ivx] - bcon_r[ivx] * ucon_r[ivy];
    flux_r[IBZ] = bcon_r[ivz] * ucon_r[ivx] - bcon_r[ivx] * ucon_r[ivz];

    // Set fluxes
    for (int n = 0; n < NHYDRO; ++n) {
      flux(n,k,j,i) =
          0.5 * (flux_l[n] + flux_r[n] - lambda * (cons_r[n] - cons_l[n]));
    }
    ey(k,j,i) =
        -0.5 * (flux_l[IBY] + flux_r[IBY] - lambda * (cons_r[IBY] - cons_l[IBY]));
    ez(k,j,i) =
        0.5 * (flux_l[IBZ] + flux_r[IBZ] - lambda * (cons_r[IBZ] - cons_l[IBZ]));

    wct(k,j,i) =
        GetWeightForCT(flux(IDN,k,j,i), prim_l(IDN,i), prim_r(IDN,i), dxw(i), dt);
  }

  return;
}
