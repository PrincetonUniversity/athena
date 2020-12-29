//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hlle_mhd_rel.cpp
//! \brief Implements HLLE Riemann solver for relativistic MHD.

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

namespace {
// Declarations
void HLLETransforming(MeshBlock *pmb, const int k, const int j,
                      const int il, const int iu, const int ivx,
                      const AthenaArray<Real> &bb, AthenaArray<Real> &bb_normal,
                      AthenaArray<Real> &g, AthenaArray<Real> &gi,
                      AthenaArray<Real> &prim_l, AthenaArray<Real> &prim_r,
                      AthenaArray<Real> &cons, AthenaArray<Real> &flux,
                      AthenaArray<Real> &ey, AthenaArray<Real> &ez);
void HLLENonTransforming(MeshBlock *pmb, const int k, const int j,
                         const int il, const int iu, const AthenaArray<Real> &bb,
                         AthenaArray<Real> &g, AthenaArray<Real> &gi,
                         AthenaArray<Real> &prim_l, AthenaArray<Real> &prim_r,
                         AthenaArray<Real> &flux,
                         AthenaArray<Real> &ey, AthenaArray<Real> &ez);
} // namespace

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
// Notes:
//   prim_l, prim_r overwritten
//   tries to implement HLLE algorithm from Mignone & Bodo 2005, MNRAS 364 126 (MB2005)
//   otherwise implements HLLE algorithm similar to that of fluxcalc() in step_ch.c in
//       Harm

void Hydro::RiemannSolver(const int k, const int j, const int il, const int iu,
                          const int ivx, const AthenaArray<Real> &bb,
                          AthenaArray<Real> &prim_l, AthenaArray<Real> &prim_r,
                          AthenaArray<Real> &flux,
                          AthenaArray<Real> &ey, AthenaArray<Real> &ez,
                          AthenaArray<Real> &wct, const AthenaArray<Real> &dxw) {
  Real dt = pmy_block->pmy_mesh->dt;
  if (GENERAL_RELATIVITY && ivx == IVY && pmy_block->pcoord->IsPole(j)) {
    HLLENonTransforming(pmy_block, k, j, il, iu, bb, g_, gi_, prim_l, prim_r, flux, ey,
                        ez);
  } else {
    HLLETransforming(pmy_block, k, j, il, iu, ivx, bb, bb_normal_,
                     g_, gi_, prim_l, prim_r, cons_, flux,
                     ey, ez);
  }
  for(int i=il; i<=iu; ++i) {
    wct(k,j,i) = GetWeightForCT(flux(IDN,k,j,i), prim_l(IDN,i), prim_r(IDN,i), dxw(i),
                                dt);
  }
  return;
}

namespace {
//----------------------------------------------------------------------------------------
// Frame-transforming HLLE implementation
// Inputs:
//   pmb: pointer to MeshBlock object
//   k,j: x3- and x2-indices
//   il,iu: lower and upper x1-indices
//   ivx: type of interface (IVX for x1, IVY for x2, IVZ for x3)
//   bb: 3D array of normal magnetic fields
//   bb_normal: 1D scratch array for normal magnetic fields
//   g,gi: 1D scratch arrays for metric coefficients
//   prim_l,prim_r: 1D arrays of left and right primitive states
//   cons: 1D scratch array for conserved quantities
// Outputs:
//   flux: 3D array of hydrodynamical fluxes across interfaces
//   ey,ez: 3D arrays of magnetic fluxes (electric fields) across interfaces
// Notes:
//   prim_l, prim_r overwritten
//   implements HLLE algorithm from Mignone & Bodo 2005, MNRAS 364 126 (MB2005)
//   references Mignone & Bodo 2006, MNRAS 368 1040 (MB2006)
//   references Mignone, Ugliano, & Bodo 2009, MNRAS 393 1141 (MUB)

void HLLETransforming(MeshBlock *pmb, const int k, const int j,
                      const int il, const int iu, const int ivx,
                      const AthenaArray<Real> &bb, AthenaArray<Real> &bb_normal,
                      AthenaArray<Real> &g, AthenaArray<Real> &gi,
                      AthenaArray<Real> &prim_l, AthenaArray<Real> &prim_r,
                      AthenaArray<Real> &cons, AthenaArray<Real> &flux,
                      AthenaArray<Real> &ey, AthenaArray<Real> &ez) {
  // Calculate metric if in GR
  int i01(0), i11(0);
#if GENERAL_RELATIVITY
  {
    switch (ivx) {
      case IVX:
        pmb->pcoord->Face1Metric(k, j, il, iu, g, gi);
        i01 = I01;
        i11 = I11;
        break;
      case IVY:
        pmb->pcoord->Face2Metric(k, j, il, iu, g, gi);
        i01 = I02;
        i11 = I22;
        break;
      case IVZ:
        pmb->pcoord->Face3Metric(k, j, il, iu, g, gi);
        i01 = I03;
        i11 = I33;
        break;
    }
  }
#endif  // GENERAL_RELATIVITY

  // Transform primitives to locally flat coordinates if in GR
#if GENERAL_RELATIVITY
  {
    switch (ivx) {
      case IVX:
        pmb->pcoord->PrimToLocal1(k, j, il, iu, bb, prim_l, prim_r, bb_normal);
        break;
      case IVY:
        pmb->pcoord->PrimToLocal2(k, j, il, iu, bb, prim_l, prim_r, bb_normal);
        break;
      case IVZ:
        pmb->pcoord->PrimToLocal3(k, j, il, iu, bb, prim_l, prim_r, bb_normal);
        break;
    }
  }
#else  // SR; need to populate 1D normal B array
  {
#pragma omp simd simdlen(SIMD_WIDTH)
    for (int i=il; i<=iu; ++i) {
      bb_normal(i) = bb(k,j,i);
    }
  }
#endif  // GENERAL_RELATIVITY

  // Calculate cyclic permutations of indices
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  // Extract ratio of specific heats
  const Real gamma_adi = pmb->peos->GetGamma();

  // Go through each interface
#pragma omp simd simdlen(SIMD_WIDTH)
  for (int i=il; i<=iu; ++i) {
    // Extract left primitives
    Real rho_l = prim_l(IDN,i);
    Real pgas_l = prim_l(IPR,i);
    Real ux_l = prim_l(ivx,i);
    Real uy_l = prim_l(ivy,i);
    Real uz_l = prim_l(ivz,i);
    Real u_l[4];
    u_l[0] = std::sqrt(1.0 + SQR(ux_l) + SQR(uy_l) + SQR(uz_l));
    u_l[1] = ux_l;
    u_l[2] = uy_l;
    u_l[3] = uz_l;
    Real bb2_l = prim_l(IBY,i);
    Real bb3_l = prim_l(IBZ,i);

    // Extract right primitives
    Real rho_r = prim_r(IDN,i);
    Real pgas_r = prim_r(IPR,i);
    Real ux_r = prim_r(ivx,i);
    Real uy_r = prim_r(ivy,i);
    Real uz_r = prim_r(ivz,i);
    Real u_r[4];
    u_r[0] = std::sqrt(1.0 + SQR(ux_r) + SQR(uy_r) + SQR(uz_r));
    u_r[1] = ux_r;
    u_r[2] = uy_r;
    u_r[3] = uz_r;
    Real bb2_r = prim_r(IBY,i);
    Real bb3_r = prim_r(IBZ,i);

    // Extract normal magnetic field
    Real bb1 = bb_normal(i);

    // Calculate 4-magnetic field in left state
    Real b_l[4];
    b_l[0] = bb1*u_l[1] + bb2_l*u_l[2] + bb3_l*u_l[3];
    b_l[1] = (bb1 + b_l[0] * u_l[1]) / u_l[0];
    b_l[2] = (bb2_l + b_l[0] * u_l[2]) / u_l[0];
    b_l[3] = (bb3_l + b_l[0] * u_l[3]) / u_l[0];
    Real b_sq_l = -SQR(b_l[0]) + SQR(b_l[1]) + SQR(b_l[2]) + SQR(b_l[3]);

    // Calculate 4-magnetic field in right state
    Real b_r[4];
    b_r[0] = bb1*u_r[1] + bb2_r*u_r[2] + bb3_r*u_r[3];
    b_r[1] = (bb1 + b_r[0] * u_r[1]) / u_r[0];
    b_r[2] = (bb2_r + b_r[0] * u_r[2]) / u_r[0];
    b_r[3] = (bb3_r + b_r[0] * u_r[3]) / u_r[0];
    Real b_sq_r = -SQR(b_r[0]) + SQR(b_r[1]) + SQR(b_r[2]) + SQR(b_r[3]);

    // Calculate left wavespeeds
    Real wgas_l = rho_l + gamma_adi / (gamma_adi - 1.0) * pgas_l;
    Real lambda_m_l, lambda_p_l;
    pmb->peos->FastMagnetosonicSpeedsGR(wgas_l, pgas_l, u_l[0], u_l[1], b_sq_l, -1.0, 0.0,
                                        1.0, &lambda_p_l, &lambda_m_l);

    // Calculate right wavespeeds
    Real wgas_r = rho_r + gamma_adi / (gamma_adi - 1.0) * pgas_r;
    Real lambda_m_r, lambda_p_r;
    pmb->peos->FastMagnetosonicSpeedsGR(wgas_r, pgas_r, u_r[0], u_r[1], b_sq_r, -1.0, 0.0,
                                        1.0, &lambda_p_r, &lambda_m_r);

    // Calculate extremal wavespeeds (MB2006 55)
    Real lambda_l = std::min(lambda_m_l, lambda_m_r);
    Real lambda_r = std::max(lambda_p_l, lambda_p_r);

    // Calculate conserved quantities in L region (MUB 8)
    Real cons_l[NWAVE];
    Real wtot_l = wgas_l + b_sq_l;
    Real ptot_l = pgas_l + 0.5*b_sq_l;
    cons_l[IDN] = rho_l * u_l[0];
    cons_l[IEN] = wtot_l * u_l[0] * u_l[0] - b_l[0] * b_l[0] - ptot_l;
    cons_l[ivx] = wtot_l * u_l[1] * u_l[0] - b_l[1] * b_l[0];
    cons_l[ivy] = wtot_l * u_l[2] * u_l[0] - b_l[2] * b_l[0];
    cons_l[ivz] = wtot_l * u_l[3] * u_l[0] - b_l[3] * b_l[0];
    cons_l[IBY] = b_l[2] * u_l[0] - b_l[0] * u_l[2];
    cons_l[IBZ] = b_l[3] * u_l[0] - b_l[0] * u_l[3];

    // Calculate fluxes in L region (MUB 15)
    Real flux_l[NWAVE];
    flux_l[IDN] = rho_l * u_l[1];
    flux_l[IEN] = wtot_l * u_l[0] * u_l[1] - b_l[0] * b_l[1];
    flux_l[ivx] = wtot_l * u_l[1] * u_l[1] - b_l[1] * b_l[1] + ptot_l;
    flux_l[ivy] = wtot_l * u_l[2] * u_l[1] - b_l[2] * b_l[1];
    flux_l[ivz] = wtot_l * u_l[3] * u_l[1] - b_l[3] * b_l[1];
    flux_l[IBY] = b_l[2] * u_l[1] - b_l[1] * u_l[2];
    flux_l[IBZ] = b_l[3] * u_l[1] - b_l[1] * u_l[3];

    // Calculate conserved quantities in R region (MUB 8)
    Real cons_r[NWAVE];
    Real wtot_r = wgas_r + b_sq_r;
    Real ptot_r = pgas_r + 0.5*b_sq_r;
    cons_r[IDN] = rho_r * u_r[0];
    cons_r[IEN] = wtot_r * u_r[0] * u_r[0] - b_r[0] * b_r[0] - ptot_r;
    cons_r[ivx] = wtot_r * u_r[1] * u_r[0] - b_r[1] * b_r[0];
    cons_r[ivy] = wtot_r * u_r[2] * u_r[0] - b_r[2] * b_r[0];
    cons_r[ivz] = wtot_r * u_r[3] * u_r[0] - b_r[3] * b_r[0];
    cons_r[IBY] = b_r[2] * u_r[0] - b_r[0] * u_r[2];
    cons_r[IBZ] = b_r[3] * u_r[0] - b_r[0] * u_r[3];

    // Calculate fluxes in R region (MUB 15)
    Real flux_r[NWAVE];
    flux_r[IDN] = rho_r * u_r[1];
    flux_r[IEN] = wtot_r * u_r[0] * u_r[1] - b_r[0] * b_r[1];
    flux_r[ivx] = wtot_r * u_r[1] * u_r[1] - b_r[1] * b_r[1] + ptot_r;
    flux_r[ivy] = wtot_r * u_r[2] * u_r[1] - b_r[2] * b_r[1];
    flux_r[ivz] = wtot_r * u_r[3] * u_r[1] - b_r[3] * b_r[1];
    flux_r[IBY] = b_r[2] * u_r[1] - b_r[1] * u_r[2];
    flux_r[IBZ] = b_r[3] * u_r[1] - b_r[1] * u_r[3];

    // Calculate conserved quantities in HLL region in GR (MB2005 9)
    Real lambda_diff_inv = 1.0 / (lambda_r-lambda_l);
    Real cons_hll[NWAVE];
    if (GENERAL_RELATIVITY) {
      for (int n = 0; n < NWAVE; ++n) {
        cons_hll[n] = (lambda_r*cons_r[n] - lambda_l*cons_l[n] + flux_l[n] - flux_r[n])
                      * lambda_diff_inv;
      }
    }

    // Calculate fluxes in HLL region (MB2005 11)
    Real flux_hll[NWAVE];
    for (int n = 0; n < NWAVE; ++n) {
      flux_hll[n] = (lambda_r*flux_l[n] - lambda_l*flux_r[n]
                     + lambda_l*lambda_r * (cons_r[n] - cons_l[n])) * lambda_diff_inv;
    }

    // Calculate interface velocity
    Real v_interface = 0.0;
    if (GENERAL_RELATIVITY) {
      v_interface = gi(i01,i) / std::sqrt(SQR(gi(i01,i)) - gi(I00,i)*gi(i11,i));
    }

    // Determine region of wavefan
    Real *cons_interface, *flux_interface;
    if (lambda_l >= v_interface) {  // L region
      cons_interface = cons_l;
      flux_interface = flux_l;
    } else if (lambda_r <= v_interface) { // R region
      cons_interface = cons_r;
      flux_interface = flux_r;
    } else {  // HLL region
      cons_interface = cons_hll;
      flux_interface = flux_hll;
    }

    // Set conserved quantities in GR
    if (GENERAL_RELATIVITY) {
      for (int n = 0; n < NWAVE; ++n) {
        cons(n,i) = cons_interface[n];
      }
    }

    // Set fluxes
    for (int n = 0; n < NHYDRO; ++n) {
      flux(n,k,j,i) = flux_interface[n];
    }
    ey(k,j,i) = -flux_interface[IBY];
    ez(k,j,i) = flux_interface[IBZ];
  }

  // Transform fluxes to global coordinates if in GR
#if GENERAL_RELATIVITY
  {
    switch (ivx) {
      case IVX:
        pmb->pcoord->FluxToGlobal1(k, j, il, iu, cons, bb_normal, flux, ey, ez);
        break;
      case IVY:
        pmb->pcoord->FluxToGlobal2(k, j, il, iu, cons, bb_normal, flux, ey, ez);
        break;
      case IVZ:
        pmb->pcoord->FluxToGlobal3(k, j, il, iu, cons, bb_normal, flux, ey, ez);
        break;
    }
  }
#endif  // GENERAL_RELATIVITY
  return;
}

//----------------------------------------------------------------------------------------
// Non-frame-transforming HLLE implementation
// Inputs:
//   pmb: pointer to MeshBlock object
//   k,j: x3- and x2-indices
//   il,iu: lower and upper x1-indices
//   bb: 3D array of normal magnetic fields
//   g,gi: 1D scratch arrays for metric coefficients
//   prim_l,prim_r: 1D arrays of left and right primitive states
// Outputs:
//   flux: 3D array of hydrodynamical fluxes across interfaces
//   ey,ez: 3D arrays of magnetic fluxes (electric fields) across interfaces
// Notes:
//   implements HLLE algorithm similar to that of fluxcalc() in step_ch.c in Harm
//   derived from RiemannSolver() in hlle_mhd_rel_no_transform.cpp assuming ivx = IVY
//   same function as in hlld_rel.cpp
void HLLENonTransforming(MeshBlock *pmb, const int k, const int j,
                         const int il, const int iu, const AthenaArray<Real> &bb,
                         AthenaArray<Real> &g, AthenaArray<Real> &gi,
                         AthenaArray<Real> &prim_l, AthenaArray<Real> &prim_r,
                         AthenaArray<Real> &flux,
                         AthenaArray<Real> &ey, AthenaArray<Real> &ez) {
#if GENERAL_RELATIVITY
  // Extract ratio of specific heats
  const Real gamma_adi = pmb->peos->GetGamma();

  // Get metric components
  pmb->pcoord->Face2Metric(k, j, il, iu, g, gi);

  // Go through each interface
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    // Extract metric
    const Real &g_00 = g(I00,i), &g_01 = g(I01,i), &g_02 = g(I02,i), &g_03 = g(I03,i),
               &g_10 = g(I01,i), &g_11 = g(I11,i), &g_12 = g(I12,i), &g_13 = g(I13,i),
               &g_20 = g(I02,i), &g_21 = g(I12,i), &g_22 = g(I22,i), &g_23 = g(I23,i),
               &g_30 = g(I03,i), &g_31 = g(I13,i), &g_32 = g(I23,i), &g_33 = g(I33,i);
    const Real &g00 = gi(I00,i), &g01 = gi(I01,i), &g02 = gi(I02,i), &g03 = gi(I03,i),
               &g10 = gi(I01,i), &g11 = gi(I11,i), &g12 = gi(I12,i), &g13 = gi(I13,i),
               &g20 = gi(I02,i), &g21 = gi(I12,i), &g22 = gi(I22,i), &g23 = gi(I23,i),
               &g30 = gi(I03,i), &g31 = gi(I13,i), &g32 = gi(I23,i), &g33 = gi(I33,i);
    Real alpha = std::sqrt(-1.0/g00);

    // Extract left primitives
    const Real &rho_l = prim_l(IDN,i);
    const Real &pgas_l = prim_l(IPR,i);
    const Real &uu1_l = prim_l(IVX,i);
    const Real &uu2_l = prim_l(IVY,i);
    const Real &uu3_l = prim_l(IVZ,i);
    const Real &bb2_l = bb(k,j,i);
    const Real &bb3_l = prim_l(IBY,i);
    const Real &bb1_l = prim_l(IBZ,i);

    // Extract right primitives
    const Real &rho_r = prim_r(IDN,i);
    const Real &pgas_r = prim_r(IPR,i);
    const Real &uu1_r = prim_r(IVX,i);
    const Real &uu2_r = prim_r(IVY,i);
    const Real &uu3_r = prim_r(IVZ,i);
    const Real &bb2_r = bb(k,j,i);
    const Real &bb3_r = prim_r(IBY,i);
    const Real &bb1_r = prim_r(IBZ,i);

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
    pmb->peos->FastMagnetosonicSpeedsGR(wgas_l, pgas_l, ucon_l[0], ucon_l[IVY], b_sq_l,
                                        g00, g02, g22, &lambda_p_l, &lambda_m_l);

    // Calculate wavespeeds in right state
    Real lambda_p_r, lambda_m_r;
    Real wgas_r = rho_r + gamma_adi/(gamma_adi-1.0) * pgas_r;
    pmb->peos->FastMagnetosonicSpeedsGR(wgas_r, pgas_r, ucon_r[0], ucon_r[IVY], b_sq_r,
                                        g00, g02, g22, &lambda_p_r, &lambda_m_r);

    // Calculate extremal wavespeeds
    Real lambda_l = std::min(lambda_m_l, lambda_m_r);
    Real lambda_r = std::max(lambda_p_l, lambda_p_r);

    // Calculate conserved quantities in L region
    // (rho u^0, T^0_\mu, and B^j = *F^{j0}, where j != IVY)
    Real cons_l[NWAVE];
    Real wtot_l = wgas_l + b_sq_l;
    Real ptot_l = pgas_l + 0.5*b_sq_l;
    cons_l[IDN] = rho_l * ucon_l[0];
    cons_l[IEN] = wtot_l * ucon_l[0] * ucov_l[0] - bcon_l[0] * bcov_l[0] + ptot_l;
    cons_l[IVX] = wtot_l * ucon_l[0] * ucov_l[1] - bcon_l[0] * bcov_l[1];
    cons_l[IVY] = wtot_l * ucon_l[0] * ucov_l[2] - bcon_l[0] * bcov_l[2];
    cons_l[IVZ] = wtot_l * ucon_l[0] * ucov_l[3] - bcon_l[0] * bcov_l[3];
    cons_l[IBY] = bcon_l[IVZ] * ucon_l[0] - bcon_l[0] * ucon_l[IVZ];
    cons_l[IBZ] = bcon_l[IVX] * ucon_l[0] - bcon_l[0] * ucon_l[IVX];

    // Calculate fluxes in L region
    // (rho u^i, T^i_\mu, and *F^{ji}, where i = IVY and j != IVY)
    Real flux_l[NWAVE];
    flux_l[IDN] = rho_l * ucon_l[IVY];
    flux_l[IEN] = wtot_l * ucon_l[IVY] * ucov_l[0] - bcon_l[IVY] * bcov_l[0];
    flux_l[IVX] = wtot_l * ucon_l[IVY] * ucov_l[1] - bcon_l[IVY] * bcov_l[1];
    flux_l[IVY] = wtot_l * ucon_l[IVY] * ucov_l[2] - bcon_l[IVY] * bcov_l[2];
    flux_l[IVZ] = wtot_l * ucon_l[IVY] * ucov_l[3] - bcon_l[IVY] * bcov_l[3];
    flux_l[IVY] += ptot_l;
    flux_l[IBY] = bcon_l[IVZ] * ucon_l[IVY] - bcon_l[IVY] * ucon_l[IVZ];
    flux_l[IBZ] = bcon_l[IVX] * ucon_l[IVY] - bcon_l[IVY] * ucon_l[IVX];

    // Calculate conserved quantities in R region
    // (rho u^0, T^0_\mu, and B^j = *F^{j0}, where j != IVY)
    Real cons_r[NWAVE];
    Real wtot_r = wgas_r + b_sq_r;
    Real ptot_r = pgas_r + 0.5*b_sq_r;
    cons_r[IDN] = rho_r * ucon_r[0];
    cons_r[IEN] = wtot_r * ucon_r[0] * ucov_r[0] - bcon_r[0] * bcov_r[0] + ptot_r;
    cons_r[IVX] = wtot_r * ucon_r[0] * ucov_r[1] - bcon_r[0] * bcov_r[1];
    cons_r[IVY] = wtot_r * ucon_r[0] * ucov_r[2] - bcon_r[0] * bcov_r[2];
    cons_r[IVZ] = wtot_r * ucon_r[0] * ucov_r[3] - bcon_r[0] * bcov_r[3];
    cons_r[IBY] = bcon_r[IVZ] * ucon_r[0] - bcon_r[0] * ucon_r[IVZ];
    cons_r[IBZ] = bcon_r[IVX] * ucon_r[0] - bcon_r[0] * ucon_r[IVX];

    // Calculate fluxes in R region
    // (rho u^i, T^i_\mu, and *F^{ji}, where i = IVY and j != IVY)
    Real flux_r[NWAVE];
    flux_r[IDN] = rho_r * ucon_r[IVY];
    flux_r[IEN] = wtot_r * ucon_r[IVY] * ucov_r[0] - bcon_r[IVY] * bcov_r[0];
    flux_r[IVX] = wtot_r * ucon_r[IVY] * ucov_r[1] - bcon_r[IVY] * bcov_r[1];
    flux_r[IVY] = wtot_r * ucon_r[IVY] * ucov_r[2] - bcon_r[IVY] * bcov_r[2];
    flux_r[IVZ] = wtot_r * ucon_r[IVY] * ucov_r[3] - bcon_r[IVY] * bcov_r[3];
    flux_r[IVY] += ptot_r;
    flux_r[IBY] = bcon_r[IVZ] * ucon_r[IVY] - bcon_r[IVY] * ucon_r[IVZ];
    flux_r[IBZ] = bcon_r[IVX] * ucon_r[IVY] - bcon_r[IVY] * ucon_r[IVX];

    // Calculate fluxes in HLL region
    Real flux_hll[NWAVE];
    for (int n = 0; n < NWAVE; ++n) {
      flux_hll[n] = (lambda_r*flux_l[n] - lambda_l*flux_r[n]
                     + lambda_r*lambda_l * (cons_r[n] - cons_l[n])) / (lambda_r-lambda_l);
    }

    // Determine region of wavefan
    Real *flux_interface;
    if (lambda_l >= 0.0) {  // L region
      flux_interface = flux_l;
    } else if (lambda_r <= 0.0) { // R region
      flux_interface = flux_r;
    } else {  // HLL region
      flux_interface = flux_hll;
    }

    // Set fluxes
    for (int n = 0; n < NHYDRO; ++n) {
      flux(n,k,j,i) = flux_interface[n];
    }
    ey(k,j,i) = -flux_interface[IBY];
    ez(k,j,i) = flux_interface[IBZ];
  }
  return;
#endif  // GENERAL_RELATIVITY
}
} // namespace
