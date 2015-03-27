// Local Lax-Friedrichs Riemann solver for relativistic hydrodynamics in pure GR

// Primary header
#include "../../fluid_integrator.hpp"

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()

// Athena headers
#include "../../../eos/eos.hpp"                     // GetGamma()
#include "../../../fluid.hpp"                       // Fluid
#include "../../../../athena.hpp"                   // enums, macros, Real
#include "../../../../athena_arrays.hpp"            // AthenaArray
#include "../../../../coordinates/coordinates.hpp"  // Coordinates
#include "../../../../mesh.hpp"                     // MeshBlock

// Declarations
static void PrimToFluxGR(Real gamma_adi_red, Real rho, Real pgas,
    const Real u_con[4], const Real u_cov[4],
    Real flux[NWAVE], int ivx);
static void PrimToConsGR(Real gamma_adi_red, Real rho, Real pgas,
    const Real u_con[4], const Real u_cov[4],
    Real cons[NWAVE], int ivx);

// Riemann solver
// Inputs:
//   k,j: x3- and x2-indices
//   il,iu: lower and upper x1-indices
//   ivx: type of interface (IVX for x1, IVY for x2, IVZ for x3)
//   b: 3D array of normal magnetic fields (not used)
//   prim_left, prim_right: left and right primitive states
// Outputs:
//   flux: fluxes across interface
// Notes:
//   prim_left, prim_right overwritten
//   implements LLF algorithm similar to that of fluxcalc() in step_ch.c in Harm
void FluidIntegrator::RiemannSolver(const int k, const int j, const int il,
    const int iu, const int ivx, const AthenaArray<Real> &b,
    AthenaArray<Real> &prim_left, AthenaArray<Real> &prim_right,
    AthenaArray<Real> &flux)
{
  // Calculate cyclic permutations of indices
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  // Extract ratio of specific heats
  const Real gamma_adi = pmy_fluid->pf_eos->GetGamma();
  const Real gamma_adi_red = gamma_adi/(gamma_adi-1.0);

  // Get metric components
  switch (ivx)
  {
    case IVX:
      pmy_fluid->pmy_block->pcoord->Face1Metric(k, j, il, iu, g_, g_inv_);
      break;
    case IVY:
      pmy_fluid->pmy_block->pcoord->Face2Metric(k, j, il, iu, g_, g_inv_);
      break;
    case IVZ:
      pmy_fluid->pmy_block->pcoord->Face3Metric(k, j, il, iu, g_, g_inv_);
      break;
  }

  // Go through each interface
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract metric
    const Real &g_00 = g_(I00,i), &g_01 = g_(I01,i), &g_02 = g_(I02,i),
          &g_03 = g_(I03,i);
    const Real &g_10 = g_(I01,i), &g_11 = g_(I11,i), &g_12 = g_(I12,i),
          &g_13 = g_(I13,i);
    const Real &g_20 = g_(I02,i), &g_21 = g_(I12,i), &g_22 = g_(I22,i),
          &g_23 = g_(I23,i);
    const Real &g_30 = g_(I03,i), &g_31 = g_(I13,i), &g_32 = g_(I23,i),
          &g_33 = g_(I33,i);

    // Extract inverse of metric
    const Real &g00 = g_inv_(I00,i), &g01 = g_inv_(I01,i), &g02 = g_inv_(I02,i),
         &g03 = g_inv_(I03,i);
    const Real &g10 = g_inv_(I01,i), &g11 = g_inv_(I11,i), &g12 = g_inv_(I12,i),
         &g13 = g_inv_(I13,i);
    const Real &g20 = g_inv_(I02,i), &g21 = g_inv_(I12,i), &g22 = g_inv_(I22,i),
         &g23 = g_inv_(I23,i);
    const Real &g30 = g_inv_(I03,i), &g31 = g_inv_(I13,i), &g32 = g_inv_(I23,i),
         &g33 = g_inv_(I33,i);
    Real g_inv_diag, g_inv_off_diag;
    switch (ivx)
    {
      case IVX:
        g_inv_diag = g11;
        g_inv_off_diag = g01;
        break;
      case IVY:
        g_inv_diag = g22;
        g_inv_off_diag = g02;
        break;
      case IVZ:
        g_inv_diag = g33;
        g_inv_off_diag = g03;
        break;
    }

    // Extract left primitives
    const Real &rho_left = prim_left(IDN,i);
    const Real &pgas_left = prim_left(IEN,i);
    const Real &v1_left = prim_left(IVX,i);
    const Real &v2_left = prim_left(IVY,i);
    const Real &v3_left = prim_left(IVZ,i);

    // Extract right primitives
    const Real &rho_right = prim_right(IDN,i);
    const Real &pgas_right = prim_right(IEN,i);
    const Real &v1_right = prim_right(IVX,i);
    const Real &v2_right = prim_right(IVY,i);
    const Real &v3_right = prim_right(IVZ,i);

    // Calculate 4-velocity for left primitives
    Real u_con_left[4], u_cov_left[4];
    Real neg_inv_u0_sq = g_00 + 2.0*g_01*v1_left + 2.0*g_02*v2_left + 2.0*g_03*v3_left
                       + g_11*SQR(v1_left) + 2.0*g_12*v1_left*v2_left
                           + 2.0*g_13*v1_left*v3_left
                       + g_22*SQR(v2_left) + 2.0*g_23*v2_left*v3_left
                       + g_33*SQR(v3_left);
    u_con_left[0] = std::sqrt(-1.0/neg_inv_u0_sq);
    u_con_left[1] = u_con_left[0] * v1_left;
    u_con_left[2] = u_con_left[0] * v2_left;
    u_con_left[3] = u_con_left[0] * v3_left;
    u_cov_left[0] = g_00*u_con_left[0] + g_01*u_con_left[1] + g_02*u_con_left[2]
        + g_03*u_con_left[3];
    u_cov_left[1] = g_10*u_con_left[0] + g_11*u_con_left[1] + g_12*u_con_left[2]
        + g_03*u_con_left[3];
    u_cov_left[2] = g_20*u_con_left[0] + g_21*u_con_left[1] + g_22*u_con_left[2]
        + g_03*u_con_left[3];
    u_cov_left[3] = g_30*u_con_left[0] + g_31*u_con_left[1] + g_32*u_con_left[2]
        + g_03*u_con_left[3];

    // Calculate 4-velocity for right primitives
    Real u_con_right[4], u_cov_right[4];
    neg_inv_u0_sq = g_00 + 2.0*g_01*v1_right + 2.0*g_02*v2_right + 2.0*g_03*v3_right
                  + g_11*SQR(v1_right) + 2.0*g_12*v1_right*v2_right
                      + 2.0*g_13*v1_right*v3_right
                  + g_22*SQR(v2_right) + 2.0*g_23*v2_right*v3_right
                  + g_33*SQR(v3_right);
    u_con_right[0] = std::sqrt(-1.0/neg_inv_u0_sq);
    u_con_right[1] = u_con_right[0] * v1_right;
    u_con_right[2] = u_con_right[0] * v2_right;
    u_con_right[3] = u_con_right[0] * v3_right;
    u_cov_right[0] = g_00*u_con_right[0] + g_01*u_con_right[1] + g_02*u_con_right[2]
        + g_03*u_con_right[3];
    u_cov_right[1] = g_10*u_con_right[0] + g_11*u_con_right[1] + g_12*u_con_right[2]
        + g_03*u_con_right[3];
    u_cov_right[2] = g_20*u_con_right[0] + g_21*u_con_right[1] + g_22*u_con_right[2]
        + g_03*u_con_right[3];
    u_cov_right[3] = g_30*u_con_right[0] + g_31*u_con_right[1] + g_32*u_con_right[2]
        + g_03*u_con_right[3];

    // Calculate wavespeeds in left region
    Real lambda_plus_left, lambda_minus_left;
    Real w_left = rho_left + gamma_adi_red * pgas_left;
    pmy_fluid->pf_eos->SoundSpeedsGR(w_left, pgas_left, u_con_left[0], u_con_left[ivx],
        g00, g_inv_off_diag, g_inv_diag,
        &lambda_plus_left, &lambda_minus_left);

    // Calculate wavespeeds in right region
    Real lambda_plus_right, lambda_minus_right;
    Real w_right = rho_right + gamma_adi_red * pgas_right;
    pmy_fluid->pf_eos->SoundSpeedsGR(
        w_right, pgas_right, u_con_right[0], u_con_right[ivx],
        g00, g_inv_off_diag, g_inv_diag,
        &lambda_plus_right, &lambda_minus_right);

    // Calculate extremal wavespeed
    Real lambda_left = std::min(lambda_minus_left, lambda_minus_right);
    Real lambda_right = std::max(lambda_plus_left, lambda_plus_right);
    Real lambda = std::max(lambda_right, -lambda_left);

    // Calculate L/R state conserved quantities and fluxes
    Real cons_left[NWAVE], cons_right[NWAVE];
    PrimToConsGR(gamma_adi_red, rho_left, pgas_left,
        u_con_left, u_cov_left,
        cons_left, ivx);
    PrimToConsGR(gamma_adi_red, rho_right, pgas_right,
        u_con_right, u_cov_right,
        cons_right, ivx);
    Real flux_left[NWAVE], flux_right[NWAVE];
    PrimToFluxGR(gamma_adi_red, rho_left, pgas_left,
        u_con_left, u_cov_left,
        flux_left, ivx);
    PrimToFluxGR(gamma_adi_red, rho_right, pgas_right,
        u_con_right, u_cov_right,
        flux_right, ivx);

    // Set fluxes in HLL state
    for (int n = 0; n < NWAVE; ++n)
      flux(n,i) = 0.5 * (flux_left[n] + flux_right[n]
          - lambda * (cons_right[n] - cons_left[n]));
  }
  return;
}

// Function for converting constant primitive state to flux state
// Notes:
//   calculates rho u^i and T^i_\mu, where i = ivx
static void PrimToFluxGR(Real gamma_adi_red, Real rho, Real pgas,
    const Real u_con[4], const Real u_cov[4],
    Real flux[NWAVE], int ivx)
{
  Real w = rho + gamma_adi_red * pgas;
  flux[IDN] = rho*u_con[ivx];
  flux[IEN] = w*u_con[ivx]*u_cov[0];
  flux[IVX] = w*u_con[ivx]*u_cov[1];
  flux[IVY] = w*u_con[ivx]*u_cov[2];
  flux[IVZ] = w*u_con[ivx]*u_cov[3];
  flux[ivx] += pgas;
  return;
}

// Function for converting primitive state to conserved state
// Notes:
//   calculates rho u^0 and T^0_\mu
static void PrimToConsGR(Real gamma_adi_red, Real rho, Real pgas,
    const Real u_con[4], const Real u_cov[4],
    Real cons[NWAVE], int ivx)
{
  Real w = rho + gamma_adi_red * pgas;
  cons[IDN] = rho*u_con[0];
  cons[IEN] = w*u_con[0]*u_cov[0] + pgas;
  cons[IVX] = w*u_con[0]*u_cov[1];
  cons[IVY] = w*u_con[0]*u_cov[2];
  cons[IVZ] = w*u_con[0]*u_cov[3];
  return;
}
