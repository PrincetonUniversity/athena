// Local Lax-Friedrichs Riemann solver for relativistic magnetohydrodynamics in pure GR

// Primary header
#include "../../fluid_integrator.hpp"

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()

// Athena headers
#include "../../../eos/eos.hpp"                     // GetGamma(),
                                                    //     FastMagnetosonicSpeedsGR()
#include "../../../fluid.hpp"                       // Fluid
#include "../../../../athena.hpp"                   // enums, macros, Real
#include "../../../../athena_arrays.hpp"            // AthenaArray
#include "../../../../coordinates/coordinates.hpp"  // Coordinates
#include "../../../../mesh.hpp"                     // MeshBlock

// Declarations
static void PrimToFluxGR(Real gamma_adi_red, Real rho, Real pgas, Real b_sq,
    const Real u_con[4], const Real u_cov[4], const Real b_con[4], const Real b_cov[4],
    Real flux[NWAVE], int ivx, int ivy, int ivz);
static void PrimToConsGR(Real gamma_adi_red, Real rho, Real pgas, Real b_sq,
    const Real u_con[4], const Real u_cov[4], const Real b_con[4], const Real b_cov[4],
    Real cons[NWAVE], int ivy, int ivz);

// Riemann solver
// Inputs:
//   il,iu: lower and upper indices for interfaces
//   prim_left, prim_right: left and right primitive states
// Outputs:
//   flux: fluxes
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
    Real b1_left, b2_left, b3_left;
    switch (ivx)
    {
      case IVX:
        b1_left = b(k,j,i);
        b2_left = prim_left(IBY,i);
        b3_left = prim_left(IBZ,i);
        break;
      case IVY:
        b2_left = b(k,j,i);
        b3_left = prim_left(IBY,i);
        b1_left = prim_left(IBZ,i);
        break;
      case IVZ:
        b3_left = b(k,j,i);
        b1_left = prim_left(IBY,i);
        b2_left = prim_left(IBZ,i);
        break;
    }

    // Extract right primitives
    const Real &rho_right = prim_right(IDN,i);
    const Real &pgas_right = prim_right(IEN,i);
    const Real &v1_right = prim_right(IVX,i);
    const Real &v2_right = prim_right(IVY,i);
    const Real &v3_right = prim_right(IVZ,i);
    Real b1_right, b2_right, b3_right;
    switch (ivx)
    {
      case IVX:
        b1_right = b(k,j,i);
        b2_right = prim_right(IBY,i);
        b3_right = prim_right(IBZ,i);
        break;
      case IVY:
        b2_right = b(k,j,i);
        b3_right = prim_right(IBY,i);
        b1_right = prim_right(IBZ,i);
        break;
      case IVZ:
        b3_right = b(k,j,i);
        b1_right = prim_right(IBY,i);
        b2_right = prim_right(IBZ,i);
        break;
    }

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
        + g_13*u_con_left[3];
    u_cov_left[2] = g_20*u_con_left[0] + g_21*u_con_left[1] + g_22*u_con_left[2]
        + g_23*u_con_left[3];
    u_cov_left[3] = g_30*u_con_left[0] + g_31*u_con_left[1] + g_32*u_con_left[2]
        + g_33*u_con_left[3];

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
        + g_13*u_con_right[3];
    u_cov_right[2] = g_20*u_con_right[0] + g_21*u_con_right[1] + g_22*u_con_right[2]
        + g_23*u_con_right[3];
    u_cov_right[3] = g_30*u_con_right[0] + g_31*u_con_right[1] + g_32*u_con_right[2]
        + g_33*u_con_right[3];

    // Calculate 4-magnetic field for left primitives
    Real b_con_left[4], b_cov_left[4];
    b_con_left[0] = u_con_left[0] * (g_01*b1_left + g_02*b2_left + g_03*b3_left)
                  + u_con_left[1] * (g_11*b1_left + g_12*b2_left + g_13*b3_left)
                  + u_con_left[2] * (g_21*b1_left + g_22*b2_left + g_23*b3_left)
                  + u_con_left[3] * (g_31*b1_left + g_32*b2_left + g_33*b3_left);
    b_con_left[1] = (b1_left + b_con_left[0] * u_con_left[1]) / u_con_left[0];
    b_con_left[2] = (b2_left + b_con_left[0] * u_con_left[2]) / u_con_left[0];
    b_con_left[3] = (b3_left + b_con_left[0] * u_con_left[3]) / u_con_left[0];
    b_cov_left[0] = g_00*b_con_left[0] + g_01*b_con_left[1] + g_02*b_con_left[2]
        + g_03*b_con_left[3];
    b_cov_left[1] = g_10*b_con_left[0] + g_11*b_con_left[1] + g_12*b_con_left[2]
        + g_13*b_con_left[3];
    b_cov_left[2] = g_20*b_con_left[0] + g_21*b_con_left[1] + g_22*b_con_left[2]
        + g_23*b_con_left[3];
    b_cov_left[3] = g_30*b_con_left[0] + g_31*b_con_left[1] + g_32*b_con_left[2]
        + g_33*b_con_left[3];
    Real b_sq_left = b_con_left[0]*b_cov_left[0] + b_con_left[1]*b_cov_left[1]
        + b_con_left[2]*b_cov_left[2] + b_con_left[3]*b_cov_left[3];

    // Calculate 4-magnetic field for right primitives
    Real b_con_right[4], b_cov_right[4];
    b_con_right[0] = u_con_right[0] * (g_01*b1_right + g_02*b2_right + g_03*b3_right)
                  + u_con_right[1] * (g_11*b1_right + g_12*b2_right + g_13*b3_right)
                  + u_con_right[2] * (g_21*b1_right + g_22*b2_right + g_23*b3_right)
                  + u_con_right[3] * (g_31*b1_right + g_32*b2_right + g_33*b3_right);
    b_con_right[1] = (b1_right + b_con_right[0] * u_con_right[1]) / u_con_right[0];
    b_con_right[2] = (b2_right + b_con_right[0] * u_con_right[2]) / u_con_right[0];
    b_con_right[3] = (b3_right + b_con_right[0] * u_con_right[3]) / u_con_right[0];
    b_cov_right[0] = g_00*b_con_right[0] + g_01*b_con_right[1] + g_02*b_con_right[2]
        + g_03*b_con_right[3];
    b_cov_right[1] = g_10*b_con_right[0] + g_11*b_con_right[1] + g_12*b_con_right[2]
        + g_13*b_con_right[3];
    b_cov_right[2] = g_20*b_con_right[0] + g_21*b_con_right[1] + g_22*b_con_right[2]
        + g_23*b_con_right[3];
    b_cov_right[3] = g_30*b_con_right[0] + g_31*b_con_right[1] + g_32*b_con_right[2]
        + g_33*b_con_right[3];
    Real b_sq_right = b_con_right[0]*b_cov_right[0] + b_con_right[1]*b_cov_right[1]
        + b_con_right[2]*b_cov_right[2] + b_con_right[3]*b_cov_right[3];

    // Calculate wavespeeds
    Real lambda_left_plus, lambda_left_minus;
    Real rho_h_left = rho_left + gamma_adi_red * pgas_left;
    pmy_fluid->pf_eos->FastMagnetosonicSpeedsGR(
        rho_h_left, pgas_left, u_con_left[0], u_con_left[ivx], b_sq_left,
        g00, g_inv_off_diag, g_inv_diag,
        &lambda_left_plus, &lambda_left_minus);
    Real lambda_right_plus, lambda_right_minus;
    Real rho_h_right = rho_right + gamma_adi_red * pgas_right;
    pmy_fluid->pf_eos->FastMagnetosonicSpeedsGR(
        rho_h_right, pgas_right, u_con_right[0], u_con_right[ivx], b_sq_right,
        g00, g_inv_off_diag, g_inv_diag,
        &lambda_right_plus, &lambda_right_minus);
    Real lambda_left = std::min(lambda_left_minus, lambda_right_minus);
    Real lambda_right = std::max(lambda_left_plus, lambda_right_plus);
    Real lambda = std::max(lambda_right, -lambda_left);

    // Calculate L/R state conserved quantities and fluxes
    Real cons_left[NWAVE], cons_right[NWAVE];
    PrimToConsGR(gamma_adi_red, rho_left, pgas_left, b_sq_left,
        u_con_left, u_cov_left, b_con_left, b_cov_left,
        cons_left, ivy, ivz);
    PrimToConsGR(gamma_adi_red, rho_right, pgas_right, b_sq_right,
        u_con_right, u_cov_right, b_con_right, b_cov_right,
        cons_right, ivy, ivz);
    Real flux_left[NWAVE], flux_right[NWAVE];
    PrimToFluxGR(gamma_adi_red, rho_left, pgas_left, b_sq_left,
        u_con_left, u_cov_left, b_con_left, b_cov_left,
        flux_left, ivx, ivy, ivz);
    PrimToFluxGR(gamma_adi_red, rho_right, pgas_right, b_sq_right,
        u_con_right, u_cov_right, b_con_right, b_cov_right,
        flux_right, ivx, ivy, ivz);

    // Set fluxes
    for (int n = 0; n < NWAVE; ++n)
      flux(n,i) = 0.5 * (flux_left[n] + flux_right[n]
          - lambda * (cons_right[n] - cons_left[n]));
  }
  return;
}

// Function for converting primitive state to flux state
// Notes:
//   calculates rho u^i and T^i_\mu, where i = ivx
static void PrimToFluxGR(Real gamma_adi_red, Real rho, Real pgas, Real b_sq,
    const Real u_con[4], const Real u_cov[4], const Real b_con[4], const Real b_cov[4],
    Real flux[NWAVE], int ivx, int ivy, int ivz)
{
  Real w = rho + gamma_adi_red * pgas + b_sq;
  Real ptot = pgas + 0.5 * b_sq;
  flux[IDN] = rho*u_con[ivx];
  flux[IEN] = w*u_con[ivx]*u_cov[0] - b_con[ivx]*b_cov[0];
  flux[IVX] = w*u_con[ivx]*u_cov[1] - b_con[ivx]*b_cov[1];
  flux[IVY] = w*u_con[ivx]*u_cov[2] - b_con[ivx]*b_cov[2];
  flux[IVZ] = w*u_con[ivx]*u_cov[3] - b_con[ivx]*b_cov[3];
  flux[IBY] = b_con[ivy]*u_con[ivx] - b_con[ivx]*u_con[ivy];
  flux[IBZ] = b_con[ivz]*u_con[ivx] - b_con[ivx]*u_con[ivz];
  flux[ivx] += ptot;
  return;
}

// Function for converting primitive state to conserved state
// Notes:
//   calculates rho u^0 and T^0_\mu
static void PrimToConsGR(Real gamma_adi_red, Real rho, Real pgas, Real b_sq,
    const Real u_con[4], const Real u_cov[4], const Real b_con[4], const Real b_cov[4],
    Real cons[NWAVE], int ivy, int ivz)
{
  Real w = rho + gamma_adi_red * pgas + b_sq;
  Real ptot = pgas + 0.5 * b_sq;
  cons[IDN] = rho*u_con[0];
  cons[IEN] = w*u_con[0]*u_cov[0] - b_con[0]*b_cov[0] + ptot;
  cons[IVX] = w*u_con[0]*u_cov[1] - b_con[0]*b_cov[1];
  cons[IVY] = w*u_con[0]*u_cov[2] - b_con[0]*b_cov[2];
  cons[IVZ] = w*u_con[0]*u_cov[3] - b_con[0]*b_cov[3];
  cons[IBY] = b_con[ivy]*u_con[0] - b_con[0]*u_con[ivy];
  cons[IBZ] = b_con[ivz]*u_con[0] - b_con[0]*u_con[ivz];
  return;
}
