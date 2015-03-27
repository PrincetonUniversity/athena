// HLLE Riemann solver for relativistic magnetohydrodynamics

// Primary header
#include "../../fluid_integrator.hpp"

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()

// Athena headers
#include "../../../eos/eos.hpp"                     // GetGamma(),
                                                    //     FastMagnetosonicSpeedsSR()
#include "../../../fluid.hpp"                       // Fluid
#include "../../../../athena.hpp"                   // enums, macros, Real
#include "../../../../athena_arrays.hpp"            // AthenaArray
#include "../../../../coordinates/coordinates.hpp"  // Coordinates
#include "../../../../mesh.hpp"                     // MeshBlock

// Declarations
static void PrimToFluxFlat(Real gamma_adi_red, Real rho, Real pgas,
    const Real u_con[4], const Real b_con[4],
    Real flux[NWAVE], int ivx, int ivy, int ivz);
static void PrimToConsFlat(Real gamma_adi_red, Real rho, Real pgas,
    const Real u_con[4], const Real b_con[4],
    Real cons[NWAVE], int ivx, int ivy, int ivz);

// Riemann solver
// Inputs:
//   il,iu: lower and upper indices for interfaces
//   prim_left, prim_right: left and right primitive states
// Outputs:
//   flux: fluxes
// Notes:
//   prim_left, prim_right overwritten
//   implements HLLE algorithm from Mignone & Bodo 2005, MNRAS 364 126 (MB2005)
//   references Mignone & Bodo 2006, MNRAS 368 1040 (MB2006)
void FluidIntegrator::RiemannSolver(const int k, const int j, const int il,
    const int iu, const int ivx, const AthenaArray<Real> &b,
    AthenaArray<Real> &prim_left, AthenaArray<Real> &prim_right,
    AthenaArray<Real> &flux)
{
  // Transform primitives to locally flat coordinates if in GR
  if (GENERAL_RELATIVITY)
    switch (ivx)
    {
      case IVX:
        pmy_fluid->pmy_block->pcoord->PrimToLocal1(k, j, il, iu, b,
            prim_left, prim_right, b_normal_);
        break;
      case IVY:
        pmy_fluid->pmy_block->pcoord->PrimToLocal2(k, j, il, iu, b,
            prim_left, prim_right, b_normal_);
        break;
      case IVZ:
        pmy_fluid->pmy_block->pcoord->PrimToLocal3(k, j, il, iu, b,
            prim_left, prim_right, b_normal_);
        break;
    }
  else  // SR; need to populate 1D normal B array
  {
    #pragma simd
    for (int i = il; i <= iu; i++)
      b_normal_(i) = b(k,j,i);
  }

  // Calculate cyclic permutations of indices
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  // Extract ratio of specific heats
  const Real gamma_adi = pmy_fluid->pf_eos->GetGamma();
  const Real gamma_adi_red = gamma_adi/(gamma_adi-1.0);

  // Go through each interface
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract left primitives
    const Real &rho_left = prim_left(IDN,i);
    const Real &pgas_left = prim_left(IEN,i);
    const Real &vx_left = prim_left(ivx,i);
    const Real &vy_left = prim_left(ivy,i);
    const Real &vz_left = prim_left(ivz,i);
    const Real &by_left = prim_left(IBY,i);
    const Real &bz_left = prim_left(IBZ,i);

    // Extract right primitives
    const Real &rho_right = prim_right(IDN,i);
    const Real &pgas_right = prim_right(IEN,i);
    const Real &vx_right = prim_right(ivx,i);
    const Real &vy_right = prim_right(ivy,i);
    const Real &vz_right = prim_right(ivz,i);
    const Real &by_right = prim_right(IBY,i);
    const Real &bz_right = prim_right(IBZ,i);

    // Extract normal magnetic field
    const Real &bx = b_normal_(i);

    // Calculate 4-vectors for left primitives
    Real u_con_left[4], b_con_left[4];
    u_con_left[0] = std::sqrt(1.0/(1.0-(SQR(vx_left)+SQR(vy_left)+SQR(vz_left))));
    u_con_left[1] = u_con_left[0] * vx_left;
    u_con_left[2] = u_con_left[0] * vy_left;
    u_con_left[3] = u_con_left[0] * vz_left;
    b_con_left[0] = bx*u_con_left[1] + by_left*u_con_left[2] + bz_left*u_con_left[3];
    b_con_left[1] = (bx + b_con_left[0] * u_con_left[1]) / u_con_left[0];
    b_con_left[2] = (by_left + b_con_left[0] * u_con_left[2]) / u_con_left[0];
    b_con_left[3] = (bz_left + b_con_left[0] * u_con_left[3]) / u_con_left[0];

    // Calculate 4-vectors for right primitives
    Real u_con_right[4], b_con_right[4];
    u_con_right[0] = std::sqrt(1.0/(1.0-(SQR(vx_right)+SQR(vy_right)+SQR(vz_right))));
    u_con_right[1] = u_con_right[0] * vx_right;
    u_con_right[2] = u_con_right[0] * vy_right;
    u_con_right[3] = u_con_right[0] * vz_right;
    b_con_right[0] = bx*u_con_right[1] + by_right*u_con_right[2]
        + bz_right*u_con_right[3];
    b_con_right[1] = (bx + b_con_right[0] * u_con_right[1]) / u_con_right[0];
    b_con_right[2] = (by_right + b_con_right[0] * u_con_right[2]) / u_con_right[0];
    b_con_right[3] = (bz_right + b_con_right[0] * u_con_right[3]) / u_con_right[0];

    // Calculate wavespeeds
    Real lambda_left_plus, lambda_left_minus;
    pmy_fluid->pf_eos->FastMagnetosonicSpeedsSR(
        rho_left, pgas_left, u_con_left, b_con_left,
        &lambda_left_plus, &lambda_left_minus);                          // (MB2006 56)
    Real lambda_right_plus, lambda_right_minus;
    pmy_fluid->pf_eos->FastMagnetosonicSpeedsSR(
        rho_right, pgas_right, u_con_right, b_con_right,
        &lambda_right_plus, &lambda_right_minus);                        // (MB2006 56)
    Real lambda_left = std::min(lambda_left_minus, lambda_right_minus);  // (MB2006 55)
    Real lambda_right = std::max(lambda_left_plus, lambda_right_plus);   // (MB2006 55)

    // Calculate L/R state fluxes
    Real flux_left[NWAVE], flux_right[NWAVE];
    PrimToFluxFlat(gamma_adi_red, rho_left, pgas_left,
        u_con_left, b_con_left,
        flux_left, ivx, ivy, ivz);
    PrimToFluxFlat(gamma_adi_red, rho_right, pgas_right,
        u_con_right, b_con_right,
        flux_right, ivx, ivy, ivz);

    // Set fluxes if in L state
    if (lambda_left >= 0.0)
    {
      for (int n = 0; n < NWAVE; ++n)
        flux(n,i) = flux_left[n];
      continue;
    }

    // Set fluxes if in R state
    if (lambda_right <= 0.0)
    {
      for (int n = 0; n < NWAVE; ++n)
        flux(n,i) = flux_right[n];
      continue;
    }

    // Set fluxes in HLL state
    Real cons_left[NWAVE], cons_right[NWAVE];
    PrimToConsFlat(gamma_adi_red, rho_left, pgas_left,
        u_con_left, b_con_left,
        cons_left, ivx, ivy, ivz);
    PrimToConsFlat(gamma_adi_red, rho_right, pgas_right,
        u_con_right, b_con_right,
        cons_right, ivx, ivy, ivz);
    for (int n = 0; n < NWAVE; ++n)
      flux(n,i) = (lambda_right*flux_left[n] - lambda_left*flux_right[n]
          + lambda_right*lambda_left * (cons_right[n] - cons_left[n]))
          / (lambda_right-lambda_left);                                   // (MB2005 11)
  }

  // Transform fluxes to global coordinates if in GR
  if (GENERAL_RELATIVITY)
    switch (ivx)
    {
      case IVX:
        pmy_fluid->pmy_block->pcoord->FluxToGlobal1(k, j, il, iu, flux);
        break;
      case IVY:
        pmy_fluid->pmy_block->pcoord->FluxToGlobal2(k, j, il, iu, flux);
        break;
      case IVZ:
        pmy_fluid->pmy_block->pcoord->FluxToGlobal3(k, j, il, iu, flux);
        break;
    }
  return;
}

// Function for converting constant primitive state to flux state in flat spacetime
// Notes:
//   same function as in hlld_mhd_rel.cpp
//   implements (15) from Mignone, Ugliano, & Bodo 2009, MNRAS 393 1141
//     note B^i v^x - B^x v^i = b^i u^x - b^x u^i
static void PrimToFluxFlat(Real gamma_adi_red, Real rho, Real pgas,
    const Real u_con[4], const Real b_con[4],
    Real flux[NWAVE], int ivx, int ivy, int ivz)
{
  Real b_sq = -SQR(b_con[0]) + SQR(b_con[1]) + SQR(b_con[2]) + SQR(b_con[3]);
  Real w = rho + gamma_adi_red * pgas + b_sq;
  Real ptot = pgas + 0.5*b_sq;
  flux[IDN] = rho*u_con[1];
  flux[IEN] = w*u_con[0]*u_con[1] - b_con[0]*b_con[1];
  flux[ivx] = w*u_con[1]*u_con[1] - b_con[1]*b_con[1] + ptot;
  flux[ivy] = w*u_con[2]*u_con[1] - b_con[2]*b_con[1];
  flux[ivz] = w*u_con[3]*u_con[1] - b_con[3]*b_con[1];
  flux[IBY] = b_con[2]*u_con[1] - b_con[1]*u_con[2];
  flux[IBZ] = b_con[3]*u_con[1] - b_con[1]*u_con[3];
  return;
}

// Function for converting primitive state to conserved state in flat spacetime
// Notes:
//   same function as in hlld_rel.cpp
//   references Mignone, Ugliano, & Bodo 2009, MNRAS 393 1141 (MUB)
static void PrimToConsFlat(Real gamma_adi_red, Real rho, Real pgas,
    const Real u_con[4], const Real b_con[4],
    Real cons[NWAVE], int ivx, int ivy, int ivz)
{
  Real b_sq = -SQR(b_con[0]) + SQR(b_con[1]) + SQR(b_con[2]) + SQR(b_con[3]);
  Real w = rho + gamma_adi_red * pgas + b_sq;
  Real ptot = pgas + 0.5*b_sq;
  cons[IDN] = rho*u_con[0];
  cons[IEN] = w*u_con[0]*u_con[0] - b_con[0]*b_con[0] - ptot;  // (MUB 8)
  cons[ivx] = w*u_con[1]*u_con[0] - b_con[1]*b_con[0];         // (MUB 8)
  cons[ivy] = w*u_con[2]*u_con[0] - b_con[2]*b_con[0];         // (MUB 8)
  cons[ivz] = w*u_con[3]*u_con[0] - b_con[3]*b_con[0];         // (MUB 8)
  cons[IBY] = b_con[2]*u_con[0] - b_con[0]*u_con[2];
  cons[IBZ] = b_con[3]*u_con[0] - b_con[0]*u_con[3];
  return;
}
