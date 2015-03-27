// HLLE Riemann solver for relativistic magnetohydrodynamics

// Primary header
#include "../../fluid_integrator.hpp"

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()

// Athena headers
#include "../../../eos/eos.hpp"                     // GetGamma(), SoundSpeedsSR()
#include "../../../fluid.hpp"                       // Fluid
#include "../../../../athena.hpp"                   // enums, macros, Real
#include "../../../../athena_arrays.hpp"            // AthenaArray
#include "../../../../coordinates/coordinates.hpp"  // Coordinates
#include "../../../../mesh.hpp"                     // MeshBlock

// Declarations
static void PrimToFluxFlat(Real gamma_adi_red, Real rho, Real pgas, const Real u_con[4],
    Real flux[NWAVE], int ivx, int ivy, int ivz);
static void PrimToConsFlat(Real gamma_adi_red, Real rho, Real pgas, const Real u_con[4],
    Real cons[NWAVE], int ivx, int ivy, int ivz);

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
//   implements HLLE algorithm from Mignone & Bodo 2005, MNRAS 364 126 (MB)
//   prim_left, prim_right overwritten
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

    // Extract right primitives
    const Real &rho_right = prim_right(IDN,i);
    const Real &pgas_right = prim_right(IEN,i);
    const Real &vx_right = prim_right(ivx,i);
    const Real &vy_right = prim_right(ivy,i);
    const Real &vz_right = prim_right(ivz,i);

    // Calculate 4-velocity for left primitives
    Real u_con_left[4];
    u_con_left[0] = std::sqrt(1.0/(1.0-(SQR(vx_left)+SQR(vy_left)+SQR(vz_left))));
    u_con_left[1] = u_con_left[0] * vx_left;
    u_con_left[2] = u_con_left[0] * vy_left;
    u_con_left[3] = u_con_left[0] * vz_left;

    // Calculate 4-velocity for right primitives
    Real u_con_right[4];
    u_con_right[0] = std::sqrt(1.0/(1.0-(SQR(vx_right)+SQR(vy_right)+SQR(vz_right))));
    u_con_right[1] = u_con_right[0] * vx_right;
    u_con_right[2] = u_con_right[0] * vy_right;
    u_con_right[3] = u_con_right[0] * vz_right;

    // Calculate wavespeeds in left region
    Real lambda_plus_left, lambda_minus_left;
    Real rho_h_left = rho_left + gamma_adi_red * pgas_left;
    Real v_sq_left = SQR(vx_left) + SQR(vy_left) + SQR(vz_left);
    Real gamma_sq_left = 1.0/(1.0-v_sq_left);
    pmy_fluid->pf_eos->SoundSpeedsSR(
        rho_h_left, pgas_left, vx_left, gamma_sq_left,
        &lambda_plus_left, &lambda_minus_left);                   // (MB 23)

    // Calculate wavespeeds in right region
    Real lambda_plus_right, lambda_minus_right;
    Real rho_h_right = rho_right + gamma_adi_red * pgas_right;
    Real v_sq_right = SQR(vx_right) + SQR(vy_right) + SQR(vz_right);
    Real gamma_sq_right = 1.0/(1.0-v_sq_right);
    pmy_fluid->pf_eos->SoundSpeedsSR(
        rho_h_right, pgas_right, vx_right, gamma_sq_right,
        &lambda_plus_right, &lambda_minus_right);                     // (MB 23)

    // Calculate extremal wavespeeds
    Real lambda_left = std::min(lambda_minus_left, lambda_minus_right);
    Real lambda_right = std::max(lambda_plus_left, lambda_plus_right);

    // Calculate L/R state fluxes
    Real flux_left[NWAVE], flux_right[NWAVE];
    PrimToFluxFlat(gamma_adi_red, rho_left, pgas_left, u_con_left,
        flux_left, ivx, ivy, ivz);
    PrimToFluxFlat(gamma_adi_red, rho_right, pgas_right, u_con_right,
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
    PrimToConsFlat(gamma_adi_red, rho_left, pgas_left, u_con_left,
        cons_left, ivx, ivy, ivz);
    PrimToConsFlat(gamma_adi_red, rho_right, pgas_right, u_con_right,
        cons_right, ivx, ivy, ivz);
    for (int n = 0; n < NWAVE; ++n)
      flux(n,i) = (lambda_right*flux_left[n] - lambda_left*flux_right[n]
          + lambda_right*lambda_left * (cons_right[n] - cons_left[n]))
          / (lambda_right-lambda_left);                                   // (MB 11)
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
//   implements (15) from Mignone, Ugliano, & Bodo 2009, MNRAS 393 1141
//   equivalent to (2) and (3) from Mignone & Bodo 2005, MNRAS 364 126
static void PrimToFluxFlat(Real gamma_adi_red, Real rho, Real pgas, const Real u_con[4],
    Real flux[NWAVE], int ivx, int ivy, int ivz)
{
  Real w = rho + gamma_adi_red * pgas;
  flux[IDN] = rho*u_con[1];
  flux[IEN] = w*u_con[0]*u_con[1];
  flux[ivx] = w*u_con[1]*u_con[1] + pgas;
  flux[ivy] = w*u_con[2]*u_con[1];
  flux[ivz] = w*u_con[3]*u_con[1];
  return;
}

// Function for converting primitive state to conserved state in flat spacetime
// Notes:
//   implements (3) from Mignone & Bodo 2005, MNRAS 364 126
static void PrimToConsFlat(Real gamma_adi_red, Real rho, Real pgas, const Real u_con[4],
    Real cons[NWAVE], int ivx, int ivy, int ivz)
{
  Real w = rho + gamma_adi_red * pgas;
  cons[IDN] = rho*u_con[0];
  cons[IEN] = w*u_con[0]*u_con[0] - pgas;
  cons[ivx] = w*u_con[1]*u_con[0];
  cons[ivy] = w*u_con[2]*u_con[0];
  cons[ivz] = w*u_con[3]*u_con[0];
  return;
}
