// HLLC Riemann solver for relativistic hydrodynamics

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

// Riemann solver
// Inputs:
//   k,j: x3- and x2-indices
//   il,iu: lower and upper x1-indices
//   ivx: type of interface (IVX for x1, IVY for x2, IVZ for x3)
//   b: 3D array of normal magnetic fields (unused)
//   prim_left, prim_right: left and right primitive states
// Outputs:
//   flux: fluxes across interface
// Notes:
//   implements HLLC algorithm from Mignone & Bodo 2005, MNRAS 364 126 (MB)
//   prim_left, prim_right overwritten
void FluidIntegrator::RiemannSolver(const int k,const int j, const int il, const int iu,
    const int ivx, const AthenaArray<Real> &b, AthenaArray<Real> &prim_left,
    AthenaArray<Real> &prim_right, AthenaArray<Real> &flux)
{
  // Calculate cyclic permutations of indices
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  // Extract ratio of specific heats
  const Real gamma_adi = pmy_fluid->pf_eos->GetGamma();
  const Real gamma_prime = gamma_adi / (gamma_adi - 1.0);

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

  // Go through each interface
  #pragma simd
  for (int i = il; i <= iu; i++)
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

    // Calculate wavespeeds in left state
    Real lambda_plus_left, lambda_minus_left;
    Real wgas_left = rho_left + gamma_prime * pgas_left;
    Real v_sq_left = SQR(vx_left) + SQR(vy_left) + SQR(vz_left);
    Real gamma_sq_left = 1.0/(1.0-v_sq_left);
    pmy_fluid->pf_eos->SoundSpeedsSR(
        wgas_left, pgas_left, vx_left, gamma_sq_left,
        &lambda_plus_left, &lambda_minus_left);                   // (MB 23)

    // Calculate wavespeeds in right state
    Real lambda_plus_right, lambda_minus_right;
    Real wgas_right = rho_right + gamma_prime * pgas_right;
    Real v_sq_right = SQR(vx_right) + SQR(vy_right) + SQR(vz_right);
    Real gamma_sq_right = 1.0/(1.0-v_sq_right);
    pmy_fluid->pf_eos->SoundSpeedsSR(
        wgas_right, pgas_right, vx_right, gamma_sq_right,
        &lambda_plus_right, &lambda_minus_right);                     // (MB 23)

    // Calculate extremal wavespeeds
    Real lambda_left = std::min(lambda_minus_left, lambda_minus_right);
    Real lambda_right = std::max(lambda_plus_left, lambda_plus_right);

    // Calculate conserved quantities in L state (MB 3)
    Real cons_left[NWAVE];
    cons_left[IDN] = std::sqrt(gamma_sq_left) * rho_left;
    cons_left[IEN] = gamma_sq_left * wgas_left - pgas_left;
    cons_left[ivx] = gamma_sq_left * wgas_left * vx_left;
    cons_left[ivy] = gamma_sq_left * wgas_left * vy_left;
    cons_left[ivz] = gamma_sq_left * wgas_left * vz_left;

    // Calculate conserved quantities in R state (MB 3)
    Real cons_right[NWAVE];
    cons_right[IDN] = std::sqrt(gamma_sq_right) * rho_right;
    cons_right[IEN] = gamma_sq_right * wgas_right - pgas_right;
    cons_right[ivx] = gamma_sq_right * wgas_right * vx_right;
    cons_right[ivy] = gamma_sq_right * wgas_right * vy_right;
    cons_right[ivz] = gamma_sq_right * wgas_right * vz_right;

    // Calculate fluxes in L state (MB 2)
    Real flux_left[NWAVE];
    flux_left[IDN] = cons_left[IDN] * vx_left;
    flux_left[IEN] = cons_left[ivx];
    flux_left[ivx] = cons_left[ivx] * vx_left + pgas_left;
    flux_left[ivy] = cons_left[ivy] * vx_left;
    flux_left[ivz] = cons_left[ivz] * vx_left;

    // Calculate fluxes in R state (MB 2)
    Real flux_right[NWAVE];
    flux_right[IDN] = cons_right[IDN] * vx_right;
    flux_right[IEN] = cons_right[ivx];
    flux_right[ivx] = cons_right[ivx] * vx_right + pgas_right;
    flux_right[ivy] = cons_right[ivy] * vx_right;
    flux_right[ivz] = cons_right[ivz] * vx_right;

    // Calculate select HLL conserved quantities (MB 9) and fluxes (MB 11)
    Real e_hll = (lambda_right*cons_right[IEN] - lambda_left*cons_left[IEN]
        + flux_left[IEN] - flux_right[IEN]) / (lambda_right-lambda_left);
    Real mx_hll = (lambda_right*cons_right[ivx] - lambda_left*cons_left[ivx]
        + flux_left[ivx] - flux_right[ivx]) / (lambda_right-lambda_left);
    Real flux_e_hll = (lambda_right*flux_left[IEN] - lambda_left*flux_right[IEN]
        + lambda_right*lambda_left * (cons_right[IEN]-cons_left[IEN]))
        / (lambda_right-lambda_left);
    Real flux_mx_hll = (lambda_right*flux_left[ivx] - lambda_left*flux_right[ivx]
        + lambda_right*lambda_left * (cons_right[ivx]-cons_left[ivx]))
        / (lambda_right-lambda_left);

    // Calculate contact wavespeed (MB 18)
    Real lambda_star;
    if (flux_e_hll > TINY_NUMBER || flux_e_hll < -TINY_NUMBER)  // use quadratic formula
    {
      // Follows algorithm in Numerical Recipes (section 5.6) for avoiding cancellations
      Real a = flux_e_hll;
      Real b = -(e_hll + flux_mx_hll);
      Real c = mx_hll;
      Real q = -0.5 * (b - std::sqrt(SQR(b) - 4.0*a*c));
      lambda_star = c / q;
    }
    else  // no quadratic term
      lambda_star = mx_hll / (e_hll + flux_mx_hll);

    // Calculate contact pressure (MB 17)
    Real a = lambda_left * cons_left[IEN] - cons_left[ivx];
    Real b = cons_left[ivx] * (lambda_left - vx_left) - pgas_left;
    Real pgas_star_left = (a * lambda_star - b) / (1.0 - lambda_left * lambda_star);
    a = lambda_right * cons_right[IEN] - cons_right[ivx];
    b = cons_right[ivx] * (lambda_right - vx_right) - pgas_right;
    Real pgas_star_right = (a * lambda_star - b) / (1.0 - lambda_right * lambda_star);
    Real pgas_star = 0.5 * (pgas_star_left + pgas_star_right);

    // Calculate conserved quantities in L* state (MB 16)
    Real cons_left_star[NWAVE];
    for (int n = 0; n < NWAVE; ++n)
      cons_left_star[n] = cons_left[n] * (lambda_left-vx_left);
    cons_left_star[IEN] += pgas_star*lambda_star - pgas_left*vx_left;
    cons_left_star[ivx] += pgas_star - pgas_left;
    for (int n = 0; n < NWAVE; ++n)
      cons_left_star[n] /= lambda_left - lambda_star;

    // Calculate conserved quantities in R* state (MB 16)
    Real cons_right_star[NWAVE];
    for (int n = 0; n < NWAVE; ++n)
      cons_right_star[n] = cons_right[n] * (lambda_right-vx_right);
    cons_right_star[IEN] += pgas_star*lambda_star - pgas_right*vx_right;
    cons_right_star[ivx] += pgas_star - pgas_right;
    for (int n = 0; n < NWAVE; ++n)
      cons_right_star[n] /= lambda_right - lambda_star;

    // Calculate fluxes in L* state
    Real flux_left_star[NWAVE];
    for (int n = 0; n < NWAVE; ++n)
      flux_left_star[n] = flux_left[n]
          + lambda_left * (cons_left_star[n] - cons_left[n]);

    // Calculate fluxes in R* state
    Real flux_right_star[NWAVE];
    for (int n = 0; n < NWAVE; ++n)
      flux_right_star[n] = flux_right[n]
          + lambda_right * (cons_right_star[n] - cons_right[n]);

    // Set fluxes
    for (int n = 0; n < NWAVE; ++n)
    {
      if (lambda_left >= 0.0)  // L state
        flux(n,i) = flux_left[n];
      else if (lambda_right <= 0.0)  // R state
        flux(n,i) = flux_right[n];
      else if (lambda_star >= 0.0)  // L* state
        flux(n,i) = flux_left_star[n];
      else  // R* state
        flux(n,i) = flux_right_star[n];
    }

    // Set conserved quantities in GR
    if (GENERAL_RELATIVITY)
      for (int n = 0; n < NWAVE; ++n)
      {
        if (lambda_left >= 0.0)  // L state
          cons_(n,i) = cons_left[n];
        else if (lambda_right <= 0.0)  // R state
          cons_(n,i) = cons_right[n];
        else if (lambda_star >= 0.0)  // L* state
          cons_(n,i) = cons_left_star[n];
        else  // R* state
          cons_(n,i) = cons_right_star[n];
      }
  }

  // Transform fluxes to global coordinates if in GR
  if (GENERAL_RELATIVITY)
    switch (ivx)
    {
      case IVX:
        pmy_fluid->pmy_block->pcoord->FluxToGlobal1(k, j, il, iu, cons_, b_normal_,
            flux);
        break;
      case IVY:
        pmy_fluid->pmy_block->pcoord->FluxToGlobal2(k, j, il, iu, cons_, b_normal_,
            flux);
        break;
      case IVZ:
        pmy_fluid->pmy_block->pcoord->FluxToGlobal3(k, j, il, iu, cons_, b_normal_,
            flux);
        break;
    }
  return;
}
