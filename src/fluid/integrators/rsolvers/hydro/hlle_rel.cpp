// HLLE Riemann solver for relativistic hydrodynamics

// Primary header
#include "../../fluid_integrator.hpp"

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()

// Athena headers
#include "../../../eos/eos.hpp"                     // GetGamma(), WavespeedsRel()
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
//   b: 3D array of normal magnetic fields (not used)
//   prim_left, prim_right: left and right primitive states
// Outputs:
//   flux: fluxes across interface
// Notes:
//   implements HLLE algorithm from Mignone & Bodo 2005, MNRAS 364 126 (MB)
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
  const Real gamma_adi_red = gamma_adi / (gamma_adi - 1.0);

  // Transform primitives to locally flat coordinates if in GR
  if (GENERAL_RELATIVITY)
    switch (ivx)
    {
      case IVX:
        pmy_fluid->pmy_block->pcoord->PrimToLocal1(k, j, b, prim_left, prim_right,
            b_normal_);
        break;
      case IVY:
        pmy_fluid->pmy_block->pcoord->PrimToLocal2(k, j, b, prim_left, prim_right,
            b_normal_);
        break;
      case IVZ:
        pmy_fluid->pmy_block->pcoord->PrimToLocal3(k, j, b, prim_left, prim_right,
            b_normal_);
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

    // Extract fluxes
    Real &flux_d = flux(IDN,i);
    Real &flux_e = flux(IEN,i);
    Real &flux_mx = flux(ivx,i);
    Real &flux_my = flux(ivy,i);
    Real &flux_mz = flux(ivz,i);

    // Calculate wavespeeds in left region
    Real lambda_plus_left, lambda_minus_left;
    Real rho_h_left = rho_left + gamma_adi_red * pgas_left;
    Real v_sq_left = SQR(vx_left) + SQR(vy_left) + SQR(vz_left);
    Real gamma_sq_left = 1.0/(1.0-v_sq_left);
    pmy_fluid->pf_eos->WavespeedsRel(
        rho_h_left, pgas_left, vx_left, gamma_sq_left,
        &lambda_plus_left, &lambda_minus_left);                   // (MB 23)

    // Calculate wavespeeds in right region
    Real lambda_plus_right, lambda_minus_right;
    Real rho_h_right = rho_right + gamma_adi_red * pgas_right;
    Real v_sq_right = SQR(vx_right) + SQR(vy_right) + SQR(vz_right);
    Real gamma_sq_right = 1.0/(1.0-v_sq_right);
    pmy_fluid->pf_eos->WavespeedsRel(
        rho_h_right, pgas_right, vx_right, gamma_sq_right,
        &lambda_plus_right, &lambda_minus_right);                     // (MB 23)

    // Calculate extremal wavespeeds
    Real lambda_left = std::min(lambda_minus_left, lambda_minus_right);
    Real lambda_right = std::max(lambda_plus_left, lambda_plus_right);

    // Set fluxes in extremal cases
    if (lambda_left >= 0.0)  // left region
    {
      // Calculate intermediate quantities (MB 3)
      Real gamma_sq_rho_h = gamma_sq_left * rho_h_left;
      Real d = std::sqrt(gamma_sq_left) * rho_left;
      Real mx = gamma_sq_rho_h * vx_left;
      Real my = gamma_sq_rho_h * vy_left;
      Real mz = gamma_sq_rho_h * vz_left;

      // Set fluxes (MB 2)
      flux_d = d * vx_left;
      flux_e = mx;
      flux_mx = mx * vx_left + pgas_left;
      flux_my = my * vx_left;
      flux_mz = mz * vx_left;
      continue;
    }
    if (lambda_right <= 0.0)  // right region
    {
      // Calculate intermediate quantities (MB 3)
      Real gamma_sq_rho_h = gamma_sq_right * rho_h_right;
      Real d = std::sqrt(gamma_sq_right) * rho_right;
      Real mx = gamma_sq_rho_h * vx_right;
      Real my = gamma_sq_rho_h * vy_right;
      Real mz = gamma_sq_rho_h * vz_right;

      // Set fluxes (MB 2)
      flux_d = d * vx_right;
      flux_e = mx;
      flux_mx = mx * vx_right + pgas_right;
      flux_my = my * vx_right;
      flux_mz = mz * vx_right;
      continue;
    }

    // Calculate left conserved quantities (MB 3)
    Real gamma_sq_rho_h_left = gamma_sq_left * rho_h_left;
    Real d_left = std::sqrt(gamma_sq_left) * rho_left;
    Real e_left = gamma_sq_rho_h_left - pgas_left;
    Real mx_left = gamma_sq_rho_h_left * vx_left;
    Real my_left = gamma_sq_rho_h_left * vy_left;
    Real mz_left = gamma_sq_rho_h_left * vz_left;

    // Calculate right conserved quantities (MB 3)
    Real gamma_sq_rho_h_right = gamma_sq_right * rho_h_right;
    Real d_right = std::sqrt(gamma_sq_right) * rho_right;
    Real e_right = gamma_sq_rho_h_right - pgas_right;
    Real mx_right = gamma_sq_rho_h_right * vx_right;
    Real my_right = gamma_sq_rho_h_right * vy_right;
    Real mz_right = gamma_sq_rho_h_right * vz_right;

    // Calculate left fluxes (MB 2)
    Real flux_d_left = d_left * vx_left;
    Real flux_e_left = mx_left;
    Real flux_mx_left = mx_left * vx_left + pgas_left;
    Real flux_my_left = my_left * vx_left;
    Real flux_mz_left = mz_left * vx_left;

    // Calculate right fluxes (MB 2)
    Real flux_d_right = d_right * vx_right;
    Real flux_e_right = mx_right;
    Real flux_mx_right = mx_right * vx_right + pgas_right;
    Real flux_my_right = my_right * vx_right;
    Real flux_mz_right = mz_right * vx_right;

    // Set fluxes (MB 11)
    Real denom_inverse = 1.0 / (lambda_right - lambda_left);
    flux_d = (lambda_right * flux_d_left - lambda_left * flux_d_right + lambda_right
        * lambda_left * (d_right - d_left)) * denom_inverse;
    flux_e = (lambda_right * flux_e_left - lambda_left * flux_e_right + lambda_right
        * lambda_left * (e_right - e_left)) * denom_inverse;
    flux_mx = (lambda_right * flux_mx_left - lambda_left * flux_mx_right + lambda_right
        * lambda_left * (mx_right - mx_left)) * denom_inverse;
    flux_my = (lambda_right * flux_my_left - lambda_left * flux_my_right + lambda_right
        * lambda_left * (my_right - my_left)) * denom_inverse;
    flux_mz = (lambda_right * flux_mz_left - lambda_left * flux_mz_right + lambda_right
        * lambda_left * (mz_right - mz_left)) * denom_inverse;
  }

  // Transform fluxes to global coordinates if in GR
  if (GENERAL_RELATIVITY)
    switch (ivx)
    {
      case IVX:
        pmy_fluid->pmy_block->pcoord->FluxToGlobal1(k, j, flux);
        break;
      case IVY:
        pmy_fluid->pmy_block->pcoord->FluxToGlobal2(k, j, flux);
        break;
      case IVZ:
        pmy_fluid->pmy_block->pcoord->FluxToGlobal3(k, j, flux);
        break;
    }
  return;
}
