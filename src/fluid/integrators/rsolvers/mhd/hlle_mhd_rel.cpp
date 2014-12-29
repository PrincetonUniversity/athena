// HLLE Riemann solver for relativistic magnetohydrodynamics

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
static void PrimToFluxFlat(Real gamma_adi_red, Real rho, Real pgas,
    Real ut, Real ux, Real uy, Real uz,
    Real bcovt, Real bcovx, Real bcovy, Real bcovz,
    int ivx, int ivy, int ivz, Real flux[NWAVE]);
static void PrimToConsFlat(Real gamma_adi_red, Real rho, Real pgas,
    Real ut, Real ux, Real uy, Real uz,
    Real bcovt, Real bcovx, Real bcovy, Real bcovz,
    int ivx, int ivy, int ivz, Real cons[NWAVE]);

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
  // Calculate cyclic permutations of indices
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

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

  // Extract ratio of specific heats
  const Real gamma_adi = pmy_fluid->pf_eos->GetGamma();
  const Real gamma_adi_red = gamma_adi/(gamma_adi-1.0);

  // Go through each interface
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract left primitives
    Real &rho_left = prim_left(IDN,i);
    Real &pgas_left = prim_left(IEN,i);
    Real &vx_left = prim_left(ivx,i);
    Real &vy_left = prim_left(ivy,i);
    Real &vz_left = prim_left(ivz,i);
    Real &by_left = prim_left(IBY,i);
    Real &bz_left = prim_left(IBZ,i);

    // Extract right primitives
    Real &rho_right = prim_right(IDN,i);
    Real &pgas_right = prim_right(IEN,i);
    Real &vx_right = prim_right(ivx,i);
    Real &vy_right = prim_right(ivy,i);
    Real &vz_right = prim_right(ivz,i);
    Real &by_right = prim_right(IBY,i);
    Real &bz_right = prim_right(IBZ,i);

    // Extract normal magnetic field
    Real &bx = b_normal_(k,j,i);

    // Calculate covariant versions of left primitives
    Real ut_left = std::sqrt(1.0/(1.0-(SQR(vx_left)+SQR(vy_left)+SQR(vz_left))));
    Real ux_left = ut_left * vx_left;
    Real uy_left = ut_left * vy_left;
    Real uz_left = ut_left * vz_left;
    Real bcovt_left = bx*ux_left + by_left*uy_left + bz_left*uz_left;
    Real bcovx_left = (bx + bcovt_left * ux_left) / ut_left;
    Real bcovy_left = (by_left + bcovt_left * uy_left) / ut_left;
    Real bcovz_left = (bz_left + bcovt_left * uz_left) / ut_left;

    // Calculate covariant versions of right primitives
    Real ut_right = std::sqrt(1.0/(1.0-(SQR(vx_right)+SQR(vy_right)+SQR(vz_right))));
    Real ux_right = ut_right * vx_right;
    Real uy_right = ut_right * vy_right;
    Real uz_right = ut_right * vz_right;
    Real bcovt_right = bx*ux_right + by_right*uy_right + bz_right*uz_right;
    Real bcovx_right = (bx + bcovt_right * ux_right) / ut_right;
    Real bcovy_right = (by_right + bcovt_right * uy_right) / ut_right;
    Real bcovz_right = (bz_right + bcovt_right * uz_right) / ut_right;

    // Calculate wavespeeds
    Real lambda_left_plus, lambda_left_minus;
    pmy_fluid->pf_eos->FastMagnetosonicSpeedsRel(
        rho_left, pgas_left,
        vx_left, vy_left, vz_left,
        ut_left, ux_left, uy_left, uz_left,
        bx, by_left, bz_left,
        bcovt_left, bcovx_left, bcovy_left, bcovz_left,
        &lambda_left_plus, &lambda_left_minus);                          // (MB2006 56)
    Real lambda_right_plus, lambda_right_minus;
    pmy_fluid->pf_eos->FastMagnetosonicSpeedsRel(
        rho_right, pgas_right,
        vx_right, vy_right, vz_right,
        ut_right, ux_right, uy_right, uz_right,
        bx, by_right, bz_right,
        bcovt_right, bcovx_right, bcovy_right, bcovz_right,
        &lambda_right_plus, &lambda_right_minus);                        // (MB2006 56)
    Real lambda_left = std::min(lambda_left_minus, lambda_right_minus);  // (MB2006 55)
    Real lambda_right = std::max(lambda_left_plus, lambda_right_plus);   // (MB2006 55)

    // Calculate L/R state fluxes
    Real flux_left[NWAVE], flux_right[NWAVE];
    PrimToFluxFlat(gamma_adi_red, rho_left, pgas_left,
        ut_left, ux_left, uy_left, uz_left,
        bcovt_left, bcovx_left, bcovy_left, bcovz_left,
        ivx, ivy, ivz, flux_left);
    PrimToFluxFlat(gamma_adi_red, rho_right, pgas_right,
        ut_right, ux_right, uy_right, uz_right,
        bcovt_right, bcovx_right, bcovy_right, bcovz_right,
        ivx, ivy, ivz, flux_right);

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
        ut_left, ux_left, uy_left, uz_left,
        bcovt_left, bcovx_left, bcovy_left, bcovz_left,
        ivx, ivy, ivz, cons_left);
    PrimToConsFlat(gamma_adi_red, rho_right, pgas_right,
        ut_right, ux_right, uy_right, uz_right,
        bcovt_right, bcovx_right, bcovy_right, bcovz_right,
        ivx, ivy, ivz, cons_right);
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

// Function for converting constant primitive state to flux state in flat spacetime
// Notes:
//   implements (15) from Mignone, Ugliano, & Bodo 2009, MNRAS 393 1141
//     note B^i v^x - B^x v^i = b^i u^x - b^x u^i
static void PrimToFluxFlat(Real gamma_adi_red, Real rho, Real pgas,
    Real ut, Real ux, Real uy, Real uz,
    Real bcovt, Real bcovx, Real bcovy, Real bcovz,
    int ivx, int ivy, int ivz, Real flux[NWAVE])
{
  Real bcov_sq = -SQR(bcovt) + SQR(bcovx) + SQR(bcovy) + SQR(bcovz);
  Real rho_h = rho + gamma_adi_red * pgas + bcov_sq;
  Real ptot = pgas + 0.5*bcov_sq;
  flux[IDN] = ux * rho;
  flux[IEN] = rho_h * ut * ux - bcovt * bcovx;
  flux[ivx] = rho_h * ux * ux - bcovx * bcovx + ptot;
  flux[ivy] = rho_h * uy * ux - bcovy * bcovx;
  flux[ivz] = rho_h * uz * ux - bcovz * bcovx;
  flux[IBY] = bcovy * ux - bcovx * uy;
  flux[IBZ] = bcovz * ux - bcovx * uz;
  return;
}

// Function for converting primitive state to conserved state in flat spacetime
// Notes:
//   references Mignone, Ugliano, & Bodo 2009, MNRAS 393 1141 (MUB)
static void PrimToConsFlat(Real gamma_adi_red, Real rho, Real pgas,
    Real ut, Real ux, Real uy, Real uz,
    Real bcovt, Real bcovx, Real bcovy, Real bcovz,
    int ivx, int ivy, int ivz, Real cons[NWAVE])
{
  Real bcov_sq = -SQR(bcovt) + SQR(bcovx) + SQR(bcovy) + SQR(bcovz);
  Real rho_h = rho + gamma_adi_red * pgas + bcov_sq;
  Real ptot = pgas + 0.5*bcov_sq;
  cons[IDN] = ut * rho;
  cons[IEN] = rho_h * ut * ut - bcovt * bcovt - ptot;  // (MUB 8)
  cons[ivx] = rho_h * ux * ut - bcovx * bcovt;         // (MUB 8)
  cons[ivy] = rho_h * uy * ut - bcovy * bcovt;         // (MUB 8)
  cons[ivz] = rho_h * uz * ut - bcovz * bcovt;         // (MUB 8)
  cons[IBY] = bcovy * ut - bcovt * uy;
  cons[IBZ] = bcovz * ut - bcovt * uz;
  return;
}
