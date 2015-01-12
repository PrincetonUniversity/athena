// Local Lax-Friedrichs Riemann solver for relativistic hydrodynamics

// Primary header
#include "../../fluid_integrator.hpp"

// C++ headers
#include <cmath>      // sqrt()

// Athena headers
#include "../../../eos/eos.hpp"                     // GetGamma()
#include "../../../fluid.hpp"                       // Fluid
#include "../../../../athena.hpp"                   // enums, macros, Real
#include "../../../../athena_arrays.hpp"            // AthenaArray
#include "../../../../coordinates/coordinates.hpp"  // Coordinates
#include "../../../../mesh.hpp"                     // MeshBlock

// Declarations
static void PrimToConsFlat(Real gamma_adi_red, Real rho, Real pgas,
    Real ut, Real ux, Real uy, Real uz,
    int ivx, int ivy, int ivz, Real cons[NWAVE]);
static void PrimToFluxFlat(Real gamma_adi_red, Real rho, Real pgas,
    Real ut, Real ux, Real uy, Real uz,
    int ivx, int ivy, int ivz, Real flux[NWAVE]);

// Riemann solver
// Inputs:
//   il,iu: lower and upper indices for interfaces
//   prim_left, prim_right: left and right primitive states
// Outputs:
//   flux: fluxes
// Notes:
//   prim_left, prim_right overwritten
//   implements LLF scheme
//   equivalent to HLLE with outer wavespeeds set to c
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

    // Calculate covariant versions of left primitives
    Real ut_left = std::sqrt(1.0/(1.0-(SQR(vx_left)+SQR(vy_left)+SQR(vz_left))));
    Real ux_left = ut_left * vx_left;
    Real uy_left = ut_left * vy_left;
    Real uz_left = ut_left * vz_left;

    // Calculate covariant versions of right primitives
    Real ut_right = std::sqrt(1.0/(1.0-(SQR(vx_right)+SQR(vy_right)+SQR(vz_right))));
    Real ux_right = ut_right * vx_right;
    Real uy_right = ut_right * vy_right;
    Real uz_right = ut_right * vz_right;

    // Calculate L/R state conserved quantities
    Real cons_left[NWAVE], cons_right[NWAVE];
    PrimToConsFlat(gamma_adi_red, rho_left, pgas_left,
        ut_left, ux_left, uy_left, uz_left,
        ivx, ivy, ivz, cons_left);
    PrimToConsFlat(gamma_adi_red, rho_right, pgas_right,
        ut_right, ux_right, uy_right, uz_right,
        ivx, ivy, ivz, cons_right);

    // Calculate L/R state fluxes
    Real flux_left[NWAVE], flux_right[NWAVE];
    PrimToFluxFlat(gamma_adi_red, rho_left, pgas_left,
        ut_left, ux_left, uy_left, uz_left,
        ivx, ivy, ivz, flux_left);
    PrimToFluxFlat(gamma_adi_red, rho_right, pgas_right,
        ut_right, ux_right, uy_right, uz_right,
        ivx, ivy, ivz, flux_right);

    // Set fluxes
    for (int n = 0; n < NWAVE; ++n)
      flux(n,i) = 0.5 * (flux_left[n] + flux_right[n] - cons_right[n] + cons_left[n]);
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

// Function for converting primitive state to conserved state in flat spacetime
// Notes:
//   implements (3) from Mignone & Bodo 2005, MNRAS 364 126 (MB)
static void PrimToConsFlat(Real gamma_adi_red, Real rho, Real pgas,
    Real ut, Real ux, Real uy, Real uz,
    int ivx, int ivy, int ivz, Real cons[NWAVE])
{
  Real rho_h = rho + gamma_adi_red * pgas;
  cons[IDN] = ut * rho;
  cons[IEN] = rho_h * ut * ut - pgas;  // (MUB 8)
  cons[ivx] = rho_h * ux * ut;         // (MUB 8)
  cons[ivy] = rho_h * uy * ut;         // (MUB 8)
  cons[ivz] = rho_h * uz * ut;         // (MUB 8)
  return;
}

// Function for converting constant primitive state to flux state in flat spacetime
// Notes:
//   implements (2) from Mignone & Bodo 2005, MNRAS 364 126 (MB)
static void PrimToFluxFlat(Real gamma_adi_red, Real rho, Real pgas,
    Real ut, Real ux, Real uy, Real uz,
    int ivx, int ivy, int ivz, Real flux[NWAVE])
{
  Real rho_h = rho + gamma_adi_red * pgas;
  flux[IDN] = ux * rho;
  flux[IEN] = rho_h * ut * ux;
  flux[ivx] = rho_h * ux * ux + pgas;
  flux[ivy] = rho_h * uy * ux;
  flux[ivz] = rho_h * uz * ux;
  return;
}
