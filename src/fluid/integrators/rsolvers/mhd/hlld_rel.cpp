// HLLD Riemann solver for relativistic magnetohydrodynamics

// Primary header
#include "../../fluid_integrator.hpp"

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // abs(), isfinite(), NAN, sqrt()

// Athena headers
#include "../../../fluid.hpp"                       // Fluid
#include "../../../eos/eos.hpp"                     // FluidEqnOfState
#include "../../../../athena.hpp"                   // enums, macros, Real
#include "../../../../athena_arrays.hpp"            // AthenaArray
#include "../../../../mesh.hpp"                     // MeshBlock
#include "../../../../coordinates/coordinates.hpp"  // Coordinates

// Declarations
static Real ConsToPFlat(const Real cons[NWAVE], Real bbx, Real gamma_adi, int ivx);
static Real EResidual(Real w_guess, Real d, Real e, Real m_sq, Real b_sq, Real s_sq,
    Real gamma_adi);
static Real EResidualPrime(Real w_guess, Real d, Real m_sq, Real b_sq, Real s_sq,
    Real gamma_adi);
static Real PResidual(Real p, Real bbx, Real lambda_l, Real lambda_r,
    const Real r_l[NWAVE], const Real r_r[NWAVE], int ivx);
static Real FindRootNR(Real w_init, Real d, Real e, Real m_sq, Real b_sq, Real s_sq,
    Real gamma_adi);
static Real FindRootSecant(Real ptot_init, Real bbx, Real lambda_l, Real lambda_r,
    const Real r_l[NWAVE], const Real r_r[NWAVE], int ivx);
static Real quadratic_root(Real a1, Real a0, bool greater_root);

// Riemann solver
// Inputs:
//   k,j: x3- and x2-indices
//   il,iu: lower and upper x1-indices
//   ivx: type of interface (IVX for x1, IVY for x2, IVZ for x3)
//   bb: 3D array of normal magnetic fields
//   prim_l, prim_r: left and right primitive states
// Outputs:
//   flux: fluxes across interface
// Notes:
//   prim_l, prim_r overwritten
//   implements HLLD solver from Mignone, Ugliano, & Bodo 2009, MNRAS 393 1141 (MUB)
//   references Mignone & Bodo 2006, MNRAS 368 1040 (MB)
//   follows Athena 4.2, hlld_sr.c, in variable choices and magic numbers
void FluidIntegrator::RiemannSolver(const int k, const int j, const int il,
    const int iu, const int ivx, const AthenaArray<Real> &bb, AthenaArray<Real> &prim_l,
    AthenaArray<Real> &prim_r, AthenaArray<Real> &flux)
{
  // Parameters
  const Real p_transition = 0.01;     // value delineating intial pressure regimes
  const Real vc_extension = 1.0e-6;   // use contact region if Alfven speeds smaller
  const Real delta_kx_aug = 1.0e-12;  // amount to add to \Delta K^x

  // Calculate metric if in GR
  int i01, i11;
  if (GENERAL_RELATIVITY)
    switch (ivx)
    {
      case IVX:
        pmy_fluid->pmy_block->pcoord->Face1Metric(k, j, il, iu, g_, gi_);
        i01 = I01;
        i11 = I11;
        break;
      case IVY:
        pmy_fluid->pmy_block->pcoord->Face2Metric(k, j, il, iu, g_, gi_);
        i01 = I02;
        i11 = I22;
        break;
      case IVZ:
        pmy_fluid->pmy_block->pcoord->Face3Metric(k, j, il, iu, g_, gi_);
        i01 = I03;
        i11 = I33;
        break;
    }

  // Transform primitives to locally flat coordinates if in GR
  if (GENERAL_RELATIVITY)
    switch (ivx)
    {
      case IVX:
        pmy_fluid->pmy_block->pcoord->PrimToLocal1(k, j, il, iu, bb, prim_l, prim_r,
            bb_normal_);
        break;
      case IVY:
        pmy_fluid->pmy_block->pcoord->PrimToLocal2(k, j, il, iu, bb, prim_l, prim_r,
            bb_normal_);
        break;
      case IVZ:
        pmy_fluid->pmy_block->pcoord->PrimToLocal3(k, j, il, iu, bb, prim_l, prim_r,
            bb_normal_);
        break;
    }
  else  // SR; need to populate 1D normal B array
  {
    #pragma simd
    for (int i = il; i <= iu; ++i)
      bb_normal_(i) = bb(k,j,i);
  }

  // Calculate cyclic permutations of indices
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  // Extract ratio of specific heats
  const Real gamma_adi = pmy_fluid->pf_eos->GetGamma();

  // Go through each interface
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Calculate interface velocity
    Real v_interface = 0.0;
    if (GENERAL_RELATIVITY)
      v_interface = gi_(i01,i) / std::sqrt(SQR(gi_(i01,i)) - gi_(I00,i)*gi_(i11,i));

    // Extract left primitives
    const Real &rho_l = prim_l(IDN,i);
    const Real &pgas_l = prim_l(IEN,i);
    Real u_l[4];
    if (GENERAL_RELATIVITY)
    {
      u_l[1] = prim_l(ivx,i);
      u_l[2] = prim_l(ivy,i);
      u_l[3] = prim_l(ivz,i);
      u_l[0] = std::sqrt(1.0 + SQR(u_l[1]) + SQR(u_l[2]) + SQR(u_l[3]));
    }
    else  // SR
    {
      const Real &vx_l = prim_l(ivx,i);
      const Real &vy_l = prim_l(ivy,i);
      const Real &vz_l = prim_l(ivz,i);
      u_l[0] = std::sqrt(1.0 / (1.0 - SQR(vx_l) - SQR(vy_l) - SQR(vz_l)));
      u_l[1] = u_l[0] * vx_l;
      u_l[2] = u_l[0] * vy_l;
      u_l[3] = u_l[0] * vz_l;
    }
    const Real &bby_l = prim_l(IBY,i);
    const Real &bbz_l = prim_l(IBZ,i);

    // Extract right primitives
    const Real &rho_r = prim_r(IDN,i);
    const Real &pgas_r = prim_r(IEN,i);
    Real u_r[4];
    if (GENERAL_RELATIVITY)
    {
      u_r[1] = prim_r(ivx,i);
      u_r[2] = prim_r(ivy,i);
      u_r[3] = prim_r(ivz,i);
      u_r[0] = std::sqrt(1.0 + SQR(u_r[1]) + SQR(u_r[2]) + SQR(u_r[3]));
    }
    else  // SR
    {
      const Real &vx_r = prim_r(ivx,i);
      const Real &vy_r = prim_r(ivy,i);
      const Real &vz_r = prim_r(ivz,i);
      u_r[0] = std::sqrt(1.0 / (1.0 - SQR(vx_r) - SQR(vy_r) - SQR(vz_r)));
      u_r[1] = u_r[0] * vx_r;
      u_r[2] = u_r[0] * vy_r;
      u_r[3] = u_r[0] * vz_r;
    }
    const Real &bby_r = prim_r(IBY,i);
    const Real &bbz_r = prim_r(IBZ,i);

    // Extract normal magnetic field
    const Real &bbx = bb(k,j,i);

    // Calculate 4-magnetic field in left state
    Real b_l[4];
    b_l[0] = bbx*u_l[1] + bby_l*u_l[2] + bbz_l*u_l[3];
    b_l[1] = (bbx + b_l[0] * u_l[1]) / u_l[0];
    b_l[2] = (bby_l + b_l[0] * u_l[2]) / u_l[0];
    b_l[3] = (bbz_l + b_l[0] * u_l[3]) / u_l[0];
    Real b_sq_l = -SQR(b_l[0]) + SQR(b_l[1]) + SQR(b_l[2]) + SQR(b_l[3]);

    // Calculate 4-magnetic field in right state
    Real b_r[4];
    b_r[0] = bbx*u_r[1] + bby_r*u_r[2] + bbz_r*u_r[3];
    b_r[1] = (bbx + b_r[0] * u_r[1]) / u_r[0];
    b_r[2] = (bby_r + b_r[0] * u_r[2]) / u_r[0];
    b_r[3] = (bbz_r + b_r[0] * u_r[3]) / u_r[0];
    Real b_sq_r = -SQR(b_r[0]) + SQR(b_r[1]) + SQR(b_r[2]) + SQR(b_r[3]);

    // Calculate wavespeeds in left state (MB 56)
    Real lambda_p_l, lambda_m_l;
    pmy_fluid->pf_eos->FastMagnetosonicSpeedsSR(rho_l, pgas_l, u_l, b_l, &lambda_p_l,
        &lambda_m_l);

    // Calculate wavespeeds in right state (MB 56)
    Real lambda_p_r, lambda_m_r;
    pmy_fluid->pf_eos->FastMagnetosonicSpeedsSR(rho_r, pgas_r, u_r, b_r, &lambda_p_r,
        &lambda_m_r);

    // Calculate extremal wavespeeds (MB 55)
    Real lambda_l = std::min(lambda_m_l, lambda_m_r);
    Real lambda_r = std::max(lambda_p_l, lambda_p_r);

    // Calculate conserved quantities in L region (MUB 8)
    Real cons_l[NWAVE];
    Real wtot_l = rho_l + gamma_adi/(gamma_adi-1.0) * pgas_l + b_sq_l;
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
    Real wtot_r = rho_r + gamma_adi/(gamma_adi-1.0) * pgas_r + b_sq_r;
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

    // Calculate jump quantities across left fast wave (MUB 12)
    Real r_l[NWAVE];
    for (int n = 0; n < NWAVE; ++n)
      r_l[n] = lambda_l * cons_l[n] - flux_l[n];

    // Calculate jump quantities across right fast wave (MUB 12)
    Real r_r[NWAVE];
    for (int n = 0; n < NWAVE; ++n)
      r_r[n] = lambda_r * cons_r[n] - flux_r[n];

    // Calculate conserved quantities in HLL region (MB 29)
    Real cons_hll[NWAVE];
    for (int n = 0; n < NWAVE; ++n)
      cons_hll[n] = (r_r[n]-r_l[n]) / (lambda_r-lambda_l);

    // Calculate fluxes in HLL region (MB 31)
    Real flux_hll[NWAVE];
    for (int n = 0; n < NWAVE; ++n)
      flux_hll[n] = (lambda_l*r_r[n] - lambda_r*r_l[n]) / (lambda_r-lambda_l);

    // Calculate total pressure in HLL region
    Real ptot_hll = ConsToPFlat(cons_hll, bbx, gamma_adi, ivx);

    // Calculate initial guess for total pressure (MUB 53)
    Real ptot_init;
    if (SQR(bbx)/ptot_hll < p_transition)  // weak magnetic field
    {
      Real a1 = cons_hll[IEN] - flux_hll[ivx];
      Real a0 = cons_hll[ivx]*flux_hll[IEN] - flux_hll[ivx]*cons_hll[IEN];
      ptot_init= quadratic_root(a1, a0, true);                              // (MUB 55)
    }
    else  // strong magnetic field
      ptot_init = ptot_hll;

    // Apply secant method to find total pressure
    Real ptot_true = FindRootSecant(ptot_init, bbx, lambda_l, lambda_r, r_l, r_r, ivx);

    // Calculate velocity in aL region
    Real al = r_l[ivx] - lambda_l*r_l[IEN] + ptot_true*(1.0-SQR(lambda_l));  // (MUB 26)
    Real gl = SQR(r_l[IBY]) + SQR(r_l[IBZ]);                                 // (MUB 27)
    Real cl = r_l[ivy]*r_l[IBY] + r_l[ivz]*r_l[IBZ];                         // (MUB 28)
    Real ql = -al - gl + SQR(bbx)*(1.0-SQR(lambda_l));                       // (MUB 29)
    Real xl = bbx * (al*lambda_l*bbx+cl)
        - (al+gl) * (lambda_l*ptot_true+r_l[IEN]);                           // (MUB 30)
    Real vx_al = (bbx * (al*bbx+lambda_l*cl)
        - (al+gl) * (ptot_true+r_l[ivx])) / xl;                              // (MUB 23)
    Real vy_al = (ql*r_l[ivy]
        + r_l[IBY] * (cl + bbx * (lambda_l*r_l[ivx]-r_l[IEN]))) / xl;        // (MUB 24)
    Real vz_al = (ql*r_l[ivz]
        + r_l[IBZ] * (cl + bbx * (lambda_l*r_l[ivx]-r_l[IEN]))) / xl;        // (MUB 24)

    // Calculate auxiliary quantites in aL region
    Real v_rm_l = vx_al*r_l[ivx] + vy_al*r_l[ivy]
        + vz_al*r_l[ivz];
    Real wtot_al = ptot_true + (r_l[IEN]-v_rm_l)
        / (lambda_l-vx_al);                                     // (MUB 31)
    Real eta_l = -copysign(std::sqrt(wtot_al), bbx);            // (MUB 35)
    Real denom_al = lambda_l*ptot_true + r_l[IEN] + bbx*eta_l;
    Real kx_l = (r_l[ivx] + ptot_true + lambda_l*bbx*eta_l)     // R_{B^x} = \lambda B^x
        / denom_al;                                             // (MUB 43)
    Real ky_l = (r_l[ivy] + r_l[IBY]*eta_l) / denom_al;         // (MUB 43)
    Real kz_l = (r_l[ivz] + r_l[IBZ]*eta_l) / denom_al;         // (MUB 43)
    Real &lambda_al = kx_l;

    // Calculate conserved quantities in aL region
    Real cons_al[NWAVE];
    cons_al[IBY] = (r_l[IBY] - bbx*vy_al) / (lambda_l-vx_al);            // (MUB 21)
    cons_al[IBZ] = (r_l[IBZ] - bbx*vz_al) / (lambda_l-vx_al);            // (MUB 21)
    Real v_bb_al = vx_al*bbx + vy_al*cons_al[IBY] + vz_al*cons_al[IBZ];
    cons_al[IDN] = r_l[IDN] / (lambda_l-vx_al);                          // (MUB 32)
    cons_al[IEN] = (r_l[IEN] + ptot_true*vx_al - v_bb_al*bbx)
        / (lambda_l-vx_al);                                              // (MUB 33)
    cons_al[ivx] = (cons_al[IEN] + ptot_true) * vx_al - v_bb_al * bbx;   // (MUB 34)
    cons_al[ivy] = (cons_al[IEN] + ptot_true) * vy_al
        - v_bb_al * cons_al[IBY];                                        // (MUB 34)
    cons_al[ivz] = (cons_al[IEN] + ptot_true) * vz_al
        - v_bb_al * cons_al[IBZ];                                        // (MUB 34)

    // Calculate fluxes in aL region (MUB 11,12)
    Real flux_al[NWAVE];
    for (int n = 0; n < NWAVE; ++n)
      flux_al[n] = lambda_l * cons_al[n] - r_l[n];

    // Calculate velocity in aR region
    Real ar = r_r[ivx] - lambda_r*r_r[IEN] + ptot_true*(1.0-SQR(lambda_r));  // (MUB 26)
    Real gr = SQR(r_r[IBY]) + SQR(r_r[IBZ]);                                 // (MUB 27)
    Real cr = r_r[ivy]*r_r[IBY] + r_r[ivz]*r_r[IBZ];                         // (MUB 28)
    Real qr = -ar - gr + SQR(bbx)*(1.0-SQR(lambda_r));                       // (MUB 29)
    Real xr = bbx * (ar*lambda_r*bbx+cr)
        - (ar+gr) * (lambda_r*ptot_true+r_r[IEN]);                           // (MUB 30)
    Real vx_ar = (bbx * (ar*bbx+lambda_r*cr)
        - (ar+gr) * (ptot_true+r_r[ivx])) / xr;                              // (MUB 23)
    Real vy_ar = (qr*r_r[ivy]
        + r_r[IBY] * (cr + bbx * (lambda_r*r_r[ivx]-r_r[IEN]))) / xr;        // (MUB 24)
    Real vz_ar = (qr*r_r[ivz]
        + r_r[IBZ] * (cr + bbx * (lambda_r*r_r[ivx]-r_r[IEN]))) / xr;        // (MUB 24)

    // Calculate auxiliary quantites in aR region
    Real v_rm_r = vx_ar*r_r[ivx] + vy_ar*r_r[ivy]
        + vz_ar*r_r[ivz];
    Real wtot_ar = ptot_true + (r_r[IEN]-v_rm_r)
        / (lambda_r-vx_ar);                                     // (MUB 31)
    Real eta_r = copysign(std::sqrt(wtot_ar), bbx);             // (MUB 35)
    Real denom_ar = lambda_r*ptot_true + r_r[IEN] + bbx*eta_r;
    Real kx_r = (r_r[ivx] + ptot_true + lambda_r*bbx*eta_r)     // R_{B^x} = \lambda B^x
        / denom_ar;                                             // (MUB 43)
    Real ky_r = (r_r[ivy] + r_r[IBY]*eta_r) / denom_ar;         // (MUB 43)
    Real kz_r = (r_r[ivz] + r_r[IBZ]*eta_r) / denom_ar;         // (MUB 43)
    Real &lambda_ar = kx_r;

    // Calculate conserved quantities in aR region
    Real cons_ar[NWAVE];
    cons_ar[IBY] = (r_r[IBY] - bbx*vy_ar) / (lambda_r-vx_ar);            // (MUB 21)
    cons_ar[IBZ] = (r_r[IBZ] - bbx*vz_ar) / (lambda_r-vx_ar);            // (MUB 21)
    Real v_bb_ar = vx_ar*bbx + vy_ar*cons_ar[IBY] + vz_ar*cons_ar[IBZ];
    cons_ar[IDN] = r_r[IDN] / (lambda_r-vx_ar);                          // (MUB 32)
    cons_ar[IEN] = (r_r[IEN] + ptot_true*vx_ar - v_bb_ar*bbx)
        / (lambda_r-vx_ar);                                              // (MUB 33)
    cons_ar[ivx] = (cons_ar[IEN] + ptot_true) * vx_ar - v_bb_ar * bbx;   // (MUB 34)
    cons_ar[ivy] = (cons_ar[IEN] + ptot_true) * vy_ar
        - v_bb_ar * cons_ar[IBY];                                        // (MUB 34)
    cons_ar[ivz] = (cons_ar[IEN] + ptot_true) * vz_ar
        - v_bb_ar * cons_ar[IBZ];                                        // (MUB 34)

    // Calculate fluxes in aR region (MUB 11,12)
    Real flux_ar[NWAVE];
    for (int n = 0; n < NWAVE; ++n)
      flux_ar[n] = lambda_r * cons_ar[n] - r_r[n];

    // Calculate B_c (MUB 45)
    Real cons_c[NWAVE];
    Real denom_c = lambda_ar - lambda_al + delta_kx_aug;
    Real numer_al = cons_al[IBY] * (lambda_al-vx_al) + bbx*vy_al;
    Real numer_ar = cons_ar[IBY] * (lambda_ar-vx_ar) + bbx*vy_ar;
    cons_c[IBY] = (numer_ar - numer_al) / denom_c;
    numer_al = cons_al[IBZ] * (lambda_al-vx_al) + bbx*vz_al;
    numer_ar = cons_ar[IBZ] * (lambda_ar-vx_ar) + bbx*vz_ar;
    cons_c[IBZ] = (numer_ar - numer_al) / denom_c;

    // Calculate v_c (MUB 47), averaging left and right values in case of disagreement
    Real k_sq = SQR(kx_l) + SQR(ky_l) + SQR(kz_l);
    Real k_bc = kx_l*bbx + ky_l*cons_c[IBY] + kz_l*cons_c[IBZ];
    Real vx_cl = kx_l - bbx * (1.0-k_sq) / (eta_l-k_bc);
    Real vy_cl = ky_l - cons_c[IBY] * (1.0-k_sq) / (eta_l-k_bc);
    Real vz_cl = kz_l - cons_c[IBZ] * (1.0-k_sq) / (eta_l-k_bc);
    k_sq = SQR(kx_r) + SQR(ky_r) + SQR(kz_r);
    k_bc = kx_r*bbx + ky_r*cons_c[IBY] + kz_r*cons_c[IBZ];
    Real vx_cr = kx_r - bbx * (1.0-k_sq) / (eta_r-k_bc);
    Real vy_cr = ky_r - cons_c[IBY] * (1.0-k_sq) / (eta_r-k_bc);
    Real vz_cr = kz_r - cons_c[IBZ] * (1.0-k_sq) / (eta_r-k_bc);
    Real vx_c = 0.5 * (vx_cl + vx_cr);
    Real vy_c = 0.5 * (vy_cl + vy_cr);
    Real vz_c = 0.5 * (vz_cl + vz_cr);

    // Calculate remaining conserved quantities in c region
    Real v_bb_c = vx_c*bbx + vy_c*cons_c[IBY] + vz_c*cons_c[IBZ];
    if (vx_c >= 0.0)  // cL region
    {
      cons_c[IDN] = cons_al[IDN] * (lambda_al-vx_al) / (lambda_al-vx_c);     // (MUB 50)
      cons_c[IEN] = (lambda_al*cons_al[IEN] - cons_al[ivx] + ptot_true*vx_c
          - v_bb_c*bbx) / (lambda_al-vx_c);                                  // (MUB 51)
    }
    else  // cR region
    {
      cons_c[IDN] = cons_ar[IDN] * (lambda_ar-vx_ar) / (lambda_ar-vx_c);     // (MUB 50)
      cons_c[IEN] = (lambda_ar*cons_ar[IEN] - cons_ar[ivx] + ptot_true*vx_c
          - v_bb_c*bbx) / (lambda_ar-vx_c);                                  // (MUB 51)
    }
    cons_c[ivx] = (cons_c[IEN] + ptot_true) * vx_c - v_bb_c * bbx;          // (MUB 52)
    cons_c[ivy] = (cons_c[IEN] + ptot_true) * vy_c - v_bb_c * cons_c[IBY];  // (MUB 52)
    cons_c[ivz] = (cons_c[IEN] + ptot_true) * vz_c - v_bb_c * cons_c[IBZ];  // (MUB 52)

    // Calculate fluxes in c region (MUB 11)
    Real flux_c[NWAVE];
    for (int n = 0; n < NWAVE; ++n)
    {
      if (vx_c >= 0.0)  // cL region
        flux_c[n] = flux_al[n] + lambda_al * (cons_c[n] - cons_al[n]);
      else  // cR region
        flux_c[n] = flux_ar[n] + lambda_ar * (cons_c[n] - cons_ar[n]);
    }

    // Set conserved quantities in GR
    if (GENERAL_RELATIVITY)
      for (int n = 0; n < NWAVE; ++n)
      {
        if (lambda_l >= v_interface)  // L region
          cons_(n,i) = cons_l[n];
        else if (lambda_r <= v_interface)  // R region
          cons_(n,i) = cons_r[n];
        else if (lambda_al >= v_interface-vc_extension)  // aL region
          cons_(n,i) = cons_al[n];
        else if (lambda_ar <= v_interface+vc_extension)  // aR region
          cons_(n,i) = cons_ar[n];
        else  // c region
          cons_(n,i) = cons_c[n];
      }

    // Set fluxes
    for (int n = 0; n < NWAVE; ++n)
    {
      if (lambda_l >= v_interface)  // L region
        flux(n,i) = flux_l[n];
      else if (lambda_r <= v_interface)  // R region
        flux(n,i) = flux_r[n];
      else if (lambda_al >= v_interface-vc_extension)  // aL region
        flux(n,i) = flux_al[n];
      else if (lambda_ar <= v_interface+vc_extension)  // aR region
        flux(n,i) = flux_ar[n];
      else  // c region
        flux(n,i) = flux_c[n];
    }
  }

  // Transform fluxes to global coordinates if in GR
  if (GENERAL_RELATIVITY)
    switch (ivx)
    {
      case IVX:
        pmy_fluid->pmy_block->pcoord->FluxToGlobal1(k, j, il, iu, cons_, bb_normal_,
            flux);
        break;
      case IVY:
        pmy_fluid->pmy_block->pcoord->FluxToGlobal2(k, j, il, iu, cons_, bb_normal_,
            flux);
        break;
      case IVZ:
        pmy_fluid->pmy_block->pcoord->FluxToGlobal3(k, j, il, iu, cons_, bb_normal_,
            flux);
        break;
    }
  return;
}

// Function for finding total pressure given conserved quantities in flat spacetime
// Inputs:
//   cons: conserved state, excluding normal magnetic field
//   bbx: normal magnetic field
//   gamma_adi: ratio of specific heats
//   ivx: type of interface (IVX for x1, IVY for x2, IVZ for x3)
// Outputs:
//   returned value: total pressure
// Notes:
//   follows Mignone & McKinney 2007, MNRAS 378 1118 (MM)
//   follows hlld_sr.c in Athena 4.2 in using W and E rather than W' and E'
static Real ConsToPFlat(const Real cons[NWAVE], Real bbx, Real gamma_adi, int ivx)
{
  // Calculate cyclic permutations of indices
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  // Calculate variations on conserved quantities
  Real m_sq = SQR(cons[ivx]) + SQR(cons[ivy]) + SQR(cons[ivz]);
  Real b_sq = SQR(bbx) + SQR(cons[IBY]) + SQR(cons[IBZ]);
  Real m_dot_b = cons[ivx]*bbx + cons[ivy]*cons[IBY] + cons[ivz]*cons[IBZ];
  Real s_sq = SQR(m_dot_b);

  // Construct initial guess for enthalpy W (MM A26-A27)
  Real a1 = 4.0/3.0 * (b_sq - cons[IEN]);
  Real a0 = 1.0/3.0 * (m_sq + b_sq * (b_sq - 2.0*cons[IEN]));
  Real w_init = quadratic_root(a1, a0, true);

  // Apply Newton-Raphson method to find new W
  Real w_true = FindRootNR(w_init, cons[IDN], cons[IEN], m_sq, b_sq, s_sq,
      gamma_adi);

  // Calculate primitives from W
  Real v_sq = (m_sq + s_sq/SQR(w_true) * (2.0*w_true + b_sq))
      / SQR(w_true + b_sq);                                              // (MM A3)
  Real gamma_lor_sq = 1.0/(1.0-v_sq);
  Real gamma_lor = std::sqrt(gamma_lor_sq);
  Real chi = (1.0 - v_sq) * (w_true - gamma_lor * cons[IDN]);            // (cf. MM A11)
  Real pgas = (gamma_adi-1.0)/gamma_adi * chi;                           // (MM A17)
  Real vx = (cons[ivx] + m_dot_b/w_true * bbx) / (w_true + b_sq);         // (MM A10)
  Real vy = (cons[ivy] + m_dot_b/w_true * cons[IBY]) / (w_true + b_sq);  // (MM A10)
  Real vz = (cons[ivz] + m_dot_b/w_true * cons[IBZ]) / (w_true + b_sq);  // (MM A10)

  // Calculate total pressure
  Real v_bb = vx*bbx + vy*cons[IBY] + vz*cons[IBZ];
  Real bcov_sq = b_sq / gamma_lor_sq + SQR(v_bb);   // (MM 2)
  Real ptot = pgas + 0.5*bcov_sq;
  return ptot;
}

// Function whose value vanishes for correct enthalpy
// Inputs:
//   w_guess: guess for total enthalpy W
//   d: relativistic density D
//   e: total energy E
//   m_sq: square magnitude of momentum \vec{m}
//   b_sq: square magnitude of magnetic field \vec{B}
//   s_sq: (\vec{m} \cdot \vec{B})^2
//   gamma_adi: ratio of specific heats
// Outputs:
//   returned value: calculated minus given value of E
// Notes:
//   follows Mignone & McKinney 2007, MNRAS 378 1118 (MM)
//   implementation follows that of hlld_sr.c in Athena 4.2
static Real EResidual(Real w_guess, Real d, Real e, Real m_sq, Real b_sq, Real s_sq,
    Real gamma_adi)
{
  Real v_sq = (m_sq + s_sq/SQR(w_guess) * (2.0*w_guess + b_sq))
      / SQR(w_guess + b_sq);                                     // (MM A3)
  Real gamma_sq = 1.0/(1.0-v_sq);
  Real gamma_lorentz = std::sqrt(gamma_sq);
  Real chi = (1.0 - v_sq) * (w_guess - gamma_lorentz * d);       // (cf. MM A11)
  Real pgas = (gamma_adi-1.0)/gamma_adi * chi;                   // (MM A17)
  Real e_calc = w_guess - pgas + 0.5 * b_sq * (1.0+v_sq)
      - s_sq / (2.0*SQR(w_guess));                               // (MM A1)
  return e_calc - e;
}

// Derivative of EResidual()
// Inputs:
//   w_guess: guess for total enthalpy W
//   d: relativistic density D
//   m_sq: square magnitude of momentum \vec{m}
//   b_sq: square magnitude of magnetic field \vec{B}
//   s_sq: (\vec{m} \cdot \vec{B})^2
//   gamma_adi: ratio of specific heats
// Outputs:
//   returned value: derivative of calculated value of E
// Notes:
//   follows Mignone & McKinney 2007, MNRAS 378 1118 (MM)
//   implementation follows that of hlld_sr.c in Athena 4.2
static Real EResidualPrime(Real w_guess, Real d, Real m_sq, Real b_sq, Real s_sq,
    Real gamma_adi)
{
  Real v_sq = (m_sq + s_sq/SQR(w_guess) * (2.0*w_guess + b_sq))
      / SQR(w_guess + b_sq);                                            // (MM A3)
  Real gamma_sq = 1.0/(1.0-v_sq);
  Real gamma_lorentz = std::sqrt(gamma_sq);
  Real chi = (1.0 - v_sq) * (w_guess - gamma_lorentz * d);              // (cf. MM A11)
  Real w_cu = SQR(w_guess) * w_guess;
  Real w_b_cu = SQR(w_guess + b_sq) * (w_guess + b_sq);
  Real dv_sq_dw = -2.0 / (w_cu*w_b_cu)
      * (s_sq * (3.0*w_guess*(w_guess+b_sq) + SQR(b_sq)) + m_sq*w_cu);  // (MM A16)
  Real dchi_dw = 1.0 - v_sq
      - gamma_lorentz/2.0 * (d + 2.0*gamma_lorentz*chi) * dv_sq_dw;     // (cf. MM A14)
  Real drho_dw = -gamma_lorentz*d/2.0 * dv_sq_dw;                       // (MM A15)
  Real dpgas_dchi = (gamma_adi-1.0)/gamma_adi;                          // (MM A18)
  Real dpgas_drho = 0.0;                                                // (MM A18)
  Real dpgas_dw = dpgas_dchi * dchi_dw + dpgas_drho * drho_dw;
  return 1.0 - dpgas_dw + 0.5*b_sq*dv_sq_dw + s_sq/w_cu;
}

// Function whose value vanishes for correct total pressure
// Inputs:
//   p: guess for total pressure
//   bbx: normal magnetic field
//   lambda_l,lambda_r: L/R fast wavespeeds
//   r_l,r_r: L/R fast wave jump quantities
//   ivx: type of interface (IVX for x1, IVY for x2, IVZ for x3)
// Outputs:
//   returned value: residual
// Notes:
//   follows Mignone, Ugliano, & Bodo 2009, MNRAS 393 1141 (MUB)
//   follows Athena 4.2, hlld_sr.c, in variable choices and magic numbers
static Real PResidual(Real p, Real bbx, Real lambda_l, Real lambda_r,
    const Real r_l[NWAVE], const Real r_r[NWAVE], int ivx)
{
  // Parameters
  const Real delta_kx_aug = 1.0e-12;

  // Calculate cyclic permutations of indices
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  // Calculate v_aL and v_aR
  Real a = r_l[ivx] - lambda_l*r_l[IEN] + p*(1.0-SQR(lambda_l));       // (MUB 26)
  Real g = SQR(r_l[IBY]) + SQR(r_l[IBZ]);                              // (MUB 27)
  Real c = r_l[ivy]*r_l[IBY] + r_l[ivz]*r_l[IBZ];                      // (MUB 28)
  Real q = -a - g + SQR(bbx)*(1.0-SQR(lambda_l));                      // (MUB 29)
  Real x = bbx * (a*lambda_l*bbx+c) - (a+g) * (lambda_l*p+r_l[IEN]);   // (MUB 30)
  Real vx_al = (bbx * (a*bbx+lambda_l*c) - (a+g) * (p+r_l[ivx])) / x;  // (MUB 23)
  Real vy_al = (q*r_l[ivy]
      + r_l[IBY] * (c + bbx * (lambda_l*r_l[ivx]-r_l[IEN]))) / x;      // (MUB 24)
  Real vz_al = (q*r_l[ivz]
      + r_l[IBZ] * (c + bbx * (lambda_l*r_l[ivx]-r_l[IEN]))) / x;      // (MUB 24)
  a = r_r[ivx] - lambda_r*r_r[IEN] + p*(1.0-SQR(lambda_r));            // (MUB 26)
  g = SQR(r_r[IBY]) + SQR(r_r[IBZ]);                                   // (MUB 27)
  c = r_r[ivy]*r_r[IBY] + r_r[ivz]*r_r[IBZ];                           // (MUB 28)
  q = -a - g + SQR(bbx)*(1.0-SQR(lambda_r));                           // (MUB 29)
  x = bbx * (a*lambda_r*bbx+c) - (a+g) * (lambda_r*p+r_r[IEN]);        // (MUB 30)
  Real vx_ar = (bbx * (a*bbx+lambda_r*c) - (a+g) * (p+r_r[ivx])) / x;  // (MUB 23)
  Real vy_ar = (q*r_r[ivy]
      + r_r[IBY] * (c + bbx * (lambda_r*r_r[ivx]-r_r[IEN]))) / x;      // (MUB 24)
  Real vz_ar = (q*r_r[ivz]
      + r_r[IBZ] * (c + bbx * (lambda_r*r_r[ivx]-r_r[IEN]))) / x;      // (MUB 24)

  // Calculate B_aL and B_aR (MUB 21)
  Real by_al = (r_l[IBY] - bbx*vy_al) / (lambda_l - vx_al);
  Real bz_al = (r_l[IBZ] - bbx*vz_al) / (lambda_l - vx_al);
  Real by_ar = (r_r[IBY] - bbx*vy_ar) / (lambda_r - vx_ar);
  Real bz_ar = (r_r[IBZ] - bbx*vz_ar) / (lambda_r - vx_ar);

  // Calculate w_aL and w_aR (MUB 31)
  Real v_dot_rm = vx_al*r_l[ivx] + vy_al*r_l[ivy] + vz_al*r_l[ivz];
  Real w_al = p + (r_l[IEN] - v_dot_rm) / (lambda_l - vx_al);
  v_dot_rm = vx_ar*r_r[ivx] + vy_ar*r_r[ivy] + vz_ar*r_r[ivz];
  Real w_ar = p + (r_r[IEN] - v_dot_rm) / (lambda_r - vx_ar);

  // Calculate eta_L and eta_R (MUB 35)
  Real eta_l = -copysign(std::sqrt(w_al), bbx);
  Real eta_r = copysign(std::sqrt(w_ar), bbx);

  // Calculate K_L and K_R (MUB 43)
  Real denom = lambda_l*p + r_l[IEN] + bbx*eta_l;
  Real kx_l = (r_l[ivx] + p + lambda_l*bbx*eta_l) / denom;  // R_{B^x} = \lambda B^x
  Real ky_l = (r_l[ivy] + r_l[IBY]*eta_l) / denom;
  Real kz_l = (r_l[ivz] + r_l[IBZ]*eta_l) / denom;
  denom = lambda_r*p + r_r[IEN] + bbx*eta_r;
  Real kx_r = (r_r[ivx] + p + lambda_r*bbx*eta_r) / denom;  // R_{B^x} = \lambda B^x
  Real ky_r = (r_r[ivy] + r_r[IBY]*eta_r) / denom;
  Real kz_r = (r_r[ivz] + r_r[IBZ]*eta_r) / denom;

  // Rename Alfven wavespeeds for what they are
  Real &lambda_al = kx_l;
  Real &lambda_ar = kx_r;

  // Calculate B_c (MUB 45)
  Real delta_kx = kx_r - kx_l + delta_kx_aug;
  Real bbx_c_delta_kx = bbx * delta_kx;
  Real by_c_delta_kx = by_ar*(lambda_ar-vx_ar) - by_al*(lambda_al-vx_al)
      + bbx*(vy_ar-vy_al);
  Real bz_c_delta_kx = bz_ar*(lambda_ar-vx_ar) - bz_al*(lambda_al-vx_al)
      + bbx*(vz_ar-vz_al);

  // Calculate residual
  Real k_sq = SQR(kx_l) + SQR(ky_l) + SQR(kz_l);
  Real k_dot_bc_delta_kx = kx_l*bbx_c_delta_kx + ky_l*by_c_delta_kx
      + kz_l*bz_c_delta_kx;
  Real y_l = (1.0-k_sq) / (eta_l*delta_kx - k_dot_bc_delta_kx);      // (MUB 49)
  k_sq = SQR(kx_r) + SQR(ky_r) + SQR(kz_r);
  k_dot_bc_delta_kx = kx_r*bbx_c_delta_kx + ky_r*by_c_delta_kx
      + kz_r*bz_c_delta_kx;
  Real y_r = (1.0-k_sq) / (eta_r*delta_kx - k_dot_bc_delta_kx);      // (MUB 49)
  return delta_kx * (1.0 - bbx * (y_r-y_l));                         // (MUB 48)
}

// Newton-Raphson root finder
// Inputs:
//   w_init: initial guess for total enthalpy W
//   d: relativistic density D
//   e: total energy E
//   m_sq: square magnitude of momentum \vec{m}
//   b_sq: square magnitude of magnetic field \vec{B}
//   s_sq: (\vec{m} \cdot \vec{B})^2
//   gamma_adi: ratio of specific heats
// Outputs:
//   returned value: total enthalpy W
// Notes:
//   returns NAN in event of failure
//   forces W to be positive
static Real FindRootNR(Real w_init, Real d, Real e, Real m_sq, Real b_sq, Real s_sq,
    Real gamma_adi)
{
  // Parameters
  const int max_iterations = 100;      // maximum number of iterations
  const Real tol_w = 1.0e-8 * w_init;  // absolute tolerance in W
  const Real tol_res = 1.0e-15;        // absolute tolerance in residual

  // Check if root has already been found
  Real new_res = EResidual(w_init, d, e, m_sq, b_sq, s_sq, gamma_adi);
  if (std::abs(new_res) < tol_res)
    return w_init;

  // Iterate to find root
  Real new_w = w_init;
  for (int i = 0; i < max_iterations; ++i)
  {
    // Prepare needed values
    Real old_w = new_w;
    Real old_res = new_res;
    Real derivative = EResidualPrime(old_w, d, m_sq, b_sq, s_sq, gamma_adi);
    Real delta = -old_res / derivative;

    // Check that update makes sense
    if (!std::isfinite(delta))
      return NAN;

    // Reduce step if root goes out of bounds
    int j;
    for (j = i; j < max_iterations; ++j)
    {
      new_w = old_w + delta;
      if (new_w > 0.0)
        break;
      else
        delta /= 2.0;
    }
    i = j;

    // Reduce step if new value is worse than old
    for (j = i; j < max_iterations; ++j)
    {
      new_res = EResidual(new_w, d, e, m_sq, b_sq, s_sq, gamma_adi);
      if (std::abs(new_res) < std::abs(old_res))
        break;
      else
      {
        delta /= 2.0;
        new_w = old_w + delta;
      }
    }
    i = j;

    // Check if root found
    if (std::abs(new_res) < tol_res || std::abs(delta) < tol_w)
      return new_w;
  }

  // Indicate failure to converge
  return NAN;
}

// Secant root finder
// Inputs:
//   ptot_init: initial guess for total pressure inside fast waves
//   bbx: normal magnetic field
//   lambda_l,lambda_r: L/R fast wavespeeds
//   r_l,r_r: L/R fast wave jump quantities
//   ivx: type of interface (IVX for x1, IVY for x2, IVZ for x3)
// Outputs:
//   returned value: total pressure
// Notes:
//   returns NAN in event of failure
static Real FindRootSecant(Real ptot_init, Real bbx, Real lambda_l, Real lambda_r,
    const Real r_l[NWAVE], const Real r_r[NWAVE], int ivx)
{
  // Parameters
  const int max_iterations = 20;             // maximum number of iterations
  const Real ptot_max = 100.0 * ptot_init;   // upper bound on root
  const Real tol_ptot = 1.0e-8 * ptot_init;  // absolute tolerance in total pressure
  const Real tol_res = 1.0e-12;              // absolute tolerance in residual
  const Real tol_ptot_factor = 1.0e3;        // scale of tol used to perturb guess

  // Check if initial guess is good enough
  Real current_p = ptot_init;
  Real current_res = PResidual(current_p, bbx, lambda_l, lambda_r, r_l, r_r, ivx);
  if (std::abs(current_res) < tol_res)
    return current_p;

  // Set other point needed to find secant
  Real new_p = current_p * (1.0 + tol_ptot_factor * tol_ptot);
  if (new_p > ptot_max)  // new_p already guaranteed to be positive if ptot_init > 0
    return NAN;

  // Check if second guess is good enough
  Real new_res = PResidual(new_p, bbx, lambda_l, lambda_r, r_l, r_r, ivx);
  if (std::abs(new_res) < tol_res)
    return new_p;

  // Iterate to find root
  for (int i = 0; i < max_iterations; ++i)
  {
    // Shift values
    Real old_p = current_p;
    Real old_res = current_res;
    current_p = new_p;
    current_res = new_res;

    // Find new abscissa and check if change is small enough
    new_p = (current_res*old_p - old_res*current_p) / (current_res-old_res);
    if (new_p <= 0.0 || new_p > ptot_max)
      return NAN;
    if (std::abs(new_p-current_p) < tol_ptot)
      return current_p;

    // Find new ordinate and check if root found
    new_res = PResidual(new_p, bbx, lambda_l, lambda_r,
        r_l, r_r, ivx);
    if (std::abs(new_res) < tol_res)
      return new_p;
  }

  // Treat non-convergence as error
  return NAN;
}

// Function for finding root of monic quadratic equation
// Inputs:
//   a1: linear coefficient
//   a0: constant coefficient
//   greater_root: flag indicating that larger root is to be returned
//     "larger" does not mean absolute value
// Outputs:
//   returned value: desired root
// Notes:
//   same function as in adiabatic_mhd_sr.cpp, adiabatic_mhd_gr.cpp, and
//       linear_wave_rel.cpp
//   solves x^2 + a_1 x + a_0 = 0 for x
//   returns abscissa of vertex if there are no real roots
//   follows advice in Numerical Recipes, 3rd ed. (5.6) for avoiding large cancellations
static Real quadratic_root(Real a1, Real a0, bool greater_root)
{
  if (a1*a1 < 4.0*a0)  // no real roots
    return -a1/2.0;
  if (greater_root)
  {
    if (a1 >= 0.0)
      return -2.0*a0 / (a1 + std::sqrt(a1*a1 - 4.0*a0));
    else
      return (-a1 + std::sqrt(a1*a1 - 4.0*a0)) / 2.0;
  }
  else
  {
    if (a1 >= 0.0)
      return (-a1 - std::sqrt(a1*a1 - 4.0*a0)) / 2.0;
    else
      return -2.0*a0 / (a1 - std::sqrt(a1*a1 - 4.0*a0));
  }
}
