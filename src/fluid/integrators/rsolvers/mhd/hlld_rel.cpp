// HLLD Riemann solver for special relativistic MHD

// Primary header
#include "../../fluid_integrator.hpp"

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()

// Athena headers
#include "../../../../athena.hpp"                   // enums, macros, Real
#include "../../../../athena_arrays.hpp"            // AthenaArray
#include "../../../eos/eos.hpp"                     // GetGamma()
#include "../../../fluid.hpp"                       // Fluid
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
static Real ConsToPFlat(const Real cons[NWAVE], Real bx, Real gamma_adi_red, int ivx);
static Real EResidual(Real w_guess, Real d, Real e, Real m_sq, Real b_sq, Real s_sq,
    Real gamma_prime);
static Real EResidualPrime(Real w_guess, Real d, Real m_sq, Real b_sq, Real s_sq,
    Real gamma_prime);
static Real PResidual(Real p, Real bx, Real l_l, Real l_r,
    const Real r_l[NWAVE], const Real r_r[NWAVE], int ivx);
static void CalculateAState(
    Real p, Real lambda, Real bx, const Real r[NWAVE], int ivx, bool left_flag,
    Real *pvx_a, Real *pvy_a, Real *pvz_a, Real *pkx, Real *pky, Real *pkz, Real *peta,
    Real cons_a[NWAVE], Real flux_a[NWAVE]);
static void CalculateCState(Real p, Real bx,
    Real l_al, Real l_ar,
    Real vx_al, Real vy_al, Real vz_al, Real vx_ar, Real vy_ar, Real vz_ar,
    Real kx_l, Real ky_l, Real kz_l, Real kx_r, Real ky_r, Real kz_r,
    Real eta_l, Real eta_r,
    const Real c_al[NWAVE], const Real c_ar[NWAVE], int ivx,
    Real *pvx_c, Real c_c[NWAVE]);
static Real FindRootNR(Real w_initial, Real d, Real e, Real m_sq, Real b_sq, Real s_sq,
    Real gamma_prime);
static Real FindRootSecant(
    Real ptot_initial, Real bx, Real lambda_left, Real lambda_right,
    const Real r_left[NWAVE], const Real r_right[NWAVE], int ivx);
static Real quadratic_root(Real a1, Real a0, bool greater_root);

// Riemann solver
// Inputs:
//   il, iu: lower and upper indices for interfaces
//   pprim_left, pprim_right: pointers to left and right primitive states
// Outputs:
//   pflux: pointer to fluxes
// Notes:
//   implements HLLD solver from Mignone, Ugliano, & Bodo 2009, MNRAS 393 1141 (MUB)
//   references Mignone & Bodo 2006, MNRAS 368 1040 (MB)
void FluidIntegrator::RiemannSolver(const int k,const int j, const int il, const int iu,
  const int ivx, const AthenaArray<Real> &b, AthenaArray<Real> &prim_left,
  AthenaArray<Real> &prim_right, AthenaArray<Real> &flux)
{
  // Parameters
  const Real p_transition = 0.1;     // value delineating intial pressure regimes
  const Real vc_extension = 1.0e-6;  // use contact region if Alfven speeds smaller

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
  const Real gamma_adi_red = gamma_adi / (gamma_adi - 1.0);

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
    const Real bx = b_normal_(i);

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
        &lambda_left_plus, &lambda_left_minus);                          // (MB 56)
    Real lambda_right_plus, lambda_right_minus;
    pmy_fluid->pf_eos->FastMagnetosonicSpeedsRel(
        rho_right, pgas_right,
        vx_right, vy_right, vz_right,
        ut_right, ux_right, uy_right, uz_right,
        bx, by_right, bz_right,
        bcovt_right, bcovx_right, bcovy_right, bcovz_right,
        &lambda_right_plus, &lambda_right_minus);                        // (MB 56)
    Real lambda_left = std::min(lambda_left_minus, lambda_right_minus);  // (MB 55)
    Real lambda_right = std::max(lambda_left_plus, lambda_right_plus);   // (MB 55)

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

    // Calculate L/R state conserved quantities
    Real cons_left[NWAVE], cons_right[NWAVE];
    PrimToConsFlat(gamma_adi_red, rho_left, pgas_left,
        ut_left, ux_left, uy_left, uz_left,
        bcovt_left, bcovx_left, bcovy_left, bcovz_left,
        ivx, ivy, ivz, cons_left);
    PrimToConsFlat(gamma_adi_red, rho_right, pgas_right,
        ut_right, ux_right, uy_right, uz_right,
        bcovt_right, bcovx_right, bcovy_right, bcovz_right,
        ivx, ivy, ivz, cons_right);

    // Calculate fast wave jump quantities and HLL state and fluxes
    Real r_left[NWAVE], r_right[NWAVE];
    Real cons_hll[NWAVE], flux_hll[NWAVE];
    for (int n = 0; n < NWAVE; ++n)
    {
      r_left[n] = lambda_left * cons_left[n] - flux_left[n];              // (MUB 12)
      r_right[n] = lambda_right * cons_right[n] - flux_right[n];          // (MUB 12)
      cons_hll[n] = (r_right[n]-r_left[n]) / (lambda_right-lambda_left);  // (MB 29)
      flux_hll[n] = (lambda_left*r_right[n] - lambda_right*r_left[n])
          / (lambda_right-lambda_left);                                   // (MB 31)
    }

    // Calculate initial guess for total pressure (MUB 53)
    Real ptot_initial;
    Real ptot_hll = ConsToPFlat(cons_hll, bx, gamma_adi_red, ivx);
    if (SQR(bx)/ptot_hll < p_transition)  // weak magnetic field
    {
      Real a1 = cons_hll[IEN] - flux_hll[ivx];
      Real a0 = cons_hll[ivx]*flux_hll[IEN] - flux_hll[ivx]*cons_hll[IEN];
      ptot_initial = quadratic_root(a1, a0, true);                          // (MUB 55)
    }
    else  // strong magnetic field
      ptot_initial = ptot_hll;

    // Apply secant method to find total pressure
    Real ptot_true = FindRootSecant(ptot_initial, bx, lambda_left, lambda_right,
        r_left, r_right, ivx);

    // Set fluxes if in aL state
    Real vx_al, vy_al, vz_al;
    Real kx_l, ky_l, kz_l, eta_l;
    Real cons_al[NWAVE], flux_al[NWAVE];
    CalculateAState(ptot_true, lambda_left, bx, r_left, ivx, true,
        &vx_al, &vy_al, &vz_al, &kx_l, &ky_l, &kz_l, &eta_l,
        cons_al, flux_al);
    Real lambda_al = kx_l;
    if (lambda_al >= -vc_extension)
    {
      for (int n = 0; n < NWAVE; ++n)
        flux(n,i) = flux_al[n];
      continue;
    }

    // Set fluxes if in aL state
    Real vx_ar, vy_ar, vz_ar;
    Real kx_r, ky_r, kz_r, eta_r;
    Real cons_ar[NWAVE], flux_ar[NWAVE];
    CalculateAState(ptot_true, lambda_right, bx, r_right, ivx, false,
        &vx_ar, &vy_ar, &vz_ar, &kx_r, &ky_r, &kz_r, &eta_r,
        cons_ar, flux_ar);
    Real lambda_ar = kx_r;
    if (lambda_ar <= vc_extension)
    {
      for (int n = 0; n < NWAVE; ++n)
        flux(n,i) = flux_ar[n];
      continue;
    }

    // Set fluxes in c state
    Real vx_c;
    Real cons_c[NWAVE];
    CalculateCState(ptot_true, bx,
        lambda_al, lambda_ar,
        vx_al, vy_al, vz_al, vx_ar, vy_ar, vz_ar,
        kx_l, ky_l, kz_l, kx_r, ky_r, kz_r,
        eta_l, eta_r,
        cons_al, cons_ar, ivx,
        &vx_c, cons_c);
    if (vx_c >= 0.0)
      for (int n = 0; n < NWAVE; ++n)
        flux(n,i) = flux_al[n] + lambda_al * (cons_c[n] - cons_al[n]);  // (MUB 11)
    else
      for (int n = 0; n < NWAVE; ++n)
        flux(n,i) = flux_ar[n] + lambda_ar * (cons_c[n] - cons_ar[n]);  // (MUB 11)
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
//   same function as in hlle_mhd_rel.cpp
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
//   same function as in hlle_mhd_rel.cpp
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

// Function for finding total pressure given conserved quantities in flat spacetime
// Notes:
//   follows Mignone & McKinney 2007, MNRAS 378 1118 (MM)
//   follows hlld_sr.c in Athena 4.2 in using W and E rather than W' and E'
static Real ConsToPFlat(const Real cons[NWAVE], Real bx, Real gamma_adi_red, int ivx)
{
  // Calculate cyclic permutations of indices
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  // Calculate variations on conserved quantities
  Real m_sq = SQR(cons[ivx]) + SQR(cons[ivy]) + SQR(cons[ivz]);
  Real b_sq = SQR(bx) + SQR(cons[IBY]) + SQR(cons[IBZ]);
  Real m_dot_b = cons[ivx]*bx + cons[ivy]*cons[IBY] + cons[ivz]*cons[IVZ];
  Real s_sq = SQR(m_dot_b);

  // Construct initial guess for enthalpy W (cf. MM A26-A27)
  Real a1 = 4.0/3.0 * (b_sq - cons[IEN]);
  Real a0 = 1.0/3.0 * (m_sq + b_sq * (b_sq - 2.0*cons[IEN]));
  Real w_initial = quadratic_root(a1, a0, true);

  // Apply Newton-Raphson method to find new W
  Real w_true = FindRootNR(w_initial, cons[IDN], cons[IEN], m_sq, b_sq, s_sq,
      gamma_adi_red);

  // Calculate primitives from W
  Real v_sq = (m_sq + s_sq/SQR(w_true) * (2.0*w_true + b_sq))
      / SQR(w_true + b_sq);                                              // (cf. MM A3)
  Real gamma_lor_sq = 1.0/(1.0-v_sq);
  Real gamma_lor = std::sqrt(gamma_lor_sq);
  Real chi = (1.0 - v_sq) * (w_true - gamma_lor * cons[IDN]);            // (cf. MM A11)
  Real pgas = chi/gamma_adi_red;                                         // (MM A17)
  Real vx = (cons[ivx] + m_dot_b/w_true * bx) / (w_true + b_sq);         // (MM A10)
  Real vy = (cons[ivy] + m_dot_b/w_true * cons[IBY]) / (w_true + b_sq);  // (MM A10)
  Real vz = (cons[ivz] + m_dot_b/w_true * cons[IBZ]) / (w_true + b_sq);  // (MM A10)

  // Calculate total pressure
  Real v_dot_b = vx*bx + vy*cons[IBY] + vz*cons[IBZ];
  Real bcov_sq = b_sq / gamma_lor_sq + SQR(v_dot_b);   // (MM 2)
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
//   gamma_prime: reduced adiabatic gas constant Gamma' = Gamma/(Gamma-1)
// Outputs:
//   returned value: calculated minus given value of E
// Notes:
//   follows Mignone & McKinney 2007, MNRAS 378 1118 (MM)
//   implementation follows that of hlld_sr.c in Athena 4.2
//   same function as in hlld_mhd_rel.cpp
static Real EResidual(Real w_guess, Real d, Real e, Real m_sq, Real b_sq, Real s_sq,
    Real gamma_prime)
{
  Real v_sq = (m_sq + s_sq/SQR(w_guess) * (2.0*w_guess + b_sq))
      / SQR(w_guess + b_sq);                                     // (cf. MM A3)
  Real gamma_sq = 1.0/(1.0-v_sq);
  Real gamma_lorentz = std::sqrt(gamma_sq);
  Real chi = (1.0 - v_sq) * (w_guess - gamma_lorentz * d);       // (cf. MM A11)
  Real pgas = chi/gamma_prime;                                   // (MM A17)
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
//   gamma_prime: reduced adiabatic gas constant Gamma' = Gamma/(Gamma-1)
// Outputs:
//   returned value: derivative of calculated value of E
// Notes:
//   follows Mignone & McKinney 2007, MNRAS 378 1118 (MM)
//   implementation follows that of hlld_sr.c in Athena 4.2
//   same function as in hlld_mhd_rel.cpp
static Real EResidualPrime(Real w_guess, Real d, Real m_sq, Real b_sq, Real s_sq,
    Real gamma_prime)
{
  Real v_sq = (m_sq + s_sq/SQR(w_guess) * (2.0*w_guess + b_sq))
      / SQR(w_guess + b_sq);                                            // (cf. MM A3)
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
  Real dpgas_dchi = 1.0/gamma_prime;                                    // (MM A18)
  Real dpgas_drho = 0.0;                                                // (MM A18)
  Real dpgas_dw = dpgas_dchi * dchi_dw + dpgas_drho * drho_dw;
  return 1.0 - dpgas_dw + 0.5*b_sq*dv_sq_dw + s_sq/w_cu;
}

// Function whose value vanishes for correct total pressure
// Notes:
//   follows Mignone, Ugliano, & Bodo 2009, MNRAS 393 1141 (MUB)
static Real PResidual(Real p, Real bx, Real l_l, Real l_r,
    const Real r_l[NWAVE], const Real r_r[NWAVE], int ivx)
{
  // Calculate cyclic permutations of indices
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  // Calculate v_aL and v_aR
  Real a = r_l[ivx] - l_l*r_l[IEN] + p*(1.0-SQR(l_l));          // (MUB 26)
  Real g = SQR(r_l[IBY]) + SQR(r_l[IBZ]);                       // (MUB 27)
  Real c = r_l[ivy]*r_l[IBY] + r_l[ivz]*r_l[IBZ];               // (MUB 28)
  Real q = -a - g + SQR(bx)*(1.0-SQR(l_l));                     // (MUB 29)
  Real x = bx * (a*l_l*bx+c) - (a+g) * (l_l*p+r_l[IEN]);        // (MUB 30)
  Real vx_al = (bx * (a*bx+l_l*c) - (a+g) * (p+r_l[ivx])) / x;  // (MUB 23)
  Real vy_al = (q*r_l[ivy]
      + r_l[IBY] * (c + bx * (l_l*r_l[ivx]-r_l[IEN]))) / x;     // (MUB 24)
  Real vz_al = (q*r_l[ivz]
      + r_l[IBZ] * (c + bx * (l_l*r_l[ivx]-r_l[IEN]))) / x;     // (MUB 24)
  a = r_r[ivx] - l_r*r_r[IEN] + p*(1.0-SQR(l_r));               // (MUB 26)
  g = SQR(r_r[IBY]) + SQR(r_r[IBZ]);                            // (MUB 27)
  c = r_r[ivy]*r_r[IBY] + r_r[ivz]*r_r[IBZ];                    // (MUB 28)
  q = -a - g + SQR(bx)*(1.0-SQR(l_r));                          // (MUB 29)
  x = bx * (a*l_r*bx+c) - (a+g) * (l_r*p+r_r[IEN]);             // (MUB 30)
  Real vx_ar = (bx * (a*bx+l_r*c) - (a+g) * (p+r_r[ivx])) / x;  // (MUB 23)
  Real vy_ar = (q*r_r[ivy]
      + r_r[IBY] * (c + bx * (l_r*r_r[ivx]-r_r[IEN]))) / x;     // (MUB 24)
  Real vz_ar = (q*r_r[ivz]
      + r_r[IBZ] * (c + bx * (l_r*r_r[ivx]-r_r[IEN]))) / x;     // (MUB 24)

  // Calculate B_aL and B_aR (MUB 21)
  Real by_al = (r_l[IBY] - bx*vy_al) / (l_l - vx_al);
  Real bz_al = (r_l[IBZ] - bx*vz_al) / (l_l - vx_al);
  Real by_ar = (r_r[IBY] - bx*vy_ar) / (l_r - vx_ar);
  Real bz_ar = (r_r[IBZ] - bx*vz_ar) / (l_r - vx_ar);

  // Calculate w_aL and w_aR (MUB 31)
  Real v_dot_rm = vx_al*r_l[ivx] + vy_al*r_l[ivy] + vz_al*r_l[ivz];
  Real w_al = p + (r_l[IEN] - v_dot_rm) / (l_l - vx_al);
  v_dot_rm = vx_ar*r_r[ivx] + vy_ar*r_r[ivy] + vz_ar*r_r[ivz];
  Real w_ar = p + (r_r[IEN] - v_dot_rm) / (l_r - vx_ar);

  // Calculate eta_L and eta_R (MUB 35)
  Real eta_l = -copysign(std::sqrt(w_al), bx);
  Real eta_r = copysign(std::sqrt(w_ar), bx);

  // Calculate K_L and K_R (MUB 43)
  Real denominator = l_l*p + r_l[IEN] + bx*eta_l;
  Real kx_l = (r_l[ivx] + p + l_l*bx*eta_l) / denominator;  // R_{B^x} = \lambda B^x
  Real ky_l = (r_l[ivy] + r_l[IBY]*eta_l) / denominator;
  Real kz_l = (r_l[ivz] + r_l[IBZ]*eta_l) / denominator;
  denominator = l_r*p + r_r[IEN] + bx*eta_r;
  Real kx_r = (r_r[ivx] + p + l_r*bx*eta_r) / denominator;  // R_{B^x} = \lambda B^x
  Real ky_r = (r_r[ivy] + r_r[IBY]*eta_r) / denominator;
  Real kz_r = (r_r[ivz] + r_r[IBZ]*eta_r) / denominator;

  // Rename Alfven wavespeeds for what they are
  Real &l_al = kx_l;
  Real &l_ar = kx_r;

  // Calculate B_c (MUB 45)
  denominator = l_ar - l_al;
  Real numerator_al = by_al * (l_al-vx_al) + bx*vy_al;
  Real numerator_ar = by_ar * (l_ar-vx_ar) + bx*vy_ar;
  Real by_c = (numerator_ar - numerator_al) / denominator;
  numerator_al = bz_al * (l_al-vx_al) + bx*vz_al;
  numerator_ar = bz_ar * (l_ar-vx_ar) + bx*vz_ar;
  Real bz_c = (numerator_ar - numerator_al) / denominator;

  // Calculate residual
  Real delta_kx = kx_r - kx_l;
  Real k_sq = SQR(kx_l) + SQR(ky_l) + SQR(kz_l);
  Real k_dot_bc = kx_l*bx + ky_l*by_c + kz_l*bz_c;
  Real y_l = (1.0-k_sq) / (delta_kx * (eta_l-k_dot_bc));  // (cf. MUB 49)
  k_sq = SQR(kx_r) + SQR(ky_r) + SQR(kz_r);
  k_dot_bc = kx_r*bx + ky_r*by_c + kz_r*bz_c;
  Real y_r = (1.0-k_sq) / (delta_kx * (eta_r-k_dot_bc));  // (cf. MUB 49)
  return delta_kx * (1.0 - bx * (y_r-y_l));               // (MUB 48)
}

// Function for calculating variables between Alfven and fast waves
// Notes:
//   follows Mignone, Ugliano, & Bodo 2009, MNRAS 393 1141 (MUB)
static void CalculateAState(
    Real p, Real lambda, Real bx, const Real r[NWAVE], int ivx, bool left_flag,
    Real *pvx_a, Real *pvy_a, Real *pvz_a, Real *pkx, Real *pky, Real *pkz, Real *peta,
    Real cons_a[NWAVE], Real flux_a[NWAVE])
{
  // Calculate cyclic permutations of indices
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  // Set v_a
  Real a = r[ivx] - lambda*r[IEN] + p*(1.0-SQR(lambda));                 // (MUB 26)
  Real g = SQR(r[IBY]) + SQR(r[IBZ]);                                    // (MUB 27)
  Real c = r[ivy]*r[IBY] + r[ivz]*r[IBZ];                                // (MUB 28)
  Real q = -a - g + SQR(bx)*(1.0-SQR(lambda));                           // (MUB 29)
  Real x = bx * (a*lambda*bx+c) - (a+g) * (lambda*p+r[IEN]);             // (MUB 30)
  *pvx_a = (bx * (a*bx+lambda*c) - (a+g) * (p+r[ivx])) / x;              // (MUB 23)
  *pvy_a = (q*r[ivy] + r[IBY] * (c + bx * (lambda*r[ivx]-r[IEN]))) / x;  // (MUB 24)
  *pvz_a = (q*r[ivz] + r[IBZ] * (c + bx * (lambda*r[ivx]-r[IEN]))) / x;  // (MUB 24)

  // Set conserved quantities
  Real denominator = lambda - (*pvx_a);
  cons_a[IBY] = (r[IBY] - bx*(*pvy_a)) / denominator;                        // (MUB 21)
  cons_a[IBZ] = (r[IBZ] - bx*(*pvz_a)) / denominator;                        // (MUB 21)
  Real v_dot_b = (*pvx_a)*bx + (*pvy_a)*cons_a[IBY] + (*pvz_a)*cons_a[IBZ];
  cons_a[IDN] = r[IDN] / denominator;                                        // (MUB 32)
  cons_a[IEN] = (r[IEN] + p*(*pvx_a) - v_dot_b*bx) / denominator;            // (MUB 33)
  cons_a[ivx] = (cons_a[IEN] + p) * (*pvx_a) - v_dot_b * bx;                 // (MUB 34)
  cons_a[ivy] = (cons_a[IEN] + p) * (*pvy_a) - v_dot_b * cons_a[IBY];        // (MUB 34)
  cons_a[ivz] = (cons_a[IEN] + p) * (*pvz_a) - v_dot_b * cons_a[IBZ];        // (MUB 34)

  // Set fluxes (cf. MUB 11,12)
  for (int n = 0; n < NWAVE; ++n)
    flux_a[n] = lambda * cons_a[n] - r[n];

  // Set eta and K
  Real v_dot_rm = (*pvx_a)*r[ivx] + (*pvy_a)*r[ivy]
    + (*pvz_a)*r[ivz];
  Real w = p + (r[IEN]-v_dot_rm) / denominator;           // (MUB 31)
  *peta = (left_flag ? -1.0 : 1.0)
      * copysign(std::sqrt(w), bx);                       // (MUB 35)
  denominator = lambda*p + r[IEN] + bx*(*peta);
  *pkx = (r[ivx] + p + lambda*bx*(*peta)) / denominator;  // R_{B^x} = \lambda B^x
  *pky = (r[ivy] + r[IBY]*(*peta)) / denominator;
  *pkz = (r[ivz] + r[IBZ]*(*peta)) / denominator;
  return;
}

// Function for calculating variables between Alfven waves
// Notes:
//   follows Mignone, Ugliano, & Bodo 2009, MNRAS 393 1141 (MUB)
static void CalculateCState(Real p, Real bx,
    Real l_al, Real l_ar,
    Real vx_al, Real vy_al, Real vz_al, Real vx_ar, Real vy_ar, Real vz_ar,
    Real kx_l, Real ky_l, Real kz_l, Real kx_r, Real ky_r, Real kz_r,
    Real eta_l, Real eta_r,
    const Real c_al[NWAVE], const Real c_ar[NWAVE], int ivx,
    Real *pvx_c, Real c_c[NWAVE])
{
  // Calculate cyclic permutations of indices
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  // Set B_c (MUB 45)
  Real denominator = l_ar - l_al;
  Real numerator_al = c_al[IBY] * (l_al-vx_al) + bx*vy_al;
  Real numerator_ar = c_ar[IBY] * (l_ar-vx_ar) + bx*vy_ar;
  c_c[IBY] = (numerator_ar - numerator_al) / denominator;
  numerator_al = c_al[IBZ] * (l_al-vx_al) + bx*vz_al;
  numerator_ar = c_ar[IBZ] * (l_ar-vx_ar) + bx*vz_ar;
  c_c[IBZ] = (numerator_ar - numerator_al) / denominator;

  // Calculate v_c (MUB 47), averaging left and right values in case of disagreement
  Real k_sq = SQR(kx_l) + SQR(ky_l) + SQR(kz_l);
  Real k_dot_bc = kx_l*bx + ky_l*c_c[IBY] + kz_l*c_c[IBZ];
  Real vx_cl = kx_l - bx * (1.0-k_sq) / (eta_l-k_dot_bc);
  Real vy_cl = ky_l - c_c[IBY] * (1.0-k_sq) / (eta_l-k_dot_bc);
  Real vz_cl = kz_l - c_c[IBZ] * (1.0-k_sq) / (eta_l-k_dot_bc);
  k_sq = SQR(kx_r) + SQR(ky_r) + SQR(kz_r);
  k_dot_bc = kx_r*bx + ky_r*c_c[IBY] + kz_r*c_c[IBZ];
  Real vx_cr = kx_r - bx * (1.0-k_sq) / (eta_r-k_dot_bc);
  Real vy_cr = ky_r - c_c[IBY] * (1.0-k_sq) / (eta_r-k_dot_bc);
  Real vz_cr = kz_r - c_c[IBZ] * (1.0-k_sq) / (eta_r-k_dot_bc);
  *pvx_c = 0.5 * (vx_cl + vx_cr);
  Real vy_c = 0.5 * (vy_cl + vy_cr);
  Real vz_c = 0.5 * (vz_cl + vz_cr);

  // Set remaining conserved quantities
  Real vc_dot_bc = (*pvx_c)*bx + vy_c*c_c[IBY] + vz_c*c_c[IVZ];
  if (*pvx_c >= 0.0)
  {
    c_c[IDN] = c_al[IDN] * (l_al-vx_al) / (l_al-(*pvx_c));               // (MUB 50)
    c_c[IEN] = (l_al*c_al[IEN] - c_al[ivx] + p*(*pvx_c) - vc_dot_bc*bx)
        / (l_al-(*pvx_c));                                               // (MUB 51)
  }
  else
  {
    c_c[IDN] = c_ar[IDN] * (l_ar-vx_ar) / (l_ar-(*pvx_c));               // (MUB 50)
    c_c[IEN] = (l_ar*c_ar[IEN] - c_ar[ivx] + p*(*pvx_c) - vc_dot_bc*bx)
        / (l_ar-(*pvx_c));                                               // (MUB 51)
  }
  c_c[ivx] = (c_c[IEN] + p) * (*pvx_c) - vc_dot_bc * bx;    // (MUB 52)
  c_c[ivy] = (c_c[IEN] + p) * vy_c - vc_dot_bc * c_c[IBY];  // (MUB 52)
  c_c[ivz] = (c_c[IEN] + p) * vz_c - vc_dot_bc * c_c[IBZ];  // (MUB 52)
  return;
}

// Newton-Raphson root finder
// Inputs:
//   w_initial: initial guess for total enthalpy W
//   d: relativistic density D
//   e: total energy E
//   m_sq: square magnitude of momentum \vec{m}
//   b_sq: square magnitude of magnetic field \vec{B}
//   s_sq: (\vec{m} \cdot \vec{B})^2
//   gamma_prime: reduced adiabatic gas constant Gamma' = Gamma/(Gamma-1)
// Outputs:
//   returned value: total enthalpy W
// Notes:
//   returns NAN in event of failure
//   forces W to be positive
//   same function as in hlld_mhd_rel.cpp
static Real FindRootNR(Real w_initial, Real d, Real e, Real m_sq, Real b_sq, Real s_sq,
    Real gamma_prime)
{
  // Parameters
  const int max_iterations = 100;         // maximum number of iterations
  const Real tol_w = 1.0e-8 * w_initial;  // absolute tolerance in W
  const Real tol_res = 1.0e-15;           // absolute tolerance in residual

  // Check if root has already been found
  Real new_res = EResidual(w_initial, d, e, m_sq, b_sq, s_sq, gamma_prime);
  if (std::abs(new_res) < tol_res)
    return w_initial;

  // Iterate to find root
  Real new_w = w_initial;
  for (int i = 0; i < max_iterations; i++)
  {
    // Prepare needed values
    Real old_w = new_w;
    Real old_res = new_res;
    Real derivative = EResidualPrime(old_w, d, m_sq, b_sq, s_sq, gamma_prime);
    Real delta = -old_res / derivative;

    // Check that update makes sense
    if (!std::isfinite(delta))
      return NAN;

    // Reduce step if root goes out of bounds
    int j;
    for (j = i; j < max_iterations; j++)
    {
      new_w = old_w + delta;
      if (new_w > 0.0)
        break;
      else
        delta /= 2.0;
    }
    i = j;

    // Reduce step if new value is worse than old
    for (j = i; j < max_iterations; j++)
    {
      new_res = EResidual(new_w, d, e, m_sq, b_sq, s_sq, gamma_prime);
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
static Real FindRootSecant(
    Real ptot_initial, Real bx, Real lambda_left, Real lambda_right,
    const Real r_left[NWAVE], const Real r_right[NWAVE], int ivx)
{
  // Parameters
  const int max_iterations = 20;                // maximum number of iterations
  const Real ptot_max = 100.0 * ptot_initial;   // upper bound on root
  const Real tol_ptot = 1.0e-8 * ptot_initial;  // absolute tolerance in total pressure
  const Real tol_res = 1.0e-12;                 // absolute tolerance in residual
  const Real tol_ptot_factor = 1.0e3;           // scale of tol used to perturb guess

  // Check if initial guess is good enough
  Real current_p = ptot_initial;
  Real current_res = PResidual(current_p, bx, lambda_left, lambda_right,
      r_left, r_right, ivx);
  if (std::abs(current_res) < tol_res)
    return current_p;

  // Set other point needed to find secant
  Real new_p = current_p * (1.0 + tol_ptot_factor * tol_ptot);
  if (new_p > ptot_max)  // new_p already guaranteed to be positive if ptot_initial > 0
    return NAN;

  // Check if second guess is good enough
  Real new_res = PResidual(new_p, bx, lambda_left, lambda_right,
      r_left, r_right, ivx);
  if (std::abs(new_res) < tol_res)
    return new_p;

  // Iterate to find root
  for (int i = 0; i < max_iterations; i++)
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
    new_res = PResidual(new_p, bx, lambda_left, lambda_right,
        r_left, r_right, ivx);
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
//   same function as in adiabatic_mhd_gr.cpp and hlld_mhd_rel.cpp
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
