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
    Real bx = b(k,j,i);

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
    pmy_fluid->pf_eos->FastMagnetosonicSpeedsRelativistic(
        rho_left, pgas_left,
        vx_left, vy_left, vz_left,
        ut_left, ux_left, uy_left, uz_left,
        bx, by_left, bz_left,
        bcovt_left, bcovx_left, bcovy_left, bcovz_left,
        lambda_left_plus, lambda_left_minus);                            // (MB 56)
    Real lambda_right_plus, lambda_right_minus;
    pmy_fluid->pf_eos->FastMagnetosonicSpeedsRelativistic(
        rho_right, pgas_right,
        vx_right, vy_right, vz_right,
        ut_right, ux_right, uy_right, uz_right,
        bx, by_right, bz_right,
        bcovt_right, bcovx_right, bcovy_right, bcovz_right,
        lambda_right_plus, lambda_right_minus);                          // (MB 56)
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
    for (int n = 0; n < NWAVE; ++n)
    {
      r_left[n] = lambda_left * cons_left - flux_left;                    // (MUB 12)
      r_right[n] = lambda_right * cons_right - flux_right;                // (MUB 12)
      cons_hll[n] = (r_right[n]-r_left[n]) / (lambda_right-lambda_left);  // (MB 29)
      flux_hll[n] = (lambda_left*r_right[n] - lambda_right*r_left[n])
          / (lambda_right-lambda_left);                                   // (MB 31)
    }

  }
  return;
}

// Function for finding total pressure given conserved quantities in flat spacetime
Real ConsToPFlat(const Real cons[NWAVE])
{
}

// Function whose value vanishes for correct enthalpy
// Inputs:
// Outputs:
// Notes:
//   follows Mignone & McKinney 2007, MNRAS 378 1118 (MM)
//   implementation follows that of hlld_sr.c in Athena 4.2
static Real Residual(Real w_guess, Real d, Real m_sq, Real b_sq, Real s_sq, Real e,
    Real gamma_adi_red)
{
  Real v_sq = (s_sq * (2.0*w_guess + b_sq) + m_sq * SQR(w_guess))
      / (SQR(w_guess+b_sq) * SQR(w_guess));                              // (MM A3)
  Real gamma_lorentz = 1.0 / std::sqrt(1.0 - v_sq);
  Real chi = (1.0 - v_sq) * (w_guess - gamma * d);                       // cf. (MM A11)
  Real p = chi / gamma_adi_red;                                      // (MM A17)
  return w_guess - p + (1.0+v_sq)/2.0*b_sq - s_sq/(2.0*SQR(w_guess))
      - e;                                                               // (MM A1)
}

// Derivative of Residual()
// Inputs:
// Outputs:
// Notes:
//   follows Mignone & McKinney 2007, MNRAS 378 1118 (MM)
//   implementation follows that of hlld_sr.c in Athena 4.2
static Real ResidualDerivative(Real w_guess, Real d, Real m_sq, Real b_sq, Real s_sq,
    Real gamma_adi_red)
{
  Real w_cu = SQR(w_guess) * w_guess;
  Real w_b_term_sq = SQR(w_guess + b_sq);
  Real w_b_term_cu = w_b_term_sq * (w_guess + b_sq);
  Real v_sq = (s_sq * (2.0*w_guess+b_sq) + m_sq*SQR(w_guess))
      / (w_b_term_sq * SQR(w_guess));                          // (MM A3)
  Real gamma_lorentz = 1.0 / std::sqrt(1.0 - v_sq);
  Real chi = (1.0 - v_sq) * (w_guess - gamma * d);             // cf. (MM A11)
  Real dv_sq_dw_a = 3.0*w_guess * (w_guess+b_sq) + SQR(b_sq);
  Real dv_sq_dw_b = s_sq * dv_sq_dw_a + m_sq * w_cu;
  Real dv_sq_dw = -2.0 * dv_sq_dw_b / (w_cu * w_b_term_cu);    // (MM A16)
  Real dchi_dw = 1.0 - v_sq
      - gamma_lorentz/2.0 * (d + 2.0*gamma_lorentz*chi)
      * dv_sq_dw;                                              // (MM A14)
  Real drho_dw = -gamma_lorentz*d/2.0 * dv_sq_dw;              // (MM A15)
  Real dp_dchi = 1.0 / gamma_adi_red;                          // (MM A18)
  Real dp_drho = 0.0;                                          // (MM A18)
  Real dp_dw = dp_dchi * dchi_dw + dp_drho * drho_dw;
  return 1.0 - dp_dw + b_sq/2.0*dv_sq_dw + s_sq/w_cu;          // derivative of (MM A1)
}

// Newton-Raphson root finder
// Inputs:
//   w_initial: initial guess for total enthalpy W
//   d_norm: D = alpha * rho * u^0
//   q_dot_n: Q_mu n^mu = -alpha^2 g^{mu 0} T^0_mu
//   q_norm_sq: \tilde{Q}^2 = alpha^2 g^{mu nu} T^0_mu T^0_nu
//                          + alpha^4 (g^{0 mu} T^0_mu)^2
//   b_norm_sq: \mathcal{B}_mu \mathcal{B}^mu = \alpha^2 g_{mu nu} B^mu B^nu
//   q_dot_b_norm_sq: (Q_mu \mathcal{B}^mu)^2 = (alpha^2 T^0_mu B^mu)^2
//   gamma_prime: reduced adiabatic gas constant Gamma' = Gamma/(Gamma-1)
// Outputs:
// Notes:
//   returns NAN in event of failure
//   forces W to be positive
static Real find_root_nr(Real w_initial, Real d, Real m_sq, Real b_sq,
    Real s_sq, Real q_dot_b_norm_sq, Real gamma_prime)
{
  // Parameters
  const int max_iterations = 100;         // maximum number of iterations
  const Real tol_w = 1.0e-8 * w_initial;  // absolute tolerance in W
  const Real tol_res = 1.0e-15;           // absolute tolerance in residual

  // Check if root has already been found
  Real new_res = residual(w_initial, d_norm, q_dot_n, q_norm_sq, b_norm_sq,
      q_dot_b_norm_sq, gamma_prime);
  if (std::abs(new_res) < tol_res)
    return w_initial;

  // Iterate to find root
  Real new_w = w_initial;
  for (int i = 0; i < max_iterations; i++)
  {
    // Prepare needed values
    Real old_w = new_w;
    Real old_res = new_res;
    Real derivative = residual_derivative(old_w, d_norm, q_norm_sq, b_norm_sq,
        q_dot_b_norm_sq, gamma_prime);
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
      new_res = residual(new_w, d_norm, q_dot_n, q_norm_sq, b_norm_sq, q_dot_b_norm_sq,
          gamma_prime);
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
