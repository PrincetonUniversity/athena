// HLLC Riemann solver for special relativistic hydro

// TODO: sort out headers
// TODO: make left and right inputs const

// Temporary includes to make compilation work
#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh.hpp"
#include "../fluid.hpp"
#include "../integrators/integrators.hpp"

/*
// Main header
#include "../integrators/integrators.hpp"

// Standard libraries
#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()

// Other headers
#include "../athena.hpp"         // array access, macros
#include "../athena_arrays.hpp"  // AthenaArray
#include "../fluid.hpp"          // GetGamma()
*/

// Declarations
double quadratic_root(double a1, double a0, bool greater_root);

// Riemann solver
// Inputs:
//   il, iu: lower and upper indices for interfaces
//   prim_left, prim_right: left and right primitive states
// Outputs:
//   flux: fluxes
// Notes:
//   implements HLLC algorithm from Mignone & Bodo 2005, MNRAS 364 126 (MB)
void FluidIntegrator::RiemannSolver(const int il, const int iu, const int ivx,
    const int ivy, const int ivz, AthenaArray<Real> &prim_left,
    AthenaArray<Real> &prim_right, AthenaArray<Real> &flux)
{
  // Extract ratio of specific heats
  const Real gamma_adi = pparent_fluid->GetGamma();
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

    // Extract right primitives
    Real &rho_right = prim_right(IDN,i);
    Real &pgas_right = prim_right(IEN,i);
    Real &vx_right = prim_right(ivx,i);
    Real &vy_right = prim_right(ivy,i);
    Real &vz_right = prim_right(ivz,i);

    // Extract fluxes
    Real &flux_d = flux(IDN,i);
    Real &flux_e = flux(IEN,i);
    Real &flux_mx = flux(ivx,i);
    Real &flux_my = flux(ivy,i);
    Real &flux_mz = flux(ivz,i);

    // Calculate wavespeeds in left region
    Real rho_h_left = rho_left + gamma_adi_red * pgas_left;
    Real cs_sq = gamma_adi * pgas_left / rho_h_left;  // (MB 4)
    Real v_sq_left = vx_left*vx_left + vy_left*vy_left + vz_left*vz_left;
    Real gamma_sq_left = 1.0 / (1.0 - v_sq_left);
    Real sigma_s = cs_sq / (gamma_sq_left * (1.0 - cs_sq));
    Real relative_speed = sqrt(sigma_s * (1.0 + sigma_s - vx_left*vx_left));
    Real lambda_plus_left = 1.0 / (1.0 + sigma_s)
        * (vx_left + relative_speed);  // (MB 23)
    Real lambda_minus_left = 1.0 / (1.0 + sigma_s)
        * (vx_left - relative_speed);  // (MB 23)

    // Calculate wavespeeds in right region
    Real rho_h_right = rho_right + gamma_adi_red * pgas_right;
    cs_sq = gamma_adi * pgas_right / rho_h_right;  // (MB 4)
    Real v_sq_right = vx_right*vx_right + vy_right*vy_right + vz_right*vz_right;
    Real gamma_sq_right = 1.0 / (1.0 - v_sq_right);
    sigma_s = cs_sq / (gamma_sq_right * (1.0 - cs_sq));
    relative_speed = sqrt(sigma_s * (1.0 + sigma_s - vx_right*vx_right));
    Real lambda_plus_right = 1.0 / (1.0 + sigma_s)
        * (vx_right + relative_speed);  // (MB 23)
    Real lambda_minus_right = 1.0 / (1.0 + sigma_s)
        * (vx_right - relative_speed);  // (MB 23)

    // Calculate extremal wavespeeds
    Real lambda_left = std::min(lambda_minus_left, lambda_minus_right);
    Real lambda_right = std::max(lambda_plus_left, lambda_plus_right);

    // Set fluxes in extremal cases
    if (lambda_left >= 0.0)  // left region
    {
      // Calculate intermediate quantities (MB 3)
      Real gamma_sq_rho_h = gamma_sq_left * rho_h_left;
      Real d = sqrt(gamma_sq_left) * rho_left;
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
      Real d = sqrt(gamma_sq_right) * rho_right;
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

    // Calculate select left conserved quantities and fluxes
    Real gamma_sq_rho_h_left = gamma_sq_left * rho_h_left;
    Real e_left = gamma_sq_rho_h_left - pgas_left;      // (MB 3)
    Real mx_left = gamma_sq_rho_h_left * vx_left;       // (MB 3)
    Real flux_e_left = mx_left;                         // (MB 2)
    Real flux_mx_left = mx_left * vx_left + pgas_left;  // (MB 2)

    // Calculate select right conserved quantities and fluxes
    Real gamma_sq_rho_h_right = gamma_sq_right * rho_h_right;
    Real e_right = gamma_sq_rho_h_right - pgas_right;       // (MB 3)
    Real mx_right = gamma_sq_rho_h_right * vx_right;        // (MB 3)
    Real flux_e_right = mx_right;                           // (MB 2)
    Real flux_mx_right = mx_right * vx_right + pgas_right;  // (MB 2)

    // Calculate select HLL conserved quantities and fluxes
    Real denom_inverse_rightleft = 1.0 / (lambda_right - lambda_left);
    Real e_hll = (lambda_right * e_right - lambda_left * e_left + flux_e_left
        - flux_e_right) * denom_inverse_rightleft;  // (MB 9)
    Real mx_hll = (lambda_right * mx_right - lambda_left * mx_left + flux_mx_left
        - flux_mx_right) * denom_inverse_rightleft;  // (MB 9)
    Real flux_e_hll = (lambda_right * flux_e_left - lambda_left * flux_e_right
        + lambda_right * lambda_left * (e_right - e_left))
        * denom_inverse_rightleft;  // (MB 11)
    Real flux_mx_hll = (lambda_right * flux_mx_left - lambda_left * flux_mx_right
        + lambda_right * lambda_left * (mx_right - mx_left))
        * denom_inverse_rightleft;  // (MB 11)

    // Calculate contact wavespeed (MB 18)
    Real lambda_star;
    if (flux_e_hll > TINY_NUMBER || flux_e_hll < -TINY_NUMBER)  // use quadratic formula
      lambda_star = quadratic_root(-(e_hll + flux_mx_hll)/flux_e_hll, mx_hll/flux_e_hll,
          false);
    else  // no quadratic term
      lambda_star = mx_hll / (e_hll + flux_mx_hll);

    // Calculate contact pressure (MB 17)
    Real a = lambda_left * e_left - mx_left;
    Real b = mx_left * (lambda_left - vx_left) - pgas_left;
    Real p_star_left = (a * lambda_star - b) / (1.0 - lambda_left * lambda_star);
    a = lambda_right * e_right - mx_right;
    b = mx_right * (lambda_right - vx_right) - pgas_right;
    Real p_star_right = (a * lambda_star - b) / (1.0 - lambda_right * lambda_star);
    Real p_star = 0.5 * (p_star_left + p_star_right);

    // Set fluxes in contact cases
    if (lambda_star >= 0.0)  // left contact region
    {
      // Calculate remaining left conserved quantities and fluxes
      Real d_left = sqrt(gamma_sq_left) * rho_left;  // (MB 3)
      Real my_left = gamma_sq_rho_h_left * vy_left;  // (MB 3)
      Real mz_left = gamma_sq_rho_h_left * vz_left;  // (MB 3)
      Real flux_d_left = d_left * vx_left;           // (MB 2)
      Real flux_my_left = my_left * vx_left;         // (MB 2)
      Real flux_mz_left = mz_left * vx_left;         // (MB 2)

      // Calculate contact conserved quantities (MB 16)
      Real denom_inverse_leftstar = 1.0 / (lambda_left - lambda_star);
      Real d_star = d_left * (lambda_left - vx_left) * denom_inverse_leftstar;
      Real e_star = (e_left * (lambda_left - vx_left) + p_star * lambda_star - pgas_left
          * vx_left) * denom_inverse_leftstar;
      Real mx_star = (mx_left * (lambda_left - vx_left) + p_star - pgas_left)
          * denom_inverse_leftstar;
      Real my_star = my_left * (lambda_left - vx_left) * denom_inverse_leftstar;
      Real mz_star = mz_left * (lambda_left - vx_left) * denom_inverse_leftstar;

      // Set fluxes (MB 14)
      flux_d = flux_d_left + lambda_left * (d_star - d_left);
      flux_e = flux_e_left + lambda_left * (e_star - e_left);
      flux_mx = flux_mx_left + lambda_left * (mx_star - mx_left);
      flux_my = flux_my_left + lambda_left * (my_star - my_left);
      flux_mz = flux_mz_left + lambda_left * (mz_star - mz_left);
    }
    else  // right contact region
    {
      // Calculate remaining right conserved quantities and fluxes
      Real d_right = sqrt(gamma_sq_right) * rho_right;  // (MB 3)
      Real my_right = gamma_sq_rho_h_right * vy_right;  // (MB 3)
      Real mz_right = gamma_sq_rho_h_right * vz_right;  // (MB 3)
      Real flux_d_right = d_right * vx_right;           // (MB 2)
      Real flux_my_right = my_right * vx_right;         // (MB 2)
      Real flux_mz_right = mz_right * vx_right;         // (MB 2)

      // Calculate contact conserved quantities (MB 16)
      Real denom_inverse_rightstar = 1.0 / (lambda_right - lambda_star);
      Real d_star = d_right * (lambda_right - vx_right) * denom_inverse_rightstar;
      Real e_star = (e_right * (lambda_right - vx_right) + p_star * lambda_star
          - pgas_right * vx_right) * denom_inverse_rightstar;
      Real mx_star = (mx_right * (lambda_right - vx_right) + p_star - pgas_right)
          * denom_inverse_rightstar;
      Real my_star = my_right * (lambda_right - vx_right) * denom_inverse_rightstar;
      Real mz_star = mz_right * (lambda_right - vx_right) * denom_inverse_rightstar;

      // Set fluxes (MB 14)
      flux_d = flux_d_right + lambda_right * (d_star - d_right);
      flux_e = flux_e_right + lambda_right * (e_star - e_right);
      flux_mx = flux_mx_right + lambda_right * (mx_star - mx_right);
      flux_my = flux_my_right + lambda_right * (my_star - my_right);
      flux_mz = flux_mz_right + lambda_right * (mz_star - mz_right);
    }
  }
  return;
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
//   solves x^2 + a_1 x + a_0 = 0 for x
//   discards imaginary parts of answers
//   follows advice in Numerical Recipes (section 5.6) for avoiding large cancellations
double quadratic_root(double a1, double a0, bool greater_root)
{
  if (greater_root)
  {
    if (a1 >= 0.0)
    {
      if (a1*a1 > 4.0*a0)
        return -2.0*a0 / (a1 + std::sqrt(a1*a1 - 4.0*a0));
      else
        return -2.0*a0/a1;
    }
    else
    {
      if (a1*a1 > 4.0*a0)
        return (-a1 + std::sqrt(a1*a1 - 4.0*a0)) / 2.0;
      else
        return -a1/2.0;
    }
  }
  else
  {
    if (a1 >= 0.0)
    {
      if (a1*a1 > 4.0*a0)
        return (-a1 - std::sqrt(a1*a1 - 4.0*a0)) / 2.0;
      else
        return -a1/2.0;
    }
    else
    {
      if (a1*a1 > 4.0*a0)
        return -2.0*a0 / (a1 - std::sqrt(a1*a1 - 4.0*a0));
      else
        return 2.0/a1;
    }
  }
}
