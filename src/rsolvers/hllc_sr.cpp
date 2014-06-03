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
//   P_L, P_R: left and right primitive states
// Outputs:
//   F: fluxes
// Notes:
//   implements HLLC algorithm from Mignone & Bodo 2005, MNRAS 364 126 (MB)
void FluidIntegrator::RiemannSolver(int il, int iu, AthenaArray<Real> &P_L,
    AthenaArray<Real> &P_R, AthenaArray<Real> &F)
{
  // Extract ratio of specific heats
  const Real Gamma = pparent_fluid->GetGamma();
  const Real Gamma_prime = Gamma / (Gamma - 1.0);

  // Go through each interface
#pragma simd
  for (int i = il; i <= iu; i++)
  {
    // Extract left primitives
    Real &rho_L = P_L(IDN,i);
    Real &pgas_L = P_L(IEN,i);
    Real &vx_L = P_L(IVX,i);
    Real &vy_L = P_L(IVY,i);
    Real &vz_L = P_L(IVZ,i);

    // Extract right primitives
    Real &rho_R = P_R(IDN,i);
    Real &pgas_R = P_R(IEN,i);
    Real &vx_R = P_R(IVX,i);
    Real &vy_R = P_R(IVY,i);
    Real &vz_R = P_R(IVZ,i);

    // Extract fluxes
    Real &F_D = F(IDN,i);
    Real &F_E = F(IEN,i);
    Real &F_Mx = F(IVX,i);
    Real &F_My = F(IVY,i);
    Real &F_Mz = F(IVZ,i);

    // Calculate wavespeeds in left region
    Real rho_h_L = rho_L + Gamma_prime * pgas_L;
    Real cs_sq = Gamma * pgas_L / rho_h_L;                                  // (MB 4)
    Real v_sq_L = vx_L*vx_L + vy_L*vy_L + vz_L*vz_L;
    Real gamma_sq_L = 1.0 / (1.0 - v_sq_L);
    Real sigma_s = cs_sq / (gamma_sq_L * (1.0 - cs_sq));
    Real relative_speed = sqrt(sigma_s * (1.0 + sigma_s - vx_L*vx_L));
    Real lambda_plus_L = 1.0 / (1.0 + sigma_s) * (vx_L + relative_speed);   // (MB 23)
    Real lambda_minus_L = 1.0 / (1.0 + sigma_s) * (vx_L - relative_speed);  // (MB 23)

    // Calculate wavespeeds in right region
    Real rho_h_R = rho_R + Gamma_prime * pgas_R;
    cs_sq = Gamma * pgas_R / rho_h_R;                                       // (MB 4)
    Real v_sq_R = vx_R*vx_R + vy_R*vy_R + vz_R*vz_R;
    Real gamma_sq_R = 1.0 / (1.0 - v_sq_R);
    sigma_s = cs_sq / (gamma_sq_R * (1.0 - cs_sq));
    relative_speed = sqrt(sigma_s * (1.0 + sigma_s - vx_R*vx_R));
    Real lambda_plus_R = 1.0 / (1.0 + sigma_s) * (vx_R + relative_speed);   // (MB 23)
    Real lambda_minus_R = 1.0 / (1.0 + sigma_s) * (vx_R - relative_speed);  // (MB 23)

    // Calculate extremal wavespeeds
    Real lambda_L = std::min(lambda_minus_L, lambda_minus_R);
    Real lambda_R = std::max(lambda_plus_L, lambda_plus_R);

    // Set fluxes in extremal cases
    if (lambda_L >= 0.0)  // left region
    {
      // Calculate intermediate quantities (MB 3)
      Real gamma_sq_rho_h = gamma_sq_L * rho_h_L;
      Real D = sqrt(gamma_sq_L) * rho_L;
      Real Mx = gamma_sq_rho_h * vx_L;
      Real My = gamma_sq_rho_h * vy_L;
      Real Mz = gamma_sq_rho_h * vz_L;

      // Set fluxes (MB 2)
      F_D = D * vx_L;
      F_E = Mx;
      F_Mx = Mx * vx_L + pgas_L;
      F_My = My * vx_L;
      F_Mz = Mz * vx_L;
      continue;
    }
    if (lambda_R <= 0.0)  // right region
    {
      // Calculate intermediate quantities (MB 3)
      Real gamma_sq_rho_h = gamma_sq_R * rho_h_R;
      Real D = sqrt(gamma_sq_R) * rho_R;
      Real Mx = gamma_sq_rho_h * vx_R;
      Real My = gamma_sq_rho_h * vy_R;
      Real Mz = gamma_sq_rho_h * vz_R;

      // Set fluxes (MB 2)
      F_D = D * vx_R;
      F_E = Mx;
      F_Mx = Mx * vx_R + pgas_R;
      F_My = My * vx_R;
      F_Mz = Mz * vx_R;
      continue;
    }

    // Calculate select left conserved quantities and fluxes
    Real gamma_sq_rho_h_L = gamma_sq_L * rho_h_L;
    Real E_L = gamma_sq_rho_h_L - pgas_L;  // (MB 3)
    Real Mx_L = gamma_sq_rho_h_L * vx_L;   // (MB 3)
    Real F_E_L = Mx_L;                     // (MB 2)
    Real F_Mx_L = Mx_L * vx_L + pgas_L;    // (MB 2)

    // Calculate select right conserved quantities and fluxes
    Real gamma_sq_rho_h_R = gamma_sq_R * rho_h_R;
    Real E_R = gamma_sq_rho_h_R - pgas_R;  // (MB 3)
    Real Mx_R = gamma_sq_rho_h_R * vx_R;   // (MB 3)
    Real F_E_R = Mx_R;                     // (MB 2)
    Real F_Mx_R = Mx_R * vx_R + pgas_R;    // (MB 2)

    // Calculate select HLL conserved quantities and fluxes
    Real denom_inverse_RL = 1.0 / (lambda_R - lambda_L);
    Real E_HLL = (lambda_R * E_R - lambda_L * E_L + F_E_L - F_E_R)
        * denom_inverse_RL;  // (MB 9)
    Real Mx_HLL = (lambda_R * Mx_R - lambda_L * Mx_L + F_Mx_L - F_Mx_R)
        * denom_inverse_RL;  // (MB 9)
    Real F_E_HLL = (lambda_R * F_E_L - lambda_L * F_E_R + lambda_R * lambda_L
        * (E_R - E_L)) * denom_inverse_RL;  // (MB 11)
    Real F_Mx_HLL = (lambda_R * F_Mx_L - lambda_L * F_Mx_R + lambda_R * lambda_L
        * (Mx_R - Mx_L)) * denom_inverse_RL;  // (MB 11)

    // Calculate contact wavespeed (MB 18)
    Real lambda_star;
    if (F_E_HLL > TINY_NUMBER || F_E_HLL < -TINY_NUMBER)  // can use quadratic formula
      lambda_star = quadratic_root(-(E_HLL + F_Mx_HLL)/F_E_HLL, Mx_HLL/F_E_HLL, false);
    else  // no quadratic term
      lambda_star = Mx_HLL / (E_HLL + F_Mx_HLL);

    // Calculate contact pressure (MB 17)
    Real A = lambda_L * E_L - Mx_L;
    Real B = Mx_L * (lambda_L - vx_L) - pgas_L;
    Real p_star_L = (A * lambda_star - B) / (1.0 - lambda_L * lambda_star);
    A = lambda_R * E_R - Mx_R;
    B = Mx_R * (lambda_R - vx_R) - pgas_R;
    Real p_star_R = (A * lambda_star - B) / (1.0 - lambda_R * lambda_star);
    Real p_star = 0.5 * (p_star_L + p_star_R);

    // Set fluxes in contact cases
    if (lambda_star >= 0.0)  // left contact region
    {
      // Calculate remaining left conserved quantities and fluxes
      Real D_L = sqrt(gamma_sq_L) * rho_L;  // (MB 3)
      Real My_L = gamma_sq_rho_h_L * vy_L;  // (MB 3)
      Real Mz_L = gamma_sq_rho_h_L * vz_L;  // (MB 3)
      Real F_D_L = D_L * vx_L;              // (MB 2)
      Real F_My_L = My_L * vx_L;            // (MB 2)
      Real F_Mz_L = Mz_L * vx_L;            // (MB 2)

      // Calculate contact conserved quantities (MB 16)
      Real denom_inverse_Lstar = 1.0 / (lambda_L - lambda_star);
      Real D_star = D_L * (lambda_L - vx_L) * denom_inverse_Lstar;
      Real E_star = (E_L * (lambda_L - vx_L) + p_star * lambda_star - pgas_L * vx_L)
          * denom_inverse_Lstar;
      Real Mx_star = (Mx_L * (lambda_L - vx_L) + p_star - pgas_L) * denom_inverse_Lstar;
      Real My_star = My_L * (lambda_L - vx_L) * denom_inverse_Lstar;
      Real Mz_star = Mz_L * (lambda_L - vx_L) * denom_inverse_Lstar;

      // Set fluxes (MB 14)
      F_D = F_D_L + lambda_L * (D_star - D_L);
      F_E = F_E_L + lambda_L * (E_star - E_L);
      F_Mx = F_Mx_L + lambda_L * (Mx_star - Mx_L);
      F_My = F_My_L + lambda_L * (My_star - My_L);
      F_Mz = F_Mz_L + lambda_L * (Mz_star - Mz_L);
    }
    else  // right contact region
    {
      // Calculate remaining right conserved quantities and fluxes
      Real D_R = sqrt(gamma_sq_R) * rho_R;  // (MB 3)
      Real My_R = gamma_sq_rho_h_R * vy_R;  // (MB 3)
      Real Mz_R = gamma_sq_rho_h_R * vz_R;  // (MB 3)
      Real F_D_R = D_R * vx_R;              // (MB 2)
      Real F_My_R = My_R * vx_R;            // (MB 2)
      Real F_Mz_R = Mz_R * vx_R;            // (MB 2)

      // Calculate contact conserved quantities (MB 16)
      Real denom_inverse_Rstar = 1.0 / (lambda_R - lambda_star);
      Real D_star = D_R * (lambda_R - vx_R) * denom_inverse_Rstar;
      Real E_star = (E_R * (lambda_R - vx_R) + p_star * lambda_star - pgas_R * vx_R)
          * denom_inverse_Rstar;
      Real Mx_star = (Mx_R * (lambda_R - vx_R) + p_star - pgas_R) * denom_inverse_Rstar;
      Real My_star = My_R * (lambda_R - vx_R) * denom_inverse_Rstar;
      Real Mz_star = Mz_R * (lambda_R - vx_R) * denom_inverse_Rstar;

      // Set fluxes (MB 14)
      F_D = F_D_R + lambda_R * (D_star - D_R);
      F_E = F_E_R + lambda_R * (E_star - E_R);
      F_Mx = F_Mx_R + lambda_R * (Mx_star - Mx_R);
      F_My = F_My_R + lambda_R * (My_star - My_R);
      F_Mz = F_Mz_R + lambda_R * (Mz_star - Mz_R);
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
