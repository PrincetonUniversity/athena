/*// HLLE Riemann solver for special relativistic hydro

// Main header
#include "../integrators/integrators.hpp"

// Standard libraries
#include <cmath>  // sqrt()

// Other headers
#include "../athena.hpp"         // array access, macros
#include "../athena_arrays.hpp"  // AthenaArray
#include "../fluid.hpp"          // GetGamma()

// Riemann solver
// Inputs:
//   il, iu: lower and upper indices for interfaces
//   P_L, P_R: left and right primitive states
// Outputs:
//   F: fluxes
// Notes:
//   implements HLLC algorithm from Mignone & Bodo 2005, MNRAS 364 126 (MB)
void FluidIntegrator::RiemannSolver(int il, int iu, const AthenaArray<Real> &P_L,
    const AthenaArray<Real> &P_R, AthenaArray<Real> &F)
{
  // Extract ratio of specific heats
  const Real Gamma = pparent_fluid->GetGamma();
  Real Gamma_prime = Gamma / (Gamma - 1.0);

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
    Real lambda_L = (lambda_minus_L <= lambda_minus_R) ?
        lambda_minus_L : lambda_minus_R;
    Real lambda_R = (lambda_plus_L >= lambda_plus_R) ? lambda_plus_L : lambda_plus_R;

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
      F_mz = Mz * vx_L;
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
      F_mz = Mz * vx_R;
      continue;
    }

    // Calculate left conserved quantities and fluxes
    Real gamma_sq_rho_h_L = gamma_sq_L * rho_h_L;
    Real D_L = sqrt(gamma_sq_L) * rho_L;   // (MB 3)
    Real E_L = gamma_sq_rho_h_L - pgas_L;  // (MB 3)
    Real Mx_L = gamma_sq_rho_h_L * vx_L;   // (MB 3)
    Real My_L = gamma_sq_rho_h_L * vy_L;   // (MB 3)
    Real Mz_L = gamma_sq_rho_h_L * vz_L;   // (MB 3)
    Real F_D_L = D_L * vx_L;               // (MB 2)
    Real F_E_L = Mx_L;                     // (MB 2)
    Real F_Mx_L = Mx_L * vx_L + pgas_L;    // (MB 2)
    Real F_My_L = My_L * vx_L;             // (MB 2)
    Real F_Mz_L = Mz_L * vx_L;             // (MB 2)

    // Calculate right conserved quantities and fluxes
    Real gamma_sq_rho_h_R = gamma_sq_R * rho_h_R;
    Real D_R = sqrt(gamma_sq_R) * rho_R;   // (MB 3)
    Real E_R = gamma_sq_rho_h_R - pgas_R;  // (MB 3)
    Real Mx_R = gamma_sq_rho_h_R * vx_R;   // (MB 3)
    Real My_R = gamma_sq_rho_h_R * vy_R;   // (MB 3)
    Real Mz_R = gamma_sq_rho_h_R * vz_R;   // (MB 3)
    Real F_D_R = D_R * vx_R;               // (MB 2)
    Real F_E_R = Mx_R;                     // (MB 2)
    Real F_Mx_R = Mx_R * vx_R + pgas_R;    // (MB 2)
    Real F_My_R = My_R * vx_R;             // (MB 2)
    Real F_Mz_R = Mz_R * vx_R;             // (MB 2)

    // Set fluxes (MB 11)
    Real denom_inverse = 1.0 / (lambda_R - lambda_L);
    F_D = (lambda_R * F_D_L - lambda_L * F_D_R + lambda_R * lambda_L * (D_R - D_L))
        * denom_inverse;
    F_E = (lambda_R * F_E_L - lambda_L * F_E_R + lambda_R * lambda_L * (E_R - E_L))
        * denom_inverse;
    F_Mx = (lambda_R * F_Mx_L - lambda_L * F_Mx_R + lambda_R * lambda_L * (Mx_R - Mx_L))
        * denom_inverse;
    F_My = (lambda_R * F_My_L - lambda_L * F_My_R + lambda_R * lambda_L * (My_R - My_L))
        * denom_inverse;
    F_Mz = (lambda_R * F_Mz_L - lambda_L * F_Mz_R + lambda_R * lambda_L * (Mz_R - Mz_L))
        * denom_inverse;
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
}*/
