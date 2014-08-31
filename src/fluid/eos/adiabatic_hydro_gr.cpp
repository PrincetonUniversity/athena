// Conserved-to-primitive inversion for adiabatic hydrodynamics in general relativity

// TODO: make conserved inputs const
// TODO: manually inline functions?

// Primary header
#include "../fluid.hpp"

// C++ headers
#include <cmath>  // NAN, sqrt(), abs(), isfinite()

// Athena headers
#include "../athena.hpp"                   // enums, macros, Real
#include "../athena_arrays.hpp"            // AthenaArray
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../mesh.hpp"                     // MeshBlock

// Declarations
Real find_root_nr(Real w_initial, Real d_norm, Real q_dot_n, Real q_norm_sq,
    Real gamma_prime);
Real residual(Real w_guess, Real d_norm, Real q_dot_n, Real q_norm_sq,
    Real gamma_prime);
Real residual_derivative(Real w_guess, Real d_norm, Real q_norm_sq, Real gamma_prime);

// Variable inverter
// Inputs:
//   cons: conserved quantities
//   prim_old: primitive quantities from previous half timestep
// Outputs:
//   prim: primitives
// Notes:
//   follows Noble et al. 2006, ApJ 641 626 (N)
//   implements formulas assuming no magnetic field
void Fluid::ConservedToPrimitive(AthenaArray<Real> &cons, AthenaArray<Real> &prim_old,
    AthenaArray<Real> &prim)
{
  // Parameters
  const Real max_velocity = 1.0 - 1.0e-15;
  const Real initial_guess_multiplier = 10.0;
  const int initial_guess_multiplications = 10;

  // Extract ratio of specific heats
  const Real gamma_adi = GetGamma();
  const Real gamma_prime = gamma_adi / (gamma_adi - 1.0);

  // Determine array bounds
  MeshBlock *pb = pmy_block;
  int is = pb->is;
  int ie = pb->ie;
  int jl = pb->js;
  int ju = pb->je;
  int kl = pb->ks;
  int ku = pb->ke;
  if (pb->block_size.nx2 > 1)
  {
    jl -= (NGHOST);
    ju += (NGHOST);
  }
  if (pb->block_size.nx3 > 1)
  {
    kl -= (NGHOST);
    ku += (NGHOST);
  }

  // Go through cells
  for (int k = kl; k <= ku; k++)
    for (int j = jl; j <= ju; j++)
    {
      pb->pcoord->CellMetric(k, j, g, g_inv);
#pragma simd
      for (int i = is-NGHOST; i <= ie+NGHOST; i++)
      {
        // Extract metric
        Real &g00 = g(I00,i), &g01 = g(I01,i), &g02 = g(I02,i), &g03 = g(I03,i);
        Real &g10 = g(I01,i), &g11 = g(I11,i), &g12 = g(I12,i), &g13 = g(I13,i);
        Real &g20 = g(I02,i), &g21 = g(I12,i), &g22 = g(I22,i), &g23 = g(I23,i);
        Real &g30 = g(I03,i), &g31 = g(I13,i), &g32 = g(I23,i), &g33 = g(I33,i);

        // Extract inverse of metric
        Real &gi00 = g_inv(I00,i), &gi01 = g_inv(I01,i), &gi02 = g_inv(I02,i),
             &gi03 = g_inv(I03,i);
        Real &gi10 = g_inv(I01,i), &gi11 = g_inv(I11,i), &gi12 = g_inv(I12,i),
             &gi13 = g_inv(I13,i);
        Real &gi20 = g_inv(I02,i), &gi21 = g_inv(I12,i), &gi22 = g_inv(I22,i),
             &gi23 = g_inv(I23,i);
        Real &gi30 = g_inv(I03,i), &gi31 = g_inv(I13,i), &gi32 = g_inv(I23,i),
             &gi33 = g_inv(I33,i);

        // Extract conserved quantities
        Real &d = cons(IDN,k,j,i);
        Real &e = cons(IEN,k,j,i);
        Real &m1 = cons(IVX,k,j,i);
        Real &m2 = cons(IVY,k,j,i);
        Real &m3 = cons(IVZ,k,j,i);

        // Extract old primitives
        Real &rho_old = prim_old(IDN,k,j,i);
        Real &pgas_old = prim_old(IEN,k,j,i);
        Real &v1_old = prim_old(IVX,k,j,i);
        Real &v2_old = prim_old(IVY,k,j,i);
        Real &v3_old = prim_old(IVZ,k,j,i);

        // Extract primitives
        Real &rho = prim(IDN,k,j,i);
        Real &pgas = prim(IEN,k,j,i);
        Real &v1 = prim(IVX,k,j,i);
        Real &v2 = prim(IVY,k,j,i);
        Real &v3 = prim(IVZ,k,j,i);

        // Calculate useful geometric quantities
        Real alpha = 1.0 / std::sqrt(-gi00);
        Real alpha_sq = alpha*alpha;
        Real beta_sq = g00 + alpha_sq;

        // Calculate useful variations on conserved quantities
        Real d_norm = alpha * d;  // (N21)
        Real q0 = alpha * e;  // (N17)
        Real q1 = alpha * m1;  // (N17)
        Real q2 = alpha * m2;  // (N17)
        Real q3 = alpha * m3;  // (N17)
        Real q_dot_n = -alpha * gi00*q0 + gi01*q1 + gi02*q2 + gi03*q3;
        Real q_norm_1 = gi10*q0 + gi11*q1 + gi12*q2 + gi13*q3 - alpha * gi01 * q_dot_n;
        Real q_norm_2 = gi20*q0 + gi21*q1 + gi22*q2 + gi23*q3 - alpha * gi02 * q_dot_n;
        Real q_norm_3 = gi30*q0 + gi31*q1 + gi32*q2 + gi33*q3 - alpha * gi03 * q_dot_n;
        Real q_norm_sq_a = gi00*q0*q0 + 2.0*gi01*q0*q1 + 2.0*gi02*q0*q2 + 2.0*gi03*q0*q3
                         + gi11*q1*q1 + 2.0*gi12*q1*q2 + 2.0*gi13*q1*q3
                         + gi22*q2*q2 + 2.0*gi23*q2*q3
                         + gi33*q3*q3;
        Real q_norm_sq_b = alpha * (gi00*q0 + gi01*q1 + gi02*q2 + gi03*q3);
        Real q_norm_sq = q_norm_sq_a + q_norm_sq_b*q_norm_sq_b;

        // Construct initial guess for enthalpy W
        Real v_sq = g11*v1_old*v1_old + g22*v2_old*v2_old + g33*v3_old*v3_old
            + 2.0 * (g12*v1_old*v2_old + g13*v1_old*v3_old + g23*v2_old*v3_old);
        Real beta_v = g01*v1_old + g02*v2_old + g03*v3_old;
        Real v_norm_sq = 1.0/alpha_sq * (v_sq + 2.0*beta_v + beta_sq);
        Real gamma_sq = 1.0 / (1.0 - v_norm_sq);
        Real w_initial = gamma_sq * (rho_old + gamma_prime * pgas_old);
        for (int count = 0; count < initial_guess_multiplications; count++)
          if (w_initial*w_initial <= q_norm_sq)  // v^2 negative according to (N28)
            w_initial *= initial_guess_multiplier;

        // Apply Newton-Raphson method to find new W
        Real w_true = find_root_nr(w_initial, d_norm, q_dot_n, q_norm_sq, gamma_prime);

        // Calculate primitives from W
        v_norm_sq = q_norm_sq / (w_true*w_true);  // (N28)
        gamma_sq = 1.0/(1.0 - v_norm_sq);
        Real gamma_rel = std::sqrt(gamma_sq);
        pgas = 1.0/gamma_prime
            * (w_true/gamma_sq - d_norm/std::sqrt(gamma_sq));  // (N32)
        rho = w_true / gamma_sq - gamma_prime * pgas;
        Real u0 = d_norm / (alpha * rho);  // (N21)
        Real u_norm_1 = gamma_rel * q_norm_1 / w_true;  // (N31)
        Real u_norm_2 = gamma_rel * q_norm_2 / w_true;  // (N31)
        Real u_norm_3 = gamma_rel * q_norm_3 / w_true;  // (N31)
        Real u1 = u_norm_1 - alpha * gamma_rel * gi01;
        Real u2 = u_norm_2 - alpha * gamma_rel * gi02;
        Real u3 = u_norm_3 - alpha * gamma_rel * gi03;
        v1 = u1/u0;
        v2 = u2/u0;
        v3 = u3/u0;
      }
    }
  return;
}

// Newton-Raphson root finder
// Inputs:
//   w_initial: initial guess for total enthalpy W
//   d_norm: D = alpha * rho * u^0
//   q_dot_n: Q_mu n^mu = -alpha^2 g^{mu 0} T^0_mu
//   q_norm_sq: \tilde{Q}^2 = alpha^2 g^{mu nu} T^0_mu T^0_nu
//                          + alpha^4 (g^{0 mu} T^0_mu)^2
//   gamma_prime: reduced adiabatic gas constant Gamma' = Gamma/(Gamma-1)
// Outputs:
// Notes:
//   returns NAN in event of failure
//   forces W to be positive
Real find_root_nr(Real w_initial, Real d_norm, Real q_dot_n, Real q_norm_sq,
    Real gamma_prime)
{
  // Parameters
  const int max_iterations = 100;         // maximum number of iterations
  const Real tol_w = 1.0e-8 * w_initial;  // absolute tolerance in W
  const Real tol_res = 1.0e-15;           // absolute tolerance in residual

  // Check if root has already been found
  Real new_res = residual(w_initial, d_norm, q_dot_n, q_norm_sq, gamma_prime);
  if (std::abs(new_res) < tol_res)
    return w_initial;

  // Iterate to find root
  Real new_w = w_initial;
  for (int i = 0; i < max_iterations; i++)
  {
    // Prepare needed values
    Real old_w = new_w;
    Real old_res = new_res;
    Real derivative = residual_derivative(old_w, d_norm, q_norm_sq, gamma_prime);
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
      new_res = residual(new_w, d_norm, q_dot_n, q_norm_sq, gamma_prime);
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

// Function whose value vanishes for correct enthalpy
// Inputs:
//   w_guess: guess for total enthalpy W
//   d_norm: D = alpha * rho * u^0
//   q_dot_n: Q_mu n^mu = -alpha^2 g^{mu 0} T^0_mu
//   q_norm_sq: \tilde{Q}^2 = alpha^2 g^{mu nu} T^0_mu T^0_nu
//                          + alpha^4 (g^{0 mu} T^0_mu)^2
//   gamma_prime: reduced adiabatic gas constant Gamma' = Gamma/(Gamma-1)
// Outputs:
//   returned value: calculated minus given value of q_dot_n
// Notes:
//   follows Noble et al. 2006, ApJ 641 626 (N)
//   implements formulas assuming no magnetic field
Real residual(Real w_guess, Real d_norm, Real q_dot_n, Real q_norm_sq, Real gamma_prime)
{
  Real v_norm_sq = q_norm_sq / (w_guess*w_guess);  // (N28)
  Real gamma_sq = 1.0/(1.0 - v_norm_sq);
  Real pgas = 1.0/gamma_prime
      * (w_guess/gamma_sq - d_norm/std::sqrt(gamma_sq));  // (N32)
  return -w_guess + pgas - q_dot_n;  // (N29)
}

// Derivative of residual()
// Inputs:
//   w_guess: guess for total enthalpy W
//   d_norm: D = alpha * rho * u^0
//   q_norm_sq: \tilde{Q}^2 = alpha^2 g^{mu nu} T^0_mu T^0_nu
//                          + alpha^4 (g^{0 mu} T^0_mu)^2
//   gamma_prime: reduced adiabatic gas constant Gamma' = Gamma/(Gamma-1)
// Outputs:
//   returned value: derivative of calculated value of Q_mu n^mu
// Notes:
//   follows Noble et al. 2006, ApJ 641 626 (N)
//   implements formulas assuming no magnetic field
Real residual_derivative(Real w_guess, Real d_norm, Real q_norm_sq, Real gamma_prime)
{
  Real v_norm_sq = q_norm_sq / (w_guess*w_guess);  // (N28)
  Real gamma_sq = 1.0/(1.0 - v_norm_sq);
  Real gamma_4 = gamma_sq*gamma_sq;
  Real d_v_norm_sq_dw = -2.0 * q_norm_sq / (w_guess*w_guess*w_guess);  // (N28)
  Real d_gamma_sq_dw = gamma_4 * d_v_norm_sq_dw;
  Real dpgas_dw = 1.0/(gamma_prime * gamma_4) * (gamma_sq
      + (0.5*d_norm*std::sqrt(gamma_sq) - w_guess) * d_gamma_sq_dw);  // (N32)
  return -1.0 + dpgas_dw;  // (N29)
}
