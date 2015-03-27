// Conserved-to-primitive inversion for adiabatic hydrodynamics in general relativity

// TODO: make inputs const?
// TODO: manually inline functions?

// Primary header
#include "eos.hpp"

// C++ headers
#include <cfloat>  // FLT_MIN
#include <cmath>   // NAN, sqrt(), abs(), isfinite(), isnan()

// Athena headers
#include "../fluid.hpp"                       // Fluid
#include "../../athena.hpp"                   // enums, macros, Real
#include "../../athena_arrays.hpp"            // AthenaArray
#include "../../coordinates/coordinates.hpp"  // Coordinates
#include "../../field/field.hpp"              // InterfaceField
#include "../../mesh.hpp"                     // MeshBlock
#include "../../parameter_input.hpp"          // GetReal()

// Declarations
static Real FindRootNR(Real w_initial, Real d_norm, Real q_dot_n, Real q_norm_sq,
    Real gamma_prime);
static Real QNResidual(Real w_guess, Real d_norm, Real q_dot_n, Real q_norm_sq,
    Real gamma_prime);
static Real QNResidualPrime(Real w_guess, Real d_norm, Real q_norm_sq,
    Real gamma_prime);

// Constructor
// Inputs:
//   pf: pointer to fluid object
//   pin: pointer to runtime inputs
FluidEqnOfState::FluidEqnOfState(Fluid *pf, ParameterInput *pin)
{
  pmy_fluid_ = pf;
  gamma_ = pin->GetReal("fluid", "gamma");
  density_floor_ = pin->GetOrAddReal("fluid", "dfloor", 1024*FLT_MIN);
  pressure_floor_ = pin->GetOrAddReal("fluid", "pfloor", 1024*FLT_MIN);
  int ncells1 = pf->pmy_block->block_size.nx1 + 2*NGHOST;
  g_.NewAthenaArray(NMETRIC,ncells1);
  g_inv_.NewAthenaArray(NMETRIC,ncells1);
}

// Destructor
FluidEqnOfState::~FluidEqnOfState()
{
  g_.DeleteAthenaArray();
  g_inv_.DeleteAthenaArray();
}

// Variable inverter
// Inputs:
//   cons: conserved quantities
//   prim_old: primitive quantities from previous half timestep
//   b: face-centered magnetic field (unused)
// Outputs:
//   prim: primitives
//   bcc: cell-centered magnetic fields (unused)
// Notes:
//   follows Noble et al. 2006, ApJ 641 626 (N)
//   implements formulas assuming no magnetic field
void FluidEqnOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
    const AthenaArray<Real> &prim_old, const InterfaceField &b, AthenaArray<Real> &prim,
    AthenaArray<Real> &bcc)
{
  // Parameters
  const Real max_velocity = 1.0 - 1.0e-15;
  const Real initial_guess_multiplier = 10.0;
  const int initial_guess_multiplications = 10;

  // Calculate reduced ratio of specific heats
  const Real gamma_prime = gamma_ / (gamma_ - 1.0);

  // Determine array bounds
  MeshBlock *pb = pmy_fluid_->pmy_block;
  int il = pb->is - NGHOST;
  int iu = pb->ie + NGHOST;
  int jl = pb->js;
  int ju = pb->je;
  int kl = pb->ks;
  int ku = pb->ke;
  if (pb->block_size.nx2 > 1)
  {
    jl -= NGHOST;
    ju += NGHOST;
  }
  if (pb->block_size.nx3 > 1)
  {
    kl -= NGHOST;
    ku += NGHOST;
  }

  // Go through cells
  for (int k = kl; k <= ku; k++)
    for (int j = jl; j <= ju; j++)
    {
      pb->pcoord->CellMetric(k, j, il, iu, g_, g_inv_);
      #pragma simd
      for (int i = il; i <= iu; i++)
      {
        // Extract metric
        Real &g00 = g_(I00,i), &g01 = g_(I01,i), &g02 = g_(I02,i), &g03 = g_(I03,i);
        Real &g10 = g_(I01,i), &g11 = g_(I11,i), &g12 = g_(I12,i), &g13 = g_(I13,i);
        Real &g20 = g_(I02,i), &g21 = g_(I12,i), &g22 = g_(I22,i), &g23 = g_(I23,i);
        Real &g30 = g_(I03,i), &g31 = g_(I13,i), &g32 = g_(I23,i), &g33 = g_(I33,i);

        // Extract inverse of metric
        Real &gi00 = g_inv_(I00,i), &gi01 = g_inv_(I01,i), &gi02 = g_inv_(I02,i),
             &gi03 = g_inv_(I03,i);
        Real &gi10 = g_inv_(I01,i), &gi11 = g_inv_(I11,i), &gi12 = g_inv_(I12,i),
             &gi13 = g_inv_(I13,i);
        Real &gi20 = g_inv_(I02,i), &gi21 = g_inv_(I12,i), &gi22 = g_inv_(I22,i),
             &gi23 = g_inv_(I23,i);
        Real &gi30 = g_inv_(I03,i), &gi31 = g_inv_(I13,i), &gi32 = g_inv_(I23,i),
             &gi33 = g_inv_(I33,i);

        // Extract conserved quantities
        Real &d = cons(IDN,k,j,i);
        Real &e = cons(IEN,k,j,i);
        const Real &m1 = cons(IVX,k,j,i);
        const Real &m2 = cons(IVY,k,j,i);
        const Real &m3 = cons(IVZ,k,j,i);

        // Extract old primitives
        const Real &rho_old = prim_old(IDN,k,j,i);
        const Real &pgas_old = prim_old(IEN,k,j,i);
        const Real &v1_old = prim_old(IVX,k,j,i);
        const Real &v2_old = prim_old(IVY,k,j,i);
        const Real &v3_old = prim_old(IVZ,k,j,i);

        // Extract primitives
        Real &rho = prim(IDN,k,j,i);
        Real &pgas = prim(IEN,k,j,i);
        Real &v1 = prim(IVX,k,j,i);
        Real &v2 = prim(IVY,k,j,i);
        Real &v3 = prim(IVZ,k,j,i);

        // Calculate geometric quantities
        Real alpha = 1.0 / std::sqrt(-gi00);
        Real alpha_sq = alpha*alpha;
        Real beta_sq = g00 + alpha_sq;

        // Calculate variations on conserved quantities
        Real d_norm = alpha * d;  // (N 21)
        Real q0 = alpha * e;  // (N 17)
        Real q1 = alpha * m1;  // (N 17)
        Real q2 = alpha * m2;  // (N 17)
        Real q3 = alpha * m3;  // (N 17)
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
        {
          if (w_initial*w_initial <= q_norm_sq)  // v^2 >= 1 according to (N 28)
            w_initial *= initial_guess_multiplier;
          else
            break;
        }

        // Apply Newton-Raphson method to find new W
        Real w_true = FindRootNR(w_initial, d_norm, q_dot_n, q_norm_sq, gamma_prime);

        // Set density, correcting only conserved density if floor applied
        v_norm_sq = q_norm_sq / (w_true*w_true);  // (N 28)
        gamma_sq = 1.0/(1.0 - v_norm_sq);
        Real gamma_rel = std::sqrt(gamma_sq);
        Real u0 = gamma_rel / alpha;              // (N 21)
        rho = d_norm / gamma_rel;                 // (N 21)
        if (rho < density_floor_ or std::isnan(rho))
        {
          rho = density_floor_;
          d = rho * u0;
        }

        // Set velocity
        Real u_norm_1 = gamma_rel * q_norm_1 / w_true;  // (N 31)
        Real u_norm_2 = gamma_rel * q_norm_2 / w_true;  // (N 31)
        Real u_norm_3 = gamma_rel * q_norm_3 / w_true;  // (N 31)
        Real u1 = u_norm_1 - alpha * gamma_rel * gi01;
        Real u2 = u_norm_2 - alpha * gamma_rel * gi02;
        Real u3 = u_norm_3 - alpha * gamma_rel * gi03;
        v1 = u1 / u0;
        v2 = u2 / u0;
        v3 = u3 / u0;

        // Set pressure, correcting only energy if floor applied
        pgas = 1.0/gamma_prime * (w_true/gamma_sq - rho);  // (N 21,32)
        if (pgas < pressure_floor_ or std::isnan(pgas))
        {
          pgas = pressure_floor_;
          Real u_0 = g00*u0 + g01*u1 + g02*u2 + g03*u3;
          e = (rho + gamma_prime * pgas) * u0 * u_0 + pgas;
        }
      }
    }
  return;
}

// Function for calculating relativistic sound speeds
// Inputs:
// Outputs:
//   plambda_plus: value set to most positive wavespeed
//   plambda_minus: value set to most negative wavespeed
// Notes:
//   same function as in adiabatic_hydro_sr.cpp
//     uses SR formula (should be called in locally flat coordinates)
//   references Mignone & Bodo 2005, MNRAS 364 126 (MB)
void FluidEqnOfState::SoundSpeedsSR(
    Real rho_h, Real pgas, Real vx, Real gamma_lorentz_sq,
    Real *plambda_plus, Real *plambda_minus)
{
  const Real gamma_adi = gamma_;
  Real cs_sq = gamma_adi * pgas / rho_h;                                 // (MB 4)
  Real sigma_s = cs_sq / (gamma_lorentz_sq * (1.0-cs_sq));
  Real relative_speed = std::sqrt(sigma_s * (1.0 + sigma_s - SQR(vx)));
  *plambda_plus = 1.0/(1.0+sigma_s) * (vx + relative_speed);             // (MB 23)
  *plambda_minus = 1.0/(1.0+sigma_s) * (vx - relative_speed);            // (MB 23)
  return;
}

// Function for calculating relativistic sound speeds in arbitrary coordinates
// Inputs:
// Outputs:
//   plambda_plus: value set to most positive wavespeed
//   plambda_minus: value set to most negative wavespeed
// Notes:
//   follows same general procedure as vchar() in phys.c in Harm
//   variables are named as though 1 is normal direction
void FluidEqnOfState::SoundSpeedsGR(
    Real rho_h, Real pgas, Real u0, Real u1,
    Real g00, Real g01, Real g11,
    Real *plambda_plus, Real *plambda_minus)
{
  // Parameters and constants
  const Real discriminant_tol = -1.0e-10;  // values between this and 0 are considered 0
  const Real gamma_adi = gamma_;

  // Calculate comoving sound speed
  Real cs_sq = gamma_adi * pgas / rho_h;

  // Set sound speeds in appropriate coordinates
  Real a = SQR(u0) - (g00 + SQR(u0)) * cs_sq;
  Real b = -2.0 * (u0*u1 - (g01 + u0*u1) * cs_sq);
  Real c = SQR(u1) - (g11 + SQR(u1)) * cs_sq;
  Real d = SQR(b) - 4.0*a*c;
  if (d < 0.0 and d > discriminant_tol)
    d = 0.0;
  Real d_sqrt = std::sqrt(d);
  Real root_1 = (-b + d_sqrt) / (2.0*a);
  Real root_2 = (-b - d_sqrt) / (2.0*a);
  if (root_1 > root_2)
  {
    *plambda_plus = root_1;
    *plambda_minus = root_2;
  }
  else
  {
    *plambda_plus = root_2;
    *plambda_minus = root_1;
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
static Real FindRootNR(Real w_initial, Real d_norm, Real q_dot_n, Real q_norm_sq,
    Real gamma_prime)
{
  // Parameters
  const int max_iterations = 100;         // maximum number of iterations
  const Real tol_w = 1.0e-8 * w_initial;  // absolute tolerance in W
  const Real tol_res = 1.0e-15;           // absolute tolerance in residual

  // Check if root has already been found
  Real new_res = QNResidual(w_initial, d_norm, q_dot_n, q_norm_sq, gamma_prime);
  if (std::abs(new_res) < tol_res)
    return w_initial;

  // Iterate to find root
  Real new_w = w_initial;
  for (int i = 0; i < max_iterations; i++)
  {
    // Prepare needed values
    Real old_w = new_w;
    Real old_res = new_res;
    Real derivative = QNResidualPrime(old_w, d_norm, q_norm_sq, gamma_prime);
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
      new_res = QNResidual(new_w, d_norm, q_dot_n, q_norm_sq, gamma_prime);
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
static Real QNResidual(Real w_guess, Real d_norm, Real q_dot_n, Real q_norm_sq,
    Real gamma_prime)
{
  Real v_norm_sq = q_norm_sq / (w_guess*w_guess);  // (N 28)
  Real gamma_sq = 1.0/(1.0 - v_norm_sq);
  Real pgas = 1.0/gamma_prime
      * (w_guess/gamma_sq - d_norm/std::sqrt(gamma_sq));  // (N 32)
  return -w_guess + pgas - q_dot_n;  // (N 29)
}

// Derivative of QNResidual()
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
static Real QNResidualPrime(Real w_guess, Real d_norm, Real q_norm_sq, Real gamma_prime)
{
  Real v_norm_sq = q_norm_sq / (w_guess*w_guess);  // (N 28)
  Real gamma_sq = 1.0/(1.0 - v_norm_sq);
  Real gamma_4 = gamma_sq*gamma_sq;
  Real d_v_norm_sq_dw = -2.0 * q_norm_sq / (w_guess*w_guess*w_guess);
  Real d_gamma_sq_dw = gamma_4 * d_v_norm_sq_dw;
  Real dpgas_dw = 1.0/(gamma_prime * gamma_4) * (gamma_sq
      + (0.5*d_norm*std::sqrt(gamma_sq) - w_guess) * d_gamma_sq_dw);
  return -1.0 + dpgas_dw;
}
