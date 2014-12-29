// Conserved-to-primitive inversion for adiabatic MHD in general relativity

// TODO: make inputs const?
// TODO: manually inline functions?

// Primary header
#include "eos.hpp"

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // NAN, sqrt(), abs(), isfinite()

// Athena headers
#include "../fluid.hpp"                       // Fluid
#include "../../athena.hpp"                   // enums, macros, Real
#include "../../athena_arrays.hpp"            // AthenaArray
#include "../../coordinates/coordinates.hpp"  // Coordinates
#include "../../field/field.hpp"              // InterfaceField
#include "../../mesh.hpp"                     // MeshBlock
#include "../../parameter_input.hpp"          // GetReal()

// Declarations
static Real residual(Real w_guess, Real d_norm, Real q_dot_n, Real q_norm_sq,
    Real b_norm_sq, Real q_dot_b_norm_sq, Real gamma_prime);
static Real residual_derivative(Real w_guess, Real d_norm, Real q_norm_sq,
    Real b_norm_sq, Real q_dot_b_norm_sq, Real gamma_prime);
static Real find_root_nr(Real w_initial, Real d_norm, Real q_dot_n, Real q_norm_sq,
    Real b_norm_sq, Real q_dot_b_norm_sq, Real gamma_prime);
static Real quadratic_root(Real a1, Real a0, bool greater_root);
static Real cubic_root_real(Real a2, Real a1, Real a0);
static void quartic_root_minmax(Real a3, Real a2, Real a1, Real a0, Real *pmin_value,
    Real *pmax_value);

// Constructor
// Inputs:
//   pf: pointer to fluid object
//   pin: pointer to runtime inputs
FluidEqnOfState::FluidEqnOfState(Fluid *pf, ParameterInput *pin)
{
  pmy_fluid_ = pf;
  gamma_ = pin->GetReal("fluid", "gamma");
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
//   b: face-centered magnetic field
// Outputs:
//   prim: primitives
//   bcc: cell-centered magnetic field
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

  static int call_number = 0;
  call_number++;
  // Go through cells
  for (int k = kl; k <= ku; k++)
    for (int j = jl; j <= ju; j++)
    {
      pb->pcoord->CellMetric(k, j, g_, g_inv_);
#pragma simd
      for (int i = is-NGHOST; i <= ie+NGHOST; i++)
      {
        // Extract metric
        const Real &g00 = g_(I00,i), &g01 = g_(I01,i), &g02 = g_(I02,i),
            &g03 = g_(I03,i);
        const Real &g10 = g_(I01,i), &g11 = g_(I11,i), &g12 = g_(I12,i),
            &g13 = g_(I13,i);
        const Real &g20 = g_(I02,i), &g21 = g_(I12,i), &g22 = g_(I22,i),
            &g23 = g_(I23,i);
        const Real &g30 = g_(I03,i), &g31 = g_(I13,i), &g32 = g_(I23,i),
            &g33 = g_(I33,i);

        // Extract inverse of metric
        const Real &gi00 = g_inv_(I00,i), &gi01 = g_inv_(I01,i), &gi02 = g_inv_(I02,i),
            &gi03 = g_inv_(I03,i);
        const Real &gi10 = g_inv_(I01,i), &gi11 = g_inv_(I11,i), &gi12 = g_inv_(I12,i),
            &gi13 = g_inv_(I13,i);
        const Real &gi20 = g_inv_(I02,i), &gi21 = g_inv_(I12,i), &gi22 = g_inv_(I22,i),
            &gi23 = g_inv_(I23,i);
        const Real &gi30 = g_inv_(I03,i), &gi31 = g_inv_(I13,i), &gi32 = g_inv_(I23,i),
            &gi33 = g_inv_(I33,i);

        // Extract conserved quantities
        Real &d = cons(IDN,k,j,i);
        Real &e = cons(IEN,k,j,i);
        Real &m1 = cons(IM1,k,j,i);
        Real &m2 = cons(IM2,k,j,i);
        Real &m3 = cons(IM3,k,j,i);

        // Extract face-centered magnetic field
        const Real &bf1_i = b.x1f(k,j,i);
        const Real &bf1_ip1 = b.x1f(k,j,i+1);
        const Real &bf2_j = b.x2f(k,j,i);
        const Real &bf2_jp1 = b.x2f(k,j+1,i);
        const Real &bf3_k = b.x3f(k,j,i);
        const Real &bf3_kp1 = b.x3f(k+1,j,i);

        // Extract cell-centered magnetic field
        Real &b1 = bcc(IB1,k,j,i);
        Real &b2 = bcc(IB2,k,j,i);
        Real &b3 = bcc(IB3,k,j,i);

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

        // Calculate cell-centered magnetic field
        Real interp_param = (pb->x1v(i) - pb->x1f(i)) / pb->dx1f(i);
        b1 = (1.0-interp_param) * bf1_i + interp_param * bf1_ip1;
        interp_param = (pb->x2v(j) - pb->x2f(j)) / pb->dx2f(j);
        b2 = (1.0-interp_param) * bf2_j + interp_param * bf2_jp1;
        interp_param = (pb->x3v(k) - pb->x3f(k)) / pb->dx3f(k);
        b3 = (1.0-interp_param) * bf3_k + interp_param * bf3_kp1;

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
        Real b_norm_1 = alpha * b1;  // (N 5)
        Real b_norm_2 = alpha * b2;  // (N 5)
        Real b_norm_3 = alpha * b3;  // (N 5)
        Real b_norm_sq_a = gi11*b1*b1 + 2.0*gi12*b1*b2 + 2.0*gi13*b1*b3
                         + gi22*b2*b2 + 2.0*gi23*b2*b3
                         + gi33*b3*b3;
        Real b_norm_sq = alpha_sq * b_norm_sq_a;
        Real q_dot_b_norm = alpha * (q1*b1 + q2*b2 + q3*b3);
        Real q_dot_b_norm_sq = SQR(q_dot_b_norm);

        // Construct initial guess for enthalpy W
        Real v_sq = g11*v1_old*v1_old + g22*v2_old*v2_old + g33*v3_old*v3_old
            + 2.0 * (g12*v1_old*v2_old + g13*v1_old*v3_old + g23*v2_old*v3_old);
        Real beta_v = g01*v1_old + g02*v2_old + g03*v3_old;
        Real v_norm_sq = 1.0/alpha_sq * (v_sq + 2.0*beta_v + beta_sq);
        Real gamma_sq = 1.0 / (1.0 - v_norm_sq);
        Real w_initial = gamma_sq * (rho_old + gamma_prime * pgas_old);
        for (int count = 0; count < initial_guess_multiplications; count++)
        {
          Real b_w_term_initial = b_norm_sq + w_initial;
          Real v_sq_initial = (q_norm_sq*SQR(w_initial)
              + q_dot_b_norm_sq*(b_norm_sq+2.0*w_initial))
              / SQR(w_initial * (b_norm_sq+w_initial));     // (N 28)
          if (v_sq_initial >= 1.0)  // guess for W leads to superluminal velocities
            w_initial *= initial_guess_multiplier;
          else
            break;
        }

        // Apply Newton-Raphson method to find new W
        Real w_true = find_root_nr(w_initial, d_norm, q_dot_n, q_norm_sq, b_norm_sq,
            q_dot_b_norm_sq, gamma_prime);

        // Calculate primitives from W
        v_norm_sq = (q_norm_sq*SQR(w_true)
            + q_dot_b_norm_sq*(b_norm_sq+2.0*w_true))
            / (SQR(w_true) * SQR(b_norm_sq + w_true));                // (N 28)
        gamma_sq = 1.0/(1.0 - v_norm_sq);
        Real gamma_rel = std::sqrt(gamma_sq);
        pgas = 1.0/gamma_prime
            * (w_true/gamma_sq - d_norm/std::sqrt(gamma_sq));         // (N 32)
        rho = w_true / gamma_sq - gamma_prime * pgas;
        Real u0 = d_norm / (alpha * rho);                             // (N 21)
        Real u_norm_a = gamma_rel / (w_true + b_norm_sq);
        Real u_norm_b = q_dot_b_norm / w_true;
        Real u_norm_1 = u_norm_a * (q_norm_1 + u_norm_b * b_norm_1);  // (N 31)
        Real u_norm_2 = u_norm_a * (q_norm_2 + u_norm_b * b_norm_2);  // (N 31)
        Real u_norm_3 = u_norm_a * (q_norm_3 + u_norm_b * b_norm_3);  // (N 31)
        Real u1 = u_norm_1 - alpha * gamma_rel * gi01;
        Real u2 = u_norm_2 - alpha * gamma_rel * gi02;
        Real u3 = u_norm_3 - alpha * gamma_rel * gi03;
        v1 = u1/u0;
        v2 = u2/u0;
        v3 = u3/u0;

        // Apply density and pressure floors
        bool floor_applied = false;
        // TODO: set floors at runtime
        const Real pgas_floor = 1.0e-9;
        const Real rho_floor = 1.0e-8;
        if (pgas < pgas_floor)
        {
          pgas = pgas_floor;
          floor_applied = true;
        }
        if (rho < rho_floor)
        {
          rho = rho_floor;
          floor_applied = true;
        }
        if (floor_applied)
        {
          // TODO: determine if momentum should be updated
          Real u_0 = g00*u0 + g01*u1 + g02*u2 + g03*u3;
          //Real u_1 = g01*u0 + g11*u1 + g12*u2 + g13*u3;
          //Real u_2 = g02*u0 + g12*u1 + g22*u2 + g23*u3;
          //Real u_3 = g03*u0 + g13*u1 + g23*u2 + g33*u3;
          Real b_cov0 = u0 * (g01*b1 + g02*b2 + g03*b3)
                      + u1 * (g11*b1 + g12*b2 + g13*b3)
                      + u2 * (g12*b1 + g22*b2 + g23*b3)
                      + u3 * (g13*b1 + g23*b2 + g33*b3);
          Real b_cov1 = (b1 + b_cov0 * u1) / u0;
          Real b_cov2 = (b2 + b_cov0 * u2) / u0;
          Real b_cov3 = (b3 + b_cov0 * u3) / u0;
          Real b_cov_0 = g00*b_cov0 + g01*b_cov1 + g02*b_cov2 + g03*b_cov3;
          Real b_cov_1 = g01*b_cov0 + g11*b_cov1 + g12*b_cov2 + g13*b_cov3;
          Real b_cov_2 = g02*b_cov0 + g12*b_cov1 + g22*b_cov2 + g23*b_cov3;
          Real b_cov_3 = g03*b_cov0 + g13*b_cov1 + g23*b_cov2 + g33*b_cov3;
          Real b_sq = b_cov_0*b_cov0 + b_cov_1*b_cov1 + b_cov_2*b_cov2 + b_cov_3*b_cov3;
          e = (rho + gamma_prime*pgas + b_sq) * u0*u_0 - b_cov0*b_cov_0
              + pgas + 0.5*b_sq;
          //m1 = (rho + gamma_prime*pgas + b_sq) * u0*u_1 - b_cov0*b_cov_1;
          //m2 = (rho + gamma_prime*pgas + b_sq) * u0*u_2 - b_cov0*b_cov_2;
          //m3 = (rho + gamma_prime*pgas + b_sq) * u0*u_3 - b_cov0*b_cov_3;
        }
      }
    }
  return;
}

// Function for calculating relativistic fast wavespeeds
// Inputs: TODO
// Outputs:
//   plambda_plus: value set to most positive wavespeed
//   plambda_minus: value set to most negative wavespeed
// Notes:
//   inputs assume x is transverse direction
//   references Mignone & Bodo 2005, MNRAS 364 126 (MB2005)
//   references Mignone & Bodo 2006, MNRAS 368 1040 (MB2006)
void FluidEqnOfState::FastMagnetosonicSpeedsRel(
    Real rho, Real pgas,
    Real vx, Real vy, Real vz,
    Real ut, Real ux, Real uy, Real uz,
    Real bx, Real by, Real bz,
    Real bcovt, Real bcovx, Real bcovy, Real bcovz,
    Real *plambda_plus, Real *plambda_minus)
{
  // Parameters
  const double v_limit = 1.0e-12;  // squared velocities less than this are considered 0
  const double b_limit = 1.0e-14;  // squared B^x less than this is considered 0

  // Calculate ratio of specific heats
  const Real gamma_adi = gamma_;
  const Real gamma_adi_red = gamma_adi/(gamma_adi-1.0);

  // Calculate intermediate quantities
  Real v_sq = SQR(vx) + SQR(vy) + SQR(vz);
  Real gamma_rel_sq = 1.0/(1.0-v_sq);
  Real rho_h = rho + gamma_adi_red * pgas;
  Real cs_sq = gamma_adi * pgas / rho_h;                              // (MB2005 4)
  Real bcov_sq = -SQR(bcovt) + SQR(bcovx) + SQR(bcovy) + SQR(bcovz);
  Real bx_sq = SQR(bx);

  // Case out on velocity and magnetic field
  if (v_sq < v_limit)  // vanishing velocity
  {
    Real denominator = rho_h + bcov_sq;
    Real a1 = -(bcov_sq + cs_sq * (rho_h + bx_sq)) / denominator;
    Real a0 = cs_sq * bx_sq / denominator;
    Real lambda_sq = quadratic_root(a1, a0, true);                 // (MB2006 57)
    *plambda_plus = std::sqrt(lambda_sq);
    *plambda_minus = -*plambda_plus;
  }
  else  // nonzero velocity
  {
    Real vx_sq = SQR(vx);
    if (bx_sq < b_limit)  // vanishing normal magnetic field
    {
      Real v_dot_b_perp = vy*by + vz*bz;
      Real q = bcov_sq - cs_sq*SQR(v_dot_b_perp);
      Real denominator = rho_h * (cs_sq + gamma_rel_sq*(1.0-cs_sq)) + q;
      Real a1 = -2.0 * rho_h * gamma_rel_sq * vx * (1.0-cs_sq)
          / denominator;
      Real a0 = (rho_h * (-cs_sq + gamma_rel_sq*vx_sq*(1.0-cs_sq)) - q)
          / denominator;
      *plambda_plus = quadratic_root(a1, a0, true);                       // (MB2006 58)
      *plambda_minus = quadratic_root(a1, a0, false);                     // (MB2006 58)
    }
    else  // nonzero normal magnetic field
    {
      Real vx_3 = vx_sq * vx;
      Real vx_4 = SQR(vx_sq);
      Real bcovt_sq = SQR(bcovt);
      Real bcovx_sq = SQR(bcovx);
      Real var_1 = SQR(gamma_rel_sq) * rho_h * (1.0-cs_sq);
      Real var_2 = gamma_rel_sq * (bcov_sq + rho_h * cs_sq);
      Real denominator = var_1 + var_2 - cs_sq * bcovt_sq;
      Real a3 = (-(4.0*var_1+2.0*var_2)*vx + 2.0*cs_sq*bcovt*bcovx)
          / denominator;
      Real a2 = (6.0*var_1*vx_sq + var_2*(vx_sq-1.0)
          + cs_sq*(bcovt_sq-bcovx_sq)) / denominator;
      Real a1 = (-4.0*var_1*vx_3 + 2.0*var_2*vx - 2.0*cs_sq*bcovt*bcovx)
          / denominator;
      Real a0 = (var_1*vx_4 - var_2*vx_sq + cs_sq*bcovx_sq)
        / denominator;
      quartic_root_minmax(a3, a2, a1, a0, plambda_minus, plambda_plus);   // (MB2006 56)
      *plambda_minus = std::max(*plambda_minus, -1.0);
      *plambda_plus = std::min(*plambda_plus, 1.0);
    }
  }
  return;
}

// Function whose value vanishes for correct enthalpy
// Inputs:
//   w_guess: guess for total enthalpy W
//   d_norm: D = alpha * rho * u^0
//   q_dot_n: Q_mu n^mu = -alpha^2 g^{mu 0} T^0_mu
//   q_norm_sq: \tilde{Q}^2 = alpha^2 g^{mu nu} T^0_mu T^0_nu
//                          + alpha^4 (g^{0 mu} T^0_mu)^2
//   b_norm_sq: \mathcal{B}_mu \mathcal{B}^mu = \alpha^2 g_{mu nu} B^mu B^nu
//   q_dot_b_norm_sq: (Q_mu \mathcal{B}^mu)^2 = (alpha^2 T^0_mu B^mu)^2
//   gamma_prime: reduced adiabatic gas constant Gamma' = Gamma/(Gamma-1)
// Outputs:
//   returned value: calculated minus given value of q_dot_n
// Notes:
//   follows Noble et al. 2006, ApJ 641 626 (N)
Real residual(Real w_guess, Real d_norm, Real q_dot_n, Real q_norm_sq, Real b_norm_sq,
    Real q_dot_b_norm_sq, Real gamma_prime)
{
  Real v_norm_sq = (q_norm_sq*SQR(w_guess)
      + q_dot_b_norm_sq*(b_norm_sq+2.0*w_guess))
      / SQR(w_guess*(b_norm_sq+w_guess));                       // (N 28)
  Real gamma_sq = 1.0/(1.0 - v_norm_sq);
  Real pgas = 1.0/gamma_prime
      * (w_guess/gamma_sq - d_norm/std::sqrt(gamma_sq));        // (N 32)
  Real q_dot_n_calc = -0.5*b_norm_sq * (1.0+v_norm_sq)
      + q_dot_b_norm_sq / (2.0*SQR(w_guess)) - w_guess + pgas;  // (N 29)
  return q_dot_n_calc - q_dot_n;
}

// Derivative of residual()
// Inputs:
//   w_guess: guess for total enthalpy W
//   d_norm: D = alpha * rho * u^0
//   q_norm_sq: \tilde{Q}^2 = alpha^2 g^{mu nu} T^0_mu T^0_nu
//                          + alpha^4 (g^{0 mu} T^0_mu)^2
//   b_norm_sq: \mathcal{B}_mu \mathcal{B}^mu = \alpha^2 g_{mu nu} B^mu B^nu
//   q_dot_b_norm_sq: (Q_mu \mathcal{B}^mu)^2 = (alpha^2 T^0_mu B^mu)^2
//   gamma_prime: reduced adiabatic gas constant Gamma' = Gamma/(Gamma-1)
// Outputs:
//   returned value: derivative of calculated value of Q_mu n^mu
// Notes:
//   follows Noble et al. 2006, ApJ 641 626 (N)
Real residual_derivative(Real w_guess, Real d_norm, Real q_norm_sq, Real b_norm_sq,
    Real q_dot_b_norm_sq, Real gamma_prime)
{
  Real w_sq = SQR(w_guess);
  Real w_cu = w_sq * w_guess;
  Real b_w_term = b_norm_sq + w_guess;
  Real b_w_term_sq = SQR(b_w_term);
  Real b_w_term_cu = b_w_term_sq * b_w_term;
  Real v_norm_sq = (q_norm_sq*w_sq
      + q_dot_b_norm_sq*(b_norm_sq+2.0*w_guess))
      / (w_sq*b_w_term_sq);                                                  // (N 28)
  Real gamma_sq = 1.0/(1.0 - v_norm_sq);
  Real gamma_4 = SQR(gamma_sq);
  Real dv_norm_sq_dw_a = 3.0*w_sq + 3.0*b_norm_sq*w_guess + SQR(b_norm_sq);
  Real dv_norm_sq_dw_b = q_norm_sq*w_cu + q_dot_b_norm_sq*dv_norm_sq_dw_a;
  Real dv_norm_sq_dw = -2.0 * dv_norm_sq_dw_b / (w_cu * b_w_term_cu);
  Real dgamma_sq_dw = gamma_4 * dv_norm_sq_dw;
  Real dpgas_dw_a = 0.5 * d_norm * std::sqrt(gamma_sq) - w_guess;
  Real dpgas_dw_b = gamma_sq + dpgas_dw_a * dgamma_sq_dw;
  Real dpgas_dw = 1.0/(gamma_prime*gamma_4) * dpgas_dw_b;
  return -0.5*b_norm_sq*dv_norm_sq_dw - q_dot_b_norm_sq/w_cu - 1.0
      + dpgas_dw;
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
static Real find_root_nr(Real w_initial, Real d_norm, Real q_dot_n, Real q_norm_sq,
    Real b_norm_sq, Real q_dot_b_norm_sq, Real gamma_prime)
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
//   returns abscissa of vertex if there are no real roots
//   follows advice in Numerical Recipes, 3rd ed. (5.6) for avoiding large cancellations
Real quadratic_root(Real a1, Real a0, bool greater_root)
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

// Function for finding real root of monic cubic equation
// Inputs:
//   a2: quadratic coefficient
//   a1: linear coefficient
//   a0: constant coefficient
// Outputs:
//   returned value: a real root
// Notes:
//   solves x^3 + a_2 x^2 + a_1 x + a_0 = 0 for x
//   references Numerical Recipes, 3rd ed. (NR)
Real cubic_root_real(Real a2, Real a1, Real a0)
{
  Real q = (a2*a2 - 3.0*a1) / 9.0;                       // (NR 5.6.10)
  Real r = (2.0*a2*a2*a2 - 9.0*a1*a2 + 27.0*a0) / 54.0;  // (NR 5.6.10)
  if (r*r - q*q*q < 0.0)
  {
    Real theta = acos(r/std::sqrt(q*q*q));                 // (NR 5.6.11)
    return -2.0 * std::sqrt(q) * cos(theta/3.0) - a2/3.0;  // (NR 5.6.12)
  }
  else
  {
    Real a = -copysign(1.0, r)
        * cbrt(std::abs(r) + std::sqrt(r*r - q*q*q));  // (NR 5.6.15)
    Real b = (a != 0.0) ? q/a : 0.0;                   // (NR 5.6.16)
    return a + b - a2/3.0;
  }
}

// Function for finding extremal real roots of monic quartic equation
// Inputs:
//   a3: cubic coefficient
//   a2: quadratic coefficient
//   a1: linear coefficient
//   a0: constant coefficient
// Outputs:
//   pmin_value: value set to least real root
//   pmax_value: value set to greatest real root
// Notes:
//   solves x^4 + a3 x^3 + a2 x^2 + a1 x + a0 = 0 for x
//   uses following procedure:
//     1) eliminate cubic term y^4 + b2 y^2 + b1 y + b0
//     2) construct resolvent cubic z^3 + c2 z^2 + c1 z + c0
//     3) find real root z0 of cubic
//     4) construct quadratics:
//          y^2 + d1 y + d0
//          y^2 + e1 y + e0
//     5) find roots of quadratics
//     6) set minimum and maximum roots of original quartic
void quartic_root_minmax(Real a3, Real a2, Real a1, Real a0, Real *pmin_value,
    Real *pmax_value)
{
  // Step 1: Find reduced quartic coefficients
  Real b2 = a2 - 3.0/8.0*a3*a3;
  Real b1 = a1 - a2*a3/2.0 + a3*a3*a3/8.0;
  Real b0 = a0 - a1*a3/4.0 + a2*a3*a3/16.0 - 3.0/256.0*a3*a3*a3*a3;

  // Step 2: Find resolvent cubic coefficients
  Real c2 = -b2;
  Real c1 = -4.0*b0;
  Real c0 = 4.0*b0*b2 - b1*b1;

  // Step 3: Solve cubic
  Real z0 = cubic_root_real(c2, c1, c0);

  // Step 4: Find quadratic coefficients
  Real d1 = std::sqrt(z0 - b2);
  Real e1 = -d1;
  Real d0, e0;
  if (b1 < 0)
  {
    d0 = z0/2.0 + std::sqrt(z0*z0/4.0 - b0);
    e0 = z0/2.0 - std::sqrt(z0*z0/4.0 - b0);
  }
  else
  {
    d0 = z0/2.0 - std::sqrt(z0*z0/4.0 - b0);
    e0 = z0/2.0 + std::sqrt(z0*z0/4.0 - b0);
  }

  // Step 5: Solve quadratics
  Real y1 = quadratic_root(d1, d0, false);
  Real y2 = quadratic_root(d1, d0, true);
  Real y3 = quadratic_root(e1, e0, false);
  Real y4 = quadratic_root(e1, e0, true);

  // Step 6: Set original quartic roots
  *pmin_value = std::min(y1, y3) - a3/4.0;
  *pmax_value = std::max(y2, y4) - a3/4.0;
  return;
}
