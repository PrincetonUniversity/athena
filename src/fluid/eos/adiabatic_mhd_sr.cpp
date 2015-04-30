// Conserved-to-primitive inversion for adiabatic MHD in special relativity

// TODO: make inputs const?
// TODO: check includes

// Primary header
#include "eos.hpp"

// C++ headers
#include <algorithm>  // max(), min()
#include <cfloat>     // FLT_MIN
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
static Real EResidual(Real w_guess, Real d, Real e, Real m_sq, Real b_sq, Real s_sq,
    Real gamma_prime);
static Real EResidualPrime(Real w_guess, Real d, Real m_sq, Real b_sq, Real s_sq,
    Real gamma_prime);
static Real FindRootNR(Real w_initial, Real d, Real e, Real m_sq, Real b_sq, Real s_sq,
    Real gamma_prime);
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
//   prim_old: primitive quantities from previous half timestep (not used)
//   b: face-centered magnetic field
// Outputs:
//   prim: primitives
//   bcc: cell-centered magnetic field
// Notes:
//   follows Mignone & McKinney 2007, MNRAS 378 1118 (MM)
//   follows hlld_sr.c in Athena 4.2 in using W and E rather than W' and E'
void FluidEqnOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
    const AthenaArray<Real> &prim_old, const InterfaceField &b, AthenaArray<Real> &prim,
    AthenaArray<Real> &bcc)
{
  // Calculate reduced ratio of specific heats
  const Real gamma_prime = gamma_/(gamma_-1.0);

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

  // Go through cells
  for (int k = kl; k <= ku; k++)
    for (int j = jl; j <= ju; j++)
    {
      #pragma simd
      for (int i = is-NGHOST; i <= ie+NGHOST; i++)
      {
        // Extract conserved quantities
        Real &d = cons(IDN,k,j,i);
        Real &e = cons(IEN,k,j,i);
        const Real &mx = cons(IM1,k,j,i);
        const Real &my = cons(IM2,k,j,i);
        const Real &mz = cons(IM3,k,j,i);

        // Extract face-centered magnetic field
        const Real &bxm = b.x1f(k,j,i);
        const Real &bxp = b.x1f(k,j,i+1);
        const Real &bym = b.x2f(k,j,i);
        const Real &byp = b.x2f(k,j+1,i);
        const Real &bzm = b.x3f(k,j,i);
        const Real &bzp = b.x3f(k+1,j,i);

        // Extract cell-centered magnetic field
        Real &bx = bcc(IB1,k,j,i);
        Real &by = bcc(IB2,k,j,i);
        Real &bz = bcc(IB3,k,j,i);

        // Calculate cell-centered magnetic field
        Real interp_param = (pb->x1v(i) - pb->x1f(i)) / pb->dx1f(i);
        bx = (1.0-interp_param) * bxm + interp_param * bxp;
        interp_param = (pb->x2v(j) - pb->x2f(j)) / pb->dx2f(j);
        by = (1.0-interp_param) * bym + interp_param * byp;
        interp_param = (pb->x3v(k) - pb->x3f(k)) / pb->dx3f(k);
        bz = (1.0-interp_param) * bzm + interp_param * bzp;

        // Calculate variations on conserved quantities
        Real m_sq = SQR(mx) + SQR(my) + SQR(mz);
        Real b_sq = SQR(bx) + SQR(by) + SQR(bz);
        Real m_dot_b = mx*bx + my*by + mz*bz;
        Real s_sq = SQR(m_dot_b);

        // Construct initial guess for enthalpy W (cf. MM A26-A27)
        Real a1 = 4.0/3.0 * (b_sq - e);
        Real a0 = 1.0/3.0 * (m_sq + b_sq * (b_sq - 2.0*e));
        Real w_initial = quadratic_root(a1, a0, true);

        // Apply Newton-Raphson method to find new W
        Real w_true = FindRootNR(w_initial, d, e, m_sq, b_sq, s_sq, gamma_prime);

        // Extract primitives
        Real &rho = prim(IDN,k,j,i);
        Real &pgas = prim(IEN,k,j,i);
        Real &vx = prim(IVX,k,j,i);
        Real &vy = prim(IVY,k,j,i);
        Real &vz = prim(IVZ,k,j,i);

        // Set density, correcting only conserved density if floor applied
        Real v_sq = (m_sq + s_sq/SQR(w_true) * (2.0*w_true + b_sq))
            / SQR(w_true + b_sq);                                    // (cf. MM A3)
        Real gamma_sq = 1.0/(1.0-v_sq);
        Real gamma_lorentz = std::sqrt(gamma_sq);
        rho = d/gamma_lorentz;                                       // (MM A12)
        if (rho < density_floor_)
        {
          rho = density_floor_;
          d = gamma_lorentz * rho;
        }

        // Set velocity
        vx = (mx + m_dot_b/w_true * bx) / (w_true + b_sq);  // (MM A10)
        vy = (my + m_dot_b/w_true * by) / (w_true + b_sq);  // (MM A10)
        vz = (mz + m_dot_b/w_true * bz) / (w_true + b_sq);  // (MM A10)

        // Set pressure, correcting only energy if floor applied
        Real chi = (1.0 - v_sq) * (w_true - gamma_lorentz * d);  // (cf. MM A11)
        pgas = chi/gamma_prime;                                  // (MM A17)
        if (pgas < pressure_floor_)
        {
          pgas = pressure_floor_;
          Real ut = gamma_lorentz;
          Real ux = ut * vx;
          Real uy = ut * vy;
          Real uz = ut * vz;
          Real bcont = bx*ux + by*uy + bz*uz;
          Real bconx = (bx + bcont * ux) / ut;
          Real bcony = (by + bcont * uy) / ut;
          Real bconz = (bz + bcont * uz) / ut;
          Real b_sq = -SQR(bcont) + SQR(bconx) + SQR(bcony) + SQR(bconz);
          Real w = rho + gamma_prime * pgas + b_sq;
          Real ptot = pgas + 0.5*b_sq;
          e = gamma_sq * w - SQR(bcont) - ptot;
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
//   inputs assume x is normal direction
//   same function as in adiabatic_mhd_gr.cpp
//   references Mignone & Bodo 2005, MNRAS 364 126 (MB2005)
//   references Mignone & Bodo 2006, MNRAS 368 1040 (MB2006)
void FluidEqnOfState::FastMagnetosonicSpeedsSR(
    Real rho, Real pgas, const Real u[4], const Real b[4],
    Real *plambda_plus, Real *plambda_minus)
{
  // Parameters
  const double v_limit = 1.0e-12;  // squared velocities less than this are considered 0
  const double b_limit = 1.0e-14;  // squared B^x less than this is considered 0

  // Calculate ratio of specific heats
  const Real gamma_adi = gamma_;
  const Real gamma_adi_red = gamma_adi/(gamma_adi-1.0);

  // Calculate 3-vector components
  Real vx = u[1]/u[0];
  Real vy = u[2]/u[0];
  Real vz = u[3]/u[0];
  Real bx = b[1]*u[0] - b[0]*u[1];
  Real by = b[2]*u[0] - b[0]*u[2];
  Real bz = b[3]*u[0] - b[0]*u[3];

  // Calculate intermediate quantities
  Real v_sq = SQR(vx) + SQR(vy) + SQR(vz);
  Real gamma_rel_sq = 1.0/(1.0-v_sq);
  Real rho_h = rho + gamma_adi_red * pgas;
  Real cs_sq = gamma_adi * pgas / rho_h;                          // (MB2005 4)
  Real bcov_sq = -SQR(b[0]) + SQR(b[1]) + SQR(b[2]) + SQR(b[3]);
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
      Real bcovt_sq = SQR(b[0]);
      Real bcovx_sq = SQR(b[1]);
      Real var_1 = SQR(gamma_rel_sq) * rho_h * (1.0-cs_sq);
      Real var_2 = gamma_rel_sq * (bcov_sq + rho_h * cs_sq);
      Real denominator = var_1 + var_2 - cs_sq * bcovt_sq;
      Real a3 = (-(4.0*var_1+2.0*var_2)*vx + 2.0*cs_sq*b[0]*b[1])
          / denominator;
      Real a2 = (6.0*var_1*vx_sq + var_2*(vx_sq-1.0)
          + cs_sq*(bcovt_sq-bcovx_sq)) / denominator;
      Real a1 = (-4.0*var_1*vx_3 + 2.0*var_2*vx - 2.0*cs_sq*b[0]*b[1])
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

// Function for finding root of monic quadratic equation
// Inputs:
//   a1: linear coefficient
//   a0: constant coefficient
//   greater_root: flag indicating that larger root is to be returned
//     "larger" does not mean absolute value
// Outputs:
//   returned value: desired root
// Notes:
//   same function as in adiabatic_mhd_gr.cpp, hlld_rel.cpp, and linear_wave_rel.cpp
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

// Function for finding real root of monic cubic equation
// Inputs:
//   a2: quadratic coefficient
//   a1: linear coefficient
//   a0: constant coefficient
// Outputs:
//   returned value: a real root
// Notes:
//   solves x^3 + a_2 x^2 + a_1 x + a_0 = 0 for x
//   same function as in adiabatic_mhd_gr.cpp and linear_wave_rel.cpp
//   references Numerical Recipes, 3rd ed. (NR)
static Real cubic_root_real(Real a2, Real a1, Real a0)
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
//   same function as in adiabatic_mhd_gr.cpp
static void quartic_root_minmax(Real a3, Real a2, Real a1, Real a0, Real *pmin_value,
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
  Real d1 = (z0 - b2 > 0.0) ? std::sqrt(z0 - b2) : 0.0;
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
