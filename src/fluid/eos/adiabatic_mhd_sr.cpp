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
#include "../../coordinates/coordinates.hpp" // Coordinates

// Declarations
static Real EResidual(Real w_guess, Real d, Real e, Real m_sq, Real b_sq, Real s_sq,
    Real gamma_prime);
static Real EResidualPrime(Real w_guess, Real d, Real m_sq, Real b_sq, Real s_sq,
    Real gamma_prime);
static Real quadratic_root(Real a1, Real a0, bool greater_root);

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
  Coordinates *pco = pb->pcoord;
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
      for (int i = is-NGHOST; i <= ie+NGHOST; ++i)
      {
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
        Real interp_param = (pco->x1v(i) - pco->x1f(i)) / pco->dx1f(i);
        bx = (1.0-interp_param) * bxm + interp_param * bxp;
        interp_param = (pco->x2v(j) - pco->x2f(j)) / pco->dx2f(j);
        by = (1.0-interp_param) * bym + interp_param * byp;
        interp_param = (pco->x3v(k) - pco->x3f(k)) / pco->dx3f(k);
        bz = (1.0-interp_param) * bzm + interp_param * bzp;
      }

      #pragma simd
      for (int i = is-NGHOST; i <= ie+NGHOST; ++i)
      {
        // Extract conserved quantities
        Real &d = cons(IDN,k,j,i);
        Real &e = cons(IEN,k,j,i);
        const Real &mx = cons(IM1,k,j,i);
        const Real &my = cons(IM2,k,j,i);
        const Real &mz = cons(IM3,k,j,i);

        // Extract cell-centered magnetic field
        const Real &bx = bcc(IB1,k,j,i);
        const Real &by = bcc(IB2,k,j,i);
        const Real &bz = bcc(IB3,k,j,i);

        // Calculate variations on conserved quantities
        Real m_sq = SQR(mx) + SQR(my) + SQR(mz);
        Real b_sq = SQR(bx) + SQR(by) + SQR(bz);
        Real m_dot_b = mx*bx + my*by + mz*bz;
        Real s_sq = SQR(m_dot_b);

        // Construct initial guess for enthalpy W (cf. MM A26-A27)
        Real a1 = 4.0/3.0 * (b_sq - e);
        Real a0 = 1.0/3.0 * (m_sq + b_sq * (b_sq - 2.0*e));
        Real s2 = SQR(a1) - 4.0*a0;
        Real s = (s2 < 0.0) ? 0.0 : std::sqrt(s2);
        Real w_init = (s2 >= 0.0 and a1 >= 0.0) ? -2.0*a0/(a1+s) : (-a1+s)/2.0;

        // Apply Newton-Raphson method to find new W
        const int num_iterations = 5;
        Real w_new = w_init;
        Real res_new = EResidual(w_new, d, e, m_sq, b_sq, s_sq, gamma_prime);
        for (int n = 0; n < num_iterations; ++n)
        {
          Real w_old = w_new;
          Real res_old = res_new;
          Real derivative = EResidualPrime(w_old, d, m_sq, b_sq, s_sq, gamma_prime);
          Real delta = -res_old / derivative;
          w_new = w_old + delta;
          res_new = EResidual(w_new, d, e, m_sq, b_sq, s_sq, gamma_prime);
        }
        Real w_true = w_new;

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

// Function for converting all primitives to conserved variables
// Inputs:
//   kl,ku,jl,ju,il,iu: index bounds of region to be updated
//   prim: 3D array of primitives
//   b: 3D array of cell-centered magnetic fields
// Outputs:
//   cons: 3D array of conserved variables
void FluidEqnOfState::PrimitiveToConserved(const int kl, const int ku, const int jl,
    const int ju, const int il, const int iu, const AthenaArray<Real> &prim,
    const AthenaArray<Real> &b, AthenaArray<Real> &cons)
{
  // Calculate reduced ratio of specific heats
  Real gamma_adi_red = gamma_ / (gamma_ - 1.0);

  // Go through all cells
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
    {
      #pragma simd
      for (int i = il; i <= iu; ++i)
      {
        // Extract primitives and magnetic fields
        const Real &rho = prim(IDN,k,j,i);
        const Real &pgas = prim(IEN,k,j,i);
        const Real &v1 = prim(IVX,k,j,i);
        const Real &v2 = prim(IVY,k,j,i);
        const Real &v3 = prim(IVZ,k,j,i);
        const Real &bb1 = b(IB1,k,j,i);
        const Real &bb2 = b(IB2,k,j,i);
        const Real &bb3 = b(IB3,k,j,i);

        // Calculate 4-velocity
        Real u0 = 1.0 / std::sqrt(1.0 - SQR(v1) - SQR(v2) - SQR(v3));
        Real u1 = u0 * v1;
        Real u2 = u0 * v2;
        Real u3 = u0 * v3;

        // Calculate 4-magnetic field
        Real b0 = bb1*u1 + bb2*u2 + bb3*u3;
        Real b1 = (bb1 + b0 * u1) / u0;
        Real b2 = (bb2 + b0 * u2) / u0;
        Real b3 = (bb3 + b0 * u3) / u0;
        Real b_sq = -SQR(b0) + SQR(b1) + SQR(b2) + SQR(b3);

        // Extract conserved quantities
        Real &d = cons(IDN,k,j,i);
        Real &e = cons(IEN,k,j,i);
        Real &m1 = cons(IM1,k,j,i);
        Real &m2 = cons(IM2,k,j,i);
        Real &m3 = cons(IM3,k,j,i);

        // Set conserved quantities
        Real wtot = rho + gamma_adi_red * pgas + b_sq;
        Real ptot = pgas + 0.5 * b_sq;
        d = rho * u0;
        e = wtot * u0 * u0 - b0 * b0 - ptot;
        m1 = wtot * u0 * u1 - b0 * b1;
        m2 = wtot * u0 * u2 - b0 * b2;
        m3 = wtot * u0 * u3 - b0 * b3;
      }
    }
  return;
}

// Function for calculating relativistic fast wavespeeds
// Inputs:
//   prim: 1D array of primitive states
//   bbx_vals: 1D array of B^x
//   il,iu: lower and upper x1-indices
//   ivx: type of interface (IVX for x1, IVY for x2, IVZ for x3)
// Outputs:
//   lambdas_p,lambdas_m: 1D arrays set to +/- wavespeeds
// Notes:
//   references Mignone & Bodo 2005, MNRAS 364 126 (MB2005)
//   references Mignone & Bodo 2006, MNRAS 368 1040 (MB2006)
//   references Numerical Recipes, 3rd ed. (NR)
//   follows advice in NR for avoiding large cancellations in solving quadratics
//   same function as in adiabatic_mhd_gr.cpp
void FluidEqnOfState::FastMagnetosonicSpeedsSR(const AthenaArray<Real> &prim,
    const AthenaArray<Real> &bbx_vals, int il, int iu, int ivx,
    AthenaArray<Real> &lambdas_p, AthenaArray<Real> &lambdas_m)
{
  // Parameters
  const double v_limit = 1.0e-12;  // squared velocities less than this are considered 0
  const double b_limit = 1.0e-14;  // squared B^x less than this is considered 0

  // Calculate cyclic permutations of indices
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  // Calculate ratio of specific heats
  const Real gamma_adi = gamma_;
  const Real gamma_adi_red = gamma_adi/(gamma_adi-1.0);

  // Go through states
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract primitives
    const Real &rho = prim(IDN,i);
    const Real &pgas = prim(IEN,i);
    Real u[4], vx, vy, vz;
    if (GENERAL_RELATIVITY)
    {
      u[1] = prim(ivx,i);
      u[2] = prim(ivy,i);
      u[3] = prim(ivz,i);
      u[0] = std::sqrt(1.0 + SQR(u[1]) + SQR(u[2]) + SQR(u[3]));
      vx = u[1]/u[0];
      vy = u[2]/u[0];
      vz = u[3]/u[0];
    }
    else  // SR
    {
      vx = prim(ivx,i);
      vy = prim(ivy,i);
      vz = prim(ivz,i);
      u[0] = std::sqrt(1.0 / (1.0 - SQR(vx) - SQR(vy) - SQR(vz)));
      u[1] = u[0]*vx;
      u[2] = u[0]*vy;
      u[3] = u[0]*vz;
    }
    const Real &bbx = bbx_vals(i);
    const Real &bby = prim(IBY,i);
    const Real &bbz = prim(IBZ,i);

    // Calculate contravariant magnetic field
    Real b[4];
    b[0] = bbx*u[1] + bby*u[2] + bbz*u[3];
    b[1] = (bbx + b[0] * u[1]) / u[0];
    b[2] = (bby + b[0] * u[2]) / u[0];
    b[3] = (bbz + b[0] * u[3]) / u[0];

    // Calculate intermediate quantities
    Real v_sq = SQR(vx) + SQR(vy) + SQR(vz);
    Real gamma_rel_sq = 1.0/(1.0-v_sq);
    Real w_gas = rho + gamma_adi_red * pgas;
    Real cs_sq = gamma_adi * pgas / w_gas;                       // (MB2005 4)
    Real b_sq = -SQR(b[0]) + SQR(b[1]) + SQR(b[2]) + SQR(b[3]);
    Real bbx_sq = SQR(bbx);

    // Calculate wavespeeds in vanishing velocity case (MB2006 57)
    Real lambda_plus_no_v, lambda_minus_no_v;
    {
      Real w_tot = w_gas + b_sq;
      Real a1 = -(b_sq + cs_sq * (w_gas + bbx_sq)) / w_tot;
      Real a0 = cs_sq * bbx_sq / w_tot;
      Real s2 = SQR(a1) - 4.0*a0;
      Real s = (s2 < 0.0) ? 0.0 : std::sqrt(s2);
      Real lambda_sq = 0.5 * (-a1 + s);
      lambda_plus_no_v = std::sqrt(lambda_sq);
      lambda_minus_no_v = -lambda_plus_no_v;
    }

    // Calculate wavespeeds in vanishing normal field case (MB2006 58)
    Real lambda_plus_no_bbx, lambda_minus_no_bbx;
    {
      Real vx_sq = SQR(vx);
      Real v_dot_bb_perp = vy*bby + vz*bbz;
      Real q = b_sq - cs_sq*SQR(v_dot_bb_perp);
      Real denominator = w_gas * (cs_sq + gamma_rel_sq*(1.0-cs_sq)) + q;
      Real a1 = -2.0 * w_gas * gamma_rel_sq * vx * (1.0-cs_sq) / denominator;
      Real a0 = (w_gas * (-cs_sq + gamma_rel_sq*vx_sq*(1.0-cs_sq)) - q) / denominator;
      Real s2 = SQR(a1) - 4.0*a0;
      Real s = (s2 < 0.0) ? 0.0 : std::sqrt(s2);
      lambda_plus_no_bbx = (s2 >= 0.0 and a1 >= 0.0) ? -2.0*a0/(a1+s) : (-a1+s)/2.0;
      lambda_minus_no_bbx = (s2 >= 0.0 and a1 < 0.0) ? -2.0*a0/(a1-s) : (-a1-s)/2.0;
    }

    // Calculate wavespeeds in general case (MB2006 56)
    Real lambda_plus, lambda_minus;
    {
      // Calculate quartic coefficients
      Real vx2 = SQR(vx);
      Real vx3 = vx2 * vx;
      Real vx4 = SQR(vx2);
      Real bt_sq = SQR(b[0]);
      Real bx_sq = SQR(b[1]);
      Real tmp1 = SQR(gamma_rel_sq) * w_gas * (1.0-cs_sq);
      Real tmp2 = gamma_rel_sq * (b_sq + w_gas * cs_sq);
      Real denominator = tmp1 + tmp2 - cs_sq * bt_sq;
      Real a3 = (-(4.0*tmp1+2.0*tmp2)*vx + 2.0*cs_sq*b[0]*b[1]) / denominator;
      Real a2 = (6.0*tmp1*vx2 + tmp2*(vx2-1.0) + cs_sq*(bt_sq-bx_sq)) / denominator;
      Real a1 = (-4.0*tmp1*vx3 + 2.0*tmp2*vx - 2.0*cs_sq*b[0]*b[1]) / denominator;
      Real a0 = (tmp1*vx4 - tmp2*vx2 + cs_sq*bx_sq) / denominator;

      // Calculate reduced quartic coefficients
      Real b2 = a2 - 3.0/8.0*SQR(a3);
      Real b1 = a1 - 1.0/2.0*a2*a3 + 1.0/8.0*a3*SQR(a3);
      Real b0 = a0 - 1.0/4.0*a1*a3 + 1.0/16.0*a2*SQR(a3) - 3.0/256.0*SQR(SQR(a3));

      // Solve reduced quartic equation
      Real y1, y2, y3, y4;
      {
        // Calculate resolvent cubic coefficients
        Real c2 = -b2;
        Real c1 = -4.0*b0;
        Real c0 = 4.0*b0*b2 - SQR(b1);

        // Solve resolvent cubic equation
        Real q = (c2*c2 - 3.0*c1) / 9.0;                       // (NR 5.6.10)
        Real r = (2.0*c2*c2*c2 - 9.0*c1*c2 + 27.0*c0) / 54.0;  // (NR 5.6.10)
        Real q3 = q*q*q;
        Real r2 = SQR(r);
        Real s2 = r2 - q3;
        Real z0;
        if (s2 < 0.0)
        {
          Real theta = std::acos(r/std::sqrt(q3));             // (NR 5.6.11)
          z0 = -2.0 * std::sqrt(q) * cos(theta/3.0) - c2/3.0;  // (NR 5.6.12)
        }
        else
        {
          Real s = std::sqrt(s2);
          Real aa = -copysign(1.0, r) * cbrt(std::abs(r) + s);  // (NR 5.6.15)
          Real bb = (aa != 0.0) ? q/aa : 0.0;                   // (NR 5.6.16)
          z0 = aa + bb - c2/3.0;
        }

        // Calculate quadratic coefficients
        Real d1 = (z0-b2 > 0.0) ? std::sqrt(z0-b2) : 0.0;
        Real e1 = -d1;
        Real s = std::sqrt(SQR(z0)/4.0 - b0);
        Real d0 = (b1 < 0) ? 0.5*z0+s : 0.5*z0-s;
        Real e0 = (b1 < 0) ? 0.5*z0-s : 0.5*z0+s;

        // Solve quadratic equations
        s2 = SQR(d1) - 4.0*d0;
        s = (s2 < 0.0) ? 0.0 : std::sqrt(s2);
        y1 = (s2 >= 0.0 and d1 < 0.0) ? -2.0*d0/(d1-s) : (-d1-s)/2.0;
        y2 = (s2 >= 0.0 and d1 >= 0.0) ? -2.0*d0/(d1+s) : (-d1+s)/2.0;
        s2 = SQR(e1) - 4.0*e0;
        s = (s2 < 0.0) ? 0.0 : std::sqrt(s2);
        y3 = (s2 >= 0.0 and e1 < 0.0) ? -2.0*e0/(e1-s) : (-e1-s)/2.0;
        y4 = (s2 >= 0.0 and e1 >= 0.0) ? -2.0*e0/(e1+s) : (-e1+s)/2.0;
      }

      // Calculate extremal original quartic roots
      lambda_minus = std::min(y1, y3) - a3/4.0;
      lambda_plus = std::max(y2, y4) - a3/4.0;

      // Ensure wavespeeds are not superluminal
      lambda_minus = std::max(lambda_minus, -1.0);
      lambda_plus = std::min(lambda_plus, 1.0);
    }

    // Set wavespeeds based on velocity and magnetic field
    if (v_sq < v_limit)
    {
      lambdas_p(i) = lambda_plus_no_v;
      lambdas_m(i) = lambda_minus_no_v;
    }
    else if (bbx_sq < b_limit)
    {
      lambdas_p(i) = lambda_plus_no_bbx;
      lambdas_m(i) = lambda_minus_no_bbx;
    }
    else
    {
      lambdas_p(i) = lambda_plus;
      lambdas_m(i) = lambda_minus;
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
