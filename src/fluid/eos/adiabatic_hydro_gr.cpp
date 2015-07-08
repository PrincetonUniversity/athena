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
static void PrimitiveToConservedSingle(const AthenaArray<Real> &prim,
    const AthenaArray<Real> &g, const AthenaArray<Real> &gi, int k, int j, int i,
    MeshBlock *pmb, Real gamma_prime, AthenaArray<Real> &cons);
static Real FindRootNR(Real w_initial, Real d_norm, Real q_dot_n, Real q_norm_sq,
    Real gamma_prime);
static void neighbor_average(AthenaArray<Real> &prim, AthenaArray<bool> &problem, int n,
    int k, int j, int i, int kl, int ku, int jl, int ju, int il, int iu);
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
  rho_min_ = pin->GetOrAddReal("fluid", "rho_min", density_floor_);
  rho_pow_ = pin->GetOrAddReal("fluid", "rho_pow", 0.0);
  u_min_ = pin->GetOrAddReal("fluid", "u_min", pressure_floor_/(gamma_-1.0));
  u_pow_ = pin->GetOrAddReal("fluid", "u_pow", 0.0);
  gamma_max_ = pin->GetOrAddReal("fluid", "gamma_max", 1000.0);
  int ncells1 = pf->pmy_block->block_size.nx1 + 2*NGHOST;
  g_.NewAthenaArray(NMETRIC, ncells1);
  g_inv_.NewAthenaArray(NMETRIC, ncells1);
  int ncells2 = (pf->pmy_block->block_size.nx2 > 1) ?
      pf->pmy_block->block_size.nx2 + 2*NGHOST : 1;
  int ncells3 = (pf->pmy_block->block_size.nx3 > 1) ?
      pf->pmy_block->block_size.nx3 + 2*NGHOST : 1;
  fixed_.NewAthenaArray(ncells3, ncells2, ncells1);
}

// Destructor
FluidEqnOfState::~FluidEqnOfState()
{
  g_.DeleteAthenaArray();
  g_inv_.DeleteAthenaArray();
  fixed_.DeleteAthenaArray();
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
  const Real max_w = 1.0e8;
  const Real initial_guess_multiplier = 10.0;
  const int initial_guess_multiplications = 10;

  // Calculate reduced ratio of specific heats
  const Real gamma_prime = gamma_ / (gamma_ - 1.0);

  // Determine array bounds
  MeshBlock *pmb = pmy_fluid_->pmy_block;
  int il = pmb->is - NGHOST;
  int iu = pmb->ie + NGHOST;
  int jl = pmb->js;
  int ju = pmb->je;
  int kl = pmb->ks;
  int ku = pmb->ke;
  if (pmb->block_size.nx2 > 1)
  {
    jl -= NGHOST;
    ju += NGHOST;
  }
  if (pmb->block_size.nx3 > 1)
  {
    kl -= NGHOST;
    ku += NGHOST;
  }

  // Go through cells
  for (int k = kl; k <= ku; k++)
    for (int j = jl; j <= ju; j++)
    {
      pmb->pcoord->CellMetric(k, j, il, iu, g_, g_inv_);
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
        const Real &unorm1_old = prim_old(IVX,k,j,i);
        const Real &unorm2_old = prim_old(IVY,k,j,i);
        const Real &unorm3_old = prim_old(IVZ,k,j,i);

        // Extract primitives
        Real &rho = prim(IDN,k,j,i);
        Real &pgas = prim(IEN,k,j,i);
        Real &unorm1 = prim(IVX,k,j,i);
        Real &unorm2 = prim(IVY,k,j,i);
        Real &unorm3 = prim(IVZ,k,j,i);

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
        Real q_dot_n = -alpha * (gi00*q0 + gi01*q1 + gi02*q2 + gi03*q3);
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
        Real tmp = g11*unorm1*unorm1 + 2.0*g12*unorm1*unorm2 + 2.0*g13*unorm1*unorm3
                 + g22*unorm2*unorm2 + 2.0*g23*unorm2*unorm3
                 + g33*unorm3*unorm3;
        Real gamma_sq = 1.0 + tmp;
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
        if (not (w_true > 0.0 and w_true < max_w))
        {
          fixed_(k,j,i) = true;
          w_true = w_initial;
        }
        Real v_norm_sq = q_norm_sq / (w_true*w_true);  // (N 28)
        gamma_sq = 1.0/(1.0 - v_norm_sq);
        Real gamma_rel = std::sqrt(gamma_sq);
        if (std::isnan(gamma_rel) or not std::isfinite(gamma_rel) or gamma_rel < 1.0)
        {
          fixed_(k,j,i) = true;
          gamma_rel = 1.0;
          gamma_sq = SQR(gamma_rel);
        }
        Real u0 = gamma_rel / alpha;              // (N 21)

        // Set density and pressure
        rho = d_norm / gamma_rel;                 // (N 21)
        pgas = 1.0/gamma_prime * (w_true/gamma_sq - rho);  // (N 21,32)

        // Set velocity
        unorm1 = gamma_rel * q_norm_1 / w_true;  // (N 31)
        unorm2 = gamma_rel * q_norm_2 / w_true;  // (N 31)
        unorm3 = gamma_rel * q_norm_3 / w_true;  // (N 31)

        // Apply floors to density and pressure
        Real density_floor_local = rho_min_ * std::pow(pmb->pcoord->x1v(i), rho_pow_);
        density_floor_local = std::max(density_floor_local, density_floor_);
        Real pressure_floor_local =
            (gamma_-1.0) * u_min_ * std::pow(pmb->pcoord->x1v(i), u_pow_);
        pressure_floor_local = std::max(pressure_floor_local, pressure_floor_);
        if (rho < density_floor_local or std::isnan(rho))
        {
          rho = density_floor_local;
          fixed_(k,j,i) = true;
        }
        if (pgas < pressure_floor_local or std::isnan(pgas))
        {
          pgas = pressure_floor_local;
          fixed_(k,j,i) = true;
        }

        // Apply ceiling to velocity
        if (gamma_rel > gamma_max_)
        {
          Real factor = std::sqrt((SQR(gamma_max_)-1.0) / (SQR(gamma_rel)-1.0));
          unorm1 *= factor;
          unorm2 *= factor;
          unorm3 *= factor;
          fixed_(k,j,i) = true;
        }
      }
    }

  // Fix corresponding conserved values if any changes made
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
    {
      pmb->pcoord->CellMetric(k, j, il, iu, g_, g_inv_);
      for (int i = il; i <= iu; ++i)
        if (fixed_(k,j,i))
        {
          PrimitiveToConservedSingle(prim, g_, g_inv_, k, j, i, pmb, gamma_prime, cons);
          fixed_(k,j,i) = false;
        }
    }

  return;
}

// Function for converting all primitives to conserved variables
// Inputs:
//   prim: 3D array of primitives
//   b: 3D array of cell-centered magnetic fields (unused)
// Outputs:
//   cons: 3D array of conserved variables
void FluidEqnOfState::PrimitiveToConserved(const AthenaArray<Real> &prim,
    const AthenaArray<Real> &b, AthenaArray<Real> &cons)
{
  // Prepare index bounds
  MeshBlock *pmb = pmy_fluid_->pmy_block;
  int il = pmb->is - NGHOST;
  int iu = pmb->ie + NGHOST;
  int jl = pmb->js;
  int ju = pmb->je;
  if (pmb->block_size.nx2 > 1)
  {
    jl -= (NGHOST);
    ju += (NGHOST);
  }
  int kl = pmb->ks;
  int ku = pmb->ke;
  if (pmb->block_size.nx3 > 1)
  {
    kl -= (NGHOST);
    ku += (NGHOST);
  }

  // Calculate reduced ratio of specific heats
  Real gamma_prime = gamma_/(gamma_-1.0);

  // Go through all cells
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
    {
      pmb->pcoord->CellMetric(k, j, il, iu, g_, g_inv_);
      #pragma simd
      for (int i = il; i <= iu; ++i)
        PrimitiveToConservedSingle(prim, g_, g_inv_, k, j, i, pmb, gamma_prime, cons);
    }
  return;
}

// Function for converting primitives to conserved variables in a single cell
// Inputs:
//   prim: 3D array of primitives
//   g,gi: 1D arrays of metric covariant and contravariant coefficients
//   k,j,i: indices of cell
//   pmb: pointer to MeshBlock
//   gamma_prime: reduced ratio of specific heats \Gamma/(\Gamma-1)
// Outputs:
//   cons: conserved variables set in desired cell
static void PrimitiveToConservedSingle(const AthenaArray<Real> &prim,
    const AthenaArray<Real> &g, const AthenaArray<Real> &gi, int k, int j, int i,
    MeshBlock *pmb, Real gamma_prime, AthenaArray<Real> &cons)
{
  // Extract primitives
  const Real &rho = prim(IDN,k,j,i);
  const Real &pgas = prim(IEN,k,j,i);
  const Real &unorm1 = prim(IVX,k,j,i);
  const Real &unorm2 = prim(IVY,k,j,i);
  const Real &unorm3 = prim(IVZ,k,j,i);

  // Calculate 4-velocity
  Real alpha = std::sqrt(-1.0/gi(I00,i));
  Real tmp = g(I11,i)*unorm1*unorm1 + 2.0*g(I12,i)*unorm1*unorm2
               + 2.0*g(I13,i)*unorm1*unorm3
           + g(I22,i)*unorm2*unorm2 + 2.0*g(I23,i)*unorm2*unorm3
           + g(I33,i)*unorm3*unorm3;
  Real gamma = std::sqrt(1.0 + tmp);
  Real u0 = gamma / alpha;
  Real u1 = unorm1 - alpha * gamma * gi(I01,i);
  Real u2 = unorm2 - alpha * gamma * gi(I02,i);
  Real u3 = unorm3 - alpha * gamma * gi(I03,i);
  Real u_0, u_1, u_2, u_3;
  pmb->pcoord->LowerVectorCell(u0, u1, u2, u3, k, j, i, &u_0, &u_1, &u_2, &u_3);

  // Extract conserved quantities
  Real &d0 = cons(IDN,k,j,i);
  Real &t0_0 = cons(IEN,k,j,i);
  Real &t0_1 = cons(IM1,k,j,i);
  Real &t0_2 = cons(IM2,k,j,i);
  Real &t0_3 = cons(IM3,k,j,i);

  // Set conserved quantities
  Real wgas = rho + gamma_prime * pgas;
  d0 = rho * u0;
  t0_0 = wgas * u0 * u_0 + pgas;
  t0_1 = wgas * u0 * u_1;
  t0_2 = wgas * u0 * u_2;
  t0_3 = wgas * u0 * u_3;
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

// Function for replacing primitive value in cell with average of neighbors
// Inputs:
//   prim: array of primitives
//   problem: array of flags of problem cells
//   n: IDN, IEN, IVX, IVY, or IVZ
//   k,j,i: indices of cell
//   kl,ku,jl,ju,il,ju: limits of array
// Outputs:
//   prim(index,k,j,i) modified
// Notes
//   average will only include in-bounds, non-NAN neighbors
//   if no such neighbors exist, value will be unmodified
//   same function as in adiabatic_mhd_gr.cpp
//   TODO: decide if function should be kept/implemented
static void neighbor_average(AthenaArray<Real> &prim, AthenaArray<bool> &problem, int n,
    int k, int j, int i, int kl, int ku, int jl, int ju, int il, int iu)
{
  Real neighbor_sum = 0.0;
  int num_neighbors = 0;
  if (i > il and not std::isnan(prim(n,k,j,i-1)) and not problem(k,j,i-1))
  {
    neighbor_sum += prim(n,k,j,i-1);
    num_neighbors += 1;
  }
  if (i < iu and not std::isnan(prim(n,k,j,i+1)) and not problem(k,j,i+1))
  {
    neighbor_sum += prim(n,k,j,i+1);
    num_neighbors += 1;
  }
  if (j > jl and not std::isnan(prim(n,k,j-1,i)) and not problem(k,j-1,i))
  {
    neighbor_sum += prim(n,k,j-1,i);
    num_neighbors += 1;
  }
  if (j < ju and not std::isnan(prim(n,k,j+1,i)) and not problem(k,j+1,i))
  {
    neighbor_sum += prim(n,k,j+1,i);
    num_neighbors += 1;
  }
  if (k > kl and not std::isnan(prim(n,k-1,j,i)) and not problem(k-1,j,i))
  {
    neighbor_sum += prim(n,k-1,j,i);
    num_neighbors += 1;
  }
  if (k < ku and not std::isnan(prim(n,k+1,j,i)) and not problem(k+1,j,i))
  {
    neighbor_sum += prim(n,k+1,j,i);
    num_neighbors += 1;
  }
  if (num_neighbors > 0)
    prim(n,k,j,i) = neighbor_sum / num_neighbors;
  return;
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
