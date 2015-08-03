// Conserved-to-primitive inversion for adiabatic hydrodynamics in general relativity

// Primary header
#include "eos.hpp"

// C++ headers
#include <algorithm>  // max()
#include <cfloat>     // FLT_MIN
#include <cmath>      // NAN, sqrt(), abs(), isfinite(), isnan(), pow()

// Athena headers
#include "../fluid.hpp"                       // Fluid
#include "../../athena.hpp"                   // enums, macros, Real
#include "../../athena_arrays.hpp"            // AthenaArray
#include "../../mesh.hpp"                     // MeshBlock
#include "../../parameter_input.hpp"          // ParameterInput
#include "../../coordinates/coordinates.hpp"  // Coordinates
#include "../../field/field.hpp"              // InterfaceField

// Declarations
static void PrimitiveToConservedSingle(const AthenaArray<Real> &prim, Real gamma_adi,
    const AthenaArray<Real> &g, const AthenaArray<Real> &gi, int k, int j, int i,
    MeshBlock *pmb, AthenaArray<Real> &cons);
static Real FindRootNR(Real w_initial, Real d, Real q_n, Real qq_sq, Real gamma_adi);
static Real QNResidual(Real w_guess, Real d, Real q_n, Real qq_sq, Real gamma_adi);
static Real QNResidualPrime(Real w_guess, Real d, Real qq_sq, Real gamma_adi);
static void neighbor_average(AthenaArray<Real> &prim, AthenaArray<bool> &problem, int n,
    int k, int j, int i, int kl, int ku, int jl, int ju, int il, int iu);

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
//   bb: face-centered magnetic field (unused)
// Outputs:
//   prim: primitives
//   bb_cc: cell-centered magnetic fields (unused)
// Notes:
//   follows Noble et al. 2006, ApJ 641 626 (N)
//       writing wgas_rel for W = \gamma^2 w
//       writing d for D
//       writing q for Q
//       writing qq for \tilde{Q}
//       writing uu for \tilde{u}
//       writing vv for v
//   implements formulas assuming no magnetic field
void FluidEqnOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
    const AthenaArray<Real> &prim_old, const InterfaceField &bb,
    AthenaArray<Real> &prim, AthenaArray<Real> &bb_cc)
{
  // Parameters
  const Real max_wgas_rel = 1.0e8;
  const Real initial_guess_multiplier = 10.0;
  const int initial_guess_multiplications = 10;

  // Extract ratio of specific heats
  const Real &gamma_adi = gamma_;

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
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
    {
      pmb->pcoord->CellMetric(k, j, il, iu, g_, g_inv_);
      #pragma simd
      for (int i = il; i <= iu; ++i)
      {
        // Extract metric
        const Real
            &g_00 = g_(I00,i), &g_01 = g_(I01,i), &g_02 = g_(I02,i), &g_03 = g_(I03,i),
            &g_10 = g_(I01,i), &g_11 = g_(I11,i), &g_12 = g_(I12,i), &g_13 = g_(I13,i),
            &g_20 = g_(I02,i), &g_21 = g_(I12,i), &g_22 = g_(I22,i), &g_23 = g_(I23,i),
            &g_30 = g_(I03,i), &g_31 = g_(I13,i), &g_32 = g_(I23,i), &g_33 = g_(I33,i);
        const Real &g00 = g_inv_(I00,i), &g01 = g_inv_(I01,i), &g02 = g_inv_(I02,i),
                   &g03 = g_inv_(I03,i), &g10 = g_inv_(I01,i), &g11 = g_inv_(I11,i),
                   &g12 = g_inv_(I12,i), &g13 = g_inv_(I13,i), &g20 = g_inv_(I02,i),
                   &g21 = g_inv_(I12,i), &g22 = g_inv_(I22,i), &g23 = g_inv_(I23,i),
                   &g30 = g_inv_(I03,i), &g31 = g_inv_(I13,i), &g32 = g_inv_(I23,i),
                   &g33 = g_inv_(I33,i);
        Real alpha = 1.0/std::sqrt(-g00);

        // Extract conserved quantities
        const Real &rho_u0 = cons(IDN,k,j,i);
        const Real &t0_0 = cons(IEN,k,j,i);
        const Real &t0_1 = cons(IVX,k,j,i);
        const Real &t0_2 = cons(IVY,k,j,i);
        const Real &t0_3 = cons(IVZ,k,j,i);

        // Calculate variations on conserved quantities
        Real d = alpha * rho_u0;                                               // (N 21)
        Real q_0 = alpha * t0_0;                                               // (N 17)
        Real q_1 = alpha * t0_1;                                               // (N 17)
        Real q_2 = alpha * t0_2;                                               // (N 17)
        Real q_3 = alpha * t0_3;                                               // (N 17)
        Real q_n = -alpha * (g00*q_0 + g01*q_1 + g02*q_2 + g03*q_3);
        Real qq1 = g10*q_0 + g11*q_1 + g12*q_2 + g13*q_3 - alpha * g01 * q_n;
        Real qq2 = g20*q_0 + g21*q_1 + g22*q_2 + g23*q_3 - alpha * g02 * q_n;
        Real qq3 = g30*q_0 + g31*q_1 + g32*q_2 + g33*q_3 - alpha * g03 * q_n;
        Real tmp1 = g00*q_0*q_0 + 2.0*g01*q_0*q_1 + 2.0*g02*q_0*q_2
                      + 2.0*g03*q_0*q_3
                  + g11*q_1*q_1 + 2.0*g12*q_1*q_2 + 2.0*g13*q_1*q_3
                  + g22*q_2*q_2 + 2.0*g23*q_2*q_3
                  + g33*q_3*q_3;
        Real tmp2 = alpha * (g00*q_0 + g01*q_1 + g02*q_2 + g03*q_3);
        Real qq_sq = tmp1 + SQR(tmp2);

        // Extract old primitives
        const Real &rho_old = prim_old(IDN,k,j,i);
        const Real &pgas_old = prim_old(IEN,k,j,i);
        const Real &uu1_old = prim_old(IVX,k,j,i);
        const Real &uu2_old = prim_old(IVY,k,j,i);
        const Real &uu3_old = prim_old(IVZ,k,j,i);

        // Construct initial guess for relativistic gas enthalpy W
        Real tmp = g_11*SQR(uu1_old) + 2.0*g_12*uu1_old*uu2_old
                     + 2.0*g_13*uu1_old*uu3_old
                 + g_22*SQR(uu2_old) + 2.0*g_23*uu2_old*uu3_old
                 + g_33*SQR(uu3_old);
        Real gamma_sq = 1.0 + tmp;
        Real wgas_rel_init =
            gamma_sq * (rho_old + gamma_adi/(gamma_adi-1.0) * pgas_old);
        for (int count = 0; count < initial_guess_multiplications; ++count)
        {
          if (SQR(wgas_rel_init) <= qq_sq)  // v^2 >= 1 according to (N 28)
            wgas_rel_init *= initial_guess_multiplier;
          else
            break;
        }

        // Apply Newton-Raphson method to find new W
        Real wgas_rel_true = FindRootNR(wgas_rel_init, d, q_n, qq_sq, gamma_adi);
        if (not (wgas_rel_true > 0.0 and wgas_rel_true < max_wgas_rel))
        {
          fixed_(k,j,i) = true;
          wgas_rel_true = wgas_rel_init;
        }
        Real vv_sq = qq_sq / SQR(wgas_rel_true);  // (N 28)
        gamma_sq = 1.0/(1.0-vv_sq);
        Real gamma_rel = std::sqrt(gamma_sq);
        if (std::isnan(gamma_rel) or not std::isfinite(gamma_rel) or gamma_rel < 1.0)
        {
          fixed_(k,j,i) = true;
          gamma_rel = 1.0;
          gamma_sq = SQR(gamma_rel);
        }
        Real u0 = gamma_rel/alpha;  // (N 21)

        // Extract primitives
        Real &rho = prim(IDN,k,j,i);
        Real &pgas = prim(IEN,k,j,i);
        Real &uu1 = prim(IVX,k,j,i);
        Real &uu2 = prim(IVY,k,j,i);
        Real &uu3 = prim(IVZ,k,j,i);

        // Set density and pressure
        rho = d/gamma_rel;                                                  // (N 21)
        pgas = (gamma_adi-1.0)/gamma_adi * (wgas_rel_true/gamma_sq - rho);  // (N 21,32)

        // Set velocity
        uu1 = gamma_rel * qq1 / wgas_rel_true;  // (N 31)
        uu2 = gamma_rel * qq2 / wgas_rel_true;  // (N 31)
        uu3 = gamma_rel * qq3 / wgas_rel_true;  // (N 31)

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
          uu1 *= factor;
          uu2 *= factor;
          uu3 *= factor;
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
          PrimitiveToConservedSingle(prim, gamma_adi, g_, g_inv_, k, j, i, pmb, cons);
          fixed_(k,j,i) = false;
        }
    }
  return;
}

// Function for converting all primitives to conserved variables
// Inputs:
//   prim: 3D array of primitives
//   bb_cc: 3D array of cell-centered magnetic fields (unused)
// Outputs:
//   cons: 3D array of conserved variables
void FluidEqnOfState::PrimitiveToConserved(const AthenaArray<Real> &prim,
    const AthenaArray<Real> &bb_cc, AthenaArray<Real> &cons)
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

  // Go through all cells
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
    {
      pmb->pcoord->CellMetric(k, j, il, iu, g_, g_inv_);
      #pragma simd
      for (int i = il; i <= iu; ++i)
        PrimitiveToConservedSingle(prim, gamma_, g_, g_inv_, k, j, i, pmb, cons);
    }
  return;
}

// Function for converting primitives to conserved variables in a single cell
// Inputs:
//   prim: 3D array of primitives
//   gamma_adi: ratio of specific heats
//   g,gi: 1D arrays of metric covariant and contravariant coefficients
//   k,j,i: indices of cell
//   pmb: pointer to MeshBlock
// Outputs:
//   cons: conserved variables set in desired cell
static void PrimitiveToConservedSingle(const AthenaArray<Real> &prim, Real gamma_adi,
    const AthenaArray<Real> &g, const AthenaArray<Real> &gi, int k, int j, int i,
    MeshBlock *pmb, AthenaArray<Real> &cons)
{
  // Extract primitives
  const Real &rho = prim(IDN,k,j,i);
  const Real &pgas = prim(IEN,k,j,i);
  const Real &uu1 = prim(IVX,k,j,i);
  const Real &uu2 = prim(IVY,k,j,i);
  const Real &uu3 = prim(IVZ,k,j,i);

  // Calculate 4-velocity
  Real alpha = std::sqrt(-1.0/gi(I00,i));
  Real tmp = g(I11,i)*uu1*uu1 + 2.0*g(I12,i)*uu1*uu2 + 2.0*g(I13,i)*uu1*uu3
           + g(I22,i)*uu2*uu2 + 2.0*g(I23,i)*uu2*uu3
           + g(I33,i)*uu3*uu3;
  Real gamma = std::sqrt(1.0 + tmp);
  Real u0 = gamma/alpha;
  Real u1 = uu1 - alpha * gamma * gi(I01,i);
  Real u2 = uu2 - alpha * gamma * gi(I02,i);
  Real u3 = uu3 - alpha * gamma * gi(I03,i);
  Real u_0, u_1, u_2, u_3;
  pmb->pcoord->LowerVectorCell(u0, u1, u2, u3, k, j, i, &u_0, &u_1, &u_2, &u_3);

  // Extract conserved quantities
  Real &rho_u0 = cons(IDN,k,j,i);
  Real &t0_0 = cons(IEN,k,j,i);
  Real &t0_1 = cons(IM1,k,j,i);
  Real &t0_2 = cons(IM2,k,j,i);
  Real &t0_3 = cons(IM3,k,j,i);

  // Set conserved quantities
  Real wgas = rho + gamma_adi/(gamma_adi-1.0) * pgas;
  rho_u0 = rho * u0;
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
void FluidEqnOfState::SoundSpeedsSR(Real rho_h, Real pgas, Real vx,
    Real gamma_lorentz_sq, Real *plambda_plus, Real *plambda_minus)
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
void FluidEqnOfState::SoundSpeedsGR(Real rho_h, Real pgas, Real u0, Real u1, Real g00,
    Real g01, Real g11, Real *plambda_plus, Real *plambda_minus)
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
//   d: D = alpha * rho * u^0
//   q_n: Q_mu n^mu = -alpha^2 g^{mu 0} T^0_mu
//   qq_sq: \tilde{Q}^2 = alpha^2 g^{mu nu} T^0_mu T^0_nu
//                          + alpha^4 (g^{0 mu} T^0_mu)^2
//   gamma_adi: ratio of specific heats
// Outputs:
// Notes:
//   returns NAN in event of failure
//   forces W to be positive
static Real FindRootNR(Real w_initial, Real d, Real q_n, Real qq_sq, Real gamma_adi)
{
  // Parameters
  const int max_iterations = 100;         // maximum number of iterations
  const Real tol_w = 1.0e-8 * w_initial;  // absolute tolerance in W
  const Real tol_res = 1.0e-15;           // absolute tolerance in residual

  // Check if root has already been found
  Real new_res = QNResidual(w_initial, d, q_n, qq_sq, gamma_adi);
  if (std::abs(new_res) < tol_res)
    return w_initial;

  // Iterate to find root
  Real new_w = w_initial;
  for (int i = 0; i < max_iterations; ++i)
  {
    // Prepare needed values
    Real old_w = new_w;
    Real old_res = new_res;
    Real derivative = QNResidualPrime(old_w, d, qq_sq, gamma_adi);
    Real delta = -old_res / derivative;

    // Check that update makes sense
    if (!std::isfinite(delta))
      return NAN;

    // Reduce step if root goes out of bounds
    int j;
    for (j = i; j < max_iterations; ++j)
    {
      new_w = old_w + delta;
      if (new_w > 0.0)
        break;
      else
        delta /= 2.0;
    }
    i = j;

    // Reduce step if new value is worse than old
    for (j = i; j < max_iterations; ++j)
    {
      new_res = QNResidual(new_w, d, q_n, qq_sq, gamma_adi);
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
//   w_guess: guess for enthalpy W
//   d: D = alpha * rho * u^0
//   q_n: Q_mu n^mu = -alpha^2 g^{mu 0} T^0_mu
//   qq_sq: \tilde{Q}^2 = alpha^2 g^{mu nu} T^0_mu T^0_nu
//                          + alpha^4 (g^{0 mu} T^0_mu)^2
//   gamma_adi: ratio of specific heats
// Outputs:
//   returned value: calculated minus given value of q_n
// Notes:
//   follows Noble et al. 2006, ApJ 641 626 (N)
//   implements formulas assuming no magnetic field
static Real QNResidual(Real w_guess, Real d, Real q_n, Real qq_sq, Real gamma_adi)
{
  Real v_norm_sq = qq_sq / (w_guess*w_guess);        // (N 28)
  Real gamma_sq = 1.0/(1.0 - v_norm_sq);
  Real pgas = (gamma_adi-1.0)/gamma_adi
      * (w_guess/gamma_sq - d/std::sqrt(gamma_sq));  // (N 32)
  return -w_guess + pgas - q_n;                      // (N 29)
}

// Derivative of QNResidual()
// Inputs:
//   w_guess: guess for enthalpy W
//   d: D = alpha * rho * u^0
//   qq_sq: \tilde{Q}^2 = alpha^2 g^{mu nu} T^0_mu T^0_nu
//                          + alpha^4 (g^{0 mu} T^0_mu)^2
//   gamma_adi: ratio of specific heats
// Outputs:
//   returned value: derivative of calculated value of Q_mu n^mu
// Notes:
//   follows Noble et al. 2006, ApJ 641 626 (N)
//   implements formulas assuming no magnetic field
static Real QNResidualPrime(Real w_guess, Real d, Real qq_sq, Real gamma_adi)
{
  Real v_norm_sq = qq_sq/SQR(w_guess);                             // (N 28)
  Real gamma_sq = 1.0/(1.0-v_norm_sq);
  Real gamma_4 = SQR(gamma_sq);
  Real d_v_norm_sq_dw = -2.0 * qq_sq / (w_guess*SQR(w_guess));
  Real d_gamma_sq_dw = gamma_4 * d_v_norm_sq_dw;
  Real dpgas_dw = (gamma_adi-1.0)/gamma_adi / gamma_4 * (gamma_sq
      + (0.5*d*std::sqrt(gamma_sq) - w_guess) * d_gamma_sq_dw);
  return -1.0 + dpgas_dw;
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
