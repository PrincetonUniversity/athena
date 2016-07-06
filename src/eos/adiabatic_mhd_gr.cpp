// Conserved-to-primitive inversion for adiabatic MHD in general relativity

// Primary header
#include "eos.hpp"

// C++ headers
#include <algorithm>  // max(), min()
#include <cfloat>     // FLT_MIN
#include <cmath>      // NAN, sqrt(), abs(), isfinite(), isnan(), pow()

// Athena headers
#include "../athena.hpp"                   // enums, macros
#include "../athena_arrays.hpp"            // AthenaArray
#include "../parameter_input.hpp"          // ParameterInput
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../field/field.hpp"              // FaceField
#include "../mesh/mesh.hpp"                // MeshBlock

// Declarations
static bool ConservedToPrimitiveSingle(const AthenaArray<Real> &cons,
    const AthenaArray<Real> &prim_old, Real gamma_adi, const AthenaArray<Real> &bb_cc,
    const AthenaArray<Real> &g, const AthenaArray<Real> &gi, int k, int j, int i,
    AthenaArray<Real> &prim);
static void PrimitiveToConservedSingle(const AthenaArray<Real> &prim, Real gamma_adi,
    const AthenaArray<Real> &bb_cc, const AthenaArray<Real> &g,
    const AthenaArray<Real> &gi, int k, int j, int i, AthenaArray<Real> &cons,
    Coordinates *pco);
static Real QNResidual(Real w_guess, Real d_norm, Real q_dot_n, Real q_norm_sq,
    Real bbb_sq, Real q_bbb_sq, Real gamma_prime);
static Real QNResidualPrime(Real w_guess, Real d_norm, Real q_norm_sq, Real bbb_sq,
    Real q_bbb_sq, Real gamma_prime);
static void neighbor_average(AthenaArray<Real> &prim, int n, int k, int j, int i,
    int kl, int ku, int jl, int ju, int il, int iu);

//--------------------------------------------------------------------------------------

// Constructor
// Inputs:
//   pmb: pointer to MeshBlock
//   pin: pointer to runtime inputs
EquationOfState::EquationOfState(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block_ = pmb;
  gamma_ = pin->GetReal("hydro", "gamma");
  density_floor_ = pin->GetOrAddReal("hydro", "dfloor", 1024*FLT_MIN);
  pressure_floor_ = pin->GetOrAddReal("hydro", "pfloor", 1024*FLT_MIN);
  rho_min_ = pin->GetOrAddReal("hydro", "rho_min", density_floor_);
  rho_pow_ = pin->GetOrAddReal("hydro", "rho_pow", 0.0);
  u_min_ = pin->GetOrAddReal("hydro", "u_min", pressure_floor_/(gamma_-1.0));
  u_pow_ = pin->GetOrAddReal("hydro", "u_pow", 0.0);
  rho_pmag_min_ = pin->GetOrAddReal("hydro", "rho_pmag_min", 0.0);
  u_pmag_min_ = pin->GetOrAddReal("hydro", "u_pmag_min", 0.0);
  gamma_max_ = pin->GetOrAddReal("hydro", "gamma_max", 1000.0);
  int ncells1 = pmb->block_size.nx1 + 2*NGHOST;
  g_.NewAthenaArray(NMETRIC, ncells1);
  g_inv_.NewAthenaArray(NMETRIC, ncells1);
  int ncells2 = (pmb->block_size.nx2 > 1) ? pmb->block_size.nx2 + 2*NGHOST : 1;
  int ncells3 = (pmb->block_size.nx3 > 1) ? pmb->block_size.nx3 + 2*NGHOST : 1;
  fixed_.NewAthenaArray(ncells3, ncells2, ncells1);
}

//--------------------------------------------------------------------------------------

// Destructor
EquationOfState::~EquationOfState()
{
  g_.DeleteAthenaArray();
  g_inv_.DeleteAthenaArray();
  fixed_.DeleteAthenaArray();
}

//--------------------------------------------------------------------------------------

// Variable inverter
// Inputs:
//   cons: conserved quantities
//   prim_old: primitive quantities from previous half timestep
//   bb: face-centered magnetic field
//   pco: pointer to Coordinates
//   is,ie,js,je,ks,ke: index bounds of region to be updated
// Outputs:
//   prim: primitives
//   bb_cc: cell-centered magnetic field
void EquationOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
    const AthenaArray<Real> &prim_old, const FaceField &bb,
    AthenaArray<Real> &prim, AthenaArray<Real> &bb_cc, Coordinates *pco, int is, int ie,
    int js, int je, int ks, int ke)
{
  // Extract ratio of specific heats
  const Real &gamma_adi = gamma_;

  // Interpolate magnetic field from faces to cell centers
  pmy_block_->pfield->CalculateCellCenteredField(bb, bb_cc, pco, is, ie, js,
      je, ks, ke);

  // Go through cells for first pass
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
    {
      pco->CellMetric(k, j, is, ie, g_, g_inv_);
      #pragma simd
      for (int i = is; i <= ie; ++i)
        fixed_(k,j,i) = ConservedToPrimitiveSingle(cons, prim_old, gamma_adi, bb_cc, g_,
            g_inv_, k, j, i, prim);
    }

  // Adjust conserved values if fix was applied
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
    {
      bool fix_in_row = false;
      for (int i = is; i <= ie; ++i)
        if (fixed_(k,j,i))
          fix_in_row = true;
      if (fix_in_row)
      {
        pco->CellMetric(k, j, is, ie, g_, g_inv_);
        for (int i = is; i <= ie; ++i)
          if (fixed_(k,j,i))
          {
            PrimitiveToConservedSingle(prim, gamma_adi, bb_cc, g_, g_inv_, k, j, i,
                cons, pco);
            fixed_(k,j,i) = false;
          }
      }
    }

  // Apply limits, adjusting conserved and then primitive quantities as necessary
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
    {
      pco->CellMetric(k, j, is, ie, g_, g_inv_);
      for (int i = is; i <= ie; ++i)
      {
        // Extract primitives
        Real &rho = prim(IDN,k,j,i);
        Real &pgas = prim(IEN,k,j,i);
        Real &uu1 = prim(IVX,k,j,i);
        Real &uu2 = prim(IVY,k,j,i);
        Real &uu3 = prim(IVZ,k,j,i);

        // Calculate 4-velocity
        Real alpha = std::sqrt(-1.0/g_inv_(I00,i));
        Real tmp = g_(I11,i)*uu1*uu1 + 2.0*g_(I12,i)*uu1*uu2 + 2.0*g_(I13,i)*uu1*uu3
                 + g_(I22,i)*uu2*uu2 + 2.0*g_(I23,i)*uu2*uu3
                 + g_(I33,i)*uu3*uu3;
        Real gamma = std::sqrt(1.0 + tmp);
        Real u0 = gamma/alpha;
        Real u1 = uu1 - alpha * gamma * g_inv_(I01,i);
        Real u2 = uu2 - alpha * gamma * g_inv_(I02,i);
        Real u3 = uu3 - alpha * gamma * g_inv_(I03,i);
        Real u_0, u_1, u_2, u_3;
        pco->LowerVectorCell(u0, u1, u2, u3, k, j, i, &u_0, &u_1, &u_2, &u_3);

        // Extract cell-centered magnetic field
        const Real &bb1 = bb_cc(IB1,k,j,i);
        const Real &bb2 = bb_cc(IB2,k,j,i);
        const Real &bb3 = bb_cc(IB3,k,j,i);

        // Calculate magnetic pressure
        Real b0 = g_(I01,i)*u0*bb1 + g_(I02,i)*u0*bb2 + g_(I03,i)*u0*bb3
                + g_(I11,i)*u1*bb1 + g_(I12,i)*u1*bb2 + g_(I13,i)*u1*bb3
                + g_(I12,i)*u2*bb1 + g_(I22,i)*u2*bb2 + g_(I23,i)*u2*bb3
                + g_(I13,i)*u3*bb1 + g_(I23,i)*u3*bb2 + g_(I33,i)*u3*bb3;
        Real b1 = (bb1 + b0 * u1) / u0;
        Real b2 = (bb2 + b0 * u2) / u0;
        Real b3 = (bb3 + b0 * u3) / u0;
        Real b_0, b_1, b_2, b_3;
        pco->LowerVectorCell(b0, b1, b2, b3, k, j, i, &b_0, &b_1, &b_2, &b_3);
        Real b_sq = b0*b_0 + b1*b_1 + b2*b_2 + b3*b_3;
        Real pmag = 0.5*b_sq;

        // Calculate floors for density and pressure
        Real density_floor_local = density_floor_;
        if (rho_pow_ != 0.0)
          density_floor_local = std::max(density_floor_local,
              rho_min_*std::pow(pco->x1v(i),rho_pow_));
        density_floor_local = std::max(density_floor_local, rho_pmag_min_*pmag);
        Real pressure_floor_local = pressure_floor_;
        if (u_pow_ != 0.0)
          pressure_floor_local = std::max(pressure_floor_local,
              (gamma_-1.0)*u_min_*std::pow(pco->x1v(i),u_pow_));
        pressure_floor_local = std::max(pressure_floor_local,
            (gamma_-1.0)*u_pmag_min_*pmag);

        // Adjust conserved quantities to add unmagnetized density and/or pressure at
        // rest in normal frame
        Real rho_add = std::max(density_floor_local-rho, 0.0);
        Real pgas_add = std::max(pressure_floor_local-pgas, 0.0);
        if (rho_add > 0.0 or pgas_add > 0.0)
        {
          Real u0_add, u1_add, u2_add, u3_add;
          pco->RaiseVectorCell(-alpha, 0.0, 0.0, 0.0, k, j, i, &u0_add, &u1_add,
              &u2_add, &u3_add);
          Real wgas_add = rho_add + gamma_adi/(gamma_adi-1.0) * pgas_add;
          cons(IDN,k,j,i) += rho_add * u0_add;
          cons(IEN,k,j,i) += wgas_add * u0_add * (-alpha) + pgas_add;
          fixed_(k,j,i) = true;
        }

        // Adjust primitive quantities to match
        if (rho_add > 0.0 or pgas_add > 0.0)
          fixed_(k,j,i) = ConservedToPrimitiveSingle(cons, prim_old, gamma_adi, bb_cc,
              g_, g_inv_, k, j, i, prim);
      }
    }

  // Readjust conserved values if fix was applied
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
    {
      bool fix_in_row = false;
      for (int i = is; i <= ie; ++i)
        if (fixed_(k,j,i))
          fix_in_row = true;
      if (fix_in_row)
      {
        pco->CellMetric(k, j, is, ie, g_, g_inv_);
        for (int i = is; i <= ie; ++i)
          if (fixed_(k,j,i))
          {
            PrimitiveToConservedSingle(prim, gamma_adi, bb_cc, g_, g_inv_, k, j, i,
                cons, pco);
            fixed_(k,j,i) = false;
          }
      }
    }

  // Apply limits, adjusting primitives as necessary
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
    {
      pco->CellMetric(k, j, is, ie, g_, g_inv_);
      for (int i = is; i <= ie; ++i)
      {
        // Extract primitives
        Real &rho = prim(IDN,k,j,i);
        Real &pgas = prim(IEN,k,j,i);
        Real &uu1 = prim(IVX,k,j,i);
        Real &uu2 = prim(IVY,k,j,i);
        Real &uu3 = prim(IVZ,k,j,i);

        // Calculate 4-velocity
        Real alpha = std::sqrt(-1.0/g_inv_(I00,i));
        Real tmp = g_(I11,i)*uu1*uu1 + 2.0*g_(I12,i)*uu1*uu2 + 2.0*g_(I13,i)*uu1*uu3
                 + g_(I22,i)*uu2*uu2 + 2.0*g_(I23,i)*uu2*uu3
                 + g_(I33,i)*uu3*uu3;
        Real gamma = std::sqrt(1.0 + tmp);
        Real u0 = gamma/alpha;
        Real u1 = uu1 - alpha * gamma * g_inv_(I01,i);
        Real u2 = uu2 - alpha * gamma * g_inv_(I02,i);
        Real u3 = uu3 - alpha * gamma * g_inv_(I03,i);
        Real u_0, u_1, u_2, u_3;
        pco->LowerVectorCell(u0, u1, u2, u3, k, j, i, &u_0, &u_1, &u_2, &u_3);

        // Extract cell-centered magnetic field
        const Real &bb1 = bb_cc(IB1,k,j,i);
        const Real &bb2 = bb_cc(IB2,k,j,i);
        const Real &bb3 = bb_cc(IB3,k,j,i);

        // Calculate magnetic pressure
        Real b0 = g_(I01,i)*u0*bb1 + g_(I02,i)*u0*bb2 + g_(I03,i)*u0*bb3
                + g_(I11,i)*u1*bb1 + g_(I12,i)*u1*bb2 + g_(I13,i)*u1*bb3
                + g_(I12,i)*u2*bb1 + g_(I22,i)*u2*bb2 + g_(I23,i)*u2*bb3
                + g_(I13,i)*u3*bb1 + g_(I23,i)*u3*bb2 + g_(I33,i)*u3*bb3;
        Real b1 = (bb1 + b0 * u1) / u0;
        Real b2 = (bb2 + b0 * u2) / u0;
        Real b3 = (bb3 + b0 * u3) / u0;
        Real b_0, b_1, b_2, b_3;
        pco->LowerVectorCell(b0, b1, b2, b3, k, j, i, &b_0, &b_1, &b_2, &b_3);
        Real b_sq = b0*b_0 + b1*b_1 + b2*b_2 + b3*b_3;
        Real pmag = 0.5*b_sq;

        // Calculate floors for density and pressure
        Real density_floor_local = density_floor_;
        if (rho_pow_ != 0.0)
          density_floor_local = std::max(density_floor_local,
              rho_min_*std::pow(pco->x1v(i),rho_pow_));
        density_floor_local = std::max(density_floor_local, rho_pmag_min_*pmag);
        Real pressure_floor_local = pressure_floor_;
        if (u_pow_ != 0.0)
          pressure_floor_local = std::max(pressure_floor_local,
              (gamma_-1.0)*u_min_*std::pow(pco->x1v(i),u_pow_));
        pressure_floor_local = std::max(pressure_floor_local,
            (gamma_-1.0)*u_pmag_min_*pmag);

        // Apply floors to density and pressure
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
        if (gamma > gamma_max_)
        {
          Real factor = std::sqrt((SQR(gamma_max_)-1.0) / (SQR(gamma)-1.0));
          uu1 *= factor;
          uu2 *= factor;
          uu3 *= factor;
          fixed_(k,j,i) = true;
        }
      }
    }

  // Finally ensure conserved quantities match primitives
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
    {
      bool fix_in_row = false;
      for (int i = is; i <= ie; ++i)
        if (fixed_(k,j,i))
          fix_in_row = true;
      if (fix_in_row)
      {
        pco->CellMetric(k, j, is, ie, g_, g_inv_);
        for (int i = is; i <= ie; ++i)
          if (fixed_(k,j,i))
          {
            PrimitiveToConservedSingle(prim, gamma_adi, bb_cc, g_, g_inv_, k, j, i,
                cons, pco);
            fixed_(k,j,i) = false;
          }
      }
    }
  return;
}

//--------------------------------------------------------------------------------------

// Function for finding primitives from conserved variables in a single cell
// Inputs:
//   cons: conserved quantities
//   prim_old: primitive quantities from previous half timestep
//   gamma_adi: ratio of specific heats
//   bb_cc: cell-centered magnetic field
//   g,gi: 1D arrays of metric covariant and contravariant coefficients
//   k,j,i: indices of cell
// Outputs:
//   prim: primitives
// Notes:
//   follows Noble et al. 2006, ApJ 641 626 (N)
//       writing uu for \tilde{u}
//       writing d for D
//       writing q for Q
//       writing qq for \tilde{Q}
//       writing vv for v
//       writing wgas_rel for W = \gamma^2 w
//       writing bb for B
//       writing bbb for \mathcal{B}
static bool ConservedToPrimitiveSingle(const AthenaArray<Real> &cons,
    const AthenaArray<Real> &prim_old, Real gamma_adi, const AthenaArray<Real> &bb_cc,
    const AthenaArray<Real> &g, const AthenaArray<Real> &gi, int k, int j, int i,
    AthenaArray<Real> &prim)
{
  // Parameters
  const Real max_wgas_rel = 1.0e8;
  const Real initial_guess_multiplier = 10.0;
  const int initial_guess_multiplications = 10;
  const int num_iterations = 5;

  // Prepare return value
  bool fixed = false;

  // Extract metric
  const Real &g_00 = g(I00,i), &g_01 = g(I01,i), &g_02 = g(I02,i), &g_03 = g(I03,i),
             &g_10 = g(I01,i), &g_11 = g(I11,i), &g_12 = g(I12,i), &g_13 = g(I13,i),
             &g_20 = g(I02,i), &g_21 = g(I12,i), &g_22 = g(I22,i), &g_23 = g(I23,i),
             &g_30 = g(I03,i), &g_31 = g(I13,i), &g_32 = g(I23,i), &g_33 = g(I33,i);
  const Real &g00 = gi(I00,i), &g01 = gi(I01,i), &g02 = gi(I02,i), &g03 = gi(I03,i),
             &g10 = gi(I01,i), &g11 = gi(I11,i), &g12 = gi(I12,i), &g13 = gi(I13,i),
             &g20 = gi(I02,i), &g21 = gi(I12,i), &g22 = gi(I22,i), &g23 = gi(I23,i),
             &g30 = gi(I03,i), &g31 = gi(I13,i), &g32 = gi(I23,i), &g33 = gi(I33,i);
  Real alpha = std::sqrt(-1.0/g00);

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

  // Extract cell-centered magnetic field
  const Real &bb1 = bb_cc(IB1,k,j,i);
  const Real &bb2 = bb_cc(IB2,k,j,i);
  const Real &bb3 = bb_cc(IB3,k,j,i);

  // Calculate variations on magnetic field
  Real bbb1 = alpha * bb1;                                                // (N 5)
  Real bbb2 = alpha * bb2;                                                // (N 5)
  Real bbb3 = alpha * bb3;                                                // (N 5)
  Real bbb_sq = g_11*bbb1*bbb1 + 2.0*g_12*bbb1*bbb2 + 2.0*g_13*bbb1*bbb3
              + g_22*bbb2*bbb2 + 2.0*g_23*bbb2*bbb3
              + g_33*bbb3*bbb3;
  Real q_bbb = q_1*bbb1 + q_2*bbb2 + q_3*bbb3;
  Real q_bbb_sq = SQR(q_bbb);

  // Extract old primitives
  const Real &rho_old = prim_old(IDN,k,j,i);
  const Real &pgas_old = prim_old(IEN,k,j,i);
  const Real &uu1_old = prim_old(IVX,k,j,i);
  const Real &uu2_old = prim_old(IVY,k,j,i);
  const Real &uu3_old = prim_old(IVZ,k,j,i);

  // Construct initial guess for relativistic gas enthalpy W
  Real tmp = g_11*SQR(uu1_old) + 2.0*g_12*uu1_old*uu2_old + 2.0*g_13*uu1_old*uu3_old
           + g_22*SQR(uu2_old) + 2.0*g_23*uu2_old*uu3_old
           + g_33*SQR(uu3_old);
  Real gamma_sq = 1.0 + tmp;
  Real wgas_rel_init = gamma_sq * (rho_old + gamma_adi/(gamma_adi-1.0) * pgas_old);
  for (int count = 0; count < initial_guess_multiplications; ++count)
  {
    Real vv_sq = (qq_sq*SQR(wgas_rel_init)
        + q_bbb_sq*(bbb_sq+2.0*wgas_rel_init))
        / SQR(wgas_rel_init*(bbb_sq+wgas_rel_init));  // (N 28)
    if (vv_sq >= 1.0)  // guess for W leads to superluminal velocities
      wgas_rel_init *= initial_guess_multiplier;
    else
      break;
  }

  // Apply Newton-Raphson method to find new W
  Real wgas_rel_new = wgas_rel_init;
  Real res_new = QNResidual(wgas_rel_new, d, q_n, qq_sq, bbb_sq, q_bbb_sq, gamma_adi);
  for (int n = 0; n < num_iterations; ++n)
  {
    Real wgas_rel_old = wgas_rel_new;
    Real res_old = res_new;
    Real derivative = QNResidualPrime(wgas_rel_old, d, qq_sq, bbb_sq, q_bbb_sq,
        gamma_adi);
    Real delta = -res_old / derivative;
    wgas_rel_new = wgas_rel_old + delta;
    res_new = QNResidual(wgas_rel_new, d, q_n, qq_sq, bbb_sq, q_bbb_sq, gamma_adi);
  }
  Real wgas_rel_true = wgas_rel_new;
  if (not (wgas_rel_true > 0.0 and wgas_rel_true < max_wgas_rel))
  {
    fixed = true;
    wgas_rel_true = wgas_rel_init;
  }
  Real vv_sq = (qq_sq*SQR(wgas_rel_true)
      + q_bbb_sq*(bbb_sq+2.0*wgas_rel_true))
      / SQR(wgas_rel_true*(bbb_sq+wgas_rel_true));  // (N 28)
  gamma_sq = 1.0/(1.0-vv_sq);
  Real gamma_rel = std::sqrt(gamma_sq);
  if (not std::isfinite(gamma_rel) or gamma_rel < 1.0)
  {
    fixed = true;
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

  // Set velocity (N 31)
  uu1 = gamma_rel/(wgas_rel_true+bbb_sq) * (qq1 + q_bbb/wgas_rel_true * bbb1);
  uu2 = gamma_rel/(wgas_rel_true+bbb_sq) * (qq2 + q_bbb/wgas_rel_true * bbb2);
  uu3 = gamma_rel/(wgas_rel_true+bbb_sq) * (qq3 + q_bbb/wgas_rel_true * bbb3);
  return fixed;
}

//--------------------------------------------------------------------------------------

// Function for converting all primitives to conserved variables
// Inputs:
//   prim: primitives
//   bb_cc: cell-centered magnetic field
//   pco: pointer to Coordinates
//   is,ie,js,je,ks,ke: index bounds of region to be updated
// Outputs:
//   cons: conserved variables
// Notes:
//   single-cell function exists for other purposes; call made to that function rather
//       than having duplicate code
void EquationOfState::PrimitiveToConserved(const AthenaArray<Real> &prim,
     const AthenaArray<Real> &bb_cc, AthenaArray<Real> &cons, Coordinates *pco,
     int is, int ie, int js, int je, int ks, int ke)
{
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
    {
      pco->CellMetric(k, j, is, ie, g_, g_inv_);
      #pragma simd
      for (int i = is; i <= ie; ++i)
        PrimitiveToConservedSingle(prim, gamma_, bb_cc, g_, g_inv_, k, j, i, cons, pco);
    }
  return;
}

//--------------------------------------------------------------------------------------

// Function for converting primitives to conserved variables in a single cell
// Inputs:
//   prim: 3D array of primitives
//   gamma_adi: ratio of specific heats
//   bb_cc: 3D array of cell-centered magnetic field
//   g,gi: 1D arrays of metric covariant and contravariant coefficients
//   k,j,i: indices of cell
//   pco: pointer to Coordinates
// Outputs:
//   cons: conserved variables set in desired cell
static void PrimitiveToConservedSingle(const AthenaArray<Real> &prim, Real gamma_adi,
    const AthenaArray<Real> &bb_cc, const AthenaArray<Real> &g,
    const AthenaArray<Real> &gi, int k, int j, int i, AthenaArray<Real> &cons,
    Coordinates *pco)
{
  // Extract primitives and magnetic fields
  const Real &rho = prim(IDN,k,j,i);
  const Real &pgas = prim(IEN,k,j,i);
  const Real &uu1 = prim(IVX,k,j,i);
  const Real &uu2 = prim(IVY,k,j,i);
  const Real &uu3 = prim(IVZ,k,j,i);
  const Real &bb1 = bb_cc(IB1,k,j,i);
  const Real &bb2 = bb_cc(IB2,k,j,i);
  const Real &bb3 = bb_cc(IB3,k,j,i);

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
  pco->LowerVectorCell(u0, u1, u2, u3, k, j, i, &u_0, &u_1, &u_2, &u_3);

  // Calculate 4-magnetic field
  Real b0 = g(I01,i)*u0*bb1 + g(I02,i)*u0*bb2 + g(I03,i)*u0*bb3
          + g(I11,i)*u1*bb1 + g(I12,i)*u1*bb2 + g(I13,i)*u1*bb3
          + g(I12,i)*u2*bb1 + g(I22,i)*u2*bb2 + g(I23,i)*u2*bb3
          + g(I13,i)*u3*bb1 + g(I23,i)*u3*bb2 + g(I33,i)*u3*bb3;
  Real b1 = (bb1 + b0 * u1) / u0;
  Real b2 = (bb2 + b0 * u2) / u0;
  Real b3 = (bb3 + b0 * u3) / u0;
  Real b_0, b_1, b_2, b_3;
  pco->LowerVectorCell(b0, b1, b2, b3, k, j, i, &b_0, &b_1, &b_2, &b_3);
  Real b_sq = b0*b_0 + b1*b_1 + b2*b_2 + b3*b_3;

  // Extract conserved quantities
  Real &rho_u0 = cons(IDN,k,j,i);
  Real &t0_0 = cons(IEN,k,j,i);
  Real &t0_1 = cons(IM1,k,j,i);
  Real &t0_2 = cons(IM2,k,j,i);
  Real &t0_3 = cons(IM3,k,j,i);

  // Set conserved quantities
  Real wtot = rho + gamma_adi/(gamma_adi-1.0) * pgas + b_sq;
  Real ptot = pgas + 0.5 * b_sq;
  rho_u0 = rho * u0;
  t0_0 = wtot * u0 * u_0 - b0 * b_0 + ptot;
  t0_1 = wtot * u0 * u_1 - b0 * b_1;
  t0_2 = wtot * u0 * u_2 - b0 * b_2;
  t0_3 = wtot * u0 * u_3 - b0 * b_3;
  return;
}

//--------------------------------------------------------------------------------------

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
//   almost same function as in adiabatic_mhd_sr.cpp
void EquationOfState::FastMagnetosonicSpeedsSR(const AthenaArray<Real> &prim,
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
    Real u[4];
    u[1] = prim(ivx,i);
    u[2] = prim(ivy,i);
    u[3] = prim(ivz,i);
    u[0] = std::sqrt(1.0 + SQR(u[1]) + SQR(u[2]) + SQR(u[3]));
    const Real &bbx = bbx_vals(i);
    const Real &bby = prim(IBY,i);
    const Real &bbz = prim(IBZ,i);

    // Calculate 3-velocity
    Real vx = u[1]/u[0];
    Real vy = u[2]/u[0];
    Real vz = u[3]/u[0];

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
      Real s = (s2 > 0.0) ? std::sqrt(s2) : 0.0;
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
      Real s = (s2 > 0.0) ? std::sqrt(s2) : 0.0;
      lambda_plus_no_bbx = (a1 >= 0.0) ? -2.0*a0/(a1+s) : (-a1+s)/2.0;
      lambda_minus_no_bbx = (a1 >= 0.0) ? (-a1-s)/2.0 : -2.0*a0/(a1-s);
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
        s2 = SQR(z0)/4.0 - b0;
        Real s = (s2 < 0.0) ? 0.0 : std::sqrt(s2);
        Real d0 = (b1 < 0) ? 0.5*z0+s : 0.5*z0-s;
        Real e0 = (b1 < 0) ? 0.5*z0-s : 0.5*z0+s;

        // Solve quadratic equations
        s2 = SQR(d1) - 4.0*d0;
        s = (s2 < 0.0) ? 0.0 : std::sqrt(s2);
        y1 = (d1 >= 0.0) ? (-d1-s)/2.0 : -2.0*d0/(d1-s);
        y2 = (d1 >= 0.0) ? -2.0*d0/(d1+s) : (-d1+s)/2.0;
        s2 = SQR(e1) - 4.0*e0;
        s = (s2 < 0.0) ? 0.0 : std::sqrt(s2);
        y3 = (e1 >= 0.0) ? (-e1-s)/2.0 : -2.0*e0/(e1-s);
        y4 = (e1 >= 0.0) ? -2.0*e0/(e1+s) : (-e1+s)/2.0;
      }

      // Calculate extremal original quartic roots
      lambda_minus = std::min(y1, y3) - a3/4.0;
      lambda_plus = std::max(y2, y4) - a3/4.0;

      // Ensure wavespeeds are not superluminal
      if (not std::isfinite(lambda_minus) or lambda_minus < -1.0)
        lambda_minus = -1.0;
      if (not std::isfinite(lambda_plus) or lambda_plus > 1.0)
        lambda_plus = 1.0;
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

//--------------------------------------------------------------------------------------

// Function for calculating relativistic fast wavespeeds in arbitrary coordinates
// Inputs:
//   rho_h: gas enthalpy
//   pgas: gas pressure
//   u0,u1: contravariant components of 4-velocity
//   b_sq: b_\mu b^\mu
//   g00,g01,g11: contravariant components of metric
// Outputs:
//   plambda_plus: value set to most positive wavespeed
//   plambda_minus: value set to most negative wavespeed
// Notes:
//   follows same general procedure as vchar() in phys.c in Harm
//   variables are named as though 1 is normal direction
void EquationOfState::FastMagnetosonicSpeedsGR(Real rho_h, Real pgas, Real u0, Real u1,
    Real b_sq, Real g00, Real g01, Real g11, Real *plambda_plus, Real *plambda_minus)
{
  // Parameters and constants
  const Real gamma_adi = gamma_;

  // Calculate comoving fast magnetosonic speed
  Real cs_sq = gamma_adi * pgas / rho_h;
  Real va_sq = b_sq / (b_sq + rho_h);
  Real cms_sq = cs_sq + va_sq - cs_sq * va_sq;

  // Set fast magnetosonic speeds in appropriate coordinates
  Real a = SQR(u0) - (g00 + SQR(u0)) * cms_sq;
  Real b = -2.0 * (u0*u1 - (g01 + u0*u1) * cms_sq);
  Real c = SQR(u1) - (g11 + SQR(u1)) * cms_sq;
  Real d = std::max(SQR(b) - 4.0*a*c, 0.0);
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

//--------------------------------------------------------------------------------------

// Function whose value vanishes for correct enthalpy
// Inputs:
//   w_guess: guess for total enthalpy W
//   d: D = alpha * rho * u^0
//   q_n: Q_mu n^mu = -alpha^2 g^{mu 0} T^0_mu
//   qq_sq: \tilde{Q}^2 = alpha^2 g^{mu nu} T^0_mu T^0_nu
//                          + alpha^4 (g^{0 mu} T^0_mu)^2
//   bbb_sq: \mathcal{B}_mu \mathcal{B}^mu = \alpha^2 g_{mu nu} B^mu B^nu
//   q_bbb_sq: (Q_mu \mathcal{B}^mu)^2 = (alpha^2 T^0_mu B^mu)^2
//   gamma_adi: ratio of specific heats
// Outputs:
//   returned value: calculated minus given value of q_n
// Notes:
//   follows Noble et al. 2006, ApJ 641 626 (N)
static Real QNResidual(Real w_guess, Real d, Real q_n, Real qq_sq, Real bbb_sq,
    Real q_bbb_sq, Real gamma_adi)
{
  Real v_norm_sq = (qq_sq*SQR(w_guess)
      + q_bbb_sq*(bbb_sq+2.0*w_guess))
      / SQR(w_guess*(bbb_sq+w_guess));                   // (N 28)
  Real gamma_sq = 1.0/(1.0 - v_norm_sq);
  Real pgas = (gamma_adi-1.0)/gamma_adi
      * (w_guess/gamma_sq - d/std::sqrt(gamma_sq));      // (N 32)
  Real q_dot_n_calc = -0.5*bbb_sq * (1.0+v_norm_sq)
      + q_bbb_sq / (2.0*SQR(w_guess)) - w_guess + pgas;  // (N 29)
  return q_dot_n_calc - q_n;
}

//--------------------------------------------------------------------------------------

// Derivative of QNResidual()
// Inputs:
//   w_guess: guess for total enthalpy W
//   d: D = alpha * rho * u^0
//   qq_sq: \tilde{Q}^2 = alpha^2 g^{mu nu} T^0_mu T^0_nu
//                          + alpha^4 (g^{0 mu} T^0_mu)^2
//   bbb_sq: \mathcal{B}_mu \mathcal{B}^mu = \alpha^2 g_{mu nu} B^mu B^nu
//   q_bbb_sq: (Q_mu \mathcal{B}^mu)^2 = (alpha^2 T^0_mu B^mu)^2
//   gamma_adi: ratio of specific heats
// Outputs:
//   returned value: derivative of calculated value of Q_mu n^mu
// Notes:
//   follows Noble et al. 2006, ApJ 641 626 (N)
static Real QNResidualPrime(Real w_guess, Real d, Real qq_sq, Real bbb_sq,
    Real q_bbb_sq, Real gamma_adi)
{
  Real w_sq = SQR(w_guess);
  Real w_cu = w_sq * w_guess;
  Real b_w_term = bbb_sq + w_guess;
  Real b_w_term_sq = SQR(b_w_term);
  Real b_w_term_cu = b_w_term_sq * b_w_term;
  Real v_norm_sq = (qq_sq*w_sq
      + q_bbb_sq*(bbb_sq+2.0*w_guess))
      / (w_sq*b_w_term_sq);                                            // (N 28)
  Real gamma_sq = 1.0/(1.0 - v_norm_sq);
  Real gamma_4 = SQR(gamma_sq);
  Real dv_norm_sq_dw_a = 3.0*w_sq + 3.0*bbb_sq*w_guess + SQR(bbb_sq);
  Real dv_norm_sq_dw_b = qq_sq*w_cu + q_bbb_sq*dv_norm_sq_dw_a;
  Real dv_norm_sq_dw = -2.0 * dv_norm_sq_dw_b / (w_cu * b_w_term_cu);
  Real dgamma_sq_dw = gamma_4 * dv_norm_sq_dw;
  Real dpgas_dw_a = 0.5 * d * std::sqrt(gamma_sq) - w_guess;
  Real dpgas_dw_b = gamma_sq + dpgas_dw_a * dgamma_sq_dw;
  Real dpgas_dw = (gamma_adi-1.0)/gamma_adi / gamma_4 * dpgas_dw_b;
  return -0.5*bbb_sq*dv_norm_sq_dw - q_bbb_sq/w_cu - 1.0
      + dpgas_dw;
}

//--------------------------------------------------------------------------------------

// Function for replacing primitive value in cell with average of neighbors
// Inputs:
//   prim: array of primitives
//   n: IDN, IEN, IVX, IVY, or IVZ
//   k,j,i: indices of cell
//   kl,ku,jl,ju,il,ju: limits of array
// Outputs:
//   prim(index,k,j,i) modified
// Notes
//   average will only include in-bounds, non-NAN neighbors
//   if no such neighbors exist, value will be unmodified
//   same function as in adiabatic_hydro_gr.cpp
//   TODO: decide if function should be kept/implemented
static void neighbor_average(AthenaArray<Real> &prim, int n, int k, int j, int i,
    int kl, int ku, int jl, int ju, int il, int iu)
{
  Real neighbor_sum = 0.0;
  int num_neighbors = 0;
  if (i > il and not std::isnan(prim(n,k,j,i-1)))
  {
    neighbor_sum += prim(n,k,j,i-1);
    num_neighbors += 1;
  }
  if (i < iu and not std::isnan(prim(n,k,j,i+1)))
  {
    neighbor_sum += prim(n,k,j,i+1);
    num_neighbors += 1;
  }
  if (j > jl and not std::isnan(prim(n,k,j-1,i)))
  {
    neighbor_sum += prim(n,k,j-1,i);
    num_neighbors += 1;
  }
  if (j < ju and not std::isnan(prim(n,k,j+1,i)))
  {
    neighbor_sum += prim(n,k,j+1,i);
    num_neighbors += 1;
  }
  if (k > kl and not std::isnan(prim(n,k-1,j,i)))
  {
    neighbor_sum += prim(n,k-1,j,i);
    num_neighbors += 1;
  }
  if (k < ku and not std::isnan(prim(n,k+1,j,i)))
  {
    neighbor_sum += prim(n,k+1,j,i);
    num_neighbors += 1;
  }
  if (num_neighbors > 0)
    prim(n,k,j,i) = neighbor_sum / num_neighbors;
  return;
}
