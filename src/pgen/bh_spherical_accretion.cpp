// General relativistic black hole accretion generator, spherically symmetric flows

// Primary header
#include "../mesh.hpp"

// C++ headers
#include <cmath>  // abs(), NAN, pow(), sqrt()

// Athena headers
#include "../athena.hpp"                   // macros, enums, Real
#include "../athena_arrays.hpp"            // AthenaArray
#include "../parameter_input.hpp"          // ParameterInput
#include "../bvals/bvals.hpp"              // BoundaryValues, InterfaceField
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../field/field.hpp"              // Field
#include "../fluid/fluid.hpp"              // Fluid
#include "../fluid/eos/eos.hpp"            // FluidEqnOfState

// Declarations
void FixedInner(MeshBlock *pmb, AthenaArray<Real> &cons,
                int is, int ie, int js, int je, int ks, int ke);
void FixedOuter(MeshBlock *pmb, AthenaArray<Real> &cons,
                int is, int ie, int js, int je, int ks, int ke);
static void set_state(Real rho, Real pgas, Real uu1, Real uu2, Real uu3, int k, int j,
    int i, AthenaArray<Real> &prim, AthenaArray<Real> &prim_half);
static Real TemperatureResidual(Real t, Real m, Real n_adi, Real r, Real c1, Real c2);
static Real TemperatureMin(Real m, Real n_adi, Real r, Real c1, Real c2, Real t_min,
    Real t_max);
static Real TemperatureBisect(Real t_min, Real t_max, Real m, Real n_adi, Real r,
    Real c1, Real c2);

// Global variables
static Real d_inner_1, e_inner_1, m1_inner_1, m2_inner_1, m3_inner_1;
static Real d_inner_2, e_inner_2, m1_inner_2, m2_inner_2, m3_inner_2;
static Real d_outer_1, e_outer_1, m1_outer_1, m2_outer_1, m3_outer_1;
static Real d_outer_2, e_outer_2, m1_outer_2, m2_outer_2, m3_outer_2;

// Function for setting initial conditions
// Inputs:
//   pfl: Fluid
//   pfd: Field (unused)
//   pin: parameters
// Outputs: (none)
// Notes:
//   sets primitive and conserved variables according to input primitives
//   references Hawley, Smarr, & Wilson 1984, ApJ 277 296 (HSW)
void Mesh::ProblemGenerator(Fluid *pfl, Field *pfd, ParameterInput *pin)
{
  // Parameters
  const Real t_min = 1.0e-2;  // lesser temperature root must be greater than this
  const Real t_max = 1.0e1;   // greater temperature root must be less than this

  // Prepare index bounds
  MeshBlock *pmb = pfl->pmy_block;
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

  // Get mass of black hole
  const Real m = pmb->pcoord->GetMass();

  // Get ratio of specific heats
  const Real gamma_adi = pfl->pf_eos->GetGamma();
  const Real gamma_adi_red = gamma_adi / (gamma_adi - 1.0);
  const Real n_adi = 1.0/(gamma_adi-1.0);

  // Read problem parameters
  const Real k_adi = pin->GetReal("fluid", "k_adi");
  const Real r_crit = pin->GetReal("problem", "r_crit");

  // Read initial magnetic field
  Real b1_flux = 0.0, b2_flux = 0.0, b3_flux = 0.0;
  if (MAGNETIC_FIELDS_ENABLED)
  {
    b1_flux = pin->GetReal("problem", "b1_flux");
    b2_flux = pin->GetReal("problem", "b2_flux");
    b3_flux = pin->GetReal("problem", "b3_flux");
  }

  // Prepare temporary arrays
  AthenaArray<Real> a1, a2m, a2p, a3m, a3p, b, g, gi;
  a1.NewAthenaArray(iu+1);
  a2m.NewAthenaArray(iu);
  a2p.NewAthenaArray(iu);
  a3m.NewAthenaArray(iu);
  a3p.NewAthenaArray(iu);
  b.NewAthenaArray(3, ku+2, ju+2, iu+2);
  g.NewAthenaArray(NMETRIC, iu+1);
  gi.NewAthenaArray(NMETRIC, iu+1);

  // Initialize magnetic field
  if (MAGNETIC_FIELDS_ENABLED)
    for (int k = kl; k <= ku+1; ++k)
    {
      Real interp_param_k = (pmb->pcoord->x3v(k) - pmb->pcoord->x3f(k))
          / pmb->pcoord->dx3f(k);
      for (int j = jl; j <= ju+1; ++j)
      {
        Real interp_param_j = (pmb->pcoord->x2v(j) - pmb->pcoord->x2f(j))
            / pmb->pcoord->dx2f(j);
        pmb->pcoord->Face1Area(k, j, il, iu+1, a1);
        pmb->pcoord->Face2Area(k, j, il, iu, a2m);
        pmb->pcoord->Face2Area(k, j+1, il, iu, a2p);
        pmb->pcoord->Face3Area(k, j, il, iu, a3m);
        pmb->pcoord->Face3Area(k+1, j, il, iu, a3p);
        for (int i = il; i <= iu+1; ++i)
        {
          Real interp_param_i = (pmb->pcoord->x1v(i) - pmb->pcoord->x1f(i))
              / pmb->pcoord->dx1f(i);
          Real b1m = b1_flux / a1(i);
          Real b1p = b1_flux / a1(i+1);
          Real b2m = b2_flux / a2m(i);
          Real b2p = b3_flux / a2p(i);
          Real b3m = b3_flux / a3m(i);
          Real b3p = b3_flux / a3p(i);
          b(IB1,k,j,i) = (1.0-interp_param_i) * b1m + interp_param_i * b1p;
          b(IB2,k,j,i) = (1.0-interp_param_j) * b2m + interp_param_j * b2p;
          b(IB3,k,j,i) = (1.0-interp_param_k) * b3m + interp_param_k * b3p;
          pfd->b.x1f(k,j,i) = b1m;
          if (i == iu)
            pfd->b.x1f(k,j,i+1) = b1p;
          pfd->b.x2f(k,j,i) = b2m;
          if (j == ju)
            pfd->b.x2f(k,j+1,i) = b2p;
          pfd->b.x3f(k,j,i) = b3m;
          if (k == ku)
            pfd->b.x3f(k+1,j,i) = b3p;
        }
      }
    }

  // Prepare various constants for determining primitives
  Real u_crit_sq = m / (2.0*r_crit);                         // (HSW 71)
  Real u_crit = -std::sqrt(u_crit_sq);
  Real t_crit = n_adi/(n_adi+1.0) * u_crit_sq
      / (1.0 - (n_adi+3.0) * u_crit_sq);
  Real c1 = std::pow(t_crit, n_adi) * u_crit * SQR(r_crit);  // (cf. HSW 69)
  Real c2 = SQR(1.0 + (n_adi+1.0) * t_crit)
      * (1.0 - 3.0*m / (2.0*r_crit));                        // (cf. HSW 68)

  // Initialize primitives
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
    {
      pmb->pcoord->CellMetric(k, j, il, iu, g, gi);
      for (int i = il; i <= iu; ++i)
      {
        // Get radius
        Real r = pmb->pcoord->x1v(i);

        // Calculate solution to (HSW 76)
        Real t_neg_res = TemperatureMin(m, n_adi, r, c1, c2, t_min, t_max);
        Real temperature;
        if (r <= r_crit)  // use lesser of two roots
          temperature = TemperatureBisect(t_min, t_neg_res, m, n_adi, r, c1, c2);
        else  // user greater of two roots
          temperature = TemperatureBisect(t_neg_res, t_max, m, n_adi, r, c1, c2);

        // Calculate primitives
        Real rho = std::pow(temperature/k_adi, n_adi);           // not same K as HSW
        Real pgas = temperature * rho;
        Real u1 = c1 / (SQR(r) * std::pow(temperature, n_adi));  // (HSW 75)
        Real u0 = std::sqrt(1.0/SQR(1.0 - 2.0*m/r) * SQR(u1)
            + 1.0/(1.0 - 2.0*m/r));
        Real u2 = 0.0;
        Real u3 = 0.0;
        Real uu1 = u1 - gi(I01,i)/gi(I00,i) * u0;
        Real uu2 = u2 - gi(I02,i)/gi(I00,i) * u0;
        Real uu3 = u3 - gi(I03,i)/gi(I00,i) * u0;
        set_state(rho, pgas, uu1, uu2, uu3, k, j, i, pfl->w,
            pfl->w1);
      }
    }

  // Initialize conserved variables
  pmb->pfluid->pf_eos->PrimitiveToConserved(pfl->w, b, pfl->u);  

  // Delete temporary arrays
  a1.DeleteAthenaArray();
  a2m.DeleteAthenaArray();
  a2p.DeleteAthenaArray();
  a3m.DeleteAthenaArray();
  a3p.DeleteAthenaArray();
  b.DeleteAthenaArray();
  g.DeleteAthenaArray();
  gi.DeleteAthenaArray();

  // Save inner boundary state
  d_inner_1 = pfl->u(IDN,pmb->ks,pmb->js,pmb->is-1);
  e_inner_1 = pfl->u(IEN,pmb->ks,pmb->js,pmb->is-1);
  m1_inner_1 = pfl->u(IM1,pmb->ks,pmb->js,pmb->is-1);
  m2_inner_1 = pfl->u(IM2,pmb->ks,pmb->js,pmb->is-1);
  m3_inner_1 = pfl->u(IM3,pmb->ks,pmb->js,pmb->is-1);
  d_inner_2 = pfl->u(IDN,pmb->ks,pmb->js,pmb->is-2);
  e_inner_2 = pfl->u(IEN,pmb->ks,pmb->js,pmb->is-2);
  m1_inner_2 = pfl->u(IM1,pmb->ks,pmb->js,pmb->is-2);
  m2_inner_2 = pfl->u(IM2,pmb->ks,pmb->js,pmb->is-2);
  m3_inner_2 = pfl->u(IM3,pmb->ks,pmb->js,pmb->is-2);

  // Save outer boundary state
  d_outer_1 = pfl->u(IDN,pmb->ks,pmb->js,pmb->ie+1);
  e_outer_1 = pfl->u(IEN,pmb->ks,pmb->js,pmb->ie+1);
  m1_outer_1 = pfl->u(IM1,pmb->ks,pmb->js,pmb->ie+1);
  m2_outer_1 = pfl->u(IM2,pmb->ks,pmb->js,pmb->ie+1);
  m3_outer_1 = pfl->u(IM3,pmb->ks,pmb->js,pmb->ie+1);
  d_outer_2 = pfl->u(IDN,pmb->ks,pmb->js,pmb->ie+2);
  e_outer_2 = pfl->u(IEN,pmb->ks,pmb->js,pmb->ie+2);
  m1_outer_2 = pfl->u(IM1,pmb->ks,pmb->js,pmb->ie+2);
  m2_outer_2 = pfl->u(IM2,pmb->ks,pmb->js,pmb->ie+2);
  m3_outer_2 = pfl->u(IM3,pmb->ks,pmb->js,pmb->ie+2);

  // Enroll boundary functions
  pmb->pbval->EnrollFluidBoundaryFunction(inner_x1, FixedInner);
  pmb->pbval->EnrollFluidBoundaryFunction(outer_x1, FixedOuter);
  return;
}

// Inner boundary condition
void FixedInner(MeshBlock *pmb, AthenaArray<Real> &cons,
                int is, int ie, int js, int je, int ks, int ke)
{
  // Set conserved values
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
    {
      cons(IDN,k,j,is-1) = d_inner_1;
      cons(IEN,k,j,is-1) = e_inner_1;
      cons(IM1,k,j,is-1) = m1_inner_1;
      cons(IM2,k,j,is-1) = m2_inner_1;
      cons(IM3,k,j,is-1) = m3_inner_1;
      cons(IDN,k,j,is-2) = d_inner_2;
      cons(IEN,k,j,is-2) = e_inner_2;
      cons(IM1,k,j,is-2) = m1_inner_2;
      cons(IM2,k,j,is-2) = m2_inner_2;
      cons(IM3,k,j,is-2) = m3_inner_2;
    }
  return;
}

// Outer boundary condition
void FixedOuter(MeshBlock *pmb, AthenaArray<Real> &cons,
                int is, int ie, int js, int je, int ks, int ke)
{
  // Set conserved values
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
    {
      cons(IDN,k,j,ie+1) = d_outer_1;
      cons(IEN,k,j,ie+1) = e_outer_1;
      cons(IM1,k,j,ie+1) = m1_outer_1;
      cons(IM2,k,j,ie+1) = m2_outer_1;
      cons(IM3,k,j,ie+1) = m3_outer_1;
      cons(IDN,k,j,ie+2) = d_outer_2;
      cons(IEN,k,j,ie+2) = e_outer_2;
      cons(IM1,k,j,ie+2) = m1_outer_2;
      cons(IM2,k,j,ie+2) = m2_outer_2;
      cons(IM3,k,j,ie+2) = m3_outer_2;
    }
  return;
}

// Function for setting conserved variables in a cell given the primitives
// Inputs:
//   rho: density
//   pgas: gas pressure
//   uu1,uu2,uu3: projected 4-velocity compontents \tilde{u}^i
//   k,j,i: indices for cell
// Outputs:
//   prim,prim_half: primitives in cell set
static void set_state(Real rho, Real pgas, Real uu1, Real uu2, Real uu3, int k, int j,
    int i, AthenaArray<Real> &prim, AthenaArray<Real> &prim_half)
{
  prim(IDN,k,j,i) = prim_half(IDN,k,j,i) = rho;
  prim(IEN,k,j,i) = prim_half(IEN,k,j,i) = pgas;
  prim(IM1,k,j,i) = prim_half(IM1,k,j,i) = uu1;
  prim(IM2,k,j,i) = prim_half(IM2,k,j,i) = uu2;
  prim(IM3,k,j,i) = prim_half(IM3,k,j,i) = uu3;
  return;
}

// Function whose value vanishes for correct temperature
// Notes:
//   implements (76) from Hawley, Smarr, & Wilson 1984, ApJ 277 296
static Real TemperatureResidual(Real t, Real m, Real n_adi, Real r, Real c1, Real c2)
{
  return SQR(1.0 + (n_adi+1.0) * t)
      * (1.0 - 2.0*m/r + SQR(c1) / (SQR(SQR(r)) * std::pow(t, 2.0*n_adi))) - c2;
}

// Function for finding temperature at which residual is minimized
// Inputs:
//   m: black hole mass
//   n_adi: polytropic index n = 1/(1-\Gamma)
//   r: Schwarzschild radius
//   c1,c2: constants as defined by (HSW 68,69)
//   t_min,t_max: bounds between which minimum must occur
// Outputs:
//   returned value: some temperature for which residual of (HSW 76) is negative
// Notes:
//   references Hawley, Smarr, & Wilson 1984, ApJ 277 296 (HSW)
//   performs golden section search (cf. Numerical Recipes, 3rd ed., 10.2)
static Real TemperatureMin(Real m, Real n_adi, Real r, Real c1, Real c2, Real t_min,
    Real t_max)
{
  // Parameters
  const Real ratio = 0.3819660112501051;  // (3+\sqrt{5})/2
  const int max_iterations = 30;          // maximum number of iterations

  // Initialize values
  Real t_mid = t_min + ratio * (t_max - t_min);
  Real res_mid = TemperatureResidual(t_mid, m, n_adi, r, c1, c2);

  // Apply golden section method
  bool larger_to_right = true;  // flag indicating larger subinterval is on right
  for (int n = 0; n < max_iterations; ++n)
  {
    if (res_mid < 0.0)
      return t_mid;
    Real t_new;
    if (larger_to_right)
    {
      t_new = t_mid + ratio * (t_max - t_mid);
      Real res_new = TemperatureResidual(t_new, m, n_adi, r, c1, c2);
      if (res_new < res_mid)
      {
        t_min = t_mid;
        t_mid = t_new;
        res_mid = res_new;
      }
      else
      {
        t_max = t_new;
        larger_to_right = false;
      }
    }
    else
    {
      t_new = t_mid - ratio * (t_mid - t_min);
      Real res_new = TemperatureResidual(t_new, m, n_adi, r, c1, c2);
      if (res_new < res_mid)
      {
        t_max = t_mid;
        t_mid = t_new;
        res_mid = res_new;
      }
      else
      {
        t_min = t_new;
        larger_to_right = true;
      }
    }
  }
  return NAN;
}

// Bisection root finder
// Inputs:
//   t_min,t_max: bounds between which root must occur
//   m: black hole mass
//   n_adi: polytropic index n = 1/(1-\Gamma)
//   r: Schwarzschild radius
//   c1,c2: constants as defined by (HSW 68,69)
// Outputs:
//   returned value: temperature that satisfies (HSW 76)
// Notes:
//   references Hawley, Smarr, & Wilson 1984, ApJ 277 296 (HSW)
//   performs bisection search
static Real TemperatureBisect(Real t_min, Real t_max, Real m, Real n_adi, Real r,
    Real c1, Real c2)
{
  // Parameters
  const int max_iterations = 20;
  const Real tol_residual = 1.0e-6;
  const Real tol_temperature = 1.0e-6;

  // Find initial residuals
  Real res_min = TemperatureResidual(t_min, m, n_adi, r, c1, c2);
  Real res_max = TemperatureResidual(t_max, m, n_adi, r, c1, c2);
  if (std::abs(res_min) < tol_residual)
    return t_min;
  if (std::abs(res_max) < tol_residual)
    return t_max;
  if ((res_min < 0.0 and res_max < 0.0) or (res_min > 0.0 and res_max > 0.0))
    return NAN;

  // Iterate to find root
  Real t_mid;
  for (int i = 0; i < max_iterations; ++i)
  {
    t_mid = (t_min + t_max) / 2.0;
    if (t_max - t_min < tol_temperature)
      return t_mid;
    Real res_mid = TemperatureResidual(t_mid, m, n_adi, r, c1, c2);
    if (std::abs(res_mid) < tol_residual)
      return t_mid;
    if ((res_mid < 0.0 and res_min < 0.0) or (res_mid > 0.0 and res_min > 0.0))
    {
      t_min = t_mid;
      res_min = res_mid;
    }
    else
    {
      t_max = t_mid;
      res_max = res_mid;
    }
  }
  return t_mid;
}
