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
#include "../hydro/hydro.hpp"
#include "../hydro/eos/eos.hpp"

// Declarations
void InnerHydro(MeshBlock *pmb, AthenaArray<Real> &cons, int is, int ie, int js, int je,
    int ks, int ke);
void OuterHydro(MeshBlock *pmb, AthenaArray<Real> &cons, int is, int ie, int js, int je,
    int ks, int ke);
void InnerField(MeshBlock *pmb, InterfaceField &bb, int is, int ie, int js, int je,
    int ks, int ke);
void OuterField(MeshBlock *pmb, InterfaceField &bb, int is, int ie, int js, int je,
    int ks, int ke);
static void CalculatePrimitives(Real r, Real temp_min, Real temp_max, Real *prho,
    Real *ppgas, Real *put, Real *pur);
static Real TemperatureMin(Real r, Real t_min, Real t_max);
static Real TemperatureBisect(Real r, Real t_min, Real t_max);
static Real TemperatureResidual(Real t, Real m, Real n_adi, Real r, Real c1, Real c2);

// Global variables
static Real m;             // black hole mass
static Real n_adi, k_adi;  // hydro parameters
static Real r_crit;        // sonic point radius
static Real c1, c2;        // useful constants
static Real bsq_over_rho;  // b^2/rho at inner radius

// Function for setting initial conditions
// Inputs:
//   phyd: Hydro
//   pfld: Field (unused)
//   pin: parameters
// Outputs: (none)
// Notes:
//   sets primitive and conserved variables according to input primitives
//   references Hawley, Smarr, & Wilson 1984, ApJ 277 296 (HSW)
void Mesh::ProblemGenerator(Hydro *phyd, Field *pfld, ParameterInput *pin)
{
  // Parameters
  const Real temp_min = 1.0e-2;  // lesser temperature root must be greater than this
  const Real temp_max = 1.0e1;   // greater temperature root must be less than this

  // Prepare index bounds
  MeshBlock *pmb = phyd->pmy_block;
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
  m = pmb->pcoord->GetMass();

  // Get ratio of specific heats
  Real gamma_adi = phyd->pf_eos->GetGamma();
  n_adi = 1.0/(gamma_adi-1.0);

  // Read problem parameters
  k_adi = pin->GetReal("hydro", "k_adi");
  r_crit = pin->GetReal("problem", "r_crit");
  bsq_over_rho = 0.0;
  if (MAGNETIC_FIELDS_ENABLED)
    bsq_over_rho = pin->GetReal("problem", "bsq_over_rho");

  // Prepare scratch arrays
  AthenaArray<Real> g, gi;
  g.NewAthenaArray(NMETRIC, iu+1);
  gi.NewAthenaArray(NMETRIC, iu+1);

  // Prepare various constants for determining primitives
  Real u_crit_sq = m/(2.0*r_crit);                                          // (HSW 71)
  Real u_crit = -std::sqrt(u_crit_sq);
  Real t_crit = n_adi/(n_adi+1.0) * u_crit_sq/(1.0-(n_adi+3.0)*u_crit_sq);  // (HSW 74)
  c1 = std::pow(t_crit, n_adi) * u_crit * SQR(r_crit);                      // (HSW 68)
  c2 = SQR(1.0 + (n_adi+1.0) * t_crit) * (1.0 - 3.0*m/(2.0*r_crit));        // (HSW 69)

  // Initialize primitive values
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
    {
      pmb->pcoord->CellMetric(k, j, il, iu, g, gi);
      for (int i = il; i <= iu; ++i)
      {
        Real r, theta, phi;
        pmb->pcoord->GetBoyerLindquistCoordinates(pmb->pcoord->x1v(i),
            pmb->pcoord->x2v(j), pmb->pcoord->x3v(k), &r, &theta, &phi);
        Real rho, pgas, ut, ur;
        CalculatePrimitives(r, temp_min, temp_max, &rho, &pgas, &ut, &ur);
        Real u0, u1, u2, u3;
        pmb->pcoord->TransformVectorCell(ut, ur, 0.0, 0.0, k, j, i, &u0, &u1, &u2, &u3);
        Real uu1 = u1 - gi(I01,i)/gi(I00,i) * u0;
        Real uu2 = u2 - gi(I02,i)/gi(I00,i) * u0;
        Real uu3 = u3 - gi(I03,i)/gi(I00,i) * u0;
        phyd->w(IDN,k,j,i) = phyd->w1(IDN,k,j,i) = rho;
        phyd->w(IEN,k,j,i) = phyd->w1(IEN,k,j,i) = pgas;
        phyd->w(IM1,k,j,i) = phyd->w1(IM1,k,j,i) = uu1;
        phyd->w(IM2,k,j,i) = phyd->w1(IM2,k,j,i) = uu2;
        phyd->w(IM3,k,j,i) = phyd->w1(IM3,k,j,i) = uu3;
      }
    }

  // Initialize magnetic field
  if (MAGNETIC_FIELDS_ENABLED)
  {
    // Find normalization
    Real r, theta, phi;
    pmb->pcoord->GetBoyerLindquistCoordinates(pmb->pcoord->x1f(pmb->is),
        pmb->pcoord->x2v((jl+ju)/2), pmb->pcoord->x3v((kl+ku)/2), &r, &theta, &phi);
    Real rho, pgas, ut, ur;
    CalculatePrimitives(r, temp_min, temp_max, &rho, &pgas, &ut, &ur);
    Real bbr = 1.0/SQR(r);
    Real bt = 1.0/(1.0-2.0*m/r) * bbr * ur;
    Real br = (bbr + bt * ur) / ut;
    Real bsq = -(1.0-2.0*m/r) * SQR(bt) + 1.0/(1.0-2.0*m/r) * SQR(br);
    Real bsq_over_rho_actual = bsq/rho;
    Real normalization = std::sqrt(bsq_over_rho/bsq_over_rho_actual);

    // Set field
    for (int k = kl; k <= ku+1; ++k)
      for (int j = jl; j <= ju+1; ++j)
        for (int i = il; i <= iu+1; ++i)
        {
          // Set B^1
          if (j != ju+1 and k != ku+1)
          {
            pmb->pcoord->GetBoyerLindquistCoordinates(pmb->pcoord->x1f(i),
                pmb->pcoord->x2v(j), pmb->pcoord->x3v(k), &r, &theta, &phi);
            CalculatePrimitives(r, temp_min, temp_max, &rho, &pgas, &ut, &ur);
            bbr = normalization/SQR(r);
            bt = 1.0/(1.0-2.0*m/r) * bbr * ur;
            br = (bbr + bt * ur) / ut;
            Real u0, u1, u2, u3;
            pmb->pcoord->TransformVectorFace1(ut, ur, 0.0, 0.0, k, j, i, &u0, &u1, &u2,
                &u3);
            Real b0, b1, b2, b3;
            pmb->pcoord->TransformVectorFace1(bt, br, 0.0, 0.0, k, j, i, &b0, &b1, &b2,
                &b3);
            pfld->b.x1f(k,j,i) = b1 * u0 - b0 * u1;
          }

          // Set B^2
          if (i != iu+1 and k != ku+1)
          {
            pmb->pcoord->GetBoyerLindquistCoordinates(pmb->pcoord->x1v(i),
                pmb->pcoord->x2f(j), pmb->pcoord->x3v(k), &r, &theta, &phi);
            CalculatePrimitives(r, temp_min, temp_max, &rho, &pgas, &ut, &ur);
            bbr = normalization/SQR(r);
            bt = 1.0/(1.0-2.0*m/r) * bbr * ur;
            br = (bbr + bt * ur) / ut;
            Real u0, u1, u2, u3;
            pmb->pcoord->TransformVectorFace2(ut, ur, 0.0, 0.0, k, j, i, &u0, &u1, &u2,
                &u3);
            Real b0, b1, b2, b3;
            pmb->pcoord->TransformVectorFace2(bt, br, 0.0, 0.0, k, j, i, &b0, &b1, &b2,
                &b3);
            pfld->b.x2f(k,j,i) = b2 * u0 - b0 * u2;
          }

          // Set B^3
          if (i != iu+1 and j != ju+1)
          {
            pmb->pcoord->GetBoyerLindquistCoordinates(pmb->pcoord->x1v(i),
                pmb->pcoord->x2v(j), pmb->pcoord->x3f(k), &r, &theta, &phi);
            CalculatePrimitives(r, temp_min, temp_max, &rho, &pgas, &ut, &ur);
            bbr = normalization/SQR(r);
            bt = 1.0/(1.0-2.0*m/r) * bbr * ur;
            br = (bbr + bt * ur) / ut;
            Real u0, u1, u2, u3;
            pmb->pcoord->TransformVectorFace3(ut, ur, 0.0, 0.0, k, j, i, &u0, &u1, &u2,
                &u3);
            Real b0, b1, b2, b3;
            pmb->pcoord->TransformVectorFace3(bt, br, 0.0, 0.0, k, j, i, &b0, &b1, &b2,
                &b3);
            pfld->b.x3f(k,j,i) = b3 * u0 - b0 * u3;
          }
        }
  }

  // Calculate cell-centered magnetic field
  AthenaArray<Real> bb;
  bb.NewAthenaArray(3, ku+1, ju+1, iu+1);
  if (MAGNETIC_FIELDS_ENABLED)
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        for (int i = il; i <= iu; ++i)
        {
          // Extract face-centered magnetic field
          const Real &bbf1m = pfld->b.x1f(k,j,i);
          const Real &bbf1p = pfld->b.x1f(k,j,i+1);
          const Real &bbf2m = pfld->b.x2f(k,j,i);
          const Real &bbf2p = pfld->b.x2f(k,j+1,i);
          const Real &bbf3m = pfld->b.x3f(k,j,i);
          const Real &bbf3p = pfld->b.x3f(k+1,j,i);

          // Calculate cell-centered magnetic field
          Real tmp = (pmb->pcoord->x1v(i) - pmb->pcoord->x1f(i)) / pmb->pcoord->dx1f(i);
          bb(IB1,k,j,i) = (1.0-tmp) * bbf1m + tmp * bbf1p;
          tmp = (pmb->pcoord->x2v(j) - pmb->pcoord->x2f(j)) / pmb->pcoord->dx2f(j);
          bb(IB2,k,j,i) = (1.0-tmp) * bbf2m + tmp * bbf2p;
          tmp = (pmb->pcoord->x3v(k) - pmb->pcoord->x3f(k)) / pmb->pcoord->dx3f(k);
          bb(IB3,k,j,i) = (1.0-tmp) * bbf3m + tmp * bbf3p;
       }

  // Initialize conserved variables
  pmb->phydro->pf_eos->PrimitiveToConserved(kl, ku, jl, ju, il, iu, phyd->w, bb, phyd->u);

  // Free scratch arrays
  g.DeleteAthenaArray();
  gi.DeleteAthenaArray();
  bb.DeleteAthenaArray();

  // Enroll boundary functions
  pmb->pbval->EnrollHydroBoundaryFunction(inner_x1, InnerHydro);
  pmb->pbval->EnrollHydroBoundaryFunction(outer_x1, OuterHydro);
  if (MAGNETIC_FIELDS_ENABLED)
  {
    pmb->pbval->EnrollFieldBoundaryFunction(inner_x1, InnerField);
    pmb->pbval->EnrollFieldBoundaryFunction(outer_x1, OuterField);
  }
  return;
}

// Inner hydro boundary condition
// TODO: change when interface changes
void InnerHydro(MeshBlock *pmb, AthenaArray<Real> &cons, int is, int ie, int js, int je,
    int ks, int ke)
{
  return;
}

// Outer hydro boundary condition
// TODO: change when interface changes
void OuterHydro(MeshBlock *pmb, AthenaArray<Real> &cons, int is, int ie, int js, int je,
    int ks, int ke)
{
  return;
}

// Inner field boundary condition
// TODO: comment
void InnerField(MeshBlock *pmb, InterfaceField &bb, int is, int ie, int js, int je,
    int ks, int ke)
{
  return;
}

// Outer field boundary condition
// TODO: comment
void OuterField(MeshBlock *pmb, InterfaceField &bb, int is, int ie, int js, int je,
    int ks, int ke)
{
  return;
}

// Function for calculating primitives given radius
// Inputs:
//   r: Schwarzschild radius
//   temp_min,temp_max: bounds on temperature
// Outputs:
//   prho: value set to density
//   ppgas: value set to gas pressure
//   put: value set to u^t in Schwarzschild coordinates
//   pur: value set to u^r in Schwarzschild coordinates
// Notes:
//   references Hawley, Smarr, & Wilson 1984, ApJ 277 296 (HSW)
static void CalculatePrimitives(Real r, Real temp_min, Real temp_max, Real *prho,
    Real *ppgas, Real *put, Real *pur)
{
  // Calculate solution to (HSW 76)
  Real temp_neg_res = TemperatureMin(r, temp_min, temp_max);
  Real temp;
  if (r <= r_crit)  // use lesser of two roots
    temp = TemperatureBisect(r, temp_min, temp_neg_res);
  else  // user greater of two roots
    temp = TemperatureBisect(r, temp_neg_res, temp_max);

  // Calculate primitives
  Real rho = std::pow(temp/k_adi, n_adi);             // not same K as HSW
  Real pgas = temp * rho;
  Real ur = c1 / (SQR(r) * std::pow(temp, n_adi));    // (HSW 75)
  Real ut = std::sqrt(1.0/SQR(1.0-2.0*m/r) * SQR(ur)
      + 1.0/(1.0-2.0*m/r));

  // Set primitives
  *prho = rho;
  *ppgas = pgas;
  *put = ut;
  *pur = ur;
  return;
}

// Function for finding temperature at which residual is minimized
// Inputs:
//   r: Schwarzschild radius
//   t_min,t_max: bounds between which minimum must occur
// Outputs:
//   returned value: some temperature for which residual of (HSW 76) is negative
// Notes:
//   references Hawley, Smarr, & Wilson 1984, ApJ 277 296 (HSW)
//   performs golden section search (cf. Numerical Recipes, 3rd ed., 10.2)
static Real TemperatureMin(Real r, Real t_min, Real t_max)
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
//   r: Schwarzschild radius
//   t_min,t_max: bounds between which root must occur
// Outputs:
//   returned value: temperature that satisfies (HSW 76)
// Notes:
//   references Hawley, Smarr, & Wilson 1984, ApJ 277 296 (HSW)
//   performs bisection search
static Real TemperatureBisect(Real r, Real t_min, Real t_max)
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

// Function whose value vanishes for correct temperature
// Notes:
//   implements (76) from Hawley, Smarr, & Wilson 1984, ApJ 277 296
static Real TemperatureResidual(Real t, Real m, Real n_adi, Real r, Real c1, Real c2)
{
  return SQR(1.0 + (n_adi+1.0) * t)
      * (1.0 - 2.0*m/r + SQR(c1) / (SQR(SQR(r)) * std::pow(t, 2.0*n_adi))) - c2;
}
