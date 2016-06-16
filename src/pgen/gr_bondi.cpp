// General relativistic black hole accretion generator, spherically symmetric flows

// Primary header
#include "../mesh/mesh.hpp"

// C++ headers
#include <cmath>  // abs(), NAN, pow(), sqrt()

// Athena headers
#include "../athena.hpp"                   // macros, enums, FaceField
#include "../athena_arrays.hpp"            // AthenaArray
#include "../parameter_input.hpp"          // ParameterInput
#include "../bvals/bvals.hpp"              // BoundaryValues
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../eos/eos.hpp"                  // EquationOfState
#include "../field/field.hpp"              // Field
#include "../hydro/hydro.hpp"              // Hydro

// Declarations
void FixedBoundary(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, int is, int ie, int js, int je, int ks, int ke);
static void CalculatePrimitives(Real r, Real temp_min, Real temp_max, Real *prho,
    Real *ppgas, Real *put, Real *pur);
static Real TemperatureMin(Real r, Real t_min, Real t_max);
static Real TemperatureBisect(Real r, Real t_min, Real t_max);
static Real TemperatureResidual(Real t, Real r);

// Global variables
static Real m;             // black hole mass
static Real n_adi, k_adi;  // hydro parameters
static Real r_crit;        // sonic point radius
static Real c1, c2;        // useful constants
static Real bsq_over_rho;  // b^2/rho at inner radius

//--------------------------------------------------------------------------------------

// Function for initializing global mesh properties
// Inputs:
//   pin: input parameters (unused)
// Outputs: (none)
void Mesh::InitUserMeshData(ParameterInput *pin)
{
  // Enroll boundary functions
  EnrollUserBoundaryFunction(INNER_X1, FixedBoundary);
  EnrollUserBoundaryFunction(OUTER_X1, FixedBoundary);
  return;
}

//--------------------------------------------------------------------------------------

// Function for setting initial conditions
// Inputs:
//   pin: parameters
// Outputs: (none)
// Notes:
//   sets primitive and conserved variables according to input primitives
//   references Hawley, Smarr, & Wilson 1984, ApJ 277 296 (HSW)
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  // Parameters
  const Real temp_min = 1.0e-2;  // lesser temperature root must be greater than this
  const Real temp_max = 1.0e1;   // greater temperature root must be less than this

  // Prepare index bounds
  int il = is - NGHOST;
  int iu = ie + NGHOST;
  int jl = js;
  int ju = je;
  if (block_size.nx2 > 1)
  {
    jl -= (NGHOST);
    ju += (NGHOST);
  }
  int kl = ks;
  int ku = ke;
  if (block_size.nx3 > 1)
  {
    kl -= (NGHOST);
    ku += (NGHOST);
  }

  // Get mass of black hole
  m = pcoord->GetMass();

  // Get ratio of specific heats
  Real gamma_adi = peos->GetGamma();
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
      pcoord->CellMetric(k, j, il, iu, g, gi);
      for (int i = il; i <= iu; ++i)
      {
        Real r, theta, phi;
        pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i),
            pcoord->x2v(j), pcoord->x3v(k), &r, &theta, &phi);
        Real rho, pgas, ut, ur;
        CalculatePrimitives(r, temp_min, temp_max, &rho, &pgas, &ut, &ur);
        Real u0, u1, u2, u3;
        pcoord->TransformVectorCell(ut, ur, 0.0, 0.0, k, j, i, &u0, &u1, &u2, &u3);
        Real uu1 = u1 - gi(I01,i)/gi(I00,i) * u0;
        Real uu2 = u2 - gi(I02,i)/gi(I00,i) * u0;
        Real uu3 = u3 - gi(I03,i)/gi(I00,i) * u0;
        phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) = rho;
        phydro->w(IEN,k,j,i) = phydro->w1(IEN,k,j,i) = pgas;
        phydro->w(IM1,k,j,i) = phydro->w1(IM1,k,j,i) = uu1;
        phydro->w(IM2,k,j,i) = phydro->w1(IM2,k,j,i) = uu2;
        phydro->w(IM3,k,j,i) = phydro->w1(IM3,k,j,i) = uu3;
      }
    }

  // Initialize magnetic field
  if (MAGNETIC_FIELDS_ENABLED)
  {
    // Find normalization
    Real r, theta, phi;
    pcoord->GetBoyerLindquistCoordinates(pcoord->x1f(is),
        pcoord->x2v((jl+ju)/2), pcoord->x3v((kl+ku)/2), &r, &theta, &phi);
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
            pcoord->GetBoyerLindquistCoordinates(pcoord->x1f(i),
                pcoord->x2v(j), pcoord->x3v(k), &r, &theta, &phi);
            CalculatePrimitives(r, temp_min, temp_max, &rho, &pgas, &ut, &ur);
            bbr = normalization/SQR(r);
            bt = 1.0/(1.0-2.0*m/r) * bbr * ur;
            br = (bbr + bt * ur) / ut;
            Real u0, u1, u2, u3;
            pcoord->TransformVectorFace1(ut, ur, 0.0, 0.0, k, j, i, &u0, &u1, &u2,
                &u3);
            Real b0, b1, b2, b3;
            pcoord->TransformVectorFace1(bt, br, 0.0, 0.0, k, j, i, &b0, &b1, &b2,
                &b3);
            pfield->b.x1f(k,j,i) = b1 * u0 - b0 * u1;
          }

          // Set B^2
          if (i != iu+1 and k != ku+1)
          {
            pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i),
                pcoord->x2f(j), pcoord->x3v(k), &r, &theta, &phi);
            CalculatePrimitives(r, temp_min, temp_max, &rho, &pgas, &ut, &ur);
            bbr = normalization/SQR(r);
            bt = 1.0/(1.0-2.0*m/r) * bbr * ur;
            br = (bbr + bt * ur) / ut;
            Real u0, u1, u2, u3;
            pcoord->TransformVectorFace2(ut, ur, 0.0, 0.0, k, j, i, &u0, &u1, &u2,
                &u3);
            Real b0, b1, b2, b3;
            pcoord->TransformVectorFace2(bt, br, 0.0, 0.0, k, j, i, &b0, &b1, &b2,
                &b3);
            pfield->b.x2f(k,j,i) = b2 * u0 - b0 * u2;
          }

          // Set B^3
          if (i != iu+1 and j != ju+1)
          {
            pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i),
                pcoord->x2v(j), pcoord->x3f(k), &r, &theta, &phi);
            CalculatePrimitives(r, temp_min, temp_max, &rho, &pgas, &ut, &ur);
            bbr = normalization/SQR(r);
            bt = 1.0/(1.0-2.0*m/r) * bbr * ur;
            br = (bbr + bt * ur) / ut;
            Real u0, u1, u2, u3;
            pcoord->TransformVectorFace3(ut, ur, 0.0, 0.0, k, j, i, &u0, &u1, &u2,
                &u3);
            Real b0, b1, b2, b3;
            pcoord->TransformVectorFace3(bt, br, 0.0, 0.0, k, j, i, &b0, &b1, &b2,
                &b3);
            pfield->b.x3f(k,j,i) = b3 * u0 - b0 * u3;
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
          const Real &bbf1m = pfield->b.x1f(k,j,i);
          const Real &bbf1p = pfield->b.x1f(k,j,i+1);
          const Real &bbf2m = pfield->b.x2f(k,j,i);
          const Real &bbf2p = pfield->b.x2f(k,j+1,i);
          const Real &bbf3m = pfield->b.x3f(k,j,i);
          const Real &bbf3p = pfield->b.x3f(k+1,j,i);

          // Calculate cell-centered magnetic field
          Real tmp = (pcoord->x1v(i) - pcoord->x1f(i)) / pcoord->dx1f(i);
          bb(IB1,k,j,i) = (1.0-tmp) * bbf1m + tmp * bbf1p;
          tmp = (pcoord->x2v(j) - pcoord->x2f(j)) / pcoord->dx2f(j);
          bb(IB2,k,j,i) = (1.0-tmp) * bbf2m + tmp * bbf2p;
          tmp = (pcoord->x3v(k) - pcoord->x3f(k)) / pcoord->dx3f(k);
          bb(IB3,k,j,i) = (1.0-tmp) * bbf3m + tmp * bbf3p;
       }

  // Initialize conserved variables
  peos->PrimitiveToConserved(phydro->w, bb, phydro->u, pcoord, il, iu, jl, ju, kl, ku);

  // Free scratch arrays
  g.DeleteAthenaArray();
  gi.DeleteAthenaArray();
  bb.DeleteAthenaArray();

  return;
}

//--------------------------------------------------------------------------------------

// Fixed boundary condition
// Inputs:
//   pmb: pointer to MeshBlock
//   pcoord: pointer to Coordinates
//   is,ie,js,je,ks,ke: indices demarkating active region
// Outputs:
//   prim: primitives set in ghost zones
//   bb: face-centered magnetic field set in ghost zones
// Notes:
//   does nothing
void FixedBoundary(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, int is, int ie, int js, int je, int ks, int ke)
{
  return;
}

//--------------------------------------------------------------------------------------

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

//--------------------------------------------------------------------------------------

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
  Real res_mid = TemperatureResidual(t_mid, r);

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
      Real res_new = TemperatureResidual(t_new, r);
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
      Real res_new = TemperatureResidual(t_new, r);
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

//--------------------------------------------------------------------------------------

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
  Real res_min = TemperatureResidual(t_min, r);
  Real res_max = TemperatureResidual(t_max, r);
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
    Real res_mid = TemperatureResidual(t_mid, r);
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

//--------------------------------------------------------------------------------------

// Function whose value vanishes for correct temperature
// Inputs:
//   t: temperature
//   r: Schwarzschild radius
// Outputs:
//   returned value: residual that should vanish for correct temperature
// Notes:
//   implements (76) from Hawley, Smarr, & Wilson 1984, ApJ 277 296
static Real TemperatureResidual(Real t, Real r);
{
  return SQR(1.0 + (n_adi+1.0) * t)
      * (1.0 - 2.0*m/r + SQR(c1) / (SQR(SQR(r)) * std::pow(t, 2.0*n_adi))) - c2;
}
