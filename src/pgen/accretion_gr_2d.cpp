// General relativistic black hole accretion generator, azimuthally symmetric flows

// Primary header
#include "../mesh.hpp"

// C++ headers
#include <cmath>  // pow(), sin(), sqrt()

// Athena headers
#include "../athena.hpp"                   // enums, Real
#include "../athena_arrays.hpp"            // AthenaArray
#include "../coordinates/coordinates.hpp"  // PrimToCons()
#include "../fluid/fluid.hpp"              // Fluid
#include "../fluid/bvals/bvals.hpp"        // EnrollBoundaryFunction()
#include "../fluid/eos/eos.hpp"            // GetGamma()
#include "../parameter_input.hpp"          // ParameterInput

// Declarations
void FixedInner(MeshBlock *pmb, AthenaArray<Real> &cons,
                int is, int ie, int js, int je, int ks, int ke);
void FixedOuter(MeshBlock *pmb, AthenaArray<Real> &cons,
                int is, int ie, int js, int je, int ks, int ke);
void FixedTop(MeshBlock *pmb, AthenaArray<Real> &cons,
              int is, int ie, int js, int je, int ks, int ke);
void FixedBottom(MeshBlock *pmb, AthenaArray<Real> &cons,
                 int is, int ie, int js, int je, int ks, int ke);
static void set_state(AthenaArray<Real> &prim, AthenaArray<Real> &prim_half, int i,
    int j, int k, Real rho, Real pgas, Real vx, Real vy, Real vz);

// Global variables
MeshBlock *pb;
static Real gamma_adi, k_adi, l, rho_floor;
// TODO: put in a better place
static Real M = 1.0;

// Function for setting initial conditions
// Inputs:
//   pfl: Fluid
//   pfd: Field (unused)
//   pin: parameters
// Outputs: (none)
// Notes:
//   sets primitive and conserved variables according to input primitives
//   calculates fat disk from Hawley, Smarr, & Wilson 1984, ApJ 277 296 (HSW)
//   TODO: assumes Schwarzschild - is this okay?
void Mesh::ProblemGenerator(Fluid *pfl, Field *pfd, ParameterInput *pin)
{
  // Prepare index bounds
  MeshBlock *pb = pfl->pmy_block;
  int il = pb->is - NGHOST;
  int iu = pb->ie + NGHOST;
  int jl = pb->js;
  int ju = pb->je;
  if (pb->block_size.nx2 > 1)
  {
    jl -= (NGHOST);
    ju += (NGHOST);
  }
  int kl = pb->ks;
  int ku = pb->ke;
  if (pb->block_size.nx3 > 1)
  {
    kl -= (NGHOST);
    ku += (NGHOST);
  }

  // Get ratio of specific heats
  gamma_adi = pf_eos->GetGamma();

  // Read other properties
  k_adi = pin->GetReal("problem", "k");
  l = pin->GetReal("problem", "l");
  rho_floor = pin->GetReal("problem", "rho_floor");

  // TODO: read and set mass

  // Initialize the material according to Hawley, Smarr, & Wilson
  AthenaArray<Real> g, g_inv;
  g.NewAthenaArray(NMETRIC,iu+1);
  g_inv.NewAthenaArray(NMETRIC,iu+1);
  for (int k = kl; k <= ku; k++)
    for (int j = jl; j <= ju; j++)
    {
      pb->pcoord->CellMetric(k, j, g, g_inv);
      Real theta = pb->x2v(j);
      for (int i = il; i <= iu; i++)
      {
        Real r = pb->x1v(i);
        Real neg_u_0_inv_sq = -g_inv(I00,i) - g_inv(I33,i) * l*l;  // (HSW 95a)
        Real u_0, rho;
        if (neg_u_0_inv_sq <= 0.0)
        {
          u_0 = -1.0;
          rho = rho_floor;
        }
        else
        {
          u_0 = -1.0 / std::sqrt(neg_u_0_inv_sq);
          if (u_0 <= -1.0)
          {
            u_0 = -1.0;
            rho = rho_floor;
          }
          else
          {
            Real epsilon = -1.0/gamma_adi * (1.0/u_0 + 1.0);  // (HSW 94a)
            rho = std::pow((gamma_adi-1.0) * epsilon / k_adi,
                1.0/(gamma_adi-1.0));  // (HSW 94c)
          }
        }
        Real pgas = k_adi * std::pow(rho, gamma_adi);  // (HSW 94b)
        Real u0 = g_inv(I00,i) * u_0;
        Real u_3 = -l * u_0;
        Real u3 = g_inv(I33,i) * u_3;
        Real v3 = u3 / u0;
        set_state(pfl->w, pfl->w1, i, j, k, rho, pgas, 0.0, 0.0, v3);
      }
    }
  g.DeleteAthenaArray();
  g_inv.DeleteAthenaArray();
  pb->pcoord->PrimToCons(pfl->w, pfl->u);

  // Enroll boundary functions
  pb->pfluid->pf_bcs->EnrollBoundaryFunction(inner_x1, FixedInner);
  pb->pfluid->pf_bcs->EnrollBoundaryFunction(outer_x1, FixedOuter);
  pb->pfluid->pf_bcs->EnrollBoundaryFunction(inner_x2, FixedTop);
  pb->pfluid->pf_bcs->EnrollBoundaryFunction(outer_x2, FixedBottom);
  return;
}

// Inner boundary condition
// TODO: only works in Schwarzschild (assumed metric)
void FixedInner(MeshBlock *pmb, AthenaArray<Real> &cons,
                int is, int ie, int js, int je, int ks, int ke)
{
  // Set conserved values
  Real gamma_adi_red = gamma_adi / (gamma_adi - 1.0);
  Real r = pb->x1v(is);
  Real g_inv_00 = -1.0 / (1.0 - M/r);
  for (int k = ks; k <= ke; k++)
    for (int j = js; j <= je; j++)
    {
      Real theta = pb->x2v(j);
      Real g_inv_33 = 1.0 / (r*r * std::sin(theta)*std::sin(theta));
      Real neg_u_0_inv_sq = -g_inv_00 - g_inv_33 * l*l;  // (HSW 95a)
      Real u_0, rho;
      if (neg_u_0_inv_sq <= 0.0)
      {
        u_0 = -1.0;
        rho = rho_floor;
      }
      else
      {
        u_0 = -1.0 / std::sqrt(neg_u_0_inv_sq);
        if (u_0 <= -1.0)
        {
          u_0 = -1.0;
          rho = rho_floor;
        }
        else
        {
          Real epsilon = -1.0/gamma_adi * (1.0/u_0 + 1.0);  // (HSW 94a)
          rho = std::pow((gamma_adi-1.0) * epsilon / k_adi,
              1.0/(gamma_adi-1.0));  // (HSW 94c)
        }
      }
      Real pgas = k_adi * std::pow(rho, gamma_adi);  // (HSW 94b)
      Real u0 = g_inv_00 * u_0;
      Real u_3 = -l * u_0;
      Real u3 = g_inv_33 * u_3;
      Real v3 = u3 / u0;
      Real d = rho * u0;
      Real rho_h = rho + gamma_adi_red * pgas;
      Real e = rho_h * u0 * u_0 + pgas;
      Real m1 = 0.0;
      Real m2 = 0.0;
      Real m3 = rho_h * u0 * u_3;
      for (int i = is-NGHOST; i <= is; i++)
      {
        cons(IDN,k,j,i) = d;
        cons(IEN,k,j,i) = e;
        cons(IM1,k,j,i) = m1;
        cons(IM2,k,j,i) = m2;
        cons(IM3,k,j,i) = m3;
      }
    }
  return;
}

// Outer boundary condition
// TODO: only works in Schwarzschild (assumed metric)
void FixedOuter(MeshBlock *pmb, AthenaArray<Real> &cons,
                int is, int ie, int js, int je, int ks, int ke)
{
  // Set conserved values
  Real gamma_adi_red = gamma_adi / (gamma_adi - 1.0);
  Real r = pb->x1v(ie);
  Real g_inv_00 = -1.0 / (1.0 - M/r);
  for (int k = ks; k <= ke; k++)
    for (int j = js; j <= je; j++)
    {
      Real theta = pb->x2v(j);
      Real g_inv_33 = 1.0 / (r*r * std::sin(theta)*std::sin(theta));
      Real neg_u_0_inv_sq = -g_inv_00 - g_inv_33 * l*l;  // (HSW 95a)
      Real u_0, rho;
      if (neg_u_0_inv_sq <= 0.0)
      {
        u_0 = -1.0;
        rho = rho_floor;
      }
      else
      {
        u_0 = -1.0 / std::sqrt(neg_u_0_inv_sq);
        if (u_0 <= -1.0)
        {
          u_0 = -1.0;
          rho = rho_floor;
        }
        else
        {
          Real epsilon = -1.0/gamma_adi * (1.0/u_0 + 1.0);  // (HSW 94a)
          rho = std::pow((gamma_adi-1.0) * epsilon / k_adi,
              1.0/(gamma_adi-1.0));  // (HSW 94c)
        }
      }
      Real pgas = k_adi * std::pow(rho, gamma_adi);  // (HSW 94b)
      Real u0 = g_inv_00 * u_0;
      Real u_3 = -l * u_0;
      Real u3 = g_inv_33 * u_3;
      Real v3 = u3 / u0;
      Real d = rho * u0;
      Real rho_h = rho + gamma_adi_red * pgas;
      Real e = rho_h * u0 * u_0 + pgas;
      Real m1 = 0.0;
      Real m2 = 0.0;
      Real m3 = rho_h * u0 * u_3;
      for (int i = ie; i <= ie+NGHOST; i++)
      {
        cons(IDN,k,j,i) = d;
        cons(IEN,k,j,i) = e;
        cons(IM1,k,j,i) = m1;
        cons(IM2,k,j,i) = m2;
        cons(IM3,k,j,i) = m3;
      }
    }
  return;
}

// Top boundary condition
// TODO: only works in Schwarzschild (assumed metric)
void FixedTop(MeshBlock *pmb, AthenaArray<Real> &cons,
              int is, int ie, int js, int je, int ks, int ke)
{
  // Set conserved values
  Real gamma_adi_red = gamma_adi / (gamma_adi - 1.0);
  Real theta = pb->x2v(js);
  for (int k = ks; k <= ke; k++)
    for (int j = js-NGHOST; j <= js; j++)
    {
      for (int i = is; i <= ie; i++)
      {
        Real r = pb->x1v(i);
        Real g_inv_00 = -1.0 / (1.0 - M/r);
        Real g_inv_33 = 1.0 / (r*r * std::sin(theta)*std::sin(theta));
        Real neg_u_0_inv_sq = -g_inv_00 - g_inv_33 * l*l;  // (HSW 95a)
        Real u_0, rho;
        if (neg_u_0_inv_sq <= 0.0)
        {
          u_0 = -1.0;
          rho = rho_floor;
        }
        else
        {
          u_0 = -1.0 / std::sqrt(neg_u_0_inv_sq);
          if (u_0 <= -1.0)
          {
            u_0 = -1.0;
            rho = rho_floor;
          }
          else
          {
            Real epsilon = -1.0/gamma_adi * (1.0/u_0 + 1.0);  // (HSW 94a)
            rho = std::pow((gamma_adi-1.0) * epsilon / k_adi,
                1.0/(gamma_adi-1.0));  // (HSW 94c)
          }
        }
        Real pgas = k_adi * std::pow(rho, gamma_adi);  // (HSW 94b)
        Real u0 = g_inv_00 * u_0;
        Real u_3 = -l * u_0;
        Real u3 = g_inv_33 * u_3;
        Real v3 = u3 / u0;
        Real rho_h = rho + gamma_adi_red * pgas;
        cons(IDN,k,j,i) = rho * u0;
        cons(IEN,k,j,i) = rho_h * u0 * u_0 + pgas;
        cons(IM1,k,j,i) = 0.0;
        cons(IM2,k,j,i) = 0.0;
        cons(IM3,k,j,i) = rho_h * u0 * u_3;
      }
    }
  return;
}

// Bottom boundary condition
// TODO: only works in Schwarzschild (assumed metric)
void FixedBottom(MeshBlock *pmb, AthenaArray<Real> &cons,
                 int is, int ie, int js, int je, int ks, int ke)
{
  // Set conserved values
  Real gamma_adi_red = gamma_adi / (gamma_adi - 1.0);
  Real theta = pb->x2v(je);
  for (int k = ks; k <= ke; k++)
    for (int j = je; j <= je+NGHOST; j++)
    {
      for (int i = is; i <= ie; i++)
      {
        Real r = pb->x1v(i);
        Real g_inv_00 = -1.0 / (1.0 - M/r);
        Real g_inv_33 = 1.0 / (r*r * std::sin(theta)*std::sin(theta));
        Real neg_u_0_inv_sq = -g_inv_00 - g_inv_33 * l*l;  // (HSW 95a)
        Real u_0, rho;
        if (neg_u_0_inv_sq <= 0.0)
        {
          u_0 = -1.0;
          rho = rho_floor;
        }
        else
        {
          u_0 = -1.0 / std::sqrt(neg_u_0_inv_sq);
          if (u_0 <= -1.0)
          {
            u_0 = -1.0;
            rho = rho_floor;
          }
          else
          {
            Real epsilon = -1.0/gamma_adi * (1.0/u_0 + 1.0);  // (HSW 94a)
            rho = std::pow((gamma_adi-1.0) * epsilon / k_adi,
                1.0/(gamma_adi-1.0));  // (HSW 94c)
          }
        }
        Real pgas = k_adi * std::pow(rho, gamma_adi);  // (HSW 94b)
        Real u0 = g_inv_00 * u_0;
        Real u_3 = -l * u_0;
        Real u3 = g_inv_33 * u_3;
        Real v3 = u3 / u0;
        Real rho_h = rho + gamma_adi_red * pgas;
        cons(IDN,k,j,i) = rho * u0;
        cons(IEN,k,j,i) = rho_h * u0 * u_0 + pgas;
        cons(IM1,k,j,i) = 0.0;
        cons(IM2,k,j,i) = 0.0;
        cons(IM3,k,j,i) = rho_h * u0 * u_3;
      }
    }
  return;
}

// Function for setting conserved variables in a cell given the primitives
static void set_state(AthenaArray<Real> &prim, AthenaArray<Real> &prim_half, int i,
    int j, int k, Real rho, Real pgas, Real vx, Real vy, Real vz)
{
  prim(IDN,k,j,i) = prim_half(IDN,k,j,i) = rho;
  prim(IEN,k,j,i) = prim_half(IEN,k,j,i) = pgas;
  prim(IM1,k,j,i) = prim_half(IM1,k,j,i) = vx;
  prim(IM2,k,j,i) = prim_half(IM2,k,j,i) = vy;
  prim(IM3,k,j,i) = prim_half(IM3,k,j,i) = vz;
  return;
}
