// General relativistic Fishbone-Moncrief torus generator

// Primary header
#include "../mesh.hpp"

// C++ headers
#include <cmath>      // exp(), pow(), sin(), sqrt()
#include <algorithm>  // max()

// Athena headers
#include "../athena.hpp"                   // enums, Real
#include "../athena_arrays.hpp"            // AthenaArray
#include "../coordinates/coordinates.hpp"  // PrimToCons()
#include "../fluid/fluid.hpp"              // Fluid
#include "../fluid/bvals/bvals.hpp"        // EnrollBoundaryFunction()
#include "../fluid/eos/eos.hpp"            // GetGamma()
#include "../parameter_input.hpp"          // ParameterInput

// Declarations
void FixedInner(MeshBlock *pmb, AthenaArray<Real> &cons);
void FixedOuter(MeshBlock *pmb, AthenaArray<Real> &cons);
void FixedTop(MeshBlock *pmb, AthenaArray<Real> &cons);
void FixedBottom(MeshBlock *pmb, AthenaArray<Real> &cons);
static Real calculate_l(Real r);
static Real log_h_aux(Real r, Real sin_theta);
static void set_state(AthenaArray<Real> &prim, AthenaArray<Real> &prim_half, int i,
    int j, int k, Real rho, Real p_gas, Real vx, Real vy, Real vz);
static void calculate_conserved(Real r, Real &d, Real &e);

// Global variables
static Real gamma_adi, k_adi, r_edge, r_peak, l, rho_min, rho_pow, eps_min, eps_pow;
// TODO: put in a better place
static Real M = 1.0;

// Function for setting initial conditions
// Inputs:
//   pfl: Fluid
//   pfd: Field (unused)
//   pin: parameters
// Outputs: (none)
// Notes:
//   initializes Fishbone-Moncrief torus
//     sets both primitive and conserved variables
//   defines and enrolls fixed r- and theta-direction boundary conditions
//   references Fishbone & Moncrief 1976, ApJ 207 962 (FM)
//              Fishbone 1977, ApJ 215 323 (F)
//              Hawley, Smarr, & Wilson 1984, ApJ 277 296 (HSW)
// TODO: only works in Schwarzschild (assumed metric)
void Mesh::ProblemGenerator(Fluid *pfl, Field *pfd, ParameterInput *pin)
{
  // Prepare index bounds
  MeshBlock *pb = pfl->pmy_block;
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

  // Get ratio of specific heats
  gamma_adi = pf_eos->GetGamma();

  // TODO: read and set mass

  // Read other properties
  k_adi = pin->GetReal("problem", "k_adi");
  r_edge = pin->GetReal("problem", "r_edge");
  r_peak = pin->GetReal("problem", "r_peak");
  l = pin->GetReal("problem", "u0_u_3");
  rho_min = pin->GetReal("problem", "rho_min");
  rho_pow = pin->GetReal("problem", "rho_pow");
  eps_min = pin->GetReal("problem", "eps_min");
  eps_pow = pin->GetReal("problem", "eps_pow");

  // Reset l if valid r_peak given
  if (r_peak >= 0.0)
    l = calculate_l(r_peak);

  // Initialize primitive values
  Real log_h_edge = log_h_aux(r_edge, 1.0);
  for (int k = kl; k <= ku; k++)
    for (int j = jl; j <= ju; j++)
    {
      Real sin_theta = std::sin(pmb->x2v(j));
      for (int i = il; i <= iu; i++)
      {
        // Determine if we are in the torus
        Real r = pmb->x1v(i);
        Real log_h;
        bool in_torus = false;
        if (r >= r_edge)
        {
          log_h = log_h_aux(r, sin_theta) - log_h_edge;  // (FM 3.6)
          if (log_h >= 0.0)
            in_torus = true;
        }

        // Calculate velocity
        Real v3, u_0;
        if (in_torus)
        {
          u_0 = -1.0 / std::exp(log_h);                           // (HSW 91,93)
          Real u0 = -1.0 / (1.0-2.0*M/r) * u_0;
          Real exp_neg_2chi =
              (1.0-2.0*M/r) / (r*r * sin_theta*sin_theta);        // (FM 2.15,3.5)
          Real u_phi_proj_a =
              std::sqrt(1.0 + 4.0 * l*l * exp_neg_2chi);
          Real u_phi_proj = std::sqrt(0.5 * (u_phi_proj_a-1.0));  // (FM 3.3)
          Real exp_psi = r * sin_theta;                           // (FM 3.5)
          Real u_3 = exp_psi * u_phi_proj;                        // (FM 2.12, F 2.5)
          Real u3 = u_3 / (r*r * sin_theta*sin_theta);
          v3 = u3/u0;
        }
        else
          v3 = 0.0;

        // Calculate thermodynamic quantities
        Real epsilon;  // internal energy 1/(Gamma-1) * p_gas/rho
        Real rho, p_gas;
        if (in_torus)
        {
          epsilon = -1.0/gamma_adi * (1.0/u_0 + 1.0);        // (HSW 94a)
          rho = std::pow((gamma_adi-1.0) * epsilon / k_adi,
              1.0/(gamma_adi-1.0));                          // (HSW 94c)
          p_gas = (gamma_adi-1.0) * rho * epsilon;           // (HSW 94b)
        }
        else
        {
          rho = rho_min * std::pow(r/r_edge, rho_pow);
          epsilon = eps_min * std::pow(r/r_edge, eps_pow);
          p_gas = (gamma_adi-1.0) * rho * epsilon;          // (HSW 94b)
          // TODO: use alternate criterion below? or something else?
          //   need exactly 2 of (rho floor, pgas floor, epsilon floor, k_adi)
          //rho = rho_min;
          //p_gas = k_adi * std::pow(rho_min, gamma_adi);
        }

        // Set primitive values
        set_state(w, w1, i, j, k, rho, p_gas, 0.0, 0.0, v3);
      }
    }

  // Normalize density and pressure
  Real rho_peak;
  if (r_peak >= 0.0)
  {
    Real log_h_peak = log_h_aux(r_peak, 1.0) - log_h_edge;       // (FM 3.6)
    Real u_0_peak = -1.0 / std::exp(log_h_peak);                 // (HSW 91,93)
    Real epsilon_peak = -1.0/gamma_adi * (1.0/u_0_peak + 1.0);   // (HSW 94a)
    rho_peak = std::pow((gamma_adi-1.0) * epsilon_peak / k_adi,
        1.0/(gamma_adi-1.0));                                    // (HSW 94c)
  }
  else
  {
    rho_peak = 0.0;
    for (int k = kl; k <= ku; k++)
      for (int j = jl; j <= ju; j++)
        for (int i = il; i <= iu; i++)
          if (w(IM3,k,j,i) > 0.0)
            rho_peak = std::max(rho_peak, w(IDN,k,j,i));
  }
  for (int k = kl; k <= ku; k++)
    for (int j = jl; j <= ju; j++)
      for (int i = il; i <= iu; i++)
        if (w(IM3,k,j,i) > 0.0)
        {
          w(IDN,k,j,i) /= rho_peak;
          // TODO: renormalize according to pgas = k_adi * rho^gamma_adi?
          w(IEN,k,j,i) /= rho_peak;
          //w(IEN,k,j,i) /= k_adi * std::pow(rho_peak, gamma_adi);
          // TODO: better way to induce splash and crash
          //w(IM3,k,j,i) /= 2.0;
        }

  // Initialize conserved values
  pmb->pcoord->PrimToCons(w, u);

  // Enroll boundary functions
  pmb->pfluid->pf_bcs->EnrollBoundaryFunction(inner_x1, FixedInner);
  pmb->pfluid->pf_bcs->EnrollBoundaryFunction(outer_x1, FixedOuter);
  pmb->pfluid->pf_bcs->EnrollBoundaryFunction(inner_x2, FixedTop);
  pmb->pfluid->pf_bcs->EnrollBoundaryFunction(outer_x2, FixedBottom);
  return;
}

// Inner boundary condition
// Inputs:
//   pmb: pointer to block
// Outputs:
//   cons: conserved quantities set along inner x1-boundary
// Notes:
//   references Hawley, Smarr, & Wilson 1984, ApJ 277 296 (HSW)
// TODO: only works in Schwarzschild (assumed metric)
void FixedInner(MeshBlock *pmb, AthenaArray<Real> &cons)
{
  // Extract boundary indices
  int il = pmb->is - NGHOST;
  int iu = pmb->is;
  int jl = pmb->js;
  int ju = pmb->je;
  int kl = pmb->ks;
  int ku = pmb->ke;

  // Set conserved values
  Real r = pmb->x1v(iu);
  Real d, e;
  calculate_conserved(r, d, e);
  for (int k = kl; k <= ku; k++)
    for (int j = jl; j <= ju; j++)
      for (int i = il; i <= iu; i++)
      {
        cons(IDN,k,j,i) = d;
        cons(IEN,k,j,i) = e;
        cons(IM1,k,j,i) = 0.0;
        cons(IM2,k,j,i) = 0.0;
        cons(IM3,k,j,i) = 0.0;
      }
  return;
}

// Outer boundary condition
// Inputs:
//   pmb: pointer to block
// Outputs:
//   cons: conserved quantities set along outer x1-boundary
// Notes:
//   references Hawley, Smarr, & Wilson 1984, ApJ 277 296 (HSW)
// TODO: only works in Schwarzschild (assumed metric)
void FixedOuter(MeshBlock *pmb, AthenaArray<Real> &cons)
{
  // Extract boundary indices
  int il = pmb->ie;
  int iu = pmb->ie + NGHOST;
  int jl = pmb->js;
  int ju = pmb->je;
  int kl = pmb->ks;
  int ku = pmb->ke;

  // Set conserved values
  Real r = pmb->x1v(il);
  Real d, e;
  calculate_conserved(r, d, e);
  for (int k = kl; k <= ku; k++)
    for (int j = jl; j <= ju; j++)
      for (int i = il; i <= iu; i++)
      {
        cons(IDN,k,j,i) = d;
        cons(IEN,k,j,i) = e;
        cons(IM1,k,j,i) = 0.0;
        cons(IM2,k,j,i) = 0.0;
        cons(IM3,k,j,i) = 0.0;
      }
  return;
}

// Top boundary condition
// Inputs:
//   pmb: pointer to block
// Outputs:
//   cons: conserved quantities set along inner x1-boundary
// Notes:
//   references Hawley, Smarr, & Wilson 1984, ApJ 277 296 (HSW)
// TODO: only works in Schwarzschild (assumed metric)
void FixedTop(MeshBlock *pmb, AthenaArray<Real> &cons)
{
  // Extract boundary indices
  int il = pmb->is;
  int iu = pmb->ie;
  int jl = pmb->js - NGHOST;
  int ju = pmb->js;
  int kl = pmb->ks;
  int ku = pmb->ke;

  // Set conserved values
  for (int k = kl; k <= ku; k++)
    for (int j = jl; j <= ju; j++)
      for (int i = il; i <= iu; i++)
      {
        Real r = pmb->x1v(i);
        Real d, e;
        calculate_conserved(r, d, e);
        cons(IDN,k,j,i) = d;
        cons(IEN,k,j,i) = e;
        cons(IM1,k,j,i) = 0.0;
        cons(IM2,k,j,i) = 0.0;
        cons(IM3,k,j,i) = 0.0;
      }
  return;
}

// Bottom boundary condition
// Inputs:
//   pmb: pointer to block
// Outputs:
//   cons: conserved quantities set along inner x1-boundary
// Notes:
//   references Hawley, Smarr, & Wilson 1984, ApJ 277 296 (HSW)
// TODO: only works in Schwarzschild (assumed metric)
void FixedBottom(MeshBlock *pmb, AthenaArray<Real> &cons)
{
  // Extract boundary indices
  int il = pmb->is;
  int iu = pmb->ie;
  int jl = pmb->je;
  int ju = pmb->je + NGHOST;
  int kl = pmb->ks;
  int ku = pmb->ke;

  // Set conserved values
  for (int k = kl; k <= ku; k++)
    for (int j = jl; j <= ju; j++)
      for (int i = il; i <= iu; i++)
      {
        Real r = pmb->x1v(i);
        Real d, e;
        calculate_conserved(r, d, e);
        cons(IDN,k,j,i) = d;
        cons(IEN,k,j,i) = e;
        cons(IM1,k,j,i) = 0.0;
        cons(IM2,k,j,i) = 0.0;
        cons(IM3,k,j,i) = 0.0;
      }
  return;
}

// Function for calculating angular momentum variable l
// Inputs:
//   r: radius of pressure maximum
// Outputs:
//   returned value: l = u^t u_phi
// Notes:
//   beware many different definitions of l abound
//     this is *not* -u_phi/u_t
//   uses simplified (Schwarzschild) version of formula from Harm
//     see lfish_calc() in init.c in Harm
static Real calculate_l(Real r)
{
  Real numerator = 1.0 + r*r*std::sqrt(r) * (r-2.0);
  Real denominator = r * (r-2.0) * (r-3.0);
  return numerator/denominator;
}

// Function for helping to calculate enthalpy
// Inputs:
//   r: radial coordinate
//   theta: polar coordinate
// Outputs:
//   returned value: log(h)
// Notes:
//   enthalpy defined here as h = p_gas/rho
//   implements first half of (3.6) in Fishbone & Moncrief 1976, ApJ 207 962
static Real log_h_aux(Real r, Real sin_theta)
{
  Real a = 1.0 - 2.0*M/r;
  Real b = (2.0*l / (r*sin_theta));
  Real c = std::sqrt(1.0 + b*b * a);
  Real d = 0.5 * std::log(1.0/a * (1.0 + c));
  Real e = 0.5 * c;
  return d - e;
}

// Function for setting conserved variables in a cell given the primitives
// Inputs:
//   i, j, k: indices for cell to be set
//   rho: density
//   p_gas: gas pressure
//   vx, vy, vz: 3-velocity
// Outputs:
//   prim, prim_half: primitive values set
static void set_state(AthenaArray<Real> &prim, AthenaArray<Real> &prim_half, int i,
    int j, int k, Real rho, Real p_gas, Real vx, Real vy, Real vz)
{
  prim(IDN,k,j,i) = prim_half(IDN,k,j,i) = rho;
  prim(IEN,k,j,i) = prim_half(IEN,k,j,i) = p_gas;
  prim(IM1,k,j,i) = prim_half(IM1,k,j,i) = vx;
  prim(IM2,k,j,i) = prim_half(IM2,k,j,i) = vy;
  prim(IM3,k,j,i) = prim_half(IM3,k,j,i) = vz;
  return;
}

// Function for calculating conserved quantities as a function of position
// Inputs:
//   r: radial coordinate
// Outputs:
//   d: conserved density rho * u^0
//   e: conserved energy T^0_0
// TODO: only works in Schwarzschild (assumed metric)
static void calculate_conserved(Real r, Real &d, Real &e)
{
  Real gamma_adi_red = gamma_adi/(gamma_adi-1.0);
  Real g_00 = -(1.0-2.0*M/r);
  Real u0 = std::sqrt(-1.0/g_00);
  Real rho = rho_min * std::pow(r/r_edge, rho_pow);
  Real epsilon = eps_min * std::pow(r/r_edge, eps_pow);
  Real p_gas = (gamma_adi-1.0) * rho * epsilon;          // (HSW 94b)
  //Real rho = rho_min;
  //Real p_gas = k_adi * std::pow(rho_min, gamma_adi);
  d = rho * u0;
  e = -(rho + gamma_adi_red * p_gas) + p_gas;
  return;
}
