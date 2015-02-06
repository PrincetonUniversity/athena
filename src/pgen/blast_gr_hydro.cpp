// Spherical blast wave generator for general relativistic hydrodynamics

// Primary header
#include "../mesh.hpp"

// C++ headers
#include <algorithm>  // min()
#include <cmath>      // sqrt()

// Athena headers
#include "../athena.hpp"                   // enums, Real
#include "../athena_arrays.hpp"            // AthenaArray
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../field/field.hpp"              // Field
#include "../fluid/fluid.hpp"              // Fluid
#include "../fluid/eos/eos.hpp"            // GetGamma()
#include "../parameter_input.hpp"          // ParameterInput

// Declarations
static void SetPrimCons(
    Real rho, Real pgas, Real gamma_adi_red, const AthenaArray<Real> &g,
    AthenaArray<Real> &prim, AthenaArray<Real> &prim_half, AthenaArray<Real> &cons,
    int i, int j, int k);

// Function for setting initial conditions
// Inputs:
//   pfl: Fluid
//   pfd: Field
//   pin: parameters
// Outputs: (none)
// Notes:
//   sets conserved variables according to input primitives
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

  // Read and set ratio of specific heats
  Real gamma_adi = pfl->pf_eos->GetGamma();
  Real gamma_adi_red = gamma_adi / (gamma_adi - 1.0);

  // Read problem parameters
  Real num_x = pin->GetReal("problem", "num_x");
  Real num_y = pin->GetReal("problem", "num_y");
  Real x_spacing = pin->GetReal("problem", "x_spacing");
  Real y_spacing = pin->GetReal("problem", "y_spacing");
  Real radius = pin->GetReal("problem", "radius");
  Real rho_inner = pin->GetReal("problem", "rho_inner");
  Real pgas_inner = pin->GetReal("problem", "pgas_inner");
  Real rho_outer = pin->GetReal("problem", "rho_outer");
  Real pgas_outer = pin->GetReal("problem", "pgas_outer");

  // Initialize the problem
  int ncells1 = pfl->pmy_block->block_size.nx1 + 2*NGHOST;
  AthenaArray<Real> g, gi;
  g.NewAthenaArray(NMETRIC,ncells1);
  gi.NewAthenaArray(NMETRIC,ncells1);
  for (int k = kl; k <= ku; k++)
    for (int j = jl; j <= ju; j++)
    {
      pfl->pmy_block->pcoord->CellMetric(k, j, g, gi);
      for (int i = il; i <= iu; i++)
      {
        Real x1 = pb->x1v(i);
        Real x2 = pb->x2v(j);
        Real x3 = pb->x3v(k);
        Real min_separation = pb->pcoord->DistanceBetweenPoints(x1, x2, x3, 0.0, 0.0,
            0.0);
        for (int x_index = -num_x; x_index <= num_x; ++x_index)
        {
          Real center_x = x_index * x_spacing;
          for (int y_index = -num_y; y_index <= num_y; ++y_index)
          {
            Real center_y = y_index * y_spacing; 
            min_separation = std::min(min_separation,
                pb->pcoord->DistanceBetweenPoints(x1, x2, x3, center_x, center_y, 0.0));
          }
        }
        if (min_separation < radius)
          SetPrimCons(rho_inner, pgas_inner, gamma_adi_red, g,
              pfl->w, pfl->w1, pfl->u, i, j, k);
        else
          SetPrimCons(rho_outer, pgas_outer, gamma_adi_red, g,
              pfl->w, pfl->w1, pfl->u, i, j, k);
      }
    }
  g.DeleteAthenaArray();
  gi.DeleteAthenaArray();
  return;
}

// Function for setting all variables in a cell given the primitives
static void SetPrimCons(
    Real rho, Real pgas, Real gamma_adi_red, const AthenaArray<Real> &g,
    AthenaArray<Real> &prim, AthenaArray<Real> &prim_half, AthenaArray<Real> &cons,
    int i, int j, int k)
{
  // Set primitives
  prim(IDN,k,j,i) = prim_half(IDN,k,j,i) = rho;
  prim(IEN,k,j,i) = prim_half(IEN,k,j,i) = pgas;
  prim(IM1,k,j,i) = prim_half(IM1,k,j,i) = 0.0;
  prim(IM2,k,j,i) = prim_half(IM2,k,j,i) = 0.0;
  prim(IM3,k,j,i) = prim_half(IM3,k,j,i) = 0.0;

  // Set conserved quantities
  Real rho_h = rho + gamma_adi_red * pgas;
  Real u0 = std::sqrt(-1.0/g(I00,i));
  cons(IDN,k,j,i) = rho * u0;
  cons(IEN,k,j,i) = -rho_h + pgas;
  cons(IM1,k,j,i) = -rho_h * g(I01,i)/g(I00,i);
  cons(IM2,k,j,i) = -rho_h * g(I02,i)/g(I00,i);
  cons(IM3,k,j,i) = -rho_h * g(I03,i)/g(I00,i);
  return;
}
