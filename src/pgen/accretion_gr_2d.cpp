// General relativistic black hole accretion generator, azimuthally symmetric flows

// Primary header
#include "../fluid/fluid.hpp"

// C++ headers
#include <cmath>  // pow(), sqrt()

// Athena headers
#include "../athena.hpp"                   // enums, Real
#include "../athena_arrays.hpp"            // AthenaArray
#include "../coordinates/coordinates.hpp"  // PrimToCons()
#include "../fluid/eos/eos.hpp"            // GetGamma()
#include "../mesh.hpp"                     // MeshBlock, MeshDomain, Mesh
#include "../parameter_input.hpp"          // ParameterInput

// Declarations
static void set_state(AthenaArray<Real> &prim, AthenaArray<Real> &prim_half, int i,
    int j, int k, Real rho, Real pgas, Real vx, Real vy, Real vz);

// Function for setting initial conditions
// Inputs:
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

  // Read and set ratio of specific heats
  Real gamma_adi = pfl->pf_eos->GetGamma();
  Real gamma_adi_red = gamma_adi / (gamma_adi - 1.0);

  // Read other properties
  Real k_adi = pin->GetReal("problem", "k");
  Real l = pin->GetReal("problem", "l");
  Real rho_floor = pin->GetReal("problem", "rho_floor");

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
        Real u_0, rho, pgas;
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
          pgas = (gamma_adi-1.0) * rho * epsilon;  // (HWS 94b)
        }
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
