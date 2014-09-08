// General relativistic black hole accretion generator

// Primary header
#include "../fluid/fluid.hpp"

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
//   sets conserved variables according to input primitives
void Fluid::InitFluid(ParameterInput *pin)
{
  // Prepare index bounds
  MeshBlock *pb = pmy_block;
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
  Real gamma_adi = pf_eos->GetGamma();
  Real gamma_adi_red = gamma_adi / (gamma_adi - 1.0);

  // TODO: read and set mass

  // Read outer boundary state
  Real rho_outer = pin->GetReal("problem", "rho_outer");
  Real pgas_outer = pin->GetReal("problem", "pgas_outer");
  Real v1_outer = pin->GetReal("problem", "v1_outer");
  Real v2_outer = pin->GetReal("problem", "v2_outer");
  Real v3_outer = pin->GetReal("problem", "v3_outer");

  // Read inner initial state
  Real rho_inner = pin->GetReal("problem", "rho_inner");
  Real pgas_inner = pin->GetReal("problem", "pgas_inner");
  Real v1_inner = pin->GetReal("problem", "v1_inner");
  Real v2_inner = pin->GetReal("problem", "v2_inner");
  Real v3_inner = pin->GetReal("problem", "v3_inner");

  // Calculate slopes
  Real rho_slope = (rho_outer - rho_inner) / (pb->x1v(iu) - pb->x1v(il));
  Real pgas_slope = (pgas_outer - pgas_inner) / (pb->x1v(iu) - pb->x1v(il));
  Real v1_slope = (v1_outer - v1_inner) / (pb->x1v(iu) - pb->x1v(il));
  Real v2_slope = (v2_outer - v2_inner) / (pb->x1v(iu) - pb->x1v(il));
  Real v3_slope = (v3_outer - v3_inner) / (pb->x1v(iu) - pb->x1v(il));

  // Initialize the infalling material
  for (int k = kl; k <= ku; k++)
    for (int j = jl; j <= ju; j++)
      for (int i = il; i <= iu; i++)
      {
        Real displacement = pb->x1v(i) - pb->x1v(il);
        Real rho_init = rho_inner + rho_slope * displacement;
        Real pgas_init = pgas_inner + pgas_slope * displacement;
        Real v1_init = v1_inner + v1_slope * displacement;
        Real v2_init = v2_inner + v2_slope * displacement;
        Real v3_init = v3_inner + v3_slope * displacement;
        set_state(w, w1, i, j, k, rho_init, pgas_init, v1_init, v2_init, v3_init);
      }
  pmy_block->pcoord->PrimToCons(w, u);
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
