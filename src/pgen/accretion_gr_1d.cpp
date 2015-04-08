// General relativistic black hole accretion generator, spherically symmetric flows

// Primary header
#include "../mesh.hpp"

// Athena headers
#include "../athena.hpp"                   // enums, Real
#include "../athena_arrays.hpp"            // AthenaArray
#include "../coordinates/coordinates.hpp"  // PrimToCons()
#include "../bvals/bvals.hpp"              // EnrollBoundaryFunction()
#include "../fluid/fluid.hpp"              // Fluid
#include "../fluid/eos/eos.hpp"            // GetGamma()
#include "../field/field.hpp"              // Field
#include "../parameter_input.hpp"          // ParameterInput

// Declarations
void FixedInner(MeshBlock *pmb, AthenaArray<Real> &cons,
                int is, int ie, int js, int je, int ks, int ke);
void FixedOuter(MeshBlock *pmb, AthenaArray<Real> &cons,
                int is, int ie, int js, int je, int ks, int ke);
static void set_state(Real rho, Real pgas, Real v1, Real v2, Real v3,
    AthenaArray<Real> &prim, AthenaArray<Real> &prim_half, int i, int j, int k);

// Global variables
static Real d_inner, e_inner, m1_inner, m2_inner, m3_inner;
static Real d_outer, e_outer, m1_outer, m2_outer, m3_outer;

// Function for setting initial conditions
// Inputs:
//   pfl: Fluid
//   pfd: Field (unused)
//   pin: parameters
// Outputs: (none)
// Notes:
//   sets primitive and conserved variables according to input primitives
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

  // TODO: read and set mass

  // Read inner initial hydro state
  Real rho_inner = pin->GetReal("problem", "rho_inner");
  Real pgas_inner = pin->GetReal("problem", "pgas_inner");
  Real v1_inner = pin->GetReal("problem", "v1_inner");
  Real v2_inner = pin->GetReal("problem", "v2_inner");
  Real v3_inner = pin->GetReal("problem", "v3_inner");

  // Read outer initial hydro state
  Real rho_outer = pin->GetReal("problem", "rho_outer");
  Real pgas_outer = pin->GetReal("problem", "pgas_outer");
  Real v1_outer = pin->GetReal("problem", "v1_outer");
  Real v2_outer = pin->GetReal("problem", "v2_outer");
  Real v3_outer = pin->GetReal("problem", "v3_outer");

  // Read initial magnetic field
  Real b1_flux, b2_flux, b3_flux;
  if (MAGNETIC_FIELDS_ENABLED)
  {
    b1_flux = pin->GetReal("problem", "b1_flux");
    b2_flux = pin->GetReal("problem", "b2_flux");
    b3_flux = pin->GetReal("problem", "b3_flux");
  }

  // Calculate hydro slopes
  Real rho_slope = (rho_outer - rho_inner) / (pb->x1v(iu) - pb->x1v(il));
  Real pgas_slope = (pgas_outer - pgas_inner) / (pb->x1v(iu) - pb->x1v(il));
  Real v1_slope = (v1_outer - v1_inner) / (pb->x1v(iu) - pb->x1v(il));
  Real v2_slope = (v2_outer - v2_inner) / (pb->x1v(iu) - pb->x1v(il));
  Real v3_slope = (v3_outer - v3_inner) / (pb->x1v(iu) - pb->x1v(il));

  // Prepare arrays for areas and magnetic fields
  AthenaArray<Real> a1, a2m, a2p, a3m, a3p, b;
  a1.NewAthenaArray(iu+1);
  a2m.NewAthenaArray(iu);
  a2p.NewAthenaArray(iu);
  a3m.NewAthenaArray(iu);
  a3p.NewAthenaArray(iu);
  b.NewAthenaArray(3,ku+1,ju+1,iu+1);

  // Initialize magnetic field
  if (MAGNETIC_FIELDS_ENABLED)
    for (int k = kl; k <= ku; k++)
    {
      Real interp_param_k = (pb->x3v(k) - pb->x3f(k)) / pb->dx3f(k);
      for (int j = jl; j <= ju; j++)
      {
        Real interp_param_j = (pb->x2v(j) - pb->x2f(j)) / pb->dx2f(j);
        pb->pcoord->Face1Area(k, j, il, iu+1, a1);
        pb->pcoord->Face2Area(k, j, il, iu, a2m);
        pb->pcoord->Face2Area(k, j+1, il, iu, a2p);
        pb->pcoord->Face3Area(k, j, il, iu, a3m);
        pb->pcoord->Face3Area(k+1, j, il, iu, a3p);
        for (int i = il; i <= iu; i++)
        {
          Real interp_param_i = (pb->x1v(i) - pb->x1f(i)) / pb->dx1f(i);
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

  // Initialize primitives
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
        set_state(rho_init, pgas_init, v1_init, v2_init, v3_init,
            pfl->w, pfl->w1, i, j, k);
      }

  // Initialize conserved variables
  pb->pcoord->PrimToCons(pfl->w, b, pfl->u);

  // Delete area and magnetic field arrays
  a1.DeleteAthenaArray();
  a2m.DeleteAthenaArray();
  a2p.DeleteAthenaArray();
  a3m.DeleteAthenaArray();
  a3p.DeleteAthenaArray();
  b.DeleteAthenaArray();

  // Read inner boundary state
  d_inner = pin->GetReal("problem", "d_inner");
  e_inner = pin->GetReal("problem", "e_inner");
  m1_inner = pin->GetReal("problem", "m1_inner");
  m2_inner = pin->GetReal("problem", "m2_inner");
  m3_inner = pin->GetReal("problem", "m3_inner");

  // Read outer boundary state
  d_outer = pin->GetReal("problem", "d_outer");
  e_outer = pin->GetReal("problem", "e_outer");
  m1_outer = pin->GetReal("problem", "m1_outer");
  m2_outer = pin->GetReal("problem", "m2_outer");
  m3_outer = pin->GetReal("problem", "m3_outer");

  // Enroll boundary functions
  pb->pbval->EnrollFluidBoundaryFunction(inner_x1, FixedInner);
  pb->pbval->EnrollFluidBoundaryFunction(outer_x1, FixedOuter);
  return;
}

// Inner boundary condition
void FixedInner(MeshBlock *pmb, AthenaArray<Real> &cons,
                int is, int ie, int js, int je, int ks, int ke)
{
  // Set conserved values
  for (int k = ks; k <= ke; k++)
    for (int j = js; j <= je; j++)
      for (int i = is-NGHOST; i <= is; i++)
      {
        cons(IDN,k,j,i) = d_inner;
        cons(IEN,k,j,i) = e_inner;
        cons(IM1,k,j,i) = m1_inner;
        cons(IM2,k,j,i) = m2_inner;
        cons(IM3,k,j,i) = m3_inner;
      }
  return;
}

// Outer boundary condition
void FixedOuter(MeshBlock *pmb, AthenaArray<Real> &cons,
                int is, int ie, int js, int je, int ks, int ke)
{
  // Set conserved values
  for (int k = ks; k <= ke; k++)
    for (int j = js; j <= je; j++)
      for (int i = ie; i <= ie+NGHOST; i++)
      {
        cons(IDN,k,j,i) = d_outer;
        cons(IEN,k,j,i) = e_outer;
        cons(IM1,k,j,i) = m1_outer;
        cons(IM2,k,j,i) = m2_outer;
        cons(IM3,k,j,i) = m3_outer;
      }
  return;
}

// Function for setting conserved variables in a cell given the primitives
static void set_state(Real rho, Real pgas, Real v1, Real v2, Real v3,
    AthenaArray<Real> &prim, AthenaArray<Real> &prim_half, int i, int j, int k)
{
  prim(IDN,k,j,i) = prim_half(IDN,k,j,i) = rho;
  prim(IEN,k,j,i) = prim_half(IEN,k,j,i) = pgas;
  prim(IM1,k,j,i) = prim_half(IM1,k,j,i) = v1;
  prim(IM2,k,j,i) = prim_half(IM2,k,j,i) = v2;
  prim(IM3,k,j,i) = prim_half(IM3,k,j,i) = v3;
  return;
}
