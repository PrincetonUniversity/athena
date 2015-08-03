// General relativistic problem generator for dust falling onto black hole

// Primary header
#include "../mesh.hpp"

// C++ headers
#include <cassert>  // assert
#include <cmath>    // pow(), sin(), sqrt()

// Athena headers
#include "../athena.hpp"                   // enums, Real
#include "../athena_arrays.hpp"            // AthenaArray
#include "../parameter_input.hpp"          // ParameterInput
#include "../bvals/bvals.hpp"              // BoundaryValues, InterfaceField
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../fluid/fluid.hpp"              // Fluid
#include "../fluid/eos/eos.hpp"            // FluidEqnOfState
#include "../field/field.hpp"              // Field

// Declarations
void OutflowPrimInnerFluid(MeshBlock *pmb, AthenaArray<Real> &cons,
    int is, int ie, int js, int je, int ks, int ke);
void FixedOuterFluid(MeshBlock *pmb, AthenaArray<Real> &cons,
    int is, int ie, int js, int je, int ks, int ke);

// Function for setting initial conditions
// Inputs:
//   pfl: Fluid
//   pfd: Field (unused)
//   pin: parameters
// Outputs: (none)
// Notes:
//   assumes x3 is axisymmetric direction
void Mesh::ProblemGenerator(Fluid *pfl, Field *pfd, ParameterInput *pin)
{
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

  // Get mass and spin of black hole
  Real m = pmb->pcoord->GetMass();
  Real a = pmb->pcoord->GetSpin();
  Real a2 = SQR(a);

  // Get ratio of specific heats
  Real gamma_adi = pfl->pf_eos->GetGamma();

  // Read other properties
  Real e = pin->GetReal("problem", "energy");
  Real lz = pin->GetReal("problem", "l_z");
  Real rho_min = pin->GetReal("fluid", "rho_min");
  Real rho_pow = pin->GetReal("fluid", "rho_pow");
  Real u_min = pin->GetReal("fluid", "u_min");
  Real u_pow = pin->GetReal("fluid", "u_pow");

  // Initialize primitive values
  AthenaArray<Real> g, gi;
  g.NewAthenaArray(NMETRIC, iu+1);
  gi.NewAthenaArray(NMETRIC, iu+1);
  for (int j = jl; j <= ju; j++)
  {
    pmb->pcoord->CellMetric(kl, j, il, iu, g, gi);
    for (int i = il; i <= iu; i++)
    {
      // Get Boyer-Lindquist coordinates of cell
      Real r, theta, phi;
      pmb->pcoord->GetBoyerLindquistCoordinates(pmb->pcoord->x1v(i),
          pmb->pcoord->x2v(j), pmb->pcoord->x3v(kl), &r, &theta, &phi);
      /*Real r2 = SQR(r);
      Real sin = std::sin(theta);
      Real sin2 = SQR(sin);
      Real cos2 = 1.0 - sin2;
      Real sigma = r2 + a2 * cos2;
      Real delta = r2 - 2.0*m*r + a2;
      Real pp = e * (r2 + a2) - lz * a;
      Real kk = SQR(a * e * sin - lz / sin) + a2 * cos2;
      Real rr = std::sqrt(SQR(pp) - delta * (r2 + kk));*/

      // Calculate primitives depending on location
      Real rho = rho_min * std::pow(r, rho_pow);
      Real pgas = (gamma_adi-1.0) * u_min * std::pow(r, u_pow);
      /*Real u0_bl = 1.0/sigma * (-a * (a * e * sin2 - lz) + (r2+a2)/delta * (rr - pp));
      Real u1_bl = rr/sigma;
      Real u2_bl = 0.0;
      Real u3_bl = 1.0/sigma * (-(a * e - lz / sin2) + a/delta * (rr - pp));
      Real u0, u1, u2, u3;
      pmb->pcoord->TransformVectorCell(u0_bl, u1_bl, u2_bl, u3_bl, pmb->ks, j, i,
          &u0, &u1, &u2, &u3);
      Real v1 = u1 / u0;
      Real v2 = u2 / u0;
      Real v3 = u3 / u0;*/
      Real uu1 = 0.0;
      Real uu2 = 0.0;
      Real uu3 = 0.0;

      // Set primitive values
      for (int k = kl; k <= ku; k++)
      {
        pfl->w(IDN,k,j,i) = pfl->w1(IDN,k,j,i) = rho;
        pfl->w(IEN,k,j,i) = pfl->w1(IEN,k,j,i) = pgas;
        pfl->w(IVX,k,j,i) = pfl->w1(IM1,k,j,i) = uu1;
        pfl->w(IVY,k,j,i) = pfl->w1(IM2,k,j,i) = uu2;
        pfl->w(IVZ,k,j,i) = pfl->w1(IM3,k,j,i) = uu3;
      }
    }
  }
  g.DeleteAthenaArray();
  gi.DeleteAthenaArray();

  // Initialize conserved values
  AthenaArray<Real> bb;
  bb.NewAthenaArray(3, ku+1, ju+1, iu+1);
  pmb->pfluid->pf_eos->PrimitiveToConserved(pfl->w, bb, pfl->u);  
  bb.DeleteAthenaArray();

  // Enroll boundary functions
  pmb->pbval->EnrollFluidBoundaryFunction(inner_x1, OutflowPrimInnerFluid);
  pmb->pbval->EnrollFluidBoundaryFunction(outer_x1, FixedOuterFluid);
  return;
}

// Inner fluid boundary condition
// Inputs:
//   pmb: pointer to block
// Outputs:
//   cons: conserved quantities set along inner x1-boundary
// Notes:
//   TODO: remove prim hack
//   TODO: note hack is wrong (assumes wrong primitives)
void OutflowPrimInnerFluid(MeshBlock *pmb, AthenaArray<Real> &cons,
    int is, int ie, int js, int je, int ks, int ke)
{
  int il = is - NGHOST;
  int iu = is;
  int jl = js;
  int ju = je;
  int kl = ks;
  int ku = ke;
  AthenaArray<Real> *pprim;
  if (&cons == &pmb->pfluid->u)
    pprim = &pmb->pfluid->w;
  else if (&cons == &pmb->pfluid->u1)
    pprim = &pmb->pfluid->w1;
  else
    assert(0);
  AthenaArray<Real> g, gi;
  g.NewAthenaArray(NMETRIC,iu+1);
  gi.NewAthenaArray(NMETRIC,iu+1);
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
    {
      pmb->pcoord->CellMetric(k, j, il, iu, g, gi);
      Real alpha = std::sqrt(-1.0/gi(I00,is));
      Real v1 = (*pprim)(IVX,k,j,is);
      Real v2 = (*pprim)(IVY,k,j,is);
      Real v3 = (*pprim)(IVZ,k,j,is);
      Real tmp = g(I00,is) + 2.0 * (g(I01,is)*v1 + g(I02,is)*v2 + g(I03,is)*v3)
          + g(I11,is)*v1*v1 + 2.0*g(I12,is)*v1*v2 + 2.0*g(I13,is)*v1*v3
          + g(I22,is)*v2*v2 + 2.0*g(I23,is)*v2*v3
          + g(I33,is)*v3*v3;
      Real u0 = std::sqrt(-1.0/tmp);
      Real u1 = u0 * v1;
      Real u2 = u0 * v2;
      Real u3 = u0 * v3;
      Real gamma = alpha * u0;
      Real utilde1 = u1 + alpha * gamma * gi(I01,is);
      Real utilde2 = u2 + alpha * gamma * gi(I02,is);
      Real utilde3 = u3 + alpha * gamma * gi(I03,is);
      for (int i = il; i <= iu; ++i)
      {
        Real &rho = (*pprim)(IDN,k,j,i);
        Real &pgas = (*pprim)(IEN,k,j,i);
        rho = (*pprim)(IDN,k,j,is);
        pgas = (*pprim)(IEN,k,j,is);
        Real alpha_new = std::sqrt(-1.0/gi(I00,i));
        Real u0_new = gamma / alpha_new;
        Real u1_new = utilde1 - alpha_new * gamma * gi(I01,i);
        Real u2_new = utilde2 - alpha_new * gamma * gi(I02,i);
        Real u3_new = utilde3 - alpha_new * gamma * gi(I03,i);
        (*pprim)(IVX,k,j,i) = u1_new / u0_new;
        (*pprim)(IVY,k,j,i) = u2_new / u0_new;
        (*pprim)(IVZ,k,j,i) = u3_new / u0_new;
        Real u_0_new, u_1_new, u_2_new, u_3_new;
        pmb->pcoord->LowerVectorCell(u0_new, u1_new, u2_new, u3_new, k, j, i, &u_0_new,
            &u_1_new, &u_2_new, &u_3_new);
        Real gamma_adi = pmb->pfluid->pf_eos->GetGamma();
        Real gamma_prime = gamma_adi/(gamma_adi-1.0);
        Real wgas = rho + gamma_prime * pgas;
        cons(IDN,k,j,i) = rho * u0_new;
        cons(IEN,k,j,i) = wgas * u0_new * u_0_new + pgas;
        cons(IM1,k,j,i) = wgas * u0_new * u_1_new;
        cons(IM2,k,j,i) = wgas * u0_new * u_2_new;
        cons(IM3,k,j,i) = wgas * u0_new * u_3_new;
      }
    }
  g.DeleteAthenaArray();
  gi.DeleteAthenaArray();
  return;
}

// Outer fluid boundary condition
// Inputs:
//   pmb: pointer to block
// Outputs:
//   cons: conserved quantities set along outer x1-boundary
// Notes:
//   remains unchanged
void FixedOuterFluid(MeshBlock *pmb, AthenaArray<Real> &cons,
    int is, int ie, int js, int je, int ks, int ke)
{
  return;
}
