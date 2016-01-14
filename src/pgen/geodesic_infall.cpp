// General relativistic problem generator for dust falling onto black hole

// Primary header
#include "../mesh.hpp"

#error "geodesic_infall.cpp is outdated and must be rewritten."

#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator does not support magnetic fields"
#endif

// C++ headers
#include <cassert>  // assert
#include <cmath>    // pow(), sin(), sqrt()

// Athena headers
#include "../athena.hpp"                   // enums, Real
#include "../athena_arrays.hpp"            // AthenaArray
#include "../parameter_input.hpp"          // ParameterInput
#include "../bvals/bvals.hpp"              // BoundaryValues, FaceField
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../hydro/hydro.hpp"
#include "../hydro/eos/eos.hpp"
#include "../field/field.hpp"              // Field

// Declarations
void OutflowInner(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &bb, int is, int ie, int js, int je, int ks, int ke);
void FixedOuter(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &bb, int is, int ie, int js, int je, int ks, int ke);


// Function for initializing global mesh properties
void Mesh::InitUserMeshProperties(ParameterInput *pin)
{
  // Enroll boundary functions
  EnrollUserBoundaryFunction(INNER_X1, OutflowInner);
  EnrollUserBoundaryFunction(OUTER_X1, FixedOuter);
  return;
}


// Function for cleaning up global mesh properties
void Mesh::TerminateUserMeshProperties(void)
{
  return;
}


// Function for setting initial conditions
// Inputs:
//   phyd: Hydro
//   pfld: Field (unused)
//   pin: parameters
// Outputs: (none)
// Notes:
//   assumes x3 is axisymmetric direction
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
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

  // Get mass and spin of black hole
  Real m = pcoord->GetMass();
  Real a = pcoord->GetSpin();
  Real a2 = SQR(a);

  // Get ratio of specific heats
  Real gamma_adi = phydro->peos->GetGamma();

  // Read other properties
  Real e = pin->GetReal("problem", "energy");
  Real lz = pin->GetReal("problem", "l_z");
  Real rho_min = pin->GetReal("hydro", "rho_min");
  Real rho_pow = pin->GetReal("hydro", "rho_pow");
  Real u_min = pin->GetReal("hydro", "u_min");
  Real u_pow = pin->GetReal("hydro", "u_pow");

  // Initialize primitive values
  AthenaArray<Real> g, gi;
  g.NewAthenaArray(NMETRIC, iu+1);
  gi.NewAthenaArray(NMETRIC, iu+1);
  for (int j = jl; j <= ju; j++)
  {
    pcoord->CellMetric(kl, j, il, iu, g, gi);
    for (int i = il; i <= iu; i++)
    {
      // Get Boyer-Lindquist coordinates of cell
      Real r, theta, phi;
      pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i),
          pcoord->x2v(j), pcoord->x3v(kl), &r, &theta, &phi);
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
      pcoord->TransformVectorCell(u0_bl, u1_bl, u2_bl, u3_bl, ks, j, i,
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
        phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) = rho;
        phydro->w(IEN,k,j,i) = phydro->w1(IEN,k,j,i) = pgas;
        phydro->w(IVX,k,j,i) = phydro->w1(IM1,k,j,i) = uu1;
        phydro->w(IVY,k,j,i) = phydro->w1(IM2,k,j,i) = uu2;
        phydro->w(IVZ,k,j,i) = phydro->w1(IM3,k,j,i) = uu3;
      }
    }
  }
  g.DeleteAthenaArray();
  gi.DeleteAthenaArray();

  // Initialize conserved values
  AthenaArray<Real> bb;
  bb.NewAthenaArray(3, ku+1, ju+1, iu+1);
  phydro->peos->PrimitiveToConserved(phydro->w, bb, phydro->u, pcoord,
                                          il, iu, jl, ju, kl, ku);
  bb.DeleteAthenaArray();

  return;
}


// User-defined work function called every time step
void MeshBlock::UserWorkInLoop(void)
{
  return;
}


// Inner boundary condition
// Inputs:
//   pmb: pointer to block
// Outputs:
//   prim: primitive quantities set along inner x1-boundary
// Notes:
//   TODO: remove prim hack
//   TODO: note hack is wrong (assumes wrong primitives)
void OutflowInner(MeshBlock *pmb,  Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &bb, int is, int ie, int js, int je, int ks, int ke)
{
  int il = is - NGHOST;
  int iu = is;
  int jl = js;
  int ju = je;
  int kl = ks;
  int ku = ke;
  AthenaArray<Real> *pprim;
  if (&cons == &pmb->phydro->u)
    pprim = &pmb->phydro->w;
  else if (&cons == &pmb->phydro->u1)
    pprim = &pmb->phydro->w1;
  else
    assert(0);
  AthenaArray<Real> g, gi;
  g.NewAthenaArray(NMETRIC,iu+1);
  gi.NewAthenaArray(NMETRIC,iu+1);
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
    {
      pco->CellMetric(k, j, il, iu, g, gi);
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
        pco->LowerVectorCell(u0_new, u1_new, u2_new, u3_new, k, j, i, &u_0_new,
            &u_1_new, &u_2_new, &u_3_new);
        Real gamma_adi = pmb->phydro->peos->GetGamma();
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

// Outer boundary condition
// Inputs:
//   pmb: pointer to block
// Outputs:
//   prim: primitive quantities set along outer x1-boundary
// Notes:
//   remains unchanged
void FixedOuter(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &bb, int is, int ie, int js, int je, int ks, int ke)
{
  return;
}
