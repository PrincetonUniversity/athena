// General relativistic problem generator for dust falling onto black hole

// Primary header
#include "../mesh.hpp"

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
void FixedOuter(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &bb,
                int is, int ie, int js, int je, int ks, int ke);

// Function for initializing global mesh properties
void Mesh::InitUserMeshData(ParameterInput *pin)
{
  // Enroll boundary functions
  EnrollUserBoundaryFunction(OUTER_X1, FixedOuter);
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
    for (int i = il; i <= iu; i++)
    {
      // Get Boyer-Lindquist coordinates of cell
      Real r, theta, phi;
      pmb->pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i), pcoord->x2v(j),
          pcoord->x3v(kl), &r, &theta, &phi);

      // Calculate primitives depending on location
      Real rho = rho_min * std::pow(r, rho_pow);
      Real pgas = (gamma_adi-1.0) * u_min * std::pow(r, u_pow);
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

  // Initialize conserved values
  AthenaArray<Real> bb;
  bb.NewAthenaArray(3, ku+1, ju+1, iu+1);
  phydro->peos->PrimitiveToConserved(phydro->w, bb, phydro->u, pcoord, il, iu, jl, ju,
      kl, ku);
  bb.DeleteAthenaArray();
  return;
}


// Outer boundary condition
// Inputs:
//   pmb: pointer to block
//   pco: pointer to coordinates
//   is,ie,js,je,ks,ke: index boundaries of active zone
// Outputs:
//   prim: primitive quantities set along outer x1-boundary
//   bb: magnetic fields set along outer x1-boundary
// Notes:
//   remains unchanged
void FixedOuter(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
    FaceField &bb, int is, int ie, int js, int je, int ks, int ke)
{
  return;
}
