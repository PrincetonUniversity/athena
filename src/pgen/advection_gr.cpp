// Uniform advection generator for GRMHD in flat spacetime

// Primary header
#include "../mesh.hpp"

// C++ headers
#include <cmath>  // sqrt()

// Athena headers
#include "../athena.hpp"                   // macros, enums, Real
#include "../athena_arrays.hpp"            // AthenaArray
#include "../parameter_input.hpp"          // ParameterInput
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../field/field.hpp"              // Field
#include "../hydro/hydro.hpp"
#include "../hydro/eos/eos.hpp"


// Function for initializing global mesh properties
void Mesh::InitUserMeshProperties(ParameterInput *pin)
{
  return;
}

// Function for cleaning up global mesh properties
void Mesh::TerminateUserMeshProperties(void)
{
  return;
}


// Function for setting initial conditions
// Inputs:
//   pin: parameters
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

  // Get ratio of specific heats
  Real gamma_adi = phydro->peos->GetGamma();
  Real gamma_adi_red = gamma_adi / (gamma_adi - 1.0);

  // Read problem parameters
  Real rho = pin->GetReal("problem", "rho");
  Real pgas = pin->GetReal("problem", "pgas");
  Real vx = pin->GetReal("problem", "vx");
  Real vy = pin->GetReal("problem", "vy");
  Real vz = pin->GetReal("problem", "vz");
  Real bx = 0.0, by = 0.0, bz = 0.0;
  if (MAGNETIC_FIELDS_ENABLED)
  {
    bx = pin->GetReal("problem", "bx");
    by = pin->GetReal("problem", "by");
    bz = pin->GetReal("problem", "bz");
  }

  // Prepare auxiliary arrays
  int ncells1 = block_size.nx1 + 2*NGHOST;
  int ncells2 = block_size.nx2;
  if (ncells2 > 1)
    ncells2 += 2*NGHOST;
  int ncells3 = block_size.nx3;
  if (ncells3 > 1)
    ncells3 += 2*NGHOST;
  AthenaArray<Real> b, g, gi;
  b.NewAthenaArray(3, ncells3, ncells2, ncells1);
  g.NewAthenaArray(NMETRIC, ncells1);
  gi.NewAthenaArray(NMETRIC, ncells1);

  // Initialize hydro variables
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j)
    {
      pcoord->CellMetric(k, j, il, iu, g, gi);
      for (int i = il; i <= iu; ++i)
      {
        // Construct 4-vectors
        Real ut = 1.0 / std::sqrt(1.0 - SQR(vx) - SQR(vy) - SQR(vz));
        Real ux = ut * vx;
        Real uy = ut * vy;
        Real uz = ut * vz;
        Real bcont = bx*ux + by*uy + bz*uz;
        Real bconx = (bx + bcont * ux) / ut;
        Real bcony = (by + bcont * uy) / ut;
        Real bconz = (bz + bcont * uz) / ut;

        // Transform 4-vectors
        Real u0, u1, u2, u3;
        Real bcon0, bcon1, bcon2, bcon3;
        pcoord->TransformVectorCell(ut, ux, uy, uz, k, j, i, &u0, &u1, &u2, &u3);
        pcoord->TransformVectorCell(bcont, bconx, bcony, bconz, k, j, i, &bcon0,
            &bcon1, &bcon2, &bcon3);

        // Set primitives
        phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) = rho;
        phydro->w(IEN,k,j,i) = phydro->w1(IEN,k,j,i) = pgas;
        phydro->w(IVX,k,j,i) = phydro->w1(IVX,k,j,i) = u1 - gi(I01,i)/gi(I00,i) * u0;
        phydro->w(IVY,k,j,i) = phydro->w1(IVY,k,j,i) = u2 - gi(I02,i)/gi(I00,i) * u0;
        phydro->w(IVZ,k,j,i) = phydro->w1(IVZ,k,j,i) = u3 - gi(I03,i)/gi(I00,i) * u0;

        // Store cell-centered magnetic fields
        b(IB1,k,j,i) = bcon1 * u0 - bcon0 * u1;
        b(IB2,k,j,i) = bcon2 * u0 - bcon0 * u2;
        b(IB3,k,j,i) = bcon3 * u0 - bcon0 * u3;
      }
    }
  }
  phydro->peos->PrimitiveToConserved(phydro->w, b, phydro->u, pcoord,
                                          il, iu, jl, ju, kl, ku);

  // Delete auxiliary arrays
  b.DeleteAthenaArray();
  g.DeleteAthenaArray();
  gi.DeleteAthenaArray();

  // Initialize magnetic field
  if (MAGNETIC_FIELDS_ENABLED)
    for (int k = kl; k <= ku+1; ++k)
      for (int j = jl; j <= ju+1; ++j)
        for (int i = il; i <= iu+1; ++i)
        {
          // Construct 4-vectors
          Real ut = 1.0 / std::sqrt(1.0 - SQR(vx) - SQR(vy) - SQR(vz));
          Real ux = ut * vx;
          Real uy = ut * vy;
          Real uz = ut * vz;
          Real bcont = bx*ux + by*uy + bz*uz;
          Real bconx = (bx + bcont * ux) / ut;
          Real bcony = (by + bcont * uy) / ut;
          Real bconz = (bz + bcont * uz) / ut;

          // Transform 4-vectors and set magnetic fields
          Real u0, u1, u2, u3;
          Real bcon0, bcon1, bcon2, bcon3;
          if (j != ju+1 && k != ku+1)
          {
            pcoord->TransformVectorFace1(ut, ux, uy, uz, k, j, i,
                &u0, &u1, &u2, &u3);
            pcoord->TransformVectorFace1(bcont, bconx, bcony, bconz, k, j, i,
                &bcon0, &bcon1, &bcon2, &bcon3);
            pfield->b.x1f(k,j,i) = bcon1 * u0 - bcon0 * u1;
          }
          if (i != iu+1 && k != ku+1)
          {
            pcoord->TransformVectorFace2(ut, ux, uy, uz, k, j, i,
                &u0, &u1, &u2, &u3);
            pcoord->TransformVectorFace2(bcont, bconx, bcony, bconz, k, j, i,
                &bcon0, &bcon1, &bcon2, &bcon3);
            pfield->b.x2f(k,j,i) = bcon2 * u0 - bcon0 * u2;
          }
          if (i != iu+1 && j != ju+1)
          {
            pcoord->TransformVectorFace3(ut, ux, uy, uz, k, j, i,
                &u0, &u1, &u2, &u3);
            pcoord->TransformVectorFace3(bcont, bconx, bcony, bconz, k, j, i,
                &bcon0, &bcon1, &bcon2, &bcon3);
            pfield->b.x3f(k,j,i) = bcon3 * u0 - bcon0 * u3;
          }
        }
  return;
}



// User-defined work function called every time step
void MeshBlock::UserWorkInLoop(void)
{
  return;
}

