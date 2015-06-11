// Uniform advection generator for GRMHD in flat spacetime

// Primary header
#include "../mesh.hpp"

// C++ headers
#include <cmath>  // sqrt()

// Athena headers
#include "../athena.hpp"                   // enums, Real
#include "../athena_arrays.hpp"            // AthenaArray
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../field/field.hpp"              // Field
#include "../fluid/fluid.hpp"              // Fluid
#include "../fluid/eos/eos.hpp"            // GetGamma()
#include "../parameter_input.hpp"          // ParameterInput

// Function for setting initial conditions
// Inputs:
//   pin: parameters
// Outputs:
//   pfl: fluid primitives, half-timestep primitives, and conserved variables set 
//   pfd: magnetic field set
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
  Real gamma_adi = pfl->pf_eos->GetGamma();
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

  // Prepare auxiliary array
  int ncells1 = pfl->pmy_block->block_size.nx1 + 2*NGHOST;
  int ncells2 = pfl->pmy_block->block_size.nx2;
  if (ncells2 > 1)
    ncells2 += 2*NGHOST;
  int ncells3 = pfl->pmy_block->block_size.nx3;
  if (ncells3 > 1)
    ncells3 += 2*NGHOST;
  AthenaArray<Real> b;
  b.NewAthenaArray(3,ncells3,ncells2,ncells1);

  // Initialize hydro variables
  for (int k = kl; k <= ku; k++)
    for (int j = jl; j <= ju; j++)
      for (int i = il; i <= iu; i++)
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
        pb->pcoord->TransformVectorCell(ut, ux, uy, uz, k, j, i, &u0, &u1, &u2, &u3);
        pb->pcoord->TransformVectorCell(bcont, bconx, bcony, bconz, k, j, i,
            &bcon0, &bcon1, &bcon2, &bcon3);

        // Set primitives
        pfl->w(IDN,k,j,i) = pfl->w1(IDN,k,j,i) = rho;
        pfl->w(IEN,k,j,i) = pfl->w1(IEN,k,j,i) = pgas;
        pfl->w(IVX,k,j,i) = pfl->w1(IM1,k,j,i) = u1 / u0;
        pfl->w(IVY,k,j,i) = pfl->w1(IM2,k,j,i) = u2 / u0;
        pfl->w(IVZ,k,j,i) = pfl->w1(IM3,k,j,i) = u3 / u0;

        // Store cell-centered magnetic fields
        b(IB1,k,j,i) = bcon1 * u0 - bcon0 * u1;
        b(IB2,k,j,i) = bcon2 * u0 - bcon0 * u2;
        b(IB3,k,j,i) = bcon3 * u0 - bcon0 * u3;
      }
  pb->pfluid->pf_eos->PrimitiveToConserved(pfl->w, b, pfl->u);  

  // Delete auxiliary array
  b.DeleteAthenaArray();

  // Initialize magnetic field
  if (MAGNETIC_FIELDS_ENABLED)
    for (int k = kl; k <= ku+1; k++)
      for (int j = jl; j <= ju+1; j++)
        for (int i = il; i <= iu+1; i++)
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
            pb->pcoord->TransformVectorFace1(ut, ux, uy, uz, k, j, i,
                &u0, &u1, &u2, &u3);
            pb->pcoord->TransformVectorFace1(bcont, bconx, bcony, bconz, k, j, i,
                &bcon0, &bcon1, &bcon2, &bcon3);
            pfd->b.x1f(k,j,i) = bcon1 * u0 - bcon0 * u1;
          }
          if (i != iu+1 && k != ku+1)
          {
            pb->pcoord->TransformVectorFace2(ut, ux, uy, uz, k, j, i,
                &u0, &u1, &u2, &u3);
            pb->pcoord->TransformVectorFace2(bcont, bconx, bcony, bconz, k, j, i,
                &bcon0, &bcon1, &bcon2, &bcon3);
            pfd->b.x2f(k,j,i) = bcon2 * u0 - bcon0 * u2;
          }
          if (i != iu+1 && j != ju+1)
          {
            pb->pcoord->TransformVectorFace3(ut, ux, uy, uz, k, j, i,
                &u0, &u1, &u2, &u3);
            pb->pcoord->TransformVectorFace3(bcont, bconx, bcony, bconz, k, j, i,
                &bcon0, &bcon1, &bcon2, &bcon3);
            pfd->b.x3f(k,j,i) = bcon3 * u0 - bcon0 * u3;
          }
        }
  return;
}
