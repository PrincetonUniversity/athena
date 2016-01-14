// Relativistic shock tube generator

// Primary header
#include "../mesh.hpp"

// C++ headers
#include <cmath>      // sqrt()
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena headers
#include "../athena.hpp"                   // macros, enums, Real
#include "../athena_arrays.hpp"            // AthenaArray
#include "../parameter_input.hpp"          // ParameterInput
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../field/field.hpp"              // Field
#include "../hydro/hydro.hpp"
#include "../hydro/eos/eos.hpp"


/// Function for initializing global mesh properties
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
// Outputs: (none)
// Notes:
//   sets conserved variables according to input primitives
//   assigns fields based on cell-center positions, rather than interface positions
//     this helps shock tube 2 from Mignone, Ugliano, & Bodo 2009, MNRAS 393 1141
//     otherwise the middle interface would go to left variables, creating a
//         particularly troublesome jump leading to NaN's
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

  // Read and set ratio of specific heats
  Real gamma_adi = phydro->peos->GetGamma();
  Real gamma_adi_red = gamma_adi / (gamma_adi - 1.0);

  // Read and check shock direction and position
  int shock_dir = pin->GetInteger("problem", "shock_dir"); 
  Real shock_pos = pin->GetReal("problem", "xshock"); 
  Real min_bound, max_bound;
  std::stringstream msg;
  switch (shock_dir)
  {
    case 1:
      min_bound = pmy_mesh->mesh_size.x1min;
      max_bound = pmy_mesh->mesh_size.x1max;
      break;
    case 2:
      min_bound = pmy_mesh->mesh_size.x2min;
      max_bound = pmy_mesh->mesh_size.x2max;
      break;
    case 3:
      min_bound = pmy_mesh->mesh_size.x3min;
      max_bound = pmy_mesh->mesh_size.x3max;
      break;
    default:
      msg << "### FATAL ERROR in Problem Generator" << std::endl
          << "shock_dir=" << shock_dir << " must be either 1, 2, or 3" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  if (shock_pos < min_bound || shock_pos > max_bound)
  {
    msg << "### FATAL ERROR in Problem Generator" << std::endl
        << "xshock=" << shock_pos << " lies outside x" << shock_dir
        << " domain for shkdir=" << shock_dir << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // Read left state
  Real rho_left = pin->GetReal("problem", "dl");
  Real pgas_left = pin->GetReal("problem", "pl");
  Real vx_left = pin->GetReal("problem", "ul");
  Real vy_left = pin->GetReal("problem", "vl");
  Real vz_left = pin->GetReal("problem", "wl");
  Real bx_left = 0.0, by_left = 0.0, bz_left = 0.0;
  if (MAGNETIC_FIELDS_ENABLED)
  {
    bx_left = pin->GetReal("problem", "bxl");
    by_left = pin->GetReal("problem", "byl");
    bz_left = pin->GetReal("problem", "bzl");
  }

  // Read right state
  Real rho_right = pin->GetReal("problem", "dr");
  Real pgas_right = pin->GetReal("problem", "pr");
  Real vx_right = pin->GetReal("problem", "ur");
  Real vy_right = pin->GetReal("problem", "vr");
  Real vz_right = pin->GetReal("problem", "wr");
  Real bx_right = 0.0, by_right = 0.0, bz_right = 0.0;
  if (MAGNETIC_FIELDS_ENABLED)
  {
    bx_right = pin->GetReal("problem", "bxr");
    by_right = pin->GetReal("problem", "byr");
    bz_right = pin->GetReal("problem", "bzr");
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
  if (GENERAL_RELATIVITY)
  {
    g.NewAthenaArray(NMETRIC, ncells1);
    gi.NewAthenaArray(NMETRIC, ncells1);
  }

  // Initialize hydro variables
  for (int k = kl; k <= ku; ++k)
  {
    for (int j = jl; j <= ju; ++j)
    {
      if (GENERAL_RELATIVITY)
        pcoord->CellMetric(k, j, il, iu, g, gi);
      for (int i = il; i <= iu; ++i) 
      {
        // Determine which variables to use
        Real rho = rho_right;
        Real pgas = pgas_right;
        Real vx = vx_right;
        Real vy = vy_right;
        Real vz = vz_right;
        Real bx = bx_right;
        Real by = by_right;
        Real bz = bz_right;
        bool left_side = false;
        switch(shock_dir)
        {
          case 1:
            left_side = pcoord->x1v(i) < shock_pos;
            break;
          case 2:
            left_side = pcoord->x2v(j) < shock_pos;
            break;
          case 3:
            left_side = pcoord->x3v(k) < shock_pos;
            break;
        }
        if (left_side)
        {
          rho = rho_left;
          pgas = pgas_left;
          vx = vx_left;
          vy = vy_left;
          vz = vz_left;
          bx = bx_left;
          by = by_left;
          bz = bz_left;
        }

        // Construct 4-vectors
        Real ut = std::sqrt(1.0 / (1.0 - (SQR(vx)+SQR(vy)+SQR(vz))));
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
        if (GENERAL_RELATIVITY)
        {
          pcoord->TransformVectorCell(ut, ux, uy, uz, k, j, i, &u0, &u1, &u2, &u3);
          pcoord->TransformVectorCell(bcont, bconx, bcony, bconz, k, j, i, &bcon0,
              &bcon1, &bcon2, &bcon3);
        }
        else
        {
          u0 = ut;
          u1 = ux;
          u2 = uy;
          u3 = uz;
          bcon0 = bcont;
          bcon1 = bconx;
          bcon2 = bcony;
          bcon3 = bconz;
        }

        // Set primitives
        phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) = rho;
        phydro->w(IEN,k,j,i) = phydro->w1(IEN,k,j,i) = pgas;
        if (GENERAL_RELATIVITY)
        {
          Real uu1 = u1 - gi(I01,i)/gi(I00,i) * u0;
          Real uu2 = u2 - gi(I02,i)/gi(I00,i) * u0;
          Real uu3 = u3 - gi(I03,i)/gi(I00,i) * u0;
          phydro->w(IVX,k,j,i) = phydro->w1(IVX,k,j,i) = uu1;
          phydro->w(IVY,k,j,i) = phydro->w1(IVY,k,j,i) = uu2;
          phydro->w(IVZ,k,j,i) = phydro->w1(IVZ,k,j,i) = uu3;
        }
        else
        {
          phydro->w(IVX,k,j,i) = phydro->w1(IM1,k,j,i) = u1 / u0;
          phydro->w(IVY,k,j,i) = phydro->w1(IM2,k,j,i) = u2 / u0;
          phydro->w(IVZ,k,j,i) = phydro->w1(IM3,k,j,i) = u3 / u0;
        }

        // Set magnetic fields
        b(IB1,k,j,i) = bcon1 * u0 - bcon0 * u1;
        b(IB2,k,j,i) = bcon2 * u0 - bcon0 * u2;
        b(IB3,k,j,i) = bcon3 * u0 - bcon0 * u3;
      }
    }
  }
  phydro->peos->PrimitiveToConserved(phydro->w, b, phydro->u, pcoord,
                                            il, iu, jl, ju, kl, ju);

  // Delete auxiliary arrays
  b.DeleteAthenaArray();
  if (GENERAL_RELATIVITY)
  {
    g.DeleteAthenaArray();
    gi.DeleteAthenaArray();
  }

  // Initialize magnetic field
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k = kl; k <= ku+1; ++k) {
      for (int j = jl; j <= ju+1; ++j) {
        for (int i = il; i <= iu+1; ++i)
        {
          // Determine which variables to use
          Real vx = vx_right;
          Real vy = vy_right;
          Real vz = vz_right;
          Real bx = bx_right;
          Real by = by_right;
          Real bz = bz_right;
          bool left_side = false;
          switch(shock_dir)
          {
            case 1:
              left_side = pcoord->x1v(i) < shock_pos;
              break;
            case 2:
              left_side = pcoord->x2v(j) < shock_pos;
              break;
            case 3:
              left_side = pcoord->x3v(k) < shock_pos;
              break;
          }
          if (left_side)
          {
            vx = vx_left;
            vy = vy_left;
            vz = vz_left;
            bx = bx_left;
            by = by_left;
            bz = bz_left;
          }

          // Construct 4-vectors
          Real ut = std::sqrt(1.0 / (1.0 - (SQR(vx)+SQR(vy)+SQR(vz))));
          Real ux = ut * vx;
          Real uy = ut * vy;
          Real uz = ut * vz;
          Real bcont = bx*ux + by*uy + bz*uz;
          Real bconx = (bx + bcont * ux) / ut;
          Real bcony = (by + bcont * uy) / ut;
          Real bconz = (bz + bcont * uz) / ut;

          // Set magnetic fields
          Real u0, u1, u2, u3;
          Real bcon0, bcon1, bcon2, bcon3;
          if (j != ju+1 && k != ku+1)
          {
            if (GENERAL_RELATIVITY)
            {
              pcoord->TransformVectorFace1(ut, ux, uy, uz, k, j, i, &u0, &u1, &u2,
                  &u3);
              pcoord->TransformVectorFace1(bcont, bconx, bcony, bconz, k, j, i,
                  &bcon0, &bcon1, &bcon2, &bcon3);
              pfield->b.x1f(k,j,i) = bcon1 * u0 - bcon0 * u1;
            }
            else
              pfield->b.x1f(k,j,i) = bx;
          }
          if (i != iu+1 && k != ku+1)
          {
            if (GENERAL_RELATIVITY)
            {
              pcoord->TransformVectorFace2(ut, ux, uy, uz, k, j, i, &u0, &u1, &u2,
                  &u3);
              pcoord->TransformVectorFace2(bcont, bconx, bcony, bconz, k, j, i,
                  &bcon0, &bcon1, &bcon2, &bcon3);
              pfield->b.x2f(k,j,i) = bcon2 * u0 - bcon0 * u2;
            }
            else
              pfield->b.x2f(k,j,i) = by;
          }
          if (i != iu+1 && j != ju+1)
          {
            if (GENERAL_RELATIVITY)
            {
              pcoord->TransformVectorFace3(ut, ux, uy, uz, k, j, i, &u0, &u1, &u2,
                  &u3);
              pcoord->TransformVectorFace3(bcont, bconx, bcony, bconz, k, j, i,
                  &bcon0, &bcon1, &bcon2, &bcon3);
              pfield->b.x3f(k,j,i) = bcon3 * u0 - bcon0 * u3;
            }
            else
              pfield->b.x3f(k,j,i) = bz;
          }
        }
      }
    }
  }
  return;
}


// User-defined work function called every time step
void MeshBlock::UserWorkInLoop(void)
{
  return;
}

