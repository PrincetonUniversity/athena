// General relativistic shock tube generator

// Primary header
#include "../mesh.hpp"

// C++ headers
#include <iostream>   // endl
#include <cmath>      // sqrt()
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena headers
#include "../athena.hpp"           // enums, Real
#include "../athena_arrays.hpp"    // AthenaArray
#include "../field/field.hpp"      // Field
#include "../fluid/fluid.hpp"      // Fluid
#include "../fluid/eos/eos.hpp"    // GetGamma()
#include "../parameter_input.hpp"  // ParameterInput

// Declarations
static void SetPrimCons(AthenaArray<Real> &prim, AthenaArray<Real> &prim_half,
    AthenaArray<Real> &cons, int i, int j, int k, Real rho, Real pgas, Real vx, Real vy,
    Real vz, Real bx, Real by, Real bz, Real gamma_adi, Real gamma_adi_red);

// Function for setting initial conditions
// Inputs:
//   pfl: Fluid
//   pfd: Field (unused)
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

  // Read and check shock direction and position
  int shock_dir = pin->GetInteger("problem", "shock_dir"); 
  Real shock_pos = pin->GetReal("problem", "xshock"); 
  Real min_bound, max_bound;
  std::stringstream msg;
  switch (shock_dir)
  {
    case 1:
      min_bound = pb->pmy_domain->pmy_mesh->mesh_size.x1min;
      max_bound = pb->pmy_domain->pmy_mesh->mesh_size.x1max;
      break;
    case 2:
      min_bound = pb->pmy_domain->pmy_mesh->mesh_size.x2min;
      max_bound = pb->pmy_domain->pmy_mesh->mesh_size.x2max;
      break;
    case 3:
      min_bound = pb->pmy_domain->pmy_mesh->mesh_size.x3min;
      max_bound = pb->pmy_domain->pmy_mesh->mesh_size.x3max;
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
  Real v1_left = pin->GetReal("problem", "ul");
  Real v2_left = pin->GetReal("problem", "vl");
  Real v3_left = pin->GetReal("problem", "wl");
  Real b1_left = 0.0, b2_left = 0.0, b3_left = 0.0;
  if (MAGNETIC_FIELDS_ENABLED)
  {
    b1_left = pin->GetReal("problem", "bxl");
    b2_left = pin->GetReal("problem", "byl");
    b3_left = pin->GetReal("problem", "bzl");
  }

  // Read right state
  Real rho_right = pin->GetReal("problem", "dr");
  Real pgas_right = pin->GetReal("problem", "pr");
  Real v1_right = pin->GetReal("problem", "ur");
  Real v2_right = pin->GetReal("problem", "vr");
  Real v3_right = pin->GetReal("problem", "wr");
  Real b1_right = 0.0, b2_right = 0.0, b3_right = 0.0;
  if (MAGNETIC_FIELDS_ENABLED)
  {
    b1_right = pin->GetReal("problem", "bxr");
    b2_right = pin->GetReal("problem", "byr");
    b3_right = pin->GetReal("problem", "bzr");
  }

  // Initialize the discontinuity
  for (int k = kl; k <= ku; k++)
    for (int j = jl; j <= ju; j++)
      for (int i = il; i <= iu; i++)
      {
        bool left_side = false;
        switch(shock_dir)
        {
          case 1:
            left_side = pb->x1v(i) < shock_pos;
            break;
          case 2:
            left_side = pb->x2v(j) < shock_pos;
            break;
          case 3:
            left_side = pb->x3v(k) < shock_pos;
            break;
        }
        // TODO: should field locations be determined by x1f(i) instead of x1v(i)?
        if (left_side)
        {
          SetPrimCons(pfl->w, pfl->w1, pfl->u, i, j, k, rho_left, pgas_left, v1_left,
              v2_left, v3_left, b1_left, b2_left, b3_left, gamma_adi, gamma_adi_red);
          if (MAGNETIC_FIELDS_ENABLED)
          {
            pfd->b.x1f(k,j,i) = b1_left;
            pfd->b.x2f(k,j,i) = b2_left;
            pfd->b.x3f(k,j,i) = b3_left;
          }
        }
        else
        {
          SetPrimCons(pfl->w, pfl->w1, pfl->u, i, j, k, rho_right, pgas_right, v1_right,
              v2_right, v3_right, b1_right, b2_right, b3_right, gamma_adi,
              gamma_adi_red);
          if (MAGNETIC_FIELDS_ENABLED)
          {
            pfd->b.x1f(k,j,i) = b1_right;
            pfd->b.x2f(k,j,i) = b2_right;
            pfd->b.x3f(k,j,i) = b3_right;
          }
        }
      }

  // Add magnetic fields at end faces
  if (MAGNETIC_FIELDS_ENABLED)
  {
    for (int k = kl; k <= ku; k++)
      for (int j = jl; j <= ju; j++)
        pfd->b.x1f(k,j,iu+1) = pfd->b.x1f(k,j,iu);
    for (int k = kl; k <= ku; k++)
      for (int i = il; i <= iu; i++)
        pfd->b.x2f(k,ju+1,i) = pfd->b.x2f(k,ju,i);
    for (int j = jl; j <= ju; j++)
      for (int i = il; i <= iu; i++)
        pfd->b.x3f(ku+1,j,i) = pfd->b.x3f(ku,j,i);
  }
  return;
}

// Function for setting conserved variables in a cell given the primitives
// TODO: only works for Minkowski Cartesian metric
static void SetPrimCons(AthenaArray<Real> &prim, AthenaArray<Real> &prim_half,
    AthenaArray<Real> &cons, int i, int j, int k, Real rho, Real pgas, Real vx, Real vy,
    Real vz, Real bx, Real by, Real bz, Real gamma_adi, Real gamma_adi_red)
{
  // Set primitives
  prim(IDN,k,j,i) = prim_half(IDN,k,j,i) = rho;
  prim(IEN,k,j,i) = prim_half(IEN,k,j,i) = pgas;
  prim(IM1,k,j,i) = prim_half(IM1,k,j,i) = vx;
  prim(IM2,k,j,i) = prim_half(IM2,k,j,i) = vy;
  prim(IM3,k,j,i) = prim_half(IM3,k,j,i) = vz;

  // Calculate intermediate quantites
  Real ut = std::sqrt(1.0 / (1.0 - (SQR(vx)+SQR(vy)+SQR(vz))));
  Real ux = ut * vx;
  Real uy = ut * vy;
  Real uz = ut * vz;
  Real bcovt = bx*ux + by*uy + bz*uz;
  Real bcovx = (bx + bcovt * ux) / ut;
  Real bcovy = (by + bcovt * uy) / ut;
  Real bcovz = (bz + bcovt * uz) / ut;
  Real bcov_sq = -SQR(bcovt) + SQR(bcovx) + SQR(bcovy) + SQR(bcovz);
  Real rho_h = rho + gamma_adi_red * pgas;
  Real ptot = pgas + 0.5*bcov_sq;

  // Set conserved quantities
  cons(IDN,k,j,i) = rho * ut;
  cons(IEN,k,j,i) = -(rho_h + bcov_sq) * ut * ut + bcovt * bcovt + ptot;
  cons(IM1,k,j,i) = (rho_h + bcov_sq) * ut * ux - bcovt * bcovx;
  cons(IM2,k,j,i) = (rho_h + bcov_sq) * ut * uy - bcovt * bcovy;
  cons(IM3,k,j,i) = (rho_h + bcov_sq) * ut * uz - bcovt * bcovz;
  return;
}
