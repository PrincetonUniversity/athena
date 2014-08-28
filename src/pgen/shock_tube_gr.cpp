// General relativistic shock tube generator

// Primary header
#include "../fluid.hpp"

// C++ headers
#include <iostream>   // endl
#include <cmath>      // sqrt()
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena headers
#include "../athena.hpp"           // enums, Real
#include "../athena_arrays.hpp"    // AthenaArray
#include "../mesh.hpp"             // MeshBlock
#include "../parameter_input.hpp"  // ParameterInput

// Declarations
static void set_state(AthenaArray<Real> &prim, AthenaArray<Real> &prim_half,
    AthenaArray<Real> &cons, int i, int j, int k, Real rho, Real pgas, Real vx, Real vy,
    Real vz, Real gamma_adi, Real gamma_adi_red);

// Function for setting initial conditions
// Inputs:
//   pin: parameters
// Outputs: (none)
// Notes:
//   sets conserved variables according to input primitives
void Fluid::InitProblem(ParameterInput *pin)
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
  gamma_ = pin->GetReal("fluid", "gamma");
  Real gamma_adi = GetGamma();
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

  // Read right state
  Real rho_right = pin->GetReal("problem", "dr");
  Real pgas_right = pin->GetReal("problem", "pr");
  Real v1_right = pin->GetReal("problem", "ur");
  Real v2_right = pin->GetReal("problem", "vr");
  Real v3_right = pin->GetReal("problem", "wr");

  // Initialize the discontinuity
  for (int k = kl; k <= ku; k++)
    for (int j = jl; j <= ju; j++)
      for (int i = il; i <= iu; i++)
      {
        bool left_side = false;
        switch(shock_dir)
        {
          case 1:
            if (pb->x1v(i) < shock_pos)
              left_side = true;
            break;
          case 2:
            if (pb->x2v(j) < shock_pos)
              left_side = true;
            break;
          case 3:
            if (pb->x3v(k) < shock_pos)
              left_side = true;
        }
        if (left_side)
          set_state(w, w1, u, i, j, k, rho_left, pgas_left, v1_left, v2_left, v3_left,
              gamma_adi, gamma_adi_red);
        else
          set_state(w, w1, u, i, j, k, rho_right, pgas_right, v1_right, v2_right, v3_right,
              gamma_adi, gamma_adi_red);
      }
  return;
}

// Function for setting conserved variables in a cell given the primitives
// TODO: only works for Minkowski Cartesian metric
static void set_state(AthenaArray<Real> &prim, AthenaArray<Real> &prim_half,
    AthenaArray<Real> &cons, int i, int j, int k, Real rho, Real pgas, Real vx, Real vy,
    Real vz, Real gamma_adi, Real gamma_adi_red)
{
  prim(IDN,k,j,i) = prim_half(IDN,k,j,i) = rho;
  prim(IEN,k,j,i) = prim_half(IEN,k,j,i) = pgas;
  prim(IM1,k,j,i) = prim_half(IM1,k,j,i) = vx;
  prim(IM2,k,j,i) = prim_half(IM2,k,j,i) = vy;
  prim(IM3,k,j,i) = prim_half(IM3,k,j,i) = vz;
  Real gamma_lor_sq = 1.0 / (1.0 - (vx*vx + vy*vy + vz*vz));
  Real gamma_lor_sq_rho_h = gamma_lor_sq * (rho + gamma_adi_red * pgas);
  cons(IDN,k,j,i) = std::sqrt(gamma_lor_sq) * rho;
  cons(IEN,k,j,i) = -gamma_lor_sq_rho_h + pgas;
  cons(IM1,k,j,i) = gamma_lor_sq_rho_h * vx;
  cons(IM2,k,j,i) = gamma_lor_sq_rho_h * vy;
  cons(IM3,k,j,i) = gamma_lor_sq_rho_h * vz;
  return;
}
