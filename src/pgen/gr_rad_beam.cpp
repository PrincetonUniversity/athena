//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gr_rad_beam.cpp
//  \brief Problem generator for radiative transport in GR with single beam

// C++ headers
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str(), string

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "../athena.hpp"                   // Real, FaceField, enums
#include "../athena_arrays.hpp"            // AthenaArray
#include "../parameter_input.hpp"          // ParameterInput
#include "../bvals/bvals.hpp"              // BoundaryValues
#include "../coordinates/coordinates.hpp"  // Coordinates

// Configuration checking
#if not RADIATION_ENABLED
#error "This problem generator must be used with radiation"
#endif
#if not GENERAL_RELATIVITY
#error "This problem generator must be used with general relativity"
#endif

// Declarations
void FixedBoundary(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke,
    int ngh);

// Global variables
static Real mass, spin;           // black hole parameters
static Real pos_1, pos_2, pos_3;  // coordinates of beam origin
static Real width;                // full proper diameter of beam
static Real dir_1, dir_2, dir_3;  // relative direction of beam center
static Real spread;               // full spread of beam in direction

//----------------------------------------------------------------------------------------
// Function for preparing Mesh
// Inputs:
//   pin: input parameters (unused)
// Outputs: (none)

void Mesh::InitUserMeshData(ParameterInput *pin) {

  // Read parameters from input file
  mass = pin->GetReal("coord", "m");
  spin = pin->GetReal("coord", "a");
  pos_1 = pin->GetReal("problem", "pos_1");
  pos_2 = pin->GetReal("problem", "pos_2");
  pos_3 = pin->GetReal("problem", "pos_3");
  width = pin->GetReal("problem", "width");
  dir_1 = pin->GetReal("problem", "dir_1");
  dir_2 = pin->GetReal("problem", "dir_2");
  dir_3 = pin->GetReal("problem", "dir_3");
  spread = pin->GetReal("problem", "spread");

  // Enroll boundary functions
  EnrollUserBoundaryFunction(INNER_X1, FixedBoundary);
  EnrollUserBoundaryFunction(OUTER_X1, FixedBoundary);
  return;
}

//----------------------------------------------------------------------------------------
// Function for preparing MeshBlock
// Inputs:
//   pin: input parameters (unused)
// Outputs: (none)

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(1);
  SetUserOutputVariableName(0, "E");
  return;
}

//----------------------------------------------------------------------------------------
// Function for setting initial conditions
// Inputs:
//   pin: parameters
// Outputs: (none)

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  return;
}

//----------------------------------------------------------------------------------------
// Fixed boundary condition
// Inputs:
//   pmb: pointer to MeshBlock (not used)
//   pcoord: pointer to Coordinates (not used)
//   is,ie,js,je,ks,ke: indices demarkating active region (not used)
// Outputs:
//   prim: primitives set in ghost zones (not used)
//   bb: face-centered magnetic field set in ghost zones (not used)
//   pmb->prad: intensity field set in ghost zones (not used)
// Notes:
//   leaves all quantities as in initial state

void FixedBoundary(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke,
    int ngh) {
  return;
}
