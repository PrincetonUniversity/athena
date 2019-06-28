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
#include <string>     // c_str, string

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "../athena.hpp"                   // Real, FaceField, enums
#include "../athena_arrays.hpp"            // AthenaArray
#include "../parameter_input.hpp"          // ParameterInput
#include "../bvals/bvals.hpp"              // BoundaryValues
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../radiation/radiation.hpp"      // Radiation

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
void Source(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
    AthenaArray<Real> &cons);

// Global variables
static Real pos_1, pos_2, pos_3;  // coordinates of beam origin
static Real width;                // full proper diameter of beam
static Real dir_1, dir_2, dir_3;  // relative direction of beam center
static Real spread;               // full spread of beam in direction
static Real zs, ze;               // index bounds on zeta
static Real ps, pe;               // index bounds on psi
static bool cylindrical;          // flag indicating cylindrical coordinates
static bool spherical;            // flag indicating spherical coordinates

//----------------------------------------------------------------------------------------
// Function for preparing Mesh
// Inputs:
//   pin: input parameters (unused)
// Outputs: (none)

void Mesh::InitUserMeshData(ParameterInput *pin) {

  // Read parameters from input file
  pos_1 = pin->GetReal("problem", "pos_1");
  pos_2 = pin->GetReal("problem", "pos_2");
  pos_3 = pin->GetReal("problem", "pos_3");
  width = pin->GetReal("problem", "width");
  dir_1 = pin->GetReal("problem", "dir_1");
  dir_2 = pin->GetReal("problem", "dir_2");
  dir_3 = pin->GetReal("problem", "dir_3");
  spread = pin->GetReal("problem", "spread");
  zs = NGHOST;
  ze = zs + pin->GetInteger("radiation", "n_polar");
  ps = NGHOST;
  pe = ps + pin->GetInteger("radiation", "n_azimuthal");

  // Determine coordinate type
  if (COORDINATE_SYSTEM == std::string("minkowski")) {
    cylindrical = false;
    spherical = false;
  } else if (COORDINATE_SYSTEM == std::string("minkowski_cyl")) {
    cylindrical = true;
    spherical = false;
  } else if (COORDINATE_SYSTEM == std::string("minkowski_sph")) {
    cylindrical = false;
    spherical = true;
  } else if (COORDINATE_SYSTEM == std::string("schwarzschild")) {
    cylindrical = false;
    spherical = true;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in problem generator\n";
    msg << "unsupported coordinate system\n";
    throw std::runtime_error(msg.str().c_str());
  }

  // Enroll boundary functions
  if (pin->GetString("mesh", "ix1_bc") == "user") {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, FixedBoundary);
  }
  if (pin->GetString("mesh", "ox1_bc") == "user") {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, FixedBoundary);
  }

  // Enroll source function
  EnrollUserExplicitRadSourceFunction(Source);
  return;
}

//----------------------------------------------------------------------------------------
// Function for preparing MeshBlock
// Inputs:
//   pin: input parameters (unused)
// Outputs: (none)

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {

  // Allocate persistent arrays
  AllocateRealUserMeshBlockDataField(2);
  ruser_meshblock_data[0].NewAthenaArray(prad->nang, ke + 1, je + 1, ie + 1);
  ruser_meshblock_data[1].NewAthenaArray(prad->nang, ke + 1, je + 1, ie + 1);

  // Prepare output variables
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

  // Calculate beam pattern
  prad->CalculateBeamSource(pos_1, pos_2, pos_3, width, dir_1, dir_2, dir_3, spread,
      ruser_meshblock_data[0], ruser_meshblock_data[1], cylindrical, spherical);

  // Set primitive and conserved values
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = prad->AngleInd(l, m);
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
          for (int i = is; i <= ie; ++i) {
            prad->prim(lm,k,j,i) = ruser_meshblock_data[0](lm,k,j,i);
            prad->cons(lm,k,j,i) = ruser_meshblock_data[1](lm,k,j,i);
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for preparing output
// Inputs:
//   pin: parameters
// Outputs: (none)
// Notes:
//   sets user_out_var array

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  prad->SetMoments(user_out_var);
  return;
}

//----------------------------------------------------------------------------------------
// Fixed boundary condition
// Inputs:
//   pmb: pointer to MeshBlock (not used)
//   pcoord: pointer to Coordinates (not used)
//   time: time of simulation
//   dt: simulation timestep
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

//----------------------------------------------------------------------------------------
// Source function
// Inputs:
//   pmb: pointer to MeshBlock (not used)
//   time: time of simulation
//   dt: simulation timestep
//   prim_in: primitive intensity
// Outputs:
//   cons_out: conserved intensity
// Notes:
//   applies pre-computed conserved intensity field
//   source is opaque; intensities are replaced, not added to

void Source(MeshBlock *pmb, const Real time, const Real dt,
    const AthenaArray<Real> &prim_in, AthenaArray<Real> &cons_out) {

  // Extract information from block
  Radiation *prad = pmb->prad;
  AthenaArray<Real> &cons = pmb->ruser_meshblock_data[1];
  int is = pmb->is;
  int ie = pmb->ie;
  int js = pmb->js;
  int je = pmb->je;
  int ks = pmb->ks;
  int ke = pmb->ke;

  // Overwrite conserved intensity with stored initial values
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = prad->AngleInd(l, m);
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
          for (int i = is; i <= ie; ++i) {
            Real val = cons(lm,k,j,i);
            if (val != 0.0) {
              cons_out(lm,k,j,i) = val;
            }
          }
        }
      }
    }
  }
  return;
}
