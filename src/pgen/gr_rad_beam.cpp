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
void Source(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
    AthenaArray<Real> &cons);

// Global variables
static Real pos_1, pos_2, pos_3;  // coordinates of beam origin
static Real width;                // full proper diameter of beam
static Real dir_1, dir_2, dir_3;  // relative direction of beam center
static Real spread;               // full spread of beam in direction
static Real dii_dt;               // injected I per unit time
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
  dii_dt = pin->GetReal("problem", "dii_dt");
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
  } else if (COORDINATE_SYSTEM == std::string("kerr-schild")) {
    cylindrical = false;
    spherical = true;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in problem generator\n";
    msg << "unsupported coordinate system\n";
    throw std::runtime_error(msg.str().c_str());
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
  AllocateRealUserMeshBlockDataField(1);
  ruser_meshblock_data[0].NewAthenaArray(prad->nang, ke + 1, je + 1, ie + 1);

  // Calculate beam pattern
  prad->CalculateBeamSource(pos_1, pos_2, pos_3, width, dir_1, dir_2, dir_3, spread,
      dii_dt, ruser_meshblock_data[0], cylindrical, spherical);

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
// Notes:
//   initializes vacuum

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = prad->AngleInd(l, m);
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
          for (int i = is; i <= ie; ++i) {
            prad->prim(lm,k,j,i) = 0.0;
            prad->cons(lm,k,j,i) = 0.0;
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
  AthenaArray<Real> &dcons_dt = pmb->ruser_meshblock_data[0];
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
            cons_out(lm,k,j,i) += dcons_dt(lm,k,j,i) * dt;
          }
        }
      }
    }
  }
  return;
}
