//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gr_rad_multibeam.cpp
//  \brief Problem generator for radiative transport in GR with two beams

// C++ headers
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str, string

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "../athena.hpp"               // Real, enums
#include "../athena_arrays.hpp"        // AthenaArray
#include "../parameter_input.hpp"      // ParameterInput
#include "../radiation/radiation.hpp"  // Radiation

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
namespace {
Real pos_a_1, pos_a_2, pos_a_3;  // coordinates of beam A origin
Real width_a;                    // full proper diameter of beam A
Real dir_a_1, dir_a_2, dir_a_3;  // relative direction of beam A center
Real spread_a;                   // full spread of beam A in direction
Real dii_dt_a;                   // injected I per unit time for beam A
Real pos_b_1, pos_b_2, pos_b_3;  // coordinates of beam B origin
Real width_b;                    // full proper diameter of beam B
Real dir_b_1, dir_b_2, dir_b_3;  // relative direction of beam B center
Real spread_b;                   // full spread of beam B in direction
Real dii_dt_b;                   // injected I per unit time for beam B
Real zs, ze;                     // index bounds on zeta
Real ps, pe;                     // index bounds on psi
bool cylindrical;                // flag indicating cylindrical coordinates
bool spherical;                  // flag indicating spherical coordinates
}  // namespace

//----------------------------------------------------------------------------------------
// Function for preparing Mesh
// Inputs:
//   pin: input parameters
// Outputs: (none)

void Mesh::InitUserMeshData(ParameterInput *pin) {

  // Read parameters from input file
  pos_a_1 = pin->GetReal("problem", "pos_a_1");
  pos_a_2 = pin->GetReal("problem", "pos_a_2");
  pos_a_3 = pin->GetReal("problem", "pos_a_3");
  width_a = pin->GetReal("problem", "width_a");
  dir_a_1 = pin->GetReal("problem", "dir_a_1");
  dir_a_2 = pin->GetReal("problem", "dir_a_2");
  dir_a_3 = pin->GetReal("problem", "dir_a_3");
  spread_a = pin->GetReal("problem", "spread_a");
  dii_dt_a = pin->GetReal("problem", "dii_dt_a");
  pos_b_1 = pin->GetReal("problem", "pos_b_1");
  pos_b_2 = pin->GetReal("problem", "pos_b_2");
  pos_b_3 = pin->GetReal("problem", "pos_b_3");
  width_b = pin->GetReal("problem", "width_b");
  dir_b_1 = pin->GetReal("problem", "dir_b_1");
  dir_b_2 = pin->GetReal("problem", "dir_b_2");
  dir_b_3 = pin->GetReal("problem", "dir_b_3");
  spread_b = pin->GetReal("problem", "spread_b");
  dii_dt_b = pin->GetReal("problem", "dii_dt_b");
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
  AllocateRealUserMeshBlockDataField(2);
  ruser_meshblock_data[0].NewAthenaArray(prad->nang, ke + 1, je + 1, ie + 1);
  ruser_meshblock_data[1].NewAthenaArray(prad->nang, ke + 1, je + 1, ie + 1);

  // Calculate beam pattern
  prad->CalculateBeamSource(pos_a_1, pos_a_2, pos_a_3, width_a, dir_a_1, dir_a_2, dir_a_3,
      spread_a, dii_dt_a, ruser_meshblock_data[0], cylindrical, spherical);
  prad->CalculateBeamSource(pos_b_1, pos_b_2, pos_b_3, width_b, dir_b_1, dir_b_2, dir_b_3,
      spread_b, dii_dt_b, ruser_meshblock_data[1], cylindrical, spherical);
  for (int lm = 0; lm < prad->nang; ++lm) {
    for (int k = 0; k < ke+1; ++k) {
      for (int j = 0; j < je+1; ++j) {
        for (int i = 0; i < ie+1; ++i) {
          ruser_meshblock_data[0](lm,k,j,i) += ruser_meshblock_data[1](lm,k,j,i);
        }
      }
    }
  }
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
// Source function
// Inputs:
//   pmb: pointer to MeshBlock
//   time: time of simulation
//   dt: simulation timestep
//   prim_rad: primitive intensity
//   cons_rad: conserved intensity
// Outputs:
//   cons_rad: conserved intensity updated
// Notes:
//   applies pre-computed conserved intensity field
//   source is transparent; intensities are added to, not replaced

void Source(MeshBlock *pmb, const Real time, const Real dt,
    const AthenaArray<Real> &prim_rad, AthenaArray<Real> &cons_rad) {

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
            cons_rad(lm,k,j,i) += dcons_dt(lm,k,j,i) * dt;
          }
        }
      }
    }
  }
  return;
}
