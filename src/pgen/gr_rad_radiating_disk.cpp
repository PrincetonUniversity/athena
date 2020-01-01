//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gr_rad_beam.cpp
//  \brief Problem generator for GR radiation from thin disk

//  C++ headers
#include <cmath>      // sqrt
#include <cstdlib>    // exit (needed for defs.hpp)
#include <iostream>   // cout (needed for defs.hpp), endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error (needed for defs.hpp)
#include <string>     // string (needed for defs.hpp)

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "../athena.hpp"                   // Real, enums, FaceField
#include "../athena_arrays.hpp"            // AthenaArray
#include "../parameter_input.hpp"          // ParameterInput
#include "../bvals/bvals.hpp"              // GetBoundaryFlag
#include "../bvals/bvals_interfaces.hpp"   // BoundaryFace, BoundaryFlag
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../radiation/radiation.hpp"      // Radiation

// Configuration checking
#if not RADIATION_ENABLED
#error "This problem generator must be used with radiation"
#endif
#if not GENERAL_RELATIVITY
#error "This problem generator must be used with general relativity"
#endif
#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator does not support magnetic fields"
#endif

// Declarations
void FixedBoundary(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, AthenaArray<Real> &prim_rad, Real time, Real dt, int il, int iu,
    int jl, int ju, int kl, int ku, int ngh);
Real ThetaGrid(Real x, RegionSize rs);
void ThinDiskSource(MeshBlock *pmb, const Real time, const Real dt,
    const AthenaArray<Real> &prim, AthenaArray<Real> &cons);

// File variables
namespace {
Real m, a;         // mass and spin of black hole
Real h_grid;       // grid compression parameter
Real m_msun;       // mass of black hole in solar masses
Real alpha;        // viscosity parameter
Real r_in, r_out;  // radial limits of radiating disk
Real e_cgs;        // code unit of energy density, in erg/cm^3
int j_min, j_max;  // limits on theta-index for disk
Real zs, ze;       // index bounds on zeta
Real ps, pe;       // index bounds on psi
}  // namespace

//----------------------------------------------------------------------------------------
// Function for preparing Mesh
// Inputs:
//   pin: input parameters
// Outputs: (none)

void Mesh::InitUserMeshData(ParameterInput *pin) {

  // Read parameters from input file
  m = pin->GetReal("coord", "m");
  a = pin->GetReal("coord", "a");
  h_grid = pin->GetOrAddReal("coord", "h", 1.0);
  m_msun = pin->GetReal("problem", "m_msun");
  alpha = pin->GetReal("problem", "alpha");
  r_in = pin->GetReal("problem", "r_in");
  r_out = pin->GetReal("problem", "r_out");
  e_cgs = pin->GetReal("problem", "e_cgs");
  int nx2 = pin->GetInteger("mesh", "nx2");
  j_min = (nx2 - 1) / 2 + NGHOST;
  j_max = nx2 / 2 + NGHOST;
  zs = NGHOST;
  ze = zs + pin->GetInteger("radiation", "n_polar");
  ps = NGHOST;
  pe = ps + pin->GetInteger("radiation", "n_azimuthal");

  // Check for single block
  if (nrbx1 != 1 or nrbx2 != 1 or nrbx3 != 1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Problem Generator\n"
        << "Must use only 1 MeshBlock" << std::endl;
    ATHENA_ERROR(msg);
    return;
  }

  // Enroll user-defined functions
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, FixedBoundary);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, FixedBoundary);
  EnrollUserMeshGenerator(X2DIR, ThetaGrid);
  EnrollUserExplicitRadSourceFunction(ThinDiskSource);
  return;
}

//----------------------------------------------------------------------------------------
// Function for preparing MeshBlock
// Inputs:
//   pin: input parameters (unused)
// Outputs: (none)
// Notes:
//   Uses inner disk formula for temperature from Novikov & Thorne 1973.

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {

  // Allocate persistent arrays
  AllocateRealUserMeshBlockDataField(2);
  ruser_meshblock_data[0].NewAthenaArray(prad->nang, ke + 1, je + 1, ie + 1);
  ruser_meshblock_data[1].NewAthenaArray(prad->nang, ke + 1, je + 1, ie + 1);

  // Calculate primitive intensity as a function of r
  for (int i = is; i <= ie; ++i) {
    Real r = pcoord->x1v(i);
    Real erad = 0.0;
    if (r >= r_in and r <= r_out) {
      Real r_1_4 = std::sqrt(std::sqrt(r));
      Real r_3_8 = r_1_4 * std::sqrt(r_1_4);
      Real aa = 1 + SQR(a) / SQR(r) + 2.0 * m * SQR(a) / (r * SQR(r));
      Real bb = 1.0 + std::sqrt(m) * a / (r * std::sqrt(r));
      Real ee = 1.0 + 4.0 * SQR(a) / SQR(r) - 4.0 * m * SQR(a) / (r * SQR(r))
          + 3.0 * SQR(SQR(a)) / SQR(SQR(r));
      Real tt_cgs = 4.0e7 / std::sqrt(std::sqrt(alpha))
          / std::sqrt(std::sqrt(m_msun / 3.0)) / r_3_8 / std::sqrt(aa) * std::sqrt(bb)
          * std::sqrt(std::sqrt(ee));
      Real erad_cgs = prad->arad_cgs * SQR(SQR(tt_cgs));
      erad = erad_cgs / e_cgs;
    }
    for (int n = 0; n < prad->nang; ++n) {
      for (int k = ks; k <= ke; ++k) {
        for (int j = j_min; j <= j_max; ++j) {
          ruser_meshblock_data[0](n,k,j,i) = erad / (4.0*PI);
        }
      }
    }
  }

  // Calculate disk source (conserved intensity as a function of r)
  prad->PrimitiveToConserved(ruser_meshblock_data[0], ruser_meshblock_data[1], pcoord, is,
      ie, j_min, j_max, ks, ke);
  return;
}

//----------------------------------------------------------------------------------------
// Function for setting initial conditions
// Inputs:
//   pin: parameters (unused)
// Outputs: (none)
// Notes:
//   Initializes constant state with no radiation.

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
// Fixed boundary
// Inputs:
//   pmb: pointer to MeshBlock (not used)
//   pcoord: pointer to Coordinates (not used)
//   time,dt: current time and timestep of simulation (not used)
//   is,ie,js,je,ks,ke: indices demarkating active region (not used)
// Outputs:
//   prim: primitives set in ghost zones
//   bb: face-centered magnetic field set in ghost zones
//   prim_rad: radiation primitives set in ghost zones
// Notes:
//   Does nothing.

void FixedBoundary(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, AthenaArray<Real> &prim_rad, Real time, Real dt, int il, int iu,
    int jl, int ju, int kl, int ku, int ngh) {
  return;
}

//----------------------------------------------------------------------------------------
// Theta-spacing function
// Inputs:
//   x2: internal theta-like coordinate, scaled from 0 to 1
//   rs: region size struct
// Outputs:
//   returned value: corresponding theta-value
// Notes:
//   Implements remapping from Gammie, McKinney, & Toth 2003, ApJ 589 444.

Real ThetaGrid(Real x2, RegionSize rs) {
  return PI * x2 + (1.0 - h_grid) / 2.0 * std::sin(2.0*PI * x2);
}

//----------------------------------------------------------------------------------------
// Source function for illuminating thin disk
// Inputs:
//   pmb: pointer to MeshBlock
//   time: time of simulation (unused)
//   dt: simulation timestep (unused)
//   prim_rad: primitive intensity (unused)
//   cons_rad: conserved intensity (unused)
// Outputs:
//   cons_rad: conserved intensity updated
// Notes:
//   Applies pre-computed conserved intensity field in midplane.

void ThinDiskSource(MeshBlock *pmb, const Real time, const Real dt,
    const AthenaArray<Real> &prim_rad, AthenaArray<Real> &cons_rad) {

  // Extract information from block
  Radiation *prad = pmb->prad;
  AthenaArray<Real> &cons = pmb->ruser_meshblock_data[1];
  int is = pmb->is;
  int ie = pmb->ie;
  int ks = pmb->ks;
  int ke = pmb->ke;

  // Overwrite conserved intensity with stored values
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = prad->AngleInd(l, m);
      for (int k = ks; k <= ke; ++k) {
        for (int j = j_min; j <= j_max; ++j) {
          for (int i = is; i <= ie; ++i) {
            cons_rad(lm,k,j,i) = cons(lm,k,j,i);
          }
        }
      }
    }
  }
  return;
}
