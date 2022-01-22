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
void ZeroInner(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, AthenaArray<Real> &prim_rad, Real time, Real dt, int il, int iu,
    int jl, int ju, int kl, int ku, int ngh);
void ZeroOuter(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, AthenaArray<Real> &prim_rad, Real time, Real dt, int il, int iu,
    int jl, int ju, int kl, int ku, int ngh);
Real ThetaGrid(Real x, RegionSize rs);
void ThinDiskSource(MeshBlock *pmb, const Real time, const Real dt,
    const AthenaArray<Real> &prim, AthenaArray<Real> &cons);

// File variables
namespace {
Real mass, spin;   // mass and spin of black hole
Real h_grid;       // grid compression parameter
Real m_msun;       // mass of black hole in solar masses
Real alpha;        // viscosity parameter
Real r_in, r_out;  // radial limits of radiating disk
Real e_cgs;        // code unit of energy density, in erg/cm^3
int root_offset;   // offset in level count
int num_root_th;   // number of blocks in theta-direction at root level
int zs, ze;        // index bounds on zeta
int ps, pe;        // index bounds on psi
}  // namespace

//----------------------------------------------------------------------------------------
// Function for preparing Mesh
// Inputs:
//   pin: input parameters
// Outputs: (none)

void Mesh::InitUserMeshData(ParameterInput *pin) {

  // Read parameters from input file
  mass = pin->GetReal("coord", "m");
  spin = pin->GetReal("coord", "a");
  h_grid = pin->GetOrAddReal("coord", "h", 1.0);
  m_msun = pin->GetReal("problem", "m_msun");
  alpha = pin->GetReal("problem", "alpha");
  r_in = pin->GetReal("problem", "r_in");
  r_out = pin->GetReal("problem", "r_out");
  e_cgs = pin->GetReal("problem", "e_cgs");
  root_offset = root_level;
  int nx2 = pin->GetInteger("mesh", "nx2");
  num_root_th = nx2 / pin->GetOrAddInteger("meshblock", "nx2", nx2);
  zs = NGHOST_RAD;
  ze = zs + pin->GetInteger("radiation", "n_polar");
  ps = NGHOST_RAD;
  pe = ps + pin->GetInteger("radiation", "n_azimuthal");

  // Enroll user-defined functions
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, ZeroInner);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, ZeroOuter);
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
//   Allocate persistent array.

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateRealUserMeshBlockDataField(1);
  ruser_meshblock_data[0].NewAthenaArray(prad->nang, ke + 1, je + 1, ie + 1);
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
          }
        }
      }
    }
  }
  prad->PrimitiveToConserved(prad->prim, prad->cons, pcoord, is, ie, js, je, ks, ke);
  return;
}

//----------------------------------------------------------------------------------------
// Vanishing radiation field inside radial coordinate
// Inputs:
//   pmb: pointer to MeshBlock
//   pcoord: pointer to Coordinates (not used)
//   time, dt: current time and timestep of simulation (not used)
//   il, iu, jl, ju, kl, ku: indices demarkating active region
// Outputs:
//   prim: primitives set in ghost zones (not used)
//   bb: face-centered magnetic field set in ghost zones (not used)
//   prim_rad: radiation primitives set in ghost zones
// Notes:
//   Sets intensity to 0.

void ZeroInner(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, AthenaArray<Real> &prim_rad, Real time, Real dt, int il, int iu,
    int jl, int ju, int kl, int ku, int ngh) {
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = pmb->prad->AngleInd(l, m);
      for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
          for (int i = il - ngh; i <= il - 1; ++i) {
            prim_rad(lm,k,j,i) = 0.0;
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Vanishing radiation field outside radial coordinate
// Inputs:
//   pmb: pointer to MeshBlock
//   pcoord: pointer to Coordinates (not used)
//   time, dt: current time and timestep of simulation (not used)
//   il, iu, jl, ju, kl, ku: indices demarkating active region
// Outputs:
//   prim: primitives set in ghost zones (not used)
//   bb: face-centered magnetic field set in ghost zones (not used)
//   prim_rad: radiation primitives set in ghost zones
// Notes:
//   Sets intensity to 0.

void ZeroOuter(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, AthenaArray<Real> &prim_rad, Real time, Real dt, int il, int iu,
    int jl, int ju, int kl, int ku, int ngh) {
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = pmb->prad->AngleInd(l, m);
      for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
          for (int i = iu + 1; i <= iu + ngh; ++i) {
            prim_rad(lm,k,j,i) = 0.0;
          }
        }
      }
    }
  }
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
  return PI * x2 + (1.0 - h_grid) / 2.0 * std::sin(2.0 * PI * x2);
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
//   Uses inner disk formula for temperature from Novikov & Thorne 1973.

void ThinDiskSource(MeshBlock *pmb, const Real time, const Real dt,
    const AthenaArray<Real> &prim_rad, AthenaArray<Real> &cons_rad) {

  // Extract information from block
  Coordinates *pcoord = pmb->pcoord;
  Radiation *prad = pmb->prad;
  int is = pmb->is;
  int ie = pmb->ie;
  int js = pmb->js;
  int je = pmb->je;
  int ks = pmb->ks;
  int ke = pmb->ke;
  int num_levels = pmb->loc.level - root_offset;

  // Calculate r-limits of disk region
  int il = 0;
  int iu = -1;
  if (pcoord->x1v(is) <= r_out and pcoord->x1v(ie) >= r_in) {
    for (il = is; pcoord->x1v(il) < r_in; ++il);
    for (iu = ie; pcoord->x1v(iu) > r_out; --iu);
  }

  // Calculate theta-limits of disk region
  int jl = 0;
  int ju = -1;
  int num_th = num_root_th;
  for (int n = 0; n < num_levels; ++n) {
    num_th *= 2;
  }
  if (num_th % 2 == 1 and pmb->loc.lx2 == num_th / 2) {
    jl = js + (je - js - 1) / 2;
    ju = js + (je - js) / 2;
  } else if (num_th % 2 == 0 and pmb->loc.lx2 == num_th / 2 - 1) {
    jl = ju = je;
  } else if (num_th % 2 == 0 and pmb->loc.lx2 == num_th / 2) {
    jl = ju = js;
  }

  // Overwrite primitive intensity in disk
  for (int i = il; i <= iu; ++i) {
    Real r = pcoord->x1v(i);
    Real r_1_4 = std::sqrt(std::sqrt(r));
    Real r_3_8 = r_1_4 * std::sqrt(r_1_4);
    Real aa = 1 + SQR(spin) / SQR(r) + 2.0 * mass * SQR(spin) / (r * SQR(r));
    Real bb = 1.0 + std::sqrt(mass) * spin / (r * std::sqrt(r));
    Real ee = 1.0 + 4.0 * SQR(spin) / SQR(r)
        - 4.0 * mass * SQR(spin) / (r * SQR(r))
        + 3.0 * SQR(SQR(spin)) / SQR(SQR(r));
    Real tt_cgs = 4.0e7 / std::sqrt(std::sqrt(alpha))
        / std::sqrt(std::sqrt(m_msun / 3.0)) / r_3_8 / std::sqrt(aa)
        * std::sqrt(bb) * std::sqrt(std::sqrt(ee));
    Real erad_cgs = prad->arad_cgs * SQR(SQR(tt_cgs));
    Real ii = erad_cgs / e_cgs / (4.0 * PI);
    for (int l = zs; l <= ze; ++l) {
      for (int m = ps; m <= pe; ++m) {
        int lm = prad->AngleInd(l, m);
        for (int k = ks; k <= ke; ++k) {
          for (int j = jl; j <= ju; ++j) {
            pmb->ruser_meshblock_data[0](lm,k,j,i) = ii;
          }
        }
      }
    }
  }

  // Overwrite conserved intensity in disk
  if (ju >= jl and iu >= il) {
    pmb->prad->PrimitiveToConserved(pmb->ruser_meshblock_data[0], cons_rad, pmb->pcoord,
        il, iu, jl, ju, ks, ke);
  }
  return;
}
