//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gr_rad_beam.cpp
//  \brief Problem generator for GR radiation in box with some hohlraum boundaries

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

// Global variables
namespace {
Real e_rad;           // initial radiation energy density
Real ii_ix1, ii_ox1;  // x1-boundary radiation intensities
Real ii_ix2, ii_ox2;  // x2-boundary radiation intensities
Real ii_ix3, ii_ox3;  // x3-boundary radiation intensities
}  // namespace

// Declarations
void HohlraumIX1(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, AthenaArray<Real> &prim_rad, Real time, Real dt, int il, int iu,
    int jl, int ju, int kl, int ku, int ngh);
void HohlraumOX1(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, AthenaArray<Real> &prim_rad, Real time, Real dt, int il, int iu,
    int jl, int ju, int kl, int ku, int ngh);
void HohlraumIX2(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, AthenaArray<Real> &prim_rad, Real time, Real dt, int il, int iu,
    int jl, int ju, int kl, int ku, int ngh);
void HohlraumOX2(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, AthenaArray<Real> &prim_rad, Real time, Real dt, int il, int iu,
    int jl, int ju, int kl, int ku, int ngh);
void HohlraumIX3(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, AthenaArray<Real> &prim_rad, Real time, Real dt, int il, int iu,
    int jl, int ju, int kl, int ku, int ngh);
void HohlraumOX3(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, AthenaArray<Real> &prim_rad, Real time, Real dt, int il, int iu,
    int jl, int ju, int kl, int ku, int ngh);

//----------------------------------------------------------------------------------------
// Function for preparing Mesh
// Inputs:
//   pin: input parameters
// Outputs: (none)

void Mesh::InitUserMeshData(ParameterInput *pin) {

  // Read parameters from input file
  e_rad = pin->GetReal("problem", "e_rad");
  ii_ix1 = pin->GetReal("problem", "e_rad_ix1") / (4.0*PI);
  ii_ox1 = pin->GetReal("problem", "e_rad_ox1") / (4.0*PI);
  ii_ix2 = pin->GetReal("problem", "e_rad_ix2") / (4.0*PI);
  ii_ox2 = pin->GetReal("problem", "e_rad_ox2") / (4.0*PI);
  ii_ix3 = pin->GetReal("problem", "e_rad_ix3") / (4.0*PI);
  ii_ox3 = pin->GetReal("problem", "e_rad_ox3") / (4.0*PI);

  // Enroll boundary functions
  if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, HohlraumIX1);
  }
  if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, HohlraumOX1);
  }
  if (mesh_bcs[BoundaryFace::inner_x2] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x2, HohlraumIX2);
  }
  if (mesh_bcs[BoundaryFace::outer_x2] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x2, HohlraumOX2);
  }
  if (mesh_bcs[BoundaryFace::inner_x3] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x3, HohlraumIX3);
  }
  if (mesh_bcs[BoundaryFace::outer_x3] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x3, HohlraumOX3);
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for setting initial conditions
// Inputs:
//   pin: parameters (unused)
// Outputs: (none)
// Notes:
//   initializes constant state with no radiation

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  prad->CalculateConstantRadiation(e_rad, 0.0, 0.0, 0.0, prad->cons);
  return;
}

//----------------------------------------------------------------------------------------
// Inner x1 hohlraum boundary condition
// Inputs:
//   pmb: pointer to MeshBlock
//   pcoord: pointer to Coordinates (not used)
//   time,dt: current time and timestep of simulation (not used)
//   is,ie,js,je,ks,ke: indices demarkating active region
// Outputs:
//   prim: primitives set in ghost zones
//   bb: face-centered magnetic field set in ghost zones
//   prim_rad: radiation primitives set in ghost zones
// Notes:
//   Sets intensities in ghost zones to be fixed value in all angles.
//   Ignores hydro and MHD variables.

void HohlraumIX1(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, AthenaArray<Real> &prim_rad, Real time, Real dt, int il, int iu,
    int jl, int ju, int kl, int ku, int ngh) {
  for (int n = 0; n < pmb->prad->nang; ++n) {
    for (int k = kl; k <= ku; ++k) {
      for (int j = jl; j <= ju; ++j) {
        for (int i = il-ngh; i <= il-1; ++i) {
          prim_rad(n,k,j,i) = ii_ix1;
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Outer x1 hohlraum boundary condition
// Inputs:
//   pmb: pointer to MeshBlock
//   pcoord: pointer to Coordinates (not used)
//   time,dt: current time and timestep of simulation (not used)
//   is,ie,js,je,ks,ke: indices demarkating active region
// Outputs:
//   prim: primitives set in ghost zones
//   bb: face-centered magnetic field set in ghost zones
//   prim_rad: radiation primitives set in ghost zones
// Notes:
//   Sets intensities in ghost zones to be fixed value in all angles.
//   Ignores hydro and MHD variables.

void HohlraumOX1(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, AthenaArray<Real> &prim_rad, Real time, Real dt, int il, int iu,
    int jl, int ju, int kl, int ku, int ngh) {
  for (int n = 0; n < pmb->prad->nang; ++n) {
    for (int k = kl; k <= ku; ++k) {
      for (int j = jl; j <= ju; ++j) {
        for (int i = iu+1; i <= iu+ngh; ++i) {
          prim_rad(n,k,j,i) = ii_ox1;
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Inner x2 hohlraum boundary condition
// Inputs:
//   pmb: pointer to MeshBlock
//   pcoord: pointer to Coordinates (not used)
//   time,dt: current time and timestep of simulation (not used)
//   is,ie,js,je,ks,ke: indices demarkating active region
// Outputs:
//   prim: primitives set in ghost zones
//   bb: face-centered magnetic field set in ghost zones
//   prim_rad: radiation primitives set in ghost zones
// Notes:
//   Sets intensities in ghost zones to be fixed value in all angles.
//   Ignores hydro and MHD variables.

void HohlraumIX2(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, AthenaArray<Real> &prim_rad, Real time, Real dt, int il, int iu,
    int jl, int ju, int kl, int ku, int ngh) {
  for (int n = 0; n < pmb->prad->nang; ++n) {
    for (int k = kl; k <= ku; ++k) {
      for (int j = jl-ngh; j <= jl-1; ++j) {
        for (int i = il; i <= iu; ++i) {
          prim_rad(n,k,j,i) = ii_ix2;
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Outer x2 hohlraum boundary condition
// Inputs:
//   pmb: pointer to MeshBlock
//   pcoord: pointer to Coordinates (not used)
//   time,dt: current time and timestep of simulation (not used)
//   is,ie,js,je,ks,ke: indices demarkating active region
// Outputs:
//   prim: primitives set in ghost zones
//   bb: face-centered magnetic field set in ghost zones
//   prim_rad: radiation primitives set in ghost zones
// Notes:
//   Sets intensities in ghost zones to be fixed value in all angles.
//   Ignores hydro and MHD variables.

void HohlraumOX2(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, AthenaArray<Real> &prim_rad, Real time, Real dt, int il, int iu,
    int jl, int ju, int kl, int ku, int ngh) {
  for (int n = 0; n < pmb->prad->nang; ++n) {
    for (int k = kl; k <= ku; ++k) {
      for (int j = ju+1; j <= ju+ngh; ++j) {
        for (int i = il; i <= iu; ++i) {
          prim_rad(n,k,j,i) = ii_ox2;
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Inner x3 hohlraum boundary condition
// Inputs:
//   pmb: pointer to MeshBlock
//   pcoord: pointer to Coordinates (not used)
//   time,dt: current time and timestep of simulation (not used)
//   is,ie,js,je,ks,ke: indices demarkating active region
// Outputs:
//   prim: primitives set in ghost zones
//   bb: face-centered magnetic field set in ghost zones
//   prim_rad: radiation primitives set in ghost zones
// Notes:
//   Sets intensities in ghost zones to be fixed value in all angles.
//   Ignores hydro and MHD variables.

void HohlraumIX3(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, AthenaArray<Real> &prim_rad, Real time, Real dt, int il, int iu,
    int jl, int ju, int kl, int ku, int ngh) {
  for (int n = 0; n < pmb->prad->nang; ++n) {
    for (int k = kl-ngh; k <= kl-1; ++k) {
      for (int j = jl; j <= ju; ++j) {
        for (int i = il; i <= iu; ++i) {
          prim_rad(n,k,j,i) = ii_ix3;
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Outer x3 hohlraum boundary condition
// Inputs:
//   pmb: pointer to MeshBlock
//   pcoord: pointer to Coordinates (not used)
//   time,dt: current time and timestep of simulation (not used)
//   is,ie,js,je,ks,ke: indices demarkating active region
// Outputs:
//   prim: primitives set in ghost zones
//   bb: face-centered magnetic field set in ghost zones
//   prim_rad: radiation primitives set in ghost zones
// Notes:
//   Sets intensities in ghost zones to be fixed value in all angles.
//   Ignores hydro and MHD variables.

void HohlraumOX3(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, AthenaArray<Real> &prim_rad, Real time, Real dt, int il, int iu,
    int jl, int ju, int kl, int ku, int ngh) {
  for (int n = 0; n < pmb->prad->nang; ++n) {
    for (int k = ku+1; k <= ku+ngh; ++k) {
      for (int j = jl; j <= ju; ++j) {
        for (int i = il; i <= iu; ++i) {
          prim_rad(n,k,j,i) = ii_ox3;
        }
      }
    }
  }
  return;
}
