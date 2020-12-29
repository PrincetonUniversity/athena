//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gr_geodesic_infall.cpp
//! \brief Problem generator for dust falling onto black hole.

// C headers

// C++ headers
#include <cassert>  // assert
#include <cmath>    // pow(), sin(), sqrt()
#include <cstring>  // strcmp()

// Athena++ headers
#include "../athena.hpp"                   // enums, Real, FaceField
#include "../athena_arrays.hpp"            // AthenaArray
#include "../bvals/bvals.hpp"              // BoundaryValues
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../eos/eos.hpp"                  // EquationOfState
#include "../field/field.hpp"              // Field
#include "../hydro/hydro.hpp"              // Hydro
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"          // ParameterInput

// Configuration checking
#if not GENERAL_RELATIVITY
#error "This problem generator must be used with general relativity"
#endif
#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator does not support magnetic fields"
#endif

// Declarations
void FixedBoundary(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
                   FaceField &bb, Real time, Real dt,
                   int il, int iu, int jl, int ju, int kl, int ku, int ngh);
namespace {
// TODO(felker): can the 4x copies of this function in pgen/ files be shared?
void GetBoyerLindquistCoordinates(Real x1, Real x2, Real x3, Real *pr,
                                  Real *ptheta, Real *pphi);
} // namespace

//----------------------------------------------------------------------------------------
// Function for initializing global mesh properties
// Inputs:
//   pin: input parameters (unused)
// Outputs: (none)

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Enroll boundary functions
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, FixedBoundary);
  return;
}

//----------------------------------------------------------------------------------------
// Function for setting initial conditions
// Inputs:
//   phyd: Hydro
//   pfld: Field (unused)
//   pin: parameters
// Outputs: (none)
// Notes:
//   assumes x3 is axisymmetric direction

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Prepare index bounds
  int il = is - NGHOST;
  int iu = ie + NGHOST;
  int jl = js;
  int ju = je;
  if (block_size.nx2 > 1) {
    jl -= NGHOST;
    ju += NGHOST;
  }
  int kl = ks;
  int ku = ke;
  if (block_size.nx3 > 1) {
    kl -= NGHOST;
    ku += NGHOST;
  }

  // Read problem properties
  Real rho_min = pin->GetReal("hydro", "rho_min");
  Real rho_pow = pin->GetReal("hydro", "rho_pow");
  Real pgas_min = pin->GetReal("hydro", "pgas_min");
  Real pgas_pow = pin->GetReal("hydro", "pgas_pow");

  // Initialize primitive values
  AthenaArray<Real> g, gi;
  g.NewAthenaArray(NMETRIC, iu+1);
  gi.NewAthenaArray(NMETRIC, iu+1);
  for (int j=jl; j<=ju; j++) {
    for (int i=il; i<=iu; i++) {
      // Get Boyer-Lindquist coordinates of cell
      Real r=0;
      Real theta=0;
      Real phi=0;
      GetBoyerLindquistCoordinates(pcoord->x1v(i), pcoord->x2v(j), pcoord->x3v(kl), &r,
                                   &theta, &phi);

      // Calculate primitives depending on location
      Real rho = rho_min * std::pow(r, rho_pow);
      Real pgas = pgas_min * std::pow(r, pgas_pow);
      Real uu1 = 0.0;
      Real uu2 = 0.0;
      Real uu3 = 0.0;

      // Set primitive values
      for (int k=kl; k<=ku; k++) {
        phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) = rho;
        phydro->w(IPR,k,j,i) = phydro->w1(IPR,k,j,i) = pgas;
        phydro->w(IVX,k,j,i) = phydro->w1(IM1,k,j,i) = uu1;
        phydro->w(IVY,k,j,i) = phydro->w1(IM2,k,j,i) = uu2;
        phydro->w(IVZ,k,j,i) = phydro->w1(IM3,k,j,i) = uu3;
      }
    }
  }

  // Initialize conserved values
  AthenaArray<Real> bb;
  bb.NewAthenaArray(3, ku+1, ju+1, iu+1);
  peos->PrimitiveToConserved(phydro->w, bb, phydro->u, pcoord, il, iu, jl, ju, kl, ku);
  return;
}

//----------------------------------------------------------------------------------------
// Fixed boundary condition
// Inputs:
//   pmb: pointer to MeshBlock
//   pcoord: pointer to Coordinates
//   time,dt: current time and timestep of simulation
//   is,ie,js,je,ks,ke: indices demarkating active region
// Outputs:
//   prim: primitives set in ghost zones
//   bb: face-centered magnetic field set in ghost zones
// Notes:
//   does nothing

void FixedBoundary(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
                   FaceField &bb, Real time, Real dt,
                   int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  return;
}

namespace {
//----------------------------------------------------------------------------------------
// Function for returning corresponding Boyer-Lindquist coordinates of point
// Inputs:
//   x1,x2,x3: global coordinates to be converted
// Outputs:
//   pr,ptheta,pphi: variables pointed to set to Boyer-Lindquist coordinates
// Notes:
//   conversion is trivial in all currently implemented coordinate systems

void GetBoyerLindquistCoordinates(Real x1, Real x2, Real x3, Real *pr,
                                  Real *ptheta, Real *pphi) {
  if (std::strcmp(COORDINATE_SYSTEM, "schwarzschild") == 0 ||
      std::strcmp(COORDINATE_SYSTEM, "kerr-schild") == 0) {
    *pr = x1;
    *ptheta = x2;
    *pphi = x3;
  }
  return;
}
} // namespace
