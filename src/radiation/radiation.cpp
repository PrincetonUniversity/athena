//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file radiation.cpp
//  \brief implementation of functions in class Radiation

// C++ headers
#include <cmath>      // acos(), cos(), sin()
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str(), string

// Athena++ headers
#include "radiation.hpp"
#include "../athena_arrays.hpp"    // AthenaArray
#include "../parameter_input.hpp"  // ParameterInput
#include "../mesh/mesh.hpp"        // MeshBlock

//----------------------------------------------------------------------------------------
// Radiation constructor
// Inputs:
//   pmb: pointer to containing MeshBlock
//   pin: pointer to runtime parameters

Radiation::Radiation(MeshBlock *pmb, ParameterInput *pin) {

  // Set object pointers
  pmy_block = pmb;

  // Set parameters
  nzeta = pin->GetInteger("radiation", "n_polar");
  npsi = pin->GetInteger("radiation", "n_azimuthal");
  nang = (nzeta + 2*NGHOST) * (npsi + 2*NGHOST);
  zs = NGHOST;
  ze = nzeta + NGHOST - 1;
  ps = NGHOST;
  pe = npsi + NGHOST - 1;
  is = pmy_block->is;
  ie = pmy_block->ie;
  js = pmy_block->js;
  je = pmy_block->je;
  ks = pmy_block->ks;
  ke = pmy_block->ke;

  // Verify numbers of angles
  std::stringstream msg;
  if (nzeta < 4) {
    msg << "### FATAL ERROR in Radiation constructor\n";
    msg << "too few polar angles\n";
    throw std::runtime_error(msg.str().c_str());
  }
  if (npsi < 4) {
    msg << "### FATAL ERROR in Radiation constructor\n";
    msg << "too few azimuthal angles\n";
    throw std::runtime_error(msg.str().c_str());
  }

  // Allocate memory for angles
  zetaf.NewAthenaArray(nzeta + 2*NGHOST + 1);
  zetav.NewAthenaArray(nzeta + 2*NGHOST);
  dzetaf.NewAthenaArray(nzeta + 2*NGHOST);
  psif.NewAthenaArray(npsi + 2*NGHOST + 1);
  psiv.NewAthenaArray(npsi + 2*NGHOST);
  dpsif.NewAthenaArray(npsi + 2*NGHOST);

  // Construct polar angles, equally spaced in cosine
  Real dczeta = -2.0 / nzeta;
  zetaf(zs) = 0.0;             // set north pole exactly
  zetaf(ze+1) = PI;            // set south pole exactly
  for (int l = zs+1; l <= (nzeta-1)/2+NGHOST; ++l) {
    Real czeta = 1.0 + (l - NGHOST) * dczeta;
    Real zeta = std::acos(czeta);
    zetaf(l) = zeta;                           // set northern active faces
    zetaf(ze+NGHOST+1-l) = PI - zeta;          // set southern active faces
  }
  if (nzeta%2 == 0) {
    zetaf(nzeta/2+NGHOST) = PI/2.0;  // set equator exactly if present
  }
  for (int l = zs-NGHOST; l <= zs-1; ++l) {
    zetaf(l) = -zetaf(2*NGHOST - l);                 // set northern ghost faces
    zetaf(ze+NGHOST+1-l) = 2.0*PI - zetaf(nzeta+l);  // set southern ghost faces
  }
  for (int l = zs-NGHOST; l <= ze+NGHOST; ++l) {
    zetav(l) = (zetaf(l+1) * std::cos(zetaf(l+1)) - std::sin(zetaf(l+1))
        - zetaf(l) * std::cos(zetaf(l)) + std::sin(zetaf(l))) / (std::cos(zetaf(l+1))
        - std::cos(zetaf(l)));
    dzetaf(l) = zetaf(l+1) - zetaf(l);
  }

  // Construct azimuthal angles, equally spaced
  Real dpsi = 2.0*PI / npsi;
  psif(ps) = 0.0;             // set origin exactly
  psif(pe+1) = 2.0*PI;        // set origin exactly
  for (int m = ps+1; m <= pe; ++m) {
    psif(m) = (m - NGHOST) * dpsi;  // set active faces
  }
  for (int m = ps-NGHOST; m <= ps-1; ++m) {
    psif(m) = psif(npsi+m) - 2.0*PI;                  // set beginning ghost faces
    psif(pe+NGHOST+1-m) = psif(2*NGHOST-m) + 2.0*PI;  // set end ghost faces
  }
  for (int m = ps-NGHOST; m <= pe+NGHOST; ++m) {
    psiv(m) = 0.5 * (psif(m) + psif(m+1));
    dpsif(m) = psif(m+1) - psif(m);
  }

  // Allocate memory for intensities
  int num_cells_1 = ie + NGHOST;
  int num_cells_2 = 1;
  if (js != je) {
    num_cells_2 = je + NGHOST;
  }
  int num_cells_3 = 1;
  if (ks != ke) {
    num_cells_3 = ke + NGHOST;
  }
  prim.NewAthenaArray(nang, num_cells_3, num_cells_2, num_cells_1);
  prim1.NewAthenaArray(nang, num_cells_3, num_cells_2, num_cells_1);
  cons.NewAthenaArray(nang, num_cells_3, num_cells_2, num_cells_1);
  cons1.NewAthenaArray(nang, num_cells_3, num_cells_2, num_cells_1);

  // Allocate memory for fluxes
  flux[X1DIR].NewAthenaArray(nang, num_cells_3, num_cells_2, num_cells_1 + 1);
  if (js != je) {
    flux[X2DIR].NewAthenaArray(nang, num_cells_3, num_cells_2 + 1, num_cells_1);
  }
  if (ks != ke) {
    flux[X3DIR].NewAthenaArray(nang, num_cells_3 + 1, num_cells_2, num_cells_1);
  }
}

//----------------------------------------------------------------------------------------
// Radiation destructor

Radiation::~Radiation() {
  zetaf.DeleteAthenaArray();
  zetav.DeleteAthenaArray();
  dzetaf.DeleteAthenaArray();
  psif.DeleteAthenaArray();
  psiv.DeleteAthenaArray();
  dpsif.DeleteAthenaArray();
  flux[X1DIR].DeleteAthenaArray();
  flux[X2DIR].DeleteAthenaArray();
  flux[X3DIR].DeleteAthenaArray();
}

//----------------------------------------------------------------------------------------
// Indexing function for intensities
// Inputs:
//   l: zeta-index
//   m: psi-index
// Outputs:
//   returned value: 1D index for both zeta and psi

int Radiation::IntInd(int l, int m) {
  return l * (pe + NGHOST + 1) + m;
}

//----------------------------------------------------------------------------------------
// Radiation conversion from primitive to conserved variables
// Inputs:
//   prim_vals: primitives
//   pcoord: pointer to Coordinates
//   il,iu,jl,ju,kl,ku: index bounds of region to be updated
// Outputs:
//   cons_vals: conserved quantities

void Radiation::PrimitiveToConserved(AthenaArray<Real> &prim_vals,
    AthenaArray<Real> &cons_vals, Coordinates *pcoord, int il, int iu, int jl, int ju,
    int kl, int ku) {
  return;
}

//----------------------------------------------------------------------------------------
// Radiation inversion from conserved to primitive variables
// Inputs:
//   cons_vals: conserved quantities
//   pcoord: pointer to Coordinates
//   il,iu,jl,ju,kl,ku: index bounds of region to be updated
// Outputs:
//   prim_vals: primitives

void Radiation::ConservedToPrimitive(AthenaArray<Real> &cons_vals,
    AthenaArray<Real> &prim_vals, Coordinates *pcoord, int il, int iu, int jl, int ju,
    int kl, int ku) {
  return;
}
