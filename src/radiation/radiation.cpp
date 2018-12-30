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

//----------------------------------------------------------------------------------------
// Radiation constructor
// Inputs:
//   pmb: pointer to containing MeshBlock
//   pin: pointer to runtime parameters

Radiation::Radiation(MeshBlock *pmb, ParameterInput *pin) {

  // Set MeshBlock pointer
  pmy_block = pmb;

  // Set numbers of angles
  nzeta = pin->GetInteger("radiation", "n_polar");
  npsi = pin->GetInteger("radiation", "n_azimuthal");

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
  zetaf(NGHOST) = 0.0;          // set north pole exactly
  zetaf(nzeta + NGHOST) = PI;   // set south pole exactly
  for (int l = NGHOST + 1; l <= (nzeta-1)/2 + NGHOST; ++l) {
    Real czeta = 1.0 + (l - NGHOST) * dczeta;
    Real zeta = std::acos(czeta);
    zetaf(l) = zeta;                          // set northern active faces
    zetaf(nzeta + 2*NGHOST - l) = PI - zeta;  // set southern active faces
  }
  if (nzeta%2 == 0) {
    zetaf(nzeta/2 + NGHOST) = PI/2.0;  // set equator exactly if present
  }
  for (int l = 0; l <= NGHOST - 1; ++l) {
    zetaf(l) = -zetaf(2*NGHOST - l);                          // set northern ghost faces
    zetaf(nzeta + 2*NGHOST - l) = 2.0*PI - zetaf(nzeta + l);  // set southern ghost faces
  }
  for (int l = 0; l <= nzeta + 2*NGHOST - 1; ++l) {
    zetav(l) = (zetaf(l+1) * std::cos(zetaf(l+1)) - std::sin(zetaf(l+1))
        - zetaf(l) * std::cos(zetaf(l)) + std::sin(zetaf(l))) / (std::cos(zetaf(l+1))
        - std::cos(zetaf(l)));
    dzetaf(l) = zetaf(l+1) - zetaf(l);
  }

  // Construct azimuthal angles, equally spaced
  Real dpsi = 2.0*PI / npsi;
  psif(NGHOST) = 0.0;            // set origin exactly
  psif(npsi + NGHOST) = 2.0*PI;  // set origin exactly
  for (int m = NGHOST + 1; m <= npsi + NGHOST - 1; ++m) {
    psif(m) = (m - NGHOST) * dpsi;  // set active faces
  }
  for (int m = 0; m <= NGHOST - 1; ++m) {
    psif(m) = psif(npsi + m) - 2.0*PI;                        // set beginning ghost faces
    psif(npsi + 2*NGHOST - m) = psif(2*NGHOST - m) + 2.0*PI;  // set end ghost faces
  }
  for (int m = 0; m <= npsi + 2*NGHOST - 1; ++m) {
    psiv(m) = 0.5 * (psif(m) + psif(m+1));
    dpsif(m) = psif(m+1) - psif(m);
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
}
