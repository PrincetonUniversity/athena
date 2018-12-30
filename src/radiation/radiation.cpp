//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file radiation.cpp
//  \brief implementation of functions in class Radiation

// C++ headers
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
  if (nzeta < 2) {
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
  zetaf.NewAthenaArray(nzeta+1);
  zetav.NewAthenaArray(nzeta);
  psif.NewAthenaArray(npsi+1);
  psiv.NewAthenaArray(npsi);
}

//----------------------------------------------------------------------------------------
// Radiation destructor

Hydro::~Hydro() {
  zetaf.DeleteAthenaArray();
  zetav.DeleteAthenaArray();
  psif.DeleteAthenaArray();
  psiv.DeleteAthenaArray();
}
