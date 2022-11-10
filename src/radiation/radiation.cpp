//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file radiation.cpp
//! \brief implementation of functions in class Radiation

//C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "integrators/rad_integrators.hpp"
#include "radiation.hpp"

//! Radiation constructor, initializes data structures and parameters

Radiation::Radiation(MeshBlock *pmb, ParameterInput *pin) {
  // read in the parameters
  integrator = RADIATION_INTEGRATOR;
  pmy_block = pmb;
  //number of frequency bands
  nfreq = pin->GetOrAddInteger("radiation","n_frequency",1);
  output_zone_sec = pin->GetOrAddBoolean("radiation", "output_zone_sec", false);

  if (integrator == "six_ray") {
    nang = 6;
  } else if (integrator == "const") {
    nang = 1;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in Radiation constructor" << std::endl
        << "integrator=" << integrator << " not valid radiation integrator, " << std::endl
        << "choose from {six_ray, const}" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  n_fre_ang = nang * nfreq;


  // allocate arrays
  int ncells1 = pmy_block->ncells1;
  int ncells2 = 1, ncells3 = 1;
  if (pmy_block->block_size.nx2 > 1) ncells2 = pmy_block->block_size.nx2 + 2*(NGHOST);
  if (pmy_block->block_size.nx3 > 1) ncells3 = pmy_block->block_size.nx3 + 2*(NGHOST);
  // store frequency and angles as [nfre][ang]
  ir.NewAthenaArray(ncells3, ncells2, ncells1, n_fre_ang);
  ir_avg.NewAthenaArray(nfreq, ncells3, ncells2, ncells1);

  //radiation integrator
  pradintegrator = new RadIntegrator(this, pin);
}

Radiation::~Radiation() {
  delete pradintegrator;
}
