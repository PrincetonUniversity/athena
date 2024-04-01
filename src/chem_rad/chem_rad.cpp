//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file chem_rad.cpp
//! \brief implementation of functions in class ChemRadiation

// C headers

// C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "chem_rad.hpp"
#include "integrators/rad_integrators.hpp"

//! ChemRadiation constructor, initializes data structures and parameters

ChemRadiation::ChemRadiation(MeshBlock *pmb, ParameterInput *pin) {
  // read in the parameters
  integrator = CHEMRADIATION_INTEGRATOR;
  pmy_block = pmb;
  //number of frequency bands
  nfreq = pin->GetOrAddInteger("chem_radiation","n_frequency",1);
  output_zone_sec = pin->GetOrAddBoolean("chem_radiation", "output_zone_sec", false);

  if (integrator == "six_ray") {
    nang = 6;
  } else if (integrator == "const") {
    nang = 1;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in ChemRadiation constructor" << std::endl
        << "integrator=" << integrator << " not valid chemradiation integrator, "
        << std::endl
        << "choose from {six_ray, const}" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  n_fre_ang = nang * nfreq;


  // allocate arrays
  int ncells1 = pmy_block->ncells1;
  int ncells2 = pmy_block->ncells2;
  int ncells3 = pmy_block->ncells3;
  // store frequency and angles as [nfre][ang]
  ir.NewAthenaArray(ncells3, ncells2, ncells1, n_fre_ang);
  ir_avg.NewAthenaArray(nfreq, ncells3, ncells2, ncells1);

  //radiation integrator
  pchemradintegrator = new ChemRadIntegrator(this, pin);
}

ChemRadiation::~ChemRadiation() {
  delete pchemradintegrator;
}
