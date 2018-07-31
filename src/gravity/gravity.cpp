//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gravity.cpp
//  \brief implementation of functions in class Gravity

// C/C++ headers
#include <iostream>
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "gravity.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../parameter_input.hpp"
#include "../bvals/bvals_grav.hpp"

// constructor, initializes data structures and parameters

Gravity::Gravity(MeshBlock *pmb, ParameterInput *pin) {
  pmy_block = pmb;
  four_pi_G=pmb->pmy_mesh->four_pi_G_; // default: 4piG=1
  if (four_pi_G==0.0) {
   std::stringstream msg;
   msg << "### FATAL ERROR in Gravity::Gravity" << std::endl
        << "Gravitational constant must be set in the Mesh::InitUserMeshData "
        << "using the SetGravitationalConstant or SetFourPiG function." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }
  grav_mean_rho=pmb->pmy_mesh->grav_mean_rho_;
  if (grav_mean_rho==-1.0) {
   std::stringstream msg;
   msg << "### FATAL ERROR in Gravity::Gravity" << std::endl
        << "Background Mean Density must be set in the Mesh::InitUserMeshData "
        << "using the SetMeanDensity function." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }

  // Allocate memory for gravitational potential, but only when needed.
  int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmb->block_size.nx2 > 1) ncells2 = pmb->block_size.nx2 + 2*(NGHOST);
  if (pmb->block_size.nx3 > 1) ncells3 = pmb->block_size.nx3 + 2*(NGHOST);

  phi.NewAthenaArray(ncells3,ncells2,ncells1);
}

// destructor

Gravity::~Gravity() {
  phi.DeleteAthenaArray();
}
