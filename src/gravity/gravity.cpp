//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gravity.cpp
//  \brief implementation of functions in class Gravity

// Athena++ headers
#include "gravity.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../parameter_input.hpp"

// constructor, initializes data structures and parameters

Gravity::Gravity(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block = pmb;
  four_pi_G=pin->GetOrAddReal("problem", "four_pi_G", 1.0); // default: 4piG=1

  // Allocate memory for gravitational potential, but only when needed.
  int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmb->block_size.nx2 > 1) ncells2 = pmb->block_size.nx2 + 2*(NGHOST);
  if (pmb->block_size.nx3 > 1) ncells3 = pmb->block_size.nx3 + 2*(NGHOST);

  phi.NewAthenaArray(ncells3,ncells2,ncells1);
  phi_old.NewAthenaArray(ncells3,ncells2,ncells1);

  Initialize(pin);
}

// destructor

Gravity::~Gravity()
{
  phi.DeleteAthenaArray();
  phi_old.DeleteAthenaArray();
}

