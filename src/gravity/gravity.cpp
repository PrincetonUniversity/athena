//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gravity.cpp
//  \brief implementation of functions in class Field

// Athena++ headers
#include "gravity.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"

// constructor, initializes data structures and parameters

Gravity::Gravity(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block = pmb;

  // Allocate memory for gravitational potential, but only when needed.
  if (SELF_GRAVITY_ENABLED) {
    int ncells1 = pmb->block_size.nx1;
    int ncells2 = 1, ncells3 = 1;
    if (pmb->block_size.nx2 > 1) ncells2 = pmb->block_size.nx2;
    if (pmb->block_size.nx3 > 1) ncells3 = pmb->block_size.nx3;

    phi.NewAthenaArray (ncells3,ncells2,ncells1);

    // Allocate memory for scratch vectors
    den_.NewAthenaArray (ncells3,ncells2,ncells1);
  }
}

// destructor

Gravity::~Gravity()
{
  phi.DeleteAthenaArray();
  den_.DeleteAthenaArray();
}

//----------------------------------------------------------------------------------------
// \! fn
// \! brief

