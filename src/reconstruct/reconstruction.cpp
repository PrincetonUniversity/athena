//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file reconstruction.cpp
//  \brief 

// Athena++ headers
#include "reconstruction.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp" 
#include "../mesh/mesh.hpp"

// constructor

Reconstruction::Reconstruction(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block_ = pmb;

  // set function pointers for reconstruction functions in each direction
  if (pmb->block_size.x1rat == 1.0) {
    ReconstructFuncX1 = PiecewiseLinearUniformX1;
  } else {
    ReconstructFuncX1 = PiecewiseLinearX1;
  }

  if (pmb->block_size.x2rat == 1.0) {
    ReconstructFuncX2 = PiecewiseLinearUniformX2;
  } else {
    ReconstructFuncX2 = PiecewiseLinearX2;
  }

  if (pmb->block_size.x3rat == 1.0) {
    ReconstructFuncX3 = PiecewiseLinearUniformX3;
  } else {
    ReconstructFuncX3 = PiecewiseLinearX3;
  }
}

// destructor

Reconstruction::~Reconstruction()
{
}
