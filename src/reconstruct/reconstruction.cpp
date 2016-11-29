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

// constructor

Reconstruction::Reconstruction(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block_ = pmb;
}

// destructor

Reconstruction::~Reconstruction()
{
}
