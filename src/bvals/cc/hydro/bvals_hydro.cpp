//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hydro_bc.cpp
//  \brief implements boundary functions for Hydro variables in a derived class of
//  the CellCenteredBoundaryFunctions abstract base class.

// Athena++ headers
#include "bvals.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"

//----------------------------------------------------------------------------------------
//! \class HydroBoundaryFunctions

HydroBoundaryFunctions::HydroBoundaryFunctions()
: CellCenteredBoundaryFunctions() {
}

// destructor
HydroBoundaryFunctions::~HydroBoundaryFunctions() {

}
