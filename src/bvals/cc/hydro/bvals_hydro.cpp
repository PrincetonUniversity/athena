//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_hydro.cpp
//  \brief implements boundary functions for Hydro variables in a derived class of
//  the CellCenteredBoundaryVariable base class.

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../mesh/mesh.hpp"
#include "../bvals_cc.hpp"

//----------------------------------------------------------------------------------------
//! \class HydroBoundaryFunctions

HydroBoundaryVariable::HydroBoundaryFunctions()
: CellCenteredBoundaryFunctions() {
}

// destructor
HydroBoundaryVariable::~HydroBoundaryFunctions() {

}
