//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================

// Primary header
#include "bfield.hpp"

// C++ headers
#include <algorithm>  // min()
#include <cfloat>     // FLT_MAX
#include <cmath>      // fabs(), sqrt()

// Athena headers
#include "../athena.hpp"                  // array access, macros, Real
#include "../athena_arrays.hpp"           // AthenaArray
#include "bvals/bfield_bvals.hpp"         // FluidBoundaryConditions
#include "../mesh.hpp"                    // MeshBlock, Mesh
#include "../coordinates/coordinates.hpp" // VolumeCenterWidth()

//======================================================================================
//! \file bfield.cpp
//  \brief implementation of functions in class BField
//======================================================================================

// constructor, initializes data structures and parameters

BField::BField(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block_ = pmb;

// Allocate memory for interface fields, but only when needed.

  if (MAGNETIC_FIELDS_ENABLED) {
    int ncells1 = pmy_block_->block_size.nx1 + 2*(NGHOST);
    int ncells2 = 1, ncells3 = 1;
    if (pmy_block_->block_size.nx2 > 1) ncells2=pmy_block_->block_size.nx2 + 2*(NGHOST);
    if (pmy_block_->block_size.nx3 > 1) ncells3=pmy_block_->block_size.nx3 + 2*(NGHOST);

//  Note the extra cell in each longitudunal dirn
    b1i.NewAthenaArray(ncells3,ncells2,(ncells1+1));
    b2i.NewAthenaArray(ncells3,(ncells2+1),ncells1);
    b3i.NewAthenaArray((ncells3+1),ncells2,ncells1);

// Allocate memory for interface fields at intermediate-time step

    b1i_1.NewAthenaArray(ncells3,ncells2,(ncells1+1));
    b2i_1.NewAthenaArray(ncells3,(ncells2+1),ncells1);
    b3i_1.NewAthenaArray((ncells3+1),ncells2,ncells1);

// Construct ptrs to objects of various classes needed to integrate B-field

//  pf_integrator = new FluidIntegrator(this);
//  pf_bcs = new FluidBoundaryConditions(this,pin);

  }
}

// destructor

BField::~BField()
{
  b1i.DeleteAthenaArray();
  b2i.DeleteAthenaArray();
  b3i.DeleteAthenaArray();
  b1i_1.DeleteAthenaArray();
  b2i_1.DeleteAthenaArray();
  b3i_1.DeleteAthenaArray();
}
