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
#include "../field.hpp"              // Field
#include "field_integrator.hpp"

// Athena headers
#include "../../athena.hpp"          // macros
#include "../../athena_arrays.hpp"   // AthenaArray
#include "../../mesh.hpp"            // MeshBlock
#include "../../parameter_input.hpp" 

//======================================================================================
//! \file field_integrator.cpp
//  \brief 
//======================================================================================

// constructor

FieldIntegrator::FieldIntegrator(Field *pfield, ParameterInput *pin)
{
  pmy_field = pfield;
  MeshBlock *pmb = pfield->pmy_mblock;

// Allocate memory for scratch vectors

  int nthreads = pmb->pmy_mesh->GetNumMeshThreads();
  int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmb->block_size.nx2 > 1) ncells2 = pmb->block_size.nx2 + 2*(NGHOST);
  if (pmb->block_size.nx3 > 1) ncells3 = pmb->block_size.nx3 + 2*(NGHOST);

  cc_e_.NewAthenaArray(ncells3,ncells2,ncells1);

  face_area_.NewAthenaArray(nthreads,ncells1);
  edge_length_.NewAthenaArray(nthreads,ncells1);
  edge_length_p1_.NewAthenaArray(nthreads,ncells1);
}

// destructor

FieldIntegrator::~FieldIntegrator()
{
  cc_e_.DeleteAthenaArray();
  face_area_.DeleteAthenaArray();
  edge_length_.DeleteAthenaArray();
  edge_length_p1_.DeleteAthenaArray();
}
