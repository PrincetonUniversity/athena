//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 *
 * This program is free software: you can redistribute and/or modify it under the terms
 * of the GNU General Public License (GPL) as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of GNU GPL in the file LICENSE included in the code
 * distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

// Primary header
#include "integrators.hpp"

// Athena headers
#include "../../athena.hpp"          // macros
#include "../../athena_arrays.hpp"   // AthenaArray
#include "../fluid.hpp"              // Fluid
#include "../../mesh.hpp"            // MeshBlock
#include "../../parameter_input.hpp" 

//======================================================================================
/*! \file integrators.cpp
 *  \brief 
 *====================================================================================*/

// constructor

FluidIntegrator::FluidIntegrator(Fluid *pf, ParameterInput *pin)
{
  pmy_fluid = pf;

// Allocate memory for scratch vectors

  int ncells1 = pf->pmy_block->block_size.nx1 + 2*(NGHOST);
  Mesh *pmm = pf->pmy_block->pmy_domain->pmy_mesh;

  wl_.NewAthenaArray(pmm->nthreads_mesh,NFLUID,ncells1);
  wr_.NewAthenaArray(pmm->nthreads_mesh,NFLUID,ncells1);
  flx_.NewAthenaArray(pmm->nthreads_mesh,NFLUID,ncells1);
  src_.NewAthenaArray(pmm->nthreads_mesh,NFLUID,ncells1);
  face_area_.NewAthenaArray(pmm->nthreads_mesh,ncells1);
  cell_volume_.NewAthenaArray(pmm->nthreads_mesh,ncells1);
}

// destructor

FluidIntegrator::~FluidIntegrator()
{
  wl_.DeleteAthenaArray();
  wr_.DeleteAthenaArray();
  flx_.DeleteAthenaArray();
  src_.DeleteAthenaArray();
  face_area_.DeleteAthenaArray();
  cell_volume_.DeleteAthenaArray();
}
