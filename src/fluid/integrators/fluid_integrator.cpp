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
#include "fluid_integrator.hpp"

// Athena headers
#include "../../athena.hpp"          // macros
#include "../../athena_arrays.hpp"   // AthenaArray
#include "../fluid.hpp"              // Fluid
#include "../../mesh.hpp"            // MeshBlock
#include "../../parameter_input.hpp" 

//======================================================================================
//! \file integrators.cpp
//  \brief 
//======================================================================================

// constructor

FluidIntegrator::FluidIntegrator(Fluid *pf, ParameterInput *pin)
{
  pmy_fluid = pf;

// Allocate memory for scratch vectors

  int nthreads = pf->pmy_block->pmy_mesh->GetNumMeshThreads();
  int ncells1 = pf->pmy_block->block_size.nx1 + 2*(NGHOST);
  int ncells2 = pf->pmy_block->block_size.nx2 + 2*(NGHOST);

  wl_.NewAthenaArray(nthreads,(NWAVE),ncells1);
  wr_.NewAthenaArray(nthreads,(NWAVE),ncells1);
  flx_.NewAthenaArray(nthreads,(NWAVE),ncells1);
  jflx_j_.NewAthenaArray(nthreads,(NWAVE),ncells1);
  kflx_k_.NewAthenaArray(nthreads,(NWAVE),ncells2,ncells1);
  face_area_.NewAthenaArray(nthreads,ncells1);
  face_area_p1_.NewAthenaArray(nthreads,ncells1);
  cell_volume_.NewAthenaArray(nthreads,ncells1);
  if (MAGNETIC_FIELDS_ENABLED && RELATIVISTIC_DYNAMICS)  // only used in SR/GRMHD
    b_normal_.NewAthenaArray(ncells1);
  if (GENERAL_RELATIVITY)  // only used in GR (and only in certain Riemann solvers)
  {
    g_.NewAthenaArray(NMETRIC,ncells1);
    g_inv_.NewAthenaArray(NMETRIC,ncells1);
    cons_.NewAthenaArray(NWAVE,ncells1);
  }
}

// destructor

FluidIntegrator::~FluidIntegrator()
{
  wl_.DeleteAthenaArray();
  wr_.DeleteAthenaArray();
  flx_.DeleteAthenaArray();
  jflx_j_.DeleteAthenaArray();
  kflx_k_.DeleteAthenaArray();
  face_area_.DeleteAthenaArray();
  face_area_p1_.DeleteAthenaArray();
  cell_volume_.DeleteAthenaArray();
  if (MAGNETIC_FIELDS_ENABLED && RELATIVISTIC_DYNAMICS)  // only used in SR/GRMHD
    b_normal_.DeleteAthenaArray();
  if (GENERAL_RELATIVITY)
  {
    g_.DeleteAthenaArray();
    g_inv_.DeleteAthenaArray();
    cons_.DeleteAthenaArray();
  }
}
