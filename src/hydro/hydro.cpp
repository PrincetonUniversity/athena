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
//! \file hydro.cpp
//  \brief implementation of functions in class Hydro
//======================================================================================

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../eos/eos.hpp"
#include "srcterms/srcterms.hpp"
#include "fluxes/fluxes.hpp"
#include "../mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../field/field.hpp"

// this class header
#include "hydro.hpp"

// constructor, initializes data structures and parameters

Hydro::Hydro(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block = pmb;

// Allocate memory for primitive/conserved variables

  int ncells1 = pmy_block->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmy_block->block_size.nx2 > 1) ncells2 = pmy_block->block_size.nx2 + 2*(NGHOST);
  if (pmy_block->block_size.nx3 > 1) ncells3 = pmy_block->block_size.nx3 + 2*(NGHOST);

  u.NewAthenaArray(NHYDRO,ncells3,ncells2,ncells1);
  w.NewAthenaArray(NHYDRO,ncells3,ncells2,ncells1);

// Allocate memory for primitive/conserved variables at intermediate-time step

  u1.NewAthenaArray(NHYDRO,ncells3,ncells2,ncells1);
  w1.NewAthenaArray(NHYDRO,ncells3,ncells2,ncells1);

  flux[x1face].NewAthenaArray(NHYDRO,ncells3,ncells2,ncells1+1);
  if (pmy_block->block_size.nx2 > 1) 
    flux[x2face].NewAthenaArray(NHYDRO,ncells3,ncells2+1,ncells1);
  if (pmy_block->block_size.nx3 > 1) 
    flux[x3face].NewAthenaArray(NHYDRO,ncells3+1,ncells2,ncells1);

  // Allocate memory for metric
  // TODO: this should only be done if we are in GR
  g.NewAthenaArray(NMETRIC, ncells1);
  g_inv.NewAthenaArray(NMETRIC, ncells1);

// Allocate memory for scratch arrays

  int nthreads = pmy_block->pmy_mesh->GetNumMeshThreads();
  dt1_.NewAthenaArray(nthreads,ncells1);
  dt2_.NewAthenaArray(nthreads,ncells1);
  dt3_.NewAthenaArray(nthreads,ncells1);

// Allocate memory for internal hydro output variables (if needed)

  ifov.NewAthenaArray(NIFOV,ncells3,ncells2,ncells1);

// Construct ptrs to objects of various classes needed to integrate hydro/MHD eqns 

  pflux = new HydroFluxes(this,pin);
  psrc  = new HydroSourceTerms(this,pin);
}

// destructor

Hydro::~Hydro()
{
  u.DeleteAthenaArray();
  w.DeleteAthenaArray();
  u1.DeleteAthenaArray();
  w1.DeleteAthenaArray();
  g.DeleteAthenaArray();
  g_inv.DeleteAthenaArray();

  flux[x1face].DeleteAthenaArray();
  if (pmy_block->block_size.nx2 > 1) flux[x2face].DeleteAthenaArray();
  if (pmy_block->block_size.nx3 > 1) flux[x3face].DeleteAthenaArray();

  dt1_.DeleteAthenaArray();
  dt2_.DeleteAthenaArray();
  dt3_.DeleteAthenaArray();

  ifov.DeleteAthenaArray();

  delete pflux;
  delete psrc;
}
