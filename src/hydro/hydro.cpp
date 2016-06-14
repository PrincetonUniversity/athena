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
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../field/field.hpp"
#include "srcterms/hydro_srcterms.hpp"

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

  flux[X1DIR].NewAthenaArray(NHYDRO,ncells3,ncells2,ncells1+1);
  if (pmy_block->block_size.nx2 > 1) 
    flux[X2DIR].NewAthenaArray(NHYDRO,ncells3,ncells2+1,ncells1);
  if (pmy_block->block_size.nx3 > 1) 
    flux[X3DIR].NewAthenaArray(NHYDRO,ncells3+1,ncells2,ncells1);

  // Allocate memory for metric
  // TODO: this should only be done if we are in GR
  g.NewAthenaArray(NMETRIC, ncells1);
  g_inv.NewAthenaArray(NMETRIC, ncells1);

// Allocate memory for internal hydro output variables (if needed)

  ifov.NewAthenaArray(NIFOV,ncells3,ncells2,ncells1);

// Allocate memory for scratch arrays

  int nthreads = pmy_block->pmy_mesh->GetNumMeshThreads();
  dt1_.NewAthenaArray(nthreads,ncells1);
  dt2_.NewAthenaArray(nthreads,ncells1);
  dt3_.NewAthenaArray(nthreads,ncells1);
  wl_.NewAthenaArray(nthreads,(NWAVE),ncells1);
  wr_.NewAthenaArray(nthreads,(NWAVE),ncells1);
  flx_.NewAthenaArray(nthreads,(NWAVE),ncells1);
  x1face_area_.NewAthenaArray(nthreads,ncells1+1);
  if(pmy_block->block_size.nx2 > 1) {
    x2face_area_.NewAthenaArray(nthreads,ncells1);
    x2face_area_p1_.NewAthenaArray(nthreads,ncells1);
  }
  if(pmy_block->block_size.nx3 > 1) {
    x3face_area_.NewAthenaArray(nthreads,ncells1);
    x3face_area_p1_.NewAthenaArray(nthreads,ncells1);
  }
  cell_volume_.NewAthenaArray(nthreads,ncells1);
  if (MAGNETIC_FIELDS_ENABLED && RELATIVISTIC_DYNAMICS)  // only used in (SR/GR)MHD
  {
    bb_normal_.NewAthenaArray(ncells1);
    lambdas_p_l_.NewAthenaArray(ncells1);
    lambdas_m_l_.NewAthenaArray(ncells1);
    lambdas_p_r_.NewAthenaArray(ncells1);
    lambdas_m_r_.NewAthenaArray(ncells1);
  }
  if (GENERAL_RELATIVITY)  // only used in GR
  {
    g_.NewAthenaArray(NMETRIC,ncells1);
    gi_.NewAthenaArray(NMETRIC,ncells1);
    cons_.NewAthenaArray(NWAVE,ncells1);
  }

// Construct ptrs to objects of various classes needed to integrate hydro/MHD eqns 

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
  ifov.DeleteAthenaArray();

  flux[X1DIR].DeleteAthenaArray();
  if (pmy_block->block_size.nx2 > 1) flux[X2DIR].DeleteAthenaArray();
  if (pmy_block->block_size.nx3 > 1) flux[X3DIR].DeleteAthenaArray();

  dt1_.DeleteAthenaArray();
  dt2_.DeleteAthenaArray();
  dt3_.DeleteAthenaArray();
  wl_.DeleteAthenaArray();
  wr_.DeleteAthenaArray();
  flx_.DeleteAthenaArray();
  x1face_area_.DeleteAthenaArray();
  if(pmy_block->block_size.nx2 > 1) {
    x2face_area_.DeleteAthenaArray();
    x2face_area_p1_.DeleteAthenaArray();
  }
  if(pmy_block->block_size.nx3 > 1) {
    x3face_area_.DeleteAthenaArray();
    x3face_area_p1_.DeleteAthenaArray();
  }
  cell_volume_.DeleteAthenaArray();
  if (MAGNETIC_FIELDS_ENABLED && RELATIVISTIC_DYNAMICS)  // only used in (SR/GR)MHD
  {
    bb_normal_.DeleteAthenaArray();
    lambdas_p_l_.DeleteAthenaArray();
    lambdas_m_l_.DeleteAthenaArray();
    lambdas_p_r_.DeleteAthenaArray();
    lambdas_m_r_.DeleteAthenaArray();
  }
  if (GENERAL_RELATIVITY)
  {
    g_.DeleteAthenaArray();
    gi_.DeleteAthenaArray();
    cons_.DeleteAthenaArray();
  }

  delete psrc;
}
