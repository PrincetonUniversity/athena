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
//! \file integrators.cpp
//  \brief 
//======================================================================================

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../hydro.hpp"
#include "../../mesh.hpp"
#include "../../parameter_input.hpp" 

// this class header
#include "hydro_integrator.hpp"

// constructor

HydroIntegrator::HydroIntegrator(Hydro *phydro, ParameterInput *pin)
{
  pmy_hydro = phydro;

// Allocate memory for scratch vectors

  int nthreads = phydro->pmy_block->pmy_mesh->GetNumMeshThreads();
  int ncells1 = phydro->pmy_block->block_size.nx1 + 2*(NGHOST);
  int ncells2 = phydro->pmy_block->block_size.nx2 + 2*(NGHOST);

  wl_.NewAthenaArray(nthreads,(NWAVE),ncells1);
  wr_.NewAthenaArray(nthreads,(NWAVE),ncells1);
  flx_.NewAthenaArray(nthreads,(NWAVE),ncells1);
  x1face_area_.NewAthenaArray(nthreads,ncells1+1);
  if(pmy_hydro->pmy_block->block_size.nx2 > 1) {
    x2face_area_.NewAthenaArray(nthreads,ncells1);
    x2face_area_p1_.NewAthenaArray(nthreads,ncells1);
  }
  if(pmy_hydro->pmy_block->block_size.nx3 > 1) {
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
}

// destructor

HydroIntegrator::~HydroIntegrator()
{
  wl_.DeleteAthenaArray();
  wr_.DeleteAthenaArray();
  flx_.DeleteAthenaArray();
  x1face_area_.DeleteAthenaArray();
  if(pmy_hydro->pmy_block->block_size.nx2 > 1) {
    x2face_area_.DeleteAthenaArray();
    x2face_area_p1_.DeleteAthenaArray();
  }
  if(pmy_hydro->pmy_block->block_size.nx3 > 1) {
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
}
