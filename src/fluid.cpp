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
 * You should have received a copy of GNU GPL in the file LICENSE included in
 * the code distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

#include <iostream>
#include <string>

#include "athena.hpp"
#include "athena_arrays.hpp"
#include "parameter_input.hpp"
#include "mesh.hpp"
#include "fluid.hpp"

//======================================================================================
//! \file fluid.cpp
//  \brief implementation of functions in class Fluid
//======================================================================================

// constructor, initializes data structures and parameters, calls problem generator

Fluid::Fluid(ParameterInput *pin, Block *pb)
{
  pmy_block = pb;
// Read adiabatic index from input file

  gamma_ = pin->GetReal("fluid","gamma");

// initialize time, timestep

  time = 0.0;
  dt   = 1.0;

// Allocate memory for primitive/conserved variables

  int ncells1 = pmy_block->block_size.nx1 + 2*(NGHOST);
  int ncells2 = pmy_block->block_size.nx2 + 2*(NGHOST);
  int ncells3 = pmy_block->block_size.nx3 + 2*(NGHOST);

  u.NewAthenaArray(NVAR,ncells3,ncells2,ncells1);
  w.NewAthenaArray(NVAR,ncells3,ncells2,ncells1);

// Allocate memory for primitive/conserved variables at half-time step, and scratch

  u1_.NewAthenaArray(NVAR,ncells3,ncells2,ncells1);
  w1_.NewAthenaArray(NVAR,ncells3,ncells2,ncells1);

  wl_.NewAthenaArray(NVAR,ncells1);
  wr_.NewAthenaArray(NVAR,ncells1);
  flx_.NewAthenaArray(NVAR,ncells1);

}

// destructor

Fluid::~Fluid()
{
  u.DeleteAthenaArray();
  w.DeleteAthenaArray();
  
  u1_.DeleteAthenaArray();
  w1_.DeleteAthenaArray();
  wl_.DeleteAthenaArray();
  wr_.DeleteAthenaArray();
  flx_.DeleteAthenaArray();
}
