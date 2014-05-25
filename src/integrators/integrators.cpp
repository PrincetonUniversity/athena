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

#include <sstream>
#include <iostream>
#include <string>
#include <stdexcept>

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh.hpp"
#include "../fluid.hpp"
#include "integrators.hpp"

//======================================================================================
/*! \file integrators.cpp
 *  \brief 
 *====================================================================================*/

// constructor

FluidIntegrator::FluidIntegrator(Fluid *pf)
{
  pparent_fluid = pf;

// Allocate memory for scratch vectors

  int ncells1 = pf->pparent_block->block_size.nx1 + 2*(NGHOST);
  wl_.NewAthenaArray(NVAR,ncells1);
  wr_.NewAthenaArray(NVAR,ncells1);
  flx_.NewAthenaArray(NVAR,ncells1);
}

// destructor

FluidIntegrator::~FluidIntegrator()
{
  wl_.DeleteAthenaArray();
  wr_.DeleteAthenaArray();
  flx_.DeleteAthenaArray();
}
