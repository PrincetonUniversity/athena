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
#include <math.h>
#include <float.h>

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh.hpp"
#include "geometry.hpp"

//======================================================================================
//! \file geometry.cpp
//  \brief implementation of functions in class Geometry
//======================================================================================

//namespace COORDINATE_SYSTEM {

// constructor

Geometry::Geometry(Block *pb)
{
  pparent_block = pb;

  int ncells1 = pb->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pb->block_size.nx2 > 1) ncells2 = pb->block_size.nx2 + 2*(NGHOST);
  if (pb->block_size.nx3 > 1) ncells3 = pb->block_size.nx3 + 2*(NGHOST);

  dx1v.NewAthenaArray(ncells1);
  dx2v.NewAthenaArray(ncells2);
  dx3v.NewAthenaArray(ncells3);

  x1v.NewAthenaArray(ncells1);
  x2v.NewAthenaArray(ncells2);
  x3v.NewAthenaArray(ncells3);

  face_area.NewAthenaArray(ncells1);
  cell_volume.NewAthenaArray(ncells1);
}

// destructor

Geometry::~Geometry()
{
}

//} // end namespace COORDINATE_SYSTEM
