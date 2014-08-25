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

// Primary header
#include "bvals.hpp"

// C++ headers
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena headers
#include "../athena.hpp"         // Real
#include "../athena_arrays.hpp"  // AthenaArray
#include "../mesh.hpp"           // Block

//======================================================================================
/*! \file bvals.cpp
 *  \brief boundary conditions for fluid (quantities in ghost zones) on each edge
 *====================================================================================*/

// constructor


FluidBoundaryConditions::FluidBoundaryConditions(Block *pb)
{
  pmy_block = pb;
  std::stringstream msg;

// Set BC function pointers for each of the 6 boundaries in turn -----------------------
// Inner x1

  switch(pb->block_bndry.ix1_bc){
    case 1:
      FluidInnerX1_ = ReflectInnerX1;
    break;
    case 2:
      FluidInnerX1_ = OutflowInnerX1;
    break;
    case 4:
      FluidInnerX1_ = PeriodicInnerX1;
    break;
    default:
      msg << "### FATAL ERROR in FluidBoundaryConditions constructor" << std::endl
          << "Flag ix1_bc=" << pb->block_bndry.ix1_bc << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
    break;
   }

// Outer x1

  switch(pb->block_bndry.ox1_bc){
    case 1:
      FluidOuterX1_ = ReflectOuterX1;
    break;
    case 2:
      FluidOuterX1_ = OutflowOuterX1;
    break;
    case 4:
      FluidOuterX1_ = PeriodicOuterX1;
    break;
    default:
      msg << "### FATAL ERROR in FluidBoundaryConditions constructor" << std::endl
          << "Flag ox1_bc=" << pb->block_bndry.ox1_bc << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
    break;
  }

// Inner x2

  if (pb->block_size.nx2 > 1) {
    switch(pb->block_bndry.ix2_bc){
      case 1:
        FluidInnerX2_ = ReflectInnerX2;
      break;
      case 2:
        FluidInnerX2_ = OutflowInnerX2;
      break;
      case 4:
        FluidInnerX2_ = PeriodicInnerX2;
      break;
      default:
        msg << "### FATAL ERROR in FluidBoundaryConditions constructor" << std::endl
            << "Flag ix2_bc=" << pb->block_bndry.ix2_bc << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
      break;
     }

// Outer x2

    switch(pb->block_bndry.ox2_bc){
      case 1:
        FluidOuterX2_ = ReflectOuterX2;
      break;
      case 2:
        FluidOuterX2_ = OutflowOuterX2;
      break;
      case 4:
        FluidOuterX2_ = PeriodicOuterX2;
      break;
      default:
        msg << "### FATAL ERROR in FluidBoundaryConditions constructor" << std::endl
            << "Flag ox2_bc=" << pb->block_bndry.ox2_bc << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
      break;
    }
  }

// Inner x3

  if (pb->block_size.nx3 > 1) {
    switch(pb->block_bndry.ix3_bc){
      case 1:
        FluidInnerX3_ = ReflectInnerX3;
      break;
      case 2:
        FluidInnerX3_ = OutflowInnerX3;
      break;
      case 4:
        FluidInnerX3_ = PeriodicInnerX3;
      break;
      default:
        msg << "### FATAL ERROR in FluidBoundaryConditions constructor" << std::endl
            << "Flag ix3_bc=" << pb->block_bndry.ix3_bc << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
      break;
     }

// Outer x3

    switch(pb->block_bndry.ox3_bc){
      case 1:
        FluidOuterX3_ = ReflectOuterX3;
      break;
      case 2:
        FluidOuterX3_ = OutflowOuterX3;
      break;
      case 4:
        FluidOuterX3_ = PeriodicOuterX3;
      break;
      default:
        msg << "### FATAL ERROR in FluidBoundaryConditions constructor" << std::endl
            << "Flag ox3_bc=" << pb->block_bndry.ox3_bc << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
      break;
    }
  }

}

// destructor

FluidBoundaryConditions::~FluidBoundaryConditions()
{
}

//--------------------------------------------------------------------------------------
/*! \fn void ApplyBoundaryConditions()
 *  \brief Calls BC functions using appropriate function pointers to set ghost zones.  
 */

void FluidBoundaryConditions::ApplyBoundaryConditions(AthenaArray<Real> &a)
{

// Boundary Conditions in x1-direction

  (*(FluidInnerX1_))(pmy_block, a);
  (*(FluidOuterX1_))(pmy_block, a);

// Boundary Conditions in x2-direction 

  if (pmy_block->block_size.nx2 > 1){

    (*(FluidInnerX2_))(pmy_block, a);
    (*(FluidOuterX2_))(pmy_block, a);

  }

// Boundary Conditions in x3-direction 

  if (pmy_block->block_size.nx3 > 1){

    (*(FluidInnerX3_))(pmy_block, a);
    (*(FluidOuterX3_))(pmy_block, a);

  }

  return;
}
