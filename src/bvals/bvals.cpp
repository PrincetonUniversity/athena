
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
#include "bvals.hpp"

// C++ headers
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena headers
#include "../athena.hpp"          // Real
#include "../athena_arrays.hpp"   // AthenaArray
#include "../mesh.hpp"            // MeshBlock
#include "../fluid/fluid.hpp"     // Fluid
#include "../parameter_input.hpp" // ParameterInput

//======================================================================================
//! \file bvals.cpp
//  \brief implements functions that initialize/apply BCs on each edge
//======================================================================================

// BoundaryValues constructor - sets functions for the appropriate
// boundary conditions at each of the 6 edges of a MeshBlock

BoundaryValues::BoundaryValues(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_mblock_ = pmb;

// Set BC functions for each of the 6 boundaries in turn -------------------------------
// Inner x1

  switch(pmb->block_bcs.ix1_bc){
    case -1:
      FluidInnerX1_ = NeighborInnerX1;
      BFieldInnerX1_ = NeighborInnerX1;
      break;
    case 1:
      FluidInnerX1_ = ReflectInnerX1;
      BFieldInnerX1_ = ReflectInnerX1;
      break;
    case 2:
      FluidInnerX1_ = OutflowInnerX1;
      BFieldInnerX1_ = OutflowInnerX1;
      break;
    case 3: // do nothing, useful for user-enrolled BCs
      break;
    case 4:
      FluidInnerX1_ = PeriodicInnerX1;
      BFieldInnerX1_ = PeriodicInnerX1;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in FluidBCs constructor" << std::endl
          << "Flag ix1_bc=" << pmb->block_bcs.ix1_bc << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
   }

// Outer x1

  switch(pmb->block_bcs.ox1_bc){
    case -1:
      FluidOuterX1_ = NeighborOuterX1;
      BFieldOuterX1_ = NeighborOuterX1;
      break;
    case 1:
      FluidOuterX1_ = ReflectOuterX1;
      BFieldOuterX1_ = ReflectOuterX1;
      break;
    case 2:
      FluidOuterX1_ = OutflowOuterX1;
      BFieldOuterX1_ = OutflowOuterX1;
      break;
    case 3: // do nothing, useful for user-enrolled BCs
      break;
    case 4:
      FluidOuterX1_ = PeriodicOuterX1;
      BFieldOuterX1_ = PeriodicOuterX1;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in FluidBCs constructor" << std::endl
          << "Flag ox1_bc=" << pmb->block_bcs.ox1_bc << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }

// Inner x2

  if (pmb->block_size.nx2 > 1) {
    switch(pmb->block_bcs.ix2_bc){
      case -1:
        FluidInnerX2_ = NeighborInnerX2;
        BFieldInnerX2_ = NeighborInnerX2;
        break;
      case 1:
        FluidInnerX2_ = ReflectInnerX2;
        BFieldInnerX2_ = ReflectInnerX2;
        break;
      case 2:
        FluidInnerX2_ = OutflowInnerX2;
        BFieldInnerX2_ = OutflowInnerX2;
        break;
      case 3: // do nothing, useful for user-enrolled BCs
        break;
      case 4:
        FluidInnerX2_ = PeriodicInnerX2;
        BFieldInnerX2_ = PeriodicInnerX2;
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in FluidBCs constructor" << std::endl
            << "Flag ix2_bc=" << pmb->block_bcs.ix2_bc << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
     }

// Outer x2

    switch(pmb->block_bcs.ox2_bc){
      case -1:
        FluidOuterX2_ = NeighborOuterX2;
        BFieldOuterX2_ = NeighborOuterX2;
        break;
      case 1:
        FluidOuterX2_ = ReflectOuterX2;
        BFieldOuterX2_ = ReflectOuterX2;
        break;
      case 2:
        FluidOuterX2_ = OutflowOuterX2;
        BFieldOuterX2_ = OutflowOuterX2;
        break;
      case 3: // do nothing, useful for user-enrolled BCs
        break;
      case 4:
        FluidOuterX2_ = PeriodicOuterX2;
        BFieldOuterX2_ = PeriodicOuterX2;
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in FluidBCs constructor" << std::endl
            << "Flag ox2_bc=" << pmb->block_bcs.ox2_bc << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
    }
  }

// Inner x3

  if (pmb->block_size.nx3 > 1) {
    switch(pmb->block_bcs.ix3_bc){
      case -1:
        FluidInnerX3_ = NeighborInnerX3;
        BFieldInnerX3_ = NeighborInnerX3;
        break;
      case 1:
        FluidInnerX3_ = ReflectInnerX3;
        BFieldInnerX3_ = ReflectInnerX3;
        break;
      case 2:
        FluidInnerX3_ = OutflowInnerX3;
        BFieldInnerX3_ = OutflowInnerX3;
        break;
      case 3: // do nothing, useful for user-enrolled BCs
        break;
      case 4:
        FluidInnerX3_ = PeriodicInnerX3;
        BFieldInnerX3_ = PeriodicInnerX3;
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in FluidBCs constructor" << std::endl
            << "Flag ix3_bc=" << pmb->block_bcs.ix3_bc << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
     }

// Outer x3

    switch(pmb->block_bcs.ox3_bc){
      case -1:
        FluidOuterX3_ = NeighborOuterX3;
        BFieldOuterX3_ = NeighborOuterX3;
        break;
      case 1:
        FluidOuterX3_ = ReflectOuterX3;
        BFieldOuterX3_ = ReflectOuterX3;
        break;
      case 2:
        FluidOuterX3_ = OutflowOuterX3;
        BFieldOuterX3_ = OutflowOuterX3;
        break;
      case 3: // do nothing, useful for user-enrolled BCs
        break;
      case 4:
        FluidOuterX3_ = PeriodicOuterX3;
        BFieldOuterX3_ = PeriodicOuterX3;
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in FluidBCs constructor" << std::endl
            << "Flag ox3_bc=" << pmb->block_bcs.ox3_bc << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
    }
  }

}

// destructor

BoundaryValues::~BoundaryValues()
{
}

//--------------------------------------------------------------------------------------
//! \fn
//  \brief

void BoundaryValues::EnrollFluidBoundaryFunction(enum EdgeNames edge, BValFluid_t my_bc)
{
  switch(edge){
  case inner_x1:
    FluidInnerX1_ = my_bc;
    break;
  case outer_x1:
    FluidOuterX1_ = my_bc;
    break;
  case inner_x2:
    FluidInnerX2_ = my_bc;
    break;
  case outer_x2:
    FluidOuterX2_ = my_bc;
    break;
  case inner_x3:
    FluidInnerX3_ = my_bc;
    break;
  case outer_x3:
    FluidOuterX3_ = my_bc;
    break;
  default:
    std::stringstream msg;
    msg << "### FATAL ERROR in EnrollFluidBoundaryCondition function" << std::endl
        << "EdgeName = " << edge << " not valid" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  return;
}

//--------------------------------------------------------------------------------------
//! \fn
//  \brief

void BoundaryValues::EnrollFieldBoundaryFunction(enum EdgeNames edge,BValBField_t my_bc)
{
  switch(edge){
  case inner_x1:
    BFieldInnerX1_ = my_bc;
    break;
  case outer_x1:
    BFieldOuterX1_ = my_bc;
    break;
  case inner_x2:
    BFieldInnerX2_ = my_bc;
    break;
  case outer_x2:
    BFieldOuterX2_ = my_bc;
    break;
  case inner_x3:
    BFieldInnerX3_ = my_bc;
    break;
  case outer_x3:
    BFieldOuterX3_ = my_bc;
    break;
  default:
    std::stringstream msg;
    msg << "### FATAL ERROR in EnrollFieldBoundaryCondition function" << std::endl
        << "EdgeName = " << edge << " not valid" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void FluidBCs::ApplyFluidBCs(AthenaArray<Real> &a)
//  \brief Calls BC functions using appropriate function pointers to set ghost zones.  

template<typename T> void BoundaryValues::ApplyBVals(T &input)
{
// Boundary Conditions in x1-direction

  BValsInnerX1_(pmy_mblock_, input);
  BValsOuterX1_(pmy_mblock_, input);

// Boundary Conditions in x2-direction 

  if (pmy_mblock_->block_size.nx2 > 1){
    BValsInnerX2_(pmy_mblock_, input);
    BValsOuterX2_(pmy_mblock_, input);
  }

// Boundary Conditions in x3-direction 

  if (pmy_mblock_->block_size.nx3 > 1){
    BValsInnerX3_(pmy_mblock_, input);
    BValsOuterX3_(pmy_mblock_, input);
  }

  return;
}

template void BoundaryValues::ApplyBVals< AthenaArray<Real> >(AthenaArray<Real> &input);
template void BoundaryValues::ApplyBVals<InterfaceField>(InterfaceField &input);

template<> void BoundaryValues::BValsInnerX1_< AthenaArray<Real> >
  (MeshBlock *pmb, AthenaArray<Real> &input) {FluidInnerX1_(pmb,input);}
template<> void BoundaryValues::BValsInnerX2_< AthenaArray<Real> >
  (MeshBlock *pmb, AthenaArray<Real> &input) {FluidInnerX2_(pmb,input);}
template<> void BoundaryValues::BValsInnerX3_< AthenaArray<Real> >
  (MeshBlock *pmb, AthenaArray<Real> &input) {FluidInnerX3_(pmb,input);}
template<> void BoundaryValues::BValsOuterX1_< AthenaArray<Real> >
  (MeshBlock *pmb, AthenaArray<Real> &input) {FluidOuterX1_(pmb,input);}
template<> void BoundaryValues::BValsOuterX2_< AthenaArray<Real> >
  (MeshBlock *pmb, AthenaArray<Real> &input) {FluidOuterX2_(pmb,input);}
template<> void BoundaryValues::BValsOuterX3_< AthenaArray<Real> >
  (MeshBlock *pmb, AthenaArray<Real> &input) {FluidOuterX3_(pmb,input);}

template<> void BoundaryValues::BValsInnerX1_< InterfaceField >
  (MeshBlock *pmb, InterfaceField &input) {BFieldInnerX1_(pmb,input);}
template<> void BoundaryValues::BValsInnerX2_< InterfaceField >
  (MeshBlock *pmb, InterfaceField &input) {BFieldInnerX2_(pmb,input);}
template<> void BoundaryValues::BValsInnerX3_< InterfaceField >
  (MeshBlock *pmb, InterfaceField &input) {BFieldInnerX3_(pmb,input);}
template<> void BoundaryValues::BValsOuterX1_< InterfaceField >
  (MeshBlock *pmb, InterfaceField &input) {BFieldOuterX1_(pmb,input);}
template<> void BoundaryValues::BValsOuterX2_< InterfaceField >
  (MeshBlock *pmb, InterfaceField &input) {BFieldOuterX2_(pmb,input);}
template<> void BoundaryValues::BValsOuterX3_< InterfaceField >
  (MeshBlock *pmb, InterfaceField &input) {BFieldOuterX3_(pmb,input);}
