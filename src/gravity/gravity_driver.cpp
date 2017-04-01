//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gravity_driver.cpp
//  \brief implementation of functions in class GravityDriver

// Athena++ headers
#include "mggravity.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../parameter_input.hpp"
#include "../multigrid/multigrid.hpp"

// constructor, initializes data structures and parameters

GravityDriver::GravityDriver(Mesh *pm, ParameterInput *pin)
{
  pblock_=pm->pblock;
  if(pblock_->block_size.nx1!=pblock_->block_size.nx2
  || pblock_->block_size.nx1!=pblock_->block_size.nx3) {
    std::stringstream msg;
    msg << "### FATAL ERROR in GravityDriver::GravityDriver" << std::endl
        << "The Multigrid solver requires logically cubic MeshBlocks." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    break;
  }
  if(pblock_->block_size.nx2==1 || pblock_->block_size.nx3==1 ) {
    std::stringstream msg;
    msg << "### FATAL ERROR in GravityDriver::GravityDriver" << std::endl
        << "Currently the Multigrid solver works only in 3D." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    break;
  }
  if(pm->use_meshgen_fn_[X1DIR]==true || pm->use_meshgen_fn_[X2DIR]==true
  || pm->use_meshgen_fn_[X3DIR]==true) {
    std::stringstream msg;
    msg << "### FATAL ERROR in GravityDriver::GravityDriver" << std::endl
        << "Non-uniform mesh spacing is not supported." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    break;
  }
  Real dx=pblock_->pcoord->dx1f(0);
  if(dx!=pblock_->pcoord->dx2f(0) || dx!=pblock_->pcoord->dx3f(0)) {
    std::stringstream msg;
    msg << "### FATAL ERROR in GravityDriver::GravityDriver" << std::endl
        << "The cell size must be cubic." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    break;
  }
  // count multigrid levels
  nlev_=0;
  for(int l=0; l<20; l++) {
    if((1<<l) == pblock_->block_size.nx1) {
      nlev_=l+1;
      break;
    }
  }
  if(nlev_==0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in GravityDriver::GravityDriver" << std::endl
        << "The MeshBlock size must be power of two." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    break;
  }
  // count multigrid levels
  int nrlx=0, nrly=0, nrlz=0;
  for(int l=0; l<20; l++) {
    if((1<<l) == pm->nrbx1)
      nrlx=l+1;
    if((1<<l) == pm->nrbx2)
      nrly=l+1;
    if((1<<l) == pm->nrbx3)
      nrlz=l+1;
  }
  if(nrlx==0 || nrly==0 || nrlz==0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in GravityDriver::GravityDriver" << std::endl
        << "The root grid size must be power of 2." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    break;
  }

  umode_=0; // SMR - 2 boundary calls per smoothing
  if(pm->multilevel==false) umode_=1; // Uniform - 1 boundary call per smoothing

  mode_=0; // V-cycle
  
  Real dx=(mesh_size.x1max-mesh_size.x1min)/mesh_size.nx1;
  mgroot_ = new MGGravity(pm->nrbx1, pm->nrbx2, pm->nrbx3, 1, dx);
}

// destructor

GravityDriver::~GravityDriver()
{
}

