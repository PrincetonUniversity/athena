//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file multigrid_driver.cpp
//  \brief implementation of functions in class MultigridDriver

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../parameter_input.hpp"
#include "multigrid.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// constructor, initializes data structures and parameters

MultigridDriver::MultigridDriver(Mesh *pm, MeshBlock *pmb, MGBoundaryFunc_t *MGBoundary,
                                 ParameterInput *pin)
{
  pmy_mesh_=pm;
  if(pblock_->block_size.nx1!=pblock_->block_size.nx2
  || pblock_->block_size.nx1!=pblock_->block_size.nx3) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MultigridDriver::MultigridDriver" << std::endl
        << "The Multigrid solver requires logically cubic MeshBlocks." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    break;
  }
  if(pblock_->block_size.nx2==1 || pblock_->block_size.nx3==1 ) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MultigridDriver::MultigridDriver" << std::endl
        << "Currently the Multigrid solver works only in 3D." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    break;
  }
  if(pm->use_meshgen_fn_[X1DIR]==true || pm->use_meshgen_fn_[X2DIR]==true
  || pm->use_meshgen_fn_[X3DIR]==true) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MultigridDriver::MultigridDriver" << std::endl
        << "Non-uniform mesh spacing is not supported." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    break;
  }
  Real dx=pblock_->pcoord->dx1f(0);
  if(dx!=pblock_->pcoord->dx2f(0) || dx!=pblock_->pcoord->dx3f(0)) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MultigridDriver::MultigridDriver" << std::endl
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
    msg << "### FATAL ERROR in MultigridDriver::MultigridDriver" << std::endl
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
    msg << "### FATAL ERROR in MultigridDriver::MultigridDriver" << std::endl
        << "The root grid size must be power of 2." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    break;
  }
  nrootlevel_=std::min(nrlx,std::min(nrly,nrlz));

  fperiodic_=false;
  for(int i=0; i<6; i++) {
    if(pm->mesh_bcs[i]==PERIODIC_BNDRY) fperidoc_=true;
    MGBoundaryFunction_[i]=pm->MGBoundaryFunction_[i];
  }

  nblocks_=0;
  MeshBlocks *pb=pmb;
  while(pb!=NULL) {
    nblocks_++;
    pb=pb->next;
  }

  mode_=0; // 0: FMG+V(1,1), 1: FMG+F(0,1), 2: V(1,1)

  mgroot_=NULL;

  // Setting up the MPI information
  // *** this part should be modified when dedicate processes are allocated ***
  pblock_=pmb;
  nranks_  = Globals::nranks;
  nslist_  = new int[nranks_];
  nblist_  = new int[nranks_];
  nvlist_  = new int[nranks_];
  nvslist_ = new int[nranks_];
#ifdef MPI_PARALLEL
  MPI_Comm_dup(MPI_COMM_WORKD, &MPI_COMM_MULTIGRID);
#endif
  for(int n=0; n<nranks_; n++) {
    nslist_[n]  = pm->nslist[n];
    nblist_[n]  = pm->nblist[n];
    nvlist_[n]  = nblist_[n]*nvar_;
    nvslist_[n] = nslist_[n]*nvar_;
  }
}

// destructor

MultigridDriver::~MultigridDriver()
{
  delete mgroot_;
  delete [] nslist_;
  delete [] nblist_;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SetupMultigrid(void)
//  \brief initialize the source assuming that the source terms are already loaded

MultigridDriver::SetupMultigrid(void)
{
  Mesh *pm=pmy_mesh_;
  MeshBlocks *pb=pblock_;
  Real *tsrc=new Real[pm->nbtotal*nvar_];

  if(fperiodic_)  // periodic boundary - calc and subtract average
    while(pb!=NULL) {
      for(int v=0; v<nvar_; v++)
      tsrc[pb->gid*nvar_+v]=pb->pmggrav->CalculateTotalSource(n);
      pb=pb->next;
    }
#ifdef MPI_PARALLEL
    MPI_Allgatherv(MPI_IN_PLACE, nblist_[Globals::my_rank]*nvar_, MPI_ATHENA_REAL,
                   tsrc, nvlist, nvslist, MPI_ATHENA_REAL, MPI_COMM_GRAVITY);
#endif
    for(int v=0; v<nvar_; v++) {
      Real total=0.0;
      for(int n=0; n<Globals::nranks; n++)
        total+=tsrc[n*nvar_+v];
      Real ave=total/(pm->mesh_size.(x1max-x1min)*pm->mesh_size.(x2max-x2min)
                                                 *pm->mesh_size.(x3max-x3min));
      for(int n=0; n<pm->nbtotal; n++)
        tsrc[n*nvar_+v]-=ave;
      pb=pblock_;
      while(pb!=NULL) {
        pb->pmggrav->SubtractAverageSource(v, ave);
        pb=pb->next;
      }
    }
  }
  if(mode_<=1) { // FMG
    pb=pblock_;
    while(pb!=NULL) {
      pb->pmggrav->RestrictFMGSource();
      pb=pb->next;
    }
    // construct the root grid
    if(pmy_mesh->multilevel) { 
      // *** implement later
    }
    else {
    }
  }
  delete [] tsrc;
  return;
}

