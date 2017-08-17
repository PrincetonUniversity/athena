//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gravity_driver.cpp
//  \brief implementation of functions in class GravityDriver

// Athena++ headers
#include "mggravity.hpp"
#include "gravity.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../coordinates/coordinates.hpp"
#include "../parameter_input.hpp"
#include "../multigrid/multigrid.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

class MeshBlock;



//----------------------------------------------------------------------------------------
//! \fn GravityDriver::GravityDriver(Mesh *pm, MGBoundaryFunc_t *MGBoundary,
//                                   ParameterInput *pin)
//  \brief GravityDriver constructor

GravityDriver::GravityDriver(Mesh *pm, MGBoundaryFunc_t *MGBoundary, ParameterInput *pin)
 : MultigridDriver(pm, MGBoundary, 1)
{
  four_pi_G_=pin->GetOrAddReal("problem", "four_pi_G", 1.0); // default: 4piG=1
  // Allocate multigrid objects
  RegionSize root_size=pmy_mesh_->mesh_size;
  root_size.nx1=pmy_mesh_->nrbx1;
  root_size.nx2=pmy_mesh_->nrbx2;
  root_size.nx3=pmy_mesh_->nrbx3;
  LogicalLocation lroot;
  lroot.lx1=0, lroot.lx2=0, lroot.lx3=0, lroot.level=0;
  mgroot_= new MGGravity(this,lroot,-1,-1,root_size,MGBoundary,pmy_mesh_->mesh_bcs,true);
  pmg_=NULL;
  Multigrid *pfirst;
  int nbs=nslist_[Globals::my_rank];
  int nbe=nbs+nblist_[Globals::my_rank]-1;
  RegionSize block_size;
  block_size.nx1=pmy_mesh_->mesh_size.nx1/pmy_mesh_->nrbx1;
  block_size.nx2=pmy_mesh_->mesh_size.nx2/pmy_mesh_->nrbx2;
  block_size.nx3=pmy_mesh_->mesh_size.nx3/pmy_mesh_->nrbx3;
  for(int i=nbs;i<=nbe;i++) {
    enum BoundaryFlag block_bcs[6];
    pmy_mesh_->SetBlockSizeAndBoundaries(pmy_mesh_->loclist[i], block_size, block_bcs);
    Multigrid *nmg=new MGGravity(this, pmy_mesh_->loclist[i], i, i-nbs, block_size,
                                 MGBoundary, block_bcs, false);
    AddMultigrid(nmg);
  }
}


//----------------------------------------------------------------------------------------
//! \fn void GravityDriver::LoadSourceAndData(void)
//  \brief load the sourterm and initial guess (if needed)

void GravityDriver::LoadSourceAndData(void)
{
  Multigrid *pmggrav=pmg_;
  while(pmggrav!=NULL) {
    MeshBlock *pmb=pmy_mesh_->FindMeshBlock(pmggrav->gid_);
    if(pmb!=NULL) {
      pmggrav->LoadSource(pmb->phydro->u, IDN, NGHOST, four_pi_G_);
      if(mode_>=2) // iterative mode - load initial guess
        pmggrav->LoadFinestData(pmb->pgrav->phi, 0, NGHOST);
    }
//    else { // on another process
//    }
    pmggrav=pmggrav->next;
  }
  return;
}

