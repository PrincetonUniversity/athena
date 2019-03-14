//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mg_gravity.cpp
//  \brief create multigrid solver for gravity

// C headers

// C++ headers
#include <iostream>
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../multigrid/multigrid.hpp"
#include "../parameter_input.hpp"
#include "gravity.hpp"
#include "mg_gravity.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

class MeshBlock;

//----------------------------------------------------------------------------------------
//! \fn MGGravityDriver::MGGravityDriver(Mesh *pm, MGBoundaryFunc *MGBoundary,
//                                   ParameterInput *pin)
//  \brief MGGravityDriver constructor

MGGravityDriver::MGGravityDriver(Mesh *pm, MGBoundaryFunc *MGBoundary,
                                 ParameterInput *pin)
    : MultigridDriver(pm, MGBoundary, 1) {
  four_pi_G_=pmy_mesh_->four_pi_G_;
  eps_=pmy_mesh_->grav_eps_;
  if (four_pi_G_==0.0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MGGravityDriver::MGGravityDriver" << std::endl
        << "Gravitational constant must be set in the Mesh::InitUserMeshData "
        << "using the SetGravitationalConstant or SetFourPiG function." << std::endl;
    ATHENA_ERROR(msg);
  }
  if (mode_>=2 && eps_<0.0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MGGravityDriver::MGGravityDriver" << std::endl
        << "Convergence threshold must be set in the Mesh::InitUserMeshData "
        << "using the SetGravitatyThreshold for the iterative mode." << std::endl
        << "Set the threshold = 0.0 for automatic convergence control." << std::endl;
    ATHENA_ERROR(msg);
  }

  // Allocate multigrid objects
  RegionSize root_size=pmy_mesh_->mesh_size;
  root_size.nx1 = pmy_mesh_->nrbx1;
  root_size.nx2 = pmy_mesh_->nrbx2;
  root_size.nx3 = pmy_mesh_->nrbx3;
  LogicalLocation lroot;
  lroot.lx1=0, lroot.lx2=0, lroot.lx3=0, lroot.level=0;
  mgroot_= new MGGravity(this,lroot,-1,-1,root_size,MGBoundary,pmy_mesh_->mesh_bcs,true);
  pmg_=nullptr;
  // Multigrid *pfirst;
  int nbs=nslist_[Globals::my_rank];
  int nbe=nbs+nblist_[Globals::my_rank]-1;
  RegionSize block_size;
  block_size.nx1 = pmy_mesh_->mesh_size.nx1/pmy_mesh_->nrbx1;
  block_size.nx2 = pmy_mesh_->mesh_size.nx2/pmy_mesh_->nrbx2;
  block_size.nx3 = pmy_mesh_->mesh_size.nx3/pmy_mesh_->nrbx3;
  for (int i=nbs; i<=nbe; i++) {
    BoundaryFlag block_bcs[6];
    pmy_mesh_->SetBlockSizeAndBoundaries(pmy_mesh_->loclist[i], block_size, block_bcs);
    Multigrid *nmg=new MGGravity(this, pmy_mesh_->loclist[i], i, i-nbs, block_size,
                                 MGBoundary, block_bcs, false);
    nmg->pmgbval->SearchAndSetNeighbors(pmy_mesh_->tree, ranklist_, nslist_);
    AddMultigrid(nmg);
  }
}


//----------------------------------------------------------------------------------------
//! \fn void MGGravityDriver::Solve(int stage)
//  \brief load the data and solve

void MGGravityDriver::Solve(int stage) {
  Multigrid *pmggrav=pmg_;
  AthenaArray<Real> in;

  // Load the source
  while (pmggrav!=nullptr) {
    MeshBlock *pmb=pmy_mesh_->FindMeshBlock(pmggrav->gid_);
    if (pmb!=nullptr) {
      in.InitWithShallowCopy(pmb->phydro->u);
      pmggrav->LoadSource(in, IDN, NGHOST, four_pi_G_);
      if (mode_>=2) // iterative mode - load initial guess
        pmggrav->LoadFinestData(pmb->pgrav->phi, 0, NGHOST);
    }
    //    else { // on another process
    //    }
    pmggrav=pmggrav->next;
  }

  SetupMultigrid();
  Real mean_rho=0.0;
  if (fperiodic_)
    mean_rho=last_ave_/four_pi_G_;

  if (mode_<=1) {
    SolveFMGCycle();
  } else {
    SolveIterative();
  }

  // Return the result
  pmggrav=pmg_;
  while (pmggrav!=nullptr) {
    MeshBlock *pmb=pmy_mesh_->FindMeshBlock(pmggrav->gid_);
    if (pmb!=nullptr) {
      pmggrav->RetrieveResult(pmb->pgrav->phi,0,NGHOST);
      pmb->pgrav->grav_mean_rho=mean_rho;
    }
    //    else { // on another process
    //    }
    pmggrav=pmggrav->next;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn  void MGGravity::Smooth(int color)
//  \brief Red-Black Gauss-Seidel Smoother
void MGGravity::Smooth(int color) {
  int c=color;
  AthenaArray<Real> &u=u_[current_level_];
  AthenaArray<Real> &src=src_[current_level_];
  int ll=nlevel_-1-current_level_;
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+(size_.nx1>>ll)-1, je=js+(size_.nx2>>ll)-1, ke=ks+(size_.nx3>>ll)-1;
  Real dx = rdx_*static_cast<Real>(1<<ll);
  Real dx2 = SQR(dx);
  Real isix=omega_/6.0;
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is+c; i<=ie; i+=2)
        u(0,k,j,i)-=((6.0*u(0,k,j,i)-u(0,k+1,j,i)-u(0,k,j+1,i)-u(0,k,j,i+1)
                      -u(0,k-1,j,i)-u(0,k,j-1,i)-u(0,k,j,i-1))+src(0,k,j,i)*dx2)*isix;
      c^=1;
    }
    c^=1;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MGGravity::CalculateDefect()
//  \brief calculate the residual

void MGGravity::CalculateDefect() {
  AthenaArray<Real> &u=u_[current_level_];
  AthenaArray<Real> &src=src_[current_level_];
  AthenaArray<Real> &def=def_[current_level_];
  int ll=nlevel_-1-current_level_;
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+(size_.nx1>>ll)-1, je=js+(size_.nx2>>ll)-1, ke=ks+(size_.nx3>>ll)-1;
  Real dx = rdx_*static_cast<Real>(1<<ll);
  Real idx2 = 1.0/SQR(dx);
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        def(0,k,j,i)=(6.0*u(0,k,j,i)-u(0,k+1,j,i)-u(0,k,j+1,i)-u(0,k,j,i+1)
                      -u(0,k-1,j,i)-u(0,k,j-1,i)-u(0,k,j,i-1))*idx2+src(0,k,j,i);
      }
    }
  }
  return;
}
