//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file multigrid_driver.cpp
//  \brief implementation of functions in class MultigridDriver

// C/C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cmath>

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
                                 int invar, ParameterInput *pin)
{
  pmy_mesh_=pm;
  pblock_=pmb;
  nvar_=invar;
  if(pblock_->block_size.nx1!=pblock_->block_size.nx2
  || pblock_->block_size.nx1!=pblock_->block_size.nx3) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MultigridDriver::MultigridDriver" << std::endl
        << "The Multigrid solver requires logically cubic MeshBlock." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }
  if(pblock_->block_size.nx2==1 || pblock_->block_size.nx3==1 ) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MultigridDriver::MultigridDriver" << std::endl
        << "Currently the Multigrid solver works only in 3D." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }
  if(pm->use_meshgen_fn_[X1DIR]==true || pm->use_meshgen_fn_[X2DIR]==true
  || pm->use_meshgen_fn_[X3DIR]==true) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MultigridDriver::MultigridDriver" << std::endl
        << "Non-uniform mesh spacing is not supported." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }
  Real dx=pblock_->pcoord->dx1f(0);
  if(dx!=pblock_->pcoord->dx2f(0) || dx!=pblock_->pcoord->dx3f(0)) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MultigridDriver::MultigridDriver" << std::endl
        << "The cell size must be cubic." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }
  // count multigrid levels
  nmblevel_=0;
  for(int l=0; l<20; l++) {
    if((1<<l) == pblock_->block_size.nx1) {
      nmblevel_=l+1;
      break;
    }
  }
  if(nmblevel_==0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MultigridDriver::MultigridDriver" << std::endl
        << "The MeshBlock size must be power of two." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
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
    return;
  }
  nrootlevel_=std::min(nrlx,std::min(nrly,nrlz));
  rootsrc_.NewAthenaArray(nvar_,pm->nrbx3,pm->nrbx2,pm->nrbx1);
  ntotallevel_=nrootlevel_+nmblevel_-1;
  current_level_=ntotallevel_-1;

  fperiodic_=false;
  for(int i=0; i<6; i++)
    MGBoundaryFunction_[i]=pm->MGBoundaryFunction_[i];
  if(MGBoundaryFunction_[INNER_X1]==MGPeriodicInnerX1
  && MGBoundaryFunction_[OUTER_X1]==MGPeriodicOuterX1
  && MGBoundaryFunction_[INNER_X2]==MGPeriodicInnerX2
  && MGBoundaryFunction_[OUTER_X2]==MGPeriodicOuterX2
  && MGBoundaryFunction_[INNER_X3]==MGPeriodicInnerX3
  && MGBoundaryFunction_[OUTER_X3]==MGPeriodicOuterX3)
    fperiodic_=true;

  mode_=0; // 0: FMG+V(1,1), 1: FMG+F(0,1), 2: V(1,1)

  mgroot_=NULL;

  // Setting up the MPI information
  // *** this part should be modified when dedicate processes are allocated ***
  pblock_=pmb;
  MeshBlock *pb=pblock_;
  nblocks_=0;
  while(pb!=NULL) {
    nblocks_++;
    pb=pb->next;
  }

  nranks_  = Globals::nranks;
  nslist_  = new int[nranks_];
  nblist_  = new int[nranks_];
  nvlist_  = new int[nranks_];
  nvslist_ = new int[nranks_];
#ifdef MPI_PARALLEL
  MPI_Comm_dup(MPI_COMM_WORLD, &MPI_COMM_MULTIGRID);
#endif
  for(int n=0; n<nranks_; n++) {
    nslist_[n]  = pm->nslist[n];
    nblist_[n]  = pm->nblist[n];
    nvslist_[n] = nslist_[n]*nvar_;
    nvlist_[n]  = nblist_[n]*nvar_;
  }
  rootbuf_=new Real[pm->nbtotal*nvar_];

  mgtlist_ = new MultigridTaskList(this);
}

// destructor

MultigridDriver::~MultigridDriver()
{
  delete mgroot_;
  delete [] nslist_;
  delete [] nblist_;
  delete [] nvslist_;
  delete [] nvlist_;
  delete [] rootbuf_;
  rootsrc_.DeleteAthenaArray();
  delete mgtlist_;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SetupMultigrid(void)
//  \brief initialize the source assuming that the source terms are already loaded

void MultigridDriver::SetupMultigrid(void)
{
  MeshBlock *pb;

  if(fperiodic_)
    SubtractAverage(0);
  if(mode_<=1) { // FMG
    pb=pblock_;
    while(pb!=NULL) {
      Multigrid *pmg=GetMultigridBlock(pb);
      pmg->RestrictFMGSource();
      pb=pb->next;
    }
    FillRootGridSource();
    mgroot_->RestrictFMGSource();
    current_level_=0;
  }
  else current_level_=ntotallevel_-1;
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SubtractAverage(int type)
//  \brief Calculate the global average and subtract it

void MultigridDriver::SubtractAverage(int type)
{
  MeshBlock *pb=pblock_;
  while(pb!=NULL) {
    Multigrid *pmg=GetMultigridBlock(pb);
    for(int v=0; v<nvar_; v++)
      rootbuf_[pb->gid*nvar_+v]=pmg->CalculateTotal(type, v);
    pb=pb->next;
  }
#ifdef MPI_PARALLEL
  MPI_Allgatherv(MPI_IN_PLACE, nblist_[Globals::my_rank]*nvar_, MPI_ATHENA_REAL,
                 rootbuf_, nvlist_, nvslist_, MPI_ATHENA_REAL, MPI_COMM_MULTIGRID);
#endif
  Real vol=(pmy_mesh_->mesh_size.x1max-pmy_mesh_->mesh_size.x1min)
          *(pmy_mesh_->mesh_size.x2max-pmy_mesh_->mesh_size.x2min)
          *(pmy_mesh_->mesh_size.x3max-pmy_mesh_->mesh_size.x3min);
  for(int v=0; v<nvar_; v++) {
    Real total=0.0;
    for(int n=0; n<Globals::nranks; n++)
      total+=rootbuf_[n*nvar_+v];
    Real ave=total/vol;
    pb=pblock_;
    while(pb!=NULL) {
      Multigrid *pmg=GetMultigridBlock(pb);
      pmg->SubtractAverage(type, v, ave);
      pb=pb->next;
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::FillRootGridSource(void)
//  \brief collect the coarsest data and fill the root grid

void MultigridDriver::FillRootGridSource(void)
{
  MeshBlock *pb=pblock_;
  while(pb!=NULL) {
    Multigrid *pmg=GetMultigridBlock(pb);
    for(int v=0; v<nvar_; v++)
      rootbuf_[pb->gid*nvar_+v]=pmg->GetRootSource(v);
    pb=pb->next;
  }
#ifdef MPI_PARALLEL
  MPI_Allgatherv(MPI_IN_PLACE, nblist_[Globals::my_rank]*nvar_, MPI_ATHENA_REAL,
                 rootbuf_, nvlist_, nvslist_, MPI_ATHENA_REAL, MPI_COMM_MULTIGRID);
#endif
  if(pmy_mesh_->multilevel) {
    // *** implement later
  }
  else { // uniform
    for(int n=0; n<pmy_mesh_->nbtotal; n++) {
      LogicalLocation &loc=pmy_mesh_->loclist[n];
      for(int v=0; v<nvar_; v++)
        rootsrc_(v,loc.lx3,loc.lx2,loc.lx1)=rootbuf_[n*nvar_+v];
    }
    mgroot_->LoadSource(rootsrc_,0,0,1.0);
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::FMGProlongate(void)
//  \brief Prolongation for FMG Cycle

void MultigridDriver::FMGProlongate(void)
{
  if(current_level_==nrootlevel_-1)
    TransferFromRootToBlocks();
  if(current_level_ >= nrootlevel_-1) {
    mgtlist_->SetMGTaskListFMGProlongate();
    mgtlist_->DoTaskListOneSubStep(this);
  }
  else { // root grid
    mgroot_->ApplyPhysicalBoundaries();
    mgroot_->FMGProlongate();
  }
  current_level_++;
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::TransferFromRootToBlocks(void)
//  \brief Transfer the data from the root grid to the coarsest level of each MeshBlock

void MultigridDriver::TransferFromRootToBlocks(void)
{
  MeshBlock *pb=pblock_;
  AthenaArray<Real> &src=mgroot_->GetCurrentData();
  if(pmy_mesh_->multilevel) {
    // *** implement later ***
  }
  else {
    while(pb!=NULL) {
      Multigrid *pmg=GetMultigridBlock(pb);
      LogicalLocation &loc=pb->loc;
      pmg->SetFromRootGrid(src, loc.lx1, loc.lx2, loc.lx3);
      pb=pb->next;
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::OneStepToFiner(int nsmooth)
//  \brief smoothing and restriction one level

void MultigridDriver::OneStepToFiner(int nsmooth)
{
  int ngh=mgroot_->ngh_;
  if(current_level_==nrootlevel_-1)
    TransferFromRootToBlocks();
  if(current_level_ >= nrootlevel_-1) {
    mgtlist_->SetMGTaskListToFiner(nsmooth, ngh);
    mgtlist_->DoTaskListOneSubStep(this);
  }
  else { // root grid
    mgroot_->ApplyPhysicalBoundaries();
    mgroot_->ProlongateAndCorrect();
    for(int n=0; n<nsmooth; n++) {
      mgroot_->ApplyPhysicalBoundaries();
      mgroot_->Smooth(0);
      mgroot_->ApplyPhysicalBoundaries();
      mgroot_->Smooth(1);
    }
  }
  current_level_++;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::OneStepToCoarser(int nsmooth)
//  \brief smoothing and restriction one level

void MultigridDriver::OneStepToCoarser(int nsmooth)
{
  int ngh=mgroot_->ngh_;
  if(current_level_ >= nrootlevel_) {
    mgtlist_->SetMGTaskListToCoarser(nsmooth, ngh);
    mgtlist_->DoTaskListOneSubStep(this);
    if(current_level_==nrootlevel_) {
      FillRootGridSource();
      mgroot_->ZeroClearData();
    }
  }
  else { // root grid
    for(int n=0; n<nsmooth; n++) {
      mgroot_->ApplyPhysicalBoundaries();
      mgroot_->Smooth(0);
      mgroot_->ApplyPhysicalBoundaries();
      mgroot_->Smooth(1);
    }
    mgroot_->ApplyPhysicalBoundaries();
    mgroot_->Restrict();
  }
  current_level_--;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SolveVCycle(int npresmooth, int npostsmooth)
//  \brief Solve the V-cycle starting from the current level

void MultigridDriver::SolveVCycle(int npresmooth, int npostsmooth)
{
  int startlevel=current_level_;
  while(current_level_>0)
    OneStepToCoarser(npostsmooth);
  SolveCoarsestGrid();
  while(current_level_<startlevel)
    OneStepToFiner(npresmooth);
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SolveFCycle(int npresmooth, int npostsmooth)
//  \brief Solve the F-cycle starting from the current level

void MultigridDriver::SolveFCycle(int npresmooth, int npostsmooth)
{
  // *** need to implement
  int startlevel=current_level_;
  while(current_level_>0)
    OneStepToCoarser(npostsmooth);
  SolveCoarsestGrid();
  while(current_level_<startlevel)
    OneStepToFiner(npresmooth);
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SolveFMGCycle(void)
//  \brief Solve the FMG Cycle using the V(1,1) or F(0,1) cycle

void MultigridDriver::SolveFMGCycle(void)
{
  for(int lev=0; lev<ntotallevel_; lev++) {
    if(mode_==0)
      SolveVCycle(1, 1);
    else if(mode_==1)
      SolveFCycle(0, 1);
    if(lev!=ntotallevel_-1) FMGProlongate();
  }
  if(fperiodic_)
    SubtractAverage(1);
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SolveCoarsestGrid(void)
//  \brief Solve the coarsest root grid

void MultigridDriver::SolveCoarsestGrid(void)
{
  Mesh *pm=pmy_mesh_;
  int ni=std::max(pm->nrbx1, std::max(pm->nrbx2, pm->nrbx3)) >> (nrootlevel_-1);
  if(fperiodic_ && ni==1) { // trivial case - all zero
    mgroot_->ZeroClearData();
    return;
  }
  else {
    for(int i=0; i<ni; ni++) { // iterate ni times
      mgroot_->ApplyPhysicalBoundaries();
      mgroot_->Smooth(0);
      mgroot_->ApplyPhysicalBoundaries();
      mgroot_->Smooth(1);
    }
    mgroot_->ApplyPhysicalBoundaries();
    if(fperiodic_) {
      Real vol=(pm->mesh_size.x1max-pm->mesh_size.x1min)
              *(pm->mesh_size.x2max-pm->mesh_size.x2min)
              *(pm->mesh_size.x3max-pm->mesh_size.x3min);
      for(int v=0; v<nvar_; v++) {
        Real ave=mgroot_->CalculateTotal(1, v)/vol;
        mgroot_->SubtractAverage(1, v, ave);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn Real MultigridDriver::CalculateDefectNorm(int n, int nrm)
//  \brief calculate the defect norm

Real MultigridDriver::CalculateDefectNorm(int n, int nrm)
{
  MeshBlock *pb=pblock_;
  Real norm=0.0;
  while(pb!=NULL) {
    Multigrid *pmg=GetMultigridBlock(pb);
    if(nrm==0)
      norm=std::max(norm, pmg->CalculateDefectNorm(n, nrm));
    else
      norm+=pmg->CalculateDefectNorm(n, nrm);
    pb=pb->next;
  }
#ifdef MPI_PARALLEL
  if(nrm==0)
    MPI_Allreduce(MPI_IN_PLACE,&norm,1,MPI_ATHENA_REAL,MPI_MAX,MPI_COMM_MULTIGRID);
  else
    MPI_Allreduce(MPI_IN_PLACE,&norm,1,MPI_ATHENA_REAL,MPI_SUM,MPI_COMM_MULTIGRID);
#endif
  if(nrm==2)
    norm=std::sqrt(norm);

  return norm;
}

