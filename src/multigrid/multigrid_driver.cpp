//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file multigrid_driver.cpp
//  \brief implementation of functions in class MultigridDriver

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdlib>    // abs
#include <iostream>   // endl
#include <iomanip>    // setprecision
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/cc/mg/bvals_mg.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "multigrid.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// constructor, initializes data structures and parameters

MultigridDriver::MultigridDriver(Mesh *pm, MGBoundaryFunc *MGBoundary, int invar) :
    nvar_(invar), maxreflevel_(pm->max_level-pm->root_level),
    mode_(2), // 0: FMG+V(1,1), 1: FMG+F(0,1), 2: V(1,1)
    nrbx1_(pm->nrbx1), nrbx2_(pm->nrbx2), nrbx3_(pm->nrbx3), pmy_mesh_(pm),
    fsubtract_average_(false), ffas_(true), eps_(-1.0), cbuf_(nvar_,3,3,3) {
  if (pmy_mesh_->mesh_size.nx2==1 || pmy_mesh_->mesh_size.nx3==1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MultigridDriver::MultigridDriver" << std::endl
        << "Currently the Multigrid solver works only in 3D." << std::endl;
    ATHENA_ERROR(msg);
    return;
  }
  if ( !(pmy_mesh_->use_uniform_meshgen_fn_[X1DIR])
    || !(pmy_mesh_->use_uniform_meshgen_fn_[X2DIR])
    || !(pmy_mesh_->use_uniform_meshgen_fn_[X3DIR])) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MultigridDriver::MultigridDriver" << std::endl
        << "Non-uniform mesh spacing is not supported." << std::endl;
    ATHENA_ERROR(msg);
    return;
  }

  for (int i=0; i<6; i++)
    MGBoundaryFunction_[i]=MGBoundary[i];

  if ( (MGBoundaryFunction_[BoundaryFace::inner_x1] == MGPeriodicInnerX1
    ||  MGBoundaryFunction_[BoundaryFace::inner_x1] == MGZeroGradientInnerX1)
    && (MGBoundaryFunction_[BoundaryFace::outer_x1] == MGPeriodicOuterX1
    ||  MGBoundaryFunction_[BoundaryFace::outer_x1] == MGZeroGradientOuterX1)
    && (MGBoundaryFunction_[BoundaryFace::inner_x2] == MGPeriodicInnerX2
    ||  MGBoundaryFunction_[BoundaryFace::inner_x2] == MGZeroGradientInnerX2)
    && (MGBoundaryFunction_[BoundaryFace::outer_x2] == MGPeriodicOuterX2
    ||  MGBoundaryFunction_[BoundaryFace::outer_x2] == MGZeroGradientOuterX2)
    && (MGBoundaryFunction_[BoundaryFace::inner_x3] == MGPeriodicInnerX3
    ||  MGBoundaryFunction_[BoundaryFace::inner_x3] == MGZeroGradientInnerX3)
    && (MGBoundaryFunction_[BoundaryFace::outer_x3] == MGPeriodicOuterX3)
    ||  MGBoundaryFunction_[BoundaryFace::outer_x3] == MGZeroGradientOuterX3)
    fsubtract_average_ = true;

  // Setting up the MPI information
  // *** this part should be modified when dedicate processes are allocated ***
  // *** we also need to construct another neighbor list for Multigrid ***
  ranklist_  = new int[pmy_mesh_->nbtotal];
  for (int n=0; n<pmy_mesh_->nbtotal; ++n)
    ranklist_[n]=pmy_mesh_->ranklist[n];
  nranks_  = Globals::nranks;
  nslist_  = new int[nranks_];
  nblist_  = new int[nranks_];
  nvlist_  = new int[nranks_];
  nvslist_ = new int[nranks_];
  nvlisti_  = new int[nranks_];
  nvslisti_ = new int[nranks_];
#ifdef MPI_PARALLEL
  MPI_Comm_dup(MPI_COMM_WORLD, &MPI_COMM_MULTIGRID);
  mg_phys_id_ = pmy_mesh_->ReserveTagPhysIDs(1);
#endif
  int nv = nvar_;
  if (ffas_) nv*=2;
  // assume the same parallelization as hydro
  for (int n=0; n<nranks_; ++n) {
    nslist_[n]  = pmy_mesh_->nslist[n];
    nblist_[n]  = pmy_mesh_->nblist[n];
    nvslist_[n] = nslist_[n]*nv;
    nvlist_[n]  = nblist_[n]*nv;
    nvslisti_[n] = nslist_[n]*nvar_;
    nvlisti_[n]  = nblist_[n]*nvar_;
  }
  rootbuf_=new Real[pm->nbtotal*nv];
  mgtlist_ = new MultigridTaskList(this);

  if (maxreflevel_ > 0) { // SMR / AMR
    octets_ = new std::vector<MGOctet>[maxreflevel_];
    octetmap_ = new std::unordered_map<LogicalLocation, int,
                                       LogicalLocationHash>[maxreflevel_];
    noctets_ = new int[maxreflevel_]();
    prevnoct_ = new int[maxreflevel_];
  }
}

// destructor

MultigridDriver::~MultigridDriver() {
  delete [] ranklist_;
  delete [] nslist_;
  delete [] nblist_;
  delete [] nvlist_;
  delete [] nvslist_;
  delete [] nvlisti_;
  delete [] nvslisti_;
  delete [] rootbuf_;
  delete mgtlist_;
  if (maxreflevel_ > 0) {
    delete [] octets_;
    delete [] octetmap_;
    delete [] noctets_;
    delete [] prevnoct_;
  }
#ifdef MPI_PARALLEL
  MPI_Comm_free(&MPI_COMM_MULTIGRID);
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SetupMultigrid()
//  \brief initialize the source assuming that the source terms are already loaded

void MultigridDriver::SetupMultigrid() {
  nrootlevel_ = mgroot_->GetNumberOfLevels();
  nmblevel_ = vmg_[0]->GetNumberOfLevels();
  nreflevel_ = pmy_mesh_->current_level - pmy_mesh_->root_level;
  ntotallevel_ = nrootlevel_ + nmblevel_ + nreflevel_ - 1;
  fmglevel_ = current_level_ = ntotallevel_ - 1;
  int ncoct = mgroot_->ngh_*2 + 2;
  os_ = mgroot_->ngh_;
  oe_ = os_+1;

  // note: the level of an Octet is one level lower than the data stored there
  if (nreflevel_ > 0 && pmy_mesh_->amr_updated) {
    for (int l=0; l<maxreflevel_; ++l) { // clear old data
      octetmap_[l].clear();
      prevnoct_[l] = noctets_[l];
      noctets_[l] = 0;
    }
    pmy_mesh_->tree.CountMGOctets(noctets_);
    for (int l=0; l<maxreflevel_; ++l) { // increase the octet array size if needed
      if (prevnoct_[l] < noctets_[l]) {
        octets_[l].resize(noctets_[l]);
        octetmap_[l].reserve(noctets_[l]);
      }
      noctets_[l] = 0;
    }
    pmy_mesh_->tree.GetMGOctetList(octets_, octetmap_, noctets_);
    for (int l=0; l<maxreflevel_; ++l) {
      for (int o=prevnoct_[l]; o<noctets_[l]; ++o) {
        octets_[l][o].u.NewAthenaArray(nvar_, ncoct, ncoct, ncoct);
        octets_[l][o].def.NewAthenaArray(nvar_, ncoct, ncoct, ncoct);
        octets_[l][o].src.NewAthenaArray(nvar_, ncoct, ncoct, ncoct);
        if (ffas_)
          octets_[l][o].uold.NewAthenaArray(nvar_, ncoct, ncoct, ncoct);
      }
    }
  }

  for (auto itr = vmg_.begin(); itr<vmg_.end(); itr++) {
    Multigrid *pmg = *itr;
    pmg->pmgbval->CopyNeighborInfoFromMeshBlock();
  }

  if (fsubtract_average_)
    SubtractAverage(MGVariable::src);
  if (mode_<=1) { // FMG
    for (auto itr = vmg_.begin(); itr<vmg_.end(); itr++) {
      Multigrid *pmg = *itr;
      pmg->RestrictFMGSource();
    }
    TransferFromBlocksToRoot(true);
    RestrictFMGSourceOctets();
    mgroot_->RestrictFMGSource();
    current_level_=0;
  } else {
    current_level_=ntotallevel_-1;
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SubtractAverage(MGVariable type)
//  \brief Calculate the global average and subtract it

void MultigridDriver::SubtractAverage(MGVariable type) {
  // ******
  for (auto itr = vmg_.begin(); itr<vmg_.end(); itr++) {
    Multigrid *pmg = *itr;
    for (int v=0; v<nvar_; ++v)
      rootbuf_[pmg->pmy_block_->gid*nvar_+v] = pmg->CalculateTotal(type, v);
  }
#ifdef MPI_PARALLEL
  MPI_Allgatherv(MPI_IN_PLACE, nblist_[Globals::my_rank]*nvar_, MPI_ATHENA_REAL,
                 rootbuf_, nvlisti_, nvslisti_, MPI_ATHENA_REAL, MPI_COMM_MULTIGRID);
#endif
  Real vol = (pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min)
           * (pmy_mesh_->mesh_size.x2max - pmy_mesh_->mesh_size.x2min)
           * (pmy_mesh_->mesh_size.x3max - pmy_mesh_->mesh_size.x3min);
  for (int v=0; v<nvar_; ++v) {
    Real total = 0.0;
    for (int n=0; n<pmy_mesh_->nbtotal; ++n)
      total += rootbuf_[n*nvar_+v];
    last_ave_ = total/vol;
    for (auto itr = vmg_.begin(); itr<vmg_.end(); itr++) {
      Multigrid *pmg = *itr;
      pmg->SubtractAverage(type, v, last_ave_);
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::TransferFromBlocksToRoot(bool initflag)
//  \brief collect the coarsest data and transfer to the root grid

void MultigridDriver::TransferFromBlocksToRoot(bool initflag) {
  int nv = nvar_, ngh = mgroot_->ngh_;
  if (ffas_ && !initflag) nv*=2;
  for (auto itr = vmg_.begin(); itr<vmg_.end(); itr++) {
    Multigrid *pmg = *itr;
    for (int v=0; v<nvar_; ++v)
      rootbuf_[pmg->pmy_block_->gid*nv+v]=pmg->GetCoarsestData(MGVariable::src, v);
    if (ffas_ && !initflag) {
      for (int v=0; v<nvar_; ++v)
        rootbuf_[pmg->pmy_block_->gid*nv+nvar_+v]=pmg->GetCoarsestData(MGVariable::u, v);
    }
  }

#ifdef MPI_PARALLEL
  if (!initflag)
    MPI_Allgatherv(MPI_IN_PLACE, nblist_[Globals::my_rank]*nv, MPI_ATHENA_REAL,
                   rootbuf_, nvlist_, nvslist_, MPI_ATHENA_REAL, MPI_COMM_MULTIGRID);
  else
    MPI_Allgatherv(MPI_IN_PLACE, nblist_[Globals::my_rank]*nvar_, MPI_ATHENA_REAL,
                   rootbuf_, nvlisti_, nvslisti_, MPI_ATHENA_REAL, MPI_COMM_MULTIGRID);
#endif

  for (int n=0; n<pmy_mesh_->nbtotal; ++n) {
    const LogicalLocation &loc=pmy_mesh_->loclist[n];
    int i = static_cast<int>(loc.lx1);
    int j = static_cast<int>(loc.lx2);
    int k = static_cast<int>(loc.lx3);
    if (loc.level == nrootlevel_-1) {
      for (int v=0; v<nvar_; ++v)
        mgroot_->SetData(MGVariable::src, v, k, j, i, rootbuf_[n*nv+v]);
      if (ffas_ && !initflag) {
        for (int v=0; v<nvar_; ++v)
          mgroot_->SetData(MGVariable::u, v, k, j, i, rootbuf_[n*nv+nvar_+v]);
      }
    } else {
      LogicalLocation oloc;
      oloc.lx1 = (loc.lx1 >> 1);
      oloc.lx2 = (loc.lx2 >> 1);
      oloc.lx3 = (loc.lx3 >> 1);
      oloc.level = loc.level - 1;
      int olev = oloc.level - (nrootlevel_ - 1);
      int oid = octetmap_[olev][oloc];
      int oi = (i&1) + ngh;
      int oj = (j&1) + ngh;
      int ok = (k&1) + ngh;
      for (int v=0; v<nvar_; ++v)
        octets_[olev][oid].src(v,ok,oj,oi) = rootbuf_[n*nv+v];
      if (ffas_ && !initflag) {
        for (int v=0; v<nvar_; ++v)
          octets_[olev][oid].u(v,ok,oj,oi) = rootbuf_[n*nv+nvar_+v];
      }
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::TransferFromRootToBlocks()
//  \brief Transfer the data from the root grid to the coarsest level of each MeshBlock

void MultigridDriver::TransferFromRootToBlocks() {
  mgroot_->pmgbval->ApplyPhysicalBoundaries();
  for (auto itr = vmg_.begin(); itr<vmg_.end(); itr++) {
    Multigrid *pmg = *itr;
    pmg->SetFromRootGrid();
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::FMGProlongate()
//  \brief Prolongation for FMG Cycle

void MultigridDriver::FMGProlongate() {
  int flag=0;
  if (current_level_ == nrootlevel_ + nreflevel_ - 1) {
    TransferFromRootToBlocks();
    flag=1;
  }
  if (current_level_ >= nrootlevel_ + nreflevel_ - 1) { // MeshBlocks
    mgtlist_->SetMGTaskListFMGProlongate(flag);
    mgtlist_->DoTaskListOneStage(this);
  } else if (current_level_ >= nrootlevel_ - 1) { // root to octets
    if (current_level_ == nrootlevel_ - 1)
      mgroot_->pmgbval->ApplyPhysicalBoundaries();
    else
      SetBoundariesOctets(true, false);
    FMGProlongateOctets();
  } else { // root grid
    mgroot_->pmgbval->ApplyPhysicalBoundaries();
    mgroot_->FMGProlongateBlock();
  }
  current_level_++;

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::OneStepToFiner(int nsmooth)
//  \brief prolongation and smoothing one level

void MultigridDriver::OneStepToFiner(int nsmooth) {
  int ngh=mgroot_->ngh_;
  int flag=0;
  std::cout << "ToFiner " << current_level_ << std::endl;
  if (current_level_ == nrootlevel_ + nreflevel_ - 1) {
    TransferFromRootToBlocks();
    SetMeshBlockCoarsestBoundaries();
    flag=1;
  }
  if (current_level_ >= nrootlevel_ + nreflevel_ - 1) { // MeshBlocks
    if (current_level_ == ntotallevel_ - 2) flag=2;
    mgtlist_->SetMGTaskListToFiner(nsmooth, ngh, flag);
    mgtlist_->DoTaskListOneStage(this);
    current_level_++;
  } else if (current_level_ >= nrootlevel_ - 1) { // non uniform octets
    std::cout << "Prolongating octets " << current_level_ << " " << nrootlevel_-1 << std::endl;
    if (current_level_ == nrootlevel_ - 1) {
      mgroot_->pmgbval->ApplyPhysicalBoundaries();
    } else {
      std::cout << "Octets boundary for prolongation" << std::endl;
      SetBoundariesOctets(true, ffas_);
    }
    ProlongateAndCorrectOctets();
    current_level_++;
    for (int n=0; n<nsmooth; ++n) {
      SetBoundariesOctets(false, false);
      SmoothOctets(0);
      SetBoundariesOctets(false, false);
      SmoothOctets(1);
    }
  } else { // root grid
    mgroot_->pmgbval->ApplyPhysicalBoundaries();
    std::cout << "Prolongating root grid " << std::endl;
    mgroot_->ProlongateAndCorrectBlock();
    current_level_++;
    for (int n=0; n<nsmooth; ++n) {
      mgroot_->pmgbval->ApplyPhysicalBoundaries();
      mgroot_->SmoothBlock(0);
      mgroot_->pmgbval->ApplyPhysicalBoundaries();
      mgroot_->SmoothBlock(1);
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::OneStepToCoarser(int nsmooth)
//  \brief smoothing and restriction one level

void MultigridDriver::OneStepToCoarser(int nsmooth) {
  int ngh=mgroot_->ngh_;
  std::cout << "ToCoarser " << current_level_ << std::endl;
  if (current_level_ >= nrootlevel_ + nreflevel_) { // MeshBlocks
    mgtlist_->SetMGTaskListToCoarser(nsmooth, ngh);
    mgtlist_->DoTaskListOneStage(this);
    if (current_level_ == nrootlevel_ + nreflevel_) {
      TransferFromBlocksToRoot();
      if (!ffas_) {
        mgroot_->ZeroClearData();
        if (nreflevel_ > 0)
          ZeroClearAllOctets();
      }
    }
  } else if (current_level_ > nrootlevel_-1) { // refined octets
    SetBoundariesOctets(false, false);
    if (ffas_ && current_level_ < fmglevel_) {
      StoreOldDataOctets();
      CalculateFASRHSOctets();
    }
    for (int n=0; n<nsmooth; ++n) {
      SmoothOctets(0);
      SetBoundariesOctets(false, false);
      SmoothOctets(1);
      SetBoundariesOctets(false, false);
    }
    RestrictOctets();
  } else { // uniform root grid
    mgroot_->pmgbval->ApplyPhysicalBoundaries();
    if (ffas_ && current_level_ < fmglevel_) {
      mgroot_->StoreOldData();
      mgroot_->CalculateFASRHSBlock();
    }
    for (int n=0; n<nsmooth; ++n) {
      mgroot_->SmoothBlock(0);
      mgroot_->pmgbval->ApplyPhysicalBoundaries();
      mgroot_->SmoothBlock(1);
      mgroot_->pmgbval->ApplyPhysicalBoundaries();
    }
    mgroot_->RestrictBlock();
  }

  current_level_--;

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SolveVCycle(int npresmooth, int npostsmooth)
//  \brief Solve the V-cycle starting from the current level

void MultigridDriver::SolveVCycle(int npresmooth, int npostsmooth) {
  int startlevel=current_level_;
  while (current_level_>0)
    OneStepToCoarser(npresmooth);
  SolveCoarsestGrid();
  while (current_level_<startlevel)
    OneStepToFiner(npostsmooth);
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SolveFCycle(int npresmooth, int npostsmooth)
//  \brief Solve the F-cycle starting from the current level

void MultigridDriver::SolveFCycle(int npresmooth, int npostsmooth) {
  int startlevel=current_level_;
  int turnlevel;
  if (startlevel == 0)
    turnlevel=0;
  else
    turnlevel=1;
  for (; turnlevel<=startlevel; turnlevel++) {
    while (current_level_>0)
      OneStepToCoarser(npresmooth);
    SolveCoarsestGrid();
    while (current_level_<turnlevel)
      OneStepToFiner(npostsmooth);
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SolveFMGCycle()
//  \brief Solve the FMG Cycle using the V(1,1) or F(0,1) cycle

void MultigridDriver::SolveFMGCycle() {
  for (fmglevel_=0; fmglevel_<ntotallevel_; fmglevel_++) {
    if (mode_ == 0)
      SolveVCycle(1, 1);
    else if (mode_ == 1)
      SolveFCycle(0, 1);
    if (fmglevel_!=ntotallevel_-1)
      FMGProlongate();
  }
  if (fsubtract_average_)
    SubtractAverage(MGVariable::u);
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SolveIterative()
//  \brief Solve iteratively until the convergence is achieved

void MultigridDriver::SolveIterative() {
  Real def=eps_+1e-10;
  int niter=0;
  std::cout << std::scientific << std::setprecision(14);
  while (def > eps_) {
    std::cout << "iter " << niter << std::endl;
    SolveVCycle(1,1);
    Real olddef=def;
    def=0.0;
    for (int v=0; v<nvar_; ++v)
      def+=CalculateDefectNorm(MGNormType::l2, v);
    std::cout << "after " << niter << " def = " << def << std::endl;
    if (niter > 0 && def/olddef > 0.9) {
      if (eps_ == 0.0) break;
      if (Globals::my_rank == 0)
        std::cout << "### Warning in MultigridDriver::SolveIterative" << std::endl
                  << "Slow multigrid convergence : defect norm = " << def
                  << ", convergence factor = " << def/olddef << "." << std::endl;
    }
    if (niter>100) {
      if (Globals::my_rank == 0) {
        std::cout
            << "### Warning in MultigridDriver::SolveIterative" << std::endl
            << "Aborting because the # iterations is too large, niter > 100." << std::endl
            << "Check the solution as it may not be accurate enough." << std::endl;
      }
      break;
    }
    niter++;
  }
  if (fsubtract_average_)
    SubtractAverage(MGVariable::u);
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SolveCoarsestGrid()
//  \brief Solve the coarsest root grid

void MultigridDriver::SolveCoarsestGrid() {
  int ni = (std::max(nrbx1_, std::max(nrbx2_, nrbx3_))
            >> (nrootlevel_-1));
  if (fsubtract_average_ && ni == 1) { // trivial case - all zero
    if (ffas_) { 
      mgroot_->pmgbval->ApplyPhysicalBoundaries();
      mgroot_->StoreOldData();
    }
    mgroot_->ZeroClearData();
    return;
  } else {
    if (fsubtract_average_) {
      Real vol=(mgroot_->size_.x1max-mgroot_->size_.x1min)
              *(mgroot_->size_.x2max-mgroot_->size_.x2min)
              *(mgroot_->size_.x3max-mgroot_->size_.x3min);
      for (int v=0; v<nvar_; ++v) {
        Real ave=mgroot_->CalculateTotal(MGVariable::u, v)/vol;
        mgroot_->SubtractAverage(MGVariable::u, v, ave);
      }
    }
    mgroot_->pmgbval->ApplyPhysicalBoundaries();
    if (ffas_) {
      mgroot_->StoreOldData();
      mgroot_->CalculateFASRHSBlock();
    }
    for (int i=0; i<ni; ++i) { // iterate ni times
      mgroot_->SmoothBlock(0);
      mgroot_->pmgbval->ApplyPhysicalBoundaries();
      mgroot_->SmoothBlock(1);
      mgroot_->pmgbval->ApplyPhysicalBoundaries();
    }
    if (fsubtract_average_) {
      Real vol=(mgroot_->size_.x1max-mgroot_->size_.x1min)
              *(mgroot_->size_.x2max-mgroot_->size_.x2min)
              *(mgroot_->size_.x3max-mgroot_->size_.x3min);
      for (int v=0; v<nvar_; ++v) {
        Real ave=mgroot_->CalculateTotal(MGVariable::u, v)/vol;
        mgroot_->SubtractAverage(MGVariable::u, v, ave);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn Real MultigridDriver::CalculateDefectNorm(MGNormType nrm, int n)
//  \brief calculate the defect norm

Real MultigridDriver::CalculateDefectNorm(MGNormType nrm, int n) {
  Real norm=0.0;
  for (auto itr = vmg_.begin(); itr < vmg_.end(); itr++) {
    Multigrid *pmg = *itr;
    if (nrm == MGNormType::max)
      norm=std::max(norm, pmg->CalculateDefectNorm(nrm, n));
    else
      norm+=pmg->CalculateDefectNorm(nrm, n);
  }
#ifdef MPI_PARALLEL
  if (nrm == MGNormType::max)
    MPI_Allreduce(MPI_IN_PLACE,&norm,1,MPI_ATHENA_REAL,MPI_MAX,MPI_COMM_MULTIGRID);
  else
    MPI_Allreduce(MPI_IN_PLACE,&norm,1,MPI_ATHENA_REAL,MPI_SUM,MPI_COMM_MULTIGRID);
#endif
  if (nrm != MGNormType::max) {
    Real vol=(mgroot_->size_.x1max-mgroot_->size_.x1min)
            *(mgroot_->size_.x2max-mgroot_->size_.x2min)
            *(mgroot_->size_.x3max-mgroot_->size_.x3min);
    norm /= vol;
  }
  if (nrm == MGNormType::l2)
    norm=std::sqrt(norm);

  return norm;
}


//----------------------------------------------------------------------------------------
//! \fn Multigrid* MultigridDriver::FindMultigrid(int tgid)
//  \brief return the Multigrid whose gid is tgid

Multigrid* MultigridDriver::FindMultigrid(int tgid) {
  for (auto itr = vmg_.begin(); itr < vmg_.end(); itr++) {
    Multigrid *pmg = *itr;
    if (pmg->pmy_block_->gid == tgid)
      return pmg;
  }

  return nullptr;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::RestrictFMGSourceOctets()
//  \brief restrict the source in octets for FMG

void MultigridDriver::RestrictFMGSourceOctets() {
  if (maxreflevel_ > 0) {
    const int &ngh = mgroot_->ngh_;
    for (int l=nreflevel_-1; l>=1; --l) {  // fine octets to coarse octets
      for (int o=0; o<noctets_[l]; ++o) {
        const LogicalLocation &loc = octets_[l][o].loc;
        LogicalLocation cloc;
        cloc.lx1 = (loc.lx1 >> 1);
        cloc.lx2 = (loc.lx2 >> 1);
        cloc.lx3 = (loc.lx3 >> 1);
        cloc.level = loc.level - 1;
        int oid = octetmap_[l-1][cloc];
        int oi = (static_cast<int>(loc.lx1) & 1) + ngh;
        int oj = (static_cast<int>(loc.lx2) & 1) + ngh;
        int ok = (static_cast<int>(loc.lx3) & 1) + ngh;
        mgroot_->Restrict(cbuf_, octets_[l][o].src, os_, oe_, os_, oe_, os_, oe_);
        for (int v=0; v<nvar_; ++v)
          octets_[l-1][oid].src(v, ok, oj, oi) = cbuf_(v, ngh, ngh, ngh);
      }
    }
    for (int o=0; o<noctets_[0]; ++o) { // octets to the root grid
      const LogicalLocation &loc = octets_[0][o].loc;
      int ri = (loc.lx1 >> 1);
      int rj = (loc.lx2 >> 1);
      int rk = (loc.lx3 >> 1);
      mgroot_->Restrict(cbuf_, octets_[0][o].src, os_, oe_, os_, oe_, os_, oe_);
      for (int v=0; v<nvar_; ++v)
        mgroot_->SetData(MGVariable::src, v, rk, rj, ri, cbuf_(v, ngh, ngh, ngh));
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::RestrictOctets()
//  \brief restrict the potential in octets

void MultigridDriver::RestrictOctets() {
  const int &ngh = mgroot_->ngh_;
  int lev = current_level_ - nrootlevel_;

  if (lev >= 1) { // fine octets to coarse octets
    for (int o=0; o<noctets_[lev]; ++o) {
      const LogicalLocation &loc = octets_[lev][o].loc;
      LogicalLocation cloc;
      cloc.lx1 = (loc.lx1 >> 1);
      cloc.lx2 = (loc.lx2 >> 1);
      cloc.lx3 = (loc.lx3 >> 1);
      cloc.level = loc.level - 1;
      int oid = octetmap_[lev-1][cloc];
      int oi = (static_cast<int>(loc.lx1) & 1) + ngh;
      int oj = (static_cast<int>(loc.lx2) & 1) + ngh;
      int ok = (static_cast<int>(loc.lx3) & 1) + ngh;
      mgroot_->CalculateDefect(octets_[lev][o].def, octets_[lev][o].u,
                               octets_[lev][o].src, lev+1, os_, oe_, os_, oe_, os_, oe_);
      mgroot_->Restrict(cbuf_, octets_[lev][o].def, os_, oe_, os_, oe_, os_, oe_);
      for (int v=0; v<nvar_; ++v)
        octets_[lev-1][oid].src(v, ok, oj, oi) = cbuf_(v, ngh, ngh, ngh);
      if (ffas_) {
        mgroot_->Restrict(cbuf_, octets_[lev][o].u, os_, oe_, os_, oe_, os_, oe_);
        for (int v=0; v<nvar_; ++v)
          octets_[lev-1][oid].u(v, ok, oj, oi) = cbuf_(v, ngh, ngh, ngh);
      }
    }
  } else { // octets to the root grid
    for (int o=0; o<noctets_[0]; ++o) {
      const LogicalLocation &loc = octets_[0][o].loc;
      int ri = loc.lx1; //(loc.lx1 >> 1);
      int rj = loc.lx2; //(loc.lx2 >> 1);
      int rk = loc.lx3; //(loc.lx3 >> 1);
      mgroot_->CalculateDefect(octets_[0][o].def, octets_[0][o].u,
                               octets_[0][o].src, 1, os_, oe_, os_, oe_, os_, oe_);
      mgroot_->Restrict(cbuf_, octets_[0][o].def, os_, oe_, os_, oe_, os_, oe_);
      for (int v=0; v<nvar_; ++v)
        mgroot_->SetData(MGVariable::src, v, rk, rj, ri, cbuf_(v, ngh, ngh, ngh));
      if (ffas_) {
        mgroot_->Restrict(cbuf_, octets_[0][o].u, os_, oe_, os_, oe_, os_, oe_);
        for (int v=0; v<nvar_; ++v)
          mgroot_->SetData(MGVariable::u, v, rk, rj, ri, cbuf_(v, ngh, ngh, ngh));
      }
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::ZeroClearAllOctets()
//  \brief zero clear the data in all the octets

void MultigridDriver::ZeroClearAllOctets() {
  // *** clear only where needed in case of FMG
  int maxlevel = std::min(fmglevel_ - nrootlevel_, nreflevel_);
  
  for (int l=0; l<nreflevel_; l++) {
    for (int o=0; o<noctets_[l]; ++o)
      octets_[l][o].u.ZeroClear();
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::StoreOldDataOctets()
//  \brief store the old u data in the uold array in octets

void MultigridDriver::StoreOldDataOctets() {
  int lev = current_level_ - nrootlevel_;

  for (int o=0; o<noctets_[lev]; ++o)
    memcpy(octets_[lev][o].uold.data(), octets_[lev][o].u.data(),
           octets_[lev][o].u.GetSizeInBytes());

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::CalculateFASRHSOctets()
//  \brief Calculate the RHS for FAS in Octets
void MultigridDriver::CalculateFASRHSOctets() {
  int lev = current_level_ - nrootlevel_;

  for (int o=0; o<noctets_[lev]; ++o)
    mgroot_->CalculateFASRHS(octets_[lev][o].src, octets_[lev][o].u,
                             lev+1, os_, oe_, os_, oe_, os_, oe_);

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SmoothOctets(int color)
//  \brief Apply the smoothing operator on octets
void MultigridDriver::SmoothOctets(int color) {
  int lev = current_level_ - nrootlevel_;

  for (int o=0; o<noctets_[lev]; ++o) {
    std::cout << "SmoothOctet before " << o  << " color " << color<< std::endl;
    if (o==0) 
    for (int k=0; k<=3; k++) 
    for (int j=0; j<=3; j++) 
    for (int i=0; i<=3; i++) 
    std::cout << k << " " << j << " " << i << " " << octets_[lev][o].u(0,k,j,i) << " " << octets_[lev][o].src(0,k,j,i) << std::endl;
    mgroot_->Smooth(octets_[lev][o].u, octets_[lev][o].src,
                    lev+1, os_, oe_, os_, oe_, os_, oe_, color);
    std::cout << "SmoothOctet after  " << o << " " << octets_[lev][o].u(0,1,1,1) << " " << octets_[lev][o].u(0,2,2,2) << std::endl;
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::ProlongateAndCorrectOctets()
//  \brief Prolongate and correct the potential in octets

void MultigridDriver::ProlongateAndCorrectOctets() {
  int clev = current_level_ - nrootlevel_;
  int flev = clev + 1;
  int ngh = mgroot_->ngh_;
  bool faceonly = false;

  if (flev == 0) {  // from root to octets
    const AthenaArray<Real> &u = mgroot_->GetCurrentData();
    const AthenaArray<Real> &uold = mgroot_->GetCurrentOldData();
    for (int o=0; o<noctets_[0]; ++o) { // octets to the root grid
      const LogicalLocation &loc = octets_[0][o].loc;
      int ri = (static_cast<int>(loc.lx1) >> 1) + ngh - 1;
      int rj = (static_cast<int>(loc.lx2) >> 1) + ngh - 1;
      int rk = (static_cast<int>(loc.lx3) >> 1) + ngh - 1;
      if (ffas_) {
        std::cout << "prolongating " << o << std::endl;
        for (int v=0; v<nvar_; ++v) {
          for (int k=0; k<=2; ++k) {
            for (int j=0; j<=2; ++j) {
              for (int i=0; i<=2; ++i) {
                cbuf_(v,k,j,i) = u(v, rk+k, rj+j, ri+i) - uold(v, rk+k, rj+j, ri+i);
                std::cout << k << " " << j << "  "<< i << " "  <<cbuf_(0,k,j,i) << std::endl;
              }
            }
          }
        }
        mgroot_->ProlongateAndCorrect(octets_[0][o].u, cbuf_,
                                      ngh, ngh, ngh, ngh, ngh, ngh, ngh, ngh, ngh);
        std::cout << "prolongated " << o << std::endl;
        for (int k=1; k<=2; ++k)
          for (int j=1; j<=2; ++j)
            for (int i=1; i<=2; ++i)
              std::cout <<k << " " << j << " " << i << " " << octets_[0][o].u(0,k,j,i)<< std::endl; 
      } else {
        mgroot_->ProlongateAndCorrect(octets_[0][o].u, u,
                                      ri, ri, rj, rj, rk, rk, ngh, ngh, ngh);
      }
    }
  } else { // from coarse octets to fine octets
    int pcid = -1;
    for (int o=0; o<noctets_[flev]; ++o) {
      const LogicalLocation &loc = octets_[flev][o].loc;
      LogicalLocation cloc;
      cloc.lx1 = (loc.lx1 >> 1);
      cloc.lx2 = (loc.lx2 >> 1);
      cloc.lx3 = (loc.lx3 >> 1);
      cloc.level = loc.level - 1;
      int cid = octetmap_[clev][cloc];
      int ci = (static_cast<int>(loc.lx1) & 1) + ngh - 1;
      int cj = (static_cast<int>(loc.lx2) & 1) + ngh - 1;
      int ck = (static_cast<int>(loc.lx3) & 1) + ngh - 1;
      const AthenaArray<Real> &uc = octets_[clev][cid].u;
      const AthenaArray<Real> &ucold = octets_[clev][cid].uold;
      if (ffas_) {
        for (int v=0; v<nvar_; ++v) {
          for (int k=0; k<=2; ++k) {
            for (int j=0; j<=2; ++j) {
              for (int i=0; i<=2; ++i)
                cbuf_(v,k,j,i) = uc(v, ck+k, cj+j, ci+i) - ucold(v, ck+k, cj+j, ci+i);
            }
          }
        }
        mgroot_->ProlongateAndCorrect(octets_[flev][o].u, cbuf_,
                                      ngh, ngh, ngh, ngh, ngh, ngh, ngh, ngh, ngh);
      } else {
        mgroot_->ProlongateAndCorrect(octets_[flev][o].u, uc,
                 ci+ngh, ci+ngh, cj+ngh, cj+ngh, ck+ngh, ck+ngh, ngh, ngh, ngh);
      }
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::FMGProlongateOctets()
//  \brief Prolongate the potential in octets for FMG

void MultigridDriver::FMGProlongateOctets() {
  int clev = current_level_ - nrootlevel_;
  int flev = clev + 1;
  int ngh = mgroot_->ngh_;

  if (flev == 0) {  // from root to octets
    for (int o=0; o<noctets_[0]; ++o) { // octets to the root grid
      const LogicalLocation &loc = octets_[0][o].loc;
      int ri = (loc.lx1 >> 1) + ngh;
      int rj = (loc.lx2 >> 1) + ngh;
      int rk = (loc.lx3 >> 1) + ngh;
      mgroot_->FMGProlongate(octets_[0][o].u, mgroot_->GetCurrentData(),
                             ri, ri, rj, rj, rk, rk, ngh, ngh, ngh);
    }
  } else { // from coarse octets to fine octets
    for (int o=0; o<noctets_[flev]; ++o) {
      const LogicalLocation &loc = octets_[flev][o].loc;
      LogicalLocation cloc;
      cloc.lx1 = (loc.lx1 >> 1);
      cloc.lx2 = (loc.lx2 >> 1);
      cloc.lx3 = (loc.lx3 >> 1);
      cloc.level = loc.level - 1;
      int cid = octetmap_[clev][cloc];
      int ci = (static_cast<int>(loc.lx1) & 1) + ngh;
      int cj = (static_cast<int>(loc.lx2) & 1) + ngh;
      int ck = (static_cast<int>(loc.lx3) & 1) + ngh;
      mgroot_->FMGProlongate(octets_[flev][o].u, octets_[clev][cid].u,
                             ci, ci, cj, cj, ck, ck, ngh, ngh, ngh);
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SetBoundariesOctets(bool fprolong, bool folddata)
//  \brief Apply boundary conditions for octets

void MultigridDriver::SetBoundariesOctets(bool fprolong, bool folddata) {
  int lev = current_level_ - nrootlevel_;
  bool faceonly = false;
  if (!fprolong && mgroot_->btypef == BoundaryQuantity::mggrav_f)
    faceonly = true;

  for (int o=0; o<noctets_[lev]; ++o) {
    if (fprolong && octets_[lev][o].fleaf == true) continue;
    const LogicalLocation &loc = octets_[lev][o].loc;
    AthenaArray<Real> &u = octets_[lev][o].u;
    AthenaArray<Real> &uold = octets_[lev][o].uold;
    LogicalLocation nloc = loc;
    for (int ox3=-1; ox3<=1; ++ox3) {
      nloc.lx3 = loc.lx3 + ox3;
      if (nloc.lx3 < 0) {
        if (mgroot_->pmgbval->MGBoundaryFunction_[BoundaryFace::inner_x3]
          == MGPeriodicInnerX3) nloc.lx3 = (nrbx3_ << lev) - 1;
        else continue;
      }
      if (nloc.lx3 >= (nrbx3_ << lev)) {
        if (mgroot_->pmgbval->MGBoundaryFunction_[BoundaryFace::outer_x3]
          == MGPeriodicOuterX3) nloc.lx3 = 0;
        else continue;
      }
      for (int ox2=-1; ox2<=1; ++ox2) {
        nloc.lx2 = loc.lx2 + ox2;
        if (nloc.lx2 < 0) {
          if (mgroot_->pmgbval->MGBoundaryFunction_[BoundaryFace::inner_x2]
            == MGPeriodicInnerX2) nloc.lx2 = (nrbx2_ << lev) - 1;
          else continue;
        }
        if (nloc.lx2 >= (nrbx2_ << lev)) {
          if (mgroot_->pmgbval->MGBoundaryFunction_[BoundaryFace::outer_x2]
            == MGPeriodicOuterX2) nloc.lx2 = 0;
          else continue;
        }
        for (int ox1=-1; ox1<=1; ++ox1) {
          int nt = std::abs(ox1)+std::abs(ox2)+std::abs(ox3);
          if (faceonly &&  nt > 1) continue;
          if (ox1 == 0 && ox2 == 0 && ox3 == 0) continue;
          // find a neighboring octet - either on the same or coarser level
          nloc.lx1 = loc.lx1 + ox1;
          if (nloc.lx1 < 0) {
            if (mgroot_->pmgbval->MGBoundaryFunction_[BoundaryFace::inner_x1]
              == MGPeriodicInnerX1) nloc.lx1 = (nrbx1_ << lev) - 1;
            else continue;
          }
          if (nloc.lx1 >= (nrbx1_ << lev)) {
            if (mgroot_->pmgbval->MGBoundaryFunction_[BoundaryFace::outer_x1]
              == MGPeriodicOuterX1) nloc.lx1 = 0;
            else continue;
          }
          std::cout << "octet boundary " << loc.lx1 << " " << loc.lx2 << " " << loc.lx3 << " " << loc.level << " " << lev << " " << ox1 << " " << ox2 << " " << ox3 << std::endl;
          std::cout << "neighbor " << nloc.lx1 << " " << nloc.lx2 << " " << nloc.lx3 << " " << nloc.level << std::endl;
          if (octetmap_[lev].count(nloc) == 1) { // on the same level
            std::cout << "neighbor on the same level" << std::endl;
            const AthenaArray<Real> &un = octets_[lev][octetmap_[lev][nloc]].u;
            const AthenaArray<Real> &unold = octets_[lev][octetmap_[lev][nloc]].uold;
            SetOctetBoundarySameLevel(u, un, uold, unold, ox1, ox2, ox3, folddata);
          } else { // on the coarser level
            if (lev > 0) {
              std::cout << "neighbor on a coarser octet" << std::endl;
              LogicalLocation cloc;
              cloc.lx1 = nloc.lx1 >> 1;
              cloc.lx2 = nloc.lx2 >> 1;
              cloc.lx3 = nloc.lx3 >> 1;
              cloc.level = nloc.level - 1;
              std::cout << "coarse neighbor " << cloc.lx1 << " " << cloc.lx2 << " " << cloc.lx3 << " " << cloc.level << std::endl;
              const AthenaArray<Real> &un = octets_[lev-1][octetmap_[lev-1][cloc]].u;
              SetOctetBoundaryFromCoarser(u, un, loc, ox1, ox2, ox3);
            } else {
              std::cout << "neighbor on the root grid" << std::endl;
              SetOctetBoundaryFromCoarser(u, mgroot_->GetCurrentData(),
                                          loc, ox1, ox2, ox3);
            }
          }
        }
      }
    }
    ApplyPhysicalBoundariesOctet(u, loc);
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SetOctetBoundarySameLevel(AthenaArray<Real> &dst,
//   const AthenaArray<Real> &un, AthenaArray<Real> &uold, const AthenaArray<Real> &unold,
//   int ox1, int ox2, int ox3, bool folddata)
//  \brief set an Octet boundary from a neighbor Octet on the same level

void MultigridDriver::SetOctetBoundarySameLevel(AthenaArray<Real> &dst,
     const AthenaArray<Real> &un, AthenaArray<Real> &uold, const AthenaArray<Real> &unold,
     int ox1, int ox2, int ox3, bool folddata) {
  int ngh = mgroot_->ngh_;
  int is, ie, js, je, ks, ke, nis, njs, nks;
  if (ox1 == 0)     is = ngh,   ie = ngh+1, nis = ngh;
  else if (ox1 < 0) is = 0,     ie = ngh-1, nis = ngh+1;
  else              is = ngh+2, ie = ngh+2, nis = ngh; 
  if (ox2 == 0)     js = ngh,   je = ngh+1, njs = ngh;
  else if (ox2 < 0) js = 0,     je = ngh-1, njs = ngh+1;
  else              js = ngh+2, je = ngh+2, njs = ngh; 
  if (ox3 == 0)     ks = ngh,   ke = ngh+1, nks = ngh;
  else if (ox3 < 0) ks = 0,     ke = ngh-1, nks = ngh+1;
  else              ks = ngh+2, ke = ngh+2, nks = ngh; 
  for (int v=0; v<nvar_; ++v) {
    for (int k=ks, nk=nks; k<=ke; ++k, ++nk) {
      for (int j=js, nj=njs; j<=je; ++j, ++nj) {
        for (int i=is, ni=nis; i<=ie; ++i, ++ni)
          dst(v, k, j, i) = un(v, nk, nj, ni);
      }
    }
  }
  if (folddata) {
    for (int v=0; v<nvar_; ++v) {
      for (int k=ks, nk=nks; k<=ke; ++k, ++nk) {
        for (int j=js, nj=njs; j<=je; ++j, ++nj) {
          for (int i=is, ni=nis; i<=ie; ++i, ++ni)
            uold(v, k, j, i) = unold(v, nk, nj, ni);
        }
      }
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::ApplyPhysicalBoundariesOctet(AthenaArray<Real> &u,
//                                                         const LogicalLocation &loc)
//  \brief Apply physical boundary conditions for an octet

void MultigridDriver::ApplyPhysicalBoundariesOctet(AthenaArray<Real> &u,
                                                   const LogicalLocation &loc) {
  int lev = current_level_ - nrootlevel_;
  int ngh = mgroot_->ngh_;
  Real time = pmy_mesh_->time;
  Real fac = 1.0/(static_cast<Real>(2<<lev));
  Real dx = mgroot_->rdx_*fac;
  Real dy = mgroot_->rdy_*fac;
  Real dz = mgroot_->rdz_*fac;
  Real x0 = mgroot_->size_.x1min + (static_cast<Real>(loc.lx1*2 - ngh) - 0.5)*dx;
  Real y0 = mgroot_->size_.x2min + (static_cast<Real>(loc.lx2*2 - ngh) - 0.5)*dy;
  Real z0 = mgroot_->size_.x3min + (static_cast<Real>(loc.lx2*2 - ngh) - 0.5)*dz;
  if (loc.lx1 == 0
    && mgroot_->pmgbval->MGBoundaryFunction_[BoundaryFace::inner_x1]
                      != MGPeriodicInnerX1
    && mgroot_->pmgbval->MGBoundaryFunction_[BoundaryFace::inner_x1] != nullptr)
      MGBoundaryFunction_[BoundaryFace::inner_x1](u, time, nvar_,
                0,     0,     ngh,   ngh+1, ngh,   ngh+1, ngh, x0, y0, z0, dx, dy, dz);
  if (loc.lx1 == (nrbx1_<<lev)-1
    && mgroot_->pmgbval->MGBoundaryFunction_[BoundaryFace::outer_x1]
                      != MGPeriodicOuterX1
    && mgroot_->pmgbval->MGBoundaryFunction_[BoundaryFace::outer_x1] != nullptr)
      MGBoundaryFunction_[BoundaryFace::outer_x1](u, time, nvar_,
                ngh+2, ngh+2, ngh,   ngh+1, ngh,   ngh+1, ngh, x0, y0, z0, dx, dy, dz);
  if (loc.lx2 == 0
    && mgroot_->pmgbval->MGBoundaryFunction_[BoundaryFace::inner_x2]
                      != MGPeriodicInnerX2
    && mgroot_->pmgbval->MGBoundaryFunction_[BoundaryFace::inner_x2] != nullptr)
      MGBoundaryFunction_[BoundaryFace::inner_x2](u, time, nvar_,
                0,     ngh+2, 0,     0,     ngh,   ngh+1, ngh, x0, y0, z0, dx, dy, dz);
  if (loc.lx2 == (nrbx2_<<lev)-1
    && mgroot_->pmgbval->MGBoundaryFunction_[BoundaryFace::outer_x2]
                      != MGPeriodicOuterX2
    && mgroot_->pmgbval->MGBoundaryFunction_[BoundaryFace::outer_x2] != nullptr)
      MGBoundaryFunction_[BoundaryFace::outer_x2](u, time, nvar_,
                0,     ngh+2, ngh+2, ngh+2, ngh,   ngh+1, ngh, x0, y0, z0, dx, dy, dz);
  if (loc.lx3 == 0
    && mgroot_->pmgbval->MGBoundaryFunction_[BoundaryFace::inner_x3]
                      != MGPeriodicInnerX3
    && mgroot_->pmgbval->MGBoundaryFunction_[BoundaryFace::inner_x3] != nullptr)
      MGBoundaryFunction_[BoundaryFace::inner_x3](u, time, nvar_,
                0,     ngh+2, 0,     ngh+2, 0,     0,     ngh, x0, y0, z0, dx, dy, dz);
  if (loc.lx3 == (nrbx3_<<lev)-1
    && mgroot_->pmgbval->MGBoundaryFunction_[BoundaryFace::outer_x3]
                      != MGPeriodicOuterX3
    && mgroot_->pmgbval->MGBoundaryFunction_[BoundaryFace::outer_x3] != nullptr)
      MGBoundaryFunction_[BoundaryFace::outer_x3](u, time, nvar_,
                0,     ngh+2, 0,     ngh+2, ngh+2, ngh+2, ngh, x0, y0, z0, dx, dy, dz);

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SetMeshBlockCoarsestBoundaries()
//  \brief Set boundaries of each MeshBlock from neighboring finer MeshBlocks

void MultigridDriver::SetMeshBlockCoarsestBoundaries() {
  int ngh = mgroot_->ngh_;
  Real vol = 0.125;
  for (auto itr = vmg_.begin(); itr<vmg_.end(); itr++) {
    Multigrid *pmg = *itr;
    AthenaArray<Real> &u = pmg->GetCurrentData();
    AthenaArray<Real> &uold = pmg->GetCurrentOldData();
    for(int ok=-1; ok<=1; ++ok) {
      for(int oj=-1; oj<=1; ++oj) {
        for(int oi=-1; oi<=1; ++oi) {
          if (pmg->pmgbval->nblevel[ok+1][oj+1][oi+1] > pmg->loc_.level) {
            LogicalLocation nloc;
            nloc.lx1 = pmg->loc_.lx1 + oi;
            nloc.lx2 = pmg->loc_.lx2 + oj;
            nloc.lx3 = pmg->loc_.lx3 + ok;
            nloc.level = pmg->loc_.level;
            for (int v=0; v<nvar_; ++v) {
              u(v, ok+ngh, oj+ngh, oi+ngh) = 0.0;
              if (ffas_)
                uold(v, ok+ngh, oj+ngh, oi+ngh) = 0.0;
            }
            GetNeighborOctetAverage(nloc, u, uold, ok+ngh, oj+ngh, oi+ngh, vol);
          }
        }
      }
    }

  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::GetNeighborOctetAverage(const LogicalLocation &nloc,
//           AthenaArray<Real> &u, AthenaArray<Real> &uold, int k, int j, int i, Real vol)
//  \brief Update the ghost cell with the average of all the octets including the leafs

void MultigridDriver::GetNeighborOctetAverage(const LogicalLocation &nloc,
     AthenaArray<Real> &u, AthenaArray<Real> &uold, int k, int j, int i, Real vol) {
  constexpr Real volfac = 0.125;
  int ngh = mgroot_->ngh_;
  const int l = ngh, r = ngh+1;
  int olev = nloc.level - (nrootlevel_ - 1);
  int oid = octetmap_[olev][nloc];
  AthenaArray<Real> &un = octets_[olev][oid].u;
  AthenaArray<Real> &unold = octets_[olev][oid].uold;
  if (octets_[olev][oid].fleaf == true) {
    for (int v=0; v<nvar_; ++v) {
      u(v, k, j, i) += vol*(un(v,l,l,l)+un(v,l,l,r)+un(v,l,r,l)+un(v,r,l,l)
                           +un(v,r,r,l)+un(v,r,l,r)+un(v,r,r,l)+un(v,r,r,r));
      if (ffas_)
        uold(v, k, j, i) += vol*
                         (unold(v,l,l,l)+unold(v,l,l,r)+unold(v,l,r,l)+unold(v,r,l,l)
                         +unold(v,r,r,l)+unold(v,r,l,r)+unold(v,r,r,l)+unold(v,r,r,r));
    }
  } else {
    for (int ok=0; ok<=1; ++ok) {
      for (int oj=0; oj<=1; ++oj) {
        for (int oi=0; oi<=1; ++oi) {
          LogicalLocation lloc;
          lloc.lx1 = (nloc.lx1 << 1) + oi;
          lloc.lx2 = (nloc.lx2 << 1) + oj;
          lloc.lx3 = (nloc.lx3 << 1) + ok;
          lloc.level = nloc.level + 1;
          int llev = olev + 1;
          if (octetmap_[llev].count(lloc) == 1) { // has leafs
            GetNeighborOctetAverage(lloc, u, uold, k, j, i, vol*volfac);
          } else {
            for (int v=0; v<nvar_; ++v) {
              u(v, k, j, i) += vol*un(v, ngh+ok, ngh+oj, ngh+oi);
              if (ffas_)
                uold(v, k, j, i) += vol*unold(v, ngh+ok, ngh+oj, ngh+oi);
            }
          }
        }
      }
    }
  }
}
