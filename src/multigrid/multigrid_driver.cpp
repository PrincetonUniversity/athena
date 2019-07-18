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
#include <iostream>   // endl
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
    mode_(0), // 0: FMG+V(1,1), 1: FMG+F(0,1), 2: V(1,1)
    nrbx1_(pm->nrbx1), nrbx2_(pm->nrbx2), nrbx3_(pm->nrbx3), pmy_mesh_(pm),
    fsubtract_average_(false), ffas_(false), eps_(-1.0), cbuf_(nvar_,3,3,3) {
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
  for (int n=0; n<pmy_mesh_->nbtotal; n++)
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
  for (int n=0; n<nranks_; n++) {
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

  if (nreflevel_ > 0 && pmy_mesh_->amr_updated) {
    for (int i=0; i<maxreflevel_; ++i) { // clear old data
      octetmap_[i].clear();
      prevnoct_[i] = noctets_[i];
      noctets_[i] = 0;
    }
    pmy_mesh_->tree.CountMGOctets(noctets_);
    for (int i=0; i<maxreflevel_; ++i) { // increase the octet array size if needed
      if (prevnoct_[i] < noctets_[i]) {
        octets_[i].resize(noctets_[i]);
        octetmap_[i].reserve(noctets_[i]);
      }
      noctets_[i] = 0;
    }
    pmy_mesh_->tree.GetMGOctetList(octets_, octetmap_, noctets_);
    for (int i=0; i<maxreflevel_; ++i) {
      for (int j=prevnoct_[i]; j<noctets_[i]; ++j) {
        octets_[i][j].u.NewAthenaArray(nvar_, ncoct, ncoct, ncoct);
        octets_[i][j].def.NewAthenaArray(nvar_, ncoct, ncoct, ncoct);
        octets_[i][j].src.NewAthenaArray(nvar_, ncoct, ncoct, ncoct);
        octets_[i][j].uold.NewAthenaArray(nvar_, ncoct, ncoct, ncoct);
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
    RestrictOctetsFMGSource();
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
  for (auto itr = vmg_.begin(); itr<vmg_.end(); itr++) {
    Multigrid *pmg = *itr;
    for (int v=0; v<nvar_; v++)
      rootbuf_[pmg->pmy_block_->gid*nvar_+v] = pmg->CalculateTotal(type, v);
  }
#ifdef MPI_PARALLEL
  MPI_Allgatherv(MPI_IN_PLACE, nblist_[Globals::my_rank]*nvar_, MPI_ATHENA_REAL,
                 rootbuf_, nvlisti_, nvslisti_, MPI_ATHENA_REAL, MPI_COMM_MULTIGRID);
#endif
  Real vol = (pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min)
           * (pmy_mesh_->mesh_size.x2max - pmy_mesh_->mesh_size.x2min)
           * (pmy_mesh_->mesh_size.x3max - pmy_mesh_->mesh_size.x3min);
  for (int v=0; v<nvar_; v++) {
    Real total = 0.0;
    for (int n=0; n<pmy_mesh_->nbtotal; n++)
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

  for (int n=0; n<pmy_mesh_->nbtotal; n++) {
    LogicalLocation &loc=pmy_mesh_->loclist[n];
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
      int olev = oloc.level-pmy_mesh_->root_level;
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
  AthenaArray<Real> &src=mgroot_->GetCurrentData();
  mgroot_->pmgbval->ApplyPhysicalBoundaries();
  if (pmy_mesh_->multilevel) {
    // *** implement later ***
  } else {
    for (auto itr = vmg_.begin(); itr<vmg_.end(); itr++) {
      Multigrid *pmg = *itr;
      LogicalLocation &loc=pmg->pmy_block_->loc;
      // KGF: possible std::int64_t overflow. Practically never occur with Multigrid
      pmg->SetFromRootGrid(src, static_cast<int>(loc.lx3), static_cast<int>(loc.lx2),
                           static_cast<int>(loc.lx1));
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::FMGProlongate()
//  \brief Prolongation for FMG Cycle

void MultigridDriver::FMGProlongate() {
  int flag=0;
  if (current_level_ == nrootlevel_-1) {
    TransferFromRootToBlocks();
    flag=1;
  }
  if (current_level_ >= nrootlevel_-1) {
    mgtlist_->SetMGTaskListFMGProlongate(flag);
    mgtlist_->DoTaskListOneStage(this);
//  } else if ( ) { // non uniform root
  } else { // root grid
    mgroot_->pmgbval->ApplyPhysicalBoundaries();
    mgroot_->FMGProlongateBlock();
  }
  current_level_++;
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::OneStepToFiner(int nsmooth)
//  \brief smoothing and restriction one level

void MultigridDriver::OneStepToFiner(int nsmooth) {
  int ngh=mgroot_->ngh_;
  int flag=0;
  if (current_level_ == nrootlevel_-1) {
    TransferFromRootToBlocks();
    flag=1;
  }
  if (current_level_ >= nrootlevel_-1) {
    if (current_level_ == ntotallevel_-2) flag=2;
    mgtlist_->SetMGTaskListToFiner(nsmooth, ngh, flag);
    mgtlist_->DoTaskListOneStage(this);
//  } else if ( ) { // non uniform root
  } else { // root grid
    mgroot_->pmgbval->ApplyPhysicalBoundaries();
    mgroot_->ProlongateAndCorrectBlock();
    for (int n=0; n<nsmooth; n++) {
      mgroot_->pmgbval->ApplyPhysicalBoundaries();
      mgroot_->SmoothBlock(0);
      mgroot_->pmgbval->ApplyPhysicalBoundaries();
      mgroot_->SmoothBlock(1);
    }
  }
  current_level_++;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::OneStepToCoarser(int nsmooth)
//  \brief smoothing and restriction one level

void MultigridDriver::OneStepToCoarser(int nsmooth) {
  int ngh=mgroot_->ngh_;
  if (current_level_ >= nrootlevel_) {
    mgtlist_->SetMGTaskListToCoarser(nsmooth, ngh);
    mgtlist_->DoTaskListOneStage(this);
    if (current_level_ == nrootlevel_) {
      TransferFromBlocksToRoot();
      if (!ffas_)
        mgroot_->ZeroClearData();
    }
//  } else if ( ) { // non uniform root
  } else { // uniform root grid
    mgroot_->pmgbval->ApplyPhysicalBoundaries();
    if (ffas_ && current_level_ < fmglevel_) {
      mgroot_->StoreOldData();
      mgroot_->CalculateFASRHSBlock();
    }
    for (int n=0; n<nsmooth; n++) {
      mgroot_->SmoothBlock(0);
      mgroot_->pmgbval->ApplyPhysicalBoundaries();
      mgroot_->SmoothBlock(1);
      mgroot_->pmgbval->ApplyPhysicalBoundaries();
    }
    mgroot_->RestrictBlock();
  }
  current_level_--;
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
    if (fmglevel_!=ntotallevel_-1) FMGProlongate();
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
  std::cout << std::scientific;
  while (def > eps_) {
    SolveVCycle(1,1);
    Real olddef=def;
    def=0.0;
    for (int n=0; n<nvar_; n++)
      def+=CalculateDefectNorm(MGNormType::l2, n);
    if (niter > 0 && def/olddef > 0.5) {
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
      for (int v=0; v<nvar_; v++) {
        Real ave=mgroot_->CalculateTotal(MGVariable::u, v)/vol;
        mgroot_->SubtractAverage(MGVariable::u, v, ave);
      }
    }
    mgroot_->pmgbval->ApplyPhysicalBoundaries();
    if (ffas_ && current_level_ < fmglevel_) {
      mgroot_->StoreOldData();
      mgroot_->CalculateFASRHSBlock();
    }
    for (int i=0; i<ni; i++) { // iterate ni times
      mgroot_->SmoothBlock(0);
      mgroot_->pmgbval->ApplyPhysicalBoundaries();
      mgroot_->SmoothBlock(1);
      mgroot_->pmgbval->ApplyPhysicalBoundaries();
    }
    if (fsubtract_average_) {
      Real vol=(mgroot_->size_.x1max-mgroot_->size_.x1min)
              *(mgroot_->size_.x2max-mgroot_->size_.x2min)
              *(mgroot_->size_.x3max-mgroot_->size_.x3min);
      for (int v=0; v<nvar_; v++) {
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
//! \fn void MultigridDriver::RestrictOctetsFMGSource()
//  \brief restrict the source in octets for FMG

void MultigridDriver::RestrictOctetsFMGSource() {
  if (maxreflevel_ > 0) {
    int ngh = mgroot_->ngh_;
    int is, ie, js, je, ks, ke;
    is = js = ks = ngh;
    ie = je = ke = ngh + 1;
    for (int l=nreflevel_-1; l>=1; --l) {
      for (int i=0; i<noctets_[l]; ++i) {
        LogicalLocation &loc = octets_[l][i].loc;
        LogicalLocation cloc;
        cloc.lx1 = (loc.lx1 >> 1);
        cloc.lx2 = (loc.lx2 >> 1);
        cloc.lx3 = (loc.lx3 >> 1);
        cloc.level = loc.level - 1;
        int clev = cloc.level - pmy_mesh_->root_level;
        int oi = (static_cast<int>(loc.lx1) & 1) + ngh;
        int oj = (static_cast<int>(loc.lx2) & 1) + ngh;
        int ok = (static_cast<int>(loc.lx3) & 1) + ngh;
        int oid = octetmap_[clev][cloc];
        mgroot_->Restrict(cbuf_, octets_[l][i].src, is, ie, js, je, ks, ke);
        for (int v=0; v<nvar_; ++v)
          octets_[clev][oid].src(v, ok, oj, oi) = cbuf_(v, ngh, ngh, ngh);
      }
    }
    for (int i=0; i<noctets_[0]; ++i) {
      LogicalLocation &loc = octets_[0][i].loc;
      int ri = (loc.lx1 >> 1);
      int rj = (loc.lx2 >> 1);
      int rk = (loc.lx3 >> 1);
      mgroot_->Restrict(cbuf_, octets_[0][i].src, is, ie, js, je, ks, ke);
      for (int v=0; v<nvar_; ++v)
        mgroot_->SetData(MGVariable::src, v, rk, rj, ri, cbuf_(v, ngh, ngh, ngh));
    }
  }

  return;
}

