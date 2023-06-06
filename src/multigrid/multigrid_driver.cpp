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
#include <iomanip>    // setprecision
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

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

// constructor, initializes data structures and parameters

MultigridDriver::MultigridDriver(Mesh *pm, MGBoundaryFunc *MGBoundary,
                                 MGSourceMaskFunc MGSourceMask, int invar) :
    nranks_(Globals::nranks), nthreads_(pm->num_mesh_threads_), nbtotal_(pm->nbtotal),
    nvar_(invar), mode_(0), // 0: FMG V(1,1) + iterative, 1: V(1,1) iterative
    maxreflevel_(pm->multilevel?pm->max_level-pm->root_level:0),
    nrbx1_(pm->nrbx1), nrbx2_(pm->nrbx2), nrbx3_(pm->nrbx3), srcmask_(MGSourceMask),
    pmy_mesh_(pm), fsubtract_average_(false), ffas_(pm->multilevel), needinit_(true),
    eps_(-1.0), niter_(-1), coffset_(0), mporder_(-1), nmpcoeff_(0), mpo_(3),
    autompo_(false), nodipole_(false), nb_rank_(0) {
  std::cout << std::scientific << std::setprecision(15);

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

  ranklist_  = new int[nbtotal_];
  rootbuf_=new Real[nbtotal_*nvar_*2];
  for (int n = 0; n < nbtotal_; ++n)
    ranklist_[n]=pmy_mesh_->ranklist[n];
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

  if (maxreflevel_ > 0) { // SMR / AMR
    octets_ = new std::vector<MGOctet>[maxreflevel_];
    octetmap_ = new std::unordered_map<LogicalLocation, int,
                                       LogicalLocationHash>[maxreflevel_];
    octetbflag_ = new std::vector<bool>[maxreflevel_];
    noctets_ = new int[maxreflevel_]();
    pmaxnoct_ = new int[maxreflevel_]();
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
  if (maxreflevel_ > 0) {
    delete [] octets_;
    delete [] octetmap_;
    delete [] octetbflag_;
    delete [] noctets_;
    delete [] pmaxnoct_;
  }
  if (mporder_ > 0)
    delete [] mpcoeff_;
#ifdef MPI_PARALLEL
  MPI_Comm_free(&MPI_COMM_MULTIGRID);
#endif
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::CheckBoundaryFunctions()
//  \brief check boundary functions and set some internal flags.

void MultigridDriver::CheckBoundaryFunctions() {
  fsubtract_average_ = true;
  switch(mg_mesh_bcs_[BoundaryFace::inner_x1]) {
    case BoundaryFlag::user:
      if (MGBoundaryFunction_[BoundaryFace::inner_x1] == nullptr) {
        std::stringstream msg;
        msg << "### FATAL ERROR in MGGravityDriver::CheckBoundaryFunctions" << std::endl
            << "A user-defined boundary condition is specified for " << std::endl
            << "inner_x1 but no function is enrolled." << std::endl;
        ATHENA_ERROR(msg);
      }
      fsubtract_average_ = false;
      break;
    case BoundaryFlag::periodic:
    case BoundaryFlag::mg_zerograd:
    case BoundaryFlag::mg_zerofixed:
      break;
    case BoundaryFlag::mg_multipole:
      mporder_ = 0;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in MGGravityDriver::CheckBoundaryFunctions" << std::endl
          << "Invalid or no boundary type is specified." << std::endl;
      ATHENA_ERROR(msg);
      break;
  }
  switch(mg_mesh_bcs_[BoundaryFace::outer_x1]) {
    case BoundaryFlag::user:
      if (MGBoundaryFunction_[BoundaryFace::outer_x1] == nullptr) {
        std::stringstream msg;
        msg << "### FATAL ERROR in MGGravityDriver::CheckBoundaryFunctions" << std::endl
            << "A user-defined boundary condition is specified for " << std::endl
            << "outer_x1 but no function is enrolled." << std::endl;
        ATHENA_ERROR(msg);
      }
      fsubtract_average_ = false;
      break;
    case BoundaryFlag::periodic:
    case BoundaryFlag::mg_zerograd:
    case BoundaryFlag::mg_zerofixed:
      break;
    case BoundaryFlag::mg_multipole:
      mporder_ = 0;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in MGGravityDriver::CheckBoundaryFunctions" << std::endl
          << "Invalid or no boundary type is specified." << std::endl;
      ATHENA_ERROR(msg);
      break;
  }
  switch(mg_mesh_bcs_[BoundaryFace::inner_x2]) {
    case BoundaryFlag::user:
      if (MGBoundaryFunction_[BoundaryFace::inner_x2] == nullptr) {
        std::stringstream msg;
        msg << "### FATAL ERROR in MGGravityDriver::CheckBoundaryFunctions" << std::endl
            << "A user-defined boundary condition is specified for " << std::endl
            << "inner_x2 but no function is enrolled." << std::endl;
        ATHENA_ERROR(msg);
      }
      fsubtract_average_ = false;
      break;
    case BoundaryFlag::periodic:
    case BoundaryFlag::mg_zerograd:
    case BoundaryFlag::mg_zerofixed:
      break;
    case BoundaryFlag::mg_multipole:
      mporder_ = 0;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in MGGravityDriver::CheckBoundaryFunctions" << std::endl
          << "Invalid or no boundary type is specified." << std::endl;
      ATHENA_ERROR(msg);
      break;
  }
  switch(mg_mesh_bcs_[BoundaryFace::outer_x2]) {
    case BoundaryFlag::user:
      if (MGBoundaryFunction_[BoundaryFace::outer_x2] == nullptr) {
        std::stringstream msg;
        msg << "### FATAL ERROR in MGGravityDriver::CheckBoundaryFunctions" << std::endl
            << "A user-defined boundary condition is specified for " << std::endl
            << "outer_x2 but no function is enrolled." << std::endl;
        ATHENA_ERROR(msg);
      }
      fsubtract_average_ = false;
      break;
    case BoundaryFlag::periodic:
    case BoundaryFlag::mg_zerograd:
    case BoundaryFlag::mg_zerofixed:
      break;
    case BoundaryFlag::mg_multipole:
      mporder_ = 0;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in MGGravityDriver::CheckBoundaryFunctions" << std::endl
          << "Invalid or no boundary type is specified." << std::endl;
      ATHENA_ERROR(msg);
      break;
  }
  switch(mg_mesh_bcs_[BoundaryFace::inner_x3]) {
    case BoundaryFlag::user:
      if (MGBoundaryFunction_[BoundaryFace::inner_x3] == nullptr) {
        std::stringstream msg;
        msg << "### FATAL ERROR in MGGravityDriver::CheckBoundaryFunctions" << std::endl
            << "A user-defined boundary condition is specified for " << std::endl
            << "inner_x3 but no function is enrolled." << std::endl;
        ATHENA_ERROR(msg);
      }
      fsubtract_average_ = false;
      break;
    case BoundaryFlag::periodic:
    case BoundaryFlag::mg_zerograd:
    case BoundaryFlag::mg_zerofixed:
      break;
    case BoundaryFlag::mg_multipole:
      mporder_ = 0;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in MGGravityDriver::CheckBoundaryFunctions" << std::endl
          << "Invalid or no boundary type is specified." << std::endl;
      ATHENA_ERROR(msg);
      break;
  }
  switch(mg_mesh_bcs_[BoundaryFace::outer_x3]) {
    case BoundaryFlag::user:
      if (MGBoundaryFunction_[BoundaryFace::outer_x3] == nullptr) {
        std::stringstream msg;
        msg << "### FATAL ERROR in MGGravityDriver::CheckBoundaryFunctions" << std::endl
            << "A user-defined boundary condition is specified for " << std::endl
            << "outer_x3 but no function is enrolled." << std::endl;
        ATHENA_ERROR(msg);
      }
      fsubtract_average_ = false;
      break;
    case BoundaryFlag::periodic:
    case BoundaryFlag::mg_zerograd:
    case BoundaryFlag::mg_zerofixed:
      break;
    case BoundaryFlag::mg_multipole:
      mporder_ = 0;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in MGGravityDriver::CheckBoundaryFunctions" << std::endl
          << "Invalid or no boundary type is specified." << std::endl;
      ATHENA_ERROR(msg);
      break;
  }

  if (mporder_ >= 0) {
    ffas_ = true;
    fsubtract_average_ = false;
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SubtractAverage(MGVariable type)
//  \brief Calculate the global average and subtract it

void MultigridDriver::SubtractAverage(MGVariable type) {
#pragma omp parallel for num_threads(nthreads_)
  for (auto itr = vmg_.begin(); itr < vmg_.end(); itr++) {
    Multigrid *pmg = *itr;
    for (int v=0; v<nvar_; ++v)
      rootbuf_[pmg->pmy_block_->gid*nvar_+v] = pmg->CalculateTotal(type, v);
  }
#ifdef MPI_PARALLEL
  if (nb_rank_ > 0)  // every rank has the same number of MeshBlocks
    MPI_Allgather(MPI_IN_PLACE, nb_rank_*nvar_, MPI_ATHENA_REAL,
                  rootbuf_, nb_rank_*nvar_, MPI_ATHENA_REAL, MPI_COMM_MULTIGRID);
  else
    MPI_Allgatherv(MPI_IN_PLACE, nblist_[Globals::my_rank]*nvar_, MPI_ATHENA_REAL,
                   rootbuf_, nvlisti_, nvslisti_, MPI_ATHENA_REAL, MPI_COMM_MULTIGRID);
#endif
  Real vol = (pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min)
           * (pmy_mesh_->mesh_size.x2max - pmy_mesh_->mesh_size.x2min)
           * (pmy_mesh_->mesh_size.x3max - pmy_mesh_->mesh_size.x3min);
  for (int v=0; v<nvar_; ++v) {
    Real total = 0.0;
    for (int n = 0; n < nbtotal_; ++n)
      total += rootbuf_[n*nvar_+v];
    last_ave_ = total/vol;
#pragma omp parallel for num_threads(nthreads_)
    for (auto itr = vmg_.begin(); itr < vmg_.end(); itr++) {
      Multigrid *pmg = *itr;
      pmg->SubtractAverage(type, v, last_ave_);
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SetupMultigrid()
//  \brief initialize the source assuming that the source terms are already loaded

void MultigridDriver::SetupMultigrid() {
  locrootlevel_ = pmy_mesh_->root_level;
  nrootlevel_ = mgroot_->GetNumberOfLevels();
  nmblevel_ = vmg_[0]->GetNumberOfLevels();
  nreflevel_ = pmy_mesh_->current_level - locrootlevel_;
  ntotallevel_ = nrootlevel_ + nmblevel_ + nreflevel_ - 1;
  fmglevel_ = current_level_ = ntotallevel_ - 1;
  int ncoct = mgroot_->ngh_*2 + 2, nccoct = mgroot_->ngh_*2 + 1;
  os_ = mgroot_->ngh_;
  oe_ = os_+1;

  if (pmy_mesh_->amr_updated)
    needinit_ = true;

  // note: the level of an Octet is one level lower than the data stored there
  if (nreflevel_ > 0 && needinit_) {
    for (int l = 0; l < nreflevel_; ++l) { // clear old data
      octetmap_[l].clear();
      pmaxnoct_[l] = std::max(pmaxnoct_[l], noctets_[l]);
      noctets_[l] = 0;
    }
    pmy_mesh_->tree.CountMGOctets(noctets_);
    for (int l = 0; l < nreflevel_; ++l) { // increase the octet array size if needed
      if (pmaxnoct_[l] < noctets_[l]) {
        octets_[l].resize(noctets_[l]);
        octetmap_[l].reserve(noctets_[l]);
        octetbflag_[l].resize(noctets_[l]);
      }
      noctets_[l] = 0;
    }
    pmy_mesh_->tree.GetMGOctetList(octets_, octetmap_, noctets_);
    for (int l = 0; l < nreflevel_; ++l) {
      for (int o = pmaxnoct_[l]; o < noctets_[l]; ++o) {
        octets_[l][o].u.NewAthenaArray(nvar_, ncoct, ncoct, ncoct);
        octets_[l][o].def.NewAthenaArray(nvar_, ncoct, ncoct, ncoct);
        octets_[l][o].src.NewAthenaArray(nvar_, ncoct, ncoct, ncoct);
        octets_[l][o].coord.AllocateMGCoordinates(ncoct, ncoct, ncoct);
        octets_[l][o].ccoord.AllocateMGCoordinates(nccoct, nccoct, nccoct);
        if (ffas_)
          octets_[l][o].uold.NewAthenaArray(nvar_, ncoct, ncoct, ncoct);
      }
    }
  }

  if (needinit_) {
    // reallocate buffers if needed
    if (nbtotal_ != pmy_mesh_->nbtotal) {
      if (nbtotal_ < pmy_mesh_->nbtotal) {
        delete [] ranklist_;
        delete [] rootbuf_;
        ranklist_ = new int[pmy_mesh_->nbtotal];
        rootbuf_ = new Real[pmy_mesh_->nbtotal*nvar_*2];
      }
      nbtotal_ = pmy_mesh_->nbtotal;
    }
    nb_rank_ = pmy_mesh_->nblist[0];
    for (int n = 1; n < nranks_; ++n) {
      if (nb_rank_ != pmy_mesh_->nblist[n]) {
        nb_rank_ = 0;
        break;
      }
    }

    // Setting up the MPI information
    // *** this part should be modified when dedicate processes are allocated ***
    // *** we also need to construct another neighbor list for Multigrid ***

    // assume the same parallelization as hydro
    for (int n = 0; n < nbtotal_; ++n)
      ranklist_[n] = pmy_mesh_->ranklist[n];
    for (int n = 0; n < nranks_; ++n) {
      nslist_[n]  = pmy_mesh_->nslist[n];
      nblist_[n]  = pmy_mesh_->nblist[n];
      nvslist_[n] = nslist_[n]*nvar_*2;
      nvlist_[n]  = nblist_[n]*nvar_*2;
      nvslisti_[n] = nslist_[n]*nvar_;
      nvlisti_[n]  = nblist_[n]*nvar_;
    }
    for (Multigrid* pmg : vmg_) {
      pmg->pmgbval->SearchAndSetNeighbors(pmy_mesh_->tree, ranklist_, nslist_);
      pmg->pmgbval->bcolor_ = 0;
    }
    if (nreflevel_ > 0)
      CalculateOctetCoordinates();
    needinit_ = false;
  }

  if (srcmask_ != nullptr) {
#pragma omp parallel for num_threads(nthreads_)
    for (auto itr = vmg_.begin(); itr < vmg_.end(); itr++) {
      Multigrid *pmg = *itr;
      pmg->ApplySourceMask();
    }
  }

  if (fsubtract_average_)
    SubtractAverage(MGVariable::src);

  if (mporder_ > 0) {
    if (autompo_)
      CalculateCenterOfMass();
    CalculateMultipoleCoefficients();
  }

  if (mode_ == 0) { // FMG
#pragma omp parallel for num_threads(nthreads_)
    for (auto itr = vmg_.begin(); itr < vmg_.end(); itr++) {
      Multigrid *pmg = *itr;
      pmg->RestrictFMGSource();
    }
    TransferFromBlocksToRoot(true);
    RestrictFMGSourceOctets();
    mgroot_->RestrictFMGSource();
    current_level_ = 0;
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::TransferFromBlocksToRoot(bool initflag)
//  \brief collect the coarsest data and transfer to the root grid

void MultigridDriver::TransferFromBlocksToRoot(bool initflag) {
  int nv = nvar_, ngh = mgroot_->ngh_;
  if (ffas_ && !initflag) nv*=2;
#pragma omp parallel for num_threads(nthreads_)
  for (auto itr = vmg_.begin(); itr < vmg_.end(); itr++) {
    Multigrid *pmg = *itr;
    for (int v = 0; v < nvar_; ++v)
      rootbuf_[pmg->pmy_block_->gid*nv+v]=pmg->GetCoarsestData(MGVariable::src, v);
    if (ffas_ && !initflag) {
      for (int v = 0; v < nvar_; ++v)
        rootbuf_[pmg->pmy_block_->gid*nv+nvar_+v]=pmg->GetCoarsestData(MGVariable::u, v);
    }
  }

#ifdef MPI_PARALLEL
  if (nb_rank_ > 0) { // every rank has the same number of MeshBlocks
    MPI_Allgather(MPI_IN_PLACE, nb_rank_*nv, MPI_ATHENA_REAL,
                  rootbuf_, nb_rank_*nv, MPI_ATHENA_REAL, MPI_COMM_MULTIGRID);
  } else {
    if (!initflag)
      MPI_Allgatherv(MPI_IN_PLACE, nblist_[Globals::my_rank]*nv, MPI_ATHENA_REAL,
                     rootbuf_, nvlist_, nvslist_, MPI_ATHENA_REAL, MPI_COMM_MULTIGRID);
    else
      MPI_Allgatherv(MPI_IN_PLACE, nblist_[Globals::my_rank]*nvar_, MPI_ATHENA_REAL,
                     rootbuf_, nvlisti_, nvslisti_, MPI_ATHENA_REAL, MPI_COMM_MULTIGRID);
  }
#endif

#pragma omp parallel for num_threads(nthreads_)
  for (int n = 0; n < nbtotal_; ++n) {
    const LogicalLocation &loc=pmy_mesh_->loclist[n];
    int i = static_cast<int>(loc.lx1);
    int j = static_cast<int>(loc.lx2);
    int k = static_cast<int>(loc.lx3);
    if (loc.level == locrootlevel_) {
      for (int v = 0; v < nvar_; ++v)
        mgroot_->SetData(MGVariable::src, v, k, j, i, rootbuf_[n*nv+v]);
      if (ffas_ && !initflag) {
        for (int v = 0; v < nvar_; ++v)
          mgroot_->SetData(MGVariable::u, v, k, j, i, rootbuf_[n*nv+nvar_+v]);
      }
    } else {
      LogicalLocation oloc;
      oloc.lx1 = (loc.lx1 >> 1);
      oloc.lx2 = (loc.lx2 >> 1);
      oloc.lx3 = (loc.lx3 >> 1);
      oloc.level = loc.level - 1;
      int olev = oloc.level - locrootlevel_;
      int oid = octetmap_[olev][oloc];
      int oi = (i&1) + ngh;
      int oj = (j&1) + ngh;
      int ok = (k&1) + ngh;
      for (int v = 0; v < nvar_; ++v)
        octets_[olev][oid].src(v,ok,oj,oi) = rootbuf_[n*nv+v];
      if (ffas_ && !initflag) {
        for (int v = 0; v < nvar_; ++v)
          octets_[olev][oid].u(v,ok,oj,oi) = rootbuf_[n*nv+nvar_+v];
      }
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::TransferFromRootToBlocks(bool folddata)
//  \brief Transfer the data from the root grid to the coarsest level of each MeshBlock

void MultigridDriver::TransferFromRootToBlocks(bool folddata) {
  if (nreflevel_ > 0) {
    RestrictOctetsBeforeTransfer();
    SetOctetBoundariesBeforeTransfer(folddata);
  }
#pragma omp parallel for num_threads(nthreads_)
  for (auto itr = vmg_.begin(); itr < vmg_.end(); itr++) {
    Multigrid *pmg = *itr;
    pmg->SetFromRootGrid(folddata);
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::FMGProlongate()
//  \brief Prolongation for FMG Cycle

void MultigridDriver::FMGProlongate() {
  int flag=0;
  if (current_level_ == nrootlevel_ + nreflevel_ - 1) {
    mgroot_->pmgbval->ApplyPhysicalBoundaries();
    TransferFromRootToBlocks(false);
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
  if (current_level_ == nrootlevel_ + nreflevel_ - 1) {
    mgroot_->pmgbval->ApplyPhysicalBoundaries();
    TransferFromRootToBlocks(ffas_);
    flag=1;
  }
  if (current_level_ >= nrootlevel_ + nreflevel_ - 1) { // MeshBlocks
    if (current_level_ == ntotallevel_ - 2) flag=2;
    mgtlist_->SetMGTaskListToFiner(nsmooth, ngh, flag);
    mgtlist_->DoTaskListOneStage(this);
    current_level_++;
  } else if (current_level_ >= nrootlevel_ - 1) { // non uniform octets
    if (current_level_ == nrootlevel_ - 1)
      mgroot_->pmgbval->ApplyPhysicalBoundaries();
    else
      SetBoundariesOctets(true, ffas_);
    ProlongateAndCorrectOctets();
    current_level_++;
    for (int n = 0; n < nsmooth; ++n) {
      SetBoundariesOctets(false, false);
      SmoothOctets(0);
      SetBoundariesOctets(false, false);
      SmoothOctets(1);
    }
  } else { // root grid
    mgroot_->pmgbval->ApplyPhysicalBoundaries();
    mgroot_->ProlongateAndCorrectBlock();
    current_level_++;
    for (int n = 0; n < nsmooth; ++n) {
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
  if (current_level_ >= nrootlevel_ + nreflevel_) { // MeshBlocks
    mgtlist_->SetMGTaskListToCoarser(nsmooth, ngh);
    mgtlist_->DoTaskListOneStage(this);
    if (current_level_ == nrootlevel_ + nreflevel_) {
      TransferFromBlocksToRoot();
      if (!ffas_) {
        mgroot_->ZeroClearData();
        if (nreflevel_ > 0)
          ZeroClearOctets();
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
    if (!ffas_ && current_level_ == fmglevel_) {
      mgroot_->ZeroClearData();
      ZeroClearOctets();
    }
  } else { // uniform root grid
    mgroot_->pmgbval->ApplyPhysicalBoundaries();
    if (ffas_ && current_level_ < fmglevel_) {
      mgroot_->StoreOldData();
      mgroot_->CalculateFASRHSBlock();
    }
    for (int n = 0; n < nsmooth; ++n) {
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
  coffset_ ^= 1;
  while (current_level_ > 0)
    OneStepToCoarser(npresmooth);
  SolveCoarsestGrid();
  while (current_level_ < startlevel)
    OneStepToFiner(npostsmooth);
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SolveFMGCycle()
//  \brief Solve the FMG Cycle using the V(1,1) or F(0,1) cycle

void MultigridDriver::SolveFMGCycle() {
  for (fmglevel_ = 0; fmglevel_ < ntotallevel_; fmglevel_++) {
    SolveVCycle(1, 1);
    if (fmglevel_ != ntotallevel_-1)
      FMGProlongate();
  }
  fmglevel_ = ntotallevel_ - 1;
  if (fsubtract_average_)
    SubtractAverage(MGVariable::u);
  if (eps_ >= 0.0)
    SolveIterative();
  else
    SolveIterativeFixedTimes();
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SolveIterative()
//  \brief Solve iteratively until the convergence is achieved

void MultigridDriver::SolveIterative() {
  int n = 0;
  Real def = 0.0;
  for (int v = 0; v < nvar_; ++v)
    def += CalculateDefectNorm(MGNormType::l2, v);
  while (def > eps_) {
    SolveVCycle(1, 1);
    Real olddef = def;
    def = 0.0;
    for (int v = 0; v < nvar_; ++v)
      def += CalculateDefectNorm(MGNormType::l2, v);
    //if (Globals::my_rank == 0)
    //  std::cout << "niter " << n << " def " << def << std::endl;
    if (def/olddef > 0.8) {
      if (eps_ == 0.0) break;
      if (Globals::my_rank == 0)
        std::cout << "### Warning in MultigridDriver::SolveIterative" << std::endl
                  << "Slow multigrid convergence : defect norm = " << def
                  << ", convergence factor = " << def/olddef << "." << std::endl;
    }
    if (n > 100) {
      if (Globals::my_rank == 0) {
        std::cout
            << "### Warning in MultigridDriver::SolveIterative" << std::endl
            << "Aborting because the # iterations is too large, n > 100." << std::endl
            << "Check the solution as it may not be accurate enough." << std::endl;
      }
      break;
    }
    n++;
  }
  if (fsubtract_average_)
    SubtractAverage(MGVariable::u);
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SolveIterativeFixedTimes()
//  \brief Solve iteratively niter_ times

void MultigridDriver::SolveIterativeFixedTimes() {
  for (int n = 0; n < niter_; ++n)
    SolveVCycle(1, 1);
  if (fsubtract_average_)
    SubtractAverage(MGVariable::u);
  for (int v = 0; v < nvar_; ++v)
    CalculateDefectNorm(MGNormType::l2, v);

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
    for (int i = 0; i < ni; ++i) { // iterate ni times
      mgroot_->SmoothBlock(0);
      mgroot_->pmgbval->ApplyPhysicalBoundaries();
      mgroot_->SmoothBlock(1);
      mgroot_->pmgbval->ApplyPhysicalBoundaries();
    }
    if (fsubtract_average_) {
      Real vol=(mgroot_->size_.x1max-mgroot_->size_.x1min)
              *(mgroot_->size_.x2max-mgroot_->size_.x2min)
              *(mgroot_->size_.x3max-mgroot_->size_.x3min);
      for (int v = 0; v < nvar_; ++v) {
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

  if (nrm == MGNormType::max) {
#pragma omp parallel for reduction(max : norm) num_threads(nthreads_)
    for (auto itr = vmg_.begin(); itr < vmg_.end(); itr++) {
      Multigrid *pmg = *itr;
      norm = std::max(norm, pmg->CalculateDefectNorm(nrm, n));
    }
  } else {
#pragma omp parallel for reduction(+ : norm) num_threads(nthreads_)
    for (auto itr = vmg_.begin(); itr < vmg_.end(); itr++) {
      Multigrid *pmg = *itr;
      norm += pmg->CalculateDefectNorm(nrm, n);
    }
  }
#ifdef MPI_PARALLEL
  if (nrm == MGNormType::max)
    MPI_Allreduce(MPI_IN_PLACE,&norm,1,MPI_ATHENA_REAL,MPI_MAX,MPI_COMM_MULTIGRID);
  else
    MPI_Allreduce(MPI_IN_PLACE,&norm,1,MPI_ATHENA_REAL,MPI_SUM,MPI_COMM_MULTIGRID);
#endif
  if (nrm != MGNormType::max) {
    Real vol = (mgroot_->size_.x1max-mgroot_->size_.x1min)
             * (mgroot_->size_.x2max-mgroot_->size_.x2min)
             * (mgroot_->size_.x3max-mgroot_->size_.x3min);
    norm /= vol;
  }
  if (nrm == MGNormType::l2)
    norm = std::sqrt(norm);

  return norm;
}


//----------------------------------------------------------------------------------------
//! \fn Multigrid* MultigridDriver::FindMultigrid(int tgid)
//  \brief return the Multigrid whose gid is tgid

Multigrid* MultigridDriver::FindMultigrid(int tgid) {
  int first = vmg_[0]->pmy_block_->gid;
  if (tgid > nblist_[Globals::my_rank] + first)
    return nullptr;
  return vmg_[tgid-first];
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::CalculateOctetCoordinates()
//  \brief calculate coordinates for Octets

void MultigridDriver::CalculateOctetCoordinates() {
  RegionSize size, csize;
  int ngh = mgroot_->ngh_;
  size.nx1  = 2, size.nx2  = 2, size.nx3  = 2;
  csize.nx1 = 1, csize.nx2 = 1, csize.nx3 = 1;
#pragma omp parallel for num_threads(nthreads_)
  for (int o = 0; o < noctets_[0]; ++o) {
    MGCoordinates &coord = mgroot_->coord_[mgroot_->nlevel_-1];
    LogicalLocation &loc = octets_[0][o].loc;
    int i = static_cast<int>(loc.lx1) + ngh;
    int j = static_cast<int>(loc.lx2) + ngh;
    int k = static_cast<int>(loc.lx3) + ngh;
    size.x1min = csize.x1min = coord.x1f(i);
    size.x1max = csize.x1max = coord.x1f(i+1);
    size.x2min = csize.x2min = coord.x2f(j);
    size.x2max = csize.x2max = coord.x2f(j+1);
    size.x3min = csize.x3min = coord.x3f(k);
    size.x3max = csize.x3max = coord.x3f(k+1);
    octets_[0][o].coord.CalculateMGCoordinates(size, 0, ngh);
    octets_[0][o].ccoord.CalculateMGCoordinates(csize, 0, ngh);
  }
  for (int l = 1; l < nreflevel_; l++) {
#pragma omp parallel for num_threads(nthreads_)
    for (int o = 0; o < noctets_[l]; ++o) {
      LogicalLocation &loc = octets_[l][o].loc;
      LogicalLocation cloc;
      cloc.lx1 = (loc.lx1 >> 1);
      cloc.lx2 = (loc.lx2 >> 1);
      cloc.lx3 = (loc.lx3 >> 1);
      cloc.level = loc.level - 1;
      int oid = octetmap_[l-1][cloc];
      MGCoordinates &coord = octets_[l-1][oid].coord;
      int i = static_cast<int>(loc.lx1&1) + ngh;
      int j = static_cast<int>(loc.lx2&1) + ngh;
      int k = static_cast<int>(loc.lx3&1) + ngh;
      size.x1min = csize.x1min = coord.x1f(i);
      size.x1max = csize.x1max = coord.x1f(i+1);
      size.x2min = csize.x2min = coord.x2f(j);
      size.x2max = csize.x2max = coord.x2f(j+1);
      size.x3min = csize.x3min = coord.x3f(k);
      size.x3max = csize.x3max = coord.x3f(k+1);
      octets_[l][o].coord.CalculateMGCoordinates(size, 0, ngh);
      octets_[l][o].ccoord.CalculateMGCoordinates(csize, 0, ngh);
    }
  }
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::RestrictFMGSourceOctets()
//  \brief restrict the source in octets for FMG

void MultigridDriver::RestrictFMGSourceOctets() {
  if (nreflevel_ > 0) {
    const int &ngh = mgroot_->ngh_;
    for (int l = nreflevel_ - 1; l >= 1; --l) {  // fine octets to coarse octets
#pragma omp parallel for num_threads(nthreads_)
      for (int o = 0; o < noctets_[l]; ++o) {
        AthenaArray<Real> cbuf;
        cbuf.NewAthenaArray(nvar_, 3, 3, 3);
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
        mgroot_->Restrict(cbuf, octets_[l][o].src, ngh, ngh, ngh, ngh, ngh, ngh, false);
        for (int v = 0; v < nvar_; ++v)
          octets_[l-1][oid].src(v, ok, oj, oi) = cbuf(v, ngh, ngh, ngh);
      }
    }
#pragma omp parallel for num_threads(nthreads_)
    for (int o = 0; o < noctets_[0]; ++o) { // octets to the root grid
      AthenaArray<Real> cbuf;
      cbuf.NewAthenaArray(nvar_, 3, 3, 3);
      const LogicalLocation &loc = octets_[0][o].loc;
      mgroot_->Restrict(cbuf, octets_[0][o].src, ngh, ngh, ngh, ngh, ngh, ngh, false);
      for (int v = 0; v < nvar_; ++v)
        mgroot_->SetData(MGVariable::src, v, static_cast<int>(loc.lx3),
                         static_cast<int>(loc.lx2), static_cast<int>(loc.lx1),
                         cbuf(v, ngh, ngh, ngh));
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
#pragma omp parallel for num_threads(nthreads_)
    for (int o = 0; o < noctets_[lev]; ++o) {
      AthenaArray<Real> cbuf;
      cbuf.NewAthenaArray(nvar_, 3, 3, 3);
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
                        octets_[lev][o].src, lev+1, os_, oe_, os_, oe_, os_, oe_, false);
      mgroot_->Restrict(cbuf, octets_[lev][o].def, ngh, ngh, ngh, ngh, ngh, ngh, false);
      for (int v = 0; v < nvar_; ++v)
        octets_[lev-1][oid].src(v, ok, oj, oi) = cbuf(v, ngh, ngh, ngh);
      if (ffas_) {
        mgroot_->Restrict(cbuf, octets_[lev][o].u, ngh, ngh, ngh, ngh, ngh, ngh, false);
        for (int v = 0; v < nvar_; ++v)
          octets_[lev-1][oid].u(v, ok, oj, oi) = cbuf(v, ngh, ngh, ngh);
      }
    }
  } else { // octets to the root grid
#pragma omp parallel for num_threads(nthreads_)
    for (int o = 0; o < noctets_[0]; ++o) {
      AthenaArray<Real> cbuf;
      cbuf.NewAthenaArray(nvar_, 3, 3, 3);
      const LogicalLocation &loc = octets_[0][o].loc;
      int ri = static_cast<int>(loc.lx1);
      int rj = static_cast<int>(loc.lx2);
      int rk = static_cast<int>(loc.lx3);
      mgroot_->CalculateDefect(octets_[0][o].def, octets_[0][o].u,
                               octets_[0][o].src, 1, os_, oe_, os_, oe_, os_, oe_, false);
      mgroot_->Restrict(cbuf, octets_[0][o].def, ngh, ngh, ngh, ngh, ngh, ngh, false);
      for (int v = 0; v < nvar_; ++v)
        mgroot_->SetData(MGVariable::src, v, rk, rj, ri, cbuf(v, ngh, ngh, ngh));
      if (ffas_) {
        mgroot_->Restrict(cbuf, octets_[0][o].u, ngh, ngh, ngh, ngh, ngh, ngh, false);
        for (int v = 0; v < nvar_; ++v)
          mgroot_->SetData(MGVariable::u, v, rk, rj, ri, cbuf(v, ngh, ngh, ngh));
      }
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::ZeroClearOctets()
//  \brief zero clear the data in all the octets

void MultigridDriver::ZeroClearOctets() {
  int maxlevel = current_level_ - 1 - nrootlevel_;
  for (int l = 0; l <= maxlevel; l++) {
#pragma omp parallel for num_threads(nthreads_)
    for (int o = 0; o < noctets_[l]; ++o)
      octets_[l][o].u.ZeroClear();
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::StoreOldDataOctets()
//  \brief store the old u data in the uold array in octets

void MultigridDriver::StoreOldDataOctets() {
  int lev = current_level_ - nrootlevel_;

#pragma omp parallel for num_threads(nthreads_)
  for (int o = 0; o < noctets_[lev]; ++o)
    memcpy(octets_[lev][o].uold.data(), octets_[lev][o].u.data(),
           octets_[lev][o].u.GetSizeInBytes());

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::CalculateFASRHSOctets()
//  \brief Calculate the RHS for FAS in Octets
void MultigridDriver::CalculateFASRHSOctets() {
  int lev = current_level_ - nrootlevel_;

#pragma omp parallel for num_threads(nthreads_)
  for (int o = 0; o < noctets_[lev]; ++o)
    mgroot_->CalculateFASRHS(octets_[lev][o].src, octets_[lev][o].u,
                             lev+1, os_, oe_, os_, oe_, os_, oe_, false);

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SmoothOctets(int color)
//  \brief Apply the smoothing operator on octets
void MultigridDriver::SmoothOctets(int color) {
  int lev = current_level_ - nrootlevel_;

#pragma omp parallel for num_threads(nthreads_)
  for (int o = 0; o < noctets_[lev]; ++o)
    mgroot_->Smooth(octets_[lev][o].u, octets_[lev][o].src,
                    lev+1, os_, oe_, os_, oe_, os_, oe_, color, false);
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
#pragma omp parallel for num_threads(nthreads_)
    for (int o = 0; o < noctets_[0]; ++o) {
      AthenaArray<Real> cbuf;
      cbuf.NewAthenaArray(nvar_, 3, 3, 3);
      const LogicalLocation &loc = octets_[0][o].loc;
      int ri = static_cast<int>(loc.lx1) + ngh - 1;
      int rj = static_cast<int>(loc.lx2) + ngh - 1;
      int rk = static_cast<int>(loc.lx3) + ngh - 1;
      if (ffas_) {
        for (int v = 0; v < nvar_; ++v) {
          for (int k = 0; k <= 2; ++k) {
            for (int j = 0; j <= 2; ++j) {
              for (int i = 0; i <= 2; ++i) {
                cbuf(v,k,j,i) = u(v, rk+k, rj+j, ri+i) - uold(v, rk+k, rj+j, ri+i);
              }
            }
          }
        }
        mgroot_->ProlongateAndCorrect(octets_[0][o].u, cbuf,
                              ngh, ngh, ngh, ngh, ngh, ngh, ngh, ngh, ngh, false);
      } else {
        mgroot_->ProlongateAndCorrect(octets_[0][o].u, u,
                              ri+1, ri+1, rj+1, rj+1, rk+1, rk+1, ngh, ngh, ngh, false);
      }
    }
  } else { // from coarse octets to fine octets
#pragma omp parallel for num_threads(nthreads_)
    for (int o = 0; o < noctets_[flev]; ++o) {
      AthenaArray<Real> cbuf;
      cbuf.NewAthenaArray(nvar_, 3, 3, 3);
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
        for (int v = 0; v < nvar_; ++v) {
          for (int k = 0; k <= 2; ++k) {
            for (int j = 0; j <= 2; ++j) {
              for (int i = 0; i <= 2; ++i)
                cbuf(v,k,j,i) = uc(v, ck+k, cj+j, ci+i) - ucold(v, ck+k, cj+j, ci+i);
            }
          }
        }
        mgroot_->ProlongateAndCorrect(octets_[flev][o].u, cbuf,
                              ngh, ngh, ngh, ngh, ngh, ngh, ngh, ngh, ngh, false);
      } else {
        mgroot_->ProlongateAndCorrect(octets_[flev][o].u, uc,
                              ci+1, ci+1, cj+1, cj+1, ck+1, ck+1, ngh, ngh, ngh, false);
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
#pragma omp parallel for num_threads(nthreads_)
    for (int o = 0; o < noctets_[0]; ++o) {
      const LogicalLocation &loc = octets_[0][o].loc;
      int ri = static_cast<int>(loc.lx1) + ngh;
      int rj = static_cast<int>(loc.lx2) + ngh;
      int rk = static_cast<int>(loc.lx3) + ngh;
      mgroot_->FMGProlongate(octets_[0][o].u, mgroot_->GetCurrentData(),
                             ri, ri, rj, rj, rk, rk, ngh, ngh, ngh, false);
    }
  } else { // from coarse octets to fine octets
#pragma omp parallel for num_threads(nthreads_)
    for (int o = 0; o < noctets_[flev]; ++o) {
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
                             ci, ci, cj, cj, ck, ck, ngh, ngh, ngh, false);
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SetBoundariesOctets(bool fprolong, bool folddata)
//  \brief Apply boundary conditions for octets

void MultigridDriver::SetBoundariesOctets(bool fprolong, bool folddata) {
  int lev = current_level_ - nrootlevel_;

#pragma omp parallel for num_threads(nthreads_) schedule(dynamic,1)
  for (int o = 0; o < noctets_[lev]; ++o) {
    if (fprolong && octets_[lev][o].fleaf == true) continue;
    AthenaArray<Real> cbuf, cbufold;
    AthenaArray<bool> ncoarse;
    cbuf.NewAthenaArray(nvar_, 3, 3, 3);
    cbufold.NewAthenaArray(nvar_, 3, 3, 3);
    ncoarse.NewAthenaArray(3, 3, 3);
    const LogicalLocation &loc = octets_[lev][o].loc;
    AthenaArray<Real> &u = octets_[lev][o].u;
    AthenaArray<Real> &uold = octets_[lev][o].uold;
    LogicalLocation nloc = loc;
    for (int ox3 = -1; ox3 <= 1; ++ox3) {
      nloc.lx3 = loc.lx3 + ox3;
      if (nloc.lx3 < 0) {
        if (mg_mesh_bcs_[BoundaryFace::inner_x3] == BoundaryFlag::periodic)
          nloc.lx3 = (nrbx3_ << lev) - 1;
        else
          continue;
      }
      if (nloc.lx3 >= (nrbx3_ << lev)) {
        if (mg_mesh_bcs_[BoundaryFace::outer_x3] == BoundaryFlag::periodic)
          nloc.lx3 = 0;
        else
          continue;
      }
      for (int ox2 = -1; ox2 <= 1; ++ox2) {
        nloc.lx2 = loc.lx2 + ox2;
        if (nloc.lx2 < 0) {
          if (mg_mesh_bcs_[BoundaryFace::inner_x2] == BoundaryFlag::periodic)
            nloc.lx2 = (nrbx2_ << lev) - 1;
          else
            continue;
        }
        if (nloc.lx2 >= (nrbx2_ << lev)) {
          if (mg_mesh_bcs_[BoundaryFace::outer_x2] == BoundaryFlag::periodic)
            nloc.lx2 = 0;
          else
            continue;
        }
        for (int ox1 = -1; ox1 <= 1; ++ox1) {
          if (ox1 == 0 && ox2 == 0 && ox3 == 0)
            continue;
          // find a neighboring octet - either on the same or coarser level
          nloc.lx1 = loc.lx1 + ox1;
          if (nloc.lx1 < 0) {
            if (mg_mesh_bcs_[BoundaryFace::inner_x1] == BoundaryFlag::periodic)
              nloc.lx1 = (nrbx1_ << lev) - 1;
            else
              continue;
          }
          if (nloc.lx1 >= (nrbx1_ << lev)) {
            if (mg_mesh_bcs_[BoundaryFace::outer_x1] == BoundaryFlag::periodic)
              nloc.lx1 = 0;
            else
              continue;
          }
          if (octetmap_[lev].count(nloc) == 1) { // on the same level
            int nid = octetmap_[lev][nloc];
            const AthenaArray<Real> &un = octets_[lev][nid].u;
            const AthenaArray<Real> &unold = octets_[lev][nid].uold;
            SetOctetBoundarySameLevel(u, un, uold, unold, cbuf, cbufold,
                                      ox1, ox2, ox3, folddata);
          } else if (!fprolong) { // on the coarser level
            // note: prolongation requires neighbors on the same level only
            ncoarse(ox3+1, ox2+1, ox1+1) = true;
            if (lev > 0) { // from octet
              LogicalLocation cloc;
              cloc.lx1 = nloc.lx1 >> 1;
              cloc.lx2 = nloc.lx2 >> 1;
              cloc.lx3 = nloc.lx3 >> 1;
              cloc.level = nloc.level - 1;
              int cid = octetmap_[lev-1][cloc];
              const AthenaArray<Real> &un = octets_[lev-1][cid].u;
              const AthenaArray<Real> &unold = octets_[lev-1][cid].uold;
              SetOctetBoundaryFromCoarser(un, unold, cbuf, cbufold,
                                          loc, ox1, ox2, ox3, false);
            } else { // from root
              const AthenaArray<Real> &un = mgroot_->GetCurrentData();
              const AthenaArray<Real> &unold = mgroot_->GetCurrentOldData();
              SetOctetBoundaryFromCoarser(un, unold, cbuf, cbufold,
                                          nloc, ox1, ox2, ox3, false);
            }
          }
          // note: finer neighbors are not needed here
        }
      }
    }
    if (!fprolong) {
      ApplyPhysicalBoundariesOctet(cbuf, loc, octets_[lev][o].ccoord, true);
      ProlongateOctetBoundariesFluxCons(u, cbuf, ncoarse);
    }
    ApplyPhysicalBoundariesOctet(u, loc, octets_[lev][o].coord, false);
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SetOctetBoundarySameLevel(AthenaArray<Real> &dst,
//   const AthenaArray<Real> &un, AthenaArray<Real> &uold, const AthenaArray<Real> &unold,
//   AthenaArray<Real> &cbuf, AthenaArray<Real> &cbufold,
//   int ox1, int ox2, int ox3, bool folddata)
//  \brief set an Octet boundary from a neighbor Octet on the same level

void MultigridDriver::SetOctetBoundarySameLevel(AthenaArray<Real> &dst,
     const AthenaArray<Real> &un, AthenaArray<Real> &uold, const AthenaArray<Real> &unold,
     AthenaArray<Real> &cbuf, AthenaArray<Real> &cbufold,
     int ox1, int ox2, int ox3, bool folddata) {
  const int ngh = mgroot_->ngh_;
  constexpr Real fac = 0.125;
  const int l = ngh, r = ngh + 1;
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
  int ci = ox1 + ngh, cj = ox2 + ngh, ck = ox3 + ngh;
  for (int v = 0; v < nvar_; ++v) {
    for (int k = ks, nk = nks; k <= ke; ++k, ++nk) {
      for (int j = js, nj = njs; j <= je; ++j, ++nj) {
        for (int i = is, ni = nis; i <= ie; ++i, ++ni)
          dst(v, k, j, i) = un(v, nk, nj, ni);
      }
    }
  }
  for (int v = 0; v < nvar_; ++v)
    cbuf(v,ck,cj,ci) = fac*(un(v,l,l,l) + un(v,l,l,r) + un(v,l,r,l) + un(v,r,l,l)
                          + un(v,r,r,l) + un(v,r,l,r) + un(v,l,r,r) + un(v,r,r,r));
  if (folddata) {
    for (int v = 0; v < nvar_; ++v) {
      for (int k = ks, nk = nks; k <= ke; ++k, ++nk) {
        for (int j = js, nj = njs; j <= je; ++j, ++nj) {
          for (int i = is, ni = nis; i <= ie; ++i, ++ni)
            uold(v, k, j, i) = unold(v, nk, nj, ni);
        }
      }
    }
    for (int v = 0; v < nvar_; ++v)
      cbufold(v,ck,cj,ci) = fac*
             (unold(v,l,l,l)+unold(v,l,l,r) + unold(v,l,r,l) + unold(v,r,l,l)
            + unold(v,r,r,l)+unold(v,r,l,r) + unold(v,l,r,r) + unold(v,r,r,r));
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SetOctetBoundaryFromCoarser(const AthenaArray<Real> &un,
//                            const AthenaArray<Real> &unold, AthenaArray<Real> &cbuf,
//                            AthenaArray<Real> &cbufold, const LogicalLocation &loc,
//                            int ox1, int ox2, int ox3, bool folddata) {
//  \brief set a boundary in the coarse buffer from a neighbor Octet on the coarser level

void MultigridDriver::SetOctetBoundaryFromCoarser(const AthenaArray<Real> &un,
                      const AthenaArray<Real> &unold, AthenaArray<Real> &cbuf,
                      AthenaArray<Real> &cbufold, const LogicalLocation &loc,
                      int ox1, int ox2, int ox3, bool folddata) {
  int ngh = mgroot_->ngh_;
  int ci, cj, ck;
  if (loc.level == locrootlevel_) { // from root
    // given loc is neighbor's location
    ci = static_cast<int>(loc.lx1) + ngh;
    cj = static_cast<int>(loc.lx2) + ngh;
    ck = static_cast<int>(loc.lx3) + ngh;
  } else { // from a neighbor octet
    // given loc is my location
    int ix1 = (static_cast<int>(loc.lx1) & 1);
    int ix2 = (static_cast<int>(loc.lx2) & 1);
    int ix3 = (static_cast<int>(loc.lx3) & 1);
    if (ox1 == 0) ci = ix1 + ngh;
    else          ci = (ix1^1) + ngh;
    if (ox2 == 0) cj = ix2 + ngh;
    else          cj = (ix2^1) + ngh;
    if (ox3 == 0) ck = ix3 + ngh;
    else          ck = (ix3^1) + ngh;
  }
  int i = ngh + ox1, j = ngh + ox2, k = ngh + ox3;
  for (int v = 0; v < nvar_; ++v)
    cbuf(v, k, j, i) = un(v, ck, cj, ci);
  if (folddata) {
    for (int v = 0; v < nvar_; ++v)
      cbufold(v, k, j, i) = unold(v, ck, cj, ci);
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::ApplyPhysicalBoundariesOctet(AthenaArray<Real> &u,
//                    const LogicalLocation &loc, const MGCoordinates &coord, bool fcbuf)
//  \brief Apply physical boundary conditions for an octet

void MultigridDriver::ApplyPhysicalBoundariesOctet(AthenaArray<Real> &u,
              const LogicalLocation &loc, const MGCoordinates &coord, bool fcbuf) {
  int lev = loc.level - locrootlevel_;
  int ngh = mgroot_->ngh_;
  int l, r;
  Real time = pmy_mesh_->time;
  AthenaArray<Real> &mpcoeff = mpcoeff_[0];
  if (fcbuf)
    l = ngh, r = ngh;
  else
    l = ngh, r = ngh+1;

  int bis = l - ngh, bie = r + ngh;
  int bjs = l,       bje = r;
  int bks = l,       bke = r;
  if (loc.lx2 != 0 || mg_mesh_bcs_[BoundaryFace::inner_x2] == BoundaryFlag::periodic)
    bjs = l - ngh;
  if (loc.lx2 != (nrbx2_<<lev)-1
    || mg_mesh_bcs_[BoundaryFace::inner_x2] == BoundaryFlag::periodic)
    bje = r + ngh;
  if (loc.lx3 != 0 || mg_mesh_bcs_[BoundaryFace::inner_x3] == BoundaryFlag::periodic)
    bks = l - ngh;
  if (loc.lx3 != (nrbx3_<<lev)-1
    || mg_mesh_bcs_[BoundaryFace::inner_x3] == BoundaryFlag::periodic)
    bke = r + ngh;

  if (loc.lx1 == 0) {
    switch (mg_mesh_bcs_[BoundaryFace::inner_x1]) {
      case BoundaryFlag::user:
        MGBoundaryFunction_[BoundaryFace::inner_x1](u, time, nvar_,
                                                    l, r, bjs, bje, bks, bke, ngh, coord);
        break;
      case BoundaryFlag::mg_zerograd:
        MGZeroGradientInnerX1(u, time, nvar_, l, r, bjs, bje, bks, bke, ngh, coord);
        break;
      case BoundaryFlag::mg_zerofixed:
        MGZeroFixedInnerX1(u, time, nvar_, l, r, bjs, bje, bks, bke, ngh, coord);
        break;
      case BoundaryFlag::mg_multipole:
        MGMultipoleInnerX1(u, time, nvar_, l, r, bjs, bje, bks, bke, ngh, coord,
                           mpcoeff, mpo_, mporder_);
        break;
      default:
        break;
    }
  }
  if (loc.lx1 == (nrbx1_<<lev)-1) {
    switch (mg_mesh_bcs_[BoundaryFace::outer_x1]) {
      case BoundaryFlag::user:
        MGBoundaryFunction_[BoundaryFace::outer_x1](u, time, nvar_,
                                                    l, r, bjs, bje, bks, bke, ngh, coord);
        break;
      case BoundaryFlag::mg_zerograd:
        MGZeroGradientOuterX1(u, time, nvar_, l, r, bjs, bje, bks, bke, ngh, coord);
        break;
      case BoundaryFlag::mg_zerofixed:
        MGZeroFixedOuterX1(u, time, nvar_, l, r, bjs, bje, bks, bke, ngh, coord);
        break;
      case BoundaryFlag::mg_multipole:
        MGMultipoleOuterX1(u, time, nvar_, l, r, bjs, bje, bks, bke, ngh, coord,
                           mpcoeff, mpo_, mporder_);
        break;
      default:
        break;
    }
  }
  if (loc.lx2 == 0) {
    switch (mg_mesh_bcs_[BoundaryFace::inner_x2]) {
      case BoundaryFlag::user:
        MGBoundaryFunction_[BoundaryFace::inner_x2](u, time, nvar_,
                                                    bis, bie, l, r, bks, bke, ngh, coord);
        break;
      case BoundaryFlag::mg_zerograd:
        MGZeroGradientInnerX2(u, time, nvar_, bis, bie, l, r, bks, bke, ngh, coord);
        break;
      case BoundaryFlag::mg_zerofixed:
        MGZeroFixedInnerX2(u, time, nvar_, bis, bie, l, r, bks, bke, ngh, coord);
        break;
      case BoundaryFlag::mg_multipole:
        MGMultipoleInnerX2(u, time, nvar_, bis, bie, l, r, bks, bke, ngh, coord,
                           mpcoeff, mpo_, mporder_);
        break;
      default:
        break;
    }
  }
  if (loc.lx2 == (nrbx2_<<lev)-1) {
    switch (mg_mesh_bcs_[BoundaryFace::outer_x2]) {
      case BoundaryFlag::user:
        MGBoundaryFunction_[BoundaryFace::outer_x2](u, time, nvar_,
                                                    bis, bie, l, r, bks, bke, ngh, coord);
        break;
      case BoundaryFlag::mg_zerograd:
        MGZeroGradientOuterX2(u, time, nvar_, bis, bie, l, r, bks, bke, ngh, coord);
        break;
      case BoundaryFlag::mg_zerofixed:
        MGZeroFixedOuterX2(u, time, nvar_, bis, bie, l, r, bks, bke, ngh, coord);
        break;
      case BoundaryFlag::mg_multipole:
        MGMultipoleOuterX2(u, time, nvar_, bis, bie, l, r, bks, bke, ngh, coord,
                           mpcoeff, mpo_, mporder_);
        break;
      default:
        break;
    }
  }
  bjs = l - ngh, bje = r + ngh;
  if (loc.lx3 == 0) {
    switch (mg_mesh_bcs_[BoundaryFace::inner_x3]) {
      case BoundaryFlag::user:
        MGBoundaryFunction_[BoundaryFace::inner_x3](u, time, nvar_,
                                                    bis, bie, bjs, bje, l, r, ngh, coord);
        break;
      case BoundaryFlag::mg_zerograd:
        MGZeroGradientInnerX3(u, time, nvar_, bis, bie, bjs, bje, l, r, ngh, coord);
        break;
      case BoundaryFlag::mg_zerofixed:
        MGZeroFixedInnerX3(u, time, nvar_, bis, bie, bjs, bje, l, r, ngh, coord);
        break;
      case BoundaryFlag::mg_multipole:
        MGMultipoleInnerX3(u, time, nvar_, bis, bie, bjs, bje, l, r, ngh, coord,
                           mpcoeff, mpo_, mporder_);
        break;
      default:
        break;
    }
  }
  if (loc.lx3 == (nrbx3_<<lev)-1) {
    switch (mg_mesh_bcs_[BoundaryFace::outer_x3]) {
      case BoundaryFlag::user:
        MGBoundaryFunction_[BoundaryFace::outer_x3](u, time, nvar_,
                                                    bis, bie, bjs, bje, l, r, ngh, coord);
        break;
      case BoundaryFlag::mg_zerograd:
        MGZeroGradientOuterX3(u, time, nvar_, bis, bie, bjs, bje, l, r, ngh, coord);
        break;
      case BoundaryFlag::mg_zerofixed:
        MGZeroFixedOuterX3(u, time, nvar_, bis, bie, bjs, bje, l, r, ngh, coord);
        break;
      case BoundaryFlag::mg_multipole:
        MGMultipoleOuterX3(u, time, nvar_, bis, bie, bjs, bje, l, r, ngh, coord,
                           mpcoeff, mpo_, mporder_);
        break;
      default:
        break;
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::RestrictOctetsBeforeTransfer()
//  \brief Restrict all the octets

void MultigridDriver::RestrictOctetsBeforeTransfer() {
  const int &ngh = mgroot_->ngh_;
  for (int l = nreflevel_ - 1; l >= 1; --l) {  // fine octets to coarse octets
#pragma omp parallel for num_threads(nthreads_)
    for (int o = 0; o < noctets_[l]; ++o) {
      AthenaArray<Real> cbuf;
      cbuf.NewAthenaArray(nvar_, 3, 3, 3);
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
      mgroot_->Restrict(cbuf, octets_[l][o].u, ngh, ngh, ngh, ngh, ngh, ngh, false);
      for (int v = 0; v < nvar_; ++v)
        octets_[l-1][oid].u(v, ok, oj, oi) = cbuf(v, ngh, ngh, ngh);
    }
  }
#pragma omp parallel for num_threads(nthreads_)
  for (int o = 0; o < noctets_[0]; ++o) { // octets to the root grid
    AthenaArray<Real> cbuf;
    cbuf.NewAthenaArray(nvar_, 3, 3, 3);
    const LogicalLocation &loc = octets_[0][o].loc;
    mgroot_->Restrict(cbuf, octets_[0][o].u, ngh, ngh, ngh, ngh, ngh, ngh, false);
    for (int v = 0; v < nvar_; ++v)
      mgroot_->SetData(MGVariable::u, v, static_cast<int>(loc.lx3),
                       static_cast<int>(loc.lx2), static_cast<int>(loc.lx1),
                       cbuf(v, ngh, ngh, ngh));
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SetOctetBoundariesBeforeTransfer(bool folddata)
//  \brief Set octet boundaries before transfer from root to blocks

void MultigridDriver::SetOctetBoundariesBeforeTransfer(bool folddata) {
  const int ngh = mgroot_->ngh_;
  Real time = pmy_mesh_->time;

  // clear octet boundary flag
  for (int l = 0; l < nreflevel_; ++l) {
    for (int o = 0; o < noctets_[l]; ++o)
      octetbflag_[l][o] = false;
  }

#pragma omp parallel for num_threads(nthreads_)
  for (auto itr = vmg_.begin(); itr < vmg_.end(); itr++) {
    Multigrid *pmg = *itr;
    LogicalLocation loc = pmg->loc_;
    if (loc.level == locrootlevel_) continue;
    AthenaArray<Real> cbuf, cbufold;
    AthenaArray<bool> ncoarse;
    cbuf.NewAthenaArray(nvar_, 3, 3, 3);
    cbufold.NewAthenaArray(nvar_, 3, 3, 3);
    ncoarse.NewAthenaArray(3, 3, 3);
    loc.lx1 = loc.lx1 >> 1;
    loc.lx2 = loc.lx2 >> 1;
    loc.lx3 = loc.lx3 >> 1;
    loc.level = loc.level - 1;
    int lev = loc.level - locrootlevel_;
    int oid = octetmap_[lev][loc];
    if (octetbflag_[lev][oid] == true) continue;
    octetbflag_[lev][oid] = true;
    AthenaArray<Real> &u = octets_[lev][oid].u;
    AthenaArray<Real> &uold = octets_[lev][oid].uold;
    LogicalLocation nloc = loc;
    for (int ox3 = -1; ox3 <= 1; ++ox3) {
      nloc.lx3 = loc.lx3 + ox3;
      if (nloc.lx3 < 0) {
        if (mg_mesh_bcs_[BoundaryFace::inner_x3] == BoundaryFlag::periodic)
          nloc.lx3 = (nrbx3_ << lev) - 1;
        else
          continue;
      }
      if (nloc.lx3 >= (nrbx3_ << lev)) {
        if (mg_mesh_bcs_[BoundaryFace::outer_x3] == BoundaryFlag::periodic)
          nloc.lx3 = 0;
        else
          continue;
      }
      for (int ox2 = -1; ox2 <= 1; ++ox2) {
        nloc.lx2 = loc.lx2 + ox2;
        if (nloc.lx2 < 0) {
          if (mg_mesh_bcs_[BoundaryFace::inner_x2] == BoundaryFlag::periodic)
            nloc.lx2 = (nrbx2_ << lev) - 1;
          else
            continue;
        }
        if (nloc.lx2 >= (nrbx2_ << lev)) {
          if (mg_mesh_bcs_[BoundaryFace::outer_x2] == BoundaryFlag::periodic)
            nloc.lx2 = 0;
          else
            continue;
        }
        for (int ox1 = -1; ox1 <= 1; ++ox1) {
          if (ox1 == 0 && ox2 == 0 && ox3 == 0) continue;
          nloc.lx1 = loc.lx1 + ox1;
          if (nloc.lx1 < 0) {
            if (mg_mesh_bcs_[BoundaryFace::inner_x1] == BoundaryFlag::periodic)
              nloc.lx1 = (nrbx1_ << lev) - 1;
            else
              continue;
          }
          if (nloc.lx1 >= (nrbx1_ << lev)) {
            if (mg_mesh_bcs_[BoundaryFace::outer_x1] == BoundaryFlag::periodic)
              nloc.lx1 = 0;
            else
              continue;
          }
          if (octetmap_[lev].count(nloc) == 1) { // same or finer
            int nid = octetmap_[lev][nloc];
            const AthenaArray<Real> &un = octets_[lev][nid].u;
            const AthenaArray<Real> &unold = octets_[lev][nid].uold;
            SetOctetBoundarySameLevel(u, un, uold, unold, cbuf, cbufold,
                                      ox1, ox2, ox3, folddata);
          } else { // coarser
            ncoarse(ox3+1, ox2+1, ox1+1) = true;
            if (lev > 0) { // from octet
              LogicalLocation cloc;
              cloc.lx1 = nloc.lx1 >> 1;
              cloc.lx2 = nloc.lx2 >> 1;
              cloc.lx3 = nloc.lx3 >> 1;
              cloc.level = nloc.level - 1;
              int cid = octetmap_[lev-1][cloc];
              const AthenaArray<Real> &un = octets_[lev-1][cid].u;
              const AthenaArray<Real> &unold = octets_[lev-1][cid].uold;
              SetOctetBoundaryFromCoarser(un, unold, cbuf, cbufold,
                                          loc, ox1, ox2, ox3, folddata);
            } else { // from root
              const AthenaArray<Real> &un = mgroot_->GetCurrentData();
              const AthenaArray<Real> &unold = mgroot_->GetCurrentOldData();
              SetOctetBoundaryFromCoarser(un, unold, cbuf, cbufold,
                                          nloc, ox1, ox2, ox3, folddata);
            }
          }
        }
      }
    }

    ApplyPhysicalBoundariesOctet(cbuf, loc, octets_[lev][oid].ccoord, true);
    if (folddata)
      ApplyPhysicalBoundariesOctet(cbufold, loc, octets_[lev][oid].ccoord, true);
    ProlongateOctetBoundaries(u, uold, cbuf, cbufold, ncoarse, folddata);
    ApplyPhysicalBoundariesOctet(u, loc, octets_[lev][oid].coord, false);
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::ProlongateOctetBoundaries(AthenaArray<Real> &u,
//           AthenaArray<Real> &uold, AthenaArray<Real> &cbuf, AthenaArray<Real> &cbufold,
//           const AthenaArray<bool> &ncoarse,bool folddata)
//  \brief Prolongate octet boundaries contacting the coarser level

void MultigridDriver::ProlongateOctetBoundaries(AthenaArray<Real> &u,
     AthenaArray<Real> &uold, AthenaArray<Real> &cbuf, AthenaArray<Real> &cbufold,
     const AthenaArray<bool> &ncoarse, bool folddata) {
  const int ngh = mgroot_->ngh_;
  const int flim = 2 + ngh;
  constexpr Real fac = 0.125;
  const int l = ngh, r = ngh + 1;

  for (int v = 0; v < nvar_; ++v)
    cbuf(v,ngh,ngh,ngh) = fac*(u(v,l,l,l)+u(v,l,l,r)+u(v,l,r,l)+u(v,r,l,l)
                              +u(v,r,r,l)+u(v,r,l,r)+u(v,l,r,r)+u(v,r,r,r));
  if (folddata) {
    for (int v = 0; v < nvar_; ++v)
      cbufold(v,ngh,ngh,ngh) = fac*
                             (uold(v,l,l,l)+uold(v,l,l,r)+uold(v,l,r,l)+uold(v,r,l,l)
                             +uold(v,r,r,l)+uold(v,r,l,r)+uold(v,l,r,r)+uold(v,r,r,r));
  }

  const AthenaArray<Real> &cb = cbuf;
  const AthenaArray<Real> &co = cbufold;
  for (int ox3 = -1; ox3 <= 1; ++ox3) {
    for (int ox2 = -1; ox2 <= 1; ++ox2) {
      for (int ox1 = -1; ox1 <= 1; ++ox1) {
        if (ncoarse(ox3+1, ox2+1, ox1+1)) { // coarser
          int i = ox1 + ngh, j = ox2 + ngh, k = ox3 + ngh;
          int fi = ox1*2 + ngh, fj = ox2*2 + ngh, fk = ox3*2 + ngh;
          for (int v = 0; v < nvar_; ++v) {
            if (fk >= 0 && fj >= 0 && fi >= 0)
              u(v, fk,   fj,   fi  ) =
                0.015625*(27.0*cb(v,k,j,i)+cb(v,k-1,j-1,i-1)
                  +9.0*(cb(v,k,j,i-1)+cb(v,k,j-1,i)+cb(v,k-1,j,i))
                  +3.0*(cb(v,k-1,j-1,i)+cb(v,k-1,j,i-1)+cb(v,k,j-1,i-1)));
            if (fk >= 0 && fj >= 0 && fi < flim)
              u(v, fk,   fj,   fi+1) =
                0.015625*(27.0*cb(v,k,j,i)+cb(v,k-1,j-1,i+1)
                  +9.0*(cb(v,k,j,i+1)+cb(v,k,j-1,i)+cb(v,k-1,j,i))
                  +3.0*(cb(v,k-1,j-1,i)+cb(v,k-1,j,i+1)+cb(v,k,j-1,i+1)));
            if (fk >= 0 && fj < flim && fi >= 0)
              u(v, fk,   fj+1, fi  ) =
                0.015625*(27.0*cb(v,k,j,i)+cb(v,k-1,j+1,i-1)
                  +9.0*(cb(v,k,j,i-1)+cb(v,k,j+1,i)+cb(v,k-1,j,i))
                  +3.0*(cb(v,k-1,j+1,i)+cb(v,k-1,j,i-1)+cb(v,k,j+1,i-1)));
            if (fk < flim && fj >= 0 && fi >= 0)
              u(v, fk+1, fj,   fi  ) =
                0.015625*(27.0*cb(v,k,j,i)+cb(v,k+1,j-1,i-1)
                  +9.0*(cb(v,k,j,i-1)+cb(v,k,j-1,i)+cb(v,k+1,j,i))
                  +3.0*(cb(v,k+1,j-1,i)+cb(v,k+1,j,i-1)+cb(v,k,j-1,i-1)));
            if (fk < flim && fj < flim && fi >= 0)
              u(v, fk+1, fj+1, fi  ) =
                0.015625*(27.0*cb(v,k,j,i)+cb(v,k+1,j+1,i-1)
                  +9.0*(cb(v,k,j,i-1)+cb(v,k,j+1,i)+cb(v,k+1,j,i))
                  +3.0*(cb(v,k+1,j+1,i)+cb(v,k+1,j,i-1)+cb(v,k,j+1,i-1)));
            if (fk < flim && fj >= 0 && fi < flim)
              u(v, fk+1, fj,   fi+1) =
                0.015625*(27.0*cb(v,k,j,i)+cb(v,k+1,j-1,i+1)
                  +9.0*(cb(v,k,j,i+1)+cb(v,k,j-1,i)+cb(v,k+1,j,i))
                  +3.0*(cb(v,k+1,j-1,i)+cb(v,k+1,j,i+1)+cb(v,k,j-1,i+1)));
            if (fk >= 0 && fj < flim && fi < flim)
              u(v, fk,  fj+1, fi+1) =
                0.015625*(27.0*cb(v,k,j,i)+cb(v,k-1,j+1,i+1)
                  +9.0*(cb(v,k,j,i+1)+cb(v,k,j+1,i)+cb(v,k-1,j,i))
                  +3.0*(cb(v,k-1,j+1,i)+cb(v,k-1,j,i+1)+cb(v,k,j+1,i+1)));
            if (fk < flim && fj < flim && fi < flim)
              u(v, fk+1, fj+1, fi+1) =
                0.015625*(27.0*cb(v,k,j,i)+cb(v,k+1,j+1,i+1)
                  +9.0*(cb(v,k,j,i+1)+cb(v,k,j+1,i)+cb(v,k+1,j,i))
                  +3.0*(cb(v,k+1,j+1,i)+cb(v,k+1,j,i+1)+cb(v,k,j+1,i+1)));
          }
          if (folddata) {
            for (int v = 0; v < nvar_; ++v) {
              if (fk >= 0 && fj >= 0 && fi >= 0)
                uold(v, fk,   fj,   fi  ) =
                  0.015625*(27.0*co(v,k,j,i)+co(v,k-1,j-1,i-1)
                    +9.0*(co(v,k,j,i-1)+co(v,k,j-1,i)+co(v,k-1,j,i))
                    +3.0*(co(v,k-1,j-1,i)+co(v,k-1,j,i-1)+co(v,k,j-1,i-1)));
              if (fk >= 0 && fj >= 0 && fi < flim)
                uold(v, fk,   fj,   fi+1) =
                  0.015625*(27.0*co(v,k,j,i)+co(v,k-1,j-1,i+1)
                    +9.0*(co(v,k,j,i+1)+co(v,k,j-1,i)+co(v,k-1,j,i))
                    +3.0*(co(v,k-1,j-1,i)+co(v,k-1,j,i+1)+co(v,k,j-1,i+1)));
              if (fk >= 0 && fj < flim && fi >= 0)
                uold(v, fk,   fj+1, fi  ) =
                  0.015625*(27.0*co(v,k,j,i)+co(v,k-1,j+1,i-1)
                    +9.0*(co(v,k,j,i-1)+co(v,k,j+1,i)+co(v,k-1,j,i))
                    +3.0*(co(v,k-1,j+1,i)+co(v,k-1,j,i-1)+co(v,k,j+1,i-1)));
              if (fk < flim && fj >= 0 && fi >= 0)
                uold(v, fk+1, fj,   fi  ) =
                  0.015625*(27.0*co(v,k,j,i)+co(v,k+1,j-1,i-1)
                    +9.0*(co(v,k,j,i-1)+co(v,k,j-1,i)+co(v,k+1,j,i))
                    +3.0*(co(v,k+1,j-1,i)+co(v,k+1,j,i-1)+co(v,k,j-1,i-1)));
              if (fk < flim && fj < flim && fi >= 0)
                uold(v, fk+1, fj+1, fi  ) =
                  0.015625*(27.0*co(v,k,j,i)+co(v,k+1,j+1,i-1)
                    +9.0*(co(v,k,j,i-1)+co(v,k,j+1,i)+co(v,k+1,j,i))
                    +3.0*(co(v,k+1,j+1,i)+co(v,k+1,j,i-1)+co(v,k,j+1,i-1)));
              if (fk < flim && fj >= 0 && fi < flim)
                uold(v, fk+1, fj,   fi+1) =
                  0.015625*(27.0*co(v,k,j,i)+co(v,k+1,j-1,i+1)
                    +9.0*(co(v,k,j,i+1)+co(v,k,j-1,i)+co(v,k+1,j,i))
                    +3.0*(co(v,k+1,j-1,i)+co(v,k+1,j,i+1)+co(v,k,j-1,i+1)));
              if (fk >= 0 && fj < flim && fi < flim)
                uold(v, fk,  fj+1, fi+1) =
                  0.015625*(27.0*co(v,k,j,i)+co(v,k-1,j+1,i+1)
                    +9.0*(co(v,k,j,i+1)+co(v,k,j+1,i)+co(v,k-1,j,i))
                    +3.0*(co(v,k-1,j+1,i)+co(v,k-1,j,i+1)+co(v,k,j+1,i+1)));
              if (fk < flim && fj < flim && fi < flim)
                uold(v, fk+1, fj+1, fi+1) =
                  0.015625*(27.0*co(v,k,j,i)+co(v,k+1,j+1,i+1)
                    +9.0*(co(v,k,j,i+1)+co(v,k,j+1,i)+co(v,k+1,j,i))
                    +3.0*(co(v,k+1,j+1,i)+co(v,k+1,j,i+1)+co(v,k,j+1,i+1)));
            }
          }
        }
      }
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::AllocateMultipoleCoefficients()
//  \brief Allocate arrays for multipole expansion coefficients

void MultigridDriver::AllocateMultipoleCoefficients() {
  nmpcoeff_ = 0;
  if (mporder_ <= 0) return;
  for (int i = 0; i <= mporder_; ++i)
    nmpcoeff_ += 2 * i + 1;

  mpcoeff_ = new AthenaArray<Real>[nthreads_];
  for (int th = 0; th < nthreads_; ++th)
    mpcoeff_[th].NewAthenaArray(nvar_, nmpcoeff_);

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::CalculateMultipoleCoefficients()
//  \brief Calculate multipole expansion coefficients

void MultigridDriver::CalculateMultipoleCoefficients() {
  for (int th = 0; th < nthreads_; ++th)
    mpcoeff_[th].ZeroClear();
#pragma omp parallel for num_threads(nthreads_)
  for (auto itr = vmg_.begin(); itr < vmg_.end(); itr++) {
    Multigrid *pmg = *itr;
    int th = 0;
#ifdef OPENMP_PARALLEL
    th = omp_get_thread_num();
#endif
    pmg->CalculateMultipoleCoefficients(mpcoeff_[th]);
  }
  for (int th = 1; th < nthreads_; ++th) {
    for (int i = 0; i < nmpcoeff_; ++i)
      mpcoeff_[0](i) += mpcoeff_[th](i);
  }
#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, mpcoeff_[0].data(), nmpcoeff_, MPI_ATHENA_REAL,
                MPI_SUM, MPI_COMM_MULTIGRID);
#endif
  ScaleMultipoleCoefficients();
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::ScaleMultipoleCoefficients()
//  \brief Scale coefficients for multipole expansion

void MultigridDriver::ScaleMultipoleCoefficients() {
  AthenaArray<Real> &mpcoeff = mpcoeff_[0];
  // constants for multipole expansion
  constexpr Real c0  = -0.25/PI;
  constexpr Real c1  = -0.25/PI;
  constexpr Real c2  = -0.0625/PI;
  constexpr Real c2a = -0.75/PI;
  constexpr Real c30 = -0.0625/PI;
  constexpr Real c31 = -0.0625*1.5/PI;
  constexpr Real c32 = -0.25*15.0/PI;
  constexpr Real c33 = -0.0625*2.5/PI;
  constexpr Real c40 = -0.0625*0.0625/PI;
  constexpr Real c41 = -0.0625*2.5/PI;
  constexpr Real c42 = -0.0625*5.0/PI;
  constexpr Real c43 = -0.0625*17.5/PI;
  constexpr Real c44 = -0.25*35.0/PI;

  mpcoeff(0) *= c0;
  mpcoeff(1) *= c1;
  mpcoeff(2) *= c1;
  mpcoeff(3) *= c1;
  mpcoeff(4) *= c2a;
  mpcoeff(5) *= c2a;
  mpcoeff(6) *= c2;
  mpcoeff(7) *= c2a;
  mpcoeff(8) *= c2a;
  if (mporder_ == 4) {
    mpcoeff(9)  *= c33;
    mpcoeff(10) *= c32;
    mpcoeff(11) *= c31;
    mpcoeff(12) *= c30;
    mpcoeff(13) *= c31;
    mpcoeff(14) *= c32;
    mpcoeff(15) *= c33;
    mpcoeff(16) *= c44;
    mpcoeff(17) *= c43;
    mpcoeff(18) *= c42;
    mpcoeff(19) *= c41;
    mpcoeff(20) *= c40;
    mpcoeff(21) *= c41;
    mpcoeff(22) *= c42;
    mpcoeff(23) *= c43;
    mpcoeff(24) *= c44;
  }
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::CalculateCenterOfMass()
//  \brief Calculate the position of the center of mass for multipole expansion

void MultigridDriver::CalculateCenterOfMass() {
  for (int th = 0; th < nthreads_; ++th) {
    for (int i = 0; i < 4; ++i)
      mpcoeff_[th](i) = 0.0;
  }
#pragma omp parallel for num_threads(nthreads_)
  for (auto itr = vmg_.begin(); itr < vmg_.end(); itr++) {
    Multigrid *pmg = *itr;
    int th = 0;
#ifdef OPENMP_PARALLEL
    th = omp_get_thread_num();
#endif
    pmg->CalculateCenterOfMass(mpcoeff_[th]);
  }
  for (int th = 1; th < nthreads_; ++th) {
    for (int i = 0; i < 4; ++i)
      mpcoeff_[0](i) += mpcoeff_[th](i);
  }
#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, mpcoeff_[0].data(), 4, MPI_ATHENA_REAL,
                MPI_SUM, MPI_COMM_MULTIGRID);
#endif

  Real im = 1.0 / mpcoeff_[0](0);
  mpo_(0) = im  * mpcoeff_[0](3);
  mpo_(1) = im  * mpcoeff_[0](1);
  mpo_(2) = im  * mpcoeff_[0](2);

  return;
}
