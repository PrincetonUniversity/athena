//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_var.cpp
//! \brief constructor/destructor and default implementations for some functions in the
//!        abstract BoundaryVariable class

// C headers

// C++ headers
#include <cstring>    // std::memcpy
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../globals.hpp"
#include "../mesh/mesh.hpp"
#include "bvals_interfaces.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

//! constructor

BoundaryVariable::BoundaryVariable(MeshBlock *pmb, bool fflux) :
                  bvar_index(), pmy_block_(pmb), pmy_mesh_(pmb->pmy_mesh),
                  pbval_(pmb->pbval), fflux_(fflux) {}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryVariable::InitBoundaryData(BoundaryData<> &bd, BoundaryQuantity type)
//! \brief Initialize BoundaryData structure

void BoundaryVariable::InitBoundaryData(BoundaryData<> &bd, BoundaryQuantity type) {
  MeshBlock *pmb = pmy_block_;
  NeighborIndexes *ni = pbval_->ni;
  int cng = pmb->cnghost;
  int size = 0;

  bd.nbmax = pbval_->maxneighbor_;
  // KGF: what is happening in the next two conditionals??
  // they are preventing the elimination of "BoundaryQuantity type" function parameter in
  // favor of a simpler boolean switch
  if (type == BoundaryQuantity::cc_flcor || type == BoundaryQuantity::fc_flcor) {
    for (bd.nbmax = 0; pbval_->ni[bd.nbmax].type == NeighborConnect::face; bd.nbmax++) {}
  }
  if (type == BoundaryQuantity::fc_flcor) {
    for (          ; pbval_->ni[bd.nbmax].type == NeighborConnect::edge; bd.nbmax++) {}
  }
  for (int n=0; n<bd.nbmax; n++) {
    // Clear flags and requests
    bd.flag[n] = BoundaryStatus::waiting;
    bd.sflag[n] = BoundaryStatus::waiting;
    bd.send[n] = nullptr;
    bd.recv[n] = nullptr;
#ifdef MPI_PARALLEL
    bd.req_send[n] = MPI_REQUEST_NULL;
    bd.req_recv[n] = MPI_REQUEST_NULL;
#endif
    // Allocate buffers, calculating the buffer size (variable vs. flux correction)
    if (type == BoundaryQuantity::cc || type == BoundaryQuantity::fc) {
      size = this->ComputeVariableBufferSize(ni[n], cng);
    } else if (type == BoundaryQuantity::cc_flcor || type == BoundaryQuantity::fc_flcor) {
      size = this->ComputeFluxCorrectionBufferSize(ni[n], cng);
    } else {
      std::stringstream msg;
      msg << "### FATAL ERROR in InitBoundaryData" << std::endl
          << "Invalid boundary type is specified." << std::endl;
      ATHENA_ERROR(msg);
    }
    bd.send[n] = new Real[size];
    bd.recv[n] = new Real[size];
  }
}


//----------------------------------------------------------------------------------------
//! \fn void BoundaryVariable::DestroyBoundaryData(BoundaryData<> &bd)
//! \brief Destroy BoundaryData structure

void BoundaryVariable::DestroyBoundaryData(BoundaryData<> &bd) {
  for (int n=0; n<bd.nbmax; n++) {
    delete [] bd.send[n];
    delete [] bd.recv[n];
#ifdef MPI_PARALLEL
    if (bd.req_send[n] != MPI_REQUEST_NULL)
      MPI_Request_free(&bd.req_send[n]);
    if (bd.req_recv[n] != MPI_REQUEST_NULL)
      MPI_Request_free(&bd.req_recv[n]);
#endif
  }
}


//----------------------------------------------------------------------------------------
//! \fn void BoundaryVariable::CopyVariableBufferSameProcess(NeighborBlock& nb, int ssize)
//! \brief Called in BoundaryVariable::SendBoundaryBuffer() and SendFluxCorrection()
//! when the destination neighbor block is on the same MPI rank as the sending MeshBlcok.
//! So std::memcpy() call requires a pointer to "void *dst" corresponding to
//! bd_var_.recv[nb.targetid] on the target block

void BoundaryVariable::CopyVariableBufferSameProcess(NeighborBlock& nb, int ssize) {
  // Locate target buffer
  // 1) which MeshBlock?
  MeshBlock *ptarget_block = pmy_mesh_->FindMeshBlock(nb.snb.gid);
  // 2) which element in vector of BoundaryVariable *?
  BoundaryData<> *ptarget_bdata = &(ptarget_block->pbval->bvars[bvar_index]->bd_var_);
  std::memcpy(ptarget_bdata->recv[nb.targetid], bd_var_.send[nb.bufid],
              ssize*sizeof(Real));
  // finally, set the BoundaryStatus flag on the destination buffer
  ptarget_bdata->flag[nb.targetid] = BoundaryStatus::arrived;
  return;
}

// KGF: change ssize to send_count


//----------------------------------------------------------------------------------------
//! \fn void BoundaryVariable::CopyFluxCorrectionBufferSameProcess(NeighborBlock& nb,
//!                                                                int ssize)
//!  \brief Same as CopyVariableBufferSameProcess but for flux correction

void BoundaryVariable::CopyFluxCorrectionBufferSameProcess(NeighborBlock& nb, int ssize) {
  // Locate target buffer
  // 1) which MeshBlock?
  MeshBlock *ptarget_block = pmy_mesh_->FindMeshBlock(nb.snb.gid);
  // 2) which element in vector of BoundaryVariable *?
  BoundaryData<> *ptarget_bdata =
      &(ptarget_block->pbval->bvars[bvar_index]->bd_var_flcor_);
  std::memcpy(ptarget_bdata->recv[nb.targetid], bd_var_flcor_.send[nb.bufid],
              ssize*sizeof(Real));
  // finally, set the BoundaryStatus flag on the destination buffer
  ptarget_bdata->flag[nb.targetid] = BoundaryStatus::arrived;
  return;
}


// no nb.targetid, nb.bufid in SimpleNeighborBlock.
// fixed "int bufid" is used for both IDs. Seems unnecessarily strict.

//----------------------------------------------------------------------------------------
//! \fn void BoundaryVariable::CopyShearBufferSameProcess(SimpleNeighborBlock& snb,
//!                                               int ssize, int bufid, bool upper)
//! \brief Same as CopyVariableBufferSameProcess but for shear boundaries

void BoundaryVariable::CopyShearBufferSameProcess(SimpleNeighborBlock& snb, int ssize,
                                                  int bufid, bool upper) {
  // Locate target buffer
  // 1) which MeshBlock?
  MeshBlock *ptarget_block = pmy_mesh_->FindMeshBlock(snb.gid);
  // 2) which element in vector of BoundaryVariable *?
  ShearingBoundaryData *ptarget_bdata =
      &(ptarget_block->pbval->bvars[bvar_index]->shear_bd_var_[upper]);
  std::memcpy(ptarget_bdata->recv[bufid], shear_bd_var_[upper].send[bufid],
              ssize*sizeof(Real));
  // finally, set the BoundaryStatus flag on the destination buffer
  ptarget_bdata->flag[bufid] = BoundaryStatus::arrived;
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void BoundaryVariable::CopyShearFluxSameProcess(SimpleNeighborBlock& snb,
//!                                                     int ssize, int bufid, bool upper)
//! \brief Same as CopyVariableBufferSameProcess but for shear flux

void BoundaryVariable::CopyShearFluxSameProcess(SimpleNeighborBlock& snb, int ssize,
                                               int bufid, bool upper) {
  // Locate target buffer
  // 1) which MeshBlock?
  MeshBlock *ptarget_block = pmy_mesh_->FindMeshBlock(snb.gid);
  // 2) which element in vector of BoundaryVariable *?
  ShearingFluxBoundaryData *ptarget_bdata =
      &(ptarget_block->pbval->bvars[bvar_index]->shear_bd_flux_[upper]);
  std::memcpy(ptarget_bdata->recv[bufid], shear_bd_flux_[upper].send[bufid],
              ssize*sizeof(Real));
  // finally, set the BoundaryStatus flag on the destination buffer
  ptarget_bdata->flag[bufid] = BoundaryStatus::arrived;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryVariable::SetCompletedFlagSameProcess(NeighborBlock& nb)
//! \brief
//!
//!  Called in CellCenteredBoundaryVariable::SendFluxCorrection() when there is no
//!  need to send any information to nb on the same process. Just set
//!  BoundaryStatus::completed in the flag.
void BoundaryVariable::SetCompletedFlagSameProcess(NeighborBlock& nb) {
  // Locate target buffer
  // 1) which MeshBlock?
  MeshBlock *ptarget_block = pmy_mesh_->FindMeshBlock(nb.snb.gid);
  // 2) which element in vector of BoundaryVariable *?
  BoundaryData<> *ptarget_bdata =
      &(ptarget_block->pbval->bvars[bvar_index]->bd_var_flcor_);
  ptarget_bdata->flag[nb.targetid] = BoundaryStatus::completed;
  return;
}

// Default / shared implementations of 4x BoundaryBuffer public functions

//----------------------------------------------------------------------------------------
//! \fn void BoundaryVariable::SendBoundaryBuffers()
//! \brief Send boundary buffers of variables

void BoundaryVariable::SendBoundaryBuffers() {
  MeshBlock *pmb = pmy_block_;
  int mylevel = pmb->loc.level;
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (bd_var_.sflag[nb.bufid] == BoundaryStatus::completed) continue;
    int ssize;
    if (nb.snb.level == mylevel)
      ssize = LoadBoundaryBufferSameLevel(bd_var_.send[nb.bufid], nb);
    else if (nb.snb.level<mylevel)
      ssize = LoadBoundaryBufferToCoarser(bd_var_.send[nb.bufid], nb);
    else
      ssize = LoadBoundaryBufferToFiner(bd_var_.send[nb.bufid], nb);
    if (nb.snb.rank == Globals::my_rank) {  // on the same process
      CopyVariableBufferSameProcess(nb, ssize);
    }
#ifdef MPI_PARALLEL
    else  // MPI
      MPI_Start(&(bd_var_.req_send[nb.bufid]));
#endif
    bd_var_.sflag[nb.bufid] = BoundaryStatus::completed;
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn bool BoundaryVariable::ReceiveBoundaryBuffers()
//! \brief receive the boundary data

bool BoundaryVariable::ReceiveBoundaryBuffers() {
  bool bflag = true;

  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (bd_var_.flag[nb.bufid] == BoundaryStatus::arrived) continue;
    if (bd_var_.flag[nb.bufid] == BoundaryStatus::waiting) {
      if (nb.snb.rank == Globals::my_rank) {  // on the same process
        bflag = false;
        continue;
      }
#ifdef MPI_PARALLEL
      else { // NOLINT // MPI boundary
        int test;
        // probe MPI communications.  This is a bit of black magic that seems to promote
        // communications to top of stack and gets them to complete more quickly
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &test, MPI_STATUS_IGNORE);
        MPI_Test(&(bd_var_.req_recv[nb.bufid]), &test, MPI_STATUS_IGNORE);
        if (!static_cast<bool>(test)) {
          bflag = false;
          continue;
        }
        bd_var_.flag[nb.bufid] = BoundaryStatus::arrived;
      }
#endif
    }
  }
  return bflag;
}


//----------------------------------------------------------------------------------------
//! \fn void BoundaryVariable::SetBoundaries()
//! \brief set the boundary data

void BoundaryVariable::SetBoundaries() {
  MeshBlock *pmb = pmy_block_;
  int mylevel = pmb->loc.level;
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (nb.snb.level == mylevel)
      SetBoundarySameLevel(bd_var_.recv[nb.bufid], nb);
    else if (nb.snb.level < mylevel) // only sets the prolongation buffer
      SetBoundaryFromCoarser(bd_var_.recv[nb.bufid], nb);
    else
      SetBoundaryFromFiner(bd_var_.recv[nb.bufid], nb);
    bd_var_.flag[nb.bufid] = BoundaryStatus::completed; // completed
  }

  if (pbval_->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar ||
      pbval_->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar)
    PolarBoundarySingleAzimuthalBlock();

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void BoundaryVariable::ReceiveAndSetBoundariesWithWait()
//! \brief receive and set the boundary data for initialization

void BoundaryVariable::ReceiveAndSetBoundariesWithWait() {
  MeshBlock *pmb = pmy_block_;
  int mylevel = pmb->loc.level;
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
#ifdef MPI_PARALLEL
    if (nb.snb.rank != Globals::my_rank)
      MPI_Wait(&(bd_var_.req_recv[nb.bufid]),MPI_STATUS_IGNORE);
#endif
    if (nb.snb.level == mylevel)
      SetBoundarySameLevel(bd_var_.recv[nb.bufid], nb);
    else if (nb.snb.level < mylevel)
      SetBoundaryFromCoarser(bd_var_.recv[nb.bufid], nb);
    else
      SetBoundaryFromFiner(bd_var_.recv[nb.bufid], nb);
    bd_var_.flag[nb.bufid] = BoundaryStatus::completed; // completed
  }

  if (pbval_->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar
      || pbval_->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar)
    PolarBoundarySingleAzimuthalBlock();

  return;
}

//PolarFieldBoundaryAverage();
