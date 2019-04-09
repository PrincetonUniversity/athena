//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_var.cpp
//  \brief constructor/destructor and default implementations for some functions in the
//         abstract BoundaryVariable class

// C headers

// C++ headers
#include <cstring>    // std::memcpy
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "bvals_interfaces.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// constructor

BoundaryVariable::BoundaryVariable(MeshBlock *pmb) {
  pmy_block_ = pmb;
  pbval_ = pmb->pbval;
  pmy_mesh_ = pmb->pmy_mesh;
  bvar_index = 0;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryVariable::InitBoundaryData(BoundaryData<> &bd, BoundaryQuantity type)
//  \brief Initialize BoundaryData structure

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
//  \brief Destroy BoundaryData structure

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
//  \brief

//  Called in BoundaryVariable::SendBoundaryBuffer(), SendFluxCorrection() calls when the
//  destination neighbor block is on the same MPI rank as the sending MeshBlcok. So
//  std::memcpy() call requires pointer to "void *dst" corresponding to
//  bd_var_.recv[nb.targetid] in separate BoundaryVariable object in separate vector in
//  separate BoundaryValues

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


void BoundaryVariable::CopyFluxCorrectionBufferSameProcess(NeighborBlock& nb, int ssize) {
  // Locate target buffer
  // 1) which MeshBlock?
  MeshBlock *ptarget_block = pmy_mesh_->FindMeshBlock(nb.snb.gid);
  // 2) which element in vector of BoundaryVariable *?
  BoundaryData<> *ptarget_bdata =
      &(ptarget_block->pbval->bvars[bvar_index]->bd_var_flcor_);
  std::memcpy(ptarget_bdata->recv[nb.targetid], bd_var_flcor_.send[nb.bufid],
              ssize*sizeof(Real));
  ptarget_bdata->flag[nb.targetid] = BoundaryStatus::arrived;
  return;
}
