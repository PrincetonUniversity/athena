//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_cc.cpp
//  \brief functions that apply BCs for VERTEX_CENTERED variables

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>    // memcpy()
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../eos/eos.hpp"
#include "../../field/field.hpp"
#include "../../globals.hpp"
#include "../../hydro/hydro.hpp"
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"
#include "../../utils/buffer_utils.hpp"
#include "../bvals.hpp"
#include "bvals_vc.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// constructor

VertexCenteredBoundaryVariable::VertexCenteredBoundaryVariable(
    MeshBlock *pmb, AthenaArray<Real> *var, AthenaArray<Real> *coarse_var,
    AthenaArray<Real> *var_flux)
    : BoundaryVariable(pmb), var_vc(var), coarse_buf(coarse_var), x1flux(var_flux[X1DIR]),
      x2flux(var_flux[X2DIR]), x3flux(var_flux[X3DIR]), nl_(0), nu_(var->GetDim4() -1),
      flip_across_pole_(nullptr) {
  // VertexCenteredBoundaryVariable should only be used w/ 4D or 3D (nx4=1) AthenaArray
  // For now, assume that full span of 4th dim of input AthenaArray should be used:
  // ---> get the index limits directly from the input AthenaArray
  // <=nu_ (inclusive), <nx4 (exclusive)
  if (nu_ < 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in VertexCenteredBoundaryVariable constructor" << std::endl
        << "An 'AthenaArray<Real> *var' of nx4_ = " << var->GetDim4() << " was passed\n"
        << "Should be nx4 >= 1 (likely uninitialized)." << std::endl;
    ATHENA_ERROR(msg);
  }

  InitBoundaryData(bd_var_, BoundaryQuantity::vc);
#ifdef MPI_PARALLEL
  // KGF: dead code, leaving for now:
  // cc_phys_id_ = pbval_->ReserveTagVariableIDs(1);
  vc_phys_id_ = pbval_->bvars_next_phys_id_;
#endif
  if (pmy_mesh_->multilevel) { // SMR or AMR
    // InitBoundaryData(bd_var_flcor_, BoundaryQuantity::cc_flcor);
#ifdef MPI_PARALLEL
    // cc_flx_phys_id_ = cc_phys_id_ + 1;
#endif
  }

  // node multiplicities-------------------------------------------------------
  AllocateNodeMult();
  //---------------------------------------------------------------------------

}

// destructor
VertexCenteredBoundaryVariable::~VertexCenteredBoundaryVariable() {
  DestroyBoundaryData(bd_var_);
 // if (pmy_mesh_->multilevel)
 //   DestroyBoundaryData(bd_var_flcor_);

  // node multiplicities-------------------------------------------------------
  node_mult.DeleteAthenaArray();
  //---------------------------------------------------------------------------
}

void VertexCenteredBoundaryVariable::ErrorIfPolarNotImplemented(
  const NeighborBlock& nb) {

  // BD: TODO implement polar coordinates
  if (nb.polar) {
  std::stringstream msg;
  msg << "### FATAL ERROR" << std::endl
      << "Polar coordinates not implemented for vertex-centered." << std::endl;
  ATHENA_ERROR(msg);
  }
  return;
}

void VertexCenteredBoundaryVariable::ErrorIfShearingBoxNotImplemented() {
  // BD: TODO implement shearing box
  if (SHEARING_BOX){
    std::stringstream msg;
    msg << "### FATAL ERROR" << std::endl
        << "Shearing box not implemented for vertex-centered." << std::endl;
    ATHENA_ERROR(msg);
  }
}

int VertexCenteredBoundaryVariable::ComputeVariableBufferSize(const NeighborIndexes& ni,
                                                              int cng) {
  // 'cng' to preserve function signature but is a dummy slot
  return NeighborVariableBufferSize(ni);
}

//----------------------------------------------------------------------------------------
//! \fn int VertexCenteredBoundaryVariable::LoadBoundaryBufferSameLevel(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set vertex-centered boundary buffers for sending to a block on the same level

int VertexCenteredBoundaryVariable::LoadBoundaryBufferSameLevel(Real *buf,
                                                                const NeighborBlock& nb) {
  //MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;
  int p = 0;
  AthenaArray<Real> &var = *var_vc;

  idxLoadSameLevelRanges(nb.ni, si, ei, sj, ej, sk, ek, false);
  BufferUtility::PackData(var, buf, nl_, nu_, si, ei, sj, ej, sk, ek, p);

  // if multilevel make use of pre-restricted internal data
  if (pmy_mesh_->multilevel) {

    // convert to coarse indices
    AthenaArray<Real> &coarse_var = *coarse_buf;

    idxLoadSameLevelRanges(nb.ni, si, ei, sj, ej, sk, ek, true);
    BufferUtility::PackData(coarse_var, buf, nl_, nu_,
                            si, ei, sj, ej, sk, ek, p);
  }

  return p;
}

//----------------------------------------------------------------------------------------
//! \fn int VertexCenteredBoundaryVariable::LoadBoundaryBufferToCoarser(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set vertex-centered boundary buffers for sending to a block on the coarser level

int VertexCenteredBoundaryVariable::LoadBoundaryBufferToCoarser(Real *buf,
                                                                const NeighborBlock& nb) {
  MeshBlock *pmb = pmy_block_;
  MeshRefinement *pmr = pmb->pmr;
  int si, sj, sk, ei, ej, ek;
  int p = 0;

  AthenaArray<Real> &var = *var_vc;
  AthenaArray<Real> &coarse_var = *coarse_buf;

  // vertices that are shared with adjacent MeshBlocks are to be copied to coarser level
  idxLoadToCoarserRanges(nb.ni, si, ei, sj, ej, sk, ek, false);
  pmr->RestrictVertexCenteredValues(var, coarse_var, nl_, nu_,
                                    si, ei, sj, ej, sk, ek);

  BufferUtility::PackData(coarse_var, buf, nl_, nu_, si, ei, sj, ej, sk, ek, p);

  if (pmy_mesh_->multilevel) {
    // double restrict required to populate coarse buffer of coarser level
    idxLoadToCoarserRanges(nb.ni, si, ei, sj, ej, sk, ek, true);
    pmr->RestrictTwiceToBufferVertexCenteredValues(var, buf, nl_, nu_,
                                                   si, ei, sj, ej, sk, ek, p);
  }

  return p;
}

//----------------------------------------------------------------------------------------
//! \fn int VertexCenteredBoundaryVariable::LoadBoundaryBufferToFiner(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set vertex-centered boundary buffers for sending to a block on the finer level

int VertexCenteredBoundaryVariable::LoadBoundaryBufferToFiner(Real *buf,
                                                              const NeighborBlock& nb) {
  AthenaArray<Real> &var = *var_vc;
  int si, sj, sk, ei, ej, ek;
  int p = 0;

  idxLoadToFinerRanges(nb.ni, si, ei, sj, ej, sk, ek);
  BufferUtility::PackData(var, buf, nl_, nu_, si, ei, sj, ej, sk, ek, p);

  return p;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::SetBoundarySameLevel(Real *buf,
//                                                              const NeighborBlock& nb)
//  \brief Set vertex-centered boundary received from a block on the same level

void VertexCenteredBoundaryVariable::SetBoundarySameLevel(Real *buf,
                                                          const NeighborBlock& nb) {
  //MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;
  AthenaArray<Real> &var = *var_vc;
  int p = 0;

  // BD: TODO implement
  ErrorIfPolarNotImplemented(nb);
  ErrorIfShearingBoxNotImplemented();

  idxSetSameLevelRanges(nb.ni, si, ei, sj, ej, sk, ek, 1);

  // vertex consistency--------------------------------------------------------
  if (!node_mult_assembled) {
    int c_si, c_ei, c_sj, c_ej, c_sk, c_ek;

    idxSetSameLevelRanges(nb.ni, c_si, c_ei, c_sj, c_ej, c_sk, c_ek, 3);

    for (int k=c_sk; k<=c_ek; ++k)
      for (int j=c_sj; j<=c_ej; ++j)
        for (int i=c_si; i<=c_ei; ++i)
          node_mult(0, k, j, i) += 1;
  }
  //---------------------------------------------------------------------------

  if (pmy_mesh_->multilevel) {
    // note: unpacked shared nodes additively unpacked-
    // consistency conditions will need to be applied to the coarse variable

    //MeshRefinement *pmr = pmb->pmr;
    AthenaArray<Real> &coarse_var = *coarse_buf;

    idxSetSameLevelRanges(nb.ni, si, ei, sj, ej, sk, ek, 2);

    BufferUtility::UnpackDataAdd(buf, coarse_var, nl_, nu_,
                                 si, ei, sj, ej, sk, ek, p);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::SetBoundaryFromCoarser(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set vertex-centered prolongation buffer received from a block on a coarser level

void VertexCenteredBoundaryVariable::SetBoundaryFromCoarser(Real *buf,
                                                            const NeighborBlock& nb) {
  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;

  AthenaArray<Real> &coarse_var = *coarse_buf;


  int p = 0;
  // BD: TODO implement
  ErrorIfPolarNotImplemented(nb);

  idxSetFromCoarserRanges(nb.ni, si, ei, sj, ej, sk, ek, false);

  // vertex consistency--------------------------------------------------------
  if (!node_mult_assembled) {
    int c_si, c_ei, c_sj, c_ej, c_sk, c_ek;

    idxSetFromCoarserRanges(nb.ni, c_si, c_ei, c_sj, c_ej, c_sk, c_ek, true);

    for (int k=c_sk; k<=c_ek; ++k)
      for (int j=c_sj; j<=c_ej; ++j)
        for (int i=c_si; i<=c_ei; ++i)
          node_mult(0, k, j, i) += 1;
  }
  //---------------------------------------------------------------------------

  BufferUtility::UnpackData(buf, coarse_var, nl_, nu_,
                            si, ei, sj, ej, sk, ek, p);

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::SetBoundaryFromFiner(Real *buf,
//                                                              const NeighborBlock& nb)
//  \brief Set vertex-centered boundary received from a block on a finer level
void VertexCenteredBoundaryVariable::SetBoundaryFromFiner(Real *buf,
                                                          const NeighborBlock& nb) {
  // populating from finer level; shared vertices are corrected

  //MeshBlock *pmb = pmy_block_;
  AthenaArray<Real> &var = *var_vc;
  // receive already restricted data
  int si, sj, sk, ei, ej, ek;
  int p = 0;

  // BD: TODO implement
  ErrorIfPolarNotImplemented(nb);

  idxSetFromFinerRanges(nb.ni, si, ei, sj, ej, sk, ek, 1);

  // vertex consistency--------------------------------------------------------
  if (!node_mult_assembled) {
    int c_si, c_ei, c_sj, c_ej, c_sk, c_ek;

    idxSetFromFinerRanges(nb.ni, c_si, c_ei, c_sj, c_ej, c_sk, c_ek, 3);

    for (int k=c_sk; k<=c_ek; ++k)
      for (int j=c_sj; j<=c_ej; ++j)
        for (int i=c_si; i<=c_ei; ++i)
          node_mult(0, k, j, i) += 1;
  }
  //---------------------------------------------------------------------------

  if (pmy_mesh_->multilevel) {
    AthenaArray<Real> &coarse_var = *coarse_buf;
    idxSetFromFinerRanges(nb.ni, si, ei, sj, ej, sk, ek, 2);
    BufferUtility::UnpackDataAdd(buf, coarse_var, nl_, nu_,
                                 si, ei, sj, ej, sk, ek, p);
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::RestrictNonGhost()
//  \brief populate coarser buffer with restricted data

void VertexCenteredBoundaryVariable::RestrictNonGhost() {
  MeshBlock *pmb = pmy_block_;
  MeshRefinement *pmr = pmb->pmr;
  int si, sj, sk, ei, ej, ek;

  AthenaArray<Real> &var = *var_vc;
  AthenaArray<Real> &coarse_var = *coarse_buf;

  si = pmb->civs; ei = pmb->cive;
  sj = pmb->cjvs; ej = pmb->cjve;
  sk = pmb->ckvs; ek = pmb->ckve;

  pmr->RestrictVertexCenteredValues(var, coarse_var, nl_, nu_,
                                    si, ei, sj, ej, sk, ek);

  return;
}


void VertexCenteredBoundaryVariable::SendBoundaryBuffers() {
  if (!node_mult_assembled)
    PrepareNodeMult();

  //MeshBlock *pmb = pmy_block_;

  // restrict all data (except ghosts) to coarse buffer
  if (pmy_mesh_->multilevel) {
    AthenaArray<Real> &coarse_var = *coarse_buf;
    coarse_var.ZeroClear();

    RestrictNonGhost();
  }
  BoundaryVariable::SendBoundaryBuffers();
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::SetBoundaries()
//  \brief set the vertex-centered boundary data

void VertexCenteredBoundaryVariable::SetBoundaries() {
  ZeroVertexGhosts();
  BoundaryVariable::SetBoundaries();
  FinalizeVertexConsistency();

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::ReceiveAndSetBoundariesWithWait()
//  \brief receive and set the vertex-centered boundary data for initialization

void VertexCenteredBoundaryVariable::ReceiveAndSetBoundariesWithWait() {
  ZeroVertexGhosts();

  BoundaryVariable::ReceiveAndSetBoundariesWithWait();
  FinalizeVertexConsistency();

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::PolarBoundarySingleAzimuthalBlock()
// \brief polar boundary edge-case: single MeshBlock spans the entire azimuthal (x3) range

void VertexCenteredBoundaryVariable::PolarBoundarySingleAzimuthalBlock() {
  return;
}

void VertexCenteredBoundaryVariable::SetupPersistentMPI() {
#ifdef MPI_PARALLEL
  MeshBlock* pmb = pmy_block_;
  int &mylevel = pmb->loc.level;

  int ssize, rsize;
  int tag;
  // Initialize non-polar neighbor communications to other ranks
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (nb.snb.rank != Globals::my_rank) {
      if (nb.snb.level == mylevel) { // same
        ssize = MPI_BufferSizeSameLevel(nb.ni, true);
        rsize = MPI_BufferSizeSameLevel(nb.ni, false);
      } else if (nb.snb.level < mylevel) { // coarser
        ssize = MPI_BufferSizeToCoarser(nb.ni);
        rsize = MPI_BufferSizeFromCoarser(nb.ni);
      } else { // finer
        ssize = MPI_BufferSizeToFiner(nb.ni);
        rsize = MPI_BufferSizeFromFiner(nb.ni);
      }
      // specify the offsets in the view point of the target block: flip ox? signs

      // Initialize persistent communication requests attached to specific BoundaryData
      // vertex-centered
      tag = pbval_->CreateBvalsMPITag(nb.snb.lid, nb.targetid, vc_phys_id_);
      if (bd_var_.req_send[nb.bufid] != MPI_REQUEST_NULL)
        MPI_Request_free(&bd_var_.req_send[nb.bufid]);
      MPI_Send_init(bd_var_.send[nb.bufid], ssize, MPI_ATHENA_REAL,
                    nb.snb.rank, tag, MPI_COMM_WORLD, &(bd_var_.req_send[nb.bufid]));
      tag = pbval_->CreateBvalsMPITag(pmb->lid, nb.bufid, vc_phys_id_);
      if (bd_var_.req_recv[nb.bufid] != MPI_REQUEST_NULL)
        MPI_Request_free(&bd_var_.req_recv[nb.bufid]);
      MPI_Recv_init(bd_var_.recv[nb.bufid], rsize, MPI_ATHENA_REAL,
                    nb.snb.rank, tag, MPI_COMM_WORLD, &(bd_var_.req_recv[nb.bufid]));

    }
  }
#endif
  return;
}

void VertexCenteredBoundaryVariable::StartReceiving(BoundaryCommSubset phase) {
#ifdef MPI_PARALLEL
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (nb.snb.rank != Globals::my_rank) {
      MPI_Start(&(bd_var_.req_recv[nb.bufid]));
    }
  }
#endif
  return;
}


void VertexCenteredBoundaryVariable::ClearBoundary(BoundaryCommSubset phase) {
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    bd_var_.flag[nb.bufid] = BoundaryStatus::waiting;
    bd_var_.sflag[nb.bufid] = BoundaryStatus::waiting;

#ifdef MPI_PARALLEL
    if (nb.snb.rank != Globals::my_rank) {
      // Wait for Isend
      MPI_Wait(&(bd_var_.req_send[nb.bufid]), MPI_STATUS_IGNORE);
    }
#endif
  }

  return;
}
