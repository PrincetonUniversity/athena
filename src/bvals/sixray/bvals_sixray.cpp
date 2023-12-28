//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_sixray.cpp
//! \brief functions that apply BCs for six-ray column density variables

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
#include "bvals_sixray.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

//! constructor
//
SixRayBoundaryVariable::SixRayBoundaryVariable(MeshBlock *pmb, AthenaArray<Real> *var)
      : BoundaryVariable(pmb, false), var(var), mu_(var->GetDim1() - 1), ml_(0) {
  //only take 0 size array (for const radiation)
  //or 5 dimention array set up for six ray
  if (var->GetSize() != 0 && var->GetDim5() != 6) {
    std::stringstream msg;
    msg << "### FATAL ERROR in SixRayBoundaryVariable constructor" << std::endl
        << "An 'AthenaArray<Real> *var' of nx5_ = " << var->GetDim5() << " was passed\n"
        << "Should be nx5 = 6 (likely uninitialized)." << std::endl;
    ATHENA_ERROR(msg);
  }

  //This will initialize the maximum neighbours similar to the cell-centered variable.
  //Leaving it for now for potential extension to mesh refinement in the future.
  InitBoundaryData(bd_var_, BoundaryQuantity::cc);
#ifdef MPI_PARALLEL
  sixray_phys_id_ = pbval_->bvars_next_phys_id_;
#endif

  if ((pmy_mesh_->multilevel)
      || (pbval_->shearing_box != 0)) { // SMR or AMR or SHEARING_BOX
    std::stringstream msg;
    msg << "### FATAL ERROR in SixRayBoundaryVariable constructor" << std::endl
        << "Not yet compatible with mesh refinement or shearingbox" << std::endl;
    ATHENA_ERROR(msg);
  }
}

SixRayBoundaryVariable::~SixRayBoundaryVariable() {
  DestroyBoundaryData(bd_var_);
}

//----------------------------------------------------------------------------------------
//! \fn int SixRayBoundaryVariable::ComputeVariableBufferSize(
//!     const NeighborIndexes& ni, int cng)
//! \brief calculate buffer size for six ray boundary
int SixRayBoundaryVariable::ComputeVariableBufferSize(const NeighborIndexes& ni,
                                                      int cng) {
  int size;
  MeshBlock *pmb = pmy_block_;

  if (ni.type == NeighborConnect::face) {
    size = ((ni.ox1 == 0) ? pmb->block_size.nx1 : 1)
           *((ni.ox2 == 0) ? pmb->block_size.nx2 : 1)
           *((ni.ox3 == 0) ? pmb->block_size.nx3 : 1);
    size *= mu_ + 1;
  } else {
    size = 0;
  }
  return size;
}


//----------------------------------------------------------------------------------------
//! \fn int SixRayBoundaryVariable::SetupPersistentMPI()
//! \brief Initialize MPI send and receive for six ray
void SixRayBoundaryVariable::SetupPersistentMPI() {
#ifdef MPI_PARALLEL
  MeshBlock* pmb = pmy_block_;
  int ssize, rsize;
  int tag;
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    //only face neighbors are used in six-ray
    if (nb.ni.type == NeighborConnect::face) {
      ssize = rsize = ((nb.ni.ox1 == 0) ? pmb->block_size.nx1 : 1)
            *((nb.ni.ox2 == 0) ? pmb->block_size.nx2 : 1)
            *((nb.ni.ox3 == 0) ? pmb->block_size.nx3 : 1);
      ssize *= (mu_ + 1); rsize *= (mu_ + 1);
      tag = pbval_->CreateBvalsMPITag(nb.snb.lid, nb.targetid, sixray_phys_id_);
      if (bd_var_.req_send[nb.bufid] != MPI_REQUEST_NULL) {
        MPI_Request_free(&bd_var_.req_send[nb.bufid]);
      }
      MPI_Send_init(bd_var_.send[nb.bufid], ssize, MPI_ATHENA_REAL,
                    nb.snb.rank, tag, MPI_COMM_WORLD, &(bd_var_.req_send[nb.bufid]));
      tag = pbval_->CreateBvalsMPITag(pmb->lid, nb.bufid, sixray_phys_id_);
      if (bd_var_.req_recv[nb.bufid] != MPI_REQUEST_NULL) {
        MPI_Request_free(&bd_var_.req_recv[nb.bufid]);
      }
      MPI_Recv_init(bd_var_.recv[nb.bufid], rsize, MPI_ATHENA_REAL,
                    nb.snb.rank, tag, MPI_COMM_WORLD, &(bd_var_.req_recv[nb.bufid]));
    }
  }
#endif
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void SixRayBoundaryVariable::StartReceiving(BoundaryCommSubset phase)
//! \brief call MPI_Start for six-ray boundary
void SixRayBoundaryVariable::StartReceiving(BoundaryCommSubset phase) {
#ifdef MPI_PARALLEL
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    // only face neighbors are used in six-ray
    if (nb.ni.type == NeighborConnect::face && nb.snb.rank != Globals::my_rank) {
      MPI_Start(&(bd_var_.req_recv[nb.bufid]));
    }
  }
#endif
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void SixRayBoundaryVariable::ClearBoundary(BoundaryCommSubset phase)
//! \brief Set six-ray boundary to waiting
void SixRayBoundaryVariable::ClearBoundary(BoundaryCommSubset phase) {
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    // only face neighbors are used in six-ray
    if (nb.ni.type == NeighborConnect::face) {
      bd_var_.flag[nb.bufid] = BoundaryStatus::waiting;
      bd_var_.sflag[nb.bufid] = BoundaryStatus::waiting;
    }
#ifdef MPI_PARALLEL
    if (nb.snb.rank != Globals::my_rank
        && nb.ni.type == NeighborConnect::face) {
      // Wait for Isend
      MPI_Wait(&(bd_var_.req_send[nb.bufid]), MPI_STATUS_IGNORE);
    }
#endif
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn int SixRayBoundaryVariable::LoadBoundaryBufferSameLevel(Real *buf,
//!                                                             const NeighborBlock& nb)
//! \brief Set six-ray boundary buffers for sending to a block on the same level
int SixRayBoundaryVariable::LoadBoundaryBufferSameLevel(Real *buf,
                                                        const NeighborBlock& nb) {
  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;
  int p = 0;
  BoundaryFace fid_opp = GetOppositeBoundaryFace(nb.fid);

  //only face neighbors are used in six ray
  if (nb.ni.type == NeighborConnect::face) {
    si = (nb.ni.ox1 > 0) ? (pmb->ie) : pmb->is;
    ei = (nb.ni.ox1 < 0) ? (pmb->is) : pmb->ie;
    sj = (nb.ni.ox2 > 0) ? (pmb->je) : pmb->js;
    ej = (nb.ni.ox2 < 0) ? (pmb->js) : pmb->je;
    sk = (nb.ni.ox3 > 0) ? (pmb->ke) : pmb->ks;
    ek = (nb.ni.ox3 < 0) ? (pmb->ks) : pmb->ke;
    //Load the opposite direction
    BufferUtility::PackData(*var, buf, fid_opp, fid_opp, ml_, mu_,
                            si, ei, sj, ej, sk, ek, p);
  }
  return p;
}


//----------------------------------------------------------------------------------------
//! \fn void SixRayBoundaryVariable::SetBoundarySameLevel(Real *buf,
//!                                                       const NeighborBlock& nb)
//! \brief
void SixRayBoundaryVariable::SetBoundarySameLevel(Real *buf, const NeighborBlock& nb) {
  std::stringstream msg;
  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;

  if (nb.ni.ox1 == 0)     si = pmb->is,        ei = pmb->ie;
  else if (nb.ni.ox1 > 0) si = pmb->ie + 1,    ei = pmb->ie + 1;
  else                    si = pmb->is - 1,    ei = pmb->is - 1;
  if (nb.ni.ox2 == 0)     sj = pmb->js,        ej = pmb->je;
  else if (nb.ni.ox2 > 0) sj = pmb->je + 1,    ej = pmb->je + 1;
  else                    sj = pmb->js - 1,    ej = pmb->js - 1;
  if (nb.ni.ox3 == 0)     sk = pmb->ks,        ek = pmb->ke;
  else if (nb.ni.ox3 > 0) sk = pmb->ke + 1,    ek = pmb->ke + 1;
  else                    sk = pmb->ks - 1,    ek = pmb->ks - 1;

  int p = 0;

  if (nb.polar) {
    msg << "### FATAL ERROR in SixRayBoundaryVariable::SetBoundarySameLevel()"
      << std::endl
      << "Chemistry BC not yet working with polar coordinates yet." << std::endl;
    ATHENA_ERROR(msg);
  } else {
    // Set the same direction
    BufferUtility::UnpackData(buf, *var, nb.fid, nb.fid, ml_, mu_,
                              si, ei, sj, ej, sk, ek, p);
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn int SixRayBoundaryVariable::LoadBoundaryBufferToCoarser(Real *buf,
//!                                                             const NeighborBlock& nb)
//! \brief
int SixRayBoundaryVariable::LoadBoundaryBufferToCoarser(Real *buf,
                                                        const NeighborBlock& nb) {
  // mesh refinement not implemented yet
  return 0;
}


//----------------------------------------------------------------------------------------
//! \fn int SixRayBoundaryVariable::LoadBoundaryBufferToFiner(Real *buf,
//!                                                           const NeighborBlock& nb)
//! \brief
int SixRayBoundaryVariable::LoadBoundaryBufferToFiner(Real *buf,
                                                      const NeighborBlock& nb) {
  // mesh refinement not implemented yet
  return 0;
}


//----------------------------------------------------------------------------------------
//! \fn void SixRayBoundaryVariable::SetBoundaryFromCoarser(Real *buf,
//!                                                         const NeighborBlock& nb)
//! \brief
void SixRayBoundaryVariable::SetBoundaryFromCoarser(Real *buf, const NeighborBlock& nb) {
  // mesh refinement not implemented yet
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void SixRayBoundaryVariable::SetBoundaryFromFiner(Real *buf,
//!                                                       const NeighborBlock& nb)
//! \brief
void SixRayBoundaryVariable::SetBoundaryFromFiner(Real *buf, const NeighborBlock& nb) {
  // mesh refinement not implemented yet
  return;
}


BoundaryFace SixRayBoundaryVariable::GetOppositeBoundaryFace(
                                       const BoundaryFace direction) {
  BoundaryFace opp_direction;

    if (direction == BoundaryFace::inner_x1) {
    opp_direction = BoundaryFace::outer_x1;
  } else if (direction == BoundaryFace::outer_x1) {
    opp_direction = BoundaryFace::inner_x1;
  } else if (direction == BoundaryFace::inner_x2) {
    opp_direction = BoundaryFace::outer_x2;
  } else if (direction == BoundaryFace::outer_x2) {
    opp_direction = BoundaryFace::inner_x2;
  } else if (direction == BoundaryFace::inner_x3) {
    opp_direction = BoundaryFace::outer_x3;
  } else if (direction == BoundaryFace::outer_x3) {
    opp_direction = BoundaryFace::inner_x3;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in BoundaryValues::GetOppositeBoundaryFace()" << std::endl
      << "BoundaryFace " << direction  << "is undefined." << std::endl;
    ATHENA_ERROR(msg);
  }

  return opp_direction;
}
NeighborBlock *SixRayBoundaryVariable::GetFaceNeighbor(const BoundaryFace direction) {
  NeighborBlock *pnb = nullptr;
  for (int n=0; n<pbval_->nneighbor; n++) {
    pnb = &pbval_->neighbor[n];
    // Only done for the first match, and should be the only match for uniform mesh
    // AMR needed to be added at a later date
    if (pnb->ni.type == NeighborConnect::face && pnb->fid == direction) {
      return pnb;
    }
  }
  return nullptr;
}


void SixRayBoundaryVariable::SendSixRayBoundaryBuffers(const BoundaryFace direction) {
  NeighborBlock *pnb = GetFaceNeighbor(direction);
  //only for face neigbor in the specified direction
  if (pnb != nullptr) {
    if (bd_var_.sflag[pnb->bufid] == BoundaryStatus::completed) {
      return;
    }
    int ssize;
    ssize = LoadBoundaryBufferSameLevel(bd_var_.send[pnb->bufid], *pnb);
    if (pnb->snb.rank == Globals::my_rank) {  // on the same process
      CopyVariableBufferSameProcess(*pnb, ssize);
    }
#ifdef MPI_PARALLEL
    else { // NOLINT // MPI
      MPI_Start(&(bd_var_.req_send[pnb->bufid]));
    }
#endif
    bd_var_.sflag[pnb->bufid] = BoundaryStatus::completed;
    return;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in SixRayBoundaryVariable::SendSixRayBoundaryBuffers()"
      << std::endl
      << "BoundaryFace " << direction  << "is undefined, or has no neighbor."
      << std::endl;
    ATHENA_ERROR(msg);
  }
  return;
}


bool SixRayBoundaryVariable::ReceiveAndSetSixRayBoundaryBuffers(
                               const BoundaryFace direction) {
  bool bflag = true;
  NeighborBlock *pnb = GetFaceNeighbor(direction);
  // only for face neigbor in the specified direction
  if (pnb != nullptr) {
    if (bd_var_.flag[pnb->bufid] == BoundaryStatus::arrived
        || bd_var_.flag[pnb->bufid] == BoundaryStatus::completed) {
      bflag = true;
    } else if (bd_var_.flag[pnb->bufid] == BoundaryStatus::waiting) {
      if (pnb->snb.rank == Globals::my_rank) {  // on the same process
        bflag = false;
      }
#ifdef MPI_PARALLEL
      else { // NOLINT // MPI boundary
        int test;
        // probe MPI communications.  This is a bit of black magic that seems to promote
        // communications to top of stack and gets them to complete more quickly
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &test, MPI_STATUS_IGNORE);
        MPI_Test(&(bd_var_.req_recv[pnb->bufid]), &test, MPI_STATUS_IGNORE);
        if (!static_cast<bool>(test)) {
          bflag = false;
        }
        bd_var_.flag[pnb->bufid] = BoundaryStatus::arrived;
      }
#endif
    }
    // set boundary
    if (bd_var_.flag[pnb->bufid] == BoundaryStatus::arrived) {
#ifdef MPI_PARALLEL
      if (pnb->snb.rank != Globals::my_rank) {
        MPI_Wait(&(bd_var_.req_recv[pnb->bufid]),MPI_STATUS_IGNORE);
      }
#endif
      SetBoundarySameLevel(bd_var_.recv[pnb->bufid], *pnb);
      bd_var_.flag[pnb->bufid] = BoundaryStatus::completed;
    }
    return bflag;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in SixRayBoundaryVariable::SendSixRayBoundaryBuffers()"
      << std::endl
      << "BoundaryFace " << direction  << " has no neighbor."
      << std::endl;
    ATHENA_ERROR(msg);
  }
}
