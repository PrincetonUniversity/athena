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
      : BoundaryVariable(pmb), var(var), mu_(var->GetDim1() - 1), ml_(0) {
  //only take 5 dimention array set up for six ray
  if (var->GetDim5() != 6) {
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
int SixRayBoundaryVariable::ComputeVariableBufferSize(const NeighborIndexes& ni, int cng) {
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
#endif //MPI_PARALLEL
  return;
}

//----------------------------------------------------------------------------------------
//! \fn int SixRayBoundaryVariable::StartReceiving(}
//! \brief call MPI_Start for six-ray boundary
void SixRayBoundaryVariable::StartReceiving(BoundaryCommSubset phase) {
  MeshBlock *pmb = pmy_block_;
#ifdef MPI_PARALLEL
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    //only face neighbors are used in six-ray
    if (nb.ni.type == NeighborConnect::face) {
      if (nb.snb.rank != Globals::my_rank) {
        MPI_Start(&(bd_var_.req_recv[nb.bufid]));
      }
    }
  }
#endif
  return;
}

//----------------------------------------------------------------------------------------
//! \fn int SixRayBoundaryVariable::ClearBoundary(}
//! \brief Set six-ray boundary to waiting
void SixRayBoundaryVariable::ClearBoundary(BoundaryCommSubset phase) {
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    //only face neighbors are used in six-ray
    if (nb.ni.type == NeighborConnect::face) {
      bd_var_.flag[nb.bufid] = BoundaryStatus::waiting;
      bd_var_.sflag[nb.bufid] = BoundaryStatus::waiting;
    }
  }
#ifdef MPI_PARALLEL
    if (nb.snb.rank != Globals::my_rank) {
      // Wait for Isend
      MPI_Wait(&(bd_var_.req_send[nb.bufid]), MPI_STATUS_IGNORE);
    }
#endif
  return;
}

//----------------------------------------------------------------------------------------
//! \fn int SixRayBoundaryVariable::LoadBoundaryBufferSameLevel(}
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
//! \fn int SixRayBoundaryVariable::SetBoundarySameLevel(}
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
    //Set the same direction
    BufferUtility::UnpackData(buf, *var, nb.fid, nb.fid, ml_, mu_,
                              si, ei, sj, ej, sk, ek, p);
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn int SixRayBoundaryVariable::LoadBoundaryBufferToCoarser(}
//! \brief
int SixRayBoundaryVariable::LoadBoundaryBufferToCoarser(Real *buf, const NeighborBlock& nb) {
  //mesh refinement not implemented yet
  return 0;
}

//----------------------------------------------------------------------------------------
//! \fn int SixRayBoundaryVariable::LoadBoundaryBufferToFiner(}
//! \brief
int SixRayBoundaryVariable::LoadBoundaryBufferToFiner(Real *buf, const NeighborBlock& nb) {
  //mesh refinement not implemented yet
  return 0;
}

//----------------------------------------------------------------------------------------
//! \fn int SixRayBoundaryVariable::SetBoundaryFromCoarser(}
//! \brief
void SixRayBoundaryVariable::SetBoundaryFromCoarser(Real *buf, const NeighborBlock& nb) {
  //mesh refinement not implemented yet
  return;
}

//----------------------------------------------------------------------------------------
//! \fn int SixRayBoundaryVariable::SetBoundaryFromFiner(}
//! \brief
void SixRayBoundaryVariable::SetBoundaryFromFiner(Real *buf, const NeighborBlock& nb) {
  //mesh refinement not implemented yet
  return;
}

//no flux correction
int SixRayBoundaryVariable::ComputeFluxCorrectionBufferSize(const NeighborIndexes& ni, int cng) {
  return 0;
}
void SixRayBoundaryVariable::SendFluxCorrection() {}
bool SixRayBoundaryVariable::ReceiveFluxCorrection() {
  return true;
}

//shearing box not implemented
void SixRayBoundaryVariable::StartReceivingShear(BoundaryCommSubset phase) {}

//polar boundary not implemented
void SixRayBoundaryVariable::PolarBoundarySingleAzimuthalBlock() {}

//physical boundary not implemented. Use user defined boundary for six-ray
void SixRayBoundaryVariable::ReflectInnerX1(Real time, Real dt,
                    int il, int jl, int ju, int kl, int ku, int ngh) {}
void SixRayBoundaryVariable::ReflectOuterX1(Real time, Real dt,
                    int iu, int jl, int ju, int kl, int ku, int ngh) {}
void SixRayBoundaryVariable::ReflectInnerX2(Real time, Real dt,
                    int il, int iu, int jl, int kl, int ku, int ngh) {}
void SixRayBoundaryVariable::ReflectOuterX2(Real time, Real dt,
                    int il, int iu, int ju, int kl, int ku, int ngh) {}
void SixRayBoundaryVariable::ReflectInnerX3(Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ngh) {}
void SixRayBoundaryVariable::ReflectOuterX3(Real time, Real dt,
                    int il, int iu, int jl, int ju, int ku, int ngh) {}
void SixRayBoundaryVariable::OutflowInnerX1(Real time, Real dt,
                    int il, int jl, int ju, int kl, int ku, int ngh) {}
void SixRayBoundaryVariable::OutflowOuterX1(Real time, Real dt,
                    int iu, int jl, int ju, int kl, int ku, int ngh) {}
void SixRayBoundaryVariable::OutflowInnerX2(Real time, Real dt,
                    int il, int iu, int jl, int kl, int ku, int ngh) {}
void SixRayBoundaryVariable::OutflowOuterX2(Real time, Real dt,
                    int il, int iu, int ju, int kl, int ku, int ngh) {}
void SixRayBoundaryVariable::OutflowInnerX3(Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ngh) {}
void SixRayBoundaryVariable::OutflowOuterX3(Real time, Real dt,
                    int il, int iu, int jl, int ju, int ku, int ngh) {}
void SixRayBoundaryVariable::PolarWedgeInnerX2(Real time, Real dt,
                       int il, int iu, int jl, int kl, int ku, int ngh) {}
void SixRayBoundaryVariable::PolarWedgeOuterX2(Real time, Real dt,
                       int il, int iu, int ju, int kl, int ku, int ngh) {}

BoundaryFace SixRayBoundaryVariable::GetOppositeBoundaryFace(const BoundaryFace direction) {
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

void SixRayBoundaryVariable::SendSixRayBoundaryBuffers(BoundaryFace direction) {
  MeshBlock *pmb = pmy_block_;
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    //only for face neigbor in the specified direction
    if (nb.ni.type == NeighborConnect::face && nb.fid == direction) {
      if (bd_var_.sflag[nb.bufid] == BoundaryStatus::completed) continue;
      int ssize;
      ssize = LoadBoundaryBufferSameLevel(bd_var_.send[nb.bufid], nb);
      if (nb.snb.rank == Globals::my_rank) {  // on the same process
        CopyVariableBufferSameProcess(nb, ssize);
      }
#ifdef MPI_PARALLEL
      else { // MPI
        MPI_Start(&(bd_var_.req_send[nb.bufid]));
      }
#endif
      bd_var_.sflag[nb.bufid] = BoundaryStatus::completed;
      //Only done for the first match, and should be the only match for uniform mesh
      //AMR needed to be add later
      return;
    }
  }
  std::stringstream msg;
  msg << "### FATAL ERROR in SixRayBoundaryVariable::SendSixRayBoundaryBuffers()"
    << std::endl
    << "BoundaryFace " << direction  << "is undefined." << std::endl;
  ATHENA_ERROR(msg);
  return;
}

bool SixRayBoundaryVariable::ReceiveSixRayBoundaryBuffers(BoundaryFace direction) {
  bool bflag = true;

  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    //only for face neigbor in the specified direction
    if (nb.ni.type == NeighborConnect::face && nb.fid == direction) {
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
      //Only done for the first match, and should be the only match for uniform mesh
      //AMR needed to be add later
      return bflag;
    }
  }
  std::stringstream msg;
  msg << "### FATAL ERROR in SixRayBoundaryVariable::SendSixRayBoundaryBuffers()"
    << std::endl
    << "BoundaryFace " << direction  << "is undefined." << std::endl;
  ATHENA_ERROR(msg);
}
