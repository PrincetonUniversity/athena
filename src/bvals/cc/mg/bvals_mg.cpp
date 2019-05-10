//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_mg.cpp
//  \brief

// C headers

// C++ headers
#include <cmath>
#include <cstdlib>
#include <cstring>    // memcpy()
#include <iomanip>
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../coordinates/coordinates.hpp"
#include "../../../globals.hpp"
#include "../../../mesh/mesh.hpp"
#include "../../../multigrid/multigrid.hpp"
#include "../../../parameter_input.hpp"
#include "../../../utils/buffer_utils.hpp"
#include "bvals_mg.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

class Multigrid;
class MultigridDriver;

//----------------------------------------------------------------------------------------
//! \fn MGBoundaryValues::MGBoundaryValues(Multigrid *pmg, BoundaryFlag *input_bcs)
//  \brief Constructor of the MGBoundaryValues class

MGBoundaryValues::MGBoundaryValues(Multigrid *pmg, BoundaryFlag *input_bcs)
    : BoundaryBase(pmg->pmy_driver_->pmy_mesh_, pmg->loc_, pmg->size_, input_bcs),
      pmy_mg_(pmg) {
#ifdef MPI_PARALLEL
  mgcomm_ = pmg->pmy_driver_->MPI_COMM_MULTIGRID;
  // currently assuming that Multigrid gravity is the only "int physid" associated with
  // the owning MultigridDriver object:
  mg_grav_phys_id_ = pmg->pmy_driver_->mg_phys_id_;
#endif
  if (pmy_mg_->pmy_block_ == nullptr) {
    for (int i=0; i<6; i++)
      MGBoundaryFunction_[i] = pmg->pmy_driver_->MGBoundaryFunction_[i];
  } else {
    for (int i=0; i<6; i++) {
      if (block_bcs[i] == BoundaryFlag::periodic || block_bcs[i] == BoundaryFlag::block)
        MGBoundaryFunction_[i] = nullptr;
      else
        MGBoundaryFunction_[i] = pmg->pmy_driver_->MGBoundaryFunction_[i];
    }
    InitBoundaryData(bd_mggrav_, BoundaryQuantity::mggrav);
  }
}


//----------------------------------------------------------------------------------------
//! \fn MGBoundaryValues::~MGBoundaryValues()
//  \brief Destructor of the MGBoundaryValues class

MGBoundaryValues::~MGBoundaryValues() {
  if (pmy_mg_->pmy_block_ != nullptr)
    DestroyBoundaryData(bd_mggrav_);
}

//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::InitBoundaryData(BoundaryData<> &bd, BoundaryQuantity type)
//  \brief Initialize BoundaryData<> structure

void MGBoundaryValues::InitBoundaryData(BoundaryData<> &bd, BoundaryQuantity type) {
  int size = 0;
  bd.nbmax = maxneighbor_;

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

    // Allocate buffers
    // calculate the buffer size
    switch (type) {
      case BoundaryQuantity::mggrav: {
        int ngh = pmy_mg_->ngh_;
        if (pmy_mesh_->multilevel) { // with refinement - NGHOST = 1
          int nc = block_size_.nx1;
          if (BoundaryValues::ni[n].type == NeighborConnect::face)
            size = SQR(nc)*ngh;
          else if (BoundaryValues::ni[n].type == NeighborConnect::edge)
            size = nc*ngh*ngh + (nc*ngh*ngh)/2;
          else if (BoundaryValues::ni[n].type == NeighborConnect::corner)
            size = ngh*ngh*ngh*2;
        } else { // uniform - NGHOST=1
          size = ((BoundaryValues::ni[n].ox1 == 0) ? block_size_.nx1 : ngh)
                *((BoundaryValues::ni[n].ox2 == 0) ? block_size_.nx2 : ngh)
                *((BoundaryValues::ni[n].ox3 == 0) ? block_size_.nx3 : ngh);
        }
      }
        break;
      default: {
        std::stringstream msg;
        msg << "### FATAL ERROR in InitBoundaryData" << std::endl
            << "Invalid boundary type is specified." << std::endl;
        ATHENA_ERROR(msg);
      }
        break;
    }
    bd.send[n] = new Real[size];
    bd.recv[n] = new Real[size];
  }
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::DestroyBoundaryData(BoundaryData<> &bd)
//  \brief Destroy BoundaryData<> structure

void MGBoundaryValues::DestroyBoundaryData(BoundaryData<> &bd) {
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
//! \fn void MGBoundaryValues::ApplyPhysicalBoundaries()
//  \brief Apply physical boundary conditions to the current Multigrid data

void MGBoundaryValues::ApplyPhysicalBoundaries() {
  AthenaArray<Real> &dst = pmy_mg_->GetCurrentData();
  int ll = pmy_mg_->nlevel_ - 1 - pmy_mg_->current_level_;
  int ngh = pmy_mg_->ngh_, nvar = pmy_mg_->nvar_;
  int ncx = block_size_.nx1 >> ll, ncy = block_size_.nx2 >> ll,
      ncz = block_size_.nx3 >> ll;
  int is = ngh, ie = ncx + ngh - 1;
  int js = ngh, je = ncy + ngh - 1;
  int ks = ngh, ke = ncz + ngh - 1;
  // cppcheck-suppress knownConditionTrueFalse
  int bis = is - ngh;
  int bie = ie + ngh;
  int bjs = js, bje = je;
  int bks = ks, bke = ke;
  Real dx = pmy_mg_->rdx_*static_cast<Real>(1<<ll);
  Real dy = pmy_mg_->rdy_*static_cast<Real>(1<<ll);
  Real dz = pmy_mg_->rdz_*static_cast<Real>(1<<ll);
  Real x0 = block_size_.x1min - (static_cast<Real>(ngh) + 0.5)*dx;
  Real y0 = block_size_.x2min - (static_cast<Real>(ngh) + 0.5)*dy;
  Real z0 = block_size_.x3min - (static_cast<Real>(ngh) + 0.5)*dz;
  Real time = pmy_mesh_->time;
  if (MGBoundaryFunction_[BoundaryFace::inner_x2] == nullptr) bjs = js - ngh;
  if (MGBoundaryFunction_[BoundaryFace::outer_x2] == nullptr) bje = je + ngh;
  if (MGBoundaryFunction_[BoundaryFace::inner_x3] == nullptr) bks = ks - ngh;
  if (MGBoundaryFunction_[BoundaryFace::outer_x3] == nullptr) bke = ke + ngh;

  // Apply boundary function on inner-x1
  if (MGBoundaryFunction_[BoundaryFace::inner_x1] != nullptr)
    MGBoundaryFunction_[BoundaryFace::inner_x1](
        dst, time, nvar, is, ie, bjs, bje, bks, bke, ngh, x0, y0, z0, dx, dy, dz);
  // Apply boundary function on outer-x1
  if (MGBoundaryFunction_[BoundaryFace::outer_x1] != nullptr)
    MGBoundaryFunction_[BoundaryFace::outer_x1](
        dst, time, nvar, is, ie, bjs, bje, bks, bke, ngh, x0, y0, z0, dx, dy, dz);

  // Apply boundary function on inner-x2
  if (MGBoundaryFunction_[BoundaryFace::inner_x2] != nullptr)
    MGBoundaryFunction_[BoundaryFace::inner_x2](
        dst, time, nvar, bis, bie, js, je, bks, bke, ngh, x0, y0, z0, dx, dy, dz);
  // Apply boundary function on outer-x2
  if (MGBoundaryFunction_[BoundaryFace::outer_x2] != nullptr)
    MGBoundaryFunction_[BoundaryFace::outer_x2](
        dst, time, nvar, bis, bie, js, je, bks, bke, ngh, x0, y0, z0, dx, dy, dz);

  bjs = js - ngh, bje = je + ngh;
  // Apply boundary function on inner-x3
  if (MGBoundaryFunction_[BoundaryFace::inner_x3] != nullptr)
    MGBoundaryFunction_[BoundaryFace::inner_x3](
        dst, time, nvar, bis, bie, bjs, bje, ks, ke, ngh, x0, y0, z0, dx, dy, dz);
  // Apply boundary function on outer-x3
  if (MGBoundaryFunction_[BoundaryFace::outer_x3] != nullptr)
    MGBoundaryFunction_[BoundaryFace::outer_x3](
        dst, time, nvar, bis, bie, bjs, bje, ks, ke, ngh, x0, y0, z0, dx, dy, dz);

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::StartReceivingMultigrid(int nc, BoundaryQuantity type)
//  \brief initiate MPI_Irecv for multigrid

void MGBoundaryValues::StartReceivingMultigrid(int nc, BoundaryQuantity type) {
  bool faceonly = false;
#ifdef MPI_PARALLEL
  int mylevel=loc.level;
  int nvar, ngh;
  int phys; // used only in MPI calculations
  BoundaryData<> *pbd;
#endif

  if (type == BoundaryQuantity::mggrav || type == BoundaryQuantity::mggrav_f) {
#ifdef MPI_PARALLEL
    pbd = &bd_mggrav_;
    nvar = 1, ngh = 1;
    phys = mg_grav_phys_id_;
#endif
  }
  if (type == BoundaryQuantity::mggrav_f)
    faceonly = true;
  for (int n=0; n<nneighbor; n++) {
    NeighborBlock& nb = neighbor[n];
    if (faceonly && nb.ni.type > NeighborConnect::face) break;
#ifdef MPI_PARALLEL
    if (nb.snb.rank!=Globals::my_rank) {
      int size = 0;
      if (pmy_mesh_->multilevel) { // with refinement - NGHOST = 1
        if (nb.snb.level == mylevel) { // same
          if (nb.ni.type == NeighborConnect::face) size = SQR(nc);
          else if (nb.ni.type == NeighborConnect::edge) size = nc + nc/2;
          else if (nb.ni.type == NeighborConnect::corner) size = 2;
        } else if (nb.snb.level < mylevel) { // coarser
          if (nb.ni.type == NeighborConnect::face) size = SQR(nc/2 + 1);
          else if (nb.ni.type == NeighborConnect::edge) size = SQR(nc/2 + 1);
          else if (nb.ni.type == NeighborConnect::corner) size = 1;
        } else { // finer
          if (nb.ni.type == NeighborConnect::face) size = SQR(nc);
          else if (nb.ni.type == NeighborConnect::edge) size = nc/2;
          else if (nb.ni.type == NeighborConnect::corner) size = 1;
        }
      } else { // no SMR/AMR
        if (nb.ni.type == NeighborConnect::face) size = nc*nc*ngh;
        else if (nb.ni.type == NeighborConnect::edge) size = nc*ngh*ngh;
        else if (nb.ni.type == NeighborConnect::corner) size = ngh*ngh*ngh;
      }
      size *= nvar;
      int tag = CreateBvalsMPITag(pmy_mg_->pmy_block_->lid, nb.bufid, phys);
      MPI_Irecv(pbd->recv[nb.bufid], size, MPI_ATHENA_REAL, nb.snb.rank, tag,
                mgcomm_, &(pbd->req_recv[nb.bufid]));
    }
#endif
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::ClearBoundaryMultigrid(BoundaryQuantity type)
//  \brief clean up the boundary flags after each loop for multigrid

void MGBoundaryValues::ClearBoundaryMultigrid(BoundaryQuantity type) {
  bool faceonly = false;
  BoundaryData<> *pbd{};

  if (type == BoundaryQuantity::mggrav || type == BoundaryQuantity::mggrav_f)
    pbd = &bd_mggrav_;
  if (type == BoundaryQuantity::mggrav_f)
    faceonly = true;

  for (int n=0; n<nneighbor; n++) {
    NeighborBlock& nb = neighbor[n];
    if (faceonly && nb.ni.type > NeighborConnect::face) break;
    pbd->flag[nb.bufid] = BoundaryStatus::waiting;
    pbd->sflag[nb.bufid] = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
    if (nb.snb.rank != Globals::my_rank)
      MPI_Wait(&(pbd->req_send[nb.bufid]),MPI_STATUS_IGNORE); // Wait for Isend
#endif
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn int MGBoundaryValues::LoadMultigridBoundaryBufferSameLevel(AthenaArray<Real> &src,
//                            int nvar, int nc, int ngh, Real *buf,
//                            const NeighborBlock& nb)
//  \brief Set hydro boundary buffers for sending to a block on the same level

int MGBoundaryValues::LoadMultigridBoundaryBufferSameLevel(
    AthenaArray<Real> &src, int nvar, int nc, int ngh, Real *buf,
    const NeighborBlock& nb) {
  int si, sj, sk, ei, ej, ek;

  si = (nb.ni.ox1 > 0) ? nc : ngh;
  ei = (nb.ni.ox1 < 0) ? (2*ngh - 1) : (ngh + nc - 1);
  sj = (nb.ni.ox2 > 0) ? nc : ngh;
  ej = (nb.ni.ox2 < 0) ? (2*ngh - 1) : (ngh + nc - 1);
  sk = (nb.ni.ox3 > 0) ? nc : ngh;
  ek = (nb.ni.ox3 < 0) ? (2*ngh - 1) : (ngh + nc - 1);
  int p = 0;
  BufferUtility::PackData(src, buf, 0, nvar-1, si, ei, sj, ej, sk, ek, p);
  return p;
}


//----------------------------------------------------------------------------------------
//! \fn bool MGBoundaryValues::SendMultigridBoundaryBuffers(AthenaArray<Real> &src,
//                                                  int nc, BoundaryQuantity type)
//  \brief Send boundary buffers

bool MGBoundaryValues::SendMultigridBoundaryBuffers(AthenaArray<Real> &src,
                                                    int nc, BoundaryQuantity type) {
  int mylevel = loc.level;
  int nvar, ngh;
#ifdef MPI_PARALLEL
  int phys;  // used only in MPI calculations
#endif
  bool faceonly = false;
  bool bflag = true;
  BoundaryData<> *pbd{}, *ptarget{};

  if (type == BoundaryQuantity::mggrav || type == BoundaryQuantity::mggrav_f) {
    pbd = &bd_mggrav_;
    nvar = 1, ngh = 1;
#ifdef MPI_PARALLEL
    phys = mg_grav_phys_id_;
#endif
  }
  if (type == BoundaryQuantity::mggrav_f)
    faceonly = true;
  for (int n=0; n<nneighbor; n++) {
    NeighborBlock& nb = neighbor[n];
    if (faceonly && nb.ni.type > NeighborConnect::face) break;
    if (pbd->sflag[nb.bufid] == BoundaryStatus::completed) continue;
    int ssize = 0;
    if (nb.snb.rank == Globals::my_rank) {
      Multigrid *pmg = pmy_mg_->pmy_driver_->FindMultigrid(nb.snb.gid);
      if (type == BoundaryQuantity::mggrav || type == BoundaryQuantity::mggrav_f)
        ptarget = &(pmg->pmgbval->bd_mggrav_);
      if (ptarget->flag[nb.targetid] != BoundaryStatus::waiting) {
        bflag = false;
        continue;
      }
    }
    if (nb.snb.level == mylevel)
      ssize = LoadMultigridBoundaryBufferSameLevel(src, nvar, nc, ngh,
                                                   pbd->send[nb.bufid], nb);
    //    else if (nb.snb.level < mylevel)
    //      ssize = LoadMultigridBoundaryBufferToCoarser(src, nvar,
    //                                                 pbd->send[nb.bufid], cbuf, nb);
    //    else
    //     ssize = LoadMultigridBoundaryBufferToFiner(src, nvar, pbd->send[nb.bufid], nb);
    if (nb.snb.rank == Globals::my_rank) {
      std::memcpy(ptarget->recv[nb.targetid], pbd->send[nb.bufid], ssize*sizeof(Real));
      ptarget->flag[nb.targetid] = BoundaryStatus::arrived;
    }
#ifdef MPI_PARALLEL
    else { // NOLINT
      int tag = CreateBvalsMPITag(nb.snb.lid, nb.targetid, phys);
      MPI_Isend(pbd->send[nb.bufid], ssize, MPI_ATHENA_REAL, nb.snb.rank, tag,
                mgcomm_, &(pbd->req_send[nb.bufid]));
    }
#endif
    pbd->sflag[nb.bufid] = BoundaryStatus::completed;
  }
  return bflag;
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::SetMultigridBoundarySameLevel(AthenaArray<Real> &dst,
//                     int nvar, int nc, int ngh, Real *buf, const NeighborBlock& nb)
//  \brief Set hydro boundary received from a block on the same level

void MGBoundaryValues::SetMultigridBoundarySameLevel(
    AthenaArray<Real> &dst, int nvar, int nc, int ngh, Real *buf,
    const NeighborBlock& nb) {
  int si, sj, sk, ei, ej, ek;

  if (nb.ni.ox1 == 0)     si = ngh,    ei = nc + ngh - 1;
  else if (nb.ni.ox1 > 0) si = nc + ngh, ei = nc + 2*ngh - 1;
  else              si = 0,      ei = ngh - 1;
  if (nb.ni.ox2 == 0)     sj = ngh,    ej = nc + ngh - 1;
  else if (nb.ni.ox2 > 0) sj = nc + ngh, ej = nc + 2*ngh - 1;
  else              sj = 0,      ej = ngh - 1;
  if (nb.ni.ox3 == 0)     sk = ngh,    ek = nc + ngh - 1;
  else if (nb.ni.ox3 > 0) sk = nc + ngh, ek = nc + 2*ngh - 1;
  else              sk = 0,      ek = ngh - 1;

  int p = 0;
  BufferUtility::UnpackData(buf, dst, 0, nvar-1, si, ei, sj, ej, sk, ek, p);
  return;
}


//----------------------------------------------------------------------------------------
//! \fn bool MGBoundaryValues::ReceiveMultigridBoundaryBuffers(AthenaArray<Real> &dst,
//                                                     int nc, BoundaryQuantity type)
//  \brief receive the boundary data

bool MGBoundaryValues::ReceiveMultigridBoundaryBuffers(AthenaArray<Real> &dst,
                                                       int nc,
                                                       BoundaryQuantity type) {
  bool bflag = true, faceonly = false;
  int nvar, ngh;
  BoundaryData<> *pbd{};

  if (type == BoundaryQuantity::mggrav || type == BoundaryQuantity::mggrav_f) {
    pbd = &bd_mggrav_;
    nvar = 1, ngh = 1;
  }
  if (type == BoundaryQuantity::mggrav_f)
    faceonly = true;

  for (int n=0; n<nneighbor; n++) {
    NeighborBlock& nb = neighbor[n];
    if (faceonly && nb.ni.type>NeighborConnect::face) break;
    if (pbd->flag[nb.bufid] == BoundaryStatus::completed) continue;
    if (pbd->flag[nb.bufid] == BoundaryStatus::waiting) {
      if (nb.snb.rank == Globals::my_rank) {// on the same process
        bflag = false;
        continue;
      }
#ifdef MPI_PARALLEL
      else { // NOLINT
        int test;
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, mgcomm_, &test, MPI_STATUS_IGNORE);
        MPI_Test(&(pbd->req_recv[nb.bufid]), &test, MPI_STATUS_IGNORE);
        if (!static_cast<bool>(test)) {
          bflag = false;
          continue;
        }
        pbd->flag[nb.bufid] = BoundaryStatus::arrived;
      }
#endif
    }
    if (nb.snb.level == loc.level)
      SetMultigridBoundarySameLevel(dst, nvar, nc, ngh, pbd->recv[nb.bufid], nb);
    //    else if (nb.snb.level<loc.level) // this set only the prolongation buffer
    //      SetMultigridBoundaryFromCoarser(nvar, nc, ngh, pbd->recv[nb.bufid], cbuf, nb);
    //    else
    //      SetMultigridBoundaryFromFiner(dst, nvar, nc, ngh, pbd->recv[nb.bufid], nb);
    pbd->flag[nb.bufid] = BoundaryStatus::completed; // completed
  }
  return bflag;
}

//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::CopyNeighborInfoFromMeshBlock()
//  \brief copy the neighbor information from the MeshBlock BoundaryValues class
void MGBoundaryValues::CopyNeighborInfoFromMeshBlock() {
  BoundaryValues *pbv = pmy_mg_->pmy_block_->pbval;
  nneighbor = pbv->nneighbor;
  for (int k=0; k<=2; k++) {
    for (int j=0; j<=2; j++) {
      for (int i=0; i<=2; i++)
        nblevel[k][j][i] = pbv->nblevel[k][j][i];
    }
  }
  for (int n=0; n<nneighbor; n++)
    neighbor[n] = pbv->neighbor[n];
}

