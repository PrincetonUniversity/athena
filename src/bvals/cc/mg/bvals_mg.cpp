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
    InitBoundaryData(BoundaryQuantity::mggrav);
    int nc = block_size_.nx1 + 2*pmy_mg_->ngh_;
    cbuf_.NewAthenaArray(pmy_mg_->nvar_, nc, nc, nc);
  }
}


//----------------------------------------------------------------------------------------
//! \fn MGBoundaryValues::~MGBoundaryValues()
//  \brief Destructor of the MGBoundaryValues class

MGBoundaryValues::~MGBoundaryValues() {
  if (pmy_mg_->pmy_block_ != nullptr)
    DestroyBoundaryData();
}

//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::InitBoundaryData(BoundaryQuantity type)
//  \brief Initialize BoundaryData<> structure

void MGBoundaryValues::InitBoundaryData(BoundaryQuantity type) {
  int size = 0;
  bdata_.nbmax = maxneighbor_;

  for (int n=0; n<bdata_.nbmax; n++) {
    // Clear flags and requests
    bdata_.flag[n] = BoundaryStatus::waiting;
    bdata_.sflag[n] = BoundaryStatus::waiting;
    bdata_.send[n] = nullptr;
    bdata_.recv[n] = nullptr;
#ifdef MPI_PARALLEL
    bdata_.req_send[n] = MPI_REQUEST_NULL;
    bdata_.req_recv[n] = MPI_REQUEST_NULL;
#endif

    // Allocate buffers
    // calculate the buffer size
    switch (type) {
      case BoundaryQuantity::mggrav: {
        int ngh = pmy_mg_->ngh_;
        int nc = block_size_.nx1;
        if (BoundaryValues::ni[n].type == NeighborConnect::face)
          size = SQR(nc)*ngh;
        else if (BoundaryValues::ni[n].type == NeighborConnect::edge)
          size = nc*ngh*ngh;
        else if (BoundaryValues::ni[n].type == NeighborConnect::corner)
          size = ngh*ngh*ngh;
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
    size *= pmy_mg_->nvar_;
    bdata_.send[n] = new Real[size];
    bdata_.recv[n] = new Real[size];
  }
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::DestroyBoundaryData()
//  \brief Destroy BoundaryData<> structure

void MGBoundaryValues::DestroyBoundaryData() {
  for (int n=0; n<bdata_.nbmax; n++) {
    delete [] bdata_.send[n];
    delete [] bdata_.recv[n];
#ifdef MPI_PARALLEL
    if (bdata_.req_send[n] != MPI_REQUEST_NULL)
      MPI_Request_free(&bdata_.req_send[n]);
    if (bdata_.req_recv[n] != MPI_REQUEST_NULL)
      MPI_Request_free(&bdata_.req_recv[n]);
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
//! \fn void MGBoundaryValues::StartReceivingMultigrid(BoundaryQuantity type)
//  \brief initiate MPI_Irecv for multigrid

void MGBoundaryValues::StartReceivingMultigrid(BoundaryQuantity type) {
#ifdef MPI_PARALLEL
  bool faceonly = false;
  int mylevel=loc.level;
  int nvar = pmy_mg_->nvar_, ngh = pmy_mg_->ngh_, cngh = ngh;
  int nc = pmy_mg_->GetCurrentNumberOfCells();

  if (type == BoundaryQuantity::mggrav_f)
    faceonly = true;

  for (int n=0; n<nneighbor; n++) {
    NeighborBlock& nb = neighbor[n];
    if (faceonly && nb.ni.type > NeighborConnect::face) break;
    if (nb.snb.rank!=Globals::my_rank) {
      int size = 0;
      if (type == BoundaryQuantity::mggrav_f) { // mass conservation formula for gravity
        if (nb.snb.level == mylevel) size = nc*nc*ngh;
        else if (nb.snb.level < mylevel) size = SQR(nc/2)*ngh;
        else size = SQR(nc/2)*ngh;
      } else { // default, full boundary communication
        if (nb.snb.level == mylevel) { // same
          if (nb.ni.type == NeighborConnect::face) size = SQR(nc)*ngh;
          else if (nb.ni.type == NeighborConnect::edge) size = nc*ngh*ngh;
          else if (nb.ni.type == NeighborConnect::corner) size = ngh*ngh*ngh;
          if (pmy_mg_->pmy_driver_->pmy_mesh_->multilevel) {
            if (nb.ni.type == NeighborConnect::face) size += SQR(nc/2)*ngh;
            else if (nb.ni.type == NeighborConnect::edge) size += nc/2*ngh*ngh;
            else if (nb.ni.type == NeighborConnect::corner) size += ngh*ngh*ngh;
          }
        } else if (nb.snb.level < mylevel) { // coarser
          if (nb.ni.type == NeighborConnect::face) size = SQR(nc/2+cngh)*cngh;
          else if (nb.ni.type == NeighborConnect::edge) size = (nc/2+cngh)*cngh*cngh;
          else if (nb.ni.type == NeighborConnect::corner) size = cngh*cngh*cngh;
        } else { // finer
          if (nb.ni.type == NeighborConnect::face) size = SQR(nc/2)*ngh;
          else if (nb.ni.type == NeighborConnect::edge) size = nc/2*ngh*ngh;
          else if (nb.ni.type == NeighborConnect::corner) size = ngh*ngh*ngh;
        }
      }
      size *= nvar;
      int tag = CreateBvalsMPITag(pmy_mg_->pmy_block_->lid, nb.bufid,
                                  pmy_mg_->pmy_driver_->mg_phys_id_);
      MPI_Irecv(bdata_.recv[nb.bufid], size, MPI_ATHENA_REAL, nb.snb.rank, tag,
                mgcomm_, &(bdata_.req_recv[nb.bufid]));
    }
  }
#endif
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::ClearBoundaryMultigrid(BoundaryQuantity type)
//  \brief clean up the boundary flags after each loop for multigrid

void MGBoundaryValues::ClearBoundaryMultigrid(BoundaryQuantity type) {
  bool faceonly = false;

  if (type == BoundaryQuantity::mggrav_f)
    faceonly = true;

  for (int n=0; n<nneighbor; n++) {
    NeighborBlock& nb = neighbor[n];
    if (faceonly && nb.ni.type > NeighborConnect::face) break;
    bdata_.flag[nb.bufid] = BoundaryStatus::waiting;
    bdata_.sflag[nb.bufid] = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
    if (nb.snb.rank != Globals::my_rank)
      MPI_Wait(&(bdata_.req_send[nb.bufid]),MPI_STATUS_IGNORE); // Wait for Isend
#endif
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn int MGBoundaryValues::LoadMultigridBoundaryBufferSameLevel(Real *buf,
//                                                           const NeighborBlock& nb)
//  \brief Set Multigrid boundary buffers for sending to a block on the same level

int MGBoundaryValues::LoadMultigridBoundaryBufferSameLevel(Real *buf,
                                                           const NeighborBlock& nb) {
  const AthenaArray<Real> &src = pmy_mg_->GetCurrentData();
  int nc = pmy_mg_->GetCurrentNumberOfCells();
  int ngh = pmy_mg_->ngh_, nvar = pmy_mg_->nvar_;
  int fs = ngh, fe = fs + nc - 1;
  int si = (nb.ni.ox1 > 0) ? (fe - ngh + 1) : fs;
  int ei = (nb.ni.ox1 < 0) ? (fs + ngh - 1) : fe;
  int sj = (nb.ni.ox2 > 0) ? (fe - ngh + 1) : fs;
  int ej = (nb.ni.ox2 < 0) ? (fs + ngh - 1) : fe;
  int sk = (nb.ni.ox3 > 0) ? (fe - ngh + 1) : fs;
  int ek = (nb.ni.ox3 < 0) ? (fs + ngh - 1) : fe;

  int p = 0;
  BufferUtility::PackData(src, buf, 0, nvar-1, si, ei, sj, ej, sk, ek, p);

  if (pmy_mg_->pmy_driver_->pmy_mesh_->multilevel) {
    int cn = ngh, fs = ngh, cs = cn, fe = fs + nc - 1, ce = cs + nc/2 - 1;
    si = (nb.ni.ox1 > 0) ? (ce - cn) : cs;
    ei = (nb.ni.ox1 < 0) ? (cs + cn) : ce;
    sj = (nb.ni.ox2 > 0) ? (ce - cn) : cs;
    ej = (nb.ni.ox2 < 0) ? (cs + cn) : ce;
    sk = (nb.ni.ox3 > 0) ? (ce - cn) : cs;
    ek = (nb.ni.ox3 < 0) ? (cs + cn) : ce;
    int fsi = (si - ngh)*2 + ngh;
    int fsj = (sj - ngh)*2 + ngh;
    int fsk = (sk - ngh)*2 + ngh;
    for (int n=0; n<nvar; ++n) {
      for (int k=sk, fk=fsk; k<=ek; ++k, fk+=2) {
        for (int j=sj, fj=fsj; j<=ej; ++j, fj+=2) {
          for (int i=si, fi=fsi; i<=ei; ++i, fi+=2)
            buf[p++] = cbuf_(n, k, j, i)
                     = 0.125*(((src(n, fk,   fj,   fi)+src(n, fk,   fj,   fi+1))
                            +  (src(n, fk,   fj+1, fi)+src(n, fk,   fj+1, fi+1)))
                            + ((src(n, fk+1, fj,   fi)+src(n, fk+1, fj,   fi+1))
                            +  (src(n, fk+1, fj+1, fi)+src(n, fk+1, fj+1, fi+1))));
        }
      }
    }
  }

  return p;
}


//----------------------------------------------------------------------------------------
//! \fn int MGBoundaryValues::LoadMultigridBoundaryBufferToCoarser(Real *buf,
//                                                         const NeighborBlock& nb)
//  \brief Set Multigrid boundary buffers for sending to a block on the coarser level

int MGBoundaryValues::LoadMultigridBoundaryBufferToCoarser(Real *buf,
                                                           const NeighborBlock& nb) {
  const AthenaArray<Real> &src = pmy_mg_->GetCurrentData();
  int ngh = pmy_mg_->ngh_, nvar = pmy_mg_->nvar_;
  int nc = pmy_mg_->GetCurrentNumberOfCells();
  int cn = ngh, fs = ngh, cs = cn, fe = fs + nc - 1, ce = cs + nc/2 - 1;
  int p=0;
  int si = (nb.ni.ox1 > 0) ? (ce - cn) : cs;
  int ei = (nb.ni.ox1 < 0) ? (cs + cn) : ce;
  int sj = (nb.ni.ox2 > 0) ? (ce - cn) : cs;
  int ej = (nb.ni.ox2 < 0) ? (cs + cn) : ce;
  int sk = (nb.ni.ox3 > 0) ? (ce - cn) : cs;
  int ek = (nb.ni.ox3 < 0) ? (cs + cn) : ce;

  // restrict and store in the coarse buffer
  int fsi = (si - ngh)*2 + ngh;
  int fsj = (sj - ngh)*2 + ngh;
  int fsk = (sk - ngh)*2 + ngh;
  for (int n=0; n<nvar; ++n) {
    for (int k=sk, fk=fsk; k<=ek; ++k, fk+=2) {
      for (int j=sj, fj=fsj; j<=ej; ++j, fj+=2) {
        for (int i=si, fi=fsi; i<=ei; ++i, fi+=2)
          buf[p++] = cbuf_(n, k, j, i)
                   = 0.125*(((src(n, fk,   fj,   fi)+src(n, fk,   fj,   fi+1))
                          +  (src(n, fk,   fj+1, fi)+src(n, fk,   fj+1, fi+1)))
                          + ((src(n, fk+1, fj,   fi)+src(n, fk+1, fj,   fi+1))
                          +  (src(n, fk+1, fj+1, fi)+src(n, fk+1, fj+1, fi+1))));
      }
    }
  }

  return p;
}


//----------------------------------------------------------------------------------------
//! \fn int MGBoundaryValues::LoadMultigridBoundaryBufferToFiner(Real *buf,
//                                                               const NeighborBlock& nb)
//  \brief Set Multigrid boundary buffers for sending to a block on the finer level
int MGBoundaryValues::LoadMultigridBoundaryBufferToFiner(Real *buf,
                                                         const NeighborBlock& nb) {
  const AthenaArray<Real> &src = pmy_mg_->GetCurrentData();
  int nc = pmy_mg_->GetCurrentNumberOfCells();
  int ngh = pmy_mg_->ngh_, nvar = pmy_mg_->nvar_;
  int cn = 1, fs = ngh, fe = fs + nc - 1;
  int si = (nb.ni.ox1 > 0) ? (fe - cn) : fs;
  int ei = (nb.ni.ox1 < 0) ? (fs + cn) : fe;
  int sj = (nb.ni.ox2 > 0) ? (fe - cn) : fs;
  int ej = (nb.ni.ox2 < 0) ? (fs + cn) : fe;
  int sk = (nb.ni.ox3 > 0) ? (fe - cn) : fs;
  int ek = (nb.ni.ox3 < 0) ? (fs + cn) : fe;

  // send the data first and later prolongate on the target block
  // need to add edges for faces, add corners for edges
  if (nb.ni.ox1 == 0) {
    if (nb.ni.fi1 == 1)   si += nc/2 - ngh;
    else                  ei -= nc/2 - ngh;
  }
  if (nb.ni.ox2 == 0) {
    if (nb.ni.ox1 != 0) {
      if (nb.ni.fi1 == 1) sj += nc/2 - ngh;
      else                ej -= nc/2 - ngh;
    } else {
      if (nb.ni.fi2 == 1) sj += nc/2 - ngh;
      else                ej -= nc/2 - ngh;
    }
  }
  if (nb.ni.ox3 == 0) {
    if (nb.ni.ox1 != 0 && nb.ni.ox2 != 0) {
      if (nb.ni.fi1 == 1) sk += nc/2 - ngh;
      else                ek -= nc/2 - ngh;
    } else {
      if (nb.ni.fi2 == 1) sk += nc/2 - ngh;
      else                ek -= nc/2 - ngh;
    }
  }

  int p = 0;
  BufferUtility::PackData(src, buf, 0, nvar-1, si, ei, sj, ej, sk, ek, p);

  return p;
}


//----------------------------------------------------------------------------------------
//! \fn int MGBoundaryValues::LoadMultigridBoundaryBufferToCoarserGravityMassCons(
//                                                Real *buf, const NeighborBlock& nb)
//  \brief Set Multigrid boundary buffers for sending to a block on the coarser level
//         using the Mass Conservation formula for gravity

int MGBoundaryValues::LoadMultigridBoundaryBufferToCoarserGravityMassCons(
                                                 Real *buf, const NeighborBlock& nb) {
  const AthenaArray<Real> &src = pmy_mg_->GetCurrentData();
  int nc = pmy_mg_->GetCurrentNumberOfCells();
  int ngh = pmy_mg_->ngh_, nvar = pmy_mg_->nvar_;
  int cn = ngh, fs = ngh, cs = cn, fe = fs + nc - 1, ce = cs + nc/2 - 1;
  int p = 0;

  // take averages of four cells contacting the interface and pack
  if (nb.ni.ox1 != 0) { // x1 face
    int fi;
    if (nb.ni.ox1 < 0) fi = fs;
    else               fi = fe;
    for (int n=0; n<nvar; ++n) {
      for (int k=cs, fk=fs; k<=ce; ++k, fk+=2) {
        for (int j=cs, fj=fs; j<=ce; ++j, fj+=2)
          buf[p++] = 0.25*((src(n, k,   j,   fi)+src(n, k,   j+1, fi))
                          +(src(n, k+1, j,   fi)+src(n, k+1, j+1, fi)));
      }
    }
  } else if (nb.ni.ox2 != 0) { // x2 face
    int fj;
    if (nb.ni.ox2 < 0) fj = fs;
    else               fj = fe;
    for (int n=0; n<nvar; ++n) {
      for (int k=cs, fk=fs; k<=ce; ++k, fk+=2) {
        for (int i=cs, fi=fs; i<=ce; ++i, fi+=2)
          buf[p++] = 0.25*((src(n, k,   fj, i)+src(n, k,   fj, i+1))
                          +(src(n, k+1, fj, i)+src(n, k+1, fj, i+1)));
      }
    }
  } else { // x3 face
    int fk;
    if (nb.ni.ox3 < 0) fk = fs;
    else               fk = fe;
    for (int n=0; n<nvar; ++n) {
      for (int j=cs, fj=fs; j<=ce; ++j, fj+=2) {
        for (int i=cs, fi=fs; i<=ce; ++i, fi+=2)
          buf[p++] = 0.25*((src(n, fk, j,   i)+src(n, fk, j,   i+1))
                          +(src(n, fk, j+1, i)+src(n, fk, j+1, i+1)));
      }
    }
  }

  return p;
}


//----------------------------------------------------------------------------------------
//! \fn bool MGBoundaryValues::SendMultigridBoundaryBuffers(BoundaryQuantity type)
//  \brief Send boundary buffers

bool MGBoundaryValues::SendMultigridBoundaryBuffers(BoundaryQuantity type) {
  int mylevel = loc.level;
  bool faceonly = false;
  bool bflag = true;

  if (type == BoundaryQuantity::mggrav_f)
    faceonly = true;
  for (int n=0; n<nneighbor; n++) {
    Multigrid *pmg;
    NeighborBlock& nb = neighbor[n];
    if (faceonly && nb.ni.type > NeighborConnect::face) break;
    if (bdata_.sflag[nb.bufid] == BoundaryStatus::completed) continue;
    int ssize = 0;
    if (nb.snb.rank == Globals::my_rank) {
      pmg = pmy_mg_->pmy_driver_->FindMultigrid(nb.snb.gid);
      if (pmg->pmgbval->bdata_.flag[nb.targetid] != BoundaryStatus::waiting) {
        bflag = false;
        continue;
      }
    }
    if (nb.snb.level == mylevel) {
      ssize = LoadMultigridBoundaryBufferSameLevel(bdata_.send[nb.bufid], nb);
    } else if (nb.snb.level < mylevel) {
      if (type == BoundaryQuantity::mggrav_f)
        ssize = LoadMultigridBoundaryBufferToCoarserGravityMassCons(bdata_.send[nb.bufid],
                                                                    nb);
      else
        ssize = LoadMultigridBoundaryBufferToCoarser(bdata_.send[nb.bufid], nb);
    } else {
      ssize = LoadMultigridBoundaryBufferToFiner(bdata_.send[nb.bufid], nb);
    }
    if (nb.snb.rank == Globals::my_rank) {
      std::memcpy(pmg->pmgbval->bdata_.recv[nb.targetid], bdata_.send[nb.bufid],
                  ssize*sizeof(Real));
      pmg->pmgbval->bdata_.flag[nb.targetid] = BoundaryStatus::arrived;
    }
#ifdef MPI_PARALLEL
    else { // NOLINT
      int tag = CreateBvalsMPITag(nb.snb.lid, nb.targetid,
                                  pmy_mg_->pmy_driver_->mg_phys_id_);
      MPI_Isend(bdata_.send[nb.bufid], ssize, MPI_ATHENA_REAL, nb.snb.rank, tag,
                mgcomm_, &(bdata_.req_send[nb.bufid]));
    }
#endif
    bdata_.sflag[nb.bufid] = BoundaryStatus::completed;
  }

  return bflag;
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::SetMultigridBoundarySameLevel(Real *buf,
//                                                           const NeighborBlock& nb)
//  \brief Set hydro boundary received from a block on the same level

void MGBoundaryValues::SetMultigridBoundarySameLevel(Real *buf, const NeighborBlock& nb) {
  AthenaArray<Real> &dst = pmy_mg_->GetCurrentData();
  int nc = pmy_mg_->GetCurrentNumberOfCells();
  int ngh = pmy_mg_->ngh_, nvar = pmy_mg_->nvar_;
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

  if (pmy_mg_->pmy_driver_->pmy_mesh_->multilevel) {
    int cng = ngh;
    int cs = cng, ce = cng + nc/2 -1;

    if (nb.ni.ox1 == 0)
      si = cs, ei = ce;
    else if (nb.ni.ox1 > 0)
      si = ce + 1,   ei = ce + cng;
    else
      si = cs - cng, ei = cs - 1;
    if (nb.ni.ox2 == 0)
      sj = cs, ej = ce;
    else if (nb.ni.ox2 > 0)
      sj = ce + 1,   ej = ce + cng;
    else
      sj = cs - cng, ej = cs - 1;
    if (nb.ni.ox3 == 0)
      sk = cs, ek = ce;
    else if (nb.ni.ox3 > 0)
      sk = ce + 1,   ek = ce + cng;
    else
      sk = cs - cng, ek = cs - 1;

    BufferUtility::UnpackData(buf, cbuf_, 0, nvar-1, si, ei, sj, ej, sk, ek, p);
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::SetMultigridBoundaryFromCoarser(Real *buf,
//                                                             const NeighborBlock& nb)
//  \brief Set hydro boundary received from a block on the same level

void MGBoundaryValues::SetMultigridBoundaryFromCoarser(Real *buf,
                                                       const NeighborBlock& nb) {
  int nc = pmy_mg_->GetCurrentNumberOfCells();
  int ngh = pmy_mg_->ngh_, nvar = pmy_mg_->nvar_;
  int si, sj, sk, ei, ej, ek;
  int cng = 1;
  int cs = cng, ce = cng + nc/2 -1;

  if (nb.ni.ox1 == 0) {
    si = cs, ei = ce;
    if ((loc.lx1 & 1LL) == 0LL) ei += cng;
    else                        si -= cng;
  } else if (nb.ni.ox1 > 0)  {
    si = ce + 1,   ei = ce + cng;
  } else {
    si = cs - cng, ei = cs - 1;
  }
  if (nb.ni.ox2 == 0) {
    sj = cs, ej = ce;
    if ((loc.lx2 & 1LL) == 0LL) ej += cng;
    else                        sj -= cng;
  } else if (nb.ni.ox2 > 0) {
    sj = ce + 1,   ej = ce + cng;
  } else {
    sj = cs - cng, ej = cs - 1;
  }
  if (nb.ni.ox3 == 0) {
    sk = cs, ek = ce;
    if ((loc.lx3 & 1LL) == 0LL) ek += cng;
    else                        sk -= cng;
  } else if (nb.ni.ox3 > 0)  {
    sk = ce + 1,   ek = ce + cng;
  } else {
    sk = cs - cng, ek = cs - 1;
  }

  int p = 0;
  BufferUtility::UnpackData(buf, cbuf_, 0, nvar-1, si, ei, sj, ej, sk, ek, p);

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::SetMultigridBoundaryFromFiner(Real *buf,
//                                                           const NeighborBlock& nb)
//  \brief Set hydro boundary received from a block on the same level

void MGBoundaryValues::SetMultigridBoundaryFromFiner(Real *buf, const NeighborBlock& nb) {
  AthenaArray<Real> &dst = pmy_mg_->GetCurrentData();
  int nc = pmy_mg_->GetCurrentNumberOfCells();
  int ngh = pmy_mg_->ngh_, nvar = pmy_mg_->nvar_;
  int fs = ngh, fe = fs + nc - 1;
  int si, sj, sk, ei, ej, ek;

  if (nb.ni.ox1 == 0) {
    si = fs, ei = fe;
    if (nb.ni.fi1 == 1)   si += nc/2;
    else                  ei -= nc/2;
  } else if (nb.ni.ox1 > 0) {
    si = fe + 1,   ei = fe + ngh;
  } else {
    si = fs - ngh, ei = fs - 1;
  }
  if (nb.ni.ox2 == 0) {
    sj = fs, ej = fe;
    if (nb.ni.ox1 != 0) {
      if (nb.ni.fi1 == 1) sj += nc/2;
      else                ej -= nc/2;
    } else {
      if (nb.ni.fi2 == 1) sj += nc/2;
      else                ej -= nc/2;
    }
  } else if (nb.ni.ox2 > 0) {
    sj = fe + 1,   ej = fe + ngh;
  } else {
    sj = fs - ngh, ej = fs - 1;
  }
  if (nb.ni.ox3 == 0) {
    sk = fs, ek = fe;
    if (nb.ni.ox1 != 0 && nb.ni.ox2 != 0) {
      if (nb.ni.fi1 == 1) sk += nc/2;
      else                ek -= nc/2;
    } else {
      if (nb.ni.fi2 == 1) sk += nc/2;
      else                ek -= nc/2;
    }
  } else if (nb.ni.ox3 > 0) {
    sk = fe + 1,   ek = fe + ngh;
  } else {
    sk = fs - ngh, ek = fs - 1;
  }

  int p = 0;
  BufferUtility::UnpackData(buf, dst, 0, nvar-1, si, ei, sj, ej, sk, ek, p);

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::SetMultigridBoundaryFromCoarserGravityMassCons(
//                                                     Real *buf, const NeighborBlock& nb)
//  \brief Set hydro boundary received from a block on the same level

void MGBoundaryValues::SetMultigridBoundaryFromCoarserGravityMassCons(Real *buf,
                                                             const NeighborBlock& nb) {
  AthenaArray<Real> &dst = pmy_mg_->GetCurrentData();
  int nc = pmy_mg_->GetCurrentNumberOfCells();
  int ngh = pmy_mg_->ngh_, nvar = pmy_mg_->nvar_;
  int fs = ngh, fe = fs + nc - 1;
  int p = 0;
  constexpr Real itw = 1.0/12.0;
  const AthenaArray<Real> &u = dst;

  if (nb.ni.ox1 != 0) { // x1 face
    int i, ig;
    if (nb.ni.ox1 == 1) i = fe, ig = i+1;
    else                i = fs, ig = i-1;
    for (int n=0; n<nvar; ++n) {
      for (int k=fs; k<=fe; k+=2) {
        for (int j=fs; j<=fe; j+=2) {
          Real cg = itw * (8.0 * buf[p++] + ((u(n, k,   j, i) + u(n, k,   j+1, i))
                                          +  (u(n, k+1, j, i) + u(n, k+1, j+1, i))));
          Real qdy = 0.25 * ((u(n, k+1, j+1, i) - u(n, k+1, j,   i))
                           + (u(n, k,   j+1, i) - u(n, k,   j,   i)));
          Real qdz = 0.25 * ((u(n, k+1, j+1, i) - u(n, k,   j+1, i))
                           + (u(n, k+1, j,   i) - u(n, k,   j,   i)));
          dst(n, k,   j,   ig) = cg - qdy - qdz;
          dst(n, k,   j+1, ig) = cg + qdy - qdz;
          dst(n, k+1, j,   ig) = cg - qdy + qdz;
          dst(n, k+1, j+1, ig) = cg + qdy + qdz;
        }
      }
    }
  } else if (nb.ni.ox2 != 0) { // x2 face
    int j, jg;
    if (nb.ni.ox2 == 1) j = fe, jg = j+1;
    else                j = fs, jg = j-1;
    for (int n=0; n<nvar; ++n) {
      for (int k=fs; k<=fe; k+=2) {
        for (int i=fs; i<=fe; i+=2) {
          Real cg = itw * (8.0 * buf[p++] + ((u(n, k,   j, i) + u(n, k,   j, i+1))
                                          +  (u(n, k+1, j, i) + u(n, k+1, j, i+1))));
          Real qdx = 0.25 * ((u(n, k+1, j, i+1) - u(n, k+1, j,   i))
                           + (u(n, k,   j, i+1) - u(n, k,   j,   i)));
          Real qdz = 0.25 * ((u(n, k+1, j, i+1) - u(n, k,   j,   i+1))
                           + (u(n, k+1, j, i  ) - u(n, k,   j,   i)));
          dst(n, k,   jg, i  ) = cg - qdx - qdz;
          dst(n, k,   jg, i+1) = cg + qdx - qdz;
          dst(n, k+1, jg, i  ) = cg - qdx + qdz;
          dst(n, k+1, jg, i+1) = cg + qdx + qdz;
        }
      }
    }
  } else { // x3 face
    int k, kg;
    if (nb.ni.ox3 == 1) k = fe, kg = k+1;
    else                k = fs, kg = k-1;
    for (int n=0; n<nvar; ++n) {
      for (int j=fs; j<=fe; j+=2) {
        for (int i=fs; i<=fe; i+=2) {
          Real cg = itw * (8.0 * buf[p++] + ((u(n, k, j,   i) + u(n, k, j,   i+1))
                                          +  (u(n, k, j+1, i) + u(n, k, j+1, i+1))));
          Real qdx = 0.25 * ((u(n, k, j+1, i+1) - u(n, k, j+1, i))
                           + (u(n, k, j,   i+1) - u(n, k, j,   i)));
          Real qdy = 0.25 * ((u(n, k, j+1, i+1) - u(n, k, j,   i+1))
                           + (u(n, k, j+1, i  ) - u(n, k, j,   i)));
          dst(n, kg, j,   i  ) = cg - qdx - qdy;
          dst(n, kg, j,   i+1) = cg + qdx - qdy;
          dst(n, kg, j+1, i  ) = cg - qdx + qdy;
          dst(n, kg, j+1, i+1) = cg + qdx + qdy;
        }
      }
    }
  }

  return;
}



//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::SetMultigridBoundaryFromFinerGravityMassCons(Real *buf,
//                                                          const NeighborBlock& nb);
void MGBoundaryValues::SetMultigridBoundaryFromFinerGravityMassCons(Real *buf,
                                                           const NeighborBlock& nb) {
  AthenaArray<Real> &dst = pmy_mg_->GetCurrentData();
  int nc = pmy_mg_->GetCurrentNumberOfCells();
  int ngh = pmy_mg_->ngh_, nvar = pmy_mg_->nvar_;
  int fs = ngh, fe = fs + nc - 1;
  int si, sj, sk, ei, ej, ek;
  int oi, oj, ok;
  oi = oj = ok = 0;

  // could reuse the index calculation for normal boundaries
  // but in general this depends on physics
  if (nb.ni.ox1 == 0) {
    si = fs, ei = fe;
    if (nb.ni.fi1 == 1)   si += nc/2;
    else                  ei -= nc/2;
  } else if (nb.ni.ox1 > 0) {
    si = fe + 1,   ei = fe + ngh, oi = -1;
  } else {
    si = fs - ngh, ei = fs - 1,   oi = 1;
  }
  if (nb.ni.ox2 == 0) {
    sj = fs, ej = fe;
    if (nb.ni.ox1 != 0) {
      if (nb.ni.fi1 == 1) sj += nc/2;
      else                ej -= nc/2;
    } else {
      if (nb.ni.fi2 == 1) sj += nc/2;
      else                ej -= nc/2;
    }
  } else if (nb.ni.ox2 > 0) {
    sj = fe + 1,   ej = fe + ngh, oj = -1;
  } else {
    sj = fs - ngh, ej = fs - 1,   oj = 1;
  }
  if (nb.ni.ox3 == 0) {
    sk = fs, ek = fe;
    if (nb.ni.ox1 != 0 && nb.ni.ox2 != 0) {
      if (nb.ni.fi1 == 1) sk += nc/2;
      else                ek -= nc/2;
    } else {
      if (nb.ni.fi2 == 1) sk += nc/2;
      else                ek -= nc/2;
    }
  } else if (nb.ni.ox3 > 0) {
    sk = fe + 1,   ek = fe + ngh, ok = -1;
  } else {
    sk = fs - ngh, ek = fs - 1,   ok = 1;
  }

  constexpr Real ot = 1.0/3.0;
  int p = 0;

  // correct the ghost values using the mass conservation formula
  for (int n=0; n<nvar; ++n) {
    for (int k=sk; k<=ek; ++k) {
      for (int j=sj; j<=ej; ++j) {
        for (int i=si; i<=ei; ++i)
          dst(n,k,j,i) = ot * (4.0*buf[p++] - dst(n,k+ok,j+oj,i+oi));
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn bool MGBoundaryValues::ReceiveMultigridBoundaryBuffers(BoundaryQuantity type)
//  \brief receive the boundary data

bool MGBoundaryValues::ReceiveMultigridBoundaryBuffers(BoundaryQuantity type) {
  bool bflag = true, faceonly = false;

  if (type == BoundaryQuantity::mggrav_f)
    faceonly = true;

  for (int n=0; n<nneighbor; n++) {
    NeighborBlock& nb = neighbor[n];
    if (faceonly && nb.ni.type>NeighborConnect::face) break;
    if (bdata_.flag[nb.bufid] == BoundaryStatus::completed) continue;
    if (bdata_.flag[nb.bufid] == BoundaryStatus::waiting) {
      if (nb.snb.rank == Globals::my_rank) {// on the same process
        bflag = false;
        continue;
      }
#ifdef MPI_PARALLEL
      else { // NOLINT
        int test;
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, mgcomm_, &test, MPI_STATUS_IGNORE);
        MPI_Test(&(bdata_.req_recv[nb.bufid]), &test, MPI_STATUS_IGNORE);
        if (!static_cast<bool>(test)) {
          bflag = false;
          continue;
        }
        bdata_.flag[nb.bufid] = BoundaryStatus::arrived;
      }
#endif
    }
    if (nb.snb.level == loc.level) {
      SetMultigridBoundarySameLevel(bdata_.recv[nb.bufid], nb);
    } else if (nb.snb.level<loc.level) { 
      if (type == BoundaryQuantity::mggrav_f)
        SetMultigridBoundaryFromCoarserGravityMassCons(bdata_.recv[nb.bufid], nb);
      else
        SetMultigridBoundaryFromCoarser(bdata_.recv[nb.bufid], nb);
    } else {
      if (type == BoundaryQuantity::mggrav_f)
        SetMultigridBoundaryFromFinerGravityMassCons(bdata_.recv[nb.bufid], nb);
      else
        SetMultigridBoundaryFromFiner(bdata_.recv[nb.bufid], nb);
    }
    bdata_.flag[nb.bufid] = BoundaryStatus::completed; // completed
  }

  return bflag;
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::ProlongateMultigridBoundaries()
//  \brief prolongate boundaries for Multigrid

void MGBoundaryValues::ProlongateMultigridBoundaries() {
  const int &mylevel = loc.level;
  AthenaArray<Real> &dst = pmy_mg_->GetCurrentData();
  const AthenaArray<Real> &src = dst;
  int nc = pmy_mg_->GetCurrentNumberOfCells();
  int ngh = pmy_mg_->ngh_, nvar = pmy_mg_->nvar_;
  Real time = pmy_mesh_->time;
  int ll = pmy_mg_->nlevel_ - 1 - pmy_mg_->current_level_;
  Real dx = pmy_mg_->rdx_*static_cast<Real>(1<<(ll-1));
  Real dy = pmy_mg_->rdy_*static_cast<Real>(1<<(ll-1));
  Real dz = pmy_mg_->rdz_*static_cast<Real>(1<<(ll-1));
  Real x0 = block_size_.x1min - (static_cast<Real>(ngh) + 0.5)*dx;
  Real y0 = block_size_.x2min - (static_cast<Real>(ngh) + 0.5)*dy;
  Real z0 = block_size_.x3min - (static_cast<Real>(ngh) + 0.5)*dz;
  const int cn = 1, cs = cn, ce = cs + nc/2 -1;

  for (int n=0; n<nneighbor; n++) {
    NeighborBlock& nb = neighbor[n];
    if (nb.snb.level >= mylevel) continue;

    // calculate the loop limits for the ghost zones
    int si, ei, sj, ej, sk, ek;
    if (nb.ni.ox1 == 0) {
      si = cs, ei = ce;
      if ((loc.lx1 & 1LL) == 0LL) ei += cn;
      else                        si -= cn;
    } else if (nb.ni.ox1 > 0) si = ce + 1,  ei = ce + cn;
    else                      si = cs-cn,   ei = cs-1;
    if (nb.ni.ox2 == 0) {
      sj = cs, ej = ce;
      if ((loc.lx2 & 1LL) == 0LL) ej += cn;
      else                        sj -= cn;
    } else if (nb.ni.ox2 > 0) sj = ce + 1,  ej = ce + cn;
    else                      sj = cs-cn,   ej = cs-1;
    if (nb.ni.ox3 == 0) {
      sk = cs, ek = ce;
      if ((loc.lx3 & 1LL) == 0LL) ek += cn;
      else                        sk -= cn;
    } else if (nb.ni.ox3 > 0) sk = ce + 1,  ek = ce + cn;
    else                      sk = cs-cn,   ek = cs-1;

    // Apply physical boundaries on the coarse boundary:
    if (nb.ni.ox1 == 0) {
      if (MGBoundaryFunction_[BoundaryFace::inner_x1] != nullptr)
        MGBoundaryFunction_[BoundaryFace::inner_x1](cbuf_, time, nvar,
                            cs, ce, sj, ej, sk, ek, ngh, x0, y0, z0, dx, dy, dz);
      if (MGBoundaryFunction_[BoundaryFace::outer_x1] != nullptr)
        MGBoundaryFunction_[BoundaryFace::outer_x1](cbuf_, time, nvar,
                            cs, ce, sj, ej, sk, ek, ngh, x0, y0, z0, dx, dy, dz);
    }
    if (nb.ni.ox2 == 0) {
      if (MGBoundaryFunction_[BoundaryFace::inner_x2] != nullptr)
        MGBoundaryFunction_[BoundaryFace::inner_x2](cbuf_, time, nvar,
                            cs, ce, sj, ej, sk, ek, ngh, x0, y0, z0, dx, dy, dz);
      if (MGBoundaryFunction_[BoundaryFace::outer_x2] != nullptr)
        MGBoundaryFunction_[BoundaryFace::outer_x2](cbuf_, time, nvar,
                            cs, ce, sj, ej, sk, ek, ngh, x0, y0, z0, dx, dy, dz);
    }
    if (nb.ni.ox3 == 0) {
      if (MGBoundaryFunction_[BoundaryFace::inner_x3] != nullptr)
        MGBoundaryFunction_[BoundaryFace::inner_x3](cbuf_, time, nvar,
                            cs, ce, sj, ej, sk, ek, ngh, x0, y0, z0, dx, dy, dz);
      if (MGBoundaryFunction_[BoundaryFace::outer_x3] != nullptr)
        MGBoundaryFunction_[BoundaryFace::outer_x3](cbuf_, time, nvar,
                            cs, ce, sj, ej, sk, ek, ngh, x0, y0, z0, dx, dy, dz);
    }

    // Prolongation using tri-linear interpolation
    int fsi = (si - ngh)*2 + ngh;
    int fsj = (sj - ngh)*2 + ngh;
    int fsk = (sk - ngh)*2 + ngh;
    for (int n=0; n<nvar; ++n) {
      for (int k=sk, fk=fsk; k<=ek; k++, fk+=2) {
        for (int j=sj, fj=fsj; j<=ej; j++, fj+=2) {
          for (int i=si, fi=fsi; i<=ei; i++, fi+=2) {
            if (nb.ni.ox1 >= 0 && nb.ni.ox2 >= 0 && nb.ni.ox3 >= 0)
              dst(n,fk  ,fj  ,fi  ) =
                0.015625*(27.0*cbuf_(n,k,j,i)+cbuf_(n,k-1,j-1,i-1)
                        +9.0*(cbuf_(n,k,j,i-1)+cbuf_(n,k,j-1,i)+cbuf_(n,k-1,j,i))
                        +3.0*(cbuf_(n,k-1,j-1,i)+cbuf_(n,k-1,j,i-1)+cbuf_(n,k,j-1,i-1)));
            if (nb.ni.ox1 <= 0 && nb.ni.ox2 >= 0 && nb.ni.ox3 >= 0)
              dst(n,fk  ,fj  ,fi+1) =
                0.015625*(27.0*cbuf_(n,k,j,i)+cbuf_(n,k-1,j-1,i+1)
                        +9.0*(cbuf_(n,k,j,i+1)+cbuf_(n,k,j-1,i)+cbuf_(n,k-1,j,i))
                        +3.0*(cbuf_(n,k-1,j-1,i)+cbuf_(n,k-1,j,i+1)+cbuf_(n,k,j-1,i+1)));
            if (nb.ni.ox1 >= 0 && nb.ni.ox2 <= 0 && nb.ni.ox3 >= 0)
              dst(n,fk  ,fj+1,fi  ) =
                0.015625*(27.0*cbuf_(n,k,j,i)+cbuf_(n,k-1,j+1,i-1)
                        +9.0*(cbuf_(n,k,j,i-1)+cbuf_(n,k,j+1,i)+cbuf_(n,k-1,j,i))
                        +3.0*(cbuf_(n,k-1,j+1,i)+cbuf_(n,k-1,j,i-1)+cbuf_(n,k,j+1,i-1)));
            if (nb.ni.ox1 >= 0 && nb.ni.ox2 >= 0 && nb.ni.ox3 <= 0)
              dst(n,fk+1,fj  ,fi  ) =
                0.015625*(27.0*cbuf_(n,k,j,i)+cbuf_(n,k+1,j-1,i-1)
                        +9.0*(cbuf_(n,k,j,i-1)+cbuf_(n,k,j-1,i)+cbuf_(n,k+1,j,i))
                        +3.0*(cbuf_(n,k+1,j-1,i)+cbuf_(n,k+1,j,i-1)+cbuf_(n,k,j-1,i-1)));
            if (nb.ni.ox1 >= 0 && nb.ni.ox2 <= 0 && nb.ni.ox3 <= 0)
              dst(n,fk+1,fj+1,fi  ) =
                0.015625*(27.0*cbuf_(n,k,j,i)+cbuf_(n,k+1,j+1,i-1)
                        +9.0*(cbuf_(n,k,j,i-1)+cbuf_(n,k,j+1,i)+cbuf_(n,k+1,j,i))
                        +3.0*(cbuf_(n,k+1,j+1,i)+cbuf_(n,k+1,j,i-1)+cbuf_(n,k,j+1,i-1)));
            if (nb.ni.ox1 <= 0 && nb.ni.ox2 >= 0 && nb.ni.ox3 <= 0)
              dst(n,fk+1,fj  ,fi+1) =
                0.015625*(27.0*cbuf_(n,k,j,i)+cbuf_(n,k+1,j-1,i+1)
                        +9.0*(cbuf_(n,k,j,i+1)+cbuf_(n,k,j-1,i)+cbuf_(n,k+1,j,i))
                        +3.0*(cbuf_(n,k+1,j-1,i)+cbuf_(n,k+1,j,i+1)+cbuf_(n,k,j-1,i+1)));
            if (nb.ni.ox1 <= 0 && nb.ni.ox2 <= 0 && nb.ni.ox3 >= 0)
              dst(n,fk  ,fj+1,fi+1) =
                0.015625*(27.0*cbuf_(n,k,j,i)+cbuf_(n,k-1,j+1,i+1)
                        +9.0*(cbuf_(n,k,j,i+1)+cbuf_(n,k,j+1,i)+cbuf_(n,k-1,j,i))
                        +3.0*(cbuf_(n,k-1,j+1,i)+cbuf_(n,k-1,j,i+1)+cbuf_(n,k,j+1,i+1)));
            if (nb.ni.ox1 <= 0 && nb.ni.ox2 <= 0 && nb.ni.ox3 <= 0)
              dst(n,fk+1,fj+1,fi+1) =
                0.015625*(27.0*cbuf_(n,k,j,i)+cbuf_(n,k+1,j+1,i+1)
                        +9.0*(cbuf_(n,k,j,i+1)+cbuf_(n,k,j+1,i)+cbuf_(n,k+1,j,i))
                        +3.0*(cbuf_(n,k+1,j+1,i)+cbuf_(n,k+1,j,i+1)+cbuf_(n,k,j+1,i+1)));
            }
        }
      }
    }
  } // end loop over nneighbor

  return;
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

  return;
}

