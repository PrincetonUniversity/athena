//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_mg.cpp
//! \brief

// C headers

// C++ headers
#include <algorithm>
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
//! \brief Constructor of the MGBoundaryValues class

MGBoundaryValues::MGBoundaryValues(Multigrid *pmg, BoundaryFlag *input_bcs)
    : BoundaryBase(pmg->pmy_driver_->pmy_mesh_, pmg->loc_, pmg->size_, input_bcs),
      pmy_mg_(pmg), bcolor_(0), triplebuf_(true) {
#ifdef MPI_PARALLEL
  mgcomm_ = pmg->pmy_driver_->MPI_COMM_MULTIGRID;
#endif
  for (int i = 0; i < 6; ++i)
    MGBoundaryFunction_[i] = pmg->pmy_driver_->MGBoundaryFunction_[i];
  for (int i = 0; i < 6; ++i)
    MGCoeffBoundaryFunction_[i] = pmg->pmy_driver_->MGCoeffBoundaryFunction_[i];
  if (pmy_mg_->pmy_block_ != nullptr) {
    InitBoundaryData(BoundaryQuantity::mg);
    int nc = block_size_.nx1 + 2*pmy_mg_->ngh_;
    int nv = std::max(pmy_mg_->nvar_, pmy_mg_->ncoeff_);
    cbuf_.NewAthenaArray(nv, nc, nc, nc);
    cbufold_.NewAthenaArray(nv, nc, nc, nc);
    for (int i = 0; i < 6; ++i) {
      if (block_bcs[i] != BoundaryFlag::block && block_bcs[i] != BoundaryFlag::periodic)
        apply_bndry_fn_[i] = true;
    }
  } else {
    for (int i = 0; i < 6; ++i)
      apply_bndry_fn_[i] = true;
  }
}


//----------------------------------------------------------------------------------------
//! \fn MGBoundaryValues::~MGBoundaryValues()
//! \brief Destructor of the MGBoundaryValues class

MGBoundaryValues::~MGBoundaryValues() {
  if (pmy_mg_->pmy_block_ != nullptr)
    DestroyBoundaryData();
}

//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::InitBoundaryData(BoundaryQuantity type)
//! \brief Initialize BoundaryData<> structure

void MGBoundaryValues::InitBoundaryData(BoundaryQuantity type) {
  int nbuf = triplebuf_?3:1;
  for (int c = 0; c < nbuf; ++c) {
    bdata_[c].nbmax = maxneighbor_;
    for (int n = 0; n < bdata_[c].nbmax; ++n) {
      // Clear flags and requests
      bdata_[c].flag[n] = BoundaryStatus::waiting;
      bdata_[c].sflag[n] = BoundaryStatus::waiting;
      bdata_[c].send[n] = nullptr;
      bdata_[c].recv[n] = nullptr;
#ifdef MPI_PARALLEL
      bdata_[c].req_send[n] = MPI_REQUEST_NULL;
      bdata_[c].req_recv[n] = MPI_REQUEST_NULL;
#endif

      // Allocate buffers
      // calculate the buffer size
      int ngh = pmy_mg_->ngh_;
      int nc = block_size_.nx1;
      int size = 0;
      if (BoundaryValues::ni[n].type == NeighborConnect::face)
        size = SQR(nc)*ngh;
      else if (BoundaryValues::ni[n].type == NeighborConnect::edge)
        size = nc*ngh*ngh;
      else if (BoundaryValues::ni[n].type == NeighborConnect::corner)
        size = ngh*ngh*ngh;
      if (pmy_mg_->pmy_driver_->pmy_mesh_->multilevel) {
        if (BoundaryValues::ni[n].type == NeighborConnect::face)
          size += SQR(nc/2)*ngh;
        else if (BoundaryValues::ni[n].type == NeighborConnect::edge)
          size += nc/2*ngh*ngh;
        else if (BoundaryValues::ni[n].type == NeighborConnect::corner)
          size += ngh*ngh*ngh;
      }

      size *= std::max(pmy_mg_->nvar_*(1+pmy_mg_->pmy_driver_->ffas_), pmy_mg_->ncoeff_);

      bdata_[c].send[n] = new Real[size];
      bdata_[c].recv[n] = new Real[size];
    }
  }
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::DestroyBoundaryData()
//! \brief Destroy BoundaryData<> structure

void MGBoundaryValues::DestroyBoundaryData() {
  int nbuf = triplebuf_?3:1;
  for (int c = 0; c < nbuf; ++c) {
    for (int n = 0; n < bdata_[c].nbmax; ++n) {
      delete [] bdata_[c].send[n];
      delete [] bdata_[c].recv[n];
#ifdef MPI_PARALLEL
      if (bdata_[c].req_send[n] != MPI_REQUEST_NULL)
        MPI_Request_free(&bdata_[c].req_send[n]);
      if (bdata_[c].req_recv[n] != MPI_REQUEST_NULL)
        MPI_Request_free(&bdata_[c].req_recv[n]);
#endif
    }
  }
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::DispatchBoundaryFunction(BoundaryFace face,
//                           AthenaArray<Real> &dst, Real time, int nvar,
//                           int is, int ie, int js, int je, int ks, int ke, int ngh,
//                           const MGCoordinates &coord, bool fcoeff)
//  \brief Call multigrid boundary function for a face

void MGBoundaryValues::DispatchBoundaryFunction(BoundaryFace face, AthenaArray<Real> &dst,
     Real time, int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
     const MGCoordinates &coord, bool fcoeff) {
  if (!fcoeff) {
    AthenaArray<Real> &mpcoeff = pmy_mg_->pmy_driver_->mpcoeff_[0];
    AthenaArray<Real> &mpo = pmy_mg_->pmy_driver_->mpo_;
    const int& mporder = pmy_mg_->pmy_driver_->mporder_;
    switch (face) {
      case BoundaryFace::inner_x1:
        switch (block_bcs[BoundaryFace::inner_x1]) {
          case BoundaryFlag::user:
            MGBoundaryFunction_[BoundaryFace::inner_x1](dst, time, nvar,
                                is, ie, js, je, ks, ke, ngh, coord);
            break;
          case BoundaryFlag::periodic:
            MGPeriodicInnerX1(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
            break;
          case BoundaryFlag::mg_zerograd:
            MGZeroGradientInnerX1(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
            break;
          case BoundaryFlag::mg_zerofixed:
            MGZeroFixedInnerX1(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
            break;
          case BoundaryFlag::mg_multipole:
            MGMultipoleInnerX1(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord,
                               mpcoeff, mpo, mporder);
            break;
          default:
            break;
        }
        break;
      case BoundaryFace::outer_x1:
        switch (block_bcs[BoundaryFace::outer_x1]) {
          case BoundaryFlag::user:
            MGBoundaryFunction_[BoundaryFace::outer_x1](dst, time, nvar,
                                is, ie, js, je, ks, ke, ngh, coord);
            break;
          case BoundaryFlag::periodic:
            MGPeriodicOuterX1(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
            break;
          case BoundaryFlag::mg_zerograd:
            MGZeroGradientOuterX1(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
            break;
          case BoundaryFlag::mg_zerofixed:
            MGZeroFixedOuterX1(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
            break;
          case BoundaryFlag::mg_multipole:
            MGMultipoleOuterX1(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord,
                               mpcoeff, mpo, mporder);
            break;
          default:
            break;
        }
        break;
      case BoundaryFace::inner_x2:
        switch (block_bcs[BoundaryFace::inner_x2]) {
          case BoundaryFlag::user:
            MGBoundaryFunction_[BoundaryFace::inner_x2](dst, time, nvar,
                                is, ie, js, je, ks, ke, ngh, coord);
            break;
          case BoundaryFlag::periodic:
            MGPeriodicInnerX2(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
            break;
          case BoundaryFlag::mg_zerograd:
            MGZeroGradientInnerX2(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
            break;
          case BoundaryFlag::mg_zerofixed:
            MGZeroFixedInnerX2(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
            break;
          case BoundaryFlag::mg_multipole:
            MGMultipoleInnerX2(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord,
                               mpcoeff, mpo, mporder);
            break;
          default:
            break;
        }
        break;
      case BoundaryFace::outer_x2:
        switch (block_bcs[BoundaryFace::outer_x2]) {
          case BoundaryFlag::user:
            MGBoundaryFunction_[BoundaryFace::outer_x2](dst, time, nvar,
                                is, ie, js, je, ks, ke, ngh, coord);
            break;
          case BoundaryFlag::periodic:
            MGPeriodicOuterX2(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
            break;
          case BoundaryFlag::mg_zerograd:
            MGZeroGradientOuterX2(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
            break;
          case BoundaryFlag::mg_zerofixed:
            MGZeroFixedOuterX2(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
            break;
          case BoundaryFlag::mg_multipole:
            MGMultipoleOuterX2(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord,
                               mpcoeff, mpo, mporder);
            break;
          default:
            break;
        }
        break;
      case BoundaryFace::inner_x3:
        switch (block_bcs[BoundaryFace::inner_x3]) {
          case BoundaryFlag::user:
            MGBoundaryFunction_[BoundaryFace::inner_x3](dst, time, nvar,
                                is, ie, js, je, ks, ke, ngh, coord);
            break;
          case BoundaryFlag::periodic:
            MGPeriodicInnerX3(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
            break;
          case BoundaryFlag::mg_zerograd:
            MGZeroGradientInnerX3(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
            break;
          case BoundaryFlag::mg_zerofixed:
            MGZeroFixedInnerX3(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
            break;
          case BoundaryFlag::mg_multipole:
            MGMultipoleInnerX3(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord,
                               mpcoeff, mpo, mporder);
            break;
          default:
            break;
        }
        break;
      case BoundaryFace::outer_x3:
        switch (block_bcs[BoundaryFace::outer_x3]) {
          case BoundaryFlag::user:
            MGBoundaryFunction_[BoundaryFace::outer_x3](dst, time, nvar,
                                is, ie, js, je, ks, ke, ngh, coord);
            break;
          case BoundaryFlag::periodic:
            MGPeriodicOuterX3(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
            break;
          case BoundaryFlag::mg_zerograd:
            MGZeroGradientOuterX3(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
            break;
          case BoundaryFlag::mg_zerofixed:
            MGZeroFixedOuterX3(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
            break;
          case BoundaryFlag::mg_multipole:
            MGMultipoleOuterX3(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord,
                               mpcoeff, mpo, mporder);
            break;
          default:
            break;
        }
        break;
      default:
        break;
    }
  } else { // coeffient
    switch (face) {
      case BoundaryFace::inner_x1:
        if (block_bcs[face] == BoundaryFlag::periodic)
          MGPeriodicInnerX1(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
        else if (MGCoeffBoundaryFunction_[face] != nullptr)
          MGCoeffBoundaryFunction_[face](dst, time, nvar, is, ie, js, je, ks, ke,
                                         ngh, coord);
        else // default: zero gradient
          MGZeroGradientInnerX1(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
        break;
      case BoundaryFace::outer_x1:
        if (block_bcs[face] == BoundaryFlag::periodic)
          MGPeriodicOuterX1(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
        else if (MGCoeffBoundaryFunction_[face] != nullptr)
          MGCoeffBoundaryFunction_[face](dst, time, nvar, is, ie, js, je, ks, ke,
                                         ngh, coord);
        else // default: zero gradient
          MGZeroGradientOuterX1(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
        break;
      case BoundaryFace::inner_x2:
        if (block_bcs[face] == BoundaryFlag::periodic)
          MGPeriodicInnerX2(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
        else if (MGCoeffBoundaryFunction_[face] != nullptr)
          MGCoeffBoundaryFunction_[face](dst, time, nvar, is, ie, js, je, ks, ke,
                                         ngh, coord);
        else // default: zero gradient
          MGZeroGradientInnerX2(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
        break;
      case BoundaryFace::outer_x2:
        if (block_bcs[face] == BoundaryFlag::periodic)
          MGPeriodicOuterX2(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
        else if (MGCoeffBoundaryFunction_[face] != nullptr)
          MGCoeffBoundaryFunction_[face](dst, time, nvar, is, ie, js, je, ks, ke,
                                         ngh, coord);
        else // default: zero gradient
          MGZeroGradientOuterX2(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
        break;
      case BoundaryFace::inner_x3:
        if (block_bcs[face] == BoundaryFlag::periodic)
          MGPeriodicInnerX3(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
        else if (MGCoeffBoundaryFunction_[face] != nullptr)
          MGCoeffBoundaryFunction_[face](dst, time, nvar, is, ie, js, je, ks, ke,
                                         ngh, coord);
        else // default: zero gradient
          MGZeroGradientInnerX3(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
        break;
      case BoundaryFace::outer_x3:
        if (block_bcs[face] == BoundaryFlag::periodic)
          MGPeriodicOuterX3(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
        else if (MGCoeffBoundaryFunction_[face] != nullptr)
          MGCoeffBoundaryFunction_[face](dst, time, nvar, is, ie, js, je, ks, ke,
                                         ngh, coord);
        else // default: zero gradient
          MGZeroGradientOuterX3(dst, time, nvar, is, ie, js, je, ks, ke, ngh, coord);
        break;
      default:
        break;
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::ApplyPhysicalBoundaries(int flag, bool fcoeff)
//! \brief Apply physical boundary conditions to the current Multigrid data

void MGBoundaryValues::ApplyPhysicalBoundaries(int flag, bool fcoeff) {
  AthenaArray<Real> *u;
  MGCoordinates *c;
  int lev = pmy_mg_->GetCurrentLevel();
  int ll = pmy_mg_->nlevel_ - 1 - lev;
  int ngh = pmy_mg_->ngh_, nvar = pmy_mg_->nvar_;
  if (!fcoeff) {
    if (flag == 0)      u = &(pmy_mg_->GetCurrentData()), c = &(pmy_mg_->coord_[lev]);
    else if (flag == 1) u = &cbuf_,                       c = &(pmy_mg_->ccoord_[lev]);
    else                u = &cbufold_,                    c = &(pmy_mg_->ccoord_[lev]);
  } else {
    nvar = pmy_mg_->ncoeff_;
    if (flag == 0)
      u = &(pmy_mg_->GetCurrentCoefficient()), c = &(pmy_mg_->coord_[lev]);
    else
      u = &cbuf_, c=&(pmy_mg_->ccoord_[lev]);
  }
  AthenaArray<Real> &dst = *u;
  MGCoordinates &coord = *c;
  if (flag > 0) ll++;
  int ncx = block_size_.nx1 >> ll, ncy = block_size_.nx2 >> ll,
      ncz = block_size_.nx3 >> ll;
  int is = ngh, ie = ncx + ngh - 1;
  int js = ngh, je = ncy + ngh - 1;
  int ks = ngh, ke = ncz + ngh - 1;
  // cppcheck-suppress knownConditionTrueFalse
  int bis = is - ngh, bie = ie + ngh;
  int bjs = js,       bje = je;
  int bks = ks,       bke = ke;
  Real time = pmy_mesh_->time;

  if (!apply_bndry_fn_[BoundaryFace::inner_x2]) bjs = js - ngh;
  if (!apply_bndry_fn_[BoundaryFace::outer_x2]) bje = je + ngh;
  if (!apply_bndry_fn_[BoundaryFace::inner_x3]) bks = ks - ngh;
  if (!apply_bndry_fn_[BoundaryFace::outer_x3]) bke = ke + ngh;

  if (apply_bndry_fn_[BoundaryFace::inner_x1])
    DispatchBoundaryFunction(BoundaryFace::inner_x1, dst, time, nvar,
                             is, ie, bjs, bje, bks, bke, ngh, coord, fcoeff);
  if (apply_bndry_fn_[BoundaryFace::outer_x1])
    DispatchBoundaryFunction(BoundaryFace::outer_x1, dst, time, nvar,
                             is, ie, bjs, bje, bks, bke, ngh, coord, fcoeff);
  if (apply_bndry_fn_[BoundaryFace::inner_x2])
    DispatchBoundaryFunction(BoundaryFace::inner_x2, dst, time, nvar,
                             bis, bie, js, je, bks, bke, ngh, coord, fcoeff);
  if (apply_bndry_fn_[BoundaryFace::outer_x2])
    DispatchBoundaryFunction(BoundaryFace::outer_x2, dst, time, nvar,
                             bis, bie, js, je, bks, bke, ngh, coord, fcoeff);
  bjs = js - ngh, bje = je + ngh;
  if (apply_bndry_fn_[BoundaryFace::inner_x3])
    DispatchBoundaryFunction(BoundaryFace::inner_x3, dst, time, nvar,
                             bis, bie, bjs, bje, ks, ke, ngh, coord, fcoeff);
  if (apply_bndry_fn_[BoundaryFace::outer_x3])
    DispatchBoundaryFunction(BoundaryFace::outer_x3, dst, time, nvar,
                             bis, bie, bjs, bje, ks, ke, ngh, coord, fcoeff);
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::StartReceivingMultigrid(BoundaryQuantity type,
//!                                                    bool folddata)
//! \brief initiate MPI_Irecv for multigrid

void MGBoundaryValues::StartReceivingMultigrid(BoundaryQuantity type, bool folddata) {
  if (triplebuf_)
    bcolor_ = (bcolor_ + 1) % 3;

#ifdef MPI_PARALLEL
  int nvar = pmy_mg_->nvar_, ngh = pmy_mg_->ngh_, cngh = ngh;
  int nc = pmy_mg_->GetCurrentNumberOfCells();
  if (type == BoundaryQuantity::mg_coeff) nvar = pmy_mg_->ncoeff_;

  for (int n = 0; n < nneighbor; ++n) {
    NeighborBlock& nb = neighbor[n];
    if (nb.snb.rank!=Globals::my_rank) {
      if (type == BoundaryQuantity::mg_faceonly && nb.snb.level > loc.level
        && nb.ni.type != NeighborConnect::face) continue;
      int size = 0;
      if (nb.snb.level == loc.level) { // same
        if (nb.ni.type == NeighborConnect::face) size = SQR(nc)*ngh;
        else if (nb.ni.type == NeighborConnect::edge) size = nc*ngh*ngh;
        else if (nb.ni.type == NeighborConnect::corner) size = ngh*ngh*ngh;
        if (pmy_mg_->pmy_driver_->nreflevel_ > 0) {
          if (nb.ni.type == NeighborConnect::face) size += SQR(nc/2)*ngh;
          else if (nb.ni.type == NeighborConnect::edge) size += nc/2*ngh*ngh;
          else if (nb.ni.type == NeighborConnect::corner) size += ngh*ngh*ngh;
        }
      } else if (nb.snb.level < loc.level) { // coarser
        if (nb.ni.type == NeighborConnect::face) size = SQR(nc/2+cngh)*cngh;
        else if (nb.ni.type == NeighborConnect::edge) size = (nc/2+cngh)*cngh*cngh;
        else if (nb.ni.type == NeighborConnect::corner) size = cngh*cngh*cngh;
      } else { // finer
        if (nb.ni.type == NeighborConnect::face) size = SQR(nc/2)*ngh;
        else if (nb.ni.type == NeighborConnect::edge) size = nc/2*ngh*ngh;
        else if (nb.ni.type == NeighborConnect::corner) size = ngh*ngh*ngh;
      }
      if (folddata) size *= 2;
      size *= nvar;
      int tag = CreateBvalsMPITag(pmy_mg_->pmy_block_->lid, nb.bufid,
                                  pmy_mg_->pmy_driver_->mg_phys_id_);
      MPI_Irecv(bdata_[bcolor_].recv[nb.bufid], size, MPI_ATHENA_REAL, nb.snb.rank, tag,
                mgcomm_, &(bdata_[bcolor_].req_recv[nb.bufid]));
    }
  }
#endif
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::ClearBoundaryMultigrid(BoundaryQuantity type)
//! \brief clean up the boundary flags after each loop for multigrid

void MGBoundaryValues::ClearBoundaryMultigrid(BoundaryQuantity type) {
  for (int n = 0; n < nneighbor; ++n) {
    NeighborBlock& nb = neighbor[n];
    bdata_[bcolor_].flag[nb.bufid] = BoundaryStatus::waiting;
    bdata_[bcolor_].sflag[nb.bufid] = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
    if (nb.snb.rank != Globals::my_rank) {
      if (!(type == BoundaryQuantity::mg_faceonly
        &&  nb.ni.type != NeighborConnect::face && nb.snb.level < loc.level))
        MPI_Wait(&(bdata_[bcolor_].req_send[nb.bufid]),MPI_STATUS_IGNORE);
    }
#endif
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn int MGBoundaryValues::LoadMultigridBoundaryBufferSameLevel(Real *buf,
//!                               const NeighborBlock& nb, bool folddata, bool fcoeff)
//! \brief Set Multigrid boundary buffers for sending to a block on the same level

int MGBoundaryValues::LoadMultigridBoundaryBufferSameLevel(Real *buf,
                                   const NeighborBlock& nb, bool folddata, bool fcoeff) {
  int nc = pmy_mg_->GetCurrentNumberOfCells();
  int ngh = pmy_mg_->ngh_, nvar = pmy_mg_->nvar_;
  int fs = ngh, fe = fs + nc - 1;
  int si = (nb.ni.ox1 > 0) ? (fe - ngh + 1) : fs;
  int ei = (nb.ni.ox1 < 0) ? (fs + ngh - 1) : fe;
  int sj = (nb.ni.ox2 > 0) ? (fe - ngh + 1) : fs;
  int ej = (nb.ni.ox2 < 0) ? (fs + ngh - 1) : fe;
  int sk = (nb.ni.ox3 > 0) ? (fe - ngh + 1) : fs;
  int ek = (nb.ni.ox3 < 0) ? (fs + ngh - 1) : fe;
  AthenaArray<Real> *t;
  if (!fcoeff) {
    t = &(pmy_mg_->GetCurrentData());
  } else {
    nvar = pmy_mg_->ncoeff_;
    t = &(pmy_mg_->GetCurrentCoefficient());
  }
  const AthenaArray<Real> &u = *t;
  const AthenaArray<Real> &old = pmy_mg_->GetCurrentOldData();

  int p = 0;
  BufferUtility::PackData(u, buf, 0, nvar-1, si, ei, sj, ej, sk, ek, p);
  if (folddata)
    BufferUtility::PackData(old, buf, 0, nvar-1, si, ei, sj, ej, sk, ek, p);

  if (pmy_mg_->pmy_driver_->nreflevel_ > 0) {
    int cn = ngh - 1, fs = ngh, cs = ngh, fe = fs + nc - 1, ce = cs + nc/2 - 1;
    si = (nb.ni.ox1 > 0) ? (ce - cn) : cs;
    ei = (nb.ni.ox1 < 0) ? (cs + cn) : ce;
    sj = (nb.ni.ox2 > 0) ? (ce - cn) : cs;
    ej = (nb.ni.ox2 < 0) ? (cs + cn) : ce;
    sk = (nb.ni.ox3 > 0) ? (ce - cn) : cs;
    ek = (nb.ni.ox3 < 0) ? (cs + cn) : ce;
    for (int n=0; n<nvar; ++n) {
      for (int k=sk; k<=ek; ++k) {
        int fk = (k - ngh) * 2 + ngh;
        for (int j=sj; j<=ej; ++j) {
          int fj = (j - ngh) * 2 + ngh;
#pragma ivdep
          for (int i=si; i<=ei; ++i) {
            int fi = (i - ngh) * 2 + ngh;
            buf[p++] = cbuf_(n, k, j, i)
                     = 0.125*(((u(n, fk,   fj,   fi)+u(n, fk,   fj,   fi+1))
                            +  (u(n, fk,   fj+1, fi)+u(n, fk,   fj+1, fi+1)))
                            + ((u(n, fk+1, fj,   fi)+u(n, fk+1, fj,   fi+1))
                            +  (u(n, fk+1, fj+1, fi)+u(n, fk+1, fj+1, fi+1))));
          }
        }
      }
    }
    if (folddata) {
      for (int n=0; n<nvar; ++n) {
        for (int k=sk; k<=ek; ++k) {
          int fk = (k - ngh) * 2 + ngh;
          for (int j=sj; j<=ej; ++j) {
            int fj = (j - ngh) * 2 + ngh;
#pragma ivdep
            for (int i=si; i<=ei; ++i) {
              int fi = (i - ngh) * 2 + ngh;
              buf[p++] = cbufold_(n, k, j, i)
                       = 0.125*(((old(n, fk,   fj,   fi)+old(n, fk,   fj,   fi+1))
                              +  (old(n, fk,   fj+1, fi)+old(n, fk,   fj+1, fi+1)))
                              + ((old(n, fk+1, fj,   fi)+old(n, fk+1, fj,   fi+1))
                              +  (old(n, fk+1, fj+1, fi)+old(n, fk+1, fj+1, fi+1))));
            }
          }
        }
      }
    }
  }

  return p;
}


//----------------------------------------------------------------------------------------
//! \fn int MGBoundaryValues::LoadMultigridBoundaryBufferToCoarser(Real *buf,
//!                               const NeighborBlock& nb, bool folddata, bool fcoeff)
//! \brief Set Multigrid boundary buffers for sending to a block on the coarser level

int MGBoundaryValues::LoadMultigridBoundaryBufferToCoarser(Real *buf,
                                   const NeighborBlock& nb, bool folddata, bool fcoeff) {
  int ngh = pmy_mg_->ngh_, nvar = pmy_mg_->nvar_;
  int nc = pmy_mg_->GetCurrentNumberOfCells();
  int cn = ngh - 1, cs = ngh, ce = cs + nc/2 - 1;
  int p=0;
  int si = (nb.ni.ox1 > 0) ? (ce - cn) : cs;
  int ei = (nb.ni.ox1 < 0) ? (cs + cn) : ce;
  int sj = (nb.ni.ox2 > 0) ? (ce - cn) : cs;
  int ej = (nb.ni.ox2 < 0) ? (cs + cn) : ce;
  int sk = (nb.ni.ox3 > 0) ? (ce - cn) : cs;
  int ek = (nb.ni.ox3 < 0) ? (cs + cn) : ce;
  AthenaArray<Real> *t;
  if (!fcoeff) {
    t = &(pmy_mg_->GetCurrentData());
  } else {
    nvar = pmy_mg_->ncoeff_;
    t = &(pmy_mg_->GetCurrentCoefficient());
  }
  const AthenaArray<Real> &u = *t;
  const AthenaArray<Real> &old = pmy_mg_->GetCurrentOldData();

  // restrict and store in the coarse buffer
  for (int n=0; n<nvar; ++n) {
    for (int k=sk; k<=ek; ++k) {
      int fk = (k - ngh) * 2 + ngh;
      for (int j=sj; j<=ej; ++j) {
        int fj = (j - ngh) * 2 + ngh;
#pragma ivdep
        for (int i=si; i<=ei; ++i) {
          int fi = (i - ngh) * 2 + ngh;
          buf[p++] = cbuf_(n, k, j, i)
                   = 0.125*(((u(n, fk,   fj,   fi)+u(n, fk,   fj,   fi+1))
                          +  (u(n, fk,   fj+1, fi)+u(n, fk,   fj+1, fi+1)))
                          + ((u(n, fk+1, fj,   fi)+u(n, fk+1, fj,   fi+1))
                          +  (u(n, fk+1, fj+1, fi)+u(n, fk+1, fj+1, fi+1))));
        }
      }
    }
  }
  if (folddata) {
    for (int n=0; n<nvar; ++n) {
      for (int k=sk; k<=ek; ++k) {
        int fk = (k - ngh) * 2 + ngh;
        for (int j=sj; j<=ej; ++j) {
          int fj = (j - ngh) * 2 + ngh;
#pragma ivdep
          for (int i=si; i<=ei; ++i) {
            int fi = (i - ngh) * 2 + ngh;
            buf[p++] = cbufold_(n, k, j, i)
                     = 0.125*(((old(n, fk,   fj,   fi)+old(n, fk,   fj,   fi+1))
                            +  (old(n, fk,   fj+1, fi)+old(n, fk,   fj+1, fi+1)))
                            + ((old(n, fk+1, fj,   fi)+old(n, fk+1, fj,   fi+1))
                            +  (old(n, fk+1, fj+1, fi)+old(n, fk+1, fj+1, fi+1))));
          }
        }
      }
    }
  }

  return p;
}


//----------------------------------------------------------------------------------------
//! \fn int MGBoundaryValues::LoadMultigridBoundaryBufferToFiner(Real *buf,
//!                               const NeighborBlock& nb, bool folddata, bool fcoeff)
//! \brief Set Multigrid boundary buffers for sending to a block on the finer level
int MGBoundaryValues::LoadMultigridBoundaryBufferToFiner(Real *buf,
                                   const NeighborBlock& nb, bool folddata, bool fcoeff) {
  int nc = pmy_mg_->GetCurrentNumberOfCells();
  int ngh = pmy_mg_->ngh_, nvar = pmy_mg_->nvar_;
  int cn = ngh - 1, fs = ngh, fe = fs + nc - 1;
  int si = (nb.ni.ox1 > 0) ? (fe - cn) : fs;
  int ei = (nb.ni.ox1 < 0) ? (fs + cn) : fe;
  int sj = (nb.ni.ox2 > 0) ? (fe - cn) : fs;
  int ej = (nb.ni.ox2 < 0) ? (fs + cn) : fe;
  int sk = (nb.ni.ox3 > 0) ? (fe - cn) : fs;
  int ek = (nb.ni.ox3 < 0) ? (fs + cn) : fe;
  AthenaArray<Real> *t;
  if (!fcoeff) {
    t = &(pmy_mg_->GetCurrentData());
  } else {
    nvar = pmy_mg_->ncoeff_;
    t = &(pmy_mg_->GetCurrentCoefficient());
  }
  const AthenaArray<Real> &u = *t;
  const AthenaArray<Real> &old = pmy_mg_->GetCurrentOldData();

  int nred = nc/2 - ngh;
  if (nb.ni.ox1 == 0) {
    if (nb.ni.fi1 == 1)   si += nred;
    else                  ei -= nred;
  }
  if (nb.ni.ox2 == 0) {
    if (nb.ni.ox1 != 0) {
      if (nb.ni.fi1 == 1) sj += nred;
      else                ej -= nred;
    } else {
      if (nb.ni.fi2 == 1) sj += nred;
      else                ej -= nred;
    }
  }
  if (nb.ni.ox3 == 0) {
    if (nb.ni.ox1 != 0 && nb.ni.ox2 != 0) {
      if (nb.ni.fi1 == 1) sk += nred;
      else                ek -= nred;
    } else {
      if (nb.ni.fi2 == 1) sk += nred;
      else                ek -= nred;
    }
  }

  int p = 0;
  BufferUtility::PackData(u, buf, 0, nvar-1, si, ei, sj, ej, sk, ek, p);
  if (folddata)
    BufferUtility::PackData(old, buf, 0, nvar-1, si, ei, sj, ej, sk, ek, p);

  return p;
}


//----------------------------------------------------------------------------------------
//! \fn bool MGBoundaryValues::SendMultigridBoundaryBuffers(BoundaryQuantity type,
//!                                                         bool folddata)
//! \brief Send boundary buffers

bool MGBoundaryValues::SendMultigridBoundaryBuffers(BoundaryQuantity type,
                                                    bool folddata) {
  bool bflag = true;
  bool fcoeff = false;
  if (type == BoundaryQuantity::mg_coeff)
    fcoeff = true;
  Multigrid *pmg = nullptr;
  for (int n = 0; n < nneighbor; ++n) {
    NeighborBlock& nb = neighbor[n];
    if (bdata_[bcolor_].sflag[nb.bufid] == BoundaryStatus::completed) continue;
    if (type == BoundaryQuantity::mg_faceonly && nb.snb.level < loc.level
      && nb.ni.type != NeighborConnect::face) continue;
    int ssize = 0;
    if (nb.snb.rank == Globals::my_rank) {
      pmg = pmy_mg_->pmy_driver_->FindMultigrid(nb.snb.gid);
      if (!triplebuf_ &&
          pmg->pmgbval->bdata_[bcolor_].flag[nb.targetid] != BoundaryStatus::waiting) {
        bflag = false;
        continue;
      }
    }
    if (nb.snb.level == loc.level) {
      ssize = LoadMultigridBoundaryBufferSameLevel(bdata_[bcolor_].send[nb.bufid],
                                                   nb, folddata, fcoeff);
    } else if (nb.snb.level < loc.level) {
      if (type == BoundaryQuantity::mg_faceonly)
        ssize = LoadMultigridBoundaryBufferToCoarserFluxCons(
                                                    bdata_[bcolor_].send[nb.bufid], nb);
      else
        ssize = LoadMultigridBoundaryBufferToCoarser(bdata_[bcolor_].send[nb.bufid],
                                                     nb, folddata, fcoeff);
    } else {
      if (type == BoundaryQuantity::mg_faceonly)
        ssize = LoadMultigridBoundaryBufferToFinerFluxCons(
                                                  bdata_[bcolor_].send[nb.bufid], nb);
      else
        ssize = LoadMultigridBoundaryBufferToFiner(bdata_[bcolor_].send[nb.bufid],
                                                   nb, folddata, fcoeff);
    }
    if (nb.snb.rank == Globals::my_rank) {
      std::memcpy(pmg->pmgbval->bdata_[bcolor_].recv[nb.targetid],
                  bdata_[bcolor_].send[nb.bufid], ssize*sizeof(Real));
      pmg->pmgbval->bdata_[bcolor_].flag[nb.targetid] = BoundaryStatus::arrived;
    }
#ifdef MPI_PARALLEL
    else { // NOLINT
      int tag = CreateBvalsMPITag(nb.snb.lid, nb.targetid,
                                  pmy_mg_->pmy_driver_->mg_phys_id_);
      MPI_Isend(bdata_[bcolor_].send[nb.bufid], ssize, MPI_ATHENA_REAL, nb.snb.rank, tag,
                mgcomm_, &(bdata_[bcolor_].req_send[nb.bufid]));
    }
#endif
    bdata_[bcolor_].sflag[nb.bufid] = BoundaryStatus::completed;
  }

  return bflag;
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::SetMultigridBoundarySameLevel(const Real *buf,
//!                               const NeighborBlock& nb, bool folddata, bool fcoeff)
//! \brief Set Multigrid boundary received from a block on the same level

void MGBoundaryValues::SetMultigridBoundarySameLevel(const Real *buf,
                                   const NeighborBlock& nb, bool folddata, bool fcoeff) {
  int nc = pmy_mg_->GetCurrentNumberOfCells();
  int ngh = pmy_mg_->ngh_, nvar = pmy_mg_->nvar_;
  int si, sj, sk, ei, ej, ek;
  AthenaArray<Real> *t;
  if (!fcoeff) {
    t = &(pmy_mg_->GetCurrentData());
  } else {
    nvar = pmy_mg_->ncoeff_;
    t = &(pmy_mg_->GetCurrentCoefficient());
  }
  AthenaArray<Real> &dst = *t;
  AthenaArray<Real> &old = pmy_mg_->GetCurrentOldData();

  if (nb.ni.ox1 == 0)     si = ngh,      ei = nc + ngh - 1;
  else if (nb.ni.ox1 > 0) si = nc + ngh, ei = nc + 2*ngh - 1;
  else                    si = 0,        ei = ngh - 1;
  if (nb.ni.ox2 == 0)     sj = ngh,      ej = nc + ngh - 1;
  else if (nb.ni.ox2 > 0) sj = nc + ngh, ej = nc + 2*ngh - 1;
  else                    sj = 0,        ej = ngh - 1;
  if (nb.ni.ox3 == 0)     sk = ngh,      ek = nc + ngh - 1;
  else if (nb.ni.ox3 > 0) sk = nc + ngh, ek = nc + 2*ngh - 1;
  else                    sk = 0,        ek = ngh - 1;

  int p = 0;
  BufferUtility::UnpackData(buf, dst, 0, nvar-1, si, ei, sj, ej, sk, ek, p);
  if (folddata)
    BufferUtility::UnpackData(buf, old, 0, nvar-1, si, ei, sj, ej, sk, ek, p);

  if (pmy_mg_->pmy_driver_->nreflevel_ > 0) {
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
    if (folddata)
      BufferUtility::UnpackData(buf, cbufold_, 0, nvar-1, si, ei, sj, ej, sk, ek, p);
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::SetMultigridBoundaryFromCoarser(const Real *buf,
//!                               const NeighborBlock& nb, boolf folddata, bool fcoeff)
//! \brief Set hydro boundary received from a block on the same level

void MGBoundaryValues::SetMultigridBoundaryFromCoarser(const Real *buf,
                                   const NeighborBlock& nb, bool folddata, bool fcoeff) {
  int nc = pmy_mg_->GetCurrentNumberOfCells();
  int ngh = pmy_mg_->ngh_, nvar = pmy_mg_->nvar_;
  if (fcoeff) nvar = pmy_mg_->ncoeff_;
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
  if (folddata)
    BufferUtility::UnpackData(buf, cbufold_, 0, nvar-1, si, ei, sj, ej, sk, ek, p);

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::SetMultigridBoundaryFromFiner(const Real *buf,
//!                               const NeighborBlock& nb, bool folddata, bool fcoeff)
//! \brief Set hydro boundary received from a block on the same level

void MGBoundaryValues::SetMultigridBoundaryFromFiner(const Real *buf,
                                   const NeighborBlock& nb, bool folddata, bool fcoeff) {
  int nc = pmy_mg_->GetCurrentNumberOfCells();
  int ngh = pmy_mg_->ngh_, nvar = pmy_mg_->nvar_;
  int fs = ngh, fe = fs + nc - 1;
  int si, sj, sk, ei, ej, ek;
  AthenaArray<Real> *t;
  if (!fcoeff) {
    t = &(pmy_mg_->GetCurrentData());
  } else {
    nvar = pmy_mg_->ncoeff_;
    t = &(pmy_mg_->GetCurrentCoefficient());
  }
  AthenaArray<Real> &dst = *t;
  AthenaArray<Real> &old = pmy_mg_->GetCurrentOldData();

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
  if (folddata)
    BufferUtility::UnpackData(buf, old, 0, nvar-1, si, ei, sj, ej, sk, ek, p);

  return;
}


//----------------------------------------------------------------------------------------
//! \fn bool MGBoundaryValues::ReceiveMultigridBoundaryBuffers(BoundaryQuantity type,
//!                                                            bool folddata)
//! \brief receive the boundary data

bool MGBoundaryValues::ReceiveMultigridBoundaryBuffers(BoundaryQuantity type,
                                                       bool folddata) {
  bool bflag = true;

  for (int n = 0; n < nneighbor; ++n) {
    NeighborBlock& nb = neighbor[n];
    if (bdata_[bcolor_].flag[nb.bufid] == BoundaryStatus::completed) continue;
    if (type == BoundaryQuantity::mg_faceonly && nb.snb.level > loc.level
      && nb.ni.type != NeighborConnect::face) continue;
    if (bdata_[bcolor_].flag[nb.bufid] == BoundaryStatus::waiting) {
      if (nb.snb.rank == Globals::my_rank) {// on the same process
        bflag = false;
        continue;
      }
#ifdef MPI_PARALLEL
      else { // NOLINT
        int test;
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, mgcomm_, &test, MPI_STATUS_IGNORE);
        MPI_Test(&(bdata_[bcolor_].req_recv[nb.bufid]), &test, MPI_STATUS_IGNORE);
        if (!static_cast<bool>(test)) {
          bflag = false;
          continue;
        }
        bdata_[bcolor_].flag[nb.bufid] = BoundaryStatus::arrived;
      }
#endif
    }
    if (nb.snb.level == loc.level) {
      SetMultigridBoundarySameLevel(bdata_[bcolor_].recv[nb.bufid], nb, folddata, false);
    } else if (nb.snb.level < loc.level) {
      if (type == BoundaryQuantity::mg_faceonly)
        SetMultigridBoundaryFromCoarserFluxCons(bdata_[bcolor_].recv[nb.bufid], nb);
      else
        SetMultigridBoundaryFromCoarser(bdata_[bcolor_].recv[nb.bufid], nb,
                                        folddata, false);
    } else {
      if (type == BoundaryQuantity::mg_faceonly)
        SetMultigridBoundaryFromFinerFluxCons(bdata_[bcolor_].recv[nb.bufid], nb);
      else
        SetMultigridBoundaryFromFiner(bdata_[bcolor_].recv[nb.bufid], nb,
                                      folddata, false);
    }
    bdata_[bcolor_].flag[nb.bufid] = BoundaryStatus::completed; // completed
  }

  return bflag;
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::ReceiveMultigridCoefficientBoundaryBuffers()
//! \brief receive the boundary data for Multigrid coefficients

void MGBoundaryValues::ReceiveMultigridCoefficientBoundaryBuffers() {
  for (int n = 0; n < nneighbor; ++n) {
    NeighborBlock& nb = neighbor[n];
#ifdef MPI_PARALLEL
    if (nb.snb.rank != Globals::my_rank)
      MPI_Wait(&(bdata_[bcolor_].req_recv[nb.bufid]),MPI_STATUS_IGNORE);
#endif
    if (nb.snb.level == loc.level)
      SetMultigridBoundarySameLevel(bdata_[bcolor_].recv[nb.bufid], nb, false, true);
    else if (nb.snb.level < loc.level)
      SetMultigridBoundaryFromCoarser(bdata_[bcolor_].recv[nb.bufid], nb, false, true);
    else
      SetMultigridBoundaryFromFiner(bdata_[bcolor_].recv[nb.bufid], nb, false, true);
    bdata_[bcolor_].flag[nb.bufid] = BoundaryStatus::completed; // completed
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::ProlongateMultigridBoundaries(bool folddata, bool fcoeff)
//! \brief prolongate boundaries for Multigrid

void MGBoundaryValues::ProlongateMultigridBoundaries(bool folddata, bool fcoeff) {
  int nc = pmy_mg_->GetCurrentNumberOfCells();
  int ngh = pmy_mg_->ngh_, nvar = pmy_mg_->nvar_;
  Real time = pmy_mesh_->time;
  int ll = pmy_mg_->nlevel_ - 1 - pmy_mg_->GetCurrentLevel();
  const int cn = 1, cs = cn, ce = cs + nc/2 -1;
  const int flim = ngh + nc;
  AthenaArray<Real> *t;
  if (!fcoeff) {
    t = &(pmy_mg_->GetCurrentData());
  } else {
    nvar = pmy_mg_->ncoeff_;
    t = &(pmy_mg_->GetCurrentCoefficient());
  }
  AthenaArray<Real> &dst = *t;
  AthenaArray<Real> &old = pmy_mg_->GetCurrentOldData();

  ApplyPhysicalBoundaries(1, fcoeff);
  if (folddata)
    ApplyPhysicalBoundaries(2, fcoeff);

  for (int n = 0; n < nneighbor; ++n) {
    NeighborBlock& nb = neighbor[n];
    if (nb.snb.level >= loc.level) continue;

    // calculate the loop limits for the ghost zones
    int si, ei, sj, ej, sk, ek;
    if (nb.ni.ox1 == 0) {
      si = cs, ei = ce;
      if ((loc.lx1 & 1LL) == 0LL) ei += cn;
      else                        si -= cn;
    } else if (nb.ni.ox1 > 0) {
      si = ce + 1,  ei = ce + cn;
    } else {
      si = cs-cn,   ei = cs-1;
    }
    if (nb.ni.ox2 == 0) {
      sj = cs, ej = ce;
      if ((loc.lx2 & 1LL) == 0LL) ej += cn;
      else                        sj -= cn;
    } else if (nb.ni.ox2 > 0) {
      sj = ce + 1,  ej = ce + cn;
    } else {
      sj = cs-cn,   ej = cs-1;
    }
    if (nb.ni.ox3 == 0) {
      sk = cs, ek = ce;
      if ((loc.lx3 & 1LL) == 0LL) ek += cn;
      else                        sk -= cn;
    } else if (nb.ni.ox3 > 0) {
      sk = ce + 1,  ek = ce + cn;
    } else {
      sk = cs-cn,   ek = cs-1;
    }

    // Prolongation using tri-linear interpolation
    for (int v=0; v<nvar; ++v) {
      for (int k=sk; k<=ek; ++k) {
        int fk = (k - cn) * 2 + ngh;
        for (int j=sj; j<=ej; ++j) {
          int fj = (j - cn) * 2 + ngh;
#pragma ivdep
          for (int i=si; i<=ei; ++i) {
            int fi = (i - cn) * 2 + ngh;
            if (fk >= 0 && fj >= 0 && fi >= 0)
              dst(v, fk,   fj,   fi  ) =
                0.015625*(27.0*cbuf_(v,k,j,i)+cbuf_(v,k-1,j-1,i-1)
                        +9.0*(cbuf_(v,k,j,i-1)+cbuf_(v,k,j-1,i)+cbuf_(v,k-1,j,i))
                        +3.0*(cbuf_(v,k-1,j-1,i)+cbuf_(v,k-1,j,i-1)+cbuf_(v,k,j-1,i-1)));
            if (fk >= 0 && fj >= 0 && fi < flim)
              dst(v, fk,   fj,   fi+1) =
                0.015625*(27.0*cbuf_(v,k,j,i)+cbuf_(v,k-1,j-1,i+1)
                        +9.0*(cbuf_(v,k,j,i+1)+cbuf_(v,k,j-1,i)+cbuf_(v,k-1,j,i))
                        +3.0*(cbuf_(v,k-1,j-1,i)+cbuf_(v,k-1,j,i+1)+cbuf_(v,k,j-1,i+1)));
            if (fk >= 0 && fj < flim && fi >= 0)
              dst(v, fk,   fj+1, fi  ) =
                0.015625*(27.0*cbuf_(v,k,j,i)+cbuf_(v,k-1,j+1,i-1)
                        +9.0*(cbuf_(v,k,j,i-1)+cbuf_(v,k,j+1,i)+cbuf_(v,k-1,j,i))
                        +3.0*(cbuf_(v,k-1,j+1,i)+cbuf_(v,k-1,j,i-1)+cbuf_(v,k,j+1,i-1)));
            if (fk < flim && fj >= 0 && fi >= 0)
              dst(v, fk+1, fj,   fi  ) =
                0.015625*(27.0*cbuf_(v,k,j,i)+cbuf_(v,k+1,j-1,i-1)
                        +9.0*(cbuf_(v,k,j,i-1)+cbuf_(v,k,j-1,i)+cbuf_(v,k+1,j,i))
                        +3.0*(cbuf_(v,k+1,j-1,i)+cbuf_(v,k+1,j,i-1)+cbuf_(v,k,j-1,i-1)));
            if (fk < flim && fj < flim && fi >= 0)
              dst(v, fk+1, fj+1, fi  ) =
                0.015625*(27.0*cbuf_(v,k,j,i)+cbuf_(v,k+1,j+1,i-1)
                        +9.0*(cbuf_(v,k,j,i-1)+cbuf_(v,k,j+1,i)+cbuf_(v,k+1,j,i))
                        +3.0*(cbuf_(v,k+1,j+1,i)+cbuf_(v,k+1,j,i-1)+cbuf_(v,k,j+1,i-1)));
            if (fk < flim && fj >= 0 && fi < flim)
              dst(v, fk+1, fj,   fi+1) =
                0.015625*(27.0*cbuf_(v,k,j,i)+cbuf_(v,k+1,j-1,i+1)
                        +9.0*(cbuf_(v,k,j,i+1)+cbuf_(v,k,j-1,i)+cbuf_(v,k+1,j,i))
                        +3.0*(cbuf_(v,k+1,j-1,i)+cbuf_(v,k+1,j,i+1)+cbuf_(v,k,j-1,i+1)));
            if (fk >= 0 && fj < flim && fi < flim)
              dst(v, fk,  fj+1, fi+1) =
                0.015625*(27.0*cbuf_(v,k,j,i)+cbuf_(v,k-1,j+1,i+1)
                        +9.0*(cbuf_(v,k,j,i+1)+cbuf_(v,k,j+1,i)+cbuf_(v,k-1,j,i))
                        +3.0*(cbuf_(v,k-1,j+1,i)+cbuf_(v,k-1,j,i+1)+cbuf_(v,k,j+1,i+1)));
            if (fk < flim && fj < flim && fi < flim)
              dst(v, fk+1, fj+1, fi+1) =
                0.015625*(27.0*cbuf_(v,k,j,i)+cbuf_(v,k+1,j+1,i+1)
                        +9.0*(cbuf_(v,k,j,i+1)+cbuf_(v,k,j+1,i)+cbuf_(v,k+1,j,i))
                        +3.0*(cbuf_(v,k+1,j+1,i)+cbuf_(v,k+1,j,i+1)+cbuf_(v,k,j+1,i+1)));
          }
        }
      }
    }
    if (folddata) {
      for (int v=0; v<nvar; ++v) {
        for (int k=sk; k<=ek; ++k) {
          int fk = (k - cn) * 2 + ngh;
          for (int j=sj; j<=ej; ++j) {
            int fj = (j - cn) * 2 + ngh;
#pragma ivdep
            for (int i=si; i<=ei; ++i) {
              int fi = (i - cn) * 2 + ngh;
              if (fk >= 0 && fj >= 0 && fi >= 0)
                old(v, fk,   fj,   fi  ) =
                0.015625*(27.0*cbufold_(v,k,j,i)+cbufold_(v,k-1,j-1,i-1)
                +9.0*(cbufold_(v,k,j,i-1)+cbufold_(v,k,j-1,i)+cbufold_(v,k-1,j,i))
                +3.0*(cbufold_(v,k-1,j-1,i)+cbufold_(v,k-1,j,i-1)+cbufold_(v,k,j-1,i-1)));
              if (fk >= 0 && fj >= 0 && fi < flim)
                old(v, fk,   fj,   fi+1) =
                0.015625*(27.0*cbufold_(v,k,j,i)+cbufold_(v,k-1,j-1,i+1)
                +9.0*(cbufold_(v,k,j,i+1)+cbufold_(v,k,j-1,i)+cbufold_(v,k-1,j,i))
                +3.0*(cbufold_(v,k-1,j-1,i)+cbufold_(v,k-1,j,i+1)+cbufold_(v,k,j-1,i+1)));
              if (fk >= 0 && fj < flim && fi >= 0)
                old(v, fk,   fj+1, fi  ) =
                0.015625*(27.0*cbufold_(v,k,j,i)+cbufold_(v,k-1,j+1,i-1)
                +9.0*(cbufold_(v,k,j,i-1)+cbufold_(v,k,j+1,i)+cbufold_(v,k-1,j,i))
                +3.0*(cbufold_(v,k-1,j+1,i)+cbufold_(v,k-1,j,i-1)+cbufold_(v,k,j+1,i-1)));
              if (fk < flim && fj >= 0 && fi >= 0)
                old(v, fk+1, fj,   fi  ) =
                0.015625*(27.0*cbufold_(v,k,j,i)+cbufold_(v,k+1,j-1,i-1)
                +9.0*(cbufold_(v,k,j,i-1)+cbufold_(v,k,j-1,i)+cbufold_(v,k+1,j,i))
                +3.0*(cbufold_(v,k+1,j-1,i)+cbufold_(v,k+1,j,i-1)+cbufold_(v,k,j-1,i-1)));
              if (fk < flim && fj < flim && fi >= 0)
                old(v, fk+1, fj+1, fi  ) =
                0.015625*(27.0*cbufold_(v,k,j,i)+cbufold_(v,k+1,j+1,i-1)
                +9.0*(cbufold_(v,k,j,i-1)+cbufold_(v,k,j+1,i)+cbufold_(v,k+1,j,i))
                +3.0*(cbufold_(v,k+1,j+1,i)+cbufold_(v,k+1,j,i-1)+cbufold_(v,k,j+1,i-1)));
              if (fk < flim && fj >= 0 && fi < flim)
                old(v, fk+1, fj,   fi+1) =
                0.015625*(27.0*cbufold_(v,k,j,i)+cbufold_(v,k+1,j-1,i+1)
                +9.0*(cbufold_(v,k,j,i+1)+cbufold_(v,k,j-1,i)+cbufold_(v,k+1,j,i))
                +3.0*(cbufold_(v,k+1,j-1,i)+cbufold_(v,k+1,j,i+1)+cbufold_(v,k,j-1,i+1)));
              if (fk >= 0 && fj < flim && fi < flim)
                old(v, fk,  fj+1, fi+1) =
                0.015625*(27.0*cbufold_(v,k,j,i)+cbufold_(v,k-1,j+1,i+1)
                +9.0*(cbufold_(v,k,j,i+1)+cbufold_(v,k,j+1,i)+cbufold_(v,k-1,j,i))
                +3.0*(cbufold_(v,k-1,j+1,i)+cbufold_(v,k-1,j,i+1)+cbufold_(v,k,j+1,i+1)));
              if (fk < flim && fj < flim && fi < flim)
                old(v, fk+1, fj+1, fi+1) =
                0.015625*(27.0*cbufold_(v,k,j,i)+cbufold_(v,k+1,j+1,i+1)
                +9.0*(cbufold_(v,k,j,i+1)+cbufold_(v,k,j+1,i)+cbufold_(v,k+1,j,i))
                +3.0*(cbufold_(v,k+1,j+1,i)+cbufold_(v,k+1,j,i+1)+cbufold_(v,k,j+1,i+1)));
            }
          }
        }
      }
    }
  } // end loop over nneighbor

  return;
}


//----------------------------------------------------------------------------------------
//! \fn int MGBoundaryValues::LoadMultigridBoundaryBufferToCoarserFluxCons(
//!                                               Real *buf, const NeighborBlock& nb)
//! \brief Set Multigrid boundary buffers for sending to a block on the coarser level
//!        with flux-conservation for the base class (just call the normal version)

int MGBoundaryValues::LoadMultigridBoundaryBufferToCoarserFluxCons(Real *buf,
                                                          const NeighborBlock& nb) {
  return LoadMultigridBoundaryBufferToCoarser(buf, nb, false, false);
}


//----------------------------------------------------------------------------------------
//! \fn int MGBoundaryValues::LoadMultigridBoundaryBufferToFinerFluxCons(Real *buf,
//!                                                             const NeighborBlock& nb)
//! \brief Set Multigrid boundary buffers for sending to a block on the finer level
//!        with flux-conservation for the base class (just call the normal version)

int MGBoundaryValues::LoadMultigridBoundaryBufferToFinerFluxCons(Real *buf,
                                                               const NeighborBlock& nb) {
  return LoadMultigridBoundaryBufferToFiner(buf, nb, false, false);
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::SetMultigridBoundaryFromCoarserFluxCons(
//!                                      const Real *buf, const NeighborBlock& nb)
//! \brief Set Multigrid boundary received from a block from the coarser level
//!        with flux-conservation for the base class (just call the normal version)

void MGBoundaryValues::SetMultigridBoundaryFromCoarserFluxCons(const Real *buf,
                                                             const NeighborBlock& nb) {
  SetMultigridBoundaryFromCoarser(buf, nb, false, false);
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::SetMultigridBoundaryFromFinerFluxCons(
//!                                      const Real *buf, const NeighborBlock& nb)
//! \brief Set Multigrid boundary received from a block from the finer level
//!        with flux-conservation for the base class (just call the normal version)

void MGBoundaryValues::SetMultigridBoundaryFromFinerFluxCons(const Real *buf,
                                                             const NeighborBlock& nb) {
  SetMultigridBoundaryFromFiner(buf, nb, false, false);
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::ProlongateMultigridBoundariesFluxCons()
//! \brief prolongate Multigrid boundaries with the flux-conservation formula
//!        for the base class (just call the normal version)

void MGBoundaryValues::ProlongateMultigridBoundariesFluxCons() {
  ProlongateMultigridBoundaries(false, false);
  return;
}
