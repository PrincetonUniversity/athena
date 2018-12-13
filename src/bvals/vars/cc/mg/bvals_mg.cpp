//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_mg.cpp
//  \brief

// C++ headers
#include <cmath>
#include <cstdlib>
#include <cstring>    // memcpy
#include <iomanip>
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ classes headers
#include "bvals_mg.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../globals.hpp"
#include "../mesh/mesh.hpp"
#include "../multigrid/multigrid.hpp"
#include "../parameter_input.hpp"
#include "../utils/buffer_utils.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

class Multigrid;
class MultigridDriver;

//----------------------------------------------------------------------------------------
//! \fn MGBoundaryValues::MGBoundaryValues(Multigrid *pmg, enum BoundaryFlag *input_bcs,
//                                         MGBoundaryFunc_t *MGBoundary)
//  \brief Constructor of the MGBoundaryValues class

MGBoundaryValues::MGBoundaryValues(Multigrid *pmg, enum BoundaryFlag *input_bcs,
                                   MGBoundaryFunc_t *MGBoundary)
  : BoundaryBase(pmg->pmy_driver_->pmy_mesh_, pmg->loc_, pmg->size_, input_bcs) {
  pmy_mg_=pmg;
#ifdef MPI_PARALLEL
  mgcomm_=pmg->pmy_driver_->MPI_COMM_MULTIGRID;
#endif
  if (pmy_mg_->root_flag_==true) {
    for (int i=0; i<6; i++)
      MGBoundaryFunction_[i]=MGBoundary[i];
  } else {
    for (int i=0; i<6; i++) {
      if (block_bcs[i]==PERIODIC_BNDRY || block_bcs[i]==BLOCK_BNDRY)
        MGBoundaryFunction_[i]=NULL;
      else
        MGBoundaryFunction_[i]=MGBoundary[i];
    }
    InitBoundaryData(bd_mggrav_, BNDRY_MGGRAV);
  }
}


//----------------------------------------------------------------------------------------
//! \fn MGBoundaryValues::~MGBoundaryValues()
//  \brief Destructor of the MGBoundaryValues class

MGBoundaryValues::~MGBoundaryValues() {
  if (pmy_mg_->root_flag_ == false)
    DestroyBoundaryData(bd_mggrav_);
}

//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::InitBoundaryData(MGBoundaryData &bd,
//                                              enum BoundaryType type)
//  \brief Initialize MGBoundaryData structure

void MGBoundaryValues::InitBoundaryData(MGBoundaryData &bd, enum BoundaryType type) {
  int size;
  bd.nbmax=maxneighbor_;
  for (int n=0;n<bd.nbmax;n++) {
    // Clear flags and requests
    bd.flag[n]=BNDRY_WAITING;
    bd.sflag[n]=BNDRY_WAITING;
    bd.send[n]=NULL;
    bd.recv[n]=NULL;
#ifdef MPI_PARALLEL
    bd.req_send[n]=MPI_REQUEST_NULL;
    bd.req_recv[n]=MPI_REQUEST_NULL;
#endif

    // Allocate buffers
    // calculate the buffer size
    switch(type) {
      case BNDRY_MGGRAV: {
        int ngh=pmy_mg_->ngh_;
        if (pmy_mesh_->multilevel) { // with refinement - NGHOST = 1
          int nc=block_size_.nx1;
          if (BoundaryValues::ni[n].type==NEIGHBOR_FACE) size=SQR(nc)*ngh;
          else if (BoundaryValues::ni[n].type==NEIGHBOR_EDGE) size=nc*ngh*ngh +
                                                                  (nc*ngh*ngh)/2;
          else if (BoundaryValues::ni[n].type==NEIGHBOR_CORNER) size=ngh*ngh*ngh*2;
        } else { // uniform - NGHOST=1
          size=((BoundaryValues::ni[n].ox1==0)?block_size_.nx1:ngh)
              *((BoundaryValues::ni[n].ox2==0)?block_size_.nx2:ngh)
              *((BoundaryValues::ni[n].ox3==0)?block_size_.nx3:ngh);
        }
      }
      break;
      default: {
        std::stringstream msg;
        msg << "### FATAL ERROR in InitBoundaryData" << std::endl
            << "Invalid boundary type is specified." << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      break;
    }
    bd.send[n]=new Real[size];
    bd.recv[n]=new Real[size];
  }
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::DestroyBoundaryData(MGBoundaryData &bd)
//  \brief Destroy MGBoundaryData structure
void MGBoundaryValues::DestroyBoundaryData(MGBoundaryData &bd) {
  for (int n=0;n<bd.nbmax;n++) {
    delete [] bd.send[n];
    delete [] bd.recv[n];
#ifdef MPI_PARALLEL
    if (bd.req_send[n]!=MPI_REQUEST_NULL)
      MPI_Request_free(&bd.req_send[n]);
    if (bd.req_recv[n]!=MPI_REQUEST_NULL)
      MPI_Request_free(&bd.req_recv[n]);
#endif
  }
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::ApplyPhysicalBoundaries(void)
//  \brief Apply physical boundary conditions to the current Multigrid data

void MGBoundaryValues::ApplyPhysicalBoundaries(void) {
  AthenaArray<Real> &dst=pmy_mg_->GetCurrentData();
  int ll=pmy_mg_->nlevel_-1-pmy_mg_->current_level_;
  int ngh=pmy_mg_->ngh_, nvar=pmy_mg_->nvar_;
  int ncx=block_size_.nx1>>ll, ncy=block_size_.nx2>>ll, ncz=block_size_.nx3>>ll;
  int is=ngh, ie=ncx+ngh-1, js=ngh, je=ncy+ngh-1, ks=ngh, ke=ncz+ngh-1;
  int bis=is-ngh, bie=ie+ngh, bjs=js, bje=je, bks=ks, bke=ke;
  Real dx=pmy_mg_->rdx_*static_cast<Real>(1<<ll);
  Real dy=pmy_mg_->rdy_*static_cast<Real>(1<<ll);
  Real dz=pmy_mg_->rdz_*static_cast<Real>(1<<ll);
  Real x0=block_size_.x1min-(static_cast<Real>(ngh)+0.5)*dx;
  Real y0=block_size_.x2min-(static_cast<Real>(ngh)+0.5)*dy;
  Real z0=block_size_.x3min-(static_cast<Real>(ngh)+0.5)*dz;
  Real time=pmy_mesh_->time;
  if (MGBoundaryFunction_[INNER_X2]==NULL) bjs=js-ngh;
  if (MGBoundaryFunction_[OUTER_X2]==NULL) bje=je+ngh;
  if (MGBoundaryFunction_[INNER_X3]==NULL) bks=ks-ngh;
  if (MGBoundaryFunction_[OUTER_X3]==NULL) bke=ke+ngh;

  // Apply boundary function on inner-x1
  if (MGBoundaryFunction_[INNER_X1] != NULL)
    MGBoundaryFunction_[INNER_X1](dst, time, nvar, is, ie, bjs, bje, bks, bke, ngh,
                                  x0, y0, z0, dx, dy, dz);
  // Apply boundary function on outer-x1
  if (MGBoundaryFunction_[OUTER_X1] != NULL)
    MGBoundaryFunction_[OUTER_X1](dst, time, nvar, is, ie, bjs, bje, bks, bke, ngh,
                                  x0, y0, z0, dx, dy, dz);

  // Apply boundary function on inner-x2
  if (MGBoundaryFunction_[INNER_X2] != NULL)
    MGBoundaryFunction_[INNER_X2](dst, time, nvar, bis, bie, js, je, bks, bke, ngh,
                                  x0, y0, z0, dx, dy, dz);
  // Apply boundary function on outer-x2
  if (MGBoundaryFunction_[OUTER_X2] != NULL)
    MGBoundaryFunction_[OUTER_X2](dst, time, nvar, bis, bie, js, je, bks, bke, ngh,
                                  x0, y0, z0, dx, dy, dz);

  bjs=js-ngh, bje=je+ngh;
  // Apply boundary function on inner-x3
  if (MGBoundaryFunction_[INNER_X3] != NULL)
    MGBoundaryFunction_[INNER_X3](dst, time, nvar, bis, bie, bjs, bje, ks, ke, ngh,
                                  x0, y0, z0, dx, dy, dz);
  // Apply boundary function on outer-x3
  if (MGBoundaryFunction_[OUTER_X3] != NULL)
    MGBoundaryFunction_[OUTER_X3](dst, time, nvar, bis, bie, bjs, bje, ks, ke, ngh,
                                  x0, y0, z0, dx, dy, dz);

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::StartReceivingMultigrid(int nc, enum BoundaryType type)
//  \brief initiate MPI_Irecv for multigrid

void MGBoundaryValues::StartReceivingMultigrid(int nc, enum BoundaryType type) {
  int mylevel=loc.level;
  int nvar, tag, phys, ngh;
  bool faceonly=false;
  MGBoundaryData *pbd;

  if (type==BNDRY_MGGRAV || type==BNDRY_MGGRAVF) {
    pbd=&bd_mggrav_;
    nvar=1, ngh=1;
    phys=TAG_MGGRAV;
  }
  if (type==BNDRY_MGGRAVF)
    faceonly=true;
  for (int n=0;n<nneighbor;n++) {
    NeighborBlock& nb = neighbor[n];
    if (faceonly && nb.type>NEIGHBOR_FACE) break;
#ifdef MPI_PARALLEL
    if (nb.rank!=Globals::my_rank) {
      int size;
      if (pmy_mesh_->multilevel==true) { // with refinement - NGHOST = 1
        if (nb.level == mylevel) { // same
          if (nb.type==NEIGHBOR_FACE) size=SQR(nc);
          else if (nb.type==NEIGHBOR_EDGE) size=nc+nc/2;
          else if (nb.type==NEIGHBOR_CORNER) size=2;
        } else if (nb.level < mylevel) { // coarser
          if (nb.type==NEIGHBOR_FACE) size=SQR(nc/2+1);
          else if (nb.type==NEIGHBOR_EDGE) size=SQR(nc/2+1);
          else if (nb.type==NEIGHBOR_CORNER) size=1;
        } else { // finer
          if (nb.type==NEIGHBOR_FACE) size=SQR(nc);
          else if (nb.type==NEIGHBOR_EDGE) size=nc/2;
          else if (nb.type==NEIGHBOR_CORNER) size=1;
        }
      } else { // no SMR/AMR - NGHOST=2 (Not necessarily! KGF)
        if (nb.type==NEIGHBOR_FACE) size=nc*nc*ngh;
        else if (nb.type==NEIGHBOR_EDGE) size=nc*ngh*ngh;
        else if (nb.type==NEIGHBOR_CORNER) size=ngh*ngh*ngh;
      }
      size*=nvar;
      tag=CreateBvalsMPITag(pmy_mg_->lid_, phys, nb.bufid);
      MPI_Irecv(pbd->recv[nb.bufid], size, MPI_ATHENA_REAL, nb.rank, tag,
                mgcomm_, &(pbd->req_recv[nb.bufid]));
    }
#endif
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::ClearBoundaryMultigrid(enum BoundaryType type)
//  \brief clean up the boundary flags after each loop for multigrid

void MGBoundaryValues::ClearBoundaryMultigrid(enum BoundaryType type) {
  bool faceonly=false;
  MGBoundaryData *pbd;

  if (type==BNDRY_MGGRAV || type==BNDRY_MGGRAVF)
    pbd=&bd_mggrav_;
  if (type==BNDRY_MGGRAVF)
    faceonly=true;

  for (int n=0;n<nneighbor;n++) {
    NeighborBlock& nb = neighbor[n];
    if (faceonly && nb.type>NEIGHBOR_FACE) break;
    pbd->flag[nb.bufid] = BNDRY_WAITING;
    pbd->sflag[nb.bufid] = BNDRY_WAITING;
#ifdef MPI_PARALLEL
    if (nb.rank!=Globals::my_rank)
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

int MGBoundaryValues::LoadMultigridBoundaryBufferSameLevel(AthenaArray<Real> &src,
                      int nvar, int nc, int ngh, Real *buf, const NeighborBlock& nb) {
  int si, sj, sk, ei, ej, ek;

  si=(nb.ox1>0)?nc:ngh;
  ei=(nb.ox1<0)?(2*ngh-1):(ngh+nc-1);
  sj=(nb.ox2>0)?nc:ngh;
  ej=(nb.ox2<0)?(2*ngh-1):(ngh+nc-1);
  sk=(nb.ox3>0)?nc:ngh;
  ek=(nb.ox3<0)?(2*ngh-1):(ngh+nc-1);
  int p=0;
  BufferUtility::Pack4DData(src, buf, 0, nvar-1, si, ei, sj, ej, sk, ek, p);
  return p;
}


//----------------------------------------------------------------------------------------
//! \fn bool MGBoundaryValues::SendMultigridBoundaryBuffers(AthenaArray<Real> &src,
//                                                  int nc, enum BoundaryType type)
//  \brief Send boundary buffers

bool MGBoundaryValues::SendMultigridBoundaryBuffers(AthenaArray<Real> &src,
                                                    int nc, enum BoundaryType type) {
  int mylevel=loc.level;
  int nvar, tag, ngh, phys;
  bool faceonly=false;
  bool bflag=true;
  MGBoundaryData *pbd, *ptarget;

  if (type==BNDRY_MGGRAV || type==BNDRY_MGGRAVF) {
    pbd=&bd_mggrav_;
    nvar=1, ngh=1;
    phys=TAG_MGGRAV;
  }
  if (type==BNDRY_MGGRAVF)
    faceonly=true;
  for (int n=0; n<nneighbor; n++) {
    NeighborBlock& nb = neighbor[n];
    if (faceonly && nb.type>NEIGHBOR_FACE) break;
    if (pbd->sflag[nb.bufid]==BNDRY_COMPLETED) continue;
    int ssize;
    if (nb.rank == Globals::my_rank) {
      Multigrid *pmg=pmy_mg_->pmy_driver_->FindMultigrid(nb.gid);
      if (type==BNDRY_MGGRAV || type==BNDRY_MGGRAVF)
        ptarget=&(pmg->pmgbval->bd_mggrav_);
      if (ptarget->flag[nb.targetid] != BNDRY_WAITING) {
        bflag=false;
        continue;
      }
    }
    if (nb.level==mylevel)
      ssize=LoadMultigridBoundaryBufferSameLevel(src, nvar, nc, ngh,
                                                 pbd->send[nb.bufid], nb);
//    else if (nb.level<mylevel)
//      ssize=LoadMultigridBoundaryBufferToCoarser(src, nvar,
//                                                 pbd->send[nb.bufid], cbuf, nb);
//    else
//      ssize=LoadMultigridBoundaryBufferToFiner(src, nvar, pbd->send[nb.bufid], nb);
    if (nb.rank == Globals::my_rank) {
      std::memcpy(ptarget->recv[nb.targetid], pbd->send[nb.bufid], ssize*sizeof(Real));
      ptarget->flag[nb.targetid] = BNDRY_ARRIVED;
#ifdef MPI_PARALLEL
    } else { // MPI
      tag=CreateBvalsMPITag(nb.lid, phys, nb.targetid);
      MPI_Isend(pbd->send[nb.bufid], ssize, MPI_ATHENA_REAL, nb.rank, tag,
                mgcomm_, &(pbd->req_send[nb.bufid]));
    }
#else
    }
#endif
    pbd->sflag[nb.bufid] = BNDRY_COMPLETED;
  }

  return bflag;
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::SetMultigridBoundarySameLevel(AthenaArray<Real> &dst,
//                     int nvar, int nc, int ngh, Real *buf, const NeighborBlock& nb)
//  \brief Set hydro boundary received from a block on the same level

void MGBoundaryValues::SetMultigridBoundarySameLevel(AthenaArray<Real> &dst,
                        int nvar, int nc, int ngh, Real *buf, const NeighborBlock& nb) {
  int si, sj, sk, ei, ej, ek;

  if (nb.ox1==0)     si=ngh,    ei=nc+ngh-1;
  else if (nb.ox1>0) si=nc+ngh, ei=nc+2*ngh-1;
  else              si=0,      ei=ngh-1;
  if (nb.ox2==0)     sj=ngh,    ej=nc+ngh-1;
  else if (nb.ox2>0) sj=nc+ngh, ej=nc+2*ngh-1;
  else              sj=0,      ej=ngh-1;
  if (nb.ox3==0)     sk=ngh,    ek=nc+ngh-1;
  else if (nb.ox3>0) sk=nc+ngh, ek=nc+2*ngh-1;
  else              sk=0,      ek=ngh-1;

  int p=0;
  BufferUtility::Unpack4DData(buf, dst, 0, nvar-1, si, ei, sj, ej, sk, ek, p);
  return;
}


//----------------------------------------------------------------------------------------
//! \fn bool MGBoundaryValues::ReceiveMultigridBoundaryBuffers(AthenaArray<Real> &dst,
//                                                     int nc, enum BoundaryType type)
//  \brief receive the boundary data

bool MGBoundaryValues::ReceiveMultigridBoundaryBuffers(AthenaArray<Real> &dst,
                                                       int nc, enum BoundaryType type) {
  bool bflag=true, faceonly=false;
  int nvar, ngh;
  MGBoundaryData *pbd;

  if (type==BNDRY_MGGRAV || type==BNDRY_MGGRAVF) {
    pbd=&bd_mggrav_;
    nvar=1, ngh=1;
  }
  if (type==BNDRY_MGGRAVF)
    faceonly=true;

  for (int n=0; n<nneighbor; n++) {
    NeighborBlock& nb = neighbor[n];
    if (faceonly && nb.type>NEIGHBOR_FACE) break;
    if (pbd->flag[nb.bufid]==BNDRY_COMPLETED) continue;
    if (pbd->flag[nb.bufid]==BNDRY_WAITING) {
      if (nb.rank==Globals::my_rank) {// on the same process
        bflag=false;
        continue;
#ifdef MPI_PARALLEL
      } else { // MPI boundary
        int test;
        MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,mgcomm_,&test,MPI_STATUS_IGNORE);
        MPI_Test(&(pbd->req_recv[nb.bufid]),&test,MPI_STATUS_IGNORE);
        if (static_cast<bool>(test) == false) {
          bflag=false;
          continue;
        }
        pbd->flag[nb.bufid] = BNDRY_ARRIVED;
      }
#else
      }
#endif
    }
    if (nb.level==loc.level)
      SetMultigridBoundarySameLevel(dst, nvar, nc, ngh, pbd->recv[nb.bufid], nb);
//    else if (nb.level<loc.level) // this set only the prolongation buffer
//      SetMultigridBoundaryFromCoarser(nvar, nc, ngh, pbd->recv[nb.bufid], cbuf, nb);
//    else
//      SetMultigridBoundaryFromFiner(dst, nvar, nc, ngh, pbd->recv[nb.bufid], nb);
    pbd->flag[nb.bufid] = BNDRY_COMPLETED; // completed
  }
  return bflag;
}
