//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_phi.cpp
//  \brief functions that apply BCs for gravitational potential 

// C++ headers
#include <iostream>   // endl
#include <iomanip>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cstring>    // memcpy
#include <cstdlib>
#include <cmath>

// Athena++ classes headers
#include "bvals.hpp"
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../coordinates/coordinates.hpp"
#include "../parameter_input.hpp"
#include "../utils/buffer_utils.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::StartReceivingMultigrid(int nc, enum MGBoundaryType type)
//  \brief initiate MPI_Irecv for multigrid

void BoundaryValues::StartReceivingMultigrid(int nc, enum MGBoundaryType type)
{
  MeshBlock *pmb=pmy_block_;
  int mylevel=pmb->loc.level;
  int nvar, tag;
  Real *rbuf;
#ifdef MPI_PARALLEL
  MPI_Request *req;
#endif

  for(int n=0;n<pmb->nneighbor;n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    switch(type) {
      case BND_MGGRAV:
        nvar=1;
        tag=CreateBvalsMPITag(pmb->lid, TAG_MGGRAV, nb.bufid);
        rbuf=mggrav_recv_[nb.bufid];
#ifdef MPI_PARALLEL
        if(nb.rank!=Globals::my_rank)
          req=&(req_mggrav_recv_[nb.bufid]);
#endif
        break;
      default:
        break;
    }
#ifdef MPI_PARALLEL
    if(nb.rank!=Globals::my_rank) {
      int size;
      if(pmb->pmy_mesh->multilevel==true) { // with refinement - NGHOST = 1
        if(nb.level == mylevel) { // same
          if(nb.type==NEIGHBOR_FACE) size=SQR(nc);
          else if(nb.type==NEIGHBOR_EDGE) size=nc+nc/2;
          else if(nb.type==NEIGHBOR_CORNER) size=2;
        }
        else if(nb.level < mylevel) { // coarser
          if(nb.type==NEIGHBOR_FACE) size=SQR(nc/2+1);
          else if(nb.type==NEIGHBOR_EDGE) size=SQR(nc/2+1);
          else if(nb.type==NEIGHBOR_CORNER) size=1;
        }
        else { // finer
          if(nb.type==NEIGHBOR_FACE) size=SQR(nc);
          else if(nb.type==NEIGHBOR_EDGE) size=nc/2;
          else if(nb.type==NEIGHBOR_CORNER) size=1;
        }
      }
      else { // uniform - NGHOST=2
        if(nb.type==NEIGHBOR_FACE) size=nc*nc*2;
        else if(nb.type==NEIGHBOR_EDGE) size=nc*4;
        else if(nb.type==NEIGHBOR_CORNER) size=8;
      }
      size*=nvar;
      MPI_Irecv(rbuf,size,MPI_ATHENA_REAL,nb.rank,tag,MPI_COMM_WORLD,req);
    }
#endif
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ClearBoundaryMultigrid(enum MGBoundaryType type)
//  \brief clean up the boundary flags after each loop for multigrid

void BoundaryValues::ClearBoundaryMultigrid(enum MGBoundaryType type)
{
  MeshBlock *pmb=pmy_block_;
  enum BoundaryStatus *flag;
#ifdef MPI_PARALLEL
  MPI_Request *req;
#endif

  for(int n=0;n<pmb->nneighbor;n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    switch(type) {
      case BND_MGGRAV:
        flag=mggrav_flag_[nb.bufid];
#ifdef MPI_PARALLEL
        if(nb.rank!=Globals::my_rank)
          req=&(req_mggrav_send_[nb.bufid]);
#endif
        break;
      default:
        break;
    }
    *flag = BNDRY_WAITING;
#ifdef MPI_PARALLEL
    if(nb.rank!=Globals::my_rank)
      MPI_Wait(req,MPI_STATUS_IGNORE); // Wait for Isend
#endif
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadMultigridBoundaryBufferSameLevel(AthenaArray<Real> &src,
//                          int nvar, int nc, int ngh, Real *buf, const NeighborBlock& nb)
//  \brief Set hydro boundary buffers for sending to a block on the same level

int BoundaryValues::LoadMultigridBoundaryBufferSameLevel(AthenaArray<Real> &src,
                    int nvar, int nc, int ngh, Real *buf, const NeighborBlock& nb)
{
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
//! \fn void BoundaryValues::SendMultigridBoundaryBuffers(AthenaArray<Real> &src,
//                                                int nc, enum MGBoundaryType type)
//  \brief Send boundary buffers

void BoundaryValues::SendMultigridBoundaryBuffers(AthenaArray<Real> &src,
                                                  int nc, enum MGBoundaryType type)
{
  MeshBlock *pmb=pmy_block_, *pbl;
  int mylevel=pmb->loc.level;
  int nvar, tag, ngh;
  Real *sbuf, *rbuf;
//  AthenaArray<Real> cbuf;
  enum BoundaryStatus *flag;
#ifdef MPI_PARALLEL
  MPI_Request *req;
#endif

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    if(nb.rank == Globals::my_rank) // on the same process
      pbl=pmb->pmy_mesh->FindMeshBlock(nb.gid);
    switch(type) {
      case BND_MGGRAV:
        sbuf=mggrav_send_[nb.bufid];
//        if(nb.level<mylevel) cbuf.InitWithShallowCopy(pmb->pmr->coarse_cons_);
        nvar=1, ngh=2;
        tag=CreateBvalsMPITag(nb.lid, TAG_MGGRAV, nb.target);
        if(nb.rank == Globals::my_rank) {
          rbuf=pbl->pbval->mggrav_recv_[nb.targetid];
          flag=&(pbl->pbval->mggrav_flag_[nb.targetid]);
        }
#ifdef MPI_PARALLEL
        else
          req=&(req_mggrav_send_[nb.bufid]);
#endif
        break;
      defualt:
        break;
    }
    int ssize;
    if(nb.level==mylevel)
      ssize=LoadMultigridBoundaryBufferSameLevel(src, nvar, nc, ngh, sbuf, nb);
//    else if(nb.level<mylevel)
//      ssize=LoadMultigridBoundaryBufferToCoarser(src, nvar, sbuf, cbuf, nb);
//    else
//      ssize=LoadMultigridBoundaryBufferToFiner(src, nvar, sbuf, nb);
    if(nb.rank == Globals::my_rank) {
      std::memcpy(rbuf, sbuf, ssize*sizeof(Real));
      *flag=BNDRY_ARRIVED;
    }
#ifdef MPI_PARALLEL
    else // MPI
      MPI_Isend(sbuf,ssize,MPI_ATHENA_REAL,nb.rank,tag,MPI_COMM_WORLD,req);
#endif
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetMultigridBoundarySameLevel(AthenaArray<Real> &dst,
//                              int nvar, int nc, int ngh, Real *buf, const NeighborBlock& nb)
//  \brief Set hydro boundary received from a block on the same level

void BoundaryValues::SetMultigridBoundarySameLevel(AthenaArray<Real> &dst,
                        int nvar, int nc, int ngh, Real *buf, const NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_block_;
  int si, sj, sk, ei, ej, ek;

  if(nb.ox1==0)     si=ngh,    ei=nc+ngh-1;
  else if(nb.ox1>0) si=nc+ngh, ei=nc+2*ngh-1;
  else              si=0,      ei=ngh-1;
  if(nb.ox2==0)     sj=ngh,    ej=nc+ngh-1;
  else if(nb.ox2>0) sj=nc+ngh, ej=nc+2*ngh-1;
  else              sj=0,      ej=ngh-1;
  if(nb.ox3==0)     sk=ngh,    ek=nc+ngh-1;
  else if(nb.ox3>0) sk=nc+ngh, ek=nc+2*ngh-1;
  else              sk=0,      ek=ngh-1;

  int p=0;
  BufferUtility::Unpack4DData(buf, dst, 0, nvar-1, si, ei, sj, ej, sk, ek, p);
  return;
}


//----------------------------------------------------------------------------------------
//! \fn bool BoundaryValues::ReceiveMultigridBoundaryBuffers(AthenaArray<Real> &dst,
//                                                   int nc, enum MGBoundaryType type)
//  \brief receive the boundary data

bool BoundaryValues::ReceiveMultigridBoundaryBuffers(AthenaArray<Real> &dst,
                                                     int nc, enum MGBoundaryType type)
{
  MeshBlock *pmb=pmy_block_;
  bool bflag=true;
  bool *flip=NULL;
  Real *rbuf;
//  AthenaArray<Real> cbuf;
  int nvar, ngh;
  enum BoundaryStatus *flag;
#ifdef MPI_PARALLEL
  MPI_Request *req;
#endif

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    switch(type) {
      case BND_MGGRAV:
        nvar=1, ngh=2;
        rbuf=mggrav_recv_[nb.bufid];
//        if(nb.level<pmb->loc.level) cbuf.InitWithShallowCopy(pmb->pmr->coarse_cons_);
        flag=&(mggrav_flag_[nb.bufid]);
#ifdef MPI_PARALLEL
        if(nb.rank!=Globals::my_rank)
          req=&(req_mggrav_recv_[nb.bufid]);
#endif
        break;
      defualt:
        break;
    }
    if(*flag==BNDRY_COMPLETED) continue;
    if(*flag==BNDRY_WAITING) {
      if(nb.rank==Globals::my_rank) {// on the same process
        bflag=false;
        continue;
      }
#ifdef MPI_PARALLEL
      else { // MPI boundary
        int test;
        MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&test,MPI_STATUS_IGNORE);
        MPI_Test(req,&test,MPI_STATUS_IGNORE);
        if(test==false) {
          bflag=false;
          continue;
        }
        *flag = BNDRY_ARRIVED;
      }
#endif
    }
    if(nb.level==pmb->loc.level)
      SetMultigridBoundarySameLevel(dst, nvar, nc, ngh, rbuf, nb);
//    else if(nb.level<pmb->loc.level) // this set only the prolongation buffer
//      SetMultigridBoundaryFromCoarser(nvar, nc, ngh, rbuf, cbuf, nb);
//    else
//      SetMultigridBoundaryFromFiner(dst, nvar, nc, ngh, rbuf, nb);
    *flag = BNDRY_COMPLETED; // completed
  }
  return bflag;
}


