//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_cc.cpp
//  \brief functions that apply BCs for CELL_CENTERED variables

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
//! \fn int BoundaryValues::LoadHydroBoundaryBufferSameLevel(AthenaArray<Real> &src,
//                                                 Real *buf, const NeighborBlock& nb)
//  \brief Set hydro boundary buffers for sending to a block on the same level

int BoundaryValues::LoadHydroBoundaryBufferSameLevel(AthenaArray<Real> &src, Real *buf,
                                                     const NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_block_;
  int si, sj, sk, ei, ej, ek;

  si=(nb.ox1>0)?(pmb->ie-NGHOST+1):pmb->is;
  ei=(nb.ox1<0)?(pmb->is+NGHOST-1):pmb->ie;
  sj=(nb.ox2>0)?(pmb->je-NGHOST+1):pmb->js;
  ej=(nb.ox2<0)?(pmb->js+NGHOST-1):pmb->je;
  sk=(nb.ox3>0)?(pmb->ke-NGHOST+1):pmb->ks;
  ek=(nb.ox3<0)?(pmb->ks+NGHOST-1):pmb->ke;
  int p=0;
  BufferUtility::Pack4DData(src, buf, 0, NHYDRO-1, si, ei, sj, ej, sk, ek, p);
  return p;
}

//----------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadHydroBoundaryBufferToCoarser(AthenaArray<Real> &src,
//                                                 Real *buf, const NeighborBlock& nb,
//                                                 bool conserved_values)
//  \brief Set hydro boundary buffers for sending to a block on the coarser level

int BoundaryValues::LoadHydroBoundaryBufferToCoarser(AthenaArray<Real> &src, Real *buf,
                                                     const NeighborBlock& nb,
                                                     bool conserved_values)
{
  MeshBlock *pmb=pmy_block_;
  MeshRefinement *pmr=pmb->pmr;
  int si, sj, sk, ei, ej, ek;
  int cn=pmb->cnghost-1;

  si=(nb.ox1>0)?(pmb->cie-cn):pmb->cis;
  ei=(nb.ox1<0)?(pmb->cis+cn):pmb->cie;
  sj=(nb.ox2>0)?(pmb->cje-cn):pmb->cjs;
  ej=(nb.ox2<0)?(pmb->cjs+cn):pmb->cje;
  sk=(nb.ox3>0)?(pmb->cke-cn):pmb->cks;
  ek=(nb.ox3<0)?(pmb->cks+cn):pmb->cke;

  int p=0;
  if (conserved_values) { // normal case
    pmr->RestrictCellCenteredValues(src, pmr->coarse_cons_, 0, NHYDRO-1,
                                    si, ei, sj, ej, sk, ek);
    BufferUtility::Pack4DData(pmr->coarse_cons_, buf, 0, NHYDRO-1,
                              si, ei, sj, ej, sk, ek, p);
  }
  else { // must be initialization of boundary primitives
    pmr->RestrictCellCenteredValues(src, pmr->coarse_prim_, 0, NHYDRO-1,
                                    si, ei, sj, ej, sk, ek);
    BufferUtility::Pack4DData(pmr->coarse_prim_, buf, 0, NHYDRO-1,
                              si, ei, sj, ej, sk, ek, p);
  }
  return p;
}

//----------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadHydroBoundaryBufferToFiner(AthenaArray<Real> &src,
//                                                 Real *buf, const NeighborBlock& nb)
//  \brief Set hydro boundary buffers for sending to a block on the finer level

int BoundaryValues::LoadHydroBoundaryBufferToFiner(AthenaArray<Real> &src, Real *buf,
                                                   const NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_block_;
  int si, sj, sk, ei, ej, ek;
  int cn=pmb->cnghost-1;

  si=(nb.ox1>0)?(pmb->ie-cn):pmb->is;
  ei=(nb.ox1<0)?(pmb->is+cn):pmb->ie;
  sj=(nb.ox2>0)?(pmb->je-cn):pmb->js;
  ej=(nb.ox2<0)?(pmb->js+cn):pmb->je;
  sk=(nb.ox3>0)?(pmb->ke-cn):pmb->ks;
  ek=(nb.ox3<0)?(pmb->ks+cn):pmb->ke;

  // send the data first and later prolongate on the target block
  // need to add edges for faces, add corners for edges
  if(nb.ox1==0) {
    if(nb.fi1==1)   si+=pmb->block_size.nx1/2-pmb->cnghost;
    else            ei-=pmb->block_size.nx1/2-pmb->cnghost;
  }
  if(nb.ox2==0 && pmb->block_size.nx2 > 1) {
    if(nb.ox1!=0) {
      if(nb.fi1==1) sj+=pmb->block_size.nx2/2-pmb->cnghost;
      else          ej-=pmb->block_size.nx2/2-pmb->cnghost;
    }
    else {
      if(nb.fi2==1) sj+=pmb->block_size.nx2/2-pmb->cnghost;
      else          ej-=pmb->block_size.nx2/2-pmb->cnghost;
    }
  }
  if(nb.ox3==0 && pmb->block_size.nx3 > 1) {
    if(nb.ox1!=0 && nb.ox2!=0) {
      if(nb.fi1==1) sk+=pmb->block_size.nx3/2-pmb->cnghost;
      else          ek-=pmb->block_size.nx3/2-pmb->cnghost;
    }
    else {
      if(nb.fi2==1) sk+=pmb->block_size.nx3/2-pmb->cnghost;
      else          ek-=pmb->block_size.nx3/2-pmb->cnghost;
    }
  }

  int p=0;
  BufferUtility::Pack4DData(src, buf, 0, NHYDRO-1, si, ei, sj, ej, sk, ek, p);
  return p;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SendHydroBoundaryBuffers(AthenaArray<Real> &src,
//                                                    bool conserved_values)
//  \brief Send boundary buffers

void BoundaryValues::SendHydroBoundaryBuffers(AthenaArray<Real> &src,
                                              bool conserved_values)
{
  MeshBlock *pmb=pmy_block_;
  int mylevel=pmb->loc.level;

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    int ssize;
    if(nb.level==mylevel)
      ssize=LoadHydroBoundaryBufferSameLevel(src, hydro_send_[nb.bufid],nb);
    else if(nb.level<mylevel)
      ssize=LoadHydroBoundaryBufferToCoarser(src, hydro_send_[nb.bufid], nb,
                                             conserved_values);
    else
      ssize=LoadHydroBoundaryBufferToFiner(src, hydro_send_[nb.bufid], nb);
    if(nb.rank == Globals::my_rank) { // on the same process
      MeshBlock *pbl=pmb->pmy_mesh->FindMeshBlock(nb.gid);
      std::memcpy(pbl->pbval->hydro_recv_[nb.targetid],
                  hydro_send_[nb.bufid], ssize*sizeof(Real));
      pbl->pbval->hydro_flag_[nb.targetid]=BNDRY_ARRIVED;
    }
#ifdef MPI_PARALLEL
    else // MPI
      MPI_Start(&req_hydro_send_[nb.bufid]);
#endif
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetHydroBoundarySameLevel(AthenaArray<Real> &dst,
//                                           Real *buf, const NeighborBlock& nb)
//  \brief Set hydro boundary received from a block on the same level

void BoundaryValues::SetHydroBoundarySameLevel(AthenaArray<Real> &dst, Real *buf,
                                               const NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_block_;
  int si, sj, sk, ei, ej, ek;

  if(nb.ox1==0)     si=pmb->is,        ei=pmb->ie;
  else if(nb.ox1>0) si=pmb->ie+1,      ei=pmb->ie+NGHOST;
  else              si=pmb->is-NGHOST, ei=pmb->is-1;
  if(nb.ox2==0)     sj=pmb->js,        ej=pmb->je;
  else if(nb.ox2>0) sj=pmb->je+1,      ej=pmb->je+NGHOST;
  else              sj=pmb->js-NGHOST, ej=pmb->js-1;
  if(nb.ox3==0)     sk=pmb->ks,        ek=pmb->ke;
  else if(nb.ox3>0) sk=pmb->ke+1,      ek=pmb->ke+NGHOST;
  else              sk=pmb->ks-NGHOST, ek=pmb->ks-1;

  int p=0;
  if (nb.polar) {
    for (int n=0; n<(NHYDRO); ++n) {
      Real sign = flip_across_pole_hydro[n] ? -1.0 : 1.0;
      for (int k=sk; k<=ek; ++k) {
        for (int j=ej; j>=sj; --j) {
#pragma simd
          for (int i=si; i<=ei; ++i)
            dst(n,k,j,i) = sign * buf[p++];
        }
      }
    }
  }
  else
    BufferUtility::Unpack4DData(buf, dst, 0, NHYDRO-1, si, ei, sj, ej, sk, ek, p);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetHydroBoundaryFromCoarser(Real *buf,
//                                                       const NeighborBlock& nb,
//                                                       bool conserved_values)
//  \brief Set hydro prolongation buffer received from a block on a coarser level

void BoundaryValues::SetHydroBoundaryFromCoarser(Real *buf, const NeighborBlock& nb,
    bool conserved_values)
{
  MeshBlock *pmb=pmy_block_;
  MeshRefinement *pmr=pmb->pmr;

  int si, sj, sk, ei, ej, ek;
  int cng=pmb->cnghost;

  if(nb.ox1==0) {
    si=pmb->cis, ei=pmb->cie;
    if((pmb->loc.lx1&1L)==0L) ei+=cng;
    else             si-=cng;
  }
  else if(nb.ox1>0)  si=pmb->cie+1,   ei=pmb->cie+cng;
  else               si=pmb->cis-cng, ei=pmb->cis-1;
  if(nb.ox2==0) {
    sj=pmb->cjs, ej=pmb->cje;
    if(pmb->block_size.nx2 > 1) {
      if((pmb->loc.lx2&1L)==0L) ej+=cng;
      else             sj-=cng;
    }
  }
  else if(nb.ox2>0)  sj=pmb->cje+1,   ej=pmb->cje+cng;
  else               sj=pmb->cjs-cng, ej=pmb->cjs-1;
  if(nb.ox3==0) {
    sk=pmb->cks, ek=pmb->cke;
    if(pmb->block_size.nx3 > 1) {
      if((pmb->loc.lx3&1L)==0L) ek+=cng;
      else             sk-=cng;
    }
  }
  else if(nb.ox3>0)  sk=pmb->cke+1,   ek=pmb->cke+cng;
  else               sk=pmb->cks-cng, ek=pmb->cks-1;

  int p=0;
  if (nb.polar) {
    for (int n=0; n<(NHYDRO); ++n) {
      Real sign = flip_across_pole_hydro[n] ? -1.0 : 1.0;
      for (int k=sk; k<=ek; ++k) {
        for (int j=ej; j>=sj; --j) {
#pragma simd
          for (int i=si; i<=ei; ++i) {
            if (conserved_values)
              pmr->coarse_cons_(n,k,j,i) = sign * buf[p++];
            else
              pmr->coarse_prim_(n,k,j,i) = sign * buf[p++];
          }
        }
      }
    }
  }
  else {
    if (conserved_values)
      BufferUtility::Unpack4DData(buf, pmr->coarse_cons_, 0, NHYDRO-1,
                                  si, ei, sj, ej, sk, ek, p);
    else
      BufferUtility::Unpack4DData(buf, pmr->coarse_prim_, 0, NHYDRO-1,
                                  si, ei, sj, ej, sk, ek, p);
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetHydroBoundaryFromFiner(AthenaArray<Real> &dst,
//                                               Real *buf, const NeighborBlock& nb)
//  \brief Set hydro boundary received from a block on a finer level

void BoundaryValues::SetHydroBoundaryFromFiner(AthenaArray<Real> &dst, Real *buf,
                                               const NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_block_;
  // receive already restricted data
  int si, sj, sk, ei, ej, ek;

  if(nb.ox1==0) {
    si=pmb->is, ei=pmb->ie;
    if(nb.fi1==1)   si+=pmb->block_size.nx1/2;
    else            ei-=pmb->block_size.nx1/2;
  }
  else if(nb.ox1>0) si=pmb->ie+1,      ei=pmb->ie+NGHOST;
  else              si=pmb->is-NGHOST, ei=pmb->is-1;
  if(nb.ox2==0) {
    sj=pmb->js, ej=pmb->je;
    if(pmb->block_size.nx2 > 1) {
      if(nb.ox1!=0) {
        if(nb.fi1==1) sj+=pmb->block_size.nx2/2;
        else          ej-=pmb->block_size.nx2/2;
      }
      else {
        if(nb.fi2==1) sj+=pmb->block_size.nx2/2;
        else          ej-=pmb->block_size.nx2/2;
      }
    }
  }
  else if(nb.ox2>0) sj=pmb->je+1,      ej=pmb->je+NGHOST;
  else              sj=pmb->js-NGHOST, ej=pmb->js-1;
  if(nb.ox3==0) {
    sk=pmb->ks, ek=pmb->ke;
    if(pmb->block_size.nx3 > 1) {
      if(nb.ox1!=0 && nb.ox2!=0) {
        if(nb.fi1==1) sk+=pmb->block_size.nx3/2;
        else          ek-=pmb->block_size.nx3/2;
      }
      else {
        if(nb.fi2==1) sk+=pmb->block_size.nx3/2;
        else          ek-=pmb->block_size.nx3/2;
      }
    }
  }
  else if(nb.ox3>0) sk=pmb->ke+1,      ek=pmb->ke+NGHOST;
  else              sk=pmb->ks-NGHOST, ek=pmb->ks-1;

  int p=0;
  if (nb.polar) {
    for (int n=0; n<(NHYDRO); ++n) {
      Real sign = flip_across_pole_hydro[n] ? -1.0 : 1.0;
      for (int k=sk; k<=ek; ++k) {
        for (int j=sj; j<=ej; ++j) {
#pragma simd
          for (int i=si; i<=ei; ++i)
            dst(n,k,j,i) = sign * buf[p++];
        }
      }
    }
  }
  else
    BufferUtility::Unpack4DData(buf, dst, 0, NHYDRO-1, si, ei, sj, ej, sk, ek, p);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn bool BoundaryValues::ReceiveHydroBoundaryBuffers(AthenaArray<Real> &dst)
//  \brief receive the boundary data

bool BoundaryValues::ReceiveHydroBoundaryBuffers(AthenaArray<Real> &dst)
{
  MeshBlock *pmb=pmy_block_;
  bool flag=true;

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    if(hydro_flag_[nb.bufid]==BNDRY_COMPLETED) continue;
    if(hydro_flag_[nb.bufid]==BNDRY_WAITING) {
      if(nb.rank==Globals::my_rank) {// on the same process
        flag=false;
        continue;
      }
#ifdef MPI_PARALLEL
      else { // MPI boundary
        int test;
        MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&test,MPI_STATUS_IGNORE);
        MPI_Test(&req_hydro_recv_[nb.bufid],&test,MPI_STATUS_IGNORE);
        if(test==false) {
          flag=false;
          continue;
        }
        hydro_flag_[nb.bufid] = BNDRY_ARRIVED;
      }
#endif
    }
    if(nb.level==pmb->loc.level)
      SetHydroBoundarySameLevel(dst, hydro_recv_[nb.bufid], nb);
    else if(nb.level<pmb->loc.level) // this set only the prolongation buffer
      SetHydroBoundaryFromCoarser(hydro_recv_[nb.bufid], nb, true);
    else
      SetHydroBoundaryFromFiner(dst, hydro_recv_[nb.bufid], nb);
    hydro_flag_[nb.bufid] = BNDRY_COMPLETED; // completed
  }

  if(flag&& (pmb->block_bcs[INNER_X2]==POLAR_BNDRY
         ||  pmb->block_bcs[OUTER_X2]==POLAR_BNDRY))
     PolarSingleHydro(dst);
  return flag;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ReceiveHydroBoundaryBuffersWithWait(AthenaArray<Real> &dst,
//                                                               bool conserved_values)
//  \brief receive the boundary data for initialization

void BoundaryValues::ReceiveHydroBoundaryBuffersWithWait(AthenaArray<Real> &dst,
                                                         bool conserved_values)
{
  MeshBlock *pmb=pmy_block_;

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb = pmb->neighbor[n];
#ifdef MPI_PARALLEL
    if(nb.rank!=Globals::my_rank)
      MPI_Wait(&req_hydro_recv_[nb.bufid],MPI_STATUS_IGNORE);
#endif
    if(nb.level==pmb->loc.level)
      SetHydroBoundarySameLevel(dst, hydro_recv_[nb.bufid], nb);
    else if(nb.level<pmb->loc.level)
      SetHydroBoundaryFromCoarser(hydro_recv_[nb.bufid], nb, conserved_values);
    else
      SetHydroBoundaryFromFiner(dst, hydro_recv_[nb.bufid], nb);
    hydro_flag_[nb.bufid] = BNDRY_COMPLETED; // completed
  }

  if (pmb->block_bcs[INNER_X2]==POLAR_BNDRY||pmb->block_bcs[OUTER_X2]==POLAR_BNDRY)
    PolarSingleHydro(dst);

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::PolarSingleHydro(AthenaArray<Real> &dst)
//
// \brief  single CPU in the azimuthal direction for the polar boundary

void BoundaryValues::PolarSingleHydro(AthenaArray<Real> &dst)
{
  MeshBlock *pmb=pmy_block_;
  if(pmb->loc.level == pmb->pmy_mesh->root_level && pmb->pmy_mesh->nrbx3 == 1){

    if(pmb->block_bcs[INNER_X2]==POLAR_BNDRY){
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      for (int n=0; n<(NHYDRO); ++n) {
        for (int j=pmb->js-NGHOST; j<=pmb->js-1; ++j) {
         for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i){
           for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
             exc_(k)=dst(n,k,j,i);
           }
           for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
             int k_shift = k;
             k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
             dst(n,k,j,i)=exc_(k_shift);
           }
         }
        }
      }
    }

    if(pmb->block_bcs[OUTER_X2]==POLAR_BNDRY){
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      for (int n=0; n<(NHYDRO); ++n) {
        for (int j=pmb->je+1; j<=pmb->je+NGHOST; ++j) {
         for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i){
           for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
             exc_(k)=dst(n,k,j,i);
           }
           for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
             int k_shift = k;
             k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
             dst(n,k,j,i)=exc_(k_shift);
           }
         }
        }
      }
    }
  }
  return;
}
