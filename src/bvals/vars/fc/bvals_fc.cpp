//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_fc.cpp
//  \brief functions that apply BCs for FACE_CENTERED variables

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
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../utils/buffer_utils.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadFieldBoundaryBufferSameLevel(FaceField &src,
//                                                 Real *buf, const NeighborBlock& nb)
//  \brief Set field boundary buffers for sending to a block on the same level

int BoundaryValues::LoadFieldBoundaryBufferSameLevel(FaceField &src, Real *buf,
                                                     const NeighborBlock& nb) {
  MeshBlock *pmb=pmy_block_;
  int si, sj, sk, ei, ej, ek;
  int p=0;

  // bx1
  if (nb.ox1==0)     si=pmb->is,          ei=pmb->ie+1;
  else if (nb.ox1>0) si=pmb->ie-NGHOST+1, ei=pmb->ie;
  else              si=pmb->is+1,        ei=pmb->is+NGHOST;
  if (nb.ox2==0)     sj=pmb->js,          ej=pmb->je;
  else if (nb.ox2>0) sj=pmb->je-NGHOST+1, ej=pmb->je;
  else              sj=pmb->js,          ej=pmb->js+NGHOST-1;
  if (nb.ox3==0)     sk=pmb->ks,          ek=pmb->ke;
  else if (nb.ox3>0) sk=pmb->ke-NGHOST+1, ek=pmb->ke;
  else              sk=pmb->ks,          ek=pmb->ks+NGHOST-1;
  // for SMR/AMR, always include the overlapping faces in edge and corner boundaries
  if (pmb->pmy_mesh->multilevel==true && nb.type != NEIGHBOR_FACE) {
    if (nb.ox1>0) ei++;
    else if (nb.ox1<0) si--;
  }
  BufferUtility::Pack3DData(src.x1f, buf, si, ei, sj, ej, sk, ek, p);

  // bx2
  if (nb.ox1==0)      si=pmb->is,          ei=pmb->ie;
  else if (nb.ox1>0)  si=pmb->ie-NGHOST+1, ei=pmb->ie;
  else               si=pmb->is,          ei=pmb->is+NGHOST-1;
  if (pmb->block_size.nx2==1) sj=pmb->js,  ej=pmb->je;
  else if (nb.ox2==0) sj=pmb->js,          ej=pmb->je+1;
  else if (nb.ox2>0)  sj=pmb->je-NGHOST+1, ej=pmb->je;
  else               sj=pmb->js+1,        ej=pmb->js+NGHOST;
  if (pmb->pmy_mesh->multilevel==true && nb.type != NEIGHBOR_FACE) {
    if (nb.ox2>0) ej++;
    else if (nb.ox2<0) sj--;
  }
  BufferUtility::Pack3DData(src.x2f, buf, si, ei, sj, ej, sk, ek, p);

  // bx3
  if (nb.ox2==0)      sj=pmb->js,          ej=pmb->je;
  else if (nb.ox2>0)  sj=pmb->je-NGHOST+1, ej=pmb->je;
  else               sj=pmb->js,          ej=pmb->js+NGHOST-1;
  if (pmb->block_size.nx3==1) sk=pmb->ks,  ek=pmb->ke;
  else if (nb.ox3==0) sk=pmb->ks,          ek=pmb->ke+1;
  else if (nb.ox3>0)  sk=pmb->ke-NGHOST+1, ek=pmb->ke;
  else               sk=pmb->ks+1,        ek=pmb->ks+NGHOST;
  if (pmb->pmy_mesh->multilevel==true && nb.type != NEIGHBOR_FACE) {
    if (nb.ox3>0) ek++;
    else if (nb.ox3<0) sk--;
  }
  BufferUtility::Pack3DData(src.x3f, buf, si, ei, sj, ej, sk, ek, p);

  return p;
}

//----------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadFieldBoundaryBufferToCoarser(FaceField &src,
//                                                 Real *buf, const NeighborBlock& nb)
//  \brief Set field boundary buffers for sending to a block on the coarser level

int BoundaryValues::LoadFieldBoundaryBufferToCoarser(FaceField &src, Real *buf,
                                                     const NeighborBlock& nb) {
  MeshBlock *pmb=pmy_block_;
  MeshRefinement *pmr=pmb->pmr;
  int si, sj, sk, ei, ej, ek;
  int cng=NGHOST;
  int p=0;

  // bx1
  if (nb.ox1==0)     si=pmb->cis,       ei=pmb->cie+1;
  else if (nb.ox1>0) si=pmb->cie-cng+1, ei=pmb->cie;
  else              si=pmb->cis+1,     ei=pmb->cis+cng;
  if (nb.ox2==0)     sj=pmb->cjs,       ej=pmb->cje;
  else if (nb.ox2>0) sj=pmb->cje-cng+1, ej=pmb->cje;
  else              sj=pmb->cjs,       ej=pmb->cjs+cng-1;
  if (nb.ox3==0)     sk=pmb->cks,       ek=pmb->cke;
  else if (nb.ox3>0) sk=pmb->cke-cng+1, ek=pmb->cke;
  else              sk=pmb->cks,       ek=pmb->cks+cng-1;
  // include the overlapping faces in edge and corner boundaries
  if (nb.type != NEIGHBOR_FACE) {
    if (nb.ox1>0) ei++;
    else if (nb.ox1<0) si--;
  }
  pmr->RestrictFieldX1(src.x1f, pmr->coarse_b_.x1f, si, ei, sj, ej, sk, ek);
  BufferUtility::Pack3DData(pmr->coarse_b_.x1f, buf, si, ei, sj, ej, sk, ek, p);

  // bx2
  if (nb.ox1==0)      si=pmb->cis,       ei=pmb->cie;
  else if (nb.ox1>0)  si=pmb->cie-cng+1, ei=pmb->cie;
  else               si=pmb->cis,       ei=pmb->cis+cng-1;
  if (pmb->block_size.nx2==1) sj=pmb->cjs, ej=pmb->cje;
  else if (nb.ox2==0) sj=pmb->cjs,       ej=pmb->cje+1;
  else if (nb.ox2>0)  sj=pmb->cje-cng+1, ej=pmb->cje;
  else               sj=pmb->cjs+1,     ej=pmb->cjs+cng;
  if (nb.type != NEIGHBOR_FACE) {
    if (nb.ox2>0) ej++;
    else if (nb.ox2<0) sj--;
  }
  pmr->RestrictFieldX2(src.x2f, pmr->coarse_b_.x2f, si, ei, sj, ej, sk, ek);
  if (pmb->block_size.nx2==1) { // 1D
    for (int i=si; i<=ei; i++)
      pmr->coarse_b_.x2f(sk,sj+1,i)=pmr->coarse_b_.x2f(sk,sj,i);
  }
  BufferUtility::Pack3DData(pmr->coarse_b_.x2f, buf, si, ei, sj, ej, sk, ek, p);

  // bx3
  if (nb.ox2==0)      sj=pmb->cjs,       ej=pmb->cje;
  else if (nb.ox2>0)  sj=pmb->cje-cng+1, ej=pmb->cje;
  else               sj=pmb->cjs,       ej=pmb->cjs+cng-1;
  if (pmb->block_size.nx3==1) sk=pmb->cks,  ek=pmb->cke;
  else if (nb.ox3==0) sk=pmb->cks,       ek=pmb->cke+1;
  else if (nb.ox3>0)  sk=pmb->cke-cng+1, ek=pmb->cke;
  else               sk=pmb->cks+1,     ek=pmb->cks+cng;
  if (nb.type != NEIGHBOR_FACE) {
    if (nb.ox3>0) ek++;
    else if (nb.ox3<0) sk--;
  }
  pmr->RestrictFieldX3(src.x3f, pmr->coarse_b_.x3f, si, ei, sj, ej, sk, ek);
  if (pmb->block_size.nx3==1) { // 1D or 2D
    for (int j=sj; j<=ej; j++) {
      for (int i=si; i<=ei; i++)
        pmr->coarse_b_.x3f(sk+1,j,i)=pmr->coarse_b_.x3f(sk,j,i);
    }
  }
  BufferUtility::Pack3DData(pmr->coarse_b_.x3f, buf, si, ei, sj, ej, sk, ek, p);

  return p;
}

//----------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadFieldBoundaryBufferToFiner(FaceField &src,
//                                                 Real *buf, const NeighborBlock& nb)
//  \brief Set field boundary buffers for sending to a block on the finer level

int BoundaryValues::LoadFieldBoundaryBufferToFiner(FaceField &src, Real *buf,
                                                   const NeighborBlock& nb) {
  MeshBlock *pmb=pmy_block_;
  int si, sj, sk, ei, ej, ek;
  int cn=pmb->cnghost-1;
  int p=0;

  // send the data first and later prolongate on the target block
  // need to add edges for faces, add corners for edges
  // bx1
  if (nb.ox1==0) {
    if (nb.fi1==1)   si=pmb->is+pmb->block_size.nx1/2-pmb->cnghost, ei=pmb->ie+1;
    else            si=pmb->is, ei=pmb->ie+1-pmb->block_size.nx1/2+pmb->cnghost;
  } else if (nb.ox1>0) { si=pmb->ie+1-pmb->cnghost, ei=pmb->ie+1;}
  else              si=pmb->is,                ei=pmb->is+pmb->cnghost;
  if (nb.ox2==0) {
    sj=pmb->js,    ej=pmb->je;
    if (pmb->block_size.nx2 > 1) {
      if (nb.ox1!=0) {
        if (nb.fi1==1) sj+=pmb->block_size.nx2/2-pmb->cnghost;
        else          ej-=pmb->block_size.nx2/2-pmb->cnghost;
      } else {
        if (nb.fi2==1) sj+=pmb->block_size.nx2/2-pmb->cnghost;
        else          ej-=pmb->block_size.nx2/2-pmb->cnghost;
      }
    }
  } else if (nb.ox2>0) { sj=pmb->je-cn, ej=pmb->je;}
  else              sj=pmb->js,    ej=pmb->js+cn;
  if (nb.ox3==0) {
    sk=pmb->ks,    ek=pmb->ke;
    if (pmb->block_size.nx3 > 1) {
      if (nb.ox1!=0 && nb.ox2!=0) {
        if (nb.fi1==1) sk+=pmb->block_size.nx3/2-pmb->cnghost;
        else          ek-=pmb->block_size.nx3/2-pmb->cnghost;
      } else {
        if (nb.fi2==1) sk+=pmb->block_size.nx3/2-pmb->cnghost;
        else          ek-=pmb->block_size.nx3/2-pmb->cnghost;
      }
    }
  } else if (nb.ox3>0) { sk=pmb->ke-cn, ek=pmb->ke;}
  else              sk=pmb->ks,    ek=pmb->ks+cn;
  BufferUtility::Pack3DData(src.x1f, buf, si, ei, sj, ej, sk, ek, p);

  // bx2
  if (nb.ox1==0) {
    if (nb.fi1==1)   si=pmb->is+pmb->block_size.nx1/2-pmb->cnghost, ei=pmb->ie;
    else            si=pmb->is, ei=pmb->ie-pmb->block_size.nx1/2+pmb->cnghost;
  } else if (nb.ox1>0) { si=pmb->ie-cn, ei=pmb->ie;}
  else              si=pmb->is,    ei=pmb->is+cn;
  if (nb.ox2==0) {
    sj=pmb->js,    ej=pmb->je;
    if (pmb->block_size.nx2 > 1) {
      ej++;
      if (nb.ox1!=0) {
        if (nb.fi1==1) sj+=pmb->block_size.nx2/2-pmb->cnghost;
        else          ej-=pmb->block_size.nx2/2-pmb->cnghost;
      } else {
        if (nb.fi2==1) sj+=pmb->block_size.nx2/2-pmb->cnghost;
        else          ej-=pmb->block_size.nx2/2-pmb->cnghost;
      }
    }
  } else if (nb.ox2>0) { sj=pmb->je+1-pmb->cnghost, ej=pmb->je+1;}
  else              sj=pmb->js,                ej=pmb->js+pmb->cnghost;
  BufferUtility::Pack3DData(src.x2f, buf, si, ei, sj, ej, sk, ek, p);

  // bx3
  if (nb.ox2==0) {
    sj=pmb->js,    ej=pmb->je;
    if (pmb->block_size.nx2 > 1) {
      if (nb.ox1!=0) {
        if (nb.fi1==1) sj+=pmb->block_size.nx2/2-pmb->cnghost;
        else          ej-=pmb->block_size.nx2/2-pmb->cnghost;
      } else {
        if (nb.fi2==1) sj+=pmb->block_size.nx2/2-pmb->cnghost;
        else          ej-=pmb->block_size.nx2/2-pmb->cnghost;
      }
    }
  } else if (nb.ox2>0) { sj=pmb->je-cn, ej=pmb->je;}
  else              sj=pmb->js,    ej=pmb->js+cn;
  if (nb.ox3==0) {
    sk=pmb->ks,    ek=pmb->ke;
    if (pmb->block_size.nx3 > 1) {
      ek++;
      if (nb.ox1!=0 && nb.ox2!=0) {
        if (nb.fi1==1) sk+=pmb->block_size.nx3/2-pmb->cnghost;
        else          ek-=pmb->block_size.nx3/2-pmb->cnghost;
      } else {
        if (nb.fi2==1) sk+=pmb->block_size.nx3/2-pmb->cnghost;
        else          ek-=pmb->block_size.nx3/2-pmb->cnghost;
      }
    }
  } else if (nb.ox3>0) { sk=pmb->ke+1-pmb->cnghost, ek=pmb->ke+1;}
  else              sk=pmb->ks,                ek=pmb->ks+pmb->cnghost;
  BufferUtility::Pack3DData(src.x3f, buf, si, ei, sj, ej, sk, ek, p);

  return p;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SendFieldBoundaryBuffers(FaceField &src)
//  \brief Send field boundary buffers

void BoundaryValues::SendFieldBoundaryBuffers(FaceField &src) {
  MeshBlock *pmb=pmy_block_;

  for (int n=0; n<nneighbor; n++) {
    NeighborBlock& nb = neighbor[n];
    int ssize;
    if (nb.level==pmb->loc.level)
      ssize=LoadFieldBoundaryBufferSameLevel(src, bd_field_.send[nb.bufid],nb);
    else if (nb.level<pmb->loc.level)
      ssize=LoadFieldBoundaryBufferToCoarser(src, bd_field_.send[nb.bufid],nb);
    else
      ssize=LoadFieldBoundaryBufferToFiner(src, bd_field_.send[nb.bufid], nb);
    if (nb.rank == Globals::my_rank) { // on the same process
      MeshBlock *pbl=pmb->pmy_mesh->FindMeshBlock(nb.gid);
      // find target buffer
      std::memcpy(pbl->pbval->bd_field_.recv[nb.targetid],
                  bd_field_.send[nb.bufid], ssize*sizeof(Real));
      pbl->pbval->bd_field_.flag[nb.targetid]=BNDRY_ARRIVED;
    }
#ifdef MPI_PARALLEL
    else // MPI
      MPI_Start(&(bd_field_.req_send[nb.bufid]));
#endif
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetFieldBoundarySameLevel(FaceField &dst,
//                                               Real *buf, const NeighborBlock& nb)
//  \brief Set field boundary received from a block on the same level

void BoundaryValues::SetFieldBoundarySameLevel(FaceField &dst, Real *buf,
                                               const NeighborBlock& nb) {
  MeshBlock *pmb=pmy_block_;
  int si, sj, sk, ei, ej, ek;

  int p=0;
  // bx1
  // for uniform grid: face-neighbors take care of the overlapping faces
  if (nb.ox1==0)     si=pmb->is,        ei=pmb->ie+1;
  else if (nb.ox1>0) si=pmb->ie+2,      ei=pmb->ie+NGHOST+1;
  else              si=pmb->is-NGHOST, ei=pmb->is-1;
  if (nb.ox2==0)     sj=pmb->js,        ej=pmb->je;
  else if (nb.ox2>0) sj=pmb->je+1,      ej=pmb->je+NGHOST;
  else              sj=pmb->js-NGHOST, ej=pmb->js-1;
  if (nb.ox3==0)     sk=pmb->ks,        ek=pmb->ke;
  else if (nb.ox3>0) sk=pmb->ke+1,      ek=pmb->ke+NGHOST;
  else              sk=pmb->ks-NGHOST, ek=pmb->ks-1;
  // for SMR/AMR, always include the overlapping faces in edge and corner boundaries
  if (pmb->pmy_mesh->multilevel==true && nb.type != NEIGHBOR_FACE) {
    if (nb.ox1>0) si--;
    else if (nb.ox1<0) ei++;
  }
  if (nb.polar) {
    Real sign = flip_across_pole_field[IB1] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma omp simd linear(p)
        for (int i=si; i<=ei; ++i)
          dst.x1f(k,j,i) = sign*buf[p++];
      }
    }
  } else {
    BufferUtility::Unpack3DData(buf, dst.x1f, si, ei, sj, ej, sk, ek, p);
  }

  // bx2
  if (nb.ox1==0)      si=pmb->is,         ei=pmb->ie;
  else if (nb.ox1>0)  si=pmb->ie+1,       ei=pmb->ie+NGHOST;
  else               si=pmb->is-NGHOST,  ei=pmb->is-1;
  if (pmb->block_size.nx2==1) sj=pmb->js, ej=pmb->je;
  else if (nb.ox2==0) sj=pmb->js,         ej=pmb->je+1;
  else if (nb.ox2>0)  sj=pmb->je+2,       ej=pmb->je+NGHOST+1;
  else               sj=pmb->js-NGHOST,  ej=pmb->js-1;
  // for SMR/AMR, always include the overlapping faces in edge and corner boundaries
  if (pmb->pmy_mesh->multilevel==true && nb.type != NEIGHBOR_FACE) {
    if (nb.ox2>0) sj--;
    else if (nb.ox2<0) ej++;
  }
  if (nb.polar) {
    Real sign = flip_across_pole_field[IB2] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma omp simd linear(p)
        for (int i=si; i<=ei; ++i)
          dst.x2f(k,j,i) = sign*buf[p++];
      }
    }
  } else {
    BufferUtility::Unpack3DData(buf, dst.x2f, si, ei, sj, ej, sk, ek, p);
  }
  if (pmb->block_size.nx2==1) { // 1D
#pragma omp simd
    for (int i=si; i<=ei; ++i)
      dst.x2f(sk,sj+1,i) = dst.x2f(sk,sj,i);
  }

  // bx3
  if (nb.ox2==0)      sj=pmb->js,         ej=pmb->je;
  else if (nb.ox2>0)  sj=pmb->je+1,       ej=pmb->je+NGHOST;
  else               sj=pmb->js-NGHOST,  ej=pmb->js-1;
  if (pmb->block_size.nx3==1) sk=pmb->ks, ek=pmb->ke;
  else if (nb.ox3==0) sk=pmb->ks,         ek=pmb->ke+1;
  else if (nb.ox3>0)  sk=pmb->ke+2,       ek=pmb->ke+NGHOST+1;
  else               sk=pmb->ks-NGHOST,  ek=pmb->ks-1;
  // for SMR/AMR, always include the overlapping faces in edge and corner boundaries
  if (pmb->pmy_mesh->multilevel==true && nb.type != NEIGHBOR_FACE) {
    if (nb.ox3>0) sk--;
    else if (nb.ox3<0) ek++;
  }
  if (nb.polar) {
    Real sign = flip_across_pole_field[IB3] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma omp simd linear(p)
        for (int i=si; i<=ei; ++i)
          dst.x3f(k,j,i) = sign*buf[p++];
      }
    }
  } else {
    BufferUtility::Unpack3DData(buf, dst.x3f, si, ei, sj, ej, sk, ek, p);
  }
  if (pmb->block_size.nx3==1) { // 1D or 2D
    for (int j=sj; j<=ej; ++j) {
#pragma omp simd
      for (int i=si; i<=ei; ++i)
        dst.x3f(sk+1,j,i) = dst.x3f(sk,j,i);
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetFieldBoundaryFromCoarser(Real *buf,
//                                                       const NeighborBlock& nb)
//  \brief Set field prolongation buffer received from a block on the same level

void BoundaryValues::SetFieldBoundaryFromCoarser(Real *buf, const NeighborBlock& nb) {
  MeshBlock *pmb=pmy_block_;
  MeshRefinement *pmr=pmb->pmr;
  int si, sj, sk, ei, ej, ek;
  int cng=pmb->cnghost;
  int p=0;

  // bx1
  if (nb.ox1==0) {
    si=pmb->cis, ei=pmb->cie+1;
    if ((pmb->loc.lx1&1L)==0L) ei+=cng;
    else             si-=cng;
  } else if (nb.ox1>0) {  si=pmb->cie+1,   ei=pmb->cie+1+cng;}
  else               si=pmb->cis-cng, ei=pmb->cis;
  if (nb.ox2==0) {
    sj=pmb->cjs, ej=pmb->cje;
    if (pmb->block_size.nx2 > 1) {
      if ((pmb->loc.lx2&1L)==0L) ej+=cng;
      else             sj-=cng;
    }
  } else if (nb.ox2>0) {  sj=pmb->cje+1,   ej=pmb->cje+cng;}
  else               sj=pmb->cjs-cng, ej=pmb->cjs-1;
  if (nb.ox3==0) {
    sk=pmb->cks, ek=pmb->cke;
    if (pmb->block_size.nx3 > 1) {
      if ((pmb->loc.lx3&1L)==0L) ek+=cng;
      else             sk-=cng;
    }
  } else if (nb.ox3>0) {  sk=pmb->cke+1,   ek=pmb->cke+cng;}
  else               sk=pmb->cks-cng, ek=pmb->cks-1;

  if (nb.polar) {
    Real sign = flip_across_pole_field[IB1] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma omp simd linear(p)
        for (int i=si; i<=ei; ++i)
          pmr->coarse_b_.x1f(k,j,i) = sign*buf[p++];
      }
    }
  } else {
    BufferUtility::Unpack3DData(buf, pmr->coarse_b_.x1f, si, ei, sj, ej, sk, ek, p);
  }

  // bx2
  if (nb.ox1==0) {
    si=pmb->cis, ei=pmb->cie;
    if ((pmb->loc.lx1&1L)==0L) ei+=cng;
    else             si-=cng;
  } else if (nb.ox1>0) {  si=pmb->cie+1,   ei=pmb->cie+cng;}
  else               si=pmb->cis-cng, ei=pmb->cis-1;
  if (nb.ox2==0) {
    sj=pmb->cjs, ej=pmb->cje;
    if (pmb->block_size.nx2 > 1) {
      ej++;
      if ((pmb->loc.lx2&1L)==0L) ej+=cng;
      else             sj-=cng;
    }
  } else if (nb.ox2>0) {  sj=pmb->cje+1,   ej=pmb->cje+1+cng;}
  else               sj=pmb->cjs-cng, ej=pmb->cjs;

  if (nb.polar) {
    Real sign = flip_across_pole_field[IB2] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma omp simd linear(p)
        for (int i=si; i<=ei; ++i)
          pmr->coarse_b_.x2f(k,j,i) = sign*buf[p++];
      }
    }
  } else {
    BufferUtility::Unpack3DData(buf, pmr->coarse_b_.x2f, si, ei, sj, ej, sk, ek, p);
    if (pmb->block_size.nx2 == 1) { // 1D
#pragma omp simd
      for (int i=si; i<=ei; ++i)
        pmr->coarse_b_.x2f(sk,sj+1,i) = pmr->coarse_b_.x2f(sk,sj,i);
    }
  }

  // bx3
  if (nb.ox2==0) {
    sj=pmb->cjs, ej=pmb->cje;
    if (pmb->block_size.nx2 > 1) {
      if ((pmb->loc.lx2&1L)==0L) ej+=cng;
      else             sj-=cng;
    }
  } else if (nb.ox2>0) {  sj=pmb->cje+1,   ej=pmb->cje+cng;}
  else               sj=pmb->cjs-cng, ej=pmb->cjs-1;
  if (nb.ox3==0) {
    sk=pmb->cks, ek=pmb->cke;
    if (pmb->block_size.nx3 > 1) {
      ek++;
      if ((pmb->loc.lx3&1L)==0L) ek+=cng;
      else             sk-=cng;
    }
  } else if (nb.ox3>0) {  sk=pmb->cke+1,   ek=pmb->cke+1+cng;}
  else               sk=pmb->cks-cng, ek=pmb->cks;

  if (nb.polar) {
    Real sign = flip_across_pole_field[IB3] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma omp simd linear(p)
        for (int i=si; i<=ei; ++i)
          pmr->coarse_b_.x3f(k,j,i) = sign*buf[p++];
      }
    }
  } else {
    BufferUtility::Unpack3DData(buf, pmr->coarse_b_.x3f, si, ei, sj, ej, sk, ek, p);
    if (pmb->block_size.nx3 == 1) { // 2D
      for (int j=sj; j<=ej; ++j) {
        for (int i=si; i<=ei; ++i)
          pmr->coarse_b_.x3f(sk+1,j,i) = pmr->coarse_b_.x3f(sk,j,i);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetFielBoundaryFromFiner(FaceField &dst,
//                                                    Real *buf, const NeighborBlock& nb)
//  \brief Set field boundary received from a block on the same level

void BoundaryValues::SetFieldBoundaryFromFiner(FaceField &dst, Real *buf,
                                               const NeighborBlock& nb) {
  MeshBlock *pmb=pmy_block_;
  // receive already restricted data
  int si, sj, sk, ei, ej, ek;
  int p=0;

  // bx1
  if (nb.ox1==0) {
    si=pmb->is, ei=pmb->ie+1;
    if (nb.fi1==1)   si+=pmb->block_size.nx1/2;
    else            ei-=pmb->block_size.nx1/2;
  } else if (nb.ox1>0) { si=pmb->ie+2,      ei=pmb->ie+NGHOST+1;}
  else              si=pmb->is-NGHOST, ei=pmb->is-1;
  // include the overlapping faces in edge and corner boundaries
  if (nb.type != NEIGHBOR_FACE) {
    if (nb.ox1>0) si--;
    else if (nb.ox1<0) ei++;
  }
  if (nb.ox2==0) {
    sj=pmb->js, ej=pmb->je;
    if (pmb->block_size.nx2 > 1) {
      if (nb.ox1!=0) {
        if (nb.fi1==1) sj+=pmb->block_size.nx2/2;
        else          ej-=pmb->block_size.nx2/2;
      } else {
        if (nb.fi2==1) sj+=pmb->block_size.nx2/2;
        else          ej-=pmb->block_size.nx2/2;
      }
    }
  } else if (nb.ox2>0) { sj=pmb->je+1,      ej=pmb->je+NGHOST;}
  else              sj=pmb->js-NGHOST, ej=pmb->js-1;
  if (nb.ox3==0) {
    sk=pmb->ks, ek=pmb->ke;
    if (pmb->block_size.nx3 > 1) {
      if (nb.ox1!=0 && nb.ox2!=0) {
        if (nb.fi1==1) sk+=pmb->block_size.nx3/2;
        else          ek-=pmb->block_size.nx3/2;
      } else {
        if (nb.fi2==1) sk+=pmb->block_size.nx3/2;
        else          ek-=pmb->block_size.nx3/2;
      }
    }
  } else if (nb.ox3>0) { sk=pmb->ke+1,      ek=pmb->ke+NGHOST;}
  else              sk=pmb->ks-NGHOST, ek=pmb->ks-1;

  if (nb.polar) {
    Real sign = flip_across_pole_field[IB1] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma omp simd linear(p)
        for (int i=si; i<=ei; ++i)
          dst.x1f(k,j,i) = sign*buf[p++];
      }
    }
  } else {
    BufferUtility::Unpack3DData(buf, dst.x1f, si, ei, sj, ej, sk, ek, p);
  }

  // bx2
  if (nb.ox1==0) {
    si=pmb->is, ei=pmb->ie;
    if (nb.fi1==1)   si+=pmb->block_size.nx1/2;
    else            ei-=pmb->block_size.nx1/2;
  } else if (nb.ox1>0) { si=pmb->ie+1,      ei=pmb->ie+NGHOST;}
  else              si=pmb->is-NGHOST, ei=pmb->is-1;
  if (nb.ox2==0) {
    sj=pmb->js, ej=pmb->je;
    if (pmb->block_size.nx2 > 1) {
      ej++;
      if (nb.ox1!=0) {
        if (nb.fi1==1) sj+=pmb->block_size.nx2/2;
        else          ej-=pmb->block_size.nx2/2;
      } else {
        if (nb.fi2==1) sj+=pmb->block_size.nx2/2;
        else          ej-=pmb->block_size.nx2/2;
      }
    }
  } else if (nb.ox2>0) { sj=pmb->je+2,      ej=pmb->je+NGHOST+1;}
  else              sj=pmb->js-NGHOST, ej=pmb->js-1;
  // include the overlapping faces in edge and corner boundaries
  if (nb.type != NEIGHBOR_FACE) {
    if (nb.ox2>0) sj--;
    else if (nb.ox2<0) ej++;
  }

  if (nb.polar) {
    Real sign = flip_across_pole_field[IB2] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma omp simd linear(p)
        for (int i=si; i<=ei; ++i)
          dst.x2f(k,j,i) = sign*buf[p++];
      }
    }
  } else {
    BufferUtility::Unpack3DData(buf, dst.x2f, si, ei, sj, ej, sk, ek, p);
  }
  if (pmb->block_size.nx2==1) { // 1D
#pragma omp simd
    for (int i=si; i<=ei; ++i)
      dst.x2f(sk,sj+1,i) = dst.x2f(sk,sj,i);
  }

  // bx3
  if (nb.ox2==0) {
    sj=pmb->js, ej=pmb->je;
    if (pmb->block_size.nx2 > 1) {
      if (nb.ox1!=0) {
        if (nb.fi1==1) sj+=pmb->block_size.nx2/2;
        else          ej-=pmb->block_size.nx2/2;
      } else {
        if (nb.fi2==1) sj+=pmb->block_size.nx2/2;
        else          ej-=pmb->block_size.nx2/2;
      }
    }
  } else if (nb.ox2>0) { sj=pmb->je+1,      ej=pmb->je+NGHOST;}
  else              sj=pmb->js-NGHOST, ej=pmb->js-1;
  if (nb.ox3==0) {
    sk=pmb->ks, ek=pmb->ke;
    if (pmb->block_size.nx3 > 1) {
      ek++;
      if (nb.ox1!=0 && nb.ox2!=0) {
        if (nb.fi1==1) sk+=pmb->block_size.nx3/2;
        else          ek-=pmb->block_size.nx3/2;
      } else {
        if (nb.fi2==1) sk+=pmb->block_size.nx3/2;
        else          ek-=pmb->block_size.nx3/2;
      }
    }
  } else if (nb.ox3>0) { sk=pmb->ke+2,      ek=pmb->ke+NGHOST+1;}
  else              sk=pmb->ks-NGHOST, ek=pmb->ks-1;
  // include the overlapping faces in edge and corner boundaries
  if (nb.type != NEIGHBOR_FACE) {
    if (nb.ox3>0) sk--;
    else if (nb.ox3<0) ek++;
  }

  if (nb.polar) {
    Real sign = flip_across_pole_field[IB3] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma omp simd linear(p)
        for (int i=si; i<=ei; ++i)
          dst.x3f(k,j,i) = sign*buf[p++];
      }
    }
  } else {
    BufferUtility::Unpack3DData(buf, dst.x3f, si, ei, sj, ej, sk, ek, p);
  }
  if (pmb->block_size.nx3==1) { // 1D or 2D
    for (int j=sj; j<=ej; ++j) {
#pragma omp simd
      for (int i=si; i<=ei; ++i)
        dst.x3f(sk+1,j,i) = dst.x3f(sk,j,i);
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn bool BoundaryValues::ReceiveFieldBoundaryBuffers(FaceField &dst)
//  \brief load boundary buffer for x1 direction into the array

bool BoundaryValues::ReceiveFieldBoundaryBuffers(FaceField &dst) {
  MeshBlock *pmb=pmy_block_;
  bool flag=true;

  for (int n=0; n<nneighbor; n++) {
    NeighborBlock& nb = neighbor[n];
    if (bd_field_.flag[nb.bufid]==BNDRY_COMPLETED) continue;
    if (bd_field_.flag[nb.bufid]==BNDRY_WAITING) {
      if (nb.rank==Globals::my_rank) {// on the same process
        flag=false;
        continue;
#ifdef MPI_PARALLEL
      } else { // MPI boundary
        int test;
        MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&test,MPI_STATUS_IGNORE);
        MPI_Test(&(bd_field_.req_recv[nb.bufid]),&test,MPI_STATUS_IGNORE);
        if (static_cast<bool>(test)==false) {
          flag=false;
          continue;
        }
        bd_field_.flag[nb.bufid] = BNDRY_ARRIVED;
      }
#else
      }
#endif
    }
    if (nb.level==pmb->loc.level)
      SetFieldBoundarySameLevel(dst, bd_field_.recv[nb.bufid], nb);
    else if (nb.level<pmb->loc.level)
      SetFieldBoundaryFromCoarser(bd_field_.recv[nb.bufid], nb);
    else
      SetFieldBoundaryFromFiner(dst, bd_field_.recv[nb.bufid], nb);
    bd_field_.flag[nb.bufid] = BNDRY_COMPLETED; // completed
  }

  if (flag
      and (block_bcs[INNER_X2] == POLAR_BNDRY or block_bcs[OUTER_X2] == POLAR_BNDRY)) {
    PolarSingleField(dst);
    PolarAxisFieldAverage(dst);
  }

  return flag;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ReceiveFieldBoundaryBuffersWithWait(FaceField &dst)
//  \brief load boundary buffer for x1 direction into the array

void BoundaryValues::ReceiveFieldBoundaryBuffersWithWait(FaceField &dst) {
  MeshBlock *pmb=pmy_block_;

  for (int n=0; n<nneighbor; n++) {
    NeighborBlock& nb = neighbor[n];
#ifdef MPI_PARALLEL
    if (nb.rank!=Globals::my_rank)
      MPI_Wait(&(bd_field_.req_recv[nb.bufid]),MPI_STATUS_IGNORE);
#endif
    if (nb.level==pmb->loc.level)
      SetFieldBoundarySameLevel(dst, bd_field_.recv[nb.bufid], nb);
    else if (nb.level<pmb->loc.level)
      SetFieldBoundaryFromCoarser(bd_field_.recv[nb.bufid], nb);
    else
      SetFieldBoundaryFromFiner(dst, bd_field_.recv[nb.bufid], nb);
    bd_field_.flag[nb.bufid] = BNDRY_COMPLETED; // completed
  }

  if (block_bcs[INNER_X2] == POLAR_BNDRY or block_bcs[OUTER_X2] == POLAR_BNDRY) {
    PolarSingleField(dst);
    PolarAxisFieldAverage(dst);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::PolarSingleField(FaceField &dst)
//
//  \brief single block in the azimuthal direction for the polar boundary

void BoundaryValues::PolarSingleField(FaceField &dst) {
  MeshBlock *pmb=pmy_block_;
  if (pmb->loc.level == pmb->pmy_mesh->root_level && pmb->pmy_mesh->nrbx3 == 1
  && pmb->block_size.nx3 > 1) {
    if (block_bcs[INNER_X2]==POLAR_BNDRY) {
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      for (int j=pmb->js-NGHOST; j<=pmb->js-1; ++j) {
       for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST+1; ++i) {
         for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k)
           exc_(k)=dst.x1f(k,j,i);
         for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
           int k_shift = k;
           k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
           dst.x1f(k,j,i)=exc_(k_shift);
         }
       }
      }
      for (int j=pmb->js-NGHOST; j<=pmb->js-1; ++j) {
       for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i) {
         for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k)
           exc_(k)=dst.x2f(k,j,i);
         for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
           int k_shift = k;
           k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
           dst.x2f(k,j,i)=exc_(k_shift);
         }
       }
      }
      for (int j=pmb->js-NGHOST; j<=pmb->js-1; ++j) {
       for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i) {
         for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST+1; ++k)
           exc_(k)=dst.x3f(k,j,i);
         for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST+1; ++k) {
           int k_shift = k;
           k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
           dst.x3f(k,j,i)=exc_(k_shift);
         }
       }
      }
    }

    if (block_bcs[OUTER_X2]==POLAR_BNDRY) {
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      for (int j=pmb->je+1; j<=pmb->je+NGHOST; ++j) {
        for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST+1; ++i) {
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k)
            exc_(k)=dst.x1f(k,j,i);
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
            int k_shift = k;
            k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
            dst.x1f(k,j,i)=exc_(k_shift);
          }
        }
      }
      for (int j=pmb->je+2; j<=pmb->je+NGHOST+1; ++j) {
        for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i) {
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k)
            exc_(k)=dst.x2f(k,j,i);
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
            int k_shift = k;
            k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
            dst.x2f(k,j,i)=exc_(k_shift);
          }
        }
      }
      for (int j=pmb->je+1; j<=pmb->je+NGHOST; ++j) {
        for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i) {
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST+1; ++k)
            exc_(k)=dst.x3f(k,j,i);
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST+1; ++k) {
            int k_shift = k;
            k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
            dst.x3f(k,j,i)=exc_(k_shift);
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::PolarAxisFieldAverage(FaceField &dst)
//
//  \brief set theta-component of field along axis

void BoundaryValues::PolarAxisFieldAverage(FaceField &dst) {
  MeshBlock *pmb = pmy_block_;
  int is = pmb->is - NGHOST;
  int ie = pmb->ie + NGHOST;
  int ks = pmb->ks;
  int ke = pmb->ke;
  if (pmb->block_size.nx3 > 1) {
    ks -= NGHOST;
    ke += NGHOST;
  }
  if (block_bcs[INNER_X2] == POLAR_BNDRY) {
    int j = pmb->js;
    for (int k = ks; k <= ke; ++k) {
      for (int i = is; i <= ie; ++i) {
        dst.x2f(k,j,i) = 0.5 * (dst.x2f(k,j-1,i) + dst.x2f(k,j+1,i));
      }
    }
  }
  if (block_bcs[OUTER_X2] == POLAR_BNDRY) {
    int j = pmb->je + 1;
    for (int k = ks; k <= ke; ++k) {
      for (int i = is; i <= ie; ++i) {
        dst.x2f(k,j,i) = 0.5 * (dst.x2f(k,j-1,i) + dst.x2f(k,j+1,i));
      }
    }
  }
  return;
}
