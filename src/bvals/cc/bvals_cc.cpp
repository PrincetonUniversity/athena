//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_cc.cpp
//  \brief functions that apply BCs for CELL_CENTERED variables

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>    // memcpy()
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
#include "bvals_cc.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// constructor

CellCenteredBoundaryVariable::CellCenteredBoundaryVariable(
    MeshBlock *pmb, AthenaArray<Real> *var, AthenaArray<Real> *var_flux)
    : BoundaryVariable(pmb) {

  // var_cc.InitWithShallowCopy(var);
  var_cc = var;
  // src.InitWithShallowCopy(var_cc);
  // dst.InitWithShallowCopy(var_cc);

  // KGF: uninitialized, for now. Relying on HydroBoundaryVariable:SelectCoarseBuffer() to
  // set it to pmr-> cons or prim, when needed. Must set it directly for other
  // CellCenteredBoundaryVarialbe objects
  // coarse_buf.InitWithShallowCopy(var_cc);

  if (var_flux) { // allow nullptr to be passed for non-AMR/SMR variables w/o seg-faulting
    x1flux.InitWithShallowCopy(var_flux[X1DIR]);
    x2flux.InitWithShallowCopy(var_flux[X2DIR]);
    x3flux.InitWithShallowCopy(var_flux[X3DIR]);
  }
  // KGF: Taken from 2x functions in flux_correction_cc.cpp
  // if (type == FluxCorrectionQuantity::hydro) {
  //   nl_=0, nu_=NHYDRO-1;
  //   x1flux.InitWithShallowCopy(pmb->phydro->flux[X1DIR]);
  //   x2flux.InitWithShallowCopy(pmb->phydro->flux[X2DIR]);
  //   x3flux.InitWithShallowCopy(pmb->phydro->flux[X3DIR]);
  // }

  // KGF: CellCenteredBoundaryVariable should only be used w/ 4D or 3D (nx4=1) AthenaArray
  // For now, assume that full span of 4th dim of input AthenaArray should be used;
  // get the index limits directly from the input AthenaArray

  // Loops over generic cc arrays are inclusive of upper limit in 4D, e.g. nl_<=n<=nu_
  // Outside of bvals/, loops over the 4th dimension are typically exclusive:
  // e.g. in calculate_fluxes.cpp, loops are n<NWAVE or n<NHYDRO
  nl_ = 0;
  nu_ = var->GetDim4() - 1;

  // class data members that are initialized only in derived class HydroBoundaryVariable()
  //coarse_buf=;  // unitialized AthenaArray
  flip_across_pole_ = nullptr;

  InitBoundaryData(bd_var_, BoundaryQuantity::cc);
#ifdef MPI_PARALLEL
  cc_phys_id_ = pbval_->ReserveTagVariableIDs(1);
#endif
  if (pmy_mesh_->multilevel == true) { // SMR or AMR
    InitBoundaryData(bd_var_flcor_, BoundaryQuantity::cc_flcor);
    // KGF: if storing pointer instead of BoundaryData in BoundaryVariable base class:
    // pbd_var_flcor_ = &(bd_var_flcor_);
#ifdef MPI_PARALLEL
    cc_flx_phys_id_ = pbval_->ReserveTagVariableIDs(1);
#endif
  }
}

// destructor

CellCenteredBoundaryVariable::~CellCenteredBoundaryVariable() {
  // KGF: similar to section in constructor, this can be automatically handled in separate
  // classes
  DestroyBoundaryData(bd_var_);
  if (pmy_mesh_->multilevel == true) // SMR or AMR
    DestroyBoundaryData(bd_var_flcor_);
}

int CellCenteredBoundaryVariable::ComputeVariableBufferSize(const NeighborIndexes& ni,
                                                            int cng) {
  MeshBlock *pmb = pmy_block_;
  int cng1, cng2, cng3;
  cng1 = cng;
  cng2 = cng*(pmb->block_size.nx2 > 1 ? 1 : 0);
  cng3 = cng*(pmb->block_size.nx3 > 1 ? 1 : 0);

  int size=((ni.ox1 == 0)?pmb->block_size.nx1:NGHOST)
           *((ni.ox2 == 0)?pmb->block_size.nx2:NGHOST)
           *((ni.ox3 == 0)?pmb->block_size.nx3:NGHOST);
  if (pmy_mesh_->multilevel) {
    int f2c=((ni.ox1 == 0) ? ((pmb->block_size.nx1+1)/2):NGHOST)
            *((ni.ox2 == 0) ? ((pmb->block_size.nx2+1)/2):NGHOST)
            *((ni.ox3 == 0) ? ((pmb->block_size.nx3+1)/2):NGHOST);
    int c2f=((ni.ox1 == 0) ?((pmb->block_size.nx1+1)/2+cng1):cng)
            *((ni.ox2 == 0)?((pmb->block_size.nx2+1)/2+cng2):cng)
            *((ni.ox3 == 0)?((pmb->block_size.nx3+1)/2+cng3):cng);
    size = std::max(size,c2f);
    size = std::max(size,f2c);
  }
  size *= nu_+1;
  return size;
}

int CellCenteredBoundaryVariable::ComputeFluxCorrectionBufferSize(
    const NeighborIndexes& ni, int cng) {
  MeshBlock *pmb = pmy_block_;
  int size = 0;
  if (ni.ox1 != 0)
    size = (pmb->block_size.nx2+1)/2*(pmb->block_size.nx3+1)/2*(nu_+1);
  if (ni.ox2 != 0)
    size = (pmb->block_size.nx1+1)/2*(pmb->block_size.nx3+1)/2*(nu_+1);
  if (ni.ox3 != 0)
    size = (pmb->block_size.nx1+1)/2*(pmb->block_size.nx2+1)/2*(nu_+1);
  return size;
}

//----------------------------------------------------------------------------------------
//! \fn int CellCenteredBoundaryVariable::LoadBoundaryBufferSameLevel(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set cell-centered boundary buffers for sending to a block on the same level

int CellCenteredBoundaryVariable::LoadBoundaryBufferSameLevel(Real *buf,
                                                              const NeighborBlock& nb) {
  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;

  si = (nb.ox1>0)?(pmb->ie-NGHOST+1):pmb->is;
  ei = (nb.ox1<0)?(pmb->is+NGHOST-1):pmb->ie;
  sj = (nb.ox2>0)?(pmb->je-NGHOST+1):pmb->js;
  ej = (nb.ox2<0)?(pmb->js+NGHOST-1):pmb->je;
  sk = (nb.ox3>0)?(pmb->ke-NGHOST+1):pmb->ks;
  ek = (nb.ox3<0)?(pmb->ks+NGHOST-1):pmb->ke;
  int p = 0;
  // KGF: 3/21/19 workarounds for handling mutable pdata ptr in AthenaArray obj. In either
  // approach, the pointer is resolved as soon as this function is called, so it should
  // contain the updated AthenaArray member pdata, unless a SwapAthenaArray() call occurs
  // while the function is running. In either Approach #1 or #2, this window of danger is
  // only if the call occurs inside Pack4DData(), since dereferencing the ptr to
  // AthenaArray when initializing the reference in Approach #1 does not imply a local,
  // independent copy of pdata, unlike the former "InitWithShallowCopy()" approach.

  // KGF: approach #1: store local reference variable to dereferenced pointer in each fn.
  AthenaArray<Real> &var = *var_cc; // KGF: need to rename from "var_cc" member....
  BufferUtility::PackData(var, buf, nl_, nu_, si, ei, sj, ej, sk, ek, p);
  // KGF: approach #2: just dereference original pointer each time it is used
  // BufferUtility::Pack4DData(*var_cc, buf, nl_, nu_, si, ei, sj, ej, sk, ek, p);
  return p;
}

//----------------------------------------------------------------------------------------
//! \fn int CellCenteredBoundaryVariable::LoadBoundaryBufferToCoarser(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set cell-centered boundary buffers for sending to a block on the coarser level

int CellCenteredBoundaryVariable::LoadBoundaryBufferToCoarser(Real *buf,
                                                              const NeighborBlock& nb) {
  MeshBlock *pmb = pmy_block_;
  MeshRefinement *pmr = pmb->pmr;
  int si, sj, sk, ei, ej, ek;
  int cn = NGHOST-1;
  AthenaArray<Real> &var = *var_cc;

  si = (nb.ox1>0)?(pmb->cie-cn):pmb->cis;
  ei = (nb.ox1<0)?(pmb->cis+cn):pmb->cie;
  sj = (nb.ox2>0)?(pmb->cje-cn):pmb->cjs;
  ej = (nb.ox2<0)?(pmb->cjs+cn):pmb->cje;
  sk = (nb.ox3>0)?(pmb->cke-cn):pmb->cks;
  ek = (nb.ox3<0)?(pmb->cks+cn):pmb->cke;

  int p = 0;
  pmr->RestrictCellCenteredValues(var, coarse_buf, nl_, nu_, si, ei, sj, ej, sk, ek);
  BufferUtility::PackData(coarse_buf, buf, nl_, nu_, si, ei, sj, ej, sk, ek, p);
  return p;
}

//----------------------------------------------------------------------------------------
//! \fn int CellCenteredBoundaryVariable::LoadBoundaryBufferToFiner(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set cell-centered boundary buffers for sending to a block on the finer level

int CellCenteredBoundaryVariable::LoadBoundaryBufferToFiner(Real *buf,
                                                            const NeighborBlock& nb) {
  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;
  int cn = pmb->cnghost-1;
  AthenaArray<Real> &var = *var_cc;

  si = (nb.ox1>0)?(pmb->ie-cn):pmb->is;
  ei = (nb.ox1<0)?(pmb->is+cn):pmb->ie;
  sj = (nb.ox2>0)?(pmb->je-cn):pmb->js;
  ej = (nb.ox2<0)?(pmb->js+cn):pmb->je;
  sk = (nb.ox3>0)?(pmb->ke-cn):pmb->ks;
  ek = (nb.ox3<0)?(pmb->ks+cn):pmb->ke;

  // send the data first and later prolongate on the target block
  // need to add edges for faces, add corners for edges
  if (nb.ox1 == 0) {
    if (nb.fi1 == 1)   si += pmb->block_size.nx1/2-pmb->cnghost;
    else            ei -= pmb->block_size.nx1/2-pmb->cnghost;
  }
  if (nb.ox2 == 0 && pmb->block_size.nx2 > 1) {
    if (nb.ox1 != 0) {
      if (nb.fi1 == 1) sj += pmb->block_size.nx2/2-pmb->cnghost;
      else          ej -= pmb->block_size.nx2/2-pmb->cnghost;
    } else {
      if (nb.fi2 == 1) sj += pmb->block_size.nx2/2-pmb->cnghost;
      else          ej -= pmb->block_size.nx2/2-pmb->cnghost;
    }
  }
  if (nb.ox3 == 0 && pmb->block_size.nx3 > 1) {
    if (nb.ox1 != 0 && nb.ox2 != 0) {
      if (nb.fi1 == 1) sk += pmb->block_size.nx3/2-pmb->cnghost;
      else          ek -= pmb->block_size.nx3/2-pmb->cnghost;
    } else {
      if (nb.fi2 == 1) sk += pmb->block_size.nx3/2-pmb->cnghost;
      else          ek -= pmb->block_size.nx3/2-pmb->cnghost;
    }
  }

  int p=0;
  BufferUtility::PackData(var, buf, nl_, nu_, si, ei, sj, ej, sk, ek, p);
  return p;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::SendBoundaryBuffers()
//  \brief Send boundary buffers of cell-centered variables

// KGF: (AthenaArray<Real> &dst, HydroBoundaryQuantity type)
void CellCenteredBoundaryVariable::SendBoundaryBuffers() {
  MeshBlock *pmb = pmy_block_;
  int mylevel = pmb->loc.level;

  // KGF: call switch over "HydroBoundaryQuantity type"

  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    int ssize;
    if (nb.level == mylevel)
      // KGF: src/var_cc
      // KGF: nl_, nu_
      ssize=LoadBoundaryBufferSameLevel(bd_var_.send[nb.bufid], nb);
    else if (nb.level<mylevel)
      // KGF: src/var_cc
      // KGF: nl_, nu_
      // KGF: coarse_buf
      ssize=LoadBoundaryBufferToCoarser(bd_var_.send[nb.bufid], nb);
    else
      // KGF: src/var_cc
      // KGF: nl_, nu_
      ssize=LoadBoundaryBufferToFiner(bd_var_.send[nb.bufid], nb);
    if (nb.rank == Globals::my_rank) {  // on the same process
      // KGF: additional "HydroBoundaryBuffer type" switch unique to
      // SendBoundaryBuffers().
      // if (type == HydroBoundaryQuantity::cons || type == HydroBoundaryQuantity::prim)
      //   ptarget=&(ptarget_block->pbval_->bd_var_);
      CopyVariableBufferSameProcess(nb, ssize);
    }
#ifdef MPI_PARALLEL
    else // MPI
      MPI_Start(&(bd_var_.req_send[nb.bufid]));
#endif
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::SetBoundarySameLevel(Real *buf,
//                                                              const NeighborBlock& nb)
//  \brief Set cell-centered boundary received from a block on the same level

void CellCenteredBoundaryVariable::SetBoundarySameLevel(Real *buf,
                                                        const NeighborBlock& nb) {
  MeshBlock *pmb=pmy_block_;
  int si, sj, sk, ei, ej, ek;
  AthenaArray<Real> &var = *var_cc;

  if (nb.ox1 == 0)     si=pmb->is,        ei=pmb->ie;
  else if (nb.ox1>0) si=pmb->ie+1,      ei=pmb->ie+NGHOST;
  else              si=pmb->is-NGHOST, ei=pmb->is-1;
  if (nb.ox2 == 0)     sj=pmb->js,        ej=pmb->je;
  else if (nb.ox2>0) sj=pmb->je+1,      ej=pmb->je+NGHOST;
  else              sj=pmb->js-NGHOST, ej=pmb->js-1;
  if (nb.ox3 == 0)     sk=pmb->ks,        ek=pmb->ke;
  else if (nb.ox3>0) sk=pmb->ke+1,      ek=pmb->ke+NGHOST;
  else              sk=pmb->ks-NGHOST, ek=pmb->ks-1;

  int p=0;

  if (nb.polar) {
    for (int n=nl_; n<=nu_; ++n) {
      Real sign = 1.0;
      if (flip_across_pole_ != nullptr) sign = flip_across_pole_[n] ? -1.0 : 1.0;
      for (int k=sk; k<=ek; ++k) {
        for (int j=ej; j>=sj; --j) {
#pragma omp simd
          for (int i=si; i<=ei; ++i) {
            // KGF: approach #2
            //(*var_cc)(n,k,j,i) = sign * buf[p++];
            var(n,k,j,i) = sign * buf[p++];
          }
        }
      }
    }
  } else {
    BufferUtility::UnpackData(buf, var, nl_, nu_, si, ei, sj, ej, sk, ek, p);
  }
  // 2d shearingbox in x-z plane: additional step to shift azimuthal velocity;
  // if (SHEARING_BOX) {
  //   if (ShBoxCoord_ == 2) {
  //     int level = pmb->loc.level - pmy_mesh_->root_level;
  //     std::int64_t nrbx1 = pmy_mesh_->nrbx1*(1L << level);
  //     Real qomL = qshear_*Omega_0_*x1size_;
  //     if ((pmb->loc.lx1 == 0) && (nb.ox1<0)) {
  //       for (int k=sk; k<=ek; ++k) {
  //         for (int j=sj; j<=ej; ++j) {
  //           for (int i=si; i<=ei; ++i) {
  //             if (NON_BAROTROPIC_EOS)
  //               dst(IEN,k,j,i) += (0.5/dst(IDN,k,j,i))
  //                                 *(SQR(dst(IM3,k,j,i)+qomL*dst(IDN,k,j,i))
  //                                   -SQR(dst(IM3,k,j,i)));
  //             dst(IM3,k,j,i) += qomL*dst(IDN,k,j,i);
  //           }
  //         }
  //       }
  //     } // inner boundary
  //     if ((pmb->loc.lx1 == (nrbx1-1)) && (nb.ox1>0)) {
  //       for (int k=sk; k<=ek; ++k) {
  //         for (int j=sj; j<=ej; ++j) {
  //           for (int i=si; i<=ei; ++i) {
  //             if (NON_BAROTROPIC_EOS)
  //               dst(IEN,k,j,i) += (0.5/dst(IDN,k,j,i))
  //                                 *(SQR(dst(IM3,k,j,i)-qomL*dst(IDN,k,j,i))
  //                                   -SQR(dst(IM3,k,j,i)));
  //             dst(IM3,k,j,i) -= qomL*dst(IDN,k,j,i);
  //           }
  //         }
  //       }
  //     } // outer boundary
  //   }
  // } // end KGF: shearing box in SetBoundarySameLevel
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::SetBoundaryFromCoarser(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set cell-centered prolongation buffer received from a block on a coarser level

void CellCenteredBoundaryVariable::SetBoundaryFromCoarser(Real *buf,
                                                          const NeighborBlock& nb) {
  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;
  int cng = pmb->cnghost;

  if (nb.ox1 == 0) {
    si = pmb->cis, ei = pmb->cie;
    if ((pmb->loc.lx1 & 1LL) == 0LL) ei += cng;
    else                             si -= cng;
  } else if (nb.ox1>0)  {
    si = pmb->cie+1,   ei = pmb->cie+cng;
  } else {
    si = pmb->cis-cng, ei = pmb->cis-1;
  }
  if (nb.ox2 == 0) {
    sj = pmb->cjs, ej = pmb->cje;
    if (pmb->block_size.nx2 > 1) {
      if ((pmb->loc.lx2 & 1LL) == 0LL) ej += cng;
      else                             sj -= cng;
    }
  } else if (nb.ox2>0) {
    sj = pmb->cje+1,   ej = pmb->cje+cng;
  } else {
    sj = pmb->cjs-cng, ej = pmb->cjs-1;
  }
  if (nb.ox3 == 0) {
    sk = pmb->cks, ek = pmb->cke;
    if (pmb->block_size.nx3 > 1) {
      if ((pmb->loc.lx3 & 1LL) == 0LL) ek += cng;
      else                             sk -= cng;
    }
  } else if (nb.ox3>0)  {
    sk = pmb->cke+1,   ek = pmb->cke+cng;
  } else {
    sk = pmb->cks-cng, ek = pmb->cks-1;
  }

  int p = 0;
  if (nb.polar) {
    for (int n=nl_; n<=nu_; ++n) {
      Real sign = 1.0;
      if (flip_across_pole_ != nullptr) sign = flip_across_pole_[n] ? -1.0 : 1.0;
      for (int k=sk; k<=ek; ++k) {
        for (int j=ej; j>=sj; --j) {
#pragma omp simd
          for (int i=si; i<=ei; ++i)
            coarse_buf(n,k,j,i) = sign * buf[p++];
        }
      }
    }
  } else {
    BufferUtility::UnpackData(buf, coarse_buf, nl_, nu_, si, ei, sj, ej, sk, ek, p);
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::SetBoundaryFromFiner(Real *buf,
//                                                              const NeighborBlock& nb)
//  \brief Set cell-centered boundary received from a block on a finer level

void CellCenteredBoundaryVariable::SetBoundaryFromFiner(Real *buf,
                                                        const NeighborBlock& nb) {
  MeshBlock *pmb = pmy_block_;
  AthenaArray<Real> &var = *var_cc;
  // receive already restricted data
  int si, sj, sk, ei, ej, ek;

  if (nb.ox1 == 0) {
    si = pmb->is, ei = pmb->ie;
    if (nb.fi1 == 1)   si += pmb->block_size.nx1/2;
    else            ei -= pmb->block_size.nx1/2;
  } else if (nb.ox1>0) {
    si = pmb->ie+1,      ei = pmb->ie+NGHOST;
  } else {
    si = pmb->is-NGHOST, ei = pmb->is-1;
  }
  if (nb.ox2 == 0) {
    sj = pmb->js, ej = pmb->je;
    if (pmb->block_size.nx2 > 1) {
      if (nb.ox1 != 0) {
        if (nb.fi1 == 1) sj += pmb->block_size.nx2/2;
        else          ej -= pmb->block_size.nx2/2;
      } else {
        if (nb.fi2 == 1) sj += pmb->block_size.nx2/2;
        else          ej -= pmb->block_size.nx2/2;
      }
    }
  } else if (nb.ox2>0) {
    sj = pmb->je+1,      ej = pmb->je+NGHOST;
  } else {
    sj = pmb->js-NGHOST, ej = pmb->js-1;
  }
  if (nb.ox3 == 0) {
    sk = pmb->ks, ek = pmb->ke;
    if (pmb->block_size.nx3 > 1) {
      if (nb.ox1 != 0 && nb.ox2 != 0) {
        if (nb.fi1 == 1) sk += pmb->block_size.nx3/2;
        else          ek -= pmb->block_size.nx3/2;
      } else {
        if (nb.fi2 == 1) sk += pmb->block_size.nx3/2;
        else          ek -= pmb->block_size.nx3/2;
      }
    }
  } else if (nb.ox3>0) {
    sk = pmb->ke+1,      ek = pmb->ke+NGHOST;
  } else {
    sk = pmb->ks-NGHOST, ek = pmb->ks-1;
  }

  int p = 0;
  if (nb.polar) {
    for (int n=nl_; n<=nu_; ++n) {
      Real sign=1.0;
      if (flip_across_pole_ != nullptr) sign = flip_across_pole_[n] ? -1.0 : 1.0;
      for (int k=sk; k<=ek; ++k) {
        for (int j=ej; j>=sj; --j) {
#pragma omp simd
          for (int i=si; i<=ei; ++i)
            var(n,k,j,i) = sign * buf[p++];
        }
      }
    }
  } else {
    BufferUtility::UnpackData(buf, var, nl_, nu_, si, ei, sj, ej, sk, ek, p);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn bool CellCenteredBoundaryVariable::ReceiveBoundaryBuffers()
//  \brief receive the cell-centered boundary data

bool CellCenteredBoundaryVariable::ReceiveBoundaryBuffers() {
  bool bflag = true;

  // KGF: call short switch over "HydroBoundaryQuantity type"

  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (bd_var_.flag[nb.bufid] == BoundaryStatus::arrived) continue;
    if (bd_var_.flag[nb.bufid] == BoundaryStatus::waiting) {
      if (nb.rank == Globals::my_rank) {// on the same process
        bflag = false;
        continue;
      }
#ifdef MPI_PARALLEL
      else { // NOLINT // MPI boundary
        int test;
        MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&test,MPI_STATUS_IGNORE);
        MPI_Test(&(bd_var_.req_recv[nb.bufid]),&test,MPI_STATUS_IGNORE);
        if (static_cast<bool>(test) == false) {
          bflag = false;
          continue;
        }
        bd_var_.flag[nb.bufid] = BoundaryStatus::arrived;
      }
#endif
    }
  }
  return bflag;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::SetBoundaries()
//  \brief set the cell-centered boundary data

void CellCenteredBoundaryVariable::SetBoundaries() {
  MeshBlock *pmb = pmy_block_;

  // KGF: call switch over "HydroBoundaryQuantity type"

  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (nb.level == pmb->loc.level)
      // KGF: dst
      // KGF: nl_, nu_
      SetBoundarySameLevel(bd_var_.recv[nb.bufid], nb);
    else if (nb.level < pmb->loc.level) // this set only the prolongation buffer
      // KGF: nl_, nu_
      // KGF: coarse_buf
      SetBoundaryFromCoarser(bd_var_.recv[nb.bufid], nb);
    else
      // KGF: dst
      // KGF: nl_, nu_
      SetBoundaryFromFiner(bd_var_.recv[nb.bufid], nb);
    bd_var_.flag[nb.bufid] = BoundaryStatus::completed; // completed
  }

  if (pbval_->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar ||
      pbval_->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar)
    PolarBoundarySingleAzimuthalBlock();

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::ReceiveAndSetBoundariesWithWait()
//  \brief receive and set the cell-centered boundary data for initialization

void CellCenteredBoundaryVariable::ReceiveAndSetBoundariesWithWait() {
  MeshBlock *pmb=pmy_block_;

  // KGF: call switch over "HydroBoundaryQuantity type"

  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
#ifdef MPI_PARALLEL
    if (nb.rank != Globals::my_rank)
      MPI_Wait(&(bd_var_.req_recv[nb.bufid]),MPI_STATUS_IGNORE);
#endif
    if (nb.level == pmb->loc.level)
      // KGF: dst
      // KGF: nl_, nu_
      SetBoundarySameLevel(bd_var_.recv[nb.bufid], nb);
    else if (nb.level<pmb->loc.level)
      // KGF: nl_, nu_
      // KGF: coarse_buf
      SetBoundaryFromCoarser(bd_var_.recv[nb.bufid], nb);
    else
      // KGF: dst
      // KGF: nl_, nu_
      SetBoundaryFromFiner(bd_var_.recv[nb.bufid], nb);
    bd_var_.flag[nb.bufid] = BoundaryStatus::completed; // completed
  }

  if (pbval_->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar
      || pbval_->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar)
    PolarBoundarySingleAzimuthalBlock();

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::PolarBoundarySingleAzimuthalBlock()
// \brief polar boundary edge-case: single MeshBlock spans the entire azimuthal (x3) range

void CellCenteredBoundaryVariable::PolarBoundarySingleAzimuthalBlock() {
  MeshBlock *pmb = pmy_block_;
  AthenaArray<Real> &var = *var_cc;

  if (pmb->loc.level  ==  pmy_mesh_->root_level && pmy_mesh_->nrbx3 == 1
      && pmb->block_size.nx3 > 1) {
    if (pbval_->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar) {
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      for (int n=nl_; n<=nu_; ++n) {
        for (int j=pmb->js-NGHOST; j<=pmb->js-1; ++j) {
          for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i) {
            for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k)
              pbval_->azimuthal_shift_(k) = var(n,k,j,i);
            for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
              int k_shift = k;
              k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
              var(n,k,j,i) = pbval_->azimuthal_shift_(k_shift);
            }
          }
        }
      }
    }

    if (pbval_->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar) {
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      for (int n=nl_; n<=nu_; ++n) {
        for (int j=pmb->je+1; j<=pmb->je+NGHOST; ++j) {
          for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i) {
            for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k)
              pbval_->azimuthal_shift_(k) = var(n,k,j,i);
            for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
              int k_shift = k;
              k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
              var(n,k,j,i) = pbval_->azimuthal_shift_(k_shift);
            }
          }
        }
      }
    }
  }
  return;
}

// ------- KGF: move to a separate file?

void CellCenteredBoundaryVariable::SetupPersistentMPI() {
#ifdef MPI_PARALLEL
  MeshBlock* pmb=pmy_block_;
  int &mylevel=pmb->loc.level;

  int f2d=0, f3d=0;
  int cng, cng1, cng2, cng3;
  if (pmb->block_size.nx2 > 1) f2d=1;
  if (pmb->block_size.nx3 > 1) f3d=1;
  cng  = cng1 =pmb->cnghost;
  cng2 = cng*f2d;
  cng3 = cng*f3d;
  int ssize, rsize;
  int tag;
  // Initialize non-polar neighbor communications to other ranks
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (nb.rank != Globals::my_rank) {
      if (nb.level == mylevel) { // same
        ssize=rsize=((nb.ox1 == 0)?pmb->block_size.nx1:NGHOST)
              *((nb.ox2 == 0)?pmb->block_size.nx2:NGHOST)
              *((nb.ox3 == 0)?pmb->block_size.nx3:NGHOST);
      } else if (nb.level<mylevel) { // coarser
        ssize=((nb.ox1 == 0) ? ((pmb->block_size.nx1+1)/2):NGHOST)
              *((nb.ox2 == 0) ? ((pmb->block_size.nx2+1)/2):NGHOST)
              *((nb.ox3 == 0) ? ((pmb->block_size.nx3+1)/2):NGHOST);
        rsize=((nb.ox1 == 0) ? ((pmb->block_size.nx1+1)/2+cng1):cng1)
              *((nb.ox2 == 0) ? ((pmb->block_size.nx2+1)/2+cng2):cng2)
              *((nb.ox3 == 0) ? ((pmb->block_size.nx3+1)/2+cng3):cng3);
      } else { // finer
        ssize=((nb.ox1 == 0) ? ((pmb->block_size.nx1+1)/2+cng1):cng1)
              *((nb.ox2 == 0) ? ((pmb->block_size.nx2+1)/2+cng2):cng2)
              *((nb.ox3 == 0) ? ((pmb->block_size.nx3+1)/2+cng3):cng3);
        rsize=((nb.ox1 == 0) ? ((pmb->block_size.nx1+1)/2):NGHOST)
              *((nb.ox2 == 0) ? ((pmb->block_size.nx2+1)/2):NGHOST)
              *((nb.ox3 == 0) ? ((pmb->block_size.nx3+1)/2):NGHOST);
      }
      ssize*=(nu_+1); rsize*=(nu_+1);
      // specify the offsets in the view point of the target block: flip ox? signs

      // Initialize persistent communication requests attached to specific BoundaryData
      // cell-centered hydro: bd_hydro_
      tag=pbval_->CreateBvalsMPITag(nb.lid, nb.targetid, cc_phys_id_);
      if (bd_var_.req_send[nb.bufid] != MPI_REQUEST_NULL)
        MPI_Request_free(&bd_var_.req_send[nb.bufid]);
      MPI_Send_init(bd_var_.send[nb.bufid],ssize,MPI_ATHENA_REAL,
                    nb.rank,tag,MPI_COMM_WORLD,&(bd_var_.req_send[nb.bufid]));
      tag=pbval_->CreateBvalsMPITag(pmb->lid, nb.bufid, cc_phys_id_);
      if (bd_var_.req_recv[nb.bufid] != MPI_REQUEST_NULL)
        MPI_Request_free(&bd_var_.req_recv[nb.bufid]);
      MPI_Recv_init(bd_var_.recv[nb.bufid],rsize,MPI_ATHENA_REAL,
                    nb.rank,tag,MPI_COMM_WORLD,&(bd_var_.req_recv[nb.bufid]));

      // hydro flux correction: bd_var_flcor_
      if (pmy_mesh_->multilevel == true && nb.type == NeighborConnect::face) {
        int size;
        if (nb.fid == 0 || nb.fid == 1)
          size=((pmb->block_size.nx2+1)/2)*((pmb->block_size.nx3+1)/2);
        else if (nb.fid == 2 || nb.fid == 3)
          size=((pmb->block_size.nx1+1)/2)*((pmb->block_size.nx3+1)/2);
        else // (nb.fid == 4 || nb.fid == 5)
          size=((pmb->block_size.nx1+1)/2)*((pmb->block_size.nx2+1)/2);
        size*=(nu_+1);
        if (nb.level<mylevel) { // send to coarser
          tag=pbval_->CreateBvalsMPITag(nb.lid, nb.targetid, cc_flx_phys_id_);
          if (bd_var_flcor_.req_send[nb.bufid] != MPI_REQUEST_NULL)
            MPI_Request_free(&bd_var_flcor_.req_send[nb.bufid]);
          MPI_Send_init(bd_var_flcor_.send[nb.bufid],size,MPI_ATHENA_REAL,
                        nb.rank,tag,MPI_COMM_WORLD,&(bd_var_flcor_.req_send[nb.bufid]));
        } else if (nb.level>mylevel) { // receive from finer
          tag=pbval_->CreateBvalsMPITag(pmb->lid, nb.bufid, cc_flx_phys_id_);
          if (bd_var_flcor_.req_recv[nb.bufid] != MPI_REQUEST_NULL)
            MPI_Request_free(&bd_var_flcor_.req_recv[nb.bufid]);
          MPI_Recv_init(bd_var_flcor_.recv[nb.bufid],size,MPI_ATHENA_REAL,
                        nb.rank,tag,MPI_COMM_WORLD,&(bd_var_flcor_.req_recv[nb.bufid]));
        }
      }
    }
  }
#endif
  return;
}

// KGF: approach #1: called inside #ifdef MPI_PARALLEL and loop:
// for (int n=0; n<nneighbor; n++) {

// -- unable to access nb without passing as function parameter
// void CellCenteredBoundaryVariable::StartReceivingAll() {
//   MPI_Start(&(bd_var_.req_recv[nb.bufid]));
//   if (nb.type == NeighborConnect::face && nb.level>mylevel)
//     MPI_Start(&(bd_var_flcor_.req_recv[nb.bufid]));
//   return;
// }

// KGF: approach #2: separate loops over nneighbor (make loop over bvar vector the
// outermost loop)
// --- performance penalty?
// --- redundant source code (could encapsulate loop over nneighbor, passing nb down to
// BoundaryVariable instances?)

// Are there any shared implemenations worth placing as the default inherited virtual
// function in BoundaryVariable
void CellCenteredBoundaryVariable::StartReceiving(BoundaryCommSubset phase) {
#ifdef MPI_PARALLEL
  MeshBlock *pmb = pmy_block_;
  int mylevel = pmb->loc.level;
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (nb.rank != Globals::my_rank) {
      MPI_Start(&(bd_var_.req_recv[nb.bufid]));
      if (phase == BoundaryCommSubset::all && nb.type == NeighborConnect::face
          && nb.level > mylevel) // opposite condition in ClearBoundary()
        MPI_Start(&(bd_var_flcor_.req_recv[nb.bufid]));
    }
  }
#endif
  return;
}

void CellCenteredBoundaryVariable::ClearBoundary(BoundaryCommSubset phase) {
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    bd_var_.flag[nb.bufid] = BoundaryStatus::waiting;
    if (nb.type == NeighborConnect::face)
      bd_var_flcor_.flag[nb.bufid] = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
    MeshBlock *pmb = pmy_block_;
    int mylevel = pmb->loc.level;
    if (nb.rank != Globals::my_rank) {
      // Wait for Isend
      MPI_Wait(&(bd_var_.req_send[nb.bufid]), MPI_STATUS_IGNORE);
      if (phase == BoundaryCommSubset::all && nb.type == NeighborConnect::face
          && nb.level < mylevel)
        MPI_Wait(&(bd_var_flcor_.req_send[nb.bufid]), MPI_STATUS_IGNORE);
    } // KGF: end block on other MPI process
#endif
  } // KGF: end loop over nneighbor
  return;
}
