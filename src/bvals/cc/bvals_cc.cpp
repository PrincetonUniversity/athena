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
#include "bvals_cc.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// constructor

CellCenteredBoundaryVariable::CellCenteredBoundaryVariable(
    MeshBlock *pmb, AthenaArray<Real> *var, AthenaArray<Real> *coarse_var,
    AthenaArray<Real> *var_flux)
    : BoundaryVariable(pmb), var_cc(var), coarse_buf(coarse_var), x1flux(var_flux[X1DIR]),
      x2flux(var_flux[X2DIR]), x3flux(var_flux[X3DIR]), nl_(0), nu_(var->GetDim4() -1),
      flip_across_pole_(nullptr) {
  // CellCenteredBoundaryVariable should only be used w/ 4D or 3D (nx4=1) AthenaArray
  // For now, assume that full span of 4th dim of input AthenaArray should be used:
  // ---> get the index limits directly from the input AthenaArray
  // <=nu_ (inclusive), <nx4 (exclusive)
  if (nu_ < 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in CellCenteredBoundaryVariable constructor" << std::endl
        << "An 'AthenaArray<Real> *var' of nx4_ = " << var->GetDim4() << " was passed\n"
        << "Should be nx4 >= 1 (likely uninitialized)." << std::endl;
    ATHENA_ERROR(msg);
  }

  InitBoundaryData(bd_var_, BoundaryQuantity::cc);
#ifdef MPI_PARALLEL
  // KGF: dead code, leaving for now:
  // cc_phys_id_ = pbval_->ReserveTagVariableIDs(1);
  cc_phys_id_ = pbval_->bvars_next_phys_id_;
#endif
  if (pmy_mesh_->multilevel) { // SMR or AMR
    InitBoundaryData(bd_var_flcor_, BoundaryQuantity::cc_flcor);
#ifdef MPI_PARALLEL
    cc_flx_phys_id_ = cc_phys_id_ + 1;
#endif
  }

  if (SHEARING_BOX) {
#ifdef MPI_PARALLEL
    shear_cc_phys_id_ = cc_phys_id_ + 2;
#endif
    if (pbval_->ShBoxCoord_ == 1) {
      int nc2 = pmb->ncells2;
      int nc3 = pmb->ncells3;
      for (int upper=0; upper<2; upper++) {
        if (pbval_->is_shear[upper]) {
          shear_cc_[upper].NewAthenaArray(nu_+1, nc3, nc2, NGHOST);
          shear_flx_cc_[upper].NewAthenaArray(nc2);

          // TODO(KGF): the rest of this should be a part of InitBoundaryData()

          // attach corner cells from L/R side
          int size = (pmb->block_size.nx2 + NGHOST)*pbval_->ssize_*(nu_ + 1);
          for (int n=0; n<2; n++) {
            shear_bd_var_[upper].send[n] = new Real[size];
            shear_bd_var_[upper].recv[n] = new Real[size];
            shear_bd_var_[upper].flag[n] = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
            shear_bd_var_[upper].req_send[n] = MPI_REQUEST_NULL;
            shear_bd_var_[upper].req_recv[n] = MPI_REQUEST_NULL;
#endif
          }
          // corner cells only
          size = NGHOST*pbval_->ssize_*(nu_ + 1);
          for (int n=2; n<4; n++) {
            shear_bd_var_[upper].send[n] = new Real[size];
            shear_bd_var_[upper].recv[n] = new Real[size];
            shear_bd_var_[upper].flag[n] = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
            shear_bd_var_[upper].req_send[n] = MPI_REQUEST_NULL;
            shear_bd_var_[upper].req_recv[n] = MPI_REQUEST_NULL;
#endif
          }
        } // end "if is a shearing boundary"
      }  // end loop over inner, outer shearing boundaries
    } // end "if (pbval_->ShBoxCoord_ == 1)"
  } // end shearing box component of ctor
}

// destructor

CellCenteredBoundaryVariable::~CellCenteredBoundaryVariable() {
  DestroyBoundaryData(bd_var_);
  if (pmy_mesh_->multilevel)
    DestroyBoundaryData(bd_var_flcor_);

  // TODO(KGF): this should be a part of DestroyBoundaryData()
  if (SHEARING_BOX) {
    for (int upper=0; upper<2; upper++) {
      if (pbval_->is_shear[upper]) { // if true for shearing inner blocks
        for (int n=0; n<4; n++) {
          delete[] shear_bd_var_[upper].send[n];
          delete[] shear_bd_var_[upper].recv[n];
        }
      }
    }
  }
}

int CellCenteredBoundaryVariable::ComputeVariableBufferSize(const NeighborIndexes& ni,
                                                            int cng) {
#ifdef DBGPR_BVALS_CC
  coutYellow("CellCenteredBoundaryVariable::ComputeVariableBufferSize\n");
#endif // DBGPR_BVALS_CC

  MeshBlock *pmb = pmy_block_;
  int cng1, cng2, cng3;
  cng1 = cng;
  cng2 = cng*(pmb->block_size.nx2 > 1 ? 1 : 0);
  cng3 = cng*(pmb->block_size.nx3 > 1 ? 1 : 0);

  int size = ((ni.ox1 == 0) ? pmb->block_size.nx1 : NGHOST)
    *((ni.ox2 == 0) ? pmb->block_size.nx2 : NGHOST)
    *((ni.ox3 == 0) ? pmb->block_size.nx3 : NGHOST);
  if (pmy_mesh_->multilevel) {
    int f2c = ((ni.ox1 == 0) ? ((pmb->block_size.nx1+1)/2) : NGHOST)
      *((ni.ox2 == 0) ? ((pmb->block_size.nx2+1)/2) : NGHOST)
      *((ni.ox3 == 0) ? ((pmb->block_size.nx3+1)/2) : NGHOST);
    int c2f = ((ni.ox1 == 0) ?((pmb->block_size.nx1+1)/2 + cng1) : cng)
      *((ni.ox2 == 0) ? ((pmb->block_size.nx2+1)/2 + cng2) : cng)
      *((ni.ox3 == 0) ? ((pmb->block_size.nx3+1)/2 + cng3) : cng);
    size = std::max(size, c2f);
    size = std::max(size, f2c);
  }
  size *= nu_ + 1;

  return size;
}

int CellCenteredBoundaryVariable::ComputeFluxCorrectionBufferSize(
    const NeighborIndexes& ni, int cng) {
#ifdef DBGPR_BVALS_CC
  coutYellow("CellCenteredBoundaryVariable::ComputeFluxCorrectionBufferSize\n");
#endif // DBGPR_BVALS_CC

  MeshBlock *pmb = pmy_block_;
  int size = 0;
  if (ni.ox1 != 0)
    size = (pmb->block_size.nx2 + 1)/2*(pmb->block_size.nx3 + 1)/2*(nu_ + 1);
  if (ni.ox2 != 0)
    size = (pmb->block_size.nx1 + 1)/2*(pmb->block_size.nx3 + 1)/2*(nu_ + 1);
  if (ni.ox3 != 0)
    size = (pmb->block_size.nx1 + 1)/2*(pmb->block_size.nx2 + 1)/2*(nu_ + 1);
  return size;
}

//----------------------------------------------------------------------------------------
//! \fn int CellCenteredBoundaryVariable::LoadBoundaryBufferSameLevel(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set cell-centered boundary buffers for sending to a block on the same level

int CellCenteredBoundaryVariable::LoadBoundaryBufferSameLevel(Real *buf,
                                                              const NeighborBlock& nb) {
#ifdef DBGPR_BVALS_CC
  coutYellow("CellCenteredBoundaryVariable::LoadBoundaryBufferSameLevel\n");
#endif // DBGPR_BVALS_CC

  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;

  si = (nb.ni.ox1 > 0) ? (pmb->ie - NGHOST + 1) : pmb->is;
  ei = (nb.ni.ox1 < 0) ? (pmb->is + NGHOST - 1) : pmb->ie;
  sj = (nb.ni.ox2 > 0) ? (pmb->je - NGHOST + 1) : pmb->js;
  ej = (nb.ni.ox2 < 0) ? (pmb->js + NGHOST - 1) : pmb->je;
  sk = (nb.ni.ox3 > 0) ? (pmb->ke - NGHOST + 1) : pmb->ks;
  ek = (nb.ni.ox3 < 0) ? (pmb->ks + NGHOST - 1) : pmb->ke;
  int p = 0;
  AthenaArray<Real> &var = *var_cc;

  BufferUtility::PackData(var, buf, nl_, nu_, si, ei, sj, ej, sk, ek, p);

  //////////////////////////////////////////////////////////////////////////////
  // inspect branching logic SL
  //
  // [complementary fcn] -> (current fcn)
  // (LoadBoundaryBufferSameLevel) -> [SetBoundarySameLevel]
  //
  // semi-circular indicate current scope
  // square indicate complementary function condition

  // modifies si, ei
  if (nb.ni.ox1 == 0) {
    // (ox1=0, ox2, ox3) [ox1=0, -ox2, -ox3]
    //
    // (si, ei) = (is, ie)
    true;
  } else if (nb.ni.ox1 > 0) {
    // (ox1>0, ox2, ox3) [ox1<0, -ox2, -ox3]
    //
    // (si, ei) = (ie - NGHOST + 1, ie)
    true;
  } else {
    // (ox1<0, ox2, ox3) [ox1>0, -ox2, -ox3]
    //
    // (si, ei) = (is, is + NGHOST - 1)
    true;
  }

  // modifies sj, ej
  if (nb.ni.ox2 == 0) {
    // (ox1, ox2=0, ox3) [-ox1, ox2=0, -ox3]
    //
    // (sj, ej) = (js, je)
    true;
  } else if (nb.ni.ox2 > 0) {
    // (ox1, ox2>0, ox3) [-ox1, ox2<0, -ox3]
    //
    // (sj, ej) = (je - NGHOST + 1, je)
    true;
  } else {
    // (ox1, ox2<0, ox3) [-ox1, ox2>0, -ox3]
    //
    // (sj, ej) = (js, js + NGHOST - 1)
    true;
  }

  // modifies sk, ek
  if (nb.ni.ox3 == 0) {
    // (ox1, ox2, ox3=0) [-ox1, -ox2, ox3=0]
    //
    // (sk, ek) = (ks, ke)
    true;
  } else if (nb.ni.ox3 > 0) {
    // (ox1, ox2, ox3>0) [-ox1, -ox2, ox3<0]
    //
    // (sk, ek) = (ke - NGHOST + 1, ke)
    true;
  } else {
    // (ox1, ox2, ox3<0) [-ox1, -ox2, ox3>0]
    //
    // (sk, ek) = (ks, ks + NGHOST - 1)
    true;
  }
  //////////////////////////////////////////////////////////////////////////////

    return p;
}

//----------------------------------------------------------------------------------------
//! \fn int CellCenteredBoundaryVariable::LoadBoundaryBufferToCoarser(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set cell-centered boundary buffers for sending to a block on the coarser level

int CellCenteredBoundaryVariable::LoadBoundaryBufferToCoarser(Real *buf,
                                                              const NeighborBlock& nb) {
#ifdef DBGPR_BVALS_CC
    coutYellow("CellCenteredBoundaryVariable::LoadBoundaryBufferToCoarser\n");
#endif // DBGPR_BVALS_CC

  MeshBlock *pmb = pmy_block_;
  MeshRefinement *pmr = pmb->pmr;
  int si, sj, sk, ei, ej, ek;
  int cn = NGHOST - 1;
  AthenaArray<Real> &var = *var_cc;
  AthenaArray<Real> &coarse_var = *coarse_buf;

  si = (nb.ni.ox1 > 0) ? (pmb->cie - cn) : pmb->cis;
  ei = (nb.ni.ox1 < 0) ? (pmb->cis + cn) : pmb->cie;
  sj = (nb.ni.ox2 > 0) ? (pmb->cje - cn) : pmb->cjs;
  ej = (nb.ni.ox2 < 0) ? (pmb->cjs + cn) : pmb->cje;
  sk = (nb.ni.ox3 > 0) ? (pmb->cke - cn) : pmb->cks;
  ek = (nb.ni.ox3 < 0) ? (pmb->cks + cn) : pmb->cke;

  int p = 0;

  pmr->RestrictCellCenteredValues(var, coarse_var, nl_, nu_, si, ei, sj, ej, sk, ek);
  BufferUtility::PackData(coarse_var, buf, nl_, nu_, si, ei, sj, ej, sk, ek, p);

  //////////////////////////////////////////////////////////////////////////////
  // inspect branching logic Fine2Coarse
  //
  // [complementary fcn] -> (current fcn)
  // (LoadBoundaryBufferToCoarser) -> [SetBoundaryFromFiner]
  //
  // semi-circular indicate current scope
  // square indicate complementary function condition

  // [1d] modifies si, ei
  if (nb.ni.ox1 == 0) {
    //
    // (si, ei) = (cis, cie)
  } else if (nb.ni.ox1 > 0)  {
    //
    // (si, ei) = (cie - NGHOST + 1, cie)
    true;
  } else {
    //
    // (si, ei) = (cis, cis + NGHOST - 1)
    true;
  }

  // [2d] modifies sj, ej
  if (nb.ni.ox2 == 0) {
    //
    // (sj, ej) = (cjs, cje)
  } else if (nb.ni.ox2 > 0) {
    //
    // (sj, ej) = (cje - NGHOST + 1, cje)
    true;
  } else {
    //
    // (sj, ej) = (cjs, cjs + NGHOST - 1)
    true;
  }

  // [3d] modifies sk, ek
  if (nb.ni.ox3 == 0) {
    //
    // (sk, ek) = (cks, cke)
  } else if (nb.ni.ox3 > 0)  {
    //
    // (sk, ek) = (cke - NGHOST + 1, cke)
    true;
  } else {
    //
    // (sk, ek) = (cks, cks + NGHOST - 1)
    true;
  }
  //////////////////////////////////////////////////////////////////////////////



  // if (nb.ni.ox1 == -1)
  //   coutBoldBlue("blk\n"), Q();

  return p;
}

//----------------------------------------------------------------------------------------
//! \fn int CellCenteredBoundaryVariable::LoadBoundaryBufferToFiner(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set cell-centered boundary buffers for sending to a block on the finer level

int CellCenteredBoundaryVariable::LoadBoundaryBufferToFiner(Real *buf,
                                                            const NeighborBlock& nb) {
#ifdef DBGPR_BVALS_CC
    coutYellow("CellCenteredBoundaryVariable::LoadBoundaryBufferToFiner\n");
#endif // DBGPR_BVALS_CC

  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;
  int cn = pmb->cnghost - 1;
  AthenaArray<Real> &var = *var_cc;

  si = (nb.ni.ox1 > 0) ? (pmb->ie - cn) : pmb->is;
  ei = (nb.ni.ox1 < 0) ? (pmb->is + cn) : pmb->ie;
  sj = (nb.ni.ox2 > 0) ? (pmb->je - cn) : pmb->js;
  ej = (nb.ni.ox2 < 0) ? (pmb->js + cn) : pmb->je;
  sk = (nb.ni.ox3 > 0) ? (pmb->ke - cn) : pmb->ks;
  ek = (nb.ni.ox3 < 0) ? (pmb->ks + cn) : pmb->ke;

  // send the data first and later prolongate on the target block
  // need to add edges for faces, add corners for edges
  if (nb.ni.ox1 == 0) {
    if (nb.ni.fi1 == 1)   si += pmb->block_size.nx1/2 - pmb->cnghost;
    else            ei -= pmb->block_size.nx1/2 - pmb->cnghost;
  }
  if (nb.ni.ox2 == 0 && pmb->block_size.nx2 > 1) {
    if (nb.ni.ox1 != 0) {
      if (nb.ni.fi1 == 1) sj += pmb->block_size.nx2/2 - pmb->cnghost;
      else          ej -= pmb->block_size.nx2/2 - pmb->cnghost;
    } else {
      if (nb.ni.fi2 == 1) sj += pmb->block_size.nx2/2 - pmb->cnghost;
      else          ej -= pmb->block_size.nx2/2 - pmb->cnghost;
    }
  }
  if (nb.ni.ox3 == 0 && pmb->block_size.nx3 > 1) {
    if (nb.ni.ox1 != 0 && nb.ni.ox2 != 0) {
      if (nb.ni.fi1 == 1) sk += pmb->block_size.nx3/2 - pmb->cnghost;
      else          ek -= pmb->block_size.nx3/2 - pmb->cnghost;
    } else {
      if (nb.ni.fi2 == 1) sk += pmb->block_size.nx3/2 - pmb->cnghost;
      else          ek -= pmb->block_size.nx3/2 - pmb->cnghost;
    }
  }

  int p = 0;

  BufferUtility::PackData(var, buf, nl_, nu_, si, ei, sj, ej, sk, ek, p);

  // if ((nb.ni.fi1 == 1) && (nb.ni.ox1 == 0) && (nb.ni.ox2 > 0))
  //   Q();
  // if ((nb.ni.ox1 == 0) && (nb.ni.ox2 < 0) && (nb.ni.fi1 == 1))
  //   Q();

  //////////////////////////////////////////////////////////////////////////////
  // inspect branching logic Coarse2Fine
  //
  // [complementary fcn] -> (current fcn)
  // (LoadBoundaryBufferToFiner) -> [SetBoundaryFromCoarser]
  //
  // semi-circular indicate current scope
  // square indicate complementary function condition

  // [1d] modifies si, ei
  if (nb.ni.ox1 == 0) {
    if (nb.ni.fi1 == 1) {
      // (ox1=0, ~(ox2=0 /\ ox3=0), fi1=1, fi2)
      //
      // (si, ei) = (is + block_size.nx1 / 2 - cnghost, ie)
      true;
    } else {
      // (ox1=0, ~(ox2=0 /\ ox3=0), fi1=0, fi2)
      //
      // (si, ei) = (is, ie - block_size.nx1 / 2 + cnghost)
      true;
    }
  } else if (nb.ni.ox1 > 0) {
    // (ox1>0, ox2, ox3, fi1, fi2)
    //
    // (si, ei) = (ie - cnghost + 1, ie)
    true;
  } else {
    // (ox1<0, ox2, ox3, fi1, fi2)
    // (si, ei) = (is, is + cnghost - 1)
    true;
  }

  // [2d] modifies sj, ej
  if (nb.ni.ox2 == 0 && pmb->block_size.nx2 > 1) {
    if (nb.ni.ox1 != 0) {
      if (nb.ni.fi1 == 1) {
        // (ox1!=0, ox2=0, ox3, fi1=1, fi2)
        //
        // (sj, ej) = (js + block_size.nx2 / 2 - cnghost, je)
        true;
      } else {
        // (ox1!=0, ox2=0, ox3, fi1=0, fi2)
        //
        // (sj, ej) = (js, je - block_size.nx2 / 2 + cnghost)
        true;
      }
    } else {
      if (nb.ni.fi2 == 1) {
        // (ox1=0, ox2=0, ox3!=0, fi1, fi2=1)
        //
        // (sj, ej) = (js + block_size.nx2 / 2 - cnghost, je)
        true;
      } else {
        // (ox1=0, ox2=0, ox3!=0, fi1, fi2=0)
        //
        // (sj, ej) = (js, je - block_size.nx2 / 2 + cnghost)
        true;
      }
    }
  } else if (nb.ni.ox2 > 0) {
    // (ox1, ox2>0, ox3, fi1, fi2)
    //
    // (sj, ej) = (je - cnghost + 1, je)
    true;
  } else {
    // (ox1, ox2<0, ox3, fi1, fi2)
    //
    // (sj, ej) = (js, js + cnghost - 1)
    true;
  }

  // [3d] modifies sj, ek
  if (nb.ni.ox3 == 0 && pmb->block_size.nx3 > 1) {
    if (nb.ni.ox1 != 0 && nb.ni.ox2 != 0) {
      if (nb.ni.fi1 == 1) {
        // (ox1!=0, ox2!=0, ox3=0, fi1=1, fi2)
        //
        // (sk, ek)= (ks + block_size.nx3 / 2 - cnghost, ke)
        true;
      } else {
        // (ox1!=0, ox2!=0, ox3=0, fi1=0, fi2)
        //
        // (sk, ek)= (ks, ke - block_size.nx3 / 2 + cnghost)
        true;
      }
    } else {
      if (nb.ni.fi2 == 1) {
        // (~(ox1!=0 /\ ox2!=0), ox3=0, fi1, fi2=1)
        //
        // (sk, ek)= (ks + block_size.nx3 / 2 - cnghost, ke)
        true;
      } else {
        // (~(ox1!=0 /\ ox2!=0), ox3=0, fi1, fi2=0)
        //
        // (sk, ek)= (ks, ke - block_size.nx3 / 2 + cnghost)
        true;
      }
    }
  } else if (nb.ni.ox3 > 0) {
    // (ox1, ox2, ox3>0, fi1, fi2)
    //
    // (sk, ek) = (ke - cnghost + 1, ke)
    true;
  } else {
    // (ox1, ox2, ox3<0, fi1, fi2)
    //
    // (sk, ek) = (ks, ks + cnghost - 1)
    true;
  }
  //////////////////////////////////////////////////////////////////////////////

  return p;
}


//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::SetBoundarySameLevel(Real *buf,
//                                                              const NeighborBlock& nb)
//  \brief Set cell-centered boundary received from a block on the same level

void CellCenteredBoundaryVariable::SetBoundarySameLevel(Real *buf,
                                                        const NeighborBlock& nb) {
#ifdef DBGPR_BVALS_CC
    coutYellow("CellCenteredBoundaryVariable::SetBoundarySameLevel\n");
#endif // DBGPR_BVALS_CC

  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;
  AthenaArray<Real> &var = *var_cc;

  if (nb.ni.ox1 == 0)     si = pmb->is,        ei = pmb->ie;
  else if (nb.ni.ox1 > 0) si = pmb->ie + 1,      ei = pmb->ie + NGHOST;
  else              si = pmb->is - NGHOST, ei = pmb->is - 1;
  if (nb.ni.ox2 == 0)     sj = pmb->js,        ej = pmb->je;
  else if (nb.ni.ox2 > 0) sj = pmb->je + 1,      ej = pmb->je + NGHOST;
  else              sj = pmb->js - NGHOST, ej = pmb->js - 1;
  if (nb.ni.ox3 == 0)     sk = pmb->ks,        ek = pmb->ke;
  else if (nb.ni.ox3 > 0) sk = pmb->ke + 1,      ek = pmb->ke + NGHOST;
  else              sk = pmb->ks - NGHOST, ek = pmb->ks - 1;

  //////////////////////////////////////////////////////////////////////////////
  // inject actual solution (see below)
  int p = 0;

  if (nb.polar) {
    for (int n=nl_; n<=nu_; ++n) {
      Real sign = 1.0;
      if (flip_across_pole_ != nullptr) sign = flip_across_pole_[n] ? -1.0 : 1.0;
      for (int k=sk; k<=ek; ++k) {
        for (int j=ej; j>=sj; --j) {
#pragma omp simd linear(p)
          for (int i=si; i<=ei; ++i) {
            var(n,k,j,i) = sign * buf[p++];
          }
        }
      }
    }
  } else {
    BufferUtility::UnpackData(buf, var, nl_, nu_, si, ei, sj, ej, sk, ek, p);
  }

  // printf("unpacking: (p, si, ei)=(%d, %d, %d)\n", p, si, ei);
  // //printf("%d, %d\n", nl_, nu_);

  // for (int kk=0; kk<NGHOST; kk++){
  //   printf("%1.3f, ", buf[kk]);
  // };
  // printf("\n");

  // var.print_data();

  if (SHEARING_BOX) {
    // 2D shearing box in x-z plane: additional step to shift azimuthal velocity
    if (pbval_->ShBoxCoord_ == 2) {
      int sign[2]{1, -1};
      Real qomL = pbval_->qomL_;
      for (int upper=0; upper<2; upper++) {
        if ((pmb->loc.lx1 == pbval_->loc_shear[upper]) && (sign[upper]*nb.ni.ox1 < 0)) {
          for (int k=sk; k<=ek; ++k) {
            for (int j=sj; j<=ej; ++j) {
              for (int i=si; i<=ei; ++i) {
                if (NON_BAROTROPIC_EOS)
                  var(IEN,k,j,i) += (0.5/var(IDN,k,j,i))*(
                      SQR(var(IM3,k,j,i) + sign[upper]*qomL*var(IDN,k,j,i))
                      - SQR(var(IM3,k,j,i)));
                var(IM3,k,j,i) += sign[upper]*qomL*var(IDN,k,j,i);
              }
            }
          }
        }
      }
    }
  } // end KGF: shearing box in SetBoundarySameLevel
  //////////////////////////////////////////////////////////////////////////////


  //////////////////////////////////////////////////////////////////////////////
  // inspect branching logic SL
  //
  // [complementary fcn] -> (current fcn)
  // [LoadBoundaryBufferSameLevel] -> (SetBoundarySameLevel)
  //
  // semi-circular indicate current scope
  // square indicate complementary function condition

  // modifies si, ei
  if (nb.ni.ox1 == 0) {
    // (ox1=0, ox2, ox3)                [ox1=0, -ox2, -ox3]
    //
    // (si, ei) = (is, ie)               [is, ie]
    true;
  } else if (nb.ni.ox1 > 0) {
    // (ox1>0, ox2, ox3)                [ox1<0, -ox2, -ox3]
    //
    // (si, ei) = (ie + 1, ie + NGHOST)  [is, is + NGHOST - 1]
    true;
  } else {
    // (ox1<0, ox2, ox3)                [ox1>0, -ox2, -ox3]
    //
    // (si, ei) = (is - NGHOST, is - 1)  [ie - NGHOST + 1, ie]
    true;
  }

  // modifies sj, ej
  if (nb.ni.ox2 == 0) {
    // (ox1, ox2=0, ox3)                [-ox1, ox2=0, -ox3]
    //
    // (sj, ej) = (js, je)               [js, je]
    true;
  } else if (nb.ni.ox2 > 0) {
    // (ox1, ox2>0, ox3)                [-ox1, ox2<0, -ox3]
    //
    // (sj, ej) = (je + 1, je + NGHOST)  [js, js + NGHOST - 1]
    true;
  } else {
    // (ox1, ox2<0, ox3)                [-ox1, ox2>0, -ox3]
    //
    // (sj, ej) = (js - NGHOST, js - 1)  [je - NGHOST + 1, je]
    true;
  }

  // modifies sk, ek
  if (nb.ni.ox3 == 0) {
    // (ox1, ox2, ox3=0)                [-ox1, -ox2, ox3=0]
    //
    // (sk, ek) = (ks, ke)               [ks, ke]
    true;
  } else if (nb.ni.ox3 > 0) {
    // (ox1, ox2, ox3>0)                [-ox1, -ox2, ox3<0]
    //
    // (sk, ek) = (ke + 1, ke + NGHOST)  [ks, ks + NGHOST - 1]
    true;
  } else {
    // (ox1, ox2, ox3<0)                [-ox1, -ox2, ox3>0]
    //
    // (sk, ek) = (ks - NGHOST, ks - 1)  [ke - NGHOST + 1, ke]
    true;
  }
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  // BD: debug - populate based on solution
#ifdef FILL_WAVE_BND_SL
  pmb->DebugWaveMeshBlock(var, si, ei, sj, ej, sk, ek, false);
#endif // FILL_WAVE_BND_SL
  //////////////////////////////////////////////////////////////////////////////
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::SetBoundaryFromCoarser(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set cell-centered prolongation buffer received from a block on a coarser level

void CellCenteredBoundaryVariable::SetBoundaryFromCoarser(Real *buf,
                                                          const NeighborBlock& nb) {
#ifdef DBGPR_BVALS_CC
    coutYellow("CellCenteredBoundaryVariable::SetBoundaryFromCoarser\n");
#endif // DBGPR_BVALS_CC

  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;
  int cng = pmb->cnghost;
  AthenaArray<Real> &coarse_var = *coarse_buf;

  if (nb.ni.ox1 == 0) {
    si = pmb->cis, ei = pmb->cie;
    if ((pmb->loc.lx1 & 1LL) == 0LL) ei += cng;
    else                             si -= cng;
  } else if (nb.ni.ox1 > 0)  {
    si = pmb->cie + 1,   ei = pmb->cie + cng;
  } else {
    si = pmb->cis - cng, ei = pmb->cis - 1;
  }
  if (nb.ni.ox2 == 0) {
    sj = pmb->cjs, ej = pmb->cje;
    if (pmb->block_size.nx2 > 1) {
      if ((pmb->loc.lx2 & 1LL) == 0LL) ej += cng;
      else                             sj -= cng;
    }
  } else if (nb.ni.ox2 > 0) {
    sj = pmb->cje + 1,   ej = pmb->cje + cng;
  } else {
    sj = pmb->cjs - cng, ej = pmb->cjs - 1;
  }
  if (nb.ni.ox3 == 0) {
    sk = pmb->cks, ek = pmb->cke;
    if (pmb->block_size.nx3 > 1) {
      if ((pmb->loc.lx3 & 1LL) == 0LL) ek += cng;
      else                             sk -= cng;
    }
  } else if (nb.ni.ox3 > 0)  {
    sk = pmb->cke + 1,   ek = pmb->cke + cng;
  } else {
    sk = pmb->cks - cng, ek = pmb->cks - 1;
  }

  //////////////////////////////////////////////////////////////////////////////
  // inject actual solution (see below)
  int p = 0;
  if (nb.polar) {
    for (int n=nl_; n<=nu_; ++n) {
      Real sign = 1.0;
      if (flip_across_pole_ != nullptr) sign = flip_across_pole_[n] ? -1.0 : 1.0;
      for (int k=sk; k<=ek; ++k) {
        for (int j=ej; j>=sj; --j) {
#pragma omp simd linear(p)
          for (int i=si; i<=ei; ++i)
            coarse_var(n,k,j,i) = sign * buf[p++];
        }
      }
    }
  } else {
    BufferUtility::UnpackData(buf, coarse_var, nl_, nu_, si, ei, sj, ej, sk, ek, p);
  }
  //////////////////////////////////////////////////////////////////////////////
  // if ((nb.ni.ox1 == 0) && (nb.ni.ox2 > 0))
  //   Q();


  //////////////////////////////////////////////////////////////////////////////
  // inspect branching logic Coarse2Fine
  //
  // [complementary fcn] -> (current fcn)
  // [LoadBoundaryBufferToFiner] -> (SetBoundaryFromCoarser)
  //
  // semi-circular indicate current scope
  // square indicate complementary function condition

  // note: non-coarse indices for fill_wave debug

  // [1d] modifies si, ei
  if (nb.ni.ox1 == 0) {
    //
    // (si, ei)= (cis, cie)
    if ((pmb->loc.lx1 & 1LL) == 0LL) {
      //
      // (si, ei) = (cis, cie + cnghost)
      si = pmb->is;
      ei = pmb->ie + NGHOST;
    } else {
      //
      // (si, ei) = (cis - cnghost, cie)
      si = pmb->is - NGHOST;
      ei = pmb->ie;
    }
  } else if (nb.ni.ox1 > 0)  {
    //
    // (si, ei) = (cie + 1, cie + cnghost)
    si = pmb->ie + 1;
    ei = pmb->ie + NGHOST;
  } else {
    //
    // (si, ei) = (cis - cnghost, cis - 1)
    si = pmb->is - NGHOST;
    ei = pmb->is - 1;
  }

  // [2d] modifies sj, ej
  if (nb.ni.ox2 == 0) {
    //
    // (sj, ej) = (cjs, cje)
    if (pmb->block_size.nx2 > 1) {
      if ((pmb->loc.lx2 & 1LL) == 0LL) {
        //
        // (sj, ej) = (cjs, cje + cnghost)
        sj = pmb->js;
        ej = pmb->je + NGHOST;
      } else {
        //
        // (sj, ej) = (cjs - cnghost, cje)
        sj = pmb->js - NGHOST;
        ej = pmb->je;
      }
    }
  } else if (nb.ni.ox2 > 0) {
    //
    // (sj, ej) = (cje + 1, cje + cnghost)
    sj = pmb->je + 1;
    ej = pmb->je + NGHOST;
  } else {
    //
    // (sj, ej) = (cjs - cnghost, cjs - 1)
    sj = pmb->js - NGHOST;
    ej = pmb->js - 1;
  }

  // [3d] modifies sk, ek
  if (nb.ni.ox3 == 0) {
    //
    // (sk, ek) = (cks, cke)
    if (pmb->block_size.nx3 > 1) {
      if ((pmb->loc.lx3 & 1LL) == 0LL) {
        //
        // (sk, ek) = (cks, cke + cnghost)
        sk = pmb->ks;
        ek = pmb->ke + NGHOST;
      } else {
        //
        // (sk, ek) = (cks - cnghost, cke)
        sk = pmb->ks - NGHOST;
        ek = pmb->ke;
      }
    }
  } else if (nb.ni.ox3 > 0)  {
    //
    // (sk, ek) = (cke + 1, cke + cnghost)
    sk = pmb->ke + 1;
    ek = pmb->ke + NGHOST;
  } else {
    //
    // (sk, ek) = (cks - cnghost, cks - 1)
    sk = pmb->ks - NGHOST;
    ek = pmb->ks - 1;
  }
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  // BD: debug - populate based on solution
  // [note we directly modify fund. not coarse]
#ifdef FILL_WAVE_BND_FRC
  AthenaArray<Real> &var = *var_cc;
  pmb->DebugWaveMeshBlock(var, si, ei, sj, ej, sk, ek, false);
#endif // FILL_WAVE_BND_FRC
  //////////////////////////////////////////////////////////////////////////////
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::SetBoundaryFromFiner(Real *buf,
//                                                              const NeighborBlock& nb)
//  \brief Set cell-centered boundary received from a block on a finer level

void CellCenteredBoundaryVariable::SetBoundaryFromFiner(Real *buf,
                                                        const NeighborBlock& nb) {
#ifdef DBGPR_BVALS_CC
    coutYellow("CellCenteredBoundaryVariable::SetBoundaryFromFiner\n");
#endif // DBGPR_BVALS_CC

  MeshBlock *pmb = pmy_block_;
  AthenaArray<Real> &var = *var_cc;
  // receive already restricted data
  int si, sj, sk, ei, ej, ek;

  if (nb.ni.ox1 == 0) {
    si = pmb->is, ei = pmb->ie;
    if (nb.ni.fi1 == 1)   si += pmb->block_size.nx1/2;
    else            ei -= pmb->block_size.nx1/2;
  } else if (nb.ni.ox1 > 0) {
    si = pmb->ie + 1,      ei = pmb->ie + NGHOST;
  } else {
    si = pmb->is - NGHOST, ei = pmb->is - 1;
  }
  if (nb.ni.ox2 == 0) {
    sj = pmb->js, ej = pmb->je;
    if (pmb->block_size.nx2 > 1) {
      if (nb.ni.ox1 != 0) {
        if (nb.ni.fi1 == 1) sj += pmb->block_size.nx2/2;
        else          ej -= pmb->block_size.nx2/2;
      } else {
        if (nb.ni.fi2 == 1) sj += pmb->block_size.nx2/2;
        else          ej -= pmb->block_size.nx2/2;
      }
    }
  } else if (nb.ni.ox2 > 0) {
    sj = pmb->je + 1,      ej = pmb->je + NGHOST;
  } else {
    sj = pmb->js - NGHOST, ej = pmb->js - 1;
  }
  if (nb.ni.ox3 == 0) {
    sk = pmb->ks, ek = pmb->ke;
    if (pmb->block_size.nx3 > 1) {
      if (nb.ni.ox1 != 0 && nb.ni.ox2 != 0) {
        if (nb.ni.fi1 == 1) sk += pmb->block_size.nx3/2;
        else          ek -= pmb->block_size.nx3/2;
      } else {
        if (nb.ni.fi2 == 1) sk += pmb->block_size.nx3/2;
        else          ek -= pmb->block_size.nx3/2;
      }
    }
  } else if (nb.ni.ox3 > 0) {
    sk = pmb->ke + 1,      ek = pmb->ke + NGHOST;
  } else {
    sk = pmb->ks - NGHOST, ek = pmb->ks - 1;
  }

  //////////////////////////////////////////////////////////////////////////////
  // inject actual solution (see below)
  int p = 0;
  if (nb.polar) {
    for (int n=nl_; n<=nu_; ++n) {
      Real sign=1.0;
      if (flip_across_pole_ != nullptr) sign = flip_across_pole_[n] ? -1.0 : 1.0;
      for (int k=sk; k<=ek; ++k) {
        for (int j=ej; j>=sj; --j) {
#pragma omp simd linear(p)
          for (int i=si; i<=ei; ++i)
            var(n,k,j,i) = sign * buf[p++];
        }
      }
    }
  } else {
    BufferUtility::UnpackData(buf, var, nl_, nu_, si, ei, sj, ej, sk, ek, p);
  }
  //////////////////////////////////////////////////////////////////////////////


  //////////////////////////////////////////////////////////////////////////////
  // inspect branching logic Fine2Coarse
  //
  // [complementary fcn] -> (current fcn)
  // [LoadBoundaryBufferToCoarser] -> (SetBoundaryFromFiner)
  //
  // semi-circular indicate current scope
  // square indicate complementary function condition

  // [1d] modifies si, ei
  if (nb.ni.ox1 == 0) {
    if (nb.ni.fi1 == 1) {
      // (ox1=0, ~(ox2=0 /\ ox3=0), fi1=1, fi2)
      //
      // (si, ei) = (is + block_size.nx1 / 2, ie)
      true;
    } else {
      // (ox1=0, ~(ox2=0 /\ ox3=0), fi1=0, fi2)
      //
      // (si, ei) = (is, ie - block_size.nx1 / 2)
      true;
    }
  } else if (nb.ni.ox1 > 0) {
    // (ox1>0, ox2, ox3, fi1, fi2)
    //
    // (si, ei) = (ie + 1, ie + NGHOST)
    true;
  } else {
    // (ox1<0, ox2, ox3, fi1, fi2)
    //
    // (si, ei) = (is - NGHOST, is - 1)
    true;
  }

  // [2d] modifies sj, ej
  if (nb.ni.ox2 == 0 && pmb->block_size.nx2 > 1) {
    if (nb.ni.ox1 != 0) {
      if (nb.ni.fi1 == 1) {
        // (ox1!=0, ox2=0, ox3, fi1=1, fi2)
        //
        // (sj, ej) = (js + block_size.nx2 / 2, je)
        true;
      } else {
        // (ox1!=0, ox2=0, ox3, fi1=0, fi2)
        //
        // (sj, ej) = (js, je - block_size.nx2 / 2)
        true;
      }
    } else {
      if (nb.ni.fi2 == 1) {
        // (ox1=0, ox2=0, ox3!=0, fi1, fi2=1)
        //
        // (sj, ej) = (js + block_size.nx2 / 2, je)
        true;
      } else {
        // (ox1=0, ox2=0, ox3!=0, fi1, fi2=0)
        //
        // (sj, ej) = (js, je - block_size.nx2 / 2)
        true;
      }
    }
  } else if (nb.ni.ox2 > 0) {
    // (ox1, ox2>0, ox3, fi1, fi2)
    //
    // (sj, ej) = (je + 1, je + NGHOST)
    true;
  } else {
    // (ox1, ox2<0, ox3, fi1, fi2)
    //
    // (sj, ej) = (js - NGHOST, js - 1)
    true;
  }

  // [3d] modifies sj, ek
  if (nb.ni.ox3 == 0 && pmb->block_size.nx3 > 1) {
    if (nb.ni.ox1 != 0 && nb.ni.ox2 != 0) {
      if (nb.ni.fi1 == 1) {
        // (ox1!=0, ox2!=0, ox3=0, fi1=1, fi2)
        //
        // (sk, ek) = (ks + block_size.nx3 / 2, ke)
        true;
      } else {
        // (ox1!=0, ox2!=0, ox3=0, fi1=0, fi2)
        //
        // (sk, ek) = (ks, ke - block_size.nx3 / 2)
        true;
      }
    } else {
      if (nb.ni.fi2 == 1) {
        // (~(ox1!=0 /\ ox2!=0), ox3=0, fi1, fi2=1)
        //
        // (sk, ek) = (ks + block_size.nx3 / 2, ke)
        true;
      } else {
        // (~(ox1!=0 /\ ox2!=0), ox3=0, fi1, fi2=0)
        //
        // (sk, ek) = (ks, ke - block_size.nx3 / 2)
        true;
      }
    }
  } else if (nb.ni.ox3 > 0) {
    // (ox1, ox2, ox3>0, fi1, fi2)
    //
    // (sk, ek) = (ke + 1, ke + NGHOST)
    true;
  } else {
    // (ox1, ox2, ox3<0, fi1, fi2)
    //
    // (sk, ek) = (ks - NGHOST, ks - 1)
    true;
  }
  //////////////////////////////////////////////////////////////////////////////


  // if (nb.ni.ox1 == -1)
  //   Q();
  // if ((nb.ni.ox1 == -1) and (nb.ni.ox2 == 0))
  //   Q();

  //////////////////////////////////////////////////////////////////////////////
  // BD: debug - populate based on solution
#ifdef FILL_WAVE_BND_FRF
  pmb->DebugWaveMeshBlock(var, si, ei, sj, ej, sk, ek, false);
#endif // FILL_WAVE_BND_FRF
  //////////////////////////////////////////////////////////////////////////////


  return;
  }


//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::PolarBoundarySingleAzimuthalBlock()
// \brief polar boundary edge-case: single MeshBlock spans the entire azimuthal (x3) range

void CellCenteredBoundaryVariable::PolarBoundarySingleAzimuthalBlock() {
#ifdef DBGPR_BVALS_CC
  coutYellow("CellCenteredBoundaryVariable::PolarBoundarySingleAzimuthalBlock\n");
#endif // DBGPR_BVALS_CC

  MeshBlock *pmb = pmy_block_;

  if (pmb->loc.level  ==  pmy_mesh_->root_level && pmy_mesh_->nrbx3 == 1
      && pmb->block_size.nx3 > 1) {
    AthenaArray<Real> &var = *var_cc;
    if (pbval_->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar) {
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      for (int n=nl_; n<=nu_; ++n) {
        for (int j=pmb->js-NGHOST; j<=pmb->js-1; ++j) {
          for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i) {
            for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k)
              pbval_->azimuthal_shift_(k) = var(n,k,j,i);
            for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
              int k_shift = k;
              k_shift += (k < (nx3_half + NGHOST) ? 1 : -1) * nx3_half;
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
              k_shift += (k < (nx3_half + NGHOST) ? 1 : -1) * nx3_half;
              var(n,k,j,i) = pbval_->azimuthal_shift_(k_shift);
            }
          }
        }
      }
    }
  }
  return;
}

void CellCenteredBoundaryVariable::SetupPersistentMPI() {
#ifdef MPI_PARALLEL
#ifdef DBGPR_BVALS_CC
  coutYellow("CellCenteredBoundaryVariable::SetupPersistentMPI\n");
#endif // DBGPR_BVALS_CC

  MeshBlock* pmb = pmy_block_;
  int &mylevel = pmb->loc.level;

  int f2 = pmy_mesh_->f2, f3 = pmy_mesh_->f3;
  int cng, cng1, cng2, cng3;
  cng  = cng1 = pmb->cnghost;
  cng2 = cng*f2;
  cng3 = cng*f3;
  int ssize, rsize;
  int tag;
  // Initialize non-polar neighbor communications to other ranks
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (nb.snb.rank != Globals::my_rank) {
      if (nb.snb.level == mylevel) { // same
        ssize = rsize = ((nb.ni.ox1 == 0) ? pmb->block_size.nx1 : NGHOST)
              *((nb.ni.ox2 == 0) ? pmb->block_size.nx2 : NGHOST)
              *((nb.ni.ox3 == 0) ? pmb->block_size.nx3 : NGHOST);
      } else if (nb.snb.level < mylevel) { // coarser
        ssize = ((nb.ni.ox1 == 0) ? ((pmb->block_size.nx1 + 1)/2) : NGHOST)
              *((nb.ni.ox2 == 0) ? ((pmb->block_size.nx2 + 1)/2) : NGHOST)
              *((nb.ni.ox3 == 0) ? ((pmb->block_size.nx3 + 1)/2) : NGHOST);
        rsize = ((nb.ni.ox1 == 0) ? ((pmb->block_size.nx1 + 1)/2 + cng1) : cng1)
              *((nb.ni.ox2 == 0) ? ((pmb->block_size.nx2 + 1)/2 + cng2) : cng2)
              *((nb.ni.ox3 == 0) ? ((pmb->block_size.nx3 + 1)/2 + cng3) : cng3);
      } else { // finer
        ssize = ((nb.ni.ox1 == 0) ? ((pmb->block_size.nx1 + 1)/2 + cng1) : cng1)
              *((nb.ni.ox2 == 0) ? ((pmb->block_size.nx2 + 1)/2 + cng2) : cng2)
              *((nb.ni.ox3 == 0) ? ((pmb->block_size.nx3 + 1)/2 + cng3) : cng3);
        rsize = ((nb.ni.ox1 == 0) ? ((pmb->block_size.nx1 + 1)/2) : NGHOST)
              *((nb.ni.ox2 == 0) ? ((pmb->block_size.nx2 + 1)/2) : NGHOST)
              *((nb.ni.ox3 == 0) ? ((pmb->block_size.nx3 + 1)/2) : NGHOST);
      }
      ssize *= (nu_ + 1); rsize *= (nu_ + 1);
      // specify the offsets in the view point of the target block: flip ox? signs

      // Initialize persistent communication requests attached to specific BoundaryData
      // cell-centered hydro: bd_hydro_
      tag = pbval_->CreateBvalsMPITag(nb.snb.lid, nb.targetid, cc_phys_id_);
      if (bd_var_.req_send[nb.bufid] != MPI_REQUEST_NULL)
        MPI_Request_free(&bd_var_.req_send[nb.bufid]);
      MPI_Send_init(bd_var_.send[nb.bufid], ssize, MPI_ATHENA_REAL,
                    nb.snb.rank, tag, MPI_COMM_WORLD, &(bd_var_.req_send[nb.bufid]));
      tag = pbval_->CreateBvalsMPITag(pmb->lid, nb.bufid, cc_phys_id_);
      if (bd_var_.req_recv[nb.bufid] != MPI_REQUEST_NULL)
        MPI_Request_free(&bd_var_.req_recv[nb.bufid]);
      MPI_Recv_init(bd_var_.recv[nb.bufid], rsize, MPI_ATHENA_REAL,
                    nb.snb.rank, tag, MPI_COMM_WORLD, &(bd_var_.req_recv[nb.bufid]));

      if (FLUID_ENABLED) {
        // hydro flux correction: bd_var_flcor_
        if (pmy_mesh_->multilevel && nb.ni.type == NeighborConnect::face) {
          int size;
          if (nb.fid == 0 || nb.fid == 1)
            size = ((pmb->block_size.nx2 + 1)/2)*((pmb->block_size.nx3 + 1)/2);
          else if (nb.fid == 2 || nb.fid == 3)
            size = ((pmb->block_size.nx1 + 1)/2)*((pmb->block_size.nx3 + 1)/2);
          else // (nb.fid == 4 || nb.fid == 5)
            size = ((pmb->block_size.nx1 + 1)/2)*((pmb->block_size.nx2 + 1)/2);
          size *= (nu_ + 1);
          if (nb.snb.level < mylevel) { // send to coarser
            tag = pbval_->CreateBvalsMPITag(nb.snb.lid, nb.targetid, cc_flx_phys_id_);
            if (bd_var_flcor_.req_send[nb.bufid] != MPI_REQUEST_NULL)
              MPI_Request_free(&bd_var_flcor_.req_send[nb.bufid]);
            MPI_Send_init(bd_var_flcor_.send[nb.bufid], size, MPI_ATHENA_REAL,
                          nb.snb.rank, tag, MPI_COMM_WORLD,
                          &(bd_var_flcor_.req_send[nb.bufid]));
          } else if (nb.snb.level > mylevel) { // receive from finer
            tag = pbval_->CreateBvalsMPITag(pmb->lid, nb.bufid, cc_flx_phys_id_);
            if (bd_var_flcor_.req_recv[nb.bufid] != MPI_REQUEST_NULL)
              MPI_Request_free(&bd_var_flcor_.req_recv[nb.bufid]);
            MPI_Recv_init(bd_var_flcor_.recv[nb.bufid], size, MPI_ATHENA_REAL,
                          nb.snb.rank, tag, MPI_COMM_WORLD,
                          &(bd_var_flcor_.req_recv[nb.bufid]));
          }
        }
      }

    }
  }
#endif
  return;
}

void CellCenteredBoundaryVariable::StartReceiving(BoundaryCommSubset phase) {
#ifdef MPI_PARALLEL
#ifdef DBGPR_BVALS_CC
  coutYellow("CellCenteredBoundaryVariable::StartReceiving\n");
#endif // DBGPR_BVALS_CC

  MeshBlock *pmb = pmy_block_;
  int mylevel = pmb->loc.level;
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (nb.snb.rank != Globals::my_rank) {
      MPI_Start(&(bd_var_.req_recv[nb.bufid]));
      if(FLUID_ENABLED) {
        if (phase == BoundaryCommSubset::all && nb.ni.type == NeighborConnect::face
            && nb.snb.level > mylevel) // opposite condition in ClearBoundary()
          MPI_Start(&(bd_var_flcor_.req_recv[nb.bufid]));
      }
    }

  }
#endif
  return;
}


void CellCenteredBoundaryVariable::ClearBoundary(BoundaryCommSubset phase) {
#ifdef DBGPR_BVALS_CC
  coutYellow("CellCenteredBoundaryVariable::ClearBoundary\n");
#endif // DBGPR_BVALS_CC

  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    bd_var_.flag[nb.bufid] = BoundaryStatus::waiting;
    bd_var_.sflag[nb.bufid] = BoundaryStatus::waiting;

    if (FLUID_ENABLED) {
      if (nb.ni.type == NeighborConnect::face) {
        bd_var_flcor_.flag[nb.bufid] = BoundaryStatus::waiting;
        bd_var_flcor_.sflag[nb.bufid] = BoundaryStatus::waiting;
      }
    }

#ifdef MPI_PARALLEL
    MeshBlock *pmb = pmy_block_;
    int mylevel = pmb->loc.level;
    if (nb.snb.rank != Globals::my_rank) {
      // Wait for Isend
      MPI_Wait(&(bd_var_.req_send[nb.bufid]), MPI_STATUS_IGNORE);

      if (FLUID_ENABLED) {
        if (phase == BoundaryCommSubset::all && nb.ni.type == NeighborConnect::face
            && nb.snb.level < mylevel)
          MPI_Wait(&(bd_var_flcor_.req_send[nb.bufid]), MPI_STATUS_IGNORE);
      }

    }

#endif
  }

  // clear shearing box boundary communications
  if (SHEARING_BOX) {
    // TODO(KGF): clear sflag arrays
    for (int upper=0; upper<2; upper++) {
      if (pbval_->is_shear[upper]) {
        for (int n=0; n<4; n++) {
          if (pbval_->shear_send_neighbor_[upper][n].rank == -1) continue;
          shear_bd_var_[upper].flag[n] = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
          if (pbval_->shear_send_neighbor_[upper][n].rank != Globals::my_rank) {
            MPI_Wait(&shear_bd_var_[upper].req_send[n], MPI_STATUS_IGNORE);
          }
#endif
        }
      }
    }
  }
  return;
}
