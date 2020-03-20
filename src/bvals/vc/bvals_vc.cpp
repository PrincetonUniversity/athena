//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_cc.cpp
//  \brief functions that apply BCs for VERTEX_CENTERED variables

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
#include "bvals_vc.hpp"

// BD: debug
#include "../../lagrange_interp.hpp"
//-

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// constructor

VertexCenteredBoundaryVariable::VertexCenteredBoundaryVariable(
    MeshBlock *pmb, AthenaArray<Real> *var, AthenaArray<Real> *coarse_var,
    AthenaArray<Real> *var_flux)
    : BoundaryVariable(pmb), var_vc(var), coarse_buf(coarse_var), x1flux(var_flux[X1DIR]),
      x2flux(var_flux[X2DIR]), x3flux(var_flux[X3DIR]), nl_(0), nu_(var->GetDim4() -1),
      flip_across_pole_(nullptr) {
  // VertexCenteredBoundaryVariable should only be used w/ 4D or 3D (nx4=1) AthenaArray
  // For now, assume that full span of 4th dim of input AthenaArray should be used:
  // ---> get the index limits directly from the input AthenaArray
  // <=nu_ (inclusive), <nx4 (exclusive)
  if (nu_ < 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in VertexCenteredBoundaryVariable constructor" << std::endl
        << "An 'AthenaArray<Real> *var' of nx4_ = " << var->GetDim4() << " was passed\n"
        << "Should be nx4 >= 1 (likely uninitialized)." << std::endl;
    ATHENA_ERROR(msg);
  }

  InitBoundaryData(bd_var_, BoundaryQuantity::vc);
#ifdef MPI_PARALLEL
  // KGF: dead code, leaving for now:
  // cc_phys_id_ = pbval_->ReserveTagVariableIDs(1);
  vc_phys_id_ = pbval_->bvars_next_phys_id_;
#endif
  if (pmy_mesh_->multilevel) { // SMR or AMR
    // InitBoundaryData(bd_var_flcor_, BoundaryQuantity::cc_flcor);
#ifdef MPI_PARALLEL
    // cc_flx_phys_id_ = cc_phys_id_ + 1;
#endif
  }

// VC
//   if (SHEARING_BOX) {
// #ifdef MPI_PARALLEL
//     shear_cc_phys_id_ = cc_phys_id_ + 2;
// #endif
//     if (pbval_->ShBoxCoord_ == 1) {
//       int nc2 = pmb->ncells2;
//       int nc3 = pmb->ncells3;
//       for (int upper=0; upper<2; upper++) {
//         if (pbval_->is_shear[upper]) {
//           shear_cc_[upper].NewAthenaArray(nu_+1, nc3, nc2, NGHOST);
//           shear_flx_cc_[upper].NewAthenaArray(nc2);

//           // TODO(KGF): the rest of this should be a part of InitBoundaryData()

//           // attach corner cells from L/R side
//           int size = (pmb->block_size.nx2 + NGHOST)*pbval_->ssize_*(nu_ + 1);
//           for (int n=0; n<2; n++) {
//             shear_bd_var_[upper].send[n] = new Real[size];
//             shear_bd_var_[upper].recv[n] = new Real[size];
//             shear_bd_var_[upper].flag[n] = BoundaryStatus::waiting;
// #ifdef MPI_PARALLEL
//             shear_bd_var_[upper].req_send[n] = MPI_REQUEST_NULL;
//             shear_bd_var_[upper].req_recv[n] = MPI_REQUEST_NULL;
// #endif
//           }
//           // corner cells only
//           size = NGHOST*pbval_->ssize_*(nu_ + 1);
//           for (int n=2; n<4; n++) {
//             shear_bd_var_[upper].send[n] = new Real[size];
//             shear_bd_var_[upper].recv[n] = new Real[size];
//             shear_bd_var_[upper].flag[n] = BoundaryStatus::waiting;
// #ifdef MPI_PARALLEL
//             shear_bd_var_[upper].req_send[n] = MPI_REQUEST_NULL;
//             shear_bd_var_[upper].req_recv[n] = MPI_REQUEST_NULL;
// #endif
//           }
//         } // end "if is a shearing boundary"
//       }  // end loop over inner, outer shearing boundaries
//     } // end "if (pbval_->ShBoxCoord_ == 1)"
//   } // end shearing box component of ctor
}

// destructor

VertexCenteredBoundaryVariable::~VertexCenteredBoundaryVariable() {
  DestroyBoundaryData(bd_var_);
  if (pmy_mesh_->multilevel)
    DestroyBoundaryData(bd_var_flcor_);

  // // TODO(KGF): this should be a part of DestroyBoundaryData()
  // if (SHEARING_BOX) {
  //   for (int upper=0; upper<2; upper++) {
  //     if (pbval_->is_shear[upper]) { // if true for shearing inner blocks
  //       for (int n=0; n<4; n++) {
  //         delete[] shear_bd_var_[upper].send[n];
  //         delete[] shear_bd_var_[upper].recv[n];
  //       }
  //     }
  //   }
  // }
}

void VertexCenteredBoundaryVariable::ErrorIfPolarNotImplemented(
  const NeighborBlock& nb) {

  // BD: TODO implement polar coordinates
  if (nb.polar) {
  std::stringstream msg;
  msg << "### FATAL ERROR" << std::endl
      << "Polar coordinates not implemented for vertex-centered." << std::endl;
  ATHENA_ERROR(msg);
  }
  return;
}

void VertexCenteredBoundaryVariable::ErrorIfShearingBoxNotImplemented() {
  // BD: TODO implement shearing box
  if (SHEARING_BOX){
    std::stringstream msg;
    msg << "### FATAL ERROR" << std::endl
        << "Shearing box not implemented for vertex-centered." << std::endl;
    ATHENA_ERROR(msg);
  }
}

int VertexCenteredBoundaryVariable::ComputeVariableBufferSize(const NeighborIndexes& ni,
                                                              int cng) {
  coutYellow("VertexCenteredBoundaryVariable::ComputeVariableBufferSize\n");
  MeshBlock *pmb = pmy_block_;
  int cng1, cng2, cng3;
  cng1 = cng;
  cng2 = cng*(pmb->block_size.nx2 > 1 ? 1 : 0);
  cng3 = cng*(pmb->block_size.nx3 > 1 ? 1 : 0);

  // BD: TODO- CALCULATE CAREFULLY
  int size = ((ni.ox1 == 0) ? pmb->nverts1 : NGHOST)
    *((ni.ox2 == 0) ? pmb->nverts2 : NGHOST)
    *((ni.ox3 == 0) ? pmb->nverts3 : NGHOST);

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

  printf("size = %d\n", size);

  // DOUBLE FOR DEBUG
  return 4 * size;
}

//----------------------------------------------------------------------------------------
//! \fn int VertexCenteredBoundaryVariable::LoadBoundaryBufferSameLevel(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set vertex-centered boundary buffers for sending to a block on the same level

int VertexCenteredBoundaryVariable::LoadBoundaryBufferSameLevel(Real *buf,
                                                                const NeighborBlock& nb) {
  coutYellow("VertexCenteredBoundaryVariable::LoadBoundaryBufferSameLevel\n");
  nb.print_all();

  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;
  int p = 0;
  AthenaArray<Real> &var = *var_vc;

  // shared vertex is packed
  si = (nb.ni.ox1 > 0) ? pmb->ige : pmb->ivs;
  ei = (nb.ni.ox1 < 0) ? pmb->igs : pmb->ive;

  sj = (nb.ni.ox2 > 0) ? pmb->jge : pmb->jvs;
  ej = (nb.ni.ox2 < 0) ? pmb->jgs : pmb->jve;

  sk = (nb.ni.ox3 > 0) ? pmb->kge : pmb->kvs;
  ek = (nb.ni.ox3 < 0) ? pmb->kgs : pmb->kve;

  BufferUtility::PackData(var, buf, nl_, nu_, si, ei, sj, ej, sk, ek, p);

  // if ((nb.ni.ox1 == -1) and (nb.ni.ox2 == -1))
  //   Q();

  // BD: debug: lagrange interp
  /*
  if (false) {
    coutBoldBlue("\nlagrange_interp scratch\n");

    Coordinates *pco = pmb->pcoord;

    // settings for interpolator
    const int interp_order = 4;
    const int interp_dim = 1;

    Real origin[1] = {pco->x1f(0)};
    Real delta[1] = {pco->dx1v(0)};
    int size[1] = {pmb->block_size.nx1 + 2*(NGHOST) + 1};

    // target coord
    Real coord[1] = {pco->x1f(2) + delta[0] / 2.};
    printf("coord: %1.18f\n", coord[0]);

    // set up n-dimensional interpolator
    LagrangeInterpND<interp_order, interp_dim> * pinterp = nullptr;

    pinterp = new LagrangeInterpND<interp_order, interp_dim>(origin, delta, size, coord);

    // test interpolation to point
    AthenaArray<Real> mySrc, myDst;

    int iv = 0; // solution component
    mySrc.InitWithShallowSlice(const_cast<AthenaArray<Real>&>(var), iv, 1);
    // myDst.InitWithShallowSlice(const_cast<AthenaArray<Real>&>(var), iv, 1);

    myDst.NewAthenaArray(var.GetDim1());

    // attempt interp.
    mySrc.Fill(200);
    mySrc(0, 0, 0, 0) = 10;
    mySrc(0, 0, 0, 1) = 11;
    mySrc(0, 0, 0, 2) = 12;
    myDst(0, 0, 0, 0) = pinterp->eval(mySrc.data());

    // inspect
    coutBoldBlue("mySrc:\n");
    mySrc.print_all();

    coutBoldBlue("myDst:\n");
    myDst.print_all();


    // clean-up
    // for (int i = 0; i < n; ++i) {
    //   delete pinterp[i];
    // }
    delete pinterp;
    //--

    Q();
  }
  */
  return p;
}


//----------------------------------------------------------------------------------------
//! \fn int VertexCenteredBoundaryVariable::LoadBoundaryBufferToCoarser(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set vertex-centered boundary buffers for sending to a block on the coarser level

int VertexCenteredBoundaryVariable::LoadBoundaryBufferToCoarser(Real *buf,
                                                                const NeighborBlock& nb) {
  coutYellow("VertexCenteredBoundaryVariable::LoadBoundaryBufferToCoarser\n");
  nb.print_all();

  MeshBlock *pmb = pmy_block_;
  MeshRefinement *pmr = pmb->pmr;
  int si, sj, sk, ei, ej, ek;

  AthenaArray<Real> &var = *var_vc;
  AthenaArray<Real> &coarse_var = *coarse_buf;

  // vertices that are shared with adjacent MeshBlocks are to be copied to coarser level
  int ng = NGHOST;
  si = (nb.ni.ox1 > 0) ? (pmb->civ - ng) : pmb->cis;
  ei = (nb.ni.ox1 < 0) ? (pmb->cis + ng) : pmb->civ;

  sj = (nb.ni.ox2 > 0) ? (pmb->cjv - ng) : pmb->cjs;
  ej = (nb.ni.ox2 < 0) ? (pmb->cjs + ng) : pmb->cjv;

  sk = (nb.ni.ox3 > 0) ? (pmb->ckv - ng) : pmb->cks;
  ek = (nb.ni.ox3 < 0) ? (pmb->cks + ng) : pmb->ckv;


  int p = 0;

  coutBoldRed("var_vc\n");
  var.print_all();

  coutBoldRed("coarse_buf\n");
  coarse_var.print_all();

  pmr->RestrictVertexCenteredValues(var, coarse_var, nl_, nu_, si, ei, sj, ej, sk, ek);
  coutBoldRed("coarse_buf, buf");
  BufferUtility::PackData(coarse_var, buf, nl_, nu_, si, ei, sj, ej, sk, ek, p);

  return p;
}

//----------------------------------------------------------------------------------------
//! \fn int VertexCenteredBoundaryVariable::LoadBoundaryBufferToFiner(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set vertex-centered boundary buffers for sending to a block on the finer level

int VertexCenteredBoundaryVariable::LoadBoundaryBufferToFiner(Real *buf,
                                                              const NeighborBlock& nb) {
  coutYellow("VertexCenteredBoundaryVariable::LoadBoundaryBufferToFiner\n");
  nb.print_all();

  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;
  int cn = pmb->cnghost - 1;
  AthenaArray<Real> &var = *var_vc;

  // si = (nb.ni.ox1 > 0) ? (pmb->ie - cn) : pmb->is;
  // ei = (nb.ni.ox1 < 0) ? (pmb->is + cn) : pmb->ie;
  // sj = (nb.ni.ox2 > 0) ? (pmb->je - cn) : pmb->js;
  // ej = (nb.ni.ox2 < 0) ? (pmb->js + cn) : pmb->je;
  // sk = (nb.ni.ox3 > 0) ? (pmb->ke - cn) : pmb->ks;
  // ek = (nb.ni.ox3 < 0) ? (pmb->ks + cn) : pmb->ke;

  si = (nb.ni.ox1 > 0) ? (pmb->iv - 1 - cn) : pmb->is;
  ei = (nb.ni.ox1 < 0) ? (pmb->is + 1 + cn) : pmb->iv;

  sj = (nb.ni.ox2 > 0) ? (pmb->je - 1 - cn) : pmb->js;
  ej = (nb.ni.ox2 < 0) ? (pmb->js + 1 + cn) : pmb->jv;

  sk = (nb.ni.ox3 > 0) ? (pmb->ke - 1 - cn) : pmb->ks;
  ek = (nb.ni.ox3 < 0) ? (pmb->ks + 1 + cn) : pmb->kv;
  // send the data first and later prolongate on the target block

  // need to add edges for faces, add corners for edges
  if (nb.ni.ox1 == 0) {
    if (nb.ni.fi1 == 1)   si += pmb->block_size.nx1/2;
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

  // modify for vc [this needs to be checked for 2/3d]
  // int io = 1; // offset for x1 to avoid edge vertices

  // if (nb.ni.ox1 > 0) {
  //   si = pmb->iv - cn - io;
  //   ei = pmb->iv - io;
  // } else {
  //   si = pmb->is + io;
  //   ei = pmb->is + cn + io;
  // }

  si = (nb.ni.ox1 > 0) ? (pmb->iv - 1 - cn) : pmb->is;
  ei = (nb.ni.ox1 < 0) ? (pmb->is + 1 + cn) : pmb->iv;

  sj = (nb.ni.ox2 > 0) ? (pmb->je - 1 - cn) : pmb->js;
  ej = (nb.ni.ox2 < 0) ? (pmb->js + 1 + cn) : pmb->jv;

  sk = (nb.ni.ox3 > 0) ? (pmb->ke - 1 - cn) : pmb->ks;
  ek = (nb.ni.ox3 < 0) ? (pmb->ks + 1 + cn) : pmb->kv;

  int p = 0;

  coutBoldRed("var_vc, buf");
  BufferUtility::PackData(var, buf, nl_, nu_, si, ei, sj, ej, sk, ek, p);

  // printf("x1f: ");
  // pmb->pcoord->x1f.print_data("%1.2f");

  coutBoldRed("MB::UWIL gid = ");
  printf("%d\n", pmb->gid);

  // if (nb.ni.ox1 > 0)
  //   Q();
  // if ((nb.ni.ox1 == 1) and (nb.ni.ox2 == -1))
  //   Q();

  return p;
}


//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::SetBoundarySameLevel(Real *buf,
//                                                              const NeighborBlock& nb)
//  \brief Set vertex-centered boundary received from a block on the same level

void VertexCenteredBoundaryVariable::SetBoundarySameLevel(Real *buf,
                                                          const NeighborBlock& nb) {

  coutYellow("VertexCenteredBoundaryVariable::SetBoundarySameLevel\n");
  nb.print_all();

  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;
  AthenaArray<Real> &var = *var_vc;
  int p = 0;

  si = ei = sj = ej = sk = ek = 0;

  // BD: TODO implement
  ErrorIfPolarNotImplemented(nb);
  ErrorIfShearingBoxNotImplemented();

  //----------------------------------------------------------------------------
  // could unpack based on dimensionality
  // shared data unpacked with add, isolated with overwrite
  // if (pmb->block_size.nx3 > 1) {
  //   // ...
  // } else if (pmb->block_size.nx2 > 1) {
  //   // ...
  // } else {
  //   // shared vertex
  //   if (nb.ni.ox1 < 0) {   // neighbor is to the left
  //     si = ei = pmb->ive;
  //   } else {            // neighbor is to the right
  //     si = ei = pmb->ivs;
  //   }
  //   BufferUtility::UnpackDataAdd(buf, var, nl_, nu_,
  //                                si, ei, sj, ej, sk, ek, p);
  //   // isolated vertices
  //   if (nb.ni.ox1 < 0) {
  //     si = pmb->ips; ei = pmb->ipe;
  //   } else {
  //     si = pmb->ims; ei = pmb->ime;
  //   }
  //   BufferUtility::UnpackData(buf, var, nl_, nu_, si, ei, sj, ej, sk, ek, p);
  // }
  //----------------------------------------------------------------------------


  // unpack all data additively
  // defer imposition (via suitable averaging) of consistency condition
  if (nb.ni.ox1 == 0) {
    si = pmb->ivs; ei = pmb->ive;
  } else if (nb.ni.ox1 > 0) {
    si = pmb->ive; ei = pmb->ipe;
  } else {
    si = pmb->ims; ei = pmb->ivs;
  }

  if (nb.ni.ox2 == 0) {
    sj = pmb->jvs; ej = pmb->jve;
  } else if (nb.ni.ox2 > 0) {
    sj = pmb->jve; ej = pmb->jpe;
  } else {
    sj = pmb->jms; ej = pmb->jvs;
  }

  if (nb.ni.ox3 == 0) {
    sk = pmb->kvs; ek = pmb->kve;
  } else if (nb.ni.ox3 > 0) {
    sk = pmb->kve; ek = pmb->kpe;
  } else {
    sk = pmb->kms; ek = pmb->kvs;
  }

  //////////////////////////////////////////////////////////////////////////////
  // BD: debug - populate based on solution
  if (FILL_WAVE_BND_SL) {
    // populate only faces [only one oxi is non-zero]
    // if (std::abs(nb.ni.ox1) + std::abs(nb.ni.ox2) + std::abs(nb.ni.ox3) == 1) {
    //pmb->DebugWaveMeshBlock(var, si, ei, sj, ej, sk, ek, true);
    // }

    // // populate only edges [only two oxi are non-zero]
    // if (std::abs(nb.ni.ox1) + std::abs(nb.ni.ox2) + std::abs(nb.ni.ox3) == 2) {
    //   //
    //   //pmb->DebugWaveMeshBlock(var, si, ei, sj, ej, sk, ek, true);
    // }

    // // populate only corners [oxi != 0]
    // if (std::abs(nb.ni.ox1) + std::abs(nb.ni.ox2) + std::abs(nb.ni.ox3) == 3) {
    // }

    pmb->DebugWaveMeshBlock(var, si, ei, sj, ej, sk, ek, true);
    var.print_all();
    // Q();
  } else {
    BufferUtility::UnpackDataAdd(buf, var, nl_, nu_, si, ei, sj, ej, sk, ek, p);
    // BufferUtility::UnpackData(buf, var, nl_, nu_, si, ei, sj, ej, sk, ek, p);
  }
  //////////////////////////////////////////////////////////////////////////////

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::SetBoundaryFromCoarser(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set vertex-centered prolongation buffer received from a block on a coarser level

void VertexCenteredBoundaryVariable::SetBoundaryFromCoarser(Real *buf,
                                                            const NeighborBlock& nb) {
  // populating from a coarser level; do not touch shared vertices
  coutYellow("VertexCenteredBoundaryVariable::SetBoundaryFromCoarser\n");
  nb.print_all();

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

  // modify for vc [this needs to be checked for 2/3d]
  // if (nb.ni.ox1 > 0)  {
  //   si = pmb->civ + 1,   ei = pmb->civ + cng;
  // } else {
  //   si = pmb->cis - cng, ei = pmb->cis - 1;
  // }

  int p = 0;
  // BD: TODO implement
  // if (nb.polar) {
  ErrorIfPolarNotImplemented(nb);
  // } else {
  //   BufferUtility::UnpackData(buf, coarse_var, nl_, nu_,
  //                             si, ei, sj, ej, sk, ek, p);
  // }

  // AthenaArray<Real> &var = *var_vc;

  // coarse_var.print_all();
  // var.print_all();
  // Q();

  // write like this to regroup idx logic for debug
  if (!FILL_WAVE_BND_FRC)
    BufferUtility::UnpackData(buf, coarse_var, nl_, nu_,
                              si, ei, sj, ej, sk, ek, p);

  //////////////////////////////////////////////////////////////////////////////
  // BD: debug - populate based on solution
  // [note we directly modify fund. not coarse]

  bool flag = false;

  // [1d] modifies si, ei
  if (nb.ni.ox1 == 0) {
    //
    // (si, ei)= (cis, cie)
    if ((pmb->loc.lx1 & 1LL) == 0LL) {
      //
      // (si, ei) = (cis, cie + cnghost)
      si = pmb->ivs;
      ei = pmb->ipe;
    } else {
      //
      // (si, ei) = (cis - cnghost, cie)
      si = pmb->ims;
      ei = pmb->ive;
    }
  } else if (nb.ni.ox1 > 0)  {
    //
    // (si, ei) = (cie + 1, cie + cnghost)
    si = pmb->ips;
    ei = pmb->ipe;
  } else {
    //
    // (si, ei) = (cis - cnghost, cis - 1)
    si = pmb->ims;
    ei = pmb->ime;
  }

  // [2d] modifies sj, ej
  if (nb.ni.ox2 == 0) {
    //
    // (sj, ej) = (cjs, cje)
    if (pmb->block_size.nx2 > 1) {
      if ((pmb->loc.lx2 & 1LL) == 0LL) {
        //
        // (sj, ej) = (cjs, cje + cnghost)
        sj = pmb->jvs;
        ej = pmb->jpe;
      } else {
        //
        // (sj, ej) = (cjs - cnghost, cje)
        sj = pmb->jms;
        ej = pmb->jve;
      }
    }
  } else if (nb.ni.ox2 > 0) {
    //
    // (sj, ej) = (cje + 1, cje + cnghost)
    sj = pmb->jps;
    ej = pmb->jpe;
  } else {
    //
    // (sj, ej) = (cjs - cnghost, cjs - 1)
    sj = pmb->jms;
    ej = pmb->jme;
  }

  // [3d] modifies sk, ek
  if (nb.ni.ox3 == 0) {
    //
    // (sk, ek) = (cks, cke)
    if (pmb->block_size.nx3 > 1) {
      if ((pmb->loc.lx3 & 1LL) == 0LL) {
        //
        // (sk, ek) = (cks, cke + cnghost)
        sk = pmb->kvs;
        ek = pmb->kpe;
      } else {
        //
        // (sk, ek) = (cks - cnghost, cke)
        sk = pmb->kms;
        ek = pmb->kve;
      }
    }
  } else if (nb.ni.ox3 > 0)  {
    //
    // (sk, ek) = (cke + 1, cke + cnghost)
    sk = pmb->kps;
    ek = pmb->kpe;
  } else {
    //
    // (sk, ek) = (cks - cnghost, cks - 1)
    sk = pmb->kms;
    ek = pmb->kme;
  }
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  // BD: debug - populate based on solution
  // [note we directly modify fund. not coarse]
  if (FILL_WAVE_BND_FRC) {
    // note that otherwise prolongation from ghost zones is performed
    AthenaArray<Real> &var = *var_vc;
    var.print_all();
    pmb->DebugWaveMeshBlock(var, si, ei, sj, ej, sk, ek, false);
    var.print_all();

    if (flag)
      Q();
  }
  //////////////////////////////////////////////////////////////////////////////

  coutBoldRed("MB::UWIL gid = ");
  printf("%d\n", pmb->gid);

  // if (nb.ni.ox1 < 0)
  //   Q();

  // printf("x1f: ");
  // pmb->pcoord->x1f.print_data("%1.2f");

  // if (nb.ni.ox1 == -1)
  //   Q();
  // if ((nb.ni.ox1 == 1) and (nb.ni.ox2 == -1))
  //   Q();

  // pmb->DebugWaveMeshBlockSolution();
  // Q();

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::SetBoundaryFromFiner(Real *buf,
//                                                              const NeighborBlock& nb)
//  \brief Set vertex-centered boundary received from a block on a finer level

void VertexCenteredBoundaryVariable::SetBoundaryFromFiner(Real *buf,
                                                          const NeighborBlock& nb) {
  // populating from finer level; shared vertices are corrected
  coutYellow("VertexCenteredBoundaryVariable::SetBoundaryFromFiner\n");
  nb.print_all();

  MeshBlock *pmb = pmy_block_;
  AthenaArray<Real> &var = *var_vc;
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


  // modify for vc [this needs to be checked for 2/3d]
  // if (nb.ni.ox1 > 0) {
  //   si = pmb->iv + 1,      ei = pmb->iv + NGHOST;
  // } else {
  //   si = pmb->is - NGHOST, ei = pmb->is - 1;
  // }

  int p = 0;

  // BD: TODO implement
  // if (nb.polar) {
  ErrorIfPolarNotImplemented(nb);
  // } else {
  //   BufferUtility::UnpackData(buf, var, nl_, nu_, si, ei, sj, ej, sk, ek, p);
  // }

  /*
  if (nb.ni.ox1 == 0) {
    si = pmb->is, ei = pmb->iv;
    // how to treat this? [overwrite shared vertex for now]
    if (nb.ni.fi1 == 1)   si += pmb->block_size.nx1/2;  // +1 to avoid vert;
    else            ei -= pmb->block_size.nx1/2;        // +1 to avoid vert
  } else if (nb.ni.ox1 > 0) {
    si = pmb->iv,      ei = pmb->iv + NGHOST;
  } else {
    si = pmb->is - NGHOST, ei = pmb->is;
  }

  if (nb.ni.ox2 == 0) {
    sj = pmb->js, ej = pmb->jv;
    // if (pmb->block_size.nx2 > 1) {
    //   if (nb.ni.ox1 != 0) {
    //     if (nb.ni.fi1 == 1) sj += pmb->block_size.nx2/2;
    //     else          ej -= pmb->block_size.nx2/2;
    //   } else {
    //     if (nb.ni.fi2 == 1) sj += pmb->block_size.nx2/2;
    //     else          ej -= pmb->block_size.nx2/2;
    //   }
    // }
  } else if (nb.ni.ox2 > 0) {
    sj = pmb->jv,      ej = pmb->jv + NGHOST;
  } else {
    sj = pmb->js - 82HOST, ej = pmb->js;
  }
  2/
  // DEBUG: prefill ghosts + shared vertices
  /*
  if (false)
    if (pmb->block_size.nx3 > 1) {
      //...
    } else if (3m2->+0oc-_size.nx2 > 1) {
      for (int j_ix=0; j_ix<=NGHOST; j_ix++){
        pco->x1f(2)0- i_ix<=NGHOST; i_ix++ / 2.){
          var(j_ix, i_ix) = 3;
          var(j_ix,
              interp_order+pmb->block_size.nx1 - i_ix) = 6;
          var(2*NGHOST+pmb->block_size.nx2 - j_ix, i_ix) = 8;
          var(2*NGHOST+pmb->blockinterp_ordersize.nx2 - j_ix,
              2*NGHOST+pmb->block_size.nx1 - i_ix) = 9;
        }
      }
    } else {
      for (int i_ix=0; i_ix<=NGHOST; i_ix++){
        var(i_ix) = 3;
        var(2*NGHOST+pmb->block_size.nx1 - i_ix) = 3;
      }
    }
  */

  bool flag = false;

  // [1d] modifies si, ei
  if (nb.ni.ox1 == 0) {
    if (nb.ni.fi1 == 1) {
      // (ox1=0, ~(ox2=0 /\ ox3=0), fi1=1, fi2)
      //
      // (si, ei) = (is + block_size.nx1 / 2, ie)
      si = pmb->ivs + pmb->block_size.nx1 / 2;
      ei = pmb->ive;
      printf("tag1\n");
    } else {
      // (ox1=0, ~(ox2=0 /\ ox3=0), fi1=0, fi2)
      //
      // (si, ei) = (is, ie - block_size.nx1 / 2)
      si = pmb->ivs;
      ei = pmb->ive - pmb->block_size.nx1 / 2;
      printf("tag2\n");
    }
  } else if (nb.ni.ox1 > 0) {
    // (ox1>0, ox2, ox3, fi1, fi2)
    //
    // (si, ei) = (ie + 1, ie + NGHOST)
    si = pmb->ive;
    ei = pmb->ipe;
    printf("tag3\n");
  } else {
    // (ox1<0, ox2, ox3, fi1, fi2)
    //
    // (si, ei) = (is - NGHOST, is - 1)
    si = pmb->ims;
    ei = pmb->ivs;
    printf("tag4\n");
  }

  // [2d] modifies sj, ej
  if (nb.ni.ox2 == 0 && pmb->block_size.nx2 > 1) {
    if (nb.ni.ox1 != 0) {
      if (nb.ni.fi1 == 1) {
        // (ox1!=0, ox2=0, ox3, fi1=1, fi2)
        //
        // (sj, ej) = (js + block_size.nx2 / 2, je)
        sj = pmb->jvs + pmb->block_size.nx2 / 2;
        ej = pmb->jve;
        printf("2tag1\n");
      } else {
        // (ox1!=0, ox2=0, ox3, fi1=0, fi2)
        //
        // (sj, ej) = (js, je - block_size.nx2 / 2)
        sj = pmb->jvs;
        ej = pmb->jve - pmb->block_size.nx2 / 2;
        printf("2tag2\n");
      }
    } else {
      if (nb.ni.fi2 == 1) {
        // (ox1=0, ox2=0, ox3!=0, fi1, fi2=1)
        //
        // (sj, ej) = (js + block_size.nx2 / 2, je)
        sj = pmb->jvs + pmb->block_size.nx2 / 2;
        ej = pmb->jve;
        printf("2tag3\n");
      } else {
        // (ox1=0, ox2=0, ox3!=0, fi1, fi2=0)
        //
        // (sj, ej) = (js, je - block_size.nx2 / 2)
        sj = pmb->jvs;
        ej = pmb->jve - pmb->block_size.nx2 / 2;
        printf("2tag4\n");
      }
    }
  } else if (nb.ni.ox2 > 0) {
    // (ox1, ox2>0, ox3, fi1, fi2)
    //
    // (sj, ej) = (je + 1, je + NGHOST)
    sj = pmb->jve;
    ej = pmb->jpe;
    printf("2tag5\n");
  } else {
    // (ox1, ox2<0, ox3, fi1, fi2)
    //
    // (sj, ej) = (js - NGHOST, js - 1)
    sj = pmb->jms;
    ej = pmb->jvs;
    printf("2tag6\n");
  }

  // [3d] modifies sj, ek
  if (nb.ni.ox3 == 0 && pmb->block_size.nx3 > 1) {
    if (nb.ni.ox1 != 0 && nb.ni.ox2 != 0) {
      if (nb.ni.fi1 == 1) {
        // (ox1!=0, ox2!=0, ox3=0, fi1=1, fi2)
        //
        // (sk, ek) = (ks + block_size.nx3 / 2, ke)
        sk = pmb->kvs + pmb->block_size.nx3 / 2;
        ek = pmb->kve;
        printf("3tag1\n");
      } else {
        // (ox1!=0, ox2!=0, ox3=0, fi1=0, fi2)
        //
        // (sk, ek) = (ks, ke - block_size.nx3 / 2)
        sk = pmb->kvs;
        ek = pmb->kve - pmb->block_size.nx3 / 2;
        printf("3tag2\n");
      }
    } else {
      if (nb.ni.fi2 == 1) {
        // (~(ox1!=0 /\ ox2!=0), ox3=0, fi1, fi2=1)
        //
        // (sk, ek) = (ks + block_size.nx3 / 2, ke)
        sk = pmb->kvs + pmb->block_size.nx3 / 2;
        ek = pmb->kve;
        printf("3tag3\n");
      } else {
        // (~(ox1!=0 /\ ox2!=0), ox3=0, fi1, fi2=0)
        //
        // (sk, ek) = (ks, ke - block_size.nx3 / 2)
        sk = pmb->kvs;
        ek = pmb->kve - pmb->block_size.nx3 / 2;
        printf("3tag4\n");
      }
    }
  } else if (nb.ni.ox3 > 0) {
    // (ox1, ox2, ox3>0, fi1, fi2)
    //
    // (sk, ek) = (ke + 1, ke + NGHOST)
    sk = pmb->kve;
    ek = pmb->kpe;
    printf("3tag5\n");
  } else {
    // (ox1, ox2, ox3<0, fi1, fi2)
    //
    // (sk, ek) = (ks - NGHOST, ks - 1)
    sk = pmb->kms;
    ek = pmb->kvs;
    printf("3tag6\n");
  }
  //////////////////////////////////////////////////////////////////////////////



  if (FILL_WAVE_BND_FRF) {
    var.print_all();
    //var.ZeroClear();
    pmb->DebugWaveMeshBlock(var, si, ei, sj, ej, sk, ek, true);
    var.print_all();
    if (flag)
      Q();
  } else {
    coutBoldRed("buf, var_vc");
    BufferUtility::UnpackDataAdd(buf, var, nl_, nu_, si, ei, sj, ej, sk, ek, p);
  }

  // if (nb.ni.ox1 == +1)
  //   Q();

  // if ((nb.ni.ox1 == -1) and (nb.ni.ox2 == 0))
  //   Q();

  // printf("x1f: ");
  // pmb->pcoord->x1f.print_data("%1.2f");

  return;
}



//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::SetBoundaries()
//  \brief set the vertex-centered boundary data

void VertexCenteredBoundaryVariable::SetBoundaries() {
  coutYellow("VertexCenteredBoundaryVariable::SetBoundaries\n");

  ZeroVertexGhosts();
  BoundaryVariable::SetBoundaries();
  FinalizeVertexConsistency();
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::ReceiveAndSetBoundariesWithWait()
//  \brief receive and set the vertex-centered boundary data for initialization

void VertexCenteredBoundaryVariable::ReceiveAndSetBoundariesWithWait() {
  coutYellow("VertexCenteredBoundaryVariable::ReceiveAndSetBoundariesWithWait\n");

  ZeroVertexGhosts();
  BoundaryVariable::ReceiveAndSetBoundariesWithWait();
  FinalizeVertexConsistency();
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::PolarBoundarySingleAzimuthalBlock()
// \brief polar boundary edge-case: single MeshBlock spans the entire azimuthal (x3) range

void VertexCenteredBoundaryVariable::PolarBoundarySingleAzimuthalBlock() {
  coutYellow("VertexCenteredBoundaryVariable::PolarBoundarySingleAzimuthalBlock\n");
  return;
}

void VertexCenteredBoundaryVariable::SetupPersistentMPI() {
#ifdef MPI_PARALLEL
  coutYellow("VertexCenteredBoundaryVariable::SetupPersistentMPI\n");
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
      // vertex-centered hydro: bd_hydro_
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
          // if (bd_var_flcor_.req_send[nb.bufid] != MPI_REQUEST_NULL)
          //   MPI_Request_free(&bd_var_flcor_.req_send[nb.bufid]);
          // MPI_Send_init(bd_var_flcor_.send[nb.bufid], size, MPI_ATHENA_REAL,
          //               nb.snb.rank, tag, MPI_COMM_WORLD,
          //               &(bd_var_flcor_.req_send[nb.bufid]));
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
#endif
  return;
}

void VertexCenteredBoundaryVariable::StartReceiving(BoundaryCommSubset phase) {
#ifdef MPI_PARALLEL
  coutYellow("VertexCenteredBoundaryVariable::StartReceiving\n");

  MeshBlock *pmb = pmy_block_;
  int mylevel = pmb->loc.level;
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (nb.snb.rank != Globals::my_rank) {
      MPI_Start(&(bd_var_.req_recv[nb.bufid]));
      // if (phase == BoundaryCommSubset::all && nb.ni.type == NeighborConnect::face
      //     && nb.snb.level > mylevel) // opposite condition in ClearBoundary()
      //   MPI_Start(&(bd_var_flcor_.req_recv[nb.bufid]));
    }
  }
#endif
  return;
}


void VertexCenteredBoundaryVariable::ClearBoundary(BoundaryCommSubset phase) {
  coutYellow("VertexCenteredBoundaryVariable::ClearBoundary\n");
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    bd_var_.flag[nb.bufid] = BoundaryStatus::waiting;
    bd_var_.sflag[nb.bufid] = BoundaryStatus::waiting;

    // if (nb.ni.type == NeighborConnect::face) {
    //   bd_var_flcor_.flag[nb.bufid] = BoundaryStatus::waiting;
    //   bd_var_flcor_.sflag[nb.bufid] = BoundaryStatus::waiting;
    // }
#ifdef MPI_PARALLEL
    MeshBlock *pmb = pmy_block_;
    int mylevel = pmb->loc.level;
    if (nb.snb.rank != Globals::my_rank) {
      // Wait for Isend
      MPI_Wait(&(bd_var_.req_send[nb.bufid]), MPI_STATUS_IGNORE);
      // if (phase == BoundaryCommSubset::all && nb.ni.type == NeighborConnect::face
      //     && nb.snb.level < mylevel)
      //   MPI_Wait(&(bd_var_flcor_.req_send[nb.bufid]), MPI_STATUS_IGNORE);
    }
#endif
  }

// VC
//   // clear shearing box boundary communications
//   if (SHEARING_BOX) {
//     // TODO(KGF): clear sflag arrays
//     for (int upper=0; upper<2; upper++) {
//       if (pbval_->is_shear[upper]) {
//         for (int n=0; n<4; n++) {
//           if (pbval_->shear_send_neighbor_[upper][n].rank == -1) continue;
//           shear_bd_var_[upper].flag[n] = BoundaryStatus::waiting;
// #ifdef MPI_PARALLEL
//           if (pbval_->shear_send_neighbor_[upper][n].rank != Globals::my_rank) {
//             MPI_Wait(&shear_bd_var_[upper].req_send[n], MPI_STATUS_IGNORE);
//           }
// #endif
//         }
//       }
//     }
//   }
  return;
}
