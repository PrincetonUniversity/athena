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

  // node multiplicities-------------------------------------------------------
  AllocateNodeMult();
  //---------------------------------------------------------------------------

}

// destructor
VertexCenteredBoundaryVariable::~VertexCenteredBoundaryVariable() {
  DestroyBoundaryData(bd_var_);
 // if (pmy_mesh_->multilevel)
 //   DestroyBoundaryData(bd_var_flcor_);

  // node multiplicities-------------------------------------------------------
  node_mult.DeleteAthenaArray();
  //---------------------------------------------------------------------------
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
  // 'cng' to preserve function signature but is a dummy slot
  return NeighborVariableBufferSize(ni);
}

//----------------------------------------------------------------------------------------
//! \fn int VertexCenteredBoundaryVariable::LoadBoundaryBufferSameLevel(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set vertex-centered boundary buffers for sending to a block on the same level

int VertexCenteredBoundaryVariable::LoadBoundaryBufferSameLevel(Real *buf,
                                                                const NeighborBlock& nb) {
  if (DBGPR_BVALS_VC) {
    coutYellow("VertexCenteredBoundaryVariable::LoadBoundaryBufferSameLevel\n");
    nb.print_all();
  }

  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;
  int p = 0;
  AthenaArray<Real> &var = *var_vc;

  idxLoadSameLevelRanges(nb.ni, si, ei, sj, ej, sk, ek, false);
  BufferUtility::PackData(var, buf, nl_, nu_, si, ei, sj, ej, sk, ek, p);

  // if multilevel make use of pre-restricted internal data
  if (pmy_mesh_->multilevel) {
    if (DBGPR_BVALS_VC)
      coutBoldRed("Packing coarse same level\n");

    // convert to coarse indices
    AthenaArray<Real> &coarse_var = *coarse_buf;

    // si = (nb.ni.ox1 > 0) ? pmb->cige : pmb->civs;
    // ei = (nb.ni.ox1 < 0) ? pmb->cigs : pmb->cive;

    // sj = (nb.ni.ox2 > 0) ? pmb->cjge : pmb->cjvs;
    // ej = (nb.ni.ox2 < 0) ? pmb->cjgs : pmb->cjve;

    // sk = (nb.ni.ox3 > 0) ? pmb->ckge : pmb->ckvs;
    // ek = (nb.ni.ox3 < 0) ? pmb->ckgs : pmb->ckve;

    // for partial restrict / partial pack
    // this would still miss a "strip"

    // if (nb.ni.ox1 == 0) {
    //   si = pmb->civs; ei = pmb->cive;
    // } else if (nb.ni.ox1 > 0) {
    //   si = pmb->cige; ei = pmb->cive-pmb->rcng-1;
    // } else {
    //   si = pmb->civs+pmb->rcng+1; ei = pmb->civs+pmb->cng;
    // }

    // if (nb.ni.ox2 == 0) {
    //   sj = pmb->cjvs; ej = pmb->cjve;
    // } else if (nb.ni.ox2 > 0) {
    //   sj = pmb->cjge; ej = pmb->cjve-pmb->rcng-1;
    // } else {
    //   sj = pmb->cjvs+pmb->rcng+1; ej = pmb->cjvs+pmb->cng;
    // }

    // if (nb.ni.ox3 == 0) {
    //   sk = pmb->ckvs; ek = pmb->ckve;
    // } else if (nb.ni.ox3 > 0) {
    //   sk = pmb->ckge; ek = pmb->ckve-pmb->rcng-1;
    // } else {
    //   sk = pmb->ckvs+pmb->rcng+1; ek = pmb->ckvs+pmb->cng;
    // }
    if (DBGPR_BVALS_VC) {
      var.print_all();
      coarse_var.print_all();
    }

    idxLoadSameLevelRanges(nb.ni, si, ei, sj, ej, sk, ek, true);
    BufferUtility::PackData(coarse_var, buf, nl_, nu_,
                            si, ei, sj, ej, sk, ek, p);
  }

  return p;
}

//----------------------------------------------------------------------------------------
//! \fn int VertexCenteredBoundaryVariable::LoadBoundaryBufferToCoarser(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set vertex-centered boundary buffers for sending to a block on the coarser level

int VertexCenteredBoundaryVariable::LoadBoundaryBufferToCoarser(Real *buf,
                                                                const NeighborBlock& nb) {

  if (DBGPR_BVALS_VC) {
    coutYellow("VertexCenteredBoundaryVariable::LoadBoundaryBufferToCoarser\n");
    nb.print_all();
  }

  MeshBlock *pmb = pmy_block_;
  MeshRefinement *pmr = pmb->pmr;
  int si, sj, sk, ei, ej, ek;
  int p = 0;

  AthenaArray<Real> &var = *var_vc;
  AthenaArray<Real> &coarse_var = *coarse_buf;

  // vertices that are shared with adjacent MeshBlocks are to be copied to coarser level
  idxLoadToCoarserRanges(nb.ni, si, ei, sj, ej, sk, ek, false);
  pmr->RestrictVertexCenteredValues(var, coarse_var, nl_, nu_,
                                    si, ei, sj, ej, sk, ek);

  if (DBGPR_BVALS_VC)
    coutBoldRed("coarse_buf, buf");


  BufferUtility::PackData(coarse_var, buf, nl_, nu_, si, ei, sj, ej, sk, ek, p);

  if (pmy_mesh_->multilevel) {
    // double restrict required to populate coarse buffer of coarser level
    idxLoadToCoarserRanges(nb.ni, si, ei, sj, ej, sk, ek, true);
    pmr->RestrictTwiceToBufferVertexCenteredValues(var, buf, nl_, nu_,
                                                   si, ei, sj, ej, sk, ek, p);
  }

  return p;
}

//----------------------------------------------------------------------------------------
//! \fn int VertexCenteredBoundaryVariable::LoadBoundaryBufferToFiner(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set vertex-centered boundary buffers for sending to a block on the finer level

int VertexCenteredBoundaryVariable::LoadBoundaryBufferToFiner(Real *buf,
                                                              const NeighborBlock& nb) {
  if (DBGPR_BVALS_VC) {
    coutYellow("VertexCenteredBoundaryVariable::LoadBoundaryBufferToFiner\n");
    nb.print_all();
  }

  MeshBlock *pmb = pmy_block_;
  AthenaArray<Real> &var = *var_vc;
  int si, sj, sk, ei, ej, ek;
  int p = 0;

  if (DBGPR_BVALS_VC)
    coutBoldRed("var_vc, buf");

  idxLoadToFinerRanges(nb.ni, si, ei, sj, ej, sk, ek);
  BufferUtility::PackData(var, buf, nl_, nu_, si, ei, sj, ej, sk, ek, p);


  if (DBGPR_BVALS_VC) {
    coutBoldRed("MB::UWIL gid = ");
    printf("%d\n", pmb->gid);
  }

  return p;
}

//----------------------------------------------------------------------------------------
//! \fn inline void MeshBlock::SetIndexRangesSBSL(...)
//  \brief Set index ranges for a given dimension
inline void VertexCenteredBoundaryVariable::SetIndexRangesSBSL(
  int ox, int &ix_s, int &ix_e,
  int ix_vs, int ix_ve, int ix_ms, int ix_pe) {

  if (ox == 0) {
    ix_s = ix_vs;
    ix_e = ix_ve;
  } else if (ox > 0) {
    ix_s = ix_ve;
    ix_e = ix_pe;
  } else {
    ix_s = ix_ms;
    ix_e = ix_vs;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::SetBoundarySameLevel(Real *buf,
//                                                              const NeighborBlock& nb)
//  \brief Set vertex-centered boundary received from a block on the same level

void VertexCenteredBoundaryVariable::SetBoundarySameLevel(Real *buf,
                                                          const NeighborBlock& nb) {
  if (DBGPR_BVALS_VC) {
    coutYellow("VertexCenteredBoundaryVariable::SetBoundarySameLevel\n");
    nb.print_all();
  }

  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;
  AthenaArray<Real> &var = *var_vc;
  int p = 0;

  // BD: TODO implement
  ErrorIfPolarNotImplemented(nb);
  ErrorIfShearingBoxNotImplemented();

  // unpack all data additively
  // defer imposition (via suitable averaging) of consistency condition
  SetIndexRangesSBSL(nb.ni.ox1, si, ei,
                     pmb->ivs, pmb->ive, pmb->ims, pmb->ipe);
  SetIndexRangesSBSL(nb.ni.ox2, sj, ej,
                     pmb->jvs, pmb->jve, pmb->jms, pmb->jpe);
  SetIndexRangesSBSL(nb.ni.ox3, sk, ek,
                     pmb->kvs, pmb->kve, pmb->kms, pmb->kpe);


  // vertex consistency--------------------------------------------------------
  if (!node_mult_assembled) {
    int c_si, c_ei, c_sj, c_ej, c_sk, c_ek;

    SetIndexRangesSBSL(nb.ni.ox1, c_si, c_ei, c_ivs, c_ive, c_ims, c_ipe);
    SetIndexRangesSBSL(nb.ni.ox2, c_sj, c_ej, c_jvs, c_jve, c_jms, c_jpe);
    SetIndexRangesSBSL(nb.ni.ox3, c_sk, c_ek, c_kvs, c_kve, c_kms, c_kpe);

    for (int k=c_sk; k<=c_ek; ++k)
      for (int j=c_sj; j<=c_ej; ++j)
        for (int i=c_si; i<=c_ei; ++i)
          node_mult(0, k, j, i) += 1;
  }
  //---------------------------------------------------------------------------

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

    if (DBGPR_BVALS_VC)
      var.print_all();

  } else {
    BufferUtility::UnpackDataAdd(buf, var, nl_, nu_,
                                 si, ei, sj, ej, sk, ek, p);
  }
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  if (pmy_mesh_->multilevel) {
    // note: unpacked shared nodes additively unpacked-
    // consistency conditions will need to be applied to the coarse variable

    MeshRefinement *pmr = pmb->pmr;
    AthenaArray<Real> &coarse_var = *coarse_buf;

    SetIndexRangesSBSL(nb.ni.ox1, si, ei,
                       pmb->civs, pmb->cive, pmb->cims, pmb->cipe);
    SetIndexRangesSBSL(nb.ni.ox2, sj, ej,
                       pmb->cjvs, pmb->cjve, pmb->cjms, pmb->cjpe);
    SetIndexRangesSBSL(nb.ni.ox3, sk, ek,
                       pmb->ckvs, pmb->ckve, pmb->ckms, pmb->ckpe);

    BufferUtility::UnpackDataAdd(buf, coarse_var, nl_, nu_,
                                 si, ei, sj, ej, sk, ek, p);


    /*
    // immediately restrict newly populated ghost-zones to coarse buffer
    // first restrict data communicated on fundamental representation
    //
    // note: unpacked shared nodes are overwritten here.
    // consistency conditions will need to be applied to the coarse variable

    coutBoldRed("Initial restriction\n");
    MeshRefinement *pmr = pmb->pmr;
    AthenaArray<Real> &coarse_var = *coarse_buf;

    if (nb.ni.ox1 == 0) {
      si = pmb->civs; ei = pmb->cive;
    } else if (nb.ni.ox1 > 0) {
      si = pmb->cive; ei = pmb->cive+pmb->rcng;
    } else {
      si = pmb->civs-pmb->rcng; ei = pmb->civs;
    }

    if (nb.ni.ox2 == 0) {
      sj = pmb->cjvs; ej = pmb->cjve;
    } else if (nb.ni.ox2 > 0) {
      sj = pmb->cjve; ej = pmb->cjve+pmb->rcng;
    } else {
      sj = pmb->cjvs-pmb->rcng; ej = pmb->cjvs;
    }

    if (nb.ni.ox3 == 0) {
      sk = pmb->ckvs; ek = pmb->ckve;
    } else if (nb.ni.ox3 > 0) {
      sk = pmb->ckve; ek = pmb->ckve+pmb->rcng;
    } else {
      sk = pmb->ckvs-pmb->rcng; ek = pmb->ckvs;
    }

    var.print_all();
    coarse_var.print_all();

    pmr->RestrictVertexCenteredValues(var, coarse_var, nl_, nu_,
                                      si, ei, sj, ej, sk, ek);
    coarse_var.print_all();

    coutBoldRed("Unpack remaining coarse\n");
    if (nb.ni.ox1 == 0) {
      si = pmb->civs; ei = pmb->cive;
    } else if (nb.ni.ox1 > 0) {
      si = pmb->cive+pmb->rcng+1; ei = pmb->cipe;
    } else {
      si = pmb->cims; ei = pmb->civs-pmb->rcng-1;
    }

    if (nb.ni.ox2 == 0) {
      sj = pmb->cjvs; ej = pmb->cjve;
    } else if (nb.ni.ox2 > 0) {
      sj = pmb->cjve+pmb->rcng+1; ej = pmb->cjpe;
    } else {
      sj = pmb->cjms; ej = pmb->cjvs-pmb->rcng-1;
    }

    if (nb.ni.ox3 == 0) {
      sk = pmb->ckvs; ek = pmb->ckve;
    } else if (nb.ni.ox3 > 0) {
      sk = pmb->ckve+pmb->rcng+1; ek = pmb->ckpe;
    } else {
      sk = pmb->ckms; ek = pmb->ckvs-pmb->rcng-1;
    }

    BufferUtility::UnpackData(buf, coarse_var, nl_, nu_,
                              si, ei, sj, ej, sk, ek, p);
    */
  }
  //////////////////////////////////////////////////////////////////////////////



  return;
}

//----------------------------------------------------------------------------------------
//! \fn inline void MeshBlock::SetIndexRangesSBFC(...)
//  \brief Set index ranges for a given dimension
inline void VertexCenteredBoundaryVariable::SetIndexRangesSBFC(
  int ox, int &ix_s, int &ix_e,
  int ix_cvs, int ix_cve, int ix_cms, int ix_cme, int ix_cps, int ix_cpe,
  bool level_flag) {

  if (ox == 0) {
    ix_s = ix_cvs; ix_e = ix_cve;
    if (level_flag) {
      ix_s = ix_cvs; ix_e = ix_cpe;
    } else {
      ix_s = ix_cms; ix_e = ix_cve;
    }
  } else if (ox > 0)  {
    ix_s = ix_cps; ix_e = ix_cpe;
  } else {
    ix_s = ix_cms; ix_e = ix_cme;
  }


  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::SetBoundaryFromCoarser(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set vertex-centered prolongation buffer received from a block on a coarser level

void VertexCenteredBoundaryVariable::SetBoundaryFromCoarser(Real *buf,
                                                            const NeighborBlock& nb) {

  // populating from a coarser level; do not touch shared vertices
  if (DBGPR_BVALS_VC) {
    coutYellow("VertexCenteredBoundaryVariable::SetBoundaryFromCoarser\n");
    nb.print_all();
  }

  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;

  AthenaArray<Real> &coarse_var = *coarse_buf;


  int p = 0;
  // BD: TODO implement
  ErrorIfPolarNotImplemented(nb);

  // // [1d] modifies si, ei
  // if (nb.ni.ox1 == 0) {
  //   si = pmb->civs; ei = pmb->cive;
  //   if ((pmb->loc.lx1 & 1LL) == 0LL) {
  //     si = pmb->civs; ei = pmb->cipe;
  //   } else {
  //     si = pmb->cims; ei = pmb->cive;
  //   }
  // } else if (nb.ni.ox1 > 0)  {
  //   si = pmb->cips; ei = pmb->cipe;
  // } else {
  //   si = pmb->cims; ei = pmb->cime;
  // }


  // [1d] modifies si, ei
  SetIndexRangesSBFC(nb.ni.ox1, si, ei,
                     pmb->civs, pmb->cive, pmb->cims, pmb->cime,
                     pmb->cips, pmb->cipe, (pmb->loc.lx1 & 1LL) == 0LL);

  // [2d] modifies sj, ej
  SetIndexRangesSBFC(nb.ni.ox2, sj, ej,
                     pmb->cjvs, pmb->cjve, pmb->cjms, pmb->cjme,
                     pmb->cjps, pmb->cjpe, (pmb->loc.lx2 & 1LL) == 0LL);

  // [3d] modifies sk, ek
  SetIndexRangesSBFC(nb.ni.ox3, sk, ek,
                     pmb->ckvs, pmb->ckve, pmb->ckms, pmb->ckme,
                     pmb->ckps, pmb->ckpe, (pmb->loc.lx3 & 1LL) == 0LL);


  //------

  // vertex consistency--------------------------------------------------------
  if (!node_mult_assembled) {
    int c_si, c_ei, c_sj, c_ej, c_sk, c_ek;

    int const c_ime = c_ims;
    int const c_ips = c_ipe;

    int const c_jme = c_jms;
    int const c_jps = c_jpe;

    int const c_kme = c_kms;
    int const c_kps = c_kpe;

    // [1d] modifies si, ei
    SetIndexRangesSBFC(nb.ni.ox1, c_si, c_ei,
                      c_ivs, c_ive, c_ims, c_ime,
                      c_ips, c_ipe, (pmb->loc.lx1 & 1LL) == 0LL);

    // [2d] modifies sj, ej
    SetIndexRangesSBFC(nb.ni.ox2, c_sj, c_ej,
                      c_jvs, c_jve, c_jms, c_jme,
                      c_jps, c_jpe, (pmb->loc.lx2 & 1LL) == 0LL);

    // [3d] modifies sk, ek
    SetIndexRangesSBFC(nb.ni.ox3, c_sk, c_ek,
                      c_kvs, c_kve, c_kms, c_kme,
                      c_kps, c_kpe, (pmb->loc.lx3 & 1LL) == 0LL);

    for (int k=c_sk; k<=c_ek; ++k)
      for (int j=c_sj; j<=c_ej; ++j)
        for (int i=c_si; i<=c_ei; ++i)
          node_mult(0, k, j, i) += 1;
  }
  //---------------------------------------------------------------------------


  // write like this to regroup idx logic for debug
  if (!FILL_WAVE_BND_FRC) {
    if (DBGPR_BVALS_VC)
      coarse_var.print_all();

    BufferUtility::UnpackData(buf, coarse_var, nl_, nu_,
                              si, ei, sj, ej, sk, ek, p);

    if (DBGPR_BVALS_VC) {
      coutBoldRed("ng, cng, pmb->cnghost=");
      printf("%d,%d,%d\n", pmb->ng, pmb->cng, pmb->cnghost);
    }

    return;
  }

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
    if (DBGPR_BVALS_VC)
      var.print_all();

    pmb->DebugWaveMeshBlock(var, si, ei, sj, ej, sk, ek, false);

    if (DBGPR_BVALS_VC)
      var.print_all();
    // pmb->DebugWaveMeshBlock(var,
    //                         pmb->ims, pmb->ipe,
    //                         pmb->jms, pmb->jpe,
    //                         pmb->kms, pmb->kpe,
    //                         false);

    if (flag)
      Q();
  }
  //////////////////////////////////////////////////////////////////////////////

  if (DBGPR_BVALS_VC) {
    coutBoldRed("MB::UWIL gid = ");
    printf("%d\n", pmb->gid);
  }

  // if (nb.ni.ox1 < 0)
  //   Q();

  // printf("x1f: ");
  // pmb->pcoord->x1f.print_data("%1.2f");

  // if (nb.ni.ox1 == -1)
  //   Q();
  // if ((nb.ni.ox1 == 1) and (nb.ni.ox2 == -1))
  //   Q();


  return;
}

//----------------------------------------------------------------------------------------
//! \fn inline void MeshBlock::SetIndexRangesSBFF(...)
//  \brief Set index ranges for a given dimension
inline void VertexCenteredBoundaryVariable::SetIndexRangesSBFF(
  int ox, int &ix_s, int &ix_e, int ix_vs, int ix_ve, int ix_ms, int ix_pe,
  int fi1, int fi2, int axis_half_size, bool size_flag, bool offset_flag) {

  if (ox == 0) {
    ix_s = ix_vs;
    ix_e = ix_ve;

    if (size_flag)
      if (offset_flag) {
        if (fi1 == 1) {
          ix_s += axis_half_size;
        } else {
          ix_e -= axis_half_size;
        }
      } else {
        if (fi2 == 1) {
          ix_s += axis_half_size;
        } else {
          ix_e -= axis_half_size;
        }
      }
  } else if (ox > 0) {
    ix_s = ix_ve;
    ix_e = ix_pe;
  } else {
    ix_s = ix_ms;
    ix_e = ix_vs;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::SetBoundaryFromFiner(Real *buf,
//                                                              const NeighborBlock& nb)
//  \brief Set vertex-centered boundary received from a block on a finer level
void VertexCenteredBoundaryVariable::SetBoundaryFromFiner(Real *buf,
                                                          const NeighborBlock& nb) {
  // populating from finer level; shared vertices are corrected
  if (DBGPR_BVALS_VC) {
    coutYellow("VertexCenteredBoundaryVariable::SetBoundaryFromFiner\n");
    nb.print_all();
  }

  MeshBlock *pmb = pmy_block_;
  AthenaArray<Real> &var = *var_vc;
  // receive already restricted data
  int si, sj, sk, ei, ej, ek;
  int p = 0;

  // BD: TODO implement
  ErrorIfPolarNotImplemented(nb);

  // refactor (see below)
  /*
  // [1d] modifies si, ei
  if (nb.ni.ox1 == 0) {
    si = pmb->ivs;
    ei = pmb->ive;
    if (nb.ni.fi1 == 1) {
      // (ox1=0, ~(ox2=0 /\ ox3=0), fi1=1, fi2)
      //
      // (si, ei) = (is + block_size.nx1 / 2, ie)
      si += pmb->block_size.nx1 / 2;
    } else {
      // (ox1=0, ~(ox2=0 /\ ox3=0), fi1=0, fi2)
      //
      // (si, ei) = (is, ie - block_size.nx1 / 2)
      ei -= pmb->block_size.nx1 / 2;
    }
  } else if (nb.ni.ox1 > 0) {
    // (ox1>0, ox2, ox3, fi1, fi2)
    //
    // (si, ei) = (ie + 1, ie + NGHOST)
    si = pmb->ive;
    ei = pmb->ipe;
  } else {
    // (ox1<0, ox2, ox3, fi1, fi2)
    //
    // (si, ei) = (is - NGHOST, is - 1)
    si = pmb->ims;
    ei = pmb->ivs;
  }

  // [2d] modifies sj, ej
  if (nb.ni.ox2 == 0) {
    sj = pmb->jvs;
    ej = pmb->jve;

    if (pmb->block_size.nx2 > 1)
      if (nb.ni.ox1 != 0) {
        if (nb.ni.fi1 == 1) {
          // (ox1!=0, ox2=0, ox3, fi1=1, fi2)
          //
          // (sj, ej) = (js + block_size.nx2 / 2, je)
          sj += pmb->block_size.nx2 / 2;
        } else {
          // (ox1!=0, ox2=0, ox3, fi1=0, fi2)
          //
          // (sj, ej) = (js, je - block_size.nx2 / 2)
          ej -= pmb->block_size.nx2 / 2;
        }
      } else {
        if (nb.ni.fi2 == 1) {
          // (ox1=0, ox2=0, ox3!=0, fi1, fi2=1)
          //
          // (sj, ej) = (js + block_size.nx2 / 2, je)
          sj += pmb->block_size.nx2 / 2;
        } else {
          // (ox1=0, ox2=0, ox3!=0, fi1, fi2=0)
          //
          // (sj, ej) = (js, je - block_size.nx2 / 2)
          ej -= pmb->block_size.nx2 / 2;
        }
      }
  } else if (nb.ni.ox2 > 0) {
    // (ox1, ox2>0, ox3, fi1, fi2)
    //
    // (sj, ej) = (je + 1, je + NGHOST)
    sj = pmb->jve;
    ej = pmb->jpe;
  } else {
    // (ox1, ox2<0, ox3, fi1, fi2)
    //
    // (sj, ej) = (js - NGHOST, js - 1)
    sj = pmb->jms;
    ej = pmb->jvs;
  }

  // [3d] modifies sj, ek
  if (nb.ni.ox3 == 0) {
    sk = pmb->kvs;
    ek = pmb->kve;

    if (pmb->block_size.nx3 > 1)
      if (nb.ni.ox1 != 0 && nb.ni.ox2 != 0) {
        if (nb.ni.fi1 == 1) {
          // (ox1!=0, ox2!=0, ox3=0, fi1=1, fi2)
          //
          // (sk, ek) = (ks + block_size.nx3 / 2, ke)
          sk += pmb->block_size.nx3 / 2;
        } else {
          // (ox1!=0, ox2!=0, ox3=0, fi1=0, fi2)
          //
          // (sk, ek) = (ks, ke - block_size.nx3 / 2)
          ek -= pmb->block_size.nx3 / 2;
        }
      } else {
        if (nb.ni.fi2 == 1) {
          // (~(ox1!=0 /\ ox2!=0), ox3=0, fi1, fi2=1)
          //
          // (sk, ek) = (ks + block_size.nx3 / 2, ke)
          sk += pmb->block_size.nx3 / 2;
        } else {
          // (~(ox1!=0 /\ ox2!=0), ox3=0, fi1, fi2=0)
          //
          // (sk, ek) = (ks, ke - block_size.nx3 / 2)
          ek -= pmb->block_size.nx3 / 2;
        }
      }
  } else if (nb.ni.ox3 > 0) {
    // (ox1, ox2, ox3>0, fi1, fi2)
    //
    // (sk, ek) = (ke + 1, ke + NGHOST)
    sk = pmb->kve;
    ek = pmb->kpe;
  } else {
    // (ox1, ox2, ox3<0, fi1, fi2)
    //
    // (sk, ek) = (ks - NGHOST, ks - 1)
    sk = pmb->kms;
    ek = pmb->kvs;
  }
  */

  SetIndexRangesSBFF(nb.ni.ox1, si, ei,
                     pmb->ivs, pmb->ive, pmb->ims, pmb->ipe,
                     nb.ni.fi1, nb.ni.fi2,
                     pmb->block_size.nx1 / 2,
                     true,
                     true);

  SetIndexRangesSBFF(nb.ni.ox2, sj, ej,
                     pmb->jvs, pmb->jve, pmb->jms, pmb->jpe,
                     nb.ni.fi1, nb.ni.fi2,
                     pmb->block_size.nx2 / 2,
                     (pmb->block_size.nx2 > 1),
                     (nb.ni.ox1 != 0));

  SetIndexRangesSBFF(nb.ni.ox3, sk, ek,
                     pmb->kvs, pmb->kve, pmb->kms, pmb->kpe,
                     nb.ni.fi1, nb.ni.fi2,
                     pmb->block_size.nx3 / 2,
                     (pmb->block_size.nx3 > 1),
                     (nb.ni.ox1 != 0 && nb.ni.ox2 != 0));

  //////////////////////////////////////////////////////////////////////////////

  if (FILL_WAVE_BND_FRF) {
    if (DBGPR_BVALS_VC)
      var.print_all();

    pmb->DebugWaveMeshBlock(var, si, ei, sj, ej, sk, ek, true);

    if (DBGPR_BVALS_VC)
      var.print_all();

  } else {

    if (DBGPR_BVALS_VC)
      coutBoldRed("buf, var_vc");

    BufferUtility::UnpackDataAdd(buf, var, nl_, nu_,
                                 si, ei, sj, ej, sk, ek, p);
  }

  // vertex consistency--------------------------------------------------------
  if (!node_mult_assembled) {
    int c_si, c_ei, c_sj, c_ej, c_sk, c_ek;

    SetIndexRangesSBFF(nb.ni.ox1, c_si, c_ei,
                      c_ivs, c_ive, c_ims, c_ipe,
                      nb.ni.fi1, nb.ni.fi2,
                      2,
                      true,
                      true);

    SetIndexRangesSBFF(nb.ni.ox2, c_sj, c_ej,
                      c_jvs, c_jve, c_jms, c_jpe,
                      nb.ni.fi1, nb.ni.fi2,
                      2,
                      (pmb->block_size.nx2 > 1),
                      (nb.ni.ox1 != 0));

    SetIndexRangesSBFF(nb.ni.ox3, c_sk, c_ek,
                      c_kvs, c_kve, c_kms, c_kpe,
                      nb.ni.fi1, nb.ni.fi2,
                      2,
                      (pmb->block_size.nx3 > 1),
                      (nb.ni.ox1 != 0 && nb.ni.ox2 != 0));

    for (int k=c_sk; k<=c_ek; ++k)
      for (int j=c_sj; j<=c_ej; ++j)
        for (int i=c_si; i<=c_ei; ++i)
          node_mult(0, k, j, i) += 1;
  }
  //---------------------------------------------------------------------------

  if (pmy_mesh_->multilevel) {
    AthenaArray<Real> &coarse_var = *coarse_buf;

    if (DBGPR_BVALS_VC) {
      coutBoldBlue("coarse_var\n");
      coarse_var.print_all();
    }

    SetIndexRangesSBFF(nb.ni.ox1, si, ei,
                       pmb->civs, pmb->cive, pmb->cims, pmb->cipe,
                       nb.ni.fi1, nb.ni.fi2,
                       pmb->block_size.nx1 / 4,
                       true,
                       true);

    SetIndexRangesSBFF(nb.ni.ox2, sj, ej,
                       pmb->cjvs, pmb->cjve, pmb->cjms, pmb->cjpe,
                       nb.ni.fi1, nb.ni.fi2,
                       pmb->block_size.nx2 / 4,
                       (pmb->block_size.nx2 > 1),
                       (nb.ni.ox1 != 0));

    SetIndexRangesSBFF(nb.ni.ox3, sk, ek,
                       pmb->ckvs, pmb->ckve, pmb->ckms, pmb->ckpe,
                       nb.ni.fi1, nb.ni.fi2,
                       pmb->block_size.nx3 / 4,
                       (pmb->block_size.nx3 > 1),
                       (nb.ni.ox1 != 0 && nb.ni.ox2 != 0));

    if (DBGPR_BVALS_VC)
      coutBoldRed("buf, coarse_var");

    BufferUtility::UnpackDataAdd(buf, coarse_var, nl_, nu_,
                                 si, ei, sj, ej, sk, ek, p);

  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::RestrictNonGhost()
//  \brief populate coarser buffer with restricted data

void VertexCenteredBoundaryVariable::RestrictNonGhost() {
  if (DBGPR_BVALS_VC)
    coutYellow("VertexCenteredBoundaryVariable::RestrictNonGhost\n");

  MeshBlock *pmb = pmy_block_;
  MeshRefinement *pmr = pmb->pmr;
  int si, sj, sk, ei, ej, ek;

  AthenaArray<Real> &var = *var_vc;
  AthenaArray<Real> &coarse_var = *coarse_buf;

  si = pmb->civs; ei = pmb->cive;
  sj = pmb->cjvs; ej = pmb->cjve;
  sk = pmb->ckvs; ek = pmb->ckve;

  pmr->RestrictVertexCenteredValues(var, coarse_var, nl_, nu_,
                                    si, ei, sj, ej, sk, ek);

  return;
}


void VertexCenteredBoundaryVariable::SendBoundaryBuffers() {
  if (!node_mult_assembled)
    PrepareNodeMult();

  MeshBlock *pmb = pmy_block_;

  // restrict all data (except ghosts) to coarse buffer
  if (pmy_mesh_->multilevel) {
    AthenaArray<Real> &coarse_var = *coarse_buf;
    coarse_var.ZeroClear();

    RestrictNonGhost();
  }
  BoundaryVariable::SendBoundaryBuffers();
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::SetBoundaries()
//  \brief set the vertex-centered boundary data

void VertexCenteredBoundaryVariable::SetBoundaries() {
  if (DBGPR_BVALS_VC)
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
  if (DBGPR_BVALS_VC)
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
  if (DBGPR_BVALS_VC)
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
      // if (pmy_mesh_->multilevel && nb.ni.type == NeighborConnect::face) {
      //   int size;
      //   if (nb.fid == 0 || nb.fid == 1)
      //     size = ((pmb->block_size.nx2 + 1)/2)*((pmb->block_size.nx3 + 1)/2);
      //   else if (nb.fid == 2 || nb.fid == 3)
      //     size = ((pmb->block_size.nx1 + 1)/2)*((pmb->block_size.nx3 + 1)/2);
      //   else // (nb.fid == 4 || nb.fid == 5)
      //     size = ((pmb->block_size.nx1 + 1)/2)*((pmb->block_size.nx2 + 1)/2);
      //   size *= (nu_ + 1);
      //   if (nb.snb.level < mylevel) { // send to coarser
      //     tag = pbval_->CreateBvalsMPITag(nb.snb.lid, nb.targetid, cc_flx_phys_id_);
      //     // if (bd_var_flcor_.req_send[nb.bufid] != MPI_REQUEST_NULL)
      //     //   MPI_Request_free(&bd_var_flcor_.req_send[nb.bufid]);
      //     // MPI_Send_init(bd_var_flcor_.send[nb.bufid], size, MPI_ATHENA_REAL,
      //     //               nb.snb.rank, tag, MPI_COMM_WORLD,
      //     //               &(bd_var_flcor_.req_send[nb.bufid]));
      //   } else if (nb.snb.level > mylevel) { // receive from finer
      //     tag = pbval_->CreateBvalsMPITag(pmb->lid, nb.bufid, cc_flx_phys_id_);
      //     if (bd_var_flcor_.req_recv[nb.bufid] != MPI_REQUEST_NULL)
      //       MPI_Request_free(&bd_var_flcor_.req_recv[nb.bufid]);
      //     MPI_Recv_init(bd_var_flcor_.recv[nb.bufid], size, MPI_ATHENA_REAL,
      //                   nb.snb.rank, tag, MPI_COMM_WORLD,
      //                   &(bd_var_flcor_.req_recv[nb.bufid]));
      //   }
      // }


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
  if (DBGPR_BVALS_VC)
    coutYellow("VertexCenteredBoundaryVariable::ClearBoundary\n");

  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    bd_var_.flag[nb.bufid] = BoundaryStatus::waiting;
    bd_var_.sflag[nb.bufid] = BoundaryStatus::waiting;

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

  return;
}
