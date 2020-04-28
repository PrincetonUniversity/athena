//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file idx_utilities_vc.cpp
//  \brief Various utitilies for indicies and buffers

// C headers

// C++ headers
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "bvals_vc.hpp"

//----------------------------------------------------------------------------------------
//! \fn inline void VertexCenteredBoundaryVariable::AccumulateBufferSize(...)
//  \brief Buffer size accumulator
inline void VertexCenteredBoundaryVariable::AccumulateBufferSize(
  int sn, int en,
  int si, int ei, int sj, int ej, int sk, int ek,
  int &offset, int ijk_step=1) {

  for (int n=sn; n<=en; ++n) {
    for (int k=sk; k<=ek; k+=ijk_step) {
      for (int j=sj; j<=ej; j+=ijk_step) {
//#pragma omp parallel for shared(offset, si, ei)
        for (int i=si; i<=ei; i+=ijk_step)
          offset++;
      }
    }
  }

}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::idxLoadSameLevelRanges(...)
//  \brief Compute indicial ranges for LoadBoundaryBufferSameLevel
void VertexCenteredBoundaryVariable::idxLoadSameLevelRanges(
  const NeighborIndexes& ni,
  int &si, int &ei, int &sj, int &ej, int &sk, int &ek,
  bool is_coarse=false) {

  MeshBlock *pmb = pmy_block_;

  if (!is_coarse) {
    // shared vertex is packed
    si = (ni.ox1 > 0) ? pmb->ige : pmb->ivs;
    ei = (ni.ox1 < 0) ? pmb->igs : pmb->ive;

    sj = (ni.ox2 > 0) ? pmb->jge : pmb->jvs;
    ej = (ni.ox2 < 0) ? pmb->jgs : pmb->jve;

    sk = (ni.ox3 > 0) ? pmb->kge : pmb->kvs;
    ek = (ni.ox3 < 0) ? pmb->kgs : pmb->kve;

  } else {
    // [coarse] shared vertex is packed
    si = (ni.ox1 > 0) ? pmb->cige : pmb->civs;
    ei = (ni.ox1 < 0) ? pmb->cigs : pmb->cive;

    sj = (ni.ox2 > 0) ? pmb->cjge : pmb->cjvs;
    ej = (ni.ox2 < 0) ? pmb->cjgs : pmb->cjve;

    sk = (ni.ox3 > 0) ? pmb->ckge : pmb->ckvs;
    ek = (ni.ox3 < 0) ? pmb->ckgs : pmb->ckve;
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::idxLoadToCoarserRanges(...)
//  \brief Compute indicial ranges for LoadBoundaryBufferToCoarser
void VertexCenteredBoundaryVariable::idxLoadToCoarserRanges(
  const NeighborIndexes& ni,
  int &si, int &ei, int &sj, int &ej, int &sk, int &ek,
  bool is_coarse=false) {

  MeshBlock *pmb = pmy_block_;

  if (!is_coarse) {
    // vertices that are shared with adjacent MeshBlocks are to be copied to
    // coarser level
    int ng = pmb->ng;

    si = (ni.ox1 > 0) ? (pmb->cive - ng) : pmb->civs;
    ei = (ni.ox1 < 0) ? (pmb->civs + ng) : pmb->cive;

    sj = (ni.ox2 > 0) ? (pmb->cjve - ng) : pmb->cjvs;
    ej = (ni.ox2 < 0) ? (pmb->cjvs + ng) : pmb->cjve;

    sk = (ni.ox3 > 0) ? (pmb->ckve - ng) : pmb->ckvs;
    ek = (ni.ox3 < 0) ? (pmb->ckvs + ng) : pmb->ckve;

  } else {
    int cng = 2 * pmb->cng;  // "2 for coarse-coarse"
    si = (ni.ox1 > 0) ? (pmb->cive - cng) : pmb->civs;
    ei = (ni.ox1 < 0) ? (pmb->civs + cng) : pmb->cive;

    sj = (ni.ox2 > 0) ? (pmb->cjve - cng) : pmb->cjvs;
    ej = (ni.ox2 < 0) ? (pmb->cjvs + cng) : pmb->cjve;

    sk = (ni.ox3 > 0) ? (pmb->ckve - cng) : pmb->ckvs;
    ek = (ni.ox3 < 0) ? (pmb->ckvs + cng) : pmb->ckve;
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::idxLoadToFinerRanges(...)
//  \brief Compute indicial ranges for LoadBoundaryBufferToFiner
void VertexCenteredBoundaryVariable::idxLoadToFinerRanges(
  const NeighborIndexes& ni,
  int &si, int &ei, int &sj, int &ej, int &sk, int &ek) {

  MeshBlock *pmb = pmy_block_;

  int ng = pmb->ng;
  int cng = pmb->cng;

  int tmp = cng;

  if (ni.ox1 > 0) {
    si = pmb->ive-cng, ei = pmb->ive-1;
  } else if (ni.ox1 < 0) {
    si = pmb->ivs+1, ei = pmb->ivs+cng;
  } else {
    si = pmb->ivs, ei = pmb->ive;
    if (ni.fi1 == 1)
      si += pmb->block_size.nx1/2 - tmp;
    else
      ei -= pmb->block_size.nx1/2 - tmp;

  }

  if (ni.ox2 > 0) {
    sj = pmb->jve-cng, ej = pmb->jve-1;
  } else if (ni.ox2 < 0) {
    sj = pmb->jvs+1, ej = pmb->jvs+cng;
  } else {
    sj = pmb->jvs, ej = pmb->jve;

    if (pmb->block_size.nx2 > 1)
      if (ni.ox1 != 0) {
        if (ni.fi1 == 1)
          sj += pmb->block_size.nx2/2 - tmp;
        else
          ej -= pmb->block_size.nx2/2 - tmp;
      } else {
        if (ni.fi2 == 1)
          sj += pmb->block_size.nx2/2 - tmp;
        else
          ej -= pmb->block_size.nx2/2 - tmp;

      }
  }

  if (ni.ox3 > 0) {
    sk = pmb->kve-cng, ek = pmb->kve-1;
  } else if (ni.ox3 < 0) {
    sk = pmb->kvs+1, ek = pmb->kvs+cng;
  } else {
    sk = pmb->kvs, ek = pmb->kve;

    if (pmb->block_size.nx3 > 1)
      if ((ni.ox1 != 0) && (ni.ox2 != 0)) {
        if (ni.fi1 == 1)
          sk += pmb->block_size.nx3/2 - tmp;
        else
          ek -= pmb->block_size.nx3/2 - tmp;
      } else {
        if (ni.fi2 == 1)
          sk += pmb->block_size.nx3/2 - tmp;
        else
          ek -= pmb->block_size.nx3/2 - tmp;

      }
  }

  return;
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
//! \fn void VertexCenteredBoundaryVariable::idxSetSameLevelRanges(...)
//  \brief Compute indicial ranges for SetBoundarySameLevel
void VertexCenteredBoundaryVariable::idxSetSameLevelRanges(
  const NeighborIndexes& ni,
  int &si, int &ei, int &sj, int &ej, int &sk, int &ek,
  int type) {
  // type = 1 for fundamental, 2 for coarse, 3 for node_mult

  MeshBlock *pmb = pmy_block_;

  if (type == 1) {
    SetIndexRangesSBSL(ni.ox1, si, ei,
                       pmb->ivs, pmb->ive, pmb->ims, pmb->ipe);
    SetIndexRangesSBSL(ni.ox2, sj, ej,
                       pmb->jvs, pmb->jve, pmb->jms, pmb->jpe);
    SetIndexRangesSBSL(ni.ox3, sk, ek,
                       pmb->kvs, pmb->kve, pmb->kms, pmb->kpe);
  } else if (type == 2) {
    SetIndexRangesSBSL(ni.ox1, si, ei,
                       pmb->civs, pmb->cive, pmb->cims, pmb->cipe);
    SetIndexRangesSBSL(ni.ox2, sj, ej,
                       pmb->cjvs, pmb->cjve, pmb->cjms, pmb->cjpe);
    SetIndexRangesSBSL(ni.ox3, sk, ek,
                       pmb->ckvs, pmb->ckve, pmb->ckms, pmb->ckpe);
  } else if (type == 3) {
    SetIndexRangesSBSL(ni.ox1, si, ei, c_ivs, c_ive, c_ims, c_ipe);
    SetIndexRangesSBSL(ni.ox2, sj, ej, c_jvs, c_jve, c_jms, c_jpe);
    SetIndexRangesSBSL(ni.ox3, sk, ek, c_kvs, c_kve, c_kms, c_kpe);
  }

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
//! \fn void VertexCenteredBoundaryVariable::idxSetFromCoarserRanges(...)
//  \brief Compute indicial ranges for SetBoundaryFromCoarser
void VertexCenteredBoundaryVariable::idxSetFromCoarserRanges(
  const NeighborIndexes& ni,
  int &si, int &ei, int &sj, int &ej, int &sk, int &ek,
  bool is_node_mult=false) {

  MeshBlock *pmb = pmy_block_;

  if (!is_node_mult) {
    // [1d] modifies si, ei
    SetIndexRangesSBFC(ni.ox1, si, ei,
                       pmb->civs, pmb->cive, pmb->cims, pmb->cime,
                       pmb->cips, pmb->cipe, (pmb->loc.lx1 & 1LL) == 0LL);

    // [2d] modifies sj, ej
    SetIndexRangesSBFC(ni.ox2, sj, ej,
                       pmb->cjvs, pmb->cjve, pmb->cjms, pmb->cjme,
                       pmb->cjps, pmb->cjpe, (pmb->loc.lx2 & 1LL) == 0LL);

    // [3d] modifies sk, ek
    SetIndexRangesSBFC(ni.ox3, sk, ek,
                       pmb->ckvs, pmb->ckve, pmb->ckms, pmb->ckme,
                       pmb->ckps, pmb->ckpe, (pmb->loc.lx3 & 1LL) == 0LL);
  } else {
    int const c_ime = c_ims;
    int const c_ips = c_ipe;

    int const c_jme = c_jms;
    int const c_jps = c_jpe;

    int const c_kme = c_kms;
    int const c_kps = c_kpe;

    // [1d] modifies si, ei
    SetIndexRangesSBFC(ni.ox1, si, ei,
                       c_ivs, c_ive, c_ims, c_ime,
                       c_ips, c_ipe, (pmb->loc.lx1 & 1LL) == 0LL);

    // [2d] modifies sj, ej
    SetIndexRangesSBFC(ni.ox2, sj, ej,
                       c_jvs, c_jve, c_jms, c_jme,
                       c_jps, c_jpe, (pmb->loc.lx2 & 1LL) == 0LL);

    // [3d] modifies sk, ek
    SetIndexRangesSBFC(ni.ox3, sk, ek,
                       c_kvs, c_kve, c_kms, c_kme,
                       c_kps, c_kpe, (pmb->loc.lx3 & 1LL) == 0LL);
  }

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
//! \fn void VertexCenteredBoundaryVariable::idxSetFromFinerRanges(...)
//  \brief Compute indicial ranges for SetBoundaryFromFiner
void VertexCenteredBoundaryVariable::idxSetFromFinerRanges(
  const NeighborIndexes& ni,
  int &si, int &ei, int &sj, int &ej, int &sk, int &ek,
  int type) {
  // type = 1 for fundamental, 2 for coarse, 3 for node_mult

  MeshBlock *pmb = pmy_block_;

  if (type == 1) {
    SetIndexRangesSBFF(ni.ox1, si, ei,
                       pmb->ivs, pmb->ive, pmb->ims, pmb->ipe,
                       ni.fi1, ni.fi2,
                       pmb->block_size.nx1 / 2,
                       true,
                       true);

    SetIndexRangesSBFF(ni.ox2, sj, ej,
                       pmb->jvs, pmb->jve, pmb->jms, pmb->jpe,
                       ni.fi1, ni.fi2,
                       pmb->block_size.nx2 / 2,
                       (pmb->block_size.nx2 > 1),
                       (ni.ox1 != 0));

    SetIndexRangesSBFF(ni.ox3, sk, ek,
                       pmb->kvs, pmb->kve, pmb->kms, pmb->kpe,
                       ni.fi1, ni.fi2,
                       pmb->block_size.nx3 / 2,
                       (pmb->block_size.nx3 > 1),
                       (ni.ox1 != 0 && ni.ox2 != 0));

  } else if (type == 2) {
    SetIndexRangesSBFF(ni.ox1, si, ei,
                       pmb->civs, pmb->cive, pmb->cims, pmb->cipe,
                       ni.fi1, ni.fi2,
                       pmb->block_size.nx1 / 4,
                       true,
                       true);

    SetIndexRangesSBFF(ni.ox2, sj, ej,
                       pmb->cjvs, pmb->cjve, pmb->cjms, pmb->cjpe,
                       ni.fi1, ni.fi2,
                       pmb->block_size.nx2 / 4,
                       (pmb->block_size.nx2 > 1),
                       (ni.ox1 != 0));

    SetIndexRangesSBFF(ni.ox3, sk, ek,
                       pmb->ckvs, pmb->ckve, pmb->ckms, pmb->ckpe,
                       ni.fi1, ni.fi2,
                       pmb->block_size.nx3 / 4,
                       (pmb->block_size.nx3 > 1),
                       (ni.ox1 != 0 && ni.ox2 != 0));

  } else if (type == 3) {
    SetIndexRangesSBFF(ni.ox1, si, ei,
                       c_ivs, c_ive, c_ims, c_ipe,
                       ni.fi1, ni.fi2,
                       2,
                       true,
                       true);

    SetIndexRangesSBFF(ni.ox2, sj, ej,
                       c_jvs, c_jve, c_jms, c_jpe,
                       ni.fi1, ni.fi2,
                       2,
                       (pmb->block_size.nx2 > 1),
                       (ni.ox1 != 0));

    SetIndexRangesSBFF(ni.ox3, sk, ek,
                       c_kvs, c_kve, c_kms, c_kpe,
                       ni.fi1, ni.fi2,
                       2,
                       (pmb->block_size.nx3 > 1),
                       (ni.ox1 != 0 && ni.ox2 != 0));

  }

  return;
}

// For calculation of buffer sizes based on neighbor index information
int VertexCenteredBoundaryVariable::NeighborVariableBufferSize(const NeighborIndexes& ni) {
  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;
  int size = 0;

  idxLoadSameLevelRanges(ni, si, ei, sj, ej, sk, ek, false);
  AccumulateBufferSize(nl_, nu_, si, ei, sj, ej, sk, ek, size);

  if (pmy_mesh_->multilevel) {
    idxLoadSameLevelRanges(ni, si, ei, sj, ej, sk, ek, true);
    AccumulateBufferSize(nl_, nu_, si, ei, sj, ej, sk, ek, size);

    // as we are multi-level we should also maximize over fine to coarse and vice versa ops.
    int sizef2c = 0, sizef2c_dr = 0;
    idxLoadToCoarserRanges(ni, si, ei, sj, ej, sk, ek, false);
    AccumulateBufferSize(nl_, nu_, si, ei, sj, ej, sk, ek, sizef2c);
    idxLoadToCoarserRanges(ni, si, ei, sj, ej, sk, ek, true);
    // double restrict means spatial indices jump by two per iterate here
    AccumulateBufferSize(nl_, nu_, si, ei, sj, ej, sk, ek, sizef2c, 2);

    int sizec2f = 0;
    idxLoadToFinerRanges(ni, si, ei, sj, ej, sk, ek);
    AccumulateBufferSize(nl_, nu_, si, ei, sj, ej, sk, ek, sizec2f);

    size = std::max(size, sizef2c);
    size = std::max(size, sizec2f);
  }

  return size;
}

//----------------------------------------------------------------------------------------
// MPI buffer sizes
#ifdef MPI_PARALLEL

int VertexCenteredBoundaryVariable::MPI_BufferSizeSameLevel(
  const NeighborIndexes& ni,
  bool is_send=true) {

  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;
  int size = 0;

  if (is_send) {
    idxLoadSameLevelRanges(ni, si, ei, sj, ej, sk, ek, false);
    AccumulateBufferSize(nl_, nu_, si, ei, sj, ej, sk, ek, size);

    if (pmy_mesh_->multilevel) {
      idxLoadSameLevelRanges(ni, si, ei, sj, ej, sk, ek, true);
      AccumulateBufferSize(nl_, nu_, si, ei, sj, ej, sk, ek, size);
    }
  } else {
    idxSetSameLevelRanges(ni, si, ei, sj, ej, sk, ek, 1);
    AccumulateBufferSize(nl_, nu_, si, ei, sj, ej, sk, ek, size);

    if (pmy_mesh_->multilevel) {
      idxSetSameLevelRanges(ni, si, ei, sj, ej, sk, ek, 2);
      AccumulateBufferSize(nl_, nu_, si, ei, sj, ej, sk, ek, size);
    }
  }

  return size;
}

int VertexCenteredBoundaryVariable::MPI_BufferSizeToCoarser(
  const NeighborIndexes& ni) {

  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;
  int size = 0;

  idxLoadToCoarserRanges(ni, si, ei, sj, ej, sk, ek, false);
  AccumulateBufferSize(nl_, nu_, si, ei, sj, ej, sk, ek, size);
  idxLoadToCoarserRanges(ni, si, ei, sj, ej, sk, ek, true);
  // double restrict means spatial indices jump by two per iterate here
  AccumulateBufferSize(nl_, nu_, si, ei, sj, ej, sk, ek, size, 2);

  return size;
}

int VertexCenteredBoundaryVariable::MPI_BufferSizeFromCoarser(
  const NeighborIndexes& ni) {

  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;
  int size = 0;

  idxSetFromCoarserRanges(ni, si, ei, sj, ej, sk, ek, false);
  AccumulateBufferSize(nl_, nu_, si, ei, sj, ej, sk, ek, size);

  return size;
}

int VertexCenteredBoundaryVariable::MPI_BufferSizeToFiner(
  const NeighborIndexes& ni) {

  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;
  int size = 0;

  idxLoadToFinerRanges(ni, si, ei, sj, ej, sk, ek);
  AccumulateBufferSize(nl_, nu_, si, ei, sj, ej, sk, ek, size);

  return size;
}

int VertexCenteredBoundaryVariable::MPI_BufferSizeFromFiner(
  const NeighborIndexes& ni) {

  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;
  int size = 0;

  idxSetFromFinerRanges(ni, si, ei, sj, ej, sk, ek, 1);
  AccumulateBufferSize(nl_, nu_, si, ei, sj, ej, sk, ek, size);

  idxSetFromFinerRanges(ni, si, ei, sj, ej, sk, ek, 2);
  AccumulateBufferSize(nl_, nu_, si, ei, sj, ej, sk, ek, size);

  return size;
}

#endif