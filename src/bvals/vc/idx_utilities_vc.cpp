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
//! \fn void VertexCenteredBoundaryVariable::AccumulateBufferSize(...)
//  \brief Buffer size accumulator
void VertexCenteredBoundaryVariable::AccumulateBufferSize(
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
