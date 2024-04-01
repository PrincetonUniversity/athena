//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_refine_gravity.cpp
//! \brief unctions for prolongation at level boundaries for gravity

// C headers

// C++ headers
#include <algorithm>  // min
#include <cmath>
#include <iterator>
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "bvals.hpp"


class CellCenteredCellCenteredBoundaryVariable;

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ProlongateBoundariesPostMG(
//!                                              CellCenteredBoundaryVariable* pbvar)
//! \brief Prolongate boundaries after Multigrid assuming physical boundaries are filled

void BoundaryValues::ProlongateBoundariesPostMG(CellCenteredBoundaryVariable* pbvar) {
  MeshBlock *pmb = pmy_block_;
  const int &mylevel = loc.level;

  // For each finer neighbor, to prolongate a boundary we need to fill one more cell
  // surrounding the boundary zone to calculate the slopes ("ghost-ghost zone"). 3x steps:
  for (int n=0; n<nneighbor; n++) {
    NeighborBlock& nb = neighbor[n];
    if (nb.snb.level >= mylevel) continue;
    // fill the required ghost-ghost zone
    int nis, nie, njs, nje, nks, nke;
    nis = std::max(nb.ni.ox1-1, -1);
    nie = std::min(nb.ni.ox1+1, 1);
    if (pmb->block_size.nx2 == 1) {
      njs = 0;
      nje = 0;
    } else {
      njs = std::max(nb.ni.ox2-1, -1);
      nje = std::min(nb.ni.ox2+1, 1);
    }

    if (pmb->block_size.nx3 == 1) {
      nks = 0;
      nke = 0;
    } else {
      nks = std::max(nb.ni.ox3-1, -1);
      nke = std::min(nb.ni.ox3+1, 1);
    }

    // Step 1. Apply necessary variable restrictions when ghost-ghost zone is on same lvl
    for (int nk=nks; nk<=nke; nk++) {
      for (int nj=njs; nj<=nje; nj++) {
        for (int ni=nis; ni<=nie; ni++) {
          int ntype = std::abs(ni) + std::abs(nj) + std::abs(nk);
          // skip myself or coarse levels; only the same level must be restricted
          if (ntype == 0 || nblevel[nk+1][nj+1][ni+1] != mylevel) continue;

          // this neighbor block is on the same level
          // and needs to be restricted for prolongation
          RestrictGhostCellsOnSameLevelPostMG(pbvar, nb, nk, nj, ni);
        }
      }
    }

    // calculate the loop limits for the ghost zones
    int cn = pmb->cnghost - 1;
    int si, ei, sj, ej, sk, ek;
    if (nb.ni.ox1 == 0) {
      std::int64_t &lx1 = loc.lx1;
      si = pmb->cis, ei = pmb->cie;
      if ((lx1 & 1LL) == 0LL) ei += cn;
      else             si -= cn;
    } else if (nb.ni.ox1 > 0) { si = pmb->cie + 1,  ei = pmb->cie + cn;}
    else              si = pmb->cis-cn, ei = pmb->cis-1;
    if (nb.ni.ox2 == 0) {
      sj = pmb->cjs, ej = pmb->cje;
      if (pmb->block_size.nx2 > 1) {
        std::int64_t &lx2 = loc.lx2;
        if ((lx2 & 1LL) == 0LL) ej += cn;
        else             sj -= cn;
      }
    } else if (nb.ni.ox2 > 0) { sj = pmb->cje + 1,  ej = pmb->cje + cn;}
    else              sj = pmb->cjs-cn, ej = pmb->cjs-1;
    if (nb.ni.ox3 == 0) {
      sk = pmb->cks, ek = pmb->cke;
      if (pmb->block_size.nx3 > 1) {
        std::int64_t &lx3 = loc.lx3;
        if ((lx3 & 1LL) == 0LL) ek += cn;
        else             sk -= cn;
      }
    } else if (nb.ni.ox3 > 0) { sk = pmb->cke + 1,  ek = pmb->cke + cn;}
    else              sk = pmb->cks-cn, ek = pmb->cks-1;


    // (Step 2. skip physical boundaries - assuming they are already filled


    // Step 3. Finally, the ghost-ghost zones are ready for prolongation:
    ProlongateGhostCellsPostMG(pbvar, si, ei, sj, ej, sk, ek);
  } // end loop over nneighbor
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::RestrictGhostCellsOnSameLevelPostMG(
//!                          CellCenteredBoundaryVariable* pbvar, const NeighborBlock& nb,
//!                          int nk, int nj, int ni)
//! \brief Restrict ghost cells on same level after Multigrid
void BoundaryValues::RestrictGhostCellsOnSameLevelPostMG(
                             CellCenteredBoundaryVariable* pbvar, const NeighborBlock& nb,
                             int nk, int nj, int ni) {
  MeshBlock *pmb = pmy_block_;
  MeshRefinement *pmr = pmb->pmr;

  int ris, rie, rjs, rje, rks, rke;
  if (ni == 0) {
    ris = pmb->cis, rie = pmb->cie;
    if (nb.ni.ox1 == 1)       ris = pmb->cie;
    else if (nb.ni.ox1 == -1) rie = pmb->cis;
  } else if (ni == 1) {
    ris = pmb->cie + 1, rie = pmb->cie + 1;
  } else { // (ni == -1)
    ris = pmb->cis - 1, rie = pmb->cis - 1;
  }
  if (nj == 0) {
    rjs = pmb->cjs, rje = pmb->cje;
    if (nb.ni.ox2 == 1)       rjs = pmb->cje;
    else if (nb.ni.ox2 == -1) rje = pmb->cjs;
  } else if (nj == 1) {
    rjs = pmb->cje + 1, rje = pmb->cje + 1;
  } else { // (nj == -1)
    rjs = pmb->cjs - 1, rje = pmb->cjs - 1;
  }
  if (nk == 0) {
    rks = pmb->cks, rke = pmb->cke;
    if (nb.ni.ox3 == 1)       rks = pmb->cke;
    else if (nb.ni.ox3 == -1) rke = pmb->cks;
  } else if (nk == 1) {
    rks = pmb->cke + 1, rke = pmb->cke + 1;
  } else { // (nk == -1)
    rks = pmb->cks - 1, rke = pmb->cks - 1;
  }

  pmb->pmr->RestrictCellCenteredValues(*(pbvar->var_cc), *(pbvar->coarse_buf), 0, 0,
                                       ris, rie, rjs, rje, rks, rke);

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ProlongateGhostCellsPostMG(
//!                                    CellCenteredBoundaryVariable* pbvar,
//!                                    int si, int ei, int sj, int ej, int sk, int ek)
//! \brief Prolongate ghost cells after Multigrid

void BoundaryValues::ProlongateGhostCellsPostMG(CellCenteredBoundaryVariable* pbvar,
                               int si, int ei, int sj, int ej, int sk, int ek) {
  MeshRefinement *pmr = pmy_block_->pmr;
  pmr->ProlongateCellCenteredValues(*(pbvar->coarse_buf), *(pbvar->var_cc), 0, 0,
                                    si, ei, sj, ej, sk, ek);
  return;
}
