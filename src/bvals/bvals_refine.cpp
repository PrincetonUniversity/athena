//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_refine.cpp
//  \brief constructor/destructor and utility functions for BoundaryValues class

// C headers

// C++ headers
#include <algorithm>  // min
#include <cmath>
#include <iterator>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "bvals.hpp"
#include "cc/hydro/bvals_hydro.hpp"
#include "fc/bvals_fc.hpp"

void BoundaryValues::ProlongateBoundaries(const Real time, const Real dt) {
  MeshBlock *pmb = pmy_block_;
  int &mylevel = pmb->loc.level;

  // TODO(KGF): temporarily hardcode Hydro and Field array access for the below switch
  // around ApplyPhysicalBoundariesOnCoarseLevel()

  // This hardcoded technique is also used to manually specify the coupling between
  // physical variables in:
  // - step 2, ApplyPhysicalBoundariesOnCoarseLevel(): calls to W(U) and user BoundaryFunc
  // - step 3, ProlongateGhostCells(): calls to calculate bcc and U(W)

  // Additionally, pmr->SetHydroRefinement() is currently used in
  // RestrictGhostCellsOnSameLevel() (GR) and ProlongateGhostCells() (always) to switch
  // between conserved and primitive tuples, but this does not require ph, pf

  // downcast BoundaryVariable pointers to known derived class pointer types:
  // RTTI via dynamic_case
  HydroBoundaryVariable *phbvar =
      dynamic_cast<HydroBoundaryVariable *>(bvars_main_int[0]);
  Hydro *ph = pmb->phydro;

  FaceCenteredBoundaryVariable *pfbvar = nullptr;
  Field *pf = nullptr;
  if (MAGNETIC_FIELDS_ENABLED) {
    pf = pmb->pfield;
    pfbvar = dynamic_cast<FaceCenteredBoundaryVariable *>(bvars_main_int[1]);
  }

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
          RestrictGhostCellsOnSameLevel(nb, nk, nj, ni);
        }
      }
    }

    // calculate the loop limits for the ghost zones
    int cn = pmb->cnghost - 1;
    int si, ei, sj, ej, sk, ek;
    if (nb.ni.ox1 == 0) {
      std::int64_t &lx1 = pmb->loc.lx1;
      si = pmb->cis, ei = pmb->cie;
      if ((lx1 & 1LL) == 0LL) ei += cn;
      else             si -= cn;
    } else if (nb.ni.ox1 > 0) { si = pmb->cie + 1,  ei = pmb->cie + cn;}
    else              si = pmb->cis-cn, ei = pmb->cis-1;
    if (nb.ni.ox2 == 0) {
      sj = pmb->cjs, ej = pmb->cje;
      if (pmb->block_size.nx2 > 1) {
        std::int64_t &lx2 = pmb->loc.lx2;
        if ((lx2 & 1LL) == 0LL) ej += cn;
        else             sj -= cn;
      }
    } else if (nb.ni.ox2 > 0) { sj = pmb->cje + 1,  ej = pmb->cje + cn;}
    else              sj = pmb->cjs-cn, ej = pmb->cjs-1;
    if (nb.ni.ox3 == 0) {
      sk = pmb->cks, ek = pmb->cke;
      if (pmb->block_size.nx3 > 1) {
        std::int64_t &lx3 = pmb->loc.lx3;
        if ((lx3 & 1LL) == 0LL) ek += cn;
        else             sk -= cn;
      }
    } else if (nb.ni.ox3 > 0) { sk = pmb->cke + 1,  ek = pmb->cke + cn;}
    else              sk = pmb->cks-cn, ek = pmb->cks-1;

    // (temp workaround) to automatically call all BoundaryFunction_[] on coarse_prim/b
    // instead of default targets var_cc=cons or prim + var_fc=b, exchange the pointers
    phbvar->var_cc = &(ph->coarse_prim_);
    if (MAGNETIC_FIELDS_ENABLED)
      pfbvar->var_fc = &(pf->coarse_b_);

    // Step 2. Re-apply physical boundaries on the coarse boundary:
    ApplyPhysicalBoundariesOnCoarseLevel(nb, time, dt, si, ei, sj, ej, sk, ek);

    // (temp workaround) swap BoundaryVariable var_cc/fc back from coarse_buf
    phbvar->var_cc = &(ph->w);
    if (MAGNETIC_FIELDS_ENABLED)
      pfbvar->var_fc = &(pf->b);

    // Step 3. Finally, the ghost-ghost zones are ready for prolongation:
    ProlongateGhostCells(nb, si, ei, sj, ej, sk, ek);
  } // end loop over nneighbor
  return;
}


void BoundaryValues::RestrictGhostCellsOnSameLevel(const NeighborBlock& nb, int nk,
                                                   int nj, int ni) {
  MeshBlock *pmb = pmy_block_;
  MeshRefinement *pmr = pmb->pmr;

  int ris, rie, rjs, rje, rks, rke;
  if (ni == 0) {
    ris = pmb->cis;
    rie = pmb->cie;
    if (nb.ni.ox1 == 1) {
      ris = pmb->cie;
    } else if (nb.ni.ox1 == -1) {
      rie = pmb->cis;
    }
  } else if (ni == 1) {
    ris = pmb->cie + 1, rie = pmb->cie + 1;
  } else { //(ni ==  - 1)
    ris = pmb->cis - 1, rie = pmb->cis - 1;
  }
  if (nj == 0) {
    rjs = pmb->cjs, rje = pmb->cje;
    if (nb.ni.ox2 == 1) rjs = pmb->cje;
    else if (nb.ni.ox2 == -1) rje = pmb->cjs;
  } else if (nj == 1) {
    rjs = pmb->cje + 1, rje = pmb->cje + 1;
  } else { //(nj == -1)
    rjs = pmb->cjs - 1, rje = pmb->cjs - 1;
  }
  if (nk == 0) {
    rks = pmb->cks, rke = pmb->cke;
    if (nb.ni.ox3 == 1) rks = pmb->cke;
    else if (nb.ni.ox3 == -1) rke = pmb->cks;
  } else if (nk == 1) {
    rks = pmb->cke + 1, rke = pmb->cke + 1;
  } else { //(nk == -1)
    rks = pmb->cks - 1, rke = pmb->cks - 1;
  }

  for (auto cc_pair : pmr->pvars_cc_) {
    AthenaArray<Real> *var_cc = std::get<0>(cc_pair);
    AthenaArray<Real> *coarse_cc = std::get<1>(cc_pair);
    int nu = var_cc->GetDim4() - 1;
    pmb->pmr->RestrictCellCenteredValues(*var_cc, *coarse_cc, 0, nu,
                                         ris, rie, rjs, rje, rks, rke);
  }
  // (unique to Hydro) also restrict primitive values in ghost zones when GR + multilevel
  if (GENERAL_RELATIVITY) {
    pmr->SetHydroRefinement(HydroBoundaryQuantity::prim);
    auto cc_pair = pmr->pvars_cc_.front();
    AthenaArray<Real> *var_cc = std::get<0>(cc_pair);
    AthenaArray<Real> *coarse_cc = std::get<1>(cc_pair);
    int nu = var_cc->GetDim4() - 1;
    pmb->pmr->RestrictCellCenteredValues(*var_cc, *coarse_cc, 0, nu,
                                         ris, rie, rjs, rje, rks, rke);
    pmr->SetHydroRefinement(HydroBoundaryQuantity::cons);
  }

  for (auto fc_pair : pmr->pvars_fc_) {
    FaceField *var_fc = std::get<0>(fc_pair);
    FaceField *coarse_fc = std::get<1>(fc_pair);
    int &mylevel = pmb->loc.level;
    int rs = ris, re = rie + 1;
    if (rs == pmb->cis   && nblevel[nk+1][nj+1][ni  ] < mylevel) rs++;
    if (re == pmb->cie+1 && nblevel[nk+1][nj+1][ni+2] < mylevel) re--;
    pmr->RestrictFieldX1((*var_fc).x1f, (*coarse_fc).x1f, rs, re, rjs, rje, rks,
                         rke);
    if (pmb->block_size.nx2 > 1) {
      rs = rjs, re = rje + 1;
      if (rs == pmb->cjs   && nblevel[nk+1][nj  ][ni+1] < mylevel) rs++;
      if (re == pmb->cje+1 && nblevel[nk+1][nj+2][ni+1] < mylevel) re--;
      pmr->RestrictFieldX2((*var_fc).x2f, (*coarse_fc).x2f, ris, rie, rs, re, rks,
                           rke);
    } else { // 1D
      pmr->RestrictFieldX2((*var_fc).x2f, (*coarse_fc).x2f, ris, rie, rjs, rje, rks,
                           rke);
      for (int i=ris; i<=rie; i++)
        (*coarse_fc).x2f(rks,rjs+1,i) = (*coarse_fc).x2f(rks,rjs,i);
    }
    if (pmb->block_size.nx3 > 1) {
      rs = rks, re =  rke + 1;
      if (rs == pmb->cks   && nblevel[nk  ][nj+1][ni+1] < mylevel) rs++;
      if (re == pmb->cke+1 && nblevel[nk+2][nj+1][ni+1] < mylevel) re--;
      pmr->RestrictFieldX3((*var_fc).x3f, (*coarse_fc).x3f, ris, rie, rjs, rje, rs,
                           re);
    } else { // 1D or 2D
      pmr->RestrictFieldX3((*var_fc).x3f, (*coarse_fc).x3f, ris, rie, rjs, rje, rks,
                           rke);
      for (int j=rjs; j<=rje; j++) {
        for (int i=ris; i<=rie; i++)
          (*coarse_fc).x3f(rks+1,j,i) = (*coarse_fc).x3f(rks,j,i);
      }
    }
  } // end loop over pvars_fc_
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ApplyPhysicalBoundariesOnCoarseLevel(
//           const NeighborBlock& nb, const Real time, const Real dt,
//           int si, int ei, int sj, int ej, int sk, int ek)
//  \brief

void BoundaryValues::ApplyPhysicalBoundariesOnCoarseLevel(
    const NeighborBlock& nb, const Real time, const Real dt,
    int si, int ei, int sj, int ej, int sk, int ek) {
  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;
  MeshRefinement *pmr = pmb->pmr;

  // temporarily hardcode Hydro and Field array access:
  Hydro *ph = pmb->phydro;
  Field *pf = nullptr;
  if (MAGNETIC_FIELDS_ENABLED) {
    pf = pmb->pfield;
  }

  // convert the ghost zone and ghost-ghost zones into primitive variables
  // this includes cell-centered field calculation
  int f1m = 0, f1p = 0, f2m = 0, f2p = 0, f3m = 0, f3p = 0;
  if (nb.ni.ox1 == 0) {
    if (nblevel[1][1][0] != -1) f1m = 1;
    if (nblevel[1][1][2] != -1) f1p = 1;
  } else {
    f1m = 1;
    f1p = 1;
  }
  if (pmb->block_size.nx2 > 1) {
    if (nb.ni.ox2 == 0) {
      if (nblevel[1][0][1] != -1) f2m = 1;
      if (nblevel[1][2][1] != -1) f2p = 1;
    } else {
      f2m = 1;
      f2p = 1;
    }
  }
  if (pmb->block_size.nx3 > 1) {
    if (nb.ni.ox3 == 0) {
      if (nblevel[0][1][1] != -1) f3m = 1;
      if (nblevel[2][1][1] != -1) f3p = 1;
    } else {
      f3m = 1;
      f3p = 1;
    }
  }

  // TODO(KGF): passing nullptrs (pf) if no MHD (coarse_* no longer in MeshRefinement)

  // KGF: COUPLING OF QUANTITIES (must be manually specified)
  pmb->peos->ConservedToPrimitive(ph->coarse_cons_, ph->coarse_prim_,
                                  pf->coarse_b_, ph->coarse_prim_,
                                  pf->coarse_bcc_, pmr->pcoarsec,
                                  si-f1m, ei+f1p, sj-f2m, ej+f2p, sk-f3m, ek+f3p);

  if (nb.ni.ox1 == 0) {
    if (apply_bndry_fn_[BoundaryFace::inner_x1]) {
      switch(block_bcs[BoundaryFace::inner_x1]) {
        case BoundaryFlag::reflect:
          for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
               ++bvars_it) {
            (*bvars_it)->ReflectInnerX1(pmb, pco, time, dt, pmb->cis, pmb->cie,
                                        sj, ej, sk, ek, 1);
          }
          break;
        case BoundaryFlag::outflow:
          for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
               ++bvars_it) {
            (*bvars_it)->OutflowInnerX1(pmb, pco, time, dt, pmb->cis, pmb->cie,
                                        sj, ej, sk, ek, 1);
          }
          break;
        case BoundaryFlag::user: // user-enrolled BCs
          pmy_mesh_->BoundaryFunction_[BoundaryFace::inner_x1](
              pmb, pmr->pcoarsec, ph->coarse_prim_, pf->coarse_b_, time, dt,
              pmb->cis, pmb->cie, sj, ej, sk, ek, 1);
          break;
        default:
          break;
      }
    }
    if (apply_bndry_fn_[BoundaryFace::outer_x1]) {
      switch(block_bcs[BoundaryFace::outer_x1]) {
        case BoundaryFlag::reflect:
          for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
               ++bvars_it) {
            (*bvars_it)->ReflectOuterX1(pmb, pco, time, dt, pmb->cis, pmb->cie,
                                        sj, ej, sk, ek, 1);
          }
          break;
        case BoundaryFlag::outflow:
          for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
               ++bvars_it) {
            (*bvars_it)->OutflowOuterX1(pmb, pco, time, dt, pmb->cis, pmb->cie,
                                        sj, ej, sk, ek, 1);
          }
          break;
        case BoundaryFlag::user: // user-enrolled BCs
          pmy_mesh_->BoundaryFunction_[BoundaryFace::outer_x1](
              pmb, pmr->pcoarsec, ph->coarse_prim_, pf->coarse_b_, time, dt,
              pmb->cis, pmb->cie, sj, ej, sk, ek, 1);
          break;
        default:
          break;
      }
    }
  }
  if (nb.ni.ox2 == 0 && pmb->block_size.nx2 > 1) {
    if (apply_bndry_fn_[BoundaryFace::inner_x2]) {
      switch(block_bcs[BoundaryFace::inner_x2]) {
        case BoundaryFlag::reflect:
          for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
               ++bvars_it) {
            (*bvars_it)->ReflectInnerX2(pmb, pco, time, dt, si, ei,
                                        pmb->cjs, pmb->cje, sk, ek, 1);
          }
          break;
        case BoundaryFlag::outflow:
          for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
               ++bvars_it) {
            (*bvars_it)->OutflowInnerX2(pmb, pco, time, dt, si, ei,
                                        pmb->cjs, pmb->cje, sk, ek, 1);
          }
          break;
        case BoundaryFlag::polar_wedge:
          for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
               ++bvars_it) {
            (*bvars_it)->PolarWedgeInnerX2(pmb, pco, time, dt, si, ei,
                                           pmb->cjs, pmb->cje, sk, ek, 1);
          }
          break;
        case BoundaryFlag::user: // user-enrolled BCs
          pmy_mesh_->BoundaryFunction_[BoundaryFace::inner_x2](
              pmb, pmr->pcoarsec, ph->coarse_prim_, pf->coarse_b_, time, dt,
              si, ei, pmb->cjs, pmb->cje, sk, ek, 1);
          break;
        default:
          break;
      }
    }
    if (apply_bndry_fn_[BoundaryFace::outer_x2]) {
      switch(block_bcs[BoundaryFace::outer_x2]) {
        case BoundaryFlag::reflect:
          for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
               ++bvars_it) {
            (*bvars_it)->ReflectOuterX2(pmb, pco, time, dt, si, ei,
                                        pmb->cjs, pmb->cje, sk, ek, 1);
          }
          break;
        case BoundaryFlag::outflow:
          for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
               ++bvars_it) {
            (*bvars_it)->OutflowOuterX2(pmb, pco, time, dt, si, ei,
                                        pmb->cjs, pmb->cje, sk, ek, 1);
          }
          break;
        case BoundaryFlag::polar_wedge:
          for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
               ++bvars_it) {
            (*bvars_it)->PolarWedgeOuterX2(pmb, pco, time, dt, si, ei,
                                           pmb->cjs, pmb->cje, sk, ek, 1);
          }
          break;
        case BoundaryFlag::user: // user-enrolled BCs
          pmy_mesh_->BoundaryFunction_[BoundaryFace::outer_x2](
              pmb, pmr->pcoarsec, ph->coarse_prim_, pf->coarse_b_, time, dt,
              si, ei, pmb->cjs, pmb->cje, sk, ek, 1);
          break;
        default:
          break;
      }
    }
  }
  if (nb.ni.ox3 == 0 && pmb->block_size.nx3 > 1) {
    if (apply_bndry_fn_[BoundaryFace::inner_x3]) {
      switch(block_bcs[BoundaryFace::inner_x3]) {
        case BoundaryFlag::reflect:
          for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
               ++bvars_it) {
            (*bvars_it)->ReflectInnerX3(pmb, pco, time, dt, si, ei,
                                        sj, ej, pmb->cks, pmb->cke, 1);
          }
          break;
        case BoundaryFlag::outflow:
          for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
               ++bvars_it) {
            (*bvars_it)->OutflowInnerX3(pmb, pco, time, dt, si, ei,
                                        sj, ej, pmb->cks, pmb->cke, 1);
          }
          break;
        case BoundaryFlag::user: // user-enrolled BCs
          pmy_mesh_->BoundaryFunction_[BoundaryFace::inner_x3](
              pmb, pmr->pcoarsec, ph->coarse_prim_, pf->coarse_b_, time, dt,
              si, ei, sj, ej, pmb->cks, pmb->cke, 1);
          break;
        default:
          break;
      }
    }
    if (apply_bndry_fn_[BoundaryFace::outer_x3]) {
      switch(block_bcs[BoundaryFace::outer_x3]) {
        case BoundaryFlag::reflect:
          for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
               ++bvars_it) {
            (*bvars_it)->ReflectOuterX3(pmb, pco, time, dt, si, ei,
                                        sj, ej, pmb->cks, pmb->cke, 1);
          }
          break;
        case BoundaryFlag::outflow:
          for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
               ++bvars_it) {
            (*bvars_it)->OutflowOuterX3(pmb, pco, time, dt, si, ei,
                                        sj, ej, pmb->cks, pmb->cke, 1);
          }
          break;
        case BoundaryFlag::user: // user-enrolled BCs
          pmy_mesh_->BoundaryFunction_[BoundaryFace::outer_x3](
              pmb, pmr->pcoarsec, ph->coarse_prim_, pf->coarse_b_, time, dt,
              si, ei, sj, ej, pmb->cks, pmb->cke, 1);
          break;
        default:
          break;
      }
    }
  }
  return;
}

void BoundaryValues::ProlongateGhostCells(const NeighborBlock& nb,
                                          int si, int ei, int sj, int ej,
                                          int sk, int ek) {
  MeshBlock *pmb = pmy_block_;
  MeshRefinement *pmr = pmb->pmr;

  // prolongate cell-centered S/AMR-enrolled quantities (hydro, passive scalars, ...)
  // (unique to Hydro) swap (u, coarse_cons_) with (w, coarse_prim)
  pmr->SetHydroRefinement(HydroBoundaryQuantity::prim);
  for (auto cc_pair : pmr->pvars_cc_) {
    AthenaArray<Real> *var_cc = std::get<0>(cc_pair);
    AthenaArray<Real> *coarse_cc = std::get<1>(cc_pair);
    int nu = var_cc->GetDim4() - 1;
    pmr->ProlongateCellCenteredValues(*coarse_cc, *var_cc, 0, nu,
                                      si, ei, sj, ej, sk, ek);
  }
  // swap back Hydro refinement quantities:
  pmr->SetHydroRefinement(HydroBoundaryQuantity::cons);


  // prolongate face-centered S/AMR-enrolled quantities (magnetic fields)
  int &mylevel = pmb->loc.level;
  int il, iu, jl, ju, kl, ku;
  il = si, iu = ei + 1;
  if ((nb.ni.ox1 >= 0) && (nblevel[nb.ni.ox3+1][nb.ni.ox2+1][nb.ni.ox1  ] >= mylevel))
    il++;
  if ((nb.ni.ox1 <= 0) && (nblevel[nb.ni.ox3+1][nb.ni.ox2+1][nb.ni.ox1+2] >= mylevel))
    iu--;
  if (pmb->block_size.nx2 > 1) {
    jl = sj, ju = ej + 1;
    if ((nb.ni.ox2 >= 0) && (nblevel[nb.ni.ox3+1][nb.ni.ox2  ][nb.ni.ox1+1] >= mylevel))
      jl++;
    if ((nb.ni.ox2 <= 0) && (nblevel[nb.ni.ox3+1][nb.ni.ox2+2][nb.ni.ox1+1] >= mylevel))
      ju--;
  } else {
    jl = sj;
    ju = ej;
  }
  if (pmb->block_size.nx3 > 1) {
    kl = sk, ku = ek + 1;
    if ((nb.ni.ox3 >= 0) && (nblevel[nb.ni.ox3  ][nb.ni.ox2+1][nb.ni.ox1+1] >= mylevel))
      kl++;
    if ((nb.ni.ox3 <= 0) && (nblevel[nb.ni.ox3+2][nb.ni.ox2+1][nb.ni.ox1+1] >= mylevel))
      ku--;
  } else {
    kl = sk;
    ku = ek;
  }
  for (auto fc_pair : pmr->pvars_fc_) {
    FaceField *var_fc = std::get<0>(fc_pair);
    FaceField *coarse_fc = std::get<1>(fc_pair);

    // step 1. calculate x1 outer surface fields and slopes
    pmr->ProlongateSharedFieldX1((*coarse_fc).x1f, (*var_fc).x1f, il, iu, sj, ej, sk, ek);
    // step 2. calculate x2 outer surface fields and slopes
    pmr->ProlongateSharedFieldX2((*coarse_fc).x2f, (*var_fc).x2f, si, ei, jl, ju, sk, ek);
    // step 3. calculate x3 outer surface fields and slopes
    pmr->ProlongateSharedFieldX3((*coarse_fc).x3f, (*var_fc).x3f, si, ei, sj, ej, kl, ku);

    // step 4. calculate the internal finer fields using the Toth & Roe method
    pmr->ProlongateInternalField((*var_fc), si, ei, sj, ej, sk, ek);
  }

  // now that the ghost-ghost zones are filled and prolongated,
  // calculate the loop limits for the finer grid
  int fsi, fei, fsj, fej, fsk, fek;
  fsi = (si - pmb->cis)*2 + pmb->is;
  fei = (ei - pmb->cis)*2 + pmb->is + 1;
  if (pmb->block_size.nx2 > 1) {
    fsj = (sj - pmb->cjs)*2 + pmb->js;
    fej = (ej - pmb->cjs)*2 + pmb->js + 1;
  } else {
    fsj = pmb->js;
    fej = pmb->je;
  }
  if (pmb->block_size.nx3 > 1) {
    fsk = (sk - pmb->cks)*2 + pmb->ks;
    fek = (ek - pmb->cks)*2 + pmb->ks + 1;
  } else {
    fsk = pmb->ks;
    fek = pmb->ke;
  }

  // temporarily hardcode Hydro and Field array access
  Hydro *ph = pmb->phydro;
  Field *pf = nullptr;

  // KGF: COUPLING OF QUANTITIES (must be manually specified)
  // Field prolongation completed, calculate cell centered fields
  if (MAGNETIC_FIELDS_ENABLED) {
    pf = pmb->pfield;
    pf->CalculateCellCenteredField(pf->b, pf->bcc, pmb->pcoord,
                                   fsi, fei, fsj, fej, fsk, fek);
  }
  // TODO(KGF): passing nullptrs (pf) if no MHD (coarse_* no longer in MeshRefinement)
  // KGF: COUPLING OF QUANTITIES (must be manually specified)
  // calculate conservative variables
  pmb->peos->PrimitiveToConserved(ph->w, pf->bcc, ph->u, pmb->pcoord,
                                  fsi, fei, fsj, fej, fsk, fek);
  return;
}
