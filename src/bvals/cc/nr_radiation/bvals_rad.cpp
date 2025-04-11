//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_rad.cpp
//  \brief implements boundary functions for Radiation variables and utilities to manage
//  primitive/conservative variable relationship in a derived class of the
//  CellCenteredBoundaryVariable base class.

// C headers

// C++ headers
#include <iostream>   // endl
#include <sstream>    // stringstream

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../coordinates/coordinates.hpp"
#include "../../../field/field.hpp"
#include "../../../globals.hpp"
#include "../../../hydro/hydro.hpp"
#include "../../../mesh/mesh.hpp"
#include "../../../nr_radiation/radiation.hpp"
#include "bvals_rad.hpp"

//----------------------------------------------------------------------------------------
//! \class RadiationBoundaryFunctions

RadBoundaryVariable::RadBoundaryVariable(MeshBlock *pmb,
    AthenaArray<Real> *var_rad, AthenaArray<Real> *coarse_var,
    AthenaArray<Real> *var_flux) :
    CellCenteredBoundaryVariable(pmb, var_rad, coarse_var, var_flux, true, 1) {
    // the radiation array is (k,j,i,n)
    // the number of variables is GetDim1
    // All the other shared functions are initialized in CellCenteredBoundaryVariable
    // radiation specific functions need to be defined here
    azimuthal_shift_rad_.NewAthenaArray(pmb->ke + NGHOST + 2,nu_+1);
    ir_cm_.NewAthenaArray(pmb->nfre_ang);
    ir_lab_.NewAthenaArray(pmb->nfre_ang);
    if (pbval_->shearing_box != 0) {
      int pnum = pmb->block_size.nx2+2*NGHOST+1;
      pflux_.NewAthenaArray(pnum,pmb->nfre_ang);
    }
}

//----------------------------------------------------------------------------------------
//! \fn int RadBoundaryVariable::LoadBoundaryBufferSameLevel(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set cell-centered boundary buffers for sending to a block on the same level

int RadBoundaryVariable::LoadBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb) {
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
  BufferUtility::PackData(var, buf, sk, ek, nl_, nu_, si, ei, sj, ej, p);

  return p;
}

//----------------------------------------------------------------------------------------
//! \fn int RadBoundaryVariable::LoadBoundaryBufferToCoarser(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set cell-centered boundary buffers for sending to a block on the coarser level

int RadBoundaryVariable::LoadBoundaryBufferToCoarser(Real *buf, const NeighborBlock& nb) {
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
  // function overload to do restriction for radiation variabbles
  pmr->RestrictCellCenteredValues(var, coarse_var, -1, nl_, nu_, si, ei, sj, ej, sk, ek);
  BufferUtility::PackData(coarse_var, buf, sk, ek, nl_, nu_, si, ei, sj, ej, p);
  return p;
}


//----------------------------------------------------------------------------------------
//! \fn int CellCenteredBoundaryVariable::LoadBoundaryBufferToFiner(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set cell-centered boundary buffers for sending to a block on the finer level

int RadBoundaryVariable::LoadBoundaryBufferToFiner(Real *buf, const NeighborBlock& nb) {
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
  BufferUtility::PackData(var, buf, sk, ek, nl_, nu_, si, ei, sj, ej, p);
  return p;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryVariable::SetBoundaries()
//  \brief set the boundary data

void RadBoundaryVariable::SetBoundaries() {
  MeshBlock *pmb = pmy_block_;
  int mylevel = pmb->loc.level;
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (nb.snb.level == mylevel)
      SetBoundarySameLevel(bd_var_.recv[nb.bufid], nb);
    else if (nb.snb.level < mylevel) // only sets the prolongation buffer
      SetBoundaryFromCoarser(bd_var_.recv[nb.bufid], nb);
    else
      SetBoundaryFromFiner(bd_var_.recv[nb.bufid], nb);
    bd_var_.flag[nb.bufid] = BoundaryStatus::completed; // completed
  }

  if (pbval_->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar ||
      pbval_->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar)
    PolarBoundarySingleAzimuthalBlock();

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::SetBoundarySameLevel(Real *buf,
//                                                              const NeighborBlock& nb)
//  \brief Set cell-centered boundary received from a block on the same level

void RadBoundaryVariable::SetBoundarySameLevel(Real *buf,
                                                        const NeighborBlock& nb) {
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

  int p = 0;
  // no need to flip for radiation
  if (nb.polar) {
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
        for (int i=si; i<=ei; ++i) {
#pragma omp simd linear(p)
          for (int n=nl_; n<=nu_; ++n) {
            var(k,j,i,n) = buf[p++];
          }
        }
      }
    }
  } else {
    BufferUtility::UnpackData(buf, var, sk, ek, nl_, nu_, si, ei, sj, ej, p);
  }
  return;
}



//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::SetBoundaryFromCoarser(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set cell-centered prolongation buffer received from a block on a coarser level

void RadBoundaryVariable::SetBoundaryFromCoarser(Real *buf,
                                                          const NeighborBlock& nb) {
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

  int p = 0;
  if (nb.polar) {
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
        for (int i=si; i<=ei; ++i) {
#pragma omp simd linear(p)
          for (int n=nl_; n<=nu_; ++n) {
            coarse_var(k,j,i,n) = buf[p++];
          }
        }
      }
    }
  } else {
    BufferUtility::UnpackData(buf, coarse_var, sk, ek, nl_, nu_, si, ei, sj, ej, p);
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::SetBoundaryFromFiner(Real *buf,
//                                                              const NeighborBlock& nb)
//  \brief Set cell-centered boundary received from a block on a finer level

void RadBoundaryVariable::SetBoundaryFromFiner(Real *buf, const NeighborBlock& nb) {
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

  int p = 0;
  if (nb.polar) {
//     Real sign=1.0;
//      if (flip_across_pole_ != nullptr) sign = flip_across_pole_[n] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
        for (int i=si; i<=ei; ++i) {
#pragma omp simd linear(p)
          for (int n=nl_; n<=nu_; ++n) {
            var(k,j,i,n) = buf[p++];
          }
        }
      }
    }
  } else {
    BufferUtility::UnpackData(buf, var,  sk, ek, nl_, nu_, si, ei, sj, ej, p);
  }
  return;
}




//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::PolarBoundarySingleAzimuthalBlock()
// \brief polar boundary edge-case: single MeshBlock spans the entire azimuthal (x3) range

void RadBoundaryVariable::PolarBoundarySingleAzimuthalBlock() {
  MeshBlock *pmb = pmy_block_;

  if (pmb->loc.level  ==  pmy_mesh_->root_level && pmy_mesh_->nrbx3 == 1
      && pmb->block_size.nx3 > 1) {
    AthenaArray<Real> &var = *var_cc;
    if (pbval_->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar) {
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      for (int j=pmb->js-NGHOST; j<=pmb->js-1; ++j) {
        for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i) {
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
            for (int n=nl_; n<=nu_; ++n) {
              azimuthal_shift_rad_(k,n) = var(k,j,i,n);
            }
          }
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
            int k_shift = k;
            k_shift += (k < (nx3_half + NGHOST) ? 1 : -1) * nx3_half;
            for (int n=nl_; n<=nu_; ++n) {
              var(k,j,i,n) = azimuthal_shift_rad_(k_shift,n);
            }
          }
        }
      }
    }

    if (pbval_->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar) {
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      for (int j=pmb->je+1; j<=pmb->je+NGHOST; ++j) {
        for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i) {
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
            for (int n=nl_; n<=nu_; ++n) {
              azimuthal_shift_rad_(k,n) = var(k,j,i,n);
            }
          }
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
            int k_shift = k;
            k_shift += (k < (nx3_half + NGHOST) ? 1 : -1) * nx3_half;
            for(int n=nl_; n<=nu_; ++n) {
              var(k,j,i,n) = azimuthal_shift_rad_(k_shift,n);
            }
          }
        }
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::ApplyRadPhysicalBoundaries()
// \brief Apply physical boundaries for specific intensities only
void RadBoundaryVariable::ApplyRadPhysicalBoundaries(const Real time, const Real dt) {
  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;
  BoundaryValues *pbval = pmb->pbval;
  int bis = pmb->is - NGHOST, bie = pmb->ie + NGHOST,
      bjs = pmb->js, bje = pmb->je,
      bks = pmb->ks, bke = pmb->ke;

  // Extend the transverse limits that correspond to periodic boundaries as they are
  // updated: x1, then x2, then x3
  if (!pbval->apply_bndry_fn_[BoundaryFace::inner_x2] && pmb->block_size.nx2 > 1)
    bjs = pmb->js - NGHOST;
  if (!pbval->apply_bndry_fn_[BoundaryFace::outer_x2] && pmb->block_size.nx2 > 1)
    bje = pmb->je + NGHOST;
  if (!pbval->apply_bndry_fn_[BoundaryFace::inner_x3] && pmb->block_size.nx3 > 1)
    bks = pmb->ks - NGHOST;
  if (!pbval->apply_bndry_fn_[BoundaryFace::outer_x3] && pmb->block_size.nx3 > 1)
    bke = pmb->ke + NGHOST;

  Hydro *ph = pmb->phydro;

  Field *pf = nullptr;
  if (MAGNETIC_FIELDS_ENABLED) {
    pf = pmb->pfield;
  }

  NRRadiation *prad = pmb->pnrrad;
  RadBoundaryVariable *pradbvar = &(prad->rad_bvar);


  // Apply boundary function on inner-x1
  if (pbval->apply_bndry_fn_[BoundaryFace::inner_x1]) {
    SetRadPhysicalFunctions(pmb, pco, time, dt,
                            pmb->is, pmb->ie, bjs, bje, bks, bke, NGHOST,
                            BoundaryFace::inner_x1, ph->w, pf->b, prad->ir);
  }

  // Apply boundary function on outer-x1
  if (pbval->apply_bndry_fn_[BoundaryFace::outer_x1]) {
    SetRadPhysicalFunctions(pmb, pco, time, dt,
                            pmb->is, pmb->ie, bjs, bje, bks, bke, NGHOST,
                            BoundaryFace::outer_x1, ph->w, pf->b, prad->ir);
  }

  if (pmb->block_size.nx2 > 1) { // 2D or 3D
    // Apply boundary function on inner-x2 and update W,bcc (if not periodic)
    if (pbval->apply_bndry_fn_[BoundaryFace::inner_x2]) {
      SetRadPhysicalFunctions(pmb, pco, time, dt,
                            bis, bie, pmb->js, pmb->je, bks, bke, NGHOST,
                            BoundaryFace::inner_x2, ph->w, pf->b, prad->ir);
    }

    if ((NR_RADIATION_ENABLED || IM_RADIATION_ENABLED) &&
           (pbval->block_bcs[BoundaryFace::inner_x2] != BoundaryFlag::block)) {
      if (prad->rotate_theta == 1) {
        pradbvar->RotateHPi_InnerX2(time, dt, bis, bie, pmb->js, bks, bke, NGHOST);
      }
    } // end radiation

    // Apply boundary function on outer-x2 and update W,bcc (if not periodic)
    if (pbval->apply_bndry_fn_[BoundaryFace::outer_x2]) {
      SetRadPhysicalFunctions(pmb, pco, time, dt,
                            bis, bie, pmb->js, pmb->je, bks, bke, NGHOST,
                            BoundaryFace::outer_x2, ph->w, pf->b, prad->ir);
    }
  }

  if (pmb->block_size.nx3 > 1) { // 3D
    bjs = pmb->js - NGHOST;
    bje = pmb->je + NGHOST;

    // Apply boundary function on inner-x3 and update W,bcc (if not periodic)
    if (pbval->apply_bndry_fn_[BoundaryFace::inner_x3]) {
      SetRadPhysicalFunctions(pmb, pco, time, dt,
                            bis, bie, bjs, bje, pmb->ks, pmb->ke, NGHOST,
                            BoundaryFace::inner_x3, ph->w, pf->b, prad->ir);
    }

    if ((NR_RADIATION_ENABLED || IM_RADIATION_ENABLED) &&
           (pbval->block_bcs[BoundaryFace::inner_x3] != BoundaryFlag::block)) {
      if (prad->rotate_phi == 1) {
        pradbvar->RotateHPi_InnerX3(time, dt, bis, bie, bjs, bje, pmb->ks, NGHOST);
      } else if (prad->rotate_phi == 2) {
        pradbvar->RotatePi_InnerX3(time, dt, bis, bie, bjs, bje, pmb->ks,NGHOST);
      }
    }

    // Apply boundary function on outer-x3 and update W,bcc (if not periodic)
    if (pbval->apply_bndry_fn_[BoundaryFace::outer_x3]) {
      SetRadPhysicalFunctions(pmb, pco, time, dt,
                            bis, bie, bjs, bje, pmb->ks, pmb->ke, NGHOST,
                            BoundaryFace::outer_x3, ph->w, pf->b, prad->ir);
    }
    if ((NR_RADIATION_ENABLED || IM_RADIATION_ENABLED) &&
           (pbval->block_bcs[BoundaryFace::outer_x3] != BoundaryFlag::block)) {
      if (prad->rotate_phi == 1) {
        pradbvar->RotateHPi_OuterX3(time, dt, bis, bie, bjs, bje, pmb->ke, NGHOST);
      } else if (prad->rotate_phi == 2) {
        pradbvar->RotatePi_OuterX3(time, dt, bis, bie, bjs, bje, pmb->ke, NGHOST);
      }
    }
  }
  return;
}




void RadBoundaryVariable::SetRadPhysicalFunctions(
    MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
    int il, int iu, int jl, int ju, int kl, int ku, int ngh,
    BoundaryFace face, AthenaArray<Real> &prim,
    FaceField &b, AthenaArray<Real> &ir) {

  NRRadiation *prad = pmb->pnrrad;
  RadBoundaryVariable *pradbvar = &(prad->rad_bvar);
  BoundaryValues *pbval = pmb->pbval;


  if (pbval->block_bcs[face] ==  BoundaryFlag::user) {  // user-enrolled BCs
      pmy_mesh_->RadBoundaryFunc_[face](pmb,pco,prad,prim,b, ir,time,dt,
                                             il,iu,jl,ju,kl,ku,NGHOST);
  }
  // KGF: this is only to silence the compiler -Wswitch warnings about not handling the
  // "undef" case when considering all possible BoundaryFace enumerator values. If "undef"
  // is actually passed to this function, it will likely die before that ATHENA_ERROR()
  // call, as it tries to access block_bcs[-1]
  std::stringstream msg;
  msg << "### FATAL ERROR in SetRadPhysicalFunctions" << std::endl
      << "face = BoundaryFace::undef passed to this function" << std::endl;


  switch (pbval->block_bcs[face]) {
    case BoundaryFlag::user: // handled above, outside loop over BoundaryVariable objs
      break;
    case BoundaryFlag::reflect:
      switch (face) {
        case BoundaryFace::undef:
          ATHENA_ERROR(msg);
        case BoundaryFace::inner_x1:
          pradbvar->ReflectInnerX1(time, dt, il, jl, ju, kl, ku, NGHOST);
          break;
        case BoundaryFace::outer_x1:
          pradbvar->ReflectOuterX1(time, dt, iu, jl, ju, kl, ku, NGHOST);
          break;
        case BoundaryFace::inner_x2:
          pradbvar->ReflectInnerX2(time, dt, il, iu, jl, kl, ku, NGHOST);
          break;
        case BoundaryFace::outer_x2:
          pradbvar->ReflectOuterX2(time, dt, il, iu, ju, kl, ku, NGHOST);
          break;
        case BoundaryFace::inner_x3:
          pradbvar->ReflectInnerX3(time, dt, il, iu, jl, ju, kl, NGHOST);
          break;
        case BoundaryFace::outer_x3:
          pradbvar->ReflectOuterX3(time, dt, il, iu, jl, ju, ku, NGHOST);
          break;
      }
      break;
    case BoundaryFlag::outflow:
      switch (face) {
        case BoundaryFace::undef:
          ATHENA_ERROR(msg);
        case BoundaryFace::inner_x1:
          pradbvar->OutflowInnerX1(time, dt, il, jl, ju, kl, ku, NGHOST);
          break;
        case BoundaryFace::outer_x1:
          pradbvar->OutflowOuterX1(time, dt, iu, jl, ju, kl, ku, NGHOST);
          break;
        case BoundaryFace::inner_x2:
          pradbvar->OutflowInnerX2(time, dt, il, iu, jl, kl, ku, NGHOST);
          break;
        case BoundaryFace::outer_x2:
          pradbvar->OutflowOuterX2(time, dt, il, iu, ju, kl, ku, NGHOST);
          break;
        case BoundaryFace::inner_x3:
          pradbvar->OutflowInnerX3(time, dt, il, iu, jl, ju, kl, NGHOST);
          break;
        case BoundaryFace::outer_x3:
          pradbvar->OutflowOuterX3(time, dt, il, iu, jl, ju, ku, NGHOST);
          break;
      }
      break;
    case BoundaryFlag::vacuum: // special boundary condition type for radiation
      switch (face) {
        case BoundaryFace::undef:
          ATHENA_ERROR(msg);
        case BoundaryFace::inner_x1:
          pradbvar->VacuumInnerX1(time, dt, il, jl, ju, kl, ku, NGHOST);
          break;
        case BoundaryFace::outer_x1:
          pradbvar->VacuumOuterX1(time, dt, iu, jl, ju, kl, ku, NGHOST);
          break;
        case BoundaryFace::inner_x2:
          pradbvar->VacuumInnerX2(time, dt, il, iu, jl, kl, ku, NGHOST);
          break;
        case BoundaryFace::outer_x2:
          pradbvar->VacuumOuterX2(time, dt, il, iu, ju, kl, ku, NGHOST);
          break;
        case BoundaryFace::inner_x3:
          pradbvar->VacuumInnerX3(time, dt, il, iu, jl, ju, kl, NGHOST);
          break;
        case BoundaryFace::outer_x3:
          pradbvar->VacuumOuterX3(time, dt, il, iu, jl, ju, ku, NGHOST);
          break;
      }
      break;
    case BoundaryFlag::polar_wedge:
      switch (face) {
        case BoundaryFace::undef:
          ATHENA_ERROR(msg);
        case BoundaryFace::inner_x2:
          pradbvar->PolarWedgeInnerX2(time, dt, il, iu, jl, kl, ku, NGHOST);
          break;
        case BoundaryFace::outer_x2:
          pradbvar->PolarWedgeOuterX2(time, dt, il, iu, ju, kl, ku, NGHOST);
          break;
        default:
          std::stringstream msg_polar;
          msg_polar << "### FATAL ERROR in SetRadPhysicalBoundary" << std::endl
                    << "Attempting to call polar wedge boundary function on \n"
                    << "MeshBlock boundary other than inner x2 or outer x2"
                    << std::endl;
          ATHENA_ERROR(msg_polar);
      }
      break;
    default:
      std::stringstream msg_flag;
      msg_flag << "### FATAL ERROR in SetRadPhysicalBoundary" << std::endl
               << "No BoundaryPhysics function associated with provided\n"
               << "block_bcs[" << face << "] = BoundaryFlag::"
               << GetBoundaryString(pbval->block_bcs[face]) << std::endl;
      ATHENA_ERROR(msg);
      break;
  } // end switch (block_bcs[face])
}
