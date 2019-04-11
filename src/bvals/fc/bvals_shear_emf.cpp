//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_shear_emf.cpp
//  \brief functions that apply BCs for face-centered flux corrections in shearing box
// calculations
//======================================================================================

// C headers

// C++ headers
#include <algorithm>  // min
#include <cmath>
#include <cstdlib>
#include <cstring>    // memcpy
#include <iomanip>
#include <iostream>   // endl
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
#include "../bvals_interfaces.hpp"
#include "../bvals.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif


//--------------------------------------------------------------------------------------
//! \fn int FaceCenteredBoundaryVariable::LoadEMFShearing(EdgeField &src,
//                                                        Real *buf, int nb)
//  \brief Load shearing box EMF boundary buffers

// KGF: EdgeField &src = shboxvar_outer_emf_, shboxvar_inner_emf_

void FaceCenteredBoundaryVariable::LoadEMFShearing(EdgeField &src,
                                                   Real *buf, const int nb) {
  MeshBlock *pmb = pmy_block_;
  int sj, sk, ej, ek;
  int psj, pej; // indices for e3
  int jo = pbval_->joverlap_;
  int nx2 = pmb->block_size.nx2 - NGHOST;
  sk = pmb->ks;        ek = pmb->ke;
  switch(nb) {
    case 0:
      sj = pmb->je - jo - (NGHOST - 1); ej = pmb->je;
      if (jo > nx2) sj = pmb->js;
      psj = sj; pej = ej;
      break;
    case 1:
      sj = pmb->js; ej = pmb->je - jo + NGHOST;
      if (jo < NGHOST) ej = pmb->je;
      psj = sj; pej = ej + 1;
      break;
    case 2:
      sj = pmb->je - (NGHOST - 1); ej = pmb->je;
      if (jo > nx2) sj = pmb->je - (jo - nx2) + 1;
      psj = sj; pej = ej;
      break;
    case 3:
      sj = pmb->js; ej = pmb->js + (NGHOST - 1);
      if (jo < NGHOST) ej = pmb->js + (NGHOST - jo) - 1;
      psj = sj + 1; pej = ej + 1;
      break;
    case 4:
      sj = pmb->js; ej = pmb->js + jo + NGHOST - 1;
      if (jo > nx2) ej = pmb->je;
      psj = sj; pej = ej + 1;
      break;
    case 5:
      sj = pmb->js + jo - NGHOST; ej = pmb->je;
      if (jo < NGHOST) sj = pmb->js;
      psj = sj; pej = ej + 1;
      break;
    case 6:
      sj = pmb->js; ej = pmb->js + (NGHOST - 1);
      if (jo > nx2) ej = pmb->js + (jo - nx2) - 1;
      psj = sj + 1; pej = ej + 1;
      break;
    case 7:
      sj = pmb->je - (NGHOST - 1); ej = pmb->je;
      if (jo < NGHOST) sj = pmb->je - (NGHOST - jo) + 1;
      psj = sj; pej = ej;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in FaceCenteredBoundaryVariable:LoadEMFShearing\n"
          << "nb = " << nb << " not valid" << std::endl;
      ATHENA_ERROR(msg);
  }
  int p = 0;
  // pack e2
  for (int k=sk; k<=ek+1; k++) {
    for (int j=sj; j<=ej; j++)
      buf[p++] = src.x2e(k,j);
  }
  // pack e3
  for (int k=sk; k<=ek; k++) {
    for (int j=psj; j<=pej; j++)
      buf[p++] = src.x3e(k,j);
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::SendEMFShearingBoxBoundaryCorrection()
//  \brief Send shearing box boundary buffers for EMF correction

void FaceCenteredBoundaryVariable::SendEMFShearingBoxBoundaryCorrection() {
  MeshBlock *pmb = pmy_block_;
  int js = pmb->js; int ks = pmb->ks;
  int je = pmb->je; int ke = pmb->ke;
  int nx2 = pmb->block_size.nx2;
  int nx3 = pmb->block_size.nx3;

  if (is_shear[0]) {
    // step 1. -- average edges of shboxvar_fc_flx_
    // average e3 for x1x2 edge
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je+1; j+=nx2)
        shboxvar_inner_fc_flx_.x3e(k,j) *= 0.5;
    }
    // average e2 for x1x3 edge
    for (int k=ks; k<=ke+1; k+=nx3) {
      for (int j=js; j<=je; j++)
        shboxvar_inner_fc_flx_.x2e(k,j) *= 0.5;
    }

    // step 2. -- load sendbuf; memcpy to recvbuf if on same rank, post
    // MPI_Isend otherwise
    for (int n=0; n<4; n++) {
      if (send_inner_rank_[n] != -1) {
        LoadEMFShearing(shboxvar_inner_fc_flx_, send_innerbuf_fc_flx_[n], n);
        if (send_inner_rank_[n] == Globals::my_rank) {// on the same process
          MeshBlock *pbl = pmb->pmy_mesh->FindMeshBlock(send_inner_gid_[n]);
          std::memcpy(pbl->pbval->recv_innerbuf_fc_flx_[n],
                      send_innerbuf_fc_flx_[n], send_innersize_fc_flx_[n]*sizeof(Real));
          pbl->pbval->shbox_inner_fc_flx_flag_[n] = BoundaryStatus::arrived;
        } else { // MPI
#ifdef MPI_PARALLEL
          int tag = CreateBvalsMPITag(send_inner_lid_[n], n,sh_fc_flx_phys_id_);
          MPI_Isend(send_innerbuf_fc_flx_[n], send_innersize_fc_flx_[n],
                    MPI_ATHENA_REAL, send_inner_rank_[n], tag,
                    MPI_COMM_WORLD, &rq_innersend_fc_flx_[n]);
#endif
        }
      }
    }
  } // inner boundaries

  if (is_shear[1]) {
    // step 1. -- average edges of shboxvar_fc_flx_
    // average e3 for x1x2 edge
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je+1; j+=nx2)
        shboxvar_outer_fc_flx_.x3e(k,j) *= 0.5;
    }
    // average e2 for x1x3 edge
    for (int k=ks; k<=ke+1; k+=nx3) {
      for (int j=js; j<=je; j++)
        shboxvar_outer_fc_flx_.x2e(k,j) *= 0.5;
    }

    // step 2. -- load sendbuf; memcpy to recvbuf if on same rank, post
    // MPI_Isend otherwise
    int offset = 4;
    for (int n=0; n<4; n++) {
      if (send_outer_rank_[n] != -1) {
        LoadEMFShearing(shboxvar_outer_fc_flx_, send_outerbuf_fc_flx_[n], n+offset);
        if (send_outer_rank_[n] == Globals::my_rank) {// on the same process
          MeshBlock *pbl = pmb->pmy_mesh->FindMeshBlock(send_outer_gid_[n]);
          std::memcpy(pbl->pbval->recv_outerbuf_fc_flx_[n],
                      send_outerbuf_fc_flx_[n], send_outersize_fc_flx_[n]*sizeof(Real));
          pbl->pbval->shbox_outer_fc_flx_flag_[n] = BoundaryStatus::arrived;
        } else { // MPI
#ifdef MPI_PARALLEL
          int tag = CreateBvalsMPITag(send_outer_lid_[n], n+offset,
                                     sh_fc_flx_phys_id_);
          MPI_Isend(send_outerbuf_fc_flx_[n],send_outersize_fc_flx_[n],
                    MPI_ATHENA_REAL, send_outer_rank_[n], tag,
                    MPI_COMM_WORLD, &rq_outersend_fc_flx_[n]);
#endif
        }
      }
    }
  } // outer boundaries
  return;
}

// --------------------------------------------------------------------------------------
// ! \fn void FaceCenteredBoundaryVariable::SetEMFShearingBoxBoundarySameLevel(
//                                   EdgeField &dst, Real *buf, const int nb)
//  \brief Set EMF shearing box boundary received from a block on the same level

// KGF: EdgeField &dst = shboxmap_outer_emf_, shboxmap_inner_emf_
void FaceCenteredBoundaryVariable::SetEMFShearingBoxBoundarySameLevel(EdgeField &dst,
                                                                      Real *buf,
                                                                      const int nb) {
  MeshBlock *pmb = pmy_block_;
  int sj, sk, ej, ek;
  int psj, pej;
  int jo = pbval_->joverlap_;
  int nx2 = pmb->block_size.nx2 - NGHOST;
  int nxo = pmb->block_size.nx2 - jo;

  sk = pmb->ks; ek = pmb->ke;
  switch(nb) {
    case 0:
      sj = pmb->js - NGHOST; ej = pmb->js + (jo - 1);
      if (jo > nx2) sj = pmb->js - nxo;
      psj = sj; pej = ej + 1;
      break;
    case 1:
      sj = pmb->js + jo; ej = pmb->je + NGHOST;
      if (jo < NGHOST) ej = pmb->je + jo;
      psj = sj; pej = ej + 1;
      break;
    case 2:
      sj = pmb->js - NGHOST; ej = pmb->js - 1;
      if (jo > nx2) ej = pmb->js - nxo - 1;
      psj = sj; pej = ej;
      break;
    case 3:
      sj = pmb->je + jo + 1; ej = pmb->je + NGHOST;
      psj = sj + 1; pej = ej + 1;
      break;
    case 4:
      sj = pmb->je - (jo - 1); ej = pmb->je + NGHOST;
      if (jo > nx2) ej = pmb->je + nxo;
      psj = sj; pej = ej + 1;
      break;
    case 5:
      sj = pmb->js - NGHOST; ej = pmb->je - jo;
      if (jo <= NGHOST) sj = pmb->js - jo;
      psj = sj; pej = ej + 1;
      break;
    case 6:
      sj = pmb->je + 1; ej = pmb->je + NGHOST;
      if (jo > nx2) sj = pmb->je + nxo + 1;
      psj = sj + 1; pej = ej + 1;
      break;
    case 7:
      sj = pmb->js - NGHOST; ej = pmb->js - jo - 1;
      psj = sj; pej = ej;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in FaceCenteredBoundaryVariable:SetFieldShearing\n"
          << "nb = " << nb << " not valid" << std::endl;
      ATHENA_ERROR(msg);
  }
  int p = 0;
  // unpack e2
  for (int k=sk; k<=ek+1; k++) {
    for (int j=sj; j<=ej; j++)
      dst.x2e(k,j) += buf[p++];
  }
  // unpack e3
  for (int k=sk; k<=ek; k++) {
    for (int j=psj; j<=pej; j++)
      dst.x3e(k,j) += buf[p++];
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn bool FaceCenteredBoundaryVariable::ReceiveEMFShearingBoxBoundaryCorrection()
//  \brief receive shearing box boundary data for EMF correction

bool FaceCenteredBoundaryVariable::ReceiveEMFShearingBoxBoundaryCorrection() {
  bool flagi = true, flago = true;

  if (is_shear[0]) { // check inner boundaries
    for (int n=0; n<4; n++) {
      if (shbox_inner_fc_flx_flag_[n] == BoundaryStatus::completed) continue;
      if (shbox_inner_fc_flx_flag_[n] == BoundaryStatus::waiting) {
        if (recv_inner_rank_[n] == Globals::my_rank) {// on the same process
          flagi = false;
          continue;
        } else { // MPI boundary
#ifdef MPI_PARALLEL
          int test;
          MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &test,
                     MPI_STATUS_IGNORE);
          MPI_Test(&rq_innerrecv_fc_flx_[n], &test, MPI_STATUS_IGNORE);
          if (static_cast<bool>(test) == false) {
            flagi = false;
            continue;
          }
          shbox_inner_fc_flx_flag_[n] = BoundaryStatus::arrived;
#endif
        }
      }
      // set dst if boundary arrived
      SetEMFShearingBoxBoundarySameLevel(shboxmap_inner_fc_flx_, recv_innerbuf_fc_flx_[n], n);
      shbox_inner_fc_flx_flag_[n] = BoundaryStatus::completed; // completed
    }
  } // inner boundary

  if (is_shear[1]) { // check outer boundaries
    int offset = 4;
    for (int n=0; n<4; n++) {
      if (shbox_outer_fc_flx_flag_[n] == BoundaryStatus::completed) continue;
      if (shbox_outer_fc_flx_flag_[n] == BoundaryStatus::waiting) {
        if (recv_outer_rank_[n] == Globals::my_rank) {  // on the same process
          flago = false;
          continue;
        } else { // MPI boundary
#ifdef MPI_PARALLEL
          int test;
          MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &test,
                     MPI_STATUS_IGNORE);
          MPI_Test(&rq_outerrecv_fc_flx_[n], &test, MPI_STATUS_IGNORE);
          if (static_cast<bool>(test) == false) {
            flago = false;
            continue;
          }
          shbox_outer_fc_flx_flag_[n] = BoundaryStatus::arrived;
#endif
        }
      }
      SetEMFShearingBoxBoundarySameLevel(shboxmap_outer_fc_flx_, recv_outerbuf_fc_flx_[n],
                                         n+offset);
      shbox_outer_fc_flx_flag_[n] = BoundaryStatus::completed; // completed
    }
  } // outer boundary
  return (flagi && flago);
}


//--------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::RemapEMFShearingBoxBoundary()
//  \brief Set EMF boundary received from a block on the finer level

void FaceCenteredBoundaryVariable::RemapEMFShearingBoxBoundary() {
  MeshBlock *pmb = pmy_block_;
  AthenaArray<Real> &e2 = pmb->pfield->e.x2e;
  // AthenaArray<Real> &e3 = pmb->pfield->e.x3e;
  int ks = pmb->ks, ke = pmb->ke;
  int js = pmb->js, je = pmb->je;
  int is = pmb->is, ie = pmb->ie;
  Real eps = pbval_->eps_;

  if (is_shear[0]) {
    ClearEMFShearing(shboxvar_inner_fc_flx_);
    // step 1.-- conservative remapping
    for (int k=ks; k<=ke+1; k++) {
      RemapFluxEMF(k, js, je+2, eps, shboxmap_inner_fc_flx_.x2e, flx_inner_fc_flx_.x2e);
      for (int j=js; j<=je; j++) {
        shboxmap_inner_fc_flx_.x2e(k,j) -= flx_inner_fc_flx_.x2e(j+1) - flx_inner_fc_flx_.x2e(j);
      }
    }
    // step 2.-- average the EMF correction
    // average e2
    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je; j++)
        e2(k,j,is) = 0.5*(e2(k,j,is) + shboxmap_inner_fc_flx_.x2e(k,j));
    }
    ClearEMFShearing(shboxmap_inner_fc_flx_);
  }

  if (is_shear[1]) {
    ClearEMFShearing(shboxvar_outer_fc_flx_);
    // step 1.-- conservative remapping
    for (int k=ks; k<=ke+1; k++) { // e2
      RemapFluxEMF(k, js-1, je+1, -eps, shboxmap_outer_fc_flx_.x2e, flx_outer_fc_flx_.x2e);
      for (int j=js; j<=je; j++)
        shboxmap_outer_fc_flx_.x2e(k,j) -= flx_outer_fc_flx_.x2e(j+1) - flx_outer_fc_flx_.x2e(j);
    }
    // step 2.-- average the EMF correction
    // average e2
    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je; j++)
        e2(k,j,ie+1) = 0.5*(e2(k,j,ie+1) + shboxmap_outer_fc_flx_.x2e(k,j));
    }
    ClearEMFShearing(shboxmap_outer_fc_flx_);
  }
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::ClearEMFShearing()
//  \brief Clear the working array for EMFs on the surface/edge contacting with
//  a shearing periodic boundary

void FaceCenteredBoundaryVariable::ClearEMFShearing(EdgeField &work) {
  MeshBlock *pmb = pmy_block_;
  AthenaArray<Real> &e2 = work.x2e;
  AthenaArray<Real> &e3 = work.x3e;
  int ks = pmb->ks, ke = pmb->ke;
  int js = pmb->js, je = pmb->je;
  for (int k=ks-NGHOST; k<=ke+NGHOST; k++) {
    for (int j=js-NGHOST; j<=je+NGHOST; j++) {
      e2(k,j) = 0.0;
      e3(k,j) = 0.0;
      if (k == ke + NGHOST) e2(k+1,j) = 0.0;
      if (j == je + NGHOST) e3(k,j+1) = 0.0;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::RemapFluxEMF(int k, int jinner, int jouter, Real
//  eps, static AthenaArray<Real> &U, AthenaArray<Real> &Flux)
//  \brief compute the flux along j indices for remapping
//  adopted from 2nd order RemapFlux of Athena4.0

void FaceCenteredBoundaryVariable::RemapFluxEMF(const int k, const int jinner,
                                                const int jouter,
                                                const Real eps,
                                                const AthenaArray<Real> &var,
                                                AthenaArray<Real> &flux) {
  int j, jl, ju;
  Real dUc, dUl, dUr, dUm, lim_slope;

  // jinner, jouter are index range over which flux must be returned.  Set loop
  // limits depending on direction of upwind differences
  if (eps > 0.0) { // eps always > 0 for inner i boundary
    jl = jinner - 1;
    ju = jouter - 1;
  } else {         // eps always < 0 for outer i boundary
    jl = jinner;
    ju = jouter;
  }

  // TODO(felker): do not reimplement PLM here; use plm.cpp.
  // TODO(felker): relax assumption that 2nd order reconstruction must be used
  for (j=jl; j<=ju; j++) {
    dUc = var(k,j+1) - var(k,j-1);
    dUl = var(k,j  ) - var(k,j-1);
    dUr = var(k,j+1) - var(k,j  );

    dUm = 0.0;
    if (dUl*dUr > 0.0) {
      lim_slope = std::min(std::fabs(dUl),std::fabs(dUr));
      dUm = SIGN(dUc)*std::min(0.5*std::fabs(dUc),2.0*lim_slope);
    }

    if (eps > 0.0) { // eps always > 0 for inner i boundary
      flux(j+1) = eps*(var(k,j) + 0.5*(1.0 - eps)*dUm);
    } else {         // eps always < 0 for outer i boundary
      flux(j  ) = eps*(var(k,j) - 0.5*(1.0 + eps)*dUm);
    }
  }
  return;
}
