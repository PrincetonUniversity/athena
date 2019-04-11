//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_shear_fc.cpp
//  \brief functions that apply shearing box BCs for face-centered quantities
//======================================================================================

// C headers

// C++ headers
#include <algorithm>  // min
#include <cmath>
#include <cstdlib>
#include <cstring>    // std::memcpy
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
//! \fn int FaceCenteredBoundaryVariable::LoadShearing(FaceField &src, Real *buf, int nb)
//  \brief Load shearing box field boundary buffers

// KGF: FaceField &src = shboxvar_inner/outer_field_, send_innerbuf/outerbuf_field_[n]
void FaceCenteredBoundaryVariable::LoadShearing(FaceField &src, Real *buf, int nb) {
  MeshBlock *pmb = pmy_block_;
  Mesh *pmesh = pmb->pmy_mesh;
  int si, sj, sk, ei, ej, ek;
  int psj, pej; // indices for bx2
  int jo = pbval_->joverlap_;
  int nx2 = pmb->block_size.nx2 - NGHOST;

  si = pmb->is - NGHOST; ei = pmb->is - 1;
  sk = pmb->ks;        ek = pmb->ke;
  if (pmesh->mesh_size.nx3 > 1)  ek += NGHOST, sk -= NGHOST;
  switch(nb) {
    case 0:
      sj = pmb->je - jo - (NGHOST - 1); ej = pmb->je;
      if (jo > nx2) sj = pmb->js;
      psj = sj; pej = ej + 1;
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
      msg << "### FATAL ERROR in FaceCenteredBoundaryVariable:LoadShearing " << std::endl
          << "nb = " << nb << " not valid" << std::endl;
      ATHENA_ERROR(msg);
  }
  int p = 0;
  BufferUtility::PackData(src.x1f, buf, si, ei, sj, ej, sk, ek, p);
  BufferUtility::PackData(src.x2f, buf, si, ei, psj, pej, sk, ek, p);
  BufferUtility::PackData(src.x3f, buf, si, ei, sj, ej, sk, ek+1, p);
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::SendShearingBoxBoundaryBuffers()
//  \brief Send shearing box boundary buffers for field variables

void FaceCenteredBoundaryVariable::SendShearingBoxBoundaryBuffers() {
  MeshBlock *pmb = pmy_block_;
  Mesh *pmesh = pmb->pmy_mesh;

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int ku, kl;
  if (pmesh->mesh_size.nx3 > 1) {
    ku = ke + NGHOST;
    kl = ks - NGHOST;
  } else {
    ku = ke;
    kl = ks;
  }
  Real eps = pbval_->eps_;
  Real qomL = pbval_->qomL_;

  if (pbval_->is_shear[0]) {
    int ib = is - NGHOST;
    int ii;
    // step 1. -- load shboxvar_fc_
    for (int k=kl; k<=ku; k++) {
      for (int j=js-NGHOST; j<=je+NGHOST; j++) {
        for (int i=0; i<NGHOST; i++) {
          ii = ib+i;
          shear_fc_[0].x1f(k,j,i) = (*var_fc).x1f(k,j,ii);
          shear_fc_[0].x2f(k,j,i) = (*var_fc).x2f(k,j,ii);
          shear_fc_[0].x3f(k,j,i) = (*var_fc).x3f(k,j,ii);
        }
      }
    }
    // fill the extra cells for B2i and B3i
    int kp = ku + 1;
    for (int j=js-NGHOST; j<=je+NGHOST; j++) {
      for (int i=0; i<NGHOST; i++) {
        ii = ib+i;
        shear_fc_[0].x3f(kp,j,i) = (*var_fc).x3f(kp,j,ii);
      }
    }
    int jp = je + NGHOST + 1;
    for (int k=kl; k<=ku; k++) {
      for (int i=0; i<NGHOST; i++) {
        ii = ib+i;
        shear_fc_[0].x2f(k,jp,i) = (*var_fc).x2f(k,jp,ii);
      }
    }

    // step 2. -- conservative remapping
    for (int k=kl; k<=ku; k++) {  // bx1
      for (int i=0; i<NGHOST; i++) {
        RemapFlux(k, js, je+2, i, eps, shear_fc_[0].x1f, shear_flx_fc_[0].x1f);
        for (int j=js; j<=je+1; j++) {
          shear_fc_[0].x1f(k,j,i) -= (shear_flx_fc_[0].x1f(j+1) -
                                               shear_flx_fc_[0].x1f(j));
        }
      }
    }
    for (int k=kl; k<=ku; k++) {  // bx2
      for (int i=0; i<NGHOST; i++) {
        RemapFlux(k, js, je+3, i, eps, shear_fc_[0].x2f, shear_flx_fc_[0].x2f);
        for (int j=js; j<=je+2; j++) {
          shear_fc_[0].x2f(k,j,i) -= (shear_flx_fc_[0].x2f(j+1) -
                                               shear_flx_fc_[0].x2f(j));
        }
      }
    }
    for (int k=kl; k<=ku+1; k++) { // bx3
      for (int i=0; i<NGHOST; i++) {
        RemapFlux(k, js, je+2, i, eps, shear_fc_[0].x3f, shear_flx_fc_[0].x3f);
        for (int j=js; j<=je+1; j++) {
          shear_fc_[0].x3f(k,j,i) -= (shear_flx_fc_[0].x3f(j+1) -
                                               shear_flx_fc_[0].x3f(j));
        }
      }
    }

    // step 3. -- load sendbuf; memcpy to recvbuf if on same rank, else post MPI_Isend
    for (int n=0; n<4; n++) {
      if (pbval_->shear_send_neighbor_[0][n].rank != -1) {
        LoadShearing(shear_fc_[0], shear_bd_fc_[0].send[n], n);
        if (pbval_->shear_send_neighbor_[0][n].rank == Globals::my_rank) {// on the same process
          MeshBlock *pbl=pmb->pmy_mesh->FindMeshBlock(pbval_->shear_send_neighbor_[0][n].gid);
          std::memcpy(pbl->pbval->shear_bd_fc_[0].recv[n],shear_bd_fc_[0].send[n],
                      shear_send_count_fc_[0][n]*sizeof(Real));
          pbl->pbval->shear_bd_fc_[0].flag[n] = BoundaryStatus::arrived;
        } else { // MPI
#ifdef MPI_PARALLEL
          // bufid = n
          int tag = pbval_->CreateBvalsMPITag(pbval_->shear_send_neighbor_[0][n].lid, n, shear_fc_phys_id_);
          MPI_Isend(shear_bd_fc_[0].send[n],shear_send_count_fc_[0][n],MPI_ATHENA_REAL,
                    pbval_->shear_send_neighbor_[0][n].rank,tag,MPI_COMM_WORLD, &shear_bd_fc_[0].req_send[n]);
#endif
        }
      }
    }
  } // inner boundaries

  if (pbval_->is_shear[1]) {
    int  ib = ie + 1;
    qomL = -qomL;
    int ii;
    // step 1. -- load shboxvar_fc_
    for (int k=kl; k<=ku; k++) {
      for (int j=js-NGHOST; j<=je+NGHOST; j++) {
        for (int i=0; i<NGHOST; i++) {
          ii = ib+i;
          shear_fc_[1].x1f(k,j,i) = (*var_fc).x1f(k,j,ii+1);
          shear_fc_[1].x2f(k,j,i) = (*var_fc).x2f(k,j,ii);
          shear_fc_[1].x3f(k,j,i) = (*var_fc).x3f(k,j,ii);
        }
      }
    }
    // fill the extra cells for B2i and B3i
    int kp = ku + 1;
    for (int j=js-NGHOST; j<=je+NGHOST; j++) {
      for (int i=0; i<NGHOST; i++) {
        ii = ib+i;
        shear_fc_[1].x3f(kp,j,i) = (*var_fc).x3f(kp,j,ii);
      }
    }
    int jp = je + NGHOST + 1;
    for (int k=kl; k<=ku; k++) {
      for (int i=0; i<NGHOST; i++) {
        ii = ib+i;
        shear_fc_[1].x2f(k,jp,i) = (*var_fc).x2f(k,jp,ii);
      }
    }

    // step 2. -- conservative remapping
    for (int k=kl; k<=ku; k++) {  // bx1
      for (int i=0; i<NGHOST; i++) {
        RemapFlux(k, js-1, je+1, i, -eps, shear_fc_[1].x1f,
                       shear_flx_fc_[1].x1f);
        for (int j=js-1; j<=je; j++) {
          shear_fc_[1].x1f(k,j,i) -= (shear_flx_fc_[1].x1f(j+1) -
                                               shear_flx_fc_[1].x1f(j));
        }
      }
    }
    for (int k=kl; k<=ku; k++) {  // bx2
      for (int i=0; i<NGHOST; i++) {
        RemapFlux(k, js-1, je+2, i, -eps, shear_fc_[1].x2f,
                       shear_flx_fc_[1].x2f);
        for (int j=js-1; j<=je+1; j++) {
          shear_fc_[1].x2f(k,j,i) -= (shear_flx_fc_[1].x2f(j+1) -
                                               shear_flx_fc_[1].x2f(j));
        }
      }
    }
    for (int k=kl; k<=ku+1; k++) {  // bx3
      for (int i=0; i<NGHOST; i++) {
        RemapFlux(k, js-1, je+1, i, -eps, shear_fc_[1].x3f,
                  shear_flx_fc_[1].x3f);
        for (int j=js-1; j<=je; j++) {
          shear_fc_[1].x3f(k,j,i) -= (shear_flx_fc_[1].x3f(j+1)
                                            - shear_flx_fc_[1].x3f(j));
        }
      }
    }

    // step 3. -- load sendbuf; memcpy to recvbuf if on same rank, else post MPI_Isend
    int offset = 4;
    for (int n=0; n<4; n++) {
      if (pbval_->shear_send_neighbor_[1][n].rank != -1) {
        LoadShearing(shear_fc_[1], shear_bd_fc_[1].send[n], n+offset);
        if (pbval_->shear_send_neighbor_[1][n].rank == Globals::my_rank) {// on the same process
          MeshBlock *pbl=pmb->pmy_mesh->FindMeshBlock(pbval_->shear_send_neighbor_[1][n].gid);
          std::memcpy(pbl->pbval->shear_bd_fc_[1].recv[n],
                      shear_bd_fc_[1].send[n], shear_send_count_fc_[1][n]*sizeof(Real));
          pbl->pbval->shear_bd_fc_[1].flag[n] = BoundaryStatus::arrived;
        } else { // MPI
#ifdef MPI_PARALLEL
          // bufid for outer(inner): 2(0) and 3(1)
          int tag = pbval_->CreateBvalsMPITag(pbval_->shear_send_neighbor_[1][n].lid, n+offset, shear_fc_phys_id_);
          MPI_Isend(shear_bd_fc_[1].send[n], shear_send_count_fc_[1][n], MPI_ATHENA_REAL,
                    pbval_->shear_send_neighbor_[1][n].rank, tag, MPI_COMM_WORLD, &shear_bd_fc_[1].req_send[n]);
#endif
        }
      }
    }
  } // outer boundaries
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::SetShearingBoxBoundarySameLevel(
//                                           Real *buf, const int nb)
//  \brief Set field shearing box boundary received from a block on the same level

// KGF: FaceField &dst = passed through pmb->pfield->b from
// ReceiveFieldShearingboxBoundaryBuffers(FaceField &dst)

void FaceCenteredBoundaryVariable::SetShearingBoxBoundarySameLevel(Real *buf,
                                                                   const int nb) {
  MeshBlock *pmb = pmy_block_;
  Mesh *pmesh = pmb->pmy_mesh;
  int si, sj, sk, ei, ej, ek;
  int psi, pei, psj, pej;
  int jo = pbval_->joverlap_;
  int nx2 = pmb->block_size.nx2 - NGHOST;
  int nxo = pmb->block_size.nx2 - jo;

  sk = pmb->ks; ek = pmb->ke;
  if (pmesh->mesh_size.nx3 > 1) ek += NGHOST, sk -= NGHOST;
  switch(nb) {
    case 0:
      si = pmb->is - NGHOST; ei = pmb->is - 1;
      sj = pmb->js - NGHOST; ej = pmb->js + (jo - 1);
      if (jo > nx2) sj = pmb->js - nxo;
      psi = si; pei = ei; psj = sj; pej = ej + 1;
      break;
    case 1:
      si = pmb->is - NGHOST; ei = pmb->is - 1;
      sj = pmb->js + jo; ej = pmb->je + NGHOST;
      if (jo < NGHOST) ej = pmb->je + jo;
      psi = si; pei = ei; psj = sj; pej = ej + 1;
      break;
    case 2:
      si = pmb->is - NGHOST; ei = pmb->is - 1;
      sj = pmb->js - NGHOST; ej = pmb->js - 1;
      if (jo > nx2) ej = pmb->js - nxo - 1;
      psi = si; pei = ei; psj = sj; pej = ej;
      break;
    case 3:
      si = pmb->is - NGHOST; ei = pmb->is - 1;
      sj = pmb->je + jo + 1; ej = pmb->je + NGHOST;
      psi = si; pei = ei; psj = sj + 1; pej = ej + 1;
      break;
    case 4:
      si = pmb->ie + 1; ei = pmb->ie + NGHOST;
      sj = pmb->je - (jo - 1); ej = pmb->je + NGHOST;
      if (jo > nx2) ej = pmb->je + nxo;
      psi = si + 1; pei = ei + 1; psj = sj; pej = ej + 1;
      break;
    case 5:
      si = pmb->ie + 1; ei = pmb->ie + NGHOST;
      sj = pmb->js - NGHOST; ej = pmb->je - jo;
      if (jo < NGHOST) sj = pmb->js - jo;
      psi = si + 1; pei = ei + 1; psj = sj; pej = ej + 1;
      break;
    case 6:
      si = pmb->ie + 1; ei = pmb->ie + NGHOST;
      sj = pmb->je + 1; ej = pmb->je + NGHOST;
      if (jo > nx2) sj = pmb->je + nxo + 1;
      psi = si + 1; pei = ei + 1; psj = sj + 1; pej = ej + 1;
      break;
    case 7:
      si = pmb->ie + 1; ei = pmb->ie + NGHOST;
      sj = pmb->js - NGHOST; ej = pmb->js - jo - 1;
      psi = si + 1; pei = ei + 1; psj = sj; pej = ej;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in FaceCenteredBoundaryVariable:SetShearing " << std::endl
          << "nb = " << nb << " not valid" << std::endl;
      ATHENA_ERROR(msg);
  }

  // set [sj:ej] of current meshblock
  int p = 0;
  BufferUtility::UnpackData(buf, (*var_fc).x1f, psi, pei, sj, ej, sk, ek, p);
  BufferUtility::UnpackData(buf, (*var_fc).x2f, si, ei, psj, pej, sk, ek, p);
  BufferUtility::UnpackData(buf, (*var_fc).x3f, si, ei, sj, ej, sk, ek+1, p);
  return;
}

//--------------------------------------------------------------------------------------
//! \fn bool FaceCenteredBoundaryVariable::ReceiveShearingBoxBoundaryBuffers()
//  \brief receive shearing box boundary data for field(face-centered) variables

bool FaceCenteredBoundaryVariable::ReceiveShearingBoxBoundaryBuffers() {
  bool flagi = true, flago = true;

  if (pbval_->is_shear[0]) { // check inner boundaries
    for (int n=0; n<4; n++) {
      if (shear_bd_fc_[0].flag[n] == BoundaryStatus::completed) continue;
      if (shear_bd_fc_[0].flag[n] == BoundaryStatus::waiting) {
        if (pbval_->shear_recv_neighbor_[0][n].rank == Globals::my_rank) {// on the same process
          flagi = false;
          continue;
        } else { // MPI boundary
#ifdef MPI_PARALLEL
          int test;
          MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &test,
                     MPI_STATUS_IGNORE);
          MPI_Test(&shear_bd_fc_[0].req_recv[n], &test, MPI_STATUS_IGNORE);
          if (static_cast<bool>(test) == false) {
            flagi = false;
            continue;
          }
          shear_bd_fc_[0].flag[n] = BoundaryStatus::arrived;
#endif
        }
      }
      // set dst if boundary arrived
      SetShearingBoxBoundarySameLevel(shear_bd_fc_[0].recv[n], n);
      shear_bd_fc_[0].flag[n] = BoundaryStatus::completed; // completed
    } // loop over recv[0] to recv[3]
  } // inner boundary

  if (pbval_->is_shear[1]) { // check outer boundaries
    int offset = 4;
    for (int n=0; n<4; n++) {
      if (shear_bd_fc_[1].flag[n] == BoundaryStatus::completed) continue;
      if (shear_bd_fc_[1].flag[n] == BoundaryStatus::waiting) {
        if (pbval_->shear_recv_neighbor_[1][n].rank == Globals::my_rank) {// on the same process
          flago = false;
          continue;
        } else { // MPI boundary
#ifdef MPI_PARALLEL
          int test;
          MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &test,
                     MPI_STATUS_IGNORE);

          MPI_Test(&shear_bd_fc_[1].req_recv[n], &test, MPI_STATUS_IGNORE);
          if (static_cast<bool>(test) == false) {
            flago = false;
            continue;
          }
          shear_bd_fc_[1].flag[n] = BoundaryStatus::arrived;
#endif
        }
      }
      SetShearingBoxBoundarySameLevel(shear_bd_fc_[1].recv[n], n+offset);
      shear_bd_fc_[1].flag[n] = BoundaryStatus::completed; // completed
    } // loop over recv[0] and recv[1]
  } // outer boundary
  return (flagi && flago);
}

//--------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::RemapFlux(const int k, const int jinner,
///                                         const int jouter, int i,
//                                          Real eps, static AthenaArray<Real> &U,
//                                          AthenaArray<Real> &Flux)
//  \brief compute the flux along j indices for remapping adopted from 2nd order RemapFlux
//         of Athena4.0

void FaceCenteredBoundaryVariable::RemapFlux(const int k, const int jinner,
                                             const int jouter,
                                             const int i, const Real eps,
                                             const AthenaArray<Real> &var,
                                             AthenaArray<Real> &flux) {
  int j, jl, ju;
  Real dUc, dUl, dUr, dUm, lim_slope;

  // jinner, jouter are index range over which flux must be returned.  Set loop
  // limits depending on direction of upwind differences

  if (eps > 0.0) { //eps always > 0 for inner i boundary
    jl = jinner - 1;
    ju = jouter - 1;
  } else {         // eps always < 0 for outer i boundary
    jl = jinner;
    ju = jouter;
  }

  for (j=jl; j<=ju; j++) {
    dUc = var(k,j+1,i) - var(k,j-1,i);
    dUl = var(k,j,  i) - var(k,j-1,i);
    dUr = var(k,j+1,i) - var(k,j,  i);

    dUm = 0.0;
    if (dUl*dUr > 0.0) {
      lim_slope = std::min(std::fabs(dUl),std::fabs(dUr));
      dUm = SIGN(dUc)*std::min(0.5*std::fabs(dUc),2.0*lim_slope);
    }

    if (eps > 0.0) { // eps always > 0 for inner i boundary
      flux(j+1) = eps*(var(k,j,i) + 0.5*(1.0 - eps)*dUm);
    } else {         // eps always < 0 for outer i boundary
      flux(j  ) = eps*(var(k,j,i) - 0.5*(1.0 + eps)*dUm);
    }
  }
  return;
}
