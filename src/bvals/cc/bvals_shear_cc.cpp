//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_shear_cc.cpp
//  \brief functions that apply shearing box BCs for cell-centered variables
//========================================================================================

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
#include "../bvals.hpp"
#include "../bvals_interfaces.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif


//--------------------------------------------------------------------------------------
//! \fn int CellCenteredBoundaryVariable::LoadShearing(AthenaArray<Real> &src, Real *buf,
//                                                     int nb)
//  \brief Load shearing box hydro boundary buffers

// KGF: AthenaArray<Real> &src = shboxvar_inner_hydro_, shboxvar_outer_hydro_
// TODO(KGF): is this identical to the copy in FaceCenteredBoundaryVariable??
void CellCenteredBoundaryVariable::LoadShearing(AthenaArray<Real> &src, Real *buf,
                                                int nb) {
  MeshBlock *pmb = pmy_block_;
  Mesh *pmesh = pmb->pmy_mesh;
  int si, sj, sk, ei, ej, ek;
  int nx2 = pmb->block_size.nx2 - NGHOST;
  int jo = pbval_->joverlap_;

  si = pmb->is - NGHOST; ei = pmb->is - 1;
  sk = pmb->ks;        ek = pmb->ke;
  if (pmesh->mesh_size.nx3 > 1)  ek += NGHOST, sk -= NGHOST;
  // nb=0-3 for inner boundary; nb=4-7 for outer boundary
  switch (nb) {
    case 0:
      sj = pmb->je - jo - (NGHOST - 1); ej = pmb->je;
      if (jo > nx2) sj = pmb->js;
      break;
    case 1:
      sj = pmb->js; ej = pmb->je - jo + NGHOST;
      if (jo < NGHOST) ej = pmb->je;
      break;
    case 2:
      sj = pmb->je - (NGHOST - 1); ej = pmb->je;
      if (jo > nx2) sj = pmb->je - (jo - nx2) + 1;
      break;
    case 3:
      sj = pmb->js; ej = pmb->js + (NGHOST - 1);
      if (jo < NGHOST) ej = pmb->js + (NGHOST - jo) - 1;
      break;
    case 4:
      sj = pmb->js; ej = pmb->js + jo + NGHOST - 1;
      if (jo > nx2) ej = pmb->je;
      break;
    case 5:
      sj = pmb->js + jo - NGHOST; ej = pmb->je;
      if (jo < NGHOST) sj = pmb->js;
      break;
    case 6:
      sj = pmb->js; ej = pmb->js + (NGHOST - 1);
      if (jo > nx2) ej = pmb->js + (jo - nx2) - 1;
      break;
    case 7:
      sj = pmb->je - (NGHOST - 1); ej = pmb->je;
      if (jo < NGHOST) sj = pmb->je - (NGHOST - jo) + 1;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in CellCenteredBoundaryVariable:LoadShearing "
          << std::endl << "nb = " << nb << " not valid" << std::endl;
      ATHENA_ERROR(msg);
  }
  int p = 0;
  BufferUtility::PackData(src, buf, 0, NHYDRO-1, si, ei, sj, ej, sk, ek, p);
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::SendShearingBoxBoundaryBuffersForInit()
//  \brief Send shearing box boundary buffers for hydro variables

void CellCenteredBoundaryVariable::SendShearingBoxBoundaryBuffersForInit() {
  MeshBlock *pmb = pmy_block_;
  Mesh *pmesh = pmb->pmy_mesh;
  AthenaArray<Real> &var = *var_cc;

  // KGF: hidden assumption that 2D?
  int jl = pmb->js - NGHOST;
  int ju = pmb->je + NGHOST;
  int kl = pmb->ks;
  int ku = pmb->ke;
  if (pmesh->mesh_size.nx3 > 1) {
    kl -= NGHOST;
    ku += NGHOST;
  }

  Real qomL = pbval_->qomL_;

  if (pbval_->is_shear[0]) {
    int ib = 0;
    int sign = +1;
    // step 1. -- add shear to the inner periodic boundary values
    for (int k=kl; k<=ku; k++) {
      for (int j=jl; j<=ju; j++) {
        for (int i=0; i<NGHOST; i++) {
          int ii = ib + i;
          // add shear to conservative
          shear_cc_[0](IM2,k,j,i) = var(IM2,k,j,ii) + sign*qomL*var(IDN,k,j,ii);
          // Update energy, then x2 momentum
          if (NON_BAROTROPIC_EOS) {
            var(IEN,k,j,ii) += (0.5/var(IDN,k,j,ii))*(SQR(shear_cc_[0](IM2,k,j,i))
                                                    - SQR(var(IM2,k,j,ii)));
          } // update energy
          var(IM2,k,j,ii) = shear_cc_[0](IM2,k,j,i);
        }
      }
    }
  }

  if (pbval_->is_shear[1]) {
    int ib = pmb->ie + 1;
    int sign = -1;
    // step 2. -- add shear to the outer periodic boundary values
    for (int k=kl; k<=ku; k++) {
      for (int j=jl; j<=ju; j++) {
        for (int i=0; i<NGHOST; i++) {
          int ii = ib + i;
          // add shear to conservative
          shear_cc_[1](IM2,k,j,i) = var(IM2,k,j,ii) + sign*qomL*var(IDN,k,j,ii);
          // Update energy, then x2 momentum
          if (NON_BAROTROPIC_EOS) {
            var(IEN,k,j,ii) += (0.5/var(IDN,k,j,ii))*(SQR(shear_cc_[1](IM2,k,j,i))
                                                      - SQR(var(IM2,k,j,ii)));
          }
          var(IM2,k,j,ii) = shear_cc_[1](IM2,k,j,i);
        }
      }
    }
  }
  return;
}
// --------------------------------------------------------------------------------------
// ! \fn void CellCenteredBoundaryVariable::SendShearingBoxBoundaryBuffers()
//  \brief Send shearing box boundary buffers for hydro variables

void CellCenteredBoundaryVariable::SendShearingBoxBoundaryBuffers() {
  MeshBlock *pmb = pmy_block_;
  Mesh *pmesh = pmb->pmy_mesh;
  AthenaArray<Real> &var = *var_cc;

  // KGF: hidden assumption that 2D?
  int jl = pmb->js - NGHOST;
  int ju = pmb->je + NGHOST;
  int kl = pmb->ks;
  int ku = pmb->ke;
  if (pmesh->mesh_size.nx3 > 1) {
    kl -= NGHOST;
    ku += NGHOST;
  }
  int js = pmb->js;
  int je = pmb->je;

  Real eps = pbval_->eps_;
  Real qomL = pbval_->qomL_;
  int ssize = NHYDRO;

  if (pbval_->is_shear[0]) {
    int ib = pmb->is - NGHOST;
    int ii;
    // step 1. -- load shboxvar_cc_
    for (int k=kl; k<=ku; k++) {
      for (int j=jl; j<=ju; j++) {
        for (int i=0; i<NGHOST; i++) {
          ii = ib + i;
          shear_cc_[0](IDN,k,j,i) = var(IDN,k,j,ii);
          shear_cc_[0](IM1,k,j,i) = var(IM1,k,j,ii);
          shear_cc_[0](IM2,k,j,i) = var(IM2,k,j,ii) + qomL*var(IDN,k,j,ii);
          shear_cc_[0](IM3,k,j,i) = var(IM3,k,j,ii);
          if (NON_BAROTROPIC_EOS) {
            shear_cc_[0](IEN,k,j,i) = var(IEN,k,j,ii) + (0.5/var(IDN,k,j,ii))
                                               *(SQR(shear_cc_[0](IM2,k,j,i))
                                                 - SQR(var(IM2,k,j,ii)));
          }
        }
      }
    }

    // step 2. -- conservative remaping
    for (int n=0; n<NHYDRO; n++) {
      for (int k=kl; k<=ku; k++) {
        for (int i=0; i<NGHOST; i++) {
          RemapFlux(n, k, js, je+2, i, eps, shear_cc_[0], shear_flx_cc_[0]);
          for (int j=js; j<=je+1; j++) {
            shear_cc_[0](n,k,j,i) -= shear_flx_cc_[0](j+1) - shear_flx_cc_[0](j);
          }
        }
      }
    }

    // step 3. -- load sendbuf; memcpy to recvbuf if on same rank, else post MPI_Isend
    for (int n=0; n<4; n++) {
      SimpleNeighborBlock& snb = pbval_->shear_send_neighbor_[0][n];
      if (snb.rank != -1) {
        LoadShearing(shear_cc_[0], shear_bd_var_[0].send[n], n);
        if (snb.rank == Globals::my_rank) {// on the same process
          CopyShearBufferSameProcess(snb, shear_send_count_cc_[0][n]*ssize, n, 0);
        } else { // MPI
#ifdef MPI_PARALLEL
          int tag = pbval_->CreateBvalsMPITag(snb.lid, n, shear_cc_phys_id_);
          MPI_Isend(shear_bd_var_[0].send[n], shear_send_count_cc_[0][n]*ssize,
                    MPI_ATHENA_REAL, snb.rank, tag, MPI_COMM_WORLD,
                    &shear_bd_var_[0].req_send[n]);
#endif
        }
      }
    }
  } // inner boundaries

  if (pbval_->is_shear[1]) {
    int  ib = pmb->ie + 1;
    qomL = -qomL;
    int ii;
    // step 1. -- load shboxvar_cc_
    for (int k=kl; k<=ku; k++) {
      for (int j=jl; j<=ju; j++) {
        for (int i=0; i<NGHOST; i++) {
          ii = ib+i;
          shear_cc_[1](IDN,k,j,i) = var(IDN,k,j,ii);
          shear_cc_[1](IM1,k,j,i) = var(IM1,k,j,ii);
          shear_cc_[1](IM2,k,j,i) = var(IM2,k,j,ii) + qomL*var(IDN,k,j,ii);
          shear_cc_[1](IM3,k,j,i) = var(IM3,k,j,ii);
          if (NON_BAROTROPIC_EOS) {
            shear_cc_[1](IEN,k,j,i) = var(IEN,k,j,ii)
                                               + (0.5/var(IDN,k,j,ii))
                                               *(SQR(shear_cc_[1](IM2,k,j,i))
                                                 - SQR(var(IM2,k,j,ii)));
          }
        }
      }
    }

    // step 2. -- conservative remaping
    for (int n=0; n<NHYDRO; n++) {
      for (int k=kl; k<=ku; k++) {
        for (int i=0; i<NGHOST; i++) {
          RemapFlux(n, k, js-1, je+1, i, -eps, shear_cc_[1], shear_flx_cc_[1]);
          for (int j=js-1; j<=je; j++) {
            shear_cc_[1](n,k,j,i) -= shear_flx_cc_[1](j+1) - shear_flx_cc_[1](j);
          }
        }
      }
    }

    // step 3. -- load sendbuf; memcpy to recvbuf if on same rank, post
    // MPI_Isend otherwise
    int offset = 4;
    for (int n=0; n<4; n++) {
      SimpleNeighborBlock& snb = pbval_->shear_send_neighbor_[1][n];
      if (snb.rank != -1) {
        LoadShearing(shear_cc_[1], shear_bd_var_[1].send[n], n+offset);
        if (snb.rank == Globals::my_rank) {
          CopyShearBufferSameProcess(snb, shear_send_count_cc_[1][n]*ssize, n, 1);
        } else { // MPI
#ifdef MPI_PARALLEL
          // bufid for outer(inner): 2(0) and 3(1)
          int tag = pbval_->CreateBvalsMPITag(snb.lid, n+offset, shear_cc_phys_id_);
          MPI_Isend(shear_bd_var_[1].send[n], shear_send_count_cc_[1][n]*ssize,
                    MPI_ATHENA_REAL, snb.rank, tag, MPI_COMM_WORLD,
                    &shear_bd_var_[1].req_send[n]);
#endif
        }
      }
    }
  } // outer boundaries
  return;
}

// --------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::SetShearingBoxBoundarySameLevel(Real *buf,
//                                                                         const int nb)
//  \brief Set hydro shearing box boundary received from a block on the same level

// KGF: AthenaArray<Real> &dst= pmb->phydro->u passed through from
// ReceiveHydroShearingboxBoundaryBuffers()

void CellCenteredBoundaryVariable::SetShearingBoxBoundarySameLevel(Real *buf,
                                                                   const int nb) {
  MeshBlock *pmb = pmy_block_;
  Mesh *pmesh = pmb->pmy_mesh;
  int si, sj, sk, ei, ej, ek;
  int jo = pbval_->joverlap_;
  int nx2 = pmb->block_size.nx2 - NGHOST;
  int nxo = pmb->block_size.nx2 - jo;

  sk = pmb->ks; ek = pmb->ke;
  if (pmesh->mesh_size.nx3 > 1) ek += NGHOST, sk -= NGHOST;
  // nb=0-3 for inner boundary; 4-7 for outer boundary.
  switch (nb) {
    case 0:
      si = pmb->is - NGHOST; ei = pmb->is - 1;
      sj = pmb->js - NGHOST; ej = pmb->js + (jo - 1);
      if (jo > nx2) sj = pmb->js - nxo;
      break;
    case 1:
      si = pmb->is - NGHOST; ei = pmb->is - 1;
      sj = pmb->js + jo; ej = pmb->je + NGHOST;
      if (jo < NGHOST) ej = pmb->je + jo;
      break;
    case 2:
      si = pmb->is - NGHOST; ei = pmb->is - 1;
      sj = pmb->js - NGHOST; ej = pmb->js - 1;
      if (jo > nx2) ej = pmb->js - nxo - 1;
      break;
    case 3:
      si = pmb->is - NGHOST; ei = pmb->is - 1;
      sj = pmb->je + jo + 1; ej = pmb->je + NGHOST;
      break;
    case 4:
      si = pmb->ie + 1; ei = pmb->ie + NGHOST;
      sj = pmb->je - (jo - 1); ej = pmb->je + NGHOST;
      if (jo > nx2) ej = pmb->je + nxo;
      break;
    case 5:
      si = pmb->ie + 1; ei = pmb->ie + NGHOST;
      sj = pmb->js - NGHOST; ej = pmb->je - jo;
      if (jo < NGHOST)   sj = pmb->js - jo;
      break;
    case 6:
      si = pmb->ie + 1; ei = pmb->ie + NGHOST;
      sj = pmb->je + 1; ej = pmb->je + NGHOST;
      if (jo > nx2) sj = pmb->je + nxo + 1;
      break;
    case 7:
      si = pmb->ie + 1; ei = pmb->ie + NGHOST;
      sj = pmb->js - NGHOST; ej = pmb->js - jo-1;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in CellCenteredBoundaryVariable:SetShearing " << std::endl
          << "nb = " << nb << " not valid" << std::endl;
      ATHENA_ERROR(msg);
  }

  // set [sj:ej] of current meshblock
  int p = 0;
  BufferUtility::UnpackData(buf, *var_cc, 0, NHYDRO-1, si, ei, sj, ej, sk, ek, p);
  return;
}


// --------------------------------------------------------------------------------------
// ! \fn bool CellCenteredBoundaryVariable::ReceiveShearingBoxBoundaryBuffers()
//  \brief receive shearing box boundary data for hydro variables

bool CellCenteredBoundaryVariable::ReceiveShearingBoxBoundaryBuffers() {
  bool flagi = true, flago = true;

  if (pbval_->is_shear[0]) { // check inner boundaries
    for (int n=0; n<4; n++) {
      if (shear_bd_var_[0].flag[n] == BoundaryStatus::completed) continue;
      if (shear_bd_var_[0].flag[n] == BoundaryStatus::waiting) {
        // on the same process
        if (pbval_->shear_recv_neighbor_[0][n].rank == Globals::my_rank) {
          flagi = false;
          continue;
        } else { // MPI boundary
#ifdef MPI_PARALLEL
          int test;
          MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &test,
                     MPI_STATUS_IGNORE);
          MPI_Test(&shear_bd_var_[0].req_recv[n], &test, MPI_STATUS_IGNORE);
          if (static_cast<bool>(test) == false) {
            flagi = false;
            continue;
          }
          shear_bd_var_[0].flag[n] = BoundaryStatus::arrived;
#endif
        }
      }
      // set var if boundary arrived
      SetShearingBoxBoundarySameLevel(shear_bd_var_[0].recv[n], n);
      shear_bd_var_[0].flag[n] = BoundaryStatus::completed; // completed
    } // loop over recv[0] to recv[3]
  } // inner boundary

  if (pbval_->is_shear[1]) { // check outer boundaries
    int offset = 4;
    for (int n=0; n<4; n++) {
      if (shear_bd_var_[1].flag[n] == BoundaryStatus::completed) continue;
      if (shear_bd_var_[1].flag[n] == BoundaryStatus::waiting) {
        // on the same process
        if (pbval_->shear_recv_neighbor_[1][n].rank == Globals::my_rank) {
          flago = false;
          continue;
        } else { // MPI boundary
#ifdef MPI_PARALLEL
          int test;
          MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &test,
                     MPI_STATUS_IGNORE);
          MPI_Test(&shear_bd_var_[1].req_recv[n], &test, MPI_STATUS_IGNORE);
          if (static_cast<bool>(test) == false) {
            flago = false;
            continue;
          }
          shear_bd_var_[1].flag[n] = BoundaryStatus::arrived;
#endif
        }
      }
      SetShearingBoxBoundarySameLevel(shear_bd_var_[1].recv[n], n+offset);
      shear_bd_var_[1].flag[n] = BoundaryStatus::completed;  // completed
    }
  } // outer boundary
  return (flagi && flago);
}


//--------------------------------------------------------------------------------------

//  \brief compute the flux along j indices for remapping adopted from 2nd
//  order RemapFlux of Athena 4.2 (C-version)

void CellCenteredBoundaryVariable::RemapFlux(const int n, const int k, const int jinner,
                               const int jouter, const int i, const Real eps,
                               const AthenaArray<Real> &var, AthenaArray<Real> &flux) {
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
    dUc = var(n,k,j+1,i) - var(n,k,j-1,i);
    dUl = var(n,k,j,  i) - var(n,k,j-1,i);
    dUr = var(n,k,j+1,i) - var(n,k,j,  i);

    dUm = 0.0;
    if (dUl*dUr > 0.0) {
      lim_slope = std::min(std::fabs(dUl), std::fabs(dUr));
      dUm = SIGN(dUc)*std::min(0.5*std::fabs(dUc), 2.0*lim_slope);
    }

    if (eps > 0.0) { // eps always > 0 for inner i boundary
      flux(j+1) = eps*(var(n,k,j,i) + 0.5*(1.0 - eps)*dUm);
    } else {         // eps always < 0 for outer i boundary
      flux(j  ) = eps*(var(n,k,j,i) - 0.5*(1.0 + eps)*dUm);
    }
  }
  return;
}
