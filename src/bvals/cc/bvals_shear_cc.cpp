//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_shear_cc.cpp
//! \brief functions that apply shearing box BCs for cell-centered variables
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
#include "../../orbital_advection/orbital_advection.hpp"
#include "../../parameter_input.hpp"
#include "../../utils/buffer_utils.hpp"
#include "../bvals.hpp"
#include "../bvals_interfaces.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif


//--------------------------------------------------------------------------------------
//! \fn int CellCenteredBoundaryVariable::LoadShearingBoxBoundarySameLevel(
//!                                       AthenaArray<Real> &src, Real *buf, int nb)
//! \brief Set CC shearing boundary buffers for sending to a block on the same level

void CellCenteredBoundaryVariable::LoadShearingBoxBoundarySameLevel(
                                 AthenaArray<Real> &src, Real *buf, int nb) {
  MeshBlock *pmb = pmy_block_;
  Mesh *pmesh = pmb->pmy_mesh;
  int si, sj, sk, ei, ej, ek;
  int jo = pbval_->joverlap_;
  int *jmin1 = pbval_->sb_data_[0].jmin_send;
  int *jmin2 = pbval_->sb_data_[1].jmin_send;
  int *jmax1 = pbval_->sb_data_[0].jmax_send;
  int *jmax2 = pbval_->sb_data_[1].jmax_send;
  sk = pmb->ks;        ek = pmb->ke;
  if (pmesh->mesh_size.nx3 > 1) ek += NGHOST, sk -= NGHOST;

  if (nb<4) { // inner boundary
    si = pmb->is-NGHOST;
    ei = pmb->is-1;
    sj = jmin1[nb]-jo-1;
    ej = jmax1[nb]-jo-1;
  } else if (nb<8) { //outer boundary
    si = pmb->ie+1;
    ei = pmb->ie+NGHOST;
    sj = jmin2[nb-4]+jo;
    ej = jmax2[nb-4]+jo;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in CellCenteredBoundaryVariable::"
        << "ShearingBoxBoundarySameLevel"<<std::endl
        << "nb = " << nb << " not valid" << std::endl;
    ATHENA_ERROR(msg);
  }
  int p = 0;
  for (int n=nl_; n<=nu_; ++n) {
    for (int k=sk; k<=ek; k++) {
      for (int i=si; i<=ei; i++) {
        for (int j=sj; j<=ej; j++) {
          buf[p++] = src(n,k,j,i);
        }
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::SendShearingBoxBoundaryBuffers()
//! \brief Send CC shearing box boundary buffers

void CellCenteredBoundaryVariable::SendShearingBoxBoundaryBuffers() {
  MeshBlock *pmb = pmy_block_;
  Mesh *pmesh = pmb->pmy_mesh;
  AthenaArray<Real> &var = *var_cc;
  int ssize = nu_ + 1;
  int offset[2]{0, 4};
  for (int upper=0; upper<2; upper++) {
    if (pbval_->is_shear[upper]) {
      for (int n=0; n<4; n++) {
        SimpleNeighborBlock& snb = pbval_->sb_data_[upper].send_neighbor[n];
        if (snb.rank != -1) {
          LoadShearingBoxBoundarySameLevel(var, shear_bd_var_[upper].send[n],
                                       n+offset[upper]);
          if (snb.rank == Globals::my_rank) {// on the same process
            CopyShearBufferSameProcess(snb, shear_send_count_cc_[upper][n]*ssize, n,
                                       upper);
          } else { // MPI
#ifdef MPI_PARALLEL
            int tag = pbval_->CreateBvalsMPITag(snb.lid, n+offset[upper],
                                                shear_cc_phys_id_);
            MPI_Isend(shear_bd_var_[upper].send[n], shear_send_count_cc_[upper][n]*ssize,
                      MPI_ATHENA_REAL, snb.rank, tag, MPI_COMM_WORLD,
                      &shear_bd_var_[upper].req_send[n]);
#endif
          }
        }
      }  // loop over recv[0] to recv[3]
    }  // if boundary is shearing
  }  // loop over inner/outer boundaries
  return;
}

// --------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::SetShearingBoxBoundarySameLevel(
//!                        AthenaArray<Real> &src, Real *buf, const int nb)
//! \brief Set CC shearing boundary data received from a block on the same level

void CellCenteredBoundaryVariable::SetShearingBoxBoundarySameLevel(
                    AthenaArray<Real> &src, Real *buf,const int nb) {
  MeshBlock *pmb = pmy_block_;
  Mesh *pmesh = pmb->pmy_mesh;
  int si, sj, sk, ei, ej, ek;
  const int& xgh = pbval_->xgh_;
  si = pmb->is-NGHOST; ei = pmb->is-1;
  sk = pmb->ks;        ek = pmb->ke;
  if (pmesh->mesh_size.nx3 > 1)  ek += NGHOST, sk -= NGHOST;

  int *jmin1 = pbval_->sb_data_[0].jmin_recv;
  int *jmin2 = pbval_->sb_data_[1].jmin_recv;
  int *jmax1 = pbval_->sb_data_[0].jmax_recv;
  int *jmax2 = pbval_->sb_data_[1].jmax_recv;

  if (nb<4) {
    sj = jmin1[nb]+xgh;     ej = jmax1[nb]+xgh;
  } else if (nb<8) {
    sj = jmin2[nb-4]+xgh;   ej = jmax2[nb-4]+xgh;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in CellCenteredBoundaryVariable::"
        << "SetShearingBoxBoundarySameLevel"<<std::endl
        << "nb = " << nb << " not valid" << std::endl;
    ATHENA_ERROR(msg);
  }
  int p = 0;
  BufferUtility::UnpackData(buf, src, nl_, nu_,
                            sj, ej, si, ei, sk, ek, p);
  return;
}


//----------------------------------------------------------------------------------------
//! \fn bool CellCenteredBoundaryVariable::ReceiveShearingBoxBoundaryBuffers()
//! \brief Receive CC shearing box boundary buffers

bool CellCenteredBoundaryVariable::ReceiveShearingBoxBoundaryBuffers() {
  bool flag[2]{true, true};
  int nb_offset[2]{0, 4};

  for (int upper=0; upper<2; upper++) {
    if (pbval_->is_shear[upper]) { // check inner boundaries
      for (int n=0; n<4; n++) {
        if (shear_bd_var_[upper].flag[n] == BoundaryStatus::completed) continue;
        if (shear_bd_var_[upper].flag[n] == BoundaryStatus::waiting) {
          // on the same process
          if (pbval_->sb_data_[upper].recv_neighbor[n].rank == Globals::my_rank) {
            flag[upper] = false;
            continue;
          } else { // MPI boundary
#ifdef MPI_PARALLEL
            int test;
            MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &test,
                       MPI_STATUS_IGNORE);
            MPI_Test(&shear_bd_var_[upper].req_recv[n], &test, MPI_STATUS_IGNORE);
            if (!static_cast<bool>(test)) {
              flag[upper] = false;
              continue;
            }
            shear_bd_var_[upper].flag[n] = BoundaryStatus::arrived;
#endif
          }
        }
        // set var if boundary arrived
        SetShearingBoxBoundarySameLevel(shear_cc_[upper], shear_bd_var_[upper].recv[n],
                                        n+nb_offset[upper]);
        shear_bd_var_[upper].flag[n] = BoundaryStatus::completed; // completed
      }  // loop over recv[0] to recv[3]
    }  // if boundary is shearing
  }  // loop over inner/outer boundaries
  return (flag[0] && flag[1]);
}

//----------------------------------------------------------------------------------------
//! \fn bool CellCenteredBoundaryVariable::SetShearingBoxBoundaryBuffers()
//! \brief Update CC shearing box boundary variables

void CellCenteredBoundaryVariable::SetShearingBoxBoundaryBuffers() {
  MeshBlock *pmb = pmy_block_;
  Mesh *pmesh = pmb->pmy_mesh;
  OrbitalAdvection *porb = pmb->porb;
  AthenaArray<Real> &var = *var_cc;
  AthenaArray<Real> &pflux = pbval_->pflux_;
  const int& xgh = pbval_->xgh_;
  const int& xorder = pbval_->xorder_;
  //int nb_offset[2]{0, 4};
  int ib[2]{pmb->is - NGHOST, pmb->ie + 1};
  int js = pmb->js, je = pmb->je;
  int kl = pmb->ks, ku = pmb->ke;
  if (pmesh->mesh_size.nx3 > 1) {
    kl -= NGHOST;
    ku += NGHOST;
  }
  int jl = js-NGHOST;
  int ju = je+NGHOST;

  for (int upper=0; upper<2; upper++) {
    if (pbval_->is_shear[upper]) { // check inner boundaries
      // step 1 -- (optionally) appy shear to shear_cc_ (does nothing by default)
      if (!porb->orbital_advection_defined)
        ShearQuantities(shear_cc_[upper], upper);  // Hydro overrides this
      // step 2. -- calculating remapping flux and update var
      Real eps = (1.0-2*upper)*pbval_->eps_;
      for (int n=nl_; n<=nu_; n++) {
        pbuf.InitWithShallowSlice(shear_cc_[upper], 4, n, 1);
        for (int k=kl; k<=ku; k++) {
          for (int i=0; i<NGHOST; i++) {
            int ii = ib[upper]+i;
            if (xorder<=2) {
              porb->RemapFluxPlm(pflux, pbuf, eps, 1-upper, k, i, jl, ju+1, xgh);
            } else {
              porb->RemapFluxPpm(pflux, pbuf, eps, 1-upper, k, i, jl, ju+1, xgh);
            }
            const int shift = xgh+1-upper;
            for (int j=jl; j<=ju; j++) {
              var(n,k,j,ii) = pbuf(k,i,j+shift) - (pflux(j+1) - pflux(j));
            }
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::StartReceivingShear(BoundaryCommSubset phase)
//! \brief
void CellCenteredBoundaryVariable::StartReceivingShear(BoundaryCommSubset phase) {
#ifdef MPI_PARALLEL
  int size, tag;
  if (phase == BoundaryCommSubset::all) {
    int tag_offset1[2]{0, 3};
    for (int upper=0; upper<2; upper++) {
      if (pbval_->is_shear[upper]) {
        for (int n=0; n<3; n++) {
          int target_rank = pbval_->sb_flux_data_[upper].recv_neighbor[n].rank;
          if ((target_rank != Globals::my_rank) && (target_rank != -1)) {
            size = (nu_ + 1)*shear_recv_count_flx_[upper][n];
            tag  = pbval_->CreateBvalsMPITag(pmy_block_->lid, n+tag_offset1[upper],
                                             shear_flx_phys_id_);
            MPI_Irecv(shear_bd_flux_[upper].recv[n], size, MPI_ATHENA_REAL,
                      target_rank, tag, MPI_COMM_WORLD,
                      &shear_bd_flux_[upper].req_recv[n]);
          }
        }
      }
    }
  }
  int tag_offset2[2]{0, 4};
  for (int upper=0; upper<2; upper++) {
    if (pbval_->is_shear[upper]) {
      for (int n=0; n<4; n++) {
        int target_rank = pbval_->sb_data_[upper].recv_neighbor[n].rank;
        if ((target_rank != Globals::my_rank) && (target_rank != -1)) {
          size = (nu_ + 1)*shear_recv_count_cc_[upper][n];
          tag  = pbval_->CreateBvalsMPITag(pmy_block_->lid, n+tag_offset2[upper],
                                           shear_cc_phys_id_);
          MPI_Irecv(shear_bd_var_[upper].recv[n], size, MPI_ATHENA_REAL,
                    target_rank, tag, MPI_COMM_WORLD,
                    &shear_bd_var_[upper].req_recv[n]);
        }
      }
    }
  }
#endif
  return;
}
