//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_shear_fc.cpp
//! \brief functions that apply shearing box BCs for face-centered quantities
//========================================================================================

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
#include "../../orbital_advection/orbital_advection.hpp"
#include "../../parameter_input.hpp"
#include "../../utils/buffer_utils.hpp"
#include "../bvals.hpp"
#include "../bvals_interfaces.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif


//----------------------------------------------------------------------------------------
//! \fn int FaceCenteredBoundaryVariable::LoadShearingBoxBoundarySameLevel(
//!                                              FaceField &src, Real *buf, int nb)
//! \brief Load shearing box field boundary buffers

void FaceCenteredBoundaryVariable::LoadShearingBoxBoundarySameLevel(
                                               FaceField &src, Real *buf, int nb) {
  // TODO(felker): deduplicate with CellCenteredBoundaryVariable::LoadShearing()
  // Only differences are the calculation of psj, pej, and 3x PackData calls
  MeshBlock *pmb = pmy_block_;
  Mesh *pmesh = pmb->pmy_mesh;
  int si, sj, sk, ei, ej, ek;
  int jo = pbval_->joverlap_;
  int *jmin1 = pbval_->sb_data_[0].jmin_send;
  int *jmin2 = pbval_->sb_data_[1].jmin_send;
  int *jmax1 = pbval_->sb_data_[0].jmax_send;
  int *jmax2 = pbval_->sb_data_[1].jmax_send;
  sk = pmb->ks;        ek = pmb->ke;
  if (pmesh->mesh_size.nx3 > 1)  ek += NGHOST, sk -= NGHOST;

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
    msg << "### FATAL ERROR in FaceCenteredBoundaryVariable::"
        << "LoadShearingBoxBoundarySameLevel"<<std::endl
        << "nb = " << nb << " not valid" << std::endl;
    ATHENA_ERROR(msg);
  }

  int p = 0;
  // bx2
  for (int k=sk; k<=ek; k++) {
    for (int i=si; i<=ei; i++) {
      for (int j=sj; j<=ej+1; j++)
        buf[p++] = src.x2f(k,j,i);
    }
  }
  // bx3
  for (int k=sk; k<=ek+1; k++) {
    for (int i=si; i<=ei; i++) {
      for (int j=sj; j<=ej; j++)
        buf[p++] = src.x3f(k,j,i);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::SendShearingBoxBoundaryBuffers()
//! \brief Send shearing box boundary buffers for field variables

void FaceCenteredBoundaryVariable::SendShearingBoxBoundaryBuffers() {
  MeshBlock *pmb = pmy_block_;
  FaceField &var = *var_fc;
  int offset[2]{0, 4};

  for (int upper=0; upper<2; upper++) {
    if (pbval_->is_shear[upper]) {
      for (int n=0; n<4; n++) {
        SimpleNeighborBlock& snb = pbval_->sb_data_[upper].send_neighbor[n];
        if (snb.rank != -1) {
          LoadShearingBoxBoundarySameLevel(var, shear_bd_var_[upper].send[n],
                                           n+offset[upper]);
          if (snb.rank == Globals::my_rank) {
            CopyShearBufferSameProcess(snb, shear_send_count_fc_[upper][n], n, upper);
          } else { // MPI
#ifdef MPI_PARALLEL
            int tag = pbval_->CreateBvalsMPITag(snb.lid, n+offset[upper],
                                                shear_fc_phys_id_);
            MPI_Isend(shear_bd_var_[upper].send[n], shear_send_count_fc_[upper][n],
                      MPI_ATHENA_REAL, snb.rank, tag, MPI_COMM_WORLD,
                      &shear_bd_var_[upper].req_send[n]);
#endif
          }
        }
      }
    }  // if boundary is shearing
  }  // loop over inner/outer boundaries
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::SetShearingBoxBoundarySameLevel(
//!                         AthenaArray<Real> &src, Real *buf, const int nb)
//! \brief Set field shearing box boundary received from a block on the same level
//!
//! TODO(felker): deduplicate with CellCenteredBoundaryVariable::SetShearingBoxBound...()
//! Only differences are the calculation of psi,pei,psj,pej, and 3x UnpackData calls

void FaceCenteredBoundaryVariable::SetShearingBoxBoundarySameLevel(
                           FaceField &dst, Real *buf, const int nb) {
  MeshBlock *pmb = pmy_block_;
  Mesh *pmesh = pmb->pmy_mesh;
  const int& xgh = pbval_->xgh_;
  int si, sj, sk, ei, ej, ek;
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
    msg << "### FATAL ERROR in FaceCenteredBoundaryVariable::"
        << "SetShearingBoxBoundarySameLevel"<<std::endl
        << "nb = " << nb << " not valid" << std::endl;
    ATHENA_ERROR(msg);
  }

  int p = 0;
  // bx2
  for (int k=sk; k<=ek; k++) {
    for (int i=si; i<=ei; i++) {
#pragma omp simd
      for (int j=sj; j<=ej+1; j++) {
        dst.x2f(k,i,j) = buf[p++];
      }
    }
  }

  // bx3
  for (int k=sk; k<=ek+1; k++) {
    for (int i=si; i<=ei; i++) {
#pragma omp simd
      for (int j=sj; j<=ej; j++) {
        dst.x3f(k,i,j) = buf[p++];
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn bool FaceCenteredBoundaryVariable::ReceiveShearingBoxBoundaryBuffers()
//! \brief receive shearing box boundary data for field(face-centered) variables
//!
//! \todo (felker):
//! - DRY. completely identical to CellCenteredBoundaryVariable implementation

bool FaceCenteredBoundaryVariable::ReceiveShearingBoxBoundaryBuffers() {
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
        // set dst if boundary arrived
        SetShearingBoxBoundarySameLevel(shear_fc_[upper],
                                        shear_bd_var_[upper].recv[n],
                                        n+nb_offset[upper]);
        shear_bd_var_[upper].flag[n] = BoundaryStatus::completed; // completed
      } // loop over recv[0] to recv[3]
    }  // if boundary is shearing
  }  // loop over inner/outer boundaries
  return (flag[0] && flag[1]);
}

//--------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::SetShearingBoxBoundaryBuffers()
//! \brief receive shearing box boundary data for field(face-centered) variables
void FaceCenteredBoundaryVariable::SetShearingBoxBoundaryBuffers() {
  MeshBlock *pmb = pmy_block_;
  Mesh *pmesh = pmb->pmy_mesh;
  const int& xgh = pbval_->xgh_;
  const int& xorder = pbval_->xorder_;
  AthenaArray<Real> &pflux = pbval_->pflux_;
  FaceField &var = *var_fc;
  Coordinates *pco = pmb->pcoord;
  OrbitalAdvection *porb = pmb->porb;
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
      // calculating remapping flux and update var
      Real eps = (1.0-2*upper)*pbval_->eps_;
      // bx2
      AthenaArray<Real> &bf2 = shear_fc_[upper].x2f;
      for (int k=kl; k<=ku; k++) {
        for (int i=0; i<NGHOST; i++) {
          int ii = ib[upper]+i;
          if (xorder<=2) {
            porb->RemapFluxPlm(pflux, bf2, eps, 1-upper, k, i, jl, ju+2, xgh);
          } else {
            porb->RemapFluxPpm(pflux, bf2, eps, 1-upper, k, i, jl, ju+2, xgh);
          }
          const int shift = xgh+1-upper;
          for (int j=jl; j<=ju+1; j++) {
            var.x2f(k,j,ii) = bf2(k,i,j+shift) - (pflux(j+1) - pflux(j));
          }
        }
      }
      // bx3
      AthenaArray<Real> &bf3 = shear_fc_[upper].x3f;
      for (int k=kl; k<=ku+1; k++) {
        for (int i=0; i<NGHOST; i++) {
          int ii = ib[upper]+i;
          if (xorder<=2) {
            porb->RemapFluxPlm(pflux, bf3, eps, 1-upper, k, i, jl, ju+1, xgh);
          } else {
            porb->RemapFluxPpm(pflux, bf3, eps, 1-upper, k, i, jl, ju+1, xgh);
          }
          const int shift = xgh+1-upper;
          for (int j=jl; j<=ju; j++) {
            var.x3f(k,j,ii) = bf3(k,i,j+shift) - (pflux(j+1) - pflux(j));
          }
        }
      }
      // bx1
      // calculating bx1 using div B = 0
      // present way is only for cartesian
      if (upper==0) {
        for (int k=kl; k<=ku; k++) {
          const Real& dx3 = pco->dx3f(k);
          for (int j=jl; j<=ju; j++) {
            const Real& dx2 = pco->dx2f(j);
            for (int i=0; i<NGHOST; i++) {
              int ii = pmb->is-1-i;
              const Real& dx1 = pco->dx1f(ii);
              var.x1f(k,j,ii) =
                 var.x1f(k,j,ii+1) +
                 dx1/dx2*(var.x2f(k,j+1,ii)-var.x2f(k,j,ii)) +
                 dx1/dx3*(var.x3f(k+1,j,ii)-var.x3f(k,j,ii));
            }
          }
        }
      } else {
        for (int k=kl; k<=ku; k++) {
          const Real& dx3 = pco->dx3f(k);
          for (int j=jl; j<=ju; j++) {
            const Real& dx2 = pco->dx2f(j);
            for (int i=0; i<NGHOST; i++) {
              int ii = pmb->ie+1+i;
              const Real& dx1 = pco->dx1f(ii);
              var.x1f(k,j,ii+1) =
                 var.x1f(k,j,ii) -
                 dx1/dx2*(var.x2f(k,j+1,ii)-var.x2f(k,j,ii)) -
                 dx1/dx3*(var.x3f(k+1,j,ii)-var.x3f(k,j,ii));
            }
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn?void FaceCenteredBoundaryVariable::StartReceivingShear(BoundaryCommSubset phase)
//! \brief initiate MPI_Irecv

void FaceCenteredBoundaryVariable::StartReceivingShear(BoundaryCommSubset phase) {
#ifdef MPI_PARALLEL
  int size, tag;
  if (phase == BoundaryCommSubset::all) {
    int tag_offset1[2]{0, 3};
    for (int upper=0; upper<2; upper++) {
      if (pbval_->is_shear[upper]) {
        for (int n=0; n<3; n++) {
          int target_rank = pbval_->sb_flux_data_[upper].recv_neighbor[n].rank;
          if ((target_rank != Globals::my_rank) && (target_rank != -1)) {
            // emf
            size = shear_recv_count_emf_[upper][n];
            tag  = pbval_->CreateBvalsMPITag(pmy_block_->lid, n+tag_offset1[upper],
                                             shear_emf_phys_id_);
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
          // var_vc
          size = shear_recv_count_fc_[upper][n];
          tag  = pbval_->CreateBvalsMPITag(pmy_block_->lid, n+tag_offset2[upper],
                                           shear_fc_phys_id_);
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
