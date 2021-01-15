//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_shear_emf.cpp
//! \brief functions that apply BCs for face-centered flux corrections in shearing box
//! calculations
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


//----------------------------------------------------------------------------------------
//! \fn int FaceCenteredBoundaryVariable::LoadEMFShearingBoxBoundarySameLevel(
//!                                               EdgeField &src, Real *buf, int nb)
//! \brief Load shearing box EMF boundary buffers

void FaceCenteredBoundaryVariable::LoadEMFShearingBoxBoundarySameLevel(
                                   EdgeField &src, Real *buf, const int nb) {
  MeshBlock *pmb = pmy_block_;
  int sj, sk, ej, ek;
  int jo = pbval_->joverlap_flux_;
  int *jmin1 = pbval_->sb_flux_data_[0].jmin_send;
  int *jmin2 = pbval_->sb_flux_data_[1].jmin_send;
  int *jmax1 = pbval_->sb_flux_data_[0].jmax_send;
  int *jmax2 = pbval_->sb_flux_data_[1].jmax_send;
  sk = pmb->ks; ek = pmb->ke;

  if (nb<3) { // inner boundary
    sj = jmin1[nb]-jo-1;
    ej = jmax1[nb]-jo-1;
  } else if (nb<6) { //outer boundary
    sj = jmin2[nb-3]+jo;
    ej = jmax2[nb-3]+jo;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in CellCenteredBoundaryVariable::"
        << "LoadShearingBoxBoundarySameLevel"<<std::endl
        << "nb = " << nb << " not valid" << std::endl;
    ATHENA_ERROR(msg);
  }

  int p = 0;
  // pack e2
  for (int k=sk; k<=ek+1; k++) {
#pragma omp simd
    for (int j=sj; j<=ej; j++) {
      buf[p++] = src.x2e(k,j);
    }
  }
  // pack e3
  for (int k=sk; k<=ek; k++) {
#pragma omp simd
    for (int j=sj; j<=ej+1; j++) {
      buf[p++] = src.x3e(k,j);
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::SendEMFShearingBoxBoundaryCorrection()
//! \brief Send shearing box boundary buffers for EMF correction

void FaceCenteredBoundaryVariable::SendEMFShearingBoxBoundaryCorrection() {
  MeshBlock *pmb = pmy_block_;
  int offset[2]{0, 3};
  for (int upper=0; upper<2; upper++) {
    if (pbval_->is_shear[upper]) {
      // step 1. -- average edges of shboxvar_fc_flx_
      // average e3 for x1x2 edge
      for (int k=pmb->ks; k<=pmb->ke; k++) {
        for (int j=pmb->js; j<=pmb->je+1; j+=pmb->block_size.nx2)
          shear_var_emf_[upper].x3e(k,j) *= 0.5;
      }
      if(pmb->block_size.nx3 > 1) {
        // average e2 for x1x3 edge
        if (pbval_->block_bcs[BoundaryFace::inner_x3] == BoundaryFlag::block
            || pbval_->block_bcs[BoundaryFace::inner_x3] == BoundaryFlag::periodic) {
          for (int j=pmb->js; j<=pmb->je; j++)
            shear_var_emf_[upper].x2e(pmb->ks,j) *= 0.5;
        }
        if (pbval_->block_bcs[BoundaryFace::outer_x3] == BoundaryFlag::block
            || pbval_->block_bcs[BoundaryFace::outer_x3] == BoundaryFlag::periodic) {
          for (int j=pmb->js; j<=pmb->je; j++)
            shear_var_emf_[upper].x2e(pmb->ke+1,j) *= 0.5;
        }
      }

      // step 2. -- load sendbuf; memcpy to recvbuf if on same rank, post
      // MPI_Isend otherwise
      for (int n=0; n<3; n++) {
        SimpleNeighborBlock& snb = pbval_->sb_flux_data_[upper].send_neighbor[n];
        if (snb.rank != -1) {
          LoadEMFShearingBoxBoundarySameLevel(shear_var_emf_[upper],
                                     shear_bd_flux_[upper].send[n], n+offset[upper]);
          if (snb.rank == Globals::my_rank) {
            CopyShearFluxSameProcess(snb, shear_send_count_emf_[upper][n], n, upper);
          } else { // MPI
#ifdef MPI_PARALLEL
            int tag = pbval_->CreateBvalsMPITag(snb.lid, n+offset[upper],
                                                shear_emf_phys_id_);
            MPI_Isend(shear_bd_flux_[upper].send[n], shear_send_count_emf_[upper][n],
                      MPI_ATHENA_REAL, snb.rank, tag,
                      MPI_COMM_WORLD, &shear_bd_flux_[upper].req_send[n]);
#endif
          }
        }
      }
    }  // if boundary is shearing
  }  // loop over inner/outer boundaries
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::SetEMFShearingBoxBoundarySameLevel(
//!                                   EdgeField &dst, Real *buf, const int nb)
//! \brief Set EMF shearing box boundary received from a block on the same level

void FaceCenteredBoundaryVariable::SetEMFShearingBoxBoundarySameLevel(
                                   EdgeField &dst, Real *buf, const int nb) {
  MeshBlock *pmb = pmy_block_;
  int &xgh = pbval_->xgh_;
  int sj, sk, ej, ek;

  sk = pmb->ks; ek = pmb->ke;

  int *jmin1 = pbval_->sb_flux_data_[0].jmin_recv;
  int *jmin2 = pbval_->sb_flux_data_[1].jmin_recv;
  int *jmax1 = pbval_->sb_flux_data_[0].jmax_recv;
  int *jmax2 = pbval_->sb_flux_data_[1].jmax_recv;

  if (nb<3) {
    sj = jmin1[nb]+xgh;     ej = jmax1[nb]+xgh;
  } else if (nb<6) {
    sj = jmin2[nb-3]+xgh;   ej = jmax2[nb-3]+xgh;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in CellCenteredBoundaryVariable::"
        << "SetEMFShearingBoxBoundarySameLevel"<<std::endl
        << "nb = " << nb << " not valid" << std::endl;
    ATHENA_ERROR(msg);
  }

  int p = 0;
  // unpack e2
  for (int k=sk; k<=ek+1; k++) {
#pragma omp simd
    for (int j=sj; j<=ej; j++) {
      dst.x2e(k,0,j) = buf[p++];
    }
  }
  // unpack e3
  for (int k=sk; k<=ek; k++) {
#pragma omp simd
    for (int j=sj; j<=ej+1; j++) {
      dst.x3e(k,0,j) = buf[p++];
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn bool FaceCenteredBoundaryVariable::ReceiveEMFShearingBoxBoundaryCorrection()
//! \brief receive shearing box boundary data for EMF correction
//!
//! \todo (felker):
//! * DRY. Identical to Face/CellCentered impl. except for "emf" identifiers

bool FaceCenteredBoundaryVariable::ReceiveEMFShearingBoxBoundaryCorrection() {
  bool flag[2]{true, true};
  int nb_offset[2]{0, 3};
  for (int upper=0; upper<2; upper++) {
    if (pbval_->is_shear[upper]) {
      for (int n=0; n<3; n++) {
        if (shear_bd_flux_[upper].flag[n] == BoundaryStatus::completed) continue;
        if (shear_bd_flux_[upper].flag[n] == BoundaryStatus::waiting) {
          if (pbval_->sb_flux_data_[upper].recv_neighbor[n].rank == Globals::my_rank) {
            flag[upper] = false;
            continue;
          } else { // MPI boundary
#ifdef MPI_PARALLEL
            int test;
            MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &test,
                       MPI_STATUS_IGNORE);
            MPI_Test(&shear_bd_flux_[upper].req_recv[n], &test, MPI_STATUS_IGNORE);
            if (!static_cast<bool>(test)) {
              flag[upper] = false;
              continue;
            }
            shear_bd_flux_[upper].flag[n] = BoundaryStatus::arrived;
#endif
          }
        }
        // set dst if boundary arrived
        SetEMFShearingBoxBoundarySameLevel(
            shear_map_emf_[upper], shear_bd_flux_[upper].recv[n], n+nb_offset[upper]);
        shear_bd_flux_[upper].flag[n] = BoundaryStatus::completed; // completed
      }
    }  // if boundary is shearing
  }  // loop over inner/outer boundaries
  return (flag[0] && flag[1]);
}


//--------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::SetEMFShearingBoxBoundaryCorrection()
//! \brief Set EMF boundary received from a block on the finer level

void FaceCenteredBoundaryVariable::SetEMFShearingBoxBoundaryCorrection() {
  MeshBlock *pmb = pmy_block_;
  OrbitalAdvection *porb = pmb->porb;
  int &xgh = pbval_->xgh_;
  int &xorder = pbval_->xorder_;
  AthenaArray<Real> &pflux = pbval_->pflux_;
  AthenaArray<Real> &e2 = pmb->pfield->e.x2e;
  AthenaArray<Real> &e3 = pmb->pfield->e.x3e;
  int ks = pmb->ks, ke = pmb->ke;
  int js = pmb->js, je = pmb->je;
  int is = pmb->is, ie = pmb->ie;

  int ib[2]{is, ie + 1};
  for (int upper=0; upper<2; upper++) {
    if (pbval_->is_shear[upper]) {
      Real eps = (1.0-2*upper)*pbval_->eps_flux_;
      int ii = ib[upper];
      // ex2
      AthenaArray<Real> &pe2 = shear_map_emf_[upper].x2e;
      for (int k=ks; k<=ke+1; k++) {
        if (xorder<=2) {
          porb->RemapFluxPlm(pflux, pe2, eps, 1-upper, k, 0, js, je+1, xgh);
        } else {
          porb->RemapFluxPpm(pflux, pe2, eps, 1-upper, k, 0, js, je+1, xgh);
        }
        const int shift = xgh+1-upper;
        for (int j=js; j<=je; j++) {
          e2(k,j,ii) = 0.5*(e2(k,j,ii)+pe2(k,0,j+shift) - (pflux(j+1) - pflux(j)));
        }
      }
      //ex3
      AthenaArray<Real> &pe3 = shear_map_emf_[upper].x3e;
      for (int k=ks; k<=ke; k++) {
        if (xorder<=2) {
          porb->RemapFluxPlm(pflux, pe3, eps, 1-upper, k, 0, js, je+2, xgh);
        } else {
          porb->RemapFluxPpm(pflux, pe3, eps, 1-upper, k, 0, js, je+2, xgh);
        }
        const int shift = xgh+1-upper;
        for (int j=js; j<=je+1; j++) {
          e3(k,j,ii) = 0.5*(e3(k,j,ii)+pe3(k,0,j+shift) - (pflux(j+1) - pflux(j)));
        }
      }
      ClearEMFShearing(shear_var_emf_[upper]);
    }  // if boundary is shearing
  }  // loop over inner/outer boundaries
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::ClearEMFShearing()
//! \brief Clear the working array for EMFs on the surface/edge contacting with
//! a shearing periodic boundary

void FaceCenteredBoundaryVariable::ClearEMFShearing(EdgeField &work) {
  AthenaArray<Real> &e2 = work.x2e;
  AthenaArray<Real> &e3 = work.x3e;
  e2.ZeroClear();
  e3.ZeroClear();
  return;
}
