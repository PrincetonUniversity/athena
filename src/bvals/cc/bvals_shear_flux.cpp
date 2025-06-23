//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_shear_flux.cpp
//  \brief functions that apply shearing box BCs for fluxes
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
//! \fn void CellCenteredBoundaryVariable::LoadFluxShearingBoxBoundarySameLevel(
//!                               AthenaArray<Real> &src, Real *buf, int nb)
//! \brief Set shearing boundary flux buffers for sending to a block on the same level

void CellCenteredBoundaryVariable::LoadFluxShearingBoxBoundarySameLevel(
                               AthenaArray<Real> &src, Real *buf, int nb) {
  MeshBlock *pmb = pmy_block_;
  Mesh *pmesh = pmb->pmy_mesh;
  int sj, sk, ej, ek;
  int jo = pbval_->joverlap_flux_;
  int *jmin1 = pbval_->sb_flux_data_[0].jmin_send;
  int *jmin2 = pbval_->sb_flux_data_[1].jmin_send;
  int *jmax1 = pbval_->sb_flux_data_[0].jmax_send;
  int *jmax2 = pbval_->sb_flux_data_[1].jmax_send;
  sk = pmb->ks;        ek = pmb->ke;

  if (nb<3) { // inner boundary
    sj = jmin1[nb]-jo-1;
    ej = jmax1[nb]-jo-1;
  } else if (nb<6) { //outer boundary
    sj = jmin2[nb-3]+jo;
    ej = jmax2[nb-3]+jo;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in CellCenteredBoundaryVariable:"
        << "LoadFluxShearingBoxBoundarySameLevel"     << std::endl
        << "nb = " << nb << " not valid" << std::endl;
    ATHENA_ERROR(msg);
  }
  int p = 0;
  BufferUtility::PackData(src, buf, sj, ej, sk, ek, nl_, nu_, p);
  return;
}


//---------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::SendShearingBoxBoundaryBuffers()
//! \brief Send shearing box boundary flux buffers

void CellCenteredBoundaryVariable::SendFluxShearingBoxBoundaryBuffers() {
  MeshBlock *pmb = pmy_block_;
  Mesh *pmesh = pmb->pmy_mesh;
  int ssize = nu_ + 1;
  int offset[2]{0, 3};
  for (int upper=0; upper<2; upper++) {
    if (pbval_->is_shear[upper]) {
      for (int n=0; n<3; n++) {
        SimpleNeighborBlock& snb = pbval_->sb_flux_data_[upper].send_neighbor[n];
        if (snb.rank != -1) {
          LoadFluxShearingBoxBoundarySameLevel(shear_var_flx_[upper],
                                   shear_bd_flux_[upper].send[n], n+offset[upper]);
          if (snb.rank == Globals::my_rank) {// on the same process
            CopyShearFluxSameProcess(snb, shear_send_count_flx_[upper][n]*ssize, n,
                                       upper);
          } else { // MPI
#ifdef MPI_PARALLEL
            int tag = pbval_->CreateBvalsMPITag(snb.lid, n+offset[upper],
                                                shear_flx_phys_id_);
            MPI_Isend(shear_bd_flux_[upper].send[n],
                      shear_send_count_flx_[upper][n]*ssize,
                      MPI_ATHENA_REAL, snb.rank, tag, MPI_COMM_WORLD,
                      &shear_bd_flux_[upper].req_send[n]);
#endif
          }
        }
      }  // loop over recv[0] to recv[3]
    }  // if boundary is shearing
  }  // loop over inner/outer boundaries
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::SetShearingBoxBoundarySameLevel(
//!                        AthenaArray<Real> &src, Real *buf, const int nb)
//! \brief Set shearing box boundary flux data received from a block on the same level

void CellCenteredBoundaryVariable::SetFluxShearingBoxBoundarySameLevel(
                       AthenaArray<Real> &src, Real *buf, const int nb) {
  MeshBlock *pmb = pmy_block_;
  int sj, sk, ej, ek;
  const int& xgh = pbval_->xgh_;
  sk = pmb->ks;        ek = pmb->ke;

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
        << "SetFluxShearingBoxBoundarySameLevel"<<std::endl
        << "nb = " << nb << " not valid" << std::endl;
    ATHENA_ERROR(msg);
  }
  int p = 0;
  for (int n=nl_; n<=nu_; ++n) {
    for (int k=sk; k<=ek; k++) {
#pragma omp simd
      for (int j=sj; j<=ej; j++) {
        src(n,k,0,j) = buf[p++];
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn bool CellCenteredBoundaryVariable::ReceiveShearingBoxBoundaryBuffers()
//! \brief Receive shearing box boundary flux data

bool CellCenteredBoundaryVariable::ReceiveFluxShearingBoxBoundaryBuffers() {
  bool flag[2]{true, true};
  int nb_offset[2]{0, 3};
  for (int upper=0; upper<2; upper++) {
    if (pbval_->is_shear[upper]) { // check inner boundaries
      for (int n=0; n<3; n++) {
        if (shear_bd_flux_[upper].flag[n] == BoundaryStatus::completed) continue;
        if (shear_bd_flux_[upper].flag[n] == BoundaryStatus::waiting) {
          // on the same process
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
        // set var if boundary arrived
        SetFluxShearingBoxBoundarySameLevel(shear_map_flx_[upper],
                                            shear_bd_flux_[upper].recv[n],
                                            n+nb_offset[upper]);
        shear_bd_flux_[upper].flag[n] = BoundaryStatus::completed; // completed
      }  // loop over recv[0] to recv[2]
    }  // if boundary is shearing
  }  // loop over inner/outer boundaries
  return (flag[0] && flag[1]);
}

//----------------------------------------------------------------------------------------
//! \fn bool CellCenteredBoundaryVariable::SetFluxShearingBoxBoundaryBuffers()
//! \brief Update shearing box boundary fluxes

void CellCenteredBoundaryVariable::SetFluxShearingBoxBoundaryBuffers() {
  MeshBlock *pmb = pmy_block_;
  OrbitalAdvection *porb = pmb->porb;
  AthenaArray<Real> &pflux = pbval_->pflux_;
  const int& xgh = pbval_->xgh_;
  const int& xorder = pbval_->xorder_;
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  int ib[2]{is, ie+1};
  for (int upper=0; upper<2; upper++) {
    if (pbval_->is_shear[upper]) { // check inner boundaries
      Real eps = (1.0-2*upper)*pbval_->eps_flux_;
      int ii = ib[upper];
      for (int n=nl_; n<=nu_; n++) {
        pbuf.InitWithShallowSlice(shear_map_flx_[upper], 4, n, 1);
        for (int k=ks; k<=ke; k++) {
          if (xorder<=2) {
            porb->RemapFluxPlm(pflux, pbuf, eps, 1-upper, k, 0, js, je+1, xgh);
          } else {
            porb->RemapFluxPpm(pflux, pbuf, eps, 1-upper, k, 0, js, je+1, xgh);
          }
          const int shift = xgh+1-upper;
          for (int j=js; j<=je; j++) {
            x1flux(n,k,j,ii) = 0.5*(x1flux(n,k,j,ii)+pbuf(k,0,j+shift)
                               - (pflux(j+1) - pflux(j)));
          }
        }
      }
    }
  }
  return;
}
