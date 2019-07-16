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
#include "../bvals.hpp"
#include "../bvals_interfaces.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif


//--------------------------------------------------------------------------------------
//! \fn int FaceCenteredBoundaryVariable::LoadEMFShearing(EdgeField &src,
//                                                        Real *buf, int nb)
//  \brief Load shearing box EMF boundary buffers

void FaceCenteredBoundaryVariable::LoadEMFShearing(EdgeField &src,
                                                   Real *buf, const int nb) {
  MeshBlock *pmb = pmy_block_;
  int sj, sk, ej, ek;
  int psj, pej; // indices for e3
  int jo = pbval_->joverlap_;
  int nx2 = pmb->block_size.nx2 - NGHOST;
  // Unlike non-emf FaceCentered::LoadShearing... these aren't extended by NGHOST in 3D:
  sk = pmb->ks; ek = pmb->ke;
  // psj, pej calculations seem to be identical to non-emf version, but no psi, pei calc.
  // also, manual packing loops in this fn versus PackData() calls in non-emf impl.
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

  int offset[2]{0, 4};
  for (int upper=0; upper<2; upper++) {
    if (pbval_->is_shear[upper]) {
      // step 1. -- average edges of shboxvar_fc_flx_
      // average e3 for x1x2 edge
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j+=nx2)
          shear_var_emf_[upper].x3e(k,j) *= 0.5;
      }
      // average e2 for x1x3 edge
      for (int k=ks; k<=ke+1; k+=nx3) {
        for (int j=js; j<=je; j++)
          shear_var_emf_[upper].x2e(k,j) *= 0.5;
      }

      // step 2. -- load sendbuf; memcpy to recvbuf if on same rank, post
      // MPI_Isend otherwise
      for (int n=0; n<4; n++) {
        SimpleNeighborBlock& snb = pbval_->shear_send_neighbor_[upper][n];
        if (snb.rank != -1) {
          LoadEMFShearing(shear_var_emf_[upper], shear_bd_emf_[upper].send[n],
                          n+offset[upper]);
          if (snb.rank == Globals::my_rank) {
            CopyShearEMFSameProcess(snb, shear_send_count_emf_[upper][n], n, upper);
          } else { // MPI
#ifdef MPI_PARALLEL
            int tag = pbval_->CreateBvalsMPITag(snb.lid, n+offset[upper],
                                                shear_emf_phys_id_);
            MPI_Isend(shear_bd_emf_[upper].send[n], shear_send_count_emf_[upper][n],
                      MPI_ATHENA_REAL, snb.rank, tag,
                      MPI_COMM_WORLD, &shear_bd_emf_[upper].req_send[n]);
#endif
          }
        }
      }
    }  // if boundary is shearing
  }  // loop over inner/outer boundaries
  return;
}

// --------------------------------------------------------------------------------------
// ! \fn void FaceCenteredBoundaryVariable::SetEMFShearingBoxBoundarySameLevel(
//                                   EdgeField &dst, Real *buf, const int nb)
//  \brief Set EMF shearing box boundary received from a block on the same level

void FaceCenteredBoundaryVariable::SetEMFShearingBoxBoundarySameLevel(EdgeField &dst,
                                                                      Real *buf,
                                                                      const int nb) {
  MeshBlock *pmb = pmy_block_;
  int sj, sk, ej, ek;
  int psj, pej;
  int jo = pbval_->joverlap_;
  int nx2 = pmb->block_size.nx2 - NGHOST;
  int nxo = pmb->block_size.nx2 - jo;

  // Unlike non-emf FaceCentered::SetShearing... these aren't extended by NGHOST in 3D:
  sk = pmb->ks; ek = pmb->ke;
  // psj, pej calculations seem to be identical to non-emf version, but no psi, pei calc.
  // also, manual unpacking loops in this fn versus UnpackData() calls in non-emf impl.
  switch (nb) {
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
      msg << "### FATAL ERROR in "
          << "FaceCenteredBoundaryVariable:SetEMFShearingBoxBoundarySameLevel\n"
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

// TODO(felker): DRY. Identical to Face/CellCentered impl. except for "emf" identifiers
bool FaceCenteredBoundaryVariable::ReceiveEMFShearingBoxBoundaryCorrection() {
  bool flag[2]{true, true};
  int nb_offset[2]{0, 4};
  for (int upper=0; upper<2; upper++) {
    if (pbval_->is_shear[upper]) {
      for (int n=0; n<4; n++) {
        if (shear_bd_emf_[upper].flag[n] == BoundaryStatus::completed) continue;
        if (shear_bd_emf_[upper].flag[n] == BoundaryStatus::waiting) {
          if (pbval_->shear_recv_neighbor_[upper][n].rank == Globals::my_rank) {
            flag[upper] = false;
            continue;
          } else { // MPI boundary
#ifdef MPI_PARALLEL
            int test;
            MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &test,
                       MPI_STATUS_IGNORE);
            MPI_Test(&shear_bd_emf_[upper].req_recv[n], &test, MPI_STATUS_IGNORE);
            if (!static_cast<bool>(test)) {
              flag[upper] = false;
              continue;
            }
            shear_bd_emf_[upper].flag[n] = BoundaryStatus::arrived;
#endif
          }
        }
        // set dst if boundary arrived
        SetEMFShearingBoxBoundarySameLevel(
            shear_map_emf_[upper], shear_bd_emf_[upper].recv[n], n+nb_offset[upper]);
        shear_bd_emf_[upper].flag[n] = BoundaryStatus::completed; // completed
      }
    }  // if boundary is shearing
  }  // loop over inner/outer boundaries
  return (flag[0] && flag[1]);
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

  int ib[2]{is, ie + 1};
  int sign[2]{1, -1};
  for (int upper=0; upper<2; upper++) {
    if (pbval_->is_shear[upper]) {
      Real eps = sign[upper]*pbval_->eps_;
      int jl_remap = js - upper;
      int ju_remap = je + 2 - upper;

      ClearEMFShearing(shear_var_emf_[upper]);
      // step 1.-- conservative remapping
      for (int k=ks; k<=ke+1; k++) {
        RemapFluxEMF(k, jl_remap, ju_remap, eps, shear_map_emf_[upper].x2e,
                     shear_flx_emf_[upper].x2e);
        for (int j=js; j<=je; j++) {
          shear_map_emf_[upper].x2e(k,j) -= shear_flx_emf_[upper].x2e(j+1)
                                        - shear_flx_emf_[upper].x2e(j);
        }
      }
      // step 2.-- average the EMF correction
      // average e2
      for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je; j++) {
          int ii = ib[upper];
          e2(k,j,ii) = 0.5*(e2(k,j,ii) + shear_map_emf_[upper].x2e(k,j));
        }
      }
      ClearEMFShearing(shear_map_emf_[upper]);
    }  // if boundary is shearing
  }  // loop over inner/outer boundaries
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
//! \fn void FaceCenteredBoundaryVariable::RemapFluxEMF(int k, int jinner, int jouter,
//                       Real eps, static AthenaArray<Real> &U, AthenaArray<Real> &Flux)
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
      lim_slope = std::min(std::abs(dUl), std::abs(dUr));
      dUm = SIGN(dUc)*std::min(0.5*std::abs(dUc), 2.0*lim_slope);
    }

    if (eps > 0.0) { // eps always > 0 for inner i boundary
      flux(j+1) = eps*(var(k,j) + 0.5*(1.0 - eps)*dUm);
    } else {         // eps always < 0 for outer i boundary
      flux(j  ) = eps*(var(k,j) - 0.5*(1.0 + eps)*dUm);
    }
  }
  return;
}
