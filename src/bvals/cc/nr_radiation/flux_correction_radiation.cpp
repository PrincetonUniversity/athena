//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file flux_correction_radiation.cpp
//! \brief
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
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../coordinates/coordinates.hpp"
#include "../../../eos/eos.hpp"
#include "../../../field/field.hpp"
#include "../../../globals.hpp"
#include "../../../hydro/hydro.hpp"
#include "../../../mesh/mesh.hpp"
#include "../../../nr_radiation/integrators/rad_integrators.hpp"
#include "../../../nr_radiation/radiation.hpp"
#include "../../../parameter_input.hpp"
#include "../../../utils/buffer_utils.hpp"
#include "../bvals_cc.hpp"
#include "./bvals_rad.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn int RadBoundaryVariable::LoadFluxBoundaryBufferSameLevel(Real *buf,
//!                                                  const NeighborBlock& nb)
//! \brief Set surface flux buffers for sending to a block on the same level

int RadBoundaryVariable::LoadFluxBoundaryBufferSameLevel(Real *buf,
                                                         const NeighborBlock& nb) {
  MeshBlock *pmb=pmy_block_;
  NRRadiation *prad = pmb->pnrrad;
  Real qomL = pbval_->qomL_;
  int p = 0;
  if (pbval_->shearing_box == 1 && nb.shear
      && (nb.fid == BoundaryFace::inner_x1 || nb.fid == BoundaryFace::outer_x1)) {
    int i;
    int sign;
    if (nb.fid == BoundaryFace::inner_x1) {
      i = pmb->is;
      sign = -1;
    } else {
      i = pmb->ie + 1;
      sign =  1;
    }
    // pack x1flux
    for (int k=pmb->ks; k<=pmb->ke; k++) {
      for (int j=pmb->js; j<=pmb->je; j++) {
        // convert flux due to velocity difference
        Real vx = pmb->phydro->w(IVX,k,j,i);
        Real vy = pmb->phydro->w(IVY,k,j,i);
        Real vz = pmb->phydro->w(IVZ,k,j,i);
        Real *mux = &(prad->mu(0,k,j,i,0));
        Real *muy = &(prad->mu(1,k,j,i,0));
        Real *muz = &(prad->mu(2,k,j,i,0));
        Real *flux_lab = &(x1flux(k,j,i,0));
        prad->pradintegrator->LabToCom(vx,vy,vz,mux,muy,muz,flux_lab,ir_cm_);
        flux_lab = &(ir_lab_(0));
        prad->pradintegrator->ComToLab(vx,vy+sign*qomL,vz,mux,muy,muz,ir_cm_,flux_lab);
        for (int nn=nl_; nn<=nu_; nn++) {
          buf[p++] = flux_lab[nn];
        }
      }
    }
  }
  return p;
}

//----------------------------------------------------------------------------------------
//! \fn int RadBoundaryVariable::LoadFluxBoundaryBufferToCoarser(Real *buf,
//!                                                        const NeighborBlock& nb)
//! \brief Set surface flux buffers for sending to a block on the coarser level

int RadBoundaryVariable::LoadFluxBoundaryBufferToCoarser(Real *buf,
                                                         const NeighborBlock& nb) {
  MeshBlock *pmb=pmy_block_;
  Coordinates *pco = pmb->pcoord;
  // cache pointers to surface area arrays (BoundaryBase protected variable)
  AthenaArray<Real> &sarea0 = pbval_->sarea_[0];
  AthenaArray<Real> &sarea1 = pbval_->sarea_[1];
  int p = 0;
  // x1 direction
  if (nb.fid == BoundaryFace::inner_x1 || nb.fid == BoundaryFace::outer_x1) {
    int i = pmb->is + (pmb->ie-pmb->is + 1)*nb.fid;
    if (pmb->block_size.nx3>1) { // 3D
      for (int k=pmb->ks; k<=pmb->ke; k+=2) {
        for (int j=pmb->js; j<=pmb->je; j+=2) {
          Real amm = pco->GetFace1Area(k,   j,   i);
          Real amp = pco->GetFace1Area(k,   j+1, i);
          Real apm = pco->GetFace1Area(k+1, j,   i);
          Real app = pco->GetFace1Area(k+1, j+1, i);
          Real tarea = amm + amp + apm + app;
          for (int nn=nl_; nn<=nu_; nn++) {
            buf[p++] = (x1flux(k  , j  , i, nn)*amm
                       + x1flux(k  , j+1, i, nn)*amp
                       + x1flux(k+1, j  , i, nn)*apm
                       + x1flux(k+1, j+1, i, nn)*app)/tarea;
          }
        }
      }
    } else if (pmb->block_size.nx2>1) { // 2D
      int k = pmb->ks;
      for (int j=pmb->js; j<=pmb->je; j+=2) {
        Real am = pco->GetFace1Area(k, j,   i);
        Real ap = pco->GetFace1Area(k, j+1, i);
        Real tarea = am + ap;
        for (int nn=nl_; nn<=nu_; nn++) {
          buf[p++] = (x1flux(k, j  , i, nn)*am + x1flux(k, j+1, i, nn)*ap)/tarea;
        }
      }
    } else { // 1D
      int k = pmb->ks, j = pmb->js;
      for (int nn=nl_; nn<=nu_; nn++)
        buf[p++] = x1flux(k, j, i, nn);
    }
    // x2 direction
  } else if (nb.fid == BoundaryFace::inner_x2 || nb.fid == BoundaryFace::outer_x2) {
    int j = pmb->js + (pmb->je-pmb->js + 1)*(nb.fid & 1);
    if (pmb->block_size.nx3>1) { // 3D
      for (int k=pmb->ks; k<=pmb->ke; k+=2) {
        pco->Face2Area(k  , j, pmb->is, pmb->ie, sarea0);
        pco->Face2Area(k+1, j, pmb->is, pmb->ie, sarea1);
        for (int i=pmb->is; i<=pmb->ie; i+=2) {
          Real tarea = sarea0(i) + sarea0(i+1) + sarea1(i) + sarea1(i+1);
          for (int nn=nl_; nn<=nu_; nn++) {
            buf[p++] = (x2flux(k  , j, i  , nn)*sarea0(i  )
                       + x2flux(k  , j, i+1, nn)*sarea0(i+1)
                       + x2flux(k+1, j, i  , nn)*sarea1(i  )
                       + x2flux(k+1, j, i+1, nn)*sarea1(i+1))/tarea;
          }
        }
      }
    } else if (pmb->block_size.nx2>1) { // 2D
      int k = pmb->ks;
      pco->Face2Area(0, j, pmb->is ,pmb->ie, sarea0);
      for (int i=pmb->is; i<=pmb->ie; i+=2) {
        Real tarea = sarea0(i) + sarea0(i+1);
        for (int nn=nl_; nn<=nu_; nn++) {
          buf[p++] = (x2flux(k, j, i  , nn)*sarea0(i  )
                     + x2flux(k, j, i+1, nn)*sarea0(i+1))/tarea;
        }
      }
    }
    // x3 direction - 3D onl_y
  } else if (nb.fid == BoundaryFace::inner_x3 || nb.fid == BoundaryFace::outer_x3) {
    int k = pmb->ks + (pmb->ke-pmb->ks + 1)*(nb.fid & 1);
      for (int j=pmb->js; j<=pmb->je; j+=2) {
      pco->Face3Area(k, j,   pmb->is, pmb->ie, sarea0);
      pco->Face3Area(k, j+1, pmb->is, pmb->ie, sarea1);
      for (int i=pmb->is; i<=pmb->ie; i+=2) {
        Real tarea = sarea0(i) + sarea0(i+1) + sarea1(i) + sarea1(i+1);
        for (int nn=nl_; nn<=nu_; nn++) {
          buf[p++] = (x3flux(k, j  , i, nn )*sarea0(i  )
                     + x3flux(k, j  , i+1, nn)*sarea0(i+1)
                     + x3flux(k, j+1, i, nn)*sarea1(i  )
                     + x3flux(k, j+1, i+1, nn)*sarea1(i+1))/tarea;
        }
      }
    }
  }
  return p;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::SendFluxCorrection()
//! \brief Send surface flux buffers

void RadBoundaryVariable::SendFluxCorrection() {
  MeshBlock *pmb=pmy_block_;
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (nb.ni.type != NeighborConnect::face) break;
    if (bd_var_flcor_.sflag[nb.bufid] == BoundaryStatus::completed) continue;
    int p = 0;
    if (nb.snb.level == pmb->loc.level) { // to same level
      p = LoadFluxBoundaryBufferSameLevel(bd_var_flcor_.send[nb.bufid],nb);
    } else if (nb.snb.level < pmb->loc.level) { // to coaser
      p = LoadFluxBoundaryBufferToCoarser(bd_var_flcor_.send[nb.bufid],nb);
    }
    // else { // to finer
    // }

    if (p>0) {
      if (nb.snb.rank == Globals::my_rank) // on the same node
        CopyFluxCorrectionBufferSameProcess(nb, p);
#ifdef MPI_PARALLEL
      else
        MPI_Start(&(bd_var_flcor_.req_send[nb.bufid]));
#endif
    } else {
      if (nb.snb.rank == Globals::my_rank) // on the same node
        SetCompletedFlagSameProcess(nb);
    }
    bd_var_flcor_.sflag[nb.bufid] = BoundaryStatus::completed;
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::SetFluxBoundarySameLevel(Real *buf,
//!                                                               const NeighborBlock& nb)
//! \brief Set surface flux data received from a block on the same level

void RadBoundaryVariable::SetFluxBoundarySameLevel(Real *buf,
                                                           const NeighborBlock& nb) {
  MeshBlock *pmb = pmy_block_;
  int p = 0;

  if (nb.fid == BoundaryFace::inner_x1) {
    for (int k=pmb->ks; k<=pmb->ke; k++) {
      for (int j=pmb->js; j<=pmb->je; j++) {
        for (int nn=nl_; nn<=nu_; nn++) {
          shear_var_flx_[0](k,j,nn) = buf[p++];
        }
      }
    }
  } else {
    for (int k=pmb->ks; k<=pmb->ke; k++) {
      for (int j=pmb->js; j<=pmb->je; j++) {
        for (int nn=nl_; nn<=nu_; nn++) {
          shear_var_flx_[1](k,j,nn) = buf[p++];
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::SetFluxBoundaryFromFiner(Real *buf,
//!                                                               const NeighborBlock& nb)
//! \brief Set surface flux data received from a block on the finer level

void RadBoundaryVariable::SetFluxBoundaryFromFiner(Real *buf,
                                                           const NeighborBlock& nb) {
  MeshBlock *pmb = pmy_block_;
  int p = 0;
  if (nb.fid == BoundaryFace::inner_x1 || nb.fid == BoundaryFace::outer_x1) {
    int il = pmb->is + (pmb->ie - pmb->is)*nb.fid+nb.fid;
    int jl = pmb->js, ju = pmb->je, kl = pmb->ks, ku = pmb->ke;
    if (nb.ni.fi1 == 0) ju -= pmb->block_size.nx2/2;
    else          jl += pmb->block_size.nx2/2;
    if (nb.ni.fi2 == 0) ku -= pmb->block_size.nx3/2;
    else          kl += pmb->block_size.nx3/2;

    for (int k=kl; k<=ku; k++) {
      for (int j=jl; j<=ju; j++)
        for (int nn=nl_; nn<=nu_; nn++) {
          x1flux(k,j,il, nn) = buf[p++];
      }
    }
  } else if (nb.fid == BoundaryFace::inner_x2 || nb.fid == BoundaryFace::outer_x2) {
    int jl = pmb->js + (pmb->je - pmb->js)*(nb.fid & 1) + (nb.fid & 1);
    int il = pmb->is, iu = pmb->ie, kl = pmb->ks, ku = pmb->ke;
    if (nb.ni.fi1 == 0) iu -= pmb->block_size.nx1/2;
    else          il += pmb->block_size.nx1/2;
    if (nb.ni.fi2 == 0) ku -= pmb->block_size.nx3/2;
    else          kl += pmb->block_size.nx3/2;

    for (int k=kl; k<=ku; k++) {
      for (int i=il; i<=iu; i++)
        for (int nn=nl_; nn<=nu_; nn++) {
          x2flux(k,jl,i, nn) = buf[p++];
      }
    }
  } else if (nb.fid == BoundaryFace::inner_x3 || nb.fid == BoundaryFace::outer_x3) {
    int kl = pmb->ks + (pmb->ke - pmb->ks)*(nb.fid & 1) + (nb.fid & 1);
    int il = pmb->is, iu = pmb->ie, jl = pmb->js, ju = pmb->je;
    if (nb.ni.fi1 == 0) iu -= pmb->block_size.nx1/2;
    else          il += pmb->block_size.nx1/2;
    if (nb.ni.fi2 == 0) ju -= pmb->block_size.nx2/2;
    else          jl += pmb->block_size.nx2/2;

    for (int j=jl; j<=ju; j++) {
      for (int i=il; i<=iu; i++)
        for (int nn=nl_; nn<=nu_; nn++) {
          x3flux(kl,j,i, nn) = buf[p++];
      }
    }
  }
  return;
}


bool RadBoundaryVariable::ReceiveFluxCorrection() {
  MeshBlock *pmb = pmy_block_;
  bool flag=true;

  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (nb.ni.type != NeighborConnect::face) break;
    if (bd_var_flcor_.flag[nb.bufid] == BoundaryStatus::completed) continue;
    // receive data
    if (bd_var_flcor_.flag[nb.bufid] == BoundaryStatus::waiting) {
      if (nb.snb.rank == Globals::my_rank) {// on the same process
        flag = false;
        continue;
      }
#ifdef MPI_PARALLEL
      else { // NOLINT
        int test;
        // probe MPI communications.  This is a bit of black magic that seems to promote
        // communications to top of stack and gets them to complete more quickly
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &test,
                   MPI_STATUS_IGNORE);
        MPI_Test(&(bd_var_flcor_.req_recv[nb.bufid]), &test, MPI_STATUS_IGNORE);
        if (!static_cast<bool>(test)) {
          flag = false;
          continue;
        }
        bd_var_flcor_.flag[nb.bufid] = BoundaryStatus::arrived;
      }
#endif
    }
    // set data
    if (bd_var_flcor_.flag[nb.bufid] == BoundaryStatus::arrived) {
      if (nb.snb.level==pmb->loc.level)      // from same level
        SetFluxBoundarySameLevel(bd_var_flcor_.recv[nb.bufid],nb);
      else if (nb.snb.level>pmb->loc.level)  // from finer
        SetFluxBoundaryFromFiner(bd_var_flcor_.recv[nb.bufid],nb);
      // else                                   // from coarser
      //   nothing to do

      bd_var_flcor_.flag[nb.bufid] = BoundaryStatus::completed;
    }
  }

  return flag;
}
