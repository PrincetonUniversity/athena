//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file particles_bvals.cpp
//! \brief implements functions for particle communications

// C++ Standard Libraries
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "particles.hpp"

// Local function prototypes
static int CheckSide(int xi, int xi1, int xi2);

//--------------------------------------------------------------------------------------
//! \fn void Particles::AMRCoarseToFine(MeshBlock* pmbc, MeshBlock* pmbf)
//! \brief load particles from a coarse meshblock to a fine meshblock.

void Particles::AMRCoarseToFine(Particles *pparc, Particles *pparf, MeshBlock* pmbf) {
  // Initialization
  const Real x1min = pmbf->block_size.x1min, x1max = pmbf->block_size.x1max;
  const Real x2min = pmbf->block_size.x2min, x2max = pmbf->block_size.x2max;
  const Real x3min = pmbf->block_size.x3min, x3max = pmbf->block_size.x3max;
  const bool active1 = pparc->active1_,
             active2 = pparc->active2_,
             active3 = pparc->active3_;
  const AthenaArray<Real> &xp = pparc->xp, &yp = pparc->yp, &zp = pparc->zp;
  const Coordinates *pcoord = pmbf->pcoord;

  // Loop over particles in the coarse meshblock.
  for (int k = 0; k < pparc->npar; ++k) {
    Real x1, x2, x3;

    if ((!active1 || (active1 && x1min <= x1 && x1 < x1max)) &&
        (!active2 || (active2 && x2min <= x2 && x2 < x2max)) &&
        (!active3 || (active3 && x3min <= x3 && x3 < x3max))) {
      // Load a particle to the fine meshblock.
      int npar = pparf->npar;
      if (npar >= pparf->nparmax) pparf->UpdateCapacity(2 * pparf->nparmax);
      for (int j = 0; j < pparf->nint; ++j)
        pparf->intprop(j,npar) = pparc->intprop(j,k);
      for (int j = 0; j < pparf->nreal; ++j)
        pparf->realprop(j,npar) = pparc->realprop(j,k);
      for (int j = 0; j < pparf->naux; ++j)
        pparf->auxprop(j,npar) = pparc->auxprop(j,k);
      ++pparf->npar;
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::AMRFineToCoarse(MeshBlock* pmbf, MeshBlock* pmbc)
//! \brief load particles from a fine meshblock to a coarse meshblock.

void Particles::AMRFineToCoarse(Particles *pparc, Particles *pparf) {
  // Check the capacity.
  int nparf = pparf->npar, nparc = pparc->npar;
  int npar_new = nparf + nparc;
  if (npar_new > pparc->nparmax) pparc->UpdateCapacity(npar_new);

  // Load the particles.
  for (int j = 0; j < pparf->nint; ++j)
    for (int k = 0; k < nparf; ++k)
      pparc->intprop(j,nparc+k) = pparf->intprop(j,k);
  for (int j = 0; j < pparf->nreal; ++j)
    for (int k = 0; k < nparf; ++k)
      pparc->realprop(j,nparc+k) = pparf->realprop(j,k);
  for (int j = 0; j < pparf->naux; ++j)
    for (int k = 0; k < nparf; ++k)
      pparc->auxprop(j,nparc+k) = pparf->auxprop(j,k);
  pparc->npar = npar_new;
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::ClearBoundary()
//! \brief resets boundary for particle transportation.

void Particles::ClearBoundary() {
  for (int i = 0; i < pbval_->nneighbor; ++i) {
    NeighborBlock& nb = pbval_->neighbor[i];
    bstatus_[nb.bufid] = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
    if (nb.snb.rank != Globals::my_rank) {
      ParticleBuffer& recv = recv_[nb.bufid];
      recv.flagn = recv.flagi = recv.flagr = 0;
      send_[nb.bufid].npar = 0;
    }
#endif
  }

  // clear boundary information for shear
  ClearBoundaryShear();
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::ClearNeighbors()
//! \brief clears links to neighbors.

void Particles::ClearNeighbors() {
  delete neighbor_[1][1][1].pnb;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k) {
        Neighbor *pn = &neighbor_[i][j][k];
        if (pn == NULL) continue;
        while (pn->next != NULL)
          pn = pn->next;
        while (pn->prev != NULL) {
          pn = pn->prev;
          delete pn->next;
          pn->next = NULL;
        }
        pn->pnb = NULL;
        pn->pmb = NULL;
      }
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::LinkNeighbors(MeshBlockTree &tree,
//!         int64_t nrbx1, int64_t nrbx2, int64_t nrbx3, int root_level)
//! \brief fetches neighbor information for later communication.

void Particles::LinkNeighbors(MeshBlockTree &tree,
    int64_t nrbx1, int64_t nrbx2, int64_t nrbx3, int root_level) {
  // Set myself as one of the neighbors.
  Neighbor *pn = &neighbor_[1][1][1];
  pn->pmb = pmy_block;
  pn->pnb = new NeighborBlock;
  pn->pnb->SetNeighbor(Globals::my_rank, pmy_block->loc.level,
      pmy_block->gid, pmy_block->lid, 0, 0, 0, NeighborConnect::none,
      -1, -1, false, false, 0, 0);

  // Save pointer to each neighbor.
  for (int i = 0; i < pbval_->nneighbor; ++i) {
    NeighborBlock& nb = pbval_->neighbor[i];
    SimpleNeighborBlock& snb = nb.snb;
    NeighborIndexes& ni = nb.ni;
    Neighbor *pn = &neighbor_[ni.ox1+1][ni.ox2+1][ni.ox3+1];
    while (pn->next != NULL)
      pn = pn->next;
    if (pn->pnb != NULL) {
      pn->next = new Neighbor;
      pn->next->prev = pn;
      pn = pn->next;
    }
    pn->pnb = &nb;
    if (snb.rank == Globals::my_rank) {
      pn->pmb = pmy_mesh->FindMeshBlock(snb.gid);
    } else {
#ifdef MPI_PARALLEL
      // assign unique tag
      // tag = local id of destination (remaining bits) + bufid (6 bits)
      // + particle container id (3 bits) + npar,intprop,realprop(2 bits)
      send_[nb.bufid].tag = (snb.lid<<11) | (nb.targetid<<5) | (my_ipar_ << 2);
      recv_[nb.bufid].tag = (pmy_block->lid<<11) | (nb.bufid<<5) | (my_ipar_ << 2);
#endif
    }
  }

  // Collect missing directions from fine to coarse level.
  if (pmy_mesh->multilevel) {
    int my_level = pbval_->loc.level;
    for (int l = 0; l < 3; l++) {
      if (!active1_ && l != 1) continue;
      for (int m = 0; m < 3; m++) {
        if (!active2_ && m != 1) continue;
        for (int n = 0; n < 3; n++) {
          if (!active3_ && n != 1) continue;
          Neighbor *pn = &neighbor_[l][m][n];
          if (pn->pnb == NULL) {
            int nblevel = pbval_->nblevel[n][m][l];
            if (0 <= nblevel && nblevel < my_level) {
              int ngid = tree.FindNeighbor(pbval_->loc, l-1, m-1, n-1)->GetGid();
              for (int i = 0; i < pbval_->nneighbor; ++i) {
                NeighborBlock& nb = pbval_->neighbor[i];
                if (nb.snb.gid == ngid) {
                  pn->pnb = &nb;
                  if (nb.snb.rank == Globals::my_rank)
                    pn->pmb = pmy_mesh->FindMeshBlock(ngid);
                  break;
                }
              }
            }
          }
        }
      }
    }
  }

  // Initiate boundary values.
  ClearBoundary();
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::SendToNeighbors()
//! \brief sends particles outside boundary to the buffers of neighboring meshblocks.

void Particles::SendToNeighbors() {
 
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::SendParticleBuffer()
//! \brief send particle buffer

#ifdef MPI_PARALLEL
void Particles::SendParticleBuffer(ParticleBuffer& send, int dst) {
  int npsend = send.npar;
  int sendtag = send.tag;
  MPI_Send(&npsend, 1, MPI_INT, dst, sendtag, my_comm);
  if (npsend > 0) {
    MPI_Request req = MPI_REQUEST_NULL;
    MPI_Isend(send.ibuf, npsend * nint_buf, MPI_INT,
              dst, sendtag + 1, my_comm, &req);
    MPI_Request_free(&req);
    MPI_Isend(send.rbuf, npsend * nreal_buf, MPI_ATHENA_REAL,
              dst, sendtag + 2, my_comm, &req);
    MPI_Request_free(&req);
  }
}
#endif


//--------------------------------------------------------------------------------------
//! \fn void Particles::ReceiveParticleBuffer()
//! \brief recv particle buffer

#ifdef MPI_PARALLEL
void Particles::ReceiveParticleBuffer(int nb_rank, ParticleBuffer& recv,
                                      enum BoundaryStatus& bstatus) {
  // Communicate with neighbor processes.
  if (nb_rank != Globals::my_rank && bstatus == BoundaryStatus::waiting) {
    if (!recv.flagn) {
      // Get the number of incoming particles.
      if (recv.reqn == MPI_REQUEST_NULL)
        MPI_Irecv(&recv.npar, 1, MPI_INT, nb_rank, recv.tag, my_comm, &recv.reqn);
      else
        MPI_Test(&recv.reqn, &recv.flagn, MPI_STATUS_IGNORE);
      if (recv.flagn) {
        if (recv.npar > 0) {
          // Check the buffer size.
          int nprecv = recv.npar;
          if (nprecv > recv.nparmax) {
            recv.npar = 0;
            recv.Reallocate(2 * nprecv - recv.nparmax, nint_buf, nreal_buf);
            recv.npar = nprecv;
          }
        } else {
          // No incoming particles.
          bstatus = BoundaryStatus::completed;
        }
      }
    }
    if (recv.flagn && recv.npar > 0) {
      // Receive data from the neighbor.
      if (!recv.flagi) {
        if (recv.reqi == MPI_REQUEST_NULL)
          MPI_Irecv(recv.ibuf, recv.npar * nint_buf, MPI_INT,
                    nb_rank, recv.tag + 1, my_comm, &recv.reqi);
        else
          MPI_Test(&recv.reqi, &recv.flagi, MPI_STATUS_IGNORE);
      }
      if (!recv.flagr) {
        if (recv.reqr == MPI_REQUEST_NULL)
          MPI_Irecv(recv.rbuf, recv.npar * nreal_buf, MPI_ATHENA_REAL,
                    nb_rank, recv.tag + 2, my_comm, &recv.reqr);
        else
          MPI_Test(&recv.reqr, &recv.flagr, MPI_STATUS_IGNORE);
      }
      if (recv.flagi && recv.flagr)
        bstatus = BoundaryStatus::arrived;
    }
  }
}
#endif

//--------------------------------------------------------------------------------------
//! \fn bool Particles::ReceiveFromNeighbors()
//! \brief receives particles from neighboring meshblocks and returns a flag indicating
//!        if all receives are completed.

bool Particles::ReceiveFromNeighbors() {
  bool flag = true;

  for (int i = 0; i < pbval_->nneighbor; ++i) {
    NeighborBlock& nb = pbval_->neighbor[i];
    enum BoundaryStatus& bstatus = bstatus_[nb.bufid];
    ParticleBuffer& recv = recv_[nb.bufid];
#ifdef MPI_PARALLEL
    ReceiveParticleBuffer(nb.snb.rank, recv, bstatus);
#endif
    switch (bstatus) {
      case BoundaryStatus::completed:
        break;

      case BoundaryStatus::waiting:
        flag = false;
        break;

      case BoundaryStatus::arrived:
        FlushReceiveBuffer(recv);
        bstatus = BoundaryStatus::completed;
        break;
    }
  }

  return flag;
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::ApplyBoundaryConditions(int k, Real &x1, Real &x2, Real &x3)
//! \brief applies boundary conditions to particle k and returns its updated mesh
//!        coordinates (x1,x2,x3).
//! \todo (ccyang):
//! - implement nonperiodic boundary conditions.

void Particles::ApplyBoundaryConditions(int k, Real &x1, Real &x2, Real &x3) {
 

}

//--------------------------------------------------------------------------------------
//! \fn MeshBlock* Particles::FindTargetNeighbor(
//!         int ox1, int ox2, int ox3, int xi1, int xi2, int xi3)
//! \brief finds the neighbor to send a particle to.

struct Neighbor* Particles::FindTargetNeighbor(
    int ox1, int ox2, int ox3, int xi1, int xi2, int xi3) {
  // Find the head of the linked list.
  Neighbor *pn = &neighbor_[ox1+1][ox2+1][ox3+1];
  return pn;
}


//--------------------------------------------------------------------------------------
//! \fn void Particles::FlushReceiveBuffer(ParticleBuffer& recv)
//! \brief adds particles from the receive buffer.

void Particles::FlushReceiveBuffer(ParticleBuffer& recv) {

}


//--------------------------------------------------------------------------------------
//! \fn void Particles::LoadParticleBuffer(ParticleBuffer& recv)
//! \brief load the k-th particle to the particle buffer.

void Particles::LoadParticleBuffer(ParticleBuffer *ppb, int k) {
  // Check the buffer size.
  if (ppb->npar >= ppb->nparmax)
    ppb->Reallocate((ppb->nparmax > 0) ? 2 * ppb->nparmax : 1, nint_buf, nreal_buf);

  // Copy the properties of the particle to the buffer.
  int *pi = ppb->ibuf + nint_buf * ppb->npar;
  for (int j = 0; j < nint; ++j)
    *pi++ = intprop(j,k);
  Real *pr = ppb->rbuf + nreal_buf * ppb->npar;
  for (int j = 0; j < nreal; ++j)
    *pr++ = realprop(j,k);
  for (int j = 0; j < naux; ++j)
    *pr++ = auxprop(j,k);
  ++ppb->npar;
}

//--------------------------------------------------------------------------------------
//! \fn LogicalLocation Particles::FindTargetGidAlongX2(Real x2)
//! \brief searching meshblock along the x2 direction due to shear
//!        (uniform mesh is assumed)

int Particles::FindTargetGidAlongX2(Real x2) {

  // failed to find corresponding meshblock
  return -1;
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::SendParticlesShear()
//! \brief send particles outside meshblock due to shear

void Particles::SendParticlesShear() {
 
}


//--------------------------------------------------------------------------------------
//! \fn bool Particles::ReceiveFromNeighborsShear()
//! \brief receives particles shifted out due to shear and returns a flag indicating
//!        if all receives are completed.

bool Particles::ReceiveFromNeighborsShear() {

  return true;
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::ApplyBoundaryConditionsShear(int k, Real &x1, Real &x2, Real &x3)
//! \brief apply shift due to shear boudnary crossings

void Particles::ApplyBoundaryConditionsShear(int k, Real &x1, Real &x2, Real &x3) {
 
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::StartReceivingParticlesShear()
//! \brief initialize boundary status and MPI tags
//!
//! Since the meshblocks to be communicated due to shear change in time,
//! communication tags, boundary status should be set after the ComputeShear call,
//! which is done in StartupTaskList. This is necessary function call at every substeps.

void Particles::StartReceivingParticlesShear() {

}
//--------------------------------------------------------------------------------------
//! \fn void Particles::ClearBoundaryShear()
//! \brief resets boundary for particle transportation due to shear.

void Particles::ClearBoundaryShear() {

}

//--------------------------------------------------------------------------------------
//! \fn int CheckSide(int xi, nx, int xi1, int xi2)
//! \brief returns -1 if xi < xi1, +1 if xi > xi2, or 0 otherwise.

inline int CheckSide(int xi, int xi1, int xi2) {
  if (xi < xi1) return -1;
  if (xi > xi2) return +1;
  return 0;
}
