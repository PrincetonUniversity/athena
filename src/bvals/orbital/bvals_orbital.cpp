//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_orbital.cpp
//! \brief functions for orbital communications

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>    // memcpy()
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../field/field.hpp"
#include "../../globals.hpp"
#include "../../hydro/hydro.hpp"
#include "../../mesh/mesh.hpp"
#include "../../orbital_advection/orbital_advection.hpp"
#include "../../scalars/scalars.hpp"
#include "../../utils/buffer_utils.hpp"
#include "../bvals.hpp"
#include "bvals_orbital.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// constructor

OrbitalBoundaryCommunication::OrbitalBoundaryCommunication(
    OrbitalAdvection *porb)
    : xgh(porb->xgh), pmy_block_(porb->pmb_), pmy_mesh_(porb->pm_),
      pbval_(porb->pbval_), pmy_orbital_(porb) {
  for (int upper=0; upper<2; upper++) {
    InitBoundaryData(orbital_bd_cc_[upper], BoundaryQuantity::orbital_cc);
  }
#ifdef MPI_PARALLEL
  orbital_advection_cc_phys_id_ = pbval_->AdvanceCounterPhysID(max_phys_id);
#endif
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int upper=0; upper<2; upper++) {
      InitBoundaryData(orbital_bd_fc_[upper], BoundaryQuantity::orbital_fc);
    }
#ifdef MPI_PARALLEL
    orbital_advection_fc_phys_id_ = orbital_advection_cc_phys_id_+1;
#endif
  }
  if(pmy_orbital_->orbital_refinement) {
    for (int n=0; n<2; n++) {
      size_cc_send[n] = new int[6];
      size_cc_recv[n] = new int[6];
      if (MAGNETIC_FIELDS_ENABLED) {
        size_fc_send[n] = new int[6];
        size_fc_recv[n] = new int[6];
      }
    }
  } else {
    for (int n=0; n<2; n++) {
      size_cc_send[n] = new int[1];
      size_cc_recv[n] = new int[1];
      if (MAGNETIC_FIELDS_ENABLED) {
        size_fc_send[n] = new int[1];
        size_fc_recv[n] = new int[1];
      }
    }
  }
}

// destructor
OrbitalBoundaryCommunication::~OrbitalBoundaryCommunication() {
  for (int upper=0; upper<2; upper++) {
    DestroyBoundaryData(orbital_bd_cc_[upper]);
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int upper=0; upper<2; upper++) {
      DestroyBoundaryData(orbital_bd_fc_[upper]);
    }
  }

  for (int n=0; n<2; n++) {
    delete[] size_cc_send[n];
    delete[] size_cc_recv[n];
    if (MAGNETIC_FIELDS_ENABLED) {
      delete[] size_fc_send[n];
      delete[] size_fc_recv[n];
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void OrbitalBoundaryCommunication::InitBoundaryData(
//!            OrbitalBoundaryData &bd, BoundaryQuantity type)
//! \brief allocate memory for orbital communicaiton
void OrbitalBoundaryCommunication::InitBoundaryData(
     OrbitalBoundaryData &bd, BoundaryQuantity type) {
  MeshBlock *pmb = pmy_block_;
  OrbitalAdvection *porb = pmy_orbital_;
  int nx1 = pmb->block_size.nx1;
  int nx2 = pmb->block_size.nx2;
  int nx3 = pmb->block_size.nx3;

  if (porb->orbital_refinement) {
    int ssize(0), lsize(0);
    if (porb->orbital_direction == 1) {
      if (nx3>1) { //3D
        bd.nbmax = 4;
        if (type == BoundaryQuantity::orbital_cc) {
          ssize = (NHYDRO+NSCALARS)*(nx1/2+2)*(nx2/2+2)*(nx3/2+2);
          lsize = (NHYDRO+NSCALARS)*nx1*nx2*nx3;
        } else if (type == BoundaryQuantity::orbital_fc) {
          ssize = (nx1/2+1)*(nx2/2+2)*(nx3/2+2)
                  +(nx1/2+2)*(nx2/2+1)*(nx3/2+2)
                  +(nx1/2+2)*(nx2/2+2)*(nx3/2+1);
          lsize = (2*nx1*nx3+nx1+nx3)*nx2;
        } else {
          std::stringstream msg;
          msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
              << "::InitBoundaryData" << std::endl
              << "Invalid boundary type is specified." << std::endl;
          ATHENA_ERROR(msg);
        }
      } else { //2D
        bd.nbmax = 2;
        if (type == BoundaryQuantity::orbital_cc) {
          ssize = (NHYDRO+NSCALARS)*(nx1/2+2)*(nx2/2+2);
          lsize = (NHYDRO+NSCALARS)*nx1*nx2;
        } else if (type == BoundaryQuantity::orbital_fc) {
          ssize = (nx1/2+1)*(nx2/2+2)
                  +(nx1/2+2)*(nx2/2+1)
                  +(nx1/2+2)*(nx2/2+2);
          lsize = (3*nx1+1)*nx2;
        } else {
          std::stringstream msg;
          msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
              << "::InitBoundaryData" << std::endl
              << "Invalid boundary type is specified." << std::endl;
          ATHENA_ERROR(msg);
        }
      }
    } else if (porb->orbital_direction == 2) {
      bd.nbmax = 4;
      if (type == BoundaryQuantity::orbital_cc) {
        ssize = (NHYDRO+NSCALARS)*(nx1/2+2)*(nx2/2+2)*(nx3/2+2);
        lsize = (NHYDRO+NSCALARS)*nx1*nx2*nx3;
      } else if (type == BoundaryQuantity::orbital_fc) {
        ssize = (nx1/2+1)*(nx2/2+2)*(nx3/2+2)
                  +(nx1/2+2)*(nx2/2+1)*(nx3/2+2)
                  +(nx1/2+2)*(nx2/2+2)*(nx3/2+1);
        lsize = (2*nx1*nx2+nx1+nx2)*nx3;
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::InitBoundaryData" << std::endl
            << "Invalid boundary type is specified." << std::endl;
        ATHENA_ERROR(msg);
      }
    }
    for (int n=0; n<bd.nbmax; n++) {
      if (n==0) {
        bd.send[n]  = new Real[lsize];
        bd.recv[n]  = new Real[lsize];
      } else {
        bd.send[n]  = new Real[ssize];
        bd.recv[n]  = new Real[ssize];
      }
      bd.flag[n]  = BoundaryStatus::waiting;
      bd.sflag[n] = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
      bd.req_send[n] = MPI_REQUEST_NULL;
      bd.req_recv[n] = MPI_REQUEST_NULL;
#endif
    }
  } else {
    bd.nbmax = 1;
    int size(0);
    if (porb->orbital_direction == 1) {
      if (type == BoundaryQuantity::orbital_cc) {
        size = (NHYDRO+NSCALARS)*nx1*nx2*nx3;
      } else if (type == BoundaryQuantity::orbital_fc) {
        size = (2*nx1*nx3+nx1+nx3)*nx2;
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::InitBoundaryData" << std::endl
            << "Invalid boundary type is specified." << std::endl;
        ATHENA_ERROR(msg);
      }
    } else if (porb->orbital_direction == 2) {
      if (type == BoundaryQuantity::orbital_cc) {
        size = (NHYDRO+NSCALARS)*nx1*nx2*nx3;
      } else if (type == BoundaryQuantity::orbital_fc) {
        size = (2*nx1*nx2+nx1+nx2)*nx3;
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::InitBoundaryData" << std::endl
            << "Invalid boundary type is specified." << std::endl;
        ATHENA_ERROR(msg);
      }
    }
    for (int n=0; n<bd.nbmax; n++) {
      bd.send[n]  = new Real[size];
      bd.recv[n]  = new Real[size];
      bd.flag[n]  = BoundaryStatus::waiting;
      bd.sflag[n] = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
      bd.req_send[n] = MPI_REQUEST_NULL;
      bd.req_recv[n] = MPI_REQUEST_NULL;
#endif
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OrbitalBoundaryCommunication::DestroyBoundaryData(OrbitalBoundaryData &bd)
//! \brief destroy BoundaryData structure

void OrbitalBoundaryCommunication::DestroyBoundaryData(OrbitalBoundaryData &bd) {
  for (int n=0; n<bd.nbmax; n++) {
    delete [] bd.send[n];
    delete [] bd.recv[n];
#ifdef MPI_PARALLEL
    if (bd.req_send[n] != MPI_REQUEST_NULL)
      MPI_Request_free(&bd.req_send[n]);
    if (bd.req_recv[n] != MPI_REQUEST_NULL)
      MPI_Request_free(&bd.req_recv[n]);
#endif
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OrbitalBoundaryCommunication::SetupPersistentMPI()
//! \brief set up for MPI in Mesh::Initiate()
void OrbitalBoundaryCommunication::SetupPersistentMPI() {
  OrbitalAdvection *porb = pmy_orbital_;
  // initialize
  for (int upper=0; upper<2; upper++) {
    for (int n=0; n<orbital_bd_cc_[upper].nbmax; n++) {
      SimpleNeighborBlock &snb1 = orbital_send_neighbor_[upper][n];
      SimpleNeighborBlock &snb2 = orbital_recv_neighbor_[upper][n];
      snb1.gid  = -1;  snb1.lid   = -1;
      snb1.rank = -1;  snb1.level = -1;
      snb2.gid  = -1;  snb2.lid   = -1;
      snb2.rank = -1;  snb2.level = -1;
    }
  }
  if (porb->orbital_direction == 1) {
    for(int n=0; n<pbval_->nneighbor; n++) {
      NeighborBlock& nb = pbval_->neighbor[n];
      if(nb.ni.ox1==0 && nb.ni.ox3==0) {
        if(nb.ni.ox2==-1) {
          SimpleNeighborBlock &snb1 = orbital_send_neighbor_[1][nb.ni.fi1+nb.ni.fi2*2];
          SimpleNeighborBlock &snb2 = orbital_recv_neighbor_[0][nb.ni.fi1+nb.ni.fi2*2];
          snb1.gid = nb.snb.gid;    snb1.lid = nb.snb.lid;
          snb1.rank = nb.snb.rank;  snb1.level = nb.snb.level;
          snb2.gid = nb.snb.gid;    snb2.lid = nb.snb.lid;
          snb2.rank = nb.snb.rank;  snb2.level = nb.snb.level;
        } else if(nb.ni.ox2== 1) {
          SimpleNeighborBlock &snb1 = orbital_send_neighbor_[0][nb.ni.fi1+nb.ni.fi2*2];
          SimpleNeighborBlock &snb2 = orbital_recv_neighbor_[1][nb.ni.fi1+nb.ni.fi2*2];
          snb1.gid = nb.snb.gid;    snb1.lid = nb.snb.lid;
          snb1.rank = nb.snb.rank;  snb1.level = nb.snb.level;
          snb2.gid = nb.snb.gid;    snb2.lid = nb.snb.lid;
          snb2.rank = nb.snb.rank;  snb2.level = nb.snb.level;
        }
      }
    }
  } else if (porb->orbital_direction == 2) {
    for(int n=0; n<pbval_->nneighbor; n++) {
      NeighborBlock& nb = pbval_->neighbor[n];
      if(nb.ni.ox1==0 && nb.ni.ox2==0) {
        if(nb.ni.ox3==-1) {
          SimpleNeighborBlock &snb1 = orbital_send_neighbor_[1][nb.ni.fi1+nb.ni.fi2*2];
          SimpleNeighborBlock &snb2 = orbital_recv_neighbor_[0][nb.ni.fi1+nb.ni.fi2*2];
          snb1.gid = nb.snb.gid;    snb1.lid = nb.snb.lid;
          snb1.rank = nb.snb.rank;  snb1.level = nb.snb.level;
          snb2.gid = nb.snb.gid;    snb2.lid = nb.snb.lid;
          snb2.rank = nb.snb.rank;  snb2.level = nb.snb.level;
        } else if(nb.ni.ox3== 1) {
          SimpleNeighborBlock &snb1 = orbital_send_neighbor_[0][nb.ni.fi1+nb.ni.fi2*2];
          SimpleNeighborBlock &snb2 = orbital_recv_neighbor_[1][nb.ni.fi1+nb.ni.fi2*2];
          snb1.gid = nb.snb.gid;    snb1.lid = nb.snb.lid;
          snb1.rank = nb.snb.rank;  snb1.level = nb.snb.level;
          snb2.gid = nb.snb.gid;    snb2.lid = nb.snb.lid;
          snb2.rank = nb.snb.rank;  snb2.level = nb.snb.level;
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OrbitalBoundaryCommunication::ComputeOrbit(const Real dt)
//! \brief Calculate number of cells for orbital communication
void OrbitalBoundaryCommunication::ComputeOrbit(const Real dt) {
  MeshBlock *pmb = pmy_block_;
  OrbitalAdvection *porb = pmy_orbital_;

  // initialize counts
  for (int upper=0; upper<2; upper++) {
    for (int n=0; n<orbital_bd_cc_[upper].nbmax; n++) {
      orbital_send_cc_count_[upper][n]   = 0;
      orbital_recv_cc_count_[upper][n]   = 0;
    }
    if (MAGNETIC_FIELDS_ENABLED) {
      for (int n=0; n<orbital_bd_fc_[upper].nbmax; n++) {
        orbital_send_fc_count_[upper][n]   = 0;
        orbital_recv_fc_count_[upper][n]   = 0;
      }
    }
  }

  //calculate offsets and count cell #
  if (porb->orbital_refinement) { // orbital refinement
    int mylevel = pmb->loc.level;
    porb->SetOrbitalEdgeCC(dt, size_cc_send, size_cc_recv);
    int coef = pmb->block_size.nx1/2+2;
    int size(0);
    if (porb->orbital_direction == 1) {
      if (pmb->block_size.nx3>1) { // 3D
        coef *= pmb->block_size.nx3/2+2;
      }
    } else if (porb->orbital_direction == 2) {
      coef *= pmb->block_size.nx2/2+2;
    }
    // upper=0: to right, upper=1: to left
    for (int upper=0; upper<2; upper++) {
      // send
      if (orbital_send_neighbor_[upper][0].level>mylevel) { // to finer
        for (int n=0; n<orbital_bd_cc_[upper].nbmax; n++) {
          size = size_cc_send[upper][2+n]*coef;
          if (size>0) {
            orbital_send_cc_count_[upper][n] = size;
          }
        }
      } else if (orbital_send_neighbor_[upper][0].level<mylevel) { // to coarser
        size = size_cc_send[upper][1];
        if (size>0) {
          orbital_send_cc_count_[upper][0] = size;
        }
      } else {
        size = size_cc_send[upper][0]; // to same level
        if (size>0) {
          orbital_send_cc_count_[upper][0] = size;
        }
      }
      //recv
      if (orbital_recv_neighbor_[upper][0].level>mylevel) {
        for (int n=0; n<orbital_bd_cc_[upper].nbmax; n++) { // from finer
          size = size_cc_recv[upper][2+n];
          if (size>0) {
            orbital_recv_cc_count_[upper][n] = size;
          }
        }
      } else if (orbital_recv_neighbor_[upper][0].level<mylevel) { // from coarser
        size = size_cc_recv[upper][1]*coef;
        if (size>0) {
          orbital_recv_cc_count_[upper][0] = size;
        }
      } else {
        size = size_cc_recv[upper][0];
        if (size>0) {
          orbital_recv_cc_count_[upper][0] = size;
        }
      }
    }

    if (MAGNETIC_FIELDS_ENABLED) {
      porb->SetOrbitalEdgeFC(dt, size_fc_send, size_fc_recv);
      int nco1=pmb->block_size.nx1/2;
      int nco2=(porb->orbital_direction==1)?
                pmb->block_size.nx3/2 : pmb->block_size.nx2/2;
      // upper=0: to right, upper=1: to left
      for (int upper=0; upper<2; upper++) {
        // send
        if (orbital_send_neighbor_[upper][0].level>mylevel) { // to finer
          for (int n=0; n<orbital_bd_fc_[upper].nbmax; n++) {
            int norb = size_fc_send[upper][2+n];
            if (pmb->block_size.nx3>1) { // 3D
              size = (nco1+1)*(nco2+2)*(norb+2)
                     +(nco1+2)*(nco2+1)*(norb+2)
                     +(nco1+2)*(nco2+2)*(norb+1);
            } else { // 2D
              size = (nco1+1)*(norb+2)
                     +(nco1+2)*(norb+1)
                     +(nco1+2)*(norb+2);
            }
            if (size>0) {
              orbital_send_fc_count_[upper][n] = size;
            }
          }
        } else if (orbital_send_neighbor_[upper][0].level<mylevel) { //to coarser
          size = size_fc_send[upper][1];
          if(size>0) {
            orbital_send_fc_count_[upper][0] = size;
          }
        } else { // to same
          size = size_fc_send[upper][0];
          if(size>0) {
            orbital_send_fc_count_[upper][0] = size;
          }
        }
        // recv
        if (orbital_recv_neighbor_[upper][0].level>mylevel) { // from finer
          for (int n=0; n<orbital_bd_fc_[upper].nbmax; n++) {
            size = size_fc_recv[upper][2+n];
            if (size>0) {
              orbital_recv_fc_count_[upper][n] = size;
            }
          }
        } else if (orbital_recv_neighbor_[upper][0].level<mylevel) { // from coarser
          int norb = size_fc_recv[upper][1];
          if (pmb->block_size.nx3>1) { // 3D
            size = (nco1+1)*(nco2+2)*(norb+2)
                   +(nco1+2)*(nco2+1)*(norb+2)
                   +(nco1+2)*(nco2+2)*(norb+1);
          } else { // 2D
            size = (nco1+1)*(norb+2)
                   +(nco1+2)*(norb+1)
                   +(nco1+2)*(norb+2);
          }
          if(size>0) {
            orbital_recv_fc_count_[upper][0] = size;
          }
        } else { // from same
          size =  size_fc_recv[upper][0];
          if(size>0) {
            orbital_recv_fc_count_[upper][0] = size;
          }
        }
      }
    }
  } else { // no orbital refinement
    porb->SetOrbitalEdgeCC(dt, size_cc_send, size_cc_recv);
    int size;
    for (int upper=0; upper<2; upper++) {
      size = size_cc_send[upper][0];
      if (size>0) {
        orbital_send_cc_count_[upper][0] = size;
      }
      size = size_cc_recv[upper][0];
      if (size>0) {
        orbital_recv_cc_count_[upper][0] = size;
      }
    }
    if (MAGNETIC_FIELDS_ENABLED) {
      porb->SetOrbitalEdgeFC(dt, size_fc_send, size_fc_recv);
      for (int upper=0; upper<2; upper++) {
        size = size_fc_send[upper][0];
        if (size>0) {
          orbital_send_fc_count_[upper][0] = size;
        }
        size = size_fc_recv[upper][0];
        if (size>0) {
          orbital_recv_fc_count_[upper][0] = size;
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OrbitalBoundaryCommunication::StartReceiving(BoundaryCommSubset phase)
//! \brief StartReceiving function for the orbital communication
void OrbitalBoundaryCommunication::StartReceiving(BoundaryCommSubset phase) {
#ifdef MPI_PARALLEL
  int tag, size;
  int tag_offset[2]{0,4};
#endif
  for (int upper=0; upper<2; upper++) {
    for (int n=0; n<orbital_bd_cc_[upper].nbmax; n++) {
      if (orbital_send_cc_count_[upper][n]>0) {
        orbital_bd_cc_[upper].sflag[n] = BoundaryStatus::waiting;
      } else {
        orbital_bd_cc_[upper].sflag[n] = BoundaryStatus::completed;
      }
      if (orbital_recv_cc_count_[upper][n]>0) {
        orbital_bd_cc_[upper].flag[n]  = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
        int target_rank = orbital_recv_neighbor_[upper][n].rank;
        if ((target_rank != Globals::my_rank) && (target_rank != -1)) {
          size = (NHYDRO+NSCALARS)*orbital_recv_cc_count_[upper][n];
          tag  = pbval_->CreateBvalsMPITag(pmy_block_->lid, n+tag_offset[upper],
                                           orbital_advection_cc_phys_id_);
          MPI_Irecv(orbital_bd_cc_[upper].recv[n], size, MPI_ATHENA_REAL,
                    target_rank, tag, MPI_COMM_WORLD, &orbital_bd_cc_[upper].req_recv[n]);
        }
#endif
      } else {
        orbital_bd_cc_[upper].flag[n]  = BoundaryStatus::completed;
      }
    }
    if (MAGNETIC_FIELDS_ENABLED) {
      for (int n=0; n<orbital_bd_fc_[upper].nbmax; n++) {
        if (orbital_send_fc_count_[upper][n]>0) {
          orbital_bd_fc_[upper].sflag[n] = BoundaryStatus::waiting;
        } else {
          orbital_bd_fc_[upper].sflag[n] = BoundaryStatus::completed;
        }
        if (orbital_recv_fc_count_[upper][n]>0) {
          orbital_bd_fc_[upper].flag[n]  = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
          int target_rank = orbital_recv_neighbor_[upper][n].rank;
          if ((target_rank != Globals::my_rank) && (target_rank != -1)) {
            size = orbital_recv_fc_count_[upper][n];
            tag  = pbval_->CreateBvalsMPITag(pmy_block_->lid, n+tag_offset[upper],
                                             orbital_advection_fc_phys_id_);
            MPI_Irecv(orbital_bd_fc_[upper].recv[n], size, MPI_ATHENA_REAL,
                      target_rank, tag, MPI_COMM_WORLD,
                      &orbital_bd_fc_[upper].req_recv[n]);
          }
#endif
        } else {
          orbital_bd_fc_[upper].flag[n]  = BoundaryStatus::completed;
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OrbitalBoundaryCommunication::ClearBoundary(BoundaryCommSubset phase)
//! \brief ClearBoundary function for the orbital communication
void OrbitalBoundaryCommunication::ClearBoundary(BoundaryCommSubset phase) {
  for (int upper=0; upper<2; upper++) {
    for (int n=0; n<orbital_bd_cc_[upper].nbmax; n++) {
      if (orbital_send_neighbor_[upper][n].rank == -1) continue;
      orbital_bd_cc_[upper].flag[n] = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
      if (orbital_send_neighbor_[upper][n].rank != Globals::my_rank) {
        MPI_Wait(&orbital_bd_cc_[upper].req_send[n], MPI_STATUS_IGNORE);
      }
#endif
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int upper=0; upper<2; upper++) {
      for (int n=0; n<orbital_bd_fc_[upper].nbmax; n++) {
        if (orbital_send_neighbor_[upper][n].rank == -1) continue;
        orbital_bd_fc_[upper].flag[n] = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
        if (orbital_send_neighbor_[upper][n].rank != Globals::my_rank) {
          MPI_Wait(&orbital_bd_fc_[upper].req_send[n], MPI_STATUS_IGNORE);
        }
#endif
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OrbitalBoundaryCommunication::SendBoundaryBuffersCC()
//! \brief load and send hydro variables and passive scalars
void OrbitalBoundaryCommunication::SendBoundaryBuffersCC() {
  MeshBlock *pmb = pmy_block_;
  OrbitalAdvection *porb = pmy_orbital_;
  int tag;
  int mylevel = pmy_block_->loc.level;
  int offset[2]{0,4};
  for (int upper=0; upper<2; upper++) {
    for (int n=0; n<orbital_bd_cc_[upper].nbmax; n++) {
      if (orbital_bd_cc_[upper].sflag[n] == BoundaryStatus::completed) continue;
      SimpleNeighborBlock& snb= orbital_send_neighbor_[upper][n];
      int p=0;
      if (snb.level == mylevel) { // to same
        LoadHydroBufferSameLevel(orbital_bd_cc_[upper].send[n], p, n+offset[upper]);
      } else if (snb.level < mylevel) { // to coarser
        LoadHydroBufferToCoarser(orbital_bd_cc_[upper].send[n], p, n+offset[upper]);
      } else { // to finer
        LoadHydroBufferToFiner(orbital_bd_cc_[upper].send[n], p, n+offset[upper]);
      }
      if (NSCALARS>0) {
        if (snb.level == mylevel) { // to same
          LoadScalarBufferSameLevel(orbital_bd_cc_[upper].send[n], p, n+offset[upper]);
        } else if (snb.level < mylevel) { // to coarser
          LoadScalarBufferToCoarser(orbital_bd_cc_[upper].send[n], p, n+offset[upper]);
        } else { // to finer
          LoadScalarBufferToFiner(orbital_bd_cc_[upper].send[n], p, n+offset[upper]);
        }
      }
      if (p != (NHYDRO+NSCALARS)*orbital_send_cc_count_[upper][n]) {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::SendBoundaryBuffersCC" << std::endl
            << "Send buffer size is incorrect." << std::endl;
        ATHENA_ERROR(msg);
      }
      if (snb.rank == Globals::my_rank) { //on the same process
        MeshBlock *tmb = pmy_mesh_->FindMeshBlock(snb.gid);
        OrbitalBoundaryData &obd = tmb->porb->orb_bc->orbital_bd_cc_[upper];
        if (snb.level == mylevel) { // to same level
          std::memcpy(obd.recv[n], orbital_bd_cc_[upper].send[n], p*sizeof(Real));
          obd.flag[n] = BoundaryStatus::arrived;
        } else if (snb.level < mylevel) { // to coarser
          int n1 = static_cast<int>(pmb->loc.lx1%2);
          if (porb->orbital_direction == 1) {
            int n3 = static_cast<int>(pmb->loc.lx3%2);
            std::memcpy(obd.recv[n1+n3*2], orbital_bd_cc_[upper].send[n], p*sizeof(Real));
            obd.flag[n1+n3*2] = BoundaryStatus::arrived;
          } else {
            int n2 = static_cast<int>(pmb->loc.lx2%2);
            std::memcpy(obd.recv[n1+n2*2], orbital_bd_cc_[upper].send[n], p*sizeof(Real));
            obd.flag[n1+n2*2] = BoundaryStatus::arrived;
          }
        } else { //to finer
          std::memcpy(obd.recv[0], orbital_bd_cc_[upper].send[n], p*sizeof(Real));
          obd.flag[0] = BoundaryStatus::arrived;
        }
      } else {
#ifdef MPI_PARALLEL
        if (snb.level == mylevel) { // to same level
          tag = pbval_->CreateBvalsMPITag(snb.lid, n+offset[upper],
                                          orbital_advection_cc_phys_id_);
        } else if (snb.level < mylevel) { // to coarser
          int n1 = static_cast<int>(pmb->loc.lx1%2);
          if (porb->orbital_direction == 1) {
            int n3 = static_cast<int>(pmb->loc.lx3%2);
            tag = pbval_->CreateBvalsMPITag(snb.lid, n1+n3*2+offset[upper],
                                            orbital_advection_cc_phys_id_);
          } else {
            int n2 = static_cast<int>(pmb->loc.lx2%2);
            tag = pbval_->CreateBvalsMPITag(snb.lid, n1+n2*2+offset[upper],
                                            orbital_advection_cc_phys_id_);
          }
        } else { //to finer
          tag = pbval_->CreateBvalsMPITag(snb.lid, offset[upper],
                                          orbital_advection_cc_phys_id_);
        }
        MPI_Isend(orbital_bd_cc_[upper].send[n], p, MPI_ATHENA_REAL, snb.rank,
                  tag, MPI_COMM_WORLD, &orbital_bd_cc_[upper].req_send[n]);
#endif //MPI
      }
      orbital_bd_cc_[upper].sflag[n] = BoundaryStatus::completed;
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn bool OrbitalBoundaryCommunication::ReceiveBoundaryBuffersCC()
//! \brief receive and set hydro variables and passive scalars
bool OrbitalBoundaryCommunication::ReceiveBoundaryBuffersCC() {
  bool flag[2]{true, true};
  int mylevel = pmy_block_->loc.level;
  int offset[2]{0,4};
  for (int upper=0; upper<2; upper++) {
    for (int n=0; n<orbital_bd_cc_[upper].nbmax; n++) {
      if (orbital_bd_cc_[upper].flag[n] == BoundaryStatus::completed) continue;
      SimpleNeighborBlock& snb= orbital_recv_neighbor_[upper][n];
      if (orbital_bd_cc_[upper].flag[n] == BoundaryStatus::waiting) {
        // on the same process
        if (snb.rank == Globals::my_rank) {
          flag[upper] = false;
          continue;
        } else { // MPI
#ifdef MPI_PARALLEL
          int test;
          MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &test,
                     MPI_STATUS_IGNORE);
          MPI_Test(&orbital_bd_cc_[upper].req_recv[n], &test, MPI_STATUS_IGNORE);
          if (!static_cast<bool>(test)) {
            flag[upper] = false;
            continue;
          }
#endif
          orbital_bd_cc_[upper].flag[n] = BoundaryStatus::arrived;
        }
      }
      // set var
      int p=0;
      if (snb.level == mylevel) { // from same level
        SetHydroBufferSameLevel(orbital_bd_cc_[upper].recv[n], p, n+offset[upper]);
      } else if (snb.level < mylevel) { // from coarser
        SetHydroBufferFromCoarser(orbital_bd_cc_[upper].recv[n], p, n+offset[upper]);
      } else { // from finer
        SetHydroBufferFromFiner(orbital_bd_cc_[upper].recv[n], p, n+offset[upper]);
      }
      if (NSCALARS>0) {
        if (snb.level == mylevel) { // from same level
          SetScalarBufferSameLevel(orbital_bd_cc_[upper].recv[n], p, n+offset[upper]);
        } else if (snb.level < mylevel) { // from coarser
          SetScalarBufferFromCoarser(orbital_bd_cc_[upper].recv[n], p, n+offset[upper]);
        } else { // from finer
          SetScalarBufferFromFiner(orbital_bd_cc_[upper].recv[n], p, n+offset[upper]);
        }
      }
      if (p != (NHYDRO+NSCALARS)*orbital_recv_cc_count_[upper][n]) {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::ReceiveBoundaryBuffersCC" << std::endl
            << "Recieve buffer size is incorrect." << std::endl;
        ATHENA_ERROR(msg);
      }
      orbital_bd_cc_[upper].flag[n] = BoundaryStatus::completed;
    } // loop over recv[0] to recv[nbmax-1]
  } // loop over comm direction
  return (flag[0] && flag[1]);
}

//----------------------------------------------------------------------------------------
//! \fn void OrbitalBoundaryCommunication::SendBoundaryBuffersFC()
//! \brief load and send magnetic fields
void OrbitalBoundaryCommunication::SendBoundaryBuffersFC() {
  MeshBlock *pmb = pmy_block_;
  OrbitalAdvection *porb = pmy_orbital_;
  int tag;
  int mylevel = pmb->loc.level;
  int offset[2]{0,4};
  for (int upper=0; upper<2; upper++) {
    for (int n=0; n<orbital_bd_fc_[upper].nbmax; n++) {
      if (orbital_bd_fc_[upper].sflag[n] == BoundaryStatus::completed) continue;
      SimpleNeighborBlock& snb= orbital_send_neighbor_[upper][n];
      int p=0;
      if (snb.level == mylevel) { // to same level
        LoadFieldBufferSameLevel(orbital_bd_fc_[upper].send[n], p, n+offset[upper]);
      } else if (snb.level < mylevel) { // to coarser
        LoadFieldBufferToCoarser(orbital_bd_fc_[upper].send[n], p, n+offset[upper]);
      } else { // to finer
        LoadFieldBufferToFiner(orbital_bd_fc_[upper].send[n], p, n+offset[upper]);
      }
      if (p != orbital_send_fc_count_[upper][n]) {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::SendBoundaryBuffersFC" << std::endl
            << "Send buffer size is incorrect." << std::endl;
        ATHENA_ERROR(msg);
      }
      if (snb.rank == Globals::my_rank) { //on the same process
        MeshBlock *tmb = pmy_mesh_->FindMeshBlock(snb.gid);
        OrbitalBoundaryData &obd = tmb->porb->orb_bc->orbital_bd_fc_[upper];
        if (snb.level == mylevel) { // to same level
          std::memcpy(obd.recv[n], orbital_bd_fc_[upper].send[n], p*sizeof(Real));
          obd.flag[n] = BoundaryStatus::arrived;
        } else if (snb.level < mylevel) { // to coarser
          int n1 = static_cast<int>(pmb->loc.lx1%2);
          if (porb->orbital_direction == 1) {
            int n3 = static_cast<int>(pmb->loc.lx3%2);
            std::memcpy(obd.recv[n1+n3*2], orbital_bd_fc_[upper].send[n], p*sizeof(Real));
            obd.flag[n1+n3*2] = BoundaryStatus::arrived;
          } else {
            int n2 = static_cast<int>(pmb->loc.lx2%2);
            std::memcpy(obd.recv[n1+n2*2], orbital_bd_fc_[upper].send[n], p*sizeof(Real));
            obd.flag[n1+n2*2] = BoundaryStatus::arrived;
          }
        } else { // to finer
          std::memcpy(obd.recv[0], orbital_bd_fc_[upper].send[n], p*sizeof(Real));
          obd.flag[0] = BoundaryStatus::arrived;
        }
      } else {
#ifdef MPI_PARALLEL
        if (snb.level == mylevel) { // to same level
          tag = pbval_->CreateBvalsMPITag(snb.lid, n+offset[upper],
                                          orbital_advection_fc_phys_id_);
        } else if (snb.level < mylevel) { // to coarser
          int n1 = static_cast<int>(pmb->loc.lx1%2);
          if (porb->orbital_direction == 1) {
            int n3 = static_cast<int>(pmb->loc.lx3%2);
            tag = pbval_->CreateBvalsMPITag(snb.lid, n1+n3*2+offset[upper],
                                            orbital_advection_fc_phys_id_);
          } else {
            int n2 = static_cast<int>(pmb->loc.lx2%2);
            tag = pbval_->CreateBvalsMPITag(snb.lid, n1+n2*2+offset[upper],
                                            orbital_advection_fc_phys_id_);
          }
        } else { // to finer
          tag = pbval_->CreateBvalsMPITag(snb.lid, offset[upper],
                                          orbital_advection_fc_phys_id_);
        }
        MPI_Isend(orbital_bd_fc_[upper].send[n], p, MPI_ATHENA_REAL, snb.rank,
                  tag, MPI_COMM_WORLD, &orbital_bd_fc_[upper].req_send[n]);
#endif //MPI
      }
      orbital_bd_fc_[upper].sflag[n] = BoundaryStatus::completed;
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn bool OrbitalBoundaryCommunication::ReceiveBoundaryBuffersFC()
//! \brief receive and set magnetic fields
bool OrbitalBoundaryCommunication::ReceiveBoundaryBuffersFC() {
  bool flag[2]{true, true};
  int mylevel = pmy_block_->loc.level;
  int offset[2]{0,4};
  for (int upper=0; upper<2; upper++) {
    for (int n=0; n<orbital_bd_fc_[upper].nbmax; n++) {
      if (orbital_bd_fc_[upper].flag[n] == BoundaryStatus::completed) continue;
      SimpleNeighborBlock& snb= orbital_recv_neighbor_[upper][n];
      if (orbital_bd_fc_[upper].flag[n] == BoundaryStatus::waiting) {
        // on the same process
        if (snb.rank == Globals::my_rank) {
          flag[upper] = false;
          continue;
        } else { // MPI
#ifdef MPI_PARALLEL
          int test;
          MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &test,
                     MPI_STATUS_IGNORE);
          MPI_Test(&orbital_bd_fc_[upper].req_recv[n], &test, MPI_STATUS_IGNORE);
          if (!static_cast<bool>(test)) {
            flag[upper] = false;
            continue;
          }
#endif
          orbital_bd_fc_[upper].flag[n] = BoundaryStatus::arrived;
        }
      }
      // set var
      int p=0;
      if (snb.level == mylevel) { // from same level
        SetFieldBufferSameLevel(orbital_bd_fc_[upper].recv[n], p, n+offset[upper]);
      } else if (snb.level < mylevel) { // from coarser
        SetFieldBufferFromCoarser(orbital_bd_fc_[upper].recv[n], p, n+offset[upper]);
      } else { // from finer
        SetFieldBufferFromFiner(orbital_bd_fc_[upper].recv[n], p, n+offset[upper]);
      }
      if (p != orbital_recv_fc_count_[upper][n]) {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::ReceiveBoundaryBuffersFC" << std::endl
            << "Recieve buffer size is incorrect." << std::endl;
        ATHENA_ERROR(msg);
      }
      orbital_bd_fc_[upper].flag[n] = BoundaryStatus::completed;
    } // loop over recv[0] to recv[nbmax-1]
  } // loop over comm direction
  return (flag[0] && flag[1]);
}

//----------------------------------------------------------------------------------------
//! \fn void OrbitalBoundaryCommunication::LoadHydroBufferSameLevel(
//!                                       Real *buf, int &p, int nb)
//! \brief load hydro variables (same level)
void OrbitalBoundaryCommunication::LoadHydroBufferSameLevel(Real *buf, int &p, int nb) {
  MeshBlock *pmb = pmy_block_;
  OrbitalAdvection *porb = pmy_orbital_;
  AthenaArray<Real> &ui = pmb->phydro->u;

  if(porb->orbital_uniform_mesh) { // uniform mesh
    if(porb->orbital_direction == 1) {
      int &onx = pmb->block_size.nx2;
      if (nb==0) {
        for(int k=pmb->ks; k<=pmb->ke; k++) {
          for(int i=pmb->is; i<=pmb->ie; i++) {
            int offset = porb->ofc(k,i);
            int xl = pmb->js-xgh-offset+onx;
            for(int nph=0 ; nph<NHYDRO; nph++) {
              for (int j=xl; j<=pmb->je; j++) {
                buf[p++] = ui(nph,k,j,i);
              }
            }
          }
        }
      } else if (nb==4) {
        for(int k=pmb->ks; k<=pmb->ke; k++) {
          for(int i=pmb->is; i<=pmb->ie; i++) {
            int offset = porb->ofc(k,i);
            int xu = pmb->je+1+xgh-offset-onx;
            for(int nph=0 ; nph<NHYDRO; nph++) {
              for (int j=pmb->js; j<=xu; j++) {
                buf[p++] = ui(nph,k,j,i);
              }
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::LoadHydroBufferSameLevel" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    } else if (porb->orbital_direction == 2) {
      int &onx = pmb->block_size.nx3;
      if (nb==0) {
        for(int j=pmb->js; j<=pmb->je; j++) {
          for(int i=pmb->is; i<=pmb->ie; i++) {
            int offset = porb->ofc(j,i);
            int xl = pmb->ks-xgh-offset+onx;
            for(int nph=0 ; nph<NHYDRO; nph++) {
              for (int k=xl; k<=pmb->ke; k++) {
                buf[p++] = ui(nph,k,j,i);
              }
            }
          }
        }
      } else if (nb==4) {
        for(int j=pmb->js; j<=pmb->je; j++) {
          for(int i=pmb->is; i<=pmb->ie; i++) {
            int offset = porb->ofc(j,i);
            int xu = pmb->ke+1+xgh-offset-onx;
            for(int nph=0 ; nph<NHYDRO; nph++) {
              for (int k=pmb->ks; k<=xu; k++) {
                buf[p++] = ui(nph,k,j,i);
              }
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::LoadHydroBufferSameLevel" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    }
  }
//  else { // non-uniform mesh
//  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OrbitalBoundaryCommunication::LoadHydroBufferToCoarser(
//!                                       Real *buf, int &p, int nb)
//! \brief load hydro variables (coarser level)
void OrbitalBoundaryCommunication::LoadHydroBufferToCoarser(Real *buf, int &p, int nb) {
  MeshBlock *pmb = pmy_block_;
  OrbitalAdvection *porb = pmy_orbital_;
  AthenaArray<Real> &ui = porb->u_coarse_send;

  if(porb->orbital_uniform_mesh) { // uniform mesh
    if(porb->orbital_direction == 1) {
      int honx = pmb->block_size.nx2/2;
      if (nb==0) {
        for(int k=pmb->cks; k<=pmb->cke; k++) {
          for(int i=pmb->cis; i<=pmb->cie; i++) {
            int offset = porb->ofc_coarse(k,i);
            int xl = pmb->cjs-xgh-offset+honx;
            for(int nph=0 ; nph<NHYDRO; nph++) {
              for (int j=xl; j<=pmb->cje; j++) {
                buf[p++] = ui(nph,k,j,i);
              }
            }
          }
        }
      } else if (nb==4) {
        for(int k=pmb->cks; k<=pmb->cke; k++) {
          for(int i=pmb->cis; i<=pmb->cie; i++) {
            int offset = porb->ofc_coarse(k,i);
            int xu = pmb->cje+1+xgh-offset-honx;
            for(int nph=0 ; nph<NHYDRO; nph++) {
              for (int j=pmb->cjs; j<=xu; j++) {
                buf[p++] = ui(nph,k,j,i);
              }
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::LoadHydroBufferToCoarser" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    } else if (porb->orbital_direction == 2) {
      int honx = pmb->block_size.nx3/2;
      if (nb==0) {
        for(int j=pmb->cjs; j<=pmb->cje; j++) {
          for(int i=pmb->cis; i<=pmb->cie; i++) {
            int offset = porb->ofc_coarse(j,i);
            int xl = pmb->cks-xgh-offset+honx;
            for(int nph=0 ; nph<NHYDRO; nph++) {
              for (int k=xl; k<=pmb->cke; k++) {
                buf[p++] = ui(nph,k,j,i);
              }
            }
          }
        }
      } else if (nb==4) {
        for(int j=pmb->cjs; j<=pmb->cje; j++) {
          for(int i=pmb->cis; i<=pmb->cie; i++) {
            int offset = porb->ofc_coarse(j,i);
            int xu = pmb->cke+1+xgh-offset-honx;
            for(int nph=0 ; nph<NHYDRO; nph++) {
              for (int k=pmb->cks; k<=xu; k++) {
                buf[p++] = ui(nph,k,j,i);
              }
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::LoadHydroBufferToCoarser" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    }
  }
//  else { // non-uniform mesh
//  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OrbitalBoundaryCommunication::LoadHydroBufferToFiner(
//!                                     Real *buf, int &p, int nb)
//! \brief load hydro variables (finer level)
void OrbitalBoundaryCommunication::LoadHydroBufferToFiner(Real *buf, int &p, int nb) {
  MeshBlock *pmb = pmy_block_;
  OrbitalAdvection *porb = pmy_orbital_;
  AthenaArray<Real> &ui = pmb->phydro->u;

  if(porb->orbital_uniform_mesh) { // uniform mesh
    if(porb->orbital_direction == 1) {
      int &onx = pmb->block_size.nx2;
      int il, iu, jl, ju, kl, ku;
      if (pmb->block_size.nx3>1) {
        if(nb%4<2) {
          kl = pmb->ks-1;
          ku = pmb->ke+1-pmb->block_size.nx3/2;
        } else {
          kl = pmb->ks-1+pmb->block_size.nx3/2;
          ku = pmb->ke+1;
        }
      } else {
          kl = pmb->ks;
          ku = pmb->ke;
      }
      if (nb%2==0) {
        il = pmb->is-1;
        iu = pmb->ie+1-pmb->block_size.nx1/2;
      } else {
        il = pmb->is-1+pmb->block_size.nx1/2;
        iu = pmb->ie+1;
      }
      if (nb<4) {
        jl = pmb->js-size_cc_send[0][2+nb]+1+onx;
        ju = pmb->je+1;
      } else if (nb<8) {
        jl = pmb->js-1;
        ju = pmb->je+size_cc_send[1][2+nb-4]-onx-1;
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::LoadHydroBufferToFiner" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
      BufferUtility::PackData(ui, buf, 0, NHYDRO-1,
                              il, iu, jl, ju, kl, ku, p);
    } else if (porb->orbital_direction == 2) {
      int &onx = pmb->block_size.nx3;
      int il, iu, jl, ju, kl, ku;
      if(nb%4<2) {
        jl = pmb->js-1;
        ju = pmb->je+1-pmb->block_size.nx2/2;
      } else {
        jl = pmb->js-1+pmb->block_size.nx2/2;
        ju = pmb->je+1;
      }
      if (nb%2==0) {
        il = pmb->is-1;
        iu = pmb->ie+1-pmb->block_size.nx1/2;
      } else {
        il = pmb->is-1+pmb->block_size.nx1/2;
        iu = pmb->ie+1;
      }
      if (nb<4) {
        kl = pmb->ks-size_cc_send[0][2+nb]+1+onx;
        ku = pmb->ke+1;
      } else if (nb<8) {
        kl = pmb->ks-1;
        ku = pmb->ke+size_cc_send[1][2+nb-4]-onx-1;
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::LoadHydroBufferToFiner" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
      BufferUtility::PackData(ui, buf, 0, NHYDRO-1,
                              il, iu, jl, ju, kl, ku, p);
    }
  }
//  else { // non-uniform mesh
//  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OrbitalBoundaryCommunication::SetHydroBufferSameLevel(
//!                                Real *buf, int &p, const int nb)
//! \brief set hydro variables (same level)
void OrbitalBoundaryCommunication::SetHydroBufferSameLevel(
                           Real *buf, int &p, const int nb) {
  MeshBlock *pmb = pmy_block_;
  OrbitalAdvection *porb = pmy_orbital_;
  AthenaArray<Real> &uo = porb->orbital_cons;

  if(porb->orbital_uniform_mesh) { // uniform mesh
    if(porb->orbital_direction == 1) {
      int &onx = pmb->block_size.nx2;
      if (nb==0) {
        for(int k=pmb->ks; k<=pmb->ke; k++) {
          for(int i=pmb->is; i<=pmb->ie; i++) {
            int offset = porb->ofc(k,i);
            int xl = pmb->js-xgh-offset+onx;
            int xu = pmb->je;
            if (offset<=0) {
              xl -= onx; xu -= onx;
            }
            for(int nph=0 ; nph<NHYDRO; nph++) {
#pragma omp simd
              for (int j=xl; j<=xu; j++) {
                uo(nph,k,i,j) = buf[p++];
              }
            }
          }
        }
      } else if (nb==4) {
        for(int k=pmb->ks; k<=pmb->ke; k++) {
          for(int i=pmb->is; i<=pmb->ie; i++) {
            int offset = porb->ofc(k,i);
            int xl = pmb->js+onx;
            int xu = pmb->je+1+xgh-offset;
            if (offset>0) {
              xl += onx; xu += onx;
            }
            for(int nph=0 ; nph<NHYDRO; nph++) {
#pragma omp simd
              for (int j=xl; j<=xu; j++) {
                uo(nph,k,i,j) = buf[p++];
              }
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::SetHydroBufferSameLevel" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    } else if(porb->orbital_direction == 2) {
      int &onx = pmb->block_size.nx3;
      if (nb==0) {
        for(int j=pmb->js; j<=pmb->je; j++) {
          for(int i=pmb->is; i<=pmb->ie; i++) {
            int offset = porb->ofc(j,i);
            int xl = pmb->ks-xgh-offset+onx;
            int xu = pmb->ke;
            if (offset<=0) {
              xl -= onx; xu -= onx;
            }
            for(int nph=0 ; nph<NHYDRO; nph++) {
#pragma omp simd
              for (int k=xl; k<=xu; k++) {
                uo(nph,j,i,k) = buf[p++];
              }
            }
          }
        }
      } else if (nb==4) {
        for(int j=pmb->js; j<=pmb->je; j++) {
          for(int i=pmb->is; i<=pmb->ie; i++) {
            int offset = porb->ofc(j,i);
            int xl = pmb->ks+onx;
            int xu = pmb->ke+1+xgh-offset;
            if (offset>0) {
              xl += onx; xu += onx;
            }
            for(int nph=0 ; nph<NHYDRO; nph++) {
#pragma omp simd
              for (int k=xl; k<=xu; k++) {
                uo(nph,j,i,k) = buf[p++];
              }
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::SetHydroBufferSameLevel" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OrbitalBoundaryCommunication::SetHydroBufferFromCoarser(
//!                                   Real *buf, int &p, const int nb)
//! \brief set hydro variables (coarser level)
void OrbitalBoundaryCommunication::SetHydroBufferFromCoarser(
                             Real *buf, int &p, const int nb) {
  MeshBlock *pmb = pmy_block_;
  OrbitalAdvection *porb = pmy_orbital_;
  AthenaArray<Real> &uo  = porb->orbital_cons;
  AthenaArray<Real> &uco = porb->u_coarse_recv;
  AthenaArray<Real> &uto = porb->u_temp;

  if(porb->orbital_uniform_mesh) { // uniform mesh
    if(porb->orbital_direction == 1) {
      int &onx = pmb->block_size.nx2;
      int il, iu, kl, ku;
      if (pmb->block_size.nx3>1) {
        kl = pmb->cks-1;
        ku = pmb->cke+1;
      } else {
        kl = pmb->cks;
        ku = pmb->cke;
      }
      il = pmb->cis-1;
      iu = pmb->cie+1;
      if (nb==0) {
        int jl = pmb->cjs-xgh-porb->max_ofc_coarse+onx/2-1;
        int ju = pmb->cje+1;
        // TODO(tomo-ono): This part has a problem with "#pragma omp simd"
        //                 when using the Intel compiler
        // BufferUtility::UnpackData(buf, uco, 0, NHYDRO-1,
        //                           il, iu, jl, ju, kl, ku, p);
        for (int n=0; n<NHYDRO; ++n) {
          for (int k=kl; k<=ku; k++) {
            for (int j=jl; j<=ju; j++) {
              for (int i=il; i<=iu; i++) {
                uco(n,k,j,i) = buf[p++];
              }
            }
          }
        }
        pmb->pmr->ProlongateCellCenteredValues(uco, uto, 0, NHYDRO-1,
                                               pmb->cis, pmb->cie, jl+1,
                                               pmb->cje, pmb->cks, pmb->cke);
        for (int k=pmb->ks; k<=pmb->ke; k++) {
          for (int i=pmb->is; i<=pmb->ie; i++) {
            int offset = porb->ofc(k,i);
            int xl = pmb->js-xgh-offset+onx;
            int xu = pmb->je;
            const int shift = (offset>0)? 0: -onx;
            for(int nph=0 ; nph<NHYDRO; nph++) {
              for (int j=xl; j<=xu; j++) {
                uo(nph,k,i,j+shift) = uto(nph,k,j,i);
              }
            }
          }
        }
      } else if (nb==4) {
        int jl = pmb->cjs-1;
        int ju = pmb->cje+2+xgh-porb->min_ofc_coarse-onx/2;
        // TODO(tomo-ono): This part has a problem with "#pragma omp simd"
        //                 when using the Intel compiler
        // BufferUtility::UnpackData(buf, uco, 0, NHYDRO-1,
        //                           il, iu, jl, ju, kl, ku, p);
        for (int n=0; n<NHYDRO; ++n) {
          for (int k=kl; k<=ku; k++) {
            for (int j=jl; j<=ju; j++) {
              for (int i=il; i<=iu; i++) {
                uco(n,k,j,i) = buf[p++];
              }
            }
          }
        }
        pmb->pmr->ProlongateCellCenteredValues(uco, uto, 0, NHYDRO-1,
                                               pmb->cis, pmb->cie, pmb->cjs,
                                               ju-1, pmb->cks, pmb->cke);
        for(int k=pmb->ks; k<=pmb->ke; k++) {
          for(int i=pmb->is; i<=pmb->ie; i++) {
            int offset = porb->ofc(k,i);
            int xl = pmb->js;
            int xu = pmb->je+1+xgh-offset-onx;
            const int shift = (offset>0)? 2*onx: onx;
            for(int nph=0 ; nph<NHYDRO; nph++) {
              for (int j=xl; j<=xu; j++) {
                uo(nph,k,i,j+shift) = uto(nph,k,j,i);
              }
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::SetHydroBufferFromCoarser" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    } else if(porb->orbital_direction == 2) {
      int &onx = pmb->block_size.nx3;
      int il, iu, jl, ju;
      il = pmb->cis-1;
      iu = pmb->cie+1;
      jl = pmb->cjs-1;
      ju = pmb->cje+1;
      if(nb==0) {
        int kl = pmb->cks-xgh-porb->max_ofc_coarse+onx/2-1;
        int ku = pmb->cke+1;
        // TODO(tomo-ono): This part has a problem with "#pragma omp simd"
        //                 when using the Intel compiler
        // BufferUtility::UnpackData(buf, uco, 0, NHYDRO-1,
        //                           il, iu, jl, ju, kl, ku, p);
        for (int n=0; n<NHYDRO; ++n) {
          for (int k=kl; k<=ku; k++) {
            for (int j=jl; j<=ju; j++) {
              for (int i=il; i<=iu; i++) {
                uco(n,k,j,i) = buf[p++];
              }
            }
          }
        }
        pmb->pmr->ProlongateCellCenteredValues(uco, uto, 0, NHYDRO-1,
                                               pmb->cis, pmb->cie, pmb->cjs,
                                               pmb->cje, kl+1, pmb->cke);
        for(int j=pmb->js; j<=pmb->je; j++) {
          for(int i=pmb->is; i<=pmb->ie; i++) {
            int offset = porb->ofc(j,i);
            int xl = pmb->ks-xgh-offset+onx;
            int xu = pmb->ke;
            const int shift = (offset>0)? 0: -onx;
            for(int nph=0 ; nph<NHYDRO; nph++) {
              for (int k=xl; k<=xu; k++) {
                uo(nph,j,i,k+shift) = uto(nph,k,j,i);
              }
            }
          }
        }
      } else if (nb==4) {
        int kl = pmb->cks-1;
        int ku = pmb->cke+2+xgh-porb->min_ofc_coarse-onx/2;
        // TODO(tomo-ono): This part has a problem with "#pragma omp simd"
        //                 when using the Intel compiler
        // BufferUtility::UnpackData(buf, uco, 0, NHYDRO-1,
        //                           il, iu, jl, ju, kl, ku, p);
        for (int n=0; n<NHYDRO; ++n) {
          for (int k=kl; k<=ku; k++) {
            for (int j=jl; j<=ju; j++) {
              for (int i=il; i<=iu; i++) {
                uco(n,k,j,i) = buf[p++];
              }
            }
          }
        }
        pmb->pmr->ProlongateCellCenteredValues(uco, uto, 0, NHYDRO-1,
                                               pmb->cis, pmb->cie, pmb->cjs,
                                               pmb->cje, pmb->cks, ku-1);
        for(int j=pmb->js; j<=pmb->je; j++) {
          for(int i=pmb->is; i<=pmb->ie; i++) {
            int offset = porb->ofc(j,i);
            int xl = pmb->ks;
            int xu = pmb->ke+1+xgh-offset-onx;
            const int shift = (offset>0)? 2*onx: onx;
            for(int nph=0 ; nph<NHYDRO; nph++) {
              for (int k=xl; k<=xu; k++) {
                uo(nph,j,i,k+shift) = uto(nph,k,j,i);
              }
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::SetHydroBufferFromCoarser" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OrbitalBoundaryCommunication::SetHydroBufferFromFiner(
//!                                Real *buf, int &p, const int nb)
//! \brief set hydro variables (finer level)
void OrbitalBoundaryCommunication::SetHydroBufferFromFiner(
                           Real *buf, int &p, const int nb) {
  MeshBlock *pmb = pmy_block_;
  OrbitalAdvection *porb = pmy_orbital_;
  AthenaArray<Real> &uo = porb->orbital_cons;

  if(porb->orbital_uniform_mesh) { // uniform mesh
    if(porb->orbital_direction == 1) {
      int &onx = pmb->block_size.nx2;
      int hnx1 = pmb->block_size.nx1/2;
      int hnx3 = pmb->block_size.nx3/2;
      int il, iu, kl, ku;
      if (pmb->block_size.nx3>1) {
        if (nb%4<2) {
          kl = pmb->ks;
          ku = pmb->ke-hnx3;
        } else {
          kl = pmb->ks+hnx3;
          ku = pmb->ke;
        }
      } else {
        if (nb%4>1) {
          std::stringstream msg;
          msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
              << "::SetHydroBufferFromFiner" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
          ATHENA_ERROR(msg);
        }
        kl = pmb->ks;
        ku = pmb->ke;
      }
      if (nb%2==0) {
        il = pmb->is;
        iu = pmb->ie-hnx1;
      } else {
        il = pmb->is+hnx1;
        iu = pmb->ie;
      }
      if (nb<4) {
        for(int k=kl; k<=ku; k++) {
          for(int i=il; i<=iu; i++) {
            int offset = porb->ofc(k,i);
            int xl = pmb->js-xgh-offset+onx;
            int xu = pmb->je;
            if (offset<=0) {
              xl -= onx; xu -= onx;
            }
            for(int nph=0 ; nph<NHYDRO; nph++) {
#pragma omp simd
              for (int j=xl; j<=xu; j++) {
                uo(nph,k,i,j) = buf[p++];
              }
            }
          }
        }
      } else if (nb<8) {
        for(int k=kl; k<=ku; k++) {
          for(int i=il; i<=iu; i++) {
            int offset = porb->ofc(k,i);
            int xl = pmb->js+onx;
            int xu = pmb->je+1+xgh-offset;
            if (offset>0) {
              xl += onx; xu += onx;
            }
            for(int nph=0 ; nph<NHYDRO; nph++) {
#pragma omp simd
              for (int j=xl; j<=xu; j++) {
                uo(nph,k,i,j) = buf[p++];
              }
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::SetHydroBufferFromFiner" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    } else if(porb->orbital_direction == 2) {
      int &onx = pmb->block_size.nx3;
      int hnx1 = pmb->block_size.nx1/2;
      int hnx2 = pmb->block_size.nx2/2;
      int il, iu, jl, ju;
      if (nb%4<2) {
        jl = pmb->js;
        ju = pmb->je-hnx2;
      } else {
        jl = pmb->js+hnx2;
        ju = pmb->je;
      }
      if (nb%2==0) {
        il = pmb->is;
        iu = pmb->ie-hnx1;
      } else {
        il = pmb->is+hnx1;
        iu = pmb->ie;
      }
      if (nb<4) {
        for(int j=jl; j<=ju; j++) {
          for(int i=il; i<=iu; i++) {
            int offset = porb->ofc(j,i);
            int xl = pmb->ks-xgh-offset+onx;
            int xu = pmb->ke;
            if (offset<=0) {
              xl -= onx; xu -= onx;
            }
            for(int nph=0 ; nph<NHYDRO; nph++) {
#pragma omp simd
              for (int k=xl; k<=xu; k++) {
                uo(nph,j,i,k) = buf[p++];
              }
            }
          }
        }
      } else if (nb<8) {
        for(int j=jl; j<=ju; j++) {
          for(int i=il; i<=iu; i++) {
            int offset = porb->ofc(j,i);
            int xl = pmb->ks+onx;
            int xu = pmb->ke+1+xgh-offset;
            if (offset>0) {
              xl += onx; xu += onx;
            }
            for(int nph=0 ; nph<NHYDRO; nph++) {
#pragma omp simd
              for (int k=xl; k<=xu; k++) {
                uo(nph,j,i,k) = buf[p++];
              }
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::SetHydroBufferFromFiner" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OrbitalBoundaryCommunication::LoadFieldBufferSameLevel(
//!                                       Real *buf, int &p, int nb)
//! \brief load magnetic fields (same level)
void OrbitalBoundaryCommunication::LoadFieldBufferSameLevel(Real *buf, int &p, int nb) {
  MeshBlock *pmb = pmy_block_;
  OrbitalAdvection *porb = pmy_orbital_;
  FaceField &bi = pmb->pfield->b;

  if(porb->orbital_uniform_mesh) { // uniform mesh
    if(porb->orbital_direction == 1) {
      int &onx = pmb->block_size.nx2;
      if (nb==0) {
        // b1
        for(int k=pmb->ks; k<=pmb->ke  ; k++) {
          for(int i=pmb->is; i<=pmb->ie+1; i++) {
            int offset = porb->off[0](k,i);
            int xl = pmb->js-xgh-offset+onx;
            for (int j=xl; j<=pmb->je; j++) {
              buf[p++] = bi.x1f(k,j,i);
            }
          }
        }
        // -b3
        for(int k=pmb->ks; k<=pmb->ke+1; k++) {
          for(int i=pmb->is; i<=pmb->ie  ; i++) {
            int offset = porb->off[1](k,i);
            int xl = pmb->js-xgh-offset+onx;
            for (int j=xl; j<=pmb->je; j++) {
              buf[p++] = -bi.x3f(k,j,i);
            }
          }
        }
      } else if (nb==4) {
        // b1
        for(int k=pmb->ks; k<=pmb->ke  ; k++) {
          for(int i=pmb->is; i<=pmb->ie+1; i++) {
            int offset = porb->off[0](k,i);
            int xu = pmb->je+1+xgh-offset-onx;
            for (int j=pmb->js; j<=xu; j++) {
              buf[p++] = bi.x1f(k,j,i);
            }
          }
        }
        // -b3
        for(int k=pmb->ks; k<=pmb->ke+1; k++) {
          for(int i=pmb->is; i<=pmb->ie  ; i++) {
            int offset = porb->off[1](k,i);
            int xu = pmb->je+1+xgh-offset-onx;
            for (int j=pmb->js; j<=xu; j++) {
              buf[p++] = -bi.x3f(k,j,i);
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::LoadFieldBufferSameLevel" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    } else if (porb->orbital_direction == 2) {
      int &onx = pmb->block_size.nx3;
      if (nb==0) {
        // -b1
        for(int j=pmb->js; j<=pmb->je  ; j++) {
          for(int i=pmb->is; i<=pmb->ie+1; i++) {
            int offset = porb->off[0](j,i);
            int xl = pmb->ks-xgh-offset+onx;
            for (int k=xl; k<=pmb->ke; k++) {
              buf[p++] = -bi.x1f(k,j,i);
            }
          }
        }
        // b2
        for(int j=pmb->js; j<=pmb->je+1; j++) {
          for(int i=pmb->is; i<=pmb->ie  ; i++) {
            int offset = porb->off[1](j,i);
            int xl = pmb->ks-xgh-offset+onx;
            for (int k=xl; k<=pmb->ke; k++) {
              buf[p++] = bi.x2f(k,j,i);
            }
          }
        }
      } else if (nb==4) {
        // -b1
        for(int j=pmb->js; j<=pmb->je  ; j++) {
          for(int i=pmb->is; i<=pmb->ie+1; i++) {
            int offset = porb->off[0](j,i);
            int xu = pmb->ke+1+xgh-offset-onx;
            for (int k=pmb->ks; k<=xu; k++) {
              buf[p++] = -bi.x1f(k,j,i);
            }
          }
        }
        // b2
        for(int j=pmb->js; j<=pmb->je+1; j++) {
          for(int i=pmb->is; i<=pmb->ie  ; i++) {
            int offset = porb->off[1](j,i);
            int xu = pmb->ke+1+xgh-offset-onx;
            for (int k=pmb->ks; k<=xu; k++) {
              buf[p++] = bi.x2f(k,j,i);
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::LoadFieldBufferSameLevel" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    }
  }
//  else { // non-uniform mesh
//  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OrbitalBoundaryCommunication::LoadFieldBufferToCoarser(
//!                                       Real *buf, int &p, int nb)
//! \brief load magnetic fields (coarser level)
void OrbitalBoundaryCommunication::LoadFieldBufferToCoarser(Real *buf, int &p, int nb) {
  MeshBlock *pmb = pmy_block_;
  OrbitalAdvection *porb = pmy_orbital_;
  FaceField &bi = pmb->pfield->b;

  if(porb->orbital_uniform_mesh) { // uniform mesh
    if(porb->orbital_direction == 1) {
      int onx = pmb->block_size.nx2/2;
      if (nb==0) {
        // b1
        for(int k=pmb->cks; k<=pmb->cke  ; k++) {
          for(int i=pmb->cis; i<=pmb->cie+1; i++) {
            int offset = porb->off_coarse[0](k,i);
            int xl = pmb->cjs-xgh-offset+onx;
            for (int j=xl; j<=pmb->cje; j++) {
              buf[p++] = porb->b1_coarse_send(k,j,i);
            }
          }
        }
        // -b3
        for(int k=pmb->cks; k<=pmb->cke+1; k++) {
          for(int i=pmb->cis; i<=pmb->cie  ; i++) {
            int offset = porb->off_coarse[1](k,i);
            int xl = pmb->cjs-xgh-offset+onx;
            for (int j=xl; j<=pmb->cje; j++) {
              buf[p++] = -porb->b2_coarse_send(k,j,i);
            }
          }
        }
      } else if (nb==4) {
        // b1
        for(int k=pmb->cks; k<=pmb->cke  ; k++) {
          for(int i=pmb->cis; i<=pmb->cie+1; i++) {
            int offset = porb->off_coarse[0](k,i);
            int xu = pmb->cje+1+xgh-offset-onx;
            for (int j=pmb->cjs; j<=xu; j++) {
              buf[p++] = porb->b1_coarse_send(k,j,i);
            }
          }
        }
        // -b3
        for(int k=pmb->cks; k<=pmb->cke+1; k++) {
          for(int i=pmb->cis; i<=pmb->cie  ; i++) {
            int offset = porb->off_coarse[1](k,i);
            int xu = pmb->cje+1+xgh-offset-onx;
            for (int j=pmb->cjs; j<=xu; j++) {
              buf[p++] = -porb->b2_coarse_send(k,j,i);
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::LoadFieldBufferToCoarser" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    } else if (porb->orbital_direction == 2) {
      int onx = pmb->block_size.nx3/2;
      if (nb==0) {
        // -b1
        for(int j=pmb->cjs; j<=pmb->cje  ; j++) {
          for(int i=pmb->cis; i<=pmb->cie+1; i++) {
            int offset = porb->off_coarse[0](j,i);
            int xl = pmb->cks-xgh-offset+onx;
            for (int k=xl; k<=pmb->cke; k++) {
              buf[p++] = -porb->b1_coarse_send(k,j,i);
            }
          }
        }
        // b2
        for(int j=pmb->cjs; j<=pmb->cje+1; j++) {
          for(int i=pmb->cis; i<=pmb->cie  ; i++) {
            int offset = porb->off_coarse[1](j,i);
            int xl = pmb->cks-xgh-offset+onx;
            for (int k=xl; k<=pmb->cke; k++) {
              buf[p++] = porb->b2_coarse_send(k,j,i);
            }
          }
        }
      } else if (nb==4) {
        // -b1
        for(int j=pmb->cjs; j<=pmb->cje  ; j++) {
          for(int i=pmb->cis; i<=pmb->cie+1; i++) {
            int offset = porb->off_coarse[0](j,i);
            int xu = pmb->cke+1+xgh-offset-onx;
            for (int k=pmb->cks; k<=xu; k++) {
              buf[p++] = -porb->b1_coarse_send(k,j,i);
            }
          }
        }
        // b2
        for(int j=pmb->cjs; j<=pmb->cje+1; j++) {
          for(int i=pmb->cis; i<=pmb->cie  ; i++) {
            int offset = porb->off_coarse[1](j,i);
            int xu = pmb->cke+1+xgh-offset-onx;
            for (int k=pmb->cks; k<=xu; k++) {
              buf[p++] = porb->b2_coarse_send(k,j,i);
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::LoadFieldBufferToCoarser" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    }
  }
//  else { // non-uniform mesh
//  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OrbitalBoundaryCommunication::LoadFieldBufferToFiner(
//!                                     Real *buf, int &p, int nb)
//! \brief load magnetic fields (finer level)
void OrbitalBoundaryCommunication::LoadFieldBufferToFiner(Real *buf, int &p, int nb) {
  MeshBlock *pmb = pmy_block_;
  OrbitalAdvection *porb = pmy_orbital_;
  FaceField &bi = pmb->pfield->b;

  if(porb->orbital_uniform_mesh) { // uniform mesh
    int onx;
    int il=pmb->is, iu=pmb->ie, jl=pmb->js, ju=pmb->je, kl=pmb->ks, ku=pmb->ke;
    if (pmb->block_size.nx3>1) { // 3D
      if(porb->orbital_direction == 1) {
        onx = pmb->block_size.nx2;
        if(nb%4<2) {
          ku -= pmb->block_size.nx3/2;
        } else {
          kl += pmb->block_size.nx3/2;
        }
        if (nb%2==0) {
          iu -= pmb->block_size.nx1/2;
        } else {
          il += pmb->block_size.nx1/2;
        }
        if (nb<4) {
          jl += -size_fc_send[0][2+nb]+onx;
        } else if (nb<8) {
          ju += size_fc_send[1][2+(nb-4)]-onx;
        } else {
          std::stringstream msg;
          msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
              << "::LoadFieldBufferToFiner" << std::endl
              << "Neighbors are read incorrectly." << std::endl;
          ATHENA_ERROR(msg);
        }
      } else if (porb->orbital_direction == 2) {
        onx = pmb->block_size.nx3;
        if(nb%4<2) {
          ju -= pmb->block_size.nx2/2;
        } else {
          jl += pmb->block_size.nx2/2;
        }
        if (nb%2==0) {
          iu -= pmb->block_size.nx1/2;
        } else {
          il += pmb->block_size.nx1/2;
        }
        if (nb<4) {
          kl += -size_fc_send[0][2+nb]+onx;
        } else if (nb<8) {
          ku += size_fc_send[1][2+(nb-4)]-onx;
        } else {
          std::stringstream msg;
          msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
              << "::LoadFieldBufferToFiner" << std::endl
              << "Neighbors are read incorrectly." << std::endl;
          ATHENA_ERROR(msg);
        }
      }

      // pack b1
      BufferUtility::PackData(bi.x1f, buf,
                              il, iu+1, jl-1, ju+1, kl-1, ku+1, p);
      // pack b2
      BufferUtility::PackData(bi.x2f, buf,
                              il-1, iu+1, jl, ju+1, kl-1, ku+1, p);
      // pack b3
      BufferUtility::PackData(bi.x3f, buf,
                              il-1, iu+1, jl-1, ju+1, kl, ku+1, p);
    } else { // 2D
      if(porb->orbital_direction == 1) {
        onx = pmb->block_size.nx2;
        if (nb%2==0) {
          iu -= pmb->block_size.nx1/2;
        } else {
          il += pmb->block_size.nx1/2;
        }
        if (nb<4) {
          jl += -size_fc_send[0][2+nb]+onx;
        } else if (nb<8) {
          ju += size_fc_send[1][2+(nb-4)]-onx;
        } else {
          std::stringstream msg;
          msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
              << "::LoadFieldBufferToFiner" << std::endl
              << "Neighbors are read incorrectly." << std::endl;
          ATHENA_ERROR(msg);
        }
      } else {
          std::stringstream msg;
          msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
              << "::LoadFieldBufferToFiner" << std::endl
              << "2D Orbital Advection is not allowed in spherical_polar." << std::endl;
          ATHENA_ERROR(msg);
      }
      // pack b1
      BufferUtility::PackData(bi.x1f, buf,
                              il, iu+1, jl-1, ju+1, kl, ku, p);
      // pack b2
      BufferUtility::PackData(bi.x2f, buf,
                              il-1, iu+1, jl, ju+1, kl, ku, p);
      // pack b3
      BufferUtility::PackData(bi.x3f, buf,
                              il-1, iu+1, jl-1, ju+1, kl, ku, p);
    }
  }
//  else { // non-uniform mesh
//  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OrbitalBoundaryCommunication::SetFieldBufferSameLevel(
//!                                Real *buf, int &p, const int nb)
//! \brief set magnetic fields (same level)
void OrbitalBoundaryCommunication::SetFieldBufferSameLevel(
                           Real *buf, int &p, const int nb) {
  MeshBlock *pmb = pmy_block_;
  OrbitalAdvection *porb = pmy_orbital_;
  AthenaArray<Real> &bo1 = porb->orbital_b1;
  AthenaArray<Real> &bo2 = porb->orbital_b2;

  if(porb->orbital_uniform_mesh) { // uniform mesh
    if(porb->orbital_direction == 1) {
      int &onx = pmb->block_size.nx2;
      if (nb==0) {
        for(int k=pmb->ks; k<=pmb->ke  ; k++) {
          for(int i=pmb->is; i<=pmb->ie+1; i++) {
            int offset = porb->off[0](k,i);
            int xl = pmb->js-xgh-offset+onx;
            int xu = pmb->je;
            if (offset<=0) {
              xl -= onx; xu -= onx;
            }
#pragma omp simd
            for (int j=xl; j<=xu; j++) {
              bo1(k,i,j) = buf[p++];
            }
          }
        }
        for(int k=pmb->ks; k<=pmb->ke+1; k++) {
          for(int i=pmb->is; i<=pmb->ie  ; i++) {
            int offset = porb->off[1](k,i);
            int xl = pmb->js-xgh-offset+onx;
            int xu = pmb->je;
            if (offset<=0) {
              xl -= onx; xu -= onx;
            }
#pragma omp simd
            for (int j=xl; j<=xu; j++) {
              bo2(k,i,j) = buf[p++];
            }
          }
        }
      } else if (nb==4) {
        for(int k=pmb->ks; k<=pmb->ke  ; k++) {
          for(int i=pmb->is; i<=pmb->ie+1; i++) {
            int offset = porb->off[0](k,i);
            int xl = pmb->js+onx;
            int xu = pmb->je+1+xgh-offset;
            if (offset>0) {
              xl += onx; xu += onx;
            }
#pragma omp simd
            for (int j=xl; j<=xu; j++) {
              bo1(k,i,j) = buf[p++];
            }
          }
        }
        for(int k=pmb->ks; k<=pmb->ke+1; k++) {
          for(int i=pmb->is; i<=pmb->ie  ; i++) {
            int offset = porb->off[1](k,i);
            int xl = pmb->js+onx;
            int xu = pmb->je+1+xgh-offset;
            if (offset>0) {
              xl += onx; xu += onx;
            }
#pragma omp simd
            for (int j=xl; j<=xu; j++) {
              bo2(k,i,j) = buf[p++];
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::SetFieldBufferSameLevel" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    } else if (porb->orbital_direction == 2) {
      int &onx = pmb->block_size.nx3;
      if (nb==0) {
        for(int j=pmb->js; j<=pmb->je  ; j++) {
          for(int i=pmb->is; i<=pmb->ie+1; i++) {
            int offset = porb->off[0](j,i);
            int xl = pmb->ks-xgh-offset+onx;
            int xu = pmb->ke;
            if (offset<=0) {
              xl -= onx; xu -= onx;
            }
#pragma omp simd
            for (int k=xl; k<=xu; k++) {
              bo1(j,i,k) = buf[p++];
            }
          }
        }
        for(int j=pmb->js; j<=pmb->je+1; j++) {
          for(int i=pmb->is; i<=pmb->ie  ; i++) {
            int offset = porb->off[1](j,i);
            int xl = pmb->ks-xgh-offset+onx;
            int xu = pmb->ke;
            if (offset<=0) {
              xl -= onx; xu -= onx;
            }
#pragma omp simd
            for (int k=xl; k<=xu; k++) {
              bo2(j,i,k) = buf[p++];
            }
          }
        }
      } else if (nb==4) {
        for(int j=pmb->js; j<=pmb->je  ; j++) {
          for(int i=pmb->is; i<=pmb->ie+1; i++) {
            int offset = porb->off[0](j,i);
            int xl = pmb->ks+onx;
            int xu = pmb->ke+1+xgh-offset;
            if (offset>0) {
              xl += onx; xu += onx;
            }
#pragma omp simd
            for (int k=xl; k<=xu; k++) {
              bo1(j,i,k) = buf[p++];
            }
          }
        }
        // b3
        for(int j=pmb->js; j<=pmb->je+1; j++) {
          for(int i=pmb->is; i<=pmb->ie  ; i++) {
            int offset = porb->off[1](j,i);
            int xl = pmb->ks+onx;
            int xu = pmb->ke+1+xgh-offset;
            if (offset>0) {
              xl += onx; xu += onx;
            }
#pragma omp simd
            for (int k=xl; k<=xu; k++) {
              bo1(j,i,k) = buf[p++];
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::SetFieldBufferSameLevel" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OrbitalBoundaryCommunication::SetFieldBufferFromCoarser(
//!                                  Real *buf, int &p, const int nb)
//! \brief set magnetic fields (coarser level)
void OrbitalBoundaryCommunication::SetFieldBufferFromCoarser(
                             Real *buf, int &p, const int nb) {
  MeshBlock *pmb = pmy_block_;
  OrbitalAdvection *porb = pmy_orbital_;
  AthenaArray<Real> &bo1 = porb->orbital_b1;
  AthenaArray<Real> &bo2 = porb->orbital_b2;
  FaceField &bto = porb->b_temp;
  FaceField &bco = porb->b_coarse_recv;

  if(porb->orbital_uniform_mesh) { // uniform mesh
    int onx = pmb->block_size.nx2;
    int il=pmb->cis, iu=pmb->cie, jl=pmb->cjs, ju=pmb->cje, kl=pmb->cks, ku=pmb->cke;
    if (pmb->block_size.nx3>1) { // 3D
      if(porb->orbital_direction == 1) {
        if(nb == 0) {
          jl += -xgh-porb->max_off_coarse+onx/2;
        } else if (nb == 4) {
          ju += 1+xgh-porb->min_off_coarse-onx/2;
        }
      } else if(porb->orbital_direction == 2) {
        onx = pmb->block_size.nx3;
        if(nb == 0) {
          kl += -xgh-porb->max_off_coarse+onx/2;
        } else if (nb == 4) {
          ku += 1+xgh-porb->min_off_coarse-onx/2;
        }
      }

      // TODO(tomo-ono): This part has a problem with "#pragma omp simd"
      //                 when using the Intel compiler
      // BufferUtility::UnpackData(buf, bco.x1f,
      //                           il, iu+1, jl-1, ju+1, kl-1, k+1u, p);
      // BufferUtility::UnpackData(buf, bco.x2f,
      //                           il-1, iu+1, jl, ju+1, kl-1, ku+1, p);
      // BufferUtility::UnpackData(buf, bco.x3f,
      //                           il-1, iu+1, jl-1, ju+1, kl, ku+1, p);
      // b1
      for (int k=kl-1; k<=ku+1; k++) {
        for (int j=jl-1; j<=ju+1; j++) {
          for (int i=il; i<=iu+1; i++) {
            bco.x1f(k,j,i) = buf[p++];
          }
        }
      }

      // b2
      for (int k=kl-1; k<=ku+1; k++) {
        for (int j=jl; j<=ju+1; j++) {
          for (int i=il-1; i<=iu+1; i++) {
            bco.x2f(k,j,i) = buf[p++];
          }
        }
      }

      // b3
      for (int k=kl; k<=ku+1; k++) {
        for (int j=jl-1; j<=ju+1; j++) {
          for (int i=il-1; i<=iu+1; i++) {
            bco.x3f(k,j,i) = buf[p++];
          }
        }
      }
    } else { // 2D
      if(porb->orbital_direction == 1) {
        onx = pmb->block_size.nx2;
        if(nb == 0) {
          jl += -xgh-porb->max_off_coarse+onx/2;
        } else if (nb == 4) {
          ju += 1+xgh-porb->min_off_coarse-onx/2;
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::LoadFieldBufferToFiner" << std::endl
            << "2D Orbital Advection is not allowed in spherical_polar." << std::endl;
        ATHENA_ERROR(msg);
      }

      // TODO(tomo-ono): This part has a problem with "#pragma omp simd"
      //                 when using the Intel compiler
      // BufferUtility::UnpackData(buf, bco.x1f,
      //                           il, iu+1, jl-1, ju+1, kl, ku, p);
      // BufferUtility::UnpackData(buf, bco.x2f,
      //                           il-1, iu+1, jl, ju+1, kl, ku, p);
      // BufferUtility::UnpackData(buf, bco.x3f,
      //                           il-1, iu+1, jl-1, ju+1, kl, ku, p);
      // b1
      for (int j=jl-1; j<=ju+1; j++) {
        for (int i=il; i<=iu+1; i++) {
          bco.x1f(kl,j,i) = buf[p++];
        }
      }

      // b2
      for (int j=jl; j<=ju+1; j++) {
        for (int i=il-1; i<=iu+1; i++) {
          bco.x2f(kl,j,i) = buf[p++];
        }
      }

      // b3
      for (int j=jl-1; j<=ju+1; j++) {
        for (int i=il-1; i<=iu+1; i++) {
          bco.x3f(kl,j,i) = buf[p++];
        }
      }
    }
    pmb->pmr->ProlongateSharedFieldX1(bco.x1f, bto.x1f, il, iu+1, jl, ju, kl, ku);
    pmb->pmr->ProlongateSharedFieldX2(bco.x2f, bto.x2f, il, iu, jl, ju+1, kl, ku);
    pmb->pmr->ProlongateSharedFieldX3(bco.x3f, bto.x3f, il, iu, jl, ju, kl, ku+1);
    pmb->pmr->ProlongateInternalField(bto, il, iu, jl, ju, kl, ku);

    if(porb->orbital_direction == 1) {
      if (nb == 0) {
        for (int k=pmb->ks; k<=pmb->ke; k++) {
          for (int i=pmb->is; i<=pmb->ie+1; i++) {
            int offset = porb->off[0](k,i);
            int xl = pmb->js-xgh-offset+onx;
            int xu = pmb->je;
            const int shift = (offset>0)? 0: -onx;
            for (int j=xl; j<=xu; j++) {
              bo1(k,i,j+shift) = bto.x1f(k,j,i);
            }
          }
        }
        if (pmb->block_size.nx3>1) { // 3D
          for (int k=pmb->ks; k<=pmb->ke+1; k++) {
            for (int i=pmb->is; i<=pmb->ie; i++) {
              int offset = porb->off[1](k,i);
              int xl = pmb->js-xgh-offset+onx;
              int xu = pmb->je;
              const int shift = (offset>0)? 0: -onx;
              for (int j=xl; j<=xu; j++) {
                bo2(k,i,j+shift) = -bto.x3f(k,j,i);
              }
            }
          }
        } else { // 2D
          int k = pmb->ks;
          for (int i=pmb->is; i<=pmb->ie; i++) {
            int offset = porb->off[1](k,i);
            int xl = pmb->js-xgh-offset+onx;
            int xu = pmb->je;
            const int shift = (offset>0)? 0: -onx;
            for (int j=xl; j<=xu; j++) {
              bo2(k,i,j+shift)   = -bto.x3f(k,j,i);
              bo2(k+1,i,j+shift) = -bto.x3f(k,j,i);
            }
          }
        }
      } else if (nb == 4) {
        for (int k=pmb->ks; k<=pmb->ke; k++) {
          for (int i=pmb->is; i<=pmb->ie+1; i++) {
            int offset = porb->off[0](k,i);
            int xl = pmb->js;
            int xu = pmb->je+1+xgh-offset-onx;
            const int shift = (offset>0)? 2*onx: onx;
            for (int j=xl; j<=xu; j++) {
              bo1(k,i,j+shift) = bto.x1f(k,j,i);
            }
          }
        }
        if (pmb->block_size.nx3>1) { // 3D
          for (int k=pmb->ks; k<=pmb->ke+1; k++) {
            for (int i=pmb->is; i<=pmb->ie; i++) {
              int offset = porb->off[1](k,i);
              int xl = pmb->js;
              int xu = pmb->je+1+xgh-offset-onx;
              const int shift = (offset>0)? 2*onx: onx;
              for (int j=xl; j<=xu; j++) {
                bo2(k,i,j+shift) = -bto.x3f(k,j,i);
              }
            }
          }
        } else { // 2D
          int k = pmb->ks;
          for (int i=pmb->is; i<=pmb->ie; i++) {
            int offset = porb->off[1](k,i);
            int xl = pmb->js;
            int xu = pmb->je+1+xgh-offset-onx;
            const int shift = (offset>0)? 2*onx: onx;
            for (int j=xl; j<=xu; j++) {
              bo2(k,i,j+shift)   = -bto.x3f(k,j,i);
              bo2(k+1,i,j+shift) = -bto.x3f(k,j,i);
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::SetFieldBufferFromCoarser" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    } else if(porb->orbital_direction == 2) {
      if (nb == 0) {
        for (int j=pmb->js; j<=pmb->je; j++) {
          for (int i=pmb->is; i<=pmb->ie+1; i++) {
            int offset = porb->off[0](j,i);
            int xl = pmb->ks-xgh-offset+onx;
            int xu = pmb->ke;
            const int shift = (offset>0)? 0: -onx;
            for (int k=xl; k<=xu; k++) {
              bo1(j,i,k+shift) = -bto.x1f(k,j,i);
            }
          }
        }
        for (int j=pmb->js; j<=pmb->je+1; j++) {
          for (int i=pmb->is; i<=pmb->ie; i++) {
            int offset = porb->off[1](j,i);
            int xl = pmb->ks-xgh-offset+onx;
            int xu = pmb->ke;
            const int shift = (offset>0)? 0: -onx;
            for (int k=xl; k<=xu; k++) {
              bo2(j,i,k+shift) = bto.x2f(k,j,i);
            }
          }
        }
      } else if (nb == 4) {
        for (int j=pmb->js; j<=pmb->je; j++) {
          for (int i=pmb->is; i<=pmb->ie+1; i++) {
            int offset = porb->off[0](j,i);
            int xl = pmb->ks;
            int xu = pmb->ke+1+xgh-offset-onx;
            const int shift = (offset>0)? 2*onx: onx;
            for (int k=xl; k<=xu; k++) {
              bo1(j,i,k+shift) = -bto.x1f(k,j,i);
            }
          }
        }
        for (int j=pmb->js; j<=pmb->je+1; j++) {
          for (int i=pmb->is; i<=pmb->ie; i++) {
            int offset = porb->off[1](j,i);
            int xl = pmb->ks;
            int xu = pmb->ke+1+xgh-offset-onx;
            const int shift = (offset>0)? 2*onx: onx;
            for (int k=xl; k<=xu; k++) {
              bo2(j,i,k+shift) = bto.x2f(k,j,i);
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::SetFieldBufferFromCoarser" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OrbitalBoundaryCommunication::SetFieldBufferFromFiner(
//!                                Real *buf, int &p, const int nb)
//! \brief set magnetic fields (finer level)
void OrbitalBoundaryCommunication::SetFieldBufferFromFiner(
                           Real *buf, int &p, const int nb) {
  MeshBlock *pmb = pmy_block_;
  OrbitalAdvection *porb = pmy_orbital_;
  AthenaArray<Real> &bo1 = porb->orbital_b1;
  AthenaArray<Real> &bo2 = porb->orbital_b2;

  if(porb->orbital_uniform_mesh) { // uniform mesh
    if(porb->orbital_direction == 1) {
      int &onx = pmb->block_size.nx2;
      int hnx1 = pmb->block_size.nx1/2;
      int hnx3 = pmb->block_size.nx3/2;
      int il, iu, kl, ku;
      if (pmb->block_size.nx3>1) {
        if (nb%4<2) {
          kl = pmb->ks;
          ku = pmb->ke-hnx3;
        } else {
          kl = pmb->ks+hnx3;
          ku = pmb->ke;
        }
      } else {
        if (nb%4>1) {
          std::stringstream msg;
          msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
              << "::SetFieldBufferFromFiner" << std::endl
              << "Neighbors are read incorrectly." << std::endl;
          ATHENA_ERROR(msg);
        }
        kl = pmb->ks;
        ku = pmb->ke;
      }
      if (nb%2==0) {
        il = pmb->is;
        iu = pmb->ie-hnx1;
      } else {
        il = pmb->is+hnx1;
        iu = pmb->ie;
      }
      if (nb<4) {
        for(int k=kl; k<=ku  ; k++) {
          for(int i=il; i<=iu+1; i++) {
            int offset = porb->off[0](k,i);
            int xl = pmb->js-xgh-offset+onx;
            int xu = pmb->je;
            if (offset<=0) {
              xl -= onx; xu -= onx;
            }
#pragma omp simd
            for (int j=xl; j<=xu; j++) {
              bo1(k,i,j) = buf[p++];
            }
          }
        }
        for(int k=kl; k<=ku+1; k++) {
          for(int i=il; i<=iu  ; i++) {
            int offset = porb->off[1](k,i);
            int xl = pmb->js-xgh-offset+onx;
            int xu = pmb->je;
            if (offset<=0) {
              xl -= onx; xu -= onx;
            }
#pragma omp simd
            for (int j=xl; j<=xu; j++) {
              bo2(k,i,j) = buf[p++];
            }
          }
        }
      } else if (nb<8) {
        for(int k=kl; k<=ku  ; k++) {
          for(int i=il; i<=iu+1; i++) {
            int offset = porb->off[0](k,i);
            int xl = pmb->js+onx;
            int xu = pmb->je+1+xgh-offset;
            if (offset>0) {
              xl += onx; xu += onx;
            }
#pragma omp simd
            for (int j=xl; j<=xu; j++) {
              bo1(k,i,j) = buf[p++];
            }
          }
        }
        for(int k=kl; k<=ku+1; k++) {
          for(int i=il; i<=iu  ; i++) {
            int offset = porb->off[1](k,i);
            int xl = pmb->js+onx;
            int xu = pmb->je+1+xgh-offset;
            if (offset>0) {
              xl += onx; xu += onx;
            }
#pragma omp simd
            for (int j=xl; j<=xu; j++) {
              bo2(k,i,j) = buf[p++];
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::SetFieldBufferFromFiner" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    } else if(porb->orbital_direction == 2) {
      int &onx = pmb->block_size.nx3;
      int hnx1 = pmb->block_size.nx1/2;
      int hnx2 = pmb->block_size.nx2/2;
      int il, iu, jl, ju;
      if (nb%4<2) {
        jl = pmb->js;
        ju = pmb->je-hnx2;
      } else {
        jl = pmb->js+hnx2;
        ju = pmb->je;
      }
      if (nb%2==0) {
        il = pmb->is;
        iu = pmb->ie-hnx1;
      } else {
        il = pmb->is+hnx1;
        iu = pmb->ie;
      }
      if (nb<4) {
        for(int j=jl; j<=ju  ; j++) {
          for(int i=il; i<=iu+1; i++) {
            int offset = porb->off[0](j,i);
            int xl = pmb->ks-xgh-offset+onx;
            int xu = pmb->ke;
            if (offset<=0) {
              xl -= onx; xu -= onx;
            }
#pragma omp simd
            for (int k=xl; k<=xu; k++) {
              bo1(j,i,k) = buf[p++];
            }
          }
        }
        for(int j=jl; j<=ju+1; j++) {
          for(int i=il; i<=iu  ; i++) {
            int offset = porb->off[1](j,i);
            int xl = pmb->ks-xgh-offset+onx;
            int xu = pmb->ke;
            if (offset<=0) {
              xl -= onx; xu -= onx;
            }
#pragma omp simd
            for (int k=xl; k<=xu; k++) {
              bo2(j,i,k) = buf[p++];
            }
          }
        }
      } else if (nb<8) {
        for(int j=jl; j<=ju  ; j++) {
          for(int i=il; i<=iu+1; i++) {
            int offset = porb->off[0](j,i);
            int xl = pmb->ks+onx;
            int xu = pmb->ke+1+xgh-offset;
            if (offset>0) {
              xl += onx; xu += onx;
            }
#pragma omp simd
            for (int k=xl; k<=xu; k++) {
              bo1(j,i,k) = buf[p++];
            }
          }
        }
        for(int j=jl; j<=ju+1; j++) {
          for(int i=il; i<=iu  ; i++) {
            int offset = porb->off[1](j,i);
            int xl = pmb->ks+onx;
            int xu = pmb->ke+1+xgh-offset;
            if (offset>0) {
              xl += onx; xu += onx;
            }
#pragma omp simd
            for (int k=xl; k<=xu; k++) {
              bo2(j,i,k) = buf[p++];
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::SetFieldBufferFromFiner" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void OrbitalBoundaryCommunication::LoadScalarBufferSameLevel(
//!                                        Real *buf, int &p, int nb)
//! \brief load passive scalars (same level)
void OrbitalBoundaryCommunication::LoadScalarBufferSameLevel(
                                   Real *buf, int &p, int nb) {
  MeshBlock *pmb = pmy_block_;
  OrbitalAdvection *porb = pmy_orbital_;
  AthenaArray<Real> &si = pmb->pscalars->s;

  if(porb->orbital_uniform_mesh) { // uniform mesh
    if(porb->orbital_direction == 1) {
      int &onx = pmb->block_size.nx2;
      if (nb==0) {
        for(int k=pmb->ks; k<=pmb->ke; k++) {
          for(int i=pmb->is; i<=pmb->ie; i++) {
            int offset = porb->ofc(k,i);
            int xl = pmb->js-xgh-offset+onx;
            for(int nsc=0 ; nsc<NSCALARS; nsc++) {
              for (int j=xl; j<=pmb->je; j++) {
                buf[p++] = si(nsc,k,j,i);
              }
            }
          }
        }
      } else if (nb==4) {
        for(int k=pmb->ks; k<=pmb->ke; k++) {
          for(int i=pmb->is; i<=pmb->ie; i++) {
            int offset = porb->ofc(k,i);
            int xu = pmb->je+1+xgh-offset-onx;
            for(int nsc=0 ; nsc<NSCALARS; nsc++) {
              for (int j=pmb->js; j<=xu; j++) {
                buf[p++] = si(nsc,k,j,i);
              }
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::LoadScalarBufferSameLevel" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    } else if (porb->orbital_direction == 2) {
      int &onx = pmb->block_size.nx3;
      if (nb==0) {
        for(int j=pmb->js; j<=pmb->je; j++) {
          for(int i=pmb->is; i<=pmb->ie; i++) {
            int offset = porb->ofc(j,i);
            int xl = pmb->ks-xgh-offset+onx;
            for(int nsc=0 ; nsc<NSCALARS; nsc++) {
              for (int k=xl; k<=pmb->ke; k++) {
                buf[p++] = si(nsc,k,j,i);
              }
            }
          }
        }
      } else if (nb==4) {
        for(int j=pmb->js; j<=pmb->je; j++) {
          for(int i=pmb->is; i<=pmb->ie; i++) {
            int offset = porb->ofc(j,i);
            int xu = pmb->ke+1+xgh-offset-onx;
            for(int nsc=0 ; nsc<NSCALARS; nsc++) {
              for (int k=pmb->ks; k<=xu; k++) {
                buf[p++] = si(nsc,k,j,i);
              }
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::LoadScalarBufferSameLevel" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    }
  }
//  else { // non-uniform mesh
//  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OrbitalBoundaryCommunication::LoadScalarBufferToCoarser(
//!                                        Real *buf, int &p, int nb)
//! \brief load passive scalars (coarser level)
void OrbitalBoundaryCommunication::LoadScalarBufferToCoarser(Real *buf, int &p, int nb) {
  MeshBlock *pmb = pmy_block_;
  OrbitalAdvection *porb = pmy_orbital_;
  AthenaArray<Real> &si = porb->s_coarse_send;

  if(porb->orbital_uniform_mesh) { // uniform mesh
    if(porb->orbital_direction == 1) {
      int honx = pmb->block_size.nx2/2;
      if (nb==0) {
        for(int k=pmb->cks; k<=pmb->cke; k++) {
          for(int i=pmb->cis; i<=pmb->cie; i++) {
            int offset = porb->ofc_coarse(k,i);
            int xl = pmb->cjs-xgh-offset+honx;
            for(int nsc=0 ; nsc<NSCALARS; nsc++) {
              for (int j=xl; j<=pmb->cje; j++) {
                buf[p++] = si(nsc,k,j,i);
              }
            }
          }
        }
      } else if (nb==4) {
        for(int k=pmb->cks; k<=pmb->cke; k++) {
          for(int i=pmb->cis; i<=pmb->cie; i++) {
            int offset = porb->ofc_coarse(k,i);
            int xu = pmb->cje+1+xgh-offset-honx;
            for(int nsc=0 ; nsc<NSCALARS; nsc++) {
              for (int j=pmb->cjs; j<=xu; j++) {
                buf[p++] = si(nsc,k,j,i);
              }
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::LoadScalarBufferToCoarser" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    } else if (porb->orbital_direction == 2) {
      int honx = pmb->block_size.nx3/2;
      if (nb==0) {
        for(int j=pmb->cjs; j<=pmb->cje; j++) {
          for(int i=pmb->cis; i<=pmb->cie; i++) {
            int offset = porb->ofc_coarse(j,i);
            int xl = pmb->cks-xgh-offset+honx;
            for(int nsc=0 ; nsc<NSCALARS; nsc++) {
              for (int k=xl; k<=pmb->cke; k++) {
                buf[p++] = si(nsc,k,j,i);
              }
            }
          }
        }
      } else if (nb==4) {
        for(int j=pmb->cjs; j<=pmb->cje; j++) {
          for(int i=pmb->cis; i<=pmb->cie; i++) {
            int offset = porb->ofc_coarse(j,i);
            int xu = pmb->cke+1+xgh-offset-honx;
            for(int nsc=0 ; nsc<NSCALARS; nsc++) {
              for (int k=pmb->cks; k<=xu; k++) {
                buf[p++] = si(nsc,k,j,i);
              }
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::LoadScalarBufferToCoarser" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    }
  }
//  else { // non-uniform mesh
//  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OrbitalBoundaryCommunication::LoadScalarBufferToFiner(
//!                                      Real *buf, int &p, int nb)
//! \brief load passive scalars (finer level)
void OrbitalBoundaryCommunication::LoadScalarBufferToFiner(Real *buf, int &p, int nb) {
  MeshBlock *pmb = pmy_block_;
  OrbitalAdvection *porb = pmy_orbital_;
  AthenaArray<Real> &si = pmb->pscalars->s;

  if(porb->orbital_uniform_mesh) { // uniform mesh
    if(porb->orbital_direction == 1) {
      int &onx = pmb->block_size.nx2;
      int il, iu, jl, ju, kl, ku;
      if (pmb->block_size.nx3>1) {
        if(nb%4<2) {
          kl = pmb->ks-1;
          ku = pmb->ke+1-pmb->block_size.nx3/2;
        } else {
          kl = pmb->ks-1+pmb->block_size.nx3/2;
          ku = pmb->ke+1;
        }
      } else {
          kl = pmb->ks;
          ku = pmb->ke;
      }
      if (nb%2==0) {
        il = pmb->is-1;
        iu = pmb->ie+1-pmb->block_size.nx1/2;
      } else {
        il = pmb->is-1+pmb->block_size.nx1/2;
        iu = pmb->ie+1;
      }
      if (nb<4) {
        jl = pmb->js-size_cc_send[0][2+nb]+1+onx;
        ju = pmb->je+1;
      } else if (nb<8) {
        jl = pmb->js-1;
        ju = pmb->je+size_cc_send[1][2+nb-4]-onx-1;
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::LoadScalarBufferToFiner" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
      BufferUtility::PackData(si, buf, 0, NSCALARS-1,
                              il, iu, jl, ju, kl, ku, p);
    } else if (porb->orbital_direction == 2) {
      int &onx = pmb->block_size.nx3;
      int il, iu, jl, ju, kl, ku;
      if (nb%4<2) {
        jl = pmb->js-1;
        ju = pmb->je+1-pmb->block_size.nx2/2;
      } else {
        jl = pmb->js-1+pmb->block_size.nx2/2;
        ju = pmb->je+1;
      }
      if (nb%2==0) {
        il = pmb->is-1;
        iu = pmb->ie+1-pmb->block_size.nx1/2;
      } else {
        il = pmb->is-1+pmb->block_size.nx1/2;
        iu = pmb->ie+1;
      }
      if (nb<4) {
        kl = pmb->ks-size_cc_send[0][2+nb]+1+onx;
        ku = pmb->ke+1;
      } else if (nb<8) {
        kl = pmb->ks-1;
        ku = pmb->ke+size_cc_send[0][2+nb-4]-onx-1;
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::LoadScalarBufferToFiner" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
      BufferUtility::PackData(si, buf, 0, NSCALARS-1,
                              il, iu, jl, ju, kl, ku, p);
    }
  }
//  else { // non-uniform mesh
//  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OrbitalBoundaryCommunication::SetScalarBufferSameLevel(
//!                                 Real *buf, int &p, const int nb)
//! \brief set passive scalars (same level)
void OrbitalBoundaryCommunication::SetScalarBufferSameLevel(
                            Real *buf, int &p, const int nb) {
  MeshBlock *pmb = pmy_block_;
  OrbitalAdvection *porb = pmy_orbital_;
  AthenaArray<Real> &so = porb->orbital_scalar;

  if(porb->orbital_uniform_mesh) { // uniform mesh
    if(porb->orbital_direction == 1) {
      int &onx = pmb->block_size.nx2;
      if(nb==0) {
        for(int k=pmb->ks; k<=pmb->ke; k++) {
          for(int i=pmb->is; i<=pmb->ie; i++) {
            int offset = porb->ofc(k,i);
            int xl = pmb->js-xgh-offset+onx;
            int xu = pmb->je;
            if (offset<=0) {
              xl -= onx; xu -= onx;
            }
            for(int nsc=0 ; nsc<NSCALARS; nsc++) {
#pragma omp simd
              for(int j=xl; j<=xu; j++) {
                so(nsc,k,i,j) = buf[p++];
              }
            }
          }
        }
      } else if(nb==4) {
        for(int k=pmb->ks; k<=pmb->ke; k++) {
          for(int i=pmb->is; i<=pmb->ie; i++) {
            int offset = porb->ofc(k,i);
            int xl = pmb->js+onx;
            int xu = pmb->je+1+xgh-offset;
            if (offset>0) {
              xl += onx; xu += onx;
            }
            for(int nsc=0 ; nsc<NSCALARS; nsc++) {
#pragma omp simd
              for(int j=xl; j<=xu; j++) {
                so(nsc,k,i,j) = buf[p++];
              }
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::SetScalarBufferSameLevel" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    } else if(porb->orbital_direction == 2) {
      int &onx = pmb->block_size.nx3;
      if (nb==0) {
        for(int j=pmb->js; j<=pmb->je; j++) {
          for(int i=pmb->is; i<=pmb->ie; i++) {
            int offset = porb->ofc(j,i);
            int xl = pmb->ks-xgh-offset+onx;
            int xu = pmb->ke;
            if (offset<=0) {
              xl -= onx; xu -= onx;
            }
            for(int nsc=0 ; nsc<NSCALARS; nsc++) {
#pragma omp simd
              for (int k=xl; k<=xu; k++) {
                so(nsc,j,i,k) = buf[p++];
              }
            }
          }
        }
      } else if(nb==4) {
        for(int j=pmb->js; j<=pmb->je; j++) {
          for(int i=pmb->is; i<=pmb->ie; i++) {
            int offset = porb->ofc(j,i);
            int xl = pmb->ks+onx;
            int xu = pmb->ke+1+xgh-offset;
            if(offset>0) {
              xl += onx; xu += onx;
            }
            for(int nsc=0 ; nsc<NSCALARS; nsc++) {
#pragma omp simd
              for(int k=xl; k<=xu; k++) {
                so(nsc,j,i,k) = buf[p++];
              }
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::SetScalarBufferSameLevel" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OrbitalBoundaryCommunication::SetScalarBufferFromCoarser(
//!                                   Real *buf, int &p, const int nb)
//! \brief set passive scalars (coarser level)
void OrbitalBoundaryCommunication::SetScalarBufferFromCoarser(
                              Real *buf, int &p, const int nb) {
  MeshBlock *pmb = pmy_block_;
  OrbitalAdvection *porb = pmy_orbital_;
  AthenaArray<Real> &so  = porb->orbital_scalar;
  AthenaArray<Real> &sco = porb->s_coarse_recv;
  AthenaArray<Real> &sto = porb->s_temp;

  if(porb->orbital_uniform_mesh) { // uniform mesh
    if(porb->orbital_direction == 1) {
      int &onx = pmb->block_size.nx2;
      int il, iu, kl, ku;
      if(pmb->block_size.nx3>1) {
        kl = pmb->cks-1;
        ku = pmb->cke+1;
      } else {
        kl = pmb->cks;
        ku = pmb->cke;
      }
      il = pmb->cis-1;
      iu = pmb->cie+1;
      if(nb==0) {
        int jl = pmb->cjs-xgh-porb->max_ofc_coarse+onx/2-1;
        int ju = pmb->cje+1;
        // TODO(tomo-ono): This part has a problem with "#pragma omp simd"
        //                 when using the Intel compiler
        // BufferUtility::UnpackData(buf, sco, 0, NSCALARS-1,
        //                           il, iu, jl, ju, kl, ku, p);
        for(int n=0; n<NSCALARS; ++n) {
          for(int k=kl; k<=ku; k++) {
            for(int j=jl; j<=ju; j++) {
              for(int i=il; i<=iu; i++) {
                sco(n,k,j,i) = buf[p++];
              }
            }
          }
        }
        pmb->pmr->ProlongateCellCenteredValues(sco, sto, 0, NSCALARS-1, pmb->cis,
                                               pmb->cie, jl+1, pmb->cje, pmb->cks,
                                               pmb->cke);
        for(int k=pmb->ks; k<=pmb->ke; k++) {
          for(int i=pmb->is; i<=pmb->ie; i++) {
            int offset = porb->ofc(k,i);
            int xl = pmb->js-xgh-offset+onx;
            int xu = pmb->je;
            const int shift = (offset>0)? 0: -onx;
            for(int nsc=0 ; nsc<NSCALARS; nsc++) {
              for(int j=xl; j<=xu; j++) {
                so(nsc,k,i,j+shift) = sto(nsc,k,j,i);
              }
            }
          }
        }
      } else if(nb==4) {
        int jl = pmb->cjs-1;
        int ju = pmb->cje+2+xgh-porb->min_ofc_coarse-onx/2;
        // TODO(tomo-ono): This part has a problem with "#pragma omp simd"
        //                 when using the Intel compiler
        // BufferUtility::UnpackData(buf, sco, 0, NSCALARS-1,
        //                           il, iu, jl, ju, kl, ku, p);
        for(int n=0; n<NSCALARS; ++n) {
          for(int k=kl; k<=ku; k++) {
            for(int j=jl; j<=ju; j++) {
              for(int i=il; i<=iu; i++) {
                sco(n,k,j,i) = buf[p++];
              }
            }
          }
        }
        pmb->pmr->ProlongateCellCenteredValues(sco, sto, 0, NSCALARS-1, pmb->cis,
                                               pmb->cie, pmb->cjs, ju-1, pmb->cks,
                                               pmb->cke);
        for(int k=pmb->ks; k<=pmb->ke; k++) {
          for(int i=pmb->is; i<=pmb->ie; i++) {
            int offset = porb->ofc(k,i);
            int xl = pmb->js;
            int xu = pmb->je+1+xgh-offset-onx;
            const int shift = (offset>0)? 2*onx: onx;
            for(int nsc=0 ; nsc<NSCALARS; nsc++) {
              for (int j=xl; j<=xu; j++) {
                so(nsc,k,i,j+shift) = sto(nsc,k,j,i);
              }
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::SetScalarBufferFromCoarser" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    } else if(porb->orbital_direction == 2) {
      int &onx = pmb->block_size.nx3;
      int il, iu, jl, ju;
      jl = pmb->cjs-1;
      ju = pmb->cje+1;
      il = pmb->cis-1;
      iu = pmb->cie+1;
      if(nb==0) {
        int kl = pmb->cks-xgh-porb->max_ofc_coarse+onx/2-1;
        int ku = pmb->cke+1;
        // TODO(tomo-ono): This part has a problem with "#pragma omp simd"
        //                 when using the Intel compiler
        // BufferUtility::UnpackData(buf, sco, 0, NSCALARS-1,
        //                           il, iu, jl, ju, kl, ku, p);
        for(int n=0; n<NSCALARS; ++n) {
          for(int k=kl; k<=ku; k++) {
            for(int j=jl; j<=ju; j++) {
              for(int i=il; i<=iu; i++) {
                sco(n,k,j,i) = buf[p++];
              }
            }
          }
        }
        pmb->pmr->ProlongateCellCenteredValues(sco, sto, 0, NSCALARS-1, pmb->cis,
                                               pmb->cie, pmb->cjs, pmb->cje,
                                               kl+1, pmb->cke);
        for(int j=pmb->js; j<=pmb->je; j++) {
          for(int i=pmb->is; i<=pmb->ie; i++) {
            int offset = porb->ofc(j,i);
            int xl = pmb->ks-xgh-offset+onx;
            int xu = pmb->ke;
            const int shift = (offset>0)? 0: -onx;
            for(int nsc=0 ; nsc<NSCALARS; nsc++) {
              for (int k=xl; k<=xu; k++) {
                so(nsc,j,i,k+shift) = sto(nsc,k,j,i);
              }
            }
          }
        }
      } else if (nb==4) {
        int kl = pmb->cks-1;
        int ku = pmb->cke+2+xgh-porb->min_ofc_coarse-onx/2;
        // TODO(tomo-ono): This part has a problem with "#pragma omp simd"
        //                 when using the Intel compiler
        // BufferUtility::UnpackData(buf, sco, 0, NSCALARS-1,
        //                           il, iu, jl, ju, kl, ku, p);
        for(int n=0; n<NSCALARS; ++n) {
          for(int k=kl; k<=ku; k++) {
            for(int j=jl; j<=ju; j++) {
              for(int i=il; i<=iu; i++) {
                sco(n,k,j,i) = buf[p++];
              }
            }
          }
        }
        pmb->pmr->ProlongateCellCenteredValues(sco, sto, 0, NSCALARS-1, pmb->cis,
                                               pmb->cie, pmb->cjs, pmb->cje,
                                               pmb->cks, ku-1);
        for(int j=pmb->js; j<=pmb->je; j++) {
          for(int i=pmb->is; i<=pmb->ie; i++) {
            int offset = porb->ofc(j,i);
            int xl = pmb->ks;
            int xu = pmb->ke+1+xgh-offset-onx;
            const int shift = (offset>0)? 2*onx: onx;
            for(int nsc=0 ; nsc<NSCALARS; nsc++) {
              for (int k=xl; k<=xu; k++) {
                so(nsc,j,i,k+shift) = sto(nsc,k,j,i);
              }
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::SetScalarBufferFromCoarser" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OrbitalBoundaryCommunication::SetScalarBufferFromFiner(
//!                                 Real *buf, int &p, const int nb)
//! \brief set passive scalars (passive level)
void OrbitalBoundaryCommunication::SetScalarBufferFromFiner(
                            Real *buf, int &p, const int nb) {
  MeshBlock *pmb = pmy_block_;
  OrbitalAdvection *porb = pmy_orbital_;
  AthenaArray<Real> &so  = porb->orbital_scalar;

  if(porb->orbital_uniform_mesh) { // uniform mesh
    if(porb->orbital_direction == 1) {
      int &onx = pmb->block_size.nx2;
      int hnx1 = pmb->block_size.nx1/2;
      int hnx3 = pmb->block_size.nx3/2;
      int il, iu, kl, ku;
      if (pmb->block_size.nx3>1)  {
        if (nb%4<2) {
          kl = pmb->ks;
          ku = pmb->ke-hnx3;
        } else {
          kl = pmb->ks+hnx3;
          ku = pmb->ke;
        }
      } else {
        if (nb%4>1) {
          std::stringstream msg;
          msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
              << "::SetScalarBufferFromFiner" << std::endl
              << "Neighbors are read incorrectly." << std::endl;
          ATHENA_ERROR(msg);
        }
        kl = pmb->ks;
        ku = pmb->ke;
      }
      if (nb%2==0) {
        il = pmb->is;
        iu = pmb->ie-hnx1;
      } else {
        il = pmb->is+hnx1;
        iu = pmb->ie;
      }
      if (nb<4) {
        for(int k=kl; k<=ku; k++) {
          for(int i=il; i<=iu; i++) {
            int offset = porb->ofc(k,i);
            int xl = pmb->js-xgh-offset+onx;
            int xu = pmb->je;
            if (offset<=0) {
              xl -= onx; xu -= onx;
            }
            for(int nsc=0 ; nsc<NSCALARS; nsc++) {
#pragma omp simd
              for (int j=xl; j<=xu; j++) {
                so(nsc,k,i,j) = buf[p++];
              }
            }
          }
        }
      } else if (nb<8) {
        for(int k=kl; k<=ku; k++) {
          for(int i=il; i<=iu; i++) {
            int offset = porb->ofc(k,i);
            int xl = pmb->js+onx;
            int xu = pmb->je+1+xgh-offset;
            if (offset>0) {
              xl += onx; xu += onx;
            }
            for(int nsc=0 ; nsc<NSCALARS; nsc++) {
#pragma omp simd
              for (int j=xl; j<=xu; j++) {
                so(nsc,k,i,j) = buf[p++];
              }
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::SetScalarBufferFromFiner" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    } else if(porb->orbital_direction == 2) {
      int &onx = pmb->block_size.nx3;
      int hnx1 = pmb->block_size.nx1/2;
      int hnx2 = pmb->block_size.nx2/2;
      int il, iu, jl, ju;
      if (nb%4<2) {
        jl = pmb->js;
        ju = pmb->je-hnx2;
      } else {
        jl = pmb->js+hnx2;
        ju = pmb->je;
      }
      if (nb%2==0) {
        il = pmb->is;
        iu = pmb->ie-hnx1;
      } else {
        il = pmb->is+hnx1;
        iu = pmb->ie;
      }
      if (nb<4) {
        for(int j=jl; j<=ju; j++) {
          for(int i=il; i<=iu; i++) {
            int offset = porb->ofc(j,i);
            int xl = pmb->ks-xgh-offset+onx;
            int xu = pmb->ke;
            if (offset<=0) {
              xl -= onx; xu -= onx;
            }
            for(int nsc=0 ; nsc<NSCALARS; nsc++) {
#pragma omp simd
              for(int k=xl; k<=xu; k++) {
                so(nsc,j,i,k) = buf[p++];
              }
            }
          }
        }
      } else if(nb<8) {
        for(int j=jl; j<=ju; j++) {
          for(int i=il; i<=iu; i++) {
            int offset = porb->ofc(j,i);
            int xl = pmb->ks+onx;
            int xu = pmb->ke+1+xgh-offset;
            if(offset>0) {
              xl += onx; xu += onx;
            }
            for(int nsc=0 ; nsc<NSCALARS; nsc++) {
#pragma omp simd
              for(int k=xl; k<=xu; k++) {
                so(nsc,j,i,k) = buf[p++];
              }
            }
          }
        }
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalBoundaryCommunication"
            << "::SetScalarBufferFromFiner" << std::endl
            << "Neighbors are read incorrectly." << std::endl;
        ATHENA_ERROR(msg);
      }
    }
  }
  return;
}
