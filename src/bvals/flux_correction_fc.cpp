//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file flux_correction_fc.cpp
//  \brief functions that perform flux correction for FACE_CENTERED variables

// C++ headers
#include <iostream>   // endl
#include <iomanip>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cstring>    // memcpy
#include <cstdlib>
#include <cmath>

// Athena++ classes headers
#include "bvals.hpp"
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../coordinates/coordinates.hpp"
#include "../parameter_input.hpp"
#include "../utils/buffer_utils.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadEMFBoundaryBufferSameLevel(Real *buf,
//                                                   const NeighborBlock& nb)
//  \brief Set EMF correction buffers for sending to a block on the same level
int BoundaryValues::LoadEMFBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb) {
  MeshBlock *pmb=pmy_block_;
  AthenaArray<Real> &e1=pmb->pfield->e.x1e;
  AthenaArray<Real> &e2=pmb->pfield->e.x2e;
  AthenaArray<Real> &e3=pmb->pfield->e.x3e;

  Real qomL = qshear_*Omega_0_*x1size_;
  AthenaArray<Real> &bx1=pmb->pfield->b.x1f;

  int p=0;
  if (nb.type==NEIGHBOR_FACE) {
    if (pmb->block_size.nx3 > 1) { // 3D
      // x1 direction
      if (nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
        int i;
        if (nb.fid==INNER_X1) i=pmb->is;
        else i=pmb->ie+1;
        // pack e2
        for (int k=pmb->ks; k<=pmb->ke+1; k++) {
          for (int j=pmb->js; j<=pmb->je; j++)
            buf[p++]=e2(k,j,i);
        }
        // pack e3
        // shift azmuthal velocity if shearing boundary blocks
        if (nb.shear && nb.fid==INNER_X1) {
          for (int k=pmb->ks; k<=pmb->ke; k++) {
            for (int j=pmb->js; j<=pmb->je+1; j++)
              buf[p++]=e3(k,j,i)-0.5*qomL*(bx1(k,j,i)+bx1(k,j-1,i));
          }
        } else if (nb.shear && nb.fid==OUTER_X1) {
          for (int k=pmb->ks; k<=pmb->ke; k++) {
            for (int j=pmb->js; j<=pmb->je+1; j++)
              buf[p++]=e3(k,j,i)+0.5*qomL*(bx1(k,j,i)+bx1(k,j-1,i));
          }
        } else {
          for (int k=pmb->ks; k<=pmb->ke; k++) {
            for (int j=pmb->js; j<=pmb->je+1; j++)
              buf[p++]=e3(k,j,i);
          }
        }
        // x2 direction
      } else if (nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
        int j;
        if (nb.fid==INNER_X2) j=pmb->js;
        else j=pmb->je+1;
        // pack e1
        for (int k=pmb->ks; k<=pmb->ke+1; k++) {
          for (int i=pmb->is; i<=pmb->ie; i++)
            buf[p++]=e1(k,j,i);
        }
        // pack e3
        for (int k=pmb->ks; k<=pmb->ke; k++) {
          for (int i=pmb->is; i<=pmb->ie+1; i++)
            buf[p++]=e3(k,j,i);
        }
      // x3 direction
      } else if (nb.fid==INNER_X3 || nb.fid==OUTER_X3) {
        int k;
        if (nb.fid==INNER_X3) k=pmb->ks;
        else k=pmb->ke+1;
        // pack e1
        for (int j=pmb->js; j<=pmb->je+1; j++) {
          for (int i=pmb->is; i<=pmb->ie; i++)
            buf[p++]=e1(k,j,i);
        }
        // pack e2
        for (int j=pmb->js; j<=pmb->je; j++) {
          for (int i=pmb->is; i<=pmb->ie+1; i++)
            buf[p++]=e2(k,j,i);
        }
      }
    } else if (pmb->block_size.nx2 > 1) { // 2D
      int k=pmb->ks;
      // x1 direction
      if (nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
        int i;
        if (nb.fid==INNER_X1) i=pmb->is;
        else i=pmb->ie+1;
        // pack e2
        // shift azimuthal velocity for x-z shearing
        if (SHEARING_BOX) {
          if (ShBoxCoord_==2 && (pmb->loc.lx1==0) && (nb.ox1==-1)) {
            for (int j=pmb->js; j<=pmb->je; j++)
              buf[p++]=e2(k,j,i)+qomL*bx1(k,j,i);
          } else if (ShBoxCoord_==2 && (pmb->loc.lx1==(pmb->pmy_mesh->nrbx1-1))
                                   && nb.ox1==1) {
            for (int j=pmb->js; j<=pmb->je; j++)
              buf[p++]=e2(k,j,i)-qomL*bx1(k,j,i);
          } else {
            for (int j=pmb->js; j<=pmb->je; j++)
              buf[p++]=e2(k,j,i);
          }
        } else {
          for (int j=pmb->js; j<=pmb->je; j++)
            buf[p++]=e2(k,j,i);
        }
        // pack e3
        for (int j=pmb->js; j<=pmb->je+1; j++)
          buf[p++]=e3(k,j,i);
      // x2 direction
      } else if (nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
        int j;
        if (nb.fid==INNER_X2) j=pmb->js;
        else j=pmb->je+1;
        // pack e1
        for (int i=pmb->is; i<=pmb->ie; i++)
          buf[p++]=e1(k,j,i);
        // pack e3
        for (int i=pmb->is; i<=pmb->ie+1; i++)
          buf[p++]=e3(k,j,i);
      }
    } else { // 1D
      int i, j=pmb->js, k=pmb->ks;
      if (nb.fid==INNER_X1) i=pmb->is;
      else i=pmb->ie+1;
      // pack e2 and e3
      buf[p++]=e2(k,j,i);
      buf[p++]=e3(k,j,i);
    }
  } else if (nb.type==NEIGHBOR_EDGE) {
    // x1x2 edge (both 2D and 3D)
    if (nb.eid>=0 && nb.eid<4) {
      int i, j;
      if ((nb.eid&1)==0) i=pmb->is;
      else i=pmb->ie+1;
      if ((nb.eid&2)==0) j=pmb->js;
      else j=pmb->je+1;
      // shift azmuthal velocity if shearing boundary blocks
      if (nb.shear and nb.ox1==-1) {
        for (int k=pmb->ks; k<=pmb->ke; k++)
          buf[p++]=e3(k,j,i)-0.5*qomL*(bx1(k,j,i)+bx1(k,j-1,i));
      } else if (nb.shear and nb.ox1==1) {
        for (int k=pmb->ks; k<=pmb->ke; k++)
          buf[p++]=e3(k,j,i)+0.5*qomL*(bx1(k,j,i)+bx1(k,j-1,i));
      } else {
      // pack e3
        for (int k=pmb->ks; k<=pmb->ke; k++)
          buf[p++]=e3(k,j,i);
      }
      // x1x3 edge
    } else if (nb.eid>=4 && nb.eid<8) {
      int i, k;
      if ((nb.eid&1)==0) i=pmb->is;
      else i=pmb->ie+1;
      if ((nb.eid&2)==0) k=pmb->ks;
      else k=pmb->ke+1;
      // pack e2
      // shift azimuthal velocity for x-z shearing
      if (SHEARING_BOX) {
        if (ShBoxCoord_==2 && (pmb->loc.lx1==0) && (nb.ox1==-1))   {
          for (int j=pmb->js; j<=pmb->je; j++)
            buf[p++]=e2(k,j,i)+qomL*bx1(k,j,i);
        } else if (ShBoxCoord_==2 && (pmb->loc.lx1==(pmb->pmy_mesh->nrbx1-1)) &&
                   nb.ox1==1) {
          for (int j=pmb->js; j<=pmb->je; j++)
            buf[p++]=e2(k,j,i)-qomL*bx1(k,j,i);
        } else {
          for (int j=pmb->js; j<=pmb->je; j++)
            buf[p++]=e2(k,j,i);
        }
      } else {
        for (int j=pmb->js; j<=pmb->je; j++)
          buf[p++]=e2(k,j,i);
      }
      // x2x3 edge
    } else if (nb.eid>=8 && nb.eid<12) {
      int j, k;
      if ((nb.eid&1)==0) j=pmb->js;
      else j=pmb->je+1;
      if ((nb.eid&2)==0) k=pmb->ks;
      else k=pmb->ke+1;
      // pack e1
      for (int i=pmb->is; i<=pmb->ie; i++)
        buf[p++]=e1(k,j,i);
    }
  }
  return p;
}


//----------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadEMFBoundaryBufferToCoarser(Real *buf,
//                                                         const NeighborBlock& nb)
//  \brief Set EMF correction buffers for sending to a block on the coarser level

int BoundaryValues::LoadEMFBoundaryBufferToCoarser(Real *buf, const NeighborBlock& nb) {
  MeshBlock *pmb=pmy_block_;
  Coordinates *pco=pmb->pcoord;
  AthenaArray<Real> &e1=pmb->pfield->e.x1e;
  AthenaArray<Real> &e2=pmb->pfield->e.x2e;
  AthenaArray<Real> &e3=pmb->pfield->e.x3e;
  // use the surface area aray as the edge length array
  AthenaArray<Real> &le1=sarea_[0];
  AthenaArray<Real> &le2=sarea_[1];
  int p=0;
  if (nb.type==NEIGHBOR_FACE) {
    if (pmb->block_size.nx3 > 1) { // 3D
      // x1 direction
      if (nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
        int i;
        if (nb.fid==INNER_X1) i=pmb->is;
        else i=pmb->ie+1;
        // restrict and pack e2
        for (int k=pmb->ks; k<=pmb->ke+1; k+=2) {
          for (int j=pmb->js; j<=pmb->je; j+=2) {
            Real el1=pco->GetEdge2Length(k,j,i);
            Real el2=pco->GetEdge2Length(k,j+1,i);
            buf[p++]=(e2(k,j,i)*el1+e2(k,j+1,i)*el2)/(el1+el2);
          }
        }
        // restrict and pack e3
        for (int k=pmb->ks; k<=pmb->ke; k+=2) {
          for (int j=pmb->js; j<=pmb->je+1; j+=2) {
            bool pole = pco->IsPole(j);
            Real el1, el2;
            if (not pole) {
              el1 = pco->GetEdge3Length(k,j,i);
              el2 = pco->GetEdge3Length(k+1,j,i);
            } else {
              el1 = pco->dx3f(k);
              el2 = pco->dx3f(k+1);
            }
            buf[p++]=(e3(k,j,i)*el1+e3(k+1,j,i)*el2)/(el1+el2);
          }
        }
        // x2 direction
      } else if (nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
        int j;
        if (nb.fid==INNER_X2) j=pmb->js;
        else j=pmb->je+1;
        bool pole = pco->IsPole(j);
        // restrict and pack e1
        for (int k=pmb->ks; k<=pmb->ke+1; k+=2) {
          if (not pole or not GENERAL_RELATIVITY) {
            pco->Edge1Length(k, j, pmb->is, pmb->ie, le1);
          } else {
            for (int i = pmb->is; i <= pmb->ie+1; i+=2) {
              le1(i) = pco->dx1f(i);
              le1(i+1) = pco->dx1f(i+1);
            }
          }
          for (int i=pmb->is; i<=pmb->ie; i+=2)
            buf[p++]=(e1(k,j,i)*le1(i)+e1(k,j,i+1)*le1(i+1))/(le1(i)+le1(i+1));
        }
        // restrict and pack e3
        for (int k=pmb->ks; k<=pmb->ke; k+=2) {
          if (not pole) {
            pco->Edge3Length(k,   j, pmb->is, pmb->ie+1, le1);
            pco->Edge3Length(k+1, j, pmb->is, pmb->ie+1, le2);
          } else {
            for (int i = pmb->is; i <= pmb->ie+1; i+=2) {
              le1(i) = pco->dx3f(k);
              le2(i) = pco->dx3f(k+1);
            }
          }
          for (int i=pmb->is; i<=pmb->ie+1; i+=2)
            buf[p++]=(e3(k,j,i)*le1(i)+e3(k+1,j,i)*le2(i))/(le1(i)+le2(i));
        }
      // x3 direction
      } else if (nb.fid==INNER_X3 || nb.fid==OUTER_X3) {
        int k;
        if (nb.fid==INNER_X3) k=pmb->ks;
        else k=pmb->ke+1;
        // restrict and pack e1
        for (int j=pmb->js; j<=pmb->je+1; j+=2) {
          bool pole = pco->IsPole(j);
          if (not pole or not GENERAL_RELATIVITY) {
            pco->Edge1Length(k, j, pmb->is, pmb->ie, le1);
          } else {
            for (int i = pmb->is; i <= pmb->ie; i+=2) {
              le1(i) = pco->dx1f(i);
              le1(i+1) = pco->dx1f(i+1);
            }
          }
          for (int i=pmb->is; i<=pmb->ie; i+=2)
            buf[p++]=(e1(k,j,i)*le1(i)+e1(k,j,i+1)*le1(i+1))/(le1(i)+le1(i+1));
        }
        // restrict and pack e2
        for (int j=pmb->js; j<=pmb->je; j+=2) {
          pco->Edge2Length(k,   j, pmb->is, pmb->ie+1, le1);
          pco->Edge2Length(k, j+1, pmb->is, pmb->ie+1, le2);
          for (int i=pmb->is; i<=pmb->ie+1; i+=2)
            buf[p++]=(e2(k,j,i)*le1(i)+e2(k,j+1,i)*le2(i))/(le1(i)+le2(i));
        }
      }
    } else if (pmb->block_size.nx2 > 1) { // 2D
      int k=pmb->ks;
      // x1 direction
      if (nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
        int i;
        if (nb.fid==INNER_X1) i=pmb->is;
        else i=pmb->ie+1;
        // restrict and pack e2
        for (int j=pmb->js; j<=pmb->je; j+=2) {
          Real el1=pco->GetEdge2Length(k,j,i);
          Real el2=pco->GetEdge2Length(k,j+1,i);
          buf[p++]=(e2(k,j,i)*el1+e2(k,j+1,i)*el2)/(el1+el2);
        }
        // pack e3
        for (int j=pmb->js; j<=pmb->je+1; j+=2)
          buf[p++]=e3(k,j,i);
        // x2 direction
      } else if (nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
        int j;
        if (nb.fid==INNER_X2) j=pmb->js;
        else j=pmb->je+1;
        bool pole = pco->IsPole(j);
        // restrict and pack e1
        if (not pole or not GENERAL_RELATIVITY) {
          pco->Edge1Length(k, j, pmb->is, pmb->ie, le1);
        } else {
          for (int i = pmb->is; i <= pmb->ie; i+=2) {
            le1(i) = pco->dx1f(i);
            le1(i+1) = pco->dx1f(i+1);
          }
        }
        for (int i=pmb->is; i<=pmb->ie; i+=2)
          buf[p++]=(e1(k,j,i)*le1(i)+e1(k,j,i+1)*le1(i+1))/(le1(i)+le1(i+1));
        // pack e3
        for (int i=pmb->is; i<=pmb->ie+1; i+=2)
          buf[p++]=e3(k,j,i);
      }
    } else { // 1D
      int i, j=pmb->js, k=pmb->ks;
      if (nb.fid==INNER_X1) i=pmb->is;
      else i=pmb->ie+1;
      // pack e2 and e3
      buf[p++]=e2(k,j,i);
      buf[p++]=e3(k,j,i);
    }
  } else if (nb.type==NEIGHBOR_EDGE) {
    if (pmb->block_size.nx3 > 1) { // 3D
      // x1x2 edge
      if (nb.eid>=0 && nb.eid<4) {
        int i, j;
        if ((nb.eid&1)==0) i=pmb->is;
        else i=pmb->ie+1;
        if ((nb.eid&2)==0) j=pmb->js;
        else j=pmb->je+1;
        bool pole = pco->IsPole(j);
        // restrict and pack e3
        for (int k=pmb->ks; k<=pmb->ke; k+=2) {
          Real el1, el2;
          if (not pole) {
            el1 = pco->GetEdge3Length(k,j,i);
            el2 = pco->GetEdge3Length(k+1,j,i);
          } else {
            el1 = pco->dx3f(k);
            el2 = pco->dx3f(k+1);
          }
          buf[p++]=(e3(k,j,i)*el1+e3(k+1,j,i)*el2)/(el1+el2);
        }
        // x1x3 edge
      } else if (nb.eid>=4 && nb.eid<8) {
        int i, k;
        if ((nb.eid&1)==0) i=pmb->is;
        else i=pmb->ie+1;
        if ((nb.eid&2)==0) k=pmb->ks;
        else k=pmb->ke+1;
        // restrict and pack e2
        for (int j=pmb->js; j<=pmb->je; j+=2) {
          Real el1=pco->GetEdge2Length(k,j,i);
          Real el2=pco->GetEdge2Length(k,j+1,i);
          buf[p++]=(e2(k,j,i)*el1+e2(k,j+1,i)*el2)/(el1+el2);
        }
        // x2x3 edge
      } else if (nb.eid>=8 && nb.eid<12) {
        int j, k;
        if ((nb.eid&1)==0) j=pmb->js;
        else j=pmb->je+1;
        bool pole = pco->IsPole(j);
        if ((nb.eid&2)==0) k=pmb->ks;
        else k=pmb->ke+1;
        // restrict and pack e1
        if (not pole or not GENERAL_RELATIVITY) {
          pco->Edge1Length(k, j, pmb->is, pmb->ie, le1);
        } else {
          for (int i = pmb->is; i <= pmb->ie; i+=2) {
            le1(i) = pco->dx1f(i);
            le1(i+1) = pco->dx1f(i+1);
          }
        }
        for (int i=pmb->is; i<=pmb->ie; i+=2)
          buf[p++]=(e1(k,j,i)*le1(i)+e1(k,j,i+1)*le1(i+1))/(le1(i)+le1(i+1));
      }
    } else if (pmb->block_size.nx2 > 1) { // 2D
      // x1x2 edge
      int i, j;
      if ((nb.eid&1)==0) i=pmb->is;
      else i=pmb->ie+1;
      if ((nb.eid&2)==0) j=pmb->js;
      else j=pmb->je+1;
      // pack e3
      buf[p++]=e3(pmb->ks,j,i);
    }
  }
  return p;
}

//----------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadEMFBoundaryPolarBuffer(Real *buf,
//          const PolarNeighborBlock &nb)
//  \brief Load EMF values along polar axis into send buffers

int BoundaryValues::LoadEMFBoundaryPolarBuffer(Real *buf, const PolarNeighborBlock &nb) {
  MeshBlock *pmb = pmy_block_;
  int count = 0;
  int j = nb.north ? pmb->js : pmb->je+1;
  for (int i = pmb->is; i <= pmb->ie; ++i) {
    Real val = 0.0;
    for (int k = pmb->ks; k <= pmb->ke; ++k) {  // avoid double counting right ends
      val += pmb->pfield->e.x1e(k, j, i);
    }
    buf[count++] = val / (pmb->ke - pmb->ks + 1);
  }
  return count;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SendEMFCorrection(void)
//  \brief Restrict, pack and send the surace EMF to the coarse neighbor(s) if
//  needed
void BoundaryValues::SendEMFCorrection(void) {
  MeshBlock *pmb=pmy_block_;

  // Send non-polar EMF values
  for (int n=0; n<nneighbor; n++) {
    NeighborBlock& nb = neighbor[n];
    if ((nb.type!=NEIGHBOR_FACE) && (nb.type!=NEIGHBOR_EDGE)) break;
    int p=0;
    if (nb.level==pmb->loc.level) {
      if ((nb.type==NEIGHBOR_FACE)
      || ((nb.type==NEIGHBOR_EDGE) && (edge_flag_[nb.eid]==true)))
        p=LoadEMFBoundaryBufferSameLevel(bd_emfcor_.send[nb.bufid], nb);
      else continue;
    } else if (nb.level==pmb->loc.level-1) {
      p=LoadEMFBoundaryBufferToCoarser(bd_emfcor_.send[nb.bufid], nb);
    } else {
      continue;
    }
    if (nb.rank==Globals::my_rank) { // on the same node
      MeshBlock *pbl=pmb->pmy_mesh->FindMeshBlock(nb.gid);
      std::memcpy(pbl->pbval->bd_emfcor_.recv[nb.targetid],
                  bd_emfcor_.send[nb.bufid], p*sizeof(Real));
      pbl->pbval->bd_emfcor_.flag[nb.targetid]=BNDRY_ARRIVED;
    }
#ifdef MPI_PARALLEL
    else
      MPI_Start(&(bd_emfcor_.req_send[nb.bufid]));
#endif
  }

  // Send polar EMF values
  for (int n = 0; n < num_north_polar_blocks_; ++n) {
    const PolarNeighborBlock &nb = polar_neighbor_north[n];
    int count = LoadEMFBoundaryPolarBuffer(emf_north_send_[n], nb);
    if (nb.rank == Globals::my_rank) { // on the same node
      MeshBlock *pbl = pmb->pmy_mesh->FindMeshBlock(nb.gid);
      std::memcpy(pbl->pbval->emf_north_recv_[pmb->loc.lx3],
          emf_north_send_[n], count * sizeof(Real));
      pbl->pbval->emf_north_flag_[pmb->loc.lx3] = BNDRY_ARRIVED;
    }
#ifdef MPI_PARALLEL
    else
      MPI_Start(&req_emf_north_send_[n]);
#endif
  }
  for (int n = 0; n < num_south_polar_blocks_; ++n) {
    const PolarNeighborBlock &nb = polar_neighbor_south[n];
    int count = LoadEMFBoundaryPolarBuffer(emf_south_send_[n], nb);
    if (nb.rank == Globals::my_rank) { // on the same node
      MeshBlock *pbl = pmb->pmy_mesh->FindMeshBlock(nb.gid);
      std::memcpy(pbl->pbval->emf_south_recv_[pmb->loc.lx3],
          emf_south_send_[n], count * sizeof(Real));
      pbl->pbval->emf_south_flag_[pmb->loc.lx3] = BNDRY_ARRIVED;
    }
#ifdef MPI_PARALLEL
    else
      MPI_Start(&req_emf_south_send_[n]);
#endif
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetEMFBoundarySameLevel(Real *buf, const NeighborBlock& nb)
//  \brief Add up the EMF received from a block on the same level
//         Later they will be divided in the AverageEMFBoundary function

void BoundaryValues::SetEMFBoundarySameLevel(Real *buf, const NeighborBlock& nb) {
  MeshBlock *pmb=pmy_block_;
  AthenaArray<Real> &e1=pmb->pfield->e.x1e;
  AthenaArray<Real> &e2=pmb->pfield->e.x2e;
  AthenaArray<Real> &e3=pmb->pfield->e.x3e;
  int p=0;
  if (nb.type==NEIGHBOR_FACE) {
    if (pmb->block_size.nx3 > 1) { // 3D
      // x1 direction
      if (nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
        int i;
        if (nb.fid==INNER_X1) i=pmb->is;
        else i=pmb->ie+1;

        if (nb.shear and nb.fid==INNER_X1) {
          // store e2 for shearing periodic bcs
          for (int k=pmb->ks; k<=pmb->ke+1; k++) {
            for (int j=pmb->js; j<=pmb->je; j++)
              shboxvar_inner_emf_.x2e(k,j) += buf[p++];
          }
          // store e3 for shearing periodic bcs
          for (int k=pmb->ks; k<=pmb->ke; k++) {
            for (int j=pmb->js; j<=pmb->je+1; j++)
              shboxvar_inner_emf_.x3e(k,j) += buf[p++];
          }
        } else if (nb.shear and nb.fid==OUTER_X1) {
          // store e2 for shearing periodic bcs
          for (int k=pmb->ks; k<=pmb->ke+1; k++) {
            for (int j=pmb->js; j<=pmb->je; j++)
              shboxvar_outer_emf_.x2e(k,j) += buf[p++];
          }
          // store e3 for shearing periodic bcs
          for (int k=pmb->ks; k<=pmb->ke; k++) {
            for (int j=pmb->js; j<=pmb->je+1; j++)
              shboxvar_outer_emf_.x3e(k,j) += buf[p++];
          }
        } else {
          // unpack e2
          for (int k=pmb->ks; k<=pmb->ke+1; k++) {
            for (int j=pmb->js; j<=pmb->je; j++)
              e2(k,j,i)+=buf[p++];
          }
          // unpack e3
          for (int k=pmb->ks; k<=pmb->ke; k++) {
            for (int j=pmb->js; j<=pmb->je+1; j++)
              e3(k,j,i)+=buf[p++];
          }
        }
        // x2 direction
      } else if (nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
        int j;
        if (nb.fid==INNER_X2) j=pmb->js;
        else j=pmb->je+1;
        // unpack e1
        Real sign = (nb.polar and flip_across_pole_field[IB1]) ? -1.0 : 1.0;
        for (int k=pmb->ks; k<=pmb->ke+1; k++) {
          for (int i=pmb->is; i<=pmb->ie; i++)
            e1(k,j,i)+=sign*buf[p++];
        }
        // unpack e3
        sign = (nb.polar and flip_across_pole_field[IB3]) ? -1.0 : 1.0;
        for (int k=pmb->ks; k<=pmb->ke; k++) {
          for (int i=pmb->is; i<=pmb->ie+1; i++)
            e3(k,j,i)+=sign*buf[p++];
        }
        // x3 direction
      } else if (nb.fid==INNER_X3 || nb.fid==OUTER_X3) {
        int k;
        if (nb.fid==INNER_X3) k=pmb->ks;
        else k=pmb->ke+1;
        // unpack e1
        for (int j=pmb->js; j<=pmb->je+1; j++) {
          for (int i=pmb->is; i<=pmb->ie; i++)
            e1(k,j,i)+=buf[p++];
        }
        // unpack e2
        for (int j=pmb->js; j<=pmb->je; j++) {
          for (int i=pmb->is; i<=pmb->ie+1; i++)
            e2(k,j,i)+=buf[p++];
        }
      }
    } else if (pmb->block_size.nx2 > 1) { // 2D
      int k=pmb->ks;
      // x1 direction
      if (nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
        int i;
        if (nb.fid==INNER_X1) i=pmb->is;
        else i=pmb->ie+1;
        // unpack e2
        for (int j=pmb->js; j<=pmb->je; j++) {
          e2(k+1,j,i)+=buf[p];
          e2(k,  j,i)+=buf[p++];
        }
        // unpack e3
        for (int j=pmb->js; j<=pmb->je+1; j++)
          e3(k,j,i)+=buf[p++];
        // x2 direction
      } else if (nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
        int j;
        if (nb.fid==INNER_X2) j=pmb->js;
        else j=pmb->je+1;
        // unpack e1
        for (int i=pmb->is; i<=pmb->ie; i++) {
          e1(k+1,j,i)+=buf[p];
          e1(k  ,j,i)+=buf[p++];
        }
        // unpack e3
        for (int i=pmb->is; i<=pmb->ie+1; i++)
          e3(k,j,i)+=buf[p++];
      }
    } else { // 1D
      int i, j=pmb->js, k=pmb->ks;
      if (nb.fid==INNER_X1) i=pmb->is;
      else i=pmb->ie+1;
      // unpack e2
      e2(k+1,j,i)+=buf[p];
      e2(k  ,j,i)+=buf[p++];
      // unpack e3
      e3(k,j+1,i)+=buf[p];
      e3(k  ,j,i)+=buf[p++];
    }
  } else if (nb.type==NEIGHBOR_EDGE) {
    // x1x2 edge (2D and 3D)
    if (nb.eid>=0 && nb.eid<4) {
      int i, j;
      if ((nb.eid&1)==0) i=pmb->is;
      else i=pmb->ie+1;
      if ((nb.eid&2)==0) j=pmb->js;
      else j=pmb->je+1;
      if (nb.shear and nb.ox1==-1) {
        // store e3 for shearing periodic bcs
        for (int k=pmb->ks; k<=pmb->ke; k++)
          shboxvar_inner_emf_.x3e(k,j) += buf[p++];
      } else if (nb.shear and nb.ox1==1) {
        // store e3 for shearing periodic bcs
        for (int k=pmb->ks; k<=pmb->ke; k++)
          shboxvar_outer_emf_.x3e(k,j) += buf[p++];
      } else {
      // unpack e3
        Real sign = (nb.polar and flip_across_pole_field[IB3]) ? -1.0 : 1.0;
        for (int k=pmb->ks; k<=pmb->ke; k++) {
          e3(k,j,i)+=sign*buf[p++];
        }
      }
      // x1x3 edge
    } else if (nb.eid>=4 && nb.eid<8) {
      int i, k;
      if ((nb.eid&1)==0) i=pmb->is;
      else i=pmb->ie+1;
      if ((nb.eid&2)==0) k=pmb->ks;
      else k=pmb->ke+1;
      if (nb.shear and nb.ox1==-1) {
        // store e2 for shearing periodic bcs
        for (int j=pmb->js; j<=pmb->je; j++)
          shboxvar_inner_emf_.x2e(k,j) += buf[p++];
      } else if (nb.shear and nb.ox1==1) {
        // store e2 for shearing periodic bcs
        for (int j=pmb->js; j<=pmb->je; j++)
          shboxvar_outer_emf_.x2e(k,j) += buf[p++];
      } else {
      // unpack e2
        for (int j=pmb->js; j<=pmb->je; j++)
          e2(k,j,i)+=buf[p++];
      }
      // x2x3 edge
    } else if (nb.eid>=8 && nb.eid<12) {
      int j, k;
      if ((nb.eid&1)==0) j=pmb->js;
      else j=pmb->je+1;
      if ((nb.eid&2)==0) k=pmb->ks;
      else k=pmb->ke+1;
      // unpack e1
      Real sign = (nb.polar and flip_across_pole_field[IB1]) ? -1.0 : 1.0;
      for (int i=pmb->is; i<=pmb->ie; i++)
        e1(k,j,i)+=sign*buf[p++];
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetEMFBoundaryFromFiner(Real *buf, const NeighborBlock& nb)
//  \brief Add up the EMF received from a block on the finer level
//         Later they will be divided in the AverageEMFBoundary function

void BoundaryValues::SetEMFBoundaryFromFiner(Real *buf, const NeighborBlock& nb) {
  MeshBlock *pmb=pmy_block_;
  AthenaArray<Real> &e1=pmb->pfield->e.x1e;
  AthenaArray<Real> &e2=pmb->pfield->e.x2e;
  AthenaArray<Real> &e3=pmb->pfield->e.x3e;
  int p=0;
  if (nb.type==NEIGHBOR_FACE) {
    if (pmb->block_size.nx3 > 1) { // 3D
      // x1 direction
      if (nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
        int i, jl=pmb->js, ju=pmb->je, kl=pmb->ks, ku=pmb->ke;
        if (nb.fid==INNER_X1) i=pmb->is;
        else i=pmb->ie+1;
        if (nb.fi1==0) ju=pmb->js+pmb->block_size.nx2/2-1;
        else jl=pmb->js+pmb->block_size.nx2/2;
        if (nb.fi2==0) ku=pmb->ks+pmb->block_size.nx3/2-1;
        else kl=pmb->ks+pmb->block_size.nx3/2;
        // unpack e2
        for (int k=kl; k<=ku+1; k++) {
          for (int j=jl; j<=ju; j++)
            e2(k,j,i)+=buf[p++];
        }
        // unpack e3
        for (int k=kl; k<=ku; k++) {
          for (int j=jl; j<=ju+1; j++)
            e3(k,j,i)+=buf[p++];
        }
        // x2 direction
      } else if (nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
        int j, il=pmb->is, iu=pmb->ie, kl=pmb->ks, ku=pmb->ke;
        if (nb.fid==INNER_X2) j=pmb->js;
        else j=pmb->je+1;
        if (nb.fi1==0) iu=pmb->is+pmb->block_size.nx1/2-1;
        else il=pmb->is+pmb->block_size.nx1/2;
        if (nb.fi2==0) ku=pmb->ks+pmb->block_size.nx3/2-1;
        else kl=pmb->ks+pmb->block_size.nx3/2;
        // unpack e1
        Real sign = (nb.polar and flip_across_pole_field[IB1]) ? -1.0 : 1.0;
        for (int k=kl; k<=ku+1; k++) {
          for (int i=il; i<=iu; i++)
            e1(k,j,i)+=sign*buf[p++];
        }
        // unpack e3
        sign = (nb.polar and flip_across_pole_field[IB3]) ? -1.0 : 1.0;
        for (int k=kl; k<=ku; k++) {
          for (int i=il; i<=iu+1; i++)
            e3(k,j,i)+=sign*buf[p++];
        }
      // x3 direction
      } else if (nb.fid==INNER_X3 || nb.fid==OUTER_X3) {
        int k, il=pmb->is, iu=pmb->ie, jl=pmb->js, ju=pmb->je;
        if (nb.fid==INNER_X3) k=pmb->ks;
        else k=pmb->ke+1;
        if (nb.fi1==0) iu=pmb->is+pmb->block_size.nx1/2-1;
        else il=pmb->is+pmb->block_size.nx1/2;
        if (nb.fi2==0) ju=pmb->js+pmb->block_size.nx2/2-1;
        else jl=pmb->js+pmb->block_size.nx2/2;
        // unpack e1
        for (int j=jl; j<=ju+1; j++) {
          for (int i=il; i<=iu; i++)
            e1(k,j,i)+=buf[p++];
        }
        // unpack e2
        for (int j=jl; j<=ju; j++) {
          for (int i=il; i<=iu+1; i++)
            e2(k,j,i)+=buf[p++];
        }
      }
    } else if (pmb->block_size.nx2 > 1) { // 2D
      int k=pmb->ks;
      // x1 direction
      if (nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
        int i, jl=pmb->js, ju=pmb->je;
        if (nb.fid==INNER_X1) i=pmb->is;
        else i=pmb->ie+1;
        if (nb.fi1==0) ju=pmb->js+pmb->block_size.nx2/2-1;
        else jl=pmb->js+pmb->block_size.nx2/2;
        // unpack e2
        for (int j=jl; j<=ju; j++) {
          e2(k+1,j,i)+=buf[p];
          e2(k,  j,i)+=buf[p++];
        }
        // unpack e3
        for (int j=jl; j<=ju+1; j++)
          e3(k,j,i)+=buf[p++];
        // x2 direction
      } else if (nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
        int j, il=pmb->is, iu=pmb->ie;
        if (nb.fid==INNER_X2) j=pmb->js;
        else j=pmb->je+1;
        if (nb.fi1==0) iu=pmb->is+pmb->block_size.nx1/2-1;
        else il=pmb->is+pmb->block_size.nx1/2;
        // unpack e1
        for (int i=il; i<=iu; i++) {
          e1(k+1,j,i)+=buf[p];
          e1(k  ,j,i)+=buf[p++];
        }
        // unpack e3
        for (int i=il; i<=iu+1; i++)
          e3(k,j,i)+=buf[p++];
      }
    } else { // 1D
      int i, j=pmb->js, k=pmb->ks;
      if (nb.fid==INNER_X1) i=pmb->is;
      else i=pmb->ie+1;
      // unpack e2
      e2(k+1,j,i)+=buf[p];
      e2(k  ,j,i)+=buf[p++];
      // unpack e3
      e3(k,j+1,i)+=buf[p];
      e3(k  ,j,i)+=buf[p++];
    }
  } else if (nb.type==NEIGHBOR_EDGE) {
    if (pmb->block_size.nx3 > 1) { // 3D
      // x1x2 edge
      if (nb.eid>=0 && nb.eid<4) {
        int i, j, kl=pmb->ks, ku=pmb->ke;
        if ((nb.eid&1)==0) i=pmb->is;
        else i=pmb->ie+1;
        if ((nb.eid&2)==0) j=pmb->js;
        else j=pmb->je+1;
        if (nb.fi1==0) ku=pmb->ks+pmb->block_size.nx3/2-1;
        else kl=pmb->ks+pmb->block_size.nx3/2;
        // unpack e3
        Real sign = (nb.polar and flip_across_pole_field[IB3]) ? -1.0 : 1.0;
        for (int k=kl; k<=ku; k++)
          e3(k,j,i)+=sign*buf[p++];
        // x1x3 edge
      } else if (nb.eid>=4 && nb.eid<8) {
        int i, k, jl=pmb->js, ju=pmb->je;
        if ((nb.eid&1)==0) i=pmb->is;
        else i=pmb->ie+1;
        if ((nb.eid&2)==0) k=pmb->ks;
        else k=pmb->ke+1;
        if (nb.fi1==0) ju=pmb->js+pmb->block_size.nx2/2-1;
        else jl=pmb->js+pmb->block_size.nx2/2;
        // unpack e2
        for (int j=jl; j<=ju; j++)
          e2(k,j,i)+=buf[p++];
        // x2x3 edge
      } else if (nb.eid>=8 && nb.eid<12) {
        int j, k, il=pmb->is, iu=pmb->ie;
        if ((nb.eid&1)==0) j=pmb->js;
        else j=pmb->je+1;
        if ((nb.eid&2)==0) k=pmb->ks;
        else k=pmb->ke+1;
        if (nb.fi1==0) iu=pmb->is+pmb->block_size.nx1/2-1;
        else il=pmb->is+pmb->block_size.nx1/2;
        // unpack e1
        Real sign = (nb.polar and flip_across_pole_field[IB1]) ? -1.0 : 1.0;
        for (int i=il; i<=iu; i++)
          e1(k,j,i)+=sign*buf[p++];
      }
    } else if (pmb->block_size.nx2 > 1) { // 2D
      int i, j, k=pmb->ks;
      if ((nb.eid&1)==0) i=pmb->is;
      else i=pmb->ie+1;
      if ((nb.eid&2)==0) j=pmb->js;
      else j=pmb->je+1;
      // unpack e3
      e3(k,j,i)+=buf[p++];
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetEMFBoundaryPolar(Real **buf_list, int num_bufs,
//          bool north)
//  \brief Overwrite EMF values along polar axis with azimuthal averages

void BoundaryValues::SetEMFBoundaryPolar(Real **buf_list, int num_bufs, bool north) {
  MeshBlock *pmb = pmy_block_;
  if (pmb->block_size.nx3 > 1) {
    int j = north ? pmb->js : pmb->je+1;
    int count = 0;
    for (int i = pmb->is; i <= pmb->ie; ++i) {
      Real val = 0.0;
      for (int n = 0; n < num_bufs; ++n)
        val += buf_list[n][count];
      for (int k = pmb->ks-NGHOST; k <= pmb->ke+NGHOST+1; ++k)
        pmb->pfield->e.x1e(k, j, i) = val / num_bufs;
      ++count;
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ClearCoarseEMFBoundary(void)
//  \brief Clear the EMFs on the surface/edge contacting with a finer block

void BoundaryValues::ClearCoarseEMFBoundary(void) {
  MeshBlock *pmb=pmy_block_;
  AthenaArray<Real> &e1=pmb->pfield->e.x1e;
  AthenaArray<Real> &e2=pmb->pfield->e.x2e;
  AthenaArray<Real> &e3=pmb->pfield->e.x3e;
  int i, j, k, nl;
  // face
  for (int n=0; n<nface_; n++) {
    if (n==INNER_X1 || n==OUTER_X1) {
      if (n==INNER_X1) i=pmb->is;
      else i=pmb->ie+1;
      nl=nblevel[1][1][2*n];
      if (nl>pmb->loc.level) { // finer
        if (pmb->block_size.nx3 > 1) { // 3D
          for (int k=pmb->ks+1; k<=pmb->ke; k++) {
            for (int j=pmb->js; j<=pmb->je; j++)
              e2(k,j,i)=0.0;
          }
          for (int k=pmb->ks; k<=pmb->ke; k++) {
            for (int j=pmb->js+1; j<=pmb->je; j++)
              e3(k,j,i)=0.0;
          }
        } else if (pmb->block_size.nx2 > 1) { // 2D
          for (int j=pmb->js; j<=pmb->je; j++)
            e2(pmb->ks,j,i)=e2(pmb->ks+1,j,i)=0.0;
          for (int j=pmb->js+1; j<=pmb->je; j++)
            e3(pmb->ks,j,i)=0.0;
        } else { // 1D
          e2(pmb->ks,pmb->js,i)=e2(pmb->ks+1,pmb->js,i)=0.0;
          e3(pmb->ks,pmb->js,i)=e3(pmb->ks,pmb->js+1,i)=0.0;
        }
      }
    }
    if (n==INNER_X2 || n==OUTER_X2) {
      if (n==INNER_X2) j=pmb->js;
      else j=pmb->je+1;
      nl=nblevel[1][2*n-4][1];
      if (nl>pmb->loc.level) { // finer
        if (pmb->block_size.nx3 > 1) { // 3D
          for (int k=pmb->ks+1; k<=pmb->ke; k++) {
            for (int i=pmb->is; i<=pmb->ie; i++)
              e1(k,j,i)=0.0;
          }
          for (int k=pmb->ks; k<=pmb->ke; k++) {
            for (int i=pmb->is+1; i<=pmb->ie; i++)
              e3(k,j,i)=0.0;
          }
        } else if (pmb->block_size.nx2 > 1) { // 2D
          for (int i=pmb->is; i<=pmb->ie; i++)
            e1(pmb->ks,j,i)=e1(pmb->ks+1,j,i)=0.0;
          for (int i=pmb->is+1; i<=pmb->ie; i++)
            e3(pmb->ks,j,i)=0.0;
        }
      }
    }
    if (n==INNER_X3 || n==OUTER_X3) {
      if (n==INNER_X3) k=pmb->ks;
      else k=pmb->ke+1;
      nl=nblevel[2*n-8][1][1];
      if (nl>pmb->loc.level) { // finer
        // this is always 3D
        for (int j=pmb->js+1; j<=pmb->je; j++) {
          for (int i=pmb->is; i<=pmb->ie; i++)
            e1(k,j,i)=0.0;
        }
        for (int j=pmb->js; j<=pmb->je; j++) {
          for (int i=pmb->is+1; i<=pmb->ie; i++)
            e2(k,j,i)=0.0;
        }
      }
    }
  }
  // edge
  for (int n=0; n<nedge_; n++) {
    if (edge_flag_[n]==true) continue;
    if (n>=0 && n<4) {
      if ((n&1)==0) i=pmb->is;
      else i=pmb->ie+1;
      if ((n&2)==0) j=pmb->js;
      else j=pmb->je+1;
      for (int k=pmb->ks; k<=pmb->ke; k++)
        e3(k,j,i)=0.0;
      // x1x3 edge
    } else if (n>=4 && n<8) {
      if ((n&1)==0) i=pmb->is;
      else i=pmb->ie+1;
      if ((n&2)==0) k=pmb->ks;
      else k=pmb->ke+1;
      for (int j=pmb->js; j<=pmb->je; j++)
        e2(k,j,i)=0.0;
      // x2x3 edge
    } else if (n>=8 && n<12) {
      if ((n&1)==0) j=pmb->js;
      else j=pmb->je+1;
      if ((n&2)==0) k=pmb->ks;
      else k=pmb->ke+1;
      for (int i=pmb->is; i<=pmb->ie; i++)
        e1(k,j,i)=0.0;
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::AverageEMFBoundary(void)
//  \brief Set EMF boundary received from a block on the finer level

void BoundaryValues::AverageEMFBoundary(void) {
  MeshBlock *pmb=pmy_block_;
  AthenaArray<Real> &e1=pmb->pfield->e.x1e;
  AthenaArray<Real> &e2=pmb->pfield->e.x2e;
  AthenaArray<Real> &e3=pmb->pfield->e.x3e;
  int i, j, k, nl;
  // face
  for (int n=0; n<nface_; n++) {
    if ((block_bcs[n] != BLOCK_BNDRY) && (block_bcs[n] != PERIODIC_BNDRY)
        && (block_bcs[n] != POLAR_BNDRY)) continue;
    if (n==INNER_X1 || n==OUTER_X1) {
      if (n==INNER_X1) i=pmb->is;
      else i=pmb->ie+1;
      nl=nblevel[1][1][2*n];
      if (nl==pmb->loc.level) { // same ; divide all the face EMFs by 2
        if (pmb->block_size.nx3 > 1) { // 3D
          for (int k=pmb->ks+1; k<=pmb->ke; k++) {
            for (int j=pmb->js; j<=pmb->je; j++)
              e2(k,j,i)*=0.5;
          }
          for (int k=pmb->ks; k<=pmb->ke; k++) {
            for (int j=pmb->js+1; j<=pmb->je; j++)
              e3(k,j,i)*=0.5;
          }
        } else if (pmb->block_size.nx2 > 1) { // 2D
          for (int j=pmb->js; j<=pmb->je; j++)
            e2(pmb->ks,j,i)*=0.5, e2(pmb->ks+1,j,i)*=0.5;
          for (int j=pmb->js+1; j<=pmb->je; j++)
            e3(pmb->ks,j,i)*=0.5;
        } else { // 1D
          e2(pmb->ks,pmb->js,i)*=0.5, e2(pmb->ks+1,pmb->js,i)*=0.5;
          e3(pmb->ks,pmb->js,i)*=0.5, e3(pmb->ks,pmb->js+1,i)*=0.5;
        }
      } else if (nl>pmb->loc.level) { // finer; divide the overlapping EMFs by 2
        if (pmb->block_size.nx3 > 1) { // 3D
          int k=pmb->ks+pmb->block_size.nx3/2;
          for (int j=pmb->js; j<=pmb->je; j++)
            e2(k,j,i)*=0.5;
        }
        if (pmb->block_size.nx2 > 1) { // 2D or 3D
          int j=pmb->js+pmb->block_size.nx2/2;
          for (int k=pmb->ks; k<=pmb->ke; k++)
            e3(k,j,i)*=0.5;
        }
      }
    }
    if (n==INNER_X2 || n==OUTER_X2) {
      if (n==INNER_X2) j=pmb->js;
      else j=pmb->je+1;
      nl=nblevel[1][2*n-4][1];
      if (nl==pmb->loc.level) { // same ; divide all the face EMFs by 2
        if (pmb->block_size.nx3 > 1) {
          for (int k=pmb->ks+1; k<=pmb->ke; k++) {
            for (int i=pmb->is; i<=pmb->ie; i++)
              e1(k,j,i)*=0.5;
          }
          for (int k=pmb->ks; k<=pmb->ke; k++) {
            for (int i=pmb->is+1; i<=pmb->ie; i++)
              e3(k,j,i)*=0.5;
          }
        } else if (pmb->block_size.nx2 > 1) {
          for (int i=pmb->is; i<=pmb->ie; i++)
            e1(pmb->ks,j,i)*=0.5, e1(pmb->ks+1,j,i)*=0.5;
          for (int i=pmb->is+1; i<=pmb->ie; i++)
            e3(pmb->ks,j,i)*=0.5;
        }
      } else if (nl>pmb->loc.level) { // finer; divide the overlapping EMFs by 2
        if (pmb->block_size.nx3 > 1) { // 3D
          int k=pmb->ks+pmb->block_size.nx3/2;
          for (int i=pmb->is; i<=pmb->ie; i++)
            e1(k,j,i)*=0.5;
        }
        if (pmb->block_size.nx2 > 1) { // 2D or 3D
          int i=pmb->is+pmb->block_size.nx1/2;
          for (int k=pmb->ks; k<=pmb->ke; k++)
            e3(k,j,i)*=0.5;
        }
      }
    }
    if (n==INNER_X3 || n==OUTER_X3) {
      if (n==INNER_X3) k=pmb->ks;
      else k=pmb->ke+1;
      nl=nblevel[2*n-8][1][1];
      if (nl==pmb->loc.level) { // same ; divide all the face EMFs by 2
        for (int j=pmb->js+1; j<=pmb->je; j++) {
          for (int i=pmb->is; i<=pmb->ie; i++)
            e1(k,j,i)*=0.5;
        }
        for (int j=pmb->js; j<=pmb->je; j++) {
          for (int i=pmb->is+1; i<=pmb->ie; i++)
            e2(k,j,i)*=0.5;
        }
      } else if (nl>pmb->loc.level) { // finer; divide the overlapping EMFs by 2
        // this is always 3D
        int j=pmb->js+pmb->block_size.nx2/2;
        for (int i=pmb->is; i<=pmb->ie; i++)
          e1(k,j,i)*=0.5;
        int i=pmb->is+pmb->block_size.nx1/2;
        for (int j=pmb->js; j<=pmb->je; j++)
          e2(k,j,i)*=0.5;
      }
    }
  }
  // edge
  for (int n=0; n<nedge_; n++) {
    if (nedge_fine_[n]==1) continue;
    Real div=1.0/static_cast<Real>(nedge_fine_[n]);
    NeighborBlock& nb=neighbor[n+6];
    Real half_div=div;
    if (nb.shear) half_div=0.5;
    // x1x2 edge (both 2D and 3D)
    if (n>=0 && n<4) {
      if ((n&1)==0) i=pmb->is;
      else i=pmb->ie+1;
      if ((n&2)==0) j=pmb->js;
      else j=pmb->je+1;
      for (int k=pmb->ks; k<=pmb->ke; k++)
        e3(k,j,i)*=half_div;
      // x1x3 edge
    } else if (n>=4 && n<8) {
      if ((n&1)==0) i=pmb->is;
      else i=pmb->ie+1;
      if ((n&2)==0) k=pmb->ks;
      else k=pmb->ke+1;
      for (int j=pmb->js; j<=pmb->je; j++)
        e2(k,j,i)*=half_div;
      // x2x3 edge
    } else if (n>=8 && n<12) {
      if ((n&1)==0) j=pmb->js;
      else j=pmb->je+1;
      if ((n&2)==0) k=pmb->ks;
      else k=pmb->ke+1;
      for (int i=pmb->is; i<=pmb->ie; i++)
        e1(k,j,i)*=div;
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::PolarSingleEMF(void)
//  \brief single CPU in the azimuthal direction for the polar boundary

void BoundaryValues::PolarSingleEMF(void) {
  MeshBlock *pmb=pmy_block_;
  AthenaArray<Real> &e1=pmb->pfield->e.x1e;
  AthenaArray<Real> &e3=pmb->pfield->e.x3e;
  int j;

  if (pmb->loc.level == pmb->pmy_mesh->root_level && pmb->pmy_mesh->nrbx3 == 1
  && pmb->block_size.nx3 > 1) {
    if (block_bcs[INNER_X2]==POLAR_BNDRY||block_bcs[INNER_X2]==POLAR_BNDRY_WEDGE) {
      j=pmb->js;
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      for (int i=pmb->is; i<=pmb->ie; i++) {
        Real tote1=0.0;
        for (int k=pmb->ks; k<=pmb->ke; k++)
          tote1+=e1(k,j,i);
        Real e1a=tote1/static_cast<double>(pmb->ke-pmb->ks+1);
        for (int k=pmb->ks; k<=pmb->ke+1; k++)
          e1(k,j,i)=e1a;
      }
      for (int i=pmb->is; i<=pmb->ie+1; i++) {
        for (int k=pmb->ks; k<=pmb->ke; k++)
          exc_(k)=e3(k,j,i);
        for (int k=pmb->ks; k<=pmb->ke; k++) {
          int k_shift = k;
          k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
          e3(k,j,i)=exc_(k_shift);
        }
      }
    }

    if (block_bcs[OUTER_X2]==POLAR_BNDRY||block_bcs[OUTER_X2]==POLAR_BNDRY_WEDGE) {
      j=pmb->je+1;
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      for (int i=pmb->is; i<=pmb->ie; i++) {
        Real tote1=0.0;
        for (int k=pmb->ks; k<=pmb->ke; ++k)
          tote1+=e1(k,j,i);
        Real e1a=tote1/static_cast<double>(pmb->ke-pmb->ks+1);
        for (int k=pmb->ks; k<=pmb->ke+1; ++k)
          e1(k,j,i)=e1a;
      }
      for (int i=pmb->is; i<=pmb->ie+1; i++) {
        for (int k=pmb->ks; k<=pmb->ke; k++)
          exc_(k)=e3(k,j,i);
        for (int k=pmb->ks; k<=pmb->ke; k++) {
          int k_shift = k;
          k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
          e3(k,j,i)=exc_(k_shift);
        }
      }
    }
  }
  return;
}



//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ReceiveEMFCorrection(void)
//  \brief Receive and Apply the surace EMF to the coarse neighbor(s) if needed

bool BoundaryValues::ReceiveEMFCorrection(void) {
  MeshBlock *pmb=pmy_block_;
  bool flag=true;

  // Receive same-level non-polar EMF values
  if (firsttime_==true) {
    for (int n=0; n<nneighbor; n++) { // first correct the same level
      NeighborBlock& nb = neighbor[n];
      if (nb.type!=NEIGHBOR_FACE && nb.type!=NEIGHBOR_EDGE) break;
      if (nb.level!=pmb->loc.level) continue;
      if ((nb.type==NEIGHBOR_FACE) || ((nb.type==NEIGHBOR_EDGE) &&
                                       (edge_flag_[nb.eid]==true))) {
        if (bd_emfcor_.flag[nb.bufid]==BNDRY_COMPLETED) continue;
        if (bd_emfcor_.flag[nb.bufid]==BNDRY_WAITING) {
          if (nb.rank==Globals::my_rank) { // on the same process
            flag=false;
            continue;
#ifdef MPI_PARALLEL
          } else { // MPI boundary
            int test;
            MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&test,MPI_STATUS_IGNORE);
            MPI_Test(&(bd_emfcor_.req_recv[nb.bufid]),&test,MPI_STATUS_IGNORE);
            if (static_cast<bool>(test)==false) {
              flag=false;
              continue;
            }
            bd_emfcor_.flag[nb.bufid] = BNDRY_ARRIVED;
          }
#else
          }
#endif
        }
        // boundary arrived; apply EMF correction
        SetEMFBoundarySameLevel(bd_emfcor_.recv[nb.bufid], nb);
        bd_emfcor_.flag[nb.bufid] = BNDRY_COMPLETED;
      }
    }

    if (flag==false) return flag;
    if (pmb->pmy_mesh->multilevel==true)
      ClearCoarseEMFBoundary();
    firsttime_=false;
  }

  // Receive finer non-polar EMF values
  if (pmb->pmy_mesh->multilevel==true) {
    for (int n=0; n<nneighbor; n++) { // then from finer
      NeighborBlock& nb = neighbor[n];
      if (nb.type!=NEIGHBOR_FACE && nb.type!=NEIGHBOR_EDGE) break;
      if (nb.level!=pmb->loc.level+1) continue;
      if (bd_emfcor_.flag[nb.bufid]==BNDRY_COMPLETED) continue;
      if (bd_emfcor_.flag[nb.bufid]==BNDRY_WAITING) {
        if (nb.rank==Globals::my_rank) {// on the same process
          flag=false;
          continue;
#ifdef MPI_PARALLEL
        } else { // MPI boundary
          int test;
          MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&test,MPI_STATUS_IGNORE);
          MPI_Test(&(bd_emfcor_.req_recv[nb.bufid]),&test,MPI_STATUS_IGNORE);
          if (static_cast<bool>(test)==false) {
            flag=false;
            continue;
          }
          bd_emfcor_.flag[nb.bufid] = BNDRY_ARRIVED;
        }
#else
        }
#endif
      }
      // boundary arrived; apply EMF correction
      SetEMFBoundaryFromFiner(bd_emfcor_.recv[nb.bufid], nb);
      bd_emfcor_.flag[nb.bufid] = BNDRY_COMPLETED;
    }
  }

  // Receive polar EMF values
  for (int n = 0; n < num_north_polar_blocks_; ++n) {
    const PolarNeighborBlock &nb = polar_neighbor_north[n];
    if (emf_north_flag_[n] == BNDRY_WAITING) {
      if (nb.rank == Globals::my_rank) { // on the same process
        flag = false;
        continue;
#ifdef MPI_PARALLEL
      } else {
        int recv_flag;
        MPI_Test(&req_emf_north_recv_[n], &recv_flag, MPI_STATUS_IGNORE);
        if (not recv_flag) {
          flag = false;
          continue;
        }
        emf_north_flag_[n] = BNDRY_ARRIVED;
      }
#else
      }
#endif
    }
  }
  for (int n = 0; n < num_south_polar_blocks_; ++n) {
    const PolarNeighborBlock &nb = polar_neighbor_south[n];
    if (emf_south_flag_[n] == BNDRY_WAITING) {
      if (nb.rank == Globals::my_rank) { // on the same process
        flag = false;
        continue;
#ifdef MPI_PARALLEL
      } else {
        int recv_flag;
        MPI_Test(&req_emf_south_recv_[n], &recv_flag, MPI_STATUS_IGNORE);
        if (not recv_flag) {
          flag = false;
          continue;
        }
        emf_south_flag_[n] = BNDRY_ARRIVED;
      }
#else
      }
#endif
    }
  }

  if (flag==true) {
    AverageEMFBoundary();
    if (num_north_polar_blocks_ > 0)
      SetEMFBoundaryPolar(emf_north_recv_, num_north_polar_blocks_, true);
    for (int n = 0; n < num_north_polar_blocks_; ++n)
      emf_north_flag_[n] = BNDRY_COMPLETED;
    if (num_south_polar_blocks_ > 0)
      SetEMFBoundaryPolar(emf_south_recv_, num_south_polar_blocks_, false);
    for (int n = 0; n < num_south_polar_blocks_; ++n)
      emf_south_flag_[n] = BNDRY_COMPLETED;
    if (block_bcs[INNER_X2]==POLAR_BNDRY||block_bcs[OUTER_X2]==POLAR_BNDRY||
        block_bcs[INNER_X2]==POLAR_BNDRY_WEDGE||block_bcs[OUTER_X2]==POLAR_BNDRY_WEDGE)
      PolarSingleEMF();
  }
  return flag;
}
