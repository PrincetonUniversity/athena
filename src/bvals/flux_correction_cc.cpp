//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file flux_correction_cc.cpp
//  \brief functions that perform flux correction for CELL_CENTERED variables

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

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../utils/buffer_utils.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SendFluxCorrection(enum FluxCorrectionType type)
//  \brief Restrict, pack and send the surace flux to the coarse neighbor(s)

void BoundaryValues::SendFluxCorrection(enum FluxCorrectionType type) {
  MeshBlock *pmb=pmy_block_, *pbl;
  Coordinates *pco=pmb->pcoord;
  AthenaArray<Real> x1flux, x2flux, x3flux;
  int ns, ne;
  BoundaryData *pbd, *ptarget;

  if (type==FLUX_HYDRO) {
    ns=0, ne=NHYDRO-1;
    x1flux.InitWithShallowCopy(pmb->phydro->flux[X1DIR]);
    x2flux.InitWithShallowCopy(pmb->phydro->flux[X2DIR]);
    x3flux.InitWithShallowCopy(pmb->phydro->flux[X3DIR]);
    pbd=&bd_flcor_;
  }

  for (int n=0; n<nneighbor; n++) {
    NeighborBlock& nb = neighbor[n];
    if (nb.type!=NEIGHBOR_FACE) break;
    if (nb.level==pmb->loc.level-1) {
      int p=0;
      Real *sbuf=pbd->send[nb.bufid];
      // x1 direction
      if (nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
        int i=pmb->is+(pmb->ie-pmb->is+1)*nb.fid;
        if (pmb->block_size.nx3>1) { // 3D
          for (int nn=ns; nn<=ne; nn++) {
            for (int k=pmb->ks; k<=pmb->ke; k+=2) {
              for (int j=pmb->js; j<=pmb->je; j+=2) {
                Real amm=pco->GetFace1Area(k,   j,   i);
                Real amp=pco->GetFace1Area(k,   j+1, i);
                Real apm=pco->GetFace1Area(k+1, j,   i);
                Real app=pco->GetFace1Area(k+1, j+1, i);
                Real tarea=amm+amp+apm+app;
                sbuf[p++]= (x1flux(nn, k  , j  , i)*amm
                           +x1flux(nn, k  , j+1, i)*amp
                           +x1flux(nn, k+1, j  , i)*apm
                           +x1flux(nn, k+1, j+1, i)*app)/tarea;
              }
            }
          }
        } else if (pmb->block_size.nx2>1) { // 2D
          int k=pmb->ks;
          for (int nn=ns; nn<=ne; nn++) {
            for (int j=pmb->js; j<=pmb->je; j+=2) {
              Real am=pco->GetFace1Area(k, j,   i);
              Real ap=pco->GetFace1Area(k, j+1, i);
              Real tarea=am+ap;
              sbuf[p++]=(x1flux(nn, k, j  , i)*am
                        +x1flux(nn, k, j+1, i)*ap)/tarea;
            }
          }
        } else { // 1D
          int k=pmb->ks, j=pmb->js;
          for (int nn=ns; nn<=ne; nn++)
            sbuf[p++]=x1flux(nn, k, j, i);
        }
        // x2 direction
      } else if (nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
        int j=pmb->js+(pmb->je-pmb->js+1)*(nb.fid&1);
        if (pmb->block_size.nx3>1) { // 3D
          for (int nn=ns; nn<=ne; nn++) {
            for (int k=pmb->ks; k<=pmb->ke; k+=2) {
              pco->Face2Area(k  , j, pmb->is, pmb->ie, sarea_[0]);
              pco->Face2Area(k+1, j, pmb->is, pmb->ie, sarea_[1]);
              for (int i=pmb->is; i<=pmb->ie; i+=2) {
                Real tarea=sarea_[0](i)+sarea_[0](i+1)+sarea_[1](i)+sarea_[1](i+1);
                sbuf[p++]=(x2flux(nn, k  , j, i  )*sarea_[0](i  )
                          +x2flux(nn, k  , j, i+1)*sarea_[0](i+1)
                          +x2flux(nn, k+1, j, i  )*sarea_[1](i  )
                          +x2flux(nn, k+1, j, i+1)*sarea_[1](i+1))/tarea;
              }
            }
          }
        } else if (pmb->block_size.nx2>1) { // 2D
          int k=pmb->ks;
          for (int nn=ns; nn<=ne; nn++) {
            pco->Face2Area(0, j, pmb->is ,pmb->ie, sarea_[0]);
            for (int i=pmb->is; i<=pmb->ie; i+=2) {
              Real tarea=sarea_[0](i)+sarea_[0](i+1);
              sbuf[p++]=(x2flux(nn, k, j, i  )*sarea_[0](i  )
                        +x2flux(nn, k, j, i+1)*sarea_[0](i+1))/tarea;
            }
          }
        }
        // x3 direction - 3D only
      } else if (nb.fid==INNER_X3 || nb.fid==OUTER_X3) {
        int k=pmb->ks+(pmb->ke-pmb->ks+1)*(nb.fid&1);
        for (int nn=ns; nn<=ne; nn++) {
          for (int j=pmb->js; j<=pmb->je; j+=2) {
            pco->Face3Area(k, j,   pmb->is, pmb->ie, sarea_[0]);
            pco->Face3Area(k, j+1, pmb->is, pmb->ie, sarea_[1]);
            for (int i=pmb->is; i<=pmb->ie; i+=2) {
              Real tarea=sarea_[0](i)+sarea_[0](i+1)+sarea_[1](i)+sarea_[1](i+1);
              sbuf[p++]=(x3flux(nn, k, j  , i  )*sarea_[0](i  )
                        +x3flux(nn, k, j  , i+1)*sarea_[0](i+1)
                        +x3flux(nn, k, j+1, i  )*sarea_[1](i  )
                        +x3flux(nn, k, j+1, i+1)*sarea_[1](i+1))/tarea;
            }
          }
        }
      }
      if (nb.rank==Globals::my_rank) { // on the same node
        pbl=pmb->pmy_mesh->FindMeshBlock(nb.gid);
        if (type==FLUX_HYDRO)
          ptarget=&(pbl->pbval->bd_flcor_);
        std::memcpy(ptarget->recv[nb.targetid], sbuf, p*sizeof(Real));
        ptarget->flag[nb.targetid]=BNDRY_ARRIVED;
      }
#ifdef MPI_PARALLEL
      else
        MPI_Start(&(pbd->req_send[nb.bufid]));
#endif
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn bool BoundaryValues::ReceiveFluxCorrection(enum FluxCorrectionType type)
//  \brief Receive and apply the surace flux from the finer neighbor(s)

bool BoundaryValues::ReceiveFluxCorrection(enum FluxCorrectionType type) {
  MeshBlock *pmb=pmy_block_;
  AthenaArray<Real> x1flux, x2flux, x3flux;
  bool bflag=true;
  int ns, ne;
  BoundaryData *pbd;

  if (type==FLUX_HYDRO) {
    pbd=&bd_flcor_;
    ns=0, ne=NHYDRO-1;
    x1flux.InitWithShallowCopy(pmb->phydro->flux[X1DIR]);
    x2flux.InitWithShallowCopy(pmb->phydro->flux[X2DIR]);
    x3flux.InitWithShallowCopy(pmb->phydro->flux[X3DIR]);
  }

  for (int n=0; n<nneighbor; n++) {
    NeighborBlock& nb = neighbor[n];
    if (nb.type!=NEIGHBOR_FACE) break;
    if (nb.level==pmb->loc.level+1) {
      if (pbd->flag[nb.bufid]==BNDRY_COMPLETED) continue;
      if (pbd->flag[nb.bufid]==BNDRY_WAITING) {
        if (nb.rank==Globals::my_rank) {// on the same process
          bflag=false;
          continue;
#ifdef MPI_PARALLEL
        } else { // MPI boundary
          int test;
          MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&test,MPI_STATUS_IGNORE);
          MPI_Test(&(pbd->req_recv[nb.bufid]),&test,MPI_STATUS_IGNORE);
          if (static_cast<bool>(test)==false) {
            bflag=false;
            continue;
          }
          pbd->flag[nb.bufid] = BNDRY_ARRIVED;
        }
#else
        }
#endif
      }
      // boundary arrived; apply flux correction
      int p=0;
      Real *rbuf=pbd->recv[nb.bufid];
      if (nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
        int is=pmb->is+(pmb->ie-pmb->is)*nb.fid+nb.fid;
        int js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
        if (nb.fi1==0) je-=pmb->block_size.nx2/2;
        else          js+=pmb->block_size.nx2/2;
        if (nb.fi2==0) ke-=pmb->block_size.nx3/2;
        else          ks+=pmb->block_size.nx3/2;
        for (int nn=ns; nn<=ne; nn++) {
          for (int k=ks; k<=ke; k++) {
            for (int j=js; j<=je; j++)
              x1flux(nn,k,j,is)=rbuf[p++];
          }
        }
      } else if (nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
        int js=pmb->js+(pmb->je-pmb->js)*(nb.fid&1)+(nb.fid&1);
        int is=pmb->is, ie=pmb->ie, ks=pmb->ks, ke=pmb->ke;
        if (nb.fi1==0) ie-=pmb->block_size.nx1/2;
        else          is+=pmb->block_size.nx1/2;
        if (nb.fi2==0) ke-=pmb->block_size.nx3/2;
        else          ks+=pmb->block_size.nx3/2;
        for (int nn=ns; nn<=ne; nn++) {
          for (int k=ks; k<=ke; k++) {
            for (int i=is; i<=ie; i++)
              x2flux(nn,k,js,i)=rbuf[p++];
          }
        }
      } else if (nb.fid==INNER_X3 || nb.fid==OUTER_X3) {
        int ks=pmb->ks+(pmb->ke-pmb->ks)*(nb.fid&1)+(nb.fid&1);
        int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je;
        if (nb.fi1==0) ie-=pmb->block_size.nx1/2;
        else          is+=pmb->block_size.nx1/2;
        if (nb.fi2==0) je-=pmb->block_size.nx2/2;
        else          js+=pmb->block_size.nx2/2;
        for (int nn=ns; nn<=ne; nn++) {
          for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++)
              x3flux(nn,ks,j,i)=rbuf[p++];
          }
        }
      }

      pbd->flag[nb.bufid] = BNDRY_COMPLETED;
    }
  }

  return bflag;
}
