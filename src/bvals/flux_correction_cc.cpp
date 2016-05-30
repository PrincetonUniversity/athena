//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file flux_correction_cc.cpp
//  \brief functions that perform flux correction for CELL_CENTERED variables
//======================================================================================

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

// this class header
#include "bvals.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SendFluxCorrection(int step)
//  \brief Restrict, pack and send the surace flux to the coarse neighbor(s)
void BoundaryValues::SendFluxCorrection(int step)
{
  MeshBlock *pmb=pmy_mblock_;
  Coordinates *pco=pmb->pcoord;
  AthenaArray<Real> &x1flux=pmb->phydro->flux[x1face];
  AthenaArray<Real> &x2flux=pmb->phydro->flux[x2face];
  AthenaArray<Real> &x3flux=pmb->phydro->flux[x3face];
  int fx1=pmb->loc.lx1&1L, fx2=pmb->loc.lx2&1L, fx3=pmb->loc.lx3&1L;
  int fi1, fi2;

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    if(nb.type!=NEIGHBOR_FACE) break;
    if(nb.level==pmb->loc.level-1) {
      int p=0;
      // x1 direction
      if(nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
        int i=pmb->is+(pmb->ie-pmb->is+1)*nb.fid;
        fi1=fx2, fi2=fx3;
        if(pmb->block_size.nx3>1) { // 3D
          for(int nn=0; nn<NHYDRO; nn++) {
            for(int k=pmb->ks; k<=pmb->ke; k+=2) {
              for(int j=pmb->js; j<=pmb->je; j+=2) {
                Real amm=pco->GetFace1Area(k,   j,   i);
                Real amp=pco->GetFace1Area(k,   j+1, i);
                Real apm=pco->GetFace1Area(k+1, j,   i);
                Real app=pco->GetFace1Area(k+1, j+1, i);
                Real tarea=amm+amp+apm+app;
                flcor_send_[step][nb.fid][p++]=
                           (x1flux(nn, k  , j  , i)*amm
                           +x1flux(nn, k  , j+1, i)*amp
                           +x1flux(nn, k+1, j  , i)*apm
                           +x1flux(nn, k+1, j+1, i)*app)/tarea;
              }
            }
          }
        }
        else if(pmb->block_size.nx2>1) { // 2D
          int k=pmb->ks;
          for(int nn=0; nn<NHYDRO; nn++) {
            for(int j=pmb->js; j<=pmb->je; j+=2) {
              Real am=pco->GetFace1Area(k, j,   i);
              Real ap=pco->GetFace1Area(k, j+1, i);
              Real tarea=am+ap;
              flcor_send_[step][nb.fid][p++]=
                         (x1flux(nn, k, j  , i)*am
                         +x1flux(nn, k, j+1, i)*ap)/tarea;
            }
          }
        }
        else { // 1D
          int k=pmb->ks, j=pmb->js;
          for(int nn=0; nn<NHYDRO; nn++)
            flcor_send_[step][nb.fid][p++]=x1flux(nn, k, j, i);
        }
      }
      // x2 direction
      else if(nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
        int j=pmb->js+(pmb->je-pmb->js+1)*(nb.fid&1);
        fi1=fx1, fi2=fx3;
        if(pmb->block_size.nx3>1) { // 3D
          for(int nn=0; nn<NHYDRO; nn++) {
            for(int k=pmb->ks; k<=pmb->ke; k+=2) {
              pco->Face2Area(k  , j, pmb->is, pmb->ie, sarea_[0]);
              pco->Face2Area(k+1, j, pmb->is, pmb->ie, sarea_[1]);
              for(int i=pmb->is; i<=pmb->ie; i+=2) {
                Real tarea=sarea_[0](i)+sarea_[0](i+1)+sarea_[1](i)+sarea_[1](i+1);
                flcor_send_[step][nb.fid][p++]=
                           (x2flux(nn, k  , j, i  )*sarea_[0](i  )
                           +x2flux(nn, k  , j, i+1)*sarea_[0](i+1)
                           +x2flux(nn, k+1, j, i  )*sarea_[1](i  )
                           +x2flux(nn, k+1, j, i+1)*sarea_[1](i+1))/tarea;
              }
            }
          }
        }
        else if(pmb->block_size.nx2>1) { // 2D
          int k=pmb->ks;
          for(int nn=0; nn<NHYDRO; nn++) {
            pco->Face2Area(0, j, pmb->is ,pmb->ie, sarea_[0]);
            for(int i=pmb->is; i<=pmb->ie; i+=2) {
              Real tarea=sarea_[0](i)+sarea_[0](i+1);
              flcor_send_[step][nb.fid][p++]=
                         (x2flux(nn, k, j, i  )*sarea_[0](i  )
                         +x2flux(nn, k, j, i+1)*sarea_[0](i+1))/tarea;
            }
          }
        }
      }
      // x3 direction - 3D only
      else if(nb.fid==INNER_X3 || nb.fid==OUTER_X3) {
        int k=pmb->ks+(pmb->ke-pmb->ks+1)*(nb.fid&1);
        fi1=fx1, fi2=fx2;
        for(int nn=0; nn<NHYDRO; nn++) {
          for(int j=pmb->js; j<=pmb->je; j+=2) {
            pco->Face3Area(k, j,   pmb->is, pmb->ie, sarea_[0]);
            pco->Face3Area(k, j+1, pmb->is, pmb->ie, sarea_[1]);
            for(int i=pmb->is; i<=pmb->ie; i+=2) {
              Real tarea=sarea_[0](i)+sarea_[0](i+1)+sarea_[1](i)+sarea_[1](i+1);
              flcor_send_[step][nb.fid][p++]=
                         (x3flux(nn, k, j  , i  )*sarea_[0](i  )
                         +x3flux(nn, k, j  , i+1)*sarea_[0](i+1)
                         +x3flux(nn, k, j+1, i  )*sarea_[1](i  )
                         +x3flux(nn, k, j+1, i+1)*sarea_[1](i+1))/tarea;
            }
          }
        }
      }
      if(nb.rank==Globals::my_rank) { // on the same node
        MeshBlock *pbl=pmb->pmy_mesh->FindMeshBlock(nb.gid);
        std::memcpy(pbl->pbval->flcor_recv_[step][(nb.fid^1)][fi2][fi1],
                    flcor_send_[step][nb.fid], p*sizeof(Real));
        pbl->pbval->flcor_flag_[step][(nb.fid^1)][fi2][fi1]=BNDRY_ARRIVED;
      }
#ifdef MPI_PARALLEL
      else
        MPI_Start(&req_flcor_send_[step][nb.fid]);
#endif
    }
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn bool BoundaryValues::ReceiveFluxCorrection(int step)
//  \brief Receive and apply the surace flux from the finer neighbor(s)
bool BoundaryValues::ReceiveFluxCorrection(int step)
{
  MeshBlock *pmb=pmy_mblock_;
  Coordinates *pco=pmb->pcoord;
  AthenaArray<Real> &x1flux=pmb->phydro->flux[x1face];
  AthenaArray<Real> &x2flux=pmb->phydro->flux[x2face];
  AthenaArray<Real> &x3flux=pmb->phydro->flux[x3face];
  bool flag=true;

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    if(nb.type!=NEIGHBOR_FACE) break;
    if(nb.level==pmb->loc.level+1) {
      if(flcor_flag_[step][nb.fid][nb.fi2][nb.fi1]==BNDRY_COMPLETED) continue;
      if(flcor_flag_[step][nb.fid][nb.fi2][nb.fi1]==BNDRY_WAITING) {
        if(nb.rank==Globals::my_rank) {// on the same process
          flag=false;
          continue;
        }
#ifdef MPI_PARALLEL
        else { // MPI boundary
          int test;
          MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&test,MPI_STATUS_IGNORE);
          MPI_Test(&req_flcor_recv_[step][nb.fid][nb.fi2][nb.fi1],&test,MPI_STATUS_IGNORE);
          if(test==false) {
            flag=false;
            continue;
          }
          flcor_flag_[step][nb.fid][nb.fi2][nb.fi1] = BNDRY_ARRIVED;
        }
#endif
      }
      // boundary arrived; apply flux correction
      Real *buf=flcor_recv_[step][nb.fid][nb.fi2][nb.fi1];
      int p=0;
      if(nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
        int is=pmb->is+(pmb->ie-pmb->is)*nb.fid+nb.fid;
        int js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
        if(nb.fi1==0) je-=pmb->block_size.nx2/2;
        else          js+=pmb->block_size.nx2/2;
        if(nb.fi2==0) ke-=pmb->block_size.nx3/2;
        else          ks+=pmb->block_size.nx3/2;
        for(int nn=0; nn<NHYDRO; nn++) {
          for(int k=ks; k<=ke; k++) {
            for(int j=js; j<=je; j++)
              x1flux(nn,k,j,is)=buf[p++];
          }
        }
      }
      else if(nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
        int js=pmb->js+(pmb->je-pmb->js)*(nb.fid&1)+(nb.fid&1);
        int is=pmb->is, ie=pmb->ie, ks=pmb->ks, ke=pmb->ke;
        if(nb.fi1==0) ie-=pmb->block_size.nx1/2;
        else          is+=pmb->block_size.nx1/2;
        if(nb.fi2==0) ke-=pmb->block_size.nx3/2;
        else          ks+=pmb->block_size.nx3/2;
        for(int nn=0; nn<NHYDRO; nn++) {
          for(int k=ks; k<=ke; k++) {
            for(int i=is; i<=ie; i++)
              x2flux(nn,k,js,i)=buf[p++];
          }
        }
      }
      else if(nb.fid==INNER_X3 || nb.fid==OUTER_X3) {
        int ks=pmb->ks+(pmb->ke-pmb->ks)*(nb.fid&1)+(nb.fid&1);
        int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je;
        if(nb.fi1==0) ie-=pmb->block_size.nx1/2;
        else          is+=pmb->block_size.nx1/2;
        if(nb.fi2==0) je-=pmb->block_size.nx2/2;
        else          js+=pmb->block_size.nx2/2;
        for(int nn=0; nn<NHYDRO; nn++) {
          for(int j=js; j<=je; j++) {
            for(int i=is; i<=ie; i++)
              x3flux(nn,ks,j,i)=buf[p++];
          }
        }
      }

      flcor_flag_[step][nb.fid][nb.fi2][nb.fi1] = BNDRY_COMPLETED;
    }
  }

  return flag;
}
