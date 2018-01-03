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
//! \file bvals_cc.cpp
//  \brief functions that apply BCs for CELL_CENTERED variables
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
//! \fn int BoundaryValues::LoadEMFShearing(FaceField &src, Real *buf, int nb)
//  \brief Load shearingbox EMF boundary buffers
void BoundaryValues::LoadEMFShearing(EdgeField &src, Real *buf, const int nb)
{
  MeshBlock *pmb=pmy_block_;
  //Mesh *pmesh=pmb->pmy_mesh;
  int si, sj, sk, ei, ej, ek;
  int psj,pej; // indices for e3
  int nx2=pmb->block_size.nx2-NGHOST;
  sk=pmb->ks;        ek=pmb->ke;
  // --- inner boundary:
  // load [je-joverlap-(NGHOST-1):je] for send_inner[0]; nb=0
  // load [je-joverlap:je]            for send_inner[0]; nb=0  if joverlap == nx2-1
  //
  // load [js:je-joverlap+NGHOST]     for send_inner[1]; nb=1
  //     ([js:je-joverlap+NGHOST+1]   for e3)
  // load [js:je]                     for send_inner[1]; nb=1  if joverlap == 0 or 1
  //     ([js:je+1]                   for e3)
  //
  // load [je-(NGHOST-1):je]          for send_inner[2]; nb=2
  // load [je-(NGHOST-2):je]          for send_inner[2]; nb=2  if joverlap == nx2-1
  //
  // load [js:js+NGHOST-1]            for send_inner[3]; nb=3 if joverlap == 0
  //     ([js+1:js+NGHOST]            for e3)
  //
  // load [js:js+NGHOST-2]            for send_inner[4]; nb=4 if joverlap == 1
  //     ([js+1:js+NGHOST-1]          for e3)

  // --- outer boundary:
  // load [js:js+joverlap+NGHOST-1]   for send_outer[0]; nb=5
  //     ([js:js+joverlap+NGHOST-1+1] for e3)
  // load [js:js+joverlap]            for send_outer[0]; nb=5  if joverlap == nx2-1
  //     ([js:js+joverlap+1]          for e3)
  //
  // load [js+joverlap-NGHOST:je]     for send_outer[1]; nb=6
  // load [js:je]                     for send_outer[1]; nb=6  if joverlap == 0
  //     ([js:je+1]                   for e3)
  // load [js:je]                     for send_outer[1]; nb=6  if joverlap == 1
  //     ([js:je]                   for e3)
  //
  // load [js:js+NGHOST-1]            for send_outer[2]; nb=7
  //     ([js+1:js+NGHOST-1+1]        for e3)
  // load [js:js]                     for send_outer[2]; nb=7  if joverlap == nx2-1
  //     ([js+1:js+1]                 for e3)
  // load [je-NGHOST+1:je]            for send_outer[3]; nb=8
  //
  // load [je-NGHOST+2:je]            for send_outer[4]; nb=9  if joverlap == 1
  switch(nb) {
    case 0:
      sj=pmb->je-joverlap_-(NGHOST-1); ej=pmb->je;
      if(joverlap_>nx2) sj=pmb->js;
      psj=sj; pej=ej;
      break;
    case 1:
      sj=pmb->js; ej=pmb->je-joverlap_+NGHOST;
      if(joverlap_<NGHOST) ej=pmb->je;
      psj=sj; pej=ej+1;
      break;
    case 2:
      sj=pmb->je-(NGHOST-1); ej=pmb->je;
      if(joverlap_>nx2) sj=pmb->je-(joverlap_-nx2)+1;
      psj=sj; pej=ej;
      break;
    case 3:
      sj=pmb->js; ej=pmb->js+(NGHOST-1);
      if(joverlap_<NGHOST) ej=pmb->js+(NGHOST-joverlap_)-1;
      psj=sj+1; pej=ej+1;
      break;
    case 4:
      sj=pmb->js; ej=pmb->js+joverlap_+NGHOST-1;
      if(joverlap_>nx2) ej=pmb->je;
      psj=sj; pej=ej+1;
      break;
    case 5:
      sj=pmb->js+joverlap_-NGHOST; ej=pmb->je;
      if(joverlap_<NGHOST) sj=pmb->js;
      psj=sj; pej=ej+1;
      break;
    case 6:
      sj=pmb->js; ej=pmb->js+(NGHOST-1);
      if(joverlap_>nx2) ej=pmb->js+(joverlap_-nx2)-1;
      psj=sj+1; pej=ej+1;
      break;
    case 7:
      sj=pmb->je-(NGHOST-1); ej=pmb->je;
      if(joverlap_<NGHOST) sj=pmb->je-(NGHOST-joverlap_)+1;
      psj=sj; pej=ej;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues:LoadEMFShearing " << std::endl
          << "nb = " << nb << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
   }

  int p=0;
  // pack e2
  for(int k=sk; k<=ek+1; k++) {
    for(int j=sj; j<=ej; j++)
      buf[p++] = src.x2e(k,j);
  }
  // pack e3
  for(int k=sk; k<=ek; k++) {
    for(int j=psj; j<=pej; j++)
      buf[p++] = src.x3e(k,j);
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SendEMFShearingboxBoundaryCorrectionForInit(void)
//  \brief Send shearingbox boundary buffers for EMF correction
//  Currently not used. As we only correct vy at t=0.

void BoundaryValues::SendEMFShearingboxBoundaryCorrectionForInit(void)
{
  return;
}
//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SendEMFShearingboxBoundaryCorrection(void)
//  \brief Send shearingbox boundary buffers for EMF correction

void BoundaryValues::SendEMFShearingboxBoundaryCorrection(void)
{
  MeshBlock *pmb=pmy_block_;
  Coordinates *pco=pmb->pcoord;
  Mesh *pmesh=pmb->pmy_mesh;

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int nx2=pmb->block_size.nx2;
  int nx3=pmb->block_size.nx3;

  Real qomL = qshear_*Omega_0_*x1size_;
  AthenaArray<Real> &bx1=pmb->pfield->b.x1f;

  if (shbb_.inner == true) {
    // step 1. -- average edges of shboxvar_emf_
    // average e3 for x1x2 edge
    for(int k=ks; k<=ke; k++) {
      for(int j=js; j<=je+1; j+=nx2)
        shboxvar_inner_emf_.x3e(k,j) *= 0.5;
    }
    // average e2 for x1x3 edge
    for(int k=ks; k<=ke+1; k+=nx3) {
      for(int j=js; j<=je; j++)
        shboxvar_inner_emf_.x2e(k,j) *= 0.5;
    }

    // step 2. -- load sendbuf; memcpy to recvbuf if on same rank, post MPI_Isend otherwise
    for(int n=0; n<4; n++) {
      if(send_inner_rank_[n] != -1) {
        LoadEMFShearing(shboxvar_inner_emf_, send_innerbuf_emf_[n], n);
        if (send_inner_rank_[n] == Globals::my_rank) {// on the same process
          MeshBlock *pbl=pmb->pmy_mesh->FindMeshBlock(send_inner_gid_[n]);
          std::memcpy(pbl->pbval->recv_innerbuf_emf_[n],
                  send_innerbuf_emf_[n], send_innersize_emf_[n]*sizeof(Real));
          pbl->pbval->shbox_inner_emf_flag_[n]=BNDRY_ARRIVED;
        } else { // MPI
#ifdef MPI_PARALLEL
          int tag=CreateBvalsMPITag(send_inner_lid_[n], TAG_SHBOX_EMF, n); //bufid = n
          MPI_Isend(send_innerbuf_emf_[n],send_innersize_emf_[n],MPI_ATHENA_REAL,
                    send_inner_rank_[n],tag,MPI_COMM_WORLD, &rq_innersend_emf_[n]);
#endif
       }
    }}
  } // inner boundaries

  if (shbb_.outer == true) {
    // step 1. -- average edges of shboxvar_emf_
    // average e3 for x1x2 edge
    for(int k=ks; k<=ke; k++) {
      for(int j=js; j<=je+1; j+=nx2)
        shboxvar_outer_emf_.x3e(k,j) *= 0.5;
    }
    // average e2 for x1x3 edge
    for(int k=ks; k<=ke+1; k+=nx3) {
      for(int j=js; j<=je; j++)
        shboxvar_outer_emf_.x2e(k,j) *= 0.5;
    }

    // step 2. -- load sendbuf; memcpy to recvbuf if on same rank, post MPI_Isend otherwise
    int offset = 4;
    for(int n=0; n<4; n++) {
      if(send_outer_rank_[n] != -1) {
        LoadEMFShearing(shboxvar_outer_emf_, send_outerbuf_emf_[n], n+offset);
        if (send_outer_rank_[n] == Globals::my_rank) {// on the same process
          MeshBlock *pbl=pmb->pmy_mesh->FindMeshBlock(send_outer_gid_[n]);
          std::memcpy(pbl->pbval->recv_outerbuf_emf_[n],
                  send_outerbuf_emf_[n], send_outersize_emf_[n]*sizeof(Real));
          pbl->pbval->shbox_outer_emf_flag_[n]=BNDRY_ARRIVED;
        } else { // MPI
#ifdef MPI_PARALLEL
          int tag=CreateBvalsMPITag(send_outer_lid_[n], TAG_SHBOX_EMF, n+offset); //bufid for outer(inner): 2(0) and 3(1)
          MPI_Isend(send_outerbuf_emf_[n],send_outersize_emf_[n],MPI_ATHENA_REAL,
                    send_outer_rank_[n],tag,MPI_COMM_WORLD, &rq_outersend_emf_[n]);
#endif
        }
    }}
  } // outer boundaries
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetEMFShearingboxBoundarySameLevel(EdgeField &dst, Real *buf, const int nb)
//  \brief Set EMF shearingbox boundary received from a block on the same level
void BoundaryValues::SetEMFShearingboxBoundarySameLevel(EdgeField &dst, Real *buf, const int nb)
{
  MeshBlock *pmb=pmy_block_;
  Mesh *pmesh=pmb->pmy_mesh;
  int si, sj, sk, ei, ej, ek;
  int psj,pej;
  int nx2=pmb->block_size.nx2-NGHOST;
  int nxo=pmb->block_size.nx2-joverlap_;

  sk = pmb->ks; ek = pmb->ke;
  //if(nb < 4) si = pmb->is;
  //else si = pmb->ie+1;
  // --- inner boundary:
  // set [js-NGHOST:js+(joverlap-1)] with recv_inner[0]; nb=0
  // set [js-1:js+(joverlap-1)]      with recv_inner[0]; nb=0 if joverlap==nx2-1
  // set [js+joverlap:je+NGHOST]     with recv_inner[1]; nb=1
  //    ([js+joverlap:je+NGHOST+1]   for e3)
  // set [js+joverlap:je]            with recv_inner[1]; nb=1 if joverlap==0
  //    ([js+joverlap:je+1]          for e3)
  // set [js-NGHOST:js-1]            with recv_inner[2]; nb=2
  // set [js-NGHOST:js-NGHOST]       with recv_inner[2]; nb=2 if joverlap==nx2-1
  // set [je+1:je+NGHOST]            with recv_inner[3]; nb=3 if joverlap==0
  //    ([je+2:je+NGHOST+1]          for e3)
  // set [je+2:je+NGHOST]            with recv_inner[4]; nb=4 if joverlap==1
  //    ([je+3:je+NGHOST+1]          for e3)
  //
  // --- outer boundary:
  // set [je-joverlap+1:je+NGHOST]   with recv_outer[0]; nb=5
  //    ([je-joverlap+1:je+NGHOST+1] for e3)
  // set [je-joverlap+1:je+1]        with recv_outer[0]; nb=5 if joverlap==nx2-1
  //    ([je-joverlap+1:je+2]        for e3)
  //
  // set [js-NGHOST:je-joverlap]     with recv_outer[1]; nb=6
  // set [js:je]                     with recv_outer[1]; nb=6 if joverlap==0 or 1
  //    ([js:je+1]                   for e3 if joverlap<=1)
  // set [je+1:je+NGHOST]            with recv_outer[2]; nb=7
  //    ([je+2:je+NGHOST+1]          for e3)
  // set [je+NGHOST:je+NGHOST]       with recv_outer[2]; nb=7 if joverlap==nx2-1
  //    ([je+NGHOST+1:je+NGHOST+1]   for e3)
  // set [js-NGHOST:js-1]            with recv_outer[3]; nb=8 if joverlap==0
  // set [js-NGHOST:js-2]            with recv_outer[4]; nb=9 if joverlap==1
  //
  switch(nb) {
    case 0:
      sj=pmb->js-NGHOST; ej=pmb->js+(joverlap_-1);
      if(joverlap_>nx2) sj=pmb->js-nxo;
      psj=sj; pej=ej+1;
      break;
    case 1:
      sj=pmb->js+joverlap_; ej=pmb->je+NGHOST;
      if(joverlap_<NGHOST) ej=pmb->je+joverlap_;
      psj=sj; pej=ej+1;
      break;
    case 2:
      sj=pmb->js-NGHOST; ej=pmb->js-1;
      if(joverlap_>nx2) ej=pmb->js-nxo-1;
      psj=sj; pej=ej;
      break;
    case 3:
      sj=pmb->je+joverlap_+1; ej=pmb->je+NGHOST;
      psj=sj+1; pej=ej+1;
      break;
    case 4:
      sj=pmb->je-(joverlap_-1); ej=pmb->je+NGHOST;
      if(joverlap_>nx2) ej=pmb->je+nxo;
      psj=sj; pej=ej+1;
      break;
    case 5:
      sj=pmb->js-NGHOST; ej=pmb->je-joverlap_;
      if(joverlap_<=NGHOST) sj=pmb->js-joverlap_;
      psj=sj; pej=ej+1;
      break;
    case 6:
      sj=pmb->je+1; ej=pmb->je+NGHOST;
      if(joverlap_>nx2) sj=pmb->je+nxo+1;
      psj=sj+1; pej=ej+1;
      break;
    case 7:
      sj=pmb->js-NGHOST; ej=pmb->js-joverlap_-1;
      psj=sj; pej=ej;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues:SetFieldShearing " << std::endl
          << "nb = " << nb << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
   }

  int p=0;
  // unpack e2
  for(int k=sk; k<=ek+1; k++) {
    for(int j=sj; j<=ej; j++)
      dst.x2e(k,j)+=buf[p++];
  }
 // unpack e3
  for(int k=sk; k<=ek; k++) {
    for(int j=psj; j<=pej; j++)
      dst.x3e(k,j)+=buf[p++];
  }

  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ReceiveEMFShearingboxBoundaryCorrectionWithWait(void)
//  \brief receive shearingbox boundary data for EMF correction for initialization
//         It is not used for now as our problem generator might start with t=0 periodic config.
void BoundaryValues::ReceiveEMFShearingboxBoundaryCorrectionWithWait(void)
{
  return;
}
//--------------------------------------------------------------------------------------
//! \fn bool BoundaryValues::ReceiveEMFShearingboxBoundaryCorrection(void)
//  \brief receive shearingbox boundary data for EMF correction
bool BoundaryValues::ReceiveEMFShearingboxBoundaryCorrection(void)
{
  MeshBlock *pmb=pmy_block_;
  Mesh *pmesh=pmb->pmy_mesh;
  bool flagi=true, flago=true;

  if(shbb_.inner == true) { // check inner boundaries
    for(int n=0; n<4; n++) {
      if(shbox_inner_emf_flag_[n]==BNDRY_COMPLETED) continue;
      if(shbox_inner_emf_flag_[n]==BNDRY_WAITING) {
        if (recv_inner_rank_[n]==Globals::my_rank) {// on the same process
          flagi=false;
          continue;
        }
        else { // MPI boundary
#ifdef MPI_PARALLEL
          int test;
          MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&test,MPI_STATUS_IGNORE);
          MPI_Test(&rq_innerrecv_emf_[n],&test,MPI_STATUS_IGNORE);
          if(test==false) {
            flagi=false;
            continue;
          }
          shbox_inner_emf_flag_[n] = BNDRY_ARRIVED;
#endif
        }
      }
      // set dst if boundary arrived
      SetEMFShearingboxBoundarySameLevel(shboxmap_inner_emf_,recv_innerbuf_emf_[n],n);
      shbox_inner_emf_flag_[n] = BNDRY_COMPLETED; // completed
    }
  } // inner boundary

  if(shbb_.outer == true) { // check outer boundaries
    int offset = 4;
    for(int n=0; n<4; n++) {
      if(shbox_outer_emf_flag_[n]==BNDRY_COMPLETED) continue;
      if(shbox_outer_emf_flag_[n]==BNDRY_WAITING) {
        if (recv_outer_rank_[n]==Globals::my_rank) {// on the same process
          flago=false;
          continue;
        }
        else { // MPI boundary
#ifdef MPI_PARALLEL
          int test;
          MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&test,MPI_STATUS_IGNORE);
          MPI_Test(&rq_outerrecv_emf_[n],&test,MPI_STATUS_IGNORE);
          if(test==false) {
            flago=false;
            continue;
          }
          shbox_outer_emf_flag_[n] = BNDRY_ARRIVED;
#endif
        }
      }
      SetEMFShearingboxBoundarySameLevel(shboxmap_outer_emf_,recv_outerbuf_emf_[n],n+offset);
      shbox_outer_emf_flag_[n] = BNDRY_COMPLETED; // completed
    }
  } // outer boundary

  return (flagi && flago);

}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::RemapEMFShearingboxBoundary(void)
//  \brief Set EMF boundary received from a block on the finer level
void BoundaryValues::RemapEMFShearingboxBoundary(void)
{
  MeshBlock *pmb=pmy_block_;
  AthenaArray<Real> &e2=pmb->pfield->e.x2e;
  AthenaArray<Real> &e3=pmb->pfield->e.x3e;
  int ks=pmb->ks, ke=pmb->ke;
  int js=pmb->js, je=pmb->je;
  int is=pmb->is, ie=pmb->ie;
  if(shbb_.inner==true) {
    ClearEMFShearing(shboxvar_inner_emf_);
    // step 1.-- conservative remapping
    /*
    for(int k=ks; k<=ke; k++) {  // e3
      RemapFluxEMF(k,js,je+3,eps_,shboxmap_inner_emf_.x3e,flx_inner_emf_.x3e);
      for(int j=js; j<=je+1; j++) {
        shboxmap_inner_emf_.x3e(k,j) -= flx_inner_emf_.x3e(j+1)-flx_inner_emf_.x3e(j);
      }
    } */
    for(int k=ks; k<=ke+1; k++) { // e2
      RemapFluxEMF(k,js,je+2,eps_,shboxmap_inner_emf_.x2e,flx_inner_emf_.x2e);
      for(int j=js; j<=je; j++) {
        shboxmap_inner_emf_.x2e(k,j) -= flx_inner_emf_.x2e(j+1)-flx_inner_emf_.x2e(j);
      }
    }
    // step 2.-- average the EMF correction
    // average e3
    /*
    for(int k=ks; k<=ke; k++) {
      for(int j=js; j<=je+1; j++)
      //for(int j=js-NGHOST; j<=je+1+NGHOST; j++)
        e3(k,j,is) = 0.5*(e3(k,j,is)+shboxmap_inner_emf_.x3e(k,j));
    } */
    // average e2
    for(int k=ks; k<=ke+1; k++) {
      for(int j=js; j<=je; j++)
      //for(int j=js-NGHOST; j<=je+NGHOST; j++)
        e2(k,j,is) = 0.5*(e2(k,j,is)+shboxmap_inner_emf_.x2e(k,j));
    }
    ClearEMFShearing(shboxmap_inner_emf_);
  }

  if(shbb_.outer==true) {
    ClearEMFShearing(shboxvar_outer_emf_);
    // step 1.-- conservative remapping
    /*
    for(int k=ks; k<=ke; k++) {  // e3
      RemapFluxEMF(k,js-1,je+2,-eps_,shboxmap_outer_emf_.x3e,flx_outer_emf_.x3e);
      for(int j=js; j<=je+1; j++)
        shboxmap_outer_emf_.x3e(k,j) -= flx_outer_emf_.x3e(j+1)-flx_outer_emf_.x3e(j);
    } */
    for(int k=ks; k<=ke+1; k++) { // e2
      RemapFluxEMF(k,js-1,je+1,-eps_,shboxmap_outer_emf_.x2e,flx_outer_emf_.x2e);
      for(int j=js; j<=je; j++)
        shboxmap_outer_emf_.x2e(k,j) -= flx_outer_emf_.x2e(j+1)-flx_outer_emf_.x2e(j);
    }
    // step 2.-- average the EMF correction
    // average e3
    /*
    for(int k=ks; k<=ke; k++) {
      for(int j=js; j<=je+1; j++)
      //for(int j=js-NGHOST; j<=je+1+NGHOST; j++)
        e3(k,j,ie+1) = 0.5*(e3(k,j,ie+1)+shboxmap_outer_emf_.x3e(k,j));
    } */
    // average e2
    for(int k=ks; k<=ke+1; k++) {
      for(int j=js; j<=je; j++)
      //for(int j=js-NGHOST; j<=je+NGHOST; j++)
        e2(k,j,ie+1) = 0.5*(e2(k,j,ie+1)+shboxmap_outer_emf_.x2e(k,j));
    }
    ClearEMFShearing(shboxmap_outer_emf_);
  }


  return;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ClearEMFShearing(void)
//  \brief Clear the working array for EMFs on the surface/edge contacting with a shearing periodic boundary
void BoundaryValues::ClearEMFShearing(EdgeField &work)
{
  MeshBlock *pmb=pmy_block_;
  AthenaArray<Real> &e2=work.x2e;
  AthenaArray<Real> &e3=work.x3e;
  int ks=pmb->ks, ke=pmb->ke;
  int js=pmb->js, je=pmb->je;
  for(int k=ks-NGHOST; k<=ke+NGHOST; k++) {
    for(int j=js-NGHOST; j<=je+NGHOST; j++) {
      e2(k,j) = 0.0;
      e3(k,j) = 0.0;
      if(k==ke+NGHOST) e2(k+1,j) = 0.0;
      if(j==je+NGHOST) e3(k,j+1) = 0.0;
  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::RemapFluxEMF(int k, int jinner, int jouter, Real eps, static AthenaArray<Real> &U, AthenaArray<Real> &Flux)
//  \brief compute the flux along j indices for remapping
//  adopted from 2nd order RemapFlux of Athena4.0
void BoundaryValues::RemapFluxEMF(const int k, const int jinner, const int jouter, const Real eps, const AthenaArray<Real> &U, AthenaArray<Real> &Flux)
{
  int j,jl,ju;
  Real dUc,dUl,dUr,dUm,lim_slope;

// jinner,jouter are index range over which flux must be returned.  Set loop
// limits depending on direction of upwind differences

  if (eps > 0.0) { // eps always > 0 for inner i boundary
    jl = jinner-1;
    ju = jouter-1;
  } else {         // eps always < 0 for outer i boundary
    jl = jinner;
    ju = jouter;
  }

  for (j=jl; j<=ju; j++) {
      dUc = U(k,j+1) - U(k,j-1);
      dUl = U(k,j  ) - U(k,j-1);
      dUr = U(k,j+1) - U(k,j  );

      dUm = 0.0;
      if (dUl*dUr > 0.0) {
        lim_slope = std::min(fabs(dUl),fabs(dUr));
        dUm = SIGN(dUc)*std::min(0.5*fabs(dUc),2.0*lim_slope);
      }

    if (eps > 0.0) { // eps always > 0 for inner i boundary
      Flux(j+1) = eps*(U(k,j) + 0.5*(1.0 - eps)*dUm);
    } else {         // eps always < 0 for outer i boundary
      Flux(j  ) = eps*(U(k,j) - 0.5*(1.0 + eps)*dUm);
    }
  }

  return;
}
