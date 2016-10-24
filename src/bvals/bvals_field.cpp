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
//! \fn int BoundaryValues::LoadFieldShearing(FaceField &src, Real *buf, int nb)
//  \brief Load shearingbox field boundary buffers
void BoundaryValues::LoadFieldShearing(FaceField &src, Real *buf, int nb)
{
  MeshBlock *pmb=pmy_mblock_;
  Mesh *pmesh=pmb->pmy_mesh;
  int si, sj, sk, ei, ej, ek;
  int psj,pej; // indices for bx2
  int nx2=pmb->block_size.nx2-1;

  si=pmb->is-NGHOST; ei=pmb->is-1;
  sk=pmb->ks;        ek=pmb->ke;
  if (pmesh->mesh_size.nx3>1)  ek += NGHOST, sk -= NGHOST;
  // --- inner boundary:
  // load [je-joverlap-(NGHOST-1):je] for send_inner[0]; nb=0
  // load [je-joverlap:je]            for send_inner[0]; nb=0  if joverlap == nx2-1
  // load [js:je-joverlap+NGHOST]     for send_inner[1]; nb=1
  //     ([js:je-joverlap+NGHOST+1]   for bx2)
  // load [js:je]                     for send_inner[1]; nb=1  if joverlap == 0
  //     ([js:je+1]                   for bx2)
  // load [je-(NGHOST-1):je]          for send_inner[2]; nb=2
  // load [je:je]                     for send_inner[2]; nb=2  if joverlap == nx2-1
  // load [js:js+NGHOST-1]            for send_inner[3]; nb=3
  //     ([js+1:js+NGHOST]            for bx2)
  // --- outer boundary:
  // load [js:js+joverlap+NGHOST-1]   for send_outer[0]; nb=4
  //     ([js:js+joverlap+NGHOST-1+1] for bx2)
  // load [js:js+joverlap]            for send_outer[0]; nb=4  if joverlap == nx2-1
  //     ([js:js+joverlap+1]          for bx2)
  // load [js+joverlap-NGHOST:je]     for send_outer[1]; nb=5
  // load [js:je]                     for send_outer[1]; nb=5  if joverlap == 0
  //     ([js:je+1]                   for bx2)
  // load [js:js+NGHOST-1]            for send_outer[2]; nb=6
  //     ([js+1:js+NGHOST-1+1]        for bx2)
  // load [js:js]                     for send_outer[2]; nb=6  if joverlap == nx2-1
  //     ([js+1:js+1]                 for bx2)
  // load [je-NGHOST+1:je]            for send_outer[3]; nb=7
  switch(nb) {
    case 0:
      sj=pmb->je-joverlap_-(NGHOST-1); ej=pmb->je;
      if(joverlap_==nx2) sj=pmb->je-joverlap_;
	  psj=sj; pej=ej;
	  break;
	case 1:
      sj=pmb->js; ej=pmb->je-joverlap_+NGHOST;
      if(joverlap_==0) ej=pmb->je;
	  psj=sj; pej=ej+1;
	  break;
	case 2:
      sj=pmb->je-(NGHOST-1); ej=pmb->je;
      if(joverlap_==nx2) sj=pmb->je;
	  psj=sj; pej=ej;
	  break;
	case 3:
	  sj=pmb->js; ej=pmb->js+(NGHOST-1);
	  psj=sj+1; pej=ej+1;
	  break;
	case 4:
	  sj=pmb->js; ej=pmb->js+joverlap_+NGHOST-1;
	  if(joverlap_==nx2) ej=pmb->js+joverlap_;
	  psj=sj; pej=ej+1;
	  break;
	case 5:
	  sj=pmb->js+joverlap_-NGHOST; ej=pmb->je;
	  psj=sj; pej=ej;
	  if(joverlap_==0) {
		sj=pmb->js;
		psj=sj;
		pej=ej+1;
	  }
      break;
	case 6:
	  sj=pmb->js; ej=pmb->js+(NGHOST-1);
	  if(joverlap_==nx2) ej=pmb->js;
	  psj=sj+1; pej=ej+1;
	  break;
	case 7:
	  sj=pmb->je-NGHOST+1; ej=pmb->je;
	  psj=sj; pej=ej;
	  break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues:LoadFieldShearing " << std::endl
          << "nb = " << nb << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
   }

  int p=0;
  BufferUtility::Pack3DData(src.x1f, buf, si, ei, sj, ej, sk, ek, p);
  BufferUtility::Pack3DData(src.x2f, buf, si, ei, psj, pej, sk, ek, p);
  BufferUtility::Pack3DData(src.x3f, buf, si, ei, sj, ej, sk, ek+1, p);

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SendFieldShearingboxBoundaryBuffersForInit(FaceField &src)
//  \brief Send shearingbox boundary buffers for field variables
//  Currently not used. As we only correct vy at t=0.

void BoundaryValues::SendFieldShearingboxBoundaryBuffersForInit(FaceField &src, bool cons)
{
  return;
}
//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SendFieldShearingboxBoundaryBuffers(FaceField &src)
//  \brief Send shearingbox boundary buffers for field variables

void BoundaryValues::SendFieldShearingboxBoundaryBuffers(FaceField &src, bool conserved_values)
{
  MeshBlock *pmb=pmy_mblock_;
  Coordinates *pco=pmb->pcoord;
  Mesh *pmesh=pmb->pmy_mesh;

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int ku,kl;
  if (pmesh->mesh_size.nx3>1) {
	ku = ke+NGHOST;
	kl = ks-NGHOST;
  } else {
	ku = ke;
	kl = ks;
  }

  Real qomL = qshear_*Omega_0_*x1size_;

  if (shbb_.inner == true) {
    int ib = is-NGHOST;
	int ii;
    // step 1. -- load shboxvar_field_
    for(int k=kl; k<=ku; k++) {
      for(int j=js-NGHOST; j<=je+NGHOST; j++) {
        for(int i=0; i<NGHOST; i++) {
          ii = ib+i;
          shboxvar_inner_field_.x1f(k,j,i) = src.x1f(k,j,ii);
          shboxvar_inner_field_.x2f(k,j,i) = src.x2f(k,j,ii);
          shboxvar_inner_field_.x3f(k,j,i) = src.x3f(k,j,ii);
	    }
    }}
	// fill the extra cells for B2i and B3i
	int kp = ku+1;
    for(int j=js-NGHOST; j<=je+NGHOST; j++) {
      for(int i=0; i<NGHOST; i++) {
        ii = ib+i;
	    shboxvar_inner_field_.x3f(kp,j,i) = src.x3f(kp,j,ii);
	}}
	int jp = je+NGHOST+1;
	for(int k=kl; k<=ku; k++) {
      for(int i=0; i<NGHOST; i++) {
        ii = ib+i;
        shboxvar_inner_field_.x2f(k,jp,i) = src.x2f(k,jp,ii);
	}}

    // step 2. -- conservative remapping
    for(int k=kl; k<=ku; k++) {  // bx1
      for(int i=0; i<NGHOST; i++) {
        RemapFluxField(k,js,je+2,i,eps_,shboxvar_inner_field_.x1f,flx_inner_field_.x1f);
	    for(int j=js; j<=je+1; j++) {
	      shboxvar_inner_field_.x1f(k,j,i) -= flx_inner_field_.x1f(j+1)-flx_inner_field_.x1f(j);
	    }
	}}
    for(int k=kl; k<=ku; k++) {  // bx2
      for(int i=0; i<NGHOST; i++) {
        RemapFluxField(k,js,je+3,i,eps_,shboxvar_inner_field_.x2f,flx_inner_field_.x2f);
	    for(int j=js; j<=je+2; j++) {
	      shboxvar_inner_field_.x2f(k,j,i) -= flx_inner_field_.x2f(j+1)-flx_inner_field_.x2f(j);
	    }
	}}
    for(int k=kl; k<=ku+1; k++) { // bx3
      for(int i=0; i<NGHOST; i++) {
        RemapFluxField(k,js,je+2,i,eps_,shboxvar_inner_field_.x3f,flx_inner_field_.x3f);
	    for(int j=js; j<=je+1; j++) {
	      shboxvar_inner_field_.x3f(k,j,i) -= flx_inner_field_.x3f(j+1)-flx_inner_field_.x3f(j);
	    }
	}}

	// step 3. -- load sendbuf; memcpy to recvbuf if on same rank, post MPI_Isend otherwise
	for(int n=0; n<4; n++) {
	  if(send_inner_rank_[n] != -1) {
        LoadFieldShearing(shboxvar_inner_field_, send_innerbuf_field_[n], n);
	    if (send_inner_rank_[n] == Globals::my_rank) {// on the same process
          MeshBlock *pbl=pmb->pmy_mesh->FindMeshBlock(send_inner_gid_[n]);
          std::memcpy(pbl->pbval->recv_innerbuf_field_[n],
                  send_innerbuf_field_[n], send_innersize_field_[n]*sizeof(Real));
	      pbl->pbval->shbox_inner_field_flag_[n]=BNDRY_ARRIVED;
		} else { // MPI
#ifdef MPI_PARALLEL
          int tag=CreateBvalsMPITag(send_inner_lid_[n], TAG_SHBOX_FIELD, n); //bufid = n
          MPI_Isend(send_innerbuf_field_[n],send_innersize_field_[n],MPI_ATHENA_REAL,send_inner_rank_[n],tag,MPI_COMM_WORLD, &rq_innersend_field_[n]);
#endif
	   }
	 }}
  } // inner boundaries

  if (shbb_.outer == true) {
	int	ib = ie+1;
	qomL = -qomL;
	int ii;
    // step 1. -- load shboxvar_field_
    for(int k=kl; k<=ku; k++) {
      for(int j=js-NGHOST; j<=je+NGHOST; j++) {
        for(int i=0; i<NGHOST; i++) {
          ii = ib+i;
          shboxvar_outer_field_.x1f(k,j,i) = src.x1f(k,j,ii+1);
          shboxvar_outer_field_.x2f(k,j,i) = src.x2f(k,j,ii);
          shboxvar_outer_field_.x3f(k,j,i) = src.x3f(k,j,ii);
	    }
    }}
	// fill the extra cells for B2i and B3i
	int kp = ku+1;
    for(int j=js-NGHOST; j<=je+NGHOST; j++) {
      for(int i=0; i<NGHOST; i++) {
        ii = ib+i;
	    shboxvar_outer_field_.x3f(kp,j,i) = src.x3f(kp,j,ii);
	}}
	int jp = je+NGHOST+1;
	for(int k=kl; k<=ku; k++) {
      for(int i=0; i<NGHOST; i++) {
        ii = ib+i;
        shboxvar_outer_field_.x2f(k,jp,i) = src.x2f(k,jp,ii);
	}}

    // step 2. -- conservative remapping
    for(int k=kl; k<=ku; k++) {  // bx1
      for(int i=0; i<NGHOST; i++) {
        RemapFluxField(k,js-1,je+1,i,-eps_,shboxvar_outer_field_.x1f,flx_outer_field_.x1f);
	    for(int j=js-1; j<=je; j++) {
	      shboxvar_outer_field_.x1f(k,j,i) -= flx_outer_field_.x1f(j+1)-flx_outer_field_.x1f(j);
	    }
	}}
    for(int k=kl; k<=ku; k++) {  // bx2
      for(int i=0; i<NGHOST; i++) {
        RemapFluxField(k,js-1,je+2,i,-eps_,shboxvar_outer_field_.x2f,flx_outer_field_.x2f);
	    for(int j=js-1; j<=je+1; j++) {
	      shboxvar_outer_field_.x2f(k,j,i) -= flx_outer_field_.x2f(j+1)-flx_outer_field_.x2f(j);
	    }
	}}
    for(int k=kl; k<=ku+1; k++) {  // bx3
      for(int i=0; i<NGHOST; i++) {
        RemapFluxField(k,js-1,je+1,i,-eps_,shboxvar_outer_field_.x3f,flx_outer_field_.x3f);
	    for(int j=js-1; j<=je; j++) {
	      shboxvar_outer_field_.x3f(k,j,i) -= flx_outer_field_.x3f(j+1)-flx_outer_field_.x3f(j);
	    }
	}}

	// step 3. -- load sendbuf; memcpy to recvbuf if on same rank, post MPI_Isend otherwise
	int offset = 4;
	for(int n=0; n<4; n++) {
	  if(send_outer_rank_[n] != -1) {
        LoadFieldShearing(shboxvar_outer_field_, send_outerbuf_field_[n], n+offset);
	    if (send_outer_rank_[n] == Globals::my_rank) {// on the same process
          MeshBlock *pbl=pmb->pmy_mesh->FindMeshBlock(send_outer_gid_[n]);
          std::memcpy(pbl->pbval->recv_outerbuf_field_[n],
                  send_outerbuf_field_[n], send_outersize_field_[n]*sizeof(Real));
	      pbl->pbval->shbox_outer_field_flag_[n]=BNDRY_ARRIVED;
	    } else { // MPI
#ifdef MPI_PARALLEL
          int tag=CreateBvalsMPITag(send_outer_lid_[n], TAG_SHBOX_FIELD, n+offset); //bufid for outer(inner): 2(0) and 3(1)
          MPI_Isend(send_outerbuf_field_[n],send_outersize_field_[n],MPI_ATHENA_REAL,send_outer_rank_[n],tag,MPI_COMM_WORLD, &rq_outersend_field_[n]);
#endif
	    }
      }}
  } // outer boundaries
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetFieldShearingboxBoundarySameLevel(FaceField &dst,
//                                           Real *buf, const int nb)
//  \brief Set field shearingbox boundary received from a block on the same level
void BoundaryValues::SetFieldShearingboxBoundarySameLevel(FaceField &dst, Real *buf,
                                               const int nb)
{
  MeshBlock *pmb=pmy_mblock_;
  Mesh *pmesh=pmb->pmy_mesh;
  int si, sj, sk, ei, ej, ek;
  int psi,pei,psj,pej;
  int nx2=pmb->block_size.nx2-1;

  sk = pmb->ks; ek = pmb->ke;
  if (pmesh->mesh_size.nx3>1) ek += NGHOST, sk -= NGHOST;
  // --- inner boundary:
  // set [js-NGHOST:js+(joverlap-1)] with recv_inner[0]; nb=0
  // set [js-1:js+(joverlap-1)]      with recv_inner[0]; nb=0 if joverlap==nx2-1
  // set [js+joverlap:je+NGHOST]     with recv_inner[1]; nb=1
  //    ([js+joverlap:je+NGHOST+1]   for bx2)
  // set [js+joverlap:je]            with recv_inner[1]; nb=1 if joverlap==0
  //    ([js+joverlap:je+1]          for bx2)
  // set [js-NGHOST:js-1]            with recv_inner[2]; nb=2
  // set [js-NGHOST:js-NGHOST]       with recv_inner[2]; nb=2 if joverlap==nx2-1
  // set [je+1:je+NGHOST]            with recv_inner[3]; nb=3
  //    ([je+2:je+NGHOST+1]          for bx2)
  // --- outer boundary:
  // set [je-joverlap+1:je+NGHOST]   with recv_outer[0]; nb=4
  //    ([je-joverlap+1:je+NGHOST+1] for bx2)
  // set [je-joverlap+1:je+1]        with recv_outer[0]; nb=4 if joverlap==nx2-1
  //    ([je-joverlap+1:je+2]        for bx2)
  // set [js-NGHOST:je-joverlap]     with recv_outer[1]; nb=5
  // set [js:je-joverlap]            with recv_outer[1]; nb=5 if joverlap==0
  //    ([js:je-joverlap+1]          for bx2 if joverlap==0)
  // set [je+1:je+NGHOST]            with recv_outer[2]; nb=6
  //    ([je+2:je+NGHOST+1]          for bx2)
  // set [je+NGHOST:je+NGHOST]       with recv_outer[2]; nb=6 if joverlap==nx2-1
  //    ([je+NGHOST+1:je+NGHOST+1]   for bx2)
  // set [js-NGHOST:js-1]            with recv_outer[3]; nb=7
  switch(nb) {
	case 0:
      si=pmb->is-NGHOST; ei=pmb->is-1; sj=pmb->js-NGHOST; ej=pmb->js+(joverlap_-1);
      if(joverlap_==nx2) sj=pmb->js-1;
	  psi=si; pei=ei; psj=sj; pej=ej;
	  break;
	case 1:
      si=pmb->is-NGHOST; ei=pmb->is-1; sj=pmb->js+joverlap_; ej=pmb->je+NGHOST;
      if(joverlap_==0)   ej=pmb->je;
	  psi=si; pei=ei; psj=sj; pej=ej+1;
	  break;
	case 2:
      si=pmb->is-NGHOST; ei=pmb->is-1; sj=pmb->js-NGHOST; ej=pmb->js-1;
      if(joverlap_==nx2) ej=pmb->js-NGHOST;
	  psi=si; pei=ei; psj=sj; pej=ej;
	  break;
	case 3:
      si=pmb->is-NGHOST; ei=pmb->is-1; sj=pmb->je+1; ej=pmb->je+NGHOST;
	  psi=si; pei=ei; psj=sj+1; pej=ej+1;
	  break;
	case 4:
      si=pmb->ie+1; ei=pmb->ie+NGHOST; sj=pmb->je-(joverlap_-1); ej=pmb->je+NGHOST;
	  if(joverlap_==nx2) ej=pmb->je+1;
	  psi=si+1; pei=ei+1; psj=sj; pej=ej+1;
	  break;
	case 5:
      si=pmb->ie+1; ei=pmb->ie+NGHOST; sj=pmb->js-NGHOST; ej=pmb->je-joverlap_;
	  psi=si+1; pei=ei+1; psj=sj; pej=ej;
	  if(joverlap_==0) {
		sj=pmb->js;
		psj=sj;
		pej=ej+1;
	  }
	  break;
	case 6:
      si=pmb->ie+1; ei=pmb->ie+NGHOST; sj=pmb->je+1; ej=pmb->je+NGHOST;
	  if(joverlap_==nx2)   sj=pmb->je+NGHOST;
	  psi=si+1; pei=ei+1; psj=sj+1; pej=ej+1;
	  break;
	case 7:
      si=pmb->ie+1; ei=pmb->ie+NGHOST; sj=pmb->js-NGHOST; ej=pmb->js-1;
	  psi=si+1; pei=ei+1; psj=sj; pej=ej;
	  break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues:SetFieldShearing " << std::endl
          << "nb = " << nb << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
   }

  // set [sj:ej] of current meshblock
  int p=0;
  BufferUtility::Unpack3DData(buf, dst.x1f, psi, pei, sj, ej, sk, ek, p);
  BufferUtility::Unpack3DData(buf, dst.x2f, si, ei, psj, pej, sk, ek, p);
  BufferUtility::Unpack3DData(buf, dst.x3f, si, ei, sj, ej, sk, ek+1, p);
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ReceiveFieldShearingboxBoundaryBuffersWithWait(FaceField &dst, bool conserved_value)
//  \brief receive shearingbox boundary data for field(face-centered) variables for initialization
//         It is not used for now as our problem generator might start with t=0 periodic config.
void BoundaryValues::ReceiveFieldShearingboxBoundaryBuffersWithWait(FaceField &dst, bool conserved_value)
{
  MeshBlock *pmb=pmy_mblock_;

  if(shbb_.inner == true) { // check inner boundaries
#ifdef MPI_PARALLEL
	if (recv_inner_rank_[1]!=Globals::my_rank)
      MPI_Wait(&rq_innerrecv_field_[1],MPI_STATUS_IGNORE);
#endif
    SetFieldShearingboxBoundarySameLevel(dst,recv_innerbuf_field_[1],1);
    shbox_inner_field_flag_[1] = BNDRY_COMPLETED;
  } // inner boundary

  if(shbb_.outer == true) { // check inner boundaries
#ifdef MPI_PARALLEL
	if (recv_outer_rank_[1]!=Globals::my_rank)
      MPI_Wait(&rq_outerrecv_field_[1],MPI_STATUS_IGNORE);
#endif
    SetFieldShearingboxBoundarySameLevel(dst,recv_outerbuf_field_[1],5);
    shbox_outer_field_flag_[1] = BNDRY_COMPLETED;
  } // outer boundary


  return;

}
//--------------------------------------------------------------------------------------
//! \fn bool BoundaryValues::ReceiveFieldShearingboxBoundaryBuffers(FaceField &dst)
//  \brief receive shearingbox boundary data for field(face-centered) variables
bool BoundaryValues::ReceiveFieldShearingboxBoundaryBuffers(FaceField &dst)
{
  MeshBlock *pmb=pmy_mblock_;
  Mesh *pmesh=pmb->pmy_mesh;
  bool flagi=true, flago=true;

  if(shbb_.inner == true) { // check inner boundaries
	for(int n=0; n<4; n++) {
      if(shbox_inner_field_flag_[n]==BNDRY_COMPLETED) continue;
	  if(shbox_inner_field_flag_[n]==BNDRY_WAITING) {
	    if (recv_inner_rank_[n]==Globals::my_rank) {// on the same process
          flagi=false;
		  continue;
		}
        else { // MPI boundary
#ifdef MPI_PARALLEL
          int test;
          MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&test,MPI_STATUS_IGNORE);
          MPI_Test(&rq_innerrecv_field_[n],&test,MPI_STATUS_IGNORE);
          if(test==false) {
            flagi=false;
            continue;
          }
		  shbox_inner_field_flag_[n] = BNDRY_ARRIVED;
#endif
	    }
      }
	  // set dst if boundary arrived
      SetFieldShearingboxBoundarySameLevel(dst,recv_innerbuf_field_[n],n);
      shbox_inner_field_flag_[n] = BNDRY_COMPLETED; // completed
    } // loop over recv[0] to recv[3]
  } // inner boundary

  if(shbb_.outer == true) { // check outer boundaries
	int offset = 4;
	for(int n=0; n<4; n++) {
      if(shbox_outer_field_flag_[n]==BNDRY_COMPLETED) continue;
	  if(shbox_outer_field_flag_[n]==BNDRY_WAITING) {
	    if (recv_outer_rank_[n]==Globals::my_rank) {// on the same process
          flago=false;
		  continue;
		}
        else { // MPI boundary
#ifdef MPI_PARALLEL
          int test;
          MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&test,MPI_STATUS_IGNORE);
          MPI_Test(&rq_outerrecv_field_[n],&test,MPI_STATUS_IGNORE);
          if(test==false) {
            flago=false;
            continue;
          }
		  shbox_outer_field_flag_[n] = BNDRY_ARRIVED;
#endif
	    }
      }
      SetFieldShearingboxBoundarySameLevel(dst,recv_outerbuf_field_[n],n+offset);
      shbox_outer_field_flag_[n] = BNDRY_COMPLETED; // completed
    } // loop over recv[0] and recv[1]
  } // outer boundary

  return (flagi && flago);

}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::RemapFluxField(int k, int jinner, int jouter, int i, Real eps, static AthenaArray<Real> &U, AthenaArray<Real> &Flux)
//  \brief compute the flux along j indices for remapping
//  adopted from 2nd order RemapFlux of Athena4.0
void BoundaryValues::RemapFluxField(const int k, const int jinner, const int jouter, const int i, const Real eps, const AthenaArray<Real> &U, AthenaArray<Real> &Flux)
{
  int j,jl,ju;
  Real dUc,dUl,dUr,dUm,lim_slope;

/* jinner,jouter are index range over which flux must be returned.  Set loop
 * limits depending on direction of upwind differences  */

  if (eps > 0.0) { /* eps always > 0 for inner i boundary */
    jl = jinner-1;
    ju = jouter-1;
  } else {         /* eps always < 0 for outer i boundary */
    jl = jinner;
    ju = jouter;
  }

  for (j=jl; j<=ju; j++) {
      dUc = U(k,j+1,i) - U(k,j-1,i);
      dUl = U(k,j,  i) - U(k,j-1,i);
      dUr = U(k,j+1,i) - U(k,j,  i);

      dUm = 0.0;
      if (dUl*dUr > 0.0) {
        lim_slope = std::min(fabs(dUl),fabs(dUr));
        dUm = SIGN(dUc)*std::min(0.5*fabs(dUc),2.0*lim_slope);
      }

    if (eps > 0.0) { /* eps always > 0 for inner i boundary */
      Flux(j+1) = eps*(U(k,j,i) + 0.5*(1.0 - eps)*dUm);
    } else {         /* eps always < 0 for outer i boundary */
      Flux(j  ) = eps*(U(k,j,i) - 0.5*(1.0 + eps)*dUm);
    }
  }

  return;
}

