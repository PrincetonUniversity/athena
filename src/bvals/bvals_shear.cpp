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
//! \fn int BoundaryValues::LoadHydroShearing(AthenaArray<Real> &src, Real *buf, int nb)
//  \brief Load shearingbox hydro boundary buffers
void BoundaryValues::LoadHydroShearing(AthenaArray<Real> &src, Real *buf, int nb)
{
  MeshBlock *pmb=pmy_mblock_;
  Mesh *pmesh=pmb->pmy_mesh;
  int si, sj, sk, ei, ej, ek;

  sk = pmb->ks;
  ek = pmb->ke;
  //if (pmesh->mesh_size.nx3>1) ek = pmb->ke+1; //one more for Bz
  if (pmesh->mesh_size.nx3>1) {
	ek += NGHOST;
	sk -= NGHOST;
  }
  // load [je-(joverlap-1):je] for send_inner[0]; nb=0
  // load [js:je-joverlap]     for send_inner[1]; nb=1
  // load [js:js+(joverlap-1)] for send_outer[0]; nb=2
  // load [js+joverlap:je]     for send_outer[1]; nb=3
  if(nb==0)  si=pmb->is-NGHOST, ei=pmb->is-1, sj=pmb->je-(joverlap_-1), ej=pmb->je;
  if(nb==1)  si=pmb->is-NGHOST, ei=pmb->is-1, sj=pmb->js,               ej=pmb->je-joverlap_;
  if(nb==2)  si=pmb->is-NGHOST, ei=pmb->is-1, sj=pmb->js,               ej=pmb->js+(joverlap_-1);
  if(nb==3)  si=pmb->is-NGHOST, ei=pmb->is-1, sj=pmb->js+joverlap_,     ej=pmb->je;
  //if(nb==2)  si=pmb->ie+1, ei=pmb->ie+NGHOST, sj=pmb->js,               ej=pmb->js+(joverlap_-1);
  //if(nb==3)  si=pmb->ie+1, ei=pmb->ie+NGHOST, sj=pmb->js+joverlap_,     ej=pmb->je;

  int p=0;
  BufferUtility::Pack4DData(src, buf, 0, NHYDRO-1, si, ei, sj, ej, sk, ek, p);

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SendHydroShearingboxBoundaryBuffersForInit(AthenaArray<Real> &src,
//                                                    bool conserved_values)
//  \brief Send shearingbox boundary buffers for hydro variables

void BoundaryValues::SendHydroShearingboxBoundaryBuffersForInit(AthenaArray<Real> &src,
                                              bool conserved_values)
{
  //std::cout << "do some SendHydroShearingboxBoundaryBuffersForInit\n" << std::endl;

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
  int ssize = ssize_*NHYDRO;

  if (shbb_.inner == true) {
    // step 1. -- add shear to the inner periodic boundary values
    for(int k=kl; k<=ku; k++) {
      for(int j=js-NGHOST; j<=je+NGHOST; j++) {
        for(int i=0; i<NGHOST; i++) {
          //shboxvar_inner_hydro_(IDN,k,j,i) = src(IDN,k,j,ii);
          //shboxvar_inner_hydro_(IM1,k,j,i) = src(IM1,k,j,ii);
          //shboxvar_inner_hydro_(IM3,k,j,i) = src(IM3,k,j,ii);
          shboxvar_inner_hydro_(IM2,k,j,i) = src(IM2,k,j,i) + qomL*src(IDN,k,j,i);//add shear to conservative
	      if (NON_BAROTROPIC_EOS) {
	        src(IEN,k,j,i) += (0.5/src(IDN,k,j,i))*
			(SQR(shboxvar_inner_hydro_(IM2,k,j,i)) - SQR(src(IM2,k,j,i)));
	      } // update energy
          src(IM2,k,j,i) = shboxvar_inner_hydro_(IM2,k,j,i);//update IM2
	    }
    }}
  }

  if (shbb_.outer == true) {
	int	ib = ie+1;
	int ii;
    // step 2. -- add shear to the outer periodic boundary values
    for(int k=kl; k<=ku; k++) {
      for(int j=js-NGHOST; j<=je+NGHOST; j++) {
        for(int i=0; i<NGHOST; i++) {
		  ii = ib+i;
          //shboxvar_inner_hydro_(IDN,k,j,i) = src(IDN,k,j,ii);
          //shboxvar_inner_hydro_(IM1,k,j,i) = src(IM1,k,j,ii);
          //shboxvar_inner_hydro_(IM3,k,j,i) = src(IM3,k,j,ii);
          shboxvar_outer_hydro_(IM2,k,j,i) = src(IM2,k,j,ii) - qomL*src(IDN,k,j,ii);//add shear to conservative
	      if (NON_BAROTROPIC_EOS) {
	        src(IEN,k,j,ii) += (0.5/src(IDN,k,j,ii))*
			(SQR(shboxvar_inner_hydro_(IM2,k,j,i)) - SQR(src(IM2,k,j,ii)));
	      } // update energy
          src(IM2,k,j,ii) = shboxvar_outer_hydro_(IM2,k,j,i);//update IM2
	    }
    }}
  }
  return;
}
//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SendHydroShearingboxBoundaryBuffers(AthenaArray<Real> &src,
//                                                    bool conserved_values)
//  \brief Send shearingbox boundary buffers for hydro variables

void BoundaryValues::SendHydroShearingboxBoundaryBuffers(AthenaArray<Real> &src,
                                              bool conserved_values)
{
  //std::cout << "do some SendHydroShearingboxBoundaryBuffers\n" << std::endl;

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
  int ssize = ssize_*NHYDRO;

  if (shbb_.inner == true) {
    int ib = is-NGHOST;
	int ii;
    // step 1. -- load shboxvar_hydro_
    for(int k=kl; k<=ku; k++) {
      for(int j=js-NGHOST; j<=je+NGHOST; j++) {
        for(int i=0; i<NGHOST; i++) {
          ii = ib+i;
          shboxvar_inner_hydro_(IDN,k,j,i) = src(IDN,k,j,ii);
          shboxvar_inner_hydro_(IM1,k,j,i) = src(IM1,k,j,ii);
          shboxvar_inner_hydro_(IM2,k,j,i) = src(IM2,k,j,ii) + qomL*src(IDN,k,j,ii);//add shear to conservative
          shboxvar_inner_hydro_(IM3,k,j,i) = src(IM3,k,j,ii);
	      if (NON_BAROTROPIC_EOS) {
	        shboxvar_inner_hydro_(IEN,k,j,i) = src(IEN,k,j,ii) + (0.5/src(IDN,k,j,ii))*
			(SQR(shboxvar_inner_hydro_(IM2,k,j,i)) - SQR(src(IM2,k,j,ii)));
	      }
	    }
    }}

    // step 2. -- conservative remaping
    for(int n=0; n<NHYDRO; n++) {
      for(int k=kl; k<=ku; k++) {
        for(int i=0; i<NGHOST; i++) {
          RemapFlux(n,k,js,je+1,i,eps_,shboxvar_inner_hydro_,flx_inner_hydro_);
	      for(int j=js; j<=je; j++) {
	        shboxvar_inner_hydro_(n,k,j,i) -= flx_inner_hydro_(j+1)-flx_inner_hydro_(j);
		  }
		}
	}}
	// debug:
	//if (i==0 &&  k==ks && n==2) {
	//  for(int j=js-NGHOST; j<=je+NGHOST; j++){
	//    std::cout << "outer step 3:[i,j,k] = "<< i+ie+1 << " " << j << " " << k << " " << "GhstZns(IM2,k,i,j), GhstZnsBuf(IM2,k,i,j)= " <<
	//	  GhstZns(n,k,i,j) << " " << GhstZnsBuf(n,k,i,j) << std::endl;
	//  }
	//}

    // step 3. load send_buf_hydro_ buffers
    LoadHydroShearing(shboxvar_inner_hydro_,send_innerbuf_hydro_[1],1);
    //LoadHydroShearing(shboxvar_inner_hydro_,send_innerbuf_hydro_[0],send_innerbuf_hydro_[1],joverlap_);

	// step 4. send the buffers: memcpy if same rank; mpi call if not.
	// step 4a. send_innerbuf_hydro_[1]
	if (send_inner_rank_[1] == Globals::my_rank) { // on the same process
      MeshBlock *pbl=pmb->pmy_mesh->FindMeshBlock(send_inner_gid_[1]);
      std::memcpy(pbl->pbval->recv_innerbuf_hydro_[1],
                  send_innerbuf_hydro_[1], send_innersize_hydro_[1]*ssize*sizeof(Real));
	  pbl->pbval->shbox_inner_hydro_flag_[1]=BNDRY_ARRIVED;
	}
	else { // MPI
#ifdef MPI_PARALLEL
      int tag=CreateBvalsMPITag(send_inner_lid_[1], TAG_SHBOX_HYDRO, 1); //bufid for outer(inner): 2(0) and 3(1)
      //int tag=CreateBvalsMPITag(pmb->lid, TAG_SHBOX_HYDRO, 1); //bufid for outer(inner): 2(0) and 3(1)
      MPI_Isend(send_innerbuf_hydro_[1],send_innersize_hydro_[1]*ssize,MPI_ATHENA_REAL,send_inner_rank_[1],tag,MPI_COMM_WORLD, &rq_innersend_hydro_[1]);
	  std::cout << "on cycle = " << pmesh->ncycle << "gid = " << pmb->gid << "send_inner[1] tag= " << tag << "(" << Globals::my_rank << "->" << send_inner_rank_[1] << ")" << std::endl;
#endif
	}
	// step 4b. send_innerbuf_hydro_[0] if overlap (rank=-1 if no overlap)
	if (send_inner_rank_[0] != -1) {
      LoadHydroShearing(shboxvar_inner_hydro_,send_innerbuf_hydro_[0],0);
	  if (send_inner_rank_[0] == Globals::my_rank) { // on the same process
         MeshBlock *pbl=pmb->pmy_mesh->FindMeshBlock(send_inner_gid_[0]);
         std::memcpy(pbl->pbval->recv_innerbuf_hydro_[0],
                     send_innerbuf_hydro_[0], send_innersize_hydro_[0]*ssize*sizeof(Real));
	    pbl->pbval->shbox_inner_hydro_flag_[0]=BNDRY_ARRIVED;
	  }
	  else { // MPI
#ifdef MPI_PARALLEL
         int tag=CreateBvalsMPITag(send_inner_lid_[0], TAG_SHBOX_HYDRO, 0);
         //int tag=CreateBvalsMPITag(pmb->lid, TAG_SHBOX_HYDRO, 0);
         MPI_Isend(send_innerbuf_hydro_[0],send_innersize_hydro_[0]*ssize,MPI_ATHENA_REAL,send_inner_rank_[0],tag,MPI_COMM_WORLD,&rq_innersend_hydro_[0]);
	  std::cout << "on cycle = " << pmesh->ncycle << "gid = " << pmb->gid << "send_inner[0] tag= " << tag << "(" << Globals::my_rank << "->" << send_inner_rank_[0] << ")" << std::endl;
#endif
	  }
    }
  } // inner boundaries

  if (shbb_.outer == true) {
	int	ib = ie+1;
	qomL = -qomL;
	int ii;
    // step 1. -- load shboxvar_hydro_
    for(int k=kl; k<=ku; k++) {
      for(int j=js-NGHOST; j<=je+NGHOST; j++) {
        for(int i=0; i<NGHOST; i++) {
          ii = ib+i;
          shboxvar_outer_hydro_(IDN,k,j,i) = src(IDN,k,j,ii);
          shboxvar_outer_hydro_(IM1,k,j,i) = src(IM1,k,j,ii);
          shboxvar_outer_hydro_(IM2,k,j,i) = src(IM2,k,j,ii) + qomL*src(IDN,k,j,ii);//add shear to conservative
          shboxvar_outer_hydro_(IM3,k,j,i) = src(IM3,k,j,ii);
	      if (NON_BAROTROPIC_EOS) {
	        shboxvar_outer_hydro_(IEN,k,j,i) = src(IEN,k,j,ii) + (0.5/src(IDN,k,j,ii))*
			(SQR(shboxvar_outer_hydro_(IM2,k,j,i)) - SQR(src(IM2,k,j,ii)));
	      }
	    }
    }}

    // step 2. -- conservative remaping
    for(int n=0; n<NHYDRO; n++) {
      for(int k=kl; k<=ku; k++) {
        for(int i=0; i<NGHOST; i++) {
          RemapFlux(n,k,js,je+1,i,-eps_,shboxvar_outer_hydro_,flx_outer_hydro_);
	      for(int j=js; j<=je; j++) {
	        shboxvar_outer_hydro_(n,k,j,i) -= flx_outer_hydro_(j+1)-flx_outer_hydro_(j);
		  }
		}
	}}

    // step 3. load send_buf_hydro_ buffers
    LoadHydroShearing(shboxvar_outer_hydro_,send_outerbuf_hydro_[1],3);
    //LoadHydroShearing(shboxvar_outer_hydro_,send_outerbuf_hydro_[0],send_outerbuf_hydro_[1],joverlap_);

	// step 4. send the buffers: memcpy if same rank; mpi call if not.
	// step 4a. send_outerbuf_hydro_[1]
	if (send_outer_rank_[1] == Globals::my_rank) { // on the same process
      MeshBlock *pbl=pmb->pmy_mesh->FindMeshBlock(send_outer_gid_[1]);
      std::memcpy(pbl->pbval->recv_outerbuf_hydro_[1],
                  send_outerbuf_hydro_[1], send_outersize_hydro_[1]*ssize*sizeof(Real));
	  pbl->pbval->shbox_outer_hydro_flag_[1]=BNDRY_ARRIVED;
	}
	else { // MPI
#ifdef MPI_PARALLEL
      int tag=CreateBvalsMPITag(send_outer_lid_[1], TAG_SHBOX_HYDRO, 3); //bufid for outer(inner): 2(0) and 3(1)
      //int tag=CreateBvalsMPITag(pmb->lid, TAG_SHBOX_HYDRO, 3); //bufid for outer(inner): 2(0) and 3(1)
      MPI_Isend(send_outerbuf_hydro_[1],send_outersize_hydro_[1]*ssize,MPI_ATHENA_REAL,send_outer_rank_[1],tag,MPI_COMM_WORLD, &rq_outersend_hydro_[1]);
	  std::cout << "on cycle = " << pmesh->ncycle << "gid = " << pmb->gid << "send_outer[1] tag= " << tag << "(" << Globals::my_rank << "->" << send_outer_rank_[1] << ")" << std::endl;
#endif
	}
	// step 4b. send_outerbuf_hydro_[0] if overlap (rank=-1 if no overlap)
	if (send_outer_rank_[0] != -1) {
      LoadHydroShearing(shboxvar_outer_hydro_,send_outerbuf_hydro_[0],2);
	  if (send_outer_rank_[0] == Globals::my_rank) { // on the same process
         MeshBlock *pbl=pmb->pmy_mesh->FindMeshBlock(send_outer_gid_[0]);
         std::memcpy(pbl->pbval->recv_outerbuf_hydro_[0],
                     send_outerbuf_hydro_[0], send_outersize_hydro_[0]*ssize*sizeof(Real));
	     pbl->pbval->shbox_outer_hydro_flag_[0]=BNDRY_ARRIVED;
	  }
	  else { // MPI
#ifdef MPI_PARALLEL
         int tag=CreateBvalsMPITag(send_outer_lid_[0], TAG_SHBOX_HYDRO, 2);
         //int tag=CreateBvalsMPITag(pmb->lid, TAG_SHBOX_HYDRO, 2);
         MPI_Isend(send_outerbuf_hydro_[0],send_outersize_hydro_[0]*ssize,MPI_ATHENA_REAL,send_outer_rank_[0],tag,MPI_COMM_WORLD,&rq_outersend_hydro_[0]);
	  std::cout << "on cycle = " << pmesh->ncycle << "gid = " << pmb->gid << "send_outer[0] tag= " << tag << "(" << Globals::my_rank << "->" << send_outer_rank_[0] << ")" << std::endl;
#endif
	  }
    }
  } // outer boundaries
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetHydroShearingboxBoundarySameLevel(AthenaArray<Real> &dst,
//                                           Real *buf, const int nb)
//  \brief Set hydro shearingbox boundary received from a block on the same level
void BoundaryValues::SetHydroShearingboxBoundarySameLevel(AthenaArray<Real> &dst, Real *buf,
                                               const int nb)
{
  MeshBlock *pmb=pmy_mblock_;
  Mesh *pmesh=pmb->pmy_mesh;
  int si, sj, sk, ei, ej, ek;

  sk = pmb->ks;
  ek = pmb->ke;
  if (pmesh->mesh_size.nx3>1) {
	ek += NGHOST;
	sk -= NGHOST;
  }
  // set [js:js+(joverlap-1)] of inner meshblock with recv_inner[0];nb=0
  // set [js+joverlap:je] of inner meshblock with     recv_inner[1];nb=1
  // set [je-(joverlap-1):je] of outer meshblock with recv_outer[0];nb=2
  // set [js:je-joverlap] of outer meshblock with     recv_outer[1];nb=3
  if(nb==0)  si=pmb->is-NGHOST, ei=pmb->is-1, sj=pmb->js,               ej=pmb->js+(joverlap_-1);
  if(nb==1)  si=pmb->is-NGHOST, ei=pmb->is-1, sj=pmb->js+joverlap_,     ej=pmb->je;
  if(nb==2)  si=pmb->ie+1, ei=pmb->ie+NGHOST, sj=pmb->je-(joverlap_-1), ej=pmb->je;
  if(nb==3)  si=pmb->ie+1, ei=pmb->ie+NGHOST, sj=pmb->js,               ej=pmb->je-joverlap_;

  // set [sj:ej] of current meshblock
  int p=0;
  BufferUtility::Unpack4DData(buf, dst, 0, NHYDRO-1, si, ei, sj, ej, sk, ek, p);
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ReceiveHydroShearingboxBoundaryBuffersWithWait(AthenaArray<Real> &dst, bool conserved_value)
//  \brief receive shearingbox boundary data for hydro variables for initialization
void BoundaryValues::ReceiveHydroShearingboxBoundaryBuffersWithWait(AthenaArray<Real> &dst, bool conserved_value)
{
  MeshBlock *pmb=pmy_mblock_;
  //std::cout << "do some ReceiveHydroShearingboxBoundaryBuffersWithWait\n" << std::endl;

  if(shbb_.inner) { // check inner boundaries
#ifdef MPI_PARALLEL
	if (recv_inner_rank_[1]!=Globals::my_rank)
      MPI_Wait(&rq_innerrecv_hydro_[1],MPI_STATUS_IGNORE);
#endif
    SetHydroShearingboxBoundarySameLevel(dst,recv_innerbuf_hydro_[1],1);
    shbox_inner_hydro_flag_[1] = BNDRY_COMPLETED;
  } // inner boundary

  if(shbb_.outer) { // check inner boundaries
#ifdef MPI_PARALLEL
	if (recv_outer_rank_[1]!=Globals::my_rank)
      MPI_Wait(&rq_outerrecv_hydro_[1],MPI_STATUS_IGNORE);
#endif
    SetHydroShearingboxBoundarySameLevel(dst,recv_outerbuf_hydro_[1],3);
    shbox_outer_hydro_flag_[1] = BNDRY_COMPLETED;
  } // outer boundary


  return;

}
//--------------------------------------------------------------------------------------
//! \fn bool BoundaryValues::ReceiveHydroShearingboxBoundaryBuffers(AthenaArray<Real> &dst)
//  \brief receive shearingbox boundary data for hydro variables
bool BoundaryValues::ReceiveHydroShearingboxBoundaryBuffers(AthenaArray<Real> &dst)
{
  MeshBlock *pmb=pmy_mblock_;
  Mesh *pmesh=pmb->pmy_mesh;
  bool flagi=true, flago=true;
  //std::cout << "do some ReceiveHydroShearingboxBoundaryBuffers\n" << std::endl;

  if(shbb_.inner) { // check inner boundaries
	for(int n=0; n<2; n++) {
      if(shbox_inner_hydro_flag_[n]==BNDRY_COMPLETED) continue;
	  if(shbox_inner_hydro_flag_[n]==BNDRY_WAITING) {
	    if (recv_inner_rank_[n]==Globals::my_rank) {// on the same process
          flagi=false;
		  continue;
		}
#ifdef MPI_PARALLEL
        else { // MPI boundary
          int test;
          MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&test,MPI_STATUS_IGNORE);
          MPI_Test(&rq_innerrecv_hydro_[n],&test,MPI_STATUS_IGNORE);
          if(test==false) {
            flagi=false;
            continue;
          }
	      std::cout << "on cycle = " << pmesh->ncycle << "gid = " << pmb->gid << " rank= " << Globals::my_rank << "rq_innerrecv[" << n << "] test= " << test << std::endl;
		  shbox_inner_hydro_flag_[n] = BNDRY_ARRIVED;
	    }
#endif
      }
	  // set dst if boundary arrived
    //if(nb.level==pmb->loc.level)
      //SetHydroShearingboxBoundarySameLevel(shboxvar_inner_hydro_[n],recv_innerbuf_hydro_[n],n);
      SetHydroShearingboxBoundarySameLevel(dst,recv_innerbuf_hydro_[n],n);
      shbox_inner_hydro_flag_[n] = BNDRY_COMPLETED; // completed
    } // loop over recv[0] and recv[1]
  } // inner boundary

  if(shbb_.outer) { // check outer boundaries
	int offset = 2;
	for(int n=0; n<2; n++) {
      if(shbox_outer_hydro_flag_[n]==BNDRY_COMPLETED) continue;
	  if(shbox_outer_hydro_flag_[n]==BNDRY_WAITING) {
	    if (recv_outer_rank_[n]==Globals::my_rank) {// on the same process
          flago=false;
		  continue;
		}
#ifdef MPI_PARALLEL
        else { // MPI boundary
          int test;
          MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&test,MPI_STATUS_IGNORE);
          MPI_Test(&rq_outerrecv_hydro_[n],&test,MPI_STATUS_IGNORE);
          if(test==false) {
            flago=false;
            continue;
          }
	      std::cout << "on cycle = " << pmesh->ncycle << "gid = " << pmb->gid << " rank= " << Globals::my_rank << "rq_outerrecv[" << n << "] test= " << test << std::endl;
		  shbox_outer_hydro_flag_[n] = BNDRY_ARRIVED;
	    }
#endif
      }
    //if(nb.level==pmb->loc.level)
      //SetHydroShearingboxBoundarySameLevel(shboxvar_outer_hydro_[n],recv_outerbuf_hydro_[n],n+offset);
      SetHydroShearingboxBoundarySameLevel(dst,recv_outerbuf_hydro_[n],n+offset);
      shbox_outer_hydro_flag_[n] = BNDRY_COMPLETED; // completed
    } // loop over recv[0] and recv[1]
  } // outer boundary

  return (flagi && flago);

}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::FindShearBlock(void)
//
//  \brief Calc the following things:
//  send_gid recv_gid send_lid recv_lid send_rank recv_rank,
//  send_size_hydro  recv_size_hydro: for MPI_Irecv
//  eps_,joverlap_: for update the conservative

void BoundaryValues::FindShearBlock(void)
{
  //std::cout << "update the send_to and get_from shearing blocks \n" << std::endl;
  MeshBlock *pmb=pmy_mblock_;
  Coordinates *pco=pmb->pcoord;
  Mesh *pmesh=pmb->pmy_mesh;

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int ku, ii,jj;

  int nrbx2 = pmesh->nrbx2; // # of meshblocks in azimuth
  int nx2   = pmb->block_size.nx2; // # of cells per meshblock

  Real qomL = qshear_*Omega_0_*x1size_;
  Real yshear = qomL*pmesh->time;
  Real deltay = fmod(yshear,x2size_);
  int joffset = (int)(deltay/pco->dx2v(js)); // this assumes uniform grid in azimuth
  int Ngrids  = (int)(joffset/nx2);
  joverlap_   = joffset - Ngrids*nx2;
  eps_ = (fmod(deltay,pco->dx2v(js)))/pco->dx2v(js);

  if (shbb_.inner == true) { // if inner block
	for (int n=0; n<2; n++){
      send_inner_gid_[n]  = -1;
	  send_inner_rank_[n] = -1;
	  send_inner_lid_[n]  = -1;
      recv_inner_gid_[n]  = -1;
	  recv_inner_rank_[n] = -1;
	  recv_inner_lid_[n]  = -1;
	  send_innersize_hydro_[n] = 0;
	  recv_innersize_hydro_[n] = 0;
	}
	int jblock;
    for (int j=0; j<nrbx2; j++) {
      if (shbb_.igidlist[j] == pmb->gid)  jblock = j;// index of current meshblock on the shearingboundary block list
    }
	// send [js:je-joverlap] of the meshblock to other
	int jtmp = jblock + Ngrids;
	if (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
    send_inner_gid_[1]  = shbb_.igidlist[jtmp];
	send_inner_rank_[1] = shbb_.irnklist[jtmp];
	send_inner_lid_[1]  = shbb_.ilidlist[jtmp];
	send_innersize_hydro_[1] = je-js-joverlap_+1;
    // recv [js+joverlap:je] from other
	jtmp = jblock - Ngrids;
	if (jtmp < 0) jtmp += nrbx2;
    recv_inner_gid_[1]  = shbb_.igidlist[jtmp];
	recv_inner_rank_[1] = shbb_.irnklist[jtmp];
	recv_inner_lid_[1]  = shbb_.ilidlist[jtmp];
	recv_innersize_hydro_[1] = je-js-joverlap_+1;
	shbox_inner_hydro_flag_[1]=BNDRY_WAITING;
	shbox_inner_hydro_flag_[0]=BNDRY_COMPLETED; // default as completed for no-overlap case
	// if there is overlap to next blocks
	if (joverlap_ != 0) {
	  // send [je-(joverlap-1):je] of the meshblock to other
	  jtmp = jblock + (Ngrids + 1);
	  if (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
      send_inner_gid_[0]  = shbb_.igidlist[jtmp];
	  send_inner_rank_[0] = shbb_.irnklist[jtmp];
	  send_inner_lid_[0]  = shbb_.ilidlist[jtmp];
	  send_innersize_hydro_[0] = joverlap_;
      // recv [js:js+(joverlap-1)] from other
	  jtmp = jblock - (Ngrids + 1);
	  if (jtmp < 0) jtmp += nrbx2;
      recv_inner_gid_[0]  = shbb_.igidlist[jtmp];
	  recv_inner_rank_[0] = shbb_.irnklist[jtmp];
	  recv_inner_lid_[0]  = shbb_.ilidlist[jtmp];
	  recv_innersize_hydro_[0] = joverlap_;
	  shbox_inner_hydro_flag_[0]=BNDRY_WAITING; // switch on the boundary status if overlap
	}
  }

  if (shbb_.outer == true) { // if outer block
	for (int n=0; n<2; n++){
      send_outer_gid_[n]  = -1;
	  send_outer_rank_[n] = -1;
	  send_outer_lid_[n]  = -1;
      recv_outer_gid_[n]  = -1;
	  recv_outer_rank_[n] = -1;
	  recv_outer_lid_[n]  = -1;
	  send_outersize_hydro_[n] = 0;
	  recv_outersize_hydro_[n] = 0;
	}
	int jblock;
    for (int j=0; j<nrbx2; j++) {
      if (shbb_.ogidlist[j] == pmb->gid) jblock = j;// index of current meshblock on the shearingboundary block list
    }
	// recv [js:je-joverlap] of the meshblock from other
	int jtmp = jblock + Ngrids;
	if (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
    recv_outer_gid_[1]  = shbb_.ogidlist[jtmp];
	recv_outer_rank_[1] = shbb_.ornklist[jtmp];
	recv_outer_lid_[1]  = shbb_.olidlist[jtmp];
    recv_outersize_hydro_[1] = je-js-joverlap_+1;
    // send [js+joverlap:je] of the meshblock to other
	jtmp = jblock - Ngrids;
	if (jtmp < 0) jtmp += nrbx2;
    send_outer_gid_[1]  = shbb_.ogidlist[jtmp];
	send_outer_rank_[1] = shbb_.ornklist[jtmp];
	send_outer_lid_[1]  = shbb_.olidlist[jtmp];
    send_outersize_hydro_[1] = je-js-joverlap_+1;
	shbox_outer_hydro_flag_[1]=BNDRY_WAITING;
	shbox_outer_hydro_flag_[0]=BNDRY_COMPLETED; // default as completed for no-overlap case
	// if there is overlap to next blocks
	if (joverlap_ != 0) {
      // recv [je-(joverlap-1):je] from other
	  jtmp = jblock + (Ngrids + 1);
	  if (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
      recv_outer_gid_[0]  = shbb_.ogidlist[jtmp];
	  recv_outer_rank_[0] = shbb_.ornklist[jtmp];
	  recv_outer_lid_[0]  = shbb_.olidlist[jtmp];
      recv_outersize_hydro_[0] = joverlap_;
	  // send [js:js+(joverlap-1)] of the meshblock to other
	  jtmp = jblock - (Ngrids + 1);
	  if (jtmp < 0) jtmp += nrbx2;
      send_outer_gid_[0]  = shbb_.ogidlist[jtmp];
	  send_outer_rank_[0] = shbb_.ornklist[jtmp];
	  send_outer_lid_[0]  = shbb_.olidlist[jtmp];
      send_outersize_hydro_[0] = joverlap_;
	  shbox_outer_hydro_flag_[0]=BNDRY_WAITING; // switch on the boundary status if overlap
	}
  }

  return;
}
//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ReceiveHydroShearingboxBoundaryBuffers(int n, int k, int jinner, int jouter, int i, Real eps, static AthenaArray<Real> &U, AthenaArray<Real> &Flux)
//  \brief compute the flux along j indices for remapping
//  adopted from 2nd order RemapFlux of Athena4.0
void BoundaryValues::RemapFlux(const int n, const int k, const int jinner, const int jouter, const int i, const Real eps, const AthenaArray<Real> &U, AthenaArray<Real> &Flux)
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
      dUc = U(n,k,j+1,i) - U(n,k,j-1,i);
      dUl = U(n,k,j,  i) - U(n,k,j-1,i);
      dUr = U(n,k,j+1,i) - U(n,k,j,  i);

      dUm = 0.0;
      if (dUl*dUr > 0.0) {
        lim_slope = std::min(fabs(dUl),fabs(dUr));
        dUm = SIGN(dUc)*std::min(0.5*fabs(dUc),2.0*lim_slope);
      }

    if (eps > 0.0) { /* eps always > 0 for inner i boundary */
      Flux(j+1) = eps*(U(n,k,j,i) + 0.5*(1.0 - eps)*dUm);
    } else {         /* eps always < 0 for outer i boundary */
      Flux(j  ) = eps*(U(n,k,j,i) - 0.5*(1.0 + eps)*dUm);
    }
  }

  return;
}

