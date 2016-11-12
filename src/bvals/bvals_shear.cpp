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
  int nx2=pmb->block_size.nx2-1;

  si=pmb->is-NGHOST; ei=pmb->is-1;
  sk=pmb->ks;        ek=pmb->ke;
  if (pmesh->mesh_size.nx3>1)  ek += NGHOST, sk -= NGHOST;
  // --- inner boundary:
  // load [je-joverlap-(NGHOST-1):je] for send_inner[0]; nb=0
  // load [je-joverlap:je]            for send_inner[0]; nb=0  if joverlap == nx2-1
  // load [js:je-joverlap+NGHOST]     for send_inner[1]; nb=1
  // load [js:je]                     for send_inner[1]; nb=1  if joverlap == 0
  // load [je-(NGHOST-1):je]          for send_inner[2]; nb=2
  // load [je:je]                     for send_inner[2]; nb=2  if joverlap == nx2-1
  // load [js:js+NGHOST-1]            for send_inner[3]; nb=3
  //
  // --- outer boundary:
  // load [js:js+joverlap+NGHOST-1]   for send_outer[0]; nb=4
  // load [js:js+joverlap]            for send_outer[0]; nb=4  if joverlap == nx2-1
  // load [js+joverlap-NGHOST:je]     for send_outer[1]; nb=5
  // load [js:je]                     for send_outer[1]; nb=5  if joverlap == 0
  // load [js:js+NGHOST-1]            for send_outer[2]; nb=6
  // load [js:js]                     for send_outer[2]; nb=6  if joverlap == nx2-1
  // load [je-NGHOST+1:je]            for send_outer[3]; nb=7
  switch(nb) {
    case 0:
      sj=pmb->je-joverlap_-(NGHOST-1); ej=pmb->je;
      if(joverlap_==nx2) sj = pmb->je-joverlap_;
      break;
    case 1:
      sj=pmb->js; ej=pmb->je-joverlap_+NGHOST;
      if(joverlap_==0) ej = pmb->je;
      break;
    case 2:
      sj=pmb->je-(NGHOST-1); ej=pmb->je;
      if(joverlap_==nx2) sj = pmb->je;
      break;
    case 3:
      sj=pmb->js; ej=pmb->js+(NGHOST-1);
      break;
    case 4:
      sj=pmb->js; ej=pmb->js+joverlap_+NGHOST-1;
      if(joverlap_==nx2) ej = pmb->js+joverlap_;
      break;
    case 5:
      sj=pmb->js+joverlap_-NGHOST; ej=pmb->je;
      if(joverlap_==0) sj = pmb->js;
      break;
    case 6:
      sj=pmb->js; ej=pmb->js+(NGHOST-1);
      if(joverlap_==nx2) ej=pmb->js;
      break;
    case 7:
      sj=pmb->je-NGHOST+1; ej=pmb->je;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues:LoadHydroShearing " << std::endl
          << "nb = " << nb << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }

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
          shboxvar_inner_hydro_(IM2,k,j,i) = src(IM2,k,j,i) + qomL*src(IDN,k,j,i);//add shear to conservative
          if (NON_BAROTROPIC_EOS) {
            src(IEN,k,j,i) += (0.5/src(IDN,k,j,i))*(SQR(shboxvar_inner_hydro_(IM2,k,j,i)) -
                                                    SQR(src(IM2,k,j,i)));
          } // update energy
          src(IM2,k,j,i) = shboxvar_inner_hydro_(IM2,k,j,i);//update IM2
        }
    }}
  }

  if (shbb_.outer == true) {
    int  ib = ie+1;
    int ii;
    // step 2. -- add shear to the outer periodic boundary values
    for(int k=kl; k<=ku; k++) {
      for(int j=js-NGHOST; j<=je+NGHOST; j++) {
        for(int i=0; i<NGHOST; i++) {
          ii = ib+i;
          shboxvar_outer_hydro_(IM2,k,j,i) = src(IM2,k,j,ii) - qomL*src(IDN,k,j,ii);//add shear to conservative
          if (NON_BAROTROPIC_EOS) {
            src(IEN,k,j,ii) += (0.5/src(IDN,k,j,ii))*(SQR(shboxvar_inner_hydro_(IM2,k,j,i)) -
                                                      SQR(src(IM2,k,j,ii)));
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
          RemapFlux(n,k,js,je+2,i,eps_,shboxvar_inner_hydro_,flx_inner_hydro_);
          for(int j=js; j<=je+1; j++) {
            shboxvar_inner_hydro_(n,k,j,i) -= flx_inner_hydro_(j+1)-flx_inner_hydro_(j);
          }
        }
    }}

  // step 3. -- load sendbuf; memcpy to recvbuf if on same rank, post MPI_Isend otherwise
    for(int n=0; n<4; n++) {
      if(send_inner_rank_[n] != -1) {
        LoadHydroShearing(shboxvar_inner_hydro_, send_innerbuf_hydro_[n], n);
        if (send_inner_rank_[n] == Globals::my_rank) {// on the same process
          MeshBlock *pbl=pmb->pmy_mesh->FindMeshBlock(send_inner_gid_[n]);
          std::memcpy(pbl->pbval->recv_innerbuf_hydro_[n],send_innerbuf_hydro_[n],
                      send_innersize_hydro_[n]*ssize*sizeof(Real));
          pbl->pbval->shbox_inner_hydro_flag_[n]=BNDRY_ARRIVED;
        } else { // MPI
#ifdef MPI_PARALLEL
          int tag=CreateBvalsMPITag(send_inner_lid_[n], TAG_SHBOX_HYDRO, n); //bufid = n
          MPI_Isend(send_innerbuf_hydro_[n],send_innersize_hydro_[n]*ssize,
                    MPI_ATHENA_REAL,send_inner_rank_[n],tag,MPI_COMM_WORLD,
                    &rq_innersend_hydro_[n]);
#endif
        }
      }}
  } // inner boundaries

  if (shbb_.outer == true) {
    int  ib = ie+1;
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
          RemapFlux(n,k,js-1,je+1,i,-eps_,shboxvar_outer_hydro_,flx_outer_hydro_);
          for(int j=js-1; j<=je; j++) {
            shboxvar_outer_hydro_(n,k,j,i) -= flx_outer_hydro_(j+1)-flx_outer_hydro_(j);
          }
        }
    }}

  // step 3. -- load sendbuf; memcpy to recvbuf if on same rank, post MPI_Isend otherwise
    int offset = 4;
    for(int n=0; n<4; n++) {
      if(send_outer_rank_[n] != -1) {
        LoadHydroShearing(shboxvar_outer_hydro_, send_outerbuf_hydro_[n], n+offset);
        if (send_outer_rank_[n] == Globals::my_rank) {// on the same process
            MeshBlock *pbl=pmb->pmy_mesh->FindMeshBlock(send_outer_gid_[n]);
            std::memcpy(pbl->pbval->recv_outerbuf_hydro_[n],send_outerbuf_hydro_[n],
                        send_outersize_hydro_[n]*ssize*sizeof(Real));
            pbl->pbval->shbox_outer_hydro_flag_[n]=BNDRY_ARRIVED;
        } else { // MPI
#ifdef MPI_PARALLEL
            int tag=CreateBvalsMPITag(send_outer_lid_[n], TAG_SHBOX_HYDRO, n+offset); //bufid for outer(inner): 2(0) and 3(1)
            MPI_Isend(send_outerbuf_hydro_[n],send_outersize_hydro_[n]*ssize,
                      MPI_ATHENA_REAL,send_outer_rank_[n],tag,MPI_COMM_WORLD,
                      &rq_outersend_hydro_[n]);
#endif
        }
    }}
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
  int nx2=pmb->block_size.nx2-1;

  sk = pmb->ks; ek = pmb->ke;
  if (pmesh->mesh_size.nx3>1) ek += NGHOST, sk -= NGHOST;
  // --- inner boundary:
  // set [js-NGHOST:js+(joverlap-1)] with recv_inner[0]; nb=0
  // set [js-1:js+(joverlap-1)]      with recv_inner[0]; nb=0 if joverlap==nx2-1
  // set [js+joverlap:je+NGHOST]     with recv_inner[1]; nb=1
  // set [js+joverlap:je]            with recv_inner[1]; nb=1 if joverlap==0
  // set [js-NGHOST:js-1]            with recv_inner[2]; nb=2
  // set [js-NGHOST:js-NGHOST]       with recv_inner[2]; nb=2 if joverlap==nx2-1
  // set [je+1:je+NGHOST]            with recv_inner[3]; nb=3
  // --- outer boundary:
  // set [je-joverlap+1:je+NGHOST]   with recv_outer[0]; nb=4
  // set [je-joverlap+1:je+1]        with recv_outer[0]; nb=4 if joverlap==nx2-1
  // set [js-NGHOST:je-joverlap]     with recv_outer[1]; nb=5
  // set [js:je-joverlap]            with recv_outer[1]; nb=5 if joverlap==0
  // set [je+1:je+NGHOST]            with recv_outer[2]; nb=6
  // set [je+NGHOST:je+NGHOST]       with recv_outer[2]; nb=6 if joverlap==nx2-1
  // set [js-NGHOST:js-1]            with recv_outer[3]; nb=7
  switch(nb) {
    case 0:
      si=pmb->is-NGHOST; ei=pmb->is-1; sj=pmb->js-NGHOST; ej=pmb->js+(joverlap_-1);
      if(joverlap_==nx2) sj=pmb->js-1;
      break;
    case 1:
      si=pmb->is-NGHOST; ei=pmb->is-1; sj=pmb->js+joverlap_; ej=pmb->je+NGHOST;
      if(joverlap_==0)   ej=pmb->je;
      break;
    case 2:
      si=pmb->is-NGHOST; ei=pmb->is-1; sj=pmb->js-NGHOST; ej=pmb->js-1;
      if(joverlap_==nx2) ej=pmb->js-NGHOST;
      break;
    case 3:
      si=pmb->is-NGHOST; ei=pmb->is-1; sj=pmb->je+1; ej=pmb->je+NGHOST;
      break;
    case 4:
      si=pmb->ie+1; ei=pmb->ie+NGHOST; sj=pmb->je-(joverlap_-1); ej=pmb->je+NGHOST;
      if(joverlap_==nx2) ej=pmb->je+1;
      break;
    case 5:
      si=pmb->ie+1; ei=pmb->ie+NGHOST; sj=pmb->js-NGHOST; ej=pmb->je-joverlap_;
      if(joverlap_==0)   sj=pmb->js;
      break;
    case 6:
      si=pmb->ie+1; ei=pmb->ie+NGHOST; sj=pmb->je+1; ej=pmb->je+NGHOST;
      if(joverlap_==nx2)   sj=pmb->je+NGHOST;
      break;
    case 7:
      si=pmb->ie+1; ei=pmb->ie+NGHOST; sj=pmb->js-NGHOST; ej=pmb->js-1;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues:SetHydroShearing " << std::endl
          << "nb = " << nb << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }

  // set [sj:ej] of current meshblock
  int p=0;
  BufferUtility::Unpack4DData(buf, dst, 0, NHYDRO-1, si, ei, sj, ej, sk, ek, p);
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ReceiveHydroShearingboxBoundaryBuffersWithWait(AthenaArray<Real> &dst, bool conserved_value)
//  \brief receive shearingbox boundary data for hydro variables for initialization
//         It is not used for now as our problem generator might start with t=0 periodic config.
void BoundaryValues::ReceiveHydroShearingboxBoundaryBuffersWithWait(AthenaArray<Real> &dst, bool conserved_value)
{
  MeshBlock *pmb=pmy_mblock_;

  if(shbb_.inner == true) { // check inner boundaries
#ifdef MPI_PARALLEL
    if (recv_inner_rank_[1]!=Globals::my_rank)
      MPI_Wait(&rq_innerrecv_hydro_[1],MPI_STATUS_IGNORE);
#endif
    SetHydroShearingboxBoundarySameLevel(dst,recv_innerbuf_hydro_[1],1);
    shbox_inner_hydro_flag_[1] = BNDRY_COMPLETED;
  } // inner boundary

  if(shbb_.outer == true) { // check inner boundaries
#ifdef MPI_PARALLEL
    if (recv_outer_rank_[1]!=Globals::my_rank)
      MPI_Wait(&rq_outerrecv_hydro_[1],MPI_STATUS_IGNORE);
#endif
    SetHydroShearingboxBoundarySameLevel(dst,recv_outerbuf_hydro_[1],5);
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

  if(shbb_.inner == true) { // check inner boundaries
    for(int n=0; n<4; n++) {
      if(shbox_inner_hydro_flag_[n]==BNDRY_COMPLETED) continue;
      if(shbox_inner_hydro_flag_[n]==BNDRY_WAITING) {
        if (recv_inner_rank_[n]==Globals::my_rank) {// on the same process
          flagi=false;
          continue;
        }
        else { // MPI boundary
#ifdef MPI_PARALLEL
          int test;
          MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&test,MPI_STATUS_IGNORE);
          MPI_Test(&rq_innerrecv_hydro_[n],&test,MPI_STATUS_IGNORE);
          if(test==false) {
            flagi=false;
            continue;
          }
          shbox_inner_hydro_flag_[n] = BNDRY_ARRIVED;
#endif
        }
      }
      // set dst if boundary arrived
      SetHydroShearingboxBoundarySameLevel(dst,recv_innerbuf_hydro_[n],n);
      shbox_inner_hydro_flag_[n] = BNDRY_COMPLETED; // completed
    } // loop over recv[0] to recv[3]
  } // inner boundary

  if(shbb_.outer == true) { // check outer boundaries
    int offset = 4;
    for(int n=0; n<4; n++) {
      if(shbox_outer_hydro_flag_[n]==BNDRY_COMPLETED) continue;
      if(shbox_outer_hydro_flag_[n]==BNDRY_WAITING) {
        if (recv_outer_rank_[n]==Globals::my_rank) {// on the same process
          flago=false;
          continue;
        }
        else { // MPI boundary
#ifdef MPI_PARALLEL
          int test;
          MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&test,MPI_STATUS_IGNORE);
          MPI_Test(&rq_outerrecv_hydro_[n],&test,MPI_STATUS_IGNORE);
          if(test==false) {
            flago=false;
            continue;
          }
          shbox_outer_hydro_flag_[n] = BNDRY_ARRIVED;
#endif
        }
      }
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

void BoundaryValues::FindShearBlock(const int step)
//void BoundaryValues::FindShearBlock(void)
{
  MeshBlock *pmb=pmy_mblock_;
  Coordinates *pco=pmb->pcoord;
  Mesh *pmesh=pmb->pmy_mesh;

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int ku, ii,jj;

  // pmesh->nrbx2 # of meshblocks in azimuth on root level
  // assuming same-level refinement across all horizontal blocks at given height
  int nrbx2 = pmesh->nrbx2<<(pmb->loc.level-pmesh->root_level);
  int nx2   = pmb->block_size.nx2; // # of cells per meshblock
  int nx3   = pmb->block_size.nx3; // # of cells per meshblock
  int ncells2 = pmb->block_size.nx2 + 2*NGHOST;
  int ncells3 = pmb->block_size.nx3;
  if (pmesh->mesh_size.nx3>1) ncells3 += 2*NGHOST;

  Real qomL = qshear_*Omega_0_*x1size_;
  Real wght;
  //(0.0,0.5,1.0) for VL2, (0.0,1.0,0.5) for RK2
  if (step == 0) wght = 0.0;
  else wght = wghts_[step-1];
  Real yshear = qomL*(pmesh->time+wght*pmesh->dt);
  //Real yshear = qomL*(pmesh->time+pmesh->dt);
  Real deltay = fmod(yshear,x2size_);
  int joffset = (int)(deltay/pco->dx2v(js)); // this assumes uniform grid in azimuth
  int Ngrids  = (int)(joffset/nx2);
  joverlap_   = joffset - Ngrids*nx2;
  eps_ = (fmod(deltay,pco->dx2v(js)))/pco->dx2v(js);

  if (shbb_.inner == true) { // if inner block
    for (int n=0; n<5; n++){
      send_inner_gid_[n]  = -1;
      send_inner_rank_[n] = -1;
      send_inner_lid_[n]  = -1;
      recv_inner_gid_[n]  = -1;
      recv_inner_rank_[n] = -1;
      recv_inner_lid_[n]  = -1;
      if (n<4) {
        send_innersize_hydro_[n] = 0;
        recv_innersize_hydro_[n] = 0;
        shbox_inner_hydro_flag_[n] = BNDRY_COMPLETED;
      }
      if (MAGNETIC_FIELDS_ENABLED) {
        if (n<4) {
          send_innersize_field_[n] = 0;
          recv_innersize_field_[n] = 0;
          shbox_inner_field_flag_[n]=BNDRY_COMPLETED;
        }
        send_innersize_emf_[n] = 0;
        recv_innersize_emf_[n] = 0;
        shbox_inner_emf_flag_[n]=BNDRY_COMPLETED;
      }
    }
    int jblock;
    for (int j=0; j<nrbx2; j++) {
      // index of current meshblock on the shearingboundary block list
      if (shbb_.igidlist[j] == pmb->gid)  jblock = j;
    }
    // send [js:je-joverlap] of the meshblock to other
    // attach [je-joverlap+1:je-joverlap+(NGHOST)] to its right end.
    int jtmp = jblock + Ngrids;
    if (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
    send_inner_gid_[1]  = shbb_.igidlist[jtmp];
    send_inner_rank_[1] = shbb_.irnklist[jtmp];
    send_inner_lid_[1]  = shbb_.ilidlist[jtmp];
    send_innersize_hydro_[1] = je-js-joverlap_+1+NGHOST;
    // recv [js+joverlap:je] from other
    // attach [je+1:je+NGHOST] to its right end.
    jtmp = jblock - Ngrids;
    if (jtmp < 0) jtmp += nrbx2;
    recv_inner_gid_[1]  = shbb_.igidlist[jtmp];
    recv_inner_rank_[1] = shbb_.irnklist[jtmp];
    recv_inner_lid_[1]  = shbb_.ilidlist[jtmp];
    recv_innersize_hydro_[1] = send_innersize_hydro_[1]; //je-js-joverlap_+1+NGHOST;
    shbox_inner_hydro_flag_[1] = BNDRY_WAITING;
    if (MAGNETIC_FIELDS_ENABLED) {
      send_innersize_field_[1] = send_innersize_hydro_[1]*NGHOST*(NFIELD*ncells3+1)
                                +NGHOST*ncells3;
      recv_innersize_field_[1] = send_innersize_field_[1];
      shbox_inner_field_flag_[1] = BNDRY_WAITING;
      send_innersize_emf_[1] = 2*send_innersize_hydro_[1]*nx3+
                              send_innersize_hydro_[1]+nx3; //(hydro_size+1)*nx3+hydro_size*(nx3+1)
      recv_innersize_emf_[1] = send_innersize_emf_[1];
      shbox_inner_emf_flag_[1] = BNDRY_WAITING;
    }


    // if there is overlap to next blocks
    if (joverlap_ != 0) {
      // send [je-joverlap-(NGHOST-1):je]
      // [je-(joverlap-1):je] + [je-joverlap-(NGHOST-1):je-joverlap]
      jtmp = jblock + (Ngrids + 1);
      if (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
      send_inner_gid_[0]  = shbb_.igidlist[jtmp];
      send_inner_rank_[0] = shbb_.irnklist[jtmp];
      send_inner_lid_[0]  = shbb_.ilidlist[jtmp];
      send_innersize_hydro_[0] = joverlap_+NGHOST;
      // recv [js-NGHOST:js+(joverlap-1)]
      // [js:js+(joverlap-1)] + [js-NGHOST:js-1] to its left end
      jtmp = jblock - (Ngrids + 1);
      if (jtmp < 0) jtmp += nrbx2;
      recv_inner_gid_[0]  = shbb_.igidlist[jtmp];
      recv_inner_rank_[0] = shbb_.irnklist[jtmp];
      recv_inner_lid_[0]  = shbb_.ilidlist[jtmp];
      recv_innersize_hydro_[0] = joverlap_+NGHOST;
      shbox_inner_hydro_flag_[0] = BNDRY_WAITING; // switch on the boundary status if overlap
      if (MAGNETIC_FIELDS_ENABLED) {
        send_innersize_field_[0] = send_innersize_hydro_[0]*NGHOST*(NFIELD*ncells3+1);
        recv_innersize_field_[0] = send_innersize_field_[0];
        shbox_inner_field_flag_[0] = BNDRY_WAITING;
        send_innersize_emf_[0] = send_innersize_hydro_[0]*(2*nx3+1);
        recv_innersize_emf_[0] = send_innersize_emf_[0];
        shbox_inner_emf_flag_[0] = BNDRY_WAITING;
      }
      // deal the left corner cells with send[2]
      if (joverlap_ == (nx2-1)) { // only pass one attached corner cell
        // send[0] sends [je-joverlap:je] now
        // recv[0] recvs [js-1:js+(joverlap-1)] now
        send_innersize_hydro_[0] -= 1;
        recv_innersize_hydro_[0] -= 1;
        if (MAGNETIC_FIELDS_ENABLED) {
          send_innersize_field_[0] = send_innersize_hydro_[0]*NGHOST*(NFIELD*ncells3+1);
          recv_innersize_field_[0] = send_innersize_field_[0];
          shbox_inner_field_flag_[0] = BNDRY_WAITING;
          send_innersize_emf_[0] = send_innersize_hydro_[0]*(2*nx3+1);
          recv_innersize_emf_[0] = send_innersize_emf_[0];
          shbox_inner_emf_flag_[0] = BNDRY_WAITING;
        }
        // send [je:je] to Right
        jtmp = jblock + (Ngrids + 2);
        while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
        send_inner_gid_[2]  = shbb_.igidlist[jtmp];
        send_inner_rank_[2] = shbb_.irnklist[jtmp];
        send_inner_lid_[2]  = shbb_.ilidlist[jtmp];
        send_innersize_hydro_[2] = 1;
        // recv [js-NGHOST:js-NGHOST] from Left
        jtmp = jblock - (Ngrids+2);
        while (jtmp < 0) jtmp += nrbx2;
        recv_inner_gid_[2]  = shbb_.igidlist[jtmp];
        recv_inner_rank_[2] = shbb_.irnklist[jtmp];
        recv_inner_lid_[2]  = shbb_.ilidlist[jtmp];
        recv_innersize_hydro_[2] = 1;
        shbox_inner_hydro_flag_[2] = BNDRY_WAITING;
        if (MAGNETIC_FIELDS_ENABLED) {
          send_innersize_field_[2] = send_innersize_hydro_[2]*NGHOST*(NFIELD*ncells3+1);
          recv_innersize_field_[2] = send_innersize_field_[2];
          shbox_inner_field_flag_[2] = BNDRY_WAITING;
          send_innersize_emf_[2] = send_innersize_hydro_[2]*(2*nx3+1);
          recv_innersize_emf_[2] = send_innersize_emf_[2];
          shbox_inner_emf_flag_[2] = BNDRY_WAITING;
        }
      }
      // deal the extra of EMF
      if (joverlap_ == 1) { // has to send[4] for EMF (similar to send[3] while joverlap=0)
        jtmp = jblock + (Ngrids - 1);
        while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
        while (jtmp < 0) jtmp += nrbx2;
        send_inner_gid_[4]  = shbb_.igidlist[jtmp];
        send_inner_rank_[4] = shbb_.irnklist[jtmp];
        send_inner_lid_[4]  = shbb_.ilidlist[jtmp];
        jtmp = jblock - (Ngrids-1);
        while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
        while (jtmp < 0) jtmp += nrbx2;
        recv_inner_gid_[4]  = shbb_.igidlist[jtmp];
        recv_inner_rank_[4] = shbb_.irnklist[jtmp];
        recv_inner_lid_[4]  = shbb_.ilidlist[jtmp];
        if (MAGNETIC_FIELDS_ENABLED) {
          send_innersize_emf_[4] = (NGHOST-1)*(2*nx3+1);
          recv_innersize_emf_[4] = send_innersize_emf_[4];
          shbox_inner_emf_flag_[4] = BNDRY_WAITING;
          send_innersize_emf_[1] = 2*nx3*nx2+nx2+nx3;
          recv_innersize_emf_[1] = send_innersize_emf_[1];
          shbox_inner_emf_flag_[1] = BNDRY_WAITING;
        }
      }
    } else { //joverlap_ == 0 need both send[2] and send[3]
      // do not pass attached L/R corner cells; use send[2][3] instead.
      // send[1] sends [js:je]
      // recv[1] recvs [js:je]
      send_innersize_hydro_[1] -= NGHOST;
      recv_innersize_hydro_[1] -= NGHOST;
      if (MAGNETIC_FIELDS_ENABLED) {
        send_innersize_field_[1] = send_innersize_hydro_[1]*NGHOST*(NFIELD*ncells3+1)
                              +NGHOST*ncells3;
        recv_innersize_field_[1] = send_innersize_field_[1];
        shbox_inner_field_flag_[1] = BNDRY_WAITING;
        send_innersize_emf_[1] = 2*send_innersize_hydro_[1]*nx3 + send_innersize_hydro_[1]
                            + nx3; //(hydro_size+1)*nx3+hydro_size*(nx3+1)
        recv_innersize_emf_[1] = send_innersize_emf_[1];
        shbox_inner_emf_flag_[1] = BNDRY_WAITING;
      }
      // send [je-(NGHOST-1):je] to Right
      jtmp = jblock + (Ngrids+1);
      while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
      send_inner_gid_[2]  = shbb_.igidlist[jtmp];
      send_inner_rank_[2] = shbb_.irnklist[jtmp];
      send_inner_lid_[2]  = shbb_.ilidlist[jtmp];
      send_innersize_hydro_[2] = NGHOST;
        // recv [js-NGHOST:js-NGHOST+1] from Left
      jtmp = jblock - (Ngrids+1);
      while (jtmp < 0) jtmp += nrbx2;
      recv_inner_gid_[2]  = shbb_.igidlist[jtmp];
      recv_inner_rank_[2] = shbb_.irnklist[jtmp];
      recv_inner_lid_[2]  = shbb_.ilidlist[jtmp];
      recv_innersize_hydro_[2] = NGHOST;
      shbox_inner_hydro_flag_[2] = BNDRY_WAITING;
      if (MAGNETIC_FIELDS_ENABLED) {
        send_innersize_field_[2] = send_innersize_hydro_[2]*NGHOST*(NFIELD*ncells3+1);
        recv_innersize_field_[2] = send_innersize_field_[2];
        shbox_inner_field_flag_[2] = BNDRY_WAITING;
        send_innersize_emf_[2] = send_innersize_hydro_[2]*(2*nx3+1);
        recv_innersize_emf_[2] = send_innersize_emf_[2];
        shbox_inner_emf_flag_[2] = BNDRY_WAITING;
      }

      // send [js:js+(NGHOST-1)] to Left
      jtmp = jblock + (Ngrids - 1);
      while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
      while (jtmp < 0) jtmp += nrbx2;
      send_inner_gid_[3]  = shbb_.igidlist[jtmp];
      send_inner_rank_[3] = shbb_.irnklist[jtmp];
      send_inner_lid_[3]  = shbb_.ilidlist[jtmp];
      send_innersize_hydro_[3] = NGHOST;
      // recv [je+1:je+(NGHOST-1)] from Right
      jtmp = jblock - (Ngrids-1);
      while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
      while (jtmp < 0) jtmp += nrbx2;
      recv_inner_gid_[3]  = shbb_.igidlist[jtmp];
      recv_inner_rank_[3] = shbb_.irnklist[jtmp];
      recv_inner_lid_[3]  = shbb_.ilidlist[jtmp];
      recv_innersize_hydro_[3] = NGHOST;
      shbox_inner_hydro_flag_[3] = BNDRY_WAITING;
      if (MAGNETIC_FIELDS_ENABLED) {
        send_innersize_field_[3] = send_innersize_hydro_[3]*NGHOST*(NFIELD*ncells3+1);
        recv_innersize_field_[3] = send_innersize_field_[3];
        shbox_inner_field_flag_[3] = BNDRY_WAITING;
        send_innersize_emf_[3] = send_innersize_hydro_[3]*(2*nx3+1);
        recv_innersize_emf_[3] = send_innersize_emf_[3];
        shbox_inner_emf_flag_[3] = BNDRY_WAITING;
      }
    }
  } // inner bc


  if (shbb_.outer == true) { // if outer block
    for (int n=0; n<5; n++){
      send_outer_gid_[n]  = -1;
      send_outer_rank_[n] = -1;
      send_outer_lid_[n]  = -1;
      recv_outer_gid_[n]  = -1;
      recv_outer_rank_[n] = -1;
      recv_outer_lid_[n]  = -1;
      if (n<4) {
        send_outersize_hydro_[n] = 0;
        recv_outersize_hydro_[n] = 0;
        shbox_outer_hydro_flag_[n] = BNDRY_COMPLETED;
      }
      if (MAGNETIC_FIELDS_ENABLED) {
        if (n<4) {
          send_outersize_field_[n] = 0;
          recv_outersize_field_[n] = 0;
          shbox_outer_field_flag_[n]=BNDRY_COMPLETED;
        }
        send_outersize_emf_[n] = 0;
        recv_outersize_emf_[n] = 0;
        shbox_outer_emf_flag_[n]=BNDRY_COMPLETED;
      }
    }
    int jblock;
    for (int j=0; j<nrbx2; j++) {
      if (shbb_.ogidlist[j] == pmb->gid) jblock = j;// index of current meshblock on the shearingboundary block list
    }
    // recv [js-NGHOST:je-joverlap] of the meshblock from other
    int jtmp = jblock + Ngrids;
    if (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
    recv_outer_gid_[1]  = shbb_.ogidlist[jtmp];
    recv_outer_rank_[1] = shbb_.ornklist[jtmp];
    recv_outer_lid_[1]  = shbb_.olidlist[jtmp];
    recv_outersize_hydro_[1] = je-js-joverlap_+1+NGHOST;
      // send [js+joverlap-NGHOST:je] of the meshblock to other
    jtmp = jblock - Ngrids;
    if (jtmp < 0) jtmp += nrbx2;
    send_outer_gid_[1]  = shbb_.ogidlist[jtmp];
    send_outer_rank_[1] = shbb_.ornklist[jtmp];
    send_outer_lid_[1]  = shbb_.olidlist[jtmp];
    send_outersize_hydro_[1] = je-js-joverlap_+1+NGHOST;
    shbox_outer_hydro_flag_[1]=BNDRY_WAITING;
    if (MAGNETIC_FIELDS_ENABLED) {
      send_outersize_field_[1] = send_outersize_hydro_[1]*NGHOST*(NFIELD*ncells3+1);
      recv_outersize_field_[1] = send_outersize_field_[1];
      shbox_outer_field_flag_[1] = BNDRY_WAITING;
      send_outersize_emf_[1] = send_outersize_hydro_[1]*(2*nx3+1);
      recv_outersize_emf_[1] = send_outersize_emf_[1];
      shbox_outer_emf_flag_[1] = BNDRY_WAITING;
    }

    // if there is overlap to next blocks
    if (joverlap_ != 0) {
      // recv [je-(joverlap-1):je+NGHOST] from other
      jtmp = jblock + (Ngrids + 1);
      if (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
      recv_outer_gid_[0]  = shbb_.ogidlist[jtmp];
      recv_outer_rank_[0] = shbb_.ornklist[jtmp];
      recv_outer_lid_[0]  = shbb_.olidlist[jtmp];
      recv_outersize_hydro_[0] = joverlap_+NGHOST;
      // send [js:js+(joverlap-1)+NGHOST] of the meshblock to other
      jtmp = jblock - (Ngrids + 1);
      if (jtmp < 0) jtmp += nrbx2;
      send_outer_gid_[0]  = shbb_.ogidlist[jtmp];
      send_outer_rank_[0] = shbb_.ornklist[jtmp];
      send_outer_lid_[0]  = shbb_.olidlist[jtmp];
      send_outersize_hydro_[0] = joverlap_+NGHOST;
      shbox_outer_hydro_flag_[0]=BNDRY_WAITING; // switch on the boundary status if overlap
      if (MAGNETIC_FIELDS_ENABLED) {
        send_outersize_field_[0] = send_outersize_hydro_[0]*NGHOST*(NFIELD*ncells3+1)
                                +NGHOST*ncells3;
        recv_outersize_field_[0] = send_outersize_field_[0];
        shbox_outer_field_flag_[0] = BNDRY_WAITING;
        send_outersize_emf_[0] = send_outersize_hydro_[0]*(2*nx3+1)+nx3;
        recv_outersize_emf_[0] = send_outersize_emf_[0];
        shbox_outer_emf_flag_[0] = BNDRY_WAITING;
      }
      // deal the left corner cells with send[2]
      if (joverlap_ == (nx2-1)) {
        // send[0] sends [js:js+(joverlap-1)+NGHOST-1] now
        // recv[0] recvs [je-(joverlap-1):je+NGHOST-1] now
        send_outersize_hydro_[0] -= 1;
        recv_outersize_hydro_[0] -= 1;
        if (MAGNETIC_FIELDS_ENABLED) {
          send_outersize_field_[0] = send_outersize_hydro_[0]*NGHOST*(NFIELD*ncells3+1)
                                    +NGHOST*ncells3;
          recv_outersize_field_[0] = send_outersize_field_[0];
          shbox_outer_field_flag_[0] = BNDRY_WAITING;
          send_outersize_emf_[0] = send_outersize_hydro_[0]*(2*nx3+1)+nx3;
          recv_outersize_emf_[0] = send_outersize_emf_[0];
          shbox_outer_emf_flag_[0] = BNDRY_WAITING;
        }
        // recv [je+NGHOST:je+NGHOST] from Left
        jtmp = jblock + (Ngrids + 2);
        while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
        recv_outer_gid_[2]  = shbb_.ogidlist[jtmp];
        recv_outer_rank_[2] = shbb_.ornklist[jtmp];
        recv_outer_lid_[2]  = shbb_.olidlist[jtmp];
        recv_outersize_hydro_[2] = 1;
        // send [js:js] to Right
        jtmp = jblock - (Ngrids+2);
        while (jtmp < 0) jtmp += nrbx2;
        send_outer_gid_[2]  = shbb_.ogidlist[jtmp];
        send_outer_rank_[2] = shbb_.ornklist[jtmp];
        send_outer_lid_[2]  = shbb_.olidlist[jtmp];
        send_outersize_hydro_[2] = 1;
        shbox_outer_hydro_flag_[2] = BNDRY_WAITING;
        if (MAGNETIC_FIELDS_ENABLED) {
          send_outersize_field_[2] = send_outersize_hydro_[2]*NGHOST*(NFIELD*ncells3+1);
          recv_outersize_field_[2] = send_outersize_field_[2];
          shbox_outer_field_flag_[2] = BNDRY_WAITING;
          send_outersize_emf_[2] = send_outersize_hydro_[2]*(2*nx3+1);
          recv_outersize_emf_[2] = send_outersize_emf_[2];
          shbox_outer_emf_flag_[2] = BNDRY_WAITING;
        }
      }
      // deal the extra of EMF
      if (joverlap_ == 1) { // has to send[4] for EMF (similar to send[3] while joverlap=0)
        jtmp = jblock + (Ngrids - 1);
        while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
        while (jtmp < 0) jtmp += nrbx2;
        recv_outer_gid_[4]  = shbb_.ogidlist[jtmp];
        recv_outer_rank_[4] = shbb_.ornklist[jtmp];
        recv_outer_lid_[4]  = shbb_.olidlist[jtmp];
        jtmp = jblock - (Ngrids-1);
        while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
        while (jtmp < 0) jtmp += nrbx2;
        send_outer_gid_[4]  = shbb_.ogidlist[jtmp];
        send_outer_rank_[4] = shbb_.ornklist[jtmp];
        send_outer_lid_[4]  = shbb_.olidlist[jtmp];
        if (MAGNETIC_FIELDS_ENABLED) {
          send_outersize_emf_[4] = (NGHOST-1)*(2*nx3+1);
          recv_outersize_emf_[4] = send_outersize_emf_[4];
          send_outersize_emf_[1] = (2*nx3+1)*nx2;
          recv_outersize_emf_[1] = send_outersize_emf_[1];
          shbox_outer_emf_flag_[1] = BNDRY_WAITING;
          shbox_outer_emf_flag_[4] = BNDRY_WAITING;
        }
      }
    } else { //joverlap_ == 0 need both send[2] and send[3]
      // send[1] sends [js:je]
      // recv[1] recvs [js:je]
      send_outersize_hydro_[1] -= NGHOST;
      recv_outersize_hydro_[1] -= NGHOST;
      if (MAGNETIC_FIELDS_ENABLED) {
        send_outersize_field_[1] = send_outersize_hydro_[1]*NGHOST*(NFIELD*ncells3+1)
                              +NGHOST*ncells3;
        recv_outersize_field_[1] = send_outersize_field_[1];
        shbox_outer_field_flag_[1] = BNDRY_WAITING;
        send_outersize_emf_[1] = send_outersize_hydro_[1]*(2*nx3+1)+nx3;
        recv_outersize_emf_[1] = send_outersize_emf_[1];
        shbox_outer_emf_flag_[1] = BNDRY_WAITING;
      }
      // recv [je+1:je+NGHOST] from Left
      jtmp = jblock + (Ngrids+1);
      while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
      recv_outer_gid_[2]  = shbb_.ogidlist[jtmp];
      recv_outer_rank_[2] = shbb_.ornklist[jtmp];
      recv_outer_lid_[2]  = shbb_.olidlist[jtmp];
      recv_outersize_hydro_[2] = NGHOST;
      // send [js:js+NGHOST-1] to Right
      jtmp = jblock - (Ngrids+1);
      while (jtmp < 0) jtmp += nrbx2;
      send_outer_gid_[2]  = shbb_.ogidlist[jtmp];
      send_outer_rank_[2] = shbb_.ornklist[jtmp];
      send_outer_lid_[2]  = shbb_.olidlist[jtmp];
      send_outersize_hydro_[2] = NGHOST;
      shbox_outer_hydro_flag_[2] = BNDRY_WAITING;
      if (MAGNETIC_FIELDS_ENABLED) {
        send_outersize_field_[2] = send_outersize_hydro_[2]*NGHOST*(NFIELD*ncells3+1);
        recv_outersize_field_[2] = send_outersize_field_[2];
        shbox_outer_field_flag_[2] = BNDRY_WAITING;
        send_outersize_emf_[2] = send_outersize_hydro_[2]*(2*nx3+1);
        recv_outersize_emf_[2] = send_outersize_emf_[2];
        shbox_outer_emf_flag_[2] = BNDRY_WAITING;
      }

      // recv [js-NGHOST:js-1] from Left
      jtmp = jblock + (Ngrids - 1);
      while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
      while (jtmp < 0) jtmp += nrbx2;
      recv_outer_gid_[3]  = shbb_.ogidlist[jtmp];
      recv_outer_rank_[3] = shbb_.ornklist[jtmp];
      recv_outer_lid_[3]  = shbb_.olidlist[jtmp];
      recv_outersize_hydro_[3] = NGHOST;
      // send [je-(NGHOST-1):je] to Right
      jtmp = jblock - (Ngrids-1);
      while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
      while (jtmp < 0) jtmp += nrbx2;
      send_outer_gid_[3]  = shbb_.ogidlist[jtmp];
      send_outer_rank_[3] = shbb_.ornklist[jtmp];
      send_outer_lid_[3]  = shbb_.olidlist[jtmp];
      send_outersize_hydro_[3] = NGHOST;
      shbox_outer_hydro_flag_[3] = BNDRY_WAITING;
      if (MAGNETIC_FIELDS_ENABLED) {
        send_outersize_field_[3] = send_outersize_hydro_[3]*NGHOST*(NFIELD*ncells3+1);
        recv_outersize_field_[3] = send_outersize_field_[3];
        shbox_outer_field_flag_[3] = BNDRY_WAITING;
        send_outersize_emf_[3] = send_outersize_hydro_[3]*(2*nx3+1);
        recv_outersize_emf_[3] = send_outersize_emf_[3];
        shbox_outer_emf_flag_[3] = BNDRY_WAITING;
      }
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

  //jinner,jouter are index range over which flux must be returned.  Set loop
  //limits depending on direction of upwind differences

  if (eps > 0.0) { // eps always > 0 for inner i boundary
    jl = jinner-1;
    ju = jouter-1;
  } else {         // eps always < 0 for outer i boundary
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

    if (eps > 0.0) { // eps always > 0 for inner i boundary
      Flux(j+1) = eps*(U(n,k,j,i) + 0.5*(1.0 - eps)*dUm);
    } else {         // eps always < 0 for outer i boundary
      Flux(j  ) = eps*(U(n,k,j,i) - 0.5*(1.0 + eps)*dUm);
    }
  }

  return;
}

