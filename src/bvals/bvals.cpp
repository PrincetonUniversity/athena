
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

// Primary header
#include "bvals.hpp"

// C++ headers
#include <iostream>   // endl
#include <iomanip>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cstring>    // memcpy

// Athena headers
#include "../athena.hpp"          // Real
#include "../athena_arrays.hpp"   // AthenaArray
#include "../mesh.hpp"            // MeshBlock
#include "../fluid/fluid.hpp"     // Fluid
#include "../parameter_input.hpp" // ParameterInput

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif


// arrays of start and end points, created in InitBoundaryBuffer
typedef struct NeighborIndexes
{
  int ox1, ox2, ox3, fi1, fi2;
} NeighborIndexes;

static int target_bufid_[56];
static NeighborIndexes ni_[56];

//======================================================================================
//! \file bvals.cpp
//  \brief implements functions that initialize/apply BCs on each dir
//======================================================================================

// BoundaryValues constructor - sets functions for the appropriate
// boundary conditions at each of the 6 dirs of a MeshBlock

BoundaryValues::BoundaryValues(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_mblock_ = pmb;

// Set BC functions for each of the 6 boundaries in turn -------------------------------
// Inner x1
  switch(pmb->block_bcs[inner_x1]){
    case 1:
      FluidBoundary_[inner_x1] = ReflectInnerX1;
      FieldBoundary_[inner_x1] = ReflectInnerX1;
      break;
    case 2:
      FluidBoundary_[inner_x1] = OutflowInnerX1;
      FieldBoundary_[inner_x1] = OutflowInnerX1;
      break;
    case -1: // block boundary
    case 3: // do nothing, useful for user-enrolled BCs
    case 4: // periodic boundary
      FluidBoundary_[inner_x1] = NULL;
      FieldBoundary_[inner_x1] = NULL;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
          << "Flag ix1_bc=" << pmb->block_bcs[inner_x1] << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
   }

// Outer x1
  switch(pmb->block_bcs[outer_x1]){
    case 1:
      FluidBoundary_[outer_x1] = ReflectOuterX1;
      FieldBoundary_[outer_x1] = ReflectOuterX1;
      break;
    case 2:
      FluidBoundary_[outer_x1] = OutflowOuterX1;
      FieldBoundary_[outer_x1] = OutflowOuterX1;
      break;
    case -1: // block boundary
    case 3: // do nothing, useful for user-enrolled BCs
    case 4: // periodic boundary
      FluidBoundary_[outer_x1] = NULL;
      FieldBoundary_[outer_x1] = NULL;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
          << "Flag ox1_bc=" << pmb->block_bcs[outer_x1] << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }

// Inner x2
  if (pmb->block_size.nx2 > 1) {
    switch(pmb->block_bcs[inner_x2]){
      case 1:
        FluidBoundary_[inner_x2] = ReflectInnerX2;
        FieldBoundary_[inner_x2] = ReflectInnerX2;
        break;
      case 2:
        FluidBoundary_[inner_x2] = OutflowInnerX2;
        FieldBoundary_[inner_x2] = OutflowInnerX2;
        break;
      case -1: // block boundary
      case 3: // do nothing, useful for user-enrolled BCs
      case 4: // periodic boundary
        FluidBoundary_[inner_x2] = NULL;
        FieldBoundary_[inner_x2] = NULL;
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
            << "Flag ix2_bc=" << pmb->block_bcs[inner_x2] << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
     }

// Outer x2
    switch(pmb->block_bcs[outer_x2]){
      case 1:
        FluidBoundary_[outer_x2] = ReflectOuterX2;
        FieldBoundary_[outer_x2] = ReflectOuterX2;
        break;
      case 2:
        FluidBoundary_[outer_x2] = OutflowOuterX2;
        FieldBoundary_[outer_x2] = OutflowOuterX2;
        break;
      case -1: // block boundary
      case 3: // do nothing, useful for user-enrolled BCs
      case 4: // periodic boundary
        FluidBoundary_[outer_x2] = NULL;
        FieldBoundary_[outer_x2] = NULL;
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
            << "Flag ox2_bc=" << pmb->block_bcs[outer_x2] << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
    }
  }

// Inner x3
  if (pmb->block_size.nx3 > 1) {
    switch(pmb->block_bcs[inner_x3]){
      case 1:
        FluidBoundary_[inner_x3] = ReflectInnerX3;
        FieldBoundary_[inner_x3] = ReflectInnerX3;
        break;
      case 2:
        FluidBoundary_[inner_x3] = OutflowInnerX3;
        FieldBoundary_[inner_x3] = OutflowInnerX3;
        break;
      case -1: // block boundary
      case 3: // do nothing, useful for user-enrolled BCs
      case 4: // periodic boundary
        FluidBoundary_[inner_x3] = NULL;
        FieldBoundary_[inner_x3] = NULL;
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
            << "Flag ix3_bc=" << pmb->block_bcs[inner_x3] << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
     }

// Outer x3
    switch(pmb->block_bcs[outer_x3]){
      case 1:
        FluidBoundary_[outer_x3] = ReflectOuterX3;
        FieldBoundary_[outer_x3] = ReflectOuterX3;
        break;
      case 2:
        FluidBoundary_[outer_x3] = OutflowOuterX3;
        FieldBoundary_[outer_x3] = OutflowOuterX3;
        break;
      case -1: // block boundary
      case 3: // do nothing, useful for user-enrolled BCs
      case 4: // periodic boundary
        FluidBoundary_[outer_x3] = NULL;
        FieldBoundary_[outer_x3] = NULL;
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
            << "Flag ox3_bc=" << pmb->block_bcs[outer_x3] << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
    }
  }

  // Clear flags and requests
  for(int l=0;l<NSTEP;l++) {
    for(int i=0;i<56;i++){
      fluid_flag_[l][i]=boundary_waiting;
      field_flag_[l][i]=boundary_waiting;
      fluid_send_[l][i]=NULL;
      fluid_recv_[l][i]=NULL;
      field_send_[l][i]=NULL;
      field_recv_[l][i]=NULL;
#ifdef MPI_PARALLEL
      req_fluid_send_[l][i]=MPI_REQUEST_NULL;
      req_fluid_recv_[l][i]=MPI_REQUEST_NULL;
#endif
    }
  }
  // Allocate Buffers
  int nface=1, nedge=1, cng=(NGHOST+1)/2+1, cng1=0, cng2=0, cng3=0;
  if(pmb->pmy_mesh->multilevel==true) {
    if(pmb->block_size.nx2>1) cng1=cng, cng2=cng;
    if(pmb->block_size.nx3>1) cng3=cng;
  }
  for(int l=0;l<NSTEP;l++) {
    for(int n=0;n<pmb->pmy_mesh->maxneighbor_;n++) {
      int size=((ni_[n].ox1==0)?pmb->block_size.nx1:NGHOST)
              *((ni_[n].ox2==0)?pmb->block_size.nx2:NGHOST)
              *((ni_[n].ox3==0)?pmb->block_size.nx3:NGHOST);
      if(pmb->pmy_mesh->multilevel==true) {
        int f2c=((ni_[n].ox1==0)?((pmb->block_size.nx1+1)/2):NGHOST)
               *((ni_[n].ox2==0)?((pmb->block_size.nx2+1)/2):NGHOST)
               *((ni_[n].ox3==0)?((pmb->block_size.nx3+1)/2):NGHOST);
        int c2f=((ni_[n].ox1==0)?((pmb->block_size.nx1+1)/2+cng1):cng)
               *((ni_[n].ox2==0)?((pmb->block_size.nx2+1)/2+cng2):cng)
               *((ni_[n].ox3==0)?((pmb->block_size.nx3+1)/2+cng3):cng);
        size=std::max(size,c2f);
        size=std::max(size,f2c);
      }
      size*=NFLUID;
      fluid_send_[l][n]=new Real[size];
      fluid_recv_[l][n]=new Real[size];
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for(int l=0;l<NSTEP;l++) {
      for(int n=0;n<pmb->pmy_mesh->maxneighbor_;n++) {
        int size1=((ni_[n].ox1==0)?(pmb->block_size.nx1+1):NGHOST)
                 *((ni_[n].ox2==0)?(pmb->block_size.nx2):NGHOST)
                 *((ni_[n].ox3==0)?(pmb->block_size.nx3):NGHOST);
        int size2=((ni_[n].ox1==0)?(pmb->block_size.nx1):NGHOST)
                 *((ni_[n].ox2==0)?(pmb->block_size.nx2+1):NGHOST)
                 *((ni_[n].ox3==0)?(pmb->block_size.nx3):NGHOST);
        int size3=((ni_[n].ox1==0)?(pmb->block_size.nx1):NGHOST)
                 *((ni_[n].ox2==0)?(pmb->block_size.nx2):NGHOST)
                 *((ni_[n].ox3==0)?(pmb->block_size.nx3+1):NGHOST);
        if(pmb->pmy_mesh->multilevel==true) {
          // *** need to implement
        }
        int size=size1+size2+size3;
        field_send_[l][n]=new Real[size];
        field_recv_[l][n]=new Real[size];
      }
    }
  }
}

// destructor

BoundaryValues::~BoundaryValues()
{
  MeshBlock *pmb=pmy_mblock_;
  int r=pmb->pmy_mesh->maxneighbor_;
  for(int l=0;l<NSTEP;l++) {
    for(int i=0;i<r;i++) {
      delete [] fluid_send_[l][i];
      delete [] fluid_recv_[l][i];
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for(int l=0;l<NSTEP;l++) {
      for(int i=0;i<r;i++) { 
        delete [] field_send_[l][i];
        delete [] field_recv_[l][i];
      }
    }
  }
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::Initialize(void)
//  \brief Initialize MPI requests
void BoundaryValues::Initialize(void)
{
#ifdef MPI_PARALLEL
  MeshBlock* pmb=pmy_mblock_;
  int mylevel=pmb->uid.GetLevel();
  int tag;
  int cng1, cng2, cng3;
  cng1=pmb->cnghost;
  cng2=(pmb->block_size.nx2>1)?cng1:0;
  cng3=(pmb->block_size.nx3>1)?cng1:0;
  
  for(int l=0;l<NSTEP;l++) {
    for(int n=0;n<pmb->nneighbor;n++) {
      NeighborBlock& nb=pmb->neighbor[n];
      if(nb.rank!=myrank) {
        int ssize, rsize;
        if(nb.level==mylevel) { // same
          ssize=rsize=((nb.ox1==0)?pmb->block_size.nx1:NGHOST)
                     *((nb.ox2==0)?pmb->block_size.nx2:NGHOST)
                     *((nb.ox3==0)?pmb->block_size.nx3:NGHOST);
        }
        else if(nb.level<mylevel) { // coarser
          ssize=((nb.ox1==0)?((pmb->block_size.nx1+1)/2):NGHOST)
               *((nb.ox2==0)?((pmb->block_size.nx2+1)/2):NGHOST)
               *((nb.ox3==0)?((pmb->block_size.nx3+1)/2):NGHOST);
          rsize=((nb.ox1==0)?((pmb->block_size.nx1+1)/2+cng1):cng1)
               *((nb.ox2==0)?((pmb->block_size.nx2+1)/2+cng2):cng2)
               *((nb.ox3==0)?((pmb->block_size.nx3+1)/2+cng3):cng3);
        }
        else { // finer
          ssize=((nb.ox1==0)?((pmb->block_size.nx1+1)/2+cng1):cng1)
               *((nb.ox2==0)?((pmb->block_size.nx2+1)/2+cng2):cng2)
               *((nb.ox3==0)?((pmb->block_size.nx3+1)/2+cng3):cng3);
          rsize=((nb.ox1==0)?((pmb->block_size.nx1+1)/2):NGHOST)
               *((nb.ox2==0)?((pmb->block_size.nx2+1)/2):NGHOST)
               *((nb.ox3==0)?((pmb->block_size.nx3+1)/2):NGHOST);
        }
        ssize*=NFLUID; rsize*=NFLUID;
        // specify the offsets in the view point of the target block: flip ox? signs
        tag=CreateMPITag(nb.lid, l, tag_fluid, -nb.ox1, -nb.ox2, -nb.ox3, 0, 0);
        MPI_Send_init(fluid_send_[l][nb.bufid],ssize,MPI_ATHENA_REAL,
                      nb.rank,tag,MPI_COMM_WORLD,&req_fluid_send_[l][nb.bufid]);
        tag=CreateMPITag(pmb->lid, l, tag_fluid, nb.ox1, nb.ox2, nb.ox3, 0, 0);
        MPI_Recv_init(fluid_recv_[l][nb.bufid],rsize,MPI_ATHENA_REAL,
                      nb.rank,tag,MPI_COMM_WORLD,&req_fluid_recv_[l][nb.bufid]);
      }
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for(int l=0;l<NSTEP;l++) {
      for(int n=0;n<pmb->nneighbor;n++) {
        NeighborBlock& nb=pmb->neighbor[n];
        if(nb.rank!=myrank) {
          int ssize, rsize;
          if(pmb->pmy_mesh->multilevel==false) { // uniform
            int size1=((nb.ox1==0)?(pmb->block_size.nx1+1):NGHOST)
                     *((nb.ox2==0)?(pmb->block_size.nx2):NGHOST)
                     *((nb.ox3==0)?(pmb->block_size.nx3):NGHOST);
            int size2=((nb.ox1==0)?(pmb->block_size.nx1):NGHOST)
                     *((nb.ox2==0)?(pmb->block_size.nx2+1):NGHOST)
                     *((nb.ox3==0)?(pmb->block_size.nx3):NGHOST);
            int size3=((nb.ox1==0)?(pmb->block_size.nx1):NGHOST)
                     *((nb.ox2==0)?(pmb->block_size.nx2):NGHOST)
                     *((nb.ox3==0)?(pmb->block_size.nx3+1):NGHOST);
            ssize=rsize=size1+size2+size3;
          }
          else {
            if(nb.level==mylevel) { // same
              // ****** need to implement ******
            }
            else if(nb.level<mylevel) { // coarser
              // ****** need to implement ******
            }
            else { // finer
              // ****** need to implement ******
            }
          }
          // specify the offsets in the view point of the target block: flip ox? signs
          tag=CreateMPITag(nb.lid, l, tag_field, -nb.ox1, -nb.ox2, -nb.ox3, 0, 0);
          MPI_Send_init(field_send_[l][nb.bufid],ssize,MPI_ATHENA_REAL,
                        nb.rank,tag,MPI_COMM_WORLD,&req_field_send_[l][nb.bufid]);
          tag=CreateMPITag(pmb->lid, l, tag_field, nb.ox1, nb.ox2, nb.ox3, 0, 0);
          MPI_Recv_init(field_recv_[l][nb.bufid],rsize,MPI_ATHENA_REAL,
                        nb.rank,tag,MPI_COMM_WORLD,&req_field_recv_[l][nb.bufid]);
        }
      }
    }
  }
#endif
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::EnrollFluidBoundaryFunction(enum direction dir,
//                                                       BValFluid_t my_bc)
//  \brief Enroll a user-defined boundary function for fluid

void BoundaryValues::EnrollFluidBoundaryFunction(enum direction dir, BValFluid_t my_bc)
{
  std::stringstream msg;
  if(dir<0 || dir>5) {
    msg << "### FATAL ERROR in EnrollFluidBoundaryCondition function" << std::endl
        << "dirName = " << dir << " not valid" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if(pmy_mblock_->block_bcs[dir]==-1) return;
  if(pmy_mblock_->block_bcs[dir]!=3) {
    msg << "### FATAL ERROR in EnrollFluidBoundaryCondition function" << std::endl
        << "A user-defined boundary condition flag (3) must be specified "
        << "in the input file to use a user-defined boundary function." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  FluidBoundary_[dir]=my_bc;
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::EnrollFieldBoundaryFunction(enum direction dir,
//                                                       BValField_t my_bc)
//  \brief Enroll a user-defined boundary function for magnetic fields

void BoundaryValues::EnrollFieldBoundaryFunction(enum direction dir,BValField_t my_bc)
{
  std::stringstream msg;
  if(dir<0 || dir>5) {
    msg << "### FATAL ERROR in EnrollFieldBoundaryCondition function" << std::endl
        << "dirName = " << dir << " is not valid" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if(pmy_mblock_->block_bcs[dir]==-1) return;
  if(pmy_mblock_->block_bcs[dir]!=3) {
    msg << "### FATAL ERROR in EnrollFieldBoundaryCondition function" << std::endl
        << "A user-defined boundary condition flag (3) must be specified "
        << "in the input file to use a user-defined boundary function." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  FieldBoundary_[dir]=my_bc;
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::CheckBoundary(void)
//  \brief checks if the boundary conditions are correctly enrolled
void BoundaryValues::CheckBoundary(void)
{
  int i, r=2;
  MeshBlock *pmb=pmy_mblock_;
  if(pmb->block_size.nx2 > 1) // 2D
    r=4;
  if(pmb->block_size.nx3 > 1) // 3D
    r=6;
  for(int i=0;i<r;i++) {
    if(pmb->block_bcs[i]==3) {
      if(FluidBoundary_[i]==NULL) {
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues::CheckBoundary" << std::endl
            << "A user-defined boundary is specified but the fluid boundary function "
            << "is not enrolled in direction " << i  << "." << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      if (MAGNETIC_FIELDS_ENABLED) {
        if(FieldBoundary_[i]==NULL) {
          std::stringstream msg;
          msg << "### FATAL ERROR in BoundaryValues::CheckBoundary" << std::endl
              << "A user-defined boundary is specified but the field boundary function "
              << "is not enrolled in direction " << i  << "." << std::endl;
          throw std::runtime_error(msg.str().c_str());
        }
      }
    }
  }
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::StartReceivingForInit(void)
//  \brief initiate MPI_Irecv for initialization
void BoundaryValues::StartReceivingForInit(void)
{
#ifdef MPI_PARALLEL
  MeshBlock *pmb=pmy_mblock_;
  for(int n=0;n<pmb->nneighbor;n++) {
    NeighborBlock& nb=pmb->neighbor[n];
    if(nb.rank!=myrank) { 
      MPI_Start(&req_fluid_recv_[0][nb.bufid]);
      if (MAGNETIC_FIELDS_ENABLED)
        MPI_Start(&req_field_recv_[0][nb.bufid]);
    }
  }
#endif
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::StartReceivingAll(void)
//  \brief initiate MPI_Irecv for all the sweeps
void BoundaryValues::StartReceivingAll(void)
{
#ifdef MPI_PARALLEL
  MeshBlock *pmb=pmy_mblock_;
  for(int l=0;l<NSTEP;l++) {
    for(int n=0;n<pmb->nneighbor;n++) {
      NeighborBlock& nb=pmb->neighbor[n];
      if(nb.rank!=myrank) { 
        MPI_Start(&req_fluid_recv_[l][nb.bufid]);
        if (MAGNETIC_FIELDS_ENABLED)
          MPI_Start(&req_field_recv_[l][nb.bufid]);
      }
    }
  }
#endif
  return;
}


//--------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadFluidBoundaryBufferSameLevel(AthenaArray<Real> &src,
//                                                 Real *buf, NeighborBlock& nb)
//  \brief Set fluid boundary buffers for sending to a block on the same level
int BoundaryValues::LoadFluidBoundaryBufferSameLevel(AthenaArray<Real> &src, Real *buf,
                                                     NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_mblock_;
  int si, sj, sk, ei, ej, ek;

  si=(nb.ox1>0)?(pmb->ie-NGHOST+1):pmb->is;
  ei=(nb.ox1<0)?(pmb->is+NGHOST-1):pmb->ie;
  sj=(nb.ox2>0)?(pmb->je-NGHOST+1):pmb->js;
  ej=(nb.ox2<0)?(pmb->js+NGHOST-1):pmb->je;
  sk=(nb.ox3>0)?(pmb->ke-NGHOST+1):pmb->ks;
  ek=(nb.ox3<0)?(pmb->ks+NGHOST-1):pmb->ke;

  // Set buffers
  int p=0;
  for (int n=0; n<(NFLUID); ++n) {
    for (int k=sk; k<=ek; ++k) {
      for (int j=sj; j<=ej; ++j) {
#pragma simd
        for (int i=si; i<=ei; ++i) {
          // buffer is always fully packed
          buf[p++]=src(n,k,j,i);
        }
      }
    }
  }
  return p;
}


//--------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadFluidBoundaryBufferToCoarser(AthenaArray<Real> &src,
//                                                 Real *buf, NeighborBlock& nb)
//  \brief Set fluid boundary buffers for sending to a block on the coarser level
int BoundaryValues::LoadFluidBoundaryBufferToCoarser(AthenaArray<Real> &src, Real *buf,
                                                     NeighborBlock& nb)
{
  int p=0;
// *** not implemented yet ***
  return p;
}


//--------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadFluidBoundaryBufferToFiner(AthenaArray<Real> &src,
//                                                 Real *buf, NeighborBlock& nb)
//  \brief Set fluid boundary buffers for sending to a block on the finer level
int BoundaryValues::LoadFluidBoundaryBufferToFiner(AthenaArray<Real> &src, Real *buf,
                                                   NeighborBlock& nb)
{
  int p=0;
// *** not implemented yet ***
  return p;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SendFluidBoundaryBuffers(AthenaArray<Real> &src, int step)
//  \brief Send boundary buffers
void BoundaryValues::SendFluidBoundaryBuffers(AthenaArray<Real> &src, int step)
{
  MeshBlock *pmb=pmy_mblock_;
  int mylevel=pmb->uid.GetLevel();

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb=pmb->neighbor[n];
    int ssize;
    if(nb.level==mylevel)
      ssize=LoadFluidBoundaryBufferSameLevel(src, fluid_send_[step][nb.bufid],nb);
    else if(nb.level<mylevel)
      ssize=LoadFluidBoundaryBufferToCoarser(src, fluid_send_[step][nb.bufid],nb);
    else
      ssize=LoadFluidBoundaryBufferToFiner(src, fluid_send_[step][nb.bufid], nb);
    if(nb.rank == myrank) { // on the same process
      MeshBlock *pbl=pmb->pmy_mesh->FindMeshBlock(nb.gid);
      // find target buffer
      int target = target_bufid_[nb.bufid];
      std::memcpy(pbl->pbval->fluid_recv_[step][target],
                  fluid_send_[step][nb.bufid], ssize*sizeof(Real));
      pbl->pbval->fluid_flag_[step][target]=boundary_arrived;
    }
#ifdef MPI_PARALLEL
    else // MPI
      MPI_Start(&req_fluid_send_[step][nb.bufid]);
#endif
  }

  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetFluidBoundarySameLevel(AthenaArray<Real> &dst,
//                                                     Real *buf, NeighborBlock& nb)
//  \brief Set fluid boundary received from a block on the same level
void BoundaryValues::SetFluidBoundarySameLevel(AthenaArray<Real> &dst, Real *buf,
                                               NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_mblock_;
  int si, sj, sk, ei, ej, ek;

  if(nb.ox1==0)     si=pmb->is,        ei=pmb->ie;
  else if(nb.ox1>0) si=pmb->ie+1,      ei=pmb->ie+NGHOST;
  else              si=pmb->is-NGHOST, ei=pmb->is-1;
  if(nb.ox2==0)     sj=pmb->js,        ej=pmb->je;
  else if(nb.ox2>0) sj=pmb->je+1,      ej=pmb->je+NGHOST;
  else              sj=pmb->js-NGHOST, ej=pmb->js-1;
  if(nb.ox3==0)     sk=pmb->ks,        ek=pmb->ke;
  else if(nb.ox3>0) sk=pmb->ke+1,      ek=pmb->ke+NGHOST;
  else              sk=pmb->ks-NGHOST, ek=pmb->ks-1;

  // Set buffers
  int p=0;
  for (int n=0; n<(NFLUID); ++n) {
    for (int k=sk; k<=ek; ++k) {
      for (int j=sj; j<=ej; ++j) {
#pragma simd
        for (int i=si; i<=ei; ++i) {
          // buffer is always fully packed
          dst(n,k,j,i) = buf[p++];
        }
      }
    }
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetFluidBoundaryFromCoarser(Real *buf, NeighborBlock& nb)
//  \brief Set fluid prolongation buffer received from a block on the same level
void BoundaryValues::SetFluidBoundaryFromCoarser(Real *buf, NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_mblock_;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetFluidBoundaryFromFiner(AthenaArray<Real> &dst,
//                                                     Real *buf, NeighborBlock& nb)
//  \brief Set fluid boundary received from a block on the same level
void BoundaryValues::SetFluidBoundaryFromFiner(AthenaArray<Real> &dst, Real *buf,
                                               NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_mblock_;
}


//--------------------------------------------------------------------------------------
//! \fn bool BoundaryValues::ReceiveFluidBoundaryBuffers(AthenaArray<Real> &dst, int step)
//  \brief receive the boundary data
bool BoundaryValues::ReceiveFluidBoundaryBuffers(AthenaArray<Real> &dst, int step)
{
  MeshBlock *pmb=pmy_mblock_;
  int mylevel=pmb->uid.GetLevel();
  int nc=0;

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb= pmb->neighbor[n];
    if(fluid_flag_[step][nb.bufid]==boundary_completed) { nc++; continue;}
    if(fluid_flag_[step][nb.bufid]==boundary_waiting) {
      if(nb.rank==myrank) // on the same process
        continue;
#ifdef MPI_PARALLEL
      else { // MPI boundary
        int test;
        MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&test,MPI_STATUS_IGNORE);
        MPI_Test(&req_fluid_recv_[step][nb.bufid],&test,MPI_STATUS_IGNORE);
        if(test==false) continue;
        fluid_flag_[step][nb.bufid] = boundary_arrived;
      }
#endif
    }
    if(nb.level==mylevel)
      SetFluidBoundarySameLevel(dst, fluid_recv_[step][nb.bufid], nb);
    else if(nb.level<mylevel)
      SetFluidBoundaryFromCoarser(fluid_recv_[step][nb.bufid], nb);
    else
      SetFluidBoundaryFromFiner(dst, fluid_recv_[step][nb.bufid], nb);

    fluid_flag_[step][nb.bufid] = boundary_completed; // completed
    nc++;
  }

  if(nc<pmb->nneighbor)
    return false;
  return true;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ReceiveFluidBoundaryBuffersWithWait(AthenaArray<Real> &dst,
//                                                               int step)
//  \brief receive the boundary data for initialization
void BoundaryValues::ReceiveFluidBoundaryBuffersWithWait(AthenaArray<Real> &dst, int step)
{
  MeshBlock *pmb=pmy_mblock_;
  int mylevel=pmb->uid.GetLevel();

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb= pmb->neighbor[n];
#ifdef MPI_PARALLEL
    if(nb.rank!=myrank)
      MPI_Wait(&req_fluid_recv_[0][nb.bufid],MPI_STATUS_IGNORE);
#endif
    if(nb.level==mylevel)
      SetFluidBoundarySameLevel(dst, fluid_recv_[0][nb.bufid], nb);
    else if(nb.level<mylevel)
      SetFluidBoundaryFromCoarser(fluid_recv_[0][nb.bufid], nb);
    else
      SetFluidBoundaryFromFiner(dst, fluid_recv_[0][nb.bufid], nb);
    fluid_flag_[0][nb.bufid] = boundary_completed; // completed
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadFieldBoundaryBufferSameLevel(InterfaceField &src,
//                                                 Real *buf, NeighborBlock& nb)
//  \brief Set field boundary buffers for sending to a block on the same level
int BoundaryValues::LoadFieldBoundaryBufferSameLevel(InterfaceField &src, Real *buf,
                                                     NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_mblock_;
  int si, sj, sk, ei, ej, ek;

  si=(nb.ox1>0)?(pmb->ie-NGHOST+1):pmb->is;
  ei=(nb.ox1<0)?(pmb->is+NGHOST-1):pmb->ie;
  sj=(nb.ox2>0)?(pmb->je-NGHOST+1):pmb->js;
  ej=(nb.ox2<0)?(pmb->js+NGHOST-1):pmb->je;
  sk=(nb.ox3>0)?(pmb->ke-NGHOST+1):pmb->ks;
  ek=(nb.ox3<0)?(pmb->ks+NGHOST-1):pmb->ke;

  // Set buffers
  int p=0;
  // bx1
  for (int n=0; n<(NFLUID); ++n) {
    for (int k=sk; k<=ek; ++k) {
      for (int j=sj; j<=ej; ++j) {
#pragma simd
        for (int i=si; i<=ei; ++i) {
          // buffer is always fully packed
          buf[p++]=src.x1f(k,j,i);
        }
      }
    }
  }
  // bx2
  for (int n=0; n<(NFLUID); ++n) {
    for (int k=sk; k<=ek; ++k) {
      for (int j=sj; j<=ej; ++j) {
#pragma simd
        for (int i=si; i<=ei; ++i) {
          // buffer is always fully packed
          buf[p++]=src.x2f(k,j,i);
        }
      }
    }
  }
  // bx3
  for (int n=0; n<(NFLUID); ++n) {
    for (int k=sk; k<=ek; ++k) {
      for (int j=sj; j<=ej; ++j) {
#pragma simd
        for (int i=si; i<=ei; ++i) {
          // buffer is always fully packed
          buf[p++]=src.x3f(k,j,i);
        }
      }
    }
  }

  return p;
}


//--------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadFieldBoundaryBufferToCoarser(InterfaceField &src,
//                                                 Real *buf, NeighborBlock& nb)
//  \brief Set field boundary buffers for sending to a block on the coarser level
int BoundaryValues::LoadFieldBoundaryBufferToCoarser(InterfaceField &src, Real *buf,
                                                     NeighborBlock& nb)
{
  return 0;
}

//--------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadFieldBoundaryBufferToFiner(InterfaceField &src, 
//                                                 Real *buf, NeighborBlock& nb)
//  \brief Set field boundary buffers for sending to a block on the finer level
int BoundaryValues::LoadFieldBoundaryBufferToFiner(InterfaceField &src, Real *buf,
                                                   NeighborBlock& nb)
{
  return 0;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SendFieldBoundaryBuffers(InterfaceField &src, int step)
//  \brief Send field boundary buffers
void BoundaryValues::SendFieldBoundaryBuffers(InterfaceField &src, int step)
{
  MeshBlock *pmb=pmy_mblock_;
  int mylevel=pmb->uid.GetLevel();

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb=pmb->neighbor[n];
    int ssize;
    if(nb.level==mylevel)
      ssize=LoadFieldBoundaryBufferSameLevel(src, field_send_[step][nb.bufid],nb);
    else if(nb.level<mylevel)
      ssize=LoadFieldBoundaryBufferToCoarser(src, field_send_[step][nb.bufid],nb);
    else
      ssize=LoadFieldBoundaryBufferToFiner(src, field_send_[step][nb.bufid], nb);
    if(nb.rank == myrank) { // on the same process
      MeshBlock *pbl=pmb->pmy_mesh->FindMeshBlock(nb.gid);
      // find target buffer
      int target = target_bufid_[nb.bufid];
      std::memcpy(pbl->pbval->field_recv_[step][target],
                  field_send_[step][nb.bufid], ssize*sizeof(Real));
      pbl->pbval->field_flag_[step][target]=boundary_arrived;
    }
#ifdef MPI_PARALLEL
    else // MPI
      MPI_Start(&req_field_send_[step][nb.bufid]);
#endif
  }

  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetFieldBoundarySameLevel(InterfaceField &dst,
//                                                     Real *buf, NeighborBlock& nb)
//  \brief Set field boundary received from a block on the same level
void BoundaryValues::SetFieldBoundarySameLevel(InterfaceField &dst, Real *buf,
                                               NeighborBlock& nb)
{
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetFieldBoundaryFromCoarser(Real *buf, NeighborBlock& nb)
//  \brief Set field prolongation buffer received from a block on the same level
void BoundaryValues::SetFieldBoundaryFromCoarser(Real *buf, NeighborBlock& nb)
{
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetFielBoundaryFromFiner(InterfaceField &dst,
//                                                     Real *buf, NeighborBlock& nb)
//  \brief Set field boundary received from a block on the same level
void BoundaryValues::SetFieldBoundaryFromFiner(InterfaceField &dst, Real *buf,
                                               NeighborBlock& nb)
{
  return;
}


//--------------------------------------------------------------------------------------
//! \fn bool BoundaryValues::ReceiveFieldBoundaryBuffers(InterfaceField &dst, int step)
//  \brief load boundary buffer for x1 direction into the array
bool BoundaryValues::ReceiveFieldBoundaryBuffers(InterfaceField &dst, int step)
{
  // implement me
  return true;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ReceiveFieldBoundaryBuffersWithWait(InterfaceField &dst,
//                                                               int step)
//  \brief load boundary buffer for x1 direction into the array
void BoundaryValues::ReceiveFieldBoundaryBuffersWithWait(InterfaceField &dst, int step)
{
  // implement me
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ClearBoundaryForInit(void)
//  \brief clean up the boundary flags for initialization
void BoundaryValues::ClearBoundaryForInit(void)
{
  MeshBlock *pmb=pmy_mblock_;
  for(int n=0;n<pmb->nneighbor;n++) {
    NeighborBlock& nb=pmb->neighbor[n];
    fluid_flag_[0][nb.bufid] = boundary_waiting;
    if (MAGNETIC_FIELDS_ENABLED)
      field_flag_[0][nb.bufid] = boundary_waiting;
#ifdef MPI_PARALLEL
    if(nb.rank!=myrank) {
      MPI_Wait(&req_fluid_send_[0][nb.bufid],MPI_STATUS_IGNORE); // Wait for Isend
      if (MAGNETIC_FIELDS_ENABLED)
        MPI_Wait(&req_field_send_[0][nb.bufid],MPI_STATUS_IGNORE); // Wait for Isend
    }
#endif
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ClearBoundaryAll(void)
//  \brief clean up the boundary flags after each loop
void BoundaryValues::ClearBoundaryAll(void)
{
  MeshBlock *pmb=pmy_mblock_;
  for(int l=0;l<NSTEP;l++) {
    for(int n=0;n<pmb->nneighbor;n++) {
      NeighborBlock& nb=pmb->neighbor[n];
      fluid_flag_[l][nb.bufid] = boundary_waiting;
      if (MAGNETIC_FIELDS_ENABLED)
        field_flag_[l][nb.bufid] = boundary_waiting;
#ifdef MPI_PARALLEL
      if(nb.rank!=myrank) {
        MPI_Wait(&req_fluid_send_[l][nb.bufid],MPI_STATUS_IGNORE); // Wait for Isend
        if (MAGNETIC_FIELDS_ENABLED)
          MPI_Wait(&req_field_send_[l][nb.bufid],MPI_STATUS_IGNORE); // Wait for Isend
      }
#endif
    }
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::FluidPhysicalBoundaries(AthenaArray<Real> &dst)
//  \brief Apply physical boundary conditions for fluid
void BoundaryValues::FluidPhysicalBoundaries(AthenaArray<Real> &dst)
{
  MeshBlock *pmb=pmy_mblock_;
  int bis=pmb->is-NGHOST, bie=pmb->ie+NGHOST;
  int bjs=pmb->js, bje=pmb->je, bks=pmb->ks, bke=pmb->ke;

  if(FluidBoundary_[inner_x2]==NULL && pmb->block_size.nx2>1) bjs=pmb->js-NGHOST;
  if(FluidBoundary_[outer_x2]==NULL && pmb->block_size.nx2>1) bje=pmb->je+NGHOST;
  if(FluidBoundary_[inner_x3]==NULL && pmb->block_size.nx3>1) bks=pmb->ks-NGHOST;
  if(FluidBoundary_[outer_x3]==NULL && pmb->block_size.nx3>1) bke=pmb->ke+NGHOST;

  if(FluidBoundary_[inner_x1]!=NULL)
    FluidBoundary_[inner_x1](pmb, dst, pmb->is, pmb->ie, bjs, bje, bks, bke);
  if(FluidBoundary_[outer_x1]!=NULL)
    FluidBoundary_[outer_x1](pmb, dst, pmb->is, pmb->ie, bjs, bje, bks, bke);
  if(pmb->block_size.nx2>1) { // 2D or 3D
    if(FluidBoundary_[inner_x2]!=NULL)
      FluidBoundary_[inner_x2](pmb, dst, bis, bie, pmb->js, pmb->je, bks, bke);
    if(FluidBoundary_[outer_x2]!=NULL)
      FluidBoundary_[outer_x2](pmb, dst, bis, bie, pmb->js, pmb->je, bks, bke);
  }
  if(pmb->block_size.nx3>1) { // 3D
    bjs=pmb->js-NGHOST;
    bje=pmb->je+NGHOST;
    if(FluidBoundary_[inner_x3]!=NULL)
      FluidBoundary_[inner_x3](pmb, dst, bis, bie, bjs, bje, pmb->ks, pmb->ke);
    if(FluidBoundary_[outer_x3]!=NULL)
      FluidBoundary_[outer_x3](pmb, dst, bis, bie, bjs, bje, pmb->ks, pmb->ke);
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::FieldPhysicalBoundaries(AthenaArray<Real> &dst)
//  \brief Apply physical boundary conditions for field
void BoundaryValues::FieldPhysicalBoundaries(InterfaceField &dst)
{
  MeshBlock *pmb=pmy_mblock_;
  int bis=pmb->is-NGHOST, bie=pmb->ie+NGHOST;
  int bjs=pmb->js, bje=pmb->je, bks=pmb->ks, bke=pmb->ke;

  if(FieldBoundary_[inner_x2]==NULL && pmb->block_size.nx2>1) bjs=pmb->js-NGHOST;
  if(FieldBoundary_[outer_x2]==NULL && pmb->block_size.nx2>1) bje=pmb->je+NGHOST;
  if(FieldBoundary_[inner_x3]==NULL && pmb->block_size.nx3>1) bks=pmb->ks-NGHOST;
  if(FieldBoundary_[outer_x3]==NULL && pmb->block_size.nx3>1) bke=pmb->ke+NGHOST;

  if(FieldBoundary_[inner_x1]!=NULL)
    FieldBoundary_[inner_x1](pmb, dst, pmb->is, pmb->ie, bjs, bje, bks, bke);
  if(FieldBoundary_[outer_x1]!=NULL)
    FieldBoundary_[outer_x1](pmb, dst, pmb->is, pmb->ie, bjs, bje, bks, bke);
  if(pmb->block_size.nx2>1) { // 2D or 3D
    if(FieldBoundary_[inner_x2]!=NULL)
      FieldBoundary_[inner_x2](pmb, dst, bis, bie, pmb->js, pmb->je, bks, bke);
    if(FieldBoundary_[outer_x2]!=NULL)
      FieldBoundary_[outer_x2](pmb, dst, bis, bie, pmb->js, pmb->je, bks, bke);
  }
  if(pmb->block_size.nx3>1) { // 3D
    bjs=pmb->js-NGHOST;
    bje=pmb->je+NGHOST;
    if(FieldBoundary_[inner_x3]!=NULL)
      FieldBoundary_[inner_x3](pmb, dst, bis, bie, bjs, bje, pmb->ks, pmb->ke);
    if(FieldBoundary_[outer_x3]!=NULL)
      FieldBoundary_[outer_x3](pmb, dst, bis, bie, bjs, bje, pmb->ks, pmb->ke);
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn int CreateBufferID(int ox1, int ox2, int ox3, int fi1, int fi2)
//  \brief calculate a buffer identifier
int CreateBufferID(int ox1, int ox2, int ox3, int fi1, int fi2)
{
  return ((((ox1+1)<<6)|((ox2+1)<<4)|(ox3+1))<<2) | (fi1<<1) | fi2;
}


//--------------------------------------------------------------------------------------
//! \fn int CreateMPITag(int lid, int flag, int phys, int ox1, int ox2, int ox3,
//                       int fi1, int fi2)
//  \brief calculate an MPI tag
int CreateMPITag(int lid, int flag, int phys, int ox1, int ox2, int ox3,
                 int fi1=0, int fi2=0)
{
  int dir=CreateBufferID(ox1, ox2, ox3, fi1, fi2);
// tag = local id of destination (17) + flag (2) + physics (4) + dir (6+2)
  return (lid<<14) | (flag<<12) | (phys<<8) | dir;
}


//--------------------------------------------------------------------------------------
//! \fn int BufferID(int dim, bool multilevel, bool face_only)
//  \brief calculate neighbor indexes and target buffer IDs
int BufferID(int dim, bool multilevel, bool face_only)
{
  int bufid[56];
  int nf1=1, nf2=1;
  if(multilevel==true) {
    if(dim>=2) nf1=2;
    if(dim>=3) nf2=2;
  }
  int b=0;
  // x1 face
  for(int n=-1; n<=1; n+=2) {
    for(int f2=0;f2<nf2;f2++) {
      for(int f1=0;f1<nf1;f1++) {
        ni_[b].ox1=n; ni_[b].ox2=0; ni_[b].ox3=0;
        ni_[b].fi1=f1; ni_[b].fi2=f2;
        b++;
      }
    }
  }
  // x2 face
  if(dim>=2) {
    for(int n=-1; n<=1; n+=2) {
      for(int f2=0;f2<nf2;f2++) {
        for(int f1=0;f1<nf1;f1++) {
          ni_[b].ox1=0; ni_[b].ox2=n; ni_[b].ox3=0;
          ni_[b].fi1=f1; ni_[b].fi2=f2;
          b++;
        }
      }
    }
  }
  if(dim==3) {
    // x3 face
    for(int n=-1; n<=1; n+=2) {
      for(int f2=0;f2<nf2;f2++) {
        for(int f1=0;f1<nf1;f1++) {
          ni_[b].ox1=0; ni_[b].ox2=0; ni_[b].ox3=n;
          ni_[b].fi1=f1; ni_[b].fi2=f2;
          b++;
        }
      }
    }
  }
  // edges
  // x1x2
  if(dim>=2) {
    for(int m=-1; m<=1; m+=2) {
      for(int n=-1; n<=1; n+=2) {
        for(int f1=0;f1<nf1;f1++) {
          ni_[b].ox1=n; ni_[b].ox2=m; ni_[b].ox3=0;
          ni_[b].fi1=f1; ni_[b].fi2=0;
          b++;
        }
      }
    }
  }
  if(dim==3) {
    // x1x3
    for(int m=-1; m<=1; m+=2) {
      for(int n=-1; n<=1; n+=2) {
        for(int f1=0;f1<nf1;f1++) {
          ni_[b].ox1=n; ni_[b].ox2=0; ni_[b].ox3=m;
          ni_[b].fi1=f1; ni_[b].fi2=0;
          b++;
        }
      }
    }
    // x2x3
    for(int m=-1; m<=1; m+=2) {
      for(int n=-1; n<=1; n+=2) {
        for(int f1=0;f1<nf1;f1++) {
          ni_[b].ox1=0; ni_[b].ox2=n; ni_[b].ox3=m;
          ni_[b].fi1=f1; ni_[b].fi2=0;
          b++;
        }
      }
    }
    // corners
    for(int l=-1; l<=1; l+=2) {
      for(int m=-1; m<=1; m+=2) {
        for(int n=-1; n<=1; n+=2) {
          ni_[b].ox1=n; ni_[b].ox2=m; ni_[b].ox3=l;
          ni_[b].fi1=0; ni_[b].fi2=0;
          b++;
        }
      }
    }
  }

  for(int n=0;n<b;n++)
    bufid[n]=CreateBufferID(ni_[n].ox1, ni_[n].ox2, ni_[n].ox3, ni_[n].fi1, ni_[n].fi2);

  // search buffer IDs
  int t=0;
  for(int n=0;n<b;n++) {
    int tid=CreateBufferID(-ni_[n].ox1, -ni_[n].ox2, -ni_[n].ox3, ni_[n].fi1, ni_[n].fi2);
    for(int i=0;i<b;i++) {
      if(tid==bufid[i]) {
        target_bufid_[n]=i;
        break;
      }
    }
  }
  return b;
}

