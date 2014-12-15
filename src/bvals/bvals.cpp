
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
int fluid_send_se_[6][6];
int fluid_recv_se_[6][6];
int field_send_se_[6][3][6];
int field_recv_se_[6][3][6];
int fluid_bufsize_[6];
int field_bufsize_[6];

//======================================================================================
//! \file bvals.cpp
//  \brief implements functions that initialize/apply BCs on each dir
//======================================================================================

// BoundaryValues constructor - sets functions for the appropriate
// boundary conditions at each of the 6 dirs of a MeshBlock

BoundaryValues::BoundaryValues(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_mblock_ = pmb;
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

// Set BC functions for each of the 6 boundaries in turn -------------------------------
// Inner x1

  switch(pmb->block_bcs.ix1_bc){
    case -1:
      FluidBoundary_[inner_x1] = NULL;
      FieldBoundary_[inner_x1] = NULL;
      break;
    case 1:
      FluidBoundary_[inner_x1] = ReflectInnerX1;
      FieldBoundary_[inner_x1] = ReflectInnerX1;
      break;
    case 2:
      FluidBoundary_[inner_x1] = OutflowInnerX1;
      FieldBoundary_[inner_x1] = OutflowInnerX1;
      break;
    case 3: // do nothing, useful for user-enrolled BCs
      FluidBoundary_[inner_x1] = NULL;
      FieldBoundary_[inner_x1] = NULL;
      break;
    case 4:
      FluidBoundary_[inner_x1] = PeriodicInnerX1;
      FieldBoundary_[inner_x1] = PeriodicInnerX1;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in FluidBCs constructor" << std::endl
          << "Flag ix1_bc=" << pmb->block_bcs.ix1_bc << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
   }

// Outer x1

  switch(pmb->block_bcs.ox1_bc){
    case -1:
      FluidBoundary_[outer_x1] = NULL;
      FieldBoundary_[outer_x1] = NULL;
      break;
    case 1:
      FluidBoundary_[outer_x1] = ReflectOuterX1;
      FieldBoundary_[outer_x1] = ReflectOuterX1;
      break;
    case 2:
      FluidBoundary_[outer_x1] = OutflowOuterX1;
      FieldBoundary_[outer_x1] = OutflowOuterX1;
      break;
    case 3: // do nothing, useful for user-enrolled BCs
      FluidBoundary_[outer_x1] = NULL;
      FieldBoundary_[outer_x1] = NULL;
      break;
    case 4:
      FluidBoundary_[outer_x1] = PeriodicOuterX1;
      FieldBoundary_[outer_x1] = PeriodicOuterX1;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in FluidBCs constructor" << std::endl
          << "Flag ox1_bc=" << pmb->block_bcs.ox1_bc << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }

// Inner x2

  if (pmb->block_size.nx2 > 1) {
    switch(pmb->block_bcs.ix2_bc){
      case -1:
        FluidBoundary_[inner_x2] = NULL;
        FieldBoundary_[inner_x2] = NULL;
        break;
      case 1:
        FluidBoundary_[inner_x2] = ReflectInnerX2;
        FieldBoundary_[inner_x2] = ReflectInnerX2;
        break;
      case 2:
        FluidBoundary_[inner_x2] = OutflowInnerX2;
        FieldBoundary_[inner_x2] = OutflowInnerX2;
        break;
      case 3: // do nothing, useful for user-enrolled BCs
        FluidBoundary_[inner_x2] = NULL;
        FieldBoundary_[inner_x2] = NULL;
        break;
      case 4:
        FluidBoundary_[inner_x2] = PeriodicInnerX2;
        FieldBoundary_[inner_x2] = PeriodicInnerX2;
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in FluidBCs constructor" << std::endl
            << "Flag ix2_bc=" << pmb->block_bcs.ix2_bc << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
     }

// Outer x2

    switch(pmb->block_bcs.ox2_bc){
      case -1:
        FluidBoundary_[outer_x2] = NULL;
        FieldBoundary_[outer_x2] = NULL;
        break;
      case 1:
        FluidBoundary_[outer_x2] = ReflectOuterX2;
        FieldBoundary_[outer_x2] = ReflectOuterX2;
        break;
      case 2:
        FluidBoundary_[outer_x2] = OutflowOuterX2;
        FieldBoundary_[outer_x2] = OutflowOuterX2;
        break;
      case 3: // do nothing, useful for user-enrolled BCs
        FluidBoundary_[outer_x2] = NULL;
        FieldBoundary_[outer_x2] = NULL;
        break;
      case 4:
        FluidBoundary_[outer_x2] = PeriodicOuterX2;
        FieldBoundary_[outer_x2] = PeriodicOuterX2;
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in FluidBCs constructor" << std::endl
            << "Flag ox2_bc=" << pmb->block_bcs.ox2_bc << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
    }
  }

// Inner x3

  if (pmb->block_size.nx3 > 1) {
    switch(pmb->block_bcs.ix3_bc){
      case -1:
        FluidBoundary_[inner_x3] = NULL;
        FieldBoundary_[inner_x3] = NULL;
        break;
      case 1:
        FluidBoundary_[inner_x3] = ReflectInnerX3;
        FieldBoundary_[inner_x3] = ReflectInnerX3;
        break;
      case 2:
        FluidBoundary_[inner_x3] = OutflowInnerX3;
        FieldBoundary_[inner_x3] = OutflowInnerX3;
        break;
      case 3: // do nothing, useful for user-enrolled BCs
        FluidBoundary_[inner_x3] = NULL;
        FieldBoundary_[inner_x3] = NULL;
        break;
      case 4:
        FluidBoundary_[inner_x3] = PeriodicInnerX3;
        FieldBoundary_[inner_x3] = PeriodicInnerX3;
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in FluidBCs constructor" << std::endl
            << "Flag ix3_bc=" << pmb->block_bcs.ix3_bc << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
     }

// Outer x3

    switch(pmb->block_bcs.ox3_bc){
      case -1:
        FluidBoundary_[outer_x3] = NULL;
        FieldBoundary_[outer_x3] = NULL;
        break;
      case 1:
        FluidBoundary_[outer_x3] = ReflectOuterX3;
        FieldBoundary_[outer_x3] = ReflectOuterX3;
        break;
      case 2:
        FluidBoundary_[outer_x3] = OutflowOuterX3;
        FieldBoundary_[outer_x3] = OutflowOuterX3;
        break;
      case 3: // do nothing, useful for user-enrolled BCs
        FluidBoundary_[outer_x3] = NULL;
        FieldBoundary_[outer_x3] = NULL;
        break;
      case 4:
        FluidBoundary_[outer_x3] = PeriodicOuterX3;
        FieldBoundary_[outer_x3] = PeriodicOuterX3;
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in FluidBCs constructor" << std::endl
            << "Flag ox3_bc=" << pmb->block_bcs.ox3_bc << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
    }
  }

  // Allocate Buffers
  int r=2;
  if(pmb->block_size.nx2 > 1) r=4;
  if(pmb->block_size.nx3 > 1) r=6;
  for(int i=0;i<r;i++)
  {
    fluid_send_[i]=new Real[fluid_bufsize_[i]];
    fluid_recv_[i]=new Real[fluid_bufsize_[i]];
  }

  if (MAGNETIC_FIELDS_ENABLED) {
    for(int i=0;i<r;i++)
    {
      field_send_[i]=new Real[field_bufsize_[i]];
      field_recv_[i]=new Real[field_bufsize_[i]];
    }
  }

  // initialize flags
  for(int i=0;i<6;i++) {
    fluid_flag_[i][0][0]=false;
    fluid_flag_[i][0][1]=false;
    fluid_flag_[i][1][0]=false;
    fluid_flag_[i][1][1]=false;
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for(int i=0;i<6;i++) {
      field_flag_[i][0][0]=false;
      field_flag_[i][0][1]=false;
      field_flag_[i][1][0]=false;
      field_flag_[i][1][1]=false;
    }
  }
}

// destructor

BoundaryValues::~BoundaryValues()
{
  int r=2;
  if(pmy_mblock_->block_size.nx2 > 1) r=4;
  if(pmy_mblock_->block_size.nx3 > 1) r=6;
  for(int i=0;i<r;i++) {
    delete [] fluid_send_[i];
    delete [] fluid_recv_[i];
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for(int i=0;i<r;i++) { 
      delete [] field_send_[i];
      delete [] field_recv_[i];
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn
//  \brief

void BoundaryValues::EnrollFluidBoundaryFunction(enum direction dir, BValFluid_t my_bc)
{
  if(dir<0 || dir>5)
  {
    std::stringstream msg;
    msg << "### FATAL ERROR in EnrollFluidBoundaryCondition function" << std::endl
        << "dirName = " << dir << " not valid" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if(pmy_mblock_->neighbor[dir][0][0].gid==-1)
    FluidBoundary_[dir]=my_bc;
  return;
}

//--------------------------------------------------------------------------------------
//! \fn
//  \brief

void BoundaryValues::EnrollFieldBoundaryFunction(enum direction dir,BValField_t my_bc)
{
  if(dir<0 || dir>5)
  {
    std::stringstream msg;
    msg << "### FATAL ERROR in EnrollFieldBoundaryCondition function" << std::endl
        << "dirName = " << dir << " not valid" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if(pmy_mblock_->neighbor[dir][0][0].gid==-1)
    FieldBoundary_[dir]=my_bc;
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::LoadAndSendFluidBoundaryBuffer
//                          (enum direction dir, AthenaArray<Real> &src)
//  \brief Set boundary buffer for x1 direction using boundary functions
//  note: some geometric boundaries (e.g. origin and pole) are not implemented yet
void BoundaryValues::LoadAndSendFluidBoundaryBuffer
                     (enum direction dir, AthenaArray<Real> &src)
{
  MeshBlock *pmb=pmy_mblock_;
  MeshBlock *pbl=pmb->pmy_mesh->pblock;
  int oside;
  std::stringstream msg;
  Real *sendbuf=fluid_send_[dir];
  int si, sj, sk, ei, ej, ek;

  si=fluid_send_se_[dir][0];
  ei=fluid_send_se_[dir][1];
  sj=fluid_send_se_[dir][2];
  ej=fluid_send_se_[dir][3];
  sk=fluid_send_se_[dir][4];
  ek=fluid_send_se_[dir][5];

  if(pmb->neighbor[dir][0][0].gid==-1)
    return; // do nothing for physical boundary

  if(dir%2==0)
    oside=dir+1;
  else
    oside=dir-1;

  // Set buffers
  int p=0;
  for (int n=0; n<(NFLUID); ++n) {
    for (int k=sk; k<=ek; ++k) {
      for (int j=sj; j<=ej; ++j) {
#pragma simd
        for (int i=si; i<=ei; ++i) {
          // buffer is always fully packed
          sendbuf[p++]=src(n,k,j,i);
        }
      }
    }
  }

  // Send the buffer; modify this for MPI and AMR
  if(pmb->neighbor[dir][0][0].rank == myrank) // myrank
  {
    while(pbl!=NULL)
    {
      if(pbl->gid==pmb->neighbor[dir][0][0].gid)
        break;
      pbl=pbl->next;
    }
    if(pbl==NULL)
    {
      msg << "### FATAL ERROR in SetFluidBoundary" << std::endl
          << "In boundary " << dir << " of block "<<  pmb->gid
          << ", neighbor block " << pmb->neighbor[dir][0][0].gid
          << " is not found on the same process." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    std::memcpy(pbl->pbval->fluid_recv_[oside], fluid_send_[dir],
                fluid_bufsize_[dir]*sizeof(Real));
    pbl->pbval->fluid_flag_[oside][0][0]=true; // the other side
  }
  else // MPI
  {
      msg << "### FATAL ERROR in SetFluidBoundary" << std::endl
          << "MPI is not implemented yet!!" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn bool BoundaryValues::ReceiveAndSetFluidBoundary(enum direction dir,
//                                                      AthenaArray<Real> &dst)
//  \brief load boundary buffer for x1 direction into the array
bool BoundaryValues::ReceiveAndSetFluidBoundary(enum direction dir, AthenaArray<Real> &dst)
{
  MeshBlock *pmb=pmy_mblock_;
  std::stringstream msg;
  Real *recvbuf=fluid_recv_[dir];
  int si, sj, sk, ei, ej, ek;

  if(pmb->neighbor[dir][0][0].gid==-1) // physical boundary
    FluidBoundary_[dir](pmb,dst);
  else // block boundary
  {
    if(fluid_flag_[dir][0][0] == false)
    {
      std::cout << "Not Ready: " << dir <<" from " << pmb->gid << std::endl;
      return false; // return if it is not ready yet
    }

    si=fluid_recv_se_[dir][0];
    ei=fluid_recv_se_[dir][1];
    sj=fluid_recv_se_[dir][2];
    ej=fluid_recv_se_[dir][3];
    sk=fluid_recv_se_[dir][4];
    ek=fluid_recv_se_[dir][5];

    int p=0;
    for (int n=0; n<(NFLUID); ++n) {
      for (int k=sk; k<=ek; ++k) {
        for (int j=sj; j<=ej; ++j) {
#pragma simd
          for (int i=si; i<=ei; ++i) {
            // buffer is always fully packed
            dst(n,k,j,i) = recvbuf[p++];
          }
        }
      }
    }
    fluid_flag_[dir][0][0] = false; // clear the flag
  }
  return true;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::LoadAndSendFieldBoundaryBuffer
//                           (enum direction dir, InterfaceField &src)
//  \brief Set boundary buffer for x1 direction using boundary functions
//  note: some geometric boundaries (e.g. origin and pole) are not implemented yet
void BoundaryValues::LoadAndSendFieldBoundaryBuffer(enum direction dir,
                                                    InterfaceField &src)
{
  MeshBlock *pmb=pmy_mblock_;
  MeshBlock *pbl=pmb->pmy_mesh->pblock;
  int oside;
  std::stringstream msg;
  Real *sendbuf=field_send_[dir];
  AthenaArray<Real>& x1src=src.x1f;
  AthenaArray<Real>& x2src=src.x2f;
  AthenaArray<Real>& x3src=src.x3f;
  int si, sj, sk, ei, ej, ek;

  if(pmb->neighbor[dir][0][0].gid==-1)
    return; // do nothing for physical boundary

  if(dir%2==0)
    oside=dir+1;
  else
    oside=dir-1;

  // Set buffers; x1f
  int p=0;
  si=field_send_se_[dir][x1face][0];
  ei=field_send_se_[dir][x1face][1];
  sj=field_send_se_[dir][x1face][2];
  ej=field_send_se_[dir][x1face][3];
  sk=field_send_se_[dir][x1face][4];
  ek=field_send_se_[dir][x1face][5];
  for (int k=sk; k<=ek; ++k) {
    for (int j=sj; j<=ej; ++j) {
#pragma simd
      for (int i=si; i<=ei; ++i) {
        // buffer is always fully packed
        sendbuf[p++]=x1src(k,j,i);
      }
    }
  }
  // Set buffers; x2f
  si=field_send_se_[dir][x2face][0];
  ei=field_send_se_[dir][x2face][1];
  sj=field_send_se_[dir][x2face][2];
  ej=field_send_se_[dir][x2face][3];
  sk=field_send_se_[dir][x2face][4];
  ek=field_send_se_[dir][x2face][5];
  for (int k=sk; k<=ek; ++k) {
    for (int j=sj; j<=ej; ++j) {
#pragma simd
      for (int i=si; i<=ei; ++i) {
        // buffer is always fully packed
        sendbuf[p++]=x2src(k,j,i);
      }
    }
  }
  // Set buffers; x3f
  si=field_send_se_[dir][x3face][0];
  ei=field_send_se_[dir][x3face][1];
  sj=field_send_se_[dir][x3face][2];
  ej=field_send_se_[dir][x3face][3];
  sk=field_send_se_[dir][x3face][4];
  ek=field_send_se_[dir][x3face][5];
  for (int k=sk; k<=ek; ++k) {
    for (int j=sj; j<=ej; ++j) {
#pragma simd
      for (int i=si; i<=ei; ++i) {
        // buffer is always fully packed
        sendbuf[p++]=x3src(k,j,i);
      }
    }
  }

  // Send the buffer; modify this for MPI and AMR
  if(pmb->neighbor[dir][0][0].rank == myrank) // myrank
  {
    while(pbl!=NULL)
    {
      if(pbl->gid==pmb->neighbor[dir][0][0].gid)
        break;
      pbl=pbl->next;
    }
    if(pbl==NULL)
    {
      msg << "### FATAL ERROR in SetFieldBoundary" << std::endl
          << "In boundary " << dir << " of block "<<  pmb->gid
          << ", neighbor block " << pmb->neighbor[dir][0][0].gid
          << " is not found on the same process." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    std::memcpy(pbl->pbval->field_recv_[oside], field_send_[dir],
                field_bufsize_[dir]*sizeof(Real));
    pbl->pbval->field_flag_[oside][0][0]=true; // the other side
  }
  else // MPI
  {
      msg << "### FATAL ERROR in SetFieldBoundary" << std::endl
          << "MPI is not implemented yet!!" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  return;
}

//--------------------------------------------------------------------------------------
//! \fn bool BoundaryValues::ReceiveAndSetFieldBoundary(enum direction dir,
//                                                      InterfaceField &dst)
//  \brief load boundary buffer for x1 direction into the array
bool BoundaryValues::ReceiveAndSetFieldBoundary(enum direction dir, InterfaceField &dst)
{
  MeshBlock *pmb=pmy_mblock_;
  std::stringstream msg;
  Real *recvbuf=field_recv_[dir];
  AthenaArray<Real>& x1dst=dst.x1f;
  AthenaArray<Real>& x2dst=dst.x2f;
  AthenaArray<Real>& x3dst=dst.x3f;
  int si, sj, sk, ei, ej, ek;

  if(pmb->neighbor[dir][0][0].gid==-1) // physical boundary
    FieldBoundary_[dir](pmb,dst);
  else // block boundary
  {
    if(field_flag_[dir][0][0] == false)
      return false; // return if it is not ready yet

    // Load buffers; x1f
    int p=0;
    si=field_recv_se_[dir][x1face][0];
    ei=field_recv_se_[dir][x1face][1];
    sj=field_recv_se_[dir][x1face][2];
    ej=field_recv_se_[dir][x1face][3];
    sk=field_recv_se_[dir][x1face][4];
    ek=field_recv_se_[dir][x1face][5];
    for (int k=sk; k<=ek; ++k) {
      for (int j=sj; j<=ej; ++j) {
#pragma simd
        for (int i=si; i<=ei; ++i) {
          // buffer is always fully packed
          x1dst(k,j,i)=recvbuf[p++];
        }
      }
    }
    // Load buffers; x2f
    si=field_recv_se_[dir][x2face][0];
    ei=field_recv_se_[dir][x2face][1];
    sj=field_recv_se_[dir][x2face][2];
    ej=field_recv_se_[dir][x2face][3];
    sk=field_recv_se_[dir][x2face][4];
    ek=field_recv_se_[dir][x2face][5];
    for (int k=sk; k<=ek; ++k) {
      for (int j=sj; j<=ej; ++j) {
#pragma simd
        for (int i=si; i<=ei; ++i) {
          // buffer is always fully packed
          x2dst(k,j,i)=recvbuf[p++];
        }
      }
    }
    // Load buffers; x3f
    si=field_recv_se_[dir][x3face][0];
    ei=field_recv_se_[dir][x3face][1];
    sj=field_recv_se_[dir][x3face][2];
    ej=field_recv_se_[dir][x3face][3];
    sk=field_recv_se_[dir][x3face][4];
    ek=field_recv_se_[dir][x3face][5];
    for (int k=sk; k<=ek; ++k) {
      for (int j=sj; j<=ej; ++j) {
#pragma simd
        for (int i=si; i<=ei; ++i) {
          // buffer is always fully packed
          x3dst(k,j,i)=recvbuf[p++];
        }
      }
    }
    field_flag_[dir][0][0] = false; // clear the flag
  }

  return true;
}


void InitBoundaryBuffer(int nx1, int nx2, int nx3)
{
  int is, ie, js, je, ks, ke;

  is = NGHOST;
  ie = is + nx1 - 1;

  if (nx2 > 1) {
    js = NGHOST;
    je = js + nx2 - 1;
  } else {
    js = je = 0;
  }

  if (nx3 > 1) {
    ks = NGHOST;
    ke = ks + nx3 - 1;
  } else {
    ks = ke = 0;
  }
  fluid_send_se_[inner_x1][0]=is;
  fluid_send_se_[inner_x1][1]=is+NGHOST-1;
  fluid_send_se_[inner_x1][2]=js;
  fluid_send_se_[inner_x1][3]=je;
  fluid_send_se_[inner_x1][4]=ks;
  fluid_send_se_[inner_x1][5]=ke;

  fluid_send_se_[outer_x1][0]=ie-NGHOST+1;
  fluid_send_se_[outer_x1][1]=ie;
  fluid_send_se_[outer_x1][2]=js;
  fluid_send_se_[outer_x1][3]=je;
  fluid_send_se_[outer_x1][4]=ks;
  fluid_send_se_[outer_x1][5]=ke;

  fluid_send_se_[inner_x2][0]=0;
  fluid_send_se_[inner_x2][1]=ie+NGHOST;
  fluid_send_se_[inner_x2][2]=js;
  fluid_send_se_[inner_x2][3]=js+NGHOST-1;
  fluid_send_se_[inner_x2][4]=ks;
  fluid_send_se_[inner_x2][5]=ke;

  fluid_send_se_[outer_x2][0]=0;
  fluid_send_se_[outer_x2][1]=ie+NGHOST;
  fluid_send_se_[outer_x2][2]=je-NGHOST+1;
  fluid_send_se_[outer_x2][3]=je;
  fluid_send_se_[outer_x2][4]=ks;
  fluid_send_se_[outer_x2][5]=ke;

  fluid_send_se_[inner_x3][0]=0;
  fluid_send_se_[inner_x3][1]=ie+NGHOST;
  fluid_send_se_[inner_x3][2]=0;
  fluid_send_se_[inner_x3][3]=je+NGHOST;
  fluid_send_se_[inner_x3][4]=ks;
  fluid_send_se_[inner_x3][5]=ks+NGHOST-1;

  fluid_send_se_[outer_x3][0]=0;
  fluid_send_se_[outer_x3][1]=ie+NGHOST;
  fluid_send_se_[outer_x3][2]=0;
  fluid_send_se_[outer_x3][3]=je+NGHOST;
  fluid_send_se_[outer_x3][4]=ke-NGHOST+1;
  fluid_send_se_[outer_x3][5]=ke;


  fluid_recv_se_[inner_x1][0]=is-NGHOST;
  fluid_recv_se_[inner_x1][1]=is-NGHOST+1;
  fluid_recv_se_[inner_x1][2]=js;
  fluid_recv_se_[inner_x1][3]=je;
  fluid_recv_se_[inner_x1][4]=ks;
  fluid_recv_se_[inner_x1][5]=ke;

  fluid_recv_se_[outer_x1][0]=ie+1;
  fluid_recv_se_[outer_x1][1]=ie+NGHOST;
  fluid_recv_se_[outer_x1][2]=js;
  fluid_recv_se_[outer_x1][3]=je;
  fluid_recv_se_[outer_x1][4]=ks;
  fluid_recv_se_[outer_x1][5]=ke;

  fluid_recv_se_[inner_x2][0]=0;
  fluid_recv_se_[inner_x2][1]=ie+NGHOST;
  fluid_recv_se_[inner_x2][2]=js-NGHOST;
  fluid_recv_se_[inner_x2][3]=js-NGHOST+1;
  fluid_recv_se_[inner_x2][4]=ks;
  fluid_recv_se_[inner_x2][5]=ke;

  fluid_recv_se_[outer_x2][0]=0;
  fluid_recv_se_[outer_x2][1]=ie+NGHOST;
  fluid_recv_se_[outer_x2][2]=je+1;
  fluid_recv_se_[outer_x2][3]=je+NGHOST;
  fluid_recv_se_[outer_x2][4]=ks;
  fluid_recv_se_[outer_x2][5]=ke;

  fluid_recv_se_[inner_x3][0]=0;
  fluid_recv_se_[inner_x3][1]=ie+NGHOST;
  fluid_recv_se_[inner_x3][2]=0;
  fluid_recv_se_[inner_x3][3]=je+NGHOST;
  fluid_recv_se_[inner_x3][4]=ks-NGHOST;
  fluid_recv_se_[inner_x3][5]=ks-NGHOST+1;

  fluid_recv_se_[outer_x3][0]=0;
  fluid_recv_se_[outer_x3][1]=ie+NGHOST;
  fluid_recv_se_[outer_x3][2]=0;
  fluid_recv_se_[outer_x3][3]=je+NGHOST;
  fluid_recv_se_[outer_x3][4]=ke+1;
  fluid_recv_se_[outer_x3][5]=ke+NGHOST;

  fluid_bufsize_[inner_x1]=NGHOST*nx2*nx3*NFLUID;
  fluid_bufsize_[outer_x1]=NGHOST*nx2*nx3*NFLUID;
  fluid_bufsize_[inner_x2]=(nx1+2*NGHOST)*NGHOST*nx3*NFLUID;
  fluid_bufsize_[outer_x2]=(nx1+2*NGHOST)*NGHOST*nx3*NFLUID;
  fluid_bufsize_[inner_x3]=(nx1+2*NGHOST)*(nx2+2*NGHOST)*NGHOST*NFLUID;
  fluid_bufsize_[outer_x3]=(nx1+2*NGHOST)*(nx2+2*NGHOST)*NGHOST*NFLUID;

  if (MAGNETIC_FIELDS_ENABLED) {
    field_send_se_[inner_x1][x1face][0]=is+1;
    field_send_se_[inner_x1][x1face][1]=is+NGHOST;
    field_send_se_[inner_x1][x1face][2]=js;
    field_send_se_[inner_x1][x1face][3]=je;
    field_send_se_[inner_x1][x1face][4]=ks;
    field_send_se_[inner_x1][x1face][5]=ke;

    field_send_se_[inner_x1][x2face][0]=is;
    field_send_se_[inner_x1][x2face][1]=is+NGHOST-1;
    field_send_se_[inner_x1][x2face][2]=js;
    field_send_se_[inner_x1][x2face][3]=je+1;
    field_send_se_[inner_x1][x2face][4]=ks;
    field_send_se_[inner_x1][x2face][5]=ke;

    field_send_se_[inner_x1][x3face][0]=is;
    field_send_se_[inner_x1][x3face][1]=is+NGHOST-1;
    field_send_se_[inner_x1][x3face][2]=js;
    field_send_se_[inner_x1][x3face][3]=je;
    field_send_se_[inner_x1][x3face][4]=ks;
    field_send_se_[inner_x1][x3face][5]=ke+1;

    field_send_se_[outer_x1][x1face][0]=ie-NGHOST+1;
    field_send_se_[outer_x1][x1face][1]=ie;
    field_send_se_[outer_x1][x1face][2]=js;
    field_send_se_[outer_x1][x1face][3]=je;
    field_send_se_[outer_x1][x1face][4]=ks;
    field_send_se_[outer_x1][x1face][5]=ke;

    field_send_se_[outer_x1][x2face][0]=ie-NGHOST+1;
    field_send_se_[outer_x1][x2face][1]=ie;
    field_send_se_[outer_x1][x2face][2]=js;
    field_send_se_[outer_x1][x2face][3]=je+1;
    field_send_se_[outer_x1][x2face][4]=ks;
    field_send_se_[outer_x1][x2face][5]=ke;

    field_send_se_[outer_x1][x3face][0]=ie-NGHOST+1;
    field_send_se_[outer_x1][x3face][1]=ie;
    field_send_se_[outer_x1][x3face][2]=js;
    field_send_se_[outer_x1][x3face][3]=je;
    field_send_se_[outer_x1][x3face][4]=ks;
    field_send_se_[outer_x1][x3face][5]=ke+1;

    field_send_se_[inner_x2][x1face][0]=0;
    field_send_se_[inner_x2][x1face][1]=ie+NGHOST+1;
    field_send_se_[inner_x2][x1face][2]=js;
    field_send_se_[inner_x2][x1face][3]=js+NGHOST-1;
    field_send_se_[inner_x2][x1face][4]=ks;
    field_send_se_[inner_x2][x1face][5]=ke;

    field_send_se_[inner_x2][x2face][0]=0;
    field_send_se_[inner_x2][x2face][1]=ie+NGHOST;
    field_send_se_[inner_x2][x2face][2]=js+1;
    field_send_se_[inner_x2][x2face][3]=js+NGHOST;
    field_send_se_[inner_x2][x2face][4]=ks;
    field_send_se_[inner_x2][x2face][5]=ke;

    field_send_se_[inner_x2][x3face][0]=0;
    field_send_se_[inner_x2][x3face][1]=ie+NGHOST;
    field_send_se_[inner_x2][x3face][2]=js;
    field_send_se_[inner_x2][x3face][3]=js+NGHOST-1;
    field_send_se_[inner_x2][x3face][4]=ks;
    field_send_se_[inner_x2][x3face][5]=ke+1;

    field_send_se_[outer_x2][x1face][0]=0;
    field_send_se_[outer_x2][x1face][1]=ie+NGHOST+1;
    field_send_se_[outer_x2][x1face][2]=je-NGHOST+1;
    field_send_se_[outer_x2][x1face][3]=je;
    field_send_se_[outer_x2][x1face][4]=ks;
    field_send_se_[outer_x2][x1face][5]=ke;

    field_send_se_[outer_x2][x2face][0]=0;
    field_send_se_[outer_x2][x2face][1]=ie+NGHOST;
    field_send_se_[outer_x2][x2face][2]=je-NGHOST+1;
    field_send_se_[outer_x2][x2face][3]=je;
    field_send_se_[outer_x2][x2face][4]=ks;
    field_send_se_[outer_x2][x2face][5]=ke;

    field_send_se_[outer_x2][x3face][0]=0;
    field_send_se_[outer_x2][x3face][1]=ie+NGHOST;
    field_send_se_[outer_x2][x3face][2]=je-NGHOST+1;
    field_send_se_[outer_x2][x3face][3]=je;
    field_send_se_[outer_x2][x3face][4]=ks;
    field_send_se_[outer_x2][x3face][5]=ke+1;

    field_send_se_[inner_x3][x1face][0]=0;
    field_send_se_[inner_x3][x1face][1]=ie+NGHOST+1;
    field_send_se_[inner_x3][x1face][2]=0;
    field_send_se_[inner_x3][x1face][3]=je+NGHOST;
    field_send_se_[inner_x3][x1face][4]=ks;
    field_send_se_[inner_x3][x1face][5]=ks+NGHOST-1;

    field_send_se_[inner_x3][x2face][0]=0;
    field_send_se_[inner_x3][x2face][1]=ie+NGHOST;
    field_send_se_[inner_x3][x2face][2]=0;
    field_send_se_[inner_x3][x2face][3]=je+NGHOST+1;
    field_send_se_[inner_x3][x2face][4]=ks;
    field_send_se_[inner_x3][x2face][5]=ks+NGHOST-1;

    field_send_se_[inner_x3][x3face][0]=0;
    field_send_se_[inner_x3][x3face][1]=ie+NGHOST;
    field_send_se_[inner_x3][x3face][2]=0;
    field_send_se_[inner_x3][x3face][3]=je+NGHOST+1;
    field_send_se_[inner_x3][x3face][4]=ks+1;
    field_send_se_[inner_x3][x3face][5]=ks+NGHOST;

    field_send_se_[outer_x3][x1face][0]=0;
    field_send_se_[outer_x3][x1face][1]=ie+NGHOST+1;
    field_send_se_[outer_x3][x1face][2]=0;
    field_send_se_[outer_x3][x1face][3]=je+NGHOST;
    field_send_se_[outer_x3][x1face][4]=ke-NGHOST+1;
    field_send_se_[outer_x3][x1face][5]=ke;

    field_send_se_[outer_x3][x2face][0]=0;
    field_send_se_[outer_x3][x2face][1]=ie+NGHOST;
    field_send_se_[outer_x3][x2face][2]=0;
    field_send_se_[outer_x3][x2face][3]=je+NGHOST+1;
    field_send_se_[outer_x3][x2face][4]=ke-NGHOST+1;
    field_send_se_[outer_x3][x2face][5]=ke;

    field_send_se_[outer_x3][x3face][0]=0;
    field_send_se_[outer_x3][x3face][1]=ie+NGHOST;
    field_send_se_[outer_x3][x3face][2]=0;
    field_send_se_[outer_x3][x3face][3]=je+NGHOST+1;
    field_send_se_[outer_x3][x3face][4]=ke-NGHOST+1;
    field_send_se_[outer_x3][x3face][5]=ke;

    field_recv_se_[inner_x1][x1face][0]=is-NGHOST;
    field_recv_se_[inner_x1][x1face][1]=is-1;
    field_recv_se_[inner_x1][x1face][2]=js;
    field_recv_se_[inner_x1][x1face][3]=je;
    field_recv_se_[inner_x1][x1face][4]=ks;
    field_recv_se_[inner_x1][x1face][5]=ke;

    field_recv_se_[inner_x1][x2face][0]=is-NGHOST;
    field_recv_se_[inner_x1][x2face][1]=is-1;
    field_recv_se_[inner_x1][x2face][2]=js;
    field_recv_se_[inner_x1][x2face][3]=je+1;
    field_recv_se_[inner_x1][x2face][4]=ks;
    field_recv_se_[inner_x1][x2face][5]=ke;

    field_recv_se_[inner_x1][x3face][0]=is-NGHOST;
    field_recv_se_[inner_x1][x3face][1]=is-1;
    field_recv_se_[inner_x1][x3face][2]=js;
    field_recv_se_[inner_x1][x3face][3]=je;
    field_recv_se_[inner_x1][x3face][4]=ks;
    field_recv_se_[inner_x1][x3face][5]=ke+1;

    field_recv_se_[outer_x1][x1face][0]=ie+2;
    field_recv_se_[outer_x1][x1face][1]=ie+NGHOST+1;
    field_recv_se_[outer_x1][x1face][2]=js;
    field_recv_se_[outer_x1][x1face][3]=je;
    field_recv_se_[outer_x1][x1face][4]=ks;
    field_recv_se_[outer_x1][x1face][5]=ke;

    field_recv_se_[outer_x1][x2face][0]=ie+1;
    field_recv_se_[outer_x1][x2face][1]=ie+NGHOST;
    field_recv_se_[outer_x1][x2face][2]=js;
    field_recv_se_[outer_x1][x2face][3]=je+1;
    field_recv_se_[outer_x1][x2face][4]=ks;
    field_recv_se_[outer_x1][x2face][5]=ke;

    field_recv_se_[outer_x1][x3face][0]=ie+1;
    field_recv_se_[outer_x1][x3face][1]=ie+NGHOST;
    field_recv_se_[outer_x1][x3face][2]=js;
    field_recv_se_[outer_x1][x3face][3]=je;
    field_recv_se_[outer_x1][x3face][4]=ks;
    field_recv_se_[outer_x1][x3face][5]=ke+1;

    field_recv_se_[inner_x2][x1face][0]=0;
    field_recv_se_[inner_x2][x1face][1]=ie+NGHOST+1;
    field_recv_se_[inner_x2][x1face][2]=js-NGHOST;
    field_recv_se_[inner_x2][x1face][3]=js-1;
    field_recv_se_[inner_x2][x1face][4]=ks;
    field_recv_se_[inner_x2][x1face][5]=ke;

    field_recv_se_[inner_x2][x2face][0]=0;
    field_recv_se_[inner_x2][x2face][1]=ie+NGHOST;
    field_recv_se_[inner_x2][x2face][2]=js-NGHOST;
    field_recv_se_[inner_x2][x2face][3]=js-1;
    field_recv_se_[inner_x2][x2face][4]=ks;
    field_recv_se_[inner_x2][x2face][5]=ke;

    field_recv_se_[inner_x2][x3face][0]=0;
    field_recv_se_[inner_x2][x3face][1]=ie+NGHOST;
    field_recv_se_[inner_x2][x3face][2]=js-NGHOST;
    field_recv_se_[inner_x2][x3face][3]=js-1;
    field_recv_se_[inner_x2][x3face][4]=ks;
    field_recv_se_[inner_x2][x3face][5]=ke+1;

    field_recv_se_[outer_x2][x1face][0]=0;
    field_recv_se_[outer_x2][x1face][1]=ie+NGHOST+1;
    field_recv_se_[outer_x2][x1face][2]=je+1;
    field_recv_se_[outer_x2][x1face][3]=je+NGHOST;
    field_recv_se_[outer_x2][x1face][4]=ks;
    field_recv_se_[outer_x2][x1face][5]=ke;

    field_recv_se_[outer_x2][x2face][0]=0;
    field_recv_se_[outer_x2][x2face][1]=ie+NGHOST;
    field_recv_se_[outer_x2][x2face][2]=je+2;
    field_recv_se_[outer_x2][x2face][3]=je+NGHOST+1;
    field_recv_se_[outer_x2][x2face][4]=ks;
    field_recv_se_[outer_x2][x2face][5]=ke;

    field_recv_se_[outer_x2][x3face][0]=0;
    field_recv_se_[outer_x2][x3face][1]=ie+NGHOST;
    field_recv_se_[outer_x2][x3face][2]=je+1;
    field_recv_se_[outer_x2][x3face][3]=je+NGHOST;
    field_recv_se_[outer_x2][x3face][4]=ks;
    field_recv_se_[outer_x2][x3face][5]=ke+1;

    field_recv_se_[inner_x3][x1face][0]=0;
    field_recv_se_[inner_x3][x1face][1]=ie+NGHOST+1;
    field_recv_se_[inner_x3][x1face][2]=0;
    field_recv_se_[inner_x3][x1face][3]=je+NGHOST;
    field_recv_se_[inner_x3][x1face][4]=ks-NGHOST;
    field_recv_se_[inner_x3][x1face][5]=ks-1;

    field_recv_se_[inner_x3][x2face][0]=0;
    field_recv_se_[inner_x3][x2face][1]=ie+NGHOST;
    field_recv_se_[inner_x3][x2face][2]=0;
    field_recv_se_[inner_x3][x2face][3]=je+NGHOST+1;
    field_recv_se_[inner_x3][x2face][4]=ks-NGHOST;
    field_recv_se_[inner_x3][x2face][5]=ks-1;

    field_recv_se_[inner_x3][x3face][0]=0;
    field_recv_se_[inner_x3][x3face][1]=ie+NGHOST;
    field_recv_se_[inner_x3][x3face][2]=0;
    field_recv_se_[inner_x3][x3face][3]=je+NGHOST+1;
    field_recv_se_[inner_x3][x3face][4]=ks-NGHOST;
    field_recv_se_[inner_x3][x3face][5]=ks-1;

    field_recv_se_[outer_x3][x1face][0]=0;
    field_recv_se_[outer_x3][x1face][1]=ie+NGHOST+1;
    field_recv_se_[outer_x3][x1face][2]=0;
    field_recv_se_[outer_x3][x1face][3]=je+NGHOST;
    field_recv_se_[outer_x3][x1face][4]=ke+1;
    field_recv_se_[outer_x3][x1face][5]=ke+NGHOST;

    field_recv_se_[outer_x3][x2face][0]=0;
    field_recv_se_[outer_x3][x2face][1]=ie+NGHOST;
    field_recv_se_[outer_x3][x2face][2]=0;
    field_recv_se_[outer_x3][x2face][3]=je+NGHOST+1;
    field_recv_se_[outer_x3][x2face][4]=ke+1;
    field_recv_se_[outer_x3][x2face][5]=ke+NGHOST;

    field_recv_se_[outer_x3][x3face][0]=0;
    field_recv_se_[outer_x3][x3face][1]=ie+NGHOST;
    field_recv_se_[outer_x3][x3face][2]=0;
    field_recv_se_[outer_x3][x3face][3]=je+NGHOST+1;
    field_recv_se_[outer_x3][x3face][4]=ke+2;
    field_recv_se_[outer_x3][x3face][5]=ke+NGHOST+1;

    field_bufsize_[inner_x1]=NGHOST*(nx2*nx3+(nx2+1)*nx3+nx2*(nx3+1));
    field_bufsize_[outer_x1]=NGHOST*(nx2*nx3+(nx2+1)*nx3+nx2*(nx3+1));
    field_bufsize_[inner_x2]=NGHOST*((nx1+2*NGHOST)*nx3+(nx1+2*NGHOST+1)*nx3
                            +(nx1+2*NGHOST)*(nx3+1));
    field_bufsize_[outer_x2]=NGHOST*((nx1+2*NGHOST)*nx3+(nx1+2*NGHOST+1)*nx3
                            +(nx1+2*NGHOST)*(nx3+1));
    field_bufsize_[inner_x3]=NGHOST*((nx1+2*NGHOST+1)*(nx2+2*NGHOST)
                      +(nx1+2*NGHOST)*(nx2+2*NGHOST+1)+(nx1+2*NGHOST)*(nx2+2*NGHOST));
    field_bufsize_[outer_x3]=NGHOST*((nx1+2*NGHOST+1)*(nx2+2*NGHOST)
                      +(nx1+2*NGHOST)*(nx2+2*NGHOST+1)+(nx1+2*NGHOST)*(nx2+2*NGHOST));
  }
  return;
}

