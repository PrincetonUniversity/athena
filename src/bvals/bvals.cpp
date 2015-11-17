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
//! \file bvals.cpp
//  \brief implements functions that initialize/apply BCs on each dir
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
#include "../mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../coordinates/coordinates.hpp"
#include "../parameter_input.hpp"

// this class header
#include "bvals.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

static NeighborIndexes ni_[56];
static int bufid_[56];

// BoundaryValues constructor - sets functions for the appropriate
// boundary conditions at each of the 6 dirs of a MeshBlock

BoundaryValues::BoundaryValues(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_mblock_ = pmb;
  int cng=pmb->cnghost, cng1=0, cng2=0, cng3=0;
  if(pmb->block_size.nx2>1) cng1=cng, cng2=cng;
  if(pmb->block_size.nx3>1) cng3=cng;
  int f2d=0, f3d=0;
  if(pmb->block_size.nx2 > 1) f2d=1;
  if(pmb->block_size.nx3 > 1) f3d=1;

// Set BC functions for each of the 6 boundaries in turn -------------------------------
// Inner x1
  nface_=2; nedge_=0;
  switch(pmb->block_bcs[inner_x1]){
    case 1:
      HydroBoundary_[inner_x1] = ReflectInnerX1;
      FieldBoundary_[inner_x1] = ReflectInnerX1;
      break;
    case 2:
      HydroBoundary_[inner_x1] = OutflowInnerX1;
      FieldBoundary_[inner_x1] = OutflowInnerX1;
      break;
    case -1: // block boundary
    case 3: // do nothing, useful for user-enrolled BCs
    case 4: // periodic boundary
      HydroBoundary_[inner_x1] = NULL;
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
      HydroBoundary_[outer_x1] = ReflectOuterX1;
      FieldBoundary_[outer_x1] = ReflectOuterX1;
      break;
    case 2:
      HydroBoundary_[outer_x1] = OutflowOuterX1;
      FieldBoundary_[outer_x1] = OutflowOuterX1;
      break;
    case -1: // block boundary
    case 3: // do nothing, useful for user-enrolled BCs
    case 4: // periodic boundary
      HydroBoundary_[outer_x1] = NULL;
      FieldBoundary_[outer_x1] = NULL;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
          << "Flag ox1_bc=" << pmb->block_bcs[outer_x1] << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }

  if (pmb->block_size.nx2 > 1) {
    nface_=4; nedge_=4;
// Inner x2
    switch(pmb->block_bcs[inner_x2]){
      case 1:
        HydroBoundary_[inner_x2] = ReflectInnerX2;
        FieldBoundary_[inner_x2] = ReflectInnerX2;
        break;
      case 2:
        HydroBoundary_[inner_x2] = OutflowInnerX2;
        FieldBoundary_[inner_x2] = OutflowInnerX2;
        break;
      case -1: // block boundary
      case 3: // do nothing, useful for user-enrolled BCs
      case 4: // periodic boundary
        HydroBoundary_[inner_x2] = NULL;
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
        HydroBoundary_[outer_x2] = ReflectOuterX2;
        FieldBoundary_[outer_x2] = ReflectOuterX2;
        break;
      case 2:
        HydroBoundary_[outer_x2] = OutflowOuterX2;
        FieldBoundary_[outer_x2] = OutflowOuterX2;
        break;
      case -1: // block boundary
      case 3: // do nothing, useful for user-enrolled BCs
      case 4: // periodic boundary
        HydroBoundary_[outer_x2] = NULL;
        FieldBoundary_[outer_x2] = NULL;
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
            << "Flag ox2_bc=" << pmb->block_bcs[outer_x2] << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
    }
  }

  if (pmb->block_size.nx3 > 1) {
    nface_=6; nedge_=12;
// Inner x3
    switch(pmb->block_bcs[inner_x3]){
      case 1:
        HydroBoundary_[inner_x3] = ReflectInnerX3;
        FieldBoundary_[inner_x3] = ReflectInnerX3;
        break;
      case 2:
        HydroBoundary_[inner_x3] = OutflowInnerX3;
        FieldBoundary_[inner_x3] = OutflowInnerX3;
        break;
      case -1: // block boundary
      case 3: // do nothing, useful for user-enrolled BCs
      case 4: // periodic boundary
        HydroBoundary_[inner_x3] = NULL;
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
        HydroBoundary_[outer_x3] = ReflectOuterX3;
        FieldBoundary_[outer_x3] = ReflectOuterX3;
        break;
      case 2:
        HydroBoundary_[outer_x3] = OutflowOuterX3;
        FieldBoundary_[outer_x3] = OutflowOuterX3;
        break;
      case -1: // block boundary
      case 3: // do nothing, useful for user-enrolled BCs
      case 4: // periodic boundary
        HydroBoundary_[outer_x3] = NULL;
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
      hydro_flag_[l][i]=boundary_waiting;
      field_flag_[l][i]=boundary_waiting;
      hydro_send_[l][i]=NULL;
      hydro_recv_[l][i]=NULL;
      field_send_[l][i]=NULL;
      field_recv_[l][i]=NULL;
#ifdef MPI_PARALLEL
      req_hydro_send_[l][i]=MPI_REQUEST_NULL;
      req_hydro_recv_[l][i]=MPI_REQUEST_NULL;
#endif
    }
    for(int i=0;i<48;i++){
      emfcor_send_[l][i]=NULL;
      emfcor_recv_[l][i]=NULL;
      emfcor_flag_[l][i]=boundary_waiting;
#ifdef MPI_PARALLEL
      req_emfcor_send_[l][i]=MPI_REQUEST_NULL;
      req_emfcor_recv_[l][i]=MPI_REQUEST_NULL;
#endif
    }
    for(int i=0;i<6;i++){
      flcor_send_[l][i]=NULL;
#ifdef MPI_PARALLEL
      req_flcor_send_[l][i]=MPI_REQUEST_NULL;
#endif
      for(int j=0;j<=1;j++) {
        for(int k=0;k<=1;k++) {
          flcor_recv_[l][i][j][k]=NULL;
          flcor_flag_[l][i][j][k]=boundary_waiting;
#ifdef MPI_PARALLEL
          req_flcor_recv_[l][i][j][k]=MPI_REQUEST_NULL;
#endif
        }
      }
    }
  }
  // Allocate Buffers
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
      size*=NHYDRO;
      hydro_send_[l][n]=new Real[size];
      hydro_recv_[l][n]=new Real[size];
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for(int l=0;l<NSTEP;l++) {
      for(int n=0;n<pmb->pmy_mesh->maxneighbor_;n++) {
        int size1=((ni_[n].ox1==0)?(pmb->block_size.nx1+1):NGHOST)
                 *((ni_[n].ox2==0)?(pmb->block_size.nx2):NGHOST)
                 *((ni_[n].ox3==0)?(pmb->block_size.nx3):NGHOST);
        int size2=((ni_[n].ox1==0)?(pmb->block_size.nx1):NGHOST)
                 *((ni_[n].ox2==0)?(pmb->block_size.nx2+f2d):NGHOST)
                 *((ni_[n].ox3==0)?(pmb->block_size.nx3):NGHOST);
        int size3=((ni_[n].ox1==0)?(pmb->block_size.nx1):NGHOST)
                 *((ni_[n].ox2==0)?(pmb->block_size.nx2):NGHOST)
                 *((ni_[n].ox3==0)?(pmb->block_size.nx3+f3d):NGHOST);
        int size=size1+size2+size3;
        if(pmb->pmy_mesh->multilevel==true) {
          if(ni_[n].type!=neighbor_face) {
            if(ni_[n].ox1!=0) size1=size1/NGHOST*(NGHOST+1);
            if(ni_[n].ox2!=0) size2=size2/NGHOST*(NGHOST+1);
            if(ni_[n].ox3!=0) size3=size3/NGHOST*(NGHOST+1);
          }
          size=size1+size2+size3;
          int f2c1=((ni_[n].ox1==0)?((pmb->block_size.nx1+1)/2+1):cng)
                  *((ni_[n].ox2==0)?((pmb->block_size.nx2+1)/2):cng)
                  *((ni_[n].ox3==0)?((pmb->block_size.nx3+1)/2):cng);
          int f2c2=((ni_[n].ox1==0)?((pmb->block_size.nx1+1)/2):cng)
                  *((ni_[n].ox2==0)?((pmb->block_size.nx2+1)/2+f2d):cng)
                  *((ni_[n].ox3==0)?((pmb->block_size.nx3+1)/2):cng);
          int f2c3=((ni_[n].ox1==0)?((pmb->block_size.nx1+1)/2):cng)
                  *((ni_[n].ox2==0)?((pmb->block_size.nx2+1)/2):cng)
                  *((ni_[n].ox3==0)?((pmb->block_size.nx3+1)/2+f3d):cng);
          if(ni_[n].type!=neighbor_face) {
            if(ni_[n].ox1!=0) f2c1=f2c1/cng*(cng+1);
            if(ni_[n].ox2!=0) f2c2=f2c2/cng*(cng+1);
            if(ni_[n].ox3!=0) f2c3=f2c3/cng*(cng+1);
          }
          int fsize=f2c1+f2c2+f2c3;
          int c2f1=((ni_[n].ox1==0)?((pmb->block_size.nx1+1)/2+1+cng):cng)
                  *((ni_[n].ox2==0)?((pmb->block_size.nx2+1)/2+cng*f2d):cng)
                  *((ni_[n].ox3==0)?((pmb->block_size.nx3+1)/2+cng*f3d):cng);
          int c2f2=((ni_[n].ox1==0)?((pmb->block_size.nx1+1)/2+cng):cng)
                  *((ni_[n].ox2==0)?((pmb->block_size.nx2+1)/2+f2d+cng*f2d):cng)
                  *((ni_[n].ox3==0)?((pmb->block_size.nx3+1)/2+cng*f3d):cng);
          int c2f3=((ni_[n].ox1==0)?((pmb->block_size.nx1+1)/2+cng):cng)
                  *((ni_[n].ox2==0)?((pmb->block_size.nx2+1)/2+f2d*cng):cng)
                  *((ni_[n].ox3==0)?((pmb->block_size.nx3+1)/2+f3d+cng*f3d):cng);
          int csize=c2f1+c2f2+c2f3;
          size=std::max(size,std::max(csize,fsize));
        }
        field_send_[l][n]=new Real[size];
        field_recv_[l][n]=new Real[size];

        // allocate EMF correction buffer
        if(ni_[n].type==neighbor_face) {
          if(pmb->block_size.nx3>1) { // 3D
            if(ni_[n].ox1!=0)
              size=(pmb->block_size.nx2+1)*(pmb->block_size.nx3)
                  +(pmb->block_size.nx2)*(pmb->block_size.nx3+1);
            else if(ni_[n].ox2!=0)
              size=(pmb->block_size.nx1+1)*(pmb->block_size.nx3)
                  +(pmb->block_size.nx1)*(pmb->block_size.nx3+1);
            else
              size=(pmb->block_size.nx1+1)*(pmb->block_size.nx2)
                  +(pmb->block_size.nx1)*(pmb->block_size.nx2+1);
          }
          else if(pmb->block_size.nx2>1) { // 2D
            if(ni_[n].ox1!=0)
              size=(pmb->block_size.nx2+1)+pmb->block_size.nx2;
            else 
              size=(pmb->block_size.nx1+1)+pmb->block_size.nx1;
          }
          else // 1D
            size=2;
        }
        else if(ni_[n].type==neighbor_edge) {
          if(pmb->block_size.nx3>1) { // 3D
            if(ni_[n].ox3==0) size=pmb->block_size.nx3;
            if(ni_[n].ox2==0) size=pmb->block_size.nx2;
            if(ni_[n].ox1==0) size=pmb->block_size.nx1;
          }
          else if(pmb->block_size.nx2>1)
            size=1;
        }
        else continue;
        emfcor_send_[l][n]=new Real[size];
        emfcor_recv_[l][n]=new Real[size];
      }
    }
  }

  if(pmb->pmy_mesh->multilevel==true) { // SMR or AMR
    // allocate arrays for volumes in the finer level
    int nc1=pmb->block_size.nx1+2*NGHOST;
    int nc2=pmb->block_size.nx2+2*NGHOST;
    int nc3=pmb->block_size.nx3+2*NGHOST;
    fvol_[0][0].NewAthenaArray(nc1+1);
    fvol_[0][1].NewAthenaArray(nc1+1);
    fvol_[1][0].NewAthenaArray(nc1+1);
    fvol_[1][1].NewAthenaArray(nc1+1);
    sarea_[0].NewAthenaArray(nc1);
    sarea_[1].NewAthenaArray(nc1);
    int size[6], im, jm, km;
    // allocate flux correction buffer
    size[0]=size[1]=(pmb->block_size.nx2+1)/2*(pmb->block_size.nx3+1)/2*NHYDRO;
    size[2]=size[3]=(pmb->block_size.nx1+1)/2*(pmb->block_size.nx3+1)/2*NHYDRO;
    size[4]=size[5]=(pmb->block_size.nx1+1)/2*(pmb->block_size.nx2+1)/2*NHYDRO;
    if(pmb->block_size.nx3>1) // 3D
      jm=2, km=2;
    else if(pmb->block_size.nx2>1) // 2D
      jm=1, km=2;
    else // 1D
      jm=1, km=1;
    for(int l=0;l<NSTEP;l++) {
      for(int i=0;i<nface_;i++){
        flcor_send_[l][i]=new Real[size[i]];
        for(int j=0;j<jm;j++) {
          for(int k=0;k<km;k++)
            flcor_recv_[l][i][j][k]=new Real[size[i]];
        }
      }
    }
    // allocate prolongation buffer
    int ncc1=pmb->block_size.nx1/2+2*cng;
    int ncc2=1;
    if(pmb->block_size.nx2>1) ncc2=pmb->block_size.nx2/2+2*cng;
    int ncc3=1;
    if(pmb->block_size.nx3>1) ncc3=pmb->block_size.nx3/2+2*cng;
    coarse_cons_.NewAthenaArray(NHYDRO,ncc3,ncc2,ncc1);

    if (MAGNETIC_FIELDS_ENABLED) {
      coarse_b_.x1f.NewAthenaArray(ncc3,ncc2,ncc1+1);
      coarse_b_.x2f.NewAthenaArray(ncc3,ncc2+1,ncc1);
      coarse_b_.x3f.NewAthenaArray(ncc3+1,ncc2,ncc1);
    }
  }
}

// destructor

BoundaryValues::~BoundaryValues()
{
  MeshBlock *pmb=pmy_mblock_;
  for(int l=0;l<NSTEP;l++) {
    for(int i=0;i<pmb->pmy_mesh->maxneighbor_;i++) {
      delete [] hydro_send_[l][i];
      delete [] hydro_recv_[l][i];
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for(int l=0;l<NSTEP;l++) {
      for(int i=0;i<pmb->pmy_mesh->maxneighbor_;i++) { 
        delete [] field_send_[l][i];
        delete [] field_recv_[l][i];
        if(ni_[i].type==neighbor_face || ni_[i].type==neighbor_edge) {
          delete [] emfcor_send_[l][i];
          delete [] emfcor_recv_[l][i];
        }
      }
    }
  }
  if(pmb->pmy_mesh->multilevel==true) {
    fvol_[0][0].DeleteAthenaArray();
    fvol_[0][1].DeleteAthenaArray();
    fvol_[1][0].DeleteAthenaArray();
    fvol_[1][1].DeleteAthenaArray();
    sarea_[0].DeleteAthenaArray();
    sarea_[1].DeleteAthenaArray();
    for(int l=0;l<NSTEP;l++) {
      for(int i=0;i<nface_;i++){
        delete [] flcor_send_[l][i];
        for(int j=0;j<2;j++) {
          for(int k=0;k<2;k++)
            delete [] flcor_recv_[l][i][j][k];
        }
      }
    }
    coarse_cons_.DeleteAthenaArray();
//  coarse_prim_.DeleteAthenaArray();
    if (MAGNETIC_FIELDS_ENABLED) {
      coarse_b_.x1f.DeleteAthenaArray();
      coarse_b_.x2f.DeleteAthenaArray();
      coarse_b_.x3f.DeleteAthenaArray();
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
  int myox1, myox2, myox3;
  int tag;
  int cng, cng1, cng2, cng3;
  int ssize, rsize;
  cng=cng1=pmb->cnghost;
  cng2=(pmb->block_size.nx2>1)?cng:0;
  cng3=(pmb->block_size.nx3>1)?cng:0;
  long int &lx1=pmb->loc.lx1;
  long int &lx2=pmb->loc.lx2;
  long int &lx3=pmb->loc.lx3;
  int &mylevel=pmb->loc.level;
  myox1=((int)(lx1&1L));
  myox2=((int)(lx2&1L));
  myox3=((int)(lx3&1L));
  int f2d=0, f3d=0;
  if(pmb->block_size.nx2 > 1) f2d=1;
  if(pmb->block_size.nx3 > 1) f3d=1;


  // count the number of the fine meshblocks contacting on each edge
  int eid=0;
  if(pmb->block_size.nx2 > 1) {
    for(int ox2=-1;ox2<=1;ox2+=2) {
      for(int ox1=-1;ox1<=1;ox1+=2) {
        int nis, nie, njs, nje;
        nis=std::max(ox1-1,-1), nie=std::min(ox1+1,1);
        njs=std::max(ox2-1,-1), nje=std::min(ox2+1,1);
        int nf=0, fl=mylevel;
        for(int nj=njs; nj<=nje; nj++) {
          for(int ni=nis; ni<=nie; ni++) {
            if(pmb->nblevel[1][nj+1][ni+1] > fl)
              fl++, nf=0;
            if(pmb->nblevel[1][nj+1][ni+1]==fl)
              nf++;
          }
        }
        edge_flag_[eid]=(fl==mylevel);
        nedge_fine_[eid++]=nf;
      }
    }
  }
  if(pmb->block_size.nx3 > 1) {
    for(int ox3=-1;ox3<=1;ox3+=2) {
      for(int ox1=-1;ox1<=1;ox1+=2) {
        int nis, nie, nks, nke;
        nis=std::max(ox1-1,-1), nie=std::min(ox1+1,1);
        nks=std::max(ox3-1,-1), nke=std::min(ox3+1,1);
        int nf=0, fl=mylevel;
        for(int nk=nks; nk<=nke; nk++) {
          for(int ni=nis; ni<=nie; ni++) {
            if(pmb->nblevel[nk+1][1][ni+1] > fl)
              fl++, nf=0;
            if(pmb->nblevel[nk+1][1][ni+1]==fl)
              nf++;
          }
        }
        edge_flag_[eid]=(fl==mylevel);
        nedge_fine_[eid++]=nf;
      }
    }
    for(int ox3=-1;ox3<=1;ox3+=2) {
      for(int ox2=-1;ox2<=1;ox2+=2) {
        int njs, nje, nks, nke;
        njs=std::max(ox2-1,-1), nje=std::min(ox2+1,1);
        nks=std::max(ox3-1,-1), nke=std::min(ox3+1,1);
        int nf=0, fl=mylevel;
        for(int nk=nks; nk<=nke; nk++) {
          for(int nj=njs; nj<=nje; nj++) {
            if(pmb->nblevel[nk+1][nj+1][1] > fl)
              fl++, nf=0;
            if(pmb->nblevel[nk+1][nj+1][1]==fl)
              nf++;
          }
        }
        edge_flag_[eid]=(fl==mylevel);
        nedge_fine_[eid++]=nf;
      }
    }
  }

  for(int l=0;l<NSTEP;l++) {
    for(int n=0;n<pmb->nneighbor;n++) {
      NeighborBlock& nb = pmb->neighbor[n];
      if(nb.rank!=Globals::my_rank) {
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
        ssize*=NHYDRO; rsize*=NHYDRO;
        // specify the offsets in the view point of the target block: flip ox? signs
        tag=CreateMPITag(nb.lid, l, tag_hydro, nb.targetid);
        MPI_Send_init(hydro_send_[l][nb.bufid],ssize,MPI_ATHENA_REAL,
                      nb.rank,tag,MPI_COMM_WORLD,&req_hydro_send_[l][nb.bufid]);
        tag=CreateMPITag(pmb->lid, l, tag_hydro, nb.bufid);
        MPI_Recv_init(hydro_recv_[l][nb.bufid],rsize,MPI_ATHENA_REAL,
                      nb.rank,tag,MPI_COMM_WORLD,&req_hydro_recv_[l][nb.bufid]);

        // flux correction
        if(pmb->pmy_mesh->multilevel==true && nb.type==neighbor_face) {
          int fi1, fi2, size;
          if(nb.fid==0 || nb.fid==1)
            fi1=myox2, fi2=myox3, size=((pmb->block_size.nx2+1)/2)*((pmb->block_size.nx3+1)/2);
          else if(nb.fid==2 || nb.fid==3)
            fi1=myox1, fi2=myox3, size=((pmb->block_size.nx1+1)/2)*((pmb->block_size.nx3+1)/2);
          else if(nb.fid==4 || nb.fid==5)
            fi1=myox1, fi2=myox2, size=((pmb->block_size.nx1+1)/2)*((pmb->block_size.nx2+1)/2);
          size*=NHYDRO;
          if(nb.level<mylevel) { // send to coarser
            tag=CreateMPITag(nb.lid, l, tag_flcor, nb.targetid);
            MPI_Send_init(flcor_send_[l][nb.fid],size,MPI_ATHENA_REAL,
                nb.rank,tag,MPI_COMM_WORLD,&req_flcor_send_[l][nb.fid]);
          }
          else if(nb.level>mylevel) { // receive from finer
            tag=CreateMPITag(pmb->lid, l, tag_flcor, nb.bufid);
            MPI_Recv_init(flcor_recv_[l][nb.fid][nb.fi2][nb.fi1],size,MPI_ATHENA_REAL,
                nb.rank,tag,MPI_COMM_WORLD,&req_flcor_recv_[l][nb.fid][nb.fi2][nb.fi1]);
          }
        }

        if (MAGNETIC_FIELDS_ENABLED) {
          int size, csize, fsize;
          int size1=((nb.ox1==0)?(pmb->block_size.nx1+1):NGHOST)
                   *((nb.ox2==0)?(pmb->block_size.nx2):NGHOST)
                   *((nb.ox3==0)?(pmb->block_size.nx3):NGHOST);
          int size2=((nb.ox1==0)?(pmb->block_size.nx1):NGHOST)
                   *((nb.ox2==0)?(pmb->block_size.nx2+f2d):NGHOST)
                   *((nb.ox3==0)?(pmb->block_size.nx3):NGHOST);
          int size3=((nb.ox1==0)?(pmb->block_size.nx1):NGHOST)
                   *((nb.ox2==0)?(pmb->block_size.nx2):NGHOST)
                   *((nb.ox3==0)?(pmb->block_size.nx3+f3d):NGHOST);
          size=size1+size2+size3;
          if(pmb->pmy_mesh->multilevel==true) {
            if(nb.type!=neighbor_face) {
              if(nb.ox1!=0) size1=size1/NGHOST*(NGHOST+1);
              if(nb.ox2!=0) size2=size2/NGHOST*(NGHOST+1);
              if(nb.ox3!=0) size3=size3/NGHOST*(NGHOST+1);
            }
            size=size1+size2+size3;
            int f2c1=((nb.ox1==0)?((pmb->block_size.nx1+1)/2+1):cng)
                    *((nb.ox2==0)?((pmb->block_size.nx2+1)/2):cng)
                    *((nb.ox3==0)?((pmb->block_size.nx3+1)/2):cng);
            int f2c2=((nb.ox1==0)?((pmb->block_size.nx1+1)/2):cng)
                    *((nb.ox2==0)?((pmb->block_size.nx2+1)/2+f2d):cng)
                    *((nb.ox3==0)?((pmb->block_size.nx3+1)/2):cng);
            int f2c3=((nb.ox1==0)?((pmb->block_size.nx1+1)/2):cng)
                    *((nb.ox2==0)?((pmb->block_size.nx2+1)/2):cng)
                    *((nb.ox3==0)?((pmb->block_size.nx3+1)/2+f3d):cng);
            if(nb.type!=neighbor_face) {
              if(nb.ox1!=0) f2c1=f2c1/cng*(cng+1);
              if(nb.ox2!=0) f2c2=f2c2/cng*(cng+1);
              if(nb.ox3!=0) f2c3=f2c3/cng*(cng+1);
            }
            fsize=f2c1+f2c2+f2c3;
            int c2f1=((nb.ox1==0)?((pmb->block_size.nx1+1)/2+1+cng):cng)
                    *((nb.ox2==0)?((pmb->block_size.nx2+1)/2+cng*f2d):cng)
                    *((nb.ox3==0)?((pmb->block_size.nx3+1)/2+cng*f3d):cng);
            int c2f2=((nb.ox1==0)?((pmb->block_size.nx1+1)/2+cng):cng)
                    *((nb.ox2==0)?((pmb->block_size.nx2+1)/2+f2d+cng*f2d):cng)
                    *((nb.ox3==0)?((pmb->block_size.nx3+1)/2+cng*f3d):cng);
            int c2f3=((nb.ox1==0)?((pmb->block_size.nx1+1)/2+cng):cng)
                    *((nb.ox2==0)?((pmb->block_size.nx2+1)/2+f2d*cng):cng)
                    *((nb.ox3==0)?((pmb->block_size.nx3+1)/2+f3d+cng*f3d):cng);
            csize=c2f1+c2f2+c2f3;
          }
          if(nb.level==mylevel) // same
            ssize=size, rsize=size;
          else if(nb.level<mylevel) // coarser
            ssize=fsize, rsize=csize;
          else // finer
            ssize=csize, rsize=fsize;

          tag=CreateMPITag(nb.lid, l, tag_field, nb.targetid);
          MPI_Send_init(field_send_[l][nb.bufid],ssize,MPI_ATHENA_REAL,
                        nb.rank,tag,MPI_COMM_WORLD,&req_field_send_[l][nb.bufid]);
          tag=CreateMPITag(pmb->lid, l, tag_field, nb.bufid);
          MPI_Recv_init(field_recv_[l][nb.bufid],rsize,MPI_ATHENA_REAL,
                        nb.rank,tag,MPI_COMM_WORLD,&req_field_recv_[l][nb.bufid]);
          // EMF correction
          int fi1, fi2, f2csize;
          if(nb.type==neighbor_face) { // face
            if(pmb->block_size.nx3 > 1) { // 3D
              if(nb.fid==inner_x1 || nb.fid==outer_x1) {
                size=(pmb->block_size.nx2+1)*(pmb->block_size.nx3)
                    +(pmb->block_size.nx2)*(pmb->block_size.nx3+1);
                f2csize=(pmb->block_size.nx2/2+1)*(pmb->block_size.nx3/2)
                    +(pmb->block_size.nx2/2)*(pmb->block_size.nx3/2+1);
              }
              else if(nb.fid==inner_x2 || nb.fid==outer_x2) {
                size=(pmb->block_size.nx1+1)*(pmb->block_size.nx3)
                    +(pmb->block_size.nx1)*(pmb->block_size.nx3+1);
                f2csize=(pmb->block_size.nx1/2+1)*(pmb->block_size.nx3/2)
                    +(pmb->block_size.nx1/2)*(pmb->block_size.nx3/2+1);
              }
              else if(nb.fid==inner_x3 || nb.fid==outer_x3) {
                size=(pmb->block_size.nx1+1)*(pmb->block_size.nx2)
                    +(pmb->block_size.nx1)*(pmb->block_size.nx2+1);
                f2csize=(pmb->block_size.nx1/2+1)*(pmb->block_size.nx2/2)
                    +(pmb->block_size.nx1/2)*(pmb->block_size.nx2/2+1);
              }
            }
            else if(pmb->block_size.nx2 > 1) { // 2D
              if(nb.fid==inner_x1 || nb.fid==outer_x1) {
                size=(pmb->block_size.nx2+1)+pmb->block_size.nx2;
                f2csize=(pmb->block_size.nx2/2+1)+pmb->block_size.nx2/2;
              }
              else if(nb.fid==inner_x2 || nb.fid==outer_x2) {
                size=(pmb->block_size.nx1+1)+pmb->block_size.nx1;
                f2csize=(pmb->block_size.nx1/2+1)+pmb->block_size.nx1/2;
              }
            }
            else // 1D
              size=f2csize=2;
          }
          else if(nb.type==neighbor_edge) { // edge
            if(pmb->block_size.nx3 > 1) { // 3D
              if(nb.eid>=0 && nb.eid<4) {
                size=pmb->block_size.nx3;
                f2csize=pmb->block_size.nx3/2;
              }
              else if(nb.eid>=4 && nb.eid<8) {
                size=pmb->block_size.nx2;
                f2csize=pmb->block_size.nx2/2;
              }
              else if(nb.eid>=8 && nb.eid<12) {
                size=pmb->block_size.nx1;
                f2csize=pmb->block_size.nx1/2;
              }
            }
            else if(pmb->block_size.nx2 > 1) // 2D
              size=f2csize=1;
          }
          else // corner
            continue;

          if(nb.level==mylevel) { // the same level
            if((nb.type==neighbor_face) || ((nb.type==neighbor_edge) && (edge_flag_[nb.eid]==true))) {
              tag=CreateMPITag(nb.lid, l, tag_emfcor, nb.targetid);
              MPI_Send_init(emfcor_send_[l][nb.bufid],size,MPI_ATHENA_REAL,
                            nb.rank,tag,MPI_COMM_WORLD,&req_emfcor_send_[l][nb.bufid]);
              tag=CreateMPITag(pmb->lid, l, tag_emfcor, nb.bufid);
              MPI_Recv_init(emfcor_recv_[l][nb.bufid],size,MPI_ATHENA_REAL,
                            nb.rank,tag,MPI_COMM_WORLD,&req_emfcor_recv_[l][nb.bufid]);
            }
          }
          if(nb.level>mylevel) { // finer neighbor
            tag=CreateMPITag(pmb->lid, l, tag_emfcor, nb.bufid);
            MPI_Recv_init(emfcor_recv_[l][nb.bufid],f2csize,MPI_ATHENA_REAL,
                          nb.rank,tag,MPI_COMM_WORLD,&req_emfcor_recv_[l][nb.bufid]);
          }
          if(nb.level<mylevel) { // coarser neighbor
            tag=CreateMPITag(nb.lid, l, tag_emfcor, nb.targetid);
            MPI_Send_init(emfcor_send_[l][nb.bufid],f2csize,MPI_ATHENA_REAL,
                          nb.rank,tag,MPI_COMM_WORLD,&req_emfcor_send_[l][nb.bufid]);
          }
        }
      }
    }
  }
#endif
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::EnrollHydroBoundaryFunction(enum direction dir,
//                                                       BValHydro_t my_bc)
//  \brief Enroll a user-defined boundary function for hydro

void BoundaryValues::EnrollHydroBoundaryFunction(enum direction dir, BValHydro_t my_bc)
{
  std::stringstream msg;
  if(dir<0 || dir>5) {
    msg << "### FATAL ERROR in EnrollHydroBoundaryCondition function" << std::endl
        << "dirName = " << dir << " not valid" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if(pmy_mblock_->block_bcs[dir]==-1) return;
  if(pmy_mblock_->block_bcs[dir]!=3) {
    msg << "### FATAL ERROR in EnrollHydroBoundaryCondition function" << std::endl
        << "A user-defined boundary condition flag (3) must be specified "
        << "in the input file to use a user-defined boundary function." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  HydroBoundary_[dir]=my_bc;
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
  int i;
  MeshBlock *pmb=pmy_mblock_;
  for(int i=0;i<nface_;i++) {
    if(pmb->block_bcs[i]==3) {
      if(HydroBoundary_[i]==NULL) {
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues::CheckBoundary" << std::endl
            << "A user-defined boundary is specified but the hydro boundary function "
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
    NeighborBlock& nb = pmb->neighbor[n];
    if(nb.rank!=Globals::my_rank) { 
      MPI_Start(&req_hydro_recv_[0][nb.bufid]);
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
  int mylevel=pmb->loc.level;
  for(int l=0;l<NSTEP;l++) {
    firsttime_[l]=true;
    for(int n=0;n<pmb->nneighbor;n++) {
      NeighborBlock& nb = pmb->neighbor[n];
      if(nb.rank!=Globals::my_rank) { 
        MPI_Start(&req_hydro_recv_[l][nb.bufid]);
        if(nb.type==neighbor_face && nb.level>mylevel)
          MPI_Start(&req_flcor_recv_[l][nb.fid][nb.fi2][nb.fi1]);
        if (MAGNETIC_FIELDS_ENABLED) {
          MPI_Start(&req_field_recv_[l][nb.bufid]);
          if(nb.type==neighbor_face || nb.type==neighbor_edge) {
            if((nb.level>mylevel) || ((nb.level==mylevel) && ((nb.type==neighbor_face)
            || ((nb.type==neighbor_edge) && (edge_flag_[nb.eid]==true)))))
              MPI_Start(&req_emfcor_recv_[l][nb.bufid]);
          }
        }
      }
    }
  }
#endif
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::RestrictHydro(AthenaArray<Real> &src,
//                           int csi, int cei, int csj, int cej, int csk, int cek)
//  \brief restrict the hydro data and set them into the coarse buffer
void BoundaryValues::RestrictHydro(AthenaArray<Real> &src, 
                             int csi, int cei, int csj, int cej, int csk, int cek)
{
  MeshBlock *pmb=pmy_mblock_;
  Coordinates *pco=pmb->pcoord;
  int si=(csi-pmb->cis)*2+pmb->is, ei=(cei-pmb->cis)*2+pmb->is+1;

  // store the restricted data in the prolongation buffer for later use
  if(pmb->block_size.nx3>1) { // 3D
    for (int n=0; n<(NHYDRO); ++n) {
      for (int ck=csk; ck<=cek; ck++) {
        int k=(ck-pmb->cks)*2+pmb->ks;
        for (int cj=csj; cj<=cej; cj++) {
          int j=(cj-pmb->cjs)*2+pmb->js;
          pco->CellVolume(k,j,si,ei,fvol_[0][0]);
          pco->CellVolume(k,j+1,si,ei,fvol_[0][1]);
          pco->CellVolume(k+1,j,si,ei,fvol_[1][0]);
          pco->CellVolume(k+1,j+1,si,ei,fvol_[1][1]);
          for (int ci=csi; ci<=cei; ci++) {
            int i=(ci-pmb->cis)*2+pmb->is;
            Real tvol=fvol_[0][0](i)+fvol_[0][0](i+1)+fvol_[0][1](i)+fvol_[0][1](i+1)
                     +fvol_[1][0](i)+fvol_[1][0](i+1)+fvol_[1][1](i)+fvol_[1][1](i+1);
            coarse_cons_(n,ck,cj,ci)=
              (src(n,k  ,j  ,i)*fvol_[0][0](i)+src(n,k  ,j  ,i+1)*fvol_[0][0](i+1)
              +src(n,k  ,j+1,i)*fvol_[0][1](i)+src(n,k  ,j+1,i+1)*fvol_[0][1](i+1)
              +src(n,k+1,j  ,i)*fvol_[1][0](i)+src(n,k+1,j  ,i+1)*fvol_[1][0](i+1)
              +src(n,k+1,j+1,i)*fvol_[1][1](i)+src(n,k+1,j+1,i+1)*fvol_[1][1](i+1))/tvol;
          }
        }
      }
    }
  }
  else if(pmb->block_size.nx2>1) { // 2D
    for (int n=0; n<(NHYDRO); ++n) {
      for (int cj=csj; cj<=cej; cj++) {
        int j=(cj-pmb->cjs)*2+pmb->js;
        pco->CellVolume(0,j  ,si,ei,fvol_[0][0]);
        pco->CellVolume(0,j+1,si,ei,fvol_[0][1]);
        for (int ci=csi; ci<=cei; ci++) {
          int i=(ci-pmb->cis)*2+pmb->is;
          Real tvol=fvol_[0][0](i)+fvol_[0][0](i+1)+fvol_[0][1](i)+fvol_[0][1](i+1);
          coarse_cons_(n,0,cj,ci)=
            (src(n,0,j  ,i)*fvol_[0][0](i)+src(n,0,j  ,i+1)*fvol_[0][0](i+1)
            +src(n,0,j+1,i)*fvol_[0][1](i)+src(n,0,j+1,i+1)*fvol_[0][1](i+1))/tvol;
        }
      }
    }
  }
  else { // 1D
    int j=pmb->js, cj=pmb->cjs, k=pmb->ks, ck=pmb->cks; 
    for (int n=0; n<(NHYDRO); ++n) {
      pco->CellVolume(k,j,si,ei,fvol_[0][0]);
      for (int ci=csi; ci<=cei; ci++) {
        int i=(ci-pmb->cis)*2+pmb->is;
        Real tvol=fvol_[0][0](i)+fvol_[0][0](i+1);
        coarse_cons_(n,ck,cj,ci)
          =(src(n,k,j,i)*fvol_[0][0](i)+src(n,k,j,i+1)*fvol_[0][0](i+1))/tvol;
      }
    }
  }
}


//--------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadHydroBoundaryBufferSameLevel(AthenaArray<Real> &src,
//                                                 Real *buf, NeighborBlock& nb)
//  \brief Set hydro boundary buffers for sending to a block on the same level
int BoundaryValues::LoadHydroBoundaryBufferSameLevel(AthenaArray<Real> &src, Real *buf,
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

  int p=0;
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=sk; k<=ek; ++k) {
      for (int j=sj; j<=ej; ++j) {
#pragma simd
        for (int i=si; i<=ei; ++i)
          buf[p++]=src(n,k,j,i);
      }
    }
  }
  return p;
}


//--------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadHydroBoundaryBufferToCoarser(AthenaArray<Real> &src,
//                                                 Real *buf, NeighborBlock& nb)
//  \brief Set hydro boundary buffers for sending to a block on the coarser level
int BoundaryValues::LoadHydroBoundaryBufferToCoarser(AthenaArray<Real> &src, Real *buf,
                                                     NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_mblock_;
  int si, sj, sk, ei, ej, ek;
  int cn=pmb->cnghost-1;

  si=(nb.ox1>0)?(pmb->cie-cn):pmb->cis;
  ei=(nb.ox1<0)?(pmb->cis+cn):pmb->cie;
  sj=(nb.ox2>0)?(pmb->cje-cn):pmb->cjs;
  ej=(nb.ox2<0)?(pmb->cjs+cn):pmb->cje;
  sk=(nb.ox3>0)?(pmb->cke-cn):pmb->cks;
  ek=(nb.ox3<0)?(pmb->cks+cn):pmb->cke;

  // restrict the data before sending
  RestrictHydro(src, si, ei, sj, ej, sk, ek);

  int p=0;
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=sk; k<=ek; k++) {
      for (int j=sj; j<=ej; j++) {
#pragma simd
        for (int i=si; i<=ei; i++)
            buf[p++]=coarse_cons_(n,k,j,i);
      }
    }
  }
  return p;
}


//--------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadHydroBoundaryBufferToFiner(AthenaArray<Real> &src,
//                                                 Real *buf, NeighborBlock& nb)
//  \brief Set hydro boundary buffers for sending to a block on the finer level
int BoundaryValues::LoadHydroBoundaryBufferToFiner(AthenaArray<Real> &src, Real *buf,
                                                   NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_mblock_;
  int si, sj, sk, ei, ej, ek;
  int cn=pmb->cnghost-1;

  si=(nb.ox1>0)?(pmb->ie-cn):pmb->is;
  ei=(nb.ox1<0)?(pmb->is+cn):pmb->ie;
  sj=(nb.ox2>0)?(pmb->je-cn):pmb->js;
  ej=(nb.ox2<0)?(pmb->js+cn):pmb->je;
  sk=(nb.ox3>0)?(pmb->ke-cn):pmb->ks;
  ek=(nb.ox3<0)?(pmb->ks+cn):pmb->ke;

  // send the data first and later prolongate on the target block
  // need to add edges for faces, add corners for edges
  if(nb.ox1==0) {
    if(nb.fi1==1)   si+=pmb->block_size.nx1/2-pmb->cnghost;
    else            ei-=pmb->block_size.nx1/2-pmb->cnghost;
  }
  if(nb.ox2==0 && pmb->block_size.nx2 > 1) {
    if(nb.ox1!=0) {
      if(nb.fi1==1) sj+=pmb->block_size.nx2/2-pmb->cnghost;
      else          ej-=pmb->block_size.nx2/2-pmb->cnghost;
    }
    else {
      if(nb.fi2==1) sj+=pmb->block_size.nx2/2-pmb->cnghost;
      else          ej-=pmb->block_size.nx2/2-pmb->cnghost;
    }
  }
  if(nb.ox3==0 && pmb->block_size.nx3 > 1) {
    if(nb.ox1!=0 && nb.ox2!=0) {
      if(nb.fi1==1) sk+=pmb->block_size.nx3/2-pmb->cnghost;
      else          ek-=pmb->block_size.nx3/2-pmb->cnghost;
    }
    else {
      if(nb.fi2==1) sk+=pmb->block_size.nx3/2-pmb->cnghost;
      else          ek-=pmb->block_size.nx3/2-pmb->cnghost;
    }
  }

  int p=0;
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=sk; k<=ek; ++k) {
      for (int j=sj; j<=ej; ++j) {
#pragma simd
        for (int i=si; i<=ei; ++i)
          buf[p++]=src(n,k,j,i);
      }
    }
  }
  return p;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SendHydroBoundaryBuffers(AthenaArray<Real> &src, int step)
//  \brief Send boundary buffers
void BoundaryValues::SendHydroBoundaryBuffers(AthenaArray<Real> &src, int step)
{
  MeshBlock *pmb=pmy_mblock_;
  int mylevel=pmb->loc.level;

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    int ssize;
    if(nb.level==mylevel)
      ssize=LoadHydroBoundaryBufferSameLevel(src, hydro_send_[step][nb.bufid],nb);
    else if(nb.level<mylevel)
      ssize=LoadHydroBoundaryBufferToCoarser(src, hydro_send_[step][nb.bufid],nb);
    else
      ssize=LoadHydroBoundaryBufferToFiner(src, hydro_send_[step][nb.bufid], nb);
    if(nb.rank == Globals::my_rank) { // on the same process
      MeshBlock *pbl=pmb->pmy_mesh->FindMeshBlock(nb.gid);
      std::memcpy(pbl->pbval->hydro_recv_[step][nb.targetid],
                  hydro_send_[step][nb.bufid], ssize*sizeof(Real));
      pbl->pbval->hydro_flag_[step][nb.targetid]=boundary_arrived;
    }
#ifdef MPI_PARALLEL
    else // MPI
      MPI_Start(&req_hydro_send_[step][nb.bufid]);
#endif
  }

  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetHydroBoundarySameLevel(AthenaArray<Real> &dst,
//                                                     Real *buf, NeighborBlock& nb)
//  \brief Set hydro boundary received from a block on the same level
void BoundaryValues::SetHydroBoundarySameLevel(AthenaArray<Real> &dst, Real *buf,
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

  int p=0;
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=sk; k<=ek; ++k) {
      for (int j=sj; j<=ej; ++j) {
#pragma simd
        for (int i=si; i<=ei; ++i)
          dst(n,k,j,i) = buf[p++];
      }
    }
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetHydroBoundaryFromCoarser(Real *buf, NeighborBlock& nb)
//  \brief Set hydro prolongation buffer received from a block on the same level
void BoundaryValues::SetHydroBoundaryFromCoarser(Real *buf, NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_mblock_;

  int si, sj, sk, ei, ej, ek;
  int cng=pmb->cnghost;

  if(nb.ox1==0) {
    si=pmb->cis, ei=pmb->cie;
    if((pmb->loc.lx1&1L)==0L) ei+=cng;
    else             si-=cng; 
  }
  else if(nb.ox1>0)  si=pmb->cie+1,   ei=pmb->cie+cng;
  else               si=pmb->cis-cng, ei=pmb->cis-1;
  if(nb.ox2==0) {
    sj=pmb->cjs, ej=pmb->cje;
    if(pmb->block_size.nx2 > 1) {
      if((pmb->loc.lx2&1L)==0L) ej+=cng;
      else             sj-=cng; 
    }
  }
  else if(nb.ox2>0)  sj=pmb->cje+1,   ej=pmb->cje+cng;
  else               sj=pmb->cjs-cng, ej=pmb->cjs-1;
  if(nb.ox3==0) {
    sk=pmb->cks, ek=pmb->cke;
    if(pmb->block_size.nx3 > 1) {
      if((pmb->loc.lx3&1L)==0L) ek+=cng;
      else             sk-=cng; 
    }
  }
  else if(nb.ox3>0)  sk=pmb->cke+1,   ek=pmb->cke+cng;
  else               sk=pmb->cks-cng, ek=pmb->cks-1;

  int p=0;
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=sk; k<=ek; ++k) {
      for (int j=sj; j<=ej; ++j) {
#pragma simd
        for (int i=si; i<=ei; ++i)
          coarse_cons_(n,k,j,i) = buf[p++];
      }
    }
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetHydroBoundaryFromFiner(AthenaArray<Real> &dst,
//                                                     Real *buf, NeighborBlock& nb)
//  \brief Set hydro boundary received from a block on the same level
void BoundaryValues::SetHydroBoundaryFromFiner(AthenaArray<Real> &dst, Real *buf,
                                               NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_mblock_;
  // receive already restricted data
  int si, sj, sk, ei, ej, ek;

  if(nb.ox1==0) {
    si=pmb->is, ei=pmb->ie;
    if(nb.fi1==1)   si+=pmb->block_size.nx1/2;
    else            ei-=pmb->block_size.nx1/2;
  }
  else if(nb.ox1>0) si=pmb->ie+1,      ei=pmb->ie+NGHOST;
  else              si=pmb->is-NGHOST, ei=pmb->is-1;
  if(nb.ox2==0) {
    sj=pmb->js, ej=pmb->je;
    if(pmb->block_size.nx2 > 1) {
      if(nb.ox1!=0) {
        if(nb.fi1==1) sj+=pmb->block_size.nx2/2;
        else          ej-=pmb->block_size.nx2/2;
      }
      else {
        if(nb.fi2==1) sj+=pmb->block_size.nx2/2;
        else          ej-=pmb->block_size.nx2/2;
      }
    }
  }
  else if(nb.ox2>0) sj=pmb->je+1,      ej=pmb->je+NGHOST;
  else              sj=pmb->js-NGHOST, ej=pmb->js-1;
  if(nb.ox3==0) {
    sk=pmb->ks, ek=pmb->ke;
    if(pmb->block_size.nx3 > 1) {
      if(nb.ox1!=0 && nb.ox2!=0) {
        if(nb.fi1==1) sk+=pmb->block_size.nx3/2;
        else          ek-=pmb->block_size.nx3/2;
      }
      else {
        if(nb.fi2==1) sk+=pmb->block_size.nx3/2;
        else          ek-=pmb->block_size.nx3/2;
      }
    }
  }
  else if(nb.ox3>0) sk=pmb->ke+1,      ek=pmb->ke+NGHOST;
  else              sk=pmb->ks-NGHOST, ek=pmb->ks-1;

  int p=0;
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=sk; k<=ek; ++k) {
      for (int j=sj; j<=ej; ++j) {
#pragma simd
        for (int i=si; i<=ei; ++i)
          dst(n,k,j,i) = buf[p++];
      }
    }
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn bool BoundaryValues::ReceiveHydroBoundaryBuffers(AthenaArray<Real> &dst, int step)
//  \brief receive the boundary data
bool BoundaryValues::ReceiveHydroBoundaryBuffers(AthenaArray<Real> &dst, int step)
{
  MeshBlock *pmb=pmy_mblock_;
  bool flag=true;

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    if(hydro_flag_[step][nb.bufid]==boundary_completed) continue;
    if(hydro_flag_[step][nb.bufid]==boundary_waiting) {
      if(nb.rank==Globals::my_rank) {// on the same process
        flag=false;
        continue;
      }
#ifdef MPI_PARALLEL
      else { // MPI boundary
        int test;
        MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&test,MPI_STATUS_IGNORE);
        MPI_Test(&req_hydro_recv_[step][nb.bufid],&test,MPI_STATUS_IGNORE);
        if(test==false) {
          flag=false;
          continue;
        }
        hydro_flag_[step][nb.bufid] = boundary_arrived;
      }
#endif
    }
    if(nb.level==pmb->loc.level)
      SetHydroBoundarySameLevel(dst, hydro_recv_[step][nb.bufid], nb);
    else if(nb.level<pmb->loc.level) // this set only the prolongation buffer
      SetHydroBoundaryFromCoarser(hydro_recv_[step][nb.bufid], nb);
    else
      SetHydroBoundaryFromFiner(dst, hydro_recv_[step][nb.bufid], nb);
    hydro_flag_[step][nb.bufid] = boundary_completed; // completed
  }
  return flag;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ReceiveHydroBoundaryBuffersWithWait(AthenaArray<Real> &dst,
//                                                               int step)
//  \brief receive the boundary data for initialization
void BoundaryValues::ReceiveHydroBoundaryBuffersWithWait(AthenaArray<Real> &dst, int step)
{
  MeshBlock *pmb=pmy_mblock_;

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb = pmb->neighbor[n];
#ifdef MPI_PARALLEL
    if(nb.rank!=Globals::my_rank)
      MPI_Wait(&req_hydro_recv_[0][nb.bufid],MPI_STATUS_IGNORE);
#endif
    if(nb.level==pmb->loc.level)
      SetHydroBoundarySameLevel(dst, hydro_recv_[0][nb.bufid], nb);
    else if(nb.level<pmb->loc.level)
      SetHydroBoundaryFromCoarser(hydro_recv_[0][nb.bufid], nb);
    else
      SetHydroBoundaryFromFiner(dst, hydro_recv_[0][nb.bufid], nb);
    hydro_flag_[0][nb.bufid] = boundary_completed; // completed
  }
  return;
}


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
    if(nb.type!=neighbor_face) break;
    if(nb.level==pmb->loc.level-1) {
      int p=0;
      // x1 direction
      if(nb.fid==inner_x1 || nb.fid==outer_x1) {
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
      else if(nb.fid==inner_x2 || nb.fid==outer_x2) {
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
      else if(nb.fid==inner_x3 || nb.fid==outer_x3) {
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
        pbl->pbval->flcor_flag_[step][(nb.fid^1)][fi2][fi1]=boundary_arrived;
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
//! \fn bool BoundaryValues::ReceiveFluxCorrection(AthenaArray<Real> &dst, int step)
//  \brief Receive and apply the surace flux from the finer neighbor(s)
bool BoundaryValues::ReceiveFluxCorrection(AthenaArray<Real> &dst, int step)
{
  MeshBlock *pmb=pmy_mblock_;
  Coordinates *pco=pmb->pcoord;
  AthenaArray<Real> &x1flux=pmb->phydro->flux[x1face];
  AthenaArray<Real> &x2flux=pmb->phydro->flux[x2face];
  AthenaArray<Real> &x3flux=pmb->phydro->flux[x3face];
  bool flag=true;

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    if(nb.type!=neighbor_face) break;
    if(nb.level==pmb->loc.level+1) {
      if(flcor_flag_[step][nb.fid][nb.fi2][nb.fi1]==boundary_completed) continue;
      if(flcor_flag_[step][nb.fid][nb.fi2][nb.fi1]==boundary_waiting) {
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
          flcor_flag_[step][nb.fid][nb.fi2][nb.fi1] = boundary_arrived;
        }
#endif
      }
      // boundary arrived; apply flux correction
      Real *buf=flcor_recv_[step][nb.fid][nb.fi2][nb.fi1];
      int p=0;
      if(nb.fid==inner_x1 || nb.fid==outer_x1) {
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
      else if(nb.fid==inner_x2 || nb.fid==outer_x2) {
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
      else if(nb.fid==inner_x3 || nb.fid==outer_x3) {
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

      flcor_flag_[step][nb.fid][nb.fi2][nb.fi1] = boundary_completed;
    }
  }

  return flag;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ProlongateHydroBoundaries(AthenaArray<Real> &dst)
//  \brief Prolongate the hydro in the ghost zones from the prolongation buffer
void BoundaryValues::ProlongateHydroBoundaries(AthenaArray<Real> &dst)
{
  MeshBlock *pmb=pmy_mblock_;
  Coordinates *pco=pmb->pcoord;
  Coordinates *pcrs=pmb->pcoarsec;
  int mox1, mox2, mox3;
  long int &lx1=pmb->loc.lx1;
  long int &lx2=pmb->loc.lx2;
  long int &lx3=pmb->loc.lx3;
  int &mylevel=pmb->loc.level;
  mox1=((int)(lx1&1L)<<1)-1;
  mox2=((int)(lx2&1L)<<1)-1;
  mox3=((int)(lx3&1L)<<1)-1;
  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    if(nb.level >= mylevel) continue;
    int mytype=std::abs(nb.ox1)+std::abs(nb.ox2)+std::abs(nb.ox3);
    // fill the required ghost-ghost zone
    int nis, nie, njs, nje, nks, nke;
    nis=std::max(nb.ox1-1,-1), nie=std::min(nb.ox1+1,1);
    if(pmb->block_size.nx2==1) njs=0, nje=0;
    else njs=std::max(nb.ox2-1,-1), nje=std::min(nb.ox2+1,1);
    if(pmb->block_size.nx3==1) nks=0, nke=0;
    else nks=std::max(nb.ox3-1,-1), nke=std::min(nb.ox3+1,1);
    for(int nk=nks; nk<=nke; nk++) {
      for(int nj=njs; nj<=nje; nj++) {
        for(int ni=nis; ni<=nie; ni++) {
          int ntype=std::abs(ni)+std::abs(nj)+std::abs(nk);
          if(ntype==0) continue; // skip myself
          if(pmb->nblevel[nk+1][nj+1][ni+1]!=mylevel
          && pmb->nblevel[nk+1][nj+1][ni+1]!=-1)
            continue; // physical boundary will also be restricted
          if(ntype>mytype) {
            if(pmb->block_size.nx3 > 1) // 3D
              if(((mox1==ni)+(mox2==nj)+(mox3==nk)) != ntype) continue;
            else if(pmb->block_size.nx2 > 1) // 2D
              if(((mox1==ni)+(mox2==nj)) != ntype) continue;
          }

          // this neighbor block is on the same level
          // and needs to be restricted for prolongation
          int ris, rie, rjs, rje, rks, rke;
          if(ni==0) {
            ris=pmb->cis, rie=pmb->cie;
            if(nb.ox1==1) ris=pmb->cie;
            else if(nb.ox1==-1) rie=pmb->cis;
          }
          else if(ni== 1) ris=pmb->cie+1, rie=pmb->cie+1;
          else if(ni==-1) ris=pmb->cis-1, rie=pmb->cis-1;
          if(nj==0) {
            rjs=pmb->cjs, rje=pmb->cje;
            if(nb.ox2==1) rjs=pmb->cje;
            else if(nb.ox2==-1) rje=pmb->cjs;
          }
          else if(nj== 1) rjs=pmb->cje+1, rje=pmb->cje+1;
          else if(nj==-1) rjs=pmb->cjs-1, rje=pmb->cjs-1;
          if(nk==0) {
            rks=pmb->cks, rke=pmb->cke;
            if(nb.ox3==1) rks=pmb->cke;
            else if(nb.ox3==-1) rke=pmb->cks;
          }
          else if(nk== 1) rks=pmb->cke+1, rke=pmb->cke+1;
          else if(nk==-1) rks=pmb->cks-1, rke=pmb->cks-1;
          RestrictHydro(dst, ris, rie, rjs, rje, rks, rke);
        }
      }
    }

    // now that the ghost-ghost zones are filled
    // calculate the slope with a limiter and interpolate the data
    int cn = (NGHOST+1)/2;
    int si, ei, sj, ej, sk, ek;
    if(nb.ox1==0) {
      si=pmb->cis, ei=pmb->cie;
      if((lx1&1L)==0L) ei++;
      else             si--;
    }
    else if(nb.ox1>0) si=pmb->cie+1,  ei=pmb->cie+cn;
    else              si=pmb->cis-cn, ei=pmb->cis-1;
    if(nb.ox2==0) {
      sj=pmb->cjs, ej=pmb->cje;
      if(pmb->block_size.nx2 > 1) {
        if((lx2&1L)==0L) ej++;
        else             sj--;
      }
    }
    else if(nb.ox2>0) sj=pmb->cje+1,  ej=pmb->cje+cn;
    else              sj=pmb->cjs-cn, ej=pmb->cjs-1;
    if(nb.ox3==0) {
      sk=pmb->cks, ek=pmb->cke;
      if(pmb->block_size.nx3 > 1) {
        if((lx3&1L)==0L) ek++;
        else             sk--;
      }
    }
    else if(nb.ox3>0) sk=pmb->cke+1,  ek=pmb->cke+cn;
    else              sk=pmb->cks-cn, ek=pmb->cks-1;

    if(pmb->block_size.nx3 > 1) { // 3D
      for(int n=0; n<NHYDRO; n++) {
        for(int k=sk; k<=ek; k++) {
          int fk=(k-pmb->cks)*2+pmb->ks;
          Real& x3m = pcrs->x3v(k-1);
          Real& x3c = pcrs->x3v(k);
          Real& x3p = pcrs->x3v(k+1);
          Real dx3m = x3c - x3m;
          Real dx3p = x3p - x3c;
          Real& fx3m = pco->x3v(fk);
          Real& fx3p = pco->x3v(fk+1);
          Real dx3fm= x3c-fx3m;
          Real dx3fp= fx3p-x3c;
          for(int j=sj; j<=ej; j++) {
            int fj=(j-pmb->cjs)*2+pmb->js;
            Real& x2m = pcrs->x2v(j-1);
            Real& x2c = pcrs->x2v(j);
            Real& x2p = pcrs->x2v(j+1);
            Real dx2m = x2c - x2m;
            Real dx2p = x2p - x2c;
            Real& fx2m = pco->x2v(fj);
            Real& fx2p = pco->x2v(fj+1);
            Real dx2fm= x2c-fx2m;
            Real dx2fp= fx2p-x2c;
            for(int i=si; i<=ei; i++) {
              int fi=(i-pmb->cis)*2+pmb->is;
              Real& x1m = pcrs->x1v(i-1);
              Real& x1c = pcrs->x1v(i);
              Real& x1p = pcrs->x1v(i+1);
              Real dx1m = x1c - x1m;
              Real dx1p = x1p - x1c;
              Real& fx1m = pco->x1v(fi);
              Real& fx1p = pco->x1v(fi+1);
              Real dx1fm= x1c-fx1m;
              Real dx1fp= fx1p-x1c;
              Real ccval=coarse_cons_(n,k,j,i);

              // calculate 3D gradients using the minmod limiter
              Real gx1m = (ccval-coarse_cons_(n,k,j,i-1))/dx1m;
              Real gx1p = (coarse_cons_(n,k,j,i+1)-ccval)/dx1p;
              Real gx1c = 0.5*(SIGN(gx1m)+SIGN(gx1p))*std::min(std::abs(gx1m),std::abs(gx1p));
              Real gx2m = (ccval-coarse_cons_(n,k,j-1,i))/dx2m;
              Real gx2p = (coarse_cons_(n,k,j+1,i)-ccval)/dx2p;
              Real gx2c = 0.5*(SIGN(gx2m)+SIGN(gx2p))*std::min(std::abs(gx2m),std::abs(gx2p));
              Real gx3m = (ccval-coarse_cons_(n,k-1,j,i))/dx3m;
              Real gx3p = (coarse_cons_(n,k+1,j,i)-ccval)/dx3p;
              Real gx3c = 0.5*(SIGN(gx3m)+SIGN(gx3p))*std::min(std::abs(gx3m),std::abs(gx3p));

              // interpolate onto the finer grid
              dst(n,fk  ,fj  ,fi  )=ccval-gx1c*dx1fm-gx2c*dx2fm-gx3c*dx3fm;
              dst(n,fk  ,fj  ,fi+1)=ccval+gx1c*dx1fp-gx2c*dx2fm-gx3c*dx3fm;
              dst(n,fk  ,fj+1,fi  )=ccval-gx1c*dx1fm+gx2c*dx2fp-gx3c*dx3fm;
              dst(n,fk  ,fj+1,fi+1)=ccval+gx1c*dx1fp+gx2c*dx2fp-gx3c*dx3fm;
              dst(n,fk+1,fj  ,fi  )=ccval-gx1c*dx1fm-gx2c*dx2fm+gx3c*dx3fp;
              dst(n,fk+1,fj  ,fi+1)=ccval+gx1c*dx1fp-gx2c*dx2fm+gx3c*dx3fp;
              dst(n,fk+1,fj+1,fi  )=ccval-gx1c*dx1fm+gx2c*dx2fp+gx3c*dx3fp;
              dst(n,fk+1,fj+1,fi+1)=ccval+gx1c*dx1fp+gx2c*dx2fp+gx3c*dx3fp;
            }
          }
        }
      }
    }
    else if(pmb->block_size.nx2 > 1) { // 2D
      int k=pmb->cks, fk=pmb->ks;
      for(int n=0; n<NHYDRO; n++) {
        for(int j=sj; j<=ej; j++) {
          int fj=(j-pmb->cjs)*2+pmb->js;
          Real& x2m = pcrs->x2v(j-1);
          Real& x2c = pcrs->x2v(j);
          Real& x2p = pcrs->x2v(j+1);
          Real dx2m = x2c - x2m;
          Real dx2p = x2p - x2c;
          Real& fx2m = pco->x2v(fj);
          Real& fx2p = pco->x2v(fj+1);
          Real dx2fm= x2c-fx2m;
          Real dx2fp= fx2p-x2c;
          for(int i=si; i<=ei; i++) {
            int fi=(i-pmb->cis)*2+pmb->is;
            Real& x1m = pcrs->x1v(i-1);
            Real& x1c = pcrs->x1v(i);
            Real& x1p = pcrs->x1v(i+1);
            Real dx1m = x1c - x1m;
            Real dx1p = x1p - x1c;
            Real& fx1m = pco->x1v(fi);
            Real& fx1p = pco->x1v(fi+1);
            Real dx1fm= x1c-fx1m;
            Real dx1fp= fx1p-x1c;
            Real ccval=coarse_cons_(n,k,j,i);

            // calculate 2D gradients using the minmod limiter
            Real gx1m = (ccval-coarse_cons_(n,k,j,i-1))/dx1m;
            Real gx1p = (coarse_cons_(n,k,j,i+1)-ccval)/dx1p;
            Real gx1c = 0.5*(SIGN(gx1m)+SIGN(gx1p))*std::min(std::abs(gx1m),std::abs(gx1p));
            Real gx2m = (ccval-coarse_cons_(n,k,j-1,i))/dx2m;
            Real gx2p = (coarse_cons_(n,k,j+1,i)-ccval)/dx2p;
            Real gx2c = 0.5*(SIGN(gx2m)+SIGN(gx2p))*std::min(std::abs(gx2m),std::abs(gx2p));

            // interpolate on to the finer grid
            dst(n,fk  ,fj  ,fi  )=ccval-gx1c*dx1fm-gx2c*dx2fm;
            dst(n,fk  ,fj  ,fi+1)=ccval+gx1c*dx1fp-gx2c*dx2fm;
            dst(n,fk  ,fj+1,fi  )=ccval-gx1c*dx1fm+gx2c*dx2fp;
            dst(n,fk  ,fj+1,fi+1)=ccval+gx1c*dx1fp+gx2c*dx2fp;
          }
        }
      }
    }
    else { // 1D
      int k=pmb->cks, fk=pmb->ks, j=pmb->cjs, fj=pmb->js;
      for(int n=0; n<NHYDRO; n++) {
        for(int i=si; i<=ei; i++) {
          int fi=(i-pmb->cis)*2+pmb->is;
          Real& x1m = pcrs->x1v(i-1);
          Real& x1c = pcrs->x1v(i);
          Real& x1p = pcrs->x1v(i+1);
          Real dx1m = x1c - x1m;
          Real dx1p = x1p - x1c;
          Real& fx1m = pco->x1v(fi);
          Real& fx1p = pco->x1v(fi+1);
          Real dx1fm= x1c-fx1m;
          Real dx1fp= fx1p-x1c;
          Real ccval=coarse_cons_(n,k,j,i);

          // calculate 1D gradient using the min-mod limiter
          Real gx1m = (ccval-coarse_cons_(n,k,j,i-1))/dx1m;
          Real gx1p = (coarse_cons_(n,k,j,i+1)-ccval)/dx1p;
          Real gx1c = 0.5*(SIGN(gx1m)+SIGN(gx1p))*std::min(std::abs(gx1m),std::abs(gx1p));

          // interpolate on to the finer grid
          dst(n,fk  ,fj  ,fi  )=ccval-gx1c*dx1fm;
          dst(n,fk  ,fj  ,fi+1)=ccval+gx1c*dx1fp;
        }
      }
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::RestrictFieldX1(AthenaArray<Real> &bx1f,
//                           int csi, int cei, int csj, int cej, int csk, int cek)
//  \brief restrict the x1 field data and set them into the coarse buffer
void BoundaryValues::RestrictFieldX1(AthenaArray<Real> &bx1f, 
                             int csi, int cei, int csj, int cej, int csk, int cek)
{
  MeshBlock *pmb=pmy_mblock_;
  Coordinates *pco=pmb->pcoord;
  int si=(csi-pmb->cis)*2+pmb->is, ei=(cei-pmb->cis)*2+pmb->is;

  // store the restricted data in the prolongation buffer for later use
  if(pmb->block_size.nx3>1) { // 3D
    for (int ck=csk; ck<=cek; ck++) {
      int k=(ck-pmb->cks)*2+pmb->ks;
      for (int cj=csj; cj<=cej; cj++) {
        int j=(cj-pmb->cjs)*2+pmb->js;
        // reuse fvol_ arrays as surface area
        pco->Face1Area(k,   j,   si, ei, fvol_[0][0]);
        pco->Face1Area(k,   j+1, si, ei, fvol_[0][1]);
        pco->Face1Area(k+1, j,   si, ei, fvol_[1][0]);
        pco->Face1Area(k+1, j+1, si, ei, fvol_[1][1]);
        for (int ci=csi; ci<=cei; ci++) {
          int i=(ci-pmb->cis)*2+pmb->is;
          Real tarea=fvol_[0][0](i)+fvol_[0][1](i)+fvol_[1][0](i)+fvol_[1][1](i);
          coarse_b_.x1f(ck,cj,ci)=
            (bx1f(k  ,j,i)*fvol_[0][0](i)+bx1f(k  ,j+1,i)*fvol_[0][1](i)
            +bx1f(k+1,j,i)*fvol_[1][0](i)+bx1f(k+1,j+1,i)*fvol_[1][1](i))/tarea;
        }
      }
    }
  }
  else if(pmb->block_size.nx2>1) { // 2D
    int k=pmb->ks;
    for (int cj=csj; cj<=cej; cj++) {
      int j=(cj-pmb->cjs)*2+pmb->js;
      // reuse fvol_ arrays as surface area
      pco->Face1Area(k,  j,   si, ei, fvol_[0][0]);
      pco->Face1Area(k,  j+1, si, ei, fvol_[0][1]);
      for (int ci=csi; ci<=cei; ci++) {
        int i=(ci-pmb->cis)*2+pmb->is;
        Real tarea=fvol_[0][0](i)+fvol_[0][1](i);
        coarse_b_.x1f(csk,cj,ci)=
          (bx1f(k,j,i)*fvol_[0][0](i)+bx1f(k,j+1,i)*fvol_[0][1](i))/tarea;
      }
    }
  }
  else { // 1D - no restriction, just copy 
    for (int ci=csi; ci<=cei; ci++) {
      int i=(ci-pmb->cis)*2+pmb->is;
      coarse_b_.x1f(csk,csj,ci)=bx1f(pmb->ks,pmb->js,i);
    }
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::RestrictFieldX2(AthenaArray<Real> &bx2f,
//                           int csi, int cei, int csj, int cej, int csk, int cek)
//  \brief restrict the x2 field data and set them into the coarse buffer
void BoundaryValues::RestrictFieldX2(AthenaArray<Real> &bx2f, 
                             int csi, int cei, int csj, int cej, int csk, int cek)
{
  MeshBlock *pmb=pmy_mblock_;
  Coordinates *pco=pmb->pcoord;
  int si=(csi-pmb->cis)*2+pmb->is, ei=(cei-pmb->cis)*2+pmb->is+1;

  // store the restricted data in the prolongation buffer for later use
  if(pmb->block_size.nx3>1) { // 3D
    for (int ck=csk; ck<=cek; ck++) {
      int k=(ck-pmb->cks)*2+pmb->ks;
      for (int cj=csj; cj<=cej; cj++) {
        int j=(cj-pmb->cjs)*2+pmb->js;
        // reuse fvol_ arrays as surface area
        pco->Face2Area(k,   j,  si, ei, fvol_[0][0]);
        pco->Face2Area(k+1, j,  si, ei, fvol_[0][1]);
        for (int ci=csi; ci<=cei; ci++) {
          int i=(ci-pmb->cis)*2+pmb->is;
          Real tarea=fvol_[0][0](i)+fvol_[0][0](i+1)+fvol_[0][1](i)+fvol_[0][1](i+1);
          coarse_b_.x2f(ck,cj,ci)=
            (bx2f(k  ,j,i)*fvol_[0][0](i)+bx2f(k  ,j,i+1)*fvol_[0][0](i+1)
            +bx2f(k+1,j,i)*fvol_[0][1](i)+bx2f(k+1,j,i+1)*fvol_[0][1](i+1))/tarea;
        }
      }
    }
  }
  else if(pmb->block_size.nx2>1) { // 2D
    int k=pmb->ks;
    for (int cj=csj; cj<=cej; cj++) {
      int j=(cj-pmb->cjs)*2+pmb->js;
      // reuse fvol_ arrays as surface area
      pco->Face2Area(k, j, si, ei, fvol_[0][0]);
      for (int ci=csi; ci<=cei; ci++) {
        int i=(ci-pmb->cis)*2+pmb->is;
        Real tarea=fvol_[0][0](i)+fvol_[0][0](i+1);
        coarse_b_.x2f(pmb->cks,cj,ci)=
          (bx2f(k,j,i)*fvol_[0][0](i)+bx2f(k,j,i+1)*fvol_[0][0](i+1))/tarea;
      }
    }
  }
  else { // 1D 
    int k=pmb->ks, j=pmb->js;
    pco->Face2Area(k, j, si, ei, fvol_[0][0]);
    for (int ci=csi; ci<=cei; ci++) {
      int i=(ci-pmb->cis)*2+pmb->is;
        Real tarea=fvol_[0][0](i)+fvol_[0][0](i+1);
        coarse_b_.x2f(pmb->cks,pmb->cjs,ci)=
          (bx2f(k,j,i)*fvol_[0][0](i)+bx2f(k,j,i+1)*fvol_[0][0](i+1))/tarea;
    }
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::RestrictFieldX3(AthenaArray<Real> &bx3f,
//                           int csi, int cei, int csj, int cej, int csk, int cek)
//  \brief restrict the x3 field data and set them into the coarse buffer
void BoundaryValues::RestrictFieldX3(AthenaArray<Real> &bx3f, 
                             int csi, int cei, int csj, int cej, int csk, int cek)
{
  MeshBlock *pmb=pmy_mblock_;
  Coordinates *pco=pmb->pcoord;
  int si=(csi-pmb->cis)*2+pmb->is, ei=(cei-pmb->cis)*2+pmb->is+1;

  // store the restricted data in the prolongation buffer for later use
  if(pmb->block_size.nx3>1) { // 3D
    for (int ck=csk; ck<=cek; ck++) {
      int k=(ck-pmb->cks)*2+pmb->ks;
      for (int cj=csj; cj<=cej; cj++) {
        int j=(cj-pmb->cjs)*2+pmb->js;
        // reuse fvol_ arrays as surface area
        pco->Face3Area(k,   j,  si, ei, fvol_[0][0]);
        pco->Face3Area(k, j+1,  si, ei, fvol_[0][1]);
        for (int ci=csi; ci<=cei; ci++) {
          int i=(ci-pmb->cis)*2+pmb->is;
          Real tarea=fvol_[0][0](i)+fvol_[0][0](i+1)+fvol_[0][1](i)+fvol_[0][1](i+1);
          coarse_b_.x3f(ck,cj,ci)=
            (bx3f(k,j  ,i)*fvol_[0][0](i)+bx3f(k,j  ,i+1)*fvol_[0][0](i+1)
            +bx3f(k,j+1,i)*fvol_[0][1](i)+bx3f(k,j+1,i+1)*fvol_[0][1](i+1))/tarea;
        }
      }
    }
  }
  else if(pmb->block_size.nx2>1) { // 2D
    int k=pmb->ks;
    for (int cj=csj; cj<=cej; cj++) {
      int j=(cj-pmb->cjs)*2+pmb->js;
      // reuse fvol_ arrays as surface area
      pco->Face3Area(k,   j, si, ei, fvol_[0][0]);
      pco->Face3Area(k, j+1, si, ei, fvol_[0][1]);
      for (int ci=csi; ci<=cei; ci++) {
        int i=(ci-pmb->cis)*2+pmb->is;
        Real tarea=fvol_[0][0](i)+fvol_[0][0](i+1)+fvol_[0][1](i)+fvol_[0][1](i+1);
        coarse_b_.x3f(pmb->cks,cj,ci)=
            (bx3f(k,j  ,i)*fvol_[0][0](i)+bx3f(k,j  ,i+1)*fvol_[0][0](i+1)
            +bx3f(k,j+1,i)*fvol_[0][1](i)+bx3f(k,j+1,i+1)*fvol_[0][1](i+1))/tarea;
      }
    }
  }
  else { // 1D 
    int k=pmb->ks, j=pmb->js;
    pco->Face3Area(k, j, si, ei, fvol_[0][0]);
    for (int ci=csi; ci<=cei; ci++) {
      int i=(ci-pmb->cis)*2+pmb->is;
        Real tarea=fvol_[0][0](i)+fvol_[0][0](i+1);
        coarse_b_.x3f(pmb->cks,pmb->cjs,ci)=
          (bx3f(k,j,i)*fvol_[0][0](i)+bx3f(k,j,i+1)*fvol_[0][0](i+1))/tarea;
    }
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
  int p=0;

  // bx1
  if(nb.ox1==0)     si=pmb->is,          ei=pmb->ie+1;
  else if(nb.ox1>0) si=pmb->ie-NGHOST+1, ei=pmb->ie;
  else              si=pmb->is+1,        ei=pmb->is+NGHOST;
  if(nb.ox2==0)     sj=pmb->js,          ej=pmb->je;
  else if(nb.ox2>0) sj=pmb->je-NGHOST+1, ej=pmb->je;
  else              sj=pmb->js,          ej=pmb->js+NGHOST-1;
  if(nb.ox3==0)     sk=pmb->ks,          ek=pmb->ke;
  else if(nb.ox3>0) sk=pmb->ke-NGHOST+1, ek=pmb->ke;
  else              sk=pmb->ks,          ek=pmb->ks+NGHOST-1;
  // for SMR/AMR, always include the overlapping faces in edge and corner boundaries
  if(pmb->pmy_mesh->multilevel==true && nb.type != neighbor_face) {
    if(nb.ox1>0) ei++;
    else if(nb.ox1<0) si--;
  }
  for (int k=sk; k<=ek; ++k) {
    for (int j=sj; j<=ej; ++j) {
#pragma simd
      for (int i=si; i<=ei; ++i)
        buf[p++]=src.x1f(k,j,i);
    }
  }

  // bx2
  if(nb.ox1==0)      si=pmb->is,          ei=pmb->ie;
  else if(nb.ox1>0)  si=pmb->ie-NGHOST+1, ei=pmb->ie;
  else               si=pmb->is,          ei=pmb->is+NGHOST-1;
  if(pmb->block_size.nx2==1) sj=pmb->js,  ej=pmb->je;
  else if(nb.ox2==0) sj=pmb->js,          ej=pmb->je+1;
  else if(nb.ox2>0)  sj=pmb->je-NGHOST+1, ej=pmb->je;
  else               sj=pmb->js+1,        ej=pmb->js+NGHOST;
  if(pmb->pmy_mesh->multilevel==true && nb.type != neighbor_face) {
    if(nb.ox2>0) ej++;
    else if(nb.ox2<0) sj--;
  }
  for (int k=sk; k<=ek; ++k) {
    for (int j=sj; j<=ej; ++j) {
#pragma simd
      for (int i=si; i<=ei; ++i)
        buf[p++]=src.x2f(k,j,i);
    }
  }

  // bx3
  if(nb.ox2==0)      sj=pmb->js,          ej=pmb->je;
  else if(nb.ox2>0)  sj=pmb->je-NGHOST+1, ej=pmb->je;
  else               sj=pmb->js,          ej=pmb->js+NGHOST-1;
  if(pmb->block_size.nx3==1) sk=pmb->ks,  ek=pmb->ke;
  else if(nb.ox3==0) sk=pmb->ks,          ek=pmb->ke+1;
  else if(nb.ox3>0)  sk=pmb->ke-NGHOST+1, ek=pmb->ke;
  else               sk=pmb->ks+1,        ek=pmb->ks+NGHOST;
  if(pmb->pmy_mesh->multilevel==true && nb.type != neighbor_face) {
    if(nb.ox3>0) ek++;
    else if(nb.ox3<0) sk--;
  }
  for (int k=sk; k<=ek; ++k) {
    for (int j=sj; j<=ej; ++j) {
#pragma simd
      for (int i=si; i<=ei; ++i)
        buf[p++]=src.x3f(k,j,i);
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
  MeshBlock *pmb=pmy_mblock_;
  int si, sj, sk, ei, ej, ek;
  int cng=pmb->cnghost;
  int p=0;

  // bx1
  if(nb.ox1==0)     si=pmb->cis,       ei=pmb->cie+1;
  else if(nb.ox1>0) si=pmb->cie-cng+1, ei=pmb->cie;
  else              si=pmb->cis+1,     ei=pmb->cis+cng;
  if(nb.ox2==0)     sj=pmb->cjs,       ej=pmb->cje;
  else if(nb.ox2>0) sj=pmb->cje-cng+1, ej=pmb->cje;
  else              sj=pmb->cjs,       ej=pmb->cjs+cng-1;
  if(nb.ox3==0)     sk=pmb->cks,       ek=pmb->cke;
  else if(nb.ox3>0) sk=pmb->cke-cng+1, ek=pmb->cke;
  else              sk=pmb->cks,       ek=pmb->cks+cng-1;
  // include the overlapping faces in edge and corner boundaries
  if(nb.type != neighbor_face) {
    if(nb.ox1>0) ei++;
    else if(nb.ox1<0) si--;
  }
  RestrictFieldX1(src.x1f, si, ei, sj, ej, sk, ek);
  for (int k=sk; k<=ek; k++) {
    for (int j=sj; j<=ej; j++) {
#pragma simd
      for (int i=si; i<=ei; i++)
        buf[p++]=coarse_b_.x1f(k,j,i);
    }
  }

  // bx2
  if(nb.ox1==0)      si=pmb->cis,       ei=pmb->cie;
  else if(nb.ox1>0)  si=pmb->cie-cng+1, ei=pmb->cie;
  else               si=pmb->cis,       ei=pmb->cis+cng-1;
  if(pmb->block_size.nx2==1) sj=pmb->cjs, ej=pmb->cje;
  else if(nb.ox2==0) sj=pmb->cjs,       ej=pmb->cje+1;
  else if(nb.ox2>0)  sj=pmb->cje-cng+1, ej=pmb->cje;
  else               sj=pmb->cjs+1,     ej=pmb->cjs+cng;
  if(nb.type != neighbor_face) {
    if(nb.ox2>0) ej++;
    else if(nb.ox2<0) sj--;
  }
  RestrictFieldX2(src.x2f, si, ei, sj, ej, sk, ek);
  for (int k=sk; k<=ek; k++) {
    for (int j=sj; j<=ej; j++) {
#pragma simd
      for (int i=si; i<=ei; i++)
        buf[p++]=coarse_b_.x2f(k,j,i);
    }
  }

  // bx3
  if(nb.ox2==0)      sj=pmb->cjs,       ej=pmb->cje;
  else if(nb.ox2>0)  sj=pmb->cje-cng+1, ej=pmb->cje;
  else               sj=pmb->cjs,       ej=pmb->cjs+cng-1;
  if(pmb->block_size.nx3==1) sk=pmb->cks,  ek=pmb->cke;
  else if(nb.ox3==0) sk=pmb->cks,       ek=pmb->cke+1;
  else if(nb.ox3>0)  sk=pmb->cke-cng+1, ek=pmb->cke;
  else               sk=pmb->cks+1,     ek=pmb->cks+cng;
  if(nb.type != neighbor_face) {
    if(nb.ox3>0) ek++;
    else if(nb.ox3<0) sk--;
  }
  RestrictFieldX3(src.x3f, si, ei, sj, ej, sk, ek);
  for (int k=sk; k<=ek; k++) {
    for (int j=sj; j<=ej; j++) {
#pragma simd
      for (int i=si; i<=ei; i++)
        buf[p++]=coarse_b_.x3f(k,j,i);
    }
  }

  return p;
}

//--------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadFieldBoundaryBufferToFiner(InterfaceField &src, 
//                                                 Real *buf, NeighborBlock& nb)
//  \brief Set field boundary buffers for sending to a block on the finer level
int BoundaryValues::LoadFieldBoundaryBufferToFiner(InterfaceField &src, Real *buf,
                                                   NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_mblock_;
  int si, sj, sk, ei, ej, ek;
  int cn=pmb->cnghost-1;
  int p=0;

  // send the data first and later prolongate on the target block
  // need to add edges for faces, add corners for edges
  // bx1
  if(nb.ox1==0) {
    if(nb.fi1==1)   si=pmb->is+pmb->block_size.nx1/2-pmb->cnghost, ei=pmb->ie+1;
    else            si=pmb->is, ei=pmb->ie+1-pmb->block_size.nx1/2+pmb->cnghost;
  }
  else if(nb.ox1>0) si=pmb->ie, ei=pmb->ie+1;
  else              si=pmb->is, ei=pmb->is+1;
  if(nb.ox2==0) {
    sj=pmb->js,    ej=pmb->je;
    if(pmb->block_size.nx2 > 1) {
      if(nb.ox1!=0) {
        if(nb.fi1==1) sj+=pmb->block_size.nx2/2-pmb->cnghost;
        else          ej-=pmb->block_size.nx2/2-pmb->cnghost;
      }
      else {
        if(nb.fi2==1) sj+=pmb->block_size.nx2/2-pmb->cnghost;
        else          ej-=pmb->block_size.nx2/2-pmb->cnghost;
      }
    }
  }
  else if(nb.ox2>0) sj=pmb->je-cn, ej=pmb->je;
  else              sj=pmb->js,    ej=pmb->js+cn;
  if(nb.ox3==0) {
    sk=pmb->ks,    ek=pmb->ke;
    if(pmb->block_size.nx3 > 1) {
      if(nb.ox1!=0 && nb.ox2!=0) {
        if(nb.fi1==1) sk+=pmb->block_size.nx3/2-pmb->cnghost;
        else          ek-=pmb->block_size.nx3/2-pmb->cnghost;
      }
      else {
        if(nb.fi2==1) sk+=pmb->block_size.nx3/2-pmb->cnghost;
        else          ek-=pmb->block_size.nx3/2-pmb->cnghost;
      }
    }
  }
  else if(nb.ox3>0) sk=pmb->ke-cn, ek=pmb->ke;
  else              sk=pmb->ks,    ek=pmb->ks+cn;
  for (int k=sk; k<=ek; ++k) {
    for (int j=sj; j<=ej; ++j) {
#pragma simd
      for (int i=si; i<=ei; ++i)
        buf[p++]=src.x1f(k,j,i);
    }
  }

  // bx2
  if(nb.ox1==0) {
    if(nb.fi1==1)   si=pmb->is+pmb->block_size.nx1/2-pmb->cnghost, ei=pmb->ie;
    else            si=pmb->is, ei=pmb->ie-pmb->block_size.nx1/2+pmb->cnghost;
  }
  else if(nb.ox1>0) si=pmb->ie-cn, ei=pmb->ie;
  else              si=pmb->is,    ei=pmb->is+cn;
  if(nb.ox2==0) {
    sj=pmb->js,    ej=pmb->je;
    if(pmb->block_size.nx2 > 1) {
      ej++;
      if(nb.ox1!=0) {
        if(nb.fi1==1) sj+=pmb->block_size.nx2/2-pmb->cnghost;
        else          ej-=pmb->block_size.nx2/2-pmb->cnghost;
      }
      else {
        if(nb.fi2==1) sj+=pmb->block_size.nx2/2-pmb->cnghost;
        else          ej-=pmb->block_size.nx2/2-pmb->cnghost;
      }
    }
  }
  else if(nb.ox2>0) sj=pmb->je, ej=pmb->je+1;
  else              sj=pmb->js, ej=pmb->js+1;
  for (int k=sk; k<=ek; ++k) {
    for (int j=sj; j<=ej; ++j) {
#pragma simd
      for (int i=si; i<=ei; ++i)
        buf[p++]=src.x2f(k,j,i);
    }
  }

  // bx3
  if(nb.ox2==0) {
    sj=pmb->js,    ej=pmb->je;
    if(pmb->block_size.nx2 > 1) {
      if(nb.ox1!=0) {
        if(nb.fi1==1) sj+=pmb->block_size.nx2/2-pmb->cnghost;
        else          ej-=pmb->block_size.nx2/2-pmb->cnghost;
      }
      else {
        if(nb.fi2==1) sj+=pmb->block_size.nx2/2-pmb->cnghost;
        else          ej-=pmb->block_size.nx2/2-pmb->cnghost;
      }
    }
  }
  else if(nb.ox2>0) sj=pmb->je-cn, ej=pmb->je;
  else              sj=pmb->js,    ej=pmb->js+cn;
  if(nb.ox3==0) {
    sk=pmb->ks,    ek=pmb->ke;
    if(pmb->block_size.nx3 > 1) {
      ek++;
      if(nb.ox1!=0 && nb.ox2!=0) {
        if(nb.fi1==1) sk+=pmb->block_size.nx3/2-pmb->cnghost;
        else          ek-=pmb->block_size.nx3/2-pmb->cnghost;
      }
      else {
        if(nb.fi2==1) sk+=pmb->block_size.nx3/2-pmb->cnghost;
        else          ek-=pmb->block_size.nx3/2-pmb->cnghost;
      }
    }
  }
  else if(nb.ox3>0) sk=pmb->ke, ek=pmb->ke+1;
  else              sk=pmb->ks, ek=pmb->ks+1;
  for (int k=sk; k<=ek; ++k) {
    for (int j=sj; j<=ej; ++j) {
#pragma simd
      for (int i=si; i<=ei; ++i)
        buf[p++]=src.x3f(k,j,i);
    }
  }

  return p;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SendFieldBoundaryBuffers(InterfaceField &src, int step)
//  \brief Send field boundary buffers
void BoundaryValues::SendFieldBoundaryBuffers(InterfaceField &src, int step)
{
  MeshBlock *pmb=pmy_mblock_;

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    int ssize;
    if(nb.level==pmb->loc.level)
      ssize=LoadFieldBoundaryBufferSameLevel(src, field_send_[step][nb.bufid],nb);
    else if(nb.level<pmb->loc.level)
      ssize=LoadFieldBoundaryBufferToCoarser(src, field_send_[step][nb.bufid],nb);
    else
      ssize=LoadFieldBoundaryBufferToFiner(src, field_send_[step][nb.bufid], nb);
    if(nb.rank == Globals::my_rank) { // on the same process
      MeshBlock *pbl=pmb->pmy_mesh->FindMeshBlock(nb.gid);
      // find target buffer
      std::memcpy(pbl->pbval->field_recv_[step][nb.targetid],
                  field_send_[step][nb.bufid], ssize*sizeof(Real));
      pbl->pbval->field_flag_[step][nb.targetid]=boundary_arrived;
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
  MeshBlock *pmb=pmy_mblock_;
  int si, sj, sk, ei, ej, ek;

  int p=0;
  // bx1
  // for uniform grid: face-neighbors take care of the overlapping faces
  if(nb.ox1==0)     si=pmb->is,        ei=pmb->ie+1;
  else if(nb.ox1>0) si=pmb->ie+2,      ei=pmb->ie+NGHOST+1;
  else              si=pmb->is-NGHOST, ei=pmb->is-1;
  if(nb.ox2==0)     sj=pmb->js,        ej=pmb->je;
  else if(nb.ox2>0) sj=pmb->je+1,      ej=pmb->je+NGHOST;
  else              sj=pmb->js-NGHOST, ej=pmb->js-1;
  if(nb.ox3==0)     sk=pmb->ks,        ek=pmb->ke;
  else if(nb.ox3>0) sk=pmb->ke+1,      ek=pmb->ke+NGHOST;
  else              sk=pmb->ks-NGHOST, ek=pmb->ks-1;
  // for SMR/AMR, always include the overlapping faces in edge and corner boundaries
  if(pmb->pmy_mesh->multilevel==true && nb.type != neighbor_face) {
    if(nb.ox1>0) si--;
    else if(nb.ox1<0) ei++;
  }
  for (int k=sk; k<=ek; ++k) {
    for (int j=sj; j<=ej; ++j) {
#pragma simd
      for (int i=si; i<=ei; ++i)
        dst.x1f(k,j,i)=buf[p++];
    }
  }
  // bx2
  if(nb.ox1==0)      si=pmb->is,         ei=pmb->ie;
  else if(nb.ox1>0)  si=pmb->ie+1,       ei=pmb->ie+NGHOST;
  else               si=pmb->is-NGHOST,  ei=pmb->is-1;
  if(pmb->block_size.nx2==1) sj=pmb->js, ej=pmb->je;
  else if(nb.ox2==0) sj=pmb->js,         ej=pmb->je+1;
  else if(nb.ox2>0)  sj=pmb->je+2,       ej=pmb->je+NGHOST+1;
  else               sj=pmb->js-NGHOST,  ej=pmb->js-1;
  // for SMR/AMR, always include the overlapping faces in edge and corner boundaries
  if(pmb->pmy_mesh->multilevel==true && nb.type != neighbor_face) {
    if(nb.ox2>0) sj--;
    else if(nb.ox2<0) ej++;
  }
  for (int k=sk; k<=ek; ++k) {
    for (int j=sj; j<=ej; ++j) {
#pragma simd
      for (int i=si; i<=ei; ++i)
        dst.x2f(k,j,i)=buf[p++];
    }
  }
  if(pmb->block_size.nx2==1) { // 1D
#pragma simd
    for (int i=si; i<=ei; ++i)
      dst.x2f(sk,sj+1,i)=dst.x2f(sk,sj,i);
  }
  // bx3
  if(nb.ox2==0)      sj=pmb->js,         ej=pmb->je;
  else if(nb.ox2>0)  sj=pmb->je+1,       ej=pmb->je+NGHOST;
  else               sj=pmb->js-NGHOST,  ej=pmb->js-1;
  if(pmb->block_size.nx3==1) sk=pmb->ks, ek=pmb->ke;
  else if(nb.ox3==0) sk=pmb->ks,         ek=pmb->ke+1;
  else if(nb.ox3>0)  sk=pmb->ke+2,       ek=pmb->ke+NGHOST+1;
  else               sk=pmb->ks-NGHOST,  ek=pmb->ks-1;
  // for SMR/AMR, always include the overlapping faces in edge and corner boundaries
  if(pmb->pmy_mesh->multilevel==true && nb.type != neighbor_face) {
    if(nb.ox3>0) sk--;
    else if(nb.ox3<0) ek++;
  }
  for (int k=sk; k<=ek; ++k) {
    for (int j=sj; j<=ej; ++j) {
#pragma simd
      for (int i=si; i<=ei; ++i)
        dst.x3f(k,j,i)=buf[p++];
    }
  }
  if(pmb->block_size.nx3==1) { // 1D or 2D
    for (int j=sj; j<=ej; ++j) {
#pragma simd
      for (int i=si; i<=ei; ++i)
        dst.x3f(sk+1,j,i)=dst.x3f(sk,j,i);
    }
  }

  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetFieldBoundaryFromCoarser(Real *buf, NeighborBlock& nb)
//  \brief Set field prolongation buffer received from a block on the same level
void BoundaryValues::SetFieldBoundaryFromCoarser(Real *buf, NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_mblock_;
  int si, sj, sk, ei, ej, ek;
  int cng=pmb->cnghost;
  int p=0;

  // bx1
  if(nb.ox1==0) {
    si=pmb->cis, ei=pmb->cie+1;
    if((pmb->loc.lx1&1L)==0L) ei+=cng;
    else             si-=cng; 
  }
  else if(nb.ox1>0)  si=pmb->cie+1, ei=pmb->cie+2;
  else               si=pmb->cis-1, ei=pmb->cis;
  if(nb.ox2==0) {
    sj=pmb->cjs, ej=pmb->cje;
    if(pmb->block_size.nx2 > 1) {
      if((pmb->loc.lx2&1L)==0L) ej+=cng;
      else             sj-=cng; 
    }
  }
  else if(nb.ox2>0)  sj=pmb->cje+1,   ej=pmb->cje+cng;
  else               sj=pmb->cjs-cng, ej=pmb->cjs-1;
  if(nb.ox3==0) {
    sk=pmb->cks, ek=pmb->cke;
    if(pmb->block_size.nx3 > 1) {
      if((pmb->loc.lx3&1L)==0L) ek+=cng;
      else             sk-=cng; 
    }
  }
  else if(nb.ox3>0)  sk=pmb->cke+1,   ek=pmb->cke+cng;
  else               sk=pmb->cks-cng, ek=pmb->cks-1;

  for (int k=sk; k<=ek; ++k) {
    for (int j=sj; j<=ej; ++j) {
#pragma simd
      for (int i=si; i<=ei; ++i)
        coarse_b_.x1f(k,j,i) = buf[p++];
    }
  }

  // bx2
  if(nb.ox1==0) {
    si=pmb->cis, ei=pmb->cie;
    if((pmb->loc.lx1&1L)==0L) ei+=cng;
    else             si-=cng; 
  }
  else if(nb.ox1>0)  si=pmb->cie+1,   ei=pmb->cie+cng;
  else               si=pmb->cis-cng, ei=pmb->cis-1;
  if(nb.ox2==0) {
    sj=pmb->cjs, ej=pmb->cje;
    if(pmb->block_size.nx2 > 1) {
      ej++;
      if((pmb->loc.lx2&1L)==0L) ej+=cng;
      else             sj-=cng; 
    }
  }
  else if(nb.ox2>0)  sj=pmb->cje+1, ej=pmb->cje+2;
  else               sj=pmb->cjs-1, ej=pmb->cjs;

  for (int k=sk; k<=ek; ++k) {
    for (int j=sj; j<=ej; ++j) {
#pragma simd
      for (int i=si; i<=ei; ++i)
        coarse_b_.x2f(k,j,i) = buf[p++];
    }
  }

  // bx3
  if(nb.ox2==0) {
    sj=pmb->cjs, ej=pmb->cje;
    if(pmb->block_size.nx2 > 1) {
      if((pmb->loc.lx2&1L)==0L) ej+=cng;
      else             sj-=cng; 
    }
  }
  else if(nb.ox2>0)  sj=pmb->cje+1,   ej=pmb->cje+cng;
  else               sj=pmb->cjs-cng, ej=pmb->cjs-1;
  if(nb.ox3==0) {
    sk=pmb->cks, ek=pmb->cke;
    if(pmb->block_size.nx3 > 1) {
      ek++;
      if((pmb->loc.lx3&1L)==0L) ek+=cng;
      else             sk-=cng; 
    }
  }
  else if(nb.ox3>0)  sk=pmb->cke+1, ek=pmb->cke+2;
  else               sk=pmb->cks-1, ek=pmb->cks;

  for (int k=sk; k<=ek; ++k) {
    for (int j=sj; j<=ej; ++j) {
#pragma simd
      for (int i=si; i<=ei; ++i)
        coarse_b_.x3f(k,j,i) = buf[p++];
    }
  }

  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetFielBoundaryFromFiner(InterfaceField &dst,
//                                                     Real *buf, NeighborBlock& nb)
//  \brief Set field boundary received from a block on the same level
void BoundaryValues::SetFieldBoundaryFromFiner(InterfaceField &dst, Real *buf,
                                               NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_mblock_;
  // receive already restricted data
  int si, sj, sk, ei, ej, ek;
  int p=0;

  // bx1
  if(nb.ox1==0) {
    si=pmb->is, ei=pmb->ie+1;
    if(nb.fi1==1)   si+=pmb->block_size.nx1/2;
    else            ei-=pmb->block_size.nx1/2;
  }
  else if(nb.ox1>0) si=pmb->ie+2,      ei=pmb->ie+NGHOST+1;
  else              si=pmb->is-NGHOST, ei=pmb->is-1;
  // include the overlapping faces in edge and corner boundaries
  if(nb.type != neighbor_face) {
    if(nb.ox1>0) si--;
    else if(nb.ox1<0) ei++;
  }
  if(nb.ox2==0) {
    sj=pmb->js, ej=pmb->je;
    if(pmb->block_size.nx2 > 1) {
      if(nb.ox1!=0) {
        if(nb.fi1==1) sj+=pmb->block_size.nx2/2;
        else          ej-=pmb->block_size.nx2/2;
      }
      else {
        if(nb.fi2==1) sj+=pmb->block_size.nx2/2;
        else          ej-=pmb->block_size.nx2/2;
      }
    }
  }
  else if(nb.ox2>0) sj=pmb->je+1,      ej=pmb->je+NGHOST;
  else              sj=pmb->js-NGHOST, ej=pmb->js-1;
  if(nb.ox3==0) {
    sk=pmb->ks, ek=pmb->ke;
    if(pmb->block_size.nx3 > 1) {
      if(nb.ox1!=0 && nb.ox2!=0) {
        if(nb.fi1==1) sk+=pmb->block_size.nx3/2;
        else          ek-=pmb->block_size.nx3/2;
      }
      else {
        if(nb.fi2==1) sk+=pmb->block_size.nx3/2;
        else          ek-=pmb->block_size.nx3/2;
      }
    }
  }
  else if(nb.ox3>0) sk=pmb->ke+1,      ek=pmb->ke+NGHOST;
  else              sk=pmb->ks-NGHOST, ek=pmb->ks-1;

  for (int k=sk; k<=ek; ++k) {
    for (int j=sj; j<=ej; ++j) {
#pragma simd
      for (int i=si; i<=ei; ++i)
        dst.x1f(k,j,i) = buf[p++];
    }
  }

  // bx2
  if(nb.ox1==0) {
    si=pmb->is, ei=pmb->ie;
    if(nb.fi1==1)   si+=pmb->block_size.nx1/2;
    else            ei-=pmb->block_size.nx1/2;
  }
  else if(nb.ox1>0) si=pmb->ie+1,      ei=pmb->ie+NGHOST;
  else              si=pmb->is-NGHOST, ei=pmb->is-1;
  if(nb.ox2==0) {
    sj=pmb->js, ej=pmb->je;
    if(pmb->block_size.nx2 > 1) {
      ej++;
      if(nb.ox1!=0) {
        if(nb.fi1==1) sj+=pmb->block_size.nx2/2;
        else          ej-=pmb->block_size.nx2/2;
      }
      else {
        if(nb.fi2==1) sj+=pmb->block_size.nx2/2;
        else          ej-=pmb->block_size.nx2/2;
      }
    }
  }
  else if(nb.ox2>0) sj=pmb->je+2,      ej=pmb->je+NGHOST+1;
  else              sj=pmb->js-NGHOST, ej=pmb->js-1;
  // include the overlapping faces in edge and corner boundaries
  if(nb.type != neighbor_face) {
    if(nb.ox2>0) sj--;
    else if(nb.ox2<0) ej++;
  }

  for (int k=sk; k<=ek; ++k) {
    for (int j=sj; j<=ej; ++j) {
#pragma simd
      for (int i=si; i<=ei; ++i)
        dst.x2f(k,j,i) = buf[p++];
    }
  }
  if(pmb->block_size.nx2==1) { // 1D
#pragma simd
    for (int i=si; i<=ei; ++i)
      dst.x2f(sk,sj+1,i)=dst.x2f(sk,sj,i);
  }

  // bx3
  if(nb.ox2==0) {
    sj=pmb->js, ej=pmb->je;
    if(pmb->block_size.nx2 > 1) {
      if(nb.ox1!=0) {
        if(nb.fi1==1) sj+=pmb->block_size.nx2/2;
        else          ej-=pmb->block_size.nx2/2;
      }
      else {
        if(nb.fi2==1) sj+=pmb->block_size.nx2/2;
        else          ej-=pmb->block_size.nx2/2;
      }
    }
  }
  else if(nb.ox2>0) sj=pmb->je+1,      ej=pmb->je+NGHOST;
  else              sj=pmb->js-NGHOST, ej=pmb->js-1;
  if(nb.ox3==0) {
    sk=pmb->ks, ek=pmb->ke;
    if(pmb->block_size.nx3 > 1) {
      ek++;
      if(nb.ox1!=0 && nb.ox2!=0) {
        if(nb.fi1==1) sk+=pmb->block_size.nx3/2;
        else          ek-=pmb->block_size.nx3/2;
      }
      else {
        if(nb.fi2==1) sk+=pmb->block_size.nx3/2;
        else          ek-=pmb->block_size.nx3/2;
      }
    }
  }
  else if(nb.ox3>0) sk=pmb->ke+2,      ek=pmb->ke+NGHOST+1;
  else              sk=pmb->ks-NGHOST, ek=pmb->ks-1;
  // include the overlapping faces in edge and corner boundaries
  if(nb.type != neighbor_face) {
    if(nb.ox3>0) sk--;
    else if(nb.ox3<0) ek++;
  }

  for (int k=sk; k<=ek; ++k) {
    for (int j=sj; j<=ej; ++j) {
#pragma simd
      for (int i=si; i<=ei; ++i)
        dst.x3f(k,j,i) = buf[p++];
    }
  }
  if(pmb->block_size.nx3==1) { // 1D or 2D
    for (int j=sj; j<=ej; ++j) {
#pragma simd
      for (int i=si; i<=ei; ++i)
        dst.x3f(sk+1,j,i)=dst.x3f(sk,j,i);
    }
  }

  return;
}


//--------------------------------------------------------------------------------------
//! \fn bool BoundaryValues::ReceiveFieldBoundaryBuffers(InterfaceField &dst, int step)
//  \brief load boundary buffer for x1 direction into the array
bool BoundaryValues::ReceiveFieldBoundaryBuffers(InterfaceField &dst, int step)
{
  MeshBlock *pmb=pmy_mblock_;
  bool flag=true;

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    if(field_flag_[step][nb.bufid]==boundary_completed) continue;
    if(field_flag_[step][nb.bufid]==boundary_waiting) {
      if(nb.rank==Globals::my_rank) {// on the same process
        flag=false;
        continue;
      }
#ifdef MPI_PARALLEL
      else { // MPI boundary
        int test;
        MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&test,MPI_STATUS_IGNORE);
        MPI_Test(&req_field_recv_[step][nb.bufid],&test,MPI_STATUS_IGNORE);
        if(test==false) {
          flag=false;
          continue;
        }
        field_flag_[step][nb.bufid] = boundary_arrived;
      }
#endif
    }
    if(nb.level==pmb->loc.level)
      SetFieldBoundarySameLevel(dst, field_recv_[step][nb.bufid], nb);
    else if(nb.level<pmb->loc.level)
      SetFieldBoundaryFromCoarser(field_recv_[step][nb.bufid], nb);
    else
      SetFieldBoundaryFromFiner(dst, field_recv_[step][nb.bufid], nb);
    field_flag_[step][nb.bufid] = boundary_completed; // completed
  }

  return flag;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ReceiveFieldBoundaryBuffersWithWait(InterfaceField &dst,
//                                                               int step)
//  \brief load boundary buffer for x1 direction into the array
void BoundaryValues::ReceiveFieldBoundaryBuffersWithWait(InterfaceField &dst, int step)
{
  MeshBlock *pmb=pmy_mblock_;

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb = pmb->neighbor[n];
#ifdef MPI_PARALLEL
    if(nb.rank!=Globals::my_rank)
      MPI_Wait(&req_field_recv_[0][nb.bufid],MPI_STATUS_IGNORE);
#endif
    if(nb.level==pmb->loc.level)
      SetFieldBoundarySameLevel(dst, field_recv_[0][nb.bufid], nb);
    else if(nb.level<pmb->loc.level)
      SetFieldBoundaryFromCoarser(field_recv_[0][nb.bufid], nb);
    else
      SetFieldBoundaryFromFiner(dst, field_recv_[0][nb.bufid], nb);
    field_flag_[0][nb.bufid] = boundary_completed; // completed
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadEMFBoundaryBufferSameLevel(Real *buf, NeighborBlock& nb)
//  \brief Set EMF correction buffers for sending to a block on the same level
int BoundaryValues::LoadEMFBoundaryBufferSameLevel(Real *buf, NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_mblock_;
  AthenaArray<Real> &e1=pmb->pfield->e.x1e;
  AthenaArray<Real> &e2=pmb->pfield->e.x2e;
  AthenaArray<Real> &e3=pmb->pfield->e.x3e;
  int p=0;
  if(nb.type==neighbor_face) {
    if(pmb->block_size.nx3 > 1) { // 3D
      // x1 direction
      if(nb.fid==inner_x1 || nb.fid==outer_x1) {
        int i;
        if(nb.fid==inner_x1) i=pmb->is;
        else i=pmb->ie+1;
        // pack e2
        for(int k=pmb->ks; k<=pmb->ke+1; k++) {
          for(int j=pmb->js; j<=pmb->je; j++)
            buf[p++]=e2(k,j,i);
        }
        // pack e3
        for(int k=pmb->ks; k<=pmb->ke; k++) {
          for(int j=pmb->js; j<=pmb->je+1; j++)
            buf[p++]=e3(k,j,i);
        }
      }
      // x2 direction
      else if(nb.fid==inner_x2 || nb.fid==outer_x2) {
        int j;
        if(nb.fid==inner_x2) j=pmb->js;
        else j=pmb->je+1;
        // pack e1
        for(int k=pmb->ks; k<=pmb->ke+1; k++) {
          for(int i=pmb->is; i<=pmb->ie; i++)
            buf[p++]=e1(k,j,i);
        }
        // pack e3
        for(int k=pmb->ks; k<=pmb->ke; k++) {
          for(int i=pmb->is; i<=pmb->ie+1; i++)
            buf[p++]=e3(k,j,i);
        }
      }
      // x3 direction
      else if(nb.fid==inner_x3 || nb.fid==outer_x3) {
        int k;
        if(nb.fid==inner_x3) k=pmb->ks;
        else k=pmb->ke+1;
        // pack e1
        for(int j=pmb->js; j<=pmb->je+1; j++) {
          for(int i=pmb->is; i<=pmb->ie; i++)
            buf[p++]=e1(k,j,i);
        }
        // pack e2
        for(int j=pmb->js; j<=pmb->je; j++) {
          for(int i=pmb->is; i<=pmb->ie+1; i++)
            buf[p++]=e2(k,j,i);
        }
      }
    }
    else if(pmb->block_size.nx2 > 1) { // 2D
      int k=pmb->ks;
      // x1 direction
      if(nb.fid==inner_x1 || nb.fid==outer_x1) {
        int i;
        if(nb.fid==inner_x1) i=pmb->is;
        else i=pmb->ie+1;
        // pack e2
        for(int j=pmb->js; j<=pmb->je; j++)
          buf[p++]=e2(k,j,i);
        // pack e3
        for(int j=pmb->js; j<=pmb->je+1; j++)
          buf[p++]=e3(k,j,i);
      }
      // x2 direction
      else if(nb.fid==inner_x2 || nb.fid==outer_x2) {
        int j;
        if(nb.fid==inner_x2) j=pmb->js;
        else j=pmb->je+1;
        // pack e1
        for(int i=pmb->is; i<=pmb->ie; i++)
          buf[p++]=e1(k,j,i);
        // pack e3
        for(int i=pmb->is; i<=pmb->ie+1; i++)
          buf[p++]=e3(k,j,i);
      }
    }
    else { // 1D
      int i, j=pmb->js, k=pmb->ks;
      if(nb.fid==inner_x1) i=pmb->is;
      else i=pmb->ie+1;
      // pack e2 and e3
      buf[p++]=e2(k,j,i);
      buf[p++]=e3(k,j,i);
    }
  }
  else if(nb.type==neighbor_edge) {
    // x1x2 edge (both 2D and 3D)
    if(nb.eid>=0 && nb.eid<4) {
      int i, j;
      if((nb.eid&1)==0) i=pmb->is;
      else i=pmb->ie+1;
      if((nb.eid&2)==0) j=pmb->js;
      else j=pmb->je+1;
      // pack e3
      for(int k=pmb->ks; k<=pmb->ke; k++)
        buf[p++]=e3(k,j,i);
    }
    // x1x3 edge
    else if(nb.eid>=4 && nb.eid<8) {
      int i, k;
      if((nb.eid&1)==0) i=pmb->is;
      else i=pmb->ie+1;
      if((nb.eid&2)==0) k=pmb->ks;
      else k=pmb->ke+1;
      // pack e2
      for(int j=pmb->js; j<=pmb->je; j++)
        buf[p++]=e2(k,j,i);
    }
    // x2x3 edge
    else if(nb.eid>=8 && nb.eid<12) {
      int j, k;
      if((nb.eid&1)==0) j=pmb->js;
      else j=pmb->je+1;
      if((nb.eid&2)==0) k=pmb->ks;
      else k=pmb->ke+1;
      // pack e1
      for(int i=pmb->is; i<=pmb->ie; i++)
        buf[p++]=e1(k,j,i);
    }
  }
  return p;
}


//--------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadEMFBoundaryBufferToCoarser(Real *buf, NeighborBlock& nb)
//  \brief Set EMF correction buffers for sending to a block on the coarser level
int BoundaryValues::LoadEMFBoundaryBufferToCoarser(Real *buf, NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_mblock_;
  Coordinates *pco=pmb->pcoord;
  AthenaArray<Real> &e1=pmb->pfield->e.x1e;
  AthenaArray<Real> &e2=pmb->pfield->e.x2e;
  AthenaArray<Real> &e3=pmb->pfield->e.x3e;
  // use the surface area aray as the edge length array
  AthenaArray<Real> &le1=sarea_[0];
  AthenaArray<Real> &le2=sarea_[1];
  int p=0;
  if(nb.type==neighbor_face) {
    if(pmb->block_size.nx3 > 1) { // 3D
      // x1 direction
      if(nb.fid==inner_x1 || nb.fid==outer_x1) {
        int i;
        if(nb.fid==inner_x1) i=pmb->is;
        else i=pmb->ie+1;
        // restrict and pack e2
        for(int k=pmb->ks; k<=pmb->ke+1; k+=2) {
          for(int j=pmb->js; j<=pmb->je; j+=2) {
            Real el1=pco->GetEdge2Length(k,j,i);
            Real el2=pco->GetEdge2Length(k,j+1,i);
            buf[p++]=(e2(k,j,i)*el1+e2(k,j+1,i)*el2)/(el1+el2);
          }
        }
        // restrict and pack e3
        for(int k=pmb->ks; k<=pmb->ke; k+=2) {
          for(int j=pmb->js; j<=pmb->je+1; j+=2) {
            Real el1=pco->GetEdge3Length(k,j,i);
            Real el2=pco->GetEdge3Length(k+1,j,i);
            buf[p++]=(e3(k,j,i)*el1+e3(k+1,j,i)*el2)/(el1+el2);
          }
        }
      }
      // x2 direction
      else if(nb.fid==inner_x2 || nb.fid==outer_x2) {
        int j;
        if(nb.fid==inner_x2) j=pmb->js;
        else j=pmb->je+1;
        // restrict and pack e1
        for(int k=pmb->ks; k<=pmb->ke+1; k+=2) {
          pco->Edge1Length(k, j, pmb->is, pmb->ie, le1);
          for(int i=pmb->is; i<=pmb->ie; i+=2)
            buf[p++]=(e1(k,j,i)*le1(i)+e1(k,j,i+1)*le1(i+1))/(le1(i)+le1(i+1));
        }
        // restrict and pack e3
        for(int k=pmb->ks; k<=pmb->ke; k+=2) {
          pco->Edge3Length(k,   j, pmb->is, pmb->ie+1, le1);
          pco->Edge3Length(k+1, j, pmb->is, pmb->ie+1, le2);
          for(int i=pmb->is; i<=pmb->ie+1; i+=2)
            buf[p++]=(e3(k,j,i)*le1(i)+e3(k+1,j,i)*le2(i))/(le1(i)+le2(i));
        }
      }
      // x3 direction
      else if(nb.fid==inner_x3 || nb.fid==outer_x3) {
        int k;
        if(nb.fid==inner_x3) k=pmb->ks;
        else k=pmb->ke+1;
        // restrict and pack e1
        for(int j=pmb->js; j<=pmb->je+1; j+=2) {
          pco->Edge1Length(k, j, pmb->is, pmb->ie, le1);
          for(int i=pmb->is; i<=pmb->ie; i+=2)
            buf[p++]=(e1(k,j,i)*le1(i)+e1(k,j,i+1)*le1(i+1))/(le1(i)+le1(i+1));
        }
        // restrict and pack e2
        for(int j=pmb->js; j<=pmb->je; j+=2) {
          pco->Edge2Length(k,   j, pmb->is, pmb->ie+1, le1);
          pco->Edge2Length(k, j+1, pmb->is, pmb->ie+1, le2);
          for(int i=pmb->is; i<=pmb->ie+1; i+=2)
            buf[p++]=(e2(k,j,i)*le1(i)+e2(k,j+1,i)*le2(i))/(le1(i)+le2(i));
        }
      }
    }
    else if(pmb->block_size.nx2 > 1) { // 2D
      int k=pmb->ks;
      // x1 direction
      if(nb.fid==inner_x1 || nb.fid==outer_x1) {
        int i;
        if(nb.fid==inner_x1) i=pmb->is;
        else i=pmb->ie+1;
        // restrict and pack e2
        for(int j=pmb->js; j<=pmb->je; j+=2) {
          Real el1=pco->GetEdge2Length(k,j,i);
          Real el2=pco->GetEdge2Length(k,j+1,i);
          buf[p++]=(e2(k,j,i)*el1+e2(k,j+1,i)*el2)/(el1+el2);
        }
        // pack e3
        for(int j=pmb->js; j<=pmb->je+1; j+=2)
          buf[p++]=e3(k,j,i);
      }
      // x2 direction
      else if(nb.fid==inner_x2 || nb.fid==outer_x2) {
        int j;
        if(nb.fid==inner_x2) j=pmb->js;
        else j=pmb->je+1;
        // restrict and pack e1
        pco->Edge1Length(k, j, pmb->is, pmb->ie, le1);
        for(int i=pmb->is; i<=pmb->ie; i+=2)
          buf[p++]=(e1(k,j,i)*le1(i)+e1(k,j,i+1)*le1(i+1))/(le1(i)+le1(i+1));
        // pack e3
        for(int i=pmb->is; i<=pmb->ie+1; i+=2)
          buf[p++]=e3(k,j,i);
      }
    }
    else { // 1D
      int i, j=pmb->js, k=pmb->ks;
      if(nb.fid==inner_x1) i=pmb->is;
      else i=pmb->ie+1;
      // pack e2 and e3
      buf[p++]=e2(k,j,i);
      buf[p++]=e3(k,j,i);
    }
  }
  else if(nb.type==neighbor_edge) {
    if(pmb->block_size.nx3 > 1) { // 3D
      // x1x2 edge
      if(nb.eid>=0 && nb.eid<4) {
        int i, j;
        if((nb.eid&1)==0) i=pmb->is;
        else i=pmb->ie+1;
        if((nb.eid&2)==0) j=pmb->js;
        else j=pmb->je+1;
        // restrict and pack e3
        for(int k=pmb->ks; k<=pmb->ke; k+=2) {
          Real el1=pco->GetEdge3Length(k,j,i);
          Real el2=pco->GetEdge3Length(k+1,j,i);
          buf[p++]=(e3(k,j,i)*el1+e3(k+1,j,i)*el2)/(el1+el2);
        }
      }
      // x1x3 edge
      else if(nb.eid>=4 && nb.eid<8) {
        int i, k;
        if((nb.eid&1)==0) i=pmb->is;
        else i=pmb->ie+1;
        if((nb.eid&2)==0) k=pmb->ks;
        else k=pmb->ke+1;
        // restrict and pack e2
        for(int j=pmb->js; j<=pmb->je; j+=2) {
          Real el1=pco->GetEdge2Length(k,j,i);
          Real el2=pco->GetEdge2Length(k,j+1,i);
          buf[p++]=(e2(k,j,i)*el1+e2(k,j+1,i)*el2)/(el1+el2);
        }
      }
      // x2x3 edge
      else if(nb.eid>=8 && nb.eid<12) {
        int j, k;
        if((nb.eid&1)==0) j=pmb->js;
        else j=pmb->je+1;
        if((nb.eid&2)==0) k=pmb->ks;
        else k=pmb->ke+1;
        // restrict and pack e1
        pco->Edge1Length(k, j, pmb->is, pmb->ie, le1);
        for(int i=pmb->is; i<=pmb->ie; i+=2)
          buf[p++]=(e1(k,j,i)*le1(i)+e1(k,j,i+1)*le1(i+1))/(le1(i)+le1(i+1));
      }
    }
    else if(pmb->block_size.nx2 > 1) { // 2D
      // x1x2 edge
      int i, j;
      if((nb.eid&1)==0) i=pmb->is;
      else i=pmb->ie+1;
      if((nb.eid&2)==0) j=pmb->js;
      else j=pmb->je+1;
      // pack e3
      buf[p++]=e3(pmb->ks,j,i);
    }
  }
  return p;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SendEMFCorrection(int step)
//  \brief Restrict, pack and send the surace EMF to the coarse neighbor(s) if needed
void BoundaryValues::SendEMFCorrection(int step)
{
  MeshBlock *pmb=pmy_mblock_;

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    if((nb.type!=neighbor_face) && (nb.type!=neighbor_edge)) break;
    int p=0;
    if(nb.level==pmb->loc.level) {
      if((nb.type==neighbor_face)
      || ((nb.type==neighbor_edge) && (edge_flag_[nb.eid]==true)))
        p=LoadEMFBoundaryBufferSameLevel(emfcor_send_[step][nb.bufid], nb);
      else continue;
    }
    else if(nb.level==pmb->loc.level-1)
      p=LoadEMFBoundaryBufferToCoarser(emfcor_send_[step][nb.bufid], nb);
    else continue;
    if(nb.rank==Globals::my_rank) { // on the same node
      MeshBlock *pbl=pmb->pmy_mesh->FindMeshBlock(nb.gid);
      std::memcpy(pbl->pbval->emfcor_recv_[step][nb.targetid],
                  emfcor_send_[step][nb.bufid], p*sizeof(Real));
      pbl->pbval->emfcor_flag_[step][nb.targetid]=boundary_arrived;
    }
#ifdef MPI_PARALLEL
    else
      MPI_Start(&req_emfcor_send_[step][nb.bufid]);
#endif
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetEMFBoundarySameLevel(Real *buf, NeighborBlock& nb)
//  \brief Add up the EMF received from a block on the same level
//         Later they will be divided in the AverageEMFBoundary function
void BoundaryValues::SetEMFBoundarySameLevel(Real *buf, NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_mblock_;
  AthenaArray<Real> &e1=pmb->pfield->e.x1e;
  AthenaArray<Real> &e2=pmb->pfield->e.x2e;
  AthenaArray<Real> &e3=pmb->pfield->e.x3e;
  int p=0;
  if(nb.type==neighbor_face) {
    if(pmb->block_size.nx3 > 1) { // 3D
      // x1 direction
      if(nb.fid==inner_x1 || nb.fid==outer_x1) {
        int i;
        if(nb.fid==inner_x1) i=pmb->is;
        else i=pmb->ie+1;
        // unpack e2
        for(int k=pmb->ks; k<=pmb->ke+1; k++) {
          for(int j=pmb->js; j<=pmb->je; j++)
            e2(k,j,i)+=buf[p++];
        }
        // unpack e3
        for(int k=pmb->ks; k<=pmb->ke; k++) {
          for(int j=pmb->js; j<=pmb->je+1; j++)
            e3(k,j,i)+=buf[p++];
        }
      }
      // x2 direction
      else if(nb.fid==inner_x2 || nb.fid==outer_x2) {
        int j;
        if(nb.fid==inner_x2) j=pmb->js;
        else j=pmb->je+1;
        // unpack e1
        for(int k=pmb->ks; k<=pmb->ke+1; k++) {
          for(int i=pmb->is; i<=pmb->ie; i++)
            e1(k,j,i)+=buf[p++];
        }
        // unpack e3
        for(int k=pmb->ks; k<=pmb->ke; k++) {
          for(int i=pmb->is; i<=pmb->ie+1; i++)
            e3(k,j,i)+=buf[p++];
        }
      }
      // x3 direction
      else if(nb.fid==inner_x3 || nb.fid==outer_x3) {
        int k;
        if(nb.fid==inner_x3) k=pmb->ks;
        else k=pmb->ke+1;
        // unpack e1
        for(int j=pmb->js; j<=pmb->je+1; j++) {
          for(int i=pmb->is; i<=pmb->ie; i++)
            e1(k,j,i)+=buf[p++];
        }
        // unpack e2
        for(int j=pmb->js; j<=pmb->je; j++) {
          for(int i=pmb->is; i<=pmb->ie+1; i++)
            e2(k,j,i)+=buf[p++];
        }
      }
    }
    else if(pmb->block_size.nx2 > 1) { // 2D
      int k=pmb->ks;
      // x1 direction
      if(nb.fid==inner_x1 || nb.fid==outer_x1) {
        int i;
        if(nb.fid==inner_x1) i=pmb->is;
        else i=pmb->ie+1;
        // unpack e2
        for(int j=pmb->js; j<=pmb->je; j++) {
          e2(k+1,j,i)+=buf[p];
          e2(k,  j,i)+=buf[p++];
        }
        // unpack e3
        for(int j=pmb->js; j<=pmb->je+1; j++)
          e3(k,j,i)+=buf[p++];
      }
      // x2 direction
      else if(nb.fid==inner_x2 || nb.fid==outer_x2) {
        int j;
        if(nb.fid==inner_x2) j=pmb->js;
        else j=pmb->je+1;
        // unpack e1
        for(int i=pmb->is; i<=pmb->ie; i++) {
          e1(k+1,j,i)+=buf[p];
          e1(k  ,j,i)+=buf[p++];
        }
        // unpack e3
        for(int i=pmb->is; i<=pmb->ie+1; i++)
          e3(k,j,i)+=buf[p++];
      }
    }
    else { // 1D
      int i, j=pmb->js, k=pmb->ks;
      if(nb.fid==inner_x1) i=pmb->is;
      else i=pmb->ie+1;
      // unpack e2
      e2(k+1,j,i)+=buf[p];
      e2(k  ,j,i)+=buf[p++];
      // unpack e3
      e3(k,j+1,i)+=buf[p];
      e3(k  ,j,i)+=buf[p++];
    }
  }
  else if(nb.type==neighbor_edge) {
    if(pmb->block_size.nx3 > 1) { // 3D
      // x1x2 edge
      if(nb.eid>=0 && nb.eid<4) {
        int i, j;
        if((nb.eid&1)==0) i=pmb->is;
        else i=pmb->ie+1;
        if((nb.eid&2)==0) j=pmb->js;
        else j=pmb->je+1;
        // unpack e3
        for(int k=pmb->ks; k<=pmb->ke; k++)
          e3(k,j,i)+=buf[p++];
      }
      // x1x3 edge
      else if(nb.eid>=4 && nb.eid<8) {
        int i, k;
        if((nb.eid&1)==0) i=pmb->is;
        else i=pmb->ie+1;
        if((nb.eid&2)==0) k=pmb->ks;
        else k=pmb->ke+1;
        // unpack e2
        for(int j=pmb->js; j<=pmb->je; j++)
          e2(k,j,i)+=buf[p++];
      }
      // x2x3 edge
      else if(nb.eid>=8 && nb.eid<12) {
        int j, k;
        if((nb.eid&1)==0) j=pmb->js;
        else j=pmb->je+1;
        if((nb.eid&2)==0) k=pmb->ks;
        else k=pmb->ke+1;
        // unpack e1
        for(int i=pmb->is; i<=pmb->ie; i++)
          e1(k,j,i)+=buf[p++];
      }
    }
    else if(pmb->block_size.nx2 > 1) { // 2D
      int i, j, k=pmb->ks;
      if((nb.eid&1)==0) i=pmb->is;
      else i=pmb->ie+1;
      if((nb.eid&2)==0) j=pmb->js;
      else j=pmb->je+1;
      // unpack e3
      e3(k,j,i)+=buf[p++];
    }
  }

  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetEMFBoundaryFromFiner(Real *buf, NeighborBlock& nb)
//  \brief Add up the EMF received from a block on the finer level
//         Later they will be divided in the AverageEMFBoundary function
void BoundaryValues::SetEMFBoundaryFromFiner(Real *buf, NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_mblock_;
  AthenaArray<Real> &e1=pmb->pfield->e.x1e;
  AthenaArray<Real> &e2=pmb->pfield->e.x2e;
  AthenaArray<Real> &e3=pmb->pfield->e.x3e;
  int p=0;
  if(nb.type==neighbor_face) {
    if(pmb->block_size.nx3 > 1) { // 3D
      // x1 direction
      if(nb.fid==inner_x1 || nb.fid==outer_x1) {
        int i, jl=pmb->js, ju=pmb->je, kl=pmb->ks, ku=pmb->ke;
        if(nb.fid==inner_x1) i=pmb->is;
        else i=pmb->ie+1;
        if(nb.fi1==0) ju=pmb->js+pmb->block_size.nx2/2-1;
        else jl=pmb->js+pmb->block_size.nx2/2;
        if(nb.fi2==0) ku=pmb->ks+pmb->block_size.nx3/2-1;
        else kl=pmb->ks+pmb->block_size.nx3/2;
        // unpack e2
        for(int k=kl; k<=ku+1; k++) {
          for(int j=jl; j<=ju; j++)
            e2(k,j,i)+=buf[p++];
        }
        // unpack e3
        for(int k=kl; k<=ku; k++) {
          for(int j=jl; j<=ju+1; j++)
            e3(k,j,i)+=buf[p++];
        }
      }
      // x2 direction
      else if(nb.fid==inner_x2 || nb.fid==outer_x2) {
        int j, il=pmb->is, iu=pmb->ie, kl=pmb->ks, ku=pmb->ke;
        if(nb.fid==inner_x2) j=pmb->js;
        else j=pmb->je+1;
        if(nb.fi1==0) iu=pmb->is+pmb->block_size.nx1/2-1;
        else il=pmb->is+pmb->block_size.nx1/2;
        if(nb.fi2==0) ku=pmb->ks+pmb->block_size.nx3/2-1;
        else kl=pmb->ks+pmb->block_size.nx3/2;
        // unpack e1
        for(int k=kl; k<=ku+1; k++) {
          for(int i=il; i<=iu; i++)
            e1(k,j,i)+=buf[p++];
        }
        // unpack e3
        for(int k=kl; k<=ku; k++) {
          for(int i=il; i<=iu+1; i++)
            e3(k,j,i)+=buf[p++];
        }
      }
      // x3 direction
      else if(nb.fid==inner_x3 || nb.fid==outer_x3) {
        int k, il=pmb->is, iu=pmb->ie, jl=pmb->js, ju=pmb->je;
        if(nb.fid==inner_x3) k=pmb->ks;
        else k=pmb->ke+1;
        if(nb.fi1==0) iu=pmb->is+pmb->block_size.nx1/2-1;
        else il=pmb->is+pmb->block_size.nx1/2;
        if(nb.fi2==0) ju=pmb->js+pmb->block_size.nx2/2-1;
        else jl=pmb->js+pmb->block_size.nx2/2;
        // unpack e1
        for(int j=jl; j<=ju+1; j++) {
          for(int i=il; i<=iu; i++)
            e1(k,j,i)+=buf[p++];
        }
        // unpack e2
        for(int j=jl; j<=ju; j++) {
          for(int i=il; i<=iu+1; i++)
            e2(k,j,i)+=buf[p++];
        }
      }
    }
    else if(pmb->block_size.nx2 > 1) { // 2D
      int k=pmb->ks;
      // x1 direction
      if(nb.fid==inner_x1 || nb.fid==outer_x1) {
        int i, jl=pmb->js, ju=pmb->je;
        if(nb.fid==inner_x1) i=pmb->is;
        else i=pmb->ie+1;
        if(nb.fi1==0) ju=pmb->js+pmb->block_size.nx2/2-1;
        else jl=pmb->js+pmb->block_size.nx2/2;
        // unpack e2
        for(int j=jl; j<=ju; j++) {
          e2(k+1,j,i)+=buf[p];
          e2(k,  j,i)+=buf[p++];
        }
        // unpack e3
        for(int j=jl; j<=ju+1; j++)
          e3(k,j,i)+=buf[p++];
      }
      // x2 direction
      else if(nb.fid==inner_x2 || nb.fid==outer_x2) {
        int j, il=pmb->is, iu=pmb->ie;
        if(nb.fid==inner_x2) j=pmb->js;
        else j=pmb->je+1;
        if(nb.fi1==0) iu=pmb->is+pmb->block_size.nx1/2-1;
        else il=pmb->is+pmb->block_size.nx1/2;
        // unpack e1
        for(int i=il; i<=iu; i++) {
          e1(k+1,j,i)+=buf[p];
          e1(k  ,j,i)+=buf[p++];
        }
        // unpack e3
        for(int i=il; i<=iu+1; i++)
          e3(k,j,i)+=buf[p++];
      }
    }
    else { // 1D
      int i, j=pmb->js, k=pmb->ks;
      if(nb.fid==inner_x1) i=pmb->is;
      else i=pmb->ie+1;
      // unpack e2
      e2(k+1,j,i)+=buf[p];
      e2(k  ,j,i)+=buf[p++];
      // unpack e3
      e3(k,j+1,i)+=buf[p];
      e3(k  ,j,i)+=buf[p++];
    }
  }
  else if(nb.type==neighbor_edge) {
    if(pmb->block_size.nx3 > 1) { // 3D
      // x1x2 edge
      if(nb.eid>=0 && nb.eid<4) {
        int i, j, kl=pmb->ks, ku=pmb->ke;
        if((nb.eid&1)==0) i=pmb->is;
        else i=pmb->ie+1;
        if((nb.eid&2)==0) j=pmb->js;
        else j=pmb->je+1;
        if(nb.fi1==0) ku=pmb->ks+pmb->block_size.nx3/2-1;
        else kl=pmb->ks+pmb->block_size.nx3/2;
        // unpack e3
        for(int k=kl; k<=ku; k++)
          e3(k,j,i)+=buf[p++];
      }
      // x1x3 edge
      else if(nb.eid>=4 && nb.eid<8) {
        int i, k, jl=pmb->js, ju=pmb->je;
        if((nb.eid&1)==0) i=pmb->is;
        else i=pmb->ie+1;
        if((nb.eid&2)==0) k=pmb->ks;
        else k=pmb->ke+1;
        if(nb.fi1==0) ju=pmb->js+pmb->block_size.nx2/2-1;
        else jl=pmb->js+pmb->block_size.nx2/2;
        // unpack e2
        for(int j=jl; j<=ju; j++)
          e2(k,j,i)+=buf[p++];
      }
      // x2x3 edge
      else if(nb.eid>=8 && nb.eid<12) {
        int j, k, il=pmb->is, iu=pmb->ie;
        if((nb.eid&1)==0) j=pmb->js;
        else j=pmb->je+1;
        if((nb.eid&2)==0) k=pmb->ks;
        else k=pmb->ke+1;
        if(nb.fi1==0) iu=pmb->is+pmb->block_size.nx1/2-1;
        else il=pmb->is+pmb->block_size.nx1/2;
        // unpack e1
        for(int i=il; i<=iu; i++)
          e1(k,j,i)+=buf[p++];
      }
    }
    else if(pmb->block_size.nx2 > 1) { // 2D
      int i, j, k=pmb->ks;
      if((nb.eid&1)==0) i=pmb->is;
      else i=pmb->ie+1;
      if((nb.eid&2)==0) j=pmb->js;
      else j=pmb->je+1;
      // unpack e3
      e3(k,j,i)+=buf[p++];
    }
  }

  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ClearCoarseEMFBoundary(void)
//  \brief Clear the EMFs on the surface/edge contacting with a finer block
void BoundaryValues::ClearCoarseEMFBoundary(void)
{
  MeshBlock *pmb=pmy_mblock_;
  AthenaArray<Real> &e1=pmb->pfield->e.x1e;
  AthenaArray<Real> &e2=pmb->pfield->e.x2e;
  AthenaArray<Real> &e3=pmb->pfield->e.x3e;
  int i, j, k, nl;
  // face
  for(int n=0; n<nface_; n++) {
    if(n==inner_x1 || n==outer_x1) {
      if(n==inner_x1) i=pmb->is;
      else i=pmb->ie+1;
      nl=pmb->nblevel[1][1][2*n];
      if(nl>pmb->loc.level) { // finer
        if(pmb->block_size.nx3 > 1) { // 3D
          for(int k=pmb->ks+1; k<=pmb->ke; k++) {
            for(int j=pmb->js; j<=pmb->je; j++)
              e2(k,j,i)=0.0;
          }
          for(int k=pmb->ks; k<=pmb->ke; k++) {
            for(int j=pmb->js+1; j<=pmb->je; j++)
              e3(k,j,i)=0.0;
          }
        }
        else if(pmb->block_size.nx2 > 1) { // 2D
          for(int j=pmb->js; j<=pmb->je; j++)
            e2(pmb->ks,j,i)=e2(pmb->ks+1,j,i)=0.0;
          for(int j=pmb->js+1; j<=pmb->je; j++)
            e3(pmb->ks,j,i)=0.0;
        }
        else { // 1D
          e2(pmb->ks,pmb->js,i)=e2(pmb->ks+1,pmb->js,i)=0.0;
          e3(pmb->ks,pmb->js,i)=e3(pmb->ks,pmb->js+1,i)=0.0;
        }
      }
    }
    if(n==inner_x2 || n==outer_x2) {
      if(n==inner_x2) j=pmb->js;
      else j=pmb->je+1;
      nl=pmb->nblevel[1][2*n-4][1];
      if(nl>pmb->loc.level) { // finer
        if(pmb->block_size.nx3 > 1) { // 3D
          for(int k=pmb->ks+1; k<=pmb->ke; k++) {
            for(int i=pmb->is; i<=pmb->ie; i++)
              e1(k,j,i)=0.0;
          }
          for(int k=pmb->ks; k<=pmb->ke; k++) {
            for(int i=pmb->is+1; i<=pmb->ie; i++)
              e3(k,j,i)=0.0;
          }
        }
        else if(pmb->block_size.nx2 > 1) { // 2D
          for(int i=pmb->is; i<=pmb->ie; i++)
            e1(pmb->ks,j,i)=e1(pmb->ks+1,j,i)=0.0;
          for(int i=pmb->is+1; i<=pmb->ie; i++)
            e3(pmb->ks,j,i)=0.0;
        }
      }
    }
    if(n==inner_x3 || n==outer_x3) {
      if(n==inner_x3) k=pmb->ks;
      else k=pmb->ke+1;
      nl=pmb->nblevel[2*n-8][1][1];
      if(nl>pmb->loc.level) { // finer
        // this is always 3D
        for(int j=pmb->js+1; j<=pmb->je; j++) {
          for(int i=pmb->is; i<=pmb->ie; i++)
            e1(k,j,i)=0.0;
        }
        for(int j=pmb->js; j<=pmb->je; j++) {
          for(int i=pmb->is+1; i<=pmb->ie; i++)
            e2(k,j,i)=0.0;
        }
      }
    }
  }
  // edge
  for(int n=0; n<nedge_; n++) {
    if(edge_flag_[n]==true) continue;
    if(n>=0 && n<4) {
      if((n&1)==0) i=pmb->is;
      else i=pmb->ie+1;
      if((n&2)==0) j=pmb->js;
      else j=pmb->je+1;
      for(int k=pmb->ks; k<=pmb->ke; k++)
        e3(k,j,i)=0.0;
    }
    // x1x3 edge
    else if(n>=4 && n<8) {
      if((n&1)==0) i=pmb->is;
      else i=pmb->ie+1;
      if((n&2)==0) k=pmb->ks;
      else k=pmb->ke+1;
      for(int j=pmb->js; j<=pmb->je; j++)
        e2(k,j,i)=0.0;
    }
    // x2x3 edge
    else if(n>=8 && n<12) {
      if((n&1)==0) j=pmb->js;
      else j=pmb->je+1;
      if((n&2)==0) k=pmb->ks;
      else k=pmb->ke+1;
      for(int i=pmb->is; i<=pmb->ie; i++)
        e1(k,j,i)=0.0;
    }
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::AverageEMFBoundary(void)
//  \brief Set EMF boundary received from a block on the finer level
void BoundaryValues::AverageEMFBoundary(void)
{
  MeshBlock *pmb=pmy_mblock_;
  AthenaArray<Real> &e1=pmb->pfield->e.x1e;
  AthenaArray<Real> &e2=pmb->pfield->e.x2e;
  AthenaArray<Real> &e3=pmb->pfield->e.x3e;
  int i, j, k, nl;
  // face
  for(int n=0; n<nface_; n++) {
    if ((pmb->block_bcs[n] != -1) && (pmb->block_bcs[n] != 4)) continue;
    if(n==inner_x1 || n==outer_x1) {
      if(n==inner_x1) i=pmb->is;
      else i=pmb->ie+1;
      nl=pmb->nblevel[1][1][2*n];
      if(nl==pmb->loc.level) { // same ; divide all the face EMFs by 2
        if(pmb->block_size.nx3 > 1) { // 3D
          for(int k=pmb->ks+1; k<=pmb->ke; k++) {
            for(int j=pmb->js; j<=pmb->je; j++)
              e2(k,j,i)*=0.5;
          }
          for(int k=pmb->ks; k<=pmb->ke; k++) {
            for(int j=pmb->js+1; j<=pmb->je; j++)
              e3(k,j,i)*=0.5;
          }
        }
        else if(pmb->block_size.nx2 > 1) { // 2D
          for(int j=pmb->js; j<=pmb->je; j++)
            e2(pmb->ks,j,i)*=0.5, e2(pmb->ks+1,j,i)*=0.5;
          for(int j=pmb->js+1; j<=pmb->je; j++)
            e3(pmb->ks,j,i)*=0.5;
        }
        else { // 1D
          e2(pmb->ks,pmb->js,i)*=0.5, e2(pmb->ks+1,pmb->js,i)*=0.5;
          e3(pmb->ks,pmb->js,i)*=0.5, e3(pmb->ks,pmb->js+1,i)*=0.5;
        }
      }
      else if(nl>pmb->loc.level) { // finer; divide the overlapping EMFs by 2
        if(pmb->block_size.nx3 > 1) { // 3D
          int k=pmb->ks+pmb->block_size.nx3/2;
          for(int j=pmb->js; j<=pmb->je; j++)
            e2(k,j,i)*=0.5;
        }
        if(pmb->block_size.nx2 > 1) { // 2D or 3D
          int j=pmb->js+pmb->block_size.nx2/2;
          for(int k=pmb->ks; k<=pmb->ke; k++)
            e3(k,j,i)*=0.5;
        }
      }
    }
    if(n==inner_x2 || n==outer_x2) {
      if(n==inner_x2) j=pmb->js;
      else j=pmb->je+1;
      nl=pmb->nblevel[1][2*n-4][1];
      if(nl==pmb->loc.level) { // same ; divide all the face EMFs by 2
        if(pmb->block_size.nx3 > 1) {
          for(int k=pmb->ks+1; k<=pmb->ke; k++) {
            for(int i=pmb->is; i<=pmb->ie; i++)
              e1(k,j,i)*=0.5;
          }
          for(int k=pmb->ks; k<=pmb->ke; k++) {
            for(int i=pmb->is+1; i<=pmb->ie; i++)
              e3(k,j,i)*=0.5;
          }
        }
        else if(pmb->block_size.nx2 > 1) {
          for(int i=pmb->is; i<=pmb->ie; i++)
            e1(pmb->ks,j,i)*=0.5, e1(pmb->ks+1,j,i)*=0.5;
          for(int i=pmb->is+1; i<=pmb->ie; i++)
            e3(pmb->ks,j,i)*=0.5;
        }
      }
      else if(nl>pmb->loc.level) { // finer; divide the overlapping EMFs by 2
        if(pmb->block_size.nx3 > 1) { // 3D
          int k=pmb->ks+pmb->block_size.nx3/2;
          for(int i=pmb->is; i<=pmb->ie; i++)
            e1(k,j,i)*=0.5;
        }
        if(pmb->block_size.nx2 > 1) { // 2D or 3D
          int i=pmb->is+pmb->block_size.nx1/2;
          for(int k=pmb->ks; k<=pmb->ke; k++)
            e3(k,j,i)*=0.5;
        }
      }
    }
    if(n==inner_x3 || n==outer_x3) {
      if(n==inner_x3) k=pmb->ks;
      else k=pmb->ke+1;
      nl=pmb->nblevel[2*n-8][1][1];
      if(nl==pmb->loc.level) { // same ; divide all the face EMFs by 2
        for(int j=pmb->js+1; j<=pmb->je; j++) {
          for(int i=pmb->is; i<=pmb->ie; i++)
            e1(k,j,i)*=0.5;
        }
        for(int j=pmb->js; j<=pmb->je; j++) {
          for(int i=pmb->is+1; i<=pmb->ie; i++)
            e2(k,j,i)*=0.5;
        }
      }
      else if(nl>pmb->loc.level) { // finer; divide the overlapping EMFs by 2
        // this is always 3D
        int j=pmb->js+pmb->block_size.nx2/2;
        for(int i=pmb->is; i<=pmb->ie; i++)
          e1(k,j,i)*=0.5;
        int i=pmb->is+pmb->block_size.nx1/2;
        for(int j=pmb->js; j<=pmb->je; j++)
          e2(k,j,i)*=0.5;
      }
    }
  }
  // edge
  for(int n=0; n<nedge_; n++) {
    if(nedge_fine_[n]==1) continue;
    Real div=1.0/(Real)nedge_fine_[n];
    if(n>=0 && n<4) {
      if((n&1)==0) i=pmb->is;
      else i=pmb->ie+1;
      if((n&2)==0) j=pmb->js;
      else j=pmb->je+1;
      for(int k=pmb->ks; k<=pmb->ke; k++)
        e3(k,j,i)*=div;
    }
    // x1x3 edge
    else if(n>=4 && n<8) {
      if((n&1)==0) i=pmb->is;
      else i=pmb->ie+1;
      if((n&2)==0) k=pmb->ks;
      else k=pmb->ke+1;
      for(int j=pmb->js; j<=pmb->je; j++)
        e2(k,j,i)*=div;
    }
    // x2x3 edge
    else if(n>=8 && n<12) {
      if((n&1)==0) j=pmb->js;
      else j=pmb->je+1;
      if((n&2)==0) k=pmb->ks;
      else k=pmb->ke+1;
      for(int i=pmb->is; i<=pmb->ie; i++)
        e1(k,j,i)*=div;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ReceiveEMFCorrection(int step)
//  \brief Receive and Apply the surace EMF to the coarse neighbor(s) if needed
bool BoundaryValues::ReceiveEMFCorrection(int step)
{
  MeshBlock *pmb=pmy_mblock_;
  bool flag=true;

  if(firsttime_[step]==true) {
    for(int n=0; n<pmb->nneighbor; n++) { // first correct the same level
      NeighborBlock& nb = pmb->neighbor[n];
      if(nb.type!=neighbor_face && nb.type!=neighbor_edge) break;
      if(nb.level!=pmb->loc.level) continue;
      if((nb.type==neighbor_face) || ((nb.type==neighbor_edge) && (edge_flag_[nb.eid]==true))) {
        if(emfcor_flag_[step][nb.bufid]==boundary_completed) continue;
        if(emfcor_flag_[step][nb.bufid]==boundary_waiting) {
          if(nb.rank==Globals::my_rank) {// on the same process
            flag=false;
            continue;
          }
#ifdef MPI_PARALLEL
          else { // MPI boundary
            int test;
            MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&test,MPI_STATUS_IGNORE);
            MPI_Test(&req_emfcor_recv_[step][nb.bufid],&test,MPI_STATUS_IGNORE);
            if(test==false) {
              flag=false;
              continue;
            }
            emfcor_flag_[step][nb.bufid] = boundary_arrived;
          }
#endif
        }
        // boundary arrived; apply EMF correction
        SetEMFBoundarySameLevel(emfcor_recv_[step][nb.bufid], nb);
        emfcor_flag_[step][nb.bufid] = boundary_completed;
      }
    }

    if(flag==false) return flag;

    ClearCoarseEMFBoundary();
    firsttime_[step]=false;
  }

  for(int n=0; n<pmb->nneighbor; n++) { // then from finer
    NeighborBlock& nb = pmb->neighbor[n];
    if(nb.type!=neighbor_face && nb.type!=neighbor_edge) break;
    if(nb.level!=pmb->loc.level+1) continue;
    if(emfcor_flag_[step][nb.bufid]==boundary_completed) continue;
    if(emfcor_flag_[step][nb.bufid]==boundary_waiting) {
      if(nb.rank==Globals::my_rank) {// on the same process
        flag=false;
        continue;
      }
#ifdef MPI_PARALLEL
      else { // MPI boundary
        int test;
        MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&test,MPI_STATUS_IGNORE);
        MPI_Test(&req_emfcor_recv_[step][nb.bufid],&test,MPI_STATUS_IGNORE);
        if(test==false) {
          flag=false;
          continue;
        }
        emfcor_flag_[step][nb.bufid] = boundary_arrived;
      }
#endif
    }
    // boundary arrived; apply EMF correction
    SetEMFBoundaryFromFiner(emfcor_recv_[step][nb.bufid], nb);
    emfcor_flag_[step][nb.bufid] = boundary_completed;
  }

  if(flag==true)
    AverageEMFBoundary();
  return flag;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ProlongateFieldBoundaries(InterfaceField &dst)
//  \brief Prolongate the fields in the ghost zones from the prolongation buffer
void BoundaryValues::ProlongateFieldBoundaries(InterfaceField &dst)
{
  MeshBlock *pmb=pmy_mblock_;
  Coordinates *pco=pmb->pcoord;
  Coordinates *pcrs=pmb->pcoarsec;
  int mox1, mox2, mox3;
  long int &lx1=pmb->loc.lx1;
  long int &lx2=pmb->loc.lx2;
  long int &lx3=pmb->loc.lx3;
  int &mylevel=pmb->loc.level;
  mox1=((int)(lx1&1L)<<1)-1;
  mox2=((int)(lx2&1L)<<1)-1;
  mox3=((int)(lx3&1L)<<1)-1;
  int cng=pmb->cnghost;

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    if(nb.level >= mylevel) continue;
    int mytype=std::abs(nb.ox1)+std::abs(nb.ox2)+std::abs(nb.ox3);
    // fill the required ghost-ghost zone
    int nis, nie, njs, nje, nks, nke;
    nis=std::max(nb.ox1-1,-1), nie=std::min(nb.ox1+1,1);
    if(pmb->block_size.nx2==1) njs=0, nje=0;
    else njs=std::max(nb.ox2-1,-1), nje=std::min(nb.ox2+1,1);
    if(pmb->block_size.nx3==1) nks=0, nke=0;
    else nks=std::max(nb.ox3-1,-1), nke=std::min(nb.ox3+1,1);
    for(int nk=nks; nk<=nke; nk++) {
      for(int nj=njs; nj<=nje; nj++) {
        for(int ni=nis; ni<=nie; ni++) {
          int ntype=std::abs(ni)+std::abs(nj)+std::abs(nk);
          if(ntype==0) continue; // skip myself
          if(pmb->nblevel[nk+1][nj+1][ni+1]!=mylevel
          && pmb->nblevel[nk+1][nj+1][ni+1]!=-1)
            continue; // physical boundary will also be restricted
          if(ntype>mytype) {
            if(pmb->block_size.nx3 > 1) // 3D
              if(((mox1==ni)+(mox2==nj)+(mox3==nk)) != ntype) continue;
            else if(pmb->block_size.nx2 > 1) // 2D
              if(((mox1==ni)+(mox2==nj)) != ntype) continue;
          }

          // this neighbor block is on the same level
          // and needs to be restricted for prolongation
          int ris, rie, rjs, rje, rks, rke;
          if(ni==0) {
            ris=pmb->cis, rie=pmb->cie;
            if(nb.ox1==1) ris=pmb->cie;
            else if(nb.ox1==-1) rie=pmb->cis;
          }
          else if(ni== 1) ris=pmb->cie+1, rie=pmb->cie+1;
          else if(ni==-1) ris=pmb->cis-1, rie=pmb->cis-1;
          if(nj==0) {
            rjs=pmb->cjs, rje=pmb->cje;
            if(nb.ox2==1) rjs=pmb->cje;
            else if(nb.ox2==-1) rje=pmb->cjs;
          }
          else if(nj== 1) rjs=pmb->cje+1, rje=pmb->cje+1;
          else if(nj==-1) rjs=pmb->cjs-1, rje=pmb->cjs-1;
          if(nk==0) {
            rks=pmb->cks, rke=pmb->cke;
            if(nb.ox3==1) rks=pmb->cke;
            else if(nb.ox3==-1) rke=pmb->cks;
          }
          else if(nk== 1) rks=pmb->cke+1, rke=pmb->cke+1;
          else if(nk==-1) rks=pmb->cks-1, rke=pmb->cks-1;

          int rs=ris, re=rie+1;
          if(rs==pmb->cis && pmb->nblevel[nk+1][nj+1][ni]<mylevel) rs++;
          if(re==pmb->cie+1 && pmb->nblevel[nk+1][nj+1][ni+2]<mylevel) re--;
          RestrictFieldX1(dst.x1f, rs, re, rjs, rje, rks, rke);
          if(pmb->block_size.nx2 > 1) {
            rs=rjs, re=rje+1;
            if(rs==pmb->cjs && pmb->nblevel[nk+1][nj][ni+1]<mylevel) rs++;
            if(re==pmb->cje+1 && pmb->nblevel[nk+1][nj+2][ni+1]<mylevel) re--;
            RestrictFieldX2(dst.x2f, ris, rie, rs, re, rks, rke);
          }
          else 
            RestrictFieldX2(dst.x2f, ris, rie, rjs, rje, rks, rke);
          if(pmb->block_size.nx3 > 1) {
            rs=rks, re=rke+1;
            if(rs==pmb->cks && pmb->nblevel[nk][nj+1][ni+1]<mylevel) rs++;
            if(re==pmb->cke+1 && pmb->nblevel[nk+2][nj+1][ni+1]<mylevel) re--;
            RestrictFieldX3(dst.x3f, ris, rie, rjs, rje, rs, re);
          }
          else
            RestrictFieldX3(dst.x3f, ris, rie, rjs, rje, rks, rke);
        }
      }
    }
    // now that the ghost-ghost zones are filled
    // reconstruct finer magnetic fields using the Li & Li 2004 method
    // caculate the loop indexes of the *cells* to be refined.
    int cn = (NGHOST+1)/2;
    int si, ei, sj, ej, sk, ek;
    if(nb.ox1==0) {
      si=pmb->cis, ei=pmb->cie;
      if((lx1&1L)==0L) ei++;
      else             si--;
    }
    else if(nb.ox1>0) si=pmb->cie+1,  ei=pmb->cie+cn;
    else              si=pmb->cis-cn, ei=pmb->cis-1;
    if(nb.ox2==0) {
      sj=pmb->cjs, ej=pmb->cje;
      if(pmb->block_size.nx2 > 1) {
        if((lx2&1L)==0L) ej++;
        else             sj--;
      }
    }
    else if(nb.ox2>0) sj=pmb->cje+1,  ej=pmb->cje+cn;
    else              sj=pmb->cjs-cn, ej=pmb->cjs-1;
    if(nb.ox3==0) {
      sk=pmb->cks, ek=pmb->cke;
      if(pmb->block_size.nx3 > 1) {
        if((lx3&1L)==0L) ek++;
        else             sk--;
      }
    }
    else if(nb.ox3>0) sk=pmb->cke+1,  ek=pmb->cke+cn;
    else              sk=pmb->cks-cn, ek=pmb->cks-1;

    if(pmb->block_size.nx3 > 1) { // 3D
      int kl=sk, ku=ek+1;
      if((nb.ox3>=0) && (pmb->nblevel[nb.ox3  ][nb.ox2+1][nb.ox1+1]>=mylevel)) kl++;
      if((nb.ox3<=0) && (pmb->nblevel[nb.ox3+2][nb.ox2+1][nb.ox1+1]>=mylevel)) ku--;
      int jl=sj, ju=ej+1;
      if((nb.ox2>=0) && (pmb->nblevel[nb.ox3+1][nb.ox2  ][nb.ox1+1]>=mylevel)) jl++;
      if((nb.ox2<=0) && (pmb->nblevel[nb.ox3+1][nb.ox2+2][nb.ox1+1]>=mylevel)) ju--;
      int il=si, iu=ei+1;
      if((nb.ox1>=0) && (pmb->nblevel[nb.ox3+1][nb.ox2+1][nb.ox1  ]>=mylevel)) il++;
      if((nb.ox1<=0) && (pmb->nblevel[nb.ox3+1][nb.ox2+1][nb.ox1+2]>=mylevel)) iu--;
      // step 1. calculate x3 outer surface fields and slopes
      for(int k=sk; k<=ek+1; k++) {
        int fk=(k-pmb->cks)*2+pmb->ks;
        for(int j=sj; j<=ej; j++) {
          int fj=(j-pmb->cjs)*2+pmb->js;
          Real& x2m = pcrs->x2s3(j-1);
          Real& x2c = pcrs->x2s3(j);
          Real& x2p = pcrs->x2s3(j+1);
          Real dx2m = x2c - x2m;
          Real dx2p = x2p - x2c;
          Real& fx2m = pco->x2s3(fj);
          Real& fx2p = pco->x2s3(fj+1);
          for(int i=si; i<=ei; i++) {
            int fi=(i-pmb->cis)*2+pmb->is;
            Real& x1m = pcrs->x1s3(i-1);
            Real& x1c = pcrs->x1s3(i);
            Real& x1p = pcrs->x1s3(i+1);
            Real dx1m = x1c - x1m;
            Real dx1p = x1p - x1c;
            Real& fx1m = pco->x1s3(fi);
            Real& fx1p = pco->x1s3(fi+1);
            Real ccval=coarse_b_.x3f(k,j,i);

            Real gx1m = (ccval-coarse_b_.x3f(k,j,i-1))/dx1m;
            Real gx1p = (coarse_b_.x3f(k,j,i+1)-ccval)/dx1p;
            Real gx1c = 0.5*(SIGN(gx1m)+SIGN(gx1p))*std::min(std::abs(gx1m),std::abs(gx1p));
            Real gx2m = (ccval-coarse_b_.x3f(k,j-1,i))/dx2m;
            Real gx2p = (coarse_b_.x3f(k,j+1,i)-ccval)/dx2p;
            Real gx2c = 0.5*(SIGN(gx2m)+SIGN(gx2p))*std::min(std::abs(gx2m),std::abs(gx2p));

            if(k>=kl && k<=ku) {
              dst.x3f(fk,fj  ,fi  )=ccval-gx1c*(x1c-fx1m)-gx2c*(x2c-fx2m);
              dst.x3f(fk,fj  ,fi+1)=ccval+gx1c*(fx1p-x1c)-gx2c*(x2c-fx2m);
              dst.x3f(fk,fj+1,fi  )=ccval-gx1c*(x1c-fx1m)+gx2c*(fx2p-x2c);
              dst.x3f(fk,fj+1,fi+1)=ccval+gx1c*(fx1p-x1c)+gx2c*(fx2p-x2c);
            }
          }
        }
      }
      // step 2. calculate x2 outer surface fields and slopes
      for(int k=sk; k<=ek; k++) {
        int fk=(k-pmb->cks)*2+pmb->ks;
        Real& x3m = pcrs->x3s2(k-1);
        Real& x3c = pcrs->x3s2(k);
        Real& x3p = pcrs->x3s2(k+1);
        Real dx3m = x3c - x3m;
        Real dx3p = x3p - x3c;
        Real& fx3m = pco->x3s2(fk);
        Real& fx3p = pco->x3s2(fk+1);
        for(int j=sj; j<=ej+1; j++) {
          int fj=(j-pmb->cjs)*2+pmb->js;
          for(int i=si; i<=ei; i++) {
            int fi=(i-pmb->cis)*2+pmb->is;
            Real& x1m = pcrs->x1s2(i-1);
            Real& x1c = pcrs->x1s2(i);
            Real& x1p = pcrs->x1s2(i+1);
            Real dx1m = x1c - x1m;
            Real dx1p = x1p - x1c;
            Real& fx1m = pco->x1s2(fi);
            Real& fx1p = pco->x1s2(fi+1);
            Real ccval=coarse_b_.x2f(k,j,i);

            Real gx1m = (ccval-coarse_b_.x2f(k,j,i-1))/dx1m;
            Real gx1p = (coarse_b_.x2f(k,j,i+1)-ccval)/dx1p;
            Real gx1c = 0.5*(SIGN(gx1m)+SIGN(gx1p))*std::min(std::abs(gx1m),std::abs(gx1p));
            Real gx3m = (ccval-coarse_b_.x2f(k-1,j,i))/dx3m;
            Real gx3p = (coarse_b_.x2f(k+1,j,i)-ccval)/dx3p;
            Real gx3c = 0.5*(SIGN(gx3m)+SIGN(gx3p))*std::min(std::abs(gx3m),std::abs(gx3p));

            if(j>=jl && j<=ju) {
              dst.x2f(fk  ,fj,fi  )=ccval-gx1c*(x1c-fx1m)-gx3c*(x3c-fx3m);
              dst.x2f(fk  ,fj,fi+1)=ccval+gx1c*(fx1p-x1c)-gx3c*(x3c-fx3m);
              dst.x2f(fk+1,fj,fi  )=ccval-gx1c*(x1c-fx1m)+gx3c*(fx3p-x3c);
              dst.x2f(fk+1,fj,fi+1)=ccval+gx1c*(fx1p-x1c)+gx3c*(fx3p-x3c);
            }
          }
        }
      }
      // step 3. calculate x1 outer surface fields and slopes
      for(int k=sk; k<=ek; k++) {
        int fk=(k-pmb->cks)*2+pmb->ks;
        Real& x3m = pcrs->x3s1(k-1);
        Real& x3c = pcrs->x3s1(k);
        Real& x3p = pcrs->x3s1(k+1);
        Real dx3m = x3c - x3m;
        Real dx3p = x3p - x3c;
        Real& fx3m = pco->x3s1(fk);
        Real& fx3p = pco->x3s1(fk+1);
        for(int j=sj; j<=ej; j++) {
          int fj=(j-pmb->cjs)*2+pmb->js;
          Real& x2m = pcrs->x2s1(j-1);
          Real& x2c = pcrs->x2s1(j);
          Real& x2p = pcrs->x2s1(j+1);
          Real dx2m = x2c - x2m;
          Real dx2p = x2p - x2c;
          Real& fx2m = pco->x2s1(fj);
          Real& fx2p = pco->x2s1(fj+1);
          for(int i=si; i<=ei+1; i++) {
            int fi=(i-pmb->cis)*2+pmb->is;
            Real ccval=coarse_b_.x1f(k,j,i);

            Real gx2m = (ccval-coarse_b_.x1f(k,j-1,i))/dx2m;
            Real gx2p = (coarse_b_.x1f(k,j+1,i)-ccval)/dx2p;
            Real gx2c = 0.5*(SIGN(gx2m)+SIGN(gx2p))*std::min(std::abs(gx2m),std::abs(gx2p));
            Real gx3m = (ccval-coarse_b_.x1f(k-1,j,i))/dx3m;
            Real gx3p = (coarse_b_.x1f(k+1,j,i)-ccval)/dx3p;
            Real gx3c = 0.5*(SIGN(gx3m)+SIGN(gx3p))*std::min(std::abs(gx3m),std::abs(gx3p));

            if(i>=il && i<=iu) {
              dst.x1f(fk  ,fj  ,fi)=ccval-gx2c*(x2c-fx2m)-gx3c*(x3c-fx3m);
              dst.x1f(fk  ,fj+1,fi)=ccval+gx2c*(fx2p-x2c)-gx3c*(x3c-fx3m);
              dst.x1f(fk+1,fj  ,fi)=ccval-gx2c*(x2c-fx2m)+gx3c*(fx3p-x3c);
              dst.x1f(fk+1,fj+1,fi)=ccval+gx2c*(fx2p-x2c)+gx3c*(fx3p-x3c);
            }
          }
        }
      }
      // step 4. calculate the internal finer fields using the Toth & Roe method
      int fsi=(si-pmb->cis)*2+pmb->is, fei=(ei-pmb->cis)*2+pmb->is+1;
      for(int k=sk; k<=ek; k++) {
        int fk=(k-pmb->cks)*2+pmb->ks;
        for(int j=sj; j<=ej; j++) {
          int fj=(j-pmb->cjs)*2+pmb->js;
          for(int i=si; i<=ei; i++) {
            int fi=(i-pmb->cis)*2+pmb->is;
            Real Uxx = 0.0, Vyy = 0.0, Wzz = 0.0;
            Real Uxyz = 0.0, Vxyz = 0.0, Wxyz = 0.0;
#pragma unroll
            for(int jj=0; jj<2; jj++){
              int js=2*jj-1, fjj=fj+jj, fjp=fj+2*jj;
#pragma unroll
              for(int ii=0; ii<2; ii++){
                int is=2*ii-1, fii=fi+ii, fip=fi+2*ii;
                Uxx += is*(js*(dst.x2f(fk  ,fjp,fii)*pco->GetFace2Area(fk,fjp,fii) + dst.x2f(fk+1,fjp,fii)*pco->GetFace2Area(fk+1,fjp,fii))
                             +(dst.x3f(fk+2,fjj,fii)*pco->GetFace3Area(fk+2,fjj,fii) - dst.x3f(fk  ,fjj,fii)*pco->GetFace3Area(fk,fjj,fii)));
                Vyy += js*(   (dst.x3f(fk+2,fjj,fii)*pco->GetFace3Area(fk+2,fjj,fii)  - dst.x3f(fk  ,fjj,fii)*pco->GetFace3Area(fk,fjj,fii))
                          +is*(dst.x1f(fk  ,fjj,fip)*pco->GetFace1Area(fk,fjj,fip) + dst.x1f(fk+1,fjj,fip)*pco->GetFace1Area(fk+1,fjj,fip)));
                Wzz +=     is*(dst.x1f(fk+1,fjj,fip)*pco->GetFace1Area(fk+1,fjj,fip) - dst.x1f(fk  ,fjj,fip)*pco->GetFace1Area(fk,fjj,fip))
                          +js*(dst.x2f(fk+1,fjp,fii)*pco->GetFace2Area(fk+1,fjp,fii) - dst.x2f(fk  ,fjp,fii)*pco->GetFace2Area(fk,fjp,fii));
                Uxyz += js*js*(dst.x1f(fk+1,fjj,fip)*pco->GetFace1Area(fk+1,fjj,fip) - dst.x1f(fk  ,fjj,fip)*pco->GetFace1Area(fk,fjj,fip));
                Vxyz += is*js*(dst.x2f(fk+1,fjp,fii)*pco->GetFace2Area(fk+1,fjp,fii) - dst.x2f(fk  ,fjp,fii)*pco->GetFace2Area(fk,fjp,fii));
                Wxyz += is*js*(dst.x3f(fk+2,fjj,fii)*pco->GetFace3Area(fk+2,fjj,fii) - dst.x3f(fk  ,fjj,fii)*pco->GetFace3Area(fk,fjj,fii));
              }
            }
	    Real Sdx1=SQR(pco->dx1f(fi)+pco->dx1f(fi+1));
            Real Sdx2=SQR(pco->GetEdge2Length(fk+1,fj,fi+1)+pco->GetEdge2Length(fk+1,fj+1,fi+1));
            Real Sdx3=SQR(pco->GetEdge3Length(fk,fj+1,fi+1)+pco->GetEdge3Length(fk+1,fj+1,fi+1));
	    Uxx *= 0.125; Vyy *= 0.125; Wzz *= 0.125;
            Uxyz *= 0.125/(Sdx2+Sdx3); Vxyz *= 0.125/(Sdx1+Sdx3); Wxyz *= 0.125/(Sdx1+Sdx2);
            dst.x1f(fk  ,fj  ,fi+1)=(0.5*(dst.x1f(fk  ,fj  ,fi  )*pco->GetFace1Area(fk  ,fj  ,fi  )+dst.x1f(fk  ,fj  ,fi+2)*pco->GetFace1Area(fk  ,fj  ,fi+2))
                                   + Uxx - Sdx3*Vxyz - Sdx2*Wxyz)/pco->GetFace1Area(fk  ,fj  ,fi+1);
            dst.x1f(fk  ,fj+1,fi+1)=(0.5*(dst.x1f(fk  ,fj+1,fi  )*pco->GetFace1Area(fk  ,fj+1,fi  )+dst.x1f(fk  ,fj+1,fi+2)*pco->GetFace1Area(fk  ,fj+1,fi+2))
                                   + Uxx - Sdx3*Vxyz + Sdx2*Wxyz)/pco->GetFace1Area(fk  ,fj+1,fi+1);
            dst.x1f(fk+1,fj  ,fi+1)=(0.5*(dst.x1f(fk+1,fj  ,fi  )*pco->GetFace1Area(fk+1,fj  ,fi  )+dst.x1f(fk+1,fj  ,fi+2)*pco->GetFace1Area(fk+1,fj  ,fi+2))
                                   + Uxx + Sdx3*Vxyz - Sdx2*Wxyz)/pco->GetFace1Area(fk+1,fj  ,fi+1);
            dst.x1f(fk+1,fj+1,fi+1)=(0.5*(dst.x1f(fk+1,fj+1,fi  )*pco->GetFace1Area(fk+1,fj+1,fi  )+dst.x1f(fk+1,fj+1,fi+2)*pco->GetFace1Area(fk+1,fj+1,fi+2))
                                   + Uxx + Sdx3*Vxyz + Sdx2*Wxyz)/pco->GetFace1Area(fk+1,fj+1,fi+1);

	    dst.x2f(fk  ,fj+1,fi  )=(0.5*(dst.x2f(fk  ,fj  ,fi  )*pco->GetFace2Area(fk  ,fj  ,fi  )+dst.x2f(fk  ,fj+2,fi  )*pco->GetFace2Area(fk  ,fj+2,fi  ))
                                   + Vyy - Sdx3*Uxyz - Sdx1*Wxyz)/pco->GetFace2Area(fk  ,fj+1,fi  );
            dst.x2f(fk  ,fj+1,fi+1)=(0.5*(dst.x2f(fk  ,fj  ,fi+1)*pco->GetFace2Area(fk  ,fj  ,fi+1)+dst.x2f(fk  ,fj+2,fi+1)*pco->GetFace2Area(fk  ,fj+2,fi+1))
                                   + Vyy - Sdx3*Uxyz + Sdx1*Wxyz)/pco->GetFace2Area(fk  ,fj+1,fi+1);
            dst.x2f(fk+1,fj+1,fi  )=(0.5*(dst.x2f(fk+1,fj  ,fi  )*pco->GetFace2Area(fk+1,fj  ,fi  )+dst.x2f(fk+1,fj+2,fi  )*pco->GetFace2Area(fk+1,fj+2,fi  ))
                                   + Vyy + Sdx3*Uxyz - Sdx1*Wxyz)/pco->GetFace2Area(fk+1,fj+1,fi  );
            dst.x2f(fk+1,fj+1,fi+1)=(0.5*(dst.x2f(fk+1,fj  ,fi+1)*pco->GetFace2Area(fk+1,fj  ,fi+1)+dst.x2f(fk+1,fj+2,fi+1)*pco->GetFace2Area(fk+1,fj+2,fi+1))
                                   + Vyy + Sdx3*Uxyz + Sdx1*Wxyz)/pco->GetFace2Area(fk+1,fj+1,fi+1);

            dst.x3f(fk+1,fj  ,fi  )=(0.5*(dst.x3f(fk+2,fj  ,fi  )*pco->GetFace3Area(fk+2,fj  ,fi  )+dst.x3f(fk  ,fj  ,fi  )*pco->GetFace3Area(fk  ,fj  ,fi  ))
                                   + Wzz - Sdx2*Uxyz - Sdx1*Vxyz)/pco->GetFace3Area(fk+1,fj  ,fi  );
            dst.x3f(fk+1,fj  ,fi+1)=(0.5*(dst.x3f(fk+2,fj  ,fi+1)*pco->GetFace3Area(fk+2,fj  ,fi+1)+dst.x3f(fk  ,fj  ,fi+1)*pco->GetFace3Area(fk  ,fj  ,fi+1))
                                   + Wzz - Sdx2*Uxyz + Sdx1*Vxyz)/pco->GetFace3Area(fk+1,fj  ,fi+1);
            dst.x3f(fk+1,fj+1,fi  )=(0.5*(dst.x3f(fk+2,fj+1,fi  )*pco->GetFace3Area(fk+2,fj+1,fi  )+dst.x3f(fk  ,fj+1,fi  )*pco->GetFace3Area(fk  ,fj+1,fi  ))
                                   + Wzz + Sdx2*Uxyz - Sdx1*Vxyz)/pco->GetFace3Area(fk+1,fj+1,fi  );
            dst.x3f(fk+1,fj+1,fi+1)=(0.5*(dst.x3f(fk+2,fj+1,fi+1)*pco->GetFace3Area(fk+2,fj+1,fi+1)+dst.x3f(fk  ,fj+1,fi+1)*pco->GetFace3Area(fk  ,fj+1,fi+1))
                                   + Wzz + Sdx2*Uxyz + Sdx1*Vxyz)/pco->GetFace3Area(fk+1,fj+1,fi+1);
          }
        }
      }
    }
    else if(pmb->block_size.nx2 > 1) { // 2D
      int k=pmb->cks, fk=pmb->ks;
      int jl=sj, ju=ej+1;
      if((nb.ox2>=0) && (pmb->nblevel[1][nb.ox2  ][nb.ox1+1]>=mylevel)) jl++;
      if((nb.ox2<=0) && (pmb->nblevel[1][nb.ox2+2][nb.ox1+1]>=mylevel)) ju--;
      int il=si, iu=ei+1;
      if((nb.ox1>=0) && (pmb->nblevel[1][nb.ox2+1][nb.ox1  ]>=mylevel)) il++;
      if((nb.ox1<=0) && (pmb->nblevel[1][nb.ox2+1][nb.ox1+2]>=mylevel)) iu--;
      // step 1. calculate x2 outer surface fields and slopes
      for(int j=sj; j<=ej+1; j++) {
        int fj=(j-pmb->cjs)*2+pmb->js;
        for(int i=si; i<=ei; i++) {
          int fi=(i-pmb->cis)*2+pmb->is;
          Real& x1m = pcrs->x1s2(i-1);
          Real& x1c = pcrs->x1s2(i);
          Real& x1p = pcrs->x1s2(i+1);
          Real& fx1m = pco->x1s2(fi);
          Real& fx1p = pco->x1s2(fi+1);
          Real ccval=coarse_b_.x2f(k,j,i);

          Real gx1m = (ccval-coarse_b_.x2f(k,j,i-1))/(x1c - x1m);
          Real gx1p = (coarse_b_.x2f(k,j,i+1)-ccval)/(x1p - x1c);
          Real gx1c = 0.5*(SIGN(gx1m)+SIGN(gx1p))*std::min(std::abs(gx1m),std::abs(gx1p));

          if(j>=jl && j<=ju) {
            dst.x2f(fk,fj,fi  )=ccval-gx1c*(x1c-fx1m);
            dst.x2f(fk,fj,fi+1)=ccval+gx1c*(fx1p-x1c);
          }
        }
      }
      // step 2. calculate x1 outer surface fields and slopes
      for(int j=sj; j<=ej; j++) {
        int fj=(j-pmb->cjs)*2+pmb->js;
        Real& x2m = pcrs->x2s1(j-1);
        Real& x2c = pcrs->x2s1(j);
        Real& x2p = pcrs->x2s1(j+1);
        Real dx2m = x2c - x2m;
        Real dx2p = x2p - x2c;
        Real& fx2m = pco->x2s1(fj);
        Real& fx2p = pco->x2s1(fj+1);
        for(int i=si; i<=ei+1; i++) {
          int fi=(i-pmb->cis)*2+pmb->is;
          Real ccval=coarse_b_.x1f(k,j,i);

          Real gx2m = (ccval-coarse_b_.x1f(k,j-1,i))/dx2m;
          Real gx2p = (coarse_b_.x1f(k,j+1,i)-ccval)/dx2p;
          Real gx2c = 0.5*(SIGN(gx2m)+SIGN(gx2p))*std::min(std::abs(gx2m),std::abs(gx2p));

          if(i>=il && i<=iu) {
            dst.x1f(fk,fj  ,fi)=ccval-gx2c*(x2c-fx2m);
            dst.x1f(fk,fj+1,fi)=ccval+gx2c*(fx2p-x2c);
          }
        }
      }
      int fsi=(si-pmb->cis)*2+pmb->is, fei=(ei-pmb->cis)*2+pmb->is+1;
      // step 3. calculate the internal finer fields using the Toth & Roe method
      for(int j=sj; j<=ej; j++) {
        int fj=(j-pmb->cjs)*2+pmb->js;
        for(int i=si; i<=ei; i++) {
          int fi=(i-pmb->cis)*2+pmb->is;
          Real tmp1=0.25*(dst.x2f(fk,fj+2,fi+1)*pco->GetFace2Area(fk,fj+2,fi+1)
                         -dst.x2f(fk,fj,fi+1)*pco->GetFace2Area(fk,fj,fi+1)
                         -dst.x2f(fk,fj+2,fi)*pco->GetFace2Area(fk,fj+2,fi)
                         +dst.x2f(fk,fj,fi)*pco->GetFace2Area(fk,fj,fi));
          Real tmp2=0.25*(dst.x1f(fk,fj,fi)*pco->GetFace1Area(fk,fj,fi)
                         -dst.x1f(fk,fj,fi+2)*pco->GetFace1Area(fk,fj,fi+2)
                         -dst.x1f(fk,fj+1,fi)*pco->GetFace1Area(fk,fj+1,fi)
                         +dst.x1f(fk,fj+1,fi+2)*pco->GetFace1Area(fk,fj+1,fi+2));
          dst.x1f(fk,fj  ,fi+1)=(0.5*(dst.x1f(fk,fj,fi)*pco->GetFace1Area(fk,fj,fi)
                                     +dst.x1f(fk,fj,fi+2)*pco->GetFace1Area(fk,fj,fi+2))+tmp1)/pco->GetFace1Area(fk,fj  ,fi+1);
          dst.x1f(fk,fj+1,fi+1)=(0.5*(dst.x1f(fk,fj+1,fi)*pco->GetFace1Area(fk,fj+1,fi)
                                     +dst.x1f(fk,fj+1,fi+2)*pco->GetFace1Area(fk,fj+1,fi+2))+tmp1)/pco->GetFace1Area(fk,fj+1,fi+1);
          dst.x2f(fk,fj+1,fi  )=(0.5*(dst.x2f(fk,fj,fi)*pco->GetFace2Area(fk,fj,fi)
                                     +dst.x2f(fk,fj+2,fi)*pco->GetFace2Area(fk,fj+2,fi))+tmp2)/pco->GetFace2Area(fk,fj+1,fi  );
          dst.x2f(fk,fj+1,fi+1)=(0.5*(dst.x2f(fk,fj,fi+1)*pco->GetFace2Area(fk,fj,fi+1)
                                     +dst.x2f(fk,fj+2,fi+1)*pco->GetFace2Area(fk,fj+2,fi+1))+tmp2)/pco->GetFace2Area(fk,fj+1,fi+1);
        }
      }
      // step 4. calculate the finer x3 fields (independent from x1 and x2)
      for(int j=sj; j<=ej; j++) {
        int fj=(j-pmb->cjs)*2+pmb->js;
        Real& x2m = pcrs->x2s3(j-1);
        Real& x2c = pcrs->x2s3(j);
        Real& x2p = pcrs->x2s3(j+1);
        Real dx2m = x2c - x2m;
        Real dx2p = x2p - x2c;
        Real& fx2m = pco->x2s3(fj);
        Real& fx2p = pco->x2s3(fj+1);
        Real dx2fm= x2c-fx2m;
        Real dx2fp= fx2p-x2c;
        for(int i=si; i<=ei; i++) {
          int fi=(i-pmb->cis)*2+pmb->is;
          Real& x1m = pcrs->x1s3(i-1);
          Real& x1c = pcrs->x1s3(i);
          Real& x1p = pcrs->x1s3(i+1);
          Real dx1m = x1c - x1m;
          Real dx1p = x1p - x1c;
          Real& fx1m = pco->x1s3(fi);
          Real& fx1p = pco->x1s3(fi+1);
          Real dx1fm= x1c-fx1m;
          Real dx1fp= fx1p-x1c;
          Real ccval=coarse_b_.x3f(k,j,i);

          // calculate 2D gradients using the minmod limiter
          Real gx1m = (ccval-coarse_b_.x3f(k,j,i-1))/dx1m;
          Real gx1p = (coarse_b_.x3f(k,j,i+1)-ccval)/dx1p;
          Real gx1c = 0.5*(SIGN(gx1m)+SIGN(gx1p))*std::min(std::abs(gx1m),std::abs(gx1p));
          Real gx2m = (ccval-coarse_b_.x3f(k,j-1,i))/dx2m;
          Real gx2p = (coarse_b_.x3f(k,j+1,i)-ccval)/dx2p;
          Real gx2c = 0.5*(SIGN(gx2m)+SIGN(gx2p))*std::min(std::abs(gx2m),std::abs(gx2p));

          // interpolate on to the finer grid
          dst.x3f(fk,fj  ,fi  )=dst.x3f(fk+1,fj  ,fi  )=ccval-gx1c*dx1fm-gx2c*dx2fm;
          dst.x3f(fk,fj  ,fi+1)=dst.x3f(fk+1,fj  ,fi+1)=ccval+gx1c*dx1fp-gx2c*dx2fm;
          dst.x3f(fk,fj+1,fi  )=dst.x3f(fk+1,fj+1,fi  )=ccval-gx1c*dx1fm+gx2c*dx2fp;
          dst.x3f(fk,fj+1,fi+1)=dst.x3f(fk+1,fj+1,fi+1)=ccval+gx1c*dx1fp+gx2c*dx2fp;
        }
      }
    }
    else { // 1D
      // bx1 - no prolongation, just divergence condition
      if(nb.ox1==1) {
        int fi=(si-pmb->cis)*2+pmb->is;
        Real ph=pco->GetFace1Area(0,0,fi)*dst.x1f(0,0,fi);
        dst.x1f(0,0,fi+1)=ph/pco->GetFace1Area(0,0,fi+1);
        dst.x1f(0,0,fi+2)=ph/pco->GetFace1Area(0,0,fi+2);
      }
      else if(nb.ox1==-1) {
        int fi=(ei+1-pmb->cis)*2+pmb->is;
        Real ph=pco->GetFace1Area(0,0,fi)*dst.x1f(0,0,fi);
        dst.x1f(0,0,fi-1)=ph/pco->GetFace1Area(0,0,fi-1);
        dst.x1f(0,0,fi-2)=ph/pco->GetFace1Area(0,0,fi-2);
      }
      int fi=(si-pmb->cis)*2+pmb->is;
      // bx2 and bx3 - interpolation
      Real gxm = (coarse_b_.x2f(0,0,si)-coarse_b_.x2f(0,0,si-1))
                 /(pcrs->x1s2(si)-pcrs->x1s2(si-1));
      Real gxp = (coarse_b_.x2f(0,0,si+1)-coarse_b_.x2f(0,0,si))
                 /(pcrs->x1s2(si+1)-pcrs->x1s2(si));
      Real gxc = 0.5*(SIGN(gxm)+SIGN(gxp))*std::min(std::abs(gxm),std::abs(gxp));
      dst.x2f(0,0,fi  )=dst.x2f(0,1,fi  )
                       =coarse_b_.x2f(0,0,si)-gxc*(pcrs->x1s2(si)-pco->x1s2(fi));
      dst.x2f(0,0,fi+1)=dst.x2f(0,1,fi+1)
                       =coarse_b_.x2f(0,0,si)+gxc*(pco->x1s2(fi+1)-pcrs->x1s2(si));
      gxm = (coarse_b_.x3f(0,0,si)-coarse_b_.x3f(0,0,si-1))
            /(pcrs->x1s3(si)-pcrs->x1s3(si-1));
      gxp = (coarse_b_.x3f(0,0,si+1)-coarse_b_.x3f(0,0,si))
            /(pcrs->x1s3(si+1)-pcrs->x1s3(si));
      gxc = 0.5*(SIGN(gxm)+SIGN(gxp))*std::min(std::abs(gxm),std::abs(gxp));
      dst.x3f(0,0,fi  )=dst.x3f(1,0,fi  )
                       =coarse_b_.x3f(0,0,si)-gxc*(pcrs->x1s3(si)-pco->x1s3(fi));
      dst.x3f(0,0,fi+1)=dst.x3f(1,0,fi+1)
                       =coarse_b_.x3f(0,0,si)+gxc*(pco->x1s3(fi+1)-pcrs->x1s3(si));
    }
  }
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ClearBoundaryForInit(void)
//  \brief clean up the boundary flags for initialization
void BoundaryValues::ClearBoundaryForInit(void)
{
  MeshBlock *pmb=pmy_mblock_;

  for(int n=0;n<pmb->nneighbor;n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    hydro_flag_[0][nb.bufid] = boundary_waiting;
    if (MAGNETIC_FIELDS_ENABLED)
      field_flag_[0][nb.bufid] = boundary_waiting;
#ifdef MPI_PARALLEL
    if(nb.rank!=Globals::my_rank) {
      MPI_Wait(&req_hydro_send_[0][nb.bufid],MPI_STATUS_IGNORE); // Wait for Isend
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
      NeighborBlock& nb = pmb->neighbor[n];
      hydro_flag_[l][nb.bufid] = boundary_waiting;
      if(nb.type==neighbor_face)
        flcor_flag_[l][nb.fid][nb.fi2][nb.fi1] = boundary_waiting;
      if (MAGNETIC_FIELDS_ENABLED) {
        field_flag_[l][nb.bufid] = boundary_waiting;
        if((nb.type==neighbor_face) || (nb.type==neighbor_edge))
          emfcor_flag_[l][nb.bufid] = boundary_waiting;
      }
#ifdef MPI_PARALLEL
      if(nb.rank!=Globals::my_rank) {
        MPI_Wait(&req_hydro_send_[l][nb.bufid],MPI_STATUS_IGNORE); // Wait for Isend
        if(nb.type==neighbor_face && nb.level<pmb->loc.level)
          MPI_Wait(&req_flcor_send_[l][nb.fid],MPI_STATUS_IGNORE); // Wait for Isend
        if (MAGNETIC_FIELDS_ENABLED) {
          MPI_Wait(&req_field_send_[l][nb.bufid],MPI_STATUS_IGNORE); // Wait for Isend
          if(nb.type==neighbor_face || nb.type==neighbor_edge) {
            if(nb.level < pmb->loc.level)
              MPI_Wait(&req_emfcor_send_[l][nb.bufid],MPI_STATUS_IGNORE); // Wait for Isend
            else if((nb.level==pmb->loc.level) && ((nb.type==neighbor_face)
                || ((nb.type==neighbor_edge) && (edge_flag_[nb.eid]==true))))
              MPI_Wait(&req_emfcor_send_[l][nb.bufid],MPI_STATUS_IGNORE); // Wait for Isend
          }
        }
      }
#endif
    }
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::HydroPhysicalBoundaries(AthenaArray<Real> &dst)
//  \brief Apply physical boundary conditions for hydro
void BoundaryValues::HydroPhysicalBoundaries(AthenaArray<Real> &dst)
{
  MeshBlock *pmb=pmy_mblock_;
  Coordinates *pco=pmb->pcoord;
  int bis=pmb->is, bie=pmb->ie, bjs=pmb->js, bje=pmb->je, bks=pmb->ks, bke=pmb->ke;

  if(pmb->pmy_mesh->face_only==false) { // extend the ghost zone
    bis=pmb->is-NGHOST;
    bie=pmb->ie+NGHOST;
    if(HydroBoundary_[inner_x2]==NULL && pmb->block_size.nx2>1) bjs=pmb->js-NGHOST;
    if(HydroBoundary_[outer_x2]==NULL && pmb->block_size.nx2>1) bje=pmb->je+NGHOST;
    if(HydroBoundary_[inner_x3]==NULL && pmb->block_size.nx3>1) bks=pmb->ks-NGHOST;
    if(HydroBoundary_[outer_x3]==NULL && pmb->block_size.nx3>1) bke=pmb->ke+NGHOST;
  }

  if(HydroBoundary_[inner_x1]!=NULL)
    HydroBoundary_[inner_x1](pmb, pco, dst, pmb->is, pmb->ie, bjs, bje, bks, bke);
  if(HydroBoundary_[outer_x1]!=NULL)
    HydroBoundary_[outer_x1](pmb, pco, dst, pmb->is, pmb->ie, bjs, bje, bks, bke);
  if(pmb->block_size.nx2>1) { // 2D or 3D
    if(HydroBoundary_[inner_x2]!=NULL)
      HydroBoundary_[inner_x2](pmb, pco, dst, bis, bie, pmb->js, pmb->je, bks, bke);
    if(HydroBoundary_[outer_x2]!=NULL)
      HydroBoundary_[outer_x2](pmb, pco, dst, bis, bie, pmb->js, pmb->je, bks, bke);
  }
  if(pmb->block_size.nx3>1) { // 3D
    if(pmb->pmy_mesh->face_only==false) {
      bjs=pmb->js-NGHOST;
      bje=pmb->je+NGHOST;
    }
    if(HydroBoundary_[inner_x3]!=NULL)
      HydroBoundary_[inner_x3](pmb, pco, dst, bis, bie, bjs, bje, pmb->ks, pmb->ke);
    if(HydroBoundary_[outer_x3]!=NULL)
      HydroBoundary_[outer_x3](pmb, pco, dst, bis, bie, bjs, bje, pmb->ks, pmb->ke);
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::FieldPhysicalBoundaries(AthenaArray<Real> &dst)
//  \brief Apply physical boundary conditions for field
void BoundaryValues::FieldPhysicalBoundaries(InterfaceField &dst)
{
  MeshBlock *pmb=pmy_mblock_;
  Coordinates *pco=pmb->pcoord;
  int bis=pmb->is-NGHOST, bie=pmb->ie+NGHOST;
  int bjs=pmb->js, bje=pmb->je, bks=pmb->ks, bke=pmb->ke;

  if(FieldBoundary_[inner_x2]==NULL && pmb->block_size.nx2>1) bjs=pmb->js-NGHOST;
  if(FieldBoundary_[outer_x2]==NULL && pmb->block_size.nx2>1) bje=pmb->je+NGHOST;
  if(FieldBoundary_[inner_x3]==NULL && pmb->block_size.nx3>1) bks=pmb->ks-NGHOST;
  if(FieldBoundary_[outer_x3]==NULL && pmb->block_size.nx3>1) bke=pmb->ke+NGHOST;

  if(FieldBoundary_[inner_x1]!=NULL)
    FieldBoundary_[inner_x1](pmb, pco, dst, pmb->is, pmb->ie, bjs, bje, bks, bke);
  if(FieldBoundary_[outer_x1]!=NULL)
    FieldBoundary_[outer_x1](pmb, pco, dst, pmb->is, pmb->ie, bjs, bje, bks, bke);
  if(pmb->block_size.nx2>1) { // 2D or 3D
    if(FieldBoundary_[inner_x2]!=NULL)
      FieldBoundary_[inner_x2](pmb, pco, dst, bis, bie, pmb->js, pmb->je, bks, bke);
    if(FieldBoundary_[outer_x2]!=NULL)
      FieldBoundary_[outer_x2](pmb, pco, dst, bis, bie, pmb->js, pmb->je, bks, bke);
  }
  if(pmb->block_size.nx3>1) { // 3D
    bjs=pmb->js-NGHOST;
    bje=pmb->je+NGHOST;
    if(FieldBoundary_[inner_x3]!=NULL)
      FieldBoundary_[inner_x3](pmb, pco, dst, bis, bie, bjs, bje, pmb->ks, pmb->ke);
    if(FieldBoundary_[outer_x3]!=NULL)
      FieldBoundary_[outer_x3](pmb, pco, dst, bis, bie, bjs, bje, pmb->ks, pmb->ke);
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn unsigned int CreateBufferID(int ox1, int ox2, int ox3, int fi1, int fi2)
//  \brief calculate a buffer identifier
unsigned int CreateBufferID(int ox1, int ox2, int ox3, int fi1, int fi2)
{
  unsigned int ux1=(unsigned)(ox1+1);
  unsigned int ux2=(unsigned)(ox2+1);
  unsigned int ux3=(unsigned)(ox3+1);
  return (ux1<<6) | (ux2<<4) | (ux3<<2) | (fi1<<1) | fi2;
}


//--------------------------------------------------------------------------------------
//! \fn unsigned int CreateMPITag(int lid, int flag, int phys, int bufid)
//  \brief calculate an MPI tag
unsigned int CreateMPITag(int lid, int flag, int phys, int bufid)
{
// tag = local id of destination (18) + flag (2) + physics (4) + bufid(7)
  return (lid<<13) | (flag<<11) | (phys<<7) | bufid;
}


//--------------------------------------------------------------------------------------
//! \fn int BufferID(int dim, bool multilevel, bool face_only)
//  \brief calculate neighbor indexes and target buffer IDs
int BufferID(int dim, bool multilevel, bool face_only)
{
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
        ni_[b].fi1=f1; ni_[b].fi2=f2; ni_[b].type=neighbor_face;
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
          ni_[b].fi1=f1; ni_[b].fi2=f2; ni_[b].type=neighbor_face;
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
          ni_[b].fi1=f1; ni_[b].fi2=f2; ni_[b].type=neighbor_face;
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
        for(int f1=0;f1<nf2;f1++) {
          ni_[b].ox1=n; ni_[b].ox2=m; ni_[b].ox3=0;
          ni_[b].fi1=f1; ni_[b].fi2=0; ni_[b].type=neighbor_edge;
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
          ni_[b].fi1=f1; ni_[b].fi2=0; ni_[b].type=neighbor_edge;
          b++;
        }
      }
    }
    // x2x3
    for(int m=-1; m<=1; m+=2) {
      for(int n=-1; n<=1; n+=2) {
        for(int f1=0;f1<nf1;f1++) {
          ni_[b].ox1=0; ni_[b].ox2=n; ni_[b].ox3=m;
          ni_[b].fi1=f1; ni_[b].fi2=0; ni_[b].type=neighbor_edge;
          b++;
        }
      }
    }
    // corners
    for(int l=-1; l<=1; l+=2) {
      for(int m=-1; m<=1; m+=2) {
        for(int n=-1; n<=1; n+=2) {
          ni_[b].ox1=n; ni_[b].ox2=m; ni_[b].ox3=l;
          ni_[b].fi1=0; ni_[b].fi2=0; ni_[b].type=neighbor_corner;
          b++;
        }
      }
    }
  }

  for(int n=0;n<b;n++)
    bufid_[n]=CreateBufferID(ni_[n].ox1, ni_[n].ox2, ni_[n].ox3, ni_[n].fi1, ni_[n].fi2);

  return b;
}

int FindBufferID(int ox1, int ox2, int ox3, int fi1, int fi2, int bmax)
{
  int bid=CreateBufferID(ox1, ox2, ox3, fi1, fi2);

  for(int i=0;i<bmax;i++) {
    if(bid==bufid_[i]) return i;
  }
  return -1;
}

