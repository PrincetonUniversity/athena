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
#include "../mesh_refinement/mesh_refinement.hpp"
#include "../mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/eos/eos.hpp"
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
  for(int i=0; i<6; i++)
    BoundaryFunction_[i]=NULL;

// Set BC functions for each of the 6 boundaries in turn -------------------------------
  // Inner x1
  nface_=2; nedge_=0;
  switch(pmb->block_bcs[INNER_X1]){
    case REFLECTING_BNDRY:
      BoundaryFunction_[INNER_X1] = ReflectInnerX1;
      break;
    case OUTFLOW_BNDRY:
      BoundaryFunction_[INNER_X1] = OutflowInnerX1;
      break;
    case BLOCK_BNDRY: // block boundary
    case PERIODIC_BNDRY: // periodic boundary
      BoundaryFunction_[INNER_X1] = NULL;
      break;
    case USER_BNDRY: // user-enrolled BCs
      BoundaryFunction_[INNER_X1] = pmb->pmy_mesh->BoundaryFunction_[INNER_X1];
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
          << "Flag ix1_bc=" << pmb->block_bcs[INNER_X1] << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
   }

  // Outer x1
  switch(pmb->block_bcs[OUTER_X1]){
    case REFLECTING_BNDRY:
      BoundaryFunction_[OUTER_X1] = ReflectOuterX1;
      break;
    case OUTFLOW_BNDRY:
      BoundaryFunction_[OUTER_X1] = OutflowOuterX1;
      break;
    case BLOCK_BNDRY: // block boundary
    case PERIODIC_BNDRY: // periodic boundary
      BoundaryFunction_[OUTER_X1] = NULL;
      break;
    case USER_BNDRY: // user-enrolled BCs
      BoundaryFunction_[OUTER_X1] = pmb->pmy_mesh->BoundaryFunction_[OUTER_X1];
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
          << "Flag ox1_bc=" << pmb->block_bcs[OUTER_X1] << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }

  if (pmb->block_size.nx2 > 1) {
    nface_=4; nedge_=4;
    // Inner x2
    switch(pmb->block_bcs[INNER_X2]){
      case REFLECTING_BNDRY:
        BoundaryFunction_[INNER_X2] = ReflectInnerX2;
        break;
      case OUTFLOW_BNDRY:
        BoundaryFunction_[INNER_X2] = OutflowInnerX2;
        break;
      case BLOCK_BNDRY: // block boundary
      case PERIODIC_BNDRY: // periodic boundary
      case POLAR_BNDRY: // polar boundary
        BoundaryFunction_[INNER_X2] = NULL;
        break;
      case USER_BNDRY: // user-enrolled BCs
        BoundaryFunction_[INNER_X2] = pmb->pmy_mesh->BoundaryFunction_[INNER_X2];
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
            << "Flag ix2_bc=" << pmb->block_bcs[INNER_X2] << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
     }

    // Outer x2
    switch(pmb->block_bcs[OUTER_X2]){
      case REFLECTING_BNDRY:
        BoundaryFunction_[OUTER_X2] = ReflectOuterX2;
        break;
      case OUTFLOW_BNDRY:
        BoundaryFunction_[OUTER_X2] = OutflowOuterX2;
        break;
      case BLOCK_BNDRY: // block boundary
      case PERIODIC_BNDRY: // periodic boundary
      case POLAR_BNDRY: // polar boundary
        BoundaryFunction_[OUTER_X2] = NULL;
        break;
      case USER_BNDRY: // user-enrolled BCs
        BoundaryFunction_[OUTER_X2] = pmb->pmy_mesh->BoundaryFunction_[OUTER_X2];
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
            << "Flag ox2_bc=" << pmb->block_bcs[OUTER_X2] << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
    }
  }

  if (pmb->block_size.nx3 > 1) {
    nface_=6; nedge_=12;
    // Inner x3
    switch(pmb->block_bcs[INNER_X3]){
      case REFLECTING_BNDRY:
        BoundaryFunction_[INNER_X3] = ReflectInnerX3;
        break;
      case OUTFLOW_BNDRY:
        BoundaryFunction_[INNER_X3] = OutflowInnerX3;
        break;
      case BLOCK_BNDRY: // block boundary
      case PERIODIC_BNDRY: // periodic boundary
        BoundaryFunction_[INNER_X3] = NULL;
        break;
      case USER_BNDRY: // user-enrolled BCs
        BoundaryFunction_[INNER_X3] = pmb->pmy_mesh->BoundaryFunction_[INNER_X3];
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
            << "Flag ix3_bc=" << pmb->block_bcs[INNER_X3] << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
     }

    // Outer x3
    switch(pmb->block_bcs[OUTER_X3]){
      case REFLECTING_BNDRY:
        BoundaryFunction_[OUTER_X3] = ReflectOuterX3;
        break;
      case OUTFLOW_BNDRY:
        BoundaryFunction_[OUTER_X3] = OutflowOuterX3;
        break;
      case BLOCK_BNDRY: // block boundary
      case PERIODIC_BNDRY: // periodic boundary
        BoundaryFunction_[OUTER_X3] = NULL;
        break;
      case USER_BNDRY: // user-enrolled BCs
        BoundaryFunction_[OUTER_X3] = pmb->pmy_mesh->BoundaryFunction_[OUTER_X3];
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
            << "Flag ox3_bc=" << pmb->block_bcs[OUTER_X3] << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
    }
  }

  // Count number of blocks wrapping around pole
  if (pmb->block_bcs[INNER_X2] == POLAR_BNDRY) {
    int level = pmb->loc.level - pmb->pmy_mesh->root_level;
    num_north_polar_blocks_ = pmb->pmy_mesh->nrbx3 * (1 << level);
  }
  else
    num_north_polar_blocks_ = 0;
  if (pmb->block_bcs[OUTER_X2] == POLAR_BNDRY) {
    int level = pmb->loc.level - pmb->pmy_mesh->root_level;
    num_south_polar_blocks_ = pmb->pmy_mesh->nrbx3 * (1 << level);
  }
  else
    num_south_polar_blocks_ = 0;

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
    if (num_north_polar_blocks_ > 0) {
      emf_north_send_[l] = new Real *[num_north_polar_blocks_];
      emf_north_recv_[l] = new Real *[num_north_polar_blocks_];
      emf_north_flag_[l] = new enum boundary_status[num_north_polar_blocks_];
#ifdef MPI_PARALLEL
      req_emf_north_send_[l] = new MPI_Request[num_north_polar_blocks_];
      req_emf_north_recv_[l] = new MPI_Request[num_north_polar_blocks_];
#endif
      for (int n = 0; n < num_north_polar_blocks_; ++n) {
        emf_north_send_[l][n] = NULL;
        emf_north_recv_[l][n] = NULL;
        emf_north_flag_[l][n] = boundary_waiting;
#ifdef MPI_PARALLEL
        req_emf_north_send_[l][n] = MPI_REQUEST_NULL;
        req_emf_north_recv_[l][n] = MPI_REQUEST_NULL;
#endif
      }
    }
    if (num_south_polar_blocks_ > 0) {
      emf_south_send_[l] = new Real *[num_south_polar_blocks_];
      emf_south_recv_[l] = new Real *[num_south_polar_blocks_];
      emf_south_flag_[l] = new enum boundary_status[num_south_polar_blocks_];
#ifdef MPI_PARALLEL
      req_emf_south_send_[l] = new MPI_Request[num_south_polar_blocks_];
      req_emf_south_recv_[l] = new MPI_Request[num_south_polar_blocks_];
#endif
      for (int n = 0; n < num_south_polar_blocks_; ++n) {
        emf_south_send_[l][n] = NULL;
        emf_south_recv_[l][n] = NULL;
        emf_south_flag_[l][n] = boundary_waiting;
#ifdef MPI_PARALLEL
        req_emf_south_send_[l][n] = MPI_REQUEST_NULL;
        req_emf_south_recv_[l][n] = MPI_REQUEST_NULL;
#endif
      }
    }
  }

  // Allocate buffers for non-polar neighbor communication
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

  // Allocate buffers for polar neighbor communication
  if (MAGNETIC_FIELDS_ENABLED) {
    if (num_north_polar_blocks_ > 0) {
      for (int l = 0; l < NSTEP; ++l) {
        for (int n = 0; n < num_north_polar_blocks_; ++n) {
          emf_north_send_[l][n] = new Real[pmb->block_size.nx1];
          emf_north_recv_[l][n] = new Real[pmb->block_size.nx1];
        }
      }
    }
    if (num_south_polar_blocks_ > 0) {
      for (int l = 0; l < NSTEP; ++l) {
        for (int n = 0; n < num_south_polar_blocks_; ++n) {
          emf_south_send_[l][n] = new Real[pmb->block_size.nx1];
          emf_south_recv_[l][n] = new Real[pmb->block_size.nx1];
        }
      }
    }
  }

  if(pmb->pmy_mesh->multilevel==true) { // SMR or AMR
    // allocate surface area array
    int nc1=pmb->block_size.nx1+2*NGHOST;
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
  }

 /* single CPU in the azimuthal direction with the polar boundary*/
  if(pmb->loc.level == pmb->pmy_mesh->root_level &&
     pmb->pmy_mesh->nrbx3 == 1 &&
     (pmb->block_bcs[INNER_X2]==POLAR_BNDRY||pmb->block_bcs[OUTER_X2]==POLAR_BNDRY))
       exc_.NewAthenaArray(pmb->ke+NGHOST+2);

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
  if (MAGNETIC_FIELDS_ENABLED) {
    if (num_north_polar_blocks_ > 0) {
      for (int l = 0; l < NSTEP; ++l) {
        for (int n = 0; n < num_north_polar_blocks_; ++n) {
          delete[] emf_north_send_[l][n];
          delete[] emf_north_recv_[l][n];
        }
      }
    }
    if (num_south_polar_blocks_ > 0) {
      for (int l = 0; l < NSTEP; ++l) {
        for (int n = 0; n < num_south_polar_blocks_; ++n) {
          delete[] emf_south_send_[l][n];
          delete[] emf_south_recv_[l][n];
        }
      }
    }
  }
  if(pmb->pmy_mesh->multilevel==true) {
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
  }
  if (num_north_polar_blocks_ > 0) {
    for (int l = 0; l < NSTEP; ++l) {
      delete[] emf_north_send_[l];
      delete[] emf_north_recv_[l];
      delete[] emf_north_flag_[l];
#ifdef MPI_PARALLEL
      delete[] req_emf_north_send_[l];
      delete[] req_emf_north_recv_[l];
#endif
    }
  }
  if (num_south_polar_blocks_ > 0) {
    for (int l = 0; l < NSTEP; ++l) {
      delete[] emf_south_send_[l];
      delete[] emf_south_recv_[l];
      delete[] emf_south_flag_[l];
#ifdef MPI_PARALLEL
      delete[] req_emf_south_send_[l];
      delete[] req_emf_south_recv_[l];
#endif
    }
  }
  if(pmb->loc.level == pmb->pmy_mesh->root_level &&
     pmb->pmy_mesh->nrbx3 == 1 &&
     (pmb->block_bcs[INNER_X2]==POLAR_BNDRY||pmb->block_bcs[OUTER_X2]==POLAR_BNDRY))
       exc_.DeleteAthenaArray();
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::Initialize(void)
//  \brief Initialize MPI requests
void BoundaryValues::Initialize(void)
{
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

#ifdef MPI_PARALLEL
  // Initialize non-polar neighbor communications to other ranks
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
              if(nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
                size=(pmb->block_size.nx2+1)*(pmb->block_size.nx3)
                    +(pmb->block_size.nx2)*(pmb->block_size.nx3+1);
                f2csize=(pmb->block_size.nx2/2+1)*(pmb->block_size.nx3/2)
                    +(pmb->block_size.nx2/2)*(pmb->block_size.nx3/2+1);
              }
              else if(nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
                size=(pmb->block_size.nx1+1)*(pmb->block_size.nx3)
                    +(pmb->block_size.nx1)*(pmb->block_size.nx3+1);
                f2csize=(pmb->block_size.nx1/2+1)*(pmb->block_size.nx3/2)
                    +(pmb->block_size.nx1/2)*(pmb->block_size.nx3/2+1);
              }
              else if(nb.fid==INNER_X3 || nb.fid==OUTER_X3) {
                size=(pmb->block_size.nx1+1)*(pmb->block_size.nx2)
                    +(pmb->block_size.nx1)*(pmb->block_size.nx2+1);
                f2csize=(pmb->block_size.nx1/2+1)*(pmb->block_size.nx2/2)
                    +(pmb->block_size.nx1/2)*(pmb->block_size.nx2/2+1);
              }
            }
            else if(pmb->block_size.nx2 > 1) { // 2D
              if(nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
                size=(pmb->block_size.nx2+1)+pmb->block_size.nx2;
                f2csize=(pmb->block_size.nx2/2+1)+pmb->block_size.nx2/2;
              }
              else if(nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
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

  // Initialize polar neighbor communications to other ranks
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int l = 0; l < NSTEP; ++l) {
      for (int n = 0; n < num_north_polar_blocks_; ++n) {
        const PolarNeighborBlock &nb = pmb->polar_neighbor_north[n];
        if(nb.rank != Globals::my_rank) {
          tag = CreateMPITag(nb.lid, l, tag_emfpole, pmb->loc.lx3);
          MPI_Send_init(emf_north_send_[l][n], pmb->block_size.nx1, MPI_ATHENA_REAL,
              nb.rank, tag, MPI_COMM_WORLD, &req_emf_north_send_[l][n]);
          tag = CreateMPITag(pmb->lid, l, tag_emfpole, n);
          MPI_Recv_init(emf_north_recv_[l][n], pmb->block_size.nx1, MPI_ATHENA_REAL,
              nb.rank, tag, MPI_COMM_WORLD, &req_emf_north_recv_[l][n]);
        }
      }
      for (int n = 0; n < num_south_polar_blocks_; ++n) {
        const PolarNeighborBlock &nb = pmb->polar_neighbor_south[n];
        if(nb.rank != Globals::my_rank) {
          tag = CreateMPITag(nb.lid, l, tag_emfpole, pmb->loc.lx3);
          MPI_Send_init(emf_south_send_[l][n], pmb->block_size.nx1, MPI_ATHENA_REAL,
              nb.rank, tag, MPI_COMM_WORLD, &req_emf_south_send_[l][n]);
          tag = CreateMPITag(pmb->lid, l, tag_emfpole, n);
          MPI_Recv_init(emf_south_recv_[l][n], pmb->block_size.nx1, MPI_ATHENA_REAL,
              nb.rank, tag, MPI_COMM_WORLD, &req_emf_south_recv_[l][n]);
        }
      }
    }
  }
#endif
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::CheckBoundary(void)
//  \brief checks if the boundary conditions are correctly enrolled
void BoundaryValues::CheckBoundary(void)
{
  MeshBlock *pmb=pmy_mblock_;
  for(int i=0;i<nface_;i++) {
    if(pmb->block_bcs[i]==USER_BNDRY) {
      if(BoundaryFunction_[i]==NULL) {
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues::CheckBoundary" << std::endl
            << "A user-defined boundary is specified but the hydro boundary function "
            << "is not enrolled in direction " << i  << "." << std::endl;
        throw std::runtime_error(msg.str().c_str());
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
      // Prep sending primitives to enable cons->prim inversion before prolongation
      if (GENERAL_RELATIVITY and pmb->pmy_mesh->multilevel)
        MPI_Start(&req_hydro_recv_[1][nb.bufid]);
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
  for(int l=0;l<NSTEP;l++)
    firsttime_[l]=true;
#ifdef MPI_PARALLEL
  MeshBlock *pmb=pmy_mblock_;
  int mylevel=pmb->loc.level;
  for(int l=0;l<NSTEP;l++) {
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
    if (MAGNETIC_FIELDS_ENABLED) {
      for (int n = 0; n < num_north_polar_blocks_; ++n) {
        const PolarNeighborBlock &nb = pmb->polar_neighbor_north[n];
        if (nb.rank != Globals::my_rank) {
          MPI_Start(&req_emf_north_recv_[l][n]);
        }
      }
      for (int n = 0; n < num_south_polar_blocks_; ++n) {
        const PolarNeighborBlock &nb = pmb->polar_neighbor_south[n];
        if (nb.rank != Globals::my_rank) {
          MPI_Start(&req_emf_south_recv_[l][n]);
        }
      }
    }
  }
#endif
  return;
}


//--------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadHydroBoundaryBufferSameLevel(AthenaArray<Real> &src,
//                                                 Real *buf, const NeighborBlock& nb)
//  \brief Set hydro boundary buffers for sending to a block on the same level
int BoundaryValues::LoadHydroBoundaryBufferSameLevel(AthenaArray<Real> &src, Real *buf,
                                                     const NeighborBlock& nb)
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
  BufferUtility::Pack4DData(src, buf, 0, NHYDRO-1, si, ei, sj, ej, sk, ek, p);
  return p;
}


//--------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadHydroBoundaryBufferToCoarser(AthenaArray<Real> &src,
//                                                 Real *buf, const NeighborBlock& nb,
//                                                 bool conserved_values)
//  \brief Set hydro boundary buffers for sending to a block on the coarser level
int BoundaryValues::LoadHydroBoundaryBufferToCoarser(AthenaArray<Real> &src, Real *buf,
                                                     const NeighborBlock& nb,
                                                     bool conserved_values)
{
  MeshBlock *pmb=pmy_mblock_;
  MeshRefinement *pmr=pmb->pmr;
  int si, sj, sk, ei, ej, ek;
  int cn=pmb->cnghost-1;

  si=(nb.ox1>0)?(pmb->cie-cn):pmb->cis;
  ei=(nb.ox1<0)?(pmb->cis+cn):pmb->cie;
  sj=(nb.ox2>0)?(pmb->cje-cn):pmb->cjs;
  ej=(nb.ox2<0)?(pmb->cjs+cn):pmb->cje;
  sk=(nb.ox3>0)?(pmb->cke-cn):pmb->cks;
  ek=(nb.ox3<0)?(pmb->cks+cn):pmb->cke;

  int p=0;
  if (conserved_values) { // normal case; restrict the data before sending
    pmr->RestrictCellCenteredValues(src, pmr->coarse_cons_, 0, NHYDRO-1,
                                    si, ei, sj, ej, sk, ek);
    BufferUtility::Pack4DData(pmr->coarse_cons_, buf, 0, NHYDRO-1,
                              si, ei, sj, ej, sk, ek, p);
  }
  else { // must be initialization; need to restrict but not send primitives
    pmr->RestrictCellCenteredValues(src, pmr->coarse_prim_, 0, NHYDRO-1,
                                    si, ei, sj, ej, sk, ek);
  }
  return p;
}


//--------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadHydroBoundaryBufferToFiner(AthenaArray<Real> &src,
//                                                 Real *buf, const NeighborBlock& nb)
//  \brief Set hydro boundary buffers for sending to a block on the finer level
int BoundaryValues::LoadHydroBoundaryBufferToFiner(AthenaArray<Real> &src, Real *buf,
                                                   const NeighborBlock& nb)
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
  BufferUtility::Pack4DData(src, buf, 0, NHYDRO-1, si, ei, sj, ej, sk, ek, p);
  return p;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SendHydroBoundaryBuffers(AthenaArray<Real> &src, int step,
//                                                    bool conserved_values)
//  \brief Send boundary buffers
void BoundaryValues::SendHydroBoundaryBuffers(AthenaArray<Real> &src, int step,
                                              bool conserved_values)
{
  MeshBlock *pmb=pmy_mblock_;
  int mylevel=pmb->loc.level;

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    int ssize;
    if(nb.level==mylevel)
      ssize=LoadHydroBoundaryBufferSameLevel(src, hydro_send_[step][nb.bufid],nb);
    else if(nb.level<mylevel)
      ssize=LoadHydroBoundaryBufferToCoarser(src, hydro_send_[step][nb.bufid], nb,
                                             conserved_values);
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
//                                           Real *buf, const NeighborBlock& nb)
//  \brief Set hydro boundary received from a block on the same level
void BoundaryValues::SetHydroBoundarySameLevel(AthenaArray<Real> &dst, Real *buf,
                                               const NeighborBlock& nb)
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
  if (nb.polar) {
    for (int n=0; n<(NHYDRO); ++n) {
      Real sign = flip_across_pole_hydro[n] ? -1.0 : 1.0;
      for (int k=sk; k<=ek; ++k) {
        for (int j=ej; j>=sj; --j) {
#pragma simd
          for (int i=si; i<=ei; ++i)
            dst(n,k,j,i) = sign * buf[p++];
        }
      }
    }
  }
  else
    BufferUtility::Unpack4DData(buf, dst, 0, NHYDRO-1, si, ei, sj, ej, sk, ek, p);
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetHydroBoundaryFromCoarser(Real *buf,
//                                                       const NeighborBlock& nb,
//                                                       bool conserved_values)
//  \brief Set hydro prolongation buffer received from a block on a coarser level
void BoundaryValues::SetHydroBoundaryFromCoarser(Real *buf, const NeighborBlock& nb,
    bool conserved_values)
{
  MeshBlock *pmb=pmy_mblock_;
  MeshRefinement *pmr=pmb->pmr;

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
  if (nb.polar) {
    for (int n=0; n<(NHYDRO); ++n) {
      Real sign = flip_across_pole_hydro[n] ? -1.0 : 1.0;
      for (int k=sk; k<=ek; ++k) {
        for (int j=ej; j>=sj; --j) {
#pragma simd
          for (int i=si; i<=ei; ++i) {
            if (conserved_values)
              pmr->coarse_cons_(n,k,j,i) = sign * buf[p++];
            else
              pmr->coarse_prim_(n,k,j,i) = sign * buf[p++];
          }
        }
      }
    }
  }
  else {
    if (conserved_values)
      BufferUtility::Unpack4DData(buf, pmr->coarse_cons_, 0, NHYDRO-1,
                                  si, ei, sj, ej, sk, ek, p);
    else
      BufferUtility::Unpack4DData(buf, pmr->coarse_prim_, 0, NHYDRO-1,
                                  si, ei, sj, ej, sk, ek, p);
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetHydroBoundaryFromFiner(AthenaArray<Real> &dst,
//                                               Real *buf, const NeighborBlock& nb)
//  \brief Set hydro boundary received from a block on a finer level
void BoundaryValues::SetHydroBoundaryFromFiner(AthenaArray<Real> &dst, Real *buf,
                                               const NeighborBlock& nb)
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
  if (nb.polar) {
    for (int n=0; n<(NHYDRO); ++n) {
      Real sign = flip_across_pole_hydro[n] ? -1.0 : 1.0;
      for (int k=sk; k<=ek; ++k) {
        for (int j=sj; j<=ej; ++j) {
#pragma simd
          for (int i=si; i<=ei; ++i)
            dst(n,k,j,i) = sign * buf[p++];
        }
      }
    }
  }
  else 
    BufferUtility::Unpack4DData(buf, dst, 0, NHYDRO-1, si, ei, sj, ej, sk, ek, p);
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
      SetHydroBoundaryFromCoarser(hydro_recv_[step][nb.bufid], nb, true);
    else
      SetHydroBoundaryFromFiner(dst, hydro_recv_[step][nb.bufid], nb);
    hydro_flag_[step][nb.bufid] = boundary_completed; // completed
  }

  if(flag&& (pmb->block_bcs[INNER_X2]==POLAR_BNDRY
         ||  pmb->block_bcs[OUTER_X2]==POLAR_BNDRY))
     PolarSingleHydro(dst);
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
      MPI_Wait(&req_hydro_recv_[step][nb.bufid],MPI_STATUS_IGNORE);
#endif
    if(nb.level==pmb->loc.level)
      SetHydroBoundarySameLevel(dst, hydro_recv_[step][nb.bufid], nb);
    else if(nb.level<pmb->loc.level)
      SetHydroBoundaryFromCoarser(hydro_recv_[step][nb.bufid], nb, step == 0);
    else
      SetHydroBoundaryFromFiner(dst, hydro_recv_[step][nb.bufid], nb);
    hydro_flag_[step][nb.bufid] = boundary_completed; // completed
  }
 
  if (pmb->block_bcs[INNER_X2]==POLAR_BNDRY||pmb->block_bcs[OUTER_X2]==POLAR_BNDRY)
    PolarSingleHydro(dst);

  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::PolarSingleHydro(AthenaArray<Real> &dst)
//
// \brief  single CPU in the azimuthal direction for the polar boundary 
void BoundaryValues::PolarSingleHydro(AthenaArray<Real> &dst)
{
  MeshBlock *pmb=pmy_mblock_;
  if(pmb->loc.level == pmb->pmy_mesh->root_level && pmb->pmy_mesh->nrbx3 == 1){

    if(pmb->block_bcs[INNER_X2]==POLAR_BNDRY){
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      for (int n=0; n<(NHYDRO); ++n) {
        for (int j=pmb->js-NGHOST; j<=pmb->js-1; ++j) {
         for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i){
           for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
             exc_(k)=dst(n,k,j,i);
           }
           for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
             int k_shift = k;
             k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
             dst(n,k,j,i)=exc_(k_shift);
           }
         }
        }
      }
    }

    if(pmb->block_bcs[OUTER_X2]==POLAR_BNDRY){
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      for (int n=0; n<(NHYDRO); ++n) {
        for (int j=pmb->je+1; j<=pmb->je+NGHOST; ++j) {
         for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i){
           for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
             exc_(k)=dst(n,k,j,i);
           }
           for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
             int k_shift = k;
             k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
             dst(n,k,j,i)=exc_(k_shift);
           }
         }
        }
      }
    }
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

      flcor_flag_[step][nb.fid][nb.fi2][nb.fi1] = boundary_completed;
    }
  }

  return flag;
}

//--------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadFieldBoundaryBufferSameLevel(FaceField &src,
//                                                 Real *buf, const NeighborBlock& nb)
//  \brief Set field boundary buffers for sending to a block on the same level
int BoundaryValues::LoadFieldBoundaryBufferSameLevel(FaceField &src, Real *buf,
                                                     const NeighborBlock& nb)
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
  BufferUtility::Pack3DData(src.x1f, buf, si, ei, sj, ej, sk, ek, p);

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
  BufferUtility::Pack3DData(src.x2f, buf, si, ei, sj, ej, sk, ek, p);

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
  BufferUtility::Pack3DData(src.x3f, buf, si, ei, sj, ej, sk, ek, p);

  return p;
}

//--------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadFieldBoundaryBufferToCoarser(FaceField &src,
//                                                 Real *buf, const NeighborBlock& nb)
//  \brief Set field boundary buffers for sending to a block on the coarser level
int BoundaryValues::LoadFieldBoundaryBufferToCoarser(FaceField &src, Real *buf,
                                                     const NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_mblock_;
  MeshRefinement *pmr=pmb->pmr;
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
  pmr->RestrictFieldX1(src.x1f, pmr->coarse_b_.x1f, si, ei, sj, ej, sk, ek);
  BufferUtility::Pack3DData(pmr->coarse_b_.x1f, buf, si, ei, sj, ej, sk, ek, p);

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
  pmr->RestrictFieldX2(src.x2f, pmr->coarse_b_.x2f, si, ei, sj, ej, sk, ek);
  BufferUtility::Pack3DData(pmr->coarse_b_.x2f, buf, si, ei, sj, ej, sk, ek, p);

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
  pmr->RestrictFieldX3(src.x3f, pmr->coarse_b_.x3f, si, ei, sj, ej, sk, ek);
  BufferUtility::Pack3DData(pmr->coarse_b_.x3f, buf, si, ei, sj, ej, sk, ek, p);

  return p;
}

//--------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadFieldBoundaryBufferToFiner(FaceField &src, 
//                                                 Real *buf, const NeighborBlock& nb)
//  \brief Set field boundary buffers for sending to a block on the finer level
int BoundaryValues::LoadFieldBoundaryBufferToFiner(FaceField &src, Real *buf,
                                                   const NeighborBlock& nb)
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
  BufferUtility::Pack3DData(src.x1f, buf, si, ei, sj, ej, sk, ek, p);

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
  BufferUtility::Pack3DData(src.x2f, buf, si, ei, sj, ej, sk, ek, p);

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
  BufferUtility::Pack3DData(src.x3f, buf, si, ei, sj, ej, sk, ek, p);

  return p;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SendFieldBoundaryBuffers(FaceField &src, int step)
//  \brief Send field boundary buffers
void BoundaryValues::SendFieldBoundaryBuffers(FaceField &src, int step)
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
//! \fn void BoundaryValues::SetFieldBoundarySameLevel(FaceField &dst,
//                                               Real *buf, const NeighborBlock& nb)
//  \brief Set field boundary received from a block on the same level
void BoundaryValues::SetFieldBoundarySameLevel(FaceField &dst, Real *buf,
                                               const NeighborBlock& nb)
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
  if (nb.polar) {
    Real sign = flip_across_pole_field[IB1] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma simd
        for (int i=si; i<=ei; ++i)
          dst.x1f(k,j,i)=sign*buf[p++];
      }
    }
  }
  else 
    BufferUtility::Unpack3DData(buf, dst.x1f, si, ei, sj, ej, sk, ek, p);

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
  if (nb.polar) {
    Real sign = flip_across_pole_field[IB2] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma simd
        for (int i=si; i<=ei; ++i)
          dst.x2f(k,j,i)=sign*buf[p++];
      }
    }
    if (nb.ox2 < 0) {  // interpolate B^theta across top pole
      for (int k=sk; k<=ek; ++k) {
        for (int i=si; i<=ei; ++i)
          dst.x2f(k,pmb->js,i) = 0.5 * (dst.x2f(k,pmb->js-1,i) + dst.x2f(k,pmb->js+1,i));
      }
    }
    if (nb.ox2 > 0) {  // interpolate B^theta across bottom pole
      for (int k=sk; k<=ek; ++k) {
        for (int i=si; i<=ei; ++i)
          dst.x2f(k,pmb->je+1,i) = 0.5 * (dst.x2f(k,pmb->je,i) + dst.x2f(k,pmb->je+2,i));
      }
    }
  }
  else
    BufferUtility::Unpack3DData(buf, dst.x2f, si, ei, sj, ej, sk, ek, p);
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
  if (nb.polar) {
    Real sign = flip_across_pole_field[IB3] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma simd
        for (int i=si; i<=ei; ++i)
          dst.x3f(k,j,i)=sign*buf[p++];
      }
    }
  }
  else 
    BufferUtility::Unpack3DData(buf, dst.x3f, si, ei, sj, ej, sk, ek, p);
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
//! \fn void BoundaryValues::SetFieldBoundaryFromCoarser(Real *buf,
//                                                       const NeighborBlock& nb)
//  \brief Set field prolongation buffer received from a block on the same level
void BoundaryValues::SetFieldBoundaryFromCoarser(Real *buf, const NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_mblock_;
  MeshRefinement *pmr=pmb->pmr;
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

  if (nb.polar) {
    Real sign = flip_across_pole_field[IB1] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma simd
        for (int i=si; i<=ei; ++i)
          pmr->coarse_b_.x1f(k,j,i) = sign*buf[p++];
      }
    }
  }
  else 
    BufferUtility::Unpack3DData(buf, pmr->coarse_b_.x1f, si, ei, sj, ej, sk, ek, p);

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

  if (nb.polar) {
    Real sign = flip_across_pole_field[IB2] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma simd
        for (int i=si; i<=ei; ++i)
          pmr->coarse_b_.x2f(k,j,i) = sign*buf[p++];
      }
    }
  }
  else
    BufferUtility::Unpack3DData(buf, pmr->coarse_b_.x2f, si, ei, sj, ej, sk, ek, p);

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

  if (nb.polar) {
    Real sign = flip_across_pole_field[IB3] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma simd
        for (int i=si; i<=ei; ++i)
          pmr->coarse_b_.x3f(k,j,i) = sign*buf[p++];
      }
    }
  }
  else
    BufferUtility::Unpack3DData(buf, pmr->coarse_b_.x3f, si, ei, sj, ej, sk, ek, p);

  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetFielBoundaryFromFiner(FaceField &dst,
//                                                    Real *buf, const NeighborBlock& nb)
//  \brief Set field boundary received from a block on the same level
void BoundaryValues::SetFieldBoundaryFromFiner(FaceField &dst, Real *buf,
                                               const NeighborBlock& nb)
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

  if (nb.polar) {
    Real sign = flip_across_pole_field[IB1] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma simd
        for (int i=si; i<=ei; ++i)
          dst.x1f(k,j,i) = sign*buf[p++];
      }
    }
  }
  else 
    BufferUtility::Unpack3DData(buf, dst.x1f, si, ei, sj, ej, sk, ek, p);

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

  if (nb.polar) {
    Real sign = flip_across_pole_field[IB2] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma simd
        for (int i=si; i<=ei; ++i)
          dst.x2f(k,j,i) = sign*buf[p++];
      }
    }
  }
  else
    BufferUtility::Unpack3DData(buf, dst.x2f, si, ei, sj, ej, sk, ek, p);
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

  if (nb.polar) {
    Real sign = flip_across_pole_field[IB3] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma simd
        for (int i=si; i<=ei; ++i)
          dst.x3f(k,j,i) = sign*buf[p++];
      }
    }
  }
  else 
    BufferUtility::Unpack3DData(buf, dst.x3f, si, ei, sj, ej, sk, ek, p);
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
//! \fn bool BoundaryValues::ReceiveFieldBoundaryBuffers(FaceField &dst, int step)
//  \brief load boundary buffer for x1 direction into the array
bool BoundaryValues::ReceiveFieldBoundaryBuffers(FaceField &dst, int step)
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

  if(flag&&(pmb->block_bcs[INNER_X2]==POLAR_BNDRY||pmb->block_bcs[OUTER_X2]==POLAR_BNDRY))
    PolarSingleField(dst);

  return flag;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ReceiveFieldBoundaryBuffersWithWait(FaceField &dst,
//                                                               int step)
//  \brief load boundary buffer for x1 direction into the array
void BoundaryValues::ReceiveFieldBoundaryBuffersWithWait(FaceField &dst, int step)
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

  if(pmb->block_bcs[INNER_X2]==POLAR_BNDRY||pmb->block_bcs[OUTER_X2]==POLAR_BNDRY)
    PolarSingleField(dst);
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::PolarSingleField(FaceField &dst)
//
//  \brief single CPU in the azimuthal direction for the polar boundary
void BoundaryValues::PolarSingleField(FaceField &dst)
{
  MeshBlock *pmb=pmy_mblock_;
  if(pmb->loc.level == pmb->pmy_mesh->root_level && pmb->pmy_mesh->nrbx3 == 1){
    if(pmb->block_bcs[INNER_X2]==POLAR_BNDRY){
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      for (int j=pmb->js-NGHOST; j<=pmb->js-1; ++j) {
       for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST+1; ++i){
         for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
           exc_(k)=dst.x1f(k,j,i);
         }
         for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
           int k_shift = k;
           k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
           dst.x1f(k,j,i)=exc_(k_shift);
         }
       }
      }
      for (int j=pmb->js-NGHOST; j<=pmb->js-1; ++j) {
       for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i){
         for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
           exc_(k)=dst.x2f(k,j,i);
         }
         for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
           int k_shift = k;
           k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
           dst.x2f(k,j,i)=exc_(k_shift);
         }
       }
      }
      for (int j=pmb->js-NGHOST; j<=pmb->js-1; ++j) {
       for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i){
         for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST+1; ++k) {
           exc_(k)=dst.x3f(k,j,i);
         }
         for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST+1; ++k) {
           int k_shift = k;
           k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
           dst.x3f(k,j,i)=exc_(k_shift);
         }
       }
      }
        /* average B2 across the pole */
      for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
        for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i){
          dst.x2f(k,pmb->js,i) = 0.5*(dst.x2f(k,pmb->js-1,i)+dst.x2f(k,pmb->js+1,i));
        }
      }
    }

    if(pmb->block_bcs[OUTER_X2]==POLAR_BNDRY){
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      for (int j=pmb->je+1; j<=pmb->je+NGHOST; ++j) {
        for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST+1; ++i){
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
            exc_(k)=dst.x1f(k,j,i);
          }
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
            int k_shift = k;
            k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
            dst.x1f(k,j,i)=exc_(k_shift);
          }
        }
      }
      for (int j=pmb->je+2; j<=pmb->je+NGHOST+1; ++j) {
        for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i){
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
            exc_(k)=dst.x2f(k,j,i);
          }
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
            int k_shift = k;
            k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
            dst.x2f(k,j,i)=exc_(k_shift);
          }
        }
      }
      for (int j=pmb->je+1; j<=pmb->je+NGHOST; ++j) {
        for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i){
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST+1; ++k) {
            exc_(k)=dst.x3f(k,j,i);
          }
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST+1; ++k) {
            int k_shift = k;
            k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
            dst.x3f(k,j,i)=exc_(k_shift);
          }
        }
      }
         /* average B2 across the pole */
      for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
        for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i){
          dst.x2f(k,pmb->je+1,i) = 0.5*(dst.x2f(k,pmb->je,i)+dst.x2f(k,pmb->je+2,i));
        }
      }

    }
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadEMFBoundaryBufferSameLevel(Real *buf,
//                                                         const NeighborBlock& nb)
//  \brief Set EMF correction buffers for sending to a block on the same level
int BoundaryValues::LoadEMFBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_mblock_;
  AthenaArray<Real> &e1=pmb->pfield->e.x1e;
  AthenaArray<Real> &e2=pmb->pfield->e.x2e;
  AthenaArray<Real> &e3=pmb->pfield->e.x3e;
  int p=0;
  if(nb.type==neighbor_face) {
    if(pmb->block_size.nx3 > 1) { // 3D
      // x1 direction
      if(nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
        int i;
        if(nb.fid==INNER_X1) i=pmb->is;
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
      else if(nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
        int j;
        if(nb.fid==INNER_X2) j=pmb->js;
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
      else if(nb.fid==INNER_X3 || nb.fid==OUTER_X3) {
        int k;
        if(nb.fid==INNER_X3) k=pmb->ks;
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
      if(nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
        int i;
        if(nb.fid==INNER_X1) i=pmb->is;
        else i=pmb->ie+1;
        // pack e2
        for(int j=pmb->js; j<=pmb->je; j++)
          buf[p++]=e2(k,j,i);
        // pack e3
        for(int j=pmb->js; j<=pmb->je+1; j++)
          buf[p++]=e3(k,j,i);
      }
      // x2 direction
      else if(nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
        int j;
        if(nb.fid==INNER_X2) j=pmb->js;
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
      if(nb.fid==INNER_X1) i=pmb->is;
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
//! \fn int BoundaryValues::LoadEMFBoundaryBufferToCoarser(Real *buf,
//                                                         const NeighborBlock& nb)
//  \brief Set EMF correction buffers for sending to a block on the coarser level
int BoundaryValues::LoadEMFBoundaryBufferToCoarser(Real *buf, const NeighborBlock& nb)
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
      if(nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
        int i;
        if(nb.fid==INNER_X1) i=pmb->is;
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
      else if(nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
        int j;
        if(nb.fid==INNER_X2) j=pmb->js;
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
      else if(nb.fid==INNER_X3 || nb.fid==OUTER_X3) {
        int k;
        if(nb.fid==INNER_X3) k=pmb->ks;
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
      if(nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
        int i;
        if(nb.fid==INNER_X1) i=pmb->is;
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
      else if(nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
        int j;
        if(nb.fid==INNER_X2) j=pmb->js;
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
      if(nb.fid==INNER_X1) i=pmb->is;
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
//! \fn int BoundaryValues::LoadEMFBoundaryPolarBuffer(Real *buf,
//          const PolarNeighborBlock &nb)
//  \brief Load EMF values along polar axis into send buffers
int BoundaryValues::LoadEMFBoundaryPolarBuffer(Real *buf, const PolarNeighborBlock &nb)
{
  MeshBlock *pmb = pmy_mblock_;
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

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SendEMFCorrection(int step)
//  \brief Restrict, pack and send the surace EMF to the coarse neighbor(s) if needed
void BoundaryValues::SendEMFCorrection(int step)
{
  MeshBlock *pmb=pmy_mblock_;

  // Send non-polar EMF values
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

  // Send polar EMF values
  for (int n = 0; n < num_north_polar_blocks_; ++n) {
    const PolarNeighborBlock &nb = pmb->polar_neighbor_north[n];
    int count = LoadEMFBoundaryPolarBuffer(emf_north_send_[step][n], nb);
    if (nb.rank == Globals::my_rank) { // on the same node
      MeshBlock *pbl = pmb->pmy_mesh->FindMeshBlock(nb.gid);
      std::memcpy(pbl->pbval->emf_north_recv_[step][pmb->loc.lx3],
          emf_north_send_[step][n], count * sizeof(Real));
      pbl->pbval->emf_north_flag_[step][pmb->loc.lx3] = boundary_arrived;
    }
#ifdef MPI_PARALLEL
    else
      MPI_Start(&req_emf_north_send_[step][n]);
#endif
  }
  for (int n = 0; n < num_south_polar_blocks_; ++n) {
    const PolarNeighborBlock &nb = pmb->polar_neighbor_south[n];
    int count = LoadEMFBoundaryPolarBuffer(emf_south_send_[step][n], nb);
    if (nb.rank == Globals::my_rank) { // on the same node
      MeshBlock *pbl = pmb->pmy_mesh->FindMeshBlock(nb.gid);
      std::memcpy(pbl->pbval->emf_south_recv_[step][pmb->loc.lx3],
          emf_south_send_[step][n], count * sizeof(Real));
      pbl->pbval->emf_south_flag_[step][pmb->loc.lx3] = boundary_arrived;
    }
#ifdef MPI_PARALLEL
    else
      MPI_Start(&req_emf_south_send_[step][n]);
#endif
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetEMFBoundarySameLevel(Real *buf, const NeighborBlock& nb)
//  \brief Add up the EMF received from a block on the same level
//         Later they will be divided in the AverageEMFBoundary function
void BoundaryValues::SetEMFBoundarySameLevel(Real *buf, const NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_mblock_;
  AthenaArray<Real> &e1=pmb->pfield->e.x1e;
  AthenaArray<Real> &e2=pmb->pfield->e.x2e;
  AthenaArray<Real> &e3=pmb->pfield->e.x3e;
  int p=0;
  if(nb.type==neighbor_face) {
    if(pmb->block_size.nx3 > 1) { // 3D
      // x1 direction
      if(nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
        int i;
        if(nb.fid==INNER_X1) i=pmb->is;
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
      else if(nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
        int j;
        if(nb.fid==INNER_X2) j=pmb->js;
        else j=pmb->je+1;
        // unpack e1
        Real sign = (nb.polar and flip_across_pole_field[IB1]) ? -1.0 : 1.0;
        for(int k=pmb->ks; k<=pmb->ke+1; k++) {
          for(int i=pmb->is; i<=pmb->ie; i++)
            e1(k,j,i)+=sign*buf[p++];
        }
        // unpack e3
        sign = (nb.polar and flip_across_pole_field[IB3]) ? -1.0 : 1.0;
        for(int k=pmb->ks; k<=pmb->ke; k++) {
          for(int i=pmb->is; i<=pmb->ie+1; i++)
            e3(k,j,i)+=sign*buf[p++];
        }
      }
      // x3 direction
      else if(nb.fid==INNER_X3 || nb.fid==OUTER_X3) {
        int k;
        if(nb.fid==INNER_X3) k=pmb->ks;
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
      if(nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
        int i;
        if(nb.fid==INNER_X1) i=pmb->is;
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
      else if(nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
        int j;
        if(nb.fid==INNER_X2) j=pmb->js;
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
      if(nb.fid==INNER_X1) i=pmb->is;
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
    // x1x2 edge (2D and 3D)
    if(nb.eid>=0 && nb.eid<4) {
      int i, j;
      if((nb.eid&1)==0) i=pmb->is;
      else i=pmb->ie+1;
      if((nb.eid&2)==0) j=pmb->js;
      else j=pmb->je+1;
      // unpack e3
      Real sign = (nb.polar and flip_across_pole_field[IB3]) ? -1.0 : 1.0;
      for(int k=pmb->ks; k<=pmb->ke; k++) {
        e3(k,j,i)+=sign*buf[p++];
      }
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
      Real sign = (nb.polar and flip_across_pole_field[IB1]) ? -1.0 : 1.0;
      for(int i=pmb->is; i<=pmb->ie; i++)
        e1(k,j,i)+=sign*buf[p++];
    }
  }

  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetEMFBoundaryFromFiner(Real *buf, const NeighborBlock& nb)
//  \brief Add up the EMF received from a block on the finer level
//         Later they will be divided in the AverageEMFBoundary function
void BoundaryValues::SetEMFBoundaryFromFiner(Real *buf, const NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_mblock_;
  AthenaArray<Real> &e1=pmb->pfield->e.x1e;
  AthenaArray<Real> &e2=pmb->pfield->e.x2e;
  AthenaArray<Real> &e3=pmb->pfield->e.x3e;
  int p=0;
  if(nb.type==neighbor_face) {
    if(pmb->block_size.nx3 > 1) { // 3D
      // x1 direction
      if(nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
        int i, jl=pmb->js, ju=pmb->je, kl=pmb->ks, ku=pmb->ke;
        if(nb.fid==INNER_X1) i=pmb->is;
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
      else if(nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
        int j, il=pmb->is, iu=pmb->ie, kl=pmb->ks, ku=pmb->ke;
        if(nb.fid==INNER_X2) j=pmb->js;
        else j=pmb->je+1;
        if(nb.fi1==0) iu=pmb->is+pmb->block_size.nx1/2-1;
        else il=pmb->is+pmb->block_size.nx1/2;
        if(nb.fi2==0) ku=pmb->ks+pmb->block_size.nx3/2-1;
        else kl=pmb->ks+pmb->block_size.nx3/2;
        // unpack e1
        Real sign = (nb.polar and flip_across_pole_field[IB1]) ? -1.0 : 1.0;
        for(int k=kl; k<=ku+1; k++) {
          for(int i=il; i<=iu; i++)
            e1(k,j,i)+=sign*buf[p++];
        }
        // unpack e3
        sign = (nb.polar and flip_across_pole_field[IB3]) ? -1.0 : 1.0;
        for(int k=kl; k<=ku; k++) {
          for(int i=il; i<=iu+1; i++)
            e3(k,j,i)+=sign*buf[p++];
        }
      }
      // x3 direction
      else if(nb.fid==INNER_X3 || nb.fid==OUTER_X3) {
        int k, il=pmb->is, iu=pmb->ie, jl=pmb->js, ju=pmb->je;
        if(nb.fid==INNER_X3) k=pmb->ks;
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
      if(nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
        int i, jl=pmb->js, ju=pmb->je;
        if(nb.fid==INNER_X1) i=pmb->is;
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
      else if(nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
        int j, il=pmb->is, iu=pmb->ie;
        if(nb.fid==INNER_X2) j=pmb->js;
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
      if(nb.fid==INNER_X1) i=pmb->is;
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
        Real sign = (nb.polar and flip_across_pole_field[IB3]) ? -1.0 : 1.0;
        for(int k=kl; k<=ku; k++)
          e3(k,j,i)+=sign*buf[p++];
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
        Real sign = (nb.polar and flip_across_pole_field[IB1]) ? -1.0 : 1.0;
        for(int i=il; i<=iu; i++)
          e1(k,j,i)+=sign*buf[p++];
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
//! \fn void BoundaryValues::SetEMFBoundaryPolar(Real **buf_list, int num_bufs,
//          bool north)
//  \brief Overwrite EMF values along polar axis with azimuthal averages
void BoundaryValues::SetEMFBoundaryPolar(Real **buf_list, int num_bufs, bool north)
{
  MeshBlock *pmb = pmy_mblock_;
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
    if(n==INNER_X1 || n==OUTER_X1) {
      if(n==INNER_X1) i=pmb->is;
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
    if(n==INNER_X2 || n==OUTER_X2) {
      if(n==INNER_X2) j=pmb->js;
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
    if(n==INNER_X3 || n==OUTER_X3) {
      if(n==INNER_X3) k=pmb->ks;
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
    if ((pmb->block_bcs[n] != BLOCK_BNDRY) && (pmb->block_bcs[n] != PERIODIC_BNDRY)
        && (pmb->block_bcs[n] != POLAR_BNDRY)) continue;
    if(n==INNER_X1 || n==OUTER_X1) {
      if(n==INNER_X1) i=pmb->is;
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
    if(n==INNER_X2 || n==OUTER_X2) {
      if(n==INNER_X2) j=pmb->js;
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
    if(n==INNER_X3 || n==OUTER_X3) {
      if(n==INNER_X3) k=pmb->ks;
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
    // x1x2 edge (both 2D and 3D)
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
//! \fn void BoundaryValues::PolarSingleEMF(void)
//  \brief single CPU in the azimuthal direction for the polar boundary
void BoundaryValues::PolarSingleEMF(void)
{
  MeshBlock *pmb=pmy_mblock_;
  AthenaArray<Real> &e1=pmb->pfield->e.x1e;
  AthenaArray<Real> &e2=pmb->pfield->e.x2e;
  AthenaArray<Real> &e3=pmb->pfield->e.x3e;

  int i, j, k, nl;
  if(pmb->loc.level == pmb->pmy_mesh->root_level && pmb->pmy_mesh->nrbx3 == 1){
    if(pmb->block_bcs[INNER_X2]==POLAR_BNDRY) {
      j=pmb->js;
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      if(pmb->block_size.nx3 > 1) {
        for(int i=pmb->is; i<=pmb->ie; i++){
          Real tote1=0.0;
          for(int k=pmb->ks; k<=pmb->ke; k++) {
            tote1+=e1(k,j,i);
          }
          Real e1a=tote1/double(pmb->ke-pmb->ks+1);
	  for(int k=pmb->ks; k<=pmb->ke+1; k++) {
            e1(k,j,i)=e1a;
          }
        }
        for(int i=pmb->is; i<=pmb->ie+1; i++){
          for(int k=pmb->ks; k<=pmb->ke; k++) {
            exc_(k)=e3(k,j,i);
          }
          for(int k=pmb->ks; k<=pmb->ke; k++) {
            int k_shift = k;
            k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
            e3(k,j,i)=exc_(k_shift);
          }
        }
      }
    }

    if(pmb->block_bcs[OUTER_X2]==POLAR_BNDRY){
      j=pmb->je+1;
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      if(pmb->block_size.nx3 > 1) {
        for(int i=pmb->is; i<=pmb->ie; i++){
          Real tote1=0.0;
          for (int k=pmb->ks; k<=pmb->ke; ++k) {
            tote1+=e1(k,j,i);
          }
          Real e1a=tote1/double(pmb->ke-pmb->ks+1);
          for (int k=pmb->ks; k<=pmb->ke+1; ++k) {
            e1(k,j,i)=e1a;
          }
        }
        for(int i=pmb->is; i<=pmb->ie+1; i++){
          for(int k=pmb->ks; k<=pmb->ke; k++) {
            exc_(k)=e3(k,j,i);
          }
          for(int k=pmb->ks; k<=pmb->ke; k++) {
            int k_shift = k;
            k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
            e3(k,j,i)=exc_(k_shift);
          }
        }
      }
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

  // Receive same-level non-polar EMF values
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
    if(pmb->pmy_mesh->multilevel==true)
      ClearCoarseEMFBoundary();
    firsttime_[step]=false;
  }

  // Receive finer non-polar EMF values
  if(pmb->pmy_mesh->multilevel==true) {
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
  }

  // Receive polar EMF values
  for (int n = 0; n < num_north_polar_blocks_; ++n) {
    const PolarNeighborBlock &nb = pmb->polar_neighbor_north[n];
    if (emf_north_flag_[step][n] == boundary_waiting) {
      if (nb.rank == Globals::my_rank) { // on the same process
        flag = false;
        continue;
      }
#ifdef MPI_PARALLEL
      else {
        int recv_flag;
        MPI_Test(&req_emf_north_recv_[step][n], &recv_flag, MPI_STATUS_IGNORE);
        if (not recv_flag) {
          flag = false;
          continue;
        }
        emf_north_flag_[step][n] = boundary_arrived;
      }
#endif
    }
  }
  for (int n = 0; n < num_south_polar_blocks_; ++n) {
    const PolarNeighborBlock &nb = pmb->polar_neighbor_south[n];
    if (emf_south_flag_[step][n] == boundary_waiting) {
      if (nb.rank == Globals::my_rank) { // on the same process
        flag = false;
        continue;
      }
#ifdef MPI_PARALLEL
      else {
        int recv_flag;
        MPI_Test(&req_emf_south_recv_[step][n], &recv_flag, MPI_STATUS_IGNORE);
        if (not recv_flag) {
          flag = false;
          continue;
        }
        emf_south_flag_[step][n] = boundary_arrived;
      }
#endif
    }
  }

  if(flag==true){
    AverageEMFBoundary();
    if (num_north_polar_blocks_ > 0)
      SetEMFBoundaryPolar(emf_north_recv_[step], num_north_polar_blocks_, true);
    for (int n = 0; n < num_north_polar_blocks_; ++n)
      emf_north_flag_[step][n] = boundary_completed;
    if (num_south_polar_blocks_ > 0)
      SetEMFBoundaryPolar(emf_south_recv_[step], num_south_polar_blocks_, false);
    for (int n = 0; n < num_south_polar_blocks_; ++n)
      emf_south_flag_[step][n] = boundary_completed;
    if (pmb->block_bcs[INNER_X2]==POLAR_BNDRY||pmb->block_bcs[OUTER_X2]==POLAR_BNDRY)
      PolarSingleEMF();
  }
  return flag;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ClearBoundaryForInit(void)
//  \brief clean up the boundary flags for initialization
void BoundaryValues::ClearBoundaryForInit(void)
{
  MeshBlock *pmb=pmy_mblock_;

  // Note step==0 corresponds to initial exchange of conserved variables, while step==1
  // corresponds to primitives sent only in the case of GR with refinement
  for(int n=0;n<pmb->nneighbor;n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    hydro_flag_[0][nb.bufid] = boundary_waiting;
    if (MAGNETIC_FIELDS_ENABLED)
      field_flag_[0][nb.bufid] = boundary_waiting;
    if (GENERAL_RELATIVITY and pmb->pmy_mesh->multilevel)
      hydro_flag_[1][nb.bufid] = boundary_waiting;
#ifdef MPI_PARALLEL
    if(nb.rank!=Globals::my_rank) {
      MPI_Wait(&req_hydro_send_[0][nb.bufid],MPI_STATUS_IGNORE); // Wait for Isend
      if (MAGNETIC_FIELDS_ENABLED)
        MPI_Wait(&req_field_send_[0][nb.bufid],MPI_STATUS_IGNORE); // Wait for Isend
      if (GENERAL_RELATIVITY and pmb->pmy_mesh->multilevel)
        MPI_Wait(&req_hydro_send_[1][nb.bufid],MPI_STATUS_IGNORE); // Wait for Isend
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

  // Clear non-polar boundary communications
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

  // Clear polar boundary communications
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int l = 0; l < NSTEP; ++l) {
      for (int n = 0; n < num_north_polar_blocks_; ++n) {
        PolarNeighborBlock &nb = pmb->polar_neighbor_north[n];
        emf_north_flag_[l][n] = boundary_waiting;
#ifdef MPI_PARALLEL
        if(nb.rank != Globals::my_rank)
          MPI_Wait(&req_emf_north_send_[l][n], MPI_STATUS_IGNORE);
#endif
      }
      for (int n = 0; n < num_south_polar_blocks_; ++n) {
        PolarNeighborBlock &nb = pmb->polar_neighbor_south[n];
        emf_south_flag_[l][n] = boundary_waiting;
#ifdef MPI_PARALLEL
        if(nb.rank != Globals::my_rank)
          MPI_Wait(&req_emf_south_send_[l][n], MPI_STATUS_IGNORE);
#endif
      }
    }
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ApplyPhysicalBoundaries(AthenaArray<Real> &pdst,
//           AthenaArray<Real> &cdst, FaceField &bfdst, AthenaArray<Real> &bcdst)
//                                                   FaceField &bdst)
//  \brief Apply all the physical boundary conditions for both hydro and field
void BoundaryValues::ApplyPhysicalBoundaries(AthenaArray<Real> &pdst,
     AthenaArray<Real> &cdst, FaceField &bfdst, AthenaArray<Real> &bcdst)
{
  MeshBlock *pmb=pmy_mblock_;
  Coordinates *pco=pmb->pcoord;
  int bis=pmb->is, bie=pmb->ie, bjs=pmb->js, bje=pmb->je, bks=pmb->ks, bke=pmb->ke;
  if(pmb->pmy_mesh->face_only==false) { // extend the ghost zone
    bis=pmb->is-NGHOST;
    bie=pmb->ie+NGHOST;
    if(BoundaryFunction_[INNER_X2]==NULL && pmb->block_size.nx2>1) bjs=pmb->js-NGHOST;
    if(BoundaryFunction_[OUTER_X2]==NULL && pmb->block_size.nx2>1) bje=pmb->je+NGHOST;
    if(BoundaryFunction_[INNER_X3]==NULL && pmb->block_size.nx3>1) bks=pmb->ks-NGHOST;
    if(BoundaryFunction_[OUTER_X3]==NULL && pmb->block_size.nx3>1) bke=pmb->ke+NGHOST;
  }
  // Apply boundary function on inner-x1
  if (BoundaryFunction_[INNER_X1] != NULL) {
    BoundaryFunction_[INNER_X1](pmb, pco, pdst, bfdst, pmb->is, pmb->ie, bjs,bje,bks,bke);
    if(MAGNETIC_FIELDS_ENABLED) {
      pmb->pfield->CalculateCellCenteredField(bfdst, bcdst, pco,
        pmb->is-NGHOST, pmb->is-1, bjs, bje, bks, bke);
    }
    pmb->phydro->peos->PrimitiveToConserved(pdst, bcdst, cdst, pco,
      pmb->is-NGHOST, pmb->is-1, bjs, bje, bks, bke);
  }

  // Apply boundary function on outer-x1
  if (BoundaryFunction_[OUTER_X1] != NULL) {
    BoundaryFunction_[OUTER_X1](pmb, pco, pdst, bfdst, pmb->is, pmb->ie, bjs,bje,bks,bke);
    if(MAGNETIC_FIELDS_ENABLED) {
      pmb->pfield->CalculateCellCenteredField(bfdst, bcdst, pco,
        pmb->ie+1, pmb->ie+NGHOST, bjs, bje, bks, bke);
    }
    pmb->phydro->peos->PrimitiveToConserved(pdst, bcdst, cdst, pco,
      pmb->ie+1, pmb->ie+NGHOST, bjs, bje, bks, bke);
  }

  if(pmb->block_size.nx2>1) { // 2D or 3D

    // Apply boundary function on inner-x2
    if (BoundaryFunction_[INNER_X2] != NULL) {
      BoundaryFunction_[INNER_X2](pmb, pco, pdst, bfdst, bis,bie, pmb->js, pmb->je, bks,bke);
      if(MAGNETIC_FIELDS_ENABLED) {
        pmb->pfield->CalculateCellCenteredField(bfdst, bcdst, pco,
          bis, bie, pmb->js-NGHOST, pmb->js-1, bks, bke);
      }
      pmb->phydro->peos->PrimitiveToConserved(pdst, bcdst, cdst, pco,
        bis, bie, pmb->js-NGHOST, pmb->js-1, bks, bke);
    }

    // Apply boundary function on outer-x2
    if (BoundaryFunction_[OUTER_X2] != NULL) {
      BoundaryFunction_[OUTER_X2](pmb, pco, pdst, bfdst, bis,bie, pmb->js, pmb->je, bks,bke);
      if(MAGNETIC_FIELDS_ENABLED) {
        pmb->pfield->CalculateCellCenteredField(bfdst, bcdst, pco,
          bis, bie, pmb->je+1, pmb->je+NGHOST, bks, bke);
      }
      pmb->phydro->peos->PrimitiveToConserved(pdst, bcdst, cdst, pco,
        bis, bie, pmb->je+1, pmb->je+NGHOST, bks, bke);
    }
  }

  if(pmb->block_size.nx3>1) { // 3D
    if(pmb->pmy_mesh->face_only==false) {
      bjs=pmb->js-NGHOST;
      bje=pmb->je+NGHOST;
    }

    // Apply boundary function on inner-x3
    if (BoundaryFunction_[INNER_X3] != NULL) {
      BoundaryFunction_[INNER_X3](pmb, pco, pdst, bfdst, bis,bie,bjs,bje, pmb->ks, pmb->ke);
      if(MAGNETIC_FIELDS_ENABLED) {
        pmb->pfield->CalculateCellCenteredField(bfdst, bcdst, pco,
          bis, bie, bjs, bje, pmb->ks-NGHOST, pmb->ks-1);
      }
      pmb->phydro->peos->PrimitiveToConserved(pdst, bcdst, cdst, pco,
        bis, bie, bjs, bje, pmb->ks-NGHOST, pmb->ks-1);
    }

    // Apply boundary function on outer-x3
    if (BoundaryFunction_[OUTER_X3] != NULL) {
      BoundaryFunction_[OUTER_X3](pmb, pco, pdst, bfdst, bis,bie,bjs,bje, pmb->ks, pmb->ke);
      if(MAGNETIC_FIELDS_ENABLED) {
        pmb->pfield->CalculateCellCenteredField(bfdst, bcdst, pco,
          bis, bie, bjs, bje, pmb->ke+1, pmb->ke+NGHOST);
      }
      pmb->phydro->peos->PrimitiveToConserved(pdst, bcdst, cdst, pco,
        bis, bie, bjs, bje, pmb->ke+1, pmb->ke+NGHOST);
    }
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

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ProlongateBoundaries(AthenaArray<Real> &pdst,
//           AthenaArray<Real> &cdst, FaceField &bdst, AthenaArray<Real> &bcdst)
//  \brief Prolongate the level boundary using the coarse data
void BoundaryValues::ProlongateBoundaries(AthenaArray<Real> &pdst,
     AthenaArray<Real> &cdst, FaceField &bfdst, AthenaArray<Real> &bcdst)
{
  MeshBlock *pmb=pmy_mblock_;
  MeshRefinement *pmr=pmb->pmr;
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
          // skip myself or coarse levels; only the same level must be restricted
          if(ntype==0 || pmb->nblevel[nk+1][nj+1][ni+1]!=mylevel) continue;

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

          pmb->pmr->RestrictCellCenteredValues(cdst, pmr->coarse_cons_, 0, NHYDRO-1,
                                               ris, rie, rjs, rje, rks, rke);
          if (GENERAL_RELATIVITY)
            pmb->pmr->RestrictCellCenteredValues(pdst, pmr->coarse_prim_, 0, NHYDRO-1,
                                                 ris, rie, rjs, rje, rks, rke);
          if (MAGNETIC_FIELDS_ENABLED) {
            int rs=ris, re=rie+1;
            if(rs==pmb->cis   && pmb->nblevel[nk+1][nj+1][ni  ]<mylevel) rs++;
            if(re==pmb->cie+1 && pmb->nblevel[nk+1][nj+1][ni+2]<mylevel) re--;
            pmr->RestrictFieldX1(bfdst.x1f, pmr->coarse_b_.x1f, rs, re, rjs, rje, rks, rke);
            if(pmb->block_size.nx2 > 1) {
              rs=rjs, re=rje+1;
              if(rs==pmb->cjs   && pmb->nblevel[nk+1][nj  ][ni+1]<mylevel) rs++;
              if(re==pmb->cje+1 && pmb->nblevel[nk+1][nj+2][ni+1]<mylevel) re--;
              pmr->RestrictFieldX2(bfdst.x2f, pmr->coarse_b_.x2f, ris, rie, rs, re, rks, rke);
            }
            else 
              pmr->RestrictFieldX2(bfdst.x2f, pmr->coarse_b_.x2f, ris, rie, rjs, rje, rks, rke);
            if(pmb->block_size.nx3 > 1) {
              rs=rks, re=rke+1;
              if(rs==pmb->cks   && pmb->nblevel[nk  ][nj+1][ni+1]<mylevel) rs++;
              if(re==pmb->cke+1 && pmb->nblevel[nk+2][nj+1][ni+1]<mylevel) re--;
              pmr->RestrictFieldX3(bfdst.x3f, pmr->coarse_b_.x3f, ris, rie, rjs, rje, rs, re);
            }
            else
              pmr->RestrictFieldX3(bfdst.x3f, pmr->coarse_b_.x3f, ris, rie, rjs, rje, rks, rke);
          }
        }
      }
    }


    // calculate the loop limits for the ghost zones
    int cn = (NGHOST+1)/2;
    int si, ei, sj, ej, sk, ek, fsi, fei, fsj, fej, fsk, fek;
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

    // convert the ghost zone and ghost-ghost zones into primitive variables
    // this includes cell-centered field calculation
    int f1m=0, f1p=0, f2m=0, f2p=0, f3m=0, f3p=0;
    if(nb.ox1==0) {
      if(pmb->nblevel[1][1][0]!=-1) f1m=1;
      if(pmb->nblevel[1][1][2]!=-1) f1p=1;
    }
    else f1m=1, f1p=1;
    if(pmb->block_size.nx2>1) {
      if(nb.ox2==0) {
        if(pmb->nblevel[1][0][1]!=-1) f2m=1;
        if(pmb->nblevel[1][2][1]!=-1) f2p=1;
      }
      else f2m=1, f2p=1;
    }
    if(pmb->block_size.nx3>1) {
      if(nb.ox3==0) {
        if(pmb->nblevel[0][1][1]!=-1) f3m=1;
        if(pmb->nblevel[2][1][1]!=-1) f3p=1;
      }
      else f3m=1, f3p=1;
    }
    pmb->phydro->peos->ConservedToPrimitive(pmr->coarse_cons_, pmr->coarse_prim_,
                 pmr->coarse_b_, pmr->coarse_prim_, pmr->coarse_bcc_, pmr->pcoarsec,
                 si-f1m, ei+f1p, sj-f2m, ej+f2p, sk-f3m, ek+f3p);

    // Apply physical boundaries
    if(nb.ox1==0) {
      if(BoundaryFunction_[INNER_X1]!=NULL) {
        BoundaryFunction_[INNER_X1](pmb, pmr->pcoarsec, pmr->coarse_prim_,
                          pmr->coarse_b_, pmb->cis, pmb->cie, sj, ej, sk, ek);
      }
      if(BoundaryFunction_[OUTER_X1]!=NULL) {
        BoundaryFunction_[OUTER_X1](pmb, pmr->pcoarsec, pmr->coarse_prim_,
                          pmr->coarse_b_, pmb->cis, pmb->cie, sj, ej, sk, ek);
      }
    }
    if(nb.ox2==0 && pmb->block_size.nx2 > 1) {
      if(BoundaryFunction_[INNER_X2]!=NULL) {
        BoundaryFunction_[INNER_X2](pmb, pmr->pcoarsec, pmr->coarse_prim_,
                          pmr->coarse_b_, si, ei, pmb->cjs, pmb->cje, sk, ek);
      }
      if(BoundaryFunction_[OUTER_X2]!=NULL) {
        BoundaryFunction_[OUTER_X2](pmb, pmr->pcoarsec, pmr->coarse_prim_,
                          pmr->coarse_b_, si, ei, pmb->cjs, pmb->cje, sk, ek);
      }
    }
    if(nb.ox3==0 && pmb->block_size.nx3 > 1) {
      if(BoundaryFunction_[INNER_X3]!=NULL) {
        BoundaryFunction_[INNER_X3](pmb, pmr->pcoarsec, pmr->coarse_prim_,
                          pmr->coarse_b_, si, ei, sj, ej, pmb->cks, pmb->cke);
      }
      if(BoundaryFunction_[OUTER_X3]!=NULL) {
        BoundaryFunction_[OUTER_X3](pmb, pmr->pcoarsec, pmr->coarse_prim_,
                          pmr->coarse_b_, si, ei, sj, ej, pmb->cks, pmb->cke);
      }
    }

    // now that the ghost-ghost zones are filled
    // calculate the loop limits for the finer grid
    fsi=(si-pmb->cis)*2+pmb->is,   fei=(ei-pmb->cis)*2+pmb->is+1;
    if(pmb->block_size.nx2 > 1)
      fsj=(sj-pmb->cjs)*2+pmb->js, fej=(ej-pmb->cjs)*2+pmb->js+1;
    else fsj=pmb->js, fej=pmb->je;
    if(pmb->block_size.nx3 > 1)
      fsk=(sk-pmb->cks)*2+pmb->ks, fek=(ek-pmb->cks)*2+pmb->ks+1;
    else fsk=pmb->ks, fek=pmb->ke;

    // prolongate hydro variables using primitive
    pmr->ProlongateCellCenteredValues(pmr->coarse_prim_, pdst, 0, NHYDRO-1,
                                      si, ei, sj, ej, sk, ek);
    // prollongate magnetic fields
    if (MAGNETIC_FIELDS_ENABLED) {
      int il, iu, jl, ju, kl, ku;
      il=si, iu=ei+1;
      if((nb.ox1>=0) && (pmb->nblevel[nb.ox3+1][nb.ox2+1][nb.ox1  ]>=mylevel)) il++;
      if((nb.ox1<=0) && (pmb->nblevel[nb.ox3+1][nb.ox2+1][nb.ox1+2]>=mylevel)) iu--;
      if(pmb->block_size.nx2 > 1) {
        jl=sj, ju=ej+1;
        if((nb.ox2>=0) && (pmb->nblevel[nb.ox3+1][nb.ox2  ][nb.ox1+1]>=mylevel)) jl++;
        if((nb.ox2<=0) && (pmb->nblevel[nb.ox3+1][nb.ox2+2][nb.ox1+1]>=mylevel)) ju--;
      }
      else jl=sj, ju=ej;
      if(pmb->block_size.nx3 > 1) {
        kl=sk, ku=ek+1;
        if((nb.ox3>=0) && (pmb->nblevel[nb.ox3  ][nb.ox2+1][nb.ox1+1]>=mylevel)) kl++;
        if((nb.ox3<=0) && (pmb->nblevel[nb.ox3+2][nb.ox2+1][nb.ox1+1]>=mylevel)) ku--;
      }
      else kl=sk, ku=ek;

      // step 1. calculate x1 outer surface fields and slopes
      pmr->ProlongateSharedFieldX1(pmr->coarse_b_.x1f, bfdst.x1f, il, iu, sj, ej, sk, ek);
      // step 2. calculate x2 outer surface fields and slopes
      pmr->ProlongateSharedFieldX2(pmr->coarse_b_.x2f, bfdst.x2f, si, ei, jl, ju, sk, ek);
      // step 3. calculate x3 outer surface fields and slopes
      pmr->ProlongateSharedFieldX3(pmr->coarse_b_.x3f, bfdst.x3f, si, ei, sj, ej, kl, ku);
      // step 4. calculate the internal finer fields using the Toth & Roe method
      pmr->ProlongateInternalField(bfdst, si, ei, sj, ej, sk, ek);

      // Field prolongation completed, calculate cell centered fields
      pmb->pfield->CalculateCellCenteredField(bfdst, bcdst, pmb->pcoord,
                                              fsi, fei, fsj, fej, fsk, fek);
    }
    // calculate conservative variables
    pmb->phydro->peos->PrimitiveToConserved(pdst, bcdst, cdst, pmb->pcoord,
                                            fsi, fei, fsj, fej, fsk, fek);
  }
}
