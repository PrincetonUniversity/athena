//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_var.cpp
//  \brief constructor/destructor and default implementations for some functions in the
//         abstract BoundaryVariable class

// C headers

// C++ headers
// #include <algorithm>  // min
// #include <cmath>
// #include <cstdlib>
// #include <cstring>    // std::memcpy
// #include <iomanip>
#include <iostream>   // endl
// #include <limits>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
// #include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
// #include "../coordinates/coordinates.hpp"
// #include "../eos/eos.hpp"
// #include "../field/field.hpp"
// #include "../globals.hpp"
// #include "../gravity/mg_gravity.hpp"
// #include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
// #include "../mesh/mesh_refinement.hpp"
// #include "../multigrid/multigrid.hpp"
// #include "../parameter_input.hpp"
// #include "../utils/buffer_utils.hpp"
#include "bvals_interfaces.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif


// constructor

BoundaryVariable::BoundaryVariable(
    MeshBlock *pmb, BoundaryValues *pbval, enum BoundaryType type)
    : BoundaryVariable() {
  // be sure to add this new BoundaryVariable object as the last node in the linked list
  // (or wherever)
}

// destructor

BoundaryVariable::~BoundaryVariable() {
  MeshBlock *pmb=pmy_block_;

  DestroyBoundaryData(bd_var_);
  if (pmb->pmy_mesh->multilevel==true) // SMR or AMR
    DestroyBoundaryData(bd_var_flcor_);
}


//----------------------------------------------------------------------------------------
//! \fn void BoundaryVariable::InitBoundaryData(BoundaryData &bd, enum BoundaryType type)
//  \brief Initialize BoundaryData structure

void BoundaryVariable::InitBoundaryData(BoundaryData &bd, enum BoundaryType type) {
  MeshBlock *pmb=pmy_block_;
  bool multilevel=pmy_mesh_->multilevel;
  NeighborIndexes *ni=pbval->ni;

  int f2d=0, f3d=0;
  int cng, cng1, cng2, cng3;
  if (pmb->block_size.nx2 > 1) f2d=1;
  if (pmb->block_size.nx3 > 1) f3d=1;
  cng=cng1=pmb->cnghost;
  cng2=cng*f2d;
  cng3=cng*f3d;
  int size=0;
  bd.nbmax=pbval->maxneighbor_;
  if (type==BNDRY_FLCOR || type==BNDRY_EMFCOR) {
    for (bd.nbmax=0; pbval->ni[bd.nbmax].type==NEIGHBOR_FACE; bd.nbmax++) {}
  }
  if (type==BNDRY_EMFCOR) {
    for (          ; pbval->ni[bd.nbmax].type==NEIGHBOR_EDGE; bd.nbmax++) {}
  }
  for (int n=0; n<bd.nbmax; n++) {
    // Clear flags and requests
    bd.flag[n]=BNDRY_WAITING;
    bd.send[n]=nullptr;
    bd.recv[n]=nullptr;
#ifdef MPI_PARALLEL
    bd.req_send[n]=MPI_REQUEST_NULL;
    bd.req_recv[n]=MPI_REQUEST_NULL;
#endif

    // Allocate buffers
    // calculate the buffer size
    switch(type) {
      case BNDRY_HYDRO: {
        size=((ni[n].ox1==0)?pmb->block_size.nx1:NGHOST)
             *((ni[n].ox2==0)?pmb->block_size.nx2:NGHOST)
             *((ni[n].ox3==0)?pmb->block_size.nx3:NGHOST);
        if (multilevel) {
          int f2c=((ni[n].ox1==0) ? ((pmb->block_size.nx1+1)/2):NGHOST)
                  *((ni[n].ox2==0) ? ((pmb->block_size.nx2+1)/2):NGHOST)
                  *((ni[n].ox3==0) ? ((pmb->block_size.nx3+1)/2):NGHOST);
          int c2f=((ni[n].ox1==0) ?((pmb->block_size.nx1+1)/2+cng1):cng)
                  *((ni[n].ox2==0)?((pmb->block_size.nx2+1)/2+cng2):cng)
                  *((ni[n].ox3==0)?((pmb->block_size.nx3+1)/2+cng3):cng);
          size=std::max(size,c2f);
          size=std::max(size,f2c);
        }
        size*=NHYDRO; // KGF: need a generic BNDRY_CC counterpart
      }
        break;
      case BNDRY_FC: {
        int size1=((ni[n].ox1==0) ? (pmb->block_size.nx1+1):NGHOST)
                  *((ni[n].ox2==0) ? (pmb->block_size.nx2):NGHOST)
                  *((ni[n].ox3==0) ? (pmb->block_size.nx3):NGHOST);
        int size2=((ni[n].ox1==0) ? (pmb->block_size.nx1):NGHOST)
                  *((ni[n].ox2==0) ? (pmb->block_size.nx2+f2d):NGHOST)
                  *((ni[n].ox3==0) ? (pmb->block_size.nx3):NGHOST);
        int size3=((ni[n].ox1==0) ? (pmb->block_size.nx1):NGHOST)
                  *((ni[n].ox2==0) ? (pmb->block_size.nx2):NGHOST)
                  *((ni[n].ox3==0) ? (pmb->block_size.nx3+f3d):NGHOST);
        size=size1+size2+size3;
        if (multilevel) {
          if (ni[n].type!=NEIGHBOR_FACE) {
            if (ni[n].ox1!=0) size1=size1/NGHOST*(NGHOST+1);
            if (ni[n].ox2!=0) size2=size2/NGHOST*(NGHOST+1);
            if (ni[n].ox3!=0) size3=size3/NGHOST*(NGHOST+1);
          }
          size=size1+size2+size3;
          int f2c1=((ni[n].ox1==0) ? ((pmb->block_size.nx1+1)/2+1):NGHOST)
                   *((ni[n].ox2==0) ? ((pmb->block_size.nx2+1)/2):NGHOST)
                   *((ni[n].ox3==0) ? ((pmb->block_size.nx3+1)/2):NGHOST);
          int f2c2=((ni[n].ox1==0) ? ((pmb->block_size.nx1+1)/2):NGHOST)
                   *((ni[n].ox2==0) ? ((pmb->block_size.nx2+1)/2+f2d)
                     : NGHOST)
                   *((ni[n].ox3==0) ? ((pmb->block_size.nx3+1)/2):NGHOST);
          int f2c3=((ni[n].ox1==0) ? ((pmb->block_size.nx1+1)/2):NGHOST)
                   *((ni[n].ox2==0) ? ((pmb->block_size.nx2+1)/2):NGHOST)
                   *((ni[n].ox3==0) ? ((pmb->block_size.nx3+1)/2+f3d)
                     : NGHOST);
          if (ni[n].type!=NEIGHBOR_FACE) {
            if (ni[n].ox1!=0) f2c1=f2c1/NGHOST*(NGHOST+1);
            if (ni[n].ox2!=0) f2c2=f2c2/NGHOST*(NGHOST+1);
            if (ni[n].ox3!=0) f2c3=f2c3/NGHOST*(NGHOST+1);
          }
          int fsize=f2c1+f2c2+f2c3;
          int c2f1=
              ((ni[n].ox1==0) ? ((pmb->block_size.nx1+1)/2+cng1+1):cng+1)
              *((ni[n].ox2==0) ? ((pmb->block_size.nx2+1)/2+cng2):cng)
              *((ni[n].ox3==0) ? ((pmb->block_size.nx3+1)/2+cng3):cng);
          int c2f2=
              ((ni[n].ox1==0) ?((pmb->block_size.nx1+1)/2+cng1):cng)
              *((ni[n].ox2==0)?((pmb->block_size.nx2+1)/2+cng2+f2d):cng+1)
              *((ni[n].ox3==0)?((pmb->block_size.nx3+1)/2+cng3):cng);
          int c2f3=
              ((ni[n].ox1==0) ? ((pmb->block_size.nx1+1)/2+cng1):cng)
              *((ni[n].ox2==0)?((pmb->block_size.nx2+1)/2+cng2):cng)
              *((ni[n].ox3==0) ?
                ((pmb->block_size.nx3+1)/2+cng3+f3d) : cng+1);
          int csize=c2f1+c2f2+c2f3;
          size=std::max(size,std::max(csize,fsize));
        }
      }
        break;
      case BNDRY_FLCOR: {
        if (ni[n].ox1!=0)
          size=(pmb->block_size.nx2+1)/2*(pmb->block_size.nx3+1)/2*NHYDRO;
        if (ni[n].ox2!=0)
          size=(pmb->block_size.nx1+1)/2*(pmb->block_size.nx3+1)/2*NHYDRO;
        if (ni[n].ox3!=0)
          size=(pmb->block_size.nx1+1)/2*(pmb->block_size.nx2+1)/2*NHYDRO;
      }
        break;
      case BNDRY_EMFCOR: {
        if (ni[n].type==NEIGHBOR_FACE) {
          if (pmb->block_size.nx3>1) { // 3D
            if (ni[n].ox1!=0)
              size=(pmb->block_size.nx2+1)*(pmb->block_size.nx3)
                   +(pmb->block_size.nx2)*(pmb->block_size.nx3+1);
            else if (ni[n].ox2!=0)
              size=(pmb->block_size.nx1+1)*(pmb->block_size.nx3)
                   +(pmb->block_size.nx1)*(pmb->block_size.nx3+1);
            else
              size=(pmb->block_size.nx1+1)*(pmb->block_size.nx2)
                   +(pmb->block_size.nx1)*(pmb->block_size.nx2+1);
          } else if (pmb->block_size.nx2>1) { // 2D
            if (ni[n].ox1!=0)
              size=(pmb->block_size.nx2+1)+pmb->block_size.nx2;
            else
              size=(pmb->block_size.nx1+1)+pmb->block_size.nx1;
          } else { // 1D
            size=2;
          }
        } else if (ni[n].type==NEIGHBOR_EDGE) {
          if (pmb->block_size.nx3>1) { // 3D
            if (ni[n].ox3==0) size=pmb->block_size.nx3;
            if (ni[n].ox2==0) size=pmb->block_size.nx2;
            if (ni[n].ox1==0) size=pmb->block_size.nx1;
          } else if (pmb->block_size.nx2>1) {
            size=1;
          }
        }
      }
        break;
      default: {
        std::stringstream msg;
        msg << "### FATAL ERROR in InitBoundaryData" << std::endl
            << "Invalid boundary type is specified." << std::endl;
        ATHENA_ERROR(msg);
      }
        break;
    }
    bd.send[n]=new Real[size];
    bd.recv[n]=new Real[size];
  }
}


//----------------------------------------------------------------------------------------
//! \fn void BoundaryVariable::DestroyBoundaryData(BoundaryData &bd)
//  \brief Destroy BoundaryData structure
void BoundaryVariable::DestroyBoundaryData(BoundaryData &bd) {
  for (int n=0; n<bd.nbmax; n++) {
    delete [] bd.send[n];
    delete [] bd.recv[n];
#ifdef MPI_PARALLEL
    if (bd.req_send[n]!=MPI_REQUEST_NULL)
      MPI_Request_free(&bd.req_send[n]);
    if (bd.req_recv[n]!=MPI_REQUEST_NULL)
      MPI_Request_free(&bd.req_recv[n]);
#endif
  }
}
