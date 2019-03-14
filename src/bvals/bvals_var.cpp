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
#include <cstring>    // std::memcpy
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

BoundaryVariable::BoundaryVariable(MeshBlock *pmb) {
  // KGF: what is the point of cyclically setting class member pmy_block_=pmb, then:
  // MeshBlock *pmb=pmy_block_;
  pmy_block_ = pmb;
  pbval_ = pmb->pbval;
  pmy_mesh_ = pmb->pmy_mesh;

  // btype_ = type;

  // KGF: should we even initialize this unsigned integer type data member?
  bvar_index = 0;

  // KGF: not sure of the value of passing around BoundaryQuantity type anymore.
  // Currently, it is always tied to the BoundaryVariable subclass type

  // Need to access both BoundaryData in generic BoundaryVariable::Copy*SameProcess()
  // KGF: option 1: storing bd_var_, bd_var_cc_ in Base class, not initializing them

  // KGF option 2: only store references to BoundaryData object
  // pbd_var_=nullptr;
  // pbd_var_flcor_=nullptr;

  // KGF option 3: initialize somewhat-generic bd_var_ in Base class, store ptr to
  // somewhat-optional flux-correction BoundaryData (unused for unrefined Hydro)
  // InitBoundaryData(bd_var_, type);   // call DestroyBoundaryData in Base dtor
  // pbd_var_flcor_=nullptr;

  // User must add this new BoundaryVariable object to the array in BoundaryValues
  // Choosing to not auto-enroll this new object in the array via its constructor
}

// destructor
// BoundaryVariable::~BoundaryVariable() {
// }

//----------------------------------------------------------------------------------------
//! \fn void BoundaryVariable::InitBoundaryData(BoundaryData &bd, BoundaryQuantity type)
//  \brief Initialize BoundaryData structure

void BoundaryVariable::InitBoundaryData(BoundaryData &bd, BoundaryQuantity type) {
  MeshBlock *pmb=pmy_block_;
  NeighborIndexes *ni=pbval_->ni;
  int cng = pmb->cnghost;
  int size=0;

  bd.nbmax=pbval_->maxneighbor_;
  // KGF: what is happening in the next two conditionals??
  // they are preventing the elimination of "BoundaryQuantity type" function parameter in
  // favor of a simpler boolean switch
  if (type==BoundaryQuantity::cc_FLCOR || type==BoundaryQuantity::fc_FLCOR) {
    for (bd.nbmax=0; pbval_->ni[bd.nbmax].type==NeighborConnect::face; bd.nbmax++) {}
  }
  if (type==BoundaryQuantity::fc_FLCOR) {
    for (          ; pbval_->ni[bd.nbmax].type==NeighborConnect::edge; bd.nbmax++) {}
  }
  for (int n=0; n<bd.nbmax; n++) {
    // Clear flags and requests
    bd.flag[n]=BoundaryStatus::waiting;
    bd.send[n]=nullptr;
    bd.recv[n]=nullptr;
#ifdef MPI_PARALLEL
    bd.req_send[n]=MPI_REQUEST_NULL;
    bd.req_recv[n]=MPI_REQUEST_NULL;
#endif
    // Allocate buffers, calculating the buffer size (variable vs. flux correction)
    if (type ==BoundaryQuantity::cc || type==BoundaryQuantity::fc) {
      size = this->ComputeVariableBufferSize(ni[n], cng);
    } else if (type ==BoundaryQuantity::cc_FLCOR || type==BoundaryQuantity::fc_FLCOR) {
      size = this->ComputeFluxCorrectionBufferSize(ni[n], cng);
    } else {
      std::stringstream msg;
      msg << "### FATAL ERROR in InitBoundaryData" << std::endl
          << "Invalid boundary type is specified." << std::endl;
      ATHENA_ERROR(msg);
    }
    // KGF: original switch statement dependencies on local variableS:
    // switch(type) { // local variables as input: ni[n],
    //case BoundaryQuantity::cc: {  // local variables as input: cng, ..., cng3
    //case BoundaryQuantity::fc: { // local variables as input: cng, ..., cng3, f2d, f3d
    //case BoundaryQuantity::cc_FLCOR: { // local variables as input: NONE
    //case BoundaryQuantity::fc_FLCOR: { // local variables as input: NONE

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



//----------------------------------------------------------------------------------------
//! \fn void BoundaryVariable::CopyVariableBufferSameProcess(NeighborBlock& nb, int ssize)
//  \brief

//  Called in BoundaryVariable::SendBoundaryBuffer(), SendFluxCorrection() calls when the
//  destination neighbor block is on the same MPI rank as the sending MeshBlcok. So
//  std::memcpy() call requires pointer to "void *dst" corresponding to
//  bd_var_.recv[nb.targetid] in separate BoundaryVariable object in separate vector in
//  separate BoundaryValues

void BoundaryVariable::CopyVariableBufferSameProcess(NeighborBlock& nb, int ssize) {
  // Locate target buffer
  // 1) which MeshBlock?
  MeshBlock *ptarget_block=pmy_mesh_->FindMeshBlock(nb.gid);
  // 2) which element in vector of BoundaryVariable *?
  BoundaryData *ptarget_bdata = &(ptarget_block->pbval->bvars[bvar_index]->bd_var_);
  // KGF: hardcoded assumptions that bvar_index is always up-to-date, and that the same
  // vector elements and ordering exist for all MeshBlocks/BoundaryValues objects

  // KGF: add check that bvar_index is initialized and valid. E.g. the mistake I made with
  // FFT self-gravity for which pgbval was never added to bvars, and bvar_index was unset,
  // breaking this function

  std::memcpy(ptarget_bdata->recv[nb.targetid], bd_var_.send[nb.bufid],
              ssize*sizeof(Real));
  // finally, set the BoundaryStatus flag on the destination buffer
  ptarget_bdata->flag[nb.targetid]=BoundaryStatus::arrived;
  return;
}

void BoundaryVariable::CopyFluxCorrectionBufferSameProcess(NeighborBlock& nb, int ssize) {
  // Locate target buffer
  // 1) which MeshBlock?
  MeshBlock *ptarget_block=pmy_mesh_->FindMeshBlock(nb.gid);
  // 2) which element in vector of BoundaryVariable *?
  BoundaryData *ptarget_bdata = &(ptarget_block->pbval->bvars[bvar_index]->bd_var_flcor_);
  std::memcpy(ptarget_bdata->recv[nb.targetid], bd_var_flcor_.send[nb.bufid],
              ssize*sizeof(Real));
  ptarget_bdata->flag[nb.targetid]=BoundaryStatus::arrived;
  return;
}
