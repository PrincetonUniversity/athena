//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_sixray.cpp
//! \brief functions that apply BCs for six-ray column density variables

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>    // memcpy()
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../eos/eos.hpp"
#include "../../field/field.hpp"
#include "../../globals.hpp"
#include "../../hydro/hydro.hpp"
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"
#include "../../utils/buffer_utils.hpp"
#include "../bvals.hpp"
#include "bvals_sixray.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

//! constructor
//
SixRayBoundaryVariable::SixRayBoundaryVariable(MeshBlock *pmb, AthenaArray<Real> *var)
      : BoundaryVariable(pmb), var(var), nu_(var->GetDim4() - 1) {
  //only take 4 dimention array
  if (nu_ < 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in SixRayBoundaryVariable constructor" << std::endl
        << "An 'AthenaArray<Real> *var' of nx4_ = " << var->GetDim4() << " was passed\n"
        << "Should be nx4 >= 1 (likely uninitialized)." << std::endl;
    ATHENA_ERROR(msg);
  }

  //This will initialize the maximum neighbours similar to the cell-centered variable.
  //Leaving it for now for potential extension to mesh refinement in the future.
  InitBoundaryData(bd_var_, BoundaryQuantity::cc);

  if ((pmy_mesh_->multilevel)
      || (pbval_->shearing_box != 0)) { // SMR or AMR or SHEARING_BOX
    std::stringstream msg;
    msg << "### FATAL ERROR in SixRayBoundaryVariable constructor" << std::endl
        << "Not yet compatible with mesh refinement or shearingbox" << std::endl;
    ATHENA_ERROR(msg);
  }
}

SixRayBoundaryVariable::~SixRayBoundaryVariable() {
  DestroyBoundaryData(bd_var_);
}

//----------------------------------------------------------------------------------------
//! \fn int SixRayBoundaryVariable::ComputeVariableBufferSize(
//!     const NeighborIndexes& ni, int cng)
//! \brief
int SixRayBoundaryVariable::ComputeVariableBufferSize(const NeighborIndexes& ni, int cng) {
  int size;
  MeshBlock *pmb = pmy_block_;

  if (ni.type == NeighborConnect::face) {
    size = ((ni.ox1 == 0) ? pmb->block_size.nx1 : 1)
           *((ni.ox2 == 0) ? pmb->block_size.nx2 : 1)
           *((ni.ox3 == 0) ? pmb->block_size.nx3 : 1);
    size *= nu_ + 1;
  } else {
    size = 0;
  }
  return size;
}


//----------------------------------------------------------------------------------------
//! \fn int SixRayBoundaryVariable::SetupPersistentMPI()
//! \brief
void SixRayBoundaryVariable::SetupPersistentMPI() {
  //TODO
  return;
}

//----------------------------------------------------------------------------------------
//! \fn int SixRayBoundaryVariable::StartReceiving(}
//! \brief
void SixRayBoundaryVariable::StartReceiving(BoundaryCommSubset phase) {
  //TODO
  return;
}

//----------------------------------------------------------------------------------------
//! \fn int SixRayBoundaryVariable::ClearBoundary(}
//! \brief
void SixRayBoundaryVariable::ClearBoundary(BoundaryCommSubset phase) {
  //TODO
  return;
}

//----------------------------------------------------------------------------------------
//! \fn int SixRayBoundaryVariable::LoadBoundaryBufferSameLevel(}
//! \brief
int SixRayBoundaryVariable::LoadBoundaryBufferSameLevel(Real *buf,
                                                        const NeighborBlock& nb) {
  //TODO
  return 0;
}

//----------------------------------------------------------------------------------------
//! \fn int SixRayBoundaryVariable::SetBoundarySameLevel(}
//! \brief
void SixRayBoundaryVariable::SetBoundarySameLevel(Real *buf, const NeighborBlock& nb) {
  //TODO
  return;
}

//----------------------------------------------------------------------------------------
//! \fn int SixRayBoundaryVariable::LoadBoundaryBufferToCoarser(}
//! \brief
int SixRayBoundaryVariable::LoadBoundaryBufferToCoarser(Real *buf, const NeighborBlock& nb) {
  //mesh refinement not implemented yet
  return 0;
}

//----------------------------------------------------------------------------------------
//! \fn int SixRayBoundaryVariable::LoadBoundaryBufferToFiner(}
//! \brief
int SixRayBoundaryVariable::LoadBoundaryBufferToFiner(Real *buf, const NeighborBlock& nb) {
  //mesh refinement not implemented yet
  return 0;
}

//----------------------------------------------------------------------------------------
//! \fn int SixRayBoundaryVariable::SetBoundaryFromCoarser(}
//! \brief
void SixRayBoundaryVariable::SetBoundaryFromCoarser(Real *buf, const NeighborBlock& nb) {
  //mesh refinement not implemented yet
  return;
}

//----------------------------------------------------------------------------------------
//! \fn int SixRayBoundaryVariable::SetBoundaryFromFiner(}
//! \brief
void SixRayBoundaryVariable::SetBoundaryFromFiner(Real *buf, const NeighborBlock& nb) {
  //mesh refinement not implemented yet
  return;
}

//no flux correction
int SixRayBoundaryVariable::ComputeFluxCorrectionBufferSize(const NeighborIndexes& ni, int cng) {
  return 0;
}
void SixRayBoundaryVariable::SendFluxCorrection() {}
bool SixRayBoundaryVariable::ReceiveFluxCorrection() {
  return true;
}

//shearing box not implemented
void SixRayBoundaryVariable::StartReceivingShear(BoundaryCommSubset phase) {}

//polar boundary not implemented
void SixRayBoundaryVariable::PolarBoundarySingleAzimuthalBlock() {}

//physical boundary not implemented. Use user defined boundary for six-ray
void SixRayBoundaryVariable::ReflectInnerX1(Real time, Real dt,
                    int il, int jl, int ju, int kl, int ku, int ngh) {}
void SixRayBoundaryVariable::ReflectOuterX1(Real time, Real dt,
                    int iu, int jl, int ju, int kl, int ku, int ngh) {}
void SixRayBoundaryVariable::ReflectInnerX2(Real time, Real dt,
                    int il, int iu, int jl, int kl, int ku, int ngh) {}
void SixRayBoundaryVariable::ReflectOuterX2(Real time, Real dt,
                    int il, int iu, int ju, int kl, int ku, int ngh) {}
void SixRayBoundaryVariable::ReflectInnerX3(Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ngh) {}
void SixRayBoundaryVariable::ReflectOuterX3(Real time, Real dt,
                    int il, int iu, int jl, int ju, int ku, int ngh) {}
void SixRayBoundaryVariable::OutflowInnerX1(Real time, Real dt,
                    int il, int jl, int ju, int kl, int ku, int ngh) {}
void SixRayBoundaryVariable::OutflowOuterX1(Real time, Real dt,
                    int iu, int jl, int ju, int kl, int ku, int ngh) {}
void SixRayBoundaryVariable::OutflowInnerX2(Real time, Real dt,
                    int il, int iu, int jl, int kl, int ku, int ngh) {}
void SixRayBoundaryVariable::OutflowOuterX2(Real time, Real dt,
                    int il, int iu, int ju, int kl, int ku, int ngh) {}
void SixRayBoundaryVariable::OutflowInnerX3(Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ngh) {}
void SixRayBoundaryVariable::OutflowOuterX3(Real time, Real dt,
                    int il, int iu, int jl, int ju, int ku, int ngh) {}
void SixRayBoundaryVariable::PolarWedgeInnerX2(Real time, Real dt,
                       int il, int iu, int jl, int kl, int ku, int ngh) {}
void SixRayBoundaryVariable::PolarWedgeOuterX2(Real time, Real dt,
                       int il, int iu, int ju, int kl, int ku, int ngh) {}


