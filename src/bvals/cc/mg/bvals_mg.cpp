//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_mg.cpp
//! \brief

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>    // memcpy()
#include <iomanip>
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../coordinates/coordinates.hpp"
#include "../../../globals.hpp"
#include "../../../mesh/mesh.hpp"
#include "../../../multigrid/multigrid.hpp"
#include "../../../parameter_input.hpp"
#include "../../../utils/buffer_utils.hpp"
#include "bvals_mg.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

class Multigrid;
class MultigridDriver;

//----------------------------------------------------------------------------------------
//! \fn MGBoundaryValues::MGBoundaryValues(Multigrid *pmg, BoundaryFlag *input_bcs)
//! \brief Constructor of the MGBoundaryValues class

MGBoundaryValues::MGBoundaryValues(Multigrid *pmg, BoundaryFlag *input_bcs)
    : BoundaryBase(pmg->pmy_driver_->pmy_mesh_, pmg->loc_, pmg->size_, input_bcs),
      pmy_mg_(pmg) {
}


//----------------------------------------------------------------------------------------
//! \fn MGBoundaryValues::~MGBoundaryValues()
//! \brief Destructor of the MGBoundaryValues class

MGBoundaryValues::~MGBoundaryValues() {
}

//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::InitBoundaryData(BoundaryQuantity type)
//! \brief Initialize BoundaryData<> structure

void MGBoundaryValues::InitBoundaryData(BoundaryQuantity type) {
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::DestroyBoundaryData()
//! \brief Destroy BoundaryData<> structure

void MGBoundaryValues::DestroyBoundaryData() {
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::ApplyPhysicalBoundaries(int flag)
//! \brief Apply physical boundary conditions to the current Multigrid data

void MGBoundaryValues::ApplyPhysicalBoundaries(int flag) {
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::StartReceivingMultigrid(BoundaryQuantity type,
//!                                                    bool folddata)
//! \brief initiate MPI_Irecv for multigrid

void MGBoundaryValues::StartReceivingMultigrid(BoundaryQuantity type, bool folddata) {
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::ClearBoundaryMultigrid(BoundaryQuantity type)
//! \brief clean up the boundary flags after each loop for multigrid

void MGBoundaryValues::ClearBoundaryMultigrid(BoundaryQuantity type) {
}


//----------------------------------------------------------------------------------------
//! \fn int MGBoundaryValues::LoadMultigridBoundaryBufferSameLevel(Real *buf,
//!                                                const NeighborBlock& nb, bool folddata)
//! \brief Set Multigrid boundary buffers for sending to a block on the same level

int MGBoundaryValues::LoadMultigridBoundaryBufferSameLevel(Real *buf,
                                           const NeighborBlock& nb, bool folddata) {
  return 0;
}


//----------------------------------------------------------------------------------------
//! \fn int MGBoundaryValues::LoadMultigridBoundaryBufferToCoarser(Real *buf,
//!                                                const NeighborBlock& nb, bool folddata)
//! \brief Set Multigrid boundary buffers for sending to a block on the coarser level

int MGBoundaryValues::LoadMultigridBoundaryBufferToCoarser(Real *buf,
                                           const NeighborBlock& nb, bool folddata) {
  return 0;
}


//----------------------------------------------------------------------------------------
//! \fn int MGBoundaryValues::LoadMultigridBoundaryBufferToFiner(Real *buf,
//!                                        const NeighborBlock& nb, bool folddata)
//! \brief Set Multigrid boundary buffers for sending to a block on the finer level
int MGBoundaryValues::LoadMultigridBoundaryBufferToFiner(Real *buf,
                                           const NeighborBlock& nb, bool folddata) {
  return 0;
}


//----------------------------------------------------------------------------------------
//! \fn bool MGBoundaryValues::SendMultigridBoundaryBuffers(BoundaryQuantity type,
//!                                                         bool folddata)
//! \brief Send boundary buffers

bool MGBoundaryValues::SendMultigridBoundaryBuffers(BoundaryQuantity type,
                                                    bool folddata) {
  return true;
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::SetMultigridBoundarySameLevel(const Real *buf,
//!                                                const NeighborBlock& nb, bool folddata)
//! \brief Set Multigrid boundary received from a block on the same level

void MGBoundaryValues::SetMultigridBoundarySameLevel(const Real *buf,
                                           const NeighborBlock& nb, bool folddata) {
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::SetMultigridBoundaryFromCoarser(const Real *buf,
//!                                        const NeighborBlock& nb, boolf folddata)
//! \brief Set hydro boundary received from a block on the same level

void MGBoundaryValues::SetMultigridBoundaryFromCoarser(const Real *buf,
                                           const NeighborBlock& nb, bool folddata) {
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::SetMultigridBoundaryFromFiner(const Real *buf,
//!                                                const NeighborBlock& nb, bool folddata)
//! \brief Set hydro boundary received from a block on the same level

void MGBoundaryValues::SetMultigridBoundaryFromFiner(const Real *buf,
                                               const NeighborBlock& nb, bool folddata) {
}


//----------------------------------------------------------------------------------------
//! \fn bool MGBoundaryValues::ReceiveMultigridBoundaryBuffers(BoundaryQuantity type,
//!                                                            bool folddata)
//! \brief receive the boundary data

bool MGBoundaryValues::ReceiveMultigridBoundaryBuffers(BoundaryQuantity type,
                                                       bool folddata) {
  return true;
}


//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::ProlongateMultigridBoundaries(bool folddata)
//! \brief prolongate boundaries for Multigrid

void MGBoundaryValues::ProlongateMultigridBoundaries(bool folddata) {
}

//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::CopyNeighborInfoFromMeshBlock()
//! \brief copy the neighbor information from the MeshBlock BoundaryValues class

void MGBoundaryValues::CopyNeighborInfoFromMeshBlock() {
}


//----------------------------------------------------------------------------------------
//! \fn int MGGravityBoundaryValues::LoadMultigridBoundaryBufferToCoarserFluxCons(
//!                                               Real *buf, const NeighborBlock& nb)
//! \brief Set Multigrid boundary buffers for sending to a block on the coarser level
//!        using the Mass Conservation formula for gravity

int MGGravityBoundaryValues::LoadMultigridBoundaryBufferToCoarserFluxCons(Real *buf,
                                                          const NeighborBlock& nb) {
  return 0;
}


//----------------------------------------------------------------------------------------
//! \fn int MGGravityBoundaryValues::LoadMultigridBoundaryBufferToFinerFluxCons(Real *buf,
//!                                                             const NeighborBlock& nb)
//! \brief Set Multigrid boundary buffers for sending to a block on the finer level
int MGGravityBoundaryValues::LoadMultigridBoundaryBufferToFinerFluxCons(Real *buf,
                                                               const NeighborBlock& nb) {
  return 0;
}


//----------------------------------------------------------------------------------------
//! \fn void MGGravityBoundaryValues::SetMultigridBoundaryFromCoarserFluxCons(
//!                                      const Real *buf, const NeighborBlock& nb)
//! \brief Set hydro boundary received from a block on the same level

void MGGravityBoundaryValues::SetMultigridBoundaryFromCoarserFluxCons(const Real *buf,
                                                             const NeighborBlock& nb) {
}


//----------------------------------------------------------------------------------------
//! \fn void MGGravityBoundaryValues::SetMultigridBoundaryFromFinerFluxCons(
//!                                      const Real *buf, const NeighborBlock& nb)

void MGGravityBoundaryValues::SetMultigridBoundaryFromFinerFluxCons(const Real *buf,
                                                             const NeighborBlock& nb) {
  return;
}



//----------------------------------------------------------------------------------------
//! \fn void MGBoundaryValues::ProlongateMultigridBoundariesFluxCons()
//! \brief prolongate boundaries for Multigrid

void MGBoundaryValues::ProlongateMultigridBoundariesFluxCons() {
}
