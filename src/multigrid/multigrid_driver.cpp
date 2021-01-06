//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file multigrid_driver.cpp
//  \brief implementation of functions in class MultigridDriver

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdlib>    // abs
#include <iomanip>    // setprecision
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/cc/mg/bvals_mg.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "multigrid.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// constructor, initializes data structures and parameters

MultigridDriver::MultigridDriver(Mesh *pm, MGBoundaryFunc *MGBoundary, int invar) :
    nvar_(invar),
    mode_(0), // 0: V(1,1) FMG one sweep, 1: FMG + iterative, 2: V(1,1) iterative
    maxreflevel_(pm->multilevel?pm->max_level-pm->root_level:0),
    nrbx1_(pm->nrbx1), nrbx2_(pm->nrbx2), nrbx3_(pm->nrbx3), pmy_mesh_(pm),
    fsubtract_average_(false), ffas_(pm->multilevel), eps_(-1.0),
    cbuf_(nvar_,3,3,3), cbufold_(nvar_,3,3,3) {
}

// destructor

MultigridDriver::~MultigridDriver() {
}

//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SetupMultigrid()
//  \brief initialize the source assuming that the source terms are already loaded

void MultigridDriver::SetupMultigrid() {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SubtractAverage(MGVariable type)
//  \brief Calculate the global average and subtract it

void MultigridDriver::SubtractAverage(MGVariable type) {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::TransferFromBlocksToRoot(bool initflag)
//  \brief collect the coarsest data and transfer to the root grid

void MultigridDriver::TransferFromBlocksToRoot(bool initflag) {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::TransferFromRootToBlocks(bool folddata)
//  \brief Transfer the data from the root grid to the coarsest level of each MeshBlock

void MultigridDriver::TransferFromRootToBlocks(bool folddata) {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::FMGProlongate()
//  \brief Prolongation for FMG Cycle

void MultigridDriver::FMGProlongate() {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::OneStepToFiner(int nsmooth)
//  \brief prolongation and smoothing one level

void MultigridDriver::OneStepToFiner(int nsmooth) {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::OneStepToCoarser(int nsmooth)
//  \brief smoothing and restriction one level

void MultigridDriver::OneStepToCoarser(int nsmooth) {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SolveVCycle(int npresmooth, int npostsmooth)
//  \brief Solve the V-cycle starting from the current level

void MultigridDriver::SolveVCycle(int npresmooth, int npostsmooth) {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SolveFMGCycle()
//  \brief Solve the FMG Cycle using the V(1,1) or F(0,1) cycle

void MultigridDriver::SolveFMGCycle() {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SolveIterative(Real inidef)
//  \brief Solve iteratively until the convergence is achieved

void MultigridDriver::SolveIterative(Real inidef) {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SolveCoarsestGrid()
//  \brief Solve the coarsest root grid

void MultigridDriver::SolveCoarsestGrid() {
}


//----------------------------------------------------------------------------------------
//! \fn Real MultigridDriver::CalculateDefectNorm(MGNormType nrm, int n)
//  \brief calculate the defect norm

Real MultigridDriver::CalculateDefectNorm(MGNormType nrm, int n) {
  return 0.0;
}


//----------------------------------------------------------------------------------------
//! \fn Multigrid* MultigridDriver::FindMultigrid(int tgid)
//  \brief return the Multigrid whose gid is tgid

Multigrid* MultigridDriver::FindMultigrid(int tgid) {
  return nullptr;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::RestrictFMGSourceOctets()
//  \brief restrict the source in octets for FMG

void MultigridDriver::RestrictFMGSourceOctets() {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::RestrictOctets()
//  \brief restrict the potential in octets

void MultigridDriver::RestrictOctets() {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::ZeroClearOctets()
//  \brief zero clear the data in all the octets

void MultigridDriver::ZeroClearOctets() {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::StoreOldDataOctets()
//  \brief store the old u data in the uold array in octets

void MultigridDriver::StoreOldDataOctets() {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::CalculateFASRHSOctets()
//  \brief Calculate the RHS for FAS in Octets
void MultigridDriver::CalculateFASRHSOctets() {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SmoothOctets(int color)
//  \brief Apply the smoothing operator on octets
void MultigridDriver::SmoothOctets(int color) {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::ProlongateAndCorrectOctets()
//  \brief Prolongate and correct the potential in octets

void MultigridDriver::ProlongateAndCorrectOctets() {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::FMGProlongateOctets()
//  \brief Prolongate the potential in octets for FMG

void MultigridDriver::FMGProlongateOctets() {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SetBoundariesOctets(bool fprolong, bool folddata)
//  \brief Apply boundary conditions for octets

void MultigridDriver::SetBoundariesOctets(bool fprolong, bool folddata) {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SetOctetBoundarySameLevel(AthenaArray<Real> &dst,
//   const AthenaArray<Real> &un, AthenaArray<Real> &uold, const AthenaArray<Real> &unold,
//   int ox1, int ox2, int ox3, bool folddata)
//  \brief set an Octet boundary from a neighbor Octet on the same level

void MultigridDriver::SetOctetBoundarySameLevel(AthenaArray<Real> &dst,
     const AthenaArray<Real> &un, AthenaArray<Real> &uold, const AthenaArray<Real> &unold,
     int ox1, int ox2, int ox3, bool folddata) {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SetOctetBoundaryFromCoarser(const AthenaArray<Real> &un,
//                            const AthenaArray<Real> &unold, const LogicalLocation &loc,
//                            int ox1, int ox2, int ox3, bool folddata) {
//  \brief set a boundary in the coarse buffer from a neighbor Octet on the coarser level

void MultigridDriver::SetOctetBoundaryFromCoarser(const AthenaArray<Real> &un,
                      const AthenaArray<Real> &unold, const LogicalLocation &loc,
                      int ox1, int ox2, int ox3, bool folddata) {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::ApplyPhysicalBoundariesOctet(AthenaArray<Real> &u,
//                                         const LogicalLocation &loc, bool fcbuf)
//  \brief Apply physical boundary conditions for an octet

void MultigridDriver::ApplyPhysicalBoundariesOctet(AthenaArray<Real> &u,
                                             const LogicalLocation &loc, bool fcbuf) {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::RestrictOctetsBeforeTransfer()
//  \brief Restrict all the octets

void MultigridDriver::RestrictOctetsBeforeTransfer() {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SetOctetBoundariesBeforeTransfer(bool folddata)
//  \brief Set octet boundaries before transfer from root to blocks

void MultigridDriver::SetOctetBoundariesBeforeTransfer(bool folddata) {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::ProlongateOctetBoundaries(AthenaArray<Real> &u,
//                                           AthenaArray<Real> &uold, bool folddata)
//  \brief prolongate octet boundaries contacting the coarser level

void MultigridDriver::ProlongateOctetBoundaries(AthenaArray<Real> &u,
                                                AthenaArray<Real> &uold, bool folddata) {
}
