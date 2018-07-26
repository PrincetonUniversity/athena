//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file multigrid_driver.cpp
//  \brief implementation of functions in class MultigridDriver

// C/C++ headers
#include <algorithm>
#include <cmath>
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals_mg.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "multigrid.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// constructor, initializes data structures and parameters

MultigridDriver::MultigridDriver(Mesh *pm, MGBoundaryFunc_t *MGBoundary, int invar) {
}

// destructor

MultigridDriver::~MultigridDriver() {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::AddMultigrid(Multigrid *nmg)
//  \brief add a Multigrid object to the linked list
void MultigridDriver::AddMultigrid(Multigrid *nmg) {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SetupMultigrid(void)
//  \brief initialize the source assuming that the source terms are already loaded

void MultigridDriver::SetupMultigrid(void) {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SubtractAverage(int type)
//  \brief Calculate the global average and subtract it

void MultigridDriver::SubtractAverage(int type) {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::FillRootGridSource(void)
//  \brief collect the coarsest data and fill the root grid

void MultigridDriver::FillRootGridSource(void) {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::FMGProlongate(void)
//  \brief Prolongation for FMG Cycle

void MultigridDriver::FMGProlongate(void) {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::TransferFromRootToBlocks(void)
//  \brief Transfer the data from the root grid to the coarsest level of each MeshBlock

void MultigridDriver::TransferFromRootToBlocks(void) {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::OneStepToFiner(int nsmooth)
//  \brief smoothing and restriction one level

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
//! \fn void MultigridDriver::SolveFCycle(int npresmooth, int npostsmooth)
//  \brief Solve the F-cycle starting from the current level

void MultigridDriver::SolveFCycle(int npresmooth, int npostsmooth) {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SolveFMGCycle(void)
//  \brief Solve the FMG Cycle using the V(1,1) or F(0,1) cycle

void MultigridDriver::SolveFMGCycle(void) {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SolveIterative(void)
//  \brief Solve iteratively until the convergence is achieved

void MultigridDriver::SolveIterative(void) {
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SolveCoarsestGrid(void)
//  \brief Solve the coarsest root grid

void MultigridDriver::SolveCoarsestGrid(void) {
}


//----------------------------------------------------------------------------------------
//! \fn Real MultigridDriver::CalculateDefectNorm(int n, int nrm)
//  \brief calculate the defect norm

Real MultigridDriver::CalculateDefectNorm(int n, int nrm) {
}


//----------------------------------------------------------------------------------------
//! \fn Multigrid* MultigridDriver::FindMultigrid(int tgid)
//  \brief return the Multigrid whose gid is tgid

Multigrid* MultigridDriver::FindMultigrid(int tgid) {
}
