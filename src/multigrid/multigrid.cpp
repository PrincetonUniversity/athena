//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file multigrid.cpp
//  \brief implementation of the functions commonly used in Multigrid

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstring>    // memset, memcpy
#include <iostream>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "multigrid.hpp"

//----------------------------------------------------------------------------------------
//! \fn Multigrid::Multigrid(MultigridDriver *pmd, MeshBlock *pmb, int invar, int nghost)
//  \brief Multigrid constructor

Multigrid::Multigrid(MultigridDriver *pmd, MeshBlock *pmb, int invar, int nghost) :
  pmy_driver_(pmd), pmy_block_(pmb), ngh_(nghost), nvar_(invar), defscale_(1.0) {
}


//----------------------------------------------------------------------------------------
//! \fn Multigrid::~Multigrid
//  \brief Multigrid destroctor

Multigrid::~Multigrid() {
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::LoadFinestData(const AthenaArray<Real> &src, int ns, int ngh)
//  \brief Fill the inital guess in the active zone of the finest level

void Multigrid::LoadFinestData(const AthenaArray<Real> &src, int ns, int ngh) {
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::LoadSource(const AthenaArray<Real> &src, int ns, int ngh,
//                                 Real fac)
//  \brief Fill the source in the active zone of the finest level

void Multigrid::LoadSource(const AthenaArray<Real> &src, int ns, int ngh, Real fac) {
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::RestrictFMGSource()
//  \brief restrict the source through all the multigrid levels

void Multigrid::RestrictFMGSource() {
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::RetrieveResult(AthenaArray<Real> &dst, int ns, int ngh)
//  \brief Set the result, including the ghost zone

void Multigrid::RetrieveResult(AthenaArray<Real> &dst, int ns, int ngh) {
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::RetrieveDefect(AthenaArray<Real> &dst, int ns, int ngh)
//  \brief Set the defect, including the ghost zone

void Multigrid::RetrieveDefect(AthenaArray<Real> &dst, int ns, int ngh) {
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::ZeroClearData()
//  \brief Clear the data array with zero

void Multigrid::ZeroClearData() {
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::RestrictBlock()
//  \brief Restrict the defect to the source

void Multigrid::RestrictBlock() {
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::ProlongateAndCorrectBlock()
//  \brief Prolongate the potential using tri-linear interpolation

void Multigrid::ProlongateAndCorrectBlock() {
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::FMGProlongateBlock()
//  \brief Prolongate the potential for Full Multigrid cycle

void Multigrid::FMGProlongateBlock() {
}


//----------------------------------------------------------------------------------------
//! \fn  void Multigrid::SmoothBlock(int color)
//  \brief Apply Smoother on the Block

void Multigrid::SmoothBlock(int color) {
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::CalculateDefectBlock()
//  \brief calculate the residual

void Multigrid::CalculateDefectBlock() {
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::CalculateFASRHSBlock()
//  \brief calculate the RHS for the Full Approximation Scheme

void Multigrid::CalculateFASRHSBlock() {
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::SetFromRootGrid(bool folddata)
//  \brief Load the data from the root grid or octets

void Multigrid::SetFromRootGrid(bool folddata) {
}


//----------------------------------------------------------------------------------------
//! \fn Real Multigrid::CalculateDefectNorm(MGNormType nrm, int n)
//  \brief calculate the residual norm

Real Multigrid::CalculateDefectNorm(MGNormType nrm, int n) {
  return 0.0;
}


//----------------------------------------------------------------------------------------
//! \fn Real Multigrid::CalculateTotal(MGVariable type, int n)
//  \brief calculate the sum of the array (type: 0=src, 1=u)

Real Multigrid::CalculateTotal(MGVariable type, int n) {
  return 0.0;
}


//----------------------------------------------------------------------------------------
//! \fn Real Multigrid::SubtractAverage(MGVariable type, int v, Real ave)
//  \brief subtract the average value (type: 0=src, 1=u)

void Multigrid::SubtractAverage(MGVariable type, int n, Real ave) {
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::StoreOldData()
//  \brief store the old u data in the uold array

void Multigrid::StoreOldData() {
}


//----------------------------------------------------------------------------------------
//! \fn Real Multigrid::GetCoarsestData(MGVariable type, int n)
//  \brief get the value on the coarsest level in the MG block (type: 0=src, 1=u)

Real Multigrid::GetCoarsestData(MGVariable type, int n) {
  return 0.0;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::SetData(MGVariable type, int n, int k, int j, int i, Real v)
//  \brief set a value to a cell on the current level

void Multigrid::SetData(MGVariable type, int n, int k, int j, int i, Real v) {
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::Restrict(AthenaArray<Real> &dst, const AthenaArray<Real> &src,
//                               int il, int iu, int jl, int ju, int kl, int ku)
//  \brief Actual implementation of prolongation and correction

void Multigrid::Restrict(AthenaArray<Real> &dst, const AthenaArray<Real> &src,
                         int il, int iu, int jl, int ju, int kl, int ku) {
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::ProlongateAndCorrect(AthenaArray<Real> &dst,
//      const AthenaArray<Real> &src, int il, int iu, int jl, int ju, int kl, int ku,
//      int fil, int fjl, int fkl)
//  \brief Actual implementation of prolongation and correction

void Multigrid::ProlongateAndCorrect(AthenaArray<Real> &dst, const AthenaArray<Real> &src,
     int il, int iu, int jl, int ju, int kl, int ku, int fil, int fjl, int fkl) {
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::FMGProlongate(AthenaArray<Real> &dst,
//           const AthenaArray<Real> &src, int il, int iu, int jl, int ju, int kl, int ku
//           int fil, int fjl, int fkl)
//  \brief Actual implementation of FMG prolongation

void Multigrid::FMGProlongate(AthenaArray<Real> &dst, const AthenaArray<Real> &src,
                              int il, int iu, int jl, int ju, int kl, int ku,
                              int fil, int fjl, int fkl) {
}

