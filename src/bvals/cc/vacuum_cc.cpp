//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file
//  \brief

// C headers

// C++ headers

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "./bvals_cc.hpp"

//----------------------------------------------------------------------------------------

void CellCenteredBoundaryVariable::VacuumInnerX1(
    Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh) {
  return;
}

//----------------------------------------------------------------------------------------

void CellCenteredBoundaryVariable::VacuumOuterX1(
    Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh) {
  return;
}

//----------------------------------------------------------------------------------------

void CellCenteredBoundaryVariable::VacuumInnerX2(
    Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh) {
  return;
}

//----------------------------------------------------------------------------------------

void CellCenteredBoundaryVariable::VacuumOuterX2(
    Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh) {
  return;
}

//----------------------------------------------------------------------------------------

void CellCenteredBoundaryVariable::VacuumInnerX3(
    Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ngh) {
  return;
}

//----------------------------------------------------------------------------------------

void CellCenteredBoundaryVariable::VacuumOuterX3(
    Real time, Real dt, int il, int iu, int jl, int ju, int ku, int ngh) {
  return;
}
