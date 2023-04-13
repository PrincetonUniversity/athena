//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file reflect.cpp
//  \brief implementation of reflecting BCs in each dimension

// C headers

// C++ headers

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "bvals_cc.hpp"

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::ReflectInnerX1(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, inner x1 boundary

void CellCenteredBoundaryVariable::VacuumInnerX1(
    Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh) {

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::ReflectOuterX1(
//          Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, outer x1 boundary

void CellCenteredBoundaryVariable::VacuumOuterX1(
    Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh) {

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::ReflectInnerX2(
//          Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, inner x2 boundary

void CellCenteredBoundaryVariable::VacuumInnerX2(
    Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh) {
                                                                 
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::ReflectOuterX2(
//          Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, outer x2 boundary

void CellCenteredBoundaryVariable::VacuumOuterX2(
    Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh) {

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::ReflectInnerX3(
//          Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ngh)
//  \brief REFLECTING boundary conditions, inner x3 boundary

void CellCenteredBoundaryVariable::VacuumInnerX3(
    Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ngh) {

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::ReflectOuterX3(
//          Real time, Real dt, int il, int iu, int jl, int ju, int ku, int ngh)
//  \brief REFLECTING boundary conditions, outer x3 boundary

void CellCenteredBoundaryVariable::VacuumOuterX3(
    Real time, Real dt, int il, int iu, int jl, int ju, int ku, int ngh) {

  return;
}
