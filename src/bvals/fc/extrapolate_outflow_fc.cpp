//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file outflow.cpp
//  \brief implementation of extrapolated outflow BCs in each dimension

// C, C++ headers
#include <iostream>

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "bvals_fc.hpp"


namespace {
  void NotImplementedError_(){
    std::cerr << "FC extrapolation not implemented." << std::endl;
    std::abort();
  };
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::ExtrapolateOutflowInnerX1(
//               Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief OUTFLOW boundary conditions with extrapolation, inner x1 boundary

void FaceCenteredBoundaryVariable::ExtrapolateOutflowInnerX1(
    Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh) {
  NotImplementedError_();
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::ExtrapolateOutflowOuterX1(
//               Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh)
//  \brief OUTFLOW boundary conditions with extrapolation, outer x1 boundary

void FaceCenteredBoundaryVariable::ExtrapolateOutflowOuterX1(
    Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh) {
  NotImplementedError_();
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::ExtrapolateOutflowInnerX2(
//               Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh)
//  \brief OUTFLOW boundary conditions with extrapolation, inner x2 boundary

void FaceCenteredBoundaryVariable::ExtrapolateOutflowInnerX2(
    Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh) {
  NotImplementedError_();
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::ExtrapolateOutflowOuterX2(
//               Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh)
//  \brief OUTFLOW boundary conditions with extrapolation, outer x2 boundary

void FaceCenteredBoundaryVariable::ExtrapolateOutflowOuterX2(
    Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh) {
  NotImplementedError_();
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::ExtrapolateOutflowInnerX3(
//               Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ngh)
//  \brief OUTFLOW boundary conditions with extrapolation, inner x3 boundary

void FaceCenteredBoundaryVariable::ExtrapolateOutflowInnerX3(
    Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ngh) {
  NotImplementedError_();
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::ExtrapolateOutflowOuterX3(
//               Real time, Real dt, int il, int iu, int jl, int ju, int ku, int ngh)
//  \brief OUTFLOW boundary conditions with extrapolation, outer x3 boundary

void FaceCenteredBoundaryVariable::ExtrapolateOutflowOuterX3(
    Real time, Real dt, int il, int iu, int jl, int ju, int ku, int ngh) {
  NotImplementedError_();
}
