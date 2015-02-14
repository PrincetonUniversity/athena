//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
//
//======================================================================================

// Primary header
#include "srcterms.hpp"

// Athena headers
#include "../../athena.hpp"          // Real
#include "../../athena_arrays.hpp"   // AthenaArray
#include "../../mesh.hpp"            // MeshBlock
#include "../../coordinates/coordinates.hpp"  // src_terms_i_
#include "../fluid.hpp"              // Fluid
#include "../../parameter_input.hpp" // ParameterInput

//======================================================================================
//! \file srcterms.cpp
//  \brief implements functions that compute physical source terms in the fluid
//======================================================================================

// FluidSourceTerms constructor - sets function pointers for each of the physical source
// terms to be included in the calculation.

FluidSourceTerms::FluidSourceTerms(Fluid *pf, ParameterInput *pin)
{
  pmy_fluid_ = pf;
  gm_ = pin->GetOrAddReal("problem","GM",0.0);
  UserSourceTerm = NULL;
}

// destructor

FluidSourceTerms::~FluidSourceTerms()
{
}

//--------------------------------------------------------------------------------------
//! \fn
//  \brief

void FluidSourceTerms::PhysicalSourceTerms(const Real time, const Real dt,
  const AthenaArray<Real> &prim, AthenaArray<Real> &cons)
{
  if (gm_ == 0.0) return;

// Source terms due to point mass gravity

  MeshBlock *pmb = pmy_fluid_->pmy_block;
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
  for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma simd
    for (int i=pmb->is; i<=pmb->ie; ++i) {
      Real src = dt*(pmb->pcoord->coord_src1_i_(i))*(gm_*prim(IDN,k,j,i)/pmb->x1v(i));
      cons(IM1,k,j,i) -= src;
      if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) -= src*prim(IVX,k,j,i);
    }
  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn
//  \brief

void FluidSourceTerms::EnrollSrcTermFunction(SrcTermFunc_t my_func)
{
  UserSourceTerm = my_func;
  return;
}
