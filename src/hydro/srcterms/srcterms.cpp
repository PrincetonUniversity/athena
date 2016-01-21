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
//======================================================================================
//! \file srcterms.cpp
//  \brief implements functions that compute physical source terms for hydro
//======================================================================================

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../hydro.hpp"
#include "../../parameter_input.hpp"

// this class header
#include "srcterms.hpp"

// HydroSourceTerms constructor - sets function pointers for each of the physical source
// terms to be included in the calculation.

HydroSourceTerms::HydroSourceTerms(Hydro *phyd, ParameterInput *pin)
{
  pmy_hydro_ = phyd;
  gm_ = pin->GetOrAddReal("problem","GM",0.0);
  g1_ = pin->GetOrAddReal("hydro","grav_acc1",0.0);
  g2_ = pin->GetOrAddReal("hydro","grav_acc2",0.0);
  g3_ = pin->GetOrAddReal("hydro","grav_acc3",0.0);
  UserSourceTerm = phyd->pmy_block->pmy_mesh->UserSourceTerm_;
}

// destructor

HydroSourceTerms::~HydroSourceTerms()
{
}

//--------------------------------------------------------------------------------------
//! \fn void HydroSourceTerms::PhysicalSourceTerms(const Real dt,
//  const AthenaArray<Real> *flux, const AthenaArray<Real> &prim, AthenaArray<Real> &cons)
//  \brief Physical (gravitational) source terms

void HydroSourceTerms::PhysicalSourceTerms(const Real dt, const AthenaArray<Real> *flux,
  const AthenaArray<Real> &prim, AthenaArray<Real> &cons)
{
  if (gm_==0.0 && g1_==0.0 && g2_==0.0 && g3_==0.0) return;

// Source terms due to point mass gravity

  MeshBlock *pmb = pmy_hydro_->pmy_block;
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
#pragma omp parallel for schedule(static)
    for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real den = prim(IDN,k,j,i);

        if (gm_!=0.0) {
          Real src = dt*den*pmb->pcoord->coord_src1_i_(i)*gm_/pmb->pcoord->x1v(i);
          cons(IM1,k,j,i) -= src;
          if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) -=
            dt*0.5*(pmb->pcoord->phy_src1_i_(i)*flux[x1face](IDN,k,j,i)*gm_
                   +pmb->pcoord->phy_src2_i_(i)*flux[x1face](IDN,k,j,i+1)*gm_);
        }

        if (g1_!=0.0) {
          Real src = dt*den*g1_;
          cons(IM1,k,j,i) += src;
          if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) += src*prim(IVX,k,j,i);
        }

        if (g2_!=0.0) {
          Real src = dt*den*g2_;
          cons(IM2,k,j,i) += src;
          if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) += src*prim(IVY,k,j,i);
        }

        if (g3_!=0.0) {
          Real src = dt*den*g3_;
          cons(IM3,k,j,i) += src;
          if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) += src*prim(IVZ,k,j,i);
        }
      }
    }
  }

  return;
}

