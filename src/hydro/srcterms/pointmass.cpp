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
//! \file pointmass.cpp
//  \brief Adds source terms due to point mass AT ORIGIN
//======================================================================================

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../hydro.hpp"

// this class header
#include "hydro_srcterms.hpp"

//--------------------------------------------------------------------------------------
//! \fn void HydroSourceTerms::PointMass
//  \brief Adds source terms due to point mass AT ORIGIN

void HydroSourceTerms::PointMass(const Real dt, const AthenaArray<Real> *flux,
  const AthenaArray<Real> &prim, AthenaArray<Real> &cons)
{

  MeshBlock *pmb = pmy_hydro_->pmy_block;
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
#pragma omp parallel for schedule(static)
    for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real den = prim(IDN,k,j,i);
        Real src = dt*den*pmb->pcoord->coord_src1_i_(i)*gm_/pmb->pcoord->x1v(i);
        cons(IM1,k,j,i) -= src;
        if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) -=
          dt*0.5*(pmb->pcoord->phy_src1_i_(i)*flux[X1DIR](IDN,k,j,i)*gm_
                 +pmb->pcoord->phy_src2_i_(i)*flux[X1DIR](IDN,k,j,i+1)*gm_);
      }
    }
  }

  return;
}
