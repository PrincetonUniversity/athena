//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file pointmass.cpp
//! \brief Adds source terms due to point mass AT ORIGIN

// C headers

// C++ headers

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../hydro.hpp"
#include "hydro_srcterms.hpp"

//----------------------------------------------------------------------------------------
//! \fn void HydroSourceTerms::PointMass
//! \brief Adds source terms due to point mass AT ORIGIN

void HydroSourceTerms::PointMass(const Real dt, const AthenaArray<Real> *flux,
                                 const AthenaArray<Real> &prim, AthenaArray<Real> &cons) {
  MeshBlock *pmb = pmy_hydro_->pmy_block;
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real den = prim(IDN,k,j,i);
        Real src = dt*den*pmb->pcoord->coord_src1_i_(i)*gm_/pmb->pcoord->x1v(i);
        cons(IM1,k,j,i) -= src;
        if (NON_BAROTROPIC_EOS) {
          cons(IEN,k,j,i) -=
              dt*0.5*(pmb->pcoord->phy_src1_i_(i)*flux[X1DIR](IDN,k,j,i)*gm_
                      +pmb->pcoord->phy_src2_i_(i)*flux[X1DIR](IDN,k,j,i+1)*gm_);
        }
      }
    }
  }

  return;
}
