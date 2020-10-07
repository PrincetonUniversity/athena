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
//! \file shearing_box.cpp
//  \brief Adds source terms due to local shearing box approximation
//======================================================================================

// C headers

// C++ headers

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../hydro.hpp"
#include "hydro_srcterms.hpp"

//--------------------------------------------------------------------------------------
//! \fn void HydroSourceTerms::ShearingBoxSourceTerms(const Real dt,
//  const AthenaArray<Real> *flux, const AthenaArray<Real> &prim, AthenaArray<Real>
//  &cons)
//  \brief Shearing Box source terms
//
//  Detailed description starts here.
//  We add shearing box source term via operator splitting method. The source terms are
//  added after the fluxes are computed in each step of the integration (in
//  FluxDivergence) to give predictions of the conservative variables for either the
//  next step or the final update.

void HydroSourceTerms::ShearingBoxSourceTerms(const Real dt,
                                              const AthenaArray<Real> *flux,
                                              const AthenaArray<Real> &prim,
                                              AthenaArray<Real> &cons) {
  MeshBlock *pmb = pmy_hydro_->pmy_block;

  // 1) Tidal force:
  //    dM1/dt = 2q\rho\Omega^2 x
  //    dE /dt = 2q\Omega^2 (\rho v_x) x
  // 2) Coriolis forces:
  //    dM1/dt = 2\Omega(\rho v_y)
  //    dM2/dt = -2\Omega(\rho v_x)
  if (pmb->block_size.nx3 > 1 || ShBoxCoord_== 1) {
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real den    = prim(IDN,k,j,i);
          Real qO2    = qshear_*SQR(Omega_0_);
          Real xc     = pmb->pcoord->x1v(i);
          cons(IM1,k,j,i)   += 2.0*den*dt*(qO2*xc
                               +Omega_0_*prim(IVY,k,j,i));
          cons(IM2,k,j,i)   -= 2.0*dt*Omega_0_*den*prim(IVX,k,j,i);
          if (NON_BAROTROPIC_EOS) {
            Real rho_v1_d    = flux[X1DIR](IDN,k,j,i)
                               +flux[X1DIR](IDN,k,j,i+1);
            cons(IEN,k,j,i) += dt*qO2*xc*rho_v1_d;
          }
        }
      }
    }
  } else if (pmb->block_size.nx3 == 1 && ShBoxCoord_ == 2) {
    int ks = pmb->ks;
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real den      = prim(IDN,ks,j,i);
        Real qO2      = qshear_*SQR(Omega_0_);
        Real xc       = pmb->pcoord->x1v(i);
        cons(IM1,ks,j,i)   += 2.0*den*dt*(qO2*xc
                              +Omega_0_*prim(IVZ,ks,j,i));
        cons(IM3,ks,j,i)   -= 2.0*dt*Omega_0_*den*prim(IVX,ks,j,i);
        if (NON_BAROTROPIC_EOS) {
          Real rho_v1_d     = flux[X1DIR](IDN,ks,j,i)
                              +flux[X1DIR](IDN,ks,j,i+1);
          cons(IEN,ks,j,i) += dt*qO2*xc*rho_v1_d;
        }
      }
    }
  } else {
    std::cout << "[ShearingBoxSourceTerms]: not compatible to 1D !!" << std::endl;
    return;
  }


  return;
}
