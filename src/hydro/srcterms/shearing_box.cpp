//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file shearing_box.cpp
//! \brief Adds source terms due to local shearing box approximation
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
//!   const AthenaArray<Real> *flux, const AthenaArray<Real> &prim,
//!   AthenaArray<Real> &cons)
//! \brief Shearing Box source terms
//!
//! We add shearing box source term via operator splitting method. The source terms are
//! added after the fluxes are computed in each step of the integration (in
//! FluxDivergence) to give predictions of the conservative variables for either the
//! next step or the final update.

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
  if (ShBoxCoord_== 1) {
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          const Real &den = prim(IDN,k,j,i);
          const Real qO2  = qshear_*SQR(Omega_0_);
          const Real mom1 = den*prim(IVX,k,j,i);
          const Real &xc  = pmb->pcoord->x1v(i);
          cons(IM1,k,j,i) += 2.0*dt*(Omega_0_*(den*prim(IVY,k,j,i))+qO2*den*xc);
          cons(IM2,k,j,i) -= 2.0*dt*Omega_0_*mom1;
          if (NON_BAROTROPIC_EOS) {
            const Real phic = qO2*SQR(xc);
            const Real phil = qO2*SQR(pmb->pcoord->x1f(i));
            const Real phir = qO2*SQR(pmb->pcoord->x1f(i+1));
            cons(IEN,k,j,i) += dt*(flux[X1DIR](IDN,k,j,i)*(phic-phil)+
                                   flux[X1DIR](IDN,k,j,i+1)*(phir-phic))
                                   /pmb->pcoord->dx1f(i);
          }
        }
      }
    }
  } else { // ShBoxCoord_== 2
    int ks = pmb->ks;
    for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        const Real &den = prim(IDN,ks,j,i);
        const Real qO2  = qshear_*SQR(Omega_0_);
        const Real mom1 = den*prim(IVX,ks,j,i);
        const Real &xc  = pmb->pcoord->x1v(i);
        cons(IM1,ks,j,i) += 2.0*dt*(Omega_0_*(den*prim(IVZ,ks,j,i))+qO2*den*xc);
        cons(IM3,ks,j,i) -= 2.0*dt*Omega_0_*mom1;
        if (NON_BAROTROPIC_EOS) {
          const Real phic = qO2*SQR(xc);
          const Real phil = qO2*SQR(pmb->pcoord->x1f(i));
          const Real phir = qO2*SQR(pmb->pcoord->x1f(i+1));
          cons(IEN,ks,j,i) += dt*(flux[X1DIR](IDN,ks,j,i)*(phic-phil)+
                                  flux[X1DIR](IDN,ks,j,i+1)*(phir-phic))
                                  /pmb->pcoord->dx1f(i);
        }
      }
    }
  }
  return;
}
