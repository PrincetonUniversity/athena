//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file self_gravity.cpp
//! \brief source terms due to self-gravity

// C headers

// C++ headers

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../gravity/gravity.hpp"
#include "../../mesh/mesh.hpp"
#include "../hydro.hpp"
#include "hydro_srcterms.hpp"

//----------------------------------------------------------------------------------------
//! \fn void HydroSourceTerms::SelfGravity
//! \brief Adds source terms for self-gravitational acceleration to conserved variables
//! \note
//! This implements the source term formula in Mullen et al. 2020, but only for
//! the momentum part. The energy source term is not conservative in this version.
//! I leave the fully conservative formula for later as it requires design consideration.

void HydroSourceTerms::SelfGravity(const Real dt,const AthenaArray<Real> *flux,
                                   const AthenaArray<Real> &prim,
                                   AthenaArray<Real> &cons) {
  MeshBlock *pmb = pmy_hydro_->pmy_block;
  Gravity *pgrav = pmb->pgrav;

  if (NON_BAROTROPIC_EOS) {
    if (GRAVITY_FLUX_ENABLED) { // energy only
      // acceleration in 1-direction
      for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            Real dx1 = pmb->pcoord->dx1v(i);
            Real dtodx1 = dt/dx1;
            Real phic = pgrav->phi(k,j,i);
            Real phil = 0.5*(pgrav->phi(k,j,i-1)+pgrav->phi(k,j,i  ));
            Real phir = 0.5*(pgrav->phi(k,j,i  )+pgrav->phi(k,j,i+1));
            cons(IEN,k,j,i) -= dtodx1*(flux[X1DIR](IDN,k,j,i  )*(phic - phil) +
                                       flux[X1DIR](IDN,k,j,i+1)*(phir - phic));
          }
        }
      }

      if (pmb->block_size.nx2 > 1) {
        // acceleration in 2-direction
        for (int k=pmb->ks; k<=pmb->ke; ++k) {
          for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
            for (int i=pmb->is; i<=pmb->ie; ++i) {
              Real dx2 = pmb->pcoord->dx2v(j);
              Real dtodx2 = dt/dx2;
              Real phic = pgrav->phi(k,j,i);
              Real phil = 0.5*(pgrav->phi(k,j-1,i)+pgrav->phi(k,j  ,i));
              Real phir = 0.5*(pgrav->phi(k,j  ,i)+pgrav->phi(k,j+1,i));
              cons(IEN,k,j,i) -= dtodx2*(flux[X2DIR](IDN,k,j  ,i)*(phic - phil) +
                                         flux[X2DIR](IDN,k,j+1,i)*(phir - phic));
            }
          }
        }
      }

      if (pmb->block_size.nx3 > 1) {
        // acceleration in 3-direction
        for (int k=pmb->ks; k<=pmb->ke; ++k) {
          for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
            for (int i=pmb->is; i<=pmb->ie; ++i) {
              Real dx3 = pmb->pcoord->dx3v(k);
              Real dtodx3 = dt/dx3;
              Real phic = pgrav->phi(k,j,i);
              Real phil = 0.5*(pgrav->phi(k-1,j,i)+pgrav->phi(k  ,j,i));
              Real phir = 0.5*(pgrav->phi(k  ,j,i)+pgrav->phi(k+1,j,i));
              cons(IEN,k,j,i) -= dtodx3*(flux[X3DIR](IDN,k  ,j,i)*(phic - phil) +
                                         flux[X3DIR](IDN,k+1,j,i)*(phir - phic));
            }
          }
        }
      }
    } else { // energy and momentum
      // acceleration in 1-direction
      for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            Real dx1 = pmb->pcoord->dx1v(i);
            Real dtodx1 = dt/dx1;
            Real phic = pgrav->phi(k,j,i);
            Real phil = 0.5*(pgrav->phi(k,j,i-1)+pgrav->phi(k,j,i  ));
            Real phir = 0.5*(pgrav->phi(k,j,i  )+pgrav->phi(k,j,i+1));
            cons(IM1,k,j,i) -= dtodx1*prim(IDN,k,j,i)*(phir-phil);
            cons(IEN,k,j,i) -= dtodx1*(flux[X1DIR](IDN,k,j,i  )*(phic - phil) +
                                       flux[X1DIR](IDN,k,j,i+1)*(phir - phic));
          }
        }
      }

      if (pmb->block_size.nx2 > 1) {
        // acceleration in 2-direction
        for (int k=pmb->ks; k<=pmb->ke; ++k) {
          for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
            for (int i=pmb->is; i<=pmb->ie; ++i) {
              Real dx2 = pmb->pcoord->dx2v(j);
              Real dtodx2 = dt/dx2;
              Real phic = pgrav->phi(k,j,i);
              Real phil = 0.5*(pgrav->phi(k,j-1,i)+pgrav->phi(k,j  ,i));
              Real phir = 0.5*(pgrav->phi(k,j  ,i)+pgrav->phi(k,j+1,i));
              cons(IM2,k,j,i) -= dtodx2*prim(IDN,k,j,i)*(phir-phil);
              cons(IEN,k,j,i) -= dtodx2*(flux[X2DIR](IDN,k,j  ,i)*(phic - phil) +
                                         flux[X2DIR](IDN,k,j+1,i)*(phir - phic));
            }
          }
        }
      }

      if (pmb->block_size.nx3 > 1) {
        // acceleration in 3-direction
        for (int k=pmb->ks; k<=pmb->ke; ++k) {
          for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
            for (int i=pmb->is; i<=pmb->ie; ++i) {
              Real dx3 = pmb->pcoord->dx3v(k);
              Real dtodx3 = dt/dx3;
              Real phic = pgrav->phi(k,j,i);
              Real phil = 0.5*(pgrav->phi(k-1,j,i)+pgrav->phi(k  ,j,i));
              Real phir = 0.5*(pgrav->phi(k  ,j,i)+pgrav->phi(k+1,j,i));
              cons(IM3,k,j,i) -= dtodx3*prim(IDN,k,j,i)*(phir-phil);
              cons(IEN,k,j,i) -= dtodx3*(flux[X3DIR](IDN,k  ,j,i)*(phic - phil) +
                                         flux[X3DIR](IDN,k+1,j,i)*(phir - phic));
            }
          }
        }
      }
    }
  } else if (!GRAVITY_FLUX_ENABLED) { // momentum only
    // acceleration in 1-direction
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real dx1 = pmb->pcoord->dx1v(i);
          Real dtodx1 = dt/dx1;
          Real phil = 0.5*(pgrav->phi(k,j,i-1)+pgrav->phi(k,j,i  ));
          Real phir = 0.5*(pgrav->phi(k,j,i  )+pgrav->phi(k,j,i+1));
          cons(IM1,k,j,i) -= dtodx1*prim(IDN,k,j,i)*(phir-phil);
        }
      }
    }

    if (pmb->block_size.nx2 > 1) {
      // acceleration in 2-direction
      for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            Real dx2 = pmb->pcoord->dx2v(j);
            Real dtodx2 = dt/dx2;
            Real phil = 0.5*(pgrav->phi(k,j-1,i)+pgrav->phi(k,j  ,i));
            Real phir = 0.5*(pgrav->phi(k,j  ,i)+pgrav->phi(k,j+1,i));
            cons(IM2,k,j,i) -= dtodx2*prim(IDN,k,j,i)*(phir-phil);
          }
        }
      }
    }

    if (pmb->block_size.nx3 > 1) {
      // acceleration in 3-direction
      for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            Real dx3 = pmb->pcoord->dx3v(k);
            Real dtodx3 = dt/dx3;
            Real phil = 0.5*(pgrav->phi(k-1,j,i)+pgrav->phi(k  ,j,i));
            Real phir = 0.5*(pgrav->phi(k  ,j,i)+pgrav->phi(k+1,j,i));
            cons(IM3,k,j,i) -= dtodx3*prim(IDN,k,j,i)*(phir-phil);
          }
        }
      }
    }
  }

  return;
}
