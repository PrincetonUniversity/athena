//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//  \brief source terms due to self-gravity

// Athena++ headers
#include "hydro_srcterms.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../hydro.hpp"
#include "../../gravity/gravity.hpp"

//----------------------------------------------------------------------------------------
//! \fn void HydroSourceTerms::SelfGravity
//  \brief Adds source terms for self-gravitational acceleration to conserved variables

void HydroSourceTerms::SelfGravity(const Real dt,const AthenaArray<Real> *flux,
  const AthenaArray<Real> &prim, AthenaArray<Real> &cons) {
  MeshBlock *pmb = pmy_hydro_->pmy_block;
  if (SELF_GRAVITY_ENABLED && NON_BAROTROPIC_EOS) {
    Gravity *pgrav = pmb->pgrav;
    Real four_pi_G = pgrav->four_pi_G;
    Real grav_mean_rho = pgrav->grav_mean_rho;
    Real phic, phir, phil;
// acceleration in 1-direction
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real dx1 = pmb->pcoord->dx1v(i);
          Real dx2 = pmb->pcoord->dx2v(j);
          Real dx3 = pmb->pcoord->dx3v(k);
          Real dtodx1 = dt/dx1;
          phic = pgrav->phi(k,j,i);
          phil = 0.5*(pgrav->phi(k,j,i-1)+pgrav->phi(k,j,i  ));
          phir = 0.5*(pgrav->phi(k,j,i  )+pgrav->phi(k,j,i+1));
          // Update momenta and energy with d/dx1 terms
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
            Real dx1 = pmb->pcoord->dx1v(i);
            Real dx2 = pmb->pcoord->dx2v(j);
            Real dx3 = pmb->pcoord->dx3v(k);
            Real dtodx2 = dt/dx2;
            phic = pgrav->phi(k,j,i);
            phil = 0.5*(pgrav->phi(k,j-1,i)+pgrav->phi(k,j  ,i));
            phir = 0.5*(pgrav->phi(k,j  ,i)+pgrav->phi(k,j+1,i));
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
            Real dx1 = pmb->pcoord->dx1v(i);
            Real dx2 = pmb->pcoord->dx2v(j);
            Real dx3 = pmb->pcoord->dx3v(k);
            Real dtodx3 = dt/dx3;
            phic = pgrav->phi(k,j,i);
            phil = 0.5*(pgrav->phi(k-1,j,i)+pgrav->phi(k  ,j,i));
            phir = 0.5*(pgrav->phi(k  ,j,i)+pgrav->phi(k+1,j,i));
            cons(IEN,k,j,i) -= dtodx3*(flux[X3DIR](IDN,k  ,j,i)*(phic - phil) +
                                           flux[X3DIR](IDN,k+1,j,i)*(phir - phic));
          }
        }
      }
    }
  }
  return;
}

