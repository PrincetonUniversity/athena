//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file add_flux_divergence.cpp
//  \brief Applies divergence of the fluxes, including geometric "source terms" added
//         by a function implemented in each Coordinate class.

// Athena++ headers
#include "hydro.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../field/field.hpp"
#include "../mesh/mesh.hpp"
#include "../bvals/bvals.hpp"
#include "../reconstruct/reconstruction.hpp"

// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn  void Hydro::AddFluxDivergenceToAverage
//  \brief Adds flux divergence to weighted average of conservative variables from
//  previous step(s) of time integrator algorithm

void Hydro::AddFluxDivergenceToAverage(AthenaArray<Real> &w, AthenaArray<Real> &bcc,
                                       const Real wght, AthenaArray<Real> &u_out) {
  MeshBlock *pmb=pmy_block;
  AthenaArray<Real> &x1flux=flux[X1DIR];
  AthenaArray<Real> &x2flux=flux[X2DIR];
  AthenaArray<Real> &x3flux=flux[X3DIR];
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  int tid=0;
  AthenaArray<Real> x1area, x2area, x2area_p1, x3area, x3area_p1, vol, dflx;
  x1area.InitWithShallowCopy(x1face_area_);
  x2area.InitWithShallowCopy(x2face_area_);
  x2area_p1.InitWithShallowCopy(x2face_area_p1_);
  x3area.InitWithShallowCopy(x3face_area_);
  x3area_p1.InitWithShallowCopy(x3face_area_p1_);
  vol.InitWithShallowCopy(cell_volume_);
  dflx.InitWithShallowCopy(dflx_);

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {

      // calculate x1-flux divergence
      pmb->pcoord->Face1Area(k,j,is,ie+1,x1area);
      for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          dflx(n,i) = (x1area(i+1) *x1flux(n,k,j,i+1) - x1area(i)*x1flux(n,k,j,i));
        }
      }

      // calculate x2-flux divergence
      if (pmb->block_size.nx2 > 1) {
        pmb->pcoord->Face2Area(k,j  ,is,ie,x2area   );
        pmb->pcoord->Face2Area(k,j+1,is,ie,x2area_p1);
        for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            dflx(n,i) += (x2area_p1(i)*x2flux(n,k,j+1,i) - x2area(i)*x2flux(n,k,j,i));
          }
        }
      }

      // calculate x3-flux divergence
      if (pmb->block_size.nx3 > 1) {
        pmb->pcoord->Face3Area(k  ,j,is,ie,x3area   );
        pmb->pcoord->Face3Area(k+1,j,is,ie,x3area_p1);
        for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            dflx(n,i) += (x3area_p1(i)*x3flux(n,k+1,j,i) - x3area(i)*x3flux(n,k,j,i));
          }
        }
      }

      // update conserved variables
      pmb->pcoord->CellVolume(k,j,is,ie,vol);
      for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          u_out(n,k,j,i) -= wght*(pmb->pmy_mesh->dt)*dflx(n,i)/vol(i);
        }
      }
    }
  }

  // add coordinate (geometric) source terms
  pmb->pcoord->CoordSrcTerms((wght*pmb->pmy_mesh->dt),pmb->phydro->flux,w,bcc,u_out);

  return;
}

//----------------------------------------------------------------------------------------
//! \fn  void Hydro::WeightedAveU
//  \brief Compute weighted average of cell-averaged U in time integrator step

void Hydro::WeightedAveU(AthenaArray<Real> &u_out, AthenaArray<Real> &u_in1,
                         AthenaArray<Real> &u_in2, const Real wght[3]) {
  MeshBlock *pmb=pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  // consider every possible simplified form of weighted sum operator:
  // U = a*U + b*U1 + c*U2
  // if c=0, c=b=0, or c=b=a=0 (in that order) to avoid extra FMA operations

  // u_in2 may be an unallocated AthenaArray if using a 2S time integrator
  if (wght[2] != 0.0) {
    for (int n=0; n<NHYDRO; ++n) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            u_out(n,k,j,i) = wght[0]*u_out(n,k,j,i) + wght[1]*u_in1(n,k,j,i)
                + wght[2]*u_in2(n,k,j,i);
          }
        }
      }
    }
  } else { // do not dereference u_in2
    if (wght[1] != 0.0) {
      for (int n=0; n<NHYDRO; ++n) {
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
#pragma omp simd
            for (int i=is; i<=ie; ++i) {
              u_out(n,k,j,i) = wght[0]*u_out(n,k,j,i) + wght[1]*u_in1(n,k,j,i);
            }
          }
        }
      }
    } else { // do not dereference u_in1
      if (wght[0] != 0.0) {
        for (int n=0; n<NHYDRO; ++n) {
          for (int k=ks; k<=ke; ++k) {
            for (int j=js; j<=je; ++j) {
#pragma omp simd
              for (int i=is; i<=ie; ++i) {
                u_out(n,k,j,i) = wght[0]*u_out(n,k,j,i);
              }
            }
          }
        }
      } else { // directly initialize u_out to 0
        for (int n=0; n<NHYDRO; ++n) {
          for (int k=ks; k<=ke; ++k) {
            for (int j=js; j<=je; ++j) {
#pragma omp simd
              for (int i=is; i<=ie; ++i) {
                u_out(n,k,j,i) = 0.0;
              }
            }
          }
        }
      }
    }
  }
  return;
}
