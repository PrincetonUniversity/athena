//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file add_flux_divergence.cpp
//! \brief Computes divergence of the Hydro fluxes and adds that to a conserved variable
//! register

// C headers

// C++ headers
#include <algorithm>  // std::binary_search
#include <vector>     // std::vector

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "hydro.hpp"

// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn  void Hydro::AddFluxDivergence
//! \brief Adds flux divergence to weighted average of conservative variables from
//! previous step(s) of time integrator algorithm

// TODO(felker): consider combining with PassiveScalars implementation + (see 57cfe28b)
// (may rename to AddPhysicalFluxDivergence or AddQuantityFluxDivergence to explicitly
// distinguish from CoordTerms)
// (may rename to AddHydroFluxDivergence and AddScalarsFluxDivergence, if
// the implementations remain completely independent / no inheritance is
// used)
void Hydro::AddFluxDivergence(const Real wght, AthenaArray<Real> &u_out) {
  MeshBlock *pmb = pmy_block;
  AthenaArray<Real> &x1flux = flux[X1DIR];
  AthenaArray<Real> &x2flux = flux[X2DIR];
  AthenaArray<Real> &x3flux = flux[X3DIR];
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  AthenaArray<Real> &x1area = x1face_area_, &x2area = x2face_area_,
                 &x2area_p1 = x2face_area_p1_, &x3area = x3face_area_,
                 &x3area_p1 = x3face_area_p1_, &vol = cell_volume_, &dflx = dflx_;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      // calculate x1-flux divergence
      pmb->pcoord->Face1Area(k, j, is, ie+1, x1area);
      for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          dflx(n,i) = (x1area(i+1)*x1flux(n,k,j,i+1) - x1area(i)*x1flux(n,k,j,i));
        }
      }

      // calculate x2-flux divergence
      if (pmb->block_size.nx2 > 1) {
        pmb->pcoord->Face2Area(k, j  , is, ie, x2area   );
        pmb->pcoord->Face2Area(k, j+1, is, ie, x2area_p1);
        for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            dflx(n,i) += (x2area_p1(i)*x2flux(n,k,j+1,i) - x2area(i)*x2flux(n,k,j,i));
          }
        }
      }

      // calculate x3-flux divergence
      if (pmb->block_size.nx3 > 1) {
        pmb->pcoord->Face3Area(k  , j, is, ie, x3area   );
        pmb->pcoord->Face3Area(k+1, j, is, ie, x3area_p1);
        for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            dflx(n,i) += (x3area_p1(i)*x3flux(n,k+1,j,i) - x3area(i)*x3flux(n,k,j,i));
          }
        }
      }

      // update conserved variables
      pmb->pcoord->CellVolume(k, j, is, ie, vol);
      for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          u_out(n,k,j,i) -= wght*dflx(n,i)/vol(i);
        }
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn  void Hydro::AddFluxDivergence_STS
//! \brief Adds diffusive flux divergence to weighted average of conservative variables
//! from previous step of time integration algorithm.  RKL2 registers are set to flux
//! divergence update if first stage of RKL2 STS.  Only a subset of ph->u indices are
//! integrated, dependent on which diffusive processes are integrated.
void Hydro::AddFluxDivergence_STS(const Real wght, int stage,
                                  AthenaArray<Real> &u_out,
                                  AthenaArray<Real> &fl_div_out,
                                  std::vector<int> idx_subset) {
  MeshBlock *pmb = pmy_block;
  AthenaArray<Real> &x1flux = flux[X1DIR];
  AthenaArray<Real> &x2flux = flux[X2DIR];
  AthenaArray<Real> &x3flux = flux[X3DIR];
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  AthenaArray<Real> &x1area = x1face_area_, &x2area = x2face_area_,
                 &x2area_p1 = x2face_area_p1_, &x3area = x3face_area_,
                 &x3area_p1 = x3face_area_p1_, &vol = cell_volume_, &dflx = dflx_;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      // calculate x1-flux divergence
      pmb->pcoord->Face1Area(k, j, is, ie+1, x1area);
      for (int n=0; n<NHYDRO; ++n) {
        if (std::binary_search(idx_subset.begin(), idx_subset.end(), n)) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            dflx(n,i) = (x1area(i+1) *x1flux(n,k,j,i+1) - x1area(i)*x1flux(n,k,j,i));
          }
        }
      }

      // calculate x2-flux divergence
      if (pmb->block_size.nx2 > 1) {
        pmb->pcoord->Face2Area(k, j  , is, ie, x2area   );
        pmb->pcoord->Face2Area(k, j+1, is, ie, x2area_p1);
        for (int n=0; n<NHYDRO; ++n) {
          if (std::binary_search(idx_subset.begin(), idx_subset.end(), n)) {
#pragma omp simd
            for (int i=is; i<=ie; ++i) {
              dflx(n,i) += (x2area_p1(i)*x2flux(n,k,j+1,i) - x2area(i)*x2flux(n,k,j,i));
            }
          }
        }
      }

      // calculate x3-flux divergence
      if (pmb->block_size.nx3 > 1) {
        pmb->pcoord->Face3Area(k  , j, is, ie, x3area   );
        pmb->pcoord->Face3Area(k+1, j, is, ie, x3area_p1);
        for (int n=0; n<NHYDRO; ++n) {
          if (std::binary_search(idx_subset.begin(), idx_subset.end(), n)) {
#pragma omp simd
            for (int i=is; i<=ie; ++i) {
              dflx(n,i) += (x3area_p1(i)*x3flux(n,k+1,j,i) - x3area(i)*x3flux(n,k,j,i));
            }
          }
        }
      }

      // update conserved variables
      pmb->pcoord->CellVolume(k, j, is, ie, vol);
      for (int n=0; n<NHYDRO; ++n) {
        if (std::binary_search(idx_subset.begin(), idx_subset.end(), n)) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            u_out(n,k,j,i) -= wght*dflx(n,i)/vol(i);
            if (stage == 1 && pmb->pmy_mesh->sts_integrator == "rkl2") {
              fl_div_out(n,k,j,i) = -0.5*pmb->pmy_mesh->dt*dflx(n,i)/vol(i);
            }
          }
        }
      }
    }
  }
  return;
}
