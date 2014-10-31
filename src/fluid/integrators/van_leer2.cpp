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
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================

#include <stdio.h>
#include <iostream>

// Primary header
#include "integrators.hpp"

// Athena headers
#include "../../athena.hpp"                  // enums, macros, Real
#include "../../athena_arrays.hpp"           // AthenaArray
#include "../../coordinates/coordinates.hpp" // Coordinates
#include "../fluid.hpp"                      // Fluid
#include "../../mesh.hpp"                    // MeshBlock
#include "../srcterms/srcterms.hpp"          // PhysicalSourceTerms()

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//======================================================================================
//! \file van_leer2.cpp
//  \brief van-Leer (MUSCL-Hancock) second-order integrator
//======================================================================================

//--------------------------------------------------------------------------------------
//! \fn  void FluidIntegrator::PredictVanLeer2
//  \brief predictor step for 2nd order VL integrator

void FluidIntegrator::Predict(MeshBlock *pmb)
{
  int tid=0;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  Real dt = pmb->pmy_domain->pmy_mesh->dt;
  int thread_max = pmb->pmy_domain->pmy_mesh->nthreads_mesh;
 
  AthenaArray<Real> u = pmb->pfluid->u.ShallowCopy();
  AthenaArray<Real> w = pmb->pfluid->w.ShallowCopy();
  AthenaArray<Real> u1 = pmb->pfluid->u1.ShallowCopy();
  AthenaArray<Real> w1 = pmb->pfluid->w1.ShallowCopy();

#pragma omp parallel default(shared) private(tid) num_threads(thread_max)
{
#ifdef OPENMP_PARALLEL
  tid=omp_get_thread_num();
#endif
  AthenaArray<Real> wl = wl_.ShallowSlice(tid,1);
  AthenaArray<Real> wr = wr_.ShallowSlice(tid,1);
  AthenaArray<Real> flx = flx_.ShallowSlice(tid,1);
  AthenaArray<Real> area = face_area_.ShallowSlice(tid,1);
  AthenaArray<Real> vol  = cell_volume_.ShallowSlice(tid,1);

//--------------------------------------------------------------------------------------
// i-direction 

  for (int k=ks; k<=ke; ++k){

#pragma omp for schedule(static)
    for (int j=js; j<=je; ++j){

      for (int n=0; n<NFLUID; ++n){
#pragma simd
        for (int i=is; i<=ie+1; ++i){
          wl(n,i) = w(n,k,j,i-1);
          wr(n,i) = w(n,k,j,i  );
        }
      }

      RiemannSolver(k,j,is,ie+1,IVX,IVY,IVZ,wl,wr,flx);

      pmb->pcoord->Face1Area(k,j,is,ie+1,area);
      pmb->pcoord->CellVolume(k,j,is,ie,vol);

      for (int n=0; n<NFLUID; ++n){
#pragma simd
        for (int i=is; i<=ie; ++i){
          Real& ui  = u (n,k,j,i);
          Real& flxi   = flx(n,  i);
          Real& flxip1 = flx(n,i+1);

          Real& area_i   = area(i);
          Real& area_ip1 = area(i+1);
          Real& dvol = vol(i);
 
          u1(n,k,j,i) = ui - 0.5*dt*(area_ip1*flxip1 - area_i*flxi)/dvol;
        }
      }
    }
  }

//--------------------------------------------------------------------------------------
// j-direction

  if (pmb->block_size.nx2 > 1) {

    for (int k=ks; k<=ke; ++k){

#pragma omp for schedule(static)
      for (int j=js; j<=je+1; ++j){

        for (int n=0; n<NFLUID; ++n){
#pragma simd
          for (int i=is; i<=ie; ++i){
            wl(n,i) = w(n,k,j-1,i);
            wr(n,i) = w(n,k,j  ,i);
          }
        }

        RiemannSolver(k,j,is,ie,IVY,IVZ,IVX,wl,wr,flx); 

        pmb->pcoord->Face2Area(k,j,is,ie,area);

        if (j>js) {
          pmb->pcoord->CellVolume(k,j-1,is,ie,vol);
          for (int n=0; n<NFLUID; ++n){
#pragma simd
            for (int i=is; i<=ie; ++i){
              Real& flxj  = flx(n,i);
              Real& area_i   = area(i);
              Real& dvol = vol(i);
  
              u1(n,k,j-1,i) -= 0.5*dt*area_i*flxj/dvol;
            }
          }
        }

        if (j<(je+1)) {
          pmb->pcoord->CellVolume(k,j,is,ie,vol);
          for (int n=0; n<NFLUID; ++n){
#pragma simd
            for (int i=is; i<=ie; ++i){
              Real& flxj  = flx(n,i);
              Real& area_i   = area(i);
              Real& dvol = vol(i);

              u1(n,k,j,i) += 0.5*dt*area_i*flxj/dvol;
            }
          }
        }

      }
    }
  }

//--------------------------------------------------------------------------------------
// k-direction 

  if (pmb->block_size.nx3 > 1) {

    for (int k=ks; k<=ke+1; ++k){

#pragma omp for schedule(static)
      for (int j=js; j<=je; ++j){

        for (int n=0; n<NFLUID; ++n){
#pragma simd
          for (int i=is; i<=ie; ++i){
            wl(n,i) = w(n,k-1,j,i);
            wr(n,i) = w(n,k  ,j,i);
          }
        }

        RiemannSolver(k,j,is,ie,IVZ,IVX,IVY,wl,wr,flx);

        pmb->pcoord->Face3Area(k,j,is,ie,area);
  
        if (k>ks) {
          pmb->pcoord->CellVolume(k-1,j,is,ie,vol);
          for (int n=0; n<NFLUID; ++n){
#pragma simd
            for (int i=is; i<=ie; ++i){
              Real& flxk = flx(n,i);
              Real& area_i   = area(i);
              Real& dvol = vol(i);

              u1(n,k-1,j,i) -= 0.5*dt*area_i*flxk/dvol;
            }
          }
        }

        if (k<(ke+1)) {
          pmb->pcoord->CellVolume(k,j,is,ie,vol);
          for (int n=0; n<NFLUID; ++n){
#pragma simd
            for (int i=is; i<=ie; ++i){
              Real& flxk = flx(n,i);
              Real& area_i   = area(i);
              Real& dvol = vol(i);

              u1(n,k,j,i) += 0.5*dt*area_i*flxk/dvol;
            }
          }
        }

      }
    }
  }

} // end of omp parallel region

//--------------------------------------------------------------------------------------
//  Add source terms for half a timestep

  pmb->pcoord->CoordinateSourceTerms(0.5*dt,w,u1);
  pmb->pfluid->pf_srcterms->PhysicalSourceTerms(pmb->pmy_domain->pmy_mesh->time,0.5*dt,w,u1);
  if (pmb->pfluid->pf_srcterms->UserSourceTerm != NULL)
    pmb->pfluid->pf_srcterms->UserSourceTerm(pmb->pmy_domain->pmy_mesh->time,0.5*dt,w,u1);

  return;
}

//--------------------------------------------------------------------------------------
// \!fn 
// \brief

void FluidIntegrator::Correct(MeshBlock *pmb)
{
  int tid=0;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  Real dt = pmb->pmy_domain->pmy_mesh->dt;
  int thread_max = pmb->pmy_domain->pmy_mesh->nthreads_mesh;

  AthenaArray<Real> u = pmb->pfluid->u.ShallowCopy();
  AthenaArray<Real> w = pmb->pfluid->w.ShallowCopy();
  AthenaArray<Real> u1 = pmb->pfluid->u1.ShallowCopy();
  AthenaArray<Real> w1 = pmb->pfluid->w1.ShallowCopy();

#pragma omp parallel default(shared) private(tid) num_threads(thread_max)
{
#ifdef OPENMP_PARALLEL
  tid=omp_get_thread_num();
#endif
  AthenaArray<Real> wl = wl_.ShallowSlice(tid,1);
  AthenaArray<Real> wr = wr_.ShallowSlice(tid,1);
  AthenaArray<Real> flx = flx_.ShallowSlice(tid,1);
  AthenaArray<Real> area = face_area_.ShallowSlice(tid,1);
  AthenaArray<Real> vol  = cell_volume_.ShallowSlice(tid,1);

//--------------------------------------------------------------------------------------
// i-direction 

  for (int k=ks; k<=ke; ++k){

#pragma omp for schedule(static)
    for (int j=js; j<=je; ++j){

      ReconstructionFuncX1(k,j,is,ie+1,w1,wl,wr);

      RiemannSolver(k,j,is,ie+1,IVX,IVY,IVZ,wl,wr,flx); 

      pmb->pcoord->Face1Area(k,j,is,ie+1,area);
      pmb->pcoord->CellVolume(k,j,is,ie,vol);

      for (int n=0; n<NFLUID; ++n){
#pragma simd
        for (int i=is; i<=ie; ++i){
          Real& flxi   = flx(n,  i);
          Real& flxip1 = flx(n,i+1);

          Real& area_i   = area(i);
          Real& area_ip1 = area(i+1);
          Real& dvol = vol(i);
   
          u(n,k,j,i) -= dt*(area_ip1*flxip1 - area_i*flxi)/dvol;
        }
      }
    }
  }

//--------------------------------------------------------------------------------------
// j-direction

  if (pmb->block_size.nx2 > 1) {

    for (int k=ks; k<=ke; ++k){

#pragma omp for schedule(static)
      for (int j=js; j<=je+1; ++j){

        ReconstructionFuncX2(k,j,is,ie,w1,wl,wr);

        RiemannSolver(k,j,is,ie,IVY,IVZ,IVX,wl,wr,flx); 

        pmb->pcoord->Face2Area(k,j,is,ie,area);

        if (j>js){
          pmb->pcoord->CellVolume(k,j-1,is,ie,vol);
          for (int n=0; n<NFLUID; ++n){
#pragma simd
            for (int i=is; i<=ie; ++i){
              Real& flxj  = flx(n,i);
              Real& area_i   = area(i);
              Real& dvol = vol(i);
  
              u(n,k,j-1,i) -= dt*area_i*flxj/dvol;
            }
          }
        }

        if (j<(je+1)){
          pmb->pcoord->CellVolume(k,j,is,ie,vol);
          for (int n=0; n<NFLUID; ++n){
#pragma simd
            for (int i=is; i<=ie; ++i){
              Real& flxj  = flx(n,i);
              Real& area_i   = area(i);
              Real& dvol = vol(i);
  
              u(n,k,j,i) += dt*area_i*flxj/dvol;
            }
          }
        }

      }
    }
  }

//--------------------------------------------------------------------------------------
// k-direction 

  if (pmb->block_size.nx3 > 1) {

    for (int k=ks; k<=ke+1; ++k){

#pragma omp for schedule(static)
      for (int j=js; j<=je; ++j){

        ReconstructionFuncX3(k,j,is,ie,w1,wl,wr);

        RiemannSolver(k,j,is,ie,IVZ,IVX,IVY,wl,wr,flx);

        pmb->pcoord->Face3Area(k,j,is,ie,area);

        if (k>ks){
          pmb->pcoord->CellVolume(k-1,j,is,ie,vol);
          for (int n=0; n<NFLUID; ++n){
#pragma simd
            for (int i=is; i<=ie; ++i){
              Real& flxk = flx(n,i);
              Real& area_i   = area(i);
              Real& dvol = vol(i);
  
              u(n,k-1,j,i) -= dt*area_i*flxk/dvol;
            }
          }
        }

        if (k<(ke+1)){
          pmb->pcoord->CellVolume(k,j,is,ie,vol);
          for (int n=0; n<NFLUID; ++n){
#pragma simd
            for (int i=is; i<=ie; ++i){
              Real& flxk = flx(n,i);
              Real& area_i   = area(i);
              Real& dvol = vol(i);

              u(n,k,j,i) += dt*area_i*flxk/dvol;
            }
          }
        }

      }
    }
  }

} // end of omp parallel region

//--------------------------------------------------------------------------------------
//  Add source terms for a full timestep

  pmb->pcoord->CoordinateSourceTerms(dt,w1,u);
  pmb->pfluid->pf_srcterms->PhysicalSourceTerms(pmb->pmy_domain->pmy_mesh->time,dt,w1,u);
  if (pmb->pfluid->pf_srcterms->UserSourceTerm != NULL)
    pmb->pfluid->pf_srcterms->UserSourceTerm(pmb->pmy_domain->pmy_mesh->time,dt,w1,u);

  return;
}
