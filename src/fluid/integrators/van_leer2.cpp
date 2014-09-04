//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 *
 * This program is free software: you can redistribute and/or modify it under the terms
 * of the GNU General Public License (GPL) as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of GNU GPL in the file LICENSE included in the code
 * distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

#include <stdio.h>

// Primary header
#include "integrators.hpp"

// Athena headers
#include "../../athena.hpp"                   // enums, macros, Real
#include "../../athena_arrays.hpp"            // AthenaArray
#include "../../coordinates/coordinates.hpp"  // Coordinates
#include "../fluid.hpp"                    // Fluid
#include "../../mesh.hpp"                     // MeshBlock

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//======================================================================================
/*! \file van_leer2.cpp
 *  \brief van-Leer (MUSCL-Hancock) second-order integrator
 *====================================================================================*/

//--------------------------------------------------------------------------------------
///*! \fn  void FluidIntegrator::PredictVanLeer2
// *  \brief predictor step for 2nd order VL integrator */

void FluidIntegrator::Predict(MeshBlock *pmb)
{
  int tid=0;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  Real dt = pmb->pmy_domain->pmy_mesh->dt;

  Real sum=0.0; Real sum_2=0.0;
  int ndata = (pmb->block_size.nx1)*(pmb->block_size.nx2)*(pmb->block_size.nx3)*(NVAR);
 
  AthenaArray<Real> u = pmb->pfluid->u.ShallowCopy();
  AthenaArray<Real> w = pmb->pfluid->w.ShallowCopy();
  AthenaArray<Real> u1 = pmb->pfluid->u1.ShallowCopy();
  AthenaArray<Real> w1 = pmb->pfluid->w1.ShallowCopy();

  AthenaArray<Real> src = src_.ShallowCopy();

//--------------------------------------------------------------------------------------
// i-direction 

#pragma omp parallel default(shared) private(tid) num_threads(ATHENA_MAX_NUM_THREADS)
{
#ifdef OPENMP_PARALLEL
  tid=omp_get_thread_num();
#endif
  AthenaArray<Real> *pwl = wl_.ShallowSlice(tid,1);
  AthenaArray<Real> *pwr = wr_.ShallowSlice(tid,1);
  AthenaArray<Real> *pflx = flx_.ShallowSlice(tid,1);
  AthenaArray<Real> *parea = pmb->pcoord->face_area.ShallowSlice(tid,1);
  AthenaArray<Real> *pvol  = pmb->pcoord->cell_volume.ShallowSlice(tid,1);

#pragma omp for schedule(dynamic,4)
  for (int k=ks; k<=ke; ++k){
  for (int j=js; j<=je; ++j){

    ReconstructionFuncX1(k,j,is,ie+1,w,pwl,pwr);

    RiemannSolver(k,j,is,ie+1,IVX,IVY,IVZ,pwl,pwr,pflx);

    pmb->pcoord->Area1Face(k,j,is,ie+1,parea);
    pmb->pcoord->CellVolume(k,j,is,ie,pvol);

    for (int n=0; n<NVAR; ++n){
#pragma simd
      for (int i=is; i<=ie; ++i){
        Real& ui  = u (n,k,j,i);
        Real& u1i = u1(n,k,j,i);
        Real& flxi   = (*pflx)(n,  i);
        Real& flxip1 = (*pflx)(n,i+1);

        Real& area_i   = (*parea)(i);
        Real& area_ip1 = (*parea)(i+1);
        Real& dvol = (*pvol)(i);
 
        u1i = ui - 0.5*dt*(area_ip1*flxip1 - area_i*flxi)/dvol;
      }
    }

  }}
}

//--------------------------------------------------------------------------------------
// j-direction

  if (pmb->block_size.nx2 > 1) {

#pragma omp parallel default(shared) private(tid) num_threads(ATHENA_MAX_NUM_THREADS)
{
#ifdef OPENMP_PARALLEL
    tid=omp_get_thread_num();
#endif
    AthenaArray<Real> *pwl = wl_.ShallowSlice(tid,1);
    AthenaArray<Real> *pwr = wr_.ShallowSlice(tid,1);
    AthenaArray<Real> *pflx = flx_.ShallowSlice(tid,1);
    AthenaArray<Real> *parea = pmb->pcoord->face_area.ShallowSlice(tid,1);
    AthenaArray<Real> *pvol  = pmb->pcoord->cell_volume.ShallowSlice(tid,1);

#pragma omp for schedule(dynamic,4)
    for (int k=ks; k<=ke; ++k){
    for (int j=js; j<=je+1; ++j){

      ReconstructionFuncX2(k,j,is,ie,w,pwl,pwr);

      RiemannSolver(k,j,is,ie,IVY,IVZ,IVX,pwl,pwr,pflx); 

      pmb->pcoord->Area2Face(k,j,is,ie,parea);

      if (j>js) {
        pmb->pcoord->CellVolume(k,j-1,is,ie,pvol);
        for (int n=0; n<NVAR; ++n){
#pragma simd
          for (int i=is; i<=ie; ++i){
            Real& u1jm1 = u1(n,k,j-1,i);
            Real& flxj  = (*pflx)(n,i);
            Real& area_i   = (*parea)(i);
            Real& dvol = (*pvol)(i);
  
            u1jm1 -= 0.5*dt*area_i*flxj/dvol;
          }
        }
      }

      if (j<(je+1)) {
        pmb->pcoord->CellVolume(k,j,is,ie,pvol);
        for (int n=0; n<NVAR; ++n){
#pragma simd
          for (int i=is; i<=ie; ++i){
            Real& u1j   = u1(n,k,j  ,i);
            Real& flxj  = (*pflx)(n,i);
            Real& area_i   = (*parea)(i);
            Real& dvol = (*pvol)(i);

            u1j   += 0.5*dt*area_i*flxj/dvol;
          }
        }
      }

    }}
}
  }

//--------------------------------------------------------------------------------------
// k-direction 

  if (pmb->block_size.nx3 > 1) {

#pragma omp parallel default(shared) private(tid) num_threads(ATHENA_MAX_NUM_THREADS)
{
#ifdef OPENMP_PARALLEL
    tid=omp_get_thread_num();
#endif
    AthenaArray<Real> *pwl = wl_.ShallowSlice(tid,1);
    AthenaArray<Real> *pwr = wr_.ShallowSlice(tid,1);
    AthenaArray<Real> *pflx = flx_.ShallowSlice(tid,1);
    AthenaArray<Real> *parea = pmb->pcoord->face_area.ShallowSlice(tid,1);
    AthenaArray<Real> *pvol  = pmb->pcoord->cell_volume.ShallowSlice(tid,1);

#pragma omp for schedule(dynamic,4)
    for (int k=ks; k<=ke+1; ++k){
    for (int j=js; j<=je; ++j){

      ReconstructionFuncX3(k,j,is,ie,w,pwl,pwr);

      RiemannSolver(k,j,is,ie,IVZ,IVX,IVY,pwl,pwr,pflx);

      pmb->pcoord->Area3Face(k,j,is,ie,parea);

      if (k>ks) {
        pmb->pcoord->CellVolume(k-1,j,is,ie,pvol);
        for (int n=0; n<NVAR; ++n){
#pragma simd
          for (int i=is; i<=ie; ++i){
            Real& u1km1 = u1(n,k-1,j,i);
            Real& flxk = (*pflx)(n,i);
            Real& area_i   = (*parea)(i);
            Real& dvol = (*pvol)(i);

            u1km1 -= 0.5*dt*area_i*flxk/dvol;
          }
        }
      }

      if (k<(ke+1)) {
        pmb->pcoord->CellVolume(k,j,is,ie,pvol);
        for (int n=0; n<NVAR; ++n){
#pragma simd
          for (int i=is; i<=ie; ++i){
            Real& u1k   = u1(n,k  ,j,i);
            Real& flxk = (*pflx)(n,i);
            Real& area_i   = (*parea)(i);
            Real& dvol = (*pvol)(i);

            u1k += 0.5*dt*area_i*flxk/dvol;
          }
        }
      }

    }}
}
  }

//--------------------------------------------------------------------------------------
//  Add source terms for half a timestep

  for (int k=ks; k<=ke; ++k){
  for (int j=js; j<=je; ++j){

    pmb->pcoord->CoordinateSourceTerms(k,j,w,src);

#pragma simd
    for (int i=is; i<=ie; ++i){
      Real& u1im1   = u1(IM1,k,j,i);
      Real& u1im2   = u1(IM2,k,j,i);
      Real& u1im3   = u1(IM3,k,j,i);

      u1im1 += 0.5*dt*src(IM1,i);
      u1im2 += 0.5*dt*src(IM2,i);
      u1im3 += 0.5*dt*src(IM3,i);
    }
  }}

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

  Real sum=0.0; Real sum_2=0.0;
  int ndata = (pmb->block_size.nx1)*(pmb->block_size.nx2)*(pmb->block_size.nx3)*(NVAR);
 
  AthenaArray<Real> u = pmb->pfluid->u.ShallowCopy();
  AthenaArray<Real> w = pmb->pfluid->w.ShallowCopy();
  AthenaArray<Real> u1 = pmb->pfluid->u1.ShallowCopy();
  AthenaArray<Real> w1 = pmb->pfluid->w1.ShallowCopy();

  AthenaArray<Real> src = src_.ShallowCopy();
 
//--------------------------------------------------------------------------------------
// i-direction 

#pragma omp parallel default(shared) private(tid) num_threads(ATHENA_MAX_NUM_THREADS)
{
#ifdef OPENMP_PARALLEL
  tid=omp_get_thread_num();;
#endif
  AthenaArray<Real> *pwl = wl_.ShallowSlice(tid,1);
  AthenaArray<Real> *pwr = wr_.ShallowSlice(tid,1);
  AthenaArray<Real> *pflx = flx_.ShallowSlice(tid,1);
  AthenaArray<Real> *parea = pmb->pcoord->face_area.ShallowSlice(tid,1);
  AthenaArray<Real> *pvol  = pmb->pcoord->cell_volume.ShallowSlice(tid,1);

#pragma omp for schedule(dynamic,4)
  for (int k=ks; k<=ke; ++k){
  for (int j=js; j<=je; ++j){

    ReconstructionFuncX1(k,j,is,ie+1,w1,pwl,pwr);

    RiemannSolver(k,j,is,ie+1,IVX,IVY,IVZ,pwl,pwr,pflx); 

    pmb->pcoord->Area1Face(k,j,is,ie+1,parea);
    pmb->pcoord->CellVolume(k,j,is,ie,pvol);

    for (int n=0; n<NVAR; ++n){
#pragma simd
      for (int i=is; i<=ie; ++i){
        Real& ui  = u (n,k,j,i);
        Real& flxi   = (*pflx)(n,  i);
        Real& flxip1 = (*pflx)(n,i+1);

        Real& area_i   = (*parea)(i);
        Real& area_ip1 = (*parea)(i+1);
        Real& dvol = (*pvol)(i);
 
        ui -= dt*(area_ip1*flxip1 - area_i*flxi)/dvol;
      }
    }

  }}
}

//--------------------------------------------------------------------------------------
// j-direction

  if (pmb->block_size.nx2 > 1) {

#pragma omp parallel default(shared) private(tid) num_threads(ATHENA_MAX_NUM_THREADS)
{
#ifdef OPENMP_PARALLEL
    tid=omp_get_thread_num();;
#endif
    AthenaArray<Real> *pwl = wl_.ShallowSlice(tid,1);
    AthenaArray<Real> *pwr = wr_.ShallowSlice(tid,1);
    AthenaArray<Real> *pflx = flx_.ShallowSlice(tid,1);
    AthenaArray<Real> *parea = pmb->pcoord->face_area.ShallowSlice(tid,1);
    AthenaArray<Real> *pvol  = pmb->pcoord->cell_volume.ShallowSlice(tid,1);

#pragma omp for schedule(dynamic,4)
    for (int k=ks; k<=ke; ++k){
    for (int j=js; j<=je+1; ++j){

      ReconstructionFuncX2(k,j,is,ie,w1,pwl,pwr);

      RiemannSolver(k,j,is,ie,IVY,IVZ,IVX,pwl,pwr,pflx); 

      pmb->pcoord->Area2Face(k,j,is,ie,parea);

      if (j>js){
        pmb->pcoord->CellVolume(k,j-1,is,ie,pvol);
        for (int n=0; n<NVAR; ++n){
#pragma simd
          for (int i=is; i<=ie; ++i){
            Real& ujm1 = u(n,k,j-1,i);
            Real& flxj  = (*pflx)(n,i);
            Real& area_i   = (*parea)(i);
            Real& dvol = (*pvol)(i);
  
            ujm1 -= dt*area_i*flxj/dvol;
          }
        }
      }

      if (j<(je+1)){
        pmb->pcoord->CellVolume(k,j,is,ie,pvol);
        for (int n=0; n<NVAR; ++n){
#pragma simd
          for (int i=is; i<=ie; ++i){
            Real& uj   = u(n,k,j  ,i);
            Real& flxj  = (*pflx)(n,i);
            Real& area_i   = (*parea)(i);
            Real& dvol = (*pvol)(i);
  
            uj += dt*area_i*flxj/dvol;
          }
        }
      }

    }}
}
  }

//--------------------------------------------------------------------------------------
// k-direction 

  if (pmb->block_size.nx3 > 1) {

#pragma omp parallel default(shared) private(tid) num_threads(ATHENA_MAX_NUM_THREADS)
{
#ifdef OPENMP_PARALLEL
    tid=omp_get_thread_num();;
#endif
    AthenaArray<Real> *pwl = wl_.ShallowSlice(tid,1);
    AthenaArray<Real> *pwr = wr_.ShallowSlice(tid,1);
    AthenaArray<Real> *pflx = flx_.ShallowSlice(tid,1);
    AthenaArray<Real> *parea = pmb->pcoord->face_area.ShallowSlice(tid,1);
    AthenaArray<Real> *pvol  = pmb->pcoord->cell_volume.ShallowSlice(tid,1);

#pragma omp for schedule(dynamic,4)
    for (int k=ks; k<=ke+1; ++k){
    for (int j=js; j<=je; ++j){

      ReconstructionFuncX3(k,j,is,ie,w1,pwl,pwr);

      RiemannSolver(k,j,is,ie,IVZ,IVX,IVY,pwl,pwr,pflx);

      pmb->pcoord->Area3Face(k,j,is,ie,parea);

      if (k>ks){
        pmb->pcoord->CellVolume(k-1,j,is,ie,pvol);
        for (int n=0; n<NVAR; ++n){
#pragma simd
          for (int i=is; i<=ie; ++i){
            Real& ukm1 = u(n,k-1,j,i);
            Real& flxk = (*pflx)(n,i);
            Real& area_i   = (*parea)(i);
            Real& dvol = (*pvol)(i);

            ukm1 -= dt*area_i*flxk/dvol;
          }
        }
      }

      if (k<(ke+1)){
        pmb->pcoord->CellVolume(k,j,is,ie,pvol);
        for (int n=0; n<NVAR; ++n){
#pragma simd
          for (int i=is; i<=ie; ++i){
            Real& uk   = u(n,k  ,j,i);
            Real& flxk = (*pflx)(n,i);
            Real& area_i   = (*parea)(i);
            Real& dvol = (*pvol)(i);

            uk   += dt*area_i*flxk/dvol;
          }
        }
      }

    }}
}
  }

//--------------------------------------------------------------------------------------
//  Add source terms for a full timestep

  for (int k=ks; k<=ke; ++k){
  for (int j=js; j<=je; ++j){

    pmb->pcoord->CoordinateSourceTerms(k,j,w1,src);

#pragma simd
    for (int i=is; i<=ie; ++i){
      Real& uim1   = u(IM1,k,j,i);
      Real& uim2   = u(IM2,k,j,i);
      Real& uim3   = u(IM3,k,j,i);

      uim1 += dt*src(IM1,i);
      uim2 += dt*src(IM2,i);
      uim3 += dt*src(IM3,i);
    }
  }}

  return;
}
