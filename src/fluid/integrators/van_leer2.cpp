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

// C++ headers
#include <algorithm>   // min,max

// Primary header
#include "fluid_integrator.hpp"

// Athena headers
#include "../../athena.hpp"                  // enums, macros, Real
#include "../../athena_arrays.hpp"           // AthenaArray
#include "../../coordinates/coordinates.hpp" // Coordinates
#include "../fluid.hpp"                      // Fluid
#include "../../field/field.hpp"             // Fields
#include "../../mesh.hpp"                    // MeshBlock
#include "../srcterms/srcterms.hpp"          // PhysicalSourceTerms()

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//======================================================================================
//! \file van_leer2.cpp
//  \brief van-Leer (MUSCL-Hancock) second-order integrator
//======================================================================================

//--------------------------------------------------------------------------------------
//! \fn  void FluidIntegrator::Predict
//  \brief predictor step for 2nd order VL integrator

void FluidIntegrator::OneStep(MeshBlock *pmb,AthenaArray<Real> &u, AthenaArray<Real> &w,
 InterfaceField &b, AthenaArray<Real> &bcc, const int step)
{
  int tid=0;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int max_nthreads = pmb->pmy_mesh->nthreads_mesh;
 
//  AthenaArray<Real> u = pmb->pfluid->u.ShallowCopy();
//  AthenaArray<Real> w = pmb->pfluid->w.ShallowCopy();
//  AthenaArray<Real> bcc = pmb->pfield->bcc.ShallowCopy();

  AthenaArray<Real> b1,b2,b3,ei_x1f,ei_x2f,ei_x3f,w_x1f,w_x2f,w_x3f;
  b1.InitWithShallowCopy(b.x1f);
  b2.InitWithShallowCopy(b.x2f);
  b3.InitWithShallowCopy(b.x3f);
  ei_x1f.InitWithShallowCopy(pmb->pfield->ei.x1f);
  ei_x2f.InitWithShallowCopy(pmb->pfield->ei.x2f);
  ei_x3f.InitWithShallowCopy(pmb->pfield->ei.x3f);
  w_x1f.InitWithShallowCopy(pmb->pfield->wght.x1f);
  w_x2f.InitWithShallowCopy(pmb->pfield->wght.x2f);
  w_x3f.InitWithShallowCopy(pmb->pfield->wght.x3f);

  Real dt;
  if (step == 1) {
    dt = 0.5*(pmb->pmy_mesh->dt);
  } else {
    dt = (pmb->pmy_mesh->dt);
  }

#pragma omp parallel default(shared) private(tid) num_threads(max_nthreads)
{
#ifdef OPENMP_PARALLEL
  tid=omp_get_thread_num();
#endif
  AthenaArray<Real> wl, wr, flx, area, vol;
  wl.InitWithShallowSlice(wl_,3,tid,1);
  wr.InitWithShallowSlice(wr_,3,tid,1);
  flx.InitWithShallowSlice(flx_,3,tid,1);
  area.InitWithShallowSlice(face_area_,2,tid,1);
  vol.InitWithShallowSlice(cell_volume_,2,tid,1);

//--------------------------------------------------------------------------------------
// i-direction

  for (int k=ks; k<=ke; ++k){ 

#pragma omp for schedule(static)
    for (int j=js; j<=je; ++j){

      if (step == 1) {  // reconstruction in predict step
        for (int n=0; n<NFLUID; ++n){
#pragma simd
        for (int i=is; i<=ie+1; ++i){
          wl(n,i) = w(n,k,j,i-1);
          wr(n,i) = w(n,k,j,i  );
        }}
        if (MAGNETIC_FIELDS_ENABLED) {
#pragma simd
          for (int i=is; i<=ie+1; ++i){
            wl(IBY,i) = bcc(IB2,k,j,i-1);
            wl(IBZ,i) = bcc(IB3,k,j,i-1);
            wr(IBY,i) = bcc(IB2,k,j,i  );
            wr(IBZ,i) = bcc(IB3,k,j,i  );
          }
        }
      } else { // reconstruction in correct step
        for (int n=0; n<NFLUID; ++n) {
          ReconstructionFuncX1(n,n,k,j,w,wl,wr);
        }
        if (MAGNETIC_FIELDS_ENABLED) {
          ReconstructionFuncX1(IB2,IBY,k,j,bcc,wl,wr);
          ReconstructionFuncX1(IB3,IBZ,k,j,bcc,wl,wr);
        }
      }

      RiemannSolver(k,j,is,ie+1,IVX,b1,wl,wr,flx);

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

      if (MAGNETIC_FIELDS_ENABLED) {
#pragma simd
        for (int i=is; i<=ie+1; ++i){
          ei_x1f(X1E3,k,j,i) = -flx(IBY,i); // flx(IBY) = (v1*b2 - v2*b1) = -EMFZ
          ei_x1f(X1E2,k,j,i) =  flx(IBZ,i); // flx(IBZ) = (v1*b3 - v3*b1) =  EMFY
          // estimate weight used to upwind electric fields in GS07 algorithm
          Real fac = (1024)*dt/pmb->pcoord->CenterWidth1(k,j,i);
          Real rat = std::min( 0.5, (fac*flx(IDN,i)/(u(IDN,k,j,i-1)+u(IDN,k,j,i))) );
          w_x1f(k,j,i) = 0.5 + std::max(-0.5,rat);
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

        if (step == 1) {  // reconstruction in predict step
          for (int n=0; n<NFLUID; ++n){
#pragma simd
          for (int i=is; i<=ie; ++i){
            wl(n,i) = w(n,k,j-1,i);
            wr(n,i) = w(n,k,j  ,i);
          }}
          if (MAGNETIC_FIELDS_ENABLED) {
#pragma simd
            for (int i=is; i<=ie; ++i){
              wl(IBY,i) = bcc(IB3,k,j-1,i);
              wl(IBZ,i) = bcc(IB1,k,j-1,i);
              wr(IBY,i) = bcc(IB3,k,j  ,i);
              wr(IBZ,i) = bcc(IB1,k,j  ,i);
            }
          }
        } else {  // reconstruction in correct step
          for (int n=0; n<NFLUID; ++n) {
            ReconstructionFuncX2(n,n,k,j,w,wl,wr);
          }
          if (MAGNETIC_FIELDS_ENABLED) {
            ReconstructionFuncX2(IB3,IBY,k,j,bcc,wl,wr);
            ReconstructionFuncX2(IB1,IBZ,k,j,bcc,wl,wr);
          }
        }

        RiemannSolver(k,j,is,ie,IVY,b2,wl,wr,flx); 

        pmb->pcoord->Face2Area(k,j,is,ie,area);

        if (j>js) {
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

        if (j<(je+1)) {
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

        if (MAGNETIC_FIELDS_ENABLED) {
#pragma simd
          for (int i=is; i<=ie; ++i){
            ei_x2f(X2E1,k,j,i) = -flx(IBY,i); // flx(IBY) = (v2*b3 - v3*b2) = -EMFX
            ei_x2f(X2E3,k,j,i) =  flx(IBZ,i); // flx(IBZ) = (v2*b1 - v1*b2) =  EMFZ
            // estimate weight used to upwind electric fields in GS07 algorithm
            Real fac = (1024)*dt/pmb->pcoord->CenterWidth2(k,j,i);
            Real rat = std::min( 0.5, (fac*flx(IDN,i)/(u(IDN,k,j-1,i)+u(IDN,k,j,i))) );
            w_x2f(k,j,i) = 0.5 + std::max(-0.5,rat);
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

        if (step == 1) {  // reconstruction in predict step
          for (int n=0; n<NFLUID; ++n){
#pragma simd
          for (int i=is; i<=ie; ++i){
            wl(n,i) = w(n,k-1,j,i);
            wr(n,i) = w(n,k  ,j,i);
          }}
          if (MAGNETIC_FIELDS_ENABLED) {
#pragma simd
            for (int i=is; i<=ie; ++i){
              wl(IBY,i) = bcc(IB1,k-1,j,i);
              wl(IBZ,i) = bcc(IB2,k-1,j,i);
              wr(IBY,i) = bcc(IB1,k  ,j,i);
              wr(IBZ,i) = bcc(IB2,k  ,j,i);
            }
          }
        } else {  // reconstruction in correct step
          for (int n=0; n<NFLUID; ++n) {
            ReconstructionFuncX3(n,n,k,j,w,wl,wr);
          }
          if (MAGNETIC_FIELDS_ENABLED) {
            ReconstructionFuncX3(IB1,IBY,k,j,bcc,wl,wr);
            ReconstructionFuncX3(IB2,IBZ,k,j,bcc,wl,wr);
          }
        }

        RiemannSolver(k,j,is,ie,IVZ,b3,wl,wr,flx);

        pmb->pcoord->Face3Area(k,j,is,ie,area);
  
        if (k>ks) {
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

        if (k<(ke+1)) {
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

        if (MAGNETIC_FIELDS_ENABLED) {
#pragma simd
          for (int i=is; i<=ie; ++i){
            ei_x3f(X3E2,k,j,i) = -flx(IBY,i); // flx(IBY) = (v3*b1 - v1*b3) = -EMFY
            ei_x3f(X3E1,k,j,i) =  flx(IBZ,i); // flx(IBZ) = (v3*b2 - v2*b3) =  EMFX
// estimate weight used to upwind electric fields in GS07 algorithm
            Real fac = (1024)*dt/pmb->pcoord->CenterWidth3(k,j,i);
            Real rat = std::min( 0.5, (fac*flx(IDN,i)/(u(IDN,k-1,j,i)+u(IDN,k,j,i))) );
            w_x3f(k,j,i) = 0.5 + std::max(-0.5,rat);
          }
        }

      }
    }
  }

} // end of omp parallel region

//--------------------------------------------------------------------------------------
//  Add source terms for half a timestep

  pmb->pcoord->CoordinateSourceTerms(dt,w,u);
  pmb->pfluid->pf_srcterms->PhysicalSourceTerms(pmb->pmy_mesh->time,dt,w,u);
  if (pmb->pfluid->pf_srcterms->UserSourceTerm != NULL)
    pmb->pfluid->pf_srcterms->UserSourceTerm(pmb->pmy_mesh->time,dt,w,u);

  return;
}
