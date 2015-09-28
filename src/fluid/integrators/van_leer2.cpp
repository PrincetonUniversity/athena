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
//! \file van_leer2.cpp
//  \brief van-Leer (MUSCL-Hancock) second-order integrator
//======================================================================================

// C/C++ headers
#include <algorithm>   // min,max

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../fluid.hpp"
#include "../../field/field.hpp"
#include "../../mesh.hpp"
#include "../srcterms/srcterms.hpp"
#include "../../bvals/bvals.hpp"
#include "../viscosity/viscosity.hpp"

// this class header
#include "fluid_integrator.hpp"

// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//--------------------------------------------------------------------------------------
//! \fn  void HydroIntegrator::Predict
//  \brief predictor step for 2nd order VL integrator

void HydroIntegrator::OneStep(MeshBlock *pmb,AthenaArray<Real> &u, AthenaArray<Real> &w,
 InterfaceField &b, AthenaArray<Real> &bcc, const int step)
{
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int il, iu, jl, ju, kl, ku;

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

  int tid=0;
  int nthreads = pmb->pmy_mesh->GetNumMeshThreads();
#pragma omp parallel default(shared) private(tid) num_threads(nthreads)
{
#ifdef OPENMP_PARALLEL
  tid=omp_get_thread_num();
#endif
  AthenaArray<Real> wl, wr, flx, jflx_j, kflx_k, area, area_p1, vol;
  wl.InitWithShallowSlice(wl_,3,tid,1);
  wr.InitWithShallowSlice(wr_,3,tid,1);
  flx.InitWithShallowSlice(flx_,3,tid,1);
  jflx_j.InitWithShallowSlice(jflx_j_,3,tid,1);
  kflx_k.InitWithShallowSlice(kflx_k_,4,tid,1);
  area.InitWithShallowSlice(face_area_,2,tid,1);
  area_p1.InitWithShallowSlice(face_area_p1_,2,tid,1);
  vol.InitWithShallowSlice(cell_volume_,2,tid,1);

//----------Viscosity update
  if(VISCOSITY) pmb->pfluid->pf_viscosity->ViscosityTerms(dt,w,u);

//--------------------------------------------------------------------------------------
// i-direction
  // set the loop limits
  jl=js, ju=je, kl=ks, ku=ke;
  if (MAGNETIC_FIELDS_ENABLED) {
    if(pmb->block_size.nx2 > 1) {
      if(pmb->block_size.nx3 == 1) // 2D
        jl=js-1, ju=je+1, kl=ks, ku=ke;
      else // 3D
        jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
    }
  }
  for (int k=kl; k<=ku; ++k){ 
#pragma omp for schedule(static)
    for (int j=jl; j<=ju; ++j){

      // reconstruct L/R states
      if (step == 1) {
        DonorCellX1(k,j,is,ie+1,w,bcc,wl,wr);
      } else {
        PiecewiseLinearX1(k,j,is,ie+1,w,bcc,wl,wr);
      }

      // compute fluxes
      RiemannSolver(k,j,is,ie+1,IVX,b1,wl,wr,flx);

      if(k>=ks && k<=ke && j>=js && j<=je) {
        // update conserved fluid variables
        pmb->pcoord->Face1Area(k,j,is,ie+1,area);
        pmb->pcoord->CellVolume(k,j,is,ie,vol);
        for (int n=0; n<NFLUID; ++n){
#pragma simd
          for (int i=is; i<=ie; ++i){
            u(n,k,j,i) -= dt*(area(i+1)*flx(n,i+1) - area(i)*flx(n,i))/vol(i);
          }
        }

        // add coordinate (geometric) source terms
        pmb->pcoord->CoordSrcTermsX1(k,j,dt,flx,w,bcc,u);
        // add physical source terms for a point mass potential
        pmb->pfluid->pf_srcterms->PhysicalSourceTermsX1(k,j,dt,flx,w,u);

        // store the surface flux for flux correction
        if(pmb->pmy_mesh->multilevel==true) {
          for(int n=0; n<NFLUID; ++n) {
            pmb->pbval->surface_flux_[inner_x1](n,k,j)=flx(n,is);
            pmb->pbval->surface_flux_[outer_x1](n,k,j)=flx(n,ie+1);
          }
        }
      }

      // store electric fields, compute weights for GS07 CT algorithm
      if (MAGNETIC_FIELDS_ENABLED) {
#pragma simd
        for (int i=is; i<=ie+1; ++i){
          ei_x1f(X1E3,k,j,i) = -flx(IBY,i); // flx(IBY) = (v1*b2 - v2*b1) = -EMFZ
          ei_x1f(X1E2,k,j,i) =  flx(IBZ,i); // flx(IBZ) = (v1*b3 - v3*b1) =  EMFY
          const Real& dx = pmb->pcoord->CenterWidth1(k,j,i);
          Real v_over_c = (1024)*dt*flx(IDN,i)/(dx*(wl(IDN,i) + wr(IDN,i)));
          Real tmp_min = std::min(0.5,v_over_c);
          w_x1f(k,j,i) = 0.5 + std::max(-0.5,tmp_min);
        }
      }

    }
  }

//--------------------------------------------------------------------------------------
// j-direction

  if (pmb->block_size.nx2 > 1) {
    // set the loop limits
    il=is, iu=ie, kl=ks, ku=ke;
    if (MAGNETIC_FIELDS_ENABLED) {
      if(pmb->block_size.nx3 == 1) // 2D
        il=is-1, iu=ie+1, kl=ks, ku=ke;
      else // 3D
        il=is-1, iu=ie+1, kl=ks-1, ku=ke+1;
    }
    for (int k=kl; k<=ku; ++k){
      bool first_time_through_loop = true;
#pragma omp for schedule(static)
      for (int j=js; j<=je; ++j){

        if (first_time_through_loop) {
          // reconstruct L/R states at j=jstart
          if (step == 1) {
            DonorCellX2(k,j,il,iu,w,bcc,wl,wr);
          } else {
            PiecewiseLinearX2(k,j,il,iu,w,bcc,wl,wr);
          }

          // compute and store fluxes at j=jstart
          RiemannSolver(k,j,il,iu,IVY,b2,wl,wr,jflx_j); 
      
          // store the surface flux for flux correction
          if(pmb->pmy_mesh->multilevel==true) {
            for(int n=0; n<NFLUID; ++n) {
              for(int i=is; i<=ie; ++i)
                pmb->pbval->surface_flux_[inner_x2](n,k,i)=jflx_j(n,i);
            }
          }

          // store electric fields, compute weights for GS07 CT algorithm at j=jstart
          if (MAGNETIC_FIELDS_ENABLED) {
#pragma simd
            for (int i=il; i<=iu; ++i){
              ei_x2f(X2E1,k,j,i) = -jflx_j(IBY,i); // flx(IBY)= (v2*b3 - v3*b2) = -EMFX
              ei_x2f(X2E3,k,j,i) =  jflx_j(IBZ,i); // flx(IBZ)= (v2*b1 - v1*b2) =  EMFZ
              const Real& dx = pmb->pcoord->CenterWidth2(k,j,i);
              Real v_over_c = (1024)*dt*jflx_j(IDN,i)/(dx*(wl(IDN,i) + wr(IDN,i)));
              Real tmp_min = std::min(0.5,v_over_c);
              w_x2f(k,j,i) = 0.5 + std::max(-0.5,tmp_min);
            }
          }

          first_time_through_loop = false;;
        }

        // reconstruct L/R states at j+1
        if (step == 1) {
          DonorCellX2(k,j+1,il,iu,w,bcc,wl,wr);
        } else {
          PiecewiseLinearX2(k,j+1,il,iu,w,bcc,wl,wr);
        }

        // compute fluxes at j+1
        RiemannSolver(k,j+1,il,iu,IVY,b2,wl,wr,flx); 

        // update conserved fluid variables
        pmb->pcoord->Face2Area(k,j  ,is,ie,area   );
        pmb->pcoord->Face2Area(k,j+1,is,ie,area_p1);
        pmb->pcoord->CellVolume(k,j,is,ie,vol);
        if(k>=ks && k<=ke) {
          for (int n=0; n<NFLUID; ++n){
#pragma simd
            for (int i=is; i<=ie; ++i){
              u(n,k,j,i) -= dt*(area_p1(i)*flx(n,i) - area(i)*jflx_j(n,i))/vol(i);
            }
          }

          // add coordinate (geometric) source terms
          pmb->pcoord->CoordSrcTermsX2(k,j,dt,jflx_j,flx,w,bcc,u);
          // add physical source terms for a point mass potential
          pmb->pfluid->pf_srcterms->PhysicalSourceTermsX2(k,j,dt,jflx_j,flx,w,u);
        }

        // store the surface flux for flux correction
        if(pmb->pmy_mesh->multilevel==true && j==je) {
          for(int n=0; n<NFLUID; ++n) {
            for(int i=is; i<=ie; ++i)
              pmb->pbval->surface_flux_[outer_x2](n,k,i)=flx(n,i);
          }
        }

        // store electric fields, compute weights for GS07 CT algorithm
        if (MAGNETIC_FIELDS_ENABLED) {
#pragma simd
          for (int i=il; i<=iu; ++i){
            ei_x2f(X2E1,k,j+1,i) = -flx(IBY,i); // flx(IBY) = (v2*b3 - v3*b2) = -EMFX
            ei_x2f(X2E3,k,j+1,i) =  flx(IBZ,i); // flx(IBZ) = (v2*b1 - v1*b2) =  EMFZ
            const Real& dx = pmb->pcoord->CenterWidth2(k,j,i);
            Real v_over_c = (1024)*dt*flx(IDN,i)/(dx*(wl(IDN,i) + wr(IDN,i)));
            Real tmp_min = std::min(0.5,v_over_c);
            w_x2f(k,j+1,i) = 0.5 + std::max(-0.5,tmp_min);
          }
        }

        // store fluxes for j=j-1 in next iteration
        if(j<je)
          jflx_j = flx; 

      }
    }
  }

//--------------------------------------------------------------------------------------
// k-direction 

  if (pmb->block_size.nx3 > 1) {
    // set the loop limits
    il=is, iu=ie, jl=js, ju=je;
    if (MAGNETIC_FIELDS_ENABLED)
      il=is-1, iu=ie+1, jl=js-1, ju=je+1;
    bool first_time_through_loop = true;
#pragma omp for schedule(static)
    for (int k=ks; k<=ke; ++k){

      if (first_time_through_loop) {
        for (int j=jl; j<=ju; ++j){

          // reconstruct L/R states at k=kstart
          if (step == 1) {
            DonorCellX3(k,j,il,iu,w,bcc,wl,wr);
          } else {
            PiecewiseLinearX3(k,j,il,iu,w,bcc,wl,wr);
          }

          // compute and store fluxes at k=kstart
          RiemannSolver(k,j,il,iu,IVZ,b3,wl,wr,flx);

          // store the surface flux for flux correction
          if(pmb->pmy_mesh->multilevel==true) {
            for(int n=0; n<NFLUID; ++n) {
              for(int i=is; i<=ie; ++i)
                pmb->pbval->surface_flux_[inner_x3](n,j,i)=flx(n,i);
            }
          }

          // store electric fields, compute weights for GS07 CT algorithm at k=kstart
          if (MAGNETIC_FIELDS_ENABLED) {
#pragma simd
            for (int i=il; i<=iu; ++i){
              ei_x3f(X3E2,k,j,i) = -flx(IBY,i); // flx(IBY) = (v3*b1 - v1*b3) = -EMFY
              ei_x3f(X3E1,k,j,i) =  flx(IBZ,i); // flx(IBZ) = (v3*b2 - v2*b3) =  EMFX
              const Real& dx = pmb->pcoord->CenterWidth3(k,j,i);
              Real v_over_c = (1024)*dt*flx(IDN,i)/(dx*(wl(IDN,i) + wr(IDN,i)));
              Real tmp_min = std::min(0.5,v_over_c);
              w_x3f(k,j,i) = 0.5 + std::max(-0.5,tmp_min);
            }
          }

          // store fluxes at k=kstart over all i,j
          for (int n=0; n<NFLUID; ++n){
#pragma simd
            for (int i=il; i<=iu; ++i){
              kflx_k(n,j,i) = flx(n,i);
            }
          }
        }
        first_time_through_loop = false;
      }

      for (int j=jl; j<=ju; ++j){

        // reconstruct L/R states at k+1
        if (step == 1) {
          DonorCellX3(k+1,j,il,iu,w,bcc,wl,wr);
        } else {
          PiecewiseLinearX3(k+1,j,il,iu,w,bcc,wl,wr);
        }

        // compute fluxes at k+1
        RiemannSolver(k+1,j,il,iu,IVZ,b3,wl,wr,flx);
        if(j>=js && j<=je) {
          // update conserved fluid variables
          pmb->pcoord->Face3Area(k  ,j,is,ie,area   );
          pmb->pcoord->Face3Area(k+1,j,is,ie,area_p1);
          pmb->pcoord->CellVolume(k,j,is,ie,vol);
          for (int n=0; n<NFLUID; ++n){
#pragma simd
            for (int i=is; i<=ie; ++i){
              u(n,k,j,i) -= dt*(area_p1(i)*flx(n,i) - area(i)*kflx_k(n,j,i))/vol(i);
            }
          }

          // add coordinate (geometric) source terms
          pmb->pcoord->CoordSrcTermsX3(k,j,dt,kflx_k,flx,w,bcc,u);
	  // add physical source terms for a point mass potential
	  pmb->pfluid->pf_srcterms->PhysicalSourceTermsX3(k,j,dt,kflx_k,flx,w,u);
        }

        // store the surface flux for flux correction
        if(pmb->pmy_mesh->multilevel==true && k==ke) {
          for(int n=0; n<NFLUID; ++n) {
            for(int i=is; i<=ie; ++i)
              pmb->pbval->surface_flux_[outer_x3](n,j,i)=flx(n,i);
          }
        }

        // store electric fields, compute weights for GS07 CT algorithm
        if (MAGNETIC_FIELDS_ENABLED) {
#pragma simd
          for (int i=il; i<=iu; ++i){
            ei_x3f(X3E2,k+1,j,i) = -flx(IBY,i); // flx(IBY) = (v3*b1 - v1*b3) = -EMFY
            ei_x3f(X3E1,k+1,j,i) =  flx(IBZ,i); // flx(IBZ) = (v3*b2 - v2*b3) =  EMFX
            const Real& dx = pmb->pcoord->CenterWidth3(k,j,i);
            Real v_over_c = (1024)*dt*flx(IDN,i)/(dx*(wl(IDN,i) + wr(IDN,i)));
            Real tmp_min = std::min(0.5,v_over_c);
            w_x3f(k+1,j,i) = 0.5 + std::max(-0.5,tmp_min);
          }
        }

        // store fluxes for k=k-1 in next iteration
        if(k<ke) {
          for (int n=0; n<NFLUID; ++n){
#pragma simd
            for (int i=is; i<=ie; ++i){
              kflx_k(n,j,i) = flx(n,i);
            }
          }
        }

      }
    }
  }

} // end of omp parallel region

//--------------------------------------------------------------------------------------
//  Add user source terms

  if (pmb->pfluid->pf_srcterms->UserSourceTerm != NULL)
    pmb->pfluid->pf_srcterms->UserSourceTerm(pmb->pmy_mesh->time,dt,w,u);

  return;
}
