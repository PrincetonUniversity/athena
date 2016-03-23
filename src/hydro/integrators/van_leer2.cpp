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
#include "../hydro.hpp"
#include "../../field/field.hpp"
#include "../../mesh.hpp"
#include "../srcterms/srcterms.hpp"
#include "../../bvals/bvals.hpp"
#include "../viscosity/viscosity.hpp"

// this class header
#include "hydro_integrator.hpp"

// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//--------------------------------------------------------------------------------------
//! \fn  void HydroIntegrator::CalculateFluxes
//  \brief Calculate Hydrodynamic Fluxes using the Riemann solver

void HydroIntegrator::CalculateFluxes(MeshBlock *pmb,AthenaArray<Real> &u,
  AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &bcc, const int step)
{
  AthenaArray<Real> &x1flux=pmb->phydro->flux[x1face];
  AthenaArray<Real> &x2flux=pmb->phydro->flux[x2face];
  AthenaArray<Real> &x3flux=pmb->phydro->flux[x3face];
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int il, iu, jl, ju, kl, ku;
  Real dt;
  if (step == 1) {
    dt = 0.5*(pmb->pmy_mesh->dt);
  } else {
    dt = (pmb->pmy_mesh->dt);
  }

  AthenaArray<Real> b1,b2,b3,ei_x1f,ei_x2f,ei_x3f,w_x1f,w_x2f,w_x3f;
  if (MAGNETIC_FIELDS_ENABLED) {
    b1.InitWithShallowCopy(b.x1f);
    b2.InitWithShallowCopy(b.x2f);
    b3.InitWithShallowCopy(b.x3f);
    ei_x1f.InitWithShallowCopy(pmb->pfield->ei.x1f);
    ei_x2f.InitWithShallowCopy(pmb->pfield->ei.x2f);
    ei_x3f.InitWithShallowCopy(pmb->pfield->ei.x3f);
    w_x1f.InitWithShallowCopy(pmb->pfield->wght.x1f);
    w_x2f.InitWithShallowCopy(pmb->pfield->wght.x2f);
    w_x3f.InitWithShallowCopy(pmb->pfield->wght.x3f);
  }

  int tid=0;
  int nthreads = pmb->pmy_mesh->GetNumMeshThreads();
#pragma omp parallel default(shared) private(tid) num_threads(nthreads)
{
#ifdef OPENMP_PARALLEL
  tid=omp_get_thread_num();
#endif

  AthenaArray<Real> wl, wr, flx;
  wl.InitWithShallowSlice(wl_,3,tid,1);
  wr.InitWithShallowSlice(wr_,3,tid,1);
  flx.InitWithShallowSlice(flx_,3,tid,1);

//----------Viscosity update
  if(VISCOSITY) pmb->phydro->pf_viscosity->ViscosityTerms(dt,w,u);

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

      // store fluxes
      if(k>=ks && k<=ke && j>=js && j<=je) {
        for(int n=0; n<NHYDRO; n++) {
#pragma simd
          for(int i=is; i<=ie+1; i++)
            x1flux(n,k,j,i)=flx(n,i);
        }
      }

      // store electric fields, compute weights for GS07 CT algorithm
      // no correction to the EMFs is required; they are corrected later
      if (MAGNETIC_FIELDS_ENABLED) {
#pragma simd
        for (int i=is; i<=ie+1; ++i){
          ei_x1f(X1E3,k,j,i) = -flx(IBY,i); // flux(IBY) = (v1*b2 - v2*b1) = -EMFZ
          ei_x1f(X1E2,k,j,i) =  flx(IBZ,i); // flux(IBZ) = (v1*b3 - v3*b1) =  EMFY
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
#pragma omp for schedule(static)
      for (int j=js; j<=je+1; ++j){

        // reconstruct L/R states at j
        if (step == 1) {
          DonorCellX2(k,j,il,iu,w,bcc,wl,wr);
        } else {
          PiecewiseLinearX2(k,j,il,iu,w,bcc,wl,wr);
        }

        // compute fluxes at j
        RiemannSolver(k,j,il,iu,IVY,b2,wl,wr,flx); 

        // store fluxes
        if(k>=ks && k<=ke) {
          for(int n=0; n<NHYDRO; n++) {
#pragma simd
            for(int i=is; i<=ie; i++)
              x2flux(n,k,j,i)=flx(n,i);
          }
        }

        // store electric fields, compute weights for GS07 CT algorithm
        // no correction to the EMFs is required; they are corrected later
        if (MAGNETIC_FIELDS_ENABLED) {
#pragma simd
          for (int i=il; i<=iu; ++i){
            ei_x2f(X2E1,k,j,i) = -flx(IBY,i); // flx(IBY) = (v2*b3 - v3*b2) = -EMFX
            ei_x2f(X2E3,k,j,i) =  flx(IBZ,i); // flx(IBZ) = (v2*b1 - v1*b2) =  EMFZ
            const Real& dx = pmb->pcoord->CenterWidth2(k,j,i);
            Real v_over_c = (1024)*dt*flx(IDN,i)/(dx*(wl(IDN,i) + wr(IDN,i)));
            Real tmp_min = std::min(0.5,v_over_c);
            w_x2f(k,j,i) = 0.5 + std::max(-0.5,tmp_min);
          }
        }
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
#pragma omp for schedule(static)
    for (int k=ks; k<=ke+1; ++k){
      for (int j=jl; j<=ju; ++j){

        // reconstruct L/R states at k
        if (step == 1) {
          DonorCellX3(k,j,il,iu,w,bcc,wl,wr);
        } else {
          PiecewiseLinearX3(k,j,il,iu,w,bcc,wl,wr);
        }

        // compute fluxes at k
        RiemannSolver(k,j,il,iu,IVZ,b3,wl,wr,flx);

        if(j>=js && j<=je) {
          for(int n=0; n<NHYDRO; n++) {
#pragma simd
            for(int i=is; i<=ie; i++)
              x3flux(n,k,j,i)=flx(n,i);
          }
        }

        // store electric fields, compute weights for GS07 CT algorithm
        // no correction to the EMFs is required; they are corrected later
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
      }
    }
  }

} // end of omp parallel region

  return;
}


//--------------------------------------------------------------------------------------
//! \fn  void HydroIntegrator::FluxDivergence
//  \brief Integrate the conservative variables using the calculated fluxes

void HydroIntegrator::FluxDivergence(MeshBlock *pmb,AthenaArray<Real> &u,
  AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &bcc, const int step)
{
  AthenaArray<Real> &x1flux=pmb->phydro->flux[x1face];
  AthenaArray<Real> &x2flux=pmb->phydro->flux[x2face];
  AthenaArray<Real> &x3flux=pmb->phydro->flux[x3face];
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
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
  AthenaArray<Real> x1area, x2area, x2area_p1, x3area, x3area_p1, vol;
  x1area.InitWithShallowSlice(x1face_area_,2,tid,1);
  x2area.InitWithShallowSlice(x2face_area_,2,tid,1);
  x2area_p1.InitWithShallowSlice(x2face_area_p1_,2,tid,1);
  x3area.InitWithShallowSlice(x3face_area_,2,tid,1);
  x3area_p1.InitWithShallowSlice(x3face_area_p1_,2,tid,1);
  vol.InitWithShallowSlice(cell_volume_,2,tid,1);

  // update conserved hydro variables
  if (pmb->block_size.nx3 > 1) {
#pragma omp for schedule(static)
    for (int k=ks; k<=ke; ++k) { 
      for (int j=js; j<=je; ++j) {
        pmb->pcoord->CellVolume(k,j,is,ie,vol);
        pmb->pcoord->Face1Area(k,j,is,ie+1,x1area);
        pmb->pcoord->Face2Area(k,j  ,is,ie,x2area   );
        pmb->pcoord->Face2Area(k,j+1,is,ie,x2area_p1);
        pmb->pcoord->Face3Area(k  ,j,is,ie,x3area   );
        pmb->pcoord->Face3Area(k+1,j,is,ie,x3area_p1);
        for (int n=0; n<NHYDRO; ++n) {
//#pragma simd
          for (int i=is; i<=ie; ++i) {
            u(n,k,j,i) -= dt*(x1area(i+1) *x1flux(n,k,j,i+1) 
                            - x1area(i)   *x1flux(n,k,j,i)
                            + x2area_p1(i)*x2flux(n,k,j+1,i)
                            - x2area(i)   *x2flux(n,k,j,i)
                            + x3area_p1(i)*x3flux(n,k+1,j,i)
                            - x3area(i)   *x3flux(n,k,j,i))/vol(i);
          }
        }
      }
    }
  }
  else if (pmb->block_size.nx2 > 1) {
    int k=pmb->ks;
#pragma omp for schedule(static)
    for (int j=js; j<=je; ++j) {
      pmb->pcoord->CellVolume(k,j,is,ie,vol);
      pmb->pcoord->Face1Area(k,j,is,ie+1,x1area);
      pmb->pcoord->Face2Area(k,j  ,is,ie,x2area   );
      pmb->pcoord->Face2Area(k,j+1,is,ie,x2area_p1);
      for (int n=0; n<NHYDRO; ++n) {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          u(n,k,j,i) -= dt*(x1area(i+1) *x1flux(n,k,j,i+1) 
                          - x1area(i)   *x1flux(n,k,j,i)
                          + x2area_p1(i)*x2flux(n,k,j+1,i)
                          - x2area(i)   *x2flux(n,k,j,i))/vol(i);
        }
      }
    }
  }
  else {
    int j=pmb->js, k=pmb->ks;
    pmb->pcoord->CellVolume(k,j,is,ie,vol);
    pmb->pcoord->Face1Area(k,j,is,ie+1,x1area);
#pragma omp for schedule(static)
    for (int n=0; n<NHYDRO; ++n) {
#pragma simd
      for (int i=is; i<=ie; ++i) {
        u(n,k,j,i) -= dt*(x1area(i+1)*x1flux(n,k,j,i+1) 
                        - x1area(i)*  x1flux(n,k,j,i))/vol(i);
      }
    }
  }

} // end of omp parallel region


  // add coordinate (geometric) source terms
  pmb->pcoord->CoordSrcTerms(dt,pmb->phydro->flux,w,bcc,u);
  // add physical source terms for a point mass potential
  pmb->phydro->pf_srcterms->PhysicalSourceTerms(dt,pmb->phydro->flux,w,u);

//--------------------------------------------------------------------------------------
//  Add user source terms

  if (pmb->phydro->pf_srcterms->UserSourceTerm != NULL)
    pmb->phydro->pf_srcterms->UserSourceTerm(pmb, pmb->pmy_mesh->time,dt,w,bcc,u);

  return;
}
