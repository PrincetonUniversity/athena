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
//! \file rad_transport.cpp
//  \brief implementation of radiation integrators
//======================================================================================


// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../parameter_input.hpp"
#include "../../mesh/mesh.hpp"
#include "../tc.hpp"
#include "../../coordinates/coordinates.hpp" //
#include "../../reconstruct/reconstruction.hpp"
#include "../../eos/eos.hpp"
#include "../../utils/utils.hpp"
#include <algorithm>   // min,max

// class header
#include "tc_integrators.hpp"

// MPI/OpenMP header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif


// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif


void TCIntegrator::CalculateFluxes(AthenaArray<Real> &w,
            AthenaArray<Real> &bcc, AthenaArray<Real> &tc, const int order)
{
  ThermalConduction *ptc=pmy_tc;
  MeshBlock *pmb=ptc->pmy_block;
  Coordinates *pco = pmb->pcoord;

  Real gamma_1=pmb->peos->GetGamma()-1.0;

  
  Real invlim = 1.0/ptc->vmax;
  
  int ncells1 = pmb->ncells1, ncells2 = pmb->ncells2, 
  ncells3 = pmb->ncells3; 

  AthenaArray<Real> &x1flux=ptc->flux[X1DIR];

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int il, iu, jl, ju, kl, ku;
  jl = js, ju=je, kl=ks, ku=ke;

  if(ncells2 > 1)
  {
    if(ncells3 == 1){
      jl=js-1, ju=je+1, kl=ks, ku=ke;
    }else{
      jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
    }
 
  }



//--------------------------------------------------------------------------------------
  // First, calculate the diffusion velocity along three coordinate system
  for (int k=0; k<ncells3; ++k){
    for(int j=0; j<ncells2; ++j){

      // diffusion velocity along the direction of sigma vector
      // We first assume B is along x coordinate
      // Then rotate according to B direction to the actual acooridnate

      for(int i=0; i<ncells1; ++i){
        Real kappa = w(IDN,k,j,i)/(gamma_1 * ptc->kappa(0,k,j,i));

        Real taux = taufact_ * ptc->vmax * kappa * pco->dx1f(i);
        taux = taux * taux/2.0;
        Real diffv = sqrt((1.0 - exp(-taux)) / taux);

        if(taux < 1.e-3)
          diffv = sqrt((1.0 - 0.5* taux));

        vdiff_(0,k,j,i) = ptc->vmax * diffv;
      }

      // y direction
      if(ncells2 > 1){
        pco->CenterWidth2(k,j,0,ncells1-1,cwidth2_);
          // get the optical depth across the cell
        for(int i=0; i<ncells1; ++i){
          Real kappa = w(IDN,k,j,i)/(gamma_1 * ptc->kappa(1,k,j,i));
          Real tauy = taufact_ * ptc->vmax * kappa * cwidth2_(i);          
          tauy = tauy * tauy/2.0;
          Real diffv = sqrt((1.0 - exp(-tauy)) / tauy);

          if(tauy < 1.e-3)
            diffv = sqrt((1.0 - 0.5* tauy));

          vdiff_(1,k,j,i) = ptc->vmax * diffv;            
        }// end i
      }else{
        for(int i=0; i<ncells1; ++i){
          vdiff_(1,k,j,i) = 0.0;
        }
      }
      // z direction
      if(ncells3 > 1){
        pco->CenterWidth3(k,j,0,ncells1-1,cwidth3_);
          // get the optical depth across the cell
        for(int i=0; i<ncells1; ++i){
          Real kappa = w(IDN,k,j,i)/(gamma_1 * ptc->kappa(2,k,j,i));
          Real tauz = taufact_ * ptc->vmax * kappa * cwidth3_(i);  
          tauz = tauz * tauz/2.0;
          Real diffv = sqrt((1.0 - exp(-tauz)) / tauz);

          if(tauz < 1.e-3)
            diffv = sqrt((1.0 - 0.5* tauz));

          vdiff_(2,k,j,i) = ptc->vmax * diffv;            
        }
      }else{
        for(int i=0; i<ncells1; ++i)
          vdiff_(2,k,j,i) = 0.0;
      }

      //rotate the v_diff vector to the local coordinate
      if(MAGNETIC_FIELDS_ENABLED){
        for(int i=0; i<ncells1; ++i){
          

          InvRotateVec(ptc->b_angle(0,k,j,i),ptc->b_angle(1,k,j,i),
                      ptc->b_angle(2,k,j,i),ptc->b_angle(3,k,j,i), 
                    vdiff_(0,k,j,i),vdiff_(1,k,j,i),vdiff_(2,k,j,i));
          // take the absolute value
          // Also add the Alfven velocity for the streaming flux
          vdiff_(0,k,j,i) = fabs(vdiff_(0,k,j,i));

          vdiff_(1,k,j,i) = fabs(vdiff_(1,k,j,i));
                                 
          vdiff_(2,k,j,i) = fabs(vdiff_(2,k,j,i));

        }

      }// end if MHD

    }//end j
  }//end k

  // prepare Array for reconstruction
  for(int n=0; n<NTC; ++n){
    for (int k=0; k<ncells3; ++k){
      for(int j=0; j<ncells2; ++j){
        for(int i=0; i<ncells1; ++i){
           utc_rho_t_(n,k,j,i) = tc(n,k,j,i);
        }// end i
      }// end j
    }// end k
  }// end n
  // add rho
  for (int k=0; k<ncells3; ++k){
    for(int j=0; j<ncells2; ++j){
      for(int i=0; i<ncells1; ++i){
         utc_rho_t_(TCRHO,k,j,i) = w(IDN,k,j,i)/gamma_1;
      }// end i
    }// end j
  }// end k
  // add Tgas
  for (int k=0; k<ncells3; ++k){
    for(int j=0; j<ncells2; ++j){
      for(int i=0; i<ncells1; ++i){
         utc_rho_t_(TCT,k,j,i) = w(IPR,k,j,i)/w(IDN,k,j,i);
      }// end i
    }// end j
  }// end k

//--------------------------------------------------------------------------------------
// i-direction


  for (int k=kl; k<=ku; ++k){
    for (int j=jl; j<=ju; ++j){
      // First, need to do reconstruction
      // to reconstruct Etc, Ftc, rho, T and 
      if (order == 1) {
        pmb->precon->DonorCellX1(k, j, is-1, ie+1, utc_rho_t_, utc_l_, utc_r_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX1(k, j, is-1, ie+1, utc_rho_t_, utc_l_, utc_r_);
      } else {
        pmb->precon->PiecewiseParabolicX1(k, j, is-1, ie+1, utc_rho_t_, utc_l_, utc_r_);
      }

      // get the optical depth across the cell
#pragma omp simd
      for(int i=is; i<=ie+1; ++i){
        vdiff_l_(i) = -0.5*(vdiff_(0,k,j,i-1)+vdiff_(0,k,j,i));
      }
#pragma omp simd
      for(int i=is; i<=ie+1; ++i){
        vdiff_r_(i) = -vdiff_l_(i);
      }      

      // calculate the flux
      TCFlux(TCF1, is, ie+1, utc_l_, utc_r_, vdiff_l_, vdiff_r_, dflx_);
      // store the flux
      for(int n=0; n<NCR; ++n){
#pragma omp simd
        for(int i=is; i<=ie+1; ++i){
          x1flux(n,k,j,i) = dflx_(n,i);
        }
      }

    }
  }



//--------------------------------------------------------------------------------------
// j-direction
  if(pmb->pmy_mesh->f2){

    AthenaArray<Real> &x2flux=ptc->flux[X2DIR];

    il=is-1; iu=ie+1; kl=ks; ku=ke;
    if (ncells3 ==  1) // 2D
      kl = ks, ku = ke;
    else // 3D
      kl = ks-1, ku = ke+1;    
    for (int k=kl; k<=ku; ++k){
      //reconstruct the first row
      if (order == 1) {
        pmb->precon->DonorCellX2(k, js-1, il, iu, utc_rho_t_, utc_l_, utc_r_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX2(k, js-1, il, iu, utc_rho_t_, utc_l_, utc_r_);
      } else {
        pmb->precon->PiecewiseParabolicX2(k, js-1, il, iu, utc_rho_t_, utc_l_, utc_r_);
      }


      for (int j=js; j<=je+1; ++j){
        if (order == 1) {
          pmb->precon->DonorCellX2(k, j, il, iu, utc_rho_t_, utc_lb_, utc_r_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX2(k, j, il, iu, utc_rho_t_, utc_lb_, utc_r_);
        } else {
          pmb->precon->PiecewiseParabolicX2(k, j, il, iu, utc_rho_t_, utc_lb_, utc_r_);
        }

        // get the optical depth across the cell
#pragma omp simd
        for(int i=il; i<=iu; ++i){
          vdiff_l_(i) = -0.5*(vdiff_(1,k,j-1,i)+vdiff_(1,k,j,i));
        }
#pragma omp simd
        for(int i=il; i<=iu; ++i){
          vdiff_r_(i) = -vdiff_l_(i);
        }
      // calculate the flux
        TCFlux(TCF2, il, iu, utc_l_, utc_r_, vdiff_l_, vdiff_r_, dflx_);        
      // store the flux
        for(int n=0; n<NTC; ++n){
#pragma omp simd
          for(int i=il; i<=iu; ++i){
            x2flux(n,k,j,i) = dflx_(n,i);
          }
        }
       // swap the array for next cycle
        utc_l_.SwapAthenaArray(utc_lb_);

      }// end j from js to je+1
    }
  }// finish j direction


//  k-direction
  if(pmb->pmy_mesh->f3){
    AthenaArray<Real> &x3flux=ptc->flux[X3DIR];
    il =is-1, iu=ie+1, jl=js-1, ju=je+1;

    for(int j=jl; j<=ju; ++j){

      if (order == 1) {
        pmb->precon->DonorCellX3(ks-1, j, il, iu, utc_rho_t_, utc_l_, utc_r_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX3(ks-1, j, il, iu, utc_rho_t_, utc_l_, utc_r_);
      } else {
        pmb->precon->PiecewiseParabolicX3(ks-1, j, il, iu, utc_rho_t_, utc_l_, utc_r_);
      }

      for (int k=ks; k<=ke+1; ++k) {
        // reconstruct L/R states at k
        if (order == 1) {
          pmb->precon->DonorCellX3(k, j, il, iu, utc_rho_t_, utc_lb_, utc_r_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX3(k, j, il, iu, utc_rho_t_, utc_lb_, utc_r_);
        } else {
          pmb->precon->PiecewiseParabolicX3(k, j, il, iu, utc_rho_t_, utc_lb_, utc_r_);
        }

#pragma omp simd
        for(int i=il; i<=iu; ++i){
          vdiff_l_(i) = -0.5*(vdiff_(2,k-1,j,i)+vdiff_(2,k,j,i));
        }
#pragma omp simd
        for(int i=il; i<=iu; ++i){
          vdiff_r_(i) = -vdiff_l_(i);
        }
      // calculate the flux
        TCFlux(TCF3, il, iu, utc_l_, utc_r_, vdiff_l_, vdiff_r_, dflx_);   
        for(int n=0; n<NTC; ++n){
#pragma omp simd
          for(int i=il; i<=iu; ++i){
            x3flux(n,k,j,i) = dflx_(n,i);
          }
        }

       // swap the array for next cycle
        utc_l_.SwapAthenaArray(utc_lb_);                

      }// end k loop
    }// end j loop 

  }// finish k direction

  //-----------------------------------------------------------------------
  // calculate coordinate source terms for Thermal conduction
  pco->AddCoordTermsDivergence(2,utc_rho_t_,coord_source_);
  
}


void TCIntegrator::FluxDivergence(const Real wght, AthenaArray<Real> &w, 
                                               AthenaArray<Real> &tc_out)
{
  ThermalConduction *ptc=pmy_tc;
  MeshBlock *pmb = ptc->pmy_block;
  Coordinates *pco = pmb->pcoord;
  Real gamma_1=pmb->peos->GetGamma()-1.0;

  AthenaArray<Real> &x1flux=ptc->flux[X1DIR];
  AthenaArray<Real> &x2flux=ptc->flux[X2DIR];
  AthenaArray<Real> &x3flux=ptc->flux[X3DIR];
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

 
  AthenaArray<Real> &x1area = x1face_area_, &x2area = x2face_area_,
                 &x2area_p1 = x2face_area_p1_, &x3area = x3face_area_,
                 &x3area_p1 = x3face_area_p1_, &vol = cell_volume_, &dflx = dflx_;

  for (int k=ks; k<=ke; ++k) { 
    for (int j=js; j<=je; ++j) {

      // calculate x1-flux divergence 
      pmb->pcoord->Face1Area(k,j,is,ie+1,x1area);

      for(int n=0; n<NTC; ++n){
 #pragma omp simd
        for(int i=is; i<=ie; ++i){
          dflx(n,i) = (x1area(i+1) *x1flux(n,k,j,i+1) - x1area(i)*x1flux(n,k,j,i));
        }// end n
      }// End i

     // calculate x2-flux
      if (pmb->block_size.nx2 > 1) {
        pmb->pcoord->Face2Area(k,j  ,is,ie,x2area   );
        pmb->pcoord->Face2Area(k,j+1,is,ie,x2area_p1);
        for(int n=0; n<NTC; ++n){
#pragma omp simd
          for(int i=is; i<=ie; ++i){
            dflx(n,i) += (x2area_p1(i)*x2flux(n,k,j+1,i) - x2area(i)*x2flux(n,k,j,i));
          }
        }
      }// end nx2

      // calculate x3-flux divergence
      if (pmb->block_size.nx3 > 1) {
        pmb->pcoord->Face3Area(k  ,j,is,ie,x3area   );
        pmb->pcoord->Face3Area(k+1,j,is,ie,x3area_p1);

        for(int n=0; n<NTC; ++n){
#pragma omp simd
          for(int i=is; i<=ie; ++i){
            dflx(n,i) += (x3area_p1(i)*x3flux(n,k+1,j,i) - x3area(i)*x3flux(n,k,j,i));
          }
        }
      }// end nx3
      // update variable with flux divergence
      // only need to update TCF1, TCF2 and TCF3
      // The evolved quantities are (T/e)F1, (T/e)F2, (T/e)F3
      // No need to apply flux divergence to u_tc(0)
      pmb->pcoord->CellVolume(k,j,is,ie,vol);
      for(int n=1; n<NTC; ++n){
#pragma omp simd
        for(int i=is; i<=ie; ++i){
          tc_out(n,k,j,i) -= w(IDN,k,j,i)*wght*dflx(n,i)/(gamma_1*vol(i));
        }
      }
      //get the energy source term
#pragma omp simd
      for(int i=is; i<=ie; ++i){
        tc_esource_(k,j,i) = -wght*(dflx(0,i)/vol(i)
                             + coord_source_(0,k,j,i));                     
      }     
    }// end j
  }// End k

  // Add coordinate source term
  for(int n=1; n<NTC; ++n)
    for(int k=ks; k<=ke; ++k)
      for(int j=js; j<=je; ++j)
#pragma omp simd
        for(int i=is; i<=ie; ++i){
          tc_out(n,k,j,i) += wght * coord_source_(n,k,j,i);
        }


}

