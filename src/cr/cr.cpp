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
//! \file radiation.cpp
//  \brief implementation of functions in class Radiation
//======================================================================================


#include <sstream>  // msg
#include <iostream>  // cout
#include <stdexcept> // runtime erro
#include <stdio.h>  // fopen and fwrite


// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp" 
#include "cr.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../globals.hpp"
#include "../coordinates/coordinates.hpp"
#include "integrators/cr_integrators.hpp"

// constructor, initializes data structures and parameters

// The default opacity function.

// This function also needs to set the streaming velocity
// This is needed to calculate the work term 
inline void DefaultOpacity(MeshBlock *pmb, AthenaArray<Real> &u_cr, 
              AthenaArray<Real> &prim, AthenaArray<Real> &bcc)
{
  // set the default opacity to be a large value in the default hydro case
  CosmicRay *pcr=pmb->pcr;
  int kl=pmb->ks, ku=pmb->ke;
  int jl=pmb->js, ju=pmb->je;
  int il=pmb->is-1, iu=pmb->ie+1;
  if(pmb->block_size.nx2 > 1){
    jl -= 1;
    ju += 1;
  }
  if(pmb->block_size.nx3 > 1){
    kl -= 1;
    ku += 1;
  }


  Real invlim = 1.0/pcr->vmax;


  for(int k=kl; k<=ku; ++k){
    for(int j=jl; j<=ju; ++j){
#pragma omp simd
      for(int i=il; i<=iu; ++i){
        pcr->sigma_diff(0,k,j,i) = pcr->max_opacity;
        pcr->sigma_diff(1,k,j,i) = pcr->max_opacity;
        pcr->sigma_diff(2,k,j,i) = pcr->max_opacity;  

      }
    }
  }

  // Need to calculate the rotation matrix 
  // We need this to determine the direction of rotation velocity


  // The information stored in the array
  // b_angle is
  // b_angle[0]=sin_theta_b
  // b_angle[1]=cos_theta_b
  // b_angle[2]=sin_phi_b
  // b_angle[3]=cos_phi_b
 

  if(MAGNETIC_FIELDS_ENABLED && (pcr->stream_flag > 0)){
    //First, calculate B_dot_grad_Pc
    for(int k=kl; k<=ku; ++k){
      for(int j=jl; j<=ju; ++j){
      // diffusion coefficient is calculated with respect to B direction
      // Use a simple estimate of Grad Pc

    // x component
        pmb->pcoord->CenterWidth1(k,j,il-1,iu+1,pcr->cwidth);
        for(int i=il; i<=iu; ++i){
          Real distance = 0.5*(pcr->cwidth(i-1) + pcr->cwidth(i+1))
                         + pcr->cwidth(i);
          Real dprdx=(u_cr(CRE,k,j,i+1) - u_cr(CRE,k,j,i-1))/3.0;
          dprdx /= distance;
          pcr->sigma_adv(0,k,j,i) = dprdx;
        }
    // y component
        pmb->pcoord->CenterWidth2(k,j-1,il,iu,pcr->cwidth1);       
        pmb->pcoord->CenterWidth2(k,j,il,iu,pcr->cwidth);
        pmb->pcoord->CenterWidth2(k,j+1,il,iu,pcr->cwidth2);

        for(int i=il; i<=iu; ++i){
          Real distance = 0.5*(pcr->cwidth1(i) + pcr->cwidth2(i))
                         + pcr->cwidth(i);
          Real dprdy=(u_cr(CRE,k,j+1,i) - u_cr(CRE,k,j-1,i))/3.0;
          dprdy /= distance;
          pcr->sigma_adv(1,k,j,i) = dprdy;

        } 
// z component
        pmb->pcoord->CenterWidth3(k-1,j,il,iu,pcr->cwidth1);       
        pmb->pcoord->CenterWidth3(k,j,il,iu,pcr->cwidth);
        pmb->pcoord->CenterWidth3(k+1,j,il,iu,pcr->cwidth2);

        for(int i=il; i<=iu; ++i){
          Real distance = 0.5*(pcr->cwidth1(i) + pcr->cwidth2(i))
                          + pcr->cwidth(i);
          Real dprdz=(u_cr(CRE,k+1,j,i) -  u_cr(CRE,k-1,j,i))/3.0;
          dprdz /= distance;
          pcr->sigma_adv(2,k,j,i) = dprdz;
        }       


        for(int i=il; i<=iu; ++i){
          // Now calculate the angles of B
          Real bxby = sqrt(bcc(IB1,k,j,i)*bcc(IB1,k,j,i) +
                           bcc(IB2,k,j,i)*bcc(IB2,k,j,i));
          Real btot = sqrt(bcc(IB1,k,j,i)*bcc(IB1,k,j,i) +
                           bcc(IB2,k,j,i)*bcc(IB2,k,j,i) + 
                           bcc(IB3,k,j,i)*bcc(IB3,k,j,i));
          
          if(btot > TINY_NUMBER){
            pcr->b_angle(0,k,j,i) = bxby/btot;
            pcr->b_angle(1,k,j,i) = bcc(IB3,k,j,i)/btot;
          }else{
            pcr->b_angle(0,k,j,i) = 1.0;
            pcr->b_angle(1,k,j,i) = 0.0;
          }
          if(bxby > TINY_NUMBER){
            pcr->b_angle(2,k,j,i) = bcc(IB2,k,j,i)/bxby;
            pcr->b_angle(3,k,j,i) = bcc(IB1,k,j,i)/bxby;
          }else{
            pcr->b_angle(2,k,j,i) = 0.0;
            pcr->b_angle(3,k,j,i) = 1.0;            
          }

          Real va = sqrt(btot*btot/prim(IDN,k,j,i));
          if(va < TINY_NUMBER){
            pcr->sigma_adv(0,k,j,i) = pcr->max_opacity;
          }else{
            Real b_grad_pc = bcc(IB1,k,j,i) * pcr->sigma_adv(0,k,j,i)
                           + bcc(IB2,k,j,i) * pcr->sigma_adv(1,k,j,i)
                           + bcc(IB3,k,j,i) * pcr->sigma_adv(2,k,j,i);
            pcr->sigma_adv(0,k,j,i) = fabs(b_grad_pc)/(btot * va * (1.0 + 1.0/3.0) 
                                               * invlim * u_cr(CRE,k,j,i));
          }
          pcr->sigma_adv(1,k,j,i) = pcr->max_opacity;
          pcr->sigma_adv(2,k,j,i) = pcr->max_opacity;

        }//end i        

      }// end j
    }// end k

  }// End MHD  
  else{

    for(int k=kl; k<=ku; ++k){
      for(int j=jl; j<=ju; ++j){
#pragma omp simd
        for(int i=il; i<=iu; ++i){

          pcr->sigma_adv(0,k,j,i) = pcr->max_opacity;
          pcr->sigma_adv(1,k,j,i) = pcr->max_opacity;
          pcr->sigma_adv(2,k,j,i) = pcr->max_opacity;  

          pcr->v_adv(0,k,j,i) = 0.0;   
          pcr->v_adv(1,k,j,i) = 0.0;
          pcr->v_adv(2,k,j,i) = 0.0;
        }
      }
    }

  }// end MHD and stream flag

}

CosmicRay::CosmicRay(MeshBlock *pmb, ParameterInput *pin):
    pmy_block(pmb), u_cr(NCR,pmb->ncells3,pmb->ncells2,pmb->ncells1),
    u_cr1(NCR,pmb->ncells3,pmb->ncells2,pmb->ncells1),
    sigma_diff(3,pmb->ncells3,pmb->ncells2,pmb->ncells1),
    sigma_adv(3,pmb->ncells3,pmb->ncells2,pmb->ncells1),
    v_adv(3,pmb->ncells3,pmb->ncells2,pmb->ncells1),
    v_diff(3,pmb->ncells3,pmb->ncells2,pmb->ncells1),
// constructor overload resolution of non-aggregate class type AthenaArray<Real>
    flux{ {NCR, pmb->ncells3, pmb->ncells2, pmb->ncells1+1},
      {NCR,pmb->ncells3, pmb->ncells2+1, pmb->ncells1,
       (pmb->pmy_mesh->f2 ? AthenaArray<Real>::DataStatus::allocated :
        AthenaArray<Real>::DataStatus::empty)},
      {NCR,pmb->ncells3+1, pmb->ncells2, pmb->ncells1,
       (pmb->pmy_mesh->f3 ? AthenaArray<Real>::DataStatus::allocated :
        AthenaArray<Real>::DataStatus::empty)}}, UserSourceTerm_{},
    coarse_cr_(NCR,pmb->ncc3, pmb->ncc2, pmb->ncc1,
             (pmb->pmy_mesh->multilevel ? AthenaArray<Real>::DataStatus::allocated :
              AthenaArray<Real>::DataStatus::empty)),
    cr_bvar(pmb, &u_cr, &coarse_cr_, flux, true){

  Mesh *pm = pmy_block->pmy_mesh;

  pmb->RegisterMeshBlockData(u_cr);
  // "Enroll" in S/AMR by adding to vector of tuples of pointers in MeshRefinement class
  if (pm->multilevel) {
    refinement_idx = pmy_block->pmr->AddToRefinement(&u_cr, &coarse_cr_);
  }

  cr_source_defined = false;

  cr_bvar.bvar_index = pmb->pbval->bvars.size();
  pmb->pbval->bvars.push_back(&cr_bvar);
  pmb->pbval->bvars_main_int.push_back(&cr_bvar);      

  vmax = pin->GetOrAddReal("cr","vmax",1.0);
  vlim = pin->GetOrAddReal("cr","vlim",0.9);
  max_opacity = pin->GetOrAddReal("cr","max_opacity",1.e10);
  stream_flag = pin->GetOrAddInteger("cr","vs_flag",1);  
  src_flag = pin->GetOrAddInteger("cr","src_flag",1);
  

  int nc1 = pmb->ncells1, nc2 = pmb->ncells2, nc3 = pmb->ncells3;

  b_grad_pc.NewAthenaArray(nc3,nc2,nc1);
  b_angle.NewAthenaArray(4,nc3,nc2,nc1);
  

  cwidth.NewAthenaArray(nc1);
  cwidth1.NewAthenaArray(nc1);
  cwidth2.NewAthenaArray(nc1);
  
  // set a default opacity function
  UpdateOpacity = DefaultOpacity;

  pcrintegrator = new CRIntegrator(this, pin);

}



//Enrol the function to update opacity

void CosmicRay::EnrollOpacityFunction(CROpacityFunc MyOpacityFunction)
{
  UpdateOpacity = MyOpacityFunction;
  
}

void CosmicRay::EnrollUserCRSource(CRSrcTermFunc my_func)
{
  UserSourceTerm_ = my_func;
  cr_source_defined = true;
  
}

