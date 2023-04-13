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
#include "../eos/eos.hpp"
#include "tc.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../globals.hpp"
#include "../coordinates/coordinates.hpp"
#include "integrators/tc_integrators.hpp"

// constructor, initializes data structures and parameters



//Default function to set the conduction coefficient
inline void DefaultOpacity(MeshBlock *pmb,  
              AthenaArray<Real> &prim, AthenaArray<Real> &bcc)
{
  // set the default opacity to be a large value in the default hydro case
  ThermalConduction *ptc=pmb->ptc;
  int nc1 = pmb->ncells1, nc2 = pmb->ncells2, nc3 = pmb->ncells3;


  for(int k=0; k<nc3; ++k){
    for(int j=0; j<nc2; ++j){
#pragma omp simd
      for(int i=0; i<nc1; ++i){

        ptc->kappa(0,k,j,i) = ptc->min_kappa;
        ptc->kappa(1,k,j,i) = ptc->min_kappa;
        ptc->kappa(2,k,j,i) = ptc->min_kappa;  

      }
    }
  }


  if(MAGNETIC_FIELDS_ENABLED){
    //First, calculate B_dot_grad_Pc
    for(int k=0; k<nc3; ++k){
      for(int j=0; j<nc2; ++j){
      // diffusion coefficient is calculated with respect to B direction
      // Use a simple estimate of Grad Pc
        for(int i=0; i<nc1; ++i){
          // Now calculate the angles of B
          Real bxby = sqrt(bcc(IB1,k,j,i)*bcc(IB1,k,j,i) +
                           bcc(IB2,k,j,i)*bcc(IB2,k,j,i));
          Real btot = sqrt(bcc(IB1,k,j,i)*bcc(IB1,k,j,i) +
                           bcc(IB2,k,j,i)*bcc(IB2,k,j,i) + 
                           bcc(IB3,k,j,i)*bcc(IB3,k,j,i));
          
          if(btot > TINY_NUMBER){
            ptc->b_angle(0,k,j,i) = bxby/btot;
            ptc->b_angle(1,k,j,i) = bcc(IB3,k,j,i)/btot;
          }else{
            ptc->b_angle(0,k,j,i) = 1.0;
            ptc->b_angle(1,k,j,i) = 0.0;
          }
          if(bxby > TINY_NUMBER){
            ptc->b_angle(2,k,j,i) = bcc(IB2,k,j,i)/bxby;
            ptc->b_angle(3,k,j,i) = bcc(IB1,k,j,i)/bxby;
          }else{
            ptc->b_angle(2,k,j,i) = 0.0;
            ptc->b_angle(3,k,j,i) = 1.0;            
          }


        }//end i        

      }// end j
    }// end k
  }// end MHD


}

ThermalConduction::ThermalConduction(MeshBlock *pmb, ParameterInput *pin):
    pmy_block(pmb), u_tc(NTC,pmb->ncells3,pmb->ncells2,pmb->ncells1),
    u_tc1(NCR,pmb->ncells3,pmb->ncells2,pmb->ncells1),
    kappa(3,pmb->ncells3,pmb->ncells2,pmb->ncells1),
    b_angle(4,pmb->ncells3,pmb->ncells2,pmb->ncells1),
// constructor overload resolution of non-aggregate class type AthenaArray<Real>
    flux{ {NTC, pmb->ncells3, pmb->ncells2, pmb->ncells1+1},
      {NTC,pmb->ncells3, pmb->ncells2+1, pmb->ncells1,
       (pmb->pmy_mesh->f2 ? AthenaArray<Real>::DataStatus::allocated :
        AthenaArray<Real>::DataStatus::empty)},
      {NTC,pmb->ncells3+1, pmb->ncells2, pmb->ncells1,
       (pmb->pmy_mesh->f3 ? AthenaArray<Real>::DataStatus::allocated :
        AthenaArray<Real>::DataStatus::empty)}},
    coarse_tc_(NTC,pmb->ncc3, pmb->ncc2, pmb->ncc1,
             (pmb->pmy_mesh->multilevel ? AthenaArray<Real>::DataStatus::allocated :
              AthenaArray<Real>::DataStatus::empty)),
    tc_bvar(pmb, &u_tc, &coarse_tc_, flux, true){

    Mesh *pm=pmy_block->pmy_mesh;

    pmb->RegisterMeshBlockData(u_tc);

  // "Enroll" in S/AMR by adding to vector of tuples of pointers in MeshRefinement class
    if (pm->multilevel) {
      refinement_idx = pmy_block->pmr->AddToRefinement(&u_tc, &coarse_tc_);
    }

    tc_bvar.bvar_index = pmb->pbval->bvars.size();
    pmb->pbval->bvars.push_back(&tc_bvar);
    pmb->pbval->bvars_main_int.push_back(&tc_bvar);      

    vmax = pin->GetOrAddReal("tc","vmax",1.0);

    min_kappa = pin->GetOrAddReal("tc","min_opacity",1.e-10);

    int nc1 = pmb->ncells1, nc2 = pmb->ncells2, nc3 = pmb->ncells3;
    
    // set a default opacity function
    UpdateOpacity = DefaultOpacity;



    ptcintegrator = new TCIntegrator(this, pin);

}



//Enrol the function to update opacity

void ThermalConduction::EnrollOpacityFunction(TCOpacityFunc MyOpacityFunction)
{
  UpdateOpacity = MyOpacityFunction;
  
}


//set the U_tc[0] to be the gas internal energy
// Also calculate the gas temperature and density

void ThermalConduction::Initialize(MeshBlock *pmb, AthenaArray<Real> &prim, 
	                                               AthenaArray<Real> &u_tc)
{
  // This is gamma - 1
  Real gamma_1=pmb->peos->GetGamma()-1.0;


  int nc1 = pmb->ncells1, nc2 = pmb->ncells2, nc3 = pmb->ncells3;


  for (int k=0; k<nc3; ++k){
    for(int j=0; j<nc2; ++j){
      for(int i=0; i<nc1; ++i){
       //set the first conserved variable to be the gas internal energy
        u_tc(0,k,j,i) = prim(IPR,k,j,i)/gamma_1;

      }
    }
  }
 
  
}


