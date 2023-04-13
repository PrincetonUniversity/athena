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
//! \file rad_source.cpp
//  \brief Add radiation source terms to both radiation and gas
//======================================================================================


// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../parameter_input.hpp"
#include "../../mesh/mesh.hpp"
#include "../tc.hpp"
#include "../../hydro/hydro.hpp"
#include "../../field/field.hpp"
#include "../../utils/utils.hpp"
#include "../../eos/eos.hpp"

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

//add the source terms implicitly



void TCIntegrator::AddSourceTerms(MeshBlock *pmb, const Real dt, 
  AthenaArray<Real> &u, AthenaArray<Real> &w, AthenaArray<Real> &u_tc)
{
  ThermalConduction *ptc=pmy_tc;

  Real gamma_1=pmb->peos->GetGamma()-1.0;


  Real vlim = ptc->vmax;
// The information stored in the array
// b_angle is
// b_angle[0]=sin_theta_b
// b_angle[1]=cos_theta_b
// b_angle[2]=sin_phi_b
// b_angle[3]=cos_phi_b
  
    
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
 
  for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je; ++j){

         Real *fc1 = &(u_tc(TCF1,k,j,0));
         Real *fc2 = &(u_tc(TCF2,k,j,0));
         Real *fc3 = &(u_tc(TCF3,k,j,0));

         // The angle of B
         Real *sint_b = &(ptc->b_angle(0,k,j,0));
         Real *cost_b = &(ptc->b_angle(1,k,j,0));
         Real *sinp_b = &(ptc->b_angle(2,k,j,0));
         Real *cosp_b = &(ptc->b_angle(3,k,j,0));
     	 

      for(int i=is; i<=ie; ++i){

         Real ftc1 = fc1[i];
         Real ftc2 = fc2[i];
         Real ftc3 = fc3[i];
         //it should be e/T
         Real rho = w(IDN,k,j,i)/gamma_1;

         // in the case with magnetic field
        // rotate the vectors to oriante to the B direction
         if(MAGNETIC_FIELDS_ENABLED){
           // Apply rotation of the vectorsw
           RotateVec(sint_b[i],cost_b[i],sinp_b[i],cosp_b[i],
           	                                  ftc1,ftc2,ftc3);

         }

         Real kappa_x = ptc->kappa(0,k,j,i);
         Real kappa_y = ptc->kappa(1,k,j,i);
         Real kappa_z = ptc->kappa(2,k,j,i);


         // Now update the momentum equation
         //\partial F/\partial t=- rho * V_m^2/kappa F
         // The solution is 
         // (1+dt * rho V_m^2/kappa)F_new=F_old


         Real dtsigma1 = dt * rho * vlim * vlim/kappa_x;
         Real dtsigma2 = dt * rho * vlim * vlim/kappa_y;
         Real dtsigma3 = dt * rho * vlim * vlim/kappa_z;


         Real newfr1 = ftc1 / (1.0 + dtsigma1);
         Real newfr2 = ftc2 / (1.0 + dtsigma2);
         Real newfr3 = ftc3 / (1.0 + dtsigma3);
   

        // Now apply the invert rotation
         if(MAGNETIC_FIELDS_ENABLED){
           // Apply rotation of the vectors
           InvRotateVec(sint_b[i],cost_b[i],sinp_b[i],cosp_b[i], 
                                           newfr1,newfr2,newfr3);
         }        

          
         u_tc(1,k,j,i) = newfr1;
         u_tc(2,k,j,i) = newfr2;
         u_tc(3,k,j,i) = newfr3;

         //add energy source term back
         if(NON_BAROTROPIC_EOS)
           u(IEN,k,j,i) += tc_esource_(k,j,i);
         
      }// end i
    }// end j
  }// end k
      
}


