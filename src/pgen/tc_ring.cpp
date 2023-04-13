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

// C++ headers
#include <iostream>   // endl
#include <fstream>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cmath>      // sqrt
#include <algorithm>  // min

// Athena++ headers
#include "../globals.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../hydro/hydro.hpp"
#include "../eos/eos.hpp"
#include "../bvals/bvals.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../field/field.hpp"
#include "../coordinates/coordinates.hpp"
#include "../thermal_conduction/tc.hpp"



//======================================================================================
/*! \file thermal_conduction.cpp
 *  \brief Test for thermal conduction
 *
 *====================================================================================*/


//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief beam test
//======================================================================================


static Real fixedkappa=0.01;
static Real fixedkappa1=1.e-8;

void FixedConductivity(MeshBlock *pmb, AthenaArray<Real> &prim, 
	                                    AthenaArray<Real> &bcc);

void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{ 
    
}



void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{

  if(TC_ENABLED){
    ptc->EnrollOpacityFunction(FixedConductivity);
  }

}



void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
    
  Real gamma = peos->GetGamma();

  Real vx=0.0;
  // Initialize hydro variable
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {

        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);

        Real radius=sqrt(x1*x1+x2*x2);
        Real tant=x2/x1;

        Real tgas = 1.0;
        Real rho = 1.0;

        if(radius < 0.7 && radius > 0.5 && fabs(tant) < 0.268 && x1 > 0.0)
        	tgas = 2.0;
      
        phydro->u(IDN,k,j,i) = rho;
        phydro->u(IM1,k,j,i) = vx;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS){
          phydro->u(IEN,k,j,i) = 0.5*vx*vx+tgas*rho/(gamma-1.0);
        }
        
        if(TC_ENABLED){
            ptc->u_tc(0,k,j,i) = tgas*rho/(gamma-1.0);
            ptc->u_tc(1,k,j,i) = 0.0;
            ptc->u_tc(2,k,j,i) = 0.0;
            ptc->u_tc(3,k,j,i) = 0.0;
        }
      }// end i
    }
  }


    // Add horizontal magnetic field lines, to show streaming and diffusion 
  // along magnetic field ines
  if(MAGNETIC_FIELDS_ENABLED){

    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          Real x1 = pcoord->x1v(i);
          Real x2 = pcoord->x2v(j);

          Real radius=sqrt(x1*x1+x2*x2);
          pfield->b.x1f(k,j,i) = -x2/radius;
        }
      }
    }

    if(block_size.nx2 > 1){

      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je+1; ++j) {
          for (int i=is; i<=ie; ++i) {

            Real x1 = pcoord->x1v(i);
            Real x2 = pcoord->x2v(j);

            Real radius=sqrt(x1*x1+x2*x2);
            pfield->b.x2f(k,j,i) = x1/radius;
          }
        }
      }

    }

    if(block_size.nx3 > 1){

      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
    }// end nx3

    // set cell centerd magnetic field
    // Add magnetic energy density to the total energy
    pfield->CalculateCellCenteredField(pfield->b,pfield->bcc,pcoord,is,ie,js,je,ks,ke);

    for(int k=ks; k<=ke; ++k){
      for(int j=js; j<=je; ++j){
        for(int i=is; i<=ie; ++i){
          phydro->u(IEN,k,j,i) +=
            0.5*(SQR((pfield->bcc(IB1,k,j,i)))
               + SQR((pfield->bcc(IB2,k,j,i)))
               + SQR((pfield->bcc(IB3,k,j,i))));
      
        }
      }
    }

  }// end MHD
  
  
  return;
}

void FixedConductivity(MeshBlock *pmb, AthenaArray<Real> &prim, 
	                                    AthenaArray<Real> &bcc)
{

  ThermalConduction *ptc=pmb->ptc;
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

  for(int k=kl; k<=ku; ++k){
    for(int j=jl; j<=ju; ++j){
#pragma omp simd
      for(int i=il; i<=iu; ++i){

        ptc->kappa(0,k,j,i) = fixedkappa;
        ptc->kappa(1,k,j,i) = fixedkappa1;
        ptc->kappa(2,k,j,i) = fixedkappa1;  

      }
    }
  }

  if(MAGNETIC_FIELDS_ENABLED){
    //First, calculate B_dot_grad_Pc
    for(int k=kl; k<=ku; ++k){
      for(int j=jl; j<=ju; ++j){
      // diffusion coefficient is calculated with respect to B direction
      // Use a simple estimate of Grad Pc
        for(int i=il; i<=iu; ++i){
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
