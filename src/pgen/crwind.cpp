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
#include "../cr/cr.hpp"
#include "../cr/integrators/cr_integrators.hpp"
#include <stdio.h>  // fopen and fwrite

// data for the cooling function




//temperature unit at 2e5K, cool_0=0.1414
//kappa_0=0.0414
//

static Real b0=1.0;
static Real ciso=1.0;
static Real rho0=1.0;
static Real pc0=1.0;
static Real ec0=3.0;
static Real isoT=1.0;
static Real vpot=10.0;

static Real rho1=1.031794457483471;
static Real pc1=1.0210855425165273;
static Real rho2=1.1000552652478548;
static Real pc2=1.0656379278445607;


//======================================================================================
/*! \file beam.cpp
 *  \brief Beam test for the radiative transfer module
 *
 *====================================================================================*/


//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief beam test
//======================================================================================

static Real sigma=1.e8;

void Diffusion(MeshBlock *pmb, AthenaArray<Real> &u_cr, 
        AthenaArray<Real> &prim, AthenaArray<Real> &bcc);

void TCKappa(MeshBlock *pmb, 
	         AthenaArray<Real> &prim, AthenaArray<Real> &bcc);

void GalacticPot(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, 
  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);

void FixMHDLeft(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, 
     int ks, int ke, int ngh);
void FixCRsourceLeft(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr, 
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie, 
    int js, int je, int ks, int ke, int ngh);

void FixMHDRight(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, 
     int ks, int ke, int ngh);
void FixCRsourceRight(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr, 
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie, 
    int js, int je, int ks, int ke, int ngh);

void Mesh::InitUserMeshData(ParameterInput *pin)
{

    
 
  EnrollUserBoundaryFunction(inner_x1, FixMHDLeft);
  EnrollUserBoundaryFunction(outer_x1, FixMHDRight);

  if(CR_ENABLED){
    EnrollUserCRBoundaryFunction(inner_x1, FixCRsourceLeft);
    EnrollUserCRBoundaryFunction(outer_x1, FixCRsourceRight);
  }

  EnrollUserExplicitSourceFunction(GalacticPot);

}


void MeshBlock::UserWorkInLoop(void)
{

/*
  int il=is, iu=ie, jl=js, ju=je, kl=ks, ku=ke;
  il -= NGHOST;
  iu += NGHOST;
  if(ju>jl){
     jl -= NGHOST;
     ju += NGHOST;
  }
  if(ku>kl){
    kl -= NGHOST;
    ku += NGHOST;
  }


  for (int k=kl; k<=ku; ++k){
    for (int j=jl; j<=ju; ++j){
     for (int i=il; i<=iu; ++i){
     
      Real& vx=phydro->w(IVX,k,j,i);

     
      Real& rho=phydro->w(IDN,k,j,i);

         if(iniflag){
          if(pcoord->x1v(i) > 50.0 && rho < 1.e-6){
              phydro->w(IDN,k,j,i) = 1.e-7;
              phydro->u(IDN,k,j,i) = 1.e-7;

          }

     }

     
      if(vx < 0.0){
        vx = 0.0;

        phydro->u(IM1,k,j,i) = 0.0;


      }

  }}}

*/
  return;

}



void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{ 



}


void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{


  AllocateUserOutputVariables(2);

  if(CR_ENABLED)
    pcr->EnrollOpacityFunction(Diffusion);


}



void MeshBlock::ProblemGenerator(ParameterInput *pin)
{

//  std::srand(gid);
//  Real amp = 0.1;

  // Initialize hydro variable
  for(int i=is; i<=ie; ++i) {
    Real &x1 = pcoord->x1v(i);
    


    for(int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {

        phydro->u(IDN,k,j,i) = rho0*pow(x1,-vpot*vpot/isoT);
        if(phydro->u(IDN,k,j,i) < 1.e-4) 
           phydro->u(IDN,k,j,i) = 1.e-4;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;

        if (NON_BAROTROPIC_EOS){

          phydro->u(IEN,k,j,i) = isoT * phydro->u(IDN,k,j,i)/(5.0/3.0-1.0);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        }
        
        if(CR_ENABLED){
            pcr->u_cr(CRE,k,j,i) = 3.0*pc0*pow(x1,-vpot*vpot/isoT);
            if(pcr->u_cr(CRE,k,j,i) < 1.e-4)
            	pcr->u_cr(CRE,k,j,i) = 1.e-4;
            pcr->u_cr(CRF1,k,j,i) = 0.0;
            pcr->u_cr(CRF2,k,j,i) = 0.0;
            pcr->u_cr(CRF3,k,j,i) = 0.0;
        }
      }// end i
    }
  }
  //Need to set opactiy sigma in the ghost zones
  if(CR_ENABLED){

  // Default values are 1/3
    int nz1 = block_size.nx1 + 2*(NGHOST);
    int nz2 = block_size.nx2;
    if(nz2 > 1) nz2 += 2*(NGHOST);
    int nz3 = block_size.nx3;
    if(nz3 > 1) nz3 += 2*(NGHOST);
    for(int k=0; k<nz3; ++k){
      for(int j=0; j<nz2; ++j){
        for(int i=0; i<nz1; ++i){
          pcr->sigma_diff(0,k,j,i) = sigma;
          pcr->sigma_diff(1,k,j,i) = sigma;
          pcr->sigma_diff(2,k,j,i) = sigma;
        }
      }
    }// end k,j,i

  }// End CR

    // Add horizontal magnetic field lines, to show streaming and diffusion 
  // along magnetic field ines
  if(MAGNETIC_FIELDS_ENABLED){

    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          pfield->b.x1f(k,j,i) = b0;
        }
      }
    }

    if(block_size.nx2 > 1){

      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je+1; ++j) {
          for (int i=is; i<=ie; ++i) {
            pfield->b.x2f(k,j,i) = 0.0;
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




void Diffusion(MeshBlock *pmb, AthenaArray<Real> &u_cr, 
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

  for(int k=kl; k<=ku; ++k){
    for(int j=jl; j<=ju; ++j){
#pragma omp simd
      for(int i=il; i<=iu; ++i){

        pcr->sigma_diff(0,k,j,i) = sigma;
        pcr->sigma_diff(1,k,j,i) = sigma;
        pcr->sigma_diff(2,k,j,i) = sigma;  

      }
    }
  }

  Real invlim=1.0/pcr->vmax;

  // The information stored in the array
  // b_angle is
  // b_angle[0]=sin_theta_b
  // b_angle[1]=cos_theta_b
  // b_angle[2]=sin_phi_b
  // b_angle[3]=cos_phi_b




  if(MAGNETIC_FIELDS_ENABLED){
    //First, calculate B_dot_grad_Pc
    for(int k=kl; k<=ku; ++k){
      for(int j=jl; j<=ju; ++j){
    // x component
        pmb->pcoord->CenterWidth1(k,j,il-1,iu+1,pcr->cwidth);
        for(int i=il; i<=iu; ++i){
          Real distance = 0.5*(pcr->cwidth(i-1) + pcr->cwidth(i+1))
                         + pcr->cwidth(i);
          Real dprdx=(u_cr(CRE,k,j,i+1)- u_cr(CRE,k,j,i-1))/3.0;
          dprdx /= distance;
          pcr->b_grad_pc(k,j,i) = bcc(IB1,k,j,i) * dprdx;
        }
    //y component
        if(ju > jl){
          pmb->pcoord->CenterWidth2(k,j-1,il,iu,pcr->cwidth1);       
          pmb->pcoord->CenterWidth2(k,j,il,iu,pcr->cwidth);
          pmb->pcoord->CenterWidth2(k,j+1,il,iu,pcr->cwidth2);

          for(int i=il; i<=iu; ++i){
            Real distance = 0.5*(pcr->cwidth1(i) + pcr->cwidth2(i))
                         + pcr->cwidth(i);
            Real dprdy=(u_cr(CRE,k,j+1,i)-  u_cr(CRE,k,j-1,i))/3.0;
            dprdy /= distance;
            pcr->b_grad_pc(k,j,i) += bcc(IB2,k,j,i) * dprdy;

          }
        }
    // z component
        if(ku > kl){
          pmb->pcoord->CenterWidth3(k-1,j,il,iu,pcr->cwidth1);       
          pmb->pcoord->CenterWidth3(k,j,il,iu,pcr->cwidth);
          pmb->pcoord->CenterWidth3(k+1,j,il,iu,pcr->cwidth2);

          for(int i=il; i<=iu; ++i){
            Real distance = 0.5*(pcr->cwidth1(i) + pcr->cwidth2(i))
                          + pcr->cwidth(i);
            Real dprdz=(u_cr(CRE,k+1,j,i) - u_cr(CRE,k-1,j,i))/3.0;
            dprdz /= distance;
            pcr->b_grad_pc(k,j,i) += bcc(IB3,k,j,i) * dprdz;

          // now only get the sign
//          if(pcr->b_grad_pc(k,j,i) > TINY_NUMBER) pcr->b_grad_pc(k,j,i) = 1.0;
//          else if(-pcr->b_grad_pc(k,j,i) > TINY_NUMBER) pcr->b_grad_pc(k,j,i) 
//            = -1.0;
//          else pcr->b_grad_pc(k,j,i) = 0.0;
          }
        }

      // now calculate the streaming velocity
      // streaming velocity is calculated with respect to the current coordinate 
      //  system
      // diffusion coefficient is calculated with respect to B direction
        for(int i=il; i<=iu; ++i){
          Real pb= bcc(IB1,k,j,i)*bcc(IB1,k,j,i)
                  +bcc(IB2,k,j,i)*bcc(IB2,k,j,i)
                  +bcc(IB3,k,j,i)*bcc(IB3,k,j,i);
          Real inv_sqrt_rho = 1.0/sqrt(prim(IDN,k,j,i));
          Real va1 = bcc(IB1,k,j,i)*inv_sqrt_rho;
          Real va2 = bcc(IB2,k,j,i)*inv_sqrt_rho;
          Real va3 = bcc(IB3,k,j,i)*inv_sqrt_rho;

          Real va = sqrt(pb/prim(IDN,k,j,i));

          Real dpc_sign = 0.0;
          if(pcr->b_grad_pc(k,j,i) > TINY_NUMBER) dpc_sign = 1.0;
          else if(-pcr->b_grad_pc(k,j,i) > TINY_NUMBER) dpc_sign = -1.0;
          
          pcr->v_adv(0,k,j,i) = -va1 * dpc_sign;
          pcr->v_adv(1,k,j,i) = -va2 * dpc_sign;
          pcr->v_adv(2,k,j,i) = -va3 * dpc_sign;

          // now the diffusion coefficient

          if(va < TINY_NUMBER){
            pcr->sigma_adv(0,k,j,i) = pcr->max_opacity;
          }else{
            pcr->sigma_adv(0,k,j,i) = fabs(pcr->b_grad_pc(k,j,i))
                                   /(sqrt(pb)* va * (4.0/3.0) 
                                    * invlim * u_cr(CRE,k,j,i)); 
          }

          pcr->sigma_adv(1,k,j,i) = pcr->max_opacity;
          pcr->sigma_adv(2,k,j,i) = pcr->max_opacity;  

          // Now calculate the angles of B
          Real bxby = sqrt(bcc(IB1,k,j,i)*bcc(IB1,k,j,i) +
                           bcc(IB2,k,j,i)*bcc(IB2,k,j,i));
          Real btot = sqrt(pb);
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

        }//        

      }// end j
    }// end k

  }// End MHD  
  else{



  for(int k=kl; k<=ku; ++k){
    for(int j=jl; j<=ju; ++j){
  // x component
      pmb->pcoord->CenterWidth1(k,j,il-1,iu+1,pcr->cwidth);
      for(int i=il; i<=iu; ++i){
         Real distance = 0.5*(pcr->cwidth(i-1) + pcr->cwidth(i+1))
                        + pcr->cwidth(i);
         Real grad_pr=(u_cr(CRE,k,j,i+1)- u_cr(CRE,k,j,i-1))/3.0;
         grad_pr /= distance;

         Real radius=pmb->pcoord->x1v(i);
         Real Bfield=b0/(radius*radius);
         Real rho=prim(IDN,k,j,i);
         Real va = sqrt(Bfield*Bfield/rho);


         if(va < TINY_NUMBER){
           pcr->sigma_adv(0,k,j,i) = sigma;
           pcr->v_adv(0,k,j,i) = 0.0;
         }else{
           Real sigma2 = fabs(grad_pr)/(va * (4.0/3.0) 
                             * invlim * u_cr(CRE,k,j,i)); 
           if(fabs(grad_pr) < TINY_NUMBER){
             pcr->sigma_adv(0,k,j,i) = 0.0;
             pcr->v_adv(0,k,j,i) = 0.0;
           }else{
             pcr->sigma_adv(0,k,j,i) = sigma2;
             pcr->v_adv(0,k,j,i) = -va * grad_pr/fabs(grad_pr);     
           }
        }

        pcr->sigma_adv(1,k,j,i) = pcr->max_opacity;
        pcr->sigma_adv(2,k,j,i) = pcr->max_opacity;
       
        pcr->v_adv(1,k,j,i) = 0.0;
        pcr->v_adv(2,k,j,i) = 0.0;




      }

    }
  }

  }
}


void FixCRsourceLeft(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr, 
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc, 
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie, 
    int js, int je, int ks, int ke, int ngh)
{

  if(CR_ENABLED){
    for (int n=0; n<(NCR); ++n) {
      if (n==(CRE)){
        for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma simd
          for (int i=1; i<=(ngh); ++i) {
            Real radius1=pmb->pcoord->x1v(is);
            Real radius2=pmb->pcoord->x1v(is-i);
//            u_cr(CRE,k,j,is-i) = 3.0*pc0*pow(radius2,-vpot*vpot/isoT);  // reflect 1-velocity
            if(i==1) u_cr(CRE,k,j,is-i)= 3.0 * pc1;
            else u_cr(CRE,k,j,is-i) = 3.0* pc2;
          }
        }}        
      }else if(n==CRF1){
        for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma simd
          for (int i=1; i<=(ngh); ++i) {
            Real radius2=pmb->pcoord->x1v(is-i);
            Real Bfield=b0/(radius2*radius2);
            Real rho=w(IDN,k,j,is-i);
            Real va = sqrt(Bfield*Bfield/rho);
            u_cr(CRF1,k,j,is-i) = (va+w(IVX,k,j,is-i)) * u_cr(CRE,k,j,is-i) 
                                    * 4.0/(3.0*pmb->pcr->vmax);  // reflect 1-velocity
          }
        }}        
      } 
      else {
        for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma simd
          for (int i=1; i<=(ngh); ++i) {
            u_cr(n,k,j,is-i) = u_cr(n,k,j,is);
          }
        }}
      }
    }
  }

}



void FixMHDLeft(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, 
     int ks, int ke, int ngh)
{

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(ngh); ++i) {
        Real radius1=pmb->pcoord->x1v(is);
        Real radius2=pmb->pcoord->x1v(is-i);
//        printf("i: %d radius2: %e\n",i,radius2);
//        prim(IDN,k,j,is-i) = prim(IDN,k,j,is) * 
//                             pow(radius2/radius1,-2.0*vpot*vpot/isoT);
//        prim(IDN,k,j,is-i) = rhobot(i-1);
//        prim(IDN,k,j,is-i) = rho0*pow(radius2,-vpot*vpot/isoT);
        if(i==1) prim(IDN,k,j,is-i) = rho1;
        else prim(IDN,k,j,is-i) = rho2;
        prim(IVX,k,j,is-i) = prim(IVX,k,j,is)*prim(IDN,k,j,is)*radius1*radius1/
                              (prim(IDN,k,j,is-i)*radius2*radius2);  // reflect 1-velocity
        prim(IVY,k,j,is-i) = prim(IVY,k,j,is);
        prim(IVZ,k,j,is-i) = prim(IVZ,k,j,is);
        if(NON_BAROTROPIC_EOS)
         prim(IEN,k,j,is-i) = prim(IEN,k,j,is);
      }
    }}

  

  // copy face-centered magnetic fields into ghost zones, reflecting b1
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) { 
    for (int j=js; j<=je; ++j) { 
#pragma simd
      for (int i=1; i<=(ngh); ++i) { 
//        b.x1f(k,j,(is-i)) = sqrt(2.0*const_pb);  // reflect 1-field
          b.x1f(k,j,(is-i)) =  b.x1f(k,j,is);
      } 
    }}
    if(je > js){ 
     for (int k=ks; k<=ke; ++k) {
     for (int j=js; j<=je+1; ++j) {
#pragma simd
      for (int i=1; i<=(ngh); ++i) {
        b.x2f(k,j,(is-i)) =  b.x2f(k,j,is);
      }
     }}  
    }
    if(ke > ks){        
     for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma simd
       for (int i=1; i<=(ngh); ++i) {
         b.x3f(k,j,(is-i)) =  b.x3f(k,j,is);
       }
      }}
    }
  }


}



void FixCRsourceRight(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr, 
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc, 
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie, 
    int js, int je, int ks, int ke, int ngh)
{


  if(CR_ENABLED){

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(ngh); ++i) {
        Real radius1=pmb->pcoord->x1v(ie);
        Real radius2=pmb->pcoord->x1v(ie+i);
        Real radius3=pmb->pcoord->x1v(ie+i-1);
        Real radius4=pmb->pcoord->x1v(ie+i-2);

        Real Bfield=b0/(radius2*radius2);
        Real rho=w(IDN,k,j,ie+i);
        Real va = sqrt(Bfield*Bfield/rho);
        Real vtot=(va+w(IVX,k,j,ie+i))/pmb->pcr->vmax;

        Real grad_Ec=(u_cr(CRE,k,j,ie+i-1) - u_cr(CRE,k,j,ie+i-2))/
                      (radius3-radius4);
        Real grad_Fc=(u_cr(CRF1,k,j,ie+i-1) - u_cr(CRF1,k,j,ie+i-2))/
                     (radius3-radius4);
        
        u_cr(CRE,k,j,ie+i) = u_cr(CRE,k,j,ie+i-1)+(radius2-radius3)*grad_Ec; 

/*        if(u_cr(CRF1,k,j,ie) > 0.0){
          u_cr(CRF1,k,j,ie+i) = u_cr(CRF1,k,j,ie+i-1)+(radius2-radius3)*grad_Fc; 
        }else{
          u_cr(CRF1,k,j,ie+i) = 0.0;
        }
*/
        u_cr(CRF1,k,j,ie+i) = u_cr(CRE,k,j,ie+i)/sqrt(3.0);
        u_cr(CRF2,k,j,ie+i) = u_cr(CRF2,k,j,ie+i);
        u_cr(CRF3,k,j,ie+i) = u_cr(CRF3,k,j,ie+i);

      }
    }}  



   }


}



void FixMHDRight(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, 
     int ks, int ke, int ngh)
{


  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(ngh); ++i) {
        Real radius1=pmb->pcoord->x1v(ie);
        Real radius2=pmb->pcoord->x1v(ie+i);
        Real radius3=pmb->pcoord->x1v(ie+i-1);
        Real radius4=pmb->pcoord->x1v(ie+i-2);
//        prim(IDN,k,j,is-i) = prim(IDN,k,j,is) * 
//                             pow(radius2/radius1,-2.0*vpot*vpot/isoT);
//        prim(IDN,k,j,is-i) = rhobot(i-1);

        Real grad_rho=(prim(IDN,k,j,ie+i-1) - prim(IDN,k,j,ie+i-2))/
                      (radius3-radius4);
//        prim(IDN,k,j,ie+i) = prim(IDN,k,j,ie+i-1)+(radius2-radius3)*grad_rho;
        
        if(prim(IVX,k,j,ie) > 0.0){

//          prim(IVX,k,j,ie+i) = prim(IVX,k,j,ie)*prim(IDN,k,j,ie)*radius1*radius1/
//                              (prim(IDN,k,j,ie+i)*radius2*radius2);
//          if(prim(IVX,k,j,ie+i) > 12.0)
//           printf("rho: %e v: %e rhoi: %e\n",prim(IDN,k,j,ie+i),prim(IVX,k,j,ie+i),prim(IDN,k,j,ie));
            prim(IVX,k,j,ie+i) = prim(IVX,k,j,ie);
            prim(IDN,k,j,ie+i) = prim(IDN,k,j,ie)*radius1*radius1/(radius2*radius2);
        }
        else{
          prim(IVX,k,j,ie+i) = 0.0;
          prim(IDN,k,j,ie+i) = prim(IDN,k,j,ie);
        }
        prim(IVY,k,j,ie+i) = prim(IVY,k,j,ie);
        prim(IVZ,k,j,ie+i) = prim(IVZ,k,j,ie);
        if(NON_BAROTROPIC_EOS)
         prim(IEN,k,j,ie+i) = prim(IEN,k,j,ie);
      }
    }}

  

  // copy face-centered magnetic fields into ghost zones, reflecting b1
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) { 
    for (int j=js; j<=je; ++j) { 
#pragma simd
      for (int i=1; i<=(ngh); ++i) { 
//        b.x1f(k,j,(is-i)) = sqrt(2.0*const_pb);  // reflect 1-field
          b.x1f(k,j,(ie+i+1)) =  b.x1f(k,j,ie+i);
      } 
    }}
    if(je > js){ 
     for (int k=ks; k<=ke; ++k) {
     for (int j=js; j<=je+1; ++j) {
#pragma simd
      for (int i=1; i<=(ngh); ++i) {
        b.x2f(k,j,(ie+i)) =  b.x2f(k,j,ie);
      }
     }}  
    }
    if(ke > ks){        
     for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma simd
       for (int i=1; i<=(ngh); ++i) {
         b.x3f(k,j,(ie+i)) =  b.x3f(k,j,ie);
       }
      }}
    }
  }

}




void GalacticPot(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, 
  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{

  AthenaArray<Real> &x1flux=pmb->phydro->flux[X1DIR];

  for(int k=pmb->ks; k<=pmb->ke; ++k){
    for(int j=pmb->js; j<=pmb->je; ++j){
      for(int i=pmb->is; i<=pmb->ie; ++i){
        Real rho = prim(IDN,k,j,i);
        Real phic = 2.0*vpot*vpot*log(pmb->pcoord->x1v(i));
        Real phil = 2.0*vpot*vpot*log(pmb->pcoord->x1f(i));
        Real phir = 2.0*vpot*vpot*log(pmb->pcoord->x1f(i+1));
        Real rr = pmb->pcoord->x1f(i+1);
        Real rl = pmb->pcoord->x1f(i);
        
        Real areal = rl * rl;
        Real arear = rr * rr;
        Real vol = (rr*rr*rr-rl*rl*rl)/3.0;
        Real src = - dt * rho * (phir - phil)/pmb->pcoord->dx1f(i);
        cons(IM1,k,j,i) += src;
        Real phidivrhov = (arear*x1flux(IDN,k,j,i+1) -
                           areal*x1flux(IDN,k,j,i))*phic/vol;
        Real divrhovphi = (arear*x1flux(IDN,k,j,i+1)*phir -
                           areal*x1flux(IDN,k,j,i)*phil)/vol;
        if(NON_BAROTROPIC_EOS)
          cons(IEN,k,j,i) += (dt*(phidivrhov - divrhovphi));
        
      }
    }
  }

}
