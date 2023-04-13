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
#include "../radiation/radiation.hpp"
#include "../radiation/integrators/rad_integrators.hpp"


//======================================================================================
/*! \file beam.cpp
 *  \brief Beam test for the radiative transfer module
 *
 *====================================================================================*/

static Real eps0 = 1.e-4;

void Inflow_X1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void Outflow_X1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void Inflow_rad_X1(MeshBlock *pmb, Coordinates *pco, Radiation *prad,
     const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void Outflow_rad_X1(MeshBlock *pmb, Coordinates *pco, Radiation *prad,
     const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);


void Atmosphereopacity(MeshBlock *pmb, AthenaArray<Real> &prim);
//======================================================================================
//! \fn void Mesh::TerminateUserMeshProperties(void)
//  \brief Clean up the Mesh properties
//======================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin)
{

  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, Inflow_X1);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, Outflow_X1);
  if(RADIATION_ENABLED || IM_RADIATION_ENABLED){
    EnrollUserRadBoundaryFunction(BoundaryFace::inner_x1, Inflow_rad_X1);
    EnrollUserRadBoundaryFunction(BoundaryFace::outer_x1, Outflow_rad_X1);
  }

}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  
  
  if(RADIATION_ENABLED || IM_RADIATION_ENABLED){
    
      prad->EnrollOpacityFunction(Atmosphereopacity);
    

  }
  

  return;
}
//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief beam test
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{

  Real tgas, er;
  int flag=1;
  if(flag==1){
    tgas=1.0;
    er =1.0;
  
  }
  
  Real rho0=1.e-3;
  Real ztop = pmy_mesh->mesh_size.x1max;
  
  Real gamma = peos->GetGamma();
  
  // Initialize hydro variable
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        Real x3 = pcoord->x1v(i);
        phydro->u(IDN,k,j,i) = rho0 * exp(fabs(x3-ztop));
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS){

          phydro->u(IEN,k,j,i) = tgas * phydro->u(IDN,k,j,i)/(gamma-1.0);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        }
      }
    }
  }
  
  //Now initialize opacity and specific intensity
  if(RADIATION_ENABLED || IM_RADIATION_ENABLED){
    int nfreq = prad->nfreq;
    int nang = prad->nang;
    AthenaArray<Real> ir_cm;
    ir_cm.NewAthenaArray(prad->n_fre_ang);

    Real *ir_lab;
    
    for(int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
         

          Real vx = phydro->u(IM1,k,j,i)/phydro->u(IDN,k,j,i);
          Real vy = phydro->u(IM2,k,j,i)/phydro->u(IDN,k,j,i);
          Real vz = phydro->u(IM3,k,j,i)/phydro->u(IDN,k,j,i);
          Real *mux = &(prad->mu(0,k,j,i,0));
          Real *muy = &(prad->mu(1,k,j,i,0));
          Real *muz = &(prad->mu(2,k,j,i,0));
          
          ir_lab = &(prad->ir(k,j,i,0));
        
//          prad->pradintegrator->ComToLab(vx,vy,vz,mux,muy,muz,ir_cm,ir_lab);
          for(int n=0; n<prad->n_fre_ang; n++){
             ir_lab[n] = er;
          }
          
          for (int ifr=0; ifr < nfreq; ++ifr){
            prad->sigma_s(k,j,i,ifr) = (1-eps0) * phydro->u(IDN,k,j,i);
            prad->sigma_a(k,j,i,ifr) = eps0 * phydro->u(IDN,k,j,i);
            prad->sigma_pe(k,j,i,ifr) = eps0 * phydro->u(IDN,k,j,i);
            prad->sigma_p(k,j,i,ifr) = eps0 * phydro->u(IDN,k,j,i);
            
          }
        }
      }
    }
    
    ir_cm.DeleteAthenaArray();
    
  }// End Rad
  
  return;
}

//======================================================================================
//! \fn void MeshBlock::UserWorkInLoop(void)
//  \brief User-defined work function for every time step
//======================================================================================

void MeshBlock::UserWorkInLoop(void)
{
  // reset the flow property to the initial condition
  Real gamma = peos->GetGamma();
  Real rho0=1.e-3;
  Real tgas = 1.0;
  Real ztop = pmy_mesh->mesh_size.x1max;
    
  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {
        Real x3 = pcoord->x1v(i);
        phydro->u(IDN,k,j,i) = rho0 * exp(fabs(x3-ztop));
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS){

          phydro->u(IEN,k,j,i) = tgas * phydro->u(IDN,k,j,i)/(gamma-1.0);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        }


      }
    }
  }

  return;
}



void Inflow_X1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{

  
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
          prim(IDN,k,j,is-i) = prim(IDN,k,j,is);
          prim(IVX,k,j,is-i) = 0.0;
          prim(IVY,k,j,is-i) = prim(IVY,k,j,is);
          prim(IVZ,k,j,is-i) = prim(IVZ,k,j,is);
          prim(IEN,k,j,is-i) = prim(IEN,k,j,is);
        
      }
    }
  }


  return;
}



void Outflow_X1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{

  
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
          prim(IDN,k,j,ie+i) = prim(IDN,k,j,ie);
          prim(IVX,k,j,ie+i) = 0.0;
          prim(IVY,k,j,ie+i) = 0.0 * prim(IVY,k,j,ie);
          prim(IVZ,k,j,ie+i) = 0.0 * prim(IVZ,ke,j,ie);
          prim(IEN,k,j,ie+i) = prim(IEN,k,j,ie);
        
      }
    }
  }


  return;
}

void Inflow_rad_X1(MeshBlock *pmb, Coordinates *pco, Radiation *prad,
     const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{


  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
 
         for(int ifr=0; ifr<prad->nfreq; ++ifr){
            for(int n=0; n<prad->nang; ++n){
              ir(k,j,is-i,ifr*prad->nang+n) = 1.0;
          }

         }   

      }//i
    }//j
  }//k


  return;
}

void Outflow_rad_X1(MeshBlock *pmb, Coordinates *pco, Radiation *prad,
     const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{


  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) { 
         for(int ifr=0; ifr<prad->nfreq; ++ifr){
            for(int n=0; n<prad->nang; ++n){
              Real miux = prad->mu(0,k,j,ie+i,ifr*prad->nang+n);
              if(miux > 0.0){
                ir(k,j,ie+i,ifr*prad->nang+n) = ir(k,j,ie,ifr*prad->nang+n);
              }else{
              	ir(k,j,ie+i,ifr*prad->nang+n) = 0.0;
              }

            }
         }   

      }//i
    }//j
  }//k


  return;
}


void Atmosphereopacity(MeshBlock *pmb, AthenaArray<Real> &prim)
{
  Radiation *prad = pmb->prad;
  Hydro *phydro=pmb->phydro;
  int il = pmb->is; int jl = pmb->js; int kl = pmb->ks;
  int iu = pmb->ie; int ju = pmb->je; int ku = pmb->ke;
  il -= NGHOST;
  iu += NGHOST;
  if(ju > jl){
    jl -= NGHOST;
    ju += NGHOST;
  }
  if(ku > kl){
    kl -= NGHOST;
    ku += NGHOST;
  }

  
  for (int k=kl; k<=ku; ++k) {
  for (int j=jl; j<=ju; ++j) {
  for (int i=il; i<=iu; ++i) {
  for (int ifr=0; ifr<prad->nfreq; ++ifr){


    prad->sigma_s(k,j,i,ifr) = (1-eps0) * phydro->u(IDN,k,j,i);
    prad->sigma_a(k,j,i,ifr) = eps0 * phydro->u(IDN,k,j,i);
    prad->sigma_ae(k,j,i,ifr) = eps0 * phydro->u(IDN,k,j,i);
    prad->sigma_planck(k,j,i,ifr) = eps0 * phydro->u(IDN,k,j,i);           
          

  }// end ifr
 }}}

}

