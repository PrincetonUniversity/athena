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


// the blackbody temperature of radiation field
const Real inject_tr = 1.0;
const Real vel_amp = 134.0;



//======================================================================================
/*! \file beam.cpp
 *  \brief Beam test for the radiative transfer module
 *
 *====================================================================================*/

void Inject(MeshBlock *pmb, Coordinates *pco, Radiation *prad, 
            const AthenaArray<Real> &w, FaceField &b,
            AthenaArray<Real> &ir,
            Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void InjectHydro(
    MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void Vacuum(MeshBlock *pmb, Coordinates *pco, Radiation *prad, 
            const AthenaArray<Real> &w, FaceField &b,
            AthenaArray<Real> &ir,
            Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void VacuumHydro(
    MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);


void Mesh::InitUserMeshData(ParameterInput *pin)
{
//  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, InjectHydro);
//  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, VacuumHydro);

  if(RADIATION_ENABLED || IM_RADIATION_ENABLED){
//    EnrollUserRadBoundaryFunction(BoundaryFace::inner_x1, Inject);
//    EnrollUserRadBoundaryFunction(BoundaryFace::outer_x1, Vacuum);
  }
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{

   AllocateUserOutputVariables(prad->n_fre_ang);
}


//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief beam test
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{


  
//  inject_tr = pin->GetOrAddReal("problem","inject_tr",1.0);

  
  Real gamma = peos->GetGamma();
  
  // Initialize hydro variable
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        phydro->u(IDN,k,j,i) = 1.0;
        Real velx = 0.0;
        Real x1 = pcoord->x1v(i);
        if(x1 >= 2 && x1 < 3.5)
          velx = vel_amp * pow(sin(2*PI*(x1-2)/6.0),2.0);
        else if(x1 >=3.5 && x1 < 6.5)
          velx = vel_amp;
        else if(x1 >=6.5 && x1 < 8)
          velx = vel_amp * pow(sin(2*PI*(x1-2)/6.0),2.0);
        else
          velx = 0.0;
        phydro->u(IM1,k,j,i) = vel_amp;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS){

          phydro->u(IEN,k,j,i) = 1.0/(gamma-1.0);
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
          for(int ifr=0; ifr<prad->nfreq-1; ++ifr){
            for(int n=0; n<nang; ++n)
              prad->ir(k,j,i,ifr*nang+n) = prad->BlackBodySpec(prad->nu_grid(ifr)/inject_tr,prad->nu_grid(ifr+1)/inject_tr);
          }// end ifr
          // the last frequency group
          for(int n=0; n<nang; ++n)
            prad->ir(k,j,i,(prad->nfreq-1)*nang+n) = 1.0-prad->FitBlackBody(prad->nu_grid(prad->nfreq-1)/inject_tr);
          
        }
      }
    }

    for(int k=0; k<ncells3; ++k)
     for(int j=0; j<ncells2; ++j)
       for(int i=0; i<ncells1; ++i){
          for (int ifr=0; ifr < nfreq; ++ifr){
            prad->sigma_s(k,j,i,ifr) = 0.0;
            prad->sigma_a(k,j,i,ifr) = 0.0;
            prad->sigma_ae(k,j,i,ifr) = 0.0;

          }

       }
    
    ir_cm.DeleteAthenaArray();
    
  }// End Rad
  
  return;
}




void Inject(MeshBlock *pmb, Coordinates *pco, Radiation *prad, 
          const AthenaArray<Real> &w, FaceField &b, 
          AthenaArray<Real> &ir,
          Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{

  int nang=prad->nang;
  int nfreq=prad->nfreq;

  Real blackbody = inject_tr * inject_tr * inject_tr * inject_tr;


  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {

        for(int ifr=0; ifr<nfreq; ++ifr){
          Real bd_emission = 0.0;
          if(ifr == nfreq - 1){
            bd_emission = 1.0 - prad->FitBlackBody(prad->nu_grid(ifr)/inject_tr);
            bd_emission *= blackbody;
          }else{
            bd_emission = prad->BlackBodySpec(prad->nu_grid(ifr)/inject_tr,
                                      prad->nu_grid(ifr+1)/inject_tr);
            bd_emission *= blackbody;
          }
          for(int n=0; n<nang; ++n){
            if(prad->mu(0,k,j,is-i,n) > 0.0){
              ir(k,j,is-i,nang*ifr+n) = bd_emission; 

            }else{
              ir(k,j,is-i,nang*ifr+n) = 0.0;
            }
          }// end angle
        }// end frequency
      }
    }
  }// end k
 

  return;
}

void InjectHydro(
    MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
  for (int n=0; n<NHYDRO; ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          prim(n,k,j,is-i) = prim(n,k,j,is);
        }
      }
    }
  }

}


void VacuumHydro(
    MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
  for (int n=0; n<NHYDRO; ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          prim(n,k,j,ie+i) = prim(n,k,j,ie);
        }
      }
    }
  }

}


void Vacuum(MeshBlock *pmb, Coordinates *pco, Radiation *prad, 
          const AthenaArray<Real> &w, FaceField &b, 
          AthenaArray<Real> &ir,
          Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{

  int nang=prad->nang;
  int nfreq=prad->nfreq;



  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {

        for(int ifr=0; ifr<nfreq; ++ifr){
          for(int n=0; n<nang; ++n){
            if(prad->mu(0,k,j,ie+i,n) > 0.0){
              ir(k,j,ie+i,nang*ifr+n) = ir(k,j,ie,nang*ifr+n); 
            }else{
              ir(k,j,ie,nang*ifr+n) = 0.0;
            }
          }// end angle
        }// end frequency
      }
    }
  }// end k
 

  return;
}



