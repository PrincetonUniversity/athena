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

const Real t_lbd=0.3;
const Real t_rbd=0.03;

static AthenaArray<Real> opacity_spectrum;


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
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, InjectHydro);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, VacuumHydro);

  if(RADIATION_ENABLED || IM_RADIATION_ENABLED){
    EnrollUserRadBoundaryFunction(BoundaryFace::inner_x1, Inject);
    EnrollUserRadBoundaryFunction(BoundaryFace::outer_x1, Vacuum);
  }


  opacity_spectrum.NewAthenaArray(52);

   FILE *fini;
   if ( (fini=fopen("./opacity.txt","r"))==NULL )
   {
     printf("Open input file error");
     return;
   }
   for(int i=0; i<52; ++i){
    fscanf(fini,"%lf",&(opacity_spectrum(i)));
   }

   fclose(fini);
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{

   AllocateUserOutputVariables(prad->n_fre_ang);
   // load the initial radiation from file

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
        Real x1 = pcoord->x1v(i);
        phydro->u(IDN,k,j,i) = 1.e5;
        Real velx = 0.0;

//        if(x1 >= 2 && x1 < 3.5)
//          velx = vel_amp * pow(sin(2*PI*(x1-2)/6.0),2.0);
//        else if(x1 >=3.5 && x1 < 6.5)
//          velx = vel_amp;
//        else if(x1 >=6.5 && x1 < 8)
//          velx = vel_amp * pow(sin(2*PI*(x1-2)/6.0),2.0);
//        else
//          velx = 0.0;



        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS){

          phydro->u(IEN,k,j,i) = phydro->u(IDN,k,j,i) * t_rbd/(gamma-1.0);
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
          for(int ifr=0; ifr<nfreq-1; ++ifr){
            for(int n=0; n<nang; ++n)
              prad->ir(k,j,i,ifr*nang+n) = pow(t_rbd,4.0)*prad->BlackBodySpec(prad->nu_grid(ifr)/t_rbd,prad->nu_grid(ifr+1)/t_rbd);
          }// end ifr
          // the last frequency group
          for(int n=0; n<nang; ++n)
            prad->ir(k,j,i,(nfreq-1)*nang+n) = pow(t_rbd,4.0)*(1.0-prad->FitBlackBody(prad->nu_grid(nfreq-1)/t_rbd));
          
        }
      }
    }

    for (int ifr=0; ifr < nfreq-1; ++ifr){
      for(int k=0; k<ncells3; ++k)
      for(int j=0; j<ncells2; ++j)
      for(int i=0; i<ncells1; ++i){
        Real x1 = pcoord->x1v(i);
        Real nu = prad->nu_cen(ifr); 
        prad->sigma_a(k,j,i,ifr) = opacity_spectrum(ifr);      
        prad->sigma_ae(k,j,i,ifr) = opacity_spectrum(ifr);
        prad->sigma_planck(k,j,i,ifr) = prad->sigma_ae(k,j,i,ifr);
        prad->sigma_s(k,j,i,ifr) = 0.0;    
      }
    }

    // the last bin
    // kappa=1
    for(int k=0; k<ncells3; ++k)
    for(int j=0; j<ncells2; ++j)
    for(int i=0; i<ncells1; ++i){
      Real x1 = pcoord->x1v(i);
      Real nu = prad->nu_grid(nfreq-1);
      prad->sigma_a(k,j,i,nfreq-1) = opacity_spectrum(nfreq-1);      
      prad->sigma_ae(k,j,i,nfreq-1) = prad->sigma_a(k,j,i,nfreq-1); 
      prad->sigma_planck(k,j,i,nfreq-1) = prad->sigma_ae(k,j,i,nfreq-1);
      prad->sigma_s(k,j,i,nfreq-1) = 0.0;  
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


  // inject radiation is gaussian distribution in nu
  // FWHM=2.0*sigma_width/3.0
  // sigma=FWHM/2.355
  // gauss=(1/\sqrt(2pi))exp(-(x-x0)^2/2\sigma^2)

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {

        for(int ifr=0; ifr<nfreq-1; ++ifr){
          for(int n=0; n<nang; ++n)
            prad->ir(k,j,is-i,ifr*nang+n) = pow(t_lbd,4.0)*prad->BlackBodySpec(prad->nu_grid(ifr)/t_lbd,prad->nu_grid(ifr+1)/t_lbd);
        }// end ifr
        // the last frequency group
        for(int n=0; n<nang; ++n)
          prad->ir(k,j,is-i,(nfreq-1)*nang+n) = pow(t_lbd,4.0)*(1.0-prad->FitBlackBody(prad->nu_grid(nfreq-1)/t_lbd));

        /*
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
              ir(k,j,is-i,nang*ifr+n) = 1.0; 

            }else{
              ir(k,j,is-i,nang*ifr+n) = 0.0;
            }
          }
          */



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
          if(n==NHYDRO-1)
           prim(n,k,j,is-i) = prim(IDN,k,j,is-i)*t_lbd;
          else
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



void MeshBlock::UserWorkInLoop(void)
{


  for(int k=ks; k<=ke; ++k)
    for(int j=js; j<=je; ++j)
      for(int i=is; i<=ie; ++i){
        for(int n=0; n<prad->n_fre_ang; ++n)
          user_out_var(n,k,j,i) = prad->ir(k,j,i,n);
      }

    Real gamma = peos->GetGamma();
    
  for(int k=ks; k<=ke; k++) {
   for(int j=js; j<=je; j++) {
     for(int i=is; i<=ie; i++) {

        phydro->u(IDN,k,j,i) = 1.e5;
        Real velx = 0.0;


        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS){

          phydro->u(IEN,k,j,i) = phydro->u(IDN,k,j,i) * t_rbd/(gamma-1.0);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        }

     }
   }
  }

  

}
