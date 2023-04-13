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
void InjectRad(MeshBlock *pmb, Coordinates *pco, Radiation *prad, 
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

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief beam test
//======================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, InjectHydro);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, VacuumHydro);

  if(RADIATION_ENABLED || IM_RADIATION_ENABLED){
    EnrollUserRadBoundaryFunction(BoundaryFace::inner_x1, InjectRad);
    EnrollUserRadBoundaryFunction(BoundaryFace::outer_x1, Vacuum);
  }
}



void MeshBlock::ProblemGenerator(ParameterInput *pin)
{

  Real tgas, er1,er2,er3, sigma1,sigma2,sigma3;

  er1 = pin->GetOrAddReal("problem","er_1",1.0);
  er2 = pin->GetOrAddReal("problem","er_2",20.0);
  er3 = pin->GetOrAddReal("problem","er_3",30.0);    
  tgas = pin->GetOrAddReal("problem","tgas",1.0);
  sigma1 = pin->GetOrAddReal("problem","sigma_1",100.0);
  sigma2 = pin->GetOrAddReal("problem","sigma_2",200.0);
  sigma3 = pin->GetOrAddReal("problem","sigma_3",300.0);

  Real tr_ini = pow(er1,0.25);
  
  Real gamma = peos->GetGamma();
  Real scale=1;
  Real tot_tau = 100.0;
  Real sum_mass = 0.0;
  // Initialize hydro variable
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {

      for (int i=is; i<=ie; ++i) {
        Real x1 = pcoord->x1v(i);
        phydro->u(IDN,k,j,i) = exp(-(x1/scale)*(x1/scale));
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS){

          phydro->u(IEN,k,j,i) = phydro->u(IDN,k,j,i) * tgas/(gamma-1.0);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        }
        sum_mass = sum_mass + phydro->u(IDN,k,j,i) * (pcoord->x1f(i+1)-pcoord->x1f(i));
      }
    }
  }
  
  //Now initialize opacity and specific intensity
  if(RADIATION_ENABLED || IM_RADIATION_ENABLED){
    int nfreq = prad->nfreq;
    int nang = prad->nang;
    AthenaArray<Real> ir_cm;
    ir_cm.NewAthenaArray(prad->n_fre_ang);

    prad->kappa_es=tot_tau/sum_mass;

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
          for(int ifr=0; ifr<prad->nfreq; ++ifr){
            ir_lab = &(prad->ir(k,j,i,ifr*nang));
            Real emission = 1.e-20;
            // Initialize with blackbody spectrum
            if(ifr == nfreq-1){
              emission *= (1.0-prad->FitBlackBody(prad->nu_grid(ifr)/tr_ini));
            }else{
              emission *= prad->BlackBodySpec(prad->nu_grid(ifr)/tr_ini,
                                      prad->nu_grid(ifr+1)/tr_ini);
            }

/*            if(ifr == 0)
              er_ini = er1;
            else if(ifr == 1)
              er_ini = er2;
            else
              er_ini = er3;          
*/        
//          prad->pradintegrator->ComToLab(vx,vy,vz,mux,muy,muz,ir_cm,ir_lab);

            for(int n=0; n<prad->nang; n++){
               ir_lab[n] = emission;
            }

          }// end ifr
          
        }
      }
    }

    for(int k=0; k<ncells3; ++k)
     for(int j=0; j<ncells2; ++j)
       for(int i=0; i<ncells1; ++i){
          for (int ifr=0; ifr < nfreq; ++ifr){
            if(ifr == 0){
              prad->sigma_s(k,j,i,ifr) = prad->kappa_es*phydro->u(IDN,k,j,i);
              prad->sigma_a(k,j,i,ifr) = 0.0;
              prad->sigma_ae(k,j,i,ifr) = 0.0;
            }else if(ifr == 1){
              prad->sigma_s(k,j,i,ifr) = prad->kappa_es*phydro->u(IDN,k,j,i);
              prad->sigma_a(k,j,i,ifr) = 0.0;
              prad->sigma_ae(k,j,i,ifr) = 0.0;
            }else{
              prad->sigma_s(k,j,i,ifr) = prad->kappa_es*phydro->u(IDN,k,j,i);
              prad->sigma_a(k,j,i,ifr) = 0.0;
              prad->sigma_ae(k,j,i,ifr) =0.0;             
            }

          }

       }
    
    ir_cm.DeleteAthenaArray();
    
  }// End Rad
  
  return;
}


void InjectRad(MeshBlock *pmb, Coordinates *pco, Radiation *prad, 
          const AthenaArray<Real> &w, FaceField &b, 
          AthenaArray<Real> &ir,
          Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{

  int nang=prad->nang;
  int noct=prad->noct;
  int nfreq=prad->nfreq;
  int ang_oct=nang/noct;

  Real tr_ini = 1.0;


  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        for(int ifr=0; ifr<nfreq; ++ifr){
          Real *ir_lab = &(prad->ir(k,j,is-i,ifr*nang));
          Real emission = 1.0;
          // Initialize with blackbody spectrum
          if(ifr == nfreq-1){
            emission *= (1.0-prad->FitBlackBody(prad->nu_grid(ifr)/tr_ini));
          }else{
            emission *= prad->BlackBodySpec(prad->nu_grid(ifr)/tr_ini,
                                    prad->nu_grid(ifr+1)/tr_ini);
          }

          for(int n=0; n<prad->nang; n++){
            ir_lab[n] = emission;
          }

  
        }

    }}
  }

 

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
