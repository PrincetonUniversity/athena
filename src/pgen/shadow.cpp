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

// File scope variables
static int ang;
static int octnum;
static Real rho_min = 1.0;
static Real rho_max = 10.0;
static Real T_min = 1.0;
static Real T_max = 6.0;
static Real lx_box = 1.0;
static Real ly_box = 1.0;
static Real sigma0 = 1.0;

//======================================================================================
/*! \file beam.cpp
 *  \brief Beam test for the radiative transfer module
 *
 *====================================================================================*/

void TwoBeams(MeshBlock *pmb, Coordinates *pco, Radiation *prad, 
            const AthenaArray<Real> &w, FaceField &b,
            AthenaArray<Real> &ir,
            Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void TwoBeamHydro(
    MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void TwoBeamHydroR(
    MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void Vacuum(MeshBlock *pmb, Coordinates *pco, Radiation *prad,
          const AthenaArray<Real> &w, FaceField &b,
          AthenaArray<Real> &ir,
          Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void FFOpacity(MeshBlock *pmb, AthenaArray<Real> &prim);

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, TwoBeamHydro);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, TwoBeamHydroR);
  if(RADIATION_ENABLED || IM_RADIATION_ENABLED){
    EnrollUserRadBoundaryFunction(BoundaryFace::inner_x1, TwoBeams);
    EnrollUserRadBoundaryFunction(BoundaryFace::outer_x1, Vacuum);
  }

  lx_box=mesh_size.x1max - mesh_size.x1min;
  ly_box=mesh_size.x2max - mesh_size.x2min;
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  ang = pin->GetOrAddInteger("problem","ang",0);
  octnum = pin->GetOrAddInteger("problem","octnum",0);

  if(RADIATION_ENABLED || IM_RADIATION_ENABLED)
    prad->EnrollOpacityFunction(FFOpacity);
  


  return;
}


//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief beam test
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  
  Real gamma = peos->GetGamma();
  Real aaxis = lx_box/10.0;
  Real baxis = ly_box/10.0;
  
  // Initialize hydro variable
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        Real xpos = pcoord->x1v(i);
        Real ypos = pcoord->x2v(j);

        Real delta = 10.0 * (xpos*xpos/(aaxis * aaxis) + ypos*ypos/(baxis*baxis) - 1.0);
        Real rho = rho_min + (rho_max - rho_min)/(1.0 + exp(delta));
        Real press = rho_min * T_min;
        phydro->u(IDN,k,j,i) = rho;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS){

          phydro->u(IEN,k,j,i) = press/(gamma-1.0);
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
    for(int k=0; k<ncells3; ++k) {
      for (int j=0; j<ncells2; ++j) {
        for (int i=0; i<ncells1; ++i) {

          for (int ifr=0; ifr < nfreq; ++ifr){
            prad->sigma_s(k,j,i,ifr) = 0.0;
            prad->sigma_a(k,j,i,ifr) = sigma0 * rho_min* rho_min* pow(T_min, -3.5);
            prad->sigma_ae(k,j,i,ifr) = prad->sigma_a(k,j,i,ifr);
            prad->sigma_planck(k,j,i,ifr) = prad->sigma_a(k,j,i,ifr);
          }
          for(int n=0; n<prad->n_fre_ang; ++n){
              prad->ir(k,j,i,n) = 0.0;
          }
        }
      }
    }
  }
  
  return;
}




void TwoBeams(MeshBlock *pmb, Coordinates *pco, Radiation *prad, 
          const AthenaArray<Real> &w, FaceField &b,
          AthenaArray<Real> &ir,
          Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{

  int nang=prad->nang;
  int noct=prad->noct;
  int nfreq=prad->nfreq;
  int ang_oct=nang/noct;

  Real jr = pow(T_max, 4.0);


  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        Real &x1 = pco->x1v(i);
        Real &x2 = pco->x2v(js-j);
        for(int ifr=0; ifr<nfreq; ++ifr){
        for(int l=0; l<noct; ++l){
        for(int n=0; n<ang_oct; ++n){
          int n_ang=l*ang_oct + n;
          if(l ==0 || l == 2){
            if(n == 3)
              ir(k,j,is-i, n_ang+ifr*nang) = 0.5*jr/prad->wmu(n_ang+ifr*nang);
            else
              ir(k,j,is-i, n_ang+ifr*nang) = 0.0;
          }else{
              ir(k,j,is-i, n_ang+ifr*nang) = 0.0;
          }

        }//
        }
        }

    }}
  }

 

  return;
}


void TwoBeamHydro(
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


  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {

        for(int ifr=0; ifr<prad->nfreq; ++ifr){
          for(int n=0; n<prad->nang; ++n){
             Real &miux = prad->mu(0,k,j,ie+i,n);
             if(miux > 0.0)
               ir(k,j,ie+i,ifr*prad->nang+n) = ir(k,j,ie,ifr*prad->nang+n); 
             else
               ir(k,j,ie+i,ifr*prad->nang+n) = 0.0;
          }

        
        }

    }}
  }

 

  return;
}


void TwoBeamHydroR(
    MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
  for (int n=0; n<NHYDRO; ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          if(n==IVX)
            prim(n,k,j,ie+i) = std::max(prim(n,k,j,ie), 0.0);
          else
            prim(n,k,j,ie+i) = prim(n,k,j,ie);
        }
      }
    }
  }

}



void FFOpacity(MeshBlock *pmb, AthenaArray<Real> &prim)
{
  Radiation *prad = pmb->prad;
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
    Real rho  = prim(IDN,k,j,i);
    Real gast = prim(IEN,k,j,i)/rho;
    Real tpower= 1.0/(gast*gast*gast*sqrt(gast));

    prad->sigma_s(k,j,i,ifr) = 0.0;
    prad->sigma_a(k,j,i,ifr) = sigma0 * rho * rho * tpower;
    prad->sigma_ae(k,j,i,ifr) = prad->sigma_a(k,j,i,ifr);
    prad->sigma_planck(k,j,i,ifr) = prad->sigma_a(k,j,i,ifr);;
  }
  }}}

}


