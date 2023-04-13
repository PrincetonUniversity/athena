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
#include "../radiation/implicit/radiation_implicit.hpp"

//void  LoadRadVariable(MeshBlock *pmb);

void RadConstantFluxInnerX1(
     MeshBlock *pmb, Coordinates *pco, Radiation *prad, 
     const AthenaArray<Real> &w, FaceField &b, 
     AthenaArray<Real> &ir,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void RadConstantFluxOuterX1(
     MeshBlock *pmb, Coordinates *pco, Radiation *prad, 
     const AthenaArray<Real> &w, FaceField &b, 
     AthenaArray<Real> &ir,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void HydroOuterX1(
    MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void HydroInnerX1(
    MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief homogeneous sphere test
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
    
  Real tgas=1.0;
  Real er=1.0;
  Real rho0=1.0;
  
  Real gamma = peos->GetGamma();
  
  // Initialize hydro variable
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        Real x1 = pcoord->x1v(i);
        if (x1 < 1.0){
            phydro->u(IDN,k,j,i) = rho0;
            phydro->u(IEN,k,j,i) = tgas * phydro->u(IDN,k,j,i)/(gamma-1.0);
        }
        else {
            phydro->u(IDN,k,j,i) = 1e-7;
            phydro->u(IEN,k,j,i) = 0.0;
        }
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
      }
    }
  }
  
  // Initialize opacity and specific intensity
  if(RADIATION_ENABLED || IM_RADIATION_ENABLED){
    int nfreq = prad->nfreq;

    Real *ir_lab;
    
    for(int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          
          ir_lab = &(prad->ir(k,j,i,0));
          for(int n=0; n<prad->n_fre_ang; n++){
             ir_lab[n] = er;
          }
          
          for (int ifr=0; ifr < nfreq; ++ifr){
            prad->sigma_s(k,j,i,ifr) = 0.0;
            if (pcoord->x1v(i) < 1.0) {
                prad->sigma_a(k,j,i,ifr) = 100.0;
                prad->sigma_pe(k,j,i,ifr) = 100.0;
                prad->sigma_p(k,j,i,ifr) = 100.0;

            }
            else {
                prad->sigma_a(k,j,i,ifr) = 0.0;
                prad->sigma_pe(k,j,i,ifr) = 0.0;
                prad->sigma_p(k,j,i,ifr) = 0.0;
            }
          }
        }
      }
    }
    
  }// End Rad
  
  return;
}



//The following returns hydro variables to their initial values at each timestep
void MeshBlock::UserWorkInLoop(void)
{
    Real rho0=1.0;
    Real tgas=1.0;
    Real gamma = peos->GetGamma();
    
    for(int k=ks; k<=ke; k++) {
        for(int j=js; j<=je; j++) {
            for(int i=is; i<=ie; i++) {
                Real x1 = pcoord->x1v(i);
                if (x1 < 1.0){
                    phydro->u(IDN,k,j,i) = rho0;
                    phydro->u(IEN,k,j,i) = tgas * phydro->u(IDN,k,j,i)/(gamma-1.0);
                }
                else {
                    phydro->u(IDN,k,j,i) = 1e-7;
                    phydro->u(IEN,k,j,i) = 0.0;
                }
                phydro->u(IM1,k,j,i) = 0.0;
                phydro->u(IM2,k,j,i) = 0.0;
                phydro->u(IM3,k,j,i) = 0.0;
            }
        }
    }
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  //ptpar = new TracePar(this, pin);

  if(RADIATION_ENABLED) {
    // Enroll Opacity Function
    //prad->EnrollOpacityFunction(Opacity);

  }

  return;
}


void Mesh::InitUserMeshData(ParameterInput *pin) {

  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, HydroInnerX1);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, HydroOuterX1);

  if(RADIATION_ENABLED || IM_RADIATION_ENABLED) {
    // Enroll radiation boundary functions
    EnrollUserRadBoundaryFunction(BoundaryFace::inner_x1, RadConstantFluxInnerX1);
    EnrollUserRadBoundaryFunction(BoundaryFace::outer_x1, RadConstantFluxOuterX1);
  }

  return;
}


void RadConstantFluxInnerX1(
     MeshBlock *pmb, Coordinates *pco, Radiation *prad, 
     const AthenaArray<Real> &w, FaceField &b, 
     AthenaArray<Real> &ir,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        for(int ifr=0; ifr<prad->nfreq; ++ifr){
          for(int n=0; n<prad->nang; ++n){
            Real &miuz = prad->mu(0,k,j,is-i,n);
            Real &weight = prad->wmu(n);
            Real ris = pco->x1v(is);
            Real ris1 = pco->x1v(is+1);
            Real r= pco->x1v(is-i);
            Real didr = (prad->ir(k,j,is+1,ifr*prad->nang+n)-prad->ir(k,j,is,ifr*prad->nang+n))/
                         (ris1 - ris);
/*            if(miuz > 0.0){
              prad->ir(k,j,is-i,ifr*prad->nang+n) = 1.;
            } else {
              prad->ir(k,j,is-i,ifr*prad->nang+n) = 1.;
            }
*/
//            prad->ir(k,j,is-i,ifr*prad->nang+n) = prad->ir(k,j,is,ifr*prad->nang+n)
//                              - didr *(ris - r);

             prad->ir(k,j,is-i,ifr*prad->nang+n) = 1.0;

            //if ((n == 0) || (n == prad->nang-1)) {
            //  std::cout << prad->ir(k,j,is-i,ifr*prad->nang+n) << " " << n << " " << i << std::endl;
            //}
          }
        }
      }
    }
  }

}


void RadConstantFluxOuterX1(
     MeshBlock *pmb, Coordinates *pco, Radiation *prad, 
     const AthenaArray<Real> &w, FaceField &b, 
     AthenaArray<Real> &ir,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        for(int ifr=0; ifr<prad->nfreq; ++ifr){
          for(int n=0; n<prad->nang; ++n){
            Real &miuz = prad->mu(0,k,j,ie+i,n);
            Real &weight = prad->wmu(n);
            if(miuz > 0.0){
              prad->ir(k,j,ie+i,ifr*prad->nang+n) = prad->ir(k,j,ie,ifr*prad->nang+n);
            } else {
              prad->ir(k,j,ie+i,ifr*prad->nang+n) = 0.;
            }
            //if ((n == 0) || (n == prad->nang-1)) {
            //  std::cout << prad->ir(k,j,is-i,ifr*prad->nang+n) << " " << n << " " << i << std::endl;
            //}
          }
        }
      }
    }
  }

}

void HydroInnerX1(
    MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{


  for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          prim(IDN,k,j,is-i) = 1.0;
          prim(IVX,k,j,is-i) = 0.0;
          prim(IVY,k,j,is-i) = 0.0;
          prim(IVZ,k,j,is-i) = 0.0;
          prim(IPR,k,j,is-i) = 1.0;
        }
      }
    }
  

}

void HydroOuterX1(
    MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          prim(IDN,k,j,ie+i) = 1.e-7;
          prim(IVX,k,j,ie+i) = 0.0;
          prim(IVY,k,j,ie+i) = 0.0;
          prim(IVZ,k,j,ie+i) = 0.0;
          prim(IPR,k,j,ie+i) = 1.e-12;

        }
      }
    }

}

