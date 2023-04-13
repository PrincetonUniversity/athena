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
#include "../radiation/implicit/radiation_implicit.hpp"
#include "../radiation/integrators/rad_integrators.hpp"

// File scope variables
// This test problem requires Pgas= rho * 0.6 * Tgas
static Real m0 = 1.2;

//======================================================================================
/*! \file radshock.cpp
 *  \brief radiation shock test for the radiative transfer module
 *
 *====================================================================================*/

void Inject(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void InjectRad(MeshBlock *pmb, Coordinates *pco, Radiation *prad, 
     const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir, 
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void InjectR(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void InjectRadR(MeshBlock *pmb, Coordinates *pco, Radiation *prad, 
     const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir, 
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  
    // Enroll boundary functions

  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, Inject);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, InjectR);
  EnrollUserRadBoundaryFunction(BoundaryFace::inner_x1, InjectRad);
  EnrollUserRadBoundaryFunction(BoundaryFace::outer_x1, InjectRadR);
  return;
}





//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief beam test
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  
  Real gamma = peos->GetGamma();
  // read input data file
  FILE *frho, *fpg, *fer, *ffr, *fvel, *fpos;
  if ((fer = fopen("data/E_4096.dat","r")) == NULL){
    printf("Open input data_E file error\n");
    return;
  }
  

  if ((fpg = fopen("data/p_4096.dat","r")) == NULL){
    printf("Open input data_Pgas file error\n");
    return;
  
  }
  
  if ((ffr = fopen("data/rF_4096.dat","r")) == NULL){
    printf("Open input data_Fr file error\n");
    return;
  
  }
  
  if ((frho = fopen("data/rho_4096.dat","r")) == NULL){
    printf("Open input data_rho file error\n");
    return;
  
  }
  
  if ((fvel = fopen("data/u_4096.dat","r")) == NULL){
    printf("Open input data_u file error\n");
    return;
  
  }
  
  if ((fpos = fopen("data/x_4096.dat","r")) == NULL){
    printf("Open input data_u file error\n");
    return;
  
  }

  Real rho, pg, er, fr, vel, pos;
  

  AthenaArray<Real> ir_cm;
  Real *ir_lab;
  
  if(RADIATION_ENABLED || IM_RADIATION_ENABLED)
    ir_cm.NewAthenaArray(prad->n_fre_ang);

  
  Real xmin = block_size.x1min;
  pos = -4.736027944533165168e-02;
  Real dx1 = pcoord->dx1v(is);
  // skip the lines not in this process
  while(pos + dx1 < xmin){
     fscanf(frho,"%lf",&(rho));
     fscanf(fpg,"%lf",&(pg));
     fscanf(fer,"%lf",&(er));
     fscanf(ffr,"%lf",&(fr));
     fscanf(fvel,"%lf",&(vel));
     fscanf(fpos,"%lf",&(pos));
  }
  
  // Initialize hydro variable
  for (int i=is; i<=ie; ++i) {
     fscanf(frho,"%lf",&(rho));
     fscanf(fpg,"%lf",&(pg));
     fscanf(fer,"%lf",&(er));
     fscanf(ffr,"%lf",&(fr));
     fscanf(fvel,"%lf",&(vel));
     fscanf(fpos,"%lf",&(pos));

     er *= 0.6*0.6*0.6*0.6;
    for(int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
      
        phydro->u(IDN,k,j,i) = rho;
        phydro->u(IM1,k,j,i) = rho * vel;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS){

          phydro->u(IEN,k,j,i) = pg/(gamma-1.0);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        }
        
        if(RADIATION_ENABLED || IM_RADIATION_ENABLED){
          for(int n=0; n<prad->n_fre_ang; ++n)
             ir_cm(n) = er;
          
          Real *mux = &(prad->mu(0,k,j,i,0));
          Real *muy = &(prad->mu(1,k,j,i,0));
          Real *muz = &(prad->mu(2,k,j,i,0));

          
          ir_lab = &(prad->ir(k,j,i,0));
          
          prad->pradintegrator->ComToLab(vel,0,0,mux,muy,muz,ir_cm,ir_lab);

        }
      }
    }
  }

  // opacity needs to include ghost zones
  int ncells1 = block_size.nx1 + 2*(NGHOST);
  int ncells2 = block_size.nx2;
  if(ncells2 > 1) ncells2 += 2*(NGHOST);
  int ncells3 = block_size.nx3;
  if(ncells3 > 1) ncells3 += 2*(NGHOST); 
  for(int k=0; k<ncells3; ++k){
    for(int j=0; j<ncells2; ++j){
      for(int i=0; i<ncells1; ++i){
  
          for(int ifr=0; ifr<prad->nfreq; ++ifr){
             prad->sigma_s(k,j,i,ifr) = 0.0;
             prad->sigma_a(k,j,i,ifr) = 577.4;
             prad->sigma_pe(k,j,i,ifr) = 577.4;
             prad->sigma_p(k,j,i,ifr) = 577.4;
          }

      }
    }
  }
  

  if(RADIATION_ENABLED || IM_RADIATION_ENABLED)
    ir_cm.DeleteAthenaArray(); 

  fclose(frho);
  fclose(fpg);
  fclose(fer);
  fclose(ffr);
  fclose(fvel);
  fclose(fpos);
  
  return;
}

//======================================================================================
//! \fn void MeshBlock::UserWorkInLoop(void)
//  \brief User-defined work function for every time step
//======================================================================================

void MeshBlock::UserWorkInLoop(void)
{
  // nothing to do
  return;
}

// This function sets boundary condition for primitive variables

void Inject(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)

{

  
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
          prim(IDN,k,j,is-i) = 1.0;
          prim(IVX,k,j,is-i) = m0;
          prim(IVY,k,j,is-i) = 0.0;
          prim(IVZ,k,j,is-i) = 0.0;
          prim(IEN,k,j,is-i) = 0.6;
        
      }
    }
  }

  return;
}

void InjectR(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)

{

  
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
          prim(IDN,k,j,ie+i) = prim(IDN,k,j,ie); 
          prim(IVX,k,j,ie+i) = prim(IVX,k,j,ie);
//         prim(IDN,k,j,ie+i) = 1.297321345222518518;
//          prim(IVX,k,j,ie+i) = 9.249828536461597306e-01;
          prim(IVY,k,j,ie+i) = 0.0;
          prim(IVZ,k,j,ie+i) = 0.0;
          prim(IEN,k,j,ie+i) = prim(IEN,k,j,ie);
//          prim(IEN,k,j,ie+i) = 0.6 * prim(IDN,k,j,ie+i)* 1.194751521050071652;
        
      }
    }
  }

  return;
}

void InjectRad(MeshBlock *pmb, Coordinates *pco, Radiation *prad, 
     const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir, 
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{

  AthenaArray<Real> ir_cm;
  Real *ir_lab;
  
  ir_cm.NewAthenaArray(prad->n_fre_ang);

  
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        for(int n=0; n<prad->n_fre_ang; ++n){
          ir_cm(n) = 1.0 * 0.6 * 0.6 * 0.6 * 0.6;
        }
        
          Real *mux = &(prad->mu(0,k,j,is-i,0));
          Real *muy = &(prad->mu(1,k,j,is-i,0));
          Real *muz = &(prad->mu(2,k,j,is-i,0));
        
          ir_lab = &(ir(k,j,is-i,0));
        
          prad->pradintegrator->ComToLab(m0,0,0,mux,muy,muz,ir_cm,ir_lab);
      }
    }
  }
  
  ir_cm.DeleteAthenaArray();

  return;
}


void InjectRadR(MeshBlock *pmb, Coordinates *pco, Radiation *prad, 
     const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir, 
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{

  AthenaArray<Real> ir_cm;
  Real *ir_lab;
  
  ir_cm.NewAthenaArray(prad->n_fre_ang);

  
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        for(int n=0; n<prad->n_fre_ang; ++n){
//          ir_cm(n) = pow(1.194751521050071652*0.6,4.0);
         ir(k,j,ie+i,n) = 2*ir(k,j,ie+i-1,n) - ir(k,j,ie+i-2,n);
        }
        
          Real *mux = &(prad->mu(0,k,j,ie+i,0));
          Real *muy = &(prad->mu(1,k,j,ie+i,0));
          Real *muz = &(prad->mu(2,k,j,ie+i,0));
        
          ir_lab = &(ir(k,j,ie+i,0));
          Real vx = w(IVX,k,j,ie);
        
//          prad->pradintegrator->ComToLab(vx,0,0,mux,muy,muz,ir_cm,ir_lab);
      }
    }
  }
  
  ir_cm.DeleteAthenaArray();

  return;
}




