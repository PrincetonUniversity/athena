///======================================================================================
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
#include <cstdlib>    // srand

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
#include "../utils/utils.hpp"

// The global parameters
// The global parameters
static Real kappaes = 9.38e3;
static Real kappaffp = 8.51705e9;
static Real kappaffr = 2.30191e8;
static Real rho0 = 0.298994; // Density normalization of torus
static Real inib0 = 0.0; // initial magnetic field strength
static Real tfloor; // temperature floor
static Real rhofloor; // density floor
static int bconf = 0; // bconf=1: pure B_phi
                      // bconf=0: vector potential proportional to density
                      // bconf=2: two loops

static Real gm1 = 1.02737e5;
static Real gm2 = 1.02737e4;
static Real qm = gm2/gm1;
static Real omega0 = 1.79925;
static Real rm2 = 32.6823;

//the initial thickness of the injected stream
static Real wheight=0.7;
//the initial injected radial velocity
static Real vr_l1=-0.447214;
static Real vphi_l1=-vr_l1*0.40035;
static Real t_l1=0.2;
static Real beta=1.0;
static Real amp=1.e-3;


//======================================================================================
/*! \file globaldisk.cpp
 *  \brief global accretion disk problem with radiation
 *
 *====================================================================================*/

void Inflow_X1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);


void Outflow_X2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);


void TidalPotential(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, 
  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);



void Mesh::InitUserMeshData(ParameterInput *pin)
{
  
    // Enroll boundary functions

  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, Outflow_X2);
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, Inflow_X1);

  tfloor = pin->GetOrAddReal("radiation", "tfloor", 0.001);
  rhofloor = pin->GetOrAddReal("hydro", "dfloor", 1.e-8);

  
  EnrollUserExplicitSourceFunction(TidalPotential);
  


  return;
}



//======================================================================================
//! \fn void Mesh::TerminateUserMeshProperties(void)
//  \brief Clean up the Mesh properties
//======================================================================================
void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{



  return;
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  


  return;
}


void MeshBlock::UserWorkInLoop(void)
{
  return;
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief beam test
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  
  Real gamma = peos->GetGamma();
  Real gm = pin->GetReal("problem", "GM");

  //initialize random number
  

  int kl=ks, ku=ke;
  if(ku > kl){
    ku += NGHOST;
    kl -= NGHOST;
  }
  int jl=js, ju=je;
  if(ju > jl){
    ju += NGHOST;
    jl -= NGHOST;
  }
  int il = is-NGHOST, iu=ie+NGHOST;

  
  // Initialize hydro variable
  for(int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {

        // Initialize the hydro quantity
        phydro->u(IDN,k,j,i) = rhofloor;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;

        if(NON_BAROTROPIC_EOS){
          phydro->u(IEN,k,j,i) = 0.1*rhofloor/(gamma-1.0);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        }
        
      }// i
    }// j
  }// k

// initialize interface B

  if (MAGNETIC_FIELDS_ENABLED) {

    for (int k=ks; k<=ke; ++k) {
      // reset loop limits for polar boundary
      jl=js; ju=je+1;
      if (pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar) jl=js+1;
      if (pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar) ju=je;
      for (int j=jl; j<=ju; ++j) {
        for (int i=is; i<=ie; ++i) {
          pfield->b.x2f(k,j,i) = 0.0;
        }
      }
    }

    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          pfield->b.x3f(k,j,i) = 0.0;
        }
      }
    }

    if (block_size.nx2 > 1) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie+1; ++i) {
            pfield->b.x1f(k,j,i) = 0.0;
          }
        }
      }
    }// nx2 
  }// End MHD

  
  return;
}

// This function sets boundary condition for primitive variables

void Inflow_X1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{

  
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
          prim(IDN,k,j,is-i) = prim(IDN,k,j,is);
          prim(IVX,k,j,is-i) = std::min(prim(IVX,k,j,is),0.0);
          prim(IVY,k,j,is-i) = prim(IVY,k,j,is);
          prim(IVZ,k,j,is-i) = prim(IVZ,k,j,is);
          prim(IEN,k,j,is-i) = prim(IEN,k,j,is);
        
      }
    }
  }
   // set magnetic field in inlet ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je; ++j){
      for(int i=1; i<=ngh; ++i){
        b.x1f(k,j,is-i) = b.x1f(k,j,is);
      }
    }}

    for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je+1; ++j){
      for(int i=1; i<=ngh; ++i){
        b.x2f(k,j,is-i) = b.x2f(k,j,is);
      }
    }}

    for(int k=ks; k<=ke+1; ++k){
    for(int j=js; j<=je; ++j){
      for(int i=1; i<=ngh; ++i){
        b.x3f(k,j,is-i) = b.x3f(k,j,is);
      }
    }}
    
  }


  return;
}

void Outflow_X2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
  //initialize random number
  std::srand(pmb->gid);
  Real Bheight = 0.5* wheight;
  
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
          Real radius=pco->x1v(ie+i);
          Real theta=pco->x2v(j);
          Real phi=pco->x3v(k);
          Real xpos=radius*sin(theta)*cos(phi);
          Real ypos=radius*sin(theta)*sin(phi);
          Real zpos=radius*cos(theta);
          Real dis=ypos*ypos+(xpos+radius)*(xpos+radius)+zpos*zpos;
          Real rhostream=rho0*exp(-dis/(wheight*wheight));

          if(rhostream < rhofloor){
            prim(IDN,k,j,ie+i) = prim(IDN,k,j,ie);
            prim(IVX,k,j,ie+i) = std::max(prim(IVX,k,j,ie),0.0);
            prim(IVY,k,j,ie+i) = prim(IVY,k,j,ie);
            prim(IVZ,k,j,ie+i) = prim(IVZ,k,j,ie);
            if(NON_BAROTROPIC_EOS) 
              prim(IEN,k,j,ie+i) = prim(IEN,k,j,ie);
          }else{
            prim(IDN,k,j,ie+i) = rhostream;
            prim(IVX,k,j,ie+i) = vr_l1;
            prim(IVY,k,j,ie+i) = 0.0;
            prim(IVZ,k,j,ie+i) = vphi_l1;
            if(NON_BAROTROPIC_EOS)
              prim(IEN,k,j,ie+i) = t_l1 * rhostream;
                
             // add random perturbation
//            if(pmb->pmy_mesh->time < 1.e-12){
//                a(IDN,k,j,ie+i) *= (1.0 + amp *
//                    ((double)rand()/(double)RAND_MAX-0.5));                
//            }
                
          }// end if
      }//i
    }// j
  }// k
   // set magnetic field in inlet ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    Real a0=0.05*sqrt(2.0*rho0*t_l1);
  
    for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je; ++j){
//#pragma simd
      for(int i=1; i<=ngh; ++i){
        Real radius=pco->x1v(ie+i);
        Real theta=pco->x2v(j);
        Real phi=pco->x3v(k);
        Real xpos=radius*sin(theta)*cos(phi);
        Real ypos=radius*sin(theta)*sin(phi);
        Real zpos=radius*cos(theta);
        Real dis=ypos*ypos+(xpos+radius)*(xpos+radius)+zpos*zpos;
        Real aphi = a0 * 2.0 * Bheight * cos(0.5 * PI * dis/(Bheight*Bheight))/PI;
        Real rhostream=rho0*exp(-dis/(wheight*wheight));
        dis=sqrt(dis);
        
        if(rhostream < rhofloor){
          b.x1f(k,j,ie+i+1) = b.x1f(k,j,ie+1);
        }else{
          b.x1f(k,j,ie+i+1) = aphi/(tan(theta)*radius)
                            - a0*sin(0.5*PI*dis/Bheight)*radius
                            *(1.0-cos(theta)*cos(phi))/dis;
        }
      }
    }}

    for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je+1; ++j){
//#pragma simd
      for(int i=1; i<=ngh; ++i){
        Real radius=pco->x1v(ie+i);
        Real theta=pco->x2v(j);
        Real phi=pco->x3v(k);
        Real xpos=radius*sin(theta)*cos(phi);
        Real ypos=radius*sin(theta)*sin(phi);
        Real zpos=radius*cos(theta);
        Real dis=ypos*ypos+(xpos+radius)*(xpos+radius)+zpos*zpos;
        Real aphi = a0 * 2.0 * Bheight * cos(0.5 * PI * dis/(Bheight*Bheight))/PI;
        Real rhostream=rho0*exp(-dis/(wheight*wheight));
        dis=sqrt(dis);
        
        if(rhostream < rhofloor){
          b.x2f(k,j,ie+i) = b.x2f(k,j,ie);
        }else{
          b.x2f(k,j,ie+i) = a0*sin(0.5*PI*dis/Bheight)*2.0*(1.0-sin(theta)*cos(phi))*radius/dis
                          - aphi/radius;
        }
      }
    }}

    for(int k=ks; k<=ke+1; ++k){
    for(int j=js; j<=je; ++j){
#pragma simd
      for(int i=1; i<=ngh; ++i){
          Real radius=pco->x1v(ie+i);
          Real theta=pco->x2v(j);
          Real phi=pco->x3v(k);
          Real xpos=radius*sin(theta)*cos(phi);
          Real ypos=radius*sin(theta)*sin(phi);
          Real zpos=radius*cos(theta);
          Real dis=ypos*ypos+(xpos+radius)*(xpos+radius)+zpos*zpos;
          Real rhostream=rho0*exp(-dis/(wheight*wheight));
          if(rhostream < rhofloor){
            b.x3f(k,j,ie+i) = b.x3f(k,j,ie);
          }else{
             b.x3f(k,j,ie+i)= sqrt(2.0*rho0*t_l1/beta);
          }
      }
    }}
    
  }
  

  return;
}





//


Real grav_pot(const Real radius, const Real theta, const Real phi)
{
  // the companion is located at \theta=90, phi=0, r=rm2
  //x=rm2, y=0, z=0
  // current point r\sin\theta \cosphi, r\sin\theta\sin phi, r\sin\theta
  Real dist_r2=sqrt(radius*radius+rm2*rm2-2.0*radius*rm2*sin(theta)*cos(phi));
  
  Real potphi=-gm1/radius-gm2/dist_r2-0.5*omega0*omega0*radius*radius
        *sin(theta)*sin(theta)+gm2*radius*sin(theta)*cos(phi)/(rm2*rm2);
  return potphi;
}


void TidalPotential(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, 
  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{

// add the effective tidal potential in the co-rotating frame
// -GM1/r-GM2/(r-R2)+GM2(r\cos\theta)/R2^2
  AthenaArray<Real> &x1flux=pmb->phydro->flux[X1DIR];
  AthenaArray<Real> &x2flux=pmb->phydro->flux[X2DIR];
  AthenaArray<Real> &x3flux=pmb->phydro->flux[X3DIR];
  
  for(int k=pmb->ks; k<=pmb->ke; ++k){
    for(int j=pmb->js; j<=pmb->je; ++j){
      for(int i=pmb->is; i<=pmb->ie; ++i){
        Real rho = prim(IDN,k,j,i);
        Real rcen = pmb->pcoord->x1v(i);
        Real rleft = pmb->pcoord->x1f(i);
        Real rright = pmb->pcoord->x1f(i+1);
        
        Real thetacen = pmb->pcoord->x2v(j);
        Real thetaleft = pmb->pcoord->x2f(j);
        Real thetaright = pmb->pcoord->x2f(j+1);
        
        Real phicen = pmb->pcoord->x3v(k);
        Real phileft = pmb->pcoord->x3f(k);
        Real phiright = pmb->pcoord->x3f(k+1);
        
        Real vol=pmb->pcoord->GetCellVolume(k,j,i);
        Real phic = grav_pot(rcen,thetacen,phicen);
        
        // radial direction
        
        Real phil = grav_pot(rleft,thetacen,phicen);
        Real phir = grav_pot(rright,thetacen,phicen);
        
        Real areal=pmb->pcoord->GetFace1Area(k,j,i);
        Real arear=pmb->pcoord->GetFace1Area(k,j,i+1);

        Real src = - dt * rho * (phir - phil)/pmb->pcoord->dx1f(i);
        cons(IM1,k,j,i) += src;
        Real phidivrhov = (arear*x1flux(IDN,k,j,i+1) -
                           areal*x1flux(IDN,k,j,i))*phic/vol;
        Real divrhovphi = (arear*x1flux(IDN,k,j,i+1)*phir -
                           areal*x1flux(IDN,k,j,i)*phil)/vol;
        if(NON_BAROTROPIC_EOS)
          cons(IEN,k,j,i) += (dt*(phidivrhov - divrhovphi));
        
        //theta direction

        phil = grav_pot(rcen,thetaleft,phicen);
        phir = grav_pot(rcen,thetaright,phicen);
        
        areal=0.5*(rright*rright-rleft*rleft)*fabs(sin(thetaleft))*
                   pmb->pcoord->dx3f(k);
        arear=0.5*(rright*rright-rleft*rleft)*fabs(sin(thetaright))*
                   pmb->pcoord->dx3f(k);
        
        src = - dt * rho * (phir - phil)/(rcen*pmb->pcoord->dx2f(j));
        cons(IM2,k,j,i) += src;
        phidivrhov = (arear*x2flux(IDN,k,j+1,i) -
                           areal*x2flux(IDN,k,j,i))*phic/vol;
        divrhovphi = (arear*x2flux(IDN,k,j+1,i)*phir -
                           areal*x2flux(IDN,k,j,i)*phil)/vol;

        if(NON_BAROTROPIC_EOS)
          cons(IEN,k,j,i) += (dt*(phidivrhov - divrhovphi));
        
        //phi direction
        
        phil = grav_pot(rcen,thetacen,phileft);
        phir = grav_pot(rcen,thetacen,phiright);
        
        areal=0.5*(rright*rright-rleft*rleft)*pmb->pcoord->dx2f(j);
        arear=areal;
        
        src = - dt * rho * (phir - phil)/(rcen*fabs(sin(thetacen))*
                                pmb->pcoord->dx3f(k));
        cons(IM3,k,j,i) += src;
        phidivrhov = (arear*x3flux(IDN,k+1,j,i) -
                           areal*x3flux(IDN,k,j,i))*phic/vol;
        divrhovphi = (arear*x3flux(IDN,k+1,j,i)*phir -
                           areal*x3flux(IDN,k,j,i)*phil)/vol;
  
        if(NON_BAROTROPIC_EOS)
          cons(IEN,k,j,i) += (dt*(phidivrhov - divrhovphi));
        
        // Add the coriolis force
       //dM/dt=-2\rho \Omega_0\times V
       // Omega_0=(\Omega_0\cos\theta,-\Omega_0\sin\theta,0)
       // because we use semi-implicit method, we need velocity
       // from conservative quantities
        rho = cons(IDN,k,j,i);
        Real vr=cons(IVX,k,j,i)/rho;
        Real vtheta=cons(IVY,k,j,i)/rho;
        Real vphi=cons(IVZ,k,j,i)/rho;
        Real dtomega = dt*omega0;
        Real sintheta=sin(thetacen);
        Real costheta=cos(thetacen);
        
        
        Real vphinew = -2.0 * sintheta*vr - 2.0*costheta*vtheta-(dtomega-1.0/dtomega)*vphi;
        vphinew /= (dtomega+1.0/dtomega);
        
        Real vrnew = dtomega * sintheta*vphinew + vr + dtomega*sintheta*vphi;
        Real vthetanew = dtomega * costheta*vphinew + vtheta + dtomega*costheta*vphi;
        
        cons(IM1,k,j,i) = vrnew * rho;
        
        cons(IM2,k,j,i) = vthetanew * rho;
        
        cons(IM3,k,j,i) = vphinew * rho;
        
      }
    }
  }

}




