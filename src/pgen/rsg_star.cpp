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


// The space for opacity table

static AthenaArray<Real> opacitytable;
static AthenaArray<Real> planckopacity;
static AthenaArray<Real> logttable;
static AthenaArray<Real> logrhottable;
static AthenaArray<Real> logttable_planck;
static AthenaArray<Real> logrhottable_planck;

static AthenaArray<Real> ini_profile;
static AthenaArray<Real> bd_data;



// The global variable

static Real consFr = 1.68545e-6;
static Real grav0 = 3.18807e-03;
static Real kappaes = 515.204;

static const Real rhounit = 2.3149e-8;
static const Real tunit = 7.0783e4;
static const Real lunit = 6.955e10;
static Real tfloor;
static Real rhofloor;

static Real lbottom=1016.3;
static int in_line=81921;

static Real rmax=1.0163e3;


//======================================================================================
/*! \file globaldisk.cpp
 *  \brief global accretion disk problem with radiation
 *
 *====================================================================================*/

void Inflow_X1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);


void Outflow_X2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void Outflow_rad_X2(MeshBlock *pmb, Coordinates *pco, Radiation *prad, 
     const AthenaArray<Real> &w, const AthenaArray<Real> &bc, AthenaArray<Real> &ir, 
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void Inflow_rad_X1(MeshBlock *pmb, Coordinates *pco, Radiation *prad, 
     const AthenaArray<Real> &w, const AthenaArray<Real> &bc, AthenaArray<Real> &ir, 
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void StarOpacity(MeshBlock *pmb, AthenaArray<Real> &prim);

void rossopacity(const Real rho, const Real tgas, Real &kappa, Real &kappa_planck);

void GravityPotential(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, 
  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);




void Mesh::InitUserMeshData(ParameterInput *pin)
{
  
    // Enroll boundary functions

  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, Outflow_X2);
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, Inflow_X1);
  
  tfloor = pin->GetOrAddReal("radiation", "tfloor", 0.001);
  rhofloor = pin->GetOrAddReal("hydro", "dfloor", 1.e-8);
  
  EnrollUserExplicitSourceFunction(GravityPotential);

  ini_profile.NewAthenaArray(in_line,5);
  FILE *fini;
  if ( (fini=fopen("./MESA_profile_combined.txt","r"))==NULL )
  {
     printf("Open input file error MESA profile");
     return;
  }  

  for(int j=0; j<in_line; j++){
    for(int i=0; i<5; i++){
      fscanf(fini,"%lf",&(ini_profile(j,i)));
    }
  }

  fclose(fini);
  bd_data.NewAthenaArray(2,NGHOST);
  

  if(RADIATION_ENABLED){
  
    EnrollUserRadBoundaryFunction(BoundaryFace::inner_x1, Inflow_rad_X1);
    EnrollUserRadBoundaryFunction(BoundaryFace::outer_x1, Outflow_rad_X2);
  
    // the opacity table
    opacitytable.NewAthenaArray(212,46);
    planckopacity.NewAthenaArray(138,37);
    logttable.NewAthenaArray(212);
    logrhottable.NewAthenaArray(46);
    logttable_planck.NewAthenaArray(138);
    logrhottable_planck.NewAthenaArray(37);
    
    // read in the opacity table
    FILE *fkappa, *flogt, *flogrhot, *fplanck, *flogt_planck, *flogrhot_planck;
      
    if ( (fkappa=fopen("./aveopacity_combined.txt","r"))==NULL )
    {
      printf("Open input file error aveopacity_combined");
      return;
    }

    if ( (fplanck=fopen("./PlanckOpacity.txt","r"))==NULL )
    {
      printf("Open input file error PlanckOpacity");
      return;
    }

    if ( (flogt=fopen("./logT.txt","r"))==NULL )
    {
      printf("Open input file error logT");
      return;
    }

    if ( (flogrhot=fopen("./logRhoT.txt","r"))==NULL )
    {
      printf("Open input file error logRhoT");
      return;
    }

    if ( (flogt_planck=fopen("./logT_planck.txt","r"))==NULL )
    {
      printf("Open input file error logT_planck");
      return;
    }

    if ( (flogrhot_planck=fopen("./logRhoT_planck.txt","r"))==NULL )
    {
      printf("Open input file error logRhoT_planck");
      return;
    }

    for(int j=0; j<212; j++){
      for(int i=0; i<46; i++){
          fscanf(fkappa,"%lf",&(opacitytable(j,i)));
      }
    }

    for(int j=0; j<138; j++){
      for(int i=0; i<37; i++){
          fscanf(fplanck,"%lf",&(planckopacity(j,i)));
      }
     }


    for(int i=0; i<46; i++){
      fscanf(flogrhot,"%lf",&(logrhottable(i)));
    }

    for(int i=0; i<212; i++){
      fscanf(flogt,"%lf",&(logttable(i)));
    }

    for(int i=0; i<37; i++){
      fscanf(flogrhot_planck,"%lf",&(logrhottable_planck(i)));
    }

    for(int i=0; i<138; i++){
      fscanf(flogt_planck,"%lf",&(logttable_planck(i)));
    }

    fclose(fkappa);
    fclose(flogt);
    fclose(flogrhot);
    fclose(fplanck);
    fclose(flogt_planck);
    fclose(flogrhot_planck);
 
  }
  

  return;
}

void MeshBlock::InitAfterRestart(ParameterInput *pin)
{

}

//======================================================================================
//! \fn void Mesh::TerminateUserMeshProperties(void)
//  \brief Clean up the Mesh properties
//======================================================================================
void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{

  ini_profile.DeleteAthenaArray();
  bd_data.DeleteAthenaArray();

  // free memory
  if(RADIATION_ENABLED){

    opacitytable.DeleteAthenaArray();
    logttable.DeleteAthenaArray();
    logrhottable.DeleteAthenaArray();
    planckopacity.DeleteAthenaArray();
    logttable_planck.DeleteAthenaArray();
    logrhottable_planck.DeleteAthenaArray();

  }


  return;
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  
  // get bottom boundary condition
  Real rbottom = pcoord->x1v(is-1);
  if(rbottom < pmy_mesh->mesh_size.x1min){
    for(int i=1; i<=NGHOST; ++i){
      Real radius = pcoord->x1v(is-i);
      int lleft=0;

      int lright=1;
      while((radius > ini_profile(lright,0)) && (lright < in_line-1)){
        lright = lright+1;
      }
      if(lright - lleft > 1) lleft = lright -1;

      Real rho = ini_profile(lleft,2) + (radius - ini_profile(lleft,0)) *
                                (ini_profile(lright,2) - ini_profile(lleft,2))
                               /(ini_profile(lright,0) - ini_profile(lleft,0));
      Real tem = ini_profile(lleft,1) + (radius - ini_profile(lleft,0)) *
                                (ini_profile(lright,1) - ini_profile(lleft,1))
                               /(ini_profile(lright,0) - ini_profile(lleft,0));

      bd_data(0,i-1) = rho;
      bd_data(1,i-1) = tem;
    }

  }
 
  
  if(RADIATION_ENABLED){
    
      prad->EnrollOpacityFunction(StarOpacity);

  }else{

  }
  

  return;
}


void MeshBlock::UserWorkInLoop(void)
{
  if(RADIATION_ENABLED){
  
    
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
    Real gamma1 = peos->GetGamma() - 1.0;
    AthenaArray<Real> ir_cm;
    ir_cm.NewAthenaArray(prad->n_fre_ang);
      
     for (int k=kl; k<=ku; ++k){
      for (int j=jl; j<=ju; ++j){
       for (int i=il; i<=iu; ++i){
         
          Real& vx=phydro->w(IVX,k,j,i);
          Real& vy=phydro->w(IVY,k,j,i);
          Real& vz=phydro->w(IVZ,k,j,i);
         
          Real& rho=phydro->w(IDN,k,j,i);
          Real& pgas=phydro->w(IEN,k,j,i);



/*         if(iniflag){
              if(pcoord->x1v(i) > 50.0 && rho < 1.e-6){
                  phydro->w(IDN,k,j,i) = 1.e-7;
                  phydro->u(IDN,k,j,i) = 1.e-7;

              }

         }
*/
         
          Real vel = sqrt(vx*vx+vy*vy+vz*vz);

          if(vel > prad->vmax * prad->crat){
            Real ratio = prad->vmax * prad->crat / vel;
            vx *= ratio;
            vy *= ratio;
            vz *= ratio;
            
            phydro->u(IM1,k,j,i) = rho*vx;
            phydro->u(IM2,k,j,i) = rho*vy;
            phydro->u(IM3,k,j,i) = rho*vz;

            Real ke = 0.5 * rho * (vx*vx+vy*vy+vz*vz);
            
            Real pb=0.0;
            if(MAGNETIC_FIELDS_ENABLED){
               pb = 0.5*(SQR(pfield->bcc(IB1,k,j,i))+SQR(pfield->bcc(IB2,k,j,i))
                     +SQR(pfield->bcc(IB3,k,j,i)));
            }
            
            Real  eint = phydro->w(IEN,k,j,i)/gamma1;
            
            phydro->u(IEN,k,j,i) = eint + ke + pb;

          }
          //replace gas pressure with the same of gas and radiation
  /*        if(iniflag > 0){
              Real temp0 = pgas/rho;
              Real coef1 = prad->prat/3.0;
              Real coef2 = rho;
              Real coef3 = -pgas;
              Real gast;
              
              gast = Rtsafe(Tequilibrium, 0.0, temp0, 1.e-12, coef1, coef2, coef3, 0.0);
              if(gast < tfloor) gast = tfloor;
              pgas = rho * gast;
              Real eint = pgas/gamma1;
              Real ke=0.5*rho*vel*vel;
              Real pb=0.0;
              if(MAGNETIC_FIELDS_ENABLED){
                  pb = 0.5*(SQR(pfield->bcc(IB1,k,j,i))+SQR(pfield->bcc(IB2,k,j,i))
                            +SQR(pfield->bcc(IB3,k,j,i)));
              }
              phydro->u(IEN,k,j,i) = eint + ke + pb;

              //initialize the radiation quantity
              for(int n=0; n<prad->n_fre_ang; ++n)
                  ir_cm(n) = gast * gast * gast * gast;
              
              Real *mux = &(prad->mu(0,k,j,i,0));
              Real *muy = &(prad->mu(1,k,j,i,0));
              Real *muz = &(prad->mu(2,k,j,i,0));
              
              Real *ir_lab = &(prad->ir(k,j,i,0));
              
              prad->pradintegrator->ComToLab(vx,vy,vz,mux,muy,muz,ir_cm,ir_lab);
            
               
          }
  */
      }}}
      
      ir_cm.DeleteAthenaArray();
      
    }else{
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
      Real gamma1 = peos->GetGamma() - 1.0;
    
      for (int k=kl; k<=ku; ++k){
       for (int j=jl; j<=ju; ++j){
        for (int i=il; i<=iu; ++i){
         
          Real& vx=phydro->w(IVX,k,j,i);
          Real& vy=phydro->w(IVY,k,j,i);
          Real& vz=phydro->w(IVZ,k,j,i);
         
          Real& rho=phydro->w(IDN,k,j,i);
          
          Real tgas=phydro->w(IEN,k,j,i)/rho;
/*
          if( pcoord->x1v(i) < 5.0){
            phydro->w(IDN,k,j,i) = rhofloor;
            phydro->u(IDN,k,j,i) = rhofloor;
            phydro->w(IEN,k,j,i) = rhofloor * tfloor;
            vx=0.0;
            vy=0.0;
            vz=0.0;
            phydro->u(IM1,k,j,i) = 0.0;
            phydro->u(IM2,k,j,i) = 0.0;
            phydro->u(IM3,k,j,i) = 0.0;
            Real eint = rhofloor * tfloor/gamma1;
            Real ke = 0.5 * rhofloor * (vx*vx+vy*vy+vz*vz);
            Real pb=0.0;
            if(MAGNETIC_FIELDS_ENABLED){
               pb = 0.5*(SQR(pfield->bcc(IB1,k,j,i))+SQR(pfield->bcc(IB2,k,j,i))
                     +SQR(pfield->bcc(IB3,k,j,i)));
            }
            
            phydro->u(IEN,k,j,i) = eint + ke + pb;
          
          }
  */  
        }
       }
      }
    
    }
  return;
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief beam test
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  
  Real gamma = peos->GetGamma();


  //initialize random number
  std::srand(gid);
 

  
  Real crat, prat;
  if(RADIATION_ENABLED){
    crat = prad->crat;
    prat = prad->prat;
  }else{
    crat = 9738.32;
    prat = 0.0;
  }
  
  Real amp = 0.0;
  
  
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

  // First, find density and temperature at rmax
  int lleft=0;
  int lright=1;
  while((rmax > ini_profile(lright,0)) && (lright < in_line-1)){
     lright = lright+1;
  }
  if(lright - lleft > 1) lleft = lright -1;
  
//  Real rho_rmax = ini_profile(lleft,2) + (rmax - ini_profile(lleft,0)) *
//                              (ini_profile(lright,2) - ini_profile(lleft,2))
//                             /(ini_profile(lright,0) - ini_profile(lleft,0));
//  Real tem_rmax = ini_profile(lleft,1) + (rmax - ini_profile(lleft,0)) *
//                              (ini_profile(lright,1) - ini_profile(lleft,1))
//                             /(ini_profile(lright,0) - ini_profile(lleft,0));

  Real rho_rmax = ini_profile(in_line-1,2);
  Real tem_rmax = ini_profile(in_line-1,1);

  Real grav_rmax = grav0*pow(lbottom/rmax,2.0);


  
  // Initialize hydro variable
  for(int i=is; i<=ie; ++i) {
    Real &x1 = pcoord->x1v(i); 
    Real tem = tem_rmax;
    Real rho = rho_rmax;   

    // get the position


    if(x1 > rmax){
      tem = tem_rmax;
      Real grav_local = grav_rmax * pow(x1/rmax,2.0);
      rho = rho_rmax * exp(-grav_local*(x1-rmax)/(tem_rmax));
      rho = std::max(rho,Real(1.e-8));
    }else{
      int lleft=0;

      int lright=1;
      while((x1 > ini_profile(lright,0)) && (lright < in_line-1)){
         lright = lright+1;
      }
      if(lright - lleft > 1) lleft = lright -1;
      
      rho = ini_profile(lleft,2) + (x1 - ini_profile(lleft,0)) *
                                  (ini_profile(lright,2) - ini_profile(lleft,2))
                                 /(ini_profile(lright,0) - ini_profile(lleft,0));
      tem = ini_profile(lleft,1) + (x1 - ini_profile(lleft,0)) *
                                  (ini_profile(lright,1) - ini_profile(lleft,1))
                                 /(ini_profile(lright,0) - ini_profile(lleft,0));

    }

    Real radflx = consFr * (lbottom/x1) * (lbottom/x1);

    
    if(rho > 0.1 ) amp = 5.e-2;
    else amp = 0.0;

    
 //   rho *= (1.0 + amp * ((double)rand()/(double)RAND_MAX-0.5));
    
    
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        phydro->u(IDN,k,j,i) = rho * (1.0 + amp * ((double)rand()/(double)RAND_MAX-0.5));
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS){

          phydro->u(IEN,k,j,i) = tem * rho/(gamma-1.0);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        }
        
        if(RADIATION_ENABLED){
          Real er = tem * tem * tem * tem;
          // geometric dialution
          
          for(int ifr=0; ifr<prad->nfreq; ++ifr){
            Real coefa = 0.0, coefb = 0.0;
            for(int n=0; n<prad->nang; ++n){
              // spherical polar coordinate
              Real &miuz = prad->mu(0,k,j,i,n);
              Real &weight = prad->wmu(n);
              if(miuz > 0.0){
                coefa += weight;
                coefb += (miuz * weight);
              }
            }
            
            for(int n=0; n<prad->nang; ++n){
              Real &miuz = prad->mu(0,k,j,i,n);
            
              if(miuz > 0.0){
                prad->ir(k,j,i,ifr*prad->nang+n) = 0.5 *
                                       (er/coefa + radflx/coefb);
              }else{
                prad->ir(k,j,i,ifr*prad->nang+n) = 0.5 *
                                       (er/coefa - radflx/coefb);
              
              }
            
            }
            
          }
        }// End Rad
 
      }// end j
    }// end k
  }// end i

  // Opacity will be set during initialization

  
  return;
}

void StarOpacity(MeshBlock *pmb, AthenaArray<Real> &prim)
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
      // electron scattering opacity
  Real kappas = 0.2 * (1.0 + 0.6);
  Real kappaa = 0.0;
  
  for (int k=kl; k<=ku; ++k) {
  for (int j=jl; j<=ju; ++j) {
  for (int i=il; i<=iu; ++i) {
  for (int ifr=0; ifr<prad->nfreq; ++ifr){
    Real rho  = prim(IDN,k,j,i);
    Real gast = std::max(prim(IEN,k,j,i)/rho,tfloor);


    Real kappa, kappa_planck;
    rossopacity(rho, gast, kappa, kappa_planck);



    if(kappa < kappas){
      if(gast < 0.14){
        kappaa = kappa;
        kappa = 0.0;
      }else{
        kappaa = 0.0;
      }
    }else{
      kappaa = kappa - kappas;
      kappa = kappas;
    }

    prad->sigma_s(k,j,i,ifr) = kappa * rho * rhounit * lunit;
    prad->sigma_a(k,j,i,ifr) = kappaa * rho * rhounit * lunit;
    prad->sigma_p(k,j,i,ifr) = kappa_planck*rho*rhounit*lunit;
    prad->sigma_pe(k,j,i,ifr) = prad->sigma_p(k,j,i,ifr);
  }    


 }}}

}



void rossopacity(const Real rho, const Real tgas, Real &kappa, Real &kappa_planck)
{
  
    
    Real logt = log10(tgas * tunit);
    Real logrhot = log10(rho* rhounit) - 3.0* logt + 18.0;
    int nrhot1_planck = 0;
    int nrhot2_planck = 0;
    
    int nrhot1 = 0;
    int nrhot2 = 0;

    while((logrhot > logrhottable_planck(nrhot2_planck)) && (nrhot2_planck < 36)){
      nrhot1_planck = nrhot2_planck;
      nrhot2_planck++;
    }
    if(nrhot2_planck==36 && (logrhot > logrhottable_planck(nrhot2_planck)))
      nrhot1_planck=nrhot2_planck;

    while((logrhot > logrhottable(nrhot2)) && (nrhot2 < 45)){
      nrhot1 = nrhot2;
      nrhot2++;
    }
    if(nrhot2==45 && (logrhot > logrhottable(nrhot2)))
      nrhot1=nrhot2;
  
  /* The data point should between NrhoT1 and NrhoT2 */
    int nt1_planck = 0;
    int nt2_planck = 0;
    int nt1 = 0;
    int nt2 = 0;
    while((logt > logttable_planck(nt2_planck)) && (nt2_planck < 137)){
      nt1_planck = nt2_planck;
      nt2_planck++;
    }
    if(nt2_planck==137 && (logt > logttable_planck(nt2_planck)))
      nt1_planck=nt2_planck;

    while((logt > logttable(nt2)) && (nt2 < 211)){
      nt1 = nt2;
      nt2++;
    }
    if(nt2==211 && (logt > logttable(nt2)))
      nt1=nt2;

  

    Real kappa_t1_rho1=opacitytable(nt1,nrhot1);
    Real kappa_t1_rho2=opacitytable(nt1,nrhot2);
    Real kappa_t2_rho1=opacitytable(nt2,nrhot1);
    Real kappa_t2_rho2=opacitytable(nt2,nrhot2);

    Real planck_t1_rho1=planckopacity(nt1_planck,nrhot1_planck);
    Real planck_t1_rho2=planckopacity(nt1_planck,nrhot2_planck);
    Real planck_t2_rho1=planckopacity(nt2_planck,nrhot1_planck);
    Real planck_t2_rho2=planckopacity(nt2_planck,nrhot2_planck);


    // in the case the temperature is out of range
    // the planck opacity should be smaller by the 
    // ratio T^-3.5
    if(nt2_planck == 137 && (logt > logttable_planck(nt2_planck))){
       Real scaling = pow(10.0, -3.5*(logt - logttable_planck(137)));
       planck_t1_rho1 *= scaling;
       planck_t1_rho2 *= scaling;
       planck_t2_rho1 *= scaling;
       planck_t2_rho2 *= scaling;
    }


    Real rho_1 = logrhottable(nrhot1);
    Real rho_2 = logrhottable(nrhot2);
    Real t_1 = logttable(nt1);
    Real t_2 = logttable(nt2);

    
    if(nrhot1 == nrhot2){
      if(nt1 == nt2){
        kappa = kappa_t1_rho1;
      }else{
        kappa = kappa_t1_rho1 + (kappa_t2_rho1 - kappa_t1_rho1) *
                                (logt - t_1)/(t_2 - t_1);
      }/* end same T*/
    }else{
      if(nt1 == nt2){
        kappa = kappa_t1_rho1 + (kappa_t1_rho2 - kappa_t1_rho1) *
                                (logrhot - rho_1)/(rho_2 - rho_1);
      }else{
        kappa = kappa_t1_rho1 * (t_2 - logt) * (rho_2 - logrhot)/
                                ((t_2 - t_1) * (rho_2 - rho_1))
              + kappa_t2_rho1 * (logt - t_1) * (rho_2 - logrhot)/
                                ((t_2 - t_1) * (rho_2 - rho_1))
              + kappa_t1_rho2 * (t_2 - logt) * (logrhot - rho_1)/
                                ((t_2 - t_1) * (rho_2 - rho_1))
              + kappa_t2_rho2 * (logt - t_1) * (logrhot - rho_1)/
                                ((t_2 - t_1) * (rho_2 - rho_1));
      }
    }/* end same rhoT */

    rho_1 = logrhottable_planck(nrhot1_planck);
    rho_2 = logrhottable_planck(nrhot2_planck);
    t_1 = logttable_planck(nt1_planck);
    t_2 = logttable_planck(nt2_planck);
 
  /* Now do the same thing for Planck mean opacity */
    if(nrhot1_planck == nrhot2_planck){
      if(nt1_planck == nt2_planck){
        kappa_planck = planck_t1_rho1;
      }else{
        kappa_planck = planck_t1_rho1 + (planck_t2_rho1 - planck_t1_rho1) *
                                (logt - t_1)/(t_2 - t_1);
      }/* end same T*/
    }else{
      if(nt1_planck == nt2_planck){
        kappa_planck = planck_t1_rho1 + (planck_t1_rho2 - planck_t1_rho1) *
                                (logrhot - rho_1)/(rho_2 - rho_1);

      }else{        
        kappa_planck = planck_t1_rho1 * (t_2 - logt) * (rho_2 - logrhot)/
                                ((t_2 - t_1) * (rho_2 - rho_1))
                     + planck_t2_rho1 * (logt - t_1) * (rho_2 - logrhot)/
                                ((t_2 - t_1) * (rho_2 - rho_1))
                     + planck_t1_rho2 * (t_2 - logt) * (logrhot - rho_1)/
                                ((t_2 - t_1) * (rho_2 - rho_1))
                     + planck_t2_rho2 * (logt - t_1) * (logrhot - rho_1)/
                                ((t_2 - t_1) * (rho_2 - rho_1));
      }
    }/* end same rhoT */

    return;

}

// This function sets boundary condition for primitive variables

void Inflow_X1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{

  
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {

        Real rho = bd_data(0,i-1);
        Real tem = bd_data(1,i-1);

        prim(IDN,k,j,is-i) = bd_data(0,i-1);
        prim(IVX,k,j,is-i) = 0.0;
        prim(IVY,k,j,is-i) = prim(IVY,k,j,is);
        prim(IVZ,k,j,is-i) = prim(IVZ,k,j,is);
        prim(IEN,k,j,is-i) = rho*tem;
        
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

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
          Real &x1g = pco->x1v(ie+i);
          Real &x1 = pco->x1v(ie+i-1);
          if(prim(IVX,k,j,ie) < 0.0){
            prim(IDN,k,j,ie+i) = prim(IDN,k,j,ie);
            prim(IVX,k,j,ie+i) = 0.0;
          }else{
            prim(IDN,k,j,ie+i) = prim(IDN,k,j,ie+i-1) * x1*x1/(x1g*x1g);
            prim(IVX,k,j,ie+i) = prim(IVX,k,j,ie);
          }
          prim(IVY,k,j,ie+i) = prim(IVY,k,j,ie);
          prim(IVZ,k,j,ie+i) = prim(IVZ,k,j,ie);
          prim(IEN,k,j,ie+i) = prim(IPR,k,j,ie);
        
      }
    }
  }
  

  return;
}


void Outflow_rad_X2(MeshBlock *pmb, Coordinates *pco, Radiation *prad, 
     const AthenaArray<Real> &w, const AthenaArray<Real> &bc, AthenaArray<Real> &ir, 
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
  
  //vacuum boundary condition at top   
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {

        for(int ifr=0; ifr<prad->nfreq; ++ifr){
          for(int n=0; n<prad->nang; ++n){
            Real miuz = prad->mu(0,k,j,ie+i,n);
            if(miuz > 0.0){
              ir(k,j,ie+i,ifr*prad->nang+n)
                            = ir(k,j,ie+i-1,ifr*prad->nang+n);
            }else{
              ir(k,j,ie+i,ifr*prad->nang+n) = 0.0;
            }
         
            
          }
        }
      }
    }
  }
  

  return;
}


//inpose a fix density, temperature, luminosity
void Inflow_rad_X1(MeshBlock *pmb, Coordinates *pco, Radiation *prad, 
     const AthenaArray<Real> &w, const AthenaArray<Real> &bc, AthenaArray<Real> &ir, 
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
  
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {

        Real &x1 = pco->x1v(is-i);    
        Real radflx = consFr * lbottom * lbottom /(x1 * x1);

        Real tem = bd_data(1,i-1);
        Real er = tem * tem * tem * tem;

        
        for(int ifr=0; ifr<prad->nfreq; ++ifr){
          Real coefa = 0.0, coefb = 0.0;
          for(int n=0; n<prad->nang; ++n){
            // spherical polar coordinate
            Real &miuz = prad->mu(0,k,j,is-i,n);
            Real &weight = prad->wmu(n);
            if(miuz > 0.0){
              coefa += weight;
              coefb += (miuz * weight);
            }
          }
          
          for(int n=0; n<prad->nang; ++n){
            Real &miuz = prad->mu(0,k,j,is-i,n);
          
            if(miuz > 0.0){
              ir(k,j,is-i,ifr*prad->nang+n) = 0.5 *
                                     (er/coefa + radflx/coefb);
            }else{
              ir(k,j,is-i,ifr*prad->nang+n) = 0.5 *
                                     (er/coefa - radflx/coefb);
            
            }
          
          }
          
        }
          
      }//i
    }//j
  }//k
  

  return;
}




//

Real grav_pot(const Real radius, const Real theta, const Real phi)
{
  // the companion is located at \theta=90, phi=0, r=rm1
  //x=rm1, y=0, z=0
  // current point r\sin\theta \cosphi, r\sin\theta\sin phi, r\cos\theta
  Real coef = -grav0*pow(lbottom,2.0);
  Real potphi = coef * pow(radius,-1.0);
  return potphi;
}


//Gravitaional acceration takes the form
// GM = grav0(r/r_0)^-2
// The potential is 
//phi = -grav0*r0 r^-1
// grav=-\partial \phi/\partial r

void GravityPotential(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, 
  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{

// add the effective tidal potential in the co-rotating frame
// -GM1/r-GM2/(r-R2)+GM2(r\cos\theta)/R2^2-0.5\Omega_0^2r^2\sin^2\theta
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
        cons(IEN,k,j,i) += (dt*(phidivrhov - divrhovphi));
        
        
      }
    }
  }

}




