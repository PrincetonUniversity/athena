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

// The space for opacity table

static AthenaArray<Real> opacitytable;
static AthenaArray<Real> planckopacity;
static AthenaArray<Real> logttable;
static AthenaArray<Real> logrhottable;

static AthenaArray<Real> ini_profile;
static AthenaArray<Real> bd_data;



// The global variable

static Real gm;
static Real consFr = 5.35879e-4;
static Real kappaes = 16.692;

static Real rhounit = 1.0e-9;
static Real tunit;
static Real lunit = 6.955e10;
static Real tfloor;

static Real lbottom=30.0;
static int ninputline = 167938;



void StarOpacity(MeshBlock *pmb, AthenaArray<Real> &prim);

//provide density and temperature, returns the opacity
void rossopacity(const Real rho, const Real tgas, Real &kappa, Real &kappa_planck);


void Inflow_X1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);


void Outflow_X2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);


void Inflow_rad_X1(MeshBlock *pmb, Coordinates *pco, Radiation *prad, 
     const AthenaArray<Real> &w, const AthenaArray<Real> &bc, AthenaArray<Real> &ir, 
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);


void Outflow_rad_X2(MeshBlock *pmb, Coordinates *pco, Radiation *prad, 
     const AthenaArray<Real> &w, const AthenaArray<Real> &bc, AthenaArray<Real> &ir, 
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);


void Mesh::InitUserMeshData(ParameterInput *pin)
{

   tfloor = pin->GetOrAddReal("radiation", "tfloor", 0.01);
   tunit = pin->GetOrAddReal("radiation","Tunit",1.55e5);
   gm = pin->GetOrAddReal("problem","GM",0.0);
  
   EnrollUserBoundaryFunction(BoundaryFace::inner_x1, Inflow_X1);
   EnrollUserBoundaryFunction(BoundaryFace::outer_x1, Outflow_X2);

   // the initial condition
   int lines=167938;

   ini_profile.NewAthenaArray(lines,3);
   FILE *fini;
   if ( (fini=fopen("./Input.txt","r"))==NULL )
   {
     printf("Open input file error MESA profile");
     return;
   }  

   for(int j=0; j<lines; j++){
     for(int i=0; i<3; i++){
       fscanf(fini,"%lf",&(ini_profile(j,i)));
     }
   }

   fclose(fini);

   bd_data.NewAthenaArray(2,NGHOST);
  

  
   if(RADIATION_ENABLED){
   
     EnrollUserRadBoundaryFunction(BoundaryFace::inner_x1, Inflow_rad_X1);
     EnrollUserRadBoundaryFunction(BoundaryFace::outer_x1, Outflow_rad_X2);
   
     
      // create the memory and read in the opacity table
     opacitytable.NewAthenaArray(138,37);
     planckopacity.NewAthenaArray(138,37);
     
     logttable.NewAthenaArray(138);
     logrhottable.NewAthenaArray(37);
     
      FILE *fkappa, *flogT, *flogrhoT, *fplanck;
     
      if ( (fkappa=fopen("./aveopacity_X0_2_solar.txt","r"))==NULL )
      {
         printf("Open input file error");
         return;
      }

      if ( (fplanck=fopen("PlanckOpacity.txt","r"))==NULL )
      {
         printf("Open input file error");
         return;
      }
  
      if ( (flogT=fopen("./logT.txt","r"))==NULL )
      {
         printf("Open input file error");
         return;
      }
  
      if ( (flogrhoT=fopen("./logRhoT.txt","r"))==NULL )
      {
         printf("Open input file error");
         return;
      }
  
     
      for(int j=0; j<138; j++){
        for(int i=0; i<37; i++){
          fscanf(fkappa,"%lf",&(opacitytable(j,i)));
        }
      }

      for(int j=0; j<138; j++){
        for(int i=0; i<37; i++){
          fscanf(fplanck,"%lf",&(planckopacity(j,i)));
        }
      }
   
      for(int i=0; i<37; i++){
        fscanf(flogrhoT,"%lf",&(logrhottable(i)));
      }
  
      for(int i=0; i<138; i++){
        fscanf(flogT,"%lf",&(logttable(i)));
      }
   
  
  
     fclose(fkappa);
     fclose(fplanck);
     fclose(flogT);
     fclose(flogrhoT);
   
   }// End Radiation

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
    Real gm1 = peos->GetGamma() - 1.0;
    
     for (int k=kl; k<=ku; ++k){
      for (int j=jl; j<=ju; ++j){
       for (int i=il; i<=iu; ++i){
         
          Real& vx=phydro->w(IVX,k,j,i);
          Real& vy=phydro->w(IVY,k,j,i);
          Real& vz=phydro->w(IVZ,k,j,i);
         
          Real& rho=phydro->w(IDN,k,j,i);
         
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
            
            Real  eint = phydro->w(IEN,k,j,i)/gm1;
            
            phydro->u(IEN,k,j,i) = eint + ke + pb;

          }
  
      }}}
    }
  return;
}






//======================================================================================
//! \fn void Mesh::TerminateUserMeshProperties(void)
//  \brief Clean up the Mesh properties
//======================================================================================
void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{
  

     ini_profile.DeleteAthenaArray();
     bd_data.DeleteAthenaArray();

     if(RADIATION_ENABLED){
     
       opacitytable.DeleteAthenaArray();
       planckopacity.DeleteAthenaArray();
       logttable.DeleteAthenaArray();
       logrhottable.DeleteAthenaArray();
    
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
      while((radius > ini_profile(lright,0)) && (lright < ninputline-1)){
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
  
  // the bottom value of z coordinate in the block
  AthenaArray<Real> height, tgas, density;
  
  Real amp = 0.0;
  
  
  height.NewAthenaArray(ninputline);
  tgas.NewAthenaArray(ninputline);
  density.NewAthenaArray(ninputline);
  
  
  // Initial profile from input data
  
  FILE *finput;
  if ( (finput=fopen("./Input.txt","r"))==NULL )
	{   
		printf("Open input file error");
		return;
	}
  for(int i=0; i<ninputline; ++i){
    fscanf(finput,"%lf",&(height(i)));
    fscanf(finput,"%lf",&(tgas(i)));
	 	fscanf(finput,"%lf",&(density(i)));
  }
  
  fclose(finput);



  //////////////////////////////////////////////
  
  // Initialize hydro variable
  for(int i=is; i<=ie; ++i) {
    Real &x1 = pcoord->x1v(i);
   
    
    // get the position
    int lleft=0;

    int lright=1;
    while((x1 > ini_profile(lright,0)) && (lright < ninputline-1)){
       lright = lright+1;
    }
    if(lright - lleft > 1) lleft = lright -1;
    
    Real rho = ini_profile(lleft,2) + (x1 - ini_profile(lleft,0)) *
                                (ini_profile(lright,2) - ini_profile(lleft,2))
                               /(ini_profile(lright,0) - ini_profile(lleft,0));
    Real tem = ini_profile(lleft,1) + (x1 - ini_profile(lleft,0)) *
                                (ini_profile(lright,1) - ini_profile(lleft,1))
                               /(ini_profile(lright,0) - ini_profile(lleft,0));
    
    if(rho > 5.e-2 && x1 > 25.0 && x1 < 35.0) amp = 5.e-2;
    else amp = 0.0;    
 //   rho *= (1.0 + amp * ((double)rand()/(double)RAND_MAX-0.5));
    
    
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        phydro->u(IDN,k,j,i) = rho;
        phydro->u(IM1,k,j,i) = rho * (amp * ((double)rand()/(double)RAND_MAX-0.5));
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
          Real radratio = (lbottom/x1) * (lbottom/x1);
          
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
                                       (er/coefa + consFr * radratio/coefb);
              }else{
                prad->ir(k,j,i,ifr*prad->nang+n) = 0.5 *
                                       (er/coefa - consFr * radratio/coefb);
              
              }
            
            }
            
          }
        }// End Rad
        
      }// end j
    }// end k
  }// end i

  
  
  
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
  Real kappas = 0.24;
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
      if(gast < 0.1){
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
    prad->sigma_ae(k,j,i,ifr) = prad->sigma_a(k,j,i,ifr);
    if(kappaa < kappa_planck)
      prad->sigma_planck(k,j,i,ifr) = (kappa_planck-kappaa)*rho*rhounit*lunit;
    else
      prad->sigma_planck(k,j,i,ifr) = 0.0;
  }
  }}}

}



void rossopacity(const Real rho, const Real tgas, Real &kappa, Real &kappa_planck)
{
  
    
    Real logt = log10(tgas * tunit);
    Real logrhot = log10(rho* rhounit) - 3.0* logt + 18.0;
    int nrhot1 = 0;
    int nrhot2 = 0;
    
    while((logrhot > logrhottable(nrhot2)) && (nrhot2 < 36)){
      nrhot1 = nrhot2;
      nrhot2++;
    }
    if(nrhot2==36 && (logrhot > logrhottable(nrhot2)))
      nrhot1=nrhot2;
  
  /* The data point should between NrhoT1 and NrhoT2 */
    int nt1 = 0;
    int nt2 = 0;
    while((logt > logttable(nt2)) && (nt2 < 137)){
      nt1 = nt2;
      nt2++;
    }
    if(nt2==137 && (logt > logttable(nt2)))
      nt1=nt2;
  

    Real kappa_t1_rho1=opacitytable(nt1,nrhot1);
    Real kappa_t1_rho2=opacitytable(nt1,nrhot2);
    Real kappa_t2_rho1=opacitytable(nt2,nrhot1);
    Real kappa_t2_rho2=opacitytable(nt2,nrhot2);

    Real planck_t1_rho1=planckopacity(nt1,nrhot1);
    Real planck_t1_rho2=planckopacity(nt1,nrhot2);
    Real planck_t2_rho1=planckopacity(nt2,nrhot1);
    Real planck_t2_rho2=planckopacity(nt2,nrhot2);


    Real rho_1 = logrhottable(nrhot1);
    Real rho_2 = logrhottable(nrhot2);
    Real t_1 = logttable(nt1);
    Real t_2 = logttable(nt2);

    
    if(nrhot1 == nrhot2){
      if(nt1 == nt2){
        kappa = kappa_t1_rho1;
        kappa_planck = planck_t1_rho1;
      }else{
        kappa = kappa_t1_rho1 + (kappa_t2_rho1 - kappa_t1_rho1) *
                                (logt - t_1)/(t_2 - t_1);
        kappa_planck = planck_t1_rho1 + (planck_t2_rho1 - planck_t1_rho1) *
                                (logt - t_1)/(t_2 - t_1);
      }/* end same T*/
    }else{
      if(nt1 == nt2){
        kappa = kappa_t1_rho1 + (kappa_t1_rho2 - kappa_t1_rho1) *
                                (logrhot - rho_1)/(rho_2 - rho_1);
        kappa_planck = planck_t1_rho1 + (planck_t1_rho2 - planck_t1_rho1) *
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
            Real miuz = prad->mu(0,k,j,ie+i,ifr*prad->nang+n);
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


