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
static Real kappaes = 0.34; //electron scattering opacity in cgs unit
static Real rho0 = 50.0; // Density normalization of the stream
static Real inib0 = 0.0; // initial magnetic field strength
static Real tfloor; // temperature floor
static Real rhofloor; // density floor
static int bconf = 0; // bconf=1: pure B_phi
                      // bconf=0: vector potential proportional to density
                      // bconf=2: two loops
static Real gm1 = 6.48607e8;
static Real gm2 = 6.48607e7;
static Real qm = gm2/gm1;
static Real omega0 = 0.732146;
static Real rm1 = 1100;

//the initial thickness of the injected stream
static Real wheight=20.0;
//the initial injected radial velocity
// set vr_l1=-0.5 * v_bin
static Real vr_l1=-539.542*0.01;
static Real vphi_l1= 539.542 * 0.1;
static Real t_l1=1.;
static Real beta=1.0;
static Real amp=1.e-3;


static Real lunit;
static Real rhounit;
static Real tunit;


// the opacity table

static AthenaArray<Real> opacitytable;
static AthenaArray<Real> planckopacity;
static AthenaArray<Real> logttable;
static AthenaArray<Real> logrhottable;

int RefinementCondition(MeshBlock *pmb);

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
     const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir, 
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void Inflow_rad_X1(MeshBlock *pmb, Coordinates *pco, Radiation *prad, 
     const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir, 
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void DiskOpacity(MeshBlock *pmb, AthenaArray<Real> &prim);

void rossopacity(const Real rho, const Real tgas, Real &kappa, Real &kappa_planck);

void TidalPotential(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,  AthenaArray<Real> &cons_scalar);

void Tequilibrium(Real temperature, Real coef1, Real coef2, Real coef3,
                  Real coef4, Real *fval, Real *dfval);

double Rtsafe(void (*funcd)(double, double, double, double,double,double *,double *),
      double x1, double x2, double xacc,
      double coef1, double coef2, double coef3, double coef4);




void Mesh::InitUserMeshData(ParameterInput *pin)
{


  if (adaptive)
    EnrollUserRefinementCondition(RefinementCondition);  
  
    // Enroll boundary functions

  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, Outflow_X2);
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, Inflow_X1);
  
  tfloor = pin->GetOrAddReal("radiation", "tfloor", 0.001);
  rhofloor = pin->GetOrAddReal("hydro", "dfloor", 1.e-7);
  rhounit = pin->GetOrAddReal("radiation", "rhounit", 1.e-9);
  tunit = pin->GetOrAddReal("radiation", "Tunit", 1.e5);
  lunit = pin->GetOrAddReal("radiation", "lunit", 1.48428e11);
  
  EnrollUserExplicitSourceFunction(TidalPotential);
  

  if(RADIATION_ENABLED){
  
    EnrollUserRadBoundaryFunction(BoundaryFace::inner_x1, Inflow_rad_X1);
    EnrollUserRadBoundaryFunction(BoundaryFace::outer_x1, Outflow_rad_X2);
  
    // the opacity table
    opacitytable.NewAthenaArray(138,37);
    planckopacity.NewAthenaArray(138,37);
    logttable.NewAthenaArray(138);
    logrhottable.NewAthenaArray(37);
    
    // read in the opacity table
    FILE *fkappa, *flogt, *flogrhot, *fplanck;
      
    if ( (fkappa=fopen("./aveopacity.txt","r"))==NULL )
    {
      printf("Open input file error");
      return;
    }

    if ( (fplanck=fopen("./PlanckOpacity.txt","r"))==NULL )
    {
      printf("Open input file error");
      return;
    }

    if ( (flogt=fopen("./logT.txt","r"))==NULL )
    {
      printf("Open input file error");
      return;
    }

    if ( (flogrhot=fopen("./logRhoT.txt","r"))==NULL )
    {
      printf("Open input file error");
      return;
    }


    int i, j;
    for(j=0; j<138; j++){
      for(i=0; i<37; i++){
          fscanf(fkappa,"%lf",&(opacitytable(j,i)));
      }
    }

    for(int j=0; j<138; j++){
      for(int i=0; i<37; i++){
          fscanf(fplanck,"%lf",&(planckopacity(j,i)));
      }
     }


    for(i=0; i<37; i++){
      fscanf(flogrhot,"%lf",&(logrhottable(i)));
    }

    for(i=0; i<138; i++){
      fscanf(flogt,"%lf",&(logttable(i)));
    }

    fclose(fkappa);
    fclose(flogt);
    fclose(flogrhot);
    fclose(fplanck);

 
  }
  

  return;
}



//======================================================================================
//! \fn void Mesh::TerminateUserMeshProperties(void)
//  \brief Clean up the Mesh properties
//======================================================================================
void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{

  // free memory
  if(RADIATION_ENABLED){

    opacitytable.DeleteAthenaArray();
    logttable.DeleteAthenaArray();
    logrhottable.DeleteAthenaArray();
    planckopacity.DeleteAthenaArray();

  }


  return;
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  
  
  if(RADIATION_ENABLED){
    
      prad->EnrollOpacityFunction(DiskOpacity);

  }else{

  }
  

  return;
}


void MeshBlock::UserWorkInLoop(void)
{


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
      Real valimit = 0.9 * 8053.61;

  for (int k=kl; k<=ku; ++k){
    for (int j=jl; j<=ju; ++j){
     for (int i=il; i<=iu; ++i){

        int flag  = 0;

        Real radius = pcoord->x1v(i);

        Real& vx=phydro->w(IVX,k,j,i);
        Real& vy=phydro->w(IVY,k,j,i);
        Real& vz=phydro->w(IVZ,k,j,i);

        Real& rho=phydro->w(IDN,k,j,i);


        Real pb=0.0;
        if(MAGNETIC_FIELDS_ENABLED){
            pb = 0.5*(SQR(pfield->bcc(IB1,k,j,i))+SQR(pfield->bcc(IB2,k,j,i))
                   +SQR(pfield->bcc(IB3,k,j,i)));
        }
        Real Va = sqrt(2.0*pb/rho);

        if(Va > valimit){
          phydro->w(IDN,k,j,i) *= (Va/valimit)*(Va/valimit);
          phydro->u(IDN,k,j,i) = phydro->w(IDN,k,j,i);
          flag = 1;
        }


        if(rho < rhofloor){
          rho = rhofloor;
          phydro->u(IDN,k,j,i) = rhofloor;
          flag = 1;
        }



        Real vel = sqrt(vx*vx+vy*vy+vz*vz);

        if(vel > valimit){
          Real ratio = valimit / vel;
          vx *= ratio;
          vy *= ratio;
          vz *= ratio;
        

          phydro->u(IM1,k,j,i) = rho*vx;
          phydro->u(IM2,k,j,i) = rho*vy;
          phydro->u(IM3,k,j,i) = rho*vz;
          flag = 1;
        }

        Real ke = 0.5 * rho * (vx*vx+vy*vy+vz*vz);
        Real tgas = phydro->w(IPR,k,j,i)/rho;
        if(tgas < 0.01)
        {
          // calculate tr
          Real sum_er=0.0;
          for(int ifr=0; ifr<prad->nfreq; ++ifr){
            for(int n=0; n<prad->nang; ++n){
              sum_er += prad->ir(k,j,i,ifr*prad->nang+n) * prad->wmu(n);
            }
          }
          Real tr = pow(sum_er,0.25);
          if(tgas < tr){
            tgas = tr;
            phydro->w(IPR,k,j,i) = tgas * rho;
            flag = 1;
          }


        }
//        if(tgas > valimit * valimit){
//          flag = 1;
//          phydro->w(IPR,k,j,i) = rho * valimit * valimit;
//        }
        // fix a locally isothermal profile
//        tgas = 1.081e5/radius;


        if(NON_BAROTROPIC_EOS){
          if(flag > 0){
            Real eint = phydro->w(IPR,k,j,i)/gamma1;
            phydro->u(IEN,k,j,i) = eint + ke + pb;
          }

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
  Real gm = pin->GetReal("problem", "GM");

  //initialize random number
  

  AthenaArray<Real> ir_cm;
  Real *ir_lab;
  
  Real crat, prat;
  if(RADIATION_ENABLED){
    ir_cm.NewAthenaArray(prad->n_fre_ang);
    crat = prad->crat;
    prat = prad->prat;
  }else{
    crat = 8053.61;
    prat = 0.0;
  }
  


  
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
        phydro->u(IEN,k,j,i) = 0.1*rhofloor/(gamma-1.0);
        phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
        phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/phydro->u(IDN,k,j,i);
        phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        
        // initialize radiation quantity
        if(RADIATION_ENABLED){
          for(int n=0; n<prad->n_fre_ang; ++n)
             prad->ir(k,j,i,n) = 0.0;
        
        }// End Rad
        
      }// i
    }// j
  }// k

  // Opacity will be set during initialization

  if(RADIATION_ENABLED)
    ir_cm.DeleteAthenaArray();
  
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
  
  
  }


      
     // Update total energy with mangefew
    if(RADIATION_ENABLED){
     // Get cell-centered magnetic field
     pfield->CalculateCellCenteredField(pfield->b,pfield->bcc,pcoord,is,ie,js,je,ks,ke);
    
    
      for(int k=ks; k<=ke; ++k){
      for(int j=js; j<=je; ++j){
      for(int i=is; i<=ie; ++i){
          phydro->u(IEN,k,j,i) +=
            0.5*(SQR((pfield->bcc(IB1,k,j,i)))
               + SQR((pfield->bcc(IB2,k,j,i)))
               + SQR((pfield->bcc(IB3,k,j,i))));
        
         }
      }}
      
    }

 
  }// End MHD

  
  return;
}




void DiskOpacity(MeshBlock *pmb, AthenaArray<Real> &prim)
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
  Real kappas = 0.2 * (1.0 + 0.7);
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

    prad->sigma_p(k,j,i,ifr) = kappa_planck*rho*rhounit*lunit;
    prad->sigma_pe(k,j,i,ifr) = prad->sigma_p(k,j,i,ifr);

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

    // in the case the temperature is out of range
    // the planck opacity should be smaller by the 
    // ratio T^-3.5
    if(nt2 == 137 && (logt > logttable(nt2))){
       Real scaling = pow(10.0, -3.5*(logt - logttable(137)));
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
  Real Bheight = 0.8* wheight;
  
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
            prim(IEN,k,j,ie+i) = prim(IEN,k,j,ie);
          }else{
            prim(IDN,k,j,ie+i) = rhostream;
            prim(IVX,k,j,ie+i) = vr_l1;
            prim(IVY,k,j,ie+i) = 0.0;
            prim(IVZ,k,j,ie+i) = vphi_l1;
            prim(IEN,k,j,ie+i) = t_l1 * rhostream;
                
             // add random perturbation
            if(pmb->pmy_mesh->time < 1.e-12){
                prim(IDN,k,j,ie+i) *= (1.0 + amp *
                    ((double)rand()/(double)RAND_MAX-0.5));                
            }
                
          }// end if
      }//i
    }// j
  }// k
   // set magnetic field in inlet ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    Real a0=0.4*sqrt(2.0*rho0*t_l1)*sqrt(1.e3);
  
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
//          b.x1f(k,j,ie+i+1) = aphi/(tan(theta)*radius)
//                            - a0*sin(0.5*PI*dis/Bheight)*radius
//                            *(1.0-cos(theta)*cos(phi))/dis;
          b.x1f(k,j,ie+i+1) = aphi/(tan(theta)*radius)+aphi*2.0*xpos*cos(theta)*cos(phi)/(wheight*wheight); 
          if(theta < 0.5*PI)
            b.x1f(k,j,ie+i+1) = -fabs(b.x1f(k,j,ie+i+1));
          else
            b.x1f(k,j,ie+i+1) = fabs(b.x1f(k,j,ie+i+1));

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
          b.x2f(k,j,ie+i)=-aphi/radius+aphi*(2*radius-2.0*xpos*sin(theta)*cos(phi))/(wheight*wheight);
//          b.x2f(k,j,ie+i) =a0*sin(0.5*PI*dis/Bheight)*2.0*(1.0-sin(theta)*cos(phi))*radius/dis
//                          - aphi/radius;
          b.x2f(k,j,ie+i)=fabs(b.x2f(k,j,ie+i));   
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
             b.x3f(k,j,ie+i)=0.0*sqrt(2.0*rho0*t_l1/beta);
          }
      }
    }}
    
  }
  
/*  if(MAGNETIC_FIELDS_ENABLED){
    AthenaArray<Real> baphi;
    baphi.NewAthenaArray(ke-ks+1,je-js+1,NGHOST);
  }
*/
  return;
}


void Outflow_rad_X2(MeshBlock *pmb, Coordinates *pco, Radiation *prad, 
     const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir, 
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
  AthenaArray<Real> ir_cm;
  Real *ir_lab;
  
  ir_cm.NewAthenaArray(prad->n_fre_ang);
  
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
         Real radius=pco->x1v(ie+i);
         Real theta=pco->x2v(j);
         Real phi=pco->x3v(k);
         Real xpos=radius*sin(theta)*cos(phi);
         Real ypos=radius*sin(theta)*sin(phi);
         Real zpos=radius*cos(theta);
         // stream injected at (-radius, 0, 0)
         Real dis=ypos*ypos+(xpos+radius)*(xpos+radius)+zpos*zpos;
         Real rhostream=rho0*exp(-dis/(wheight*wheight));

         if(rhostream > rhofloor){
           
           for(int n=0; n<prad->n_fre_ang; ++n)
             ir_cm(n) = t_l1*t_l1*t_l1*t_l1;
           
           Real vr=w(IVX,k,j,ie+i);
           Real vtheta=w(IVY,k,j,ie+i);
           Real vphi=w(IVZ,k,j,ie+i);
           
           Real *mux = &(prad->mu(0,k,j,ie+i,0));
           Real *muy = &(prad->mu(1,k,j,ie+i,0));
           Real *muz = &(prad->mu(2,k,j,ie+i,0));

           ir_lab = &(ir(k,j,ie+i,0));
          
           prad->pradintegrator->ComToLab(vr,vtheta,vphi,mux,muy,muz,ir_cm,ir_lab);
              
         }else{
           for(int ifr=0; ifr<prad->nfreq; ++ifr){
              for(int n=0; n<prad->nang; ++n){
                Real miux = prad->mu(0,k,j,ie+i,ifr*prad->nang+n);
                if(miux > 0.0){
                  ir(k,j,ie+i,ifr*prad->nang+n)
                                = ir(k,j,ie+i-1,ifr*prad->nang+n);
                }else{
                  ir(k,j,ie+i,ifr*prad->nang+n) = 0.0;
                }
         
            }
         
           }
         }
        
      }//i
    }//j
  }//k
  
  ir_cm.DeleteAthenaArray();

  return;
}



void Inflow_rad_X1(MeshBlock *pmb, Coordinates *pco, Radiation *prad, 
     const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir, 
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
  AthenaArray<Real> ir_cm;
  Real *ir_lab;
  
  ir_cm.NewAthenaArray(prad->n_fre_ang);
  
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
         Real radius=pco->x1v(is-i);
         Real theta=pco->x2v(j);
         Real phi=pco->x3v(k);
         Real xpos=radius*sin(theta)*cos(phi);
         Real ypos=radius*sin(theta)*sin(phi);
         Real zpos=radius*cos(theta);
         // stream injected at (-radius, 0, 0)
         Real dis=ypos*ypos+(xpos+radius)*(xpos+radius)+zpos*zpos;
         Real rhostream=rho0*exp(-dis/(wheight*wheight));


         for(int ifr=0; ifr<prad->nfreq; ++ifr){
            for(int n=0; n<prad->nang; ++n){
              Real miux = prad->mu(0,k,j,is-i,ifr*prad->nang+n);
              if(miux < 0.0){
                ir(k,j,is-i,ifr*prad->nang+n)
                              = ir(k,j,is-i+1,ifr*prad->nang+n);
              }else{
                ir(k,j,is-i,ifr*prad->nang+n) = 0.0;
              }
         
            }
         
           }
     
        
      }//i
    }//j
  }//k
  
  ir_cm.DeleteAthenaArray();
  return;
}




//

Real grav_pot(const Real radius, const Real theta, const Real phi)
{
  // the companion is located at \theta=90, phi=0, r=rm1
  //x=rm1, y=0, z=0
  // current point r\sin\theta \cosphi, r\sin\theta\sin phi, r\cos\theta
  Real dist_r1=sqrt(radius*radius+rm1*rm1-2.0*radius*rm1*sin(theta)*cos(phi));
  Real rs_m2 = 2.0; //Schwarzschild radius for m2
  Real rs_m1 = rs_m2/qm; // mass ratio qm = M2/M1
  Real potphi=-gm2/(radius-rs_m2)-gm1/(dist_r1-rs_m1)-0.5*omega0*omega0*radius*radius
        *sin(theta)*sin(theta)+gm1*radius*sin(theta)*cos(phi)/(rm1*rm1);
  return potphi;
}


void TidalPotential(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,  AthenaArray<Real> &cons_scalar)
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


/*
void MeshBlock::InitAfterRestart(ParameterInput *pin)
{ 

  Real gamma1 = peos->GetGamma() - 1.0;
  Real prat = prad->prat;
  AthenaArray<Real> ir_cm;
  Real *ir_lab;
  ir_cm.NewAthenaArray(prad->n_fre_ang);
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

  if (MAGNETIC_FIELDS_ENABLED) {

    int nx1 = ie-is+1+2*NGHOST;
    int nx2 = 1;
    if(je > js) nx2 = je-js+1+2*NGHOST;
    int nx3 = 1;
    if(ke > ks) nx3 = ke-ks+1+2*NGHOST;

    pfield->CalculateCellCenteredField(pfield->b,pfield->bcc,pcoord,is,ie,js,je,ks,ke);

    peos->ConservedToPrimitive(phydro->u, phydro->w1, pfield->b,
                                        phydro->w, pfield->bcc, pcoord,
                                        il, iu, jl, ju, kl, ku);

  }


  for(int k=ks; k<=ke; ++k)
    for(int j=js; j<=je; ++j)
      for(int i=is; i<=ie; ++i){
        Real rho = phydro->w(IDN,k,j,i);
        Real press0 = phydro->w(IPR,k,j,i);
        Real temp0 = press0/rho;
        Real coef1 = prat/3.0;
        Real coef2 = rho;
        Real coef3 = -press0;
        Real gast = press0/rho;
        if(gast > 1.e-1)
          gast = Rtsafe(Tequilibrium, 0.0, temp0, 1.e-12, coef1, coef2, coef3, 0.0);

        Real pb = 0.0;
        if(MAGNETIC_FIELDS_ENABLED)
          pb = 0.5 *(pfield->bcc(IB1,k,j,i) * pfield->bcc(IB1,k,j,i)
                    +pfield->bcc(IB2,k,j,i) * pfield->bcc(IB2,k,j,i)
                    +pfield->bcc(IB3,k,j,i) * pfield->bcc(IB3,k,j,i));

        Real ekin = 0.5 * (phydro->u(IM1,k,j,i) * phydro->u(IM1,k,j,i)
                          +phydro->u(IM2,k,j,i) * phydro->u(IM2,k,j,i)
                          +phydro->u(IM3,k,j,i) * phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        phydro->u(IEN,k,j,i) = ekin + rho * gast/gamma1 + pb;

        for(int n=0; n<prad->n_fre_ang; ++n){
          ir_cm(n) = gast * gast *gast * gast;
        }

        Real *mux = &(prad->mu(0,k,j,i,0));
        Real *muy = &(prad->mu(1,k,j,i,0));
        Real *muz = &(prad->mu(2,k,j,i,0));

        ir_lab = &(prad->ir(k,j,i,0));
        prad->pradintegrator->ComToLab(phydro->w(IVX,k,j,i),phydro->w(IVY,k,j,i),phydro->w(IVZ,k,j,i),mux,muy,muz,ir_cm,ir_lab);

      }


      peos->ConservedToPrimitive(phydro->u, phydro->w1, pfield->b,
                                        phydro->w, pfield->bcc, pcoord,
                                        is, ie, js, je, ks, ke);


}

*/

void Tequilibrium(Real temperature, Real coef1, Real coef2, Real coef3,
      Real coef4, Real *fval, Real *dfval)
{
  // function is
  //  coef1 * T^4 + coef2 * T + coef3 == 0
  Real temp3 = temperature * temperature * temperature;
  Real temp4 = temp3 * temperature;


  *fval = coef1 * temp4 + coef2 * temperature + coef3;
  *dfval = 4.0 * coef1 * temp3 + coef2;

  return;

}


double Rtsafe(void (*funcd)(double, double, double, double,double,double *,double *),
      double x1, double x2, double xacc,
      double coef1, double coef2, double coef3, double coef4)
{
  int j;
  double df,dx,dxold,f,fh,fl;
  double temp,xh,xl,rts;
  std::stringstream msg;

  int maxit = 400;

  (*funcd)(x1,coef1, coef2, coef3,coef4, &fl,&df);
  (*funcd)(x2,coef1, coef2, coef3,coef4, &fh,&df);
  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)){
    std::cout << "[rtsafe]:Root must be bracketed in rtsafe: Tl:" << x1 << "Th: "
    << x2 << "\n fl: " << fl << "\n fh: " << fh << "coef1: " << coef1 << "coef2: "
    << coef2 <<  "coef3: " << coef3 << "coef4: " << coef4 << "\n" << std::endl;

    msg << "### FATAL ERROR in function [rtsafe]" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (fl == 0.0) return x1;
  if (fh == 0.0) return x2;
  if (fl < 0.0) {
    xl=x1;
    xh=x2;
  } else {
    xh=x1;
    xl=x2;
  }
  rts=0.5*(x1+x2);
  dxold=fabs(x2-x1);
  dx=dxold;
  (*funcd)(rts,coef1, coef2, coef3,coef4, &f,&df);
  for (j=1;j<=maxit;j++) {
    if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0)
        || (fabs(2.0*f) > fabs(dxold*df))) {
      dxold=dx;
      dx=0.5*(xh-xl);
      rts=xl+dx;
      if (xl == rts) return rts;
    } else {
      dxold=dx;
      dx=f/df;
      temp=rts;
      rts -= dx;
      if (temp == rts) return rts;
    }
    if (fabs(dx) < xacc) return rts;
    (*funcd)(rts,coef1, coef2, coef3,coef4, &f,&df);
    if (f < 0.0)
      xl=rts;
    else
      xh=rts;
  }

  std::cout << "[rtsafe]:Maximum number of iterations exceeded in rtsafe: Tl:"
    << x1 << "Th: "
    << x2 << "\n fl: " << fl << "\n fh: " << fh << "coef1: " << coef1 << "coef2: "
    << coef2 <<  "coef3: " << coef3 << "coef4: " << coef4 << "\n" << std::endl;

  msg << "### FATAL ERROR in function [rtsafe]" << std::endl;
    throw std::runtime_error(msg.str().c_str());


  return 0.0;
}


            

// refinement condition: density curvature
int RefinementCondition(MeshBlock *pmb) {

  Coordinates *pco = pmb->pcoord;

  Real xmin = pco->x1f(pmb->is);
  Real xmax = pco->x1f(pmb->ie+1);
  Real ymin = pco->x2f(pmb->js);
  Real ymax = pco->x2f(pmb->je+1);
  Real zmin = pco->x3f(pmb->ks);
  Real zmax = pco->x3f(pmb->ke+1);
  // refine : delta rho > 0.9*amp
  if (((xmin >= 8 && xmin <= 80) ||
     (xmax >=8 && xmax <= 80)) && (
     (ymin >= 1.3962634016 && ymin <= 1.7453292519) ||
     (ymax >= 1.3962634016 && ymax <= 1.7453292519) ))
    return 1;
  return 0;
}



