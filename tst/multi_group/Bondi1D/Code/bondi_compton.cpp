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

static Real lunit;
static Real rhounit;
static Real tunit;
static Real rhoinf;
static Real csinf;
static Real gm;
static Real mdot;
static Real rsonic;
static Real csonic;
static Real p_rho_coef; // this is coefficient K for P=K\rho^gamma
static Real bondi_gamma;
static Real fgam;
static Real rho_outbd;
static Real tg_outbd;

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

//void CentralGravity( MeshBlock *pmb, const Real time, const Real dt,
//  const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
//  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons, AthenaArray<Real> &cons_scalar);

void DiskOpacity(MeshBlock *pmb, AthenaArray<Real> &prim);


void BondiVel(double vel, double coef1, double coef2, double coef3,
    double gamma, double * fval, double *dfval);

double Rtsafe(void (*funcd)(double, double, double, double,double,double *,double *),
      double x1, double x2, double xacc,
      double coef1, double coef2, double coef3, double coef4);

void Tequilibrium(Real temperature, Real coef1, Real coef2, Real coef3,
                  Real coef4, Real *fval, Real *dfval);

// AMR refinement condition
int RefinementCondition(MeshBlock *pmb);

void Mesh::InitUserMeshData(ParameterInput *pin)
{

  if (adaptive)
    EnrollUserRefinementCondition(RefinementCondition);

//  Real gamma = my_blocks(0)->peos->GetGamma();
  bondi_gamma = pin->GetOrAddReal("problem", "bondi_gamma", 1.4);

  rhoinf = pin->GetOrAddReal("problem", "density_inf", 1.e-7);
  rho_outbd = rhoinf * 2.25606;


  fgam = pow(2.0/(5.0 - 3.0 * bondi_gamma), 
                  (5.0 - 3.0 * bondi_gamma)/(2*(bondi_gamma - 1.0)));

  csinf = sqrt(bondi_gamma);
  tg_outbd = 0.0185423; // fix outer boundary to 1.e4 K


  gm = csinf * csinf; 
  mdot = PI * csinf * rhoinf * fgam;

  csonic = csinf * sqrt(2.0/(5.0-3.0*bondi_gamma));
  rsonic = 0.5 * gm/(csonic*csonic);

  // cs^2 = gamma P/rho, P=rho cs^2/gamma
  p_rho_coef = csinf*csinf/(bondi_gamma * pow(rhoinf, bondi_gamma-1.0));

  rhounit = pin->GetOrAddReal("radiation", "density_unit", 1.e-10);
  tunit = pin->GetOrAddReal("radiation", "T_unit", 1.e5);
  lunit = pin->GetOrAddReal("radiation", "length_unit", 7.07766e16);

//  EnrollUserExplicitSourceFunction(CentralGravity);

  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, InjectHydro);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, VacuumHydro);

  if(RADIATION_ENABLED || IM_RADIATION_ENABLED){
    EnrollUserRadBoundaryFunction(BoundaryFace::inner_x1, Inject);
    EnrollUserRadBoundaryFunction(BoundaryFace::outer_x1, Vacuum);
  }
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  
  if(RADIATION_ENABLED || IM_RADIATION_ENABLED){
    
    prad->EnrollOpacityFunction(DiskOpacity);

  }





}


//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief beam test
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  
  //initialize random number
  std::srand(gid);

  
//  inject_tr = pin->GetOrAddReal("problem","inject_tr",1.0);
  Real prat = prad->prat;

  int nfreq = prad->nfreq;
  int nang = prad->nang;
  AthenaArray<Real> ir_cm;
  ir_cm.NewAthenaArray(prad->n_fre_ang);
  Real *ir_lab;
  
  Real gamma = peos->GetGamma();
  Real amp = 0.0;

  // solve the equation for density
  // coef1 * rho^(\gamma+1) + coef2 * rho^2 + coef3 = 0

      // electron scattering opacity
  Real kappas = 0.2 * (1.0 + 0.7) * rhounit * lunit;
  
  // Initialize hydro variable
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        Real radius = pcoord->x1v(i);

        Real coef1 =  0.5;
        Real coef2 = -(1/radius+1/(bondi_gamma-1.0));
        Real coef3 = (1/(bondi_gamma-1.0))*pow(0.25/(radius*radius),bondi_gamma-1.0)*
                    pow(2.0/(5-3*bondi_gamma),0.5*(5.0-3*bondi_gamma));

        Real vsonic = csonic/csinf;
        Real vfrefall = sqrt(2*(1./radius+1./(bondi_gamma-1.0)));

        Real vel = 0.0;


        if(radius > rsonic)
          vel = Rtsafe(BondiVel, 1.e-6, vsonic, 1.e-12, coef1, coef2, coef3, bondi_gamma);
        else
          vel = Rtsafe(BondiVel, vsonic, vfrefall, 1.e-12, coef1, coef2, coef3, bondi_gamma);


  
        vel *= -csinf;
        Real rho = -mdot/(4.0*PI*vel*radius*radius);
        Real cssq = csinf*csinf + gm*(bondi_gamma-1)/radius-0.5*vel*vel*(bondi_gamma-1.0);

        Real gast =  100*cssq/bondi_gamma;

        Real totp = rho * gast;


        coef1 = prat/3.0;
        coef2 = rho;
        coef3 = -totp;

        gast = Rtsafe(Tequilibrium, 0.0, gast, 1.e-12, coef1, coef2, coef3, 0.0);


        // assuming sum of gas and radiation pressure = P
        // assuming T_r=T_g, blackbody spectrum


        phydro->u(IDN,k,j,i) = rho;

        phydro->u(IM1,k,j,i) = rho * vel*10.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS){

          phydro->u(IEN,k,j,i) = rho * gast/(gamma-1.0);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        }
          
        for(int ifr=0; ifr<nfreq; ++ifr){
          Real bd_emission = 0.0;
          if(ifr == nfreq - 1){
            bd_emission = 1.0 - prad->FitBlackBody(prad->nu_grid(ifr)/gast);
          }else{
            bd_emission = prad->BlackBodySpec(prad->nu_grid(ifr)/gast,prad->nu_grid(ifr+1)/gast);
          }
          for(int n=0; n<nang; ++n){
            prad->ir(k,j,i,nang*ifr+n) = bd_emission * gast*gast*gast*gast;
          }// end angle
        }// end frequency

        // convert co-moving frame to lab frame
  //      ir_lab = &(prad->ir(k,j,i,0));
  //      prad->pradintegrator->ComToLabMultiGroup(vel, 0, 0, &(prad->mu(0,k,j,i,0)), 
  //               &(prad->mu(1,k,j,i,0)), &(prad->mu(2,k,j,i,0)),ir_cm, ir_lab);

        Real opacity_t = std::min(gast,0.1);

        Real kappa_ff = 91170.5 * pow(opacity_t,-0.5) *rho;

        for (int ifr=0; ifr<prad->nfreq; ++ifr){

          prad->sigma_s(k,j,i,ifr) = rho * kappas;
          if(ifr < prad->nfreq-1){
            Real &nu_cen = prad->nu_cen(ifr);
            prad->sigma_a(k,j,i,ifr) = rho * kappa_ff * pow(nu_cen,-3.0) * (1.0-exp(-nu_cen/gast));

          }else{
            Real &nu_cen = prad->nu_cen(prad->nfreq-2);
            prad->sigma_a(k,j,i,ifr) = rho * kappa_ff * pow(nu_cen,-3.0) * (1.0-exp(-nu_cen/gast));      
          }

          prad->sigma_ae(k,j,i,ifr) = prad->sigma_a(k,j,i,ifr);
          prad->sigma_planck(k,j,i,ifr) = prad->sigma_a(k,j,i,ifr);

        }  


      }
    }
  }
  
  //Now initialize opacity and specific intensity


  
  return;
}


void Inject(MeshBlock *pmb, Coordinates *pco, Radiation *prad, 
          const AthenaArray<Real> &w, FaceField &b, 
          AthenaArray<Real> &ir,
          Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{

  int nang=prad->nang;
  int nfreq=prad->nfreq;


  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        Real rlocal = pmb->pcoord->x1v(is-i);
        Real router = pmb->pcoord->x1v(is);
//        Real router2=pmb->pcoord->x1v(is-i+2);

        for(int ifr=0; ifr<nfreq; ++ifr){
          for(int n=0; n<nang; ++n){
//            Real gradient = (ir(k,j,is-i+2,nang*ifr+n) - ir(k,j,is-i+1,nang*ifr+n))/
//                            (router2-router);
//            ir(k,j,is-i,nang*ifr+n) = std::max(ir(k,j,is-i+1,nang*ifr+n) 
//                                - gradient * (router-rlocal), TINY_NUMBER); 
            ir(k,j,is-i,nang*ifr+n) = ir(k,j,is,nang*ifr+n) * router *router/(rlocal*rlocal);
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

    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=1; i<=ngh; ++i) {
          Real rin = pco->x1v(is);
          Real rlocal = pco->x1v(is-i);
          Real router = pco->x1v(is-i+1);
          Real router2= pco->x1v(is-i+2);
          Real gradient = (prim(IDN,k,j,is-i+2) - prim(IDN,k,j,is-i+1))/
                            (router2-router);
          Real mdot = prim(IVX,k,j,is) * rin * rin * prim(IDN,k,j,is);

          prim(IDN,k,j,is-i) = prim(IDN,k,j,is-i+1) 
                              - gradient * (router-rlocal); 
          prim(IVX,k,j,is-i) = mdot/(rlocal*rlocal*prim(IDN,k,j,is-i));

          prim(IVY,k,j,is-i) = prim(IVY,k,j,is);
          prim(IVZ,k,j,is-i) = prim(IVY,k,j,is);

          gradient = (prim(IPR,k,j,is-i+2) - prim(IPR,k,j,is-i+1))/
                            (router2-router);
          prim(IPR,k,j,is-i) = prim(IPR,k,j,is-i+1) 
                              - gradient * (router-rlocal); 


        }
      }
    }


}


void VacuumHydro(
    MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{

    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=1; i<=ngh; ++i) {

          Real rout = pco->x1v(ie);
          Real rlocal = pco->x1v(ie+i);
          Real router = pco->x1v(ie+i-1);
          Real router2= pco->x1v(ie+i-2);
          Real gradient = (prim(IDN,k,j,ie+i-2) - prim(IDN,k,j,ie+i-1))/
                            (router2-router);
          Real mdot = prim(IVX,k,j,ie) * rout * rout * prim(IDN,k,j,ie);

          prim(IDN,k,j,ie+i) = prim(IDN,k,j,ie+i-1) 
                              - gradient * (router-rlocal); 
          prim(IVX,k,j,ie+i) = mdot/(rlocal*rlocal*prim(IDN,k,j,ie+i));

          prim(IVY,k,j,ie+i) = prim(IVY,k,j,ie);
          prim(IVZ,k,j,ie+i) = prim(IVZ,k,j,ie);

          gradient = (prim(IPR,k,j,ie+i-2) - prim(IPR,k,j,ie+i-1))/
                            (router2-router);
          prim(IPR,k,j,ie+i) = prim(IPR,k,j,ie+i-1) 
                              - gradient * (router-rlocal); 

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

        Real rlocal = pmb->pcoord->x1v(ie+i);
        Real router = pmb->pcoord->x1v(ie);

        for(int ifr=0; ifr<nfreq; ++ifr){
          for(int n=0; n<nang; ++n){
            ir(k,j,ie+i,nang*ifr+n) = ir(k,j,ie,nang*ifr+n) * router *router/(rlocal*rlocal);

          }// end angle
        }// end frequency
      }
    }
  }// end k
 

  return;
}


// Add gravity -omega^2 x
/*
void CentralGravity( MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons, AthenaArray<Real> &cons_scalar)
{

  for(int k=pmb->ks; k<=pmb->ke; ++k){
    for(int j=pmb->js; j<=pmb->je; ++j){
      for(int i=pmb->is; i<=pmb->ie; ++i){
        Real rho = prim(IDN,k,j,i);
        Real radius = pmb->pcoord->x1v(i);
        

        Real src = - dt * rho * ;
        cons(IM1,k,j,i) += src;

        cons(IEN,k,j,i) += prim(IVX,k,j,i) * src;
        
      }
    }
  }

}
*/

// electron scattering
// free-free opacity
// 3.68e56 T^-1/2 \rho \nu^-3 (1-exp(-h\nu/kT))
// 3.68e56 T0^-1/2 \rho_0 \nu_0^-3, \nu_0=kT_0/h
//2883.07 T^-1/2 \rho \nu^-3(1-exp(-\nu/\nu_0))
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
  Real kappas = 0.2 * (1.0 + 0.7) * rhounit * lunit;

  
  for (int k=kl; k<=ku; ++k) {
  for (int j=jl; j<=ju; ++j) {
  for (int i=il; i<=iu; ++i) {
    Real rho  = prim(IDN,k,j,i);
    Real gast = prim(IEN,k,j,i)/rho;
    gast = std::max(gast,0.1);
    Real kappa_ff = 91170.5 * pow(gast,-0.5) *rho;

  for (int ifr=0; ifr<prad->nfreq; ++ifr){

    prad->sigma_s(k,j,i,ifr) = rho * kappas;
    if(ifr < prad->nfreq-1){
      Real &nu_cen = prad->nu_cen(ifr);
      prad->sigma_a(k,j,i,ifr) = rho * kappa_ff * pow(nu_cen,-3.0) * (1.0-exp(-nu_cen/gast));

    }else{
      Real &nu_cen = prad->nu_cen(prad->nfreq-2);
      prad->sigma_a(k,j,i,ifr) = rho * kappa_ff * pow(nu_cen,-3.0) * (1.0-exp(-nu_cen/gast));      
    }

    prad->sigma_ae(k,j,i,ifr) = prad->sigma_a(k,j,i,ifr);
    prad->sigma_planck(k,j,i,ifr) = prad->sigma_a(k,j,i,ifr);

  }    


 }}}

}// end opacity function

// refinement condition: density gradient
int RefinementCondition(MeshBlock *pmb) 
{

  AthenaArray<Real> &w = pmb->phydro->w;
  Real maxeps=0.0;
  for(int k=pmb->ks; k<=pmb->ke; k++){
    for(int j=pmb->js; j<=pmb->je; j++) {
      for(int i=pmb->is; i<=pmb->ie; i++) {
//        Real eps= (std::abs(w(IDN,k,j,i+1)-w(IDN,k,j,i-1)))/w(IDN,k,j,i);
//        if(pmb->je > pmb->js)
//          eps += (std::abs(w(IDN,k,j+1,i)-w(IDN,k,j-1,i)))/w(IDN,k,j,i);
        maxeps = std::max(maxeps, w(IDN,k,j,i));
      }// end i
    }// end j
  }// end k
  Real x1max = pmb->pcoord->x1v(pmb->ie);
  if(maxeps < 50) return -1;
  if(maxeps > 200 ) return 1;
  return 0;

}

/* Function to find the density for Bondi solution */
void BondiVel(double vel, double coef1, double coef2, double coef3,
    double gamma, double * fval, double *dfval)
{

/* function is
 *  coef1 * vel^(gamma+1) + coef2 * vel^(gamma-1) + coef3 == 0 *
 */

  *fval = coef1 * pow(vel, gamma + 1.0) + coef2 * pow(vel,gamma-1.0) + coef3;
  *dfval = (gamma + 1.0) * coef1 * pow(vel, gamma) + (gamma-1.0) * coef2 * pow(vel,gamma-2.0);

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

