//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================

// Athena++ headers
#include "diffusion.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../hydro.hpp"
#include "../../eos/eos.hpp"

//----------------------------------------------------------------------------------------
//! \fn void HydroDiffusion::Viscosity
//  \brief Adds viscous stress as flux to the hydro flux

void HydroDiffusion::Viscosity(const AthenaArray<Real> &prim,
    const AthenaArray<Real> &cons, AthenaArray<Real> *diflx)
{
  AthenaArray<Real> &x1flux=diflx[X1DIR];
  AthenaArray<Real> &x2flux=diflx[X2DIR];
  AthenaArray<Real> &x3flux=diflx[X3DIR];
  int il, iu, jl, ju, kl, ku;
  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;
  Real denf;
  Real nuiso2 = -(2./3.);


  Divv(prim,divv_);

// step-2. calculate the flux across each face
// i-direction
// set the loop limits
  jl=js, ju=je, kl=ks, ku=ke;
  if (MAGNETIC_FIELDS_ENABLED) {
    if(pmb_->block_size.nx2 > 1) {
      if(pmb_->block_size.nx3 == 1) // 2D
        jl=js-1, ju=je+1, kl=ks, ku=ke;
      else // 3D
        jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
    }
  }
  for (int k=kl; k<=ku; ++k){
    for (int j=jl; j<=ju; ++j){
      // compute fluxes
      FaceXdx(k,j,is,ie+1,prim,fx_);
      FaceXdy(k,j,is,ie+1,prim,fy_);
      FaceXdz(k,j,is,ie+1,prim,fz_);
      // store fluxes
      for (int i=is; i<=ie+1; ++i){
        Real nu1 = nuiso1(prim,IM1,k,j,i);
        denf = 0.5*(prim(IDN,k,j,i)+prim(IDN,k,j,i-1));
        x1flux(IM1,k,j,i) = -denf*nu1*(fx_(i)+nuiso2*0.5*(divv_(k,j,i)+divv_(k,j,i-1)));
        x1flux(IM2,k,j,i) = -denf*nu1*fy_(i);
        x1flux(IM3,k,j,i) = -denf*nu1*fz_(i);
        if(NON_BAROTROPIC_EOS)
          x1flux(IEN,k,j,i) = 0.5*((prim(IM1,k,j,i-1)+prim(IM1,k,j,i))*x1flux(IM1,k,j,i) +
                                (prim(IM2,k,j,i-1)+prim(IM2,k,j,i))*x1flux(IM2,k,j,i) +
                                (prim(IM3,k,j,i-1)+prim(IM3,k,j,i))*x1flux(IM3,k,j,i));
      }
  }}

// j-direction
// set the loop limits
  il=is, iu=ie, kl=ks, ku=ke;
  if (MAGNETIC_FIELDS_ENABLED) {
    if(pmb_->block_size.nx3 == 1) // 2D
      il=is-1, iu=ie+1, kl=ks, ku=ke;
    else // 3D
      il=is-1, iu=ie+1, kl=ks-1, ku=ke+1;
  }
  if(pmb_->block_size.nx2 > 1) { //2D or 3D
    for (int k=kl; k<=ku; ++k){
      for (int j=js; j<=je+1; ++j){
        // compute fluxes
        FaceYdx(k,j,is,ie,prim,fx_);
        FaceYdy(k,j,is,ie,prim,fy_);
        FaceYdz(k,j,is,ie,prim,fz_);
        // store fluxes
        for(int i=il; i<=iu; i++) {
          Real nu1 = nuiso1(prim,IM1,k,j,i);
          denf = 0.5*(prim(IDN,k,j-1,i)+prim(IDN,k,j,i));
          x2flux(IM1,k,j,i) = -denf*nu1*fx_(i);
          x2flux(IM2,k,j,i) = -denf*nu1*
                           (fy_(i)+nuiso2*0.5*(divv_(k,j-1,i)+divv_(k,j,i)));
          x2flux(IM3,k,j,i) = -denf*nu1*fz_(i);
          if (NON_BAROTROPIC_EOS)
            x2flux(IEN,k,j,i) = 0.5*((prim(IM1,k,j,i)+prim(IM1,k,j-1,i))*x2flux(IM1,k,j,i) +
                                  (prim(IM2,k,j,i)+prim(IM2,k,j-1,i))*x2flux(IM2,k,j,i) +
                                  (prim(IM3,k,j,i)+prim(IM3,k,j-1,i))*x2flux(IM3,k,j,i));
        }
    }}
  } else { // 1D
    // compute fluxes
    FaceYdx(ks,js,is,ie,prim,fx_);
    FaceYdy(ks,js,is,ie,prim,fy_);
    FaceYdz(ks,js,is,ie,prim,fz_);
    // store fluxes
    for(int i=il; i<=iu; i++) {
      Real nu1 = nuiso1(prim,IM1,ks,js,i);
      denf = prim(IDN,ks,js,i);
      x2flux(IM1,ks,js,i) = -denf*nu1*fx_(i);
      x2flux(IM2,ks,js,i) = -denf*nu1*(fy_(i)+nuiso2*divv_(ks,js,i));
      x2flux(IM3,ks,js,i) = -denf*nu1*fz_(i);
      if (NON_BAROTROPIC_EOS)
        x2flux(IEN,ks,js,i) = prim(IM1,ks,js,i)*x2flux(IM1,ks,js,i) +
                              prim(IM2,ks,js,i)*x2flux(IM2,ks,js,i) +
                              prim(IM3,ks,js,i)*x2flux(IM3,ks,js,i);
    }
    for(int i=il; i<=iu; i++) {
      x2flux(IM1,ks,je+1,i) = x2flux(IM1,ks,js,i);
      x2flux(IM2,ks,je+1,i) = x2flux(IM2,ks,js,i);
      x2flux(IM3,ks,je+1,i) = x2flux(IM3,ks,js,i);
      if (NON_BAROTROPIC_EOS)
        x2flux(IEN,ks,je+1,i) =x2flux(IEN,ks,js,i);
    }
  }
// k-direction
// set the loop limits
  il=is, iu=ie, jl=js, ju=je;
  if (MAGNETIC_FIELDS_ENABLED) {
    if(pmb_->block_size.nx2 > 1) // 2D or 3D
      il=is-1, iu=ie+1, jl=js-1, ju=je+1;
    else // 1D
      il=is-1, iu=ie+1;
  }
  if(pmb_->block_size.nx3 > 1) { //3D
    for (int k=ks; k<=ke+1; ++k){
      for (int j=jl; j<=ju; ++j){
        // compute fluxes
        FaceZdx(k,j,is,ie,prim,fx_);
        FaceZdy(k,j,is,ie,prim,fy_);
        FaceZdz(k,j,is,ie,prim,fz_);
        // store fluxes
        for(int i=il; i<=iu; i++) {
          Real nu1 = nuiso1(prim,IM1,k,j,i);
          denf = 0.5*(prim(IDN,k-1,j,i)+prim(IDN,k,j,i));
          x3flux(IM1,k,j,i) = -denf*nu1*fx_(i);
          x3flux(IM2,k,j,i) = -denf*nu1*fy_(i);
          x3flux(IM3,k,j,i) = -denf*nu1*(fz_(i)+nuiso2*0.5*(divv_(k-1,j,i)+divv_(k,j,i)));
          if (NON_BAROTROPIC_EOS)
            x3flux(IEN,k,j,i) = 0.5*((prim(IM1,k,j,i)+prim(IM1,k-1,j,i))*x3flux(IM1,k,j,i) +
                                    (prim(IM2,k,j,i)+prim(IM2,k-1,j,i))*x3flux(IM2,k,j,i) +
                                    (prim(IM3,k,j,i)+prim(IM3,k-1,j,i))*x3flux(IM3,k,j,i));
        }
    }}
  } else {
    for (int j=jl; j<=ju; ++j){
      // compute fluxes
      FaceZdx(ks,j,is,ie,prim,fx_);
      FaceZdy(ks,j,is,ie,prim,fy_);
      FaceZdz(ks,j,is,ie,prim,fz_);
      // store fluxes
      for(int i=il; i<=iu; i++) {
        Real nu1 = nuiso1(prim,IM1,ks,j,i);
        denf = prim(IDN,ks,j,i);
        x3flux(IM1,ks,j,i) = -denf*nu1*fx_(i);
        x3flux(IM2,ks,j,i) = -denf*nu1*fy_(i);
        x3flux(IM3,ks,j,i) = -denf*nu1*(fz_(i)+nuiso2*divv_(ks,j,i));
        x3flux(IM1,ke+1,j,i) = x3flux(IM1,ks,j,i);
        x3flux(IM2,ke+1,j,i) = x3flux(IM2,ks,j,i);
        x3flux(IM3,ke+1,j,i) = x3flux(IM3,ks,j,i);
        if (NON_BAROTROPIC_EOS){
          x3flux(IEN,ks,j,i) = prim(IM1,ks,j,i)*x3flux(IM1,ks,j,i) +
                               prim(IM2,ks,j,i)*x3flux(IM2,ks,j,i) +
                               prim(IM3,ks,j,i)*x3flux(IM3,ks,j,i);
          x3flux(IEN,ke+1,j,i) = x3flux(IEN,ks,j,i);
        }
      }
    }}

  return;
}

//-------------------------------------------------------------------------------------
// Get the first coefficient nuiso1
Real HydroDiffusion::nuiso1(const AthenaArray<Real> &prim, const int n, const int k, const int j, const int i)
{
  if(inu_==0) return (nuiso_);
  else {
  // if inu != 0, prescribe viscosity according to nu = alpha cs^2/Omega
  // where alpha = nuiso_ and Omega = sqrt(GM/r^3) the Keplerian orb freq
  if(COORDINATE_SYSTEM == "cartesian" && NON_BAROTROPIC_EOS) {
    Real x1 = pmb_->pcoord->x1v(i);
    Real x2 = pmb_->pcoord->x2v(j);
    Real x3 = pmb_->pcoord->x3v(k);
    Real rad = sqrt(SQR(x1)+SQR(x2)); //+SQR(x3)+SQR(0.05));// softening by 0.05
    Real gamma = pmb_->peos->GetGamma();
    // 1) calc visc coef using updated prim variable (unstable)
    //Real cs2  = gamma*prim(IPR,k,j,i)/prim(IDN,k,j,i);
    //Real omg = 1.0/sqrt(SQR(rad)*rad); // Kepler freq with GM=1.0
    //Real nuprof = cs2/omg; //
    // 2) calc visc coef using not updated cons variable
    Real dens  = pmb_->phydro->u(IDN,k,j,i);
    Real m1 = pmb_->phydro->u(IM1,k,j,i);
    Real m2 = pmb_->phydro->u(IM2,k,j,i);
    Real m3 = pmb_->phydro->u(IM3,k,j,i);
    Real pres = (gamma-1.0)*(pmb_->phydro->u(IEN,k,j,i)-0.5*(SQR(m1)+SQR(m2)+SQR(m3))/dens);
    Real omg = 1.0/sqrt(SQR(rad)*rad); // Kepler freq with GM=1.0
    Real nuprof = gamma*pres/dens/omg;
    // 3) calc visc coef using pure power regardless of binary location
    //Real nuprof = gamma*0.01*sqrt(rad); // power law
    if (nuprof <= 0.0 || isnan(nuprof) || isinf(nuprof)) nuprof = SQR(0.10)/omg;
      return (nuiso_*nuprof);
  } else {
      std::cout << "inu_= " << inu_ << "for other coord and eos are not implemented yet" << std::endl;
    return 0.0;
  }
  }
}

//-------------------------------------------------------------------------------------
// Calculate Divv
void HydroDiffusion::Divv(const AthenaArray<Real> &prim, AthenaArray<Real> &divv)
{

  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;
  int il = is-1; int iu = ie+1;
  int jl, ju, kl, ku;
  Real area_p1, area;
  Real vel_p1, vel;

  if(pmb_->block_size.nx2 == 1) // 1D
    jl=js, ju=je, kl=ks, ku=ke;
  else if(pmb_->block_size.nx3 == 1) // 2D
    jl=js-1, ju=je+1, kl=ks, ku=ke;
  else // 3D
    jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;

  for (int k=kl; k<=ku; ++k){
    for (int j=jl; j<=ju; ++j){
      // calculate x1-flux divergence
      pmb_->pcoord->Face1Area(k,j,il,iu+1,x1area_);
      for (int i=il; i<=iu; ++i){
        area_p1 = x1area_(i+1);
        area    = x1area_(i);
        vel_p1  = 0.5*(prim(IM1,k,j,i+1)+prim(IM1,k,j,i));
        vel     = 0.5*(prim(IM1,k,j,i)+prim(IM1,k,j,i-1));
        divv(k,j,i) = area_p1*vel_p1 - area*vel;
      }
      // calculate x2-flux divergnece
      if (pmb_->block_size.nx2 > 1) {
        pmb_->pcoord->Face2Area(k,j  ,il,iu,x2area_);
        pmb_->pcoord->Face2Area(k,j+1,il,iu,x2area_p1_);
        for (int i=il; i<=iu; ++i){
          area_p1 = x2area_p1_(i);
          area    = x2area_(i);
          vel_p1  = 0.5*(prim(IM2,k,j+1,i)+prim(IM2,k,j,i));
          vel     = 0.5*(prim(IM2,k,j,i)+prim(IM2,k,j-1,i));
          divv(k,j,i) += area_p1*vel_p1 - area*vel;
        }
      }
      if (pmb_->block_size.nx3 > 1) {
        pmb_->pcoord->Face3Area(k  ,j,il,iu,x3area_);
        pmb_->pcoord->Face3Area(k+1,j,il,iu,x3area_p1_);
        for (int i=il; i<=iu; ++i){
          area_p1 = x3area_p1_(i);
          area    = x3area_(i);
          vel_p1  = 0.5*(prim(IM3,k+1,j,i)+prim(IM3,k,j,i));
          vel     = 0.5*(prim(IM3,k,j,i)+prim(IM3,k-1,j,i));
          divv(k,j,i) += area_p1*vel_p1 - area*vel;
        }
      }
      pmb_->pcoord->CellVolume(k,j,il,iu,vol_);
      for (int i=il; i<=iu; ++i){
        divv(k,j,i) = divv(k,j,i)/vol_(i);
      }
    }
  }

  return;
}

// v_{x1;x1}  covariant derivative at x1 interface
void HydroDiffusion::FaceXdx(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
  for (int i=il; i<=iu; ++i){
    len(i) = 2.0*(prim(IM1,k,j,i)-prim(IM1,k,j,i-1))/pco_->dx1v(i-1);
  }
  return;
}

// v_{x2;x1}+v_{x1;x2}  covariant derivative at x1 interface
void HydroDiffusion::FaceXdy(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
  if(pmb_->block_size.nx2 > 1) {
    for (int i=il; i<=iu; ++i){
      len(i) = pco_->h2f(i)*(prim(IM2,k,j,i)/pco_->h2v(i)
                      -prim(IM2,k,j,i-1)/pco_->h2v(i-1))/pco_->dx1v(i-1)
               +0.5*(prim(IM1,k,j+1,i)+prim(IM1,k,j+1,i-1)-prim(IM1,k,j-1,i)-prim(IM1,k,j-1,i-1))
                /pco_->h2f(i)/(pco_->dx2v(j-1)+pco_->dx2v(j));
    }
  } else {
    for (int i=il; i<=iu; ++i){
      len(i) = pco_->h2f(i)*(prim(IM2,k,j,i)/pco_->h2v(i)
                      -prim(IM2,k,j,i-1)/pco_->h2v(i-1))/pco_->dx1v(i-1);
    }
  }
  return;
}

// v_{x3;x1}+v_{x1;x3}  covariant derivative at x1 interface
void HydroDiffusion::FaceXdz(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
  if(pmb_->block_size.nx3 > 1) {
    for (int i=il; i<=iu; ++i){
      len(i) = pco_->h31f(i)*(prim(IM3,k,j,i)/pco_->h31v(i)-prim(IM3,k,j,i-1)/pco_->h31v(i-1))/pco_->dx1v(i-1)
               +0.5*(prim(IM1,k+1,j,i)+prim(IM1,k+1,j,i-1)-prim(IM1,k-1,j,i)-prim(IM1,k-1,j,i-1))
               /pco_->h31f(i)/pco_->h32v(j)/(pco_->dx3v(k-1)+pco_->dx3v(k));
    }
  } else {
    for (int i=il; i<=iu; ++i)
      len(i) = pco_->h31f(i)*(prim(IM3,k,j,i)/pco_->h31v(i)-prim(IM3,k,j,i-1)/pco_->h31v(i-1))/pco_->dx1v(i-1);

  }
  return;
}

// v_{x1;x2}+v_{x2;x1}  covariant derivative at x2 interface
void HydroDiffusion::FaceYdx(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
  if(pmb_->block_size.nx2 > 1) {
    for (int i=il; i<=iu; ++i){
      len(i) = (prim(IM1,k,j,i)-prim(IM1,k,j-1,i))/pco_->h2v(i)/pco_->dx2v(j-1)
               +pco_->h2v(i)*0.5*((prim(IM2,k,j,i+1)+prim(IM2,k,j-1,i+1))/pco_->h2v(i+1)
                           -(prim(IM2,k,j,i-1)+prim(IM2,k,j-1,i-1))/pco_->h2v(i-1))
               /(pco_->dx1v(i-1)+pco_->dx1v(i));
    }
  } else {
    for (int i=il; i<=iu; ++i){
      len(i) = pco_->h2v(i)*(prim(IM2,k,j,i+1)/pco_->h2v(i+1)
                            -prim(IM2,k,j,i-1)/pco_->h2v(i-1))
               /(pco_->dx1v(i-1)+pco_->dx1v(i));
    }
  }
  return;
}

// v_{x2;x2}  covariant derivative at x2 interface
void HydroDiffusion::FaceYdy(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
  if(pmb_->block_size.nx2 > 1) {
    for (int i=il; i<=iu; ++i){
      len(i) = 2.0*(prim(IM2,k,j,i)-prim(IM2,k,j-1,i))/pco_->h2v(i)/pco_->dx2v(j-1)
                +(prim(IM1,k,j,i)+prim(IM1,k,j-1,i))/pco_->h2v(i)*pco_->dh2vd1(i);
    }
  } else {
    for (int i=il; i<=iu; ++i)
      len(i) = 2.0*prim(IM1,k,j,i)/pco_->h2v(i)*pco_->dh2vd1(i);
  }
  return;
}

// v_{x3;x2}+v_{x2;x3}  covariant derivative at x2 interface
void HydroDiffusion::FaceYdz(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
  if(pmb_->block_size.nx3 > 1) {
    for (int i=il; i<=iu; ++i){
      len(i) = pco_->h32f(j)*(prim(IM3,k,j,i)/pco_->h32f(j)
               -prim(IM3,k,j-1,i)/pco_->h32f(j-1))/pco_->h2v(i)/pco_->dx2v(j-1)
               +0.5*(prim(IM2,k+1,j,i)+prim(IM2,k+1,j-1,i)
               -prim(IM2,k-1,j,i)-prim(IM2,k-1,j-1,i))/pco_->h31v(i)
               /pco_->h32f(j)/(pco_->dx3v(k-1)+pco_->dx3v(k));
    }
  } else if (pmb_->block_size.nx2 > 1) {
    for (int i=il; i<=iu; ++i){
      len(i) = pco_->h32f(j)*(prim(IM3,k,j,i)/pco_->h32f(j)-prim(IM3,k,j-1,i)/pco_->h32f(j-1))
               /pco_->h2v(i)/pco_->dx2v(j-1);
    }
  } else {
    for (int i=il; i<=iu; ++i) len(i) = 0.0;
  }
  return;
}

// v_{x1;x3}+v_{x3;x1}  covariant derivative at x3 interface
void HydroDiffusion::FaceZdx(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
  if(pmb_->block_size.nx3 > 1) {
    for (int i=il; i<=iu; ++i){
      len(i) = (prim(IM1,k,j,i)-prim(IM1,k-1,j,i))/pco_->dx3v(k-1)
               +0.5*pco_->h31v(i)*((prim(IM3,k,j,i+1)+prim(IM3,k-1,j,i+1))/pco_->h31v(i+1)
               -(prim(IM3,k,j,i-1)+prim(IM3,k-1,j,i-1))/pco_->h31v(i-1))
               /(pco_->dx1v(i-1)+pco_->dx1v(i));
    }
  } else {
    for (int i=il; i<=iu; ++i){
      len(i) = pco_->h31v(i)*(prim(IM3,k,j,i+1)/pco_->h31v(i+1)
               -prim(IM3,k,j,i-1)/pco_->h31v(i-1))
               /(pco_->dx1v(i-1)+pco_->dx1v(i));
    }
  }
  return;
}

// v_{x2;x3}+v_{x3;x2}  covariant derivative at x3 interface
void HydroDiffusion::FaceZdy(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
  if(pmb_->block_size.nx3 > 1) {
    for (int i=il; i<=iu; ++i){
      len(i) = (prim(IM2,k,j,i)-prim(IM2,k-1,j,i))/pco_->h31v(i)/pco_->h32v(j)/pco_->dx3v(k-1)
               +0.5*pco_->h32v(j)*((prim(IM3,k,j+1,i)+prim(IM3,k-1,j+1,i))/pco_->h32v(j+1)
                                  -(prim(IM3,k,j-1,i)+prim(IM3,k-1,j-1,i))/pco_->h32v(j-1))
               /pco_->h2v(i)/(pco_->dx2v(j-1)+pco_->dx2v(j));
    }
  } else if(pmb_->block_size.nx2 > 1) {
    for (int i=il; i<=iu; ++i){
      len(i) = pco_->h32v(j)*(prim(IM3,k,j+1,i)/pco_->h32v(j+1)
                            -prim(IM3,k,j-1,i)/pco_->h32v(j-1))
               /pco_->h2v(i)/(pco_->dx2v(j-1)+pco_->dx2v(j));
    }
  } else {
    for (int i=il; i<=iu; ++i) len(i) = 0.0;
  }
  return;
}

// v_{x3;x3}  covariant derivative at x3 interface
void HydroDiffusion::FaceZdz(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
  if(pmb_->block_size.nx3 > 1) {
    for (int i=il; i<=iu; ++i){
      len(i) = 2.0*(prim(IM3,k,j,i)-prim(IM3,k-1,j,i))/pco_->dx3v(k-1)/pco_->h31v(i)/pco_->h32v(j)
               + (prim(IM1,k,j,i)+prim(IM1,k-1,j,i))*pco_->dh31vd1(i)/pco_->h31v(i) +
                 (prim(IM2,k,j,i)+prim(IM2,k-1,j,i))*pco_->dh32vd2(j)/pco_->h32v(j);
    }
  } else {
    for (int i=il; i<=iu; ++i){
      len(i) = 2.0*prim(IM1,k,j,i)*pco_->dh31vd1(i)/pco_->h31v(i) +
               2.0*prim(IM2,k,j,i)*pco_->dh32vd2(j)/pco_->h32v(j);
    }
  }
  return;
}


