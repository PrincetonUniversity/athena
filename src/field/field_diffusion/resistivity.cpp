//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================

// Athena++ headers
#include "field_diffusion.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../field.hpp"

//----------------------------------------------------------------------------------------
//! \fn void FieldDiffusion::Resistivity
//  \brief Adds Ohmic dissipation EMF to the field EMF

void FieldDiffusion::Resistivity(const FaceField &bi,
    const AthenaArray<Real> &bc, EdgeField &emf)
{
  AthenaArray<Real> &e1=emf.x1e;
  AthenaArray<Real> &e2=emf.x2e;
  AthenaArray<Real> &e3=emf.x3e;
  AthenaArray<Real> &j1=j_.x1e;
  AthenaArray<Real> &j2=j_.x2e;
  AthenaArray<Real> &j3=j_.x3e;
  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;


// step-1. get current as curlB
  CurlB(bi,bc,j_); //calc j_ with bi and bcc

// step-2. calculate the EMF= eta_Ohmic*J
//------ 1D update:
  if (pmb_->block_size.nx2 == 1) {
    for (int i=is; i<=ie+1; ++i) {
      Real eta_O = 0.5*(Eta_Ohmic(ks,js,i-1)+Eta_Ohmic(ks,js,i));
      e2(ks  ,js  ,i) = eta_O*j2(ks,js,i);
      e2(ke+1,js  ,i) = e2(ks,js,i);
      e3(ks  ,js  ,i) = eta_O*j3(ks,js,i);
      e3(ks  ,je+1,i) = e3(ks,js,i);
    }
    return;
  }

//------ 2D update:
  if (pmb_->block_size.nx3 == 1) {
    for (int j=js; j<=je+1; ++j) {
      for (int i=is; i<=ie+1; ++i) {
        Real eta_O = 0.5*(Eta_Ohmic(ks,j,i)+Eta_Ohmic(ks,j-1,i));
        e1(ks  ,j,i) = eta_O*j1(ks,j,i);
        e1(ke+1,j,i) = e1(ks,j,i);

        eta_O = 0.5*(Eta_Ohmic(ks,j,i)+Eta_Ohmic(ks,j,i-1));
        e2(ks  ,j,i) = eta_O*j2(ks,j,i);
        e2(ke+1,j,i) = e2(ks,j,i);

        eta_O = 0.25*(Eta_Ohmic(ks,j,i)+Eta_Ohmic(ks,j,i-1)
                     +Eta_Ohmic(ks,j-1,i)+Eta_Ohmic(ks,j-1,i-1));
        e3(ks,  j,i) = eta_O*j3(ks,j,i);
      }
    }
    return;
  }

//------ 3D update:

  for (int k=ks; k<=ke+1; ++k) {
  for (int j=js; j<=je+1; ++j) {
    for (int i=is; i<=ie+1; ++i) {
      Real eta_O = 0.25*(Eta_Ohmic(k,j,i)+Eta_Ohmic(k,j-1,i)
                        +Eta_Ohmic(k-1,j,i)+Eta_Ohmic(k-1,j-1,i));

      e1(k,j,i) = eta_O*j1(k,j,i);

      eta_O = 0.25*(Eta_Ohmic(k,j,i)+Eta_Ohmic(k,j,i-1)
                   +Eta_Ohmic(k-1,j,i)+Eta_Ohmic(k-1,j,i-1));

      e2(k,j,i) = eta_O*j2(k,j,i);

      eta_O = 0.25*(Eta_Ohmic(k,j,i)+Eta_Ohmic(k,j,i-1)
                   +Eta_Ohmic(k,j-1,i)+Eta_Ohmic(k,j-1,i-1));

      e3(k,j,i) = eta_O*j3(k,j,i);
    }
  }}

  return;
}

//-------------------------------------------------------------------------------------
// Calculate current density
void FieldDiffusion::CurlB(const FaceField &bi, const AthenaArray<Real> &bc, EdgeField &crnt)
{

  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;

  AthenaArray<Real> &crnt1=crnt.x1e;
  AthenaArray<Real> &crnt2=crnt.x2e;
  AthenaArray<Real> &crnt3=crnt.x3e;
  const AthenaArray<Real> &bf1=bi.x1f;
  const AthenaArray<Real> &bf2=bi.x2f;
  const AthenaArray<Real> &bf3=bi.x3f;

  // 1-dim case
  if(pmb_->block_size.nx2==1) {
    for (int k=ks;k<=ke;k++) {
    for (int j=js;j<=je;j++) {
    for (int i=is-1;i<=ie+2;i++) {
      //crnt1(k,j,i) = 0.0;
      crnt2(k,j,i) = (pco_->h31v(i-1)*bc(IB3,k,j,i-1)
                     -pco_->h31v(i)*bc(IB3,k,j,i))/pco_->dx1v(i-1)/pco_->h31f(i);
      crnt3(k,j,i) = (pco_->h2v(i)*bc(IB2,k,j,i)
                     -pco_->h2v(i-1)*bc(IB2,k,j,i-1))/pco_->dx1v(i-1)/pco_->h2f(i);
    }}}
    return;
  }
  // 2-dim case
  if(pmb_->block_size.nx3==1) {
    for (int k=ks;k<=ke;k++) {
    for (int j=js-1;j<=je+2;j++) {
    for (int i=is-2;i<=ie+2;i++) {
      crnt1(k,j,i)=(pco_->h32v(j)*bc(IB3,k,j,i)
                   -pco_->h32v(j-1)*bc(IB3,k,j-1,i))
                   /pco_->dx2v(j-1)/(pco_->h2v(i)*pco_->h32f(j));
    }}}
    for (int k=ks;k<=ke;k++) {
    for (int j=js-2;j<=je+2;j++) {
    for (int i=is-1;i<=ie+2;i++) {
      crnt2(k,j,i)=(pco_->h31v(i-1)*bc(IB3,k,j,i-1)
                   -pco_->h31v(i)*bc(IB3,k,j,i))
                   /pco_->dx1v(i-1)/pco_->h31f(i);
    }}}
    for (int k=ks;k<=ke;k++) {
    for (int j=js-1;j<=je+2;j++) {
    for (int i=is-1;i<=ie+2;i++) {
      crnt3(k,j,i)=((pco_->h2v(i)*bf2(k,j,i)
                    -pco_->h2v(i-1)*bf2(k,j,i-1))/pco_->dx1v(i-1)-
                    (bf1(k,j,i)-bf1(k,j-1,i))/pco_->dx2v(j-1))/pco_->h2f(i);
    }}}
    return;
  }

  // 3-dim case
  for (int k=ks-1;k<=ke+2;k++){
  for (int j=js-1;j<=je+2;j++) {
  for (int i=is-2;i<=ie+2;i++) {
    const Real& h2h3=pco_->h2v(i)*pco_->h31v(i)*pco_->h32f(j);
    const Real& h2=pco_->h2v(i);
    const Real& h2p1=pco_->h2v(i);
    const Real& h3=pco_->h31v(i)*pco_->h32v(j-1);
    const Real& h3p1=pco_->h31v(i)*pco_->h32v(j);
    crnt1(k,j,i)=((h3p1*bf3(k,j,i)-h3*bf3(k,j-1,i))/pco_->dx2v(j-1) -
                  (h2p1*bf2(k,j,i)-h2*bf2(k-1,j,i))/pco_->dx3v(k-1))/h2h3;
  }}}
  for (int k=ks-1;k<=ke+2;k++){
  for (int j=js-2;j<=je+2;j++) {
  for (int i=is-1;i<=ie+2;i++) {
    const Real& h1h3=pco_->h31f(i)*pco_->h32v(j);
    const Real& h3=pco_->h31v(i-1)*pco_->h32v(j);
    const Real& h3p1=pco_->h31v(i)*pco_->h32v(j);
    crnt2(k,j,i)=((bf1(k,j,i)-bf1(k-1,j,i))/pco_->dx3v(k-1) -
                  (h3p1*bf3(k,j,i)-h3*bf3(k,j,i-1))/pco_->dx1v(i-1))/h1h3;
  }}}
  for (int k=ks-2;k<=ke+1;k++){
  for (int j=js-1;j<=je+2;j++) {
  for (int i=is-1;i<=ie+2;i++) {
    const Real& h1h2=pco_->h2f(i);
    const Real& h2=pco_->h2v(i-1);
    const Real& h2p1=pco_->h2v(i);
    crnt3(k,j,i)=((h2p1*bf2(k,j,i)-h2*bf2(k,j,i-1))/pco_->dx1v(i-1) -
                  (bf1(k,j,i)-bf1(k,j-1,i))/pco_->dx2v(j-1))/h1h2;
  }}}
  return;
}

//-------------------------------------------------------------------------------------
// Get the spatial variable eta_Ohmic
Real FieldDiffusion::Eta_Ohmic(const int k, const int j, const int i)
{
  if(ieta_==0) return etaO_;
  else {
    std::cout << "ieta_= " << ieta_ << " not implemented yet" << std::endl;
    return 0.0;
  }
}
