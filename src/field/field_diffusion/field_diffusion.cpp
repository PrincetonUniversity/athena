//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//  \brief Class to implement diffusion processes in the induction equations

// Athena++ headers
#include "field_diffusion.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../field.hpp"
#include "../../parameter_input.hpp"

// FieldDiffusion constructor

FieldDiffusion::FieldDiffusion(Field *pfld, ParameterInput *pin)
{
  pmy_field_ = pfld;
  pmb_ = pmy_field_->pmy_block;
  //pmy_hydro_ = pmb_->phydro;
  pco_ = pmb_->pcoord;
  field_diffusion_defined = false;

  // read parameters such as resistivity from input block
  etaO_ = pin->GetOrAddReal("problem","eta_O",0.0); //ohmic dissipation
  ieta_ = pin->GetOrAddInteger("problem","ieta",0); // default constant eta_O

  if (etaO_ != 0.0) {
    field_diffusion_defined = true;
    // Allocate memory for EMF
    int ncells1 = pmb_->block_size.nx1 + 2*(NGHOST);
    int ncells2 = 1, ncells3 = 1;
    if (pmb_->block_size.nx2 > 1) ncells2 = pmb_->block_size.nx2 + 2*(NGHOST);
    if (pmb_->block_size.nx3 > 1) ncells3 = pmb_->block_size.nx3 + 2*(NGHOST);
    emf.x1e.NewAthenaArray((ncells3+1),(ncells2+1),ncells1    );
    emf.x2e.NewAthenaArray((ncells3+1), ncells2   ,(ncells1+1));
    emf.x3e.NewAthenaArray( ncells3   ,(ncells2+1),(ncells1+1));

    j_.x1e.NewAthenaArray((ncells3+1),(ncells2+1),ncells1);
    j_.x2e.NewAthenaArray((ncells3+1),ncells2,(ncells1+1));
    j_.x3e.NewAthenaArray(ncells3,(ncells2+1),(ncells1+1));
  }
  //else {
  //  std::cout << "[FieldDiffusion]: unable to construct fielddiffusion class" << std::endl;
  //}
}

// destructor

FieldDiffusion::~FieldDiffusion()
{
  if (etaO_ != 0.0){
    emf.x1e.DeleteAthenaArray();
    emf.x2e.DeleteAthenaArray();
    emf.x3e.DeleteAthenaArray();
    j_.x1e.DeleteAthenaArray();
    j_.x2e.DeleteAthenaArray();
    j_.x3e.DeleteAthenaArray();
  }
}

//----------------------------------------------------------------------------------------
//! \fn void FieldDiffusion::CalcFieldDiffusionEMF
//  \brief Calculate diffusion EMF

void FieldDiffusion::CalcFieldDiffusionEMF(const FaceField &bi,
     const AthenaArray<Real> &bc, EdgeField &e)
{

  //int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  //int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;

  //AthenaArray<Real> &e1=e.x1e;
  //AthenaArray<Real> &e2=e.x2e;
  //AthenaArray<Real> &e3=e.x3e;
  //AthenaArray<Real> &emf1=emf.x1e;
  //AthenaArray<Real> &emf2=emf.x2e;
  //AthenaArray<Real> &emf3=emf.x3e;

  // Ohmic dissipation: to calc and add the EMF =\eta_Ohmic*J
  if (etaO_ != 0.0) {
    Resistivity(bi, bc, emf);

  //  for (int k=ks; k<=ke+1; ++k){
  //  for (int j=js; j<=je+1; ++j){
  //  for (int i=is; i<=ie; ++i){
  //    e1(k,j,i) += emf1(k,j,i);
  //  }}}

  //  for (int k=ks; k<=ke+1; ++k){
  //  for (int j=js; j<=je; ++j){
  //  for (int i=is; i<=ie+1; ++i){
  //    e2(k,j,i) += emf2(k,j,i);
  //  }}}

  //  for (int k=ks; k<=ke; ++k){
  //  for (int j=js; j<=je+1; ++j){
  //  for (int i=is; i<=ie+1; ++i){
  //    e3(k,j,i) += emf3(k,j,i);
  //  }}}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FieldDiffusion::AddFieldDiffusionEMF
//  \brief Add diffusion EMF

void FieldDiffusion::AddFieldDiffusionEMF(EdgeField &e)
{

  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;

  AthenaArray<Real> &e1=e.x1e;
  AthenaArray<Real> &e2=e.x2e;
  AthenaArray<Real> &e3=e.x3e;
  AthenaArray<Real> &emf1=emf.x1e;
  AthenaArray<Real> &emf2=emf.x2e;
  AthenaArray<Real> &emf3=emf.x3e;

  for (int k=ks; k<=ke+1; ++k){
  for (int j=js; j<=je+1; ++j){
  for (int i=is; i<=ie; ++i){
    e1(k,j,i) += emf1(k,j,i);
  }}}

  for (int k=ks; k<=ke+1; ++k){
  for (int j=js; j<=je; ++j){
  for (int i=is; i<=ie+1; ++i){
    e2(k,j,i) += emf2(k,j,i);
  }}}

  for (int k=ks; k<=ke; ++k){
  for (int j=js; j<=je+1; ++j){
  for (int i=is; i<=ie+1; ++i){
    e3(k,j,i) += emf3(k,j,i);
  }}}

  return;
}
//----------------------------------------------------------------------------------------
//! \fn void FieldDiffusion::AddEnergyFlux
//  \brief Adds diffusion flux (due to field dissipation) to hydro energy flux

void FieldDiffusion::AddEnergyFlux(const AthenaArray<Real> &bc, AthenaArray<Real> *flux)
{

  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;

  AthenaArray<Real> &x1flux=flux[X1DIR];
  AthenaArray<Real> &x2flux=flux[X2DIR];
  AthenaArray<Real> &x3flux=flux[X3DIR];
  //AthenaArray<Real> &e1=pmb_->pfield->e.x1e;
  //AthenaArray<Real> &e2=pmb_->pfield->e.x2e;
  //AthenaArray<Real> &e3=pmb_->pfield->e.x3e;
  AthenaArray<Real> &e1=emf.x1e;
  AthenaArray<Real> &e2=emf.x2e;
  AthenaArray<Real> &e3=emf.x3e;

  Real exb;
  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
  for (int i=is; i<=ie; ++i) {
    // x1flux
    exb = 0.25*(e2(k,j,i)+e2(k+1,j,i))*(bc(IB3,k,j,i)+bc(IB3,k,j,i-1))-
          0.25*(e3(k,j,i)+e3(k,j+1,i))*(bc(IB2,k,j,i)+bc(IB2,k,j,i-1));
    x1flux(IEN,k,j,i) += exb;
    if(i==ie) {
      exb = 0.25*(e2(k,j,i+1)+e2(k+1,j,i+1))*(bc(IB3,k,j,i+1)+bc(IB3,k,j,i))-
            0.25*(e3(k,j,i+1)+e3(k,j+1,i+1))*(bc(IB2,k,j,i+1)+bc(IB2,k,j,i));
      x1flux(IEN,k,j,i+1) += exb;
    }
    // x2flux
    if(pmb_->block_size.nx2 > 1) {
      exb = 0.25*(e3(k,j,i)+e3(k,j,i+1))*(bc(IB1,k,j,i)+bc(IB1,k,j-1,i))-
            0.25*(e1(k,j,i)+e1(k+1,j,i))*(bc(IB3,k,j,i)+bc(IB3,k,j-1,i));
      x2flux(IEN,k,j,i) += exb;
      if(j==je) {
        exb = 0.25*(e3(k,j+1,i)+e3(k,j+1,i+1))*(bc(IB1,k,j+1,i)+bc(IB1,k,j,i))-
              0.25*(e1(k,j+1,i)+e1(k+1,j+1,i))*(bc(IB3,k,j+1,i)+bc(IB3,k,j,i));
        x2flux(IEN,k,j+1,i) += exb;
	    }
    }
    // x3flux
    if(pmb_->block_size.nx3 > 1) {
      exb = 0.25*(e1(k,j,i)+e1(k,j+1,i))*(bc(IB2,k,j,i)+bc(IB2,k-1,j,i))-
            0.25*(e2(k,j,i)+e2(k,j,i+1))*(bc(IB1,k,j,i)+bc(IB1,k-1,j,i));
      x3flux(IEN,k,j,i) += exb;
      if(k==ke) {
        exb = 0.25*(e1(k+1,j,i)+e1(k+1,j+1,i))*(bc(IB2,k+1,j,i)+bc(IB2,k,j,i))-
              0.25*(e2(k+1,j,i)+e2(k+1,j,i+1))*(bc(IB1,k+1,j,i)+bc(IB1,k,j,i));
        x3flux(IEN,k+1,j,i) += exb;
      }
    }
  }}}
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FieldDiffusion::NewDtFldDiff(Real len, int k, int j, int i)
//  \brief return the time step constraints due to explicit diffusion processes
Real FieldDiffusion::NewDtFldDiff(Real len, int k, int j, int i)
{
  Real diff_dt;
  if(pmb_->block_size.nx3>1){
    diff_dt = SQR(len)/6.0/etaO_;
  } else {
    if(pmb_->block_size.nx2>1)
      diff_dt = SQR(len)/8.0/etaO_;
    else
      diff_dt = SQR(len)/4.0/etaO_;
  }
  return diff_dt;
}

