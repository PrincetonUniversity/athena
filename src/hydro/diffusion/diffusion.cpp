//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//  \brief Class to implement diffusion processes in the hydro equations

// Athena++ headers
#include "diffusion.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../hydro.hpp"
#include "../../parameter_input.hpp"

// HydroDiffusion constructor

HydroDiffusion::HydroDiffusion(Hydro *phyd, ParameterInput *pin)
{
  pmy_hydro_ = phyd;
  pmb_ = pmy_hydro_->pmy_block;
  pco_ = pmb_->pcoord;
  hydro_diffusion_defined = false;

  // read parameters such as viscosity from input block
  nuiso_ = pin->GetOrAddReal("problem","nuiso",0.0);
  if (nuiso_ != 0.0) {
    hydro_diffusion_defined = true;
    // Allocate memory for fluxes
    int ncells1 = pmb_->block_size.nx1 + 2*(NGHOST);
    int ncells2 = 1, ncells3 = 1;
    if (pmb_->block_size.nx2 > 1) ncells2 = pmb_->block_size.nx2 + 2*(NGHOST);
    if (pmb_->block_size.nx3 > 1) ncells3 = pmb_->block_size.nx3 + 2*(NGHOST);
    diflx[X1DIR].NewAthenaArray(NHYDRO,ncells3,ncells2,ncells1+1);
    diflx[X2DIR].NewAthenaArray(NHYDRO,ncells3,ncells2+1,ncells1);
    diflx[X3DIR].NewAthenaArray(NHYDRO,ncells3+1,ncells2,ncells1);

    x1area_.NewAthenaArray(ncells1+1);
    x2area_.NewAthenaArray(ncells1);
    x3area_.NewAthenaArray(ncells1);
    x2area_p1_.NewAthenaArray(ncells1);
    x3area_p1_.NewAthenaArray(ncells1);
    vol_.NewAthenaArray(ncells1);
    fx_.NewAthenaArray(ncells1);
    fy_.NewAthenaArray(ncells1);
    fz_.NewAthenaArray(ncells1);
    divv_.NewAthenaArray(ncells3,ncells2,ncells1);
  } else {
    std::cout << "[HydroDiffusion]: unable to construct hydrodiffusion class" << std::endl;
  }
}

// destructor

HydroDiffusion::~HydroDiffusion()
{
  if (nuiso_ != 0.0){
    diflx[X1DIR].DeleteAthenaArray();
    diflx[X2DIR].DeleteAthenaArray();
    diflx[X3DIR].DeleteAthenaArray();
    x1area_.DeleteAthenaArray();
    x2area_.DeleteAthenaArray();
    x3area_.DeleteAthenaArray();
    x2area_p1_.DeleteAthenaArray();
    x3area_p1_.DeleteAthenaArray();
    vol_.DeleteAthenaArray();
    fx_.DeleteAthenaArray();
    fy_.DeleteAthenaArray();
    fz_.DeleteAthenaArray();
  }
}

//----------------------------------------------------------------------------------------
//! \fn void HydroDiffusion::AddHydroDiffusionFlux
//  \brief Adds diffusion flux to hydro flux

void HydroDiffusion::AddHydroDiffusionFlux(const AthenaArray<Real> &prim,
     const AthenaArray<Real> &cons, AthenaArray<Real> *flux)
{

  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;

  AthenaArray<Real> &x1flux=flux[X1DIR];
  AthenaArray<Real> &x2flux=flux[X2DIR];
  AthenaArray<Real> &x3flux=flux[X3DIR];
  AthenaArray<Real> &x1diflx=diflx[X1DIR];
  AthenaArray<Real> &x2diflx=diflx[X2DIR];
  AthenaArray<Real> &x3diflx=diflx[X3DIR];

  // isotropic viscosity: to calc and add the flux due to viscous stress
  if (nuiso_ != 0.0) Viscosity(prim, cons, diflx);

  for (int n=0; n<(NHYDRO); ++n){
  for (int k=ks; k<=ke; ++k){
  for (int j=js; j<=je; ++j){
  for (int i=is; i<=ie; ++i){
    x1flux(n,k,j,i) += x1diflx(n,k,j,i);
    x2flux(n,k,j,i) += x2diflx(n,k,j,i);
    x3flux(n,k,j,i) += x3diflx(n,k,j,i);
    if(i==ie) x1flux(n,k,j,i+1) += x1diflx(n,k,j,i+1);
    if(j==je) x2flux(n,k,j+1,i) += x2diflx(n,k,j+1,i);
    if(k==ke) x3flux(n,k+1,j,i) += x3diflx(n,k+1,j,i);
  }}}}

  // template to add other diffusion processes
  //if (input_parameter_ != 0.0) OtherDiffusionProcess(dt, flux, prim,cons);

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void HydroDiffusion::NewDtDiff(Real len, int k, int j, int i)
//  \brief return the time step constraints due to explicit diffusion processes
Real HydroDiffusion::NewDtDiff(Real len, int k, int j, int i)
{
  Real diff_dt;
  if(pmb_->block_size.nx3>1){
    diff_dt = SQR(len)/6.0/nuiso_;
  } else {
    if(pmb_->block_size.nx2>1)
      diff_dt = SQR(len)/8.0/nuiso_;
    else
      diff_dt = SQR(len)/4.0/nuiso_;
  }
  return diff_dt;
}

