//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//  \brief Class to implement diffusion processes in the hydro equations

// Athena++ headers
#include <cfloat>     // FLT_MAX
#include "hydro_diffusion.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../hydro.hpp"
#include "../../parameter_input.hpp"
#include "../../field/field.hpp"

// HydroDiffusion constructor

HydroDiffusion::HydroDiffusion(Hydro *phyd, ParameterInput *pin)
{
  pmy_hydro_ = phyd;
  pmb_ = pmy_hydro_->pmy_block;
  pco_ = pmb_->pcoord;
  hydro_diffusion_defined = false;

  int nthreads = pmb_->pmy_mesh->GetNumMeshThreads();
  int ncells1 = pmb_->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmb_->block_size.nx2 > 1) ncells2 = pmb_->block_size.nx2 + 2*(NGHOST);
  if (pmb_->block_size.nx3 > 1) ncells3 = pmb_->block_size.nx3 + 2*(NGHOST);

  // Check if viscous process.
  coeff_nuiso = pin->GetOrAddReal("problem","nuiso",0.0); // iso viscosity
  coeff_nuani = pin->GetOrAddReal("problem","nuani",0.0); // aniso viscosity
  if (coeff_nuiso > 0.0 || coeff_kani  > 0.0) {
    hydro_diffusion_defined = true;
    // Allocate memory for fluxes.
    visflx[X1DIR].NewAthenaArray(NHYDRO,ncells3,ncells2,ncells1+1);
    visflx[X2DIR].NewAthenaArray(NHYDRO,ncells3,ncells2+1,ncells1);
    visflx[X3DIR].NewAthenaArray(NHYDRO,ncells3+1,ncells2,ncells1);
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

    nu.NewAthenaArray(2,ncells3,ncells2,ncells1);
    if(pmb_->pmy_mesh->ViscosityCoeff_==NULL)
      CalcViscCoeff_ = ConstViscosity;
    else
      CalcViscCoeff_ = pmb_->pmy_mesh->ViscosityCoeff_;
  }

  // Check if thermal conduction.
  if (NON_BAROTROPIC_EOS) {
    coeff_kiso  = pin->GetOrAddReal("problem","kiso",0.0); // iso thermal conduction
    coeff_kani  = pin->GetOrAddReal("problem","kani",0.0); // aniso thermal conduction
    if (coeff_kiso > 0.0 || coeff_kani > 0.0) {
      hydro_diffusion_defined = true;
      cndflx[X1DIR].NewAthenaArray(ncells3,ncells2,ncells1+1);
      cndflx[X2DIR].NewAthenaArray(ncells3,ncells2+1,ncells1);
      cndflx[X3DIR].NewAthenaArray(ncells3+1,ncells2,ncells1);

      kappa.NewAthenaArray(2,ncells3,ncells2,ncells1);
      if(pmb_->pmy_mesh->ConductionCoeff_==NULL)
        CalcCondCoeff_ = ConstConduction;
      else
        CalcCondCoeff_ = pmb_->pmy_mesh->ConductionCoeff_;
    }
  }

  if (hydro_diffusion_defined) {
    dx1_.NewAthenaArray(nthreads,ncells1);
    dx2_.NewAthenaArray(nthreads,ncells1);
    dx3_.NewAthenaArray(nthreads,ncells1);
    nu_tot_.NewAthenaArray(nthreads,ncells1);
    kappa_tot_.NewAthenaArray(nthreads,ncells1);
  }
}

// destructor

HydroDiffusion::~HydroDiffusion()
{
  if (coeff_nuiso > 0.0 || coeff_nuani > 0.0) {
    visflx[X1DIR].DeleteAthenaArray();
    visflx[X2DIR].DeleteAthenaArray();
    visflx[X3DIR].DeleteAthenaArray();
    x1area_.DeleteAthenaArray();
    x2area_.DeleteAthenaArray();
    x3area_.DeleteAthenaArray();
    x2area_p1_.DeleteAthenaArray();
    x3area_p1_.DeleteAthenaArray();
    vol_.DeleteAthenaArray();
    fx_.DeleteAthenaArray();
    fy_.DeleteAthenaArray();
    fz_.DeleteAthenaArray();
    divv_.DeleteAthenaArray();
  }
  if (coeff_kiso  > 0.0 || coeff_kani > 0.0) {
    cndflx[X1DIR].DeleteAthenaArray();
    cndflx[X2DIR].DeleteAthenaArray();
    cndflx[X3DIR].DeleteAthenaArray();
  }
  if (hydro_diffusion_defined) {
    dx1_.DeleteAthenaArray();
    dx2_.DeleteAthenaArray();
    dx3_.DeleteAthenaArray();
    nu_tot_.DeleteAthenaArray();
    kappa_tot_.DeleteAthenaArray();
  }
}

//----------------------------------------------------------------------------------------
//! \fn void HydroDiffusion::CalcHydroDiffusionFlux
//  \brief Calculate diffusion flux for hydro flux

void HydroDiffusion::CalcHydroDiffusionFlux(const AthenaArray<Real> &prim,
     const AthenaArray<Real> &cons, AthenaArray<Real> *flux)
{

  Hydro *ph=pmb_->phydro;
  Field *pf=pmb_->pfield;

  SetHydroDiffusivity(ph->w,pf->bcc);

  if (coeff_nuiso > 0.0 || coeff_nuani > 0.0) ClearHydroFlux(visflx);
  if (coeff_nuiso > 0.0) ViscousFlux_iso(prim, cons, visflx);
  if (coeff_nuani > 0.0) ViscousFlux_aniso(prim, cons, visflx);

  if (coeff_kiso > 0.0 || coeff_kani > 0.0) ClearHydroFlux(cndflx);
  if (coeff_kiso > 0.0) ThermalFlux_iso(prim, cons, cndflx);
  if (coeff_kani > 0.0) ThermalFlux_aniso(prim, cons, cndflx);


  return;
}

//----------------------------------------------------------------------------------------
//! \fn void HydroDiffusion::AddHydroDiffusionEnergyFlux
//  \brief Adds only diffusion energy flux to hydro flux

void HydroDiffusion::AddHydroDiffusionEnergyFlux(AthenaArray<Real> *flux_src, AthenaArray<Real> *flux_des)
{

  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;

  AthenaArray<Real> &x1flux=flux_des[X1DIR];
  AthenaArray<Real> &x2flux=flux_des[X2DIR];
  AthenaArray<Real> &x3flux=flux_des[X3DIR];
  AthenaArray<Real> &x1diflx=flux_src[X1DIR];
  AthenaArray<Real> &x2diflx=flux_src[X2DIR];
  AthenaArray<Real> &x3diflx=flux_src[X3DIR];

  for (int k=ks; k<=ke; ++k){
  for (int j=js; j<=je; ++j){
  for (int i=is; i<=ie; ++i){
    x1flux(IEN,k,j,i) += x1diflx(k,j,i);
    if(i==ie) x1flux(IEN,k,j,i+1) += x1diflx(k,j,i+1);
    if (pmb_->block_size.nx2 > 1) {
      x2flux(IEN,k,j,i) += x2diflx(k,j,i);
      if(j==je) x2flux(IEN,k,j+1,i) += x2diflx(k,j+1,i);
	  }
    if (pmb_->block_size.nx3 > 1) {
      x3flux(IEN,k,j,i) += x3diflx(k,j,i);
      if(k==ke) x3flux(IEN,k+1,j,i) += x3diflx(k+1,j,i);
	  }
  }}}

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void HydroDiffusion::AddHydroDiffusionFlux
//  \brief Adds all componenets of diffusion flux to hydro flux

void HydroDiffusion::AddHydroDiffusionFlux(AthenaArray<Real> *flux_src, AthenaArray<Real> *flux_des)
{

  int size1 = flux_des[X1DIR].GetSize();
#pragma simd
  for (int i=0; i<size1; ++i)
    flux_des[X1DIR](i) += flux_src[X1DIR](i);

  if (pmb_->block_size.nx2 > 1) {
    int size2 = flux_des[X2DIR].GetSize();
#pragma simd
    for (int i=0; i<size2; ++i)
      flux_des[X2DIR](i) += flux_src[X2DIR](i);
  }
  if (pmb_->block_size.nx3 > 1) {
    int size3 = flux_des[X3DIR].GetSize();
#pragma simd
    for (int i=0; i<size3; ++i)
      flux_des[X3DIR](i) += flux_src[X3DIR](i);
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void HydroDiffusion::ClearHydroDiffusionFlux
//  \brief Reset diffusion flux back to zeros

void HydroDiffusion::ClearHydroFlux(AthenaArray<Real> *flux)
{
  int size1 = flux[X1DIR].GetSize();
  int size2 = flux[X2DIR].GetSize();
  int size3 = flux[X3DIR].GetSize();

#pragma simd
  for (int i=0; i<size1; ++i)
    flux[X1DIR](i) = 0.0;

#pragma simd
  for (int i=0; i<size2; ++i)
    flux[X2DIR](i) = 0.0;

#pragma simd
  for (int i=0; i<size3; ++i)
    flux[X3DIR](i) = 0.0;

  return;
}

//---------------------------------------------------------------------------------------
//! \fn void HydroDiffusion::SetHydroDiffusivity(Athena<Real> &w, AthenaArray<Real> &bc)
// Set hydro diffusion coefficients; called in CALC_DIFFUSIVITY in tasklist

void HydroDiffusion::SetHydroDiffusivity(AthenaArray<Real> &w, AthenaArray<Real> &bc)
{
  int il = pmb_->is-NGHOST; int jl = pmb_->js; int kl = pmb_->ks;
  int iu = pmb_->ie+NGHOST; int ju = pmb_->je; int ku = pmb_->ke;
  if (pmb_->block_size.nx2 > 1) {
    jl -= NGHOST; ju += NGHOST;
  }
  if (pmb_->block_size.nx3 > 1) {
    kl -= NGHOST; ku += NGHOST;
  }

  // set viscosity using func ptr
  if (coeff_nuiso !=0.0 || coeff_nuani !=0.0)
    CalcViscCoeff_(this, w, bc, il, iu, jl, ju, kl, ku);
  // set thermal conduction using func ptr
  if (coeff_kiso !=0.0 || coeff_kani !=0.0)
    CalcCondCoeff_(this, w, bc, il, iu, jl, ju, kl, ku);

  return;
}

// Get the hydro diffusion timestep
// currently return dt for viscous and conduction processes
void HydroDiffusion::NewHydroDiffusionDt(Real &dt_vis, Real &dt_cnd)
{
  int tid=0;
  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;

  Real fac;
  if(pmb_->block_size.nx3>1) fac = 1.0/6.0;
  else if(pmb_->block_size.nx2>1) fac = 0.25;
  else fac = 0.5;

  int nthreads = pmb_->pmy_mesh->GetNumMeshThreads();
  Real *ptd_mindt_vis;
  Real *ptd_mindt_cnd;
  ptd_mindt_vis = new Real [nthreads];
  ptd_mindt_cnd = new Real [nthreads];

  for (int n=0; n<nthreads; ++n) {
    ptd_mindt_vis[n] = (FLT_MAX);
    ptd_mindt_cnd[n] = (FLT_MAX);
  }

#pragma omp parallel default(shared) private(tid) num_threads(nthreads)
{
#ifdef OPENMP_PARALLEL
  tid=omp_get_thread_num();
#endif
  AthenaArray<Real> nu_t;
  nu_t.InitWithShallowSlice(nu_tot_,2,tid,1);
  AthenaArray<Real> kappa_t;
  kappa_t.InitWithShallowSlice(kappa_tot_,2,tid,1);
  AthenaArray<Real> len, dx2, dx3;
  len.InitWithShallowSlice(dx1_,2,tid,1);
  dx2.InitWithShallowSlice(dx2_,2,tid,1);
  dx3.InitWithShallowSlice(dx3_,2,tid,1);

  for (int k=ks; k<=ke; ++k){
#pragma omp for schedule(static)
    for (int j=js; j<=je; ++j){
#pragma simd
      for (int i=is; i<=ie; ++i){
        nu_t(i) = 0.0;
        kappa_t(i) = 0.0;
      }
      if (coeff_nuiso > 0.0){
#pragma simd
        for (int i=is; i<=ie; ++i) nu_t(i) += nu(ISO,k,j,i);
      }
      if (coeff_nuani > 0.0){
#pragma simd
        for (int i=is; i<=ie; ++i) nu_t(i) += nu(ANI,k,j,i);
      }
      if (coeff_kiso > 0.0){
#pragma simd
        for (int i=is; i<=ie; ++i) kappa_t(i) += kappa(ISO,k,j,i);
      }
      if (coeff_kani > 0.0){
#pragma simd
        for (int i=is; i<=ie; ++i) kappa_t(i) += kappa(ANI,k,j,i);
      }
      pmb_->pcoord->CenterWidth1(k,j,is,ie,len);
      pmb_->pcoord->CenterWidth2(k,j,is,ie,dx2);
      pmb_->pcoord->CenterWidth3(k,j,is,ie,dx3);
#pragma simd
      for (int i=is; i<=ie; ++i){
        len(i) = (pmb_->block_size.nx2 > 1) ? std::min(len(i),dx2(i)):len(i);
        len(i) = (pmb_->block_size.nx3 > 1) ? std::min(len(i),dx3(i)):len(i);
      }
      if ((coeff_nuiso > 0.0) || (coeff_nuani > 0.0)) {
#pragma simd
        for (int i=is; i<=ie; ++i)
          ptd_mindt_vis[tid] = std::min(ptd_mindt_vis[tid], SQR(len(i))*fac/(nu_t(i)+TINY_NUMBER));
      }
      if ((coeff_kiso > 0.0) || (coeff_kani > 0.0)) {
#pragma simd
        for (int i=is; i<=ie; ++i)
          ptd_mindt_cnd[tid]= std::min(ptd_mindt_cnd[tid], SQR(len(i))*fac/(kappa_t(i)+TINY_NUMBER));
      }
    }
  }
} // end of omp parallel region

  dt_vis = ptd_mindt_vis[0];
  dt_cnd = ptd_mindt_cnd[0];
  for (int n=1; n<nthreads; ++n) {
    dt_vis = std::min(dt_vis,ptd_mindt_vis[n]);
    dt_cnd = std::min(dt_cnd,ptd_mindt_cnd[n]);
  }

  delete[] ptd_mindt_vis;
  delete[] ptd_mindt_cnd;

  return;
}
