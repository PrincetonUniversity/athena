//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hydro_diffusion.cpp
//! \brief Class to implement diffusion processes in the hydro equations

// C headers

// C++ headers
#include <algorithm>  // min()
#include <cstring>    // strcmp()
#include <iostream>   // endl
#include <limits>
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../field/field.hpp"
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"
#include "../hydro.hpp"
#include "hydro_diffusion.hpp"

//! HydroDiffusion constructor

HydroDiffusion::HydroDiffusion(Hydro *phyd, ParameterInput *pin) :
    hydro_diffusion_defined(false),
    nu_iso{pin->GetOrAddReal("problem", "nu_iso", 0.0)},
    nu_aniso{pin->GetOrAddReal("problem", "nu_aniso", 0.0)},
    kappa_iso{}, kappa_aniso{},
    pmy_hydro_(phyd), pmb_(pmy_hydro_->pmy_block), pco_(pmb_->pcoord) {
  int nc1 = pmb_->ncells1, nc2 = pmb_->ncells2, nc3 = pmb_->ncells3;

  // Check if viscous process are active
  if (nu_iso > 0.0 || nu_aniso  > 0.0) {
    hydro_diffusion_defined = true;
    // Allocate memory for fluxes
    visflx[X1DIR].NewAthenaArray(NHYDRO, nc3, nc2, nc1+1);
    visflx[X2DIR].NewAthenaArray(NHYDRO, nc3, nc2+1, nc1);
    visflx[X3DIR].NewAthenaArray(NHYDRO, nc3+1, nc2, nc1);
    x1area_.NewAthenaArray(nc1+1);
    x2area_.NewAthenaArray(nc1);
    x3area_.NewAthenaArray(nc1);
    x2area_p1_.NewAthenaArray(nc1);
    x3area_p1_.NewAthenaArray(nc1);
    vol_.NewAthenaArray(nc1);
    fx_.NewAthenaArray(nc1);
    fy_.NewAthenaArray(nc1);
    fz_.NewAthenaArray(nc1);
    div_vel_.NewAthenaArray(nc3, nc2, nc1);

    nu.NewAthenaArray(2, nc3, nc2, nc1);
    if (pmb_->pmy_mesh->ViscosityCoeff_ == nullptr)
      CalcViscCoeff_ = ConstViscosity;
    else
      CalcViscCoeff_ = pmb_->pmy_mesh->ViscosityCoeff_;
  }

  // Check if thermal conduction is active
  if (NON_BAROTROPIC_EOS) {
    kappa_iso  = pin->GetOrAddReal("problem", "kappa_iso", 0.0); // iso thermal conduction
    kappa_aniso  = pin->GetOrAddReal("problem", "kappa_aniso", 0.0); // aniso conduction
    if (kappa_iso > 0.0 || kappa_aniso > 0.0) {
      hydro_diffusion_defined = true;
      cndflx[X1DIR].NewAthenaArray(nc3, nc2, nc1+1);
      cndflx[X2DIR].NewAthenaArray(nc3, nc2+1, nc1);
      cndflx[X3DIR].NewAthenaArray(nc3+1, nc2, nc1);

      kappa.NewAthenaArray(2, nc3, nc2, nc1);
      if (pmb_->pmy_mesh->ConductionCoeff_ == nullptr)
        CalcCondCoeff_ = ConstConduction;
      else
        CalcCondCoeff_ = pmb_->pmy_mesh->ConductionCoeff_;
    }
  }

  if (hydro_diffusion_defined && RELATIVISTIC_DYNAMICS) {
    std::stringstream msg;
    msg << "### FATAL ERROR in HydroDiffusion" << std::endl
        << "Diffusion is incompatibile with relativistic dynamics" << std::endl;
    ATHENA_ERROR(msg);
  }

  if (hydro_diffusion_defined) {
    dx1_.NewAthenaArray(nc1);
    dx2_.NewAthenaArray(nc1);
    dx3_.NewAthenaArray(nc1);
    nu_tot_.NewAthenaArray(nc1);
    kappa_tot_.NewAthenaArray(nc1);
  }
}


//----------------------------------------------------------------------------------------
//! \fn void HydroDiffusion::CalcDiffusionFlux
//! \brief Wrapper function for calculating local diffusivity coefficients and then
//!  diffusion fluxes (of various flavors) for overall hydro flux.

void HydroDiffusion::CalcDiffusionFlux(const AthenaArray<Real> &prim,
                                       const AthenaArray<Real> &iprim,
                                       const AthenaArray<Real> &bcc) {
  SetDiffusivity(prim, bcc);

  if (nu_iso > 0.0 || nu_aniso > 0.0) ClearFlux(visflx);
  if (nu_iso > 0.0) ViscousFluxIso(prim, iprim, visflx);
  if (nu_aniso > 0.0) ViscousFluxAniso(prim, iprim, visflx);

  if (kappa_iso > 0.0 || kappa_aniso > 0.0) ClearFlux(cndflx);
  if (kappa_iso > 0.0) ThermalFluxIso(prim, cndflx);
  if (kappa_aniso > 0.0) ThermalFluxAniso(prim, cndflx);

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void HydroDiffusion::AddDiffusionEnergyFlux
//! \brief Adds only diffusion energy flux to hydro flux
//
// TODO(felker): completely general operation for any 2x AthenaArray<Real>
// flux[3]. Replace calls with below "AddDiffusionFlux" on shallow-sliced inputs.

void HydroDiffusion::AddDiffusionEnergyFlux(AthenaArray<Real> *flux_src,
                                            AthenaArray<Real> *flux_des) {
  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;

  AthenaArray<Real> &x1flux = flux_des[X1DIR];
  AthenaArray<Real> &x2flux = flux_des[X2DIR];
  AthenaArray<Real> &x3flux = flux_des[X3DIR];
  AthenaArray<Real> &x1diflx = flux_src[X1DIR];
  AthenaArray<Real> &x2diflx = flux_src[X2DIR];
  AthenaArray<Real> &x3diflx = flux_src[X3DIR];

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        x1flux(IEN,k,j,i) += x1diflx(k,j,i);
        if (i == ie) x1flux(IEN,k,j,i+1) += x1diflx(k,j,i+1);
        if (pmb_->block_size.nx2 > 1) {
          x2flux(IEN,k,j,i) += x2diflx(k,j,i);
          if (j == je) x2flux(IEN,k,j+1,i) += x2diflx(k,j+1,i);
        }
        if (pmb_->block_size.nx3 > 1) {
          x3flux(IEN,k,j,i) += x3diflx(k,j,i);
          if (k == ke) x3flux(IEN,k+1,j,i) += x3diflx(k+1,j,i);
        }
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void HydroDiffusion::AddDiffusionFlux
//!  \brief Adds all componenets of diffusion flux to hydro flux
//
// TODO(felker): completely general operation for any 2x AthenaArray<Real> flux[3]. Move
// from this class, and remove "Diffusion" from name.

void HydroDiffusion::AddDiffusionFlux(AthenaArray<Real> *flux_src,
                                      AthenaArray<Real> *flux_des) {
  int size1 = flux_des[X1DIR].GetSize();
#pragma omp simd
  for (int i=0; i<size1; ++i)
    flux_des[X1DIR](i) += flux_src[X1DIR](i);

  if (pmb_->block_size.nx2 > 1) {
    int size2 = flux_des[X2DIR].GetSize();
#pragma omp simd
    for (int i=0; i<size2; ++i)
      flux_des[X2DIR](i) += flux_src[X2DIR](i);
  }
  if (pmb_->block_size.nx3 > 1) {
    int size3 = flux_des[X3DIR].GetSize();
#pragma omp simd
    for (int i=0; i<size3; ++i)
      flux_des[X3DIR](i) += flux_src[X3DIR](i);
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void HydroDiffusion::ClearFlux
//! \brief Reset diffusion flux back to zeros
//
// TODO(felker): completely general operation for any 1x AthenaArray<Real> flux[3]. Move
// from this class.

void HydroDiffusion::ClearFlux(AthenaArray<Real> *flux) {
  flux[X1DIR].ZeroClear();
  flux[X2DIR].ZeroClear();
  flux[X3DIR].ZeroClear();
  return;
}

//---------------------------------------------------------------------------------------
//! \fn void HydroDiffusion::SetDiffusivity(const Athena<Real> &w,
//                                          const AthenaArray<Real> &bc)
//! \brief  Set local hydro diffusion coefficients based on the fluid and
//! magnetic field variables across the mesh.
//! Called within HydroDiffusion::CalcDiffusionFlux wrapper function.

void HydroDiffusion::SetDiffusivity(const AthenaArray<Real> &w,
                                    const AthenaArray<Real> &bc) {
  int il = pmb_->is - NGHOST; int jl = pmb_->js; int kl = pmb_->ks;
  int iu = pmb_->ie + NGHOST; int ju = pmb_->je; int ku = pmb_->ke;
  if (pmb_->block_size.nx2 > 1) {
    jl -= NGHOST; ju += NGHOST;
  }
  if (pmb_->block_size.nx3 > 1) {
    kl -= NGHOST; ku += NGHOST;
  }

  // set viscosity using func ptr
  if (nu_iso > 0.0 || nu_aniso > 0.0)
    CalcViscCoeff_(this, pmb_, w, bc, il, iu, jl, ju, kl, ku);
  // set thermal conduction using func ptr
  if (kappa_iso > 0.0 || kappa_aniso > 0.0)
    CalcCondCoeff_(this, pmb_, w, bc, il, iu, jl, ju, kl, ku);

  return;
}

//---------------------------------------------------------------------------------------
//! \fn void HydroDiffusion::NewDiffusionDt(Real &dt_vis, Real &dt_cnd)
//! \brief Get the hydro diffusion timestep
//!
//! currently return dt for viscous and conduction processes

void HydroDiffusion::NewDiffusionDt(Real &dt_vis, Real &dt_cnd) {
  Real real_max = std::numeric_limits<Real>::max();
  const bool f2 = pmb_->pmy_mesh->f2;
  const bool f3 = pmb_->pmy_mesh->f3;
  int il = pmb_->is - NGHOST; int jl = pmb_->js; int kl = pmb_->ks;
  int iu = pmb_->ie + NGHOST; int ju = pmb_->je; int ku = pmb_->ke;
  Real fac;
  if (f3)
    fac = 1.0/6.0;
  else if (f2)
    fac = 0.25;
  else
    fac = 0.5;

  dt_vis = real_max;
  dt_cnd = real_max;

  AthenaArray<Real> &nu_t = nu_tot_;
  AthenaArray<Real> &kappa_t = kappa_tot_;
  AthenaArray<Real> &len = dx1_, &dx2 = dx2_, &dx3 = dx3_;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        nu_t(i) = 0.0;
        kappa_t(i) = 0.0;
      }
      if (nu_iso > 0.0) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) nu_t(i) += nu(DiffProcess::iso,k,j,i);
      }
      if (nu_aniso > 0.0) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) nu_t(i) += nu(DiffProcess::aniso,k,j,i);
      }
      if (kappa_iso > 0.0) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) kappa_t(i) += kappa(DiffProcess::iso,k,j,i);
      }
      if (kappa_aniso > 0.0) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) kappa_t(i) += kappa(DiffProcess::aniso,k,j,i);
      }
      pmb_->pcoord->CenterWidth1(k, j, il, iu, len);
      pmb_->pcoord->CenterWidth2(k, j, il, iu, dx2);
      pmb_->pcoord->CenterWidth3(k, j, il, iu, dx3);
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        len(i) = (f2) ? std::min(len(i), dx2(i)) : len(i);
        len(i) = (f3) ? std::min(len(i), dx3(i)) : len(i);
      }
      if ((nu_iso > 0.0) || (nu_aniso > 0.0)) {
        for (int i=il; i<=iu; ++i)
          dt_vis = std::min(dt_vis, static_cast<Real>(
              SQR(len(i))*fac/(nu_t(i) + TINY_NUMBER)));
      }
      if ((kappa_iso > 0.0) || (kappa_aniso > 0.0)) {
        for (int i=il; i<=iu; ++i)
          dt_cnd = std::min(dt_cnd, static_cast<Real>(
              SQR(len(i))*fac/(kappa_t(i) + TINY_NUMBER)));
      }
    }
  }
  return;
}
