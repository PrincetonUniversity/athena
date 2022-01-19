//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file viscosity.cpp
//! \brief functions to calculate viscous stresses

// C headers

// C++ headers

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../eos/eos.hpp"
#include "../../mesh/mesh.hpp"
#include "../hydro.hpp"
#include "hydro_diffusion.hpp"

//----------------------------------------------------------------------------------------
//! \fn void HydroDiffusion::ViscousFluxIso
//! \brief Calculate isotropic viscous stress as fluxes

void HydroDiffusion::ViscousFluxIso(const AthenaArray<Real> &p,
                     const AthenaArray<Real> &p_i, AthenaArray<Real> *flx) {
  Hydro *ph = pmb_->phydro;
  const bool f2 = pmb_->pmy_mesh->f2;
  const bool f3 = pmb_->pmy_mesh->f3;
  AthenaArray<Real> &x1flux = flx[X1DIR];
  AthenaArray<Real> &x2flux = flx[X2DIR];
  AthenaArray<Real> &x3flux = flx[X3DIR];
  int il, iu, jl, ju, kl, ku;
  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;
  Real nu1, denf, flx1, flx2, flx3;
  Real nuiso2 = - TWO_3RD;

  DivVelocity(p_i, div_vel_);

  // Calculate the flux across each face.
  // i-direction
  jl = js, ju = je, kl = ks, ku = ke;
  if (MAGNETIC_FIELDS_ENABLED) {
    if (f2) {
      if (!f3) // 2D MHD limits
        jl = js-1, ju = je+1, kl = ks, ku = ke;
      else // 3D MHD limits
        jl = js-1, ju = je+1, kl = ks-1, ku = ke+1;
    }
  }
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      FaceXdx(k, j, is, ie+1, p_i, fx_);
      FaceXdy(k, j, is, ie+1, p_i, fy_);
      FaceXdz(k, j, is, ie+1, p_i, fz_);
#pragma omp simd private(nu1, denf, flx1, flx2, flx3)
      for (int i=is; i<=ie+1; ++i) {
        nu1  = 0.5*(nu(DiffProcess::iso,k,j,i)   + nu(DiffProcess::iso,k,j,i-1));
        denf = 0.5*(p_i(IDN,k,j,i) + p_i(IDN,k,j,i-1));
        flx1 = -denf*nu1*(fx_(i) + nuiso2*0.5*(div_vel_(k,j,i) + div_vel_(k,j,i-1)));
        flx2 = -denf*nu1*fy_(i);
        flx3 = -denf*nu1*fz_(i);
        x1flux(IM1,k,j,i) += flx1;
        x1flux(IM2,k,j,i) += flx2;
        x1flux(IM3,k,j,i) += flx3;
        if (NON_BAROTROPIC_EOS)
          x1flux(IEN,k,j,i) += 0.5*((p(IM1,k,j,i-1) + p(IM1,k,j,i))*flx1 +
                                    (p(IM2,k,j,i-1) + p(IM2,k,j,i))*flx2 +
                                    (p(IM3,k,j,i-1) + p(IM3,k,j,i))*flx3);
      }
    }
  }

  // j-direction
  il = is, iu = ie, kl = ks, ku = ke;
  if (MAGNETIC_FIELDS_ENABLED) {
    if (!f3) // 2D MHD limits
      il = is-1, iu = ie+1, kl = ks, ku = ke;
    else // 3D MHD limits
      il = is-1, iu = ie+1, kl = ks-1, ku = ke+1;
  }
  if (f2) { // modify x2flux for 2D or 3D
    for (int k=kl; k<=ku; ++k) {
      for (int j=js; j<=je+1; ++j) {
        // compute fluxes
        FaceYdx(k, j, is, ie, p_i, fx_);
        FaceYdy(k, j, is, ie, p_i, fy_);
        FaceYdz(k, j, is, ie, p_i, fz_);
        // store fluxes
#pragma omp simd private(nu1, denf, flx1, flx2, flx3)
        for (int i=il; i<=iu; i++) {
          nu1  = 0.5*(nu(DiffProcess::iso,k,j,i)    + nu(DiffProcess::iso,k,j-1,i));
          denf = 0.5*(p_i(IDN,k,j-1,i)+ p_i(IDN,k,j,i));
          flx1 = -denf*nu1*fx_(i);
          flx2 = -denf*nu1*(fy_(i) + nuiso2*0.5*(div_vel_(k,j-1,i) + div_vel_(k,j,i)));
          flx3 = -denf*nu1*fz_(i);
          x2flux(IM1,k,j,i) += flx1;
          x2flux(IM2,k,j,i) += flx2;
          x2flux(IM3,k,j,i) += flx3;
          if (NON_BAROTROPIC_EOS)
            x2flux(IEN,k,j,i) += 0.5*((p(IM1,k,j,i) + p(IM1,k,j-1,i))*flx1 +
                                      (p(IM2,k,j,i) + p(IM2,k,j-1,i))*flx2 +
                                      (p(IM3,k,j,i) + p(IM3,k,j-1,i))*flx3);
        }
      }
    }
  } else { // modify x2flux for 1D
    // compute fluxes
    FaceYdx(ks, js, is, ie, p_i, fx_);
    FaceYdy(ks, js, is, ie, p_i, fy_);
    FaceYdz(ks, js, is, ie, p_i, fz_);
    // store fluxes
#pragma omp simd private(nu1, denf, flx1, flx2, flx3)
    for (int i=il; i<=iu; i++) {
      nu1  = nu(DiffProcess::iso,ks,js,i);
      denf = p_i(IDN,ks,js,i);
      flx1 = -denf*nu1*fx_(i);
      flx2 = -denf*nu1*(fy_(i) + nuiso2*div_vel_(ks,js,i));
      flx3 = -denf*nu1*fz_(i);
      x2flux(IM1,ks,js,i) += flx1;
      x2flux(IM2,ks,js,i) += flx2;
      x2flux(IM3,ks,js,i) += flx3;
      if (NON_BAROTROPIC_EOS)
        x2flux(IEN,ks,js,i) += p(IM1,ks,js,i)*flx1 +
                               p(IM2,ks,js,i)*flx2 +
                               p(IM3,ks,js,i)*flx3;
    }
#pragma omp simd
    for (int i=il; i<=iu; i++) {
      x2flux(IM1,ks,je+1,i) = x2flux(IM1,ks,js,i);
      x2flux(IM2,ks,je+1,i) = x2flux(IM2,ks,js,i);
      x2flux(IM3,ks,je+1,i) = x2flux(IM3,ks,js,i);
      if (NON_BAROTROPIC_EOS)
        x2flux(IEN,ks,je+1,i) = x2flux(IEN,ks,js,i);
    }
  }
  // k-direction
  // set the loop limits
  il = is, iu = ie, jl = js, ju = je;
  if (MAGNETIC_FIELDS_ENABLED) {
    if (f2) // 2D or 3D MHD limits
      il = is-1, iu = ie+1, jl = js-1, ju = je+1;
    else // 1D MHD limits
      il = is-1, iu = ie+1;
  }
  if (pmb_->block_size.nx3 > 1) { // modify x3flux for 3D
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=jl; j<=ju; ++j) {
        // compute fluxes
        FaceZdx(k, j, is, ie, p_i, fx_);
        FaceZdy(k, j, is, ie, p_i, fy_);
        FaceZdz(k, j, is, ie, p_i, fz_);
        // store fluxes
#pragma omp simd private(nu1, denf, flx1, flx2, flx3)
        for (int i=il; i<=iu; i++) {
          nu1  = 0.5*(nu(DiffProcess::iso,k,j,i)     + nu(DiffProcess::iso,k-1,j,i));
          denf = 0.5*(p_i(IDN,k-1,j,i) + p_i(IDN,k,j,i));
          flx1 = -denf*nu1*fx_(i);
          flx2 = -denf*nu1*fy_(i);
          flx3 = -denf*nu1*(fz_(i) + nuiso2*0.5*(div_vel_(k-1,j,i) + div_vel_(k,j,i)));
          x3flux(IM1,k,j,i) += flx1;
          x3flux(IM2,k,j,i) += flx2;
          x3flux(IM3,k,j,i) += flx3;
          if (NON_BAROTROPIC_EOS)
            x3flux(IEN,k,j,i) += 0.5*((p(IM1,k,j,i) + p(IM1,k-1,j,i))*flx1 +
                                      (p(IM2,k,j,i) + p(IM2,k-1,j,i))*flx2 +
                                      (p(IM3,k,j,i) + p(IM3,k-1,j,i))*flx3);
        }
      }
    }
  } else { // modify x2flux for 1D or 2D
    for (int j=jl; j<=ju; ++j) {
      // compute fluxes
      FaceZdx(ks, j, is, ie, p_i, fx_);
      FaceZdy(ks, j, is, ie, p_i, fy_);
      FaceZdz(ks, j, is, ie, p_i, fz_);
      // store fluxes
#pragma omp simd private(nu1, denf, flx1, flx2, flx3)
      for (int i=il; i<=iu; i++) {
        nu1 = nu(DiffProcess::iso,ks,j,i);
        denf = p_i(IDN,ks,j,i);
        flx1 = -denf*nu1*fx_(i);
        flx2 = -denf*nu1*fy_(i);
        flx3 = -denf*nu1*(fz_(i) + nuiso2*div_vel_(ks,j,i));
        x3flux(IM1,ks,j,i) += flx1;
        x3flux(IM2,ks,j,i) += flx2;
        x3flux(IM3,ks,j,i) += flx3;
        x3flux(IM1,ke+1,j,i) = x3flux(IM1,ks,j,i);
        x3flux(IM2,ke+1,j,i) = x3flux(IM2,ks,j,i);
        x3flux(IM3,ke+1,j,i) = x3flux(IM3,ks,j,i);
        if (NON_BAROTROPIC_EOS) {
          x3flux(IEN,ks,j,i) += p(IM1,ks,j,i)*flx1 +
                                p(IM2,ks,j,i)*flx2 +
                                p(IM3,ks,j,i)*flx3;
          x3flux(IEN,ke+1,j,i) = x3flux(IEN,ks,j,i);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void HydroDiffusion::ViscousFluxAniso
//! \brief Calculate anisotropic viscous stress as fluxes

void HydroDiffusion::ViscousFluxAniso(const AthenaArray<Real> &p,
                            const AthenaArray<Real> &p_i, AthenaArray<Real> *flx) {
  return;
}

//-------------------------------------------------------------------------------------
//! \fn void HydroDiffusion::DivVelocity(const AthenaArray<Real> &prim,
//!                                 AthenaArray<Real> &div_vel)
//! \brief Calculate divergence of momenta

void HydroDiffusion::DivVelocity(const AthenaArray<Real> &prim,
                                 AthenaArray<Real> &div_vel) {
  const bool f2 = pmb_->pmy_mesh->f2;
  const bool f3 = pmb_->pmy_mesh->f3;

  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;
  int il = is - 1; int iu = ie + 1;
  int jl, ju, kl, ku;
  Real area_p1, area;
  Real vel_p1, vel;

  if (!f2) // 1D
    jl = js, ju = je, kl = ks, ku = ke;
  else if (!f3) // 2D
    jl = js-1, ju = je+1, kl = ks, ku = ke;
  else // 3D
    jl = js-1, ju = je+1, kl = ks-1, ku = ke+1;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      // calculate x1-flux divergence
      pmb_->pcoord->Face1Area(k, j, il, iu+1, x1area_);
#pragma omp simd private(area_p1, area, vel_p1, vel)
      for (int i=il; i<=iu; ++i) {
        area_p1 = x1area_(i+1);
        area    = x1area_(i);
        vel_p1  = 0.5*(prim(IM1,k,j,i+1) + prim(IM1,k,j,i  ));
        vel     = 0.5*(prim(IM1,k,j,i  ) + prim(IM1,k,j,i-1));
        div_vel(k,j,i) = area_p1*vel_p1 - area*vel;
      }
      // calculate x2-flux divergnece
      if (f2) {
        pmb_->pcoord->Face2Area(k, j  , il, iu, x2area_);
        pmb_->pcoord->Face2Area(k, j+1, il, iu, x2area_p1_);
#pragma omp simd private(area_p1, area, vel_p1, vel)
        for (int i=il; i<=iu; ++i) {
          area_p1 = x2area_p1_(i);
          area    = x2area_(i);
          vel_p1  = 0.5*(prim(IM2,k,j+1,i) + prim(IM2,k,j  ,i));
          vel     = 0.5*(prim(IM2,k,j  ,i) + prim(IM2,k,j-1,i));
          div_vel(k,j,i) += area_p1*vel_p1 - area*vel;
        }
      }
      if (f3) {
        pmb_->pcoord->Face3Area(k  , j, il, iu, x3area_);
        pmb_->pcoord->Face3Area(k+1, j, il, iu, x3area_p1_);
#pragma omp simd private(area_p1, area, vel_p1, vel)
        for (int i=il; i<=iu; ++i) {
          area_p1 = x3area_p1_(i);
          area    = x3area_(i);
          vel_p1  = 0.5*(prim(IM3,k+1,j,i) + prim(IM3, k  ,j,i));
          vel     = 0.5*(prim(IM3,k  ,j,i) + prim(IM3, k-1,j,i));
          div_vel(k,j,i) += area_p1*vel_p1 - area*vel;
        }
      }
      pmb_->pcoord->CellVolume(k,j,il,iu,vol_);
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        div_vel(k,j,i) = div_vel(k,j,i)/vol_(i);
      }
    }
  }
  return;
}

// v_{x1;x1}  covariant derivative at x1 interface
void HydroDiffusion::FaceXdx(const int k, const int j, const int il, const int iu,
                             const AthenaArray<Real> &prim, AthenaArray<Real> &len) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    len(i) = 2.0*(prim(IM1,k,j,i) - prim(IM1,k,j,i-1)) / pco_->dx1v(i-1);
  }
  return;
}

// v_{x2;x1}+v_{x1;x2}  covariant derivative at x1 interface
void HydroDiffusion::FaceXdy(const int k, const int j, const int il, const int iu,
                             const AthenaArray<Real> &prim, AthenaArray<Real> &len) {
  if (pmb_->pmy_mesh->f2) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      len(i) = pco_->h2f(i)
               * (prim(IM2,k,j,i)/pco_->h2v(i) - prim(IM2,k,j,i-1)/pco_->h2v(i-1))
               / pco_->dx1v(i-1)
               // KGF: add the off-centered quantities first to preserve FP symmetry
               + 0.5*(   (prim(IM1,k,j+1,i) + prim(IM1,k,j+1,i-1))
                         - (prim(IM1,k,j-1,i) + prim(IM1,k,j-1,i-1)) )
               / pco_->h2f(i)
               / (pco_->dx2v(j-1) + pco_->dx2v(j));
    }
  } else {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      len(i) = pco_->h2f(i)
               * ( prim(IM2,k,j,i)/pco_->h2v(i) - prim(IM2,k,j,i-1)/pco_->h2v(i-1) )
               / pco_->dx1v(i-1);
    }
  }
  return;
}

// v_{x3;x1}+v_{x1;x3}  covariant derivative at x1 interface
void HydroDiffusion::FaceXdz(const int k, const int j, const int il, const int iu,
                             const AthenaArray<Real> &prim, AthenaArray<Real> &len) {
  if (pmb_->pmy_mesh->f3) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      len(i) = pco_->h31f(i)
               * (prim(IM3,k,j,i)/pco_->h31v(i) - prim(IM3,k,j,i-1)/pco_->h31v(i-1))
               / pco_->dx1v(i-1)
               // KGF: add the off-centered quantities first to preserve FP symmetry
               + 0.5*(   (prim(IM1,k+1,j,i) + prim(IM1,k+1,j,i-1))
                         - (prim(IM1,k-1,j,i) + prim(IM1,k-1,j,i-1)) )
               / pco_->h31f(i)/pco_->h32v(j) // note, more terms than FaceXdy() line
               / (pco_->dx3v(k-1) + pco_->dx3v(k));
    }
  } else {
#pragma omp simd
    for (int i=il; i<=iu; ++i)
      len(i) = pco_->h31f(i)
               * ( prim(IM3,k,j,i)/pco_->h31v(i) - prim(IM3,k,j,i-1)/pco_->h31v(i-1) )
               / pco_->dx1v(i-1);
  }
  return;
}

// v_{x1;x2}+v_{x2;x1}  covariant derivative at x2 interface
void HydroDiffusion::FaceYdx(const int k, const int j, const int il, const int iu,
                             const AthenaArray<Real> &prim, AthenaArray<Real> &len) {
  if (pmb_->pmy_mesh->f2) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      len(i) = (prim(IM1,k,j,i) - prim(IM1,k,j-1,i)) / pco_->h2v(i) / pco_->dx2v(j-1)
               + pco_->h2v(i)*0.5*
               (  (prim(IM2,k,j,i+1) + prim(IM2,k,j-1,i+1)) /pco_->h2v(i+1)
                  - (prim(IM2,k,j,i-1) + prim(IM2,k,j-1,i-1)) /pco_->h2v(i-1)
                  ) / (pco_->dx1v(i-1) + pco_->dx1v(i));
    }
  } else {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      len(i) = pco_->h2v(i)
               * ( prim(IM2,k,j,i+1)/pco_->h2v(i+1) - prim(IM2,k,j,i-1)/pco_->h2v(i-1) )
               / (pco_->dx1v(i-1) + pco_->dx1v(i));
    }
  }
  return;
}

// v_{x2;x2}  covariant derivative at x2 interface
void HydroDiffusion::FaceYdy(const int k, const int j, const int il, const int iu,
                             const AthenaArray<Real> &prim, AthenaArray<Real> &len) {
  if (pmb_->pmy_mesh->f2) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      len(i) = 2.0*(prim(IM2,k,j,i) - prim(IM2,k,j-1,i)) / pco_->h2v(i) / pco_->dx2v(j-1)
               + (prim(IM1,k,j,i) + prim(IM1,k,j-1,i)) / pco_->h2v(i) * pco_->dh2vd1(i);
    }
  } else {
#pragma omp simd
    for (int i=il; i<=iu; ++i)
      len(i) = 2.0*prim(IM1,k,j,i) / pco_->h2v(i) * pco_->dh2vd1(i);
  }
  return;
}

// v_{x3;x2}+v_{x2;x3}  covariant derivative at x2 interface
void HydroDiffusion::FaceYdz(const int k, const int j, const int il, const int iu,
                             const AthenaArray<Real> &prim, AthenaArray<Real> &len) {
  if (pmb_->pmy_mesh->f3) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      len(i) = pco_->h32f(j)
               * ( prim(IM3,k,j,i)/pco_->h32v(j) - prim(IM3,k,j-1,i)/pco_->h32v(j-1) )
               / pco_->h2v(i) / pco_->dx2v(j-1)
               // KGF: add the off-centered quantities first to preserve FP symmetry
               + 0.5*(    (prim(IM2,k+1,j,i) + prim(IM2,k+1,j-1,i))
                          - (prim(IM2,k-1,j,i) + prim(IM2,k-1,j-1,i)) )
               / pco_->h31v(i)
               / pco_->h32f(j) / (pco_->dx3v(k-1) + pco_->dx3v(k));
    }
  } else if (pmb_->pmy_mesh->f2) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      len(i) = pco_->h32f(j)
               * ( prim(IM3,k,j,i)/pco_->h32v(j) - prim(IM3,k,j-1,i)/pco_->h32v(j-1) )
               / pco_->h2v(i) / pco_->dx2v(j-1);
    }
  } else {
#pragma omp simd
    for (int i=il; i<=iu; ++i)
      len(i) = 0.0;
  }
  return;
}

// v_{x1;x3}+v_{x3;x1}  covariant derivative at x3 interface
void HydroDiffusion::FaceZdx(const int k, const int j, const int il, const int iu,
                             const AthenaArray<Real> &prim, AthenaArray<Real> &len) {
  if (pmb_->pmy_mesh->f3) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      len(i) = (prim(IM1,k,j,i) - prim(IM1,k-1,j,i))/pco_->dx3v(k-1)
               + 0.5*pco_->h31v(i)*(
                   (prim(IM3,k,j,i+1) + prim(IM3,k-1,j,i+1))/pco_->h31v(i+1)
                   -(prim(IM3,k,j,i-1) + prim(IM3,k-1,j,i-1))/pco_->h31v(i-1)
                                    ) / (pco_->dx1v(i-1) + pco_->dx1v(i));
    }
  } else {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      len(i) = pco_->h31v(i)
               * ( prim(IM3,k,j,i+1)/pco_->h31v(i+1) - prim(IM3,k,j,i-1)/pco_->h31v(i-1) )
               / (pco_->dx1v(i-1) + pco_->dx1v(i));
    }
  }
  return;
}

// v_{x2;x3}+v_{x3;x2}  covariant derivative at x3 interface
void HydroDiffusion::FaceZdy(const int k, const int j, const int il, const int iu,
                             const AthenaArray<Real> &prim, AthenaArray<Real> &len) {
  if (pmb_->pmy_mesh->f3) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      len(i) = (prim(IM2,k,j,i) - prim(IM2,k-1,j,i))
               / pco_->h31v(i) / pco_->h32v(j) / pco_->dx3v(k-1)
               + 0.5*pco_->h32v(j)
               * ( (prim(IM3,k,j+1,i) + prim(IM3,k-1,j+1,i))/pco_->h32v(j+1)
                   -(prim(IM3,k,j-1,i) + prim(IM3,k-1,j-1,i))/pco_->h32v(j-1) )
               / pco_->h2v(i) / (pco_->dx2v(j-1) + pco_->dx2v(j));
    }
  } else if (pmb_->pmy_mesh->f2) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      len(i) = pco_->h32v(j)
               * ( prim(IM3,k,j+1,i)/pco_->h32v(j+1) - prim(IM3,k,j-1,i)/pco_->h32v(j-1) )
               / pco_->h2v(i) / (pco_->dx2v(j-1) + pco_->dx2v(j));
    }
  } else {
#pragma omp simd
    for (int i=il; i<=iu; ++i)
      len(i) = 0.0;
  }
  return;
}

// v_{x3;x3}  covariant derivative at x3 interface
void HydroDiffusion::FaceZdz(const int k, const int j, const int il, const int iu,
                             const AthenaArray<Real> &prim, AthenaArray<Real> &len) {
  if (pmb_->pmy_mesh->f3) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      len(i) = 2.0*(prim(IM3,k,j,i) - prim(IM3,k-1,j,i))
               / pco_->dx3v(k-1) / pco_->h31v(i) / pco_->h32v(j)
               + ((prim(IM1,k,j,i) + prim(IM1,k-1,j,i))
                  * pco_->dh31vd1(i)/pco_->h31v(i))
               + ((prim(IM2,k,j,i) + prim(IM2,k-1,j,i))
                  * pco_->dh32vd2(j)/pco_->h32v(j)/pco_->h2v(i));
    }
  } else {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      len(i) = 2.0*prim(IM1,k,j,i)*pco_->dh31vd1(i)/pco_->h31v(i)
               + 2.0*prim(IM2,k,j,i)*pco_->dh32vd2(j)/pco_->h32v(j)/pco_->h2v(i);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! constant viscosity

void ConstViscosity(HydroDiffusion *phdif, MeshBlock *pmb, const AthenaArray<Real> &prim,
                    const AthenaArray<Real> &bcc, int is, int ie, int js, int je,
                    int ks, int ke) {
  if (phdif->nu_iso > 0.0) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i)
          phdif->nu(HydroDiffusion::DiffProcess::iso,k,j,i) = phdif->nu_iso;
      }
    }
  }
  if (phdif->nu_aniso > 0.0) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i)
          phdif->nu(HydroDiffusion::DiffProcess::aniso,k,j,i) = phdif->nu_aniso;
      }
    }
  }
  return;
}
