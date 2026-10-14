//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hydro_diffusion_fourth.cpp
//! \brief fourth-order accurate isotropic viscous and conductive fluxes
//!
//! Companions to the second-order operators in viscosity.cpp/conduction.cpp for the
//! fourth-order finite-volume scheme (time/xorder=4) on uniform Cartesian grids.
//! Strategy (all steps are O(h^4) accurate):
//!  1. deconvolve the cell-averaged primitives to cell-centered point values,
//!     wc = <w> - h^2/24 Lap(<w>);
//!  2. evaluate the diffusive flux at face-centered points using fourth-order
//!     finite-difference/interpolation stencils of the point values
//!     (normal derivative, transverse derivative, face interpolation);
//!  3. convert the face-centered point flux to a face-averaged flux with the
//!     transverse Laplacian correction, <F> = F + h^2/24 Lap_perp(F).
//! The fluxes stored in visflx[]/cndflx[] are therefore face-averaged and are added to
//! the (face-averaged) hydro fluxes unchanged in Hydro::AddDiffusionFluxes().
//!
//! Variable diffusion coefficients (nu, kappa arrays) are interpolated to faces at
//! fourth order but are treated as point values as provided by the coefficient
//! functions; this is exact for constant coefficients.

// C headers

// C++ headers

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../hydro.hpp"
#include "hydro_diffusion.hpp"

namespace {
// fourth-order stencil coefficients (uniform grid spacing h)
// interpolation of point values to the face between cells m (=f-1) and f:
//   I4 = (9(a_m + a_f) - (a_{m-1} + a_{f+1}))/16
// first derivative at the same face:
//   D4 = (27(a_f - a_m) - (a_{f+1} - a_{m-1}))/(24 h)
// first derivative at a cell center:
//   C4 = (8(a_{c+1} - a_{c-1}) - (a_{c+2} - a_{c-2}))/(12 h)
constexpr Real ONE_24TH = 1.0/24.0;
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void HydroDiffusion::DeconvolvePrimitivesFourth
//! \brief compute point-valued (cell-centered) primitives wc_ (and temperature tc_)
//! from the cell-averaged primitives via <w> - h^2/24 Lap(<w>)

void HydroDiffusion::DeconvolvePrimitivesFourth(const AthenaArray<Real> &p) {
  const bool f2 = pmb_->pmy_mesh->f2;
  const bool f3 = pmb_->pmy_mesh->f3;
  const int is = pmb_->is, ie = pmb_->ie, js = pmb_->js, je = pmb_->je,
            ks = pmb_->ks, ke = pmb_->ke;
  const int il = is - NGHOST + 1, iu = ie + NGHOST - 1;
  const int jl = f2 ? js - NGHOST + 1 : js, ju = f2 ? je + NGHOST - 1 : je;
  const int kl = f3 ? ks - NGHOST + 1 : ks, ku = f3 ? ke + NGHOST - 1 : ke;

  for (int n=0; n<NHYDRO; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          Real lap = (p(n,k,j,i-1) - 2.0*p(n,k,j,i) + p(n,k,j,i+1));
          if (f2) lap += (p(n,k,j-1,i) - 2.0*p(n,k,j,i) + p(n,k,j+1,i));
          if (f3) lap += (p(n,k-1,j,i) - 2.0*p(n,k,j,i) + p(n,k+1,j,i));
          wc_(n,k,j,i) = p(n,k,j,i) - ONE_24TH*lap;
        }
      }
    }
  }

  if (kappa_iso > 0.0 && NON_BAROTROPIC_EOS) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i)
          tc_(k,j,i) = wc_(IPR,k,j,i)/wc_(IDN,k,j,i);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void HydroDiffusion::ViscousFluxIsoFourth
//! \brief fourth-order isotropic viscous fluxes on a uniform Cartesian grid

void HydroDiffusion::ViscousFluxIsoFourth(const AthenaArray<Real> &p,
                                          AthenaArray<Real> *flx) {
  const bool f2 = pmb_->pmy_mesh->f2;
  const bool f3 = pmb_->pmy_mesh->f3;
  const int is = pmb_->is, ie = pmb_->ie, js = pmb_->js, je = pmb_->je,
            ks = pmb_->ks, ke = pmb_->ke;
  const Real h1 = pco_->dx1f(is), h2 = pco_->dx2f(js), h3 = pco_->dx3f(ks);
  const Real nuiso2 = -TWO_3RD;
  AthenaArray<Real> &nu_d = nu;
  const int iso = DiffProcess::iso;

  //--- x1-fluxes ----------------------------------------------------------------------
  {
    // face-averaged output rows (transverse extension matches the 2nd-order operators)
    int jl = js, ju = je, kl = ks, ku = ke;
    if (MAGNETIC_FIELDS_ENABLED && f2) {
      jl = js-1, ju = je+1;
      if (f3) kl = ks-1, ku = ke+1;
    }
    // face-centered point-value rows (+1 in each active transverse dim for averaging)
    const int pjl = jl - (f2 ? 1 : 0), pju = ju + (f2 ? 1 : 0);
    const int pkl = kl - (f3 ? 1 : 0), pku = ku + (f3 ? 1 : 0);

    // cell-centered transverse gradients: gc1_=dv1/dx2, gc2_=dv2/dx2 (f2);
    //                                     gc3_=dv1/dx3, gc4_=dv3/dx3 (f3)
    for (int k=pkl; k<=pku; ++k) {
      for (int j=pjl; j<=pju; ++j) {
        if (f2) {
#pragma omp simd
          for (int i=is-2; i<=ie+2; ++i) {
            gc1_(k,j,i) = (8.0*(wc_(IVX,k,j+1,i) - wc_(IVX,k,j-1,i))
                           - (wc_(IVX,k,j+2,i) - wc_(IVX,k,j-2,i)))/(12.0*h2);
            gc2_(k,j,i) = (8.0*(wc_(IVY,k,j+1,i) - wc_(IVY,k,j-1,i))
                           - (wc_(IVY,k,j+2,i) - wc_(IVY,k,j-2,i)))/(12.0*h2);
          }
        }
        if (f3) {
#pragma omp simd
          for (int i=is-2; i<=ie+2; ++i) {
            gc3_(k,j,i) = (8.0*(wc_(IVX,k+1,j,i) - wc_(IVX,k-1,j,i))
                           - (wc_(IVX,k+2,j,i) - wc_(IVX,k-2,j,i)))/(12.0*h3);
            gc4_(k,j,i) = (8.0*(wc_(IVZ,k+1,j,i) - wc_(IVZ,k-1,j,i))
                           - (wc_(IVZ,k+2,j,i) - wc_(IVZ,k-2,j,i)))/(12.0*h3);
          }
        }
      }
    }

    // face-centered point fluxes
    for (int k=pkl; k<=pku; ++k) {
      for (int j=pjl; j<=pju; ++j) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          Real dv1dx = (27.0*(wc_(IVX,k,j,i) - wc_(IVX,k,j,i-1))
                        - (wc_(IVX,k,j,i+1) - wc_(IVX,k,j,i-2)))/(24.0*h1);
          Real dv2dx = (27.0*(wc_(IVY,k,j,i) - wc_(IVY,k,j,i-1))
                        - (wc_(IVY,k,j,i+1) - wc_(IVY,k,j,i-2)))/(24.0*h1);
          Real dv3dx = (27.0*(wc_(IVZ,k,j,i) - wc_(IVZ,k,j,i-1))
                        - (wc_(IVZ,k,j,i+1) - wc_(IVZ,k,j,i-2)))/(24.0*h1);
          Real dv1dy = 0.0, dv2dy = 0.0, dv1dz = 0.0, dv3dz = 0.0;
          if (f2) {
            dv1dy = (9.0*(gc1_(k,j,i-1) + gc1_(k,j,i))
                     - (gc1_(k,j,i-2) + gc1_(k,j,i+1)))/16.0;
            dv2dy = (9.0*(gc2_(k,j,i-1) + gc2_(k,j,i))
                     - (gc2_(k,j,i-2) + gc2_(k,j,i+1)))/16.0;
          }
          if (f3) {
            dv1dz = (9.0*(gc3_(k,j,i-1) + gc3_(k,j,i))
                     - (gc3_(k,j,i-2) + gc3_(k,j,i+1)))/16.0;
            dv3dz = (9.0*(gc4_(k,j,i-1) + gc4_(k,j,i))
                     - (gc4_(k,j,i-2) + gc4_(k,j,i+1)))/16.0;
          }
          Real divv = dv1dx + dv2dy + dv3dz;
          Real rhof = (9.0*(wc_(IDN,k,j,i-1) + wc_(IDN,k,j,i))
                       - (wc_(IDN,k,j,i-2) + wc_(IDN,k,j,i+1)))/16.0;
          Real nuf  = (9.0*(nu_d(iso,k,j,i-1) + nu_d(iso,k,j,i))
                       - (nu_d(iso,k,j,i-2) + nu_d(iso,k,j,i+1)))/16.0;
          Real mu = rhof*nuf;
          Real flx1 = -mu*(2.0*dv1dx + nuiso2*divv);
          Real flx2 = -mu*(dv2dx + dv1dy);
          Real flx3 = -mu*(dv3dx + dv1dz);
          fpt_(IM1,k,j,i) = flx1;
          fpt_(IM2,k,j,i) = flx2;
          fpt_(IM3,k,j,i) = flx3;
          if (NON_BAROTROPIC_EOS) {
            Real v1f = (9.0*(wc_(IVX,k,j,i-1) + wc_(IVX,k,j,i))
                        - (wc_(IVX,k,j,i-2) + wc_(IVX,k,j,i+1)))/16.0;
            Real v2f = (9.0*(wc_(IVY,k,j,i-1) + wc_(IVY,k,j,i))
                        - (wc_(IVY,k,j,i-2) + wc_(IVY,k,j,i+1)))/16.0;
            Real v3f = (9.0*(wc_(IVZ,k,j,i-1) + wc_(IVZ,k,j,i))
                        - (wc_(IVZ,k,j,i-2) + wc_(IVZ,k,j,i+1)))/16.0;
            fpt_(IEN,k,j,i) = v1f*flx1 + v2f*flx2 + v3f*flx3;
          }
        }
      }
    }

    // face-averaging via the transverse Laplacian correction
    AthenaArray<Real> &x1flux = flx[X1DIR];
    for (int n=IM1; n<=(NON_BAROTROPIC_EOS ? IEN : IM3); ++n) {
      for (int k=kl; k<=ku; ++k) {
        for (int j=jl; j<=ju; ++j) {
#pragma omp simd
          for (int i=is; i<=ie+1; ++i) {
            Real corr = 0.0;
            if (f2) corr += (fpt_(n,k,j-1,i) - 2.0*fpt_(n,k,j,i) + fpt_(n,k,j+1,i));
            if (f3) corr += (fpt_(n,k-1,j,i) - 2.0*fpt_(n,k,j,i) + fpt_(n,k+1,j,i));
            x1flux(n,k,j,i) += fpt_(n,k,j,i) + ONE_24TH*corr;
          }
        }
      }
    }
  }

  //--- x2-fluxes ----------------------------------------------------------------------
  if (f2) {
    int il = is, iu = ie, kl = ks, ku = ke;
    if (MAGNETIC_FIELDS_ENABLED) {
      il = is-1, iu = ie+1;
      if (f3) kl = ks-1, ku = ke+1;
    }
    const int pil = il - 1, piu = iu + 1;
    const int pkl = kl - (f3 ? 1 : 0), pku = ku + (f3 ? 1 : 0);

    // gc1_=dv1/dx1, gc2_=dv2/dx1; gc3_=dv2/dx3, gc4_=dv3/dx3 (f3)
    for (int k=pkl; k<=pku; ++k) {
      for (int j=js-2; j<=je+2; ++j) {
#pragma omp simd
        for (int i=pil; i<=piu; ++i) {
          gc1_(k,j,i) = (8.0*(wc_(IVX,k,j,i+1) - wc_(IVX,k,j,i-1))
                         - (wc_(IVX,k,j,i+2) - wc_(IVX,k,j,i-2)))/(12.0*h1);
          gc2_(k,j,i) = (8.0*(wc_(IVY,k,j,i+1) - wc_(IVY,k,j,i-1))
                         - (wc_(IVY,k,j,i+2) - wc_(IVY,k,j,i-2)))/(12.0*h1);
        }
        if (f3) {
#pragma omp simd
          for (int i=pil; i<=piu; ++i) {
            gc3_(k,j,i) = (8.0*(wc_(IVY,k+1,j,i) - wc_(IVY,k-1,j,i))
                           - (wc_(IVY,k+2,j,i) - wc_(IVY,k-2,j,i)))/(12.0*h3);
            gc4_(k,j,i) = (8.0*(wc_(IVZ,k+1,j,i) - wc_(IVZ,k-1,j,i))
                           - (wc_(IVZ,k+2,j,i) - wc_(IVZ,k-2,j,i)))/(12.0*h3);
          }
        }
      }
    }

    for (int k=pkl; k<=pku; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=pil; i<=piu; ++i) {
          Real dv1dy = (27.0*(wc_(IVX,k,j,i) - wc_(IVX,k,j-1,i))
                        - (wc_(IVX,k,j+1,i) - wc_(IVX,k,j-2,i)))/(24.0*h2);
          Real dv2dy = (27.0*(wc_(IVY,k,j,i) - wc_(IVY,k,j-1,i))
                        - (wc_(IVY,k,j+1,i) - wc_(IVY,k,j-2,i)))/(24.0*h2);
          Real dv3dy = (27.0*(wc_(IVZ,k,j,i) - wc_(IVZ,k,j-1,i))
                        - (wc_(IVZ,k,j+1,i) - wc_(IVZ,k,j-2,i)))/(24.0*h2);
          Real dv1dx = (9.0*(gc1_(k,j-1,i) + gc1_(k,j,i))
                        - (gc1_(k,j-2,i) + gc1_(k,j+1,i)))/16.0;
          Real dv2dx = (9.0*(gc2_(k,j-1,i) + gc2_(k,j,i))
                        - (gc2_(k,j-2,i) + gc2_(k,j+1,i)))/16.0;
          Real dv2dz = 0.0, dv3dz = 0.0;
          if (f3) {
            dv2dz = (9.0*(gc3_(k,j-1,i) + gc3_(k,j,i))
                     - (gc3_(k,j-2,i) + gc3_(k,j+1,i)))/16.0;
            dv3dz = (9.0*(gc4_(k,j-1,i) + gc4_(k,j,i))
                     - (gc4_(k,j-2,i) + gc4_(k,j+1,i)))/16.0;
          }
          Real divv = dv1dx + dv2dy + dv3dz;
          Real rhof = (9.0*(wc_(IDN,k,j-1,i) + wc_(IDN,k,j,i))
                       - (wc_(IDN,k,j-2,i) + wc_(IDN,k,j+1,i)))/16.0;
          Real nuf  = (9.0*(nu_d(iso,k,j-1,i) + nu_d(iso,k,j,i))
                       - (nu_d(iso,k,j-2,i) + nu_d(iso,k,j+1,i)))/16.0;
          Real mu = rhof*nuf;
          Real flx1 = -mu*(dv1dy + dv2dx);
          Real flx2 = -mu*(2.0*dv2dy + nuiso2*divv);
          Real flx3 = -mu*(dv3dy + dv2dz);
          fpt_(IM1,k,j,i) = flx1;
          fpt_(IM2,k,j,i) = flx2;
          fpt_(IM3,k,j,i) = flx3;
          if (NON_BAROTROPIC_EOS) {
            Real v1f = (9.0*(wc_(IVX,k,j-1,i) + wc_(IVX,k,j,i))
                        - (wc_(IVX,k,j-2,i) + wc_(IVX,k,j+1,i)))/16.0;
            Real v2f = (9.0*(wc_(IVY,k,j-1,i) + wc_(IVY,k,j,i))
                        - (wc_(IVY,k,j-2,i) + wc_(IVY,k,j+1,i)))/16.0;
            Real v3f = (9.0*(wc_(IVZ,k,j-1,i) + wc_(IVZ,k,j,i))
                        - (wc_(IVZ,k,j-2,i) + wc_(IVZ,k,j+1,i)))/16.0;
            fpt_(IEN,k,j,i) = v1f*flx1 + v2f*flx2 + v3f*flx3;
          }
        }
      }
    }

    AthenaArray<Real> &x2flux = flx[X2DIR];
    for (int n=IM1; n<=(NON_BAROTROPIC_EOS ? IEN : IM3); ++n) {
      for (int k=kl; k<=ku; ++k) {
        for (int j=js; j<=je+1; ++j) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            Real corr = (fpt_(n,k,j,i-1) - 2.0*fpt_(n,k,j,i) + fpt_(n,k,j,i+1));
            if (f3) corr += (fpt_(n,k-1,j,i) - 2.0*fpt_(n,k,j,i) + fpt_(n,k+1,j,i));
            x2flux(n,k,j,i) += fpt_(n,k,j,i) + ONE_24TH*corr;
          }
        }
      }
    }
  }

  //--- x3-fluxes ----------------------------------------------------------------------
  if (f3) {
    int il = is, iu = ie, jl = js, ju = je;
    if (MAGNETIC_FIELDS_ENABLED) {
      il = is-1, iu = ie+1, jl = js-1, ju = je+1;
    }
    const int pil = il - 1, piu = iu + 1;
    const int pjl = jl - 1, pju = ju + 1;

    // gc1_=dv1/dx1, gc2_=dv3/dx1; gc3_=dv2/dx2, gc4_=dv3/dx2
    for (int k=ks-2; k<=ke+2; ++k) {
      for (int j=pjl; j<=pju; ++j) {
#pragma omp simd
        for (int i=pil; i<=piu; ++i) {
          gc1_(k,j,i) = (8.0*(wc_(IVX,k,j,i+1) - wc_(IVX,k,j,i-1))
                         - (wc_(IVX,k,j,i+2) - wc_(IVX,k,j,i-2)))/(12.0*h1);
          gc2_(k,j,i) = (8.0*(wc_(IVZ,k,j,i+1) - wc_(IVZ,k,j,i-1))
                         - (wc_(IVZ,k,j,i+2) - wc_(IVZ,k,j,i-2)))/(12.0*h1);
          gc3_(k,j,i) = (8.0*(wc_(IVY,k,j+1,i) - wc_(IVY,k,j-1,i))
                         - (wc_(IVY,k,j+2,i) - wc_(IVY,k,j-2,i)))/(12.0*h2);
          gc4_(k,j,i) = (8.0*(wc_(IVZ,k,j+1,i) - wc_(IVZ,k,j-1,i))
                         - (wc_(IVZ,k,j+2,i) - wc_(IVZ,k,j-2,i)))/(12.0*h2);
        }
      }
    }

    for (int k=ks; k<=ke+1; ++k) {
      for (int j=pjl; j<=pju; ++j) {
#pragma omp simd
        for (int i=pil; i<=piu; ++i) {
          Real dv1dz = (27.0*(wc_(IVX,k,j,i) - wc_(IVX,k-1,j,i))
                        - (wc_(IVX,k+1,j,i) - wc_(IVX,k-2,j,i)))/(24.0*h3);
          Real dv2dz = (27.0*(wc_(IVY,k,j,i) - wc_(IVY,k-1,j,i))
                        - (wc_(IVY,k+1,j,i) - wc_(IVY,k-2,j,i)))/(24.0*h3);
          Real dv3dz = (27.0*(wc_(IVZ,k,j,i) - wc_(IVZ,k-1,j,i))
                        - (wc_(IVZ,k+1,j,i) - wc_(IVZ,k-2,j,i)))/(24.0*h3);
          Real dv1dx = (9.0*(gc1_(k-1,j,i) + gc1_(k,j,i))
                        - (gc1_(k-2,j,i) + gc1_(k+1,j,i)))/16.0;
          Real dv3dx = (9.0*(gc2_(k-1,j,i) + gc2_(k,j,i))
                        - (gc2_(k-2,j,i) + gc2_(k+1,j,i)))/16.0;
          Real dv2dy = (9.0*(gc3_(k-1,j,i) + gc3_(k,j,i))
                        - (gc3_(k-2,j,i) + gc3_(k+1,j,i)))/16.0;
          Real dv3dy = (9.0*(gc4_(k-1,j,i) + gc4_(k,j,i))
                        - (gc4_(k-2,j,i) + gc4_(k+1,j,i)))/16.0;
          Real divv = dv1dx + dv2dy + dv3dz;
          Real rhof = (9.0*(wc_(IDN,k-1,j,i) + wc_(IDN,k,j,i))
                       - (wc_(IDN,k-2,j,i) + wc_(IDN,k+1,j,i)))/16.0;
          Real nuf  = (9.0*(nu_d(iso,k-1,j,i) + nu_d(iso,k,j,i))
                       - (nu_d(iso,k-2,j,i) + nu_d(iso,k+1,j,i)))/16.0;
          Real mu = rhof*nuf;
          Real flx1 = -mu*(dv1dz + dv3dx);
          Real flx2 = -mu*(dv2dz + dv3dy);
          Real flx3 = -mu*(2.0*dv3dz + nuiso2*divv);
          fpt_(IM1,k,j,i) = flx1;
          fpt_(IM2,k,j,i) = flx2;
          fpt_(IM3,k,j,i) = flx3;
          if (NON_BAROTROPIC_EOS) {
            Real v1f = (9.0*(wc_(IVX,k-1,j,i) + wc_(IVX,k,j,i))
                        - (wc_(IVX,k-2,j,i) + wc_(IVX,k+1,j,i)))/16.0;
            Real v2f = (9.0*(wc_(IVY,k-1,j,i) + wc_(IVY,k,j,i))
                        - (wc_(IVY,k-2,j,i) + wc_(IVY,k+1,j,i)))/16.0;
            Real v3f = (9.0*(wc_(IVZ,k-1,j,i) + wc_(IVZ,k,j,i))
                        - (wc_(IVZ,k-2,j,i) + wc_(IVZ,k+1,j,i)))/16.0;
            fpt_(IEN,k,j,i) = v1f*flx1 + v2f*flx2 + v3f*flx3;
          }
        }
      }
    }

    AthenaArray<Real> &x3flux = flx[X3DIR];
    for (int n=IM1; n<=(NON_BAROTROPIC_EOS ? IEN : IM3); ++n) {
      for (int k=ks; k<=ke+1; ++k) {
        for (int j=jl; j<=ju; ++j) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            Real corr = (fpt_(n,k,j,i-1) - 2.0*fpt_(n,k,j,i) + fpt_(n,k,j,i+1))
                        + (fpt_(n,k,j-1,i) - 2.0*fpt_(n,k,j,i) + fpt_(n,k,j+1,i));
            x3flux(n,k,j,i) += fpt_(n,k,j,i) + ONE_24TH*corr;
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void HydroDiffusion::ThermalFluxIsoFourth
//! \brief fourth-order isotropic conductive flux, F = -kappa dT/dx, on a uniform
//! Cartesian grid (temperature point values tc_ from DeconvolvePrimitivesFourth)

void HydroDiffusion::ThermalFluxIsoFourth(const AthenaArray<Real> &p,
                                          AthenaArray<Real> *flx) {
  const bool f2 = pmb_->pmy_mesh->f2;
  const bool f3 = pmb_->pmy_mesh->f3;
  const int is = pmb_->is, ie = pmb_->ie, js = pmb_->js, je = pmb_->je,
            ks = pmb_->ks, ke = pmb_->ke;
  const Real h1 = pco_->dx1f(is), h2 = pco_->dx2f(js), h3 = pco_->dx3f(ks);
  AthenaArray<Real> &kappa_d = kappa;
  const int iso = DiffProcess::iso;

  //--- x1-fluxes
  {
    int jl = js, ju = je, kl = ks, ku = ke;
    if (MAGNETIC_FIELDS_ENABLED && f2) {
      jl = js-1, ju = je+1;
      if (f3) kl = ks-1, ku = ke+1;
    }
    const int pjl = jl - (f2 ? 1 : 0), pju = ju + (f2 ? 1 : 0);
    const int pkl = kl - (f3 ? 1 : 0), pku = ku + (f3 ? 1 : 0);
    for (int k=pkl; k<=pku; ++k) {
      for (int j=pjl; j<=pju; ++j) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          Real dTdx = (27.0*(tc_(k,j,i) - tc_(k,j,i-1))
                       - (tc_(k,j,i+1) - tc_(k,j,i-2)))/(24.0*h1);
          Real kappaf = (9.0*(kappa_d(iso,k,j,i-1) + kappa_d(iso,k,j,i))
                         - (kappa_d(iso,k,j,i-2) + kappa_d(iso,k,j,i+1)))/16.0;
          fpt_(0,k,j,i) = -kappaf*dTdx;
        }
      }
    }
    AthenaArray<Real> &x1flux = flx[X1DIR];
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          Real corr = 0.0;
          if (f2) corr += (fpt_(0,k,j-1,i) - 2.0*fpt_(0,k,j,i) + fpt_(0,k,j+1,i));
          if (f3) corr += (fpt_(0,k-1,j,i) - 2.0*fpt_(0,k,j,i) + fpt_(0,k+1,j,i));
          x1flux(k,j,i) += fpt_(0,k,j,i) + ONE_24TH*corr;
        }
      }
    }
  }

  //--- x2-fluxes
  if (f2) {
    int il = is, iu = ie, kl = ks, ku = ke;
    if (MAGNETIC_FIELDS_ENABLED) {
      il = is-1, iu = ie+1;
      if (f3) kl = ks-1, ku = ke+1;
    }
    const int pil = il - 1, piu = iu + 1;
    const int pkl = kl - (f3 ? 1 : 0), pku = ku + (f3 ? 1 : 0);
    for (int k=pkl; k<=pku; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=pil; i<=piu; ++i) {
          Real dTdy = (27.0*(tc_(k,j,i) - tc_(k,j-1,i))
                       - (tc_(k,j+1,i) - tc_(k,j-2,i)))/(24.0*h2);
          Real kappaf = (9.0*(kappa_d(iso,k,j-1,i) + kappa_d(iso,k,j,i))
                         - (kappa_d(iso,k,j-2,i) + kappa_d(iso,k,j+1,i)))/16.0;
          fpt_(0,k,j,i) = -kappaf*dTdy;
        }
      }
    }
    AthenaArray<Real> &x2flux = flx[X2DIR];
    for (int k=kl; k<=ku; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          Real corr = (fpt_(0,k,j,i-1) - 2.0*fpt_(0,k,j,i) + fpt_(0,k,j,i+1));
          if (f3) corr += (fpt_(0,k-1,j,i) - 2.0*fpt_(0,k,j,i) + fpt_(0,k+1,j,i));
          x2flux(k,j,i) += fpt_(0,k,j,i) + ONE_24TH*corr;
        }
      }
    }
  }

  //--- x3-fluxes
  if (f3) {
    int il = is, iu = ie, jl = js, ju = je;
    if (MAGNETIC_FIELDS_ENABLED) {
      il = is-1, iu = ie+1, jl = js-1, ju = je+1;
    }
    const int pil = il - 1, piu = iu + 1;
    const int pjl = jl - 1, pju = ju + 1;
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=pjl; j<=pju; ++j) {
#pragma omp simd
        for (int i=pil; i<=piu; ++i) {
          Real dTdz = (27.0*(tc_(k,j,i) - tc_(k-1,j,i))
                       - (tc_(k+1,j,i) - tc_(k-2,j,i)))/(24.0*h3);
          Real kappaf = (9.0*(kappa_d(iso,k-1,j,i) + kappa_d(iso,k,j,i))
                         - (kappa_d(iso,k-2,j,i) + kappa_d(iso,k+1,j,i)))/16.0;
          fpt_(0,k,j,i) = -kappaf*dTdz;
        }
      }
    }
    AthenaArray<Real> &x3flux = flx[X3DIR];
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          Real corr = (fpt_(0,k,j,i-1) - 2.0*fpt_(0,k,j,i) + fpt_(0,k,j,i+1))
                      + (fpt_(0,k,j-1,i) - 2.0*fpt_(0,k,j,i) + fpt_(0,k,j+1,i));
          x3flux(k,j,i) += fpt_(0,k,j,i) + ONE_24TH*corr;
        }
      }
    }
  }
  return;
}
