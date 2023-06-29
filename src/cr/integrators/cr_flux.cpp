//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file dc.cpp
//  \brief piecewise constant (donor cell) reconstruction

// Athena++ headers
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../cr.hpp"
#include "./cr_integrators.hpp"

//----------------------------------------------------------------------------------------
//! \fn ReIntegrator::CRFlux()
//  \brief HLLE flux for Cosmic Ray Transport

void CRIntegrator::CRFlux(int fdir, int il, int iu,
      AthenaArray<Real> &w_l, AthenaArray<Real> &w_r,
      AthenaArray<Real> &vdiff_l, AthenaArray<Real> &vdiff_r,
                                      AthenaArray<Real> &flx) {
  // First, get the photon diffusion speed
  Real vmax = pmy_cr->vmax;
  Real eddl = 1.0/3.0;
  Real eddr = 1.0/3.0;
  Real edd12 = 0.0;
  Real edd13 = 0.0;
  Real edd23 = 0.0;
  // use sigma_l, sigma_r to get diffusion speed
  for (int i=il; i<=iu; ++i) {
    Real vl = w_l(NCR,i);
    Real vr = w_r(NCR,i);

    Real meanadv=0.5*(vl + vr);
    Real meandiffv = 0.5*(vdiff_l(i)+vdiff_r(i));

    Real al = std::min((meanadv - meandiffv),(vl-vdiff_l(i)));
    Real ar = std::max((meanadv + meandiffv),(vr+vdiff_r(i)));
    ar = std::min(ar,vmax * std::sqrt(eddr));
    al = std::max(al,-vmax * std::sqrt(eddl));

    Real bp = ar > 0.0 ? ar : 0.0;
    Real bm = al < 0.0 ? al : 0.0;

    // computer L/R fluxes along lines
    // F_L - (S_L)U_L
    // F_R - (S_R)U_R

    Real fl_e = vmax * w_l(fdir,i) - bm * w_l(CRE,i);
    Real fr_e = vmax * w_r(fdir,i) - bp * w_r(CRE,i);
    Real fl_f1, fr_f1, fl_f2, fr_f2, fl_f3, fr_f3;
    if (fdir == CRF1) {
      fl_f1 = vmax * eddl * w_l(CRE,i) - bm * w_l(CRF1,i);
      fr_f1 = vmax * eddr * w_r(CRE,i) - bp * w_r(CRF1,i);

      fl_f2 = vmax * edd12 * w_l(CRE,i) - bm * w_l(CRF2,i);
      fr_f2 = vmax * edd12 * w_r(CRE,i) - bp * w_r(CRF2,i);

      fl_f3 = vmax * edd13 * w_l(CRE,i) - bm * w_l(CRF3,i);
      fr_f3 = vmax * edd13 * w_r(CRE,i) - bp * w_r(CRF3,i);
    } else if (fdir == CRF2) {
      fl_f1 = vmax * edd12 * w_l(CRE,i) - bm * w_l(CRF1,i);
      fr_f1 = vmax * edd12 * w_r(CRE,i) - bp * w_r(CRF1,i);

      fl_f2 = vmax * eddl * w_l(CRE,i) - bm * w_l(CRF2,i);
      fr_f2 = vmax * eddr * w_r(CRE,i) - bp * w_r(CRF2,i);

      fl_f3 = vmax * edd23 * w_l(CRE,i) - bm * w_l(CRF3,i);
      fr_f3 = vmax * edd23 * w_r(CRE,i) - bp * w_r(CRF3,i);
    } else {
      fl_f1 = vmax * edd13 * w_l(CRE,i) - bm * w_l(CRF1,i);
      fr_f1 = vmax * edd13 * w_r(CRE,i) - bp * w_r(CRF1,i);

      fl_f2 = vmax * edd23 * w_l(CRE,i) - bm * w_l(CRF2,i);
      fr_f2 = vmax * edd23 * w_r(CRE,i) - bp * w_r(CRF2,i);

      fl_f3 = vmax * eddl * w_l(CRE,i) - bm * w_l(CRF3,i);
      fr_f3 = vmax * eddr * w_r(CRE,i) - bp * w_r(CRF3,i);
    }
    // calculate the HLLE flux
    Real tmp = 0.0;
    if (std::abs(bm-bp) > TINY_NUMBER) {
      tmp = 0.5*(bp + bm)/(bp - bm);
    }
    flx(CRE,i) = 0.5*(fl_e + fr_e) + (fl_e - fr_e) * tmp;
    flx(CRF1,i) = 0.5*(fl_f1 + fr_f1) + (fl_f1 - fr_f1) * tmp;
    flx(CRF2,i) = 0.5*(fl_f2 + fr_f2) + (fl_f2 - fr_f2) * tmp;
    flx(CRF3,i) = 0.5*(fl_f3 + fr_f3) + (fl_f3 - fr_f3) * tmp;
  }
  return;
}
