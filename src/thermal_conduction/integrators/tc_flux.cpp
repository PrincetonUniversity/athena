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
#include "../../mesh/mesh.hpp"
#include "../tc.hpp"
#include "tc_integrators.hpp"


//----------------------------------------------------------------------------------------
//! \fn Reconstruction::NTCFlux()
//  \brief HLLE flux for Thermal Conduction

//The four independent variables are:
// e, (T/e)F_1, (T/e)F_2, (T/e)F_3
// But the stored variables are
// e, F_1, F_2, F_3
//the fluxes are:
// V_m (F_1, F_2, F_3)
// V_m T * (Unity Tensor), all off-diagonal components are zero

void TCIntegrator::TCFlux(int fdir, int il, int iu, 
      AthenaArray<Real> &w_l, AthenaArray<Real> &w_r,  
      AthenaArray<Real> &vdiff_l, AthenaArray<Real> &vdiff_r,      
                                      AthenaArray<Real> &flx)  
{

  Real vmax = pmy_tc->vmax;
  for (int i=il; i<=iu; ++i){
    
    Real meandiffv = 0.5*(vdiff_l(i)+vdiff_r(i));
    Real rho_l = w_l(TCRHO,i);
    Real rho_r = w_r(TCRHO,i);
    Real t_l = w_l(TCT,i);
    Real t_r = w_r(TCT,i);
    Real meanrho = 0.5*(rho_l + rho_r);
    Real rhoratiol=meanrho/rho_l;
    Real rhoratior=meanrho/rho_r;


    Real al = std::min(meandiffv,vdiff_l(i));
    Real ar = std::max(meandiffv,vdiff_r(i));


    Real bp = ar > 0.0 ? ar : 0.0;
    Real bm = al < 0.0 ? al : 0.0;

    // computer L/R fluxes along lines
    // F_L - (S_L)U_L
    // F_R - (S_R)U_R

    Real fl_e = vmax * w_l(fdir,i) - bm * w_l(0,i)*rhoratiol;
    Real fr_e = vmax * w_r(fdir,i) - bp * w_r(0,i)*rhoratior;
    Real fl_f1, fr_f1, fl_f2, fr_f2, fl_f3, fr_f3;
    if(fdir == 1){

      fl_f1 = vmax * t_l - bm * w_l(1,i)/meanrho;
      fr_f1 = vmax * t_r - bp * w_r(1,i)/meanrho;

//      fl_f2 = -bm * w_l(2,i);
//      fr_f2 = -bp * w_r(2,i);

//      fl_f3 = -bm * w_l(3,i);
//      fr_f3 = -bp * w_r(3,i);

      fl_f2 = 0.0;
      fr_f2 = 0.0;

      fl_f3 = 0.0;
      fr_f3 = 0.0;

    } else if(fdir == 2){
 
//      fl_f1 = -bm * w_l(CRF1,i);
//      fr_f1 = -bp * w_r(CRF1,i);
      fl_f1 = 0.0;
      fr_f1 = 0.0;

      fl_f2 = vmax * t_l - bm * w_l(2,i)/meanrho;
      fr_f2 = vmax * t_r - bp * w_r(2,i)/meanrho;

//      fl_f3 = -bm * w_l(CRF3,i);
//      fr_f3 = -bp * w_r(CRF3,i);

      fl_f3 = 0.0;
      fr_f3 = 0.0;

    }else if(fdir == 3){

//      fl_f1 = vmax * eddl(PC13,i) * w_l(CRE,i) - bm * w_l(CRF1,i);
//      fr_f1 = vmax * eddr(PC13,i) * w_r(CRE,i) - bp * w_r(CRF1,i);

//      fl_f2 = vmax * eddl(PC23,i) * w_l(CRE,i) - bm * w_l(CRF2,i);
//      fr_f2 = vmax * eddr(PC23,i) * w_r(CRE,i) - bp * w_r(CRF2,i);
      fl_f1 = 0.0;
      fr_f1 = 0.0;

      fl_f2 = 0.0;
      fr_f2 = 0.0;

      fl_f3 = vmax * t_l - bm * w_l(3,i)/meanrho;
      fr_f3 = vmax * t_r - bp * w_r(3,i)/meanrho;

    }
  
    //calculate the HLLE flux
    Real tmp = 0.0;
    if(fabs(bm-bp) > TINY_NUMBER)
    	tmp = 0.5*(bp + bm)/(bp - bm);
    flx(0,i) = 0.5*(fl_e + fr_e) + (fl_e - fr_e) * tmp;
    flx(1,i) = 0.5*(fl_f1 + fr_f1) + (fl_f1 - fr_f1) * tmp;
    flx(2,i) = 0.5*(fl_f2 + fr_f2) + (fl_f2 - fr_f2) * tmp;
    flx(3,i) = 0.5*(fl_f3 + fr_f3) + (fl_f3 - fr_f3) * tmp;       

  }


  return;
}


