//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file emission.cpp
//  \brief coupling of radiation to matter via optically thin emission

// Athena++ headers
#include "../radiation.hpp"
#include "../../athena_arrays.hpp"  // AthenaArray
#include "../../mesh/mesh.hpp"
#include "../../eos/eos.hpp"
#include "../../athena.hpp"

//----------------------------------------------------------------------------------------
// Function for calculating effect on radiation from optically thin emission
// Inputs:
//   prim_hydro: primitive hydro variables:
//     index 0: density (IDN) or pressure (IPR)
//     indices 1-3: k, j, i
//   n: unit normal null vector:
//     index 0: fluid-frame component (0 through 3)
//     index 1: angle (0 through nzeta * npsi - 1)
//     index 2: i
//   omega: fractional solid angle (normalized to 1) in fluid frame:
//     index 0: angle (0 through nzeta * npsi - 1)
//     index 1: i
//   dt: time interval over which coupling should be applied
//   k, j: x3- and x2-indices
//   intensity: fluid-frame I:
//     index 0: angle (0 through nzeta * npsi - 1)
//     index 1: i
// Outputs:
//   intensity: fluid-frame I updated
// Notes:
//   Implements emission with no absorption or scattering.

int FouthPolyRoot(const Real coef4, const Real tconst, Real &root);

void Radiation::Coupling(const AthenaArray<Real> &prim_hydro, const AthenaArray<Real> &n,
    const AthenaArray<Real> &omega, Real dt, int k, int j, AthenaArray<Real> &intensity) {

  Real gamma = pmy_block->peos->GetGamma();
  int nang = nzeta * npsi;


  Real coef[2];
  for (int i=0; i<2; ++i)
    coef[i] = 0.0;

  //First, swap the order of intensity(n,i) to (i,n)
  for(int i=is; i<=ie; ++i){
    Real rho = prim_hydro(IDN,k,j,i);
    Real pgas = prim_hydro(IPR,k,j,i);
    Real tgas = pgas/rho;
    Real tgasnew = tgas;
  
    for(int m=0; m< nang; ++m){
      intensity_scr_(m) = intensity(m,i);
      tran_coef_(m) = n(0,m,i);
      weight_(m) = omega(m,i);
    }

    Real suma1=0.0, suma2=0.0, suma3=0.0;
    Real jr_cm=0.0;

    Real sigma_s=opacity(OPAS,k,j,i);
    Real sigma_a=opacity(OPAA,k,j,i);
    Real sigma_p=0.0;
    
    if(using_planck_mean){
      sigma_p = opacity(OPAP,k,j,i);
    }

    Real dtcsigmaa = dt * sigma_a;
    Real dtcsigmas = dt * sigma_s;
    Real dtcsigmap = dt * sigma_p;

    bool badcell = 0;
    
    for(int m=0; m<nang; m++){
       Real vncsigma = 1.0/(1.0 + (dtcsigmaa + dtcsigmas) * tran_coef_(m));
       vncsigma2_(m) = tran_coef_(m) * vncsigma;
       Real ir_weight = intensity_scr_(m) * weight_(m);
       jr_cm += ir_weight;
       suma1 += (weight_(m) * vncsigma2_(m));
       suma2 += (ir_weight * vncsigma);
    }
    suma3 = suma1 * (dtcsigmas - dtcsigmap);
    suma1 *= (dtcsigmaa + dtcsigmap);

    coef[1] = (dtcsigmaa + dtcsigmap - (dtcsigmaa + dtcsigmap) * suma1/(1.0-suma3))
                   * (gamma - 1.0)/rho;
    coef[0] = -tgas - (dtcsigmaa + dtcsigmap) * suma2 * (gamma - 1.0)/(rho*(1.0-suma3));
    
    if(fabs(coef[1]) > TINY_NUMBER){
      int flag = FouthPolyRoot(coef[1], coef[0], tgasnew);
      if(flag == -1 || tgasnew != tgasnew){
        badcell = 1;
        tgasnew = tgas;
      }
    }else{
      tgasnew = -coef[0];
    }
    // even if tr=told, there can be change for intensity, making them isotropic
    if(!badcell){
    
      Real emission = arad* tgasnew * tgasnew * tgasnew * tgasnew;
      
      // get update jr_cm
      jr_cm = (suma1 * emission + suma2)/(1.0-suma3);
    
    // Update the co-moving frame specific intensity

      for(int m=0; m<nang; m++){
        intensity_scr_(m) +=
                         ((dtcsigmas - dtcsigmap) * jr_cm + (dtcsigmaa + dtcsigmap) * emission
                            - (dtcsigmas + dtcsigmaa) * intensity_scr_(m)) * vncsigma2_(m);
      }

      // copy intensity back 
      for(int m=0; m< nang; ++m){
        intensity(m,i) = intensity_scr_(m);
      }
    }// finish badcell

  }// finish i loop



  return;
}




// Exact solution for fourth order polynomical with the format
// coef4 * x^4 + x + tconst == 0

int FouthPolyRoot(const Real coef4, const Real tconst, Real &root)
{

// First, get the real root of
// z^3-4*tconst/coef4 * z - 1/coef4^2==0
  Real asquar = coef4 * coef4;
  Real acubic = coef4 * asquar;
  Real ccubic = tconst * tconst * tconst;
  Real delta1 = 0.25 - 64.0 * ccubic * coef4/27.0;
  if(delta1 < 0.0) return -1;
  else delta1 = sqrt(delta1);
  if(delta1 < 0.5) return -1;
  Real zroot = 0.0;
  if(delta1 > 1.e11){
    // to avoid small number cancellation
    zroot = pow(delta1,-2.0/3.0)/3.0;
  }else{
    zroot = (pow(0.5 + delta1, 1.0/3.0) -
               pow(-0.5 + delta1, 1.0/3.0));
  }

  if(zroot < 0.0) return -1;

  zroot *= pow(coef4,-2.0/3.0);
  
  Real rcoef = sqrt(zroot);
  Real delta2 = -zroot + 2.0/(coef4*rcoef);
  if(delta2 < 0.0) return -1;
  else delta2 = sqrt(delta2);
  root = 0.5 * (delta2 - rcoef);
  if(root < 0.0) return -1;

  return 0;
}
