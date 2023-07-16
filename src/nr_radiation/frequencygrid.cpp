//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file frequencygrid.cpp
//  \brief implementation of frequency grid in class Radiation
//======================================================================================

// C headers

// C++ headers
#include <cstdio>  // fopen and fwrite
#include <sstream>  // msg

// Athena++ headers
#include "../utils/utils.hpp"
#include "./radiation.hpp"

//--------------------------------------------------------------------------------------
// \!fn void FrequencyGrid()

// \brief function to create the frequency grid
// specific intensities are still defined as frequency integrated over each

void NRRadiation::FrequencyGrid() {
  Real h_planck = 6.6260755e-27; // Planck constant
  Real k_b = 1.380649e-16;   // Boltzman constant

  // convert frequency to unit of kT_unit/h
  if (nu_min < 0.0) {
    nu_min = -nu_min;
  } else {
    nu_min = nu_min * h_planck/(k_b * tunit);
  }
  if (nu_max < 0.0) {
    nu_max = -nu_max;
  } else {
    nu_max = nu_max * h_planck/(k_b * tunit);
  }

  if ((nu_max <= nu_min) && (nfreq > 2)) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Radiation Class" << std::endl
        << "frequency_max needs to be larger than frequency_min!";
    throw std::runtime_error(msg.str().c_str());
  }

  if (log_fre_ == 1) {
    if (fre_ratio > 1) {
      if (nfreq > 2) {
        nu_max = nu_min * std::pow(fre_ratio,nfreq-2);
      } else {
        nfreq = std::log10(nu_max/nu_min)/std::log10(fre_ratio)+2;
      }
    } else {
      if (nfreq > 2) {
        fre_ratio = std::log10(nu_max/nu_min)/(nfreq - 2);
        fre_ratio = std::pow(10.0,fre_ratio);
      }
    }
  }

  if (nfreq > 1) {
    nu_grid.NewAthenaArray(nfreq);
    nu_grid(0) = 0.0;
    nu_grid(1) = nu_min;

    if (nfreq > 2) {
      if (log_fre_ == 1) {
        for(int n=2; n<nfreq; ++n)
          nu_grid(n) = nu_grid(n-1) * fre_ratio;
      } else {
        Real d_nu = (nu_max - nu_min)/(nfreq-2);

        for(int n=2; n<nfreq; ++n)
          nu_grid(n) = nu_grid(n-1) + d_nu;
      }
    }// end nfreq > 2

    nu_cen.NewAthenaArray(nfreq);

    for(int n=0; n<nfreq-1; ++n)
      nu_cen(n) = 0.5*(nu_grid(n)+nu_grid(n+1));

    delta_nu.NewAthenaArray(nfreq);

    for(int n=0; n<nfreq-1; ++n) {
      delta_nu(n) = nu_grid(n+1)-nu_grid(n);
    }
  }
  emission_spec.NewAthenaArray(nfreq);

  // initialize with default emission spectrum assuming tgas=1
  if (nfreq == 1) {
    emission_spec(0) = 1.0;
  } else {
    for (int ifr=0; ifr<nfreq-1; ++ifr) {
      emission_spec(ifr) =  IntPlanckFunc(nu_grid(ifr), nu_grid(ifr+1));
    }
      emission_spec(nfreq-1) = 1.0 - FitIntPlanckFunc(nu_grid(nfreq-1));
  }
}

// return the integral (15/pi^4)\int_0^{nu/T} x^3 dx/(exp(x)-1)
// frequency is scaled with kT_0/h
// using fitting formula to return \int_0^nu_min and \int_0^nu_max
// The fitting formula is based on asymptotic expression of the integral,
// as well as polynomial fitting to the numerical values calculated by
// directly integrating the function in mathematica
Real NRRadiation::FitIntPlanckFunc(Real nu_t) {
  // the integral at nu_t=1.5 is 0.6154949828394710
  Real integral = 0.0;
  Real nu_2 = nu_t * nu_t;
  Real nu_3 = nu_t * nu_2;
  Real nu_7 = nu_2 * nu_2 * nu_3;
  if (nu_t < 1.9434) {
    integral = 0.051329911273422 * nu_3 -0.019248716727533 * nu_t * nu_3
               + 0.002566495563671 * nu_2 * nu_3
               -3.055351861513195*1.e-5*nu_7;
  } else if (nu_t < 5.23754) {
    Real logfre = std::log10(nu_t);
    integral = -0.6874*logfre*logfre*logfre
               -0.5698*logfre*logfre+2.671*logfre-1.476;
    integral = std::pow(10.0,integral);
  } else if (nu_t < 50) {
    integral = 1.0-15*ONE_PI_FOUR_POWER*std::exp(-nu_t) *
               (6.0+6.0*nu_t+3.0*nu_t*nu_t+nu_t*nu_t*nu_t);
  } else {
    integral = 1.0;
  }
  return integral;
}


Real NRRadiation::IntPlanckFunc(Real nu_min, Real nu_max) {
  return (FitIntPlanckFunc(nu_max) - FitIntPlanckFunc(nu_min));
}

// In the last frequency bin, [nu, infty]
// we assume the spectrum is blackbody with effective temeprature Tr
// so that intensity=T_r^4 (15/pi^4) int_{nu/T_r}^{infty} x^3dx/(exp(x)-1)
// This is rearranged to
// intensity/nu^4=A=(1/y)^4(15/pi^4)\int_y^{infty} x^3dx/(exp(x)-1)
// we use a fitting formula to get y
Real NRRadiation::EffectiveBlackBody(Real intensity, Real nu) {
  Real ir = std::max(intensity,static_cast<Real>(TINY_NUMBER));
  Real a_nu = ir/(nu*nu*nu*nu); // I/nu^4
  Real nu_tr = 1.0;
  if (a_nu > 0.184077200146896) {
    // solve the fourth order polynomial (1/y)^4-(5/pi^4)(1/y)-a_nu=0
    // -(pi^4/5) (1/y)^4 + (1/y) + a_nu(pi^4/5)==0
    Real coef4=-PI_FOUR_POWER*0.2;
    Real coef = a_nu*0.2*PI_FOUR_POWER;
    int flag=FouthPolyRoot(coef4, coef, nu_tr);
    if (flag == -1)
      nu_tr = std::pow((1.0/a_nu),0.25);
    else
      nu_tr = 1.0/nu_tr;
  } else {
    Real loganu = -std::log(a_nu);
    if (loganu < 40.0) {
      nu_tr = -0.000525 * loganu * loganu * loganu + 0.03138 * loganu * loganu
              + 0.3223 * loganu + 0.8278;
    } else {
      nu_tr = 30.3278;
    }

    if (a_nu < 1.e-10) {
      // improve the accuracy with iteration
      Real yini = nu_tr;
      int count = 0;
      Real residual=1.0;
      while((count < 6) && (residual > 1.e-6)) {
        Real ratio = yini*yini*yini*yini/(6.0+6.0*yini+3.0*yini*yini+yini*yini*yini);
        nu_tr = -std::log((PI_FOUR_POWER*a_nu/15.0)*ratio);
        residual = std::abs((nu_tr-yini)/nu_tr);
        count++;
        yini=nu_tr;
      }
    }
  }
  return nu_tr;
}

// In the last frequency bin, [nu, infty]
// we assum, n(nu)=1/(exp(nu/T)-1), the input value
// n_nu2 = \int_{nu_f}^{infty} n\nu^2 d\nu
// we fit the formula n_nu2/nu^3=A=(1/y)^3\int_y^{infty} x^2/(exp(x)-1) dx

Real NRRadiation::EffectiveBlackBodyNNu2(Real n_nu2, Real nu) {
  Real fit_a = n_nu2/(nu*nu*nu);
  Real log_fit_a=std::log(fit_a);
  Real nu_tr = 1.0;
  if (fit_a < 0.001) {
    nu_tr = 0.001177*log_fit_a*log_fit_a-0.8812*log_fit_a-0.7428;
  } else if (fit_a < 5.0) {
    Real exp_nu = (log_fit_a+13.28)/9.551;
    nu_tr=8.647*std::exp(-exp_nu*exp_nu);
  } else {
    nu_tr=std::pow(2.404113806319301/fit_a,1.0/3.0);
  }
  return nu_tr;
}

// return the integral (15/pi^4)\int_{\nu/T}^{\infty} \nu J_nu d\nu
// =(15/pi^4)\int_{\nu/T}^{\infty} x^4/(exp(x)-1) dx
// the input is nu_t=nu_f/T
Real NRRadiation::IntegrateBBNuJ(Real nu_t) {
  Real nu_sq = nu_t*nu_t;
  Real nu_four = nu_sq * nu_sq;
  Real nu_three = nu_t * nu_sq;
  Real nu_five = nu_t * nu_four;
  Real nu_six = nu_sq * nu_four;
  Real exp_nu = std::exp(-nu_t);
  Real integral = 0.0;
  if (nu_t < 2.596) {
    integral = 3.832229496128511
               - 3.75 * ONE_PI_FOUR_POWER * nu_four
               + 1.5 * ONE_PI_FOUR_POWER * nu_five
               - (5.0/24.0) * ONE_PI_FOUR_POWER * nu_six;
  } else {
    integral = 15.0*ONE_PI_FOUR_POWER*exp_nu*(24.0
               +24.0*nu_t+12.0*nu_sq+4.0*nu_three+nu_four);
  }
  return integral;
}

// return the integral (15/pi^4)\int_{\nu/T}^{\infty} (J_nu/nu)^2 d\nu
// =(15/pi^4)\int_{\nu/T}^{\infty} x^4/(exp(x)-1)^2 dx
// the input is nu_t=nu_f/T
Real NRRadiation::IntegrateBBJONuSq(Real nu_t) {
  Real nu_sq = nu_t*nu_t;
  Real nu_four = nu_sq * nu_sq;
  Real nu_three = nu_t * nu_sq;
  Real exp_nu = std::exp(-2.0*nu_t);
  Real integral = 0.0;
  if (nu_t < 10.0) {
    Real top = 0.01872 * nu_sq - 0.2732 * nu_t + 0.9735;
    Real bottom = nu_sq - 1.854 * nu_t + 5.828;
    integral = top/bottom;
  } else {
    integral = 3.75*ONE_PI_FOUR_POWER*exp_nu*(3.0
               +6.0*nu_t+6.0*nu_sq+4.0*nu_three+2.0*nu_four);
  }
  return integral;
}

//For a given value J=(15/pi^4)T_r^4\int x^3/(exp(x)-1) dx,
// get the integral \int (J/\nu)^2 dx=
// (15/pi^4)^2 T_r^5\int x^4/(exp(x)-1)^2 dx
// The input is J and nu_f

Real NRRadiation::BBJToJONuSq(Real &bb_j, Real &nu_f) {
  Real j_nu4 = bb_j/(nu_f*nu_f*nu_f*nu_f);
  Real log_j_nu4 = std::log10(j_nu4);
  Real log_jonu2 = 0.0;
  if (log_j_nu4 < -2.71777) {
    log_jonu2 = (1.987*log_j_nu4*log_j_nu4-2.123*log_j_nu4+4.532)/
       (log_j_nu4-1.733);
  } else if ((log_j_nu4 >= -2.71777) && (log_j_nu4 < 2.66367)) {
    log_jonu2 = 0.006977*log_j_nu4*log_j_nu4*log_j_nu4
         -0.03641*log_j_nu4*log_j_nu4+1.314*log_j_nu4-1.632;
  } else {
    log_jonu2 = 1.25*log_j_nu4-1.588;
  }

  Real jonu2 = std::pow(10.0, log_jonu2);
  jonu2 = jonu2 * nu_f*nu_f*nu_f*nu_f*nu_f;
  return jonu2;
}

// For a given value J=(15/pi^4)T_r^4\int x^3/(exp(x)-1) dx,
// return the value n(nu_f)=1/(exp(x)-1)
// The input is J and nu_f

Real NRRadiation::BBJtoNnu(Real &bb_j, Real &nu_f) {
  Real j_nu4 = bb_j/(nu_f*nu_f*nu_f*nu_f);
  Real log_j_nu4 = std::log10(j_nu4)+3.0;
  Real log_nnu = 0.0;

  if (log_j_nu4 < 0.0) {
    log_nnu = -0.003321*log_j_nu4*log_j_nu4+0.8649*log_j_nu4-1.869;
  } else if (log_j_nu4 >= 0.0 && log_j_nu4 < 7.07128) {
    log_nnu = 0.004199*log_j_nu4*log_j_nu4*log_j_nu4
            -0.07701*log_j_nu4*log_j_nu4+0.7398*log_j_nu4-1.869;
  } else {
    log_nnu = -0.002325*log_j_nu4*log_j_nu4+0.2955*log_j_nu4-0.977;
  }

  Real nnu = std::pow(10.0, log_nnu);
  return nnu;
}

// For a given value J=(15/pi^4)T_r^4\int x^3/(exp(x)-1) dx,
// return the value Jnu=(15/pi^4)T_r^5\int x^4/(exp(x)-1) dx
// The input is J and nu_f

Real NRRadiation::BBJtoJnu(Real &bb_j, Real &nu_f) {
  Real j_nu4 = bb_j/(nu_f*nu_f*nu_f*nu_f);
  Real log_j_nu4 = std::log10(j_nu4);
  Real log_jnu = 0.0;
  if (log_j_nu4 < -3.1936) {
    log_jnu = 6.363e-5 * log_j_nu4*log_j_nu4*log_j_nu4
              +0.002858*log_j_nu4*log_j_nu4
              +1.043*log_j_nu4+0.2264;
  } else if ((log_j_nu4 >= -3.1936) && (log_j_nu4 < 2.6701)) {
    log_jnu = -0.003183*log_j_nu4*log_j_nu4*log_j_nu4
              +0.0158*log_j_nu4*log_j_nu4
              +1.225*log_j_nu4+0.5983;
  } else {
    log_jnu = 1.25*log_j_nu4+0.5836;
  }

  Real jnu = std::pow(10.0, log_jnu);
  jnu *= (nu_f*nu_f*nu_f*nu_f*nu_f);

  return jnu;
}

Real NRRadiation::DBBjDNNu2(Real &bb_j, Real &nu_f) {
  Real j_nu4 = bb_j/(nu_f*nu_f*nu_f*nu_f);
  Real log_j_nu4 = std::log10(j_nu4);
  Real log_ratio=0.0;
  if (log_j_nu4 < -3.06556) {
    log_ratio=(-0.8125*log_j_nu4*log_j_nu4-1.888*log_j_nu4-2.331)/
              (log_j_nu4*log_j_nu4+2.091*log_j_nu4+3.25);
  } else if (log_j_nu4 >= -3.06556 && log_j_nu4 < 2.87328) {
    log_ratio=(0.2166*log_j_nu4*log_j_nu4*log_j_nu4
              +1.491*log_j_nu4*log_j_nu4+4.095*log_j_nu4-5.699)/
              (log_j_nu4*log_j_nu4+6.111*log_j_nu4+24.98);
  } else {
    log_ratio=0.2499*log_j_nu4-0.2551;
  }

  Real ratio = std::pow(10.0, log_ratio);
  ratio *= nu_f;
  return ratio;
}

// return the integral \int_{nu/T}^{\infty} n \nu^2 d\nu
// = \int_{nu/T}^{\infty} x^2/(exp(x)-1) dx
// the input is nu_t=nu_f/T_r
Real NRRadiation::IntegrateBBNNu2(Real nu_t) {
  Real nu_sq = nu_t*nu_t;
  Real nu_three = nu_t * nu_sq;
  Real nu_four = nu_sq * nu_sq;
  Real nu_six = nu_three * nu_three;
  Real nu_eight = nu_four * nu_four;
  Real exp_nu = std::exp(-nu_t);
  Real integral = 0.0;
  if (nu_t < 3.297) {
    integral = 2.404113806319301 - 0.5 * nu_sq + (1.0/6.0) * nu_three
             - (1.0/48.0) * nu_four + (1.0/4320.0) * nu_six
             - (1.0/241920.0) * nu_eight;
  } else {
    integral = exp_nu * (2.0 + 2.0 * nu_t + nu_sq);
  }
  return integral;
}

// We assume blackbody spectrum in the frequency range nu_f to infty
// For a given value of J=(15/\pi^4)\int_{\nu_f}^infty nu^3/(exp(nu/Tr)-1) dnu
// and frequency nu_f, we return the integral
// nnu2=\int_{nu_f}^{infty} nu^2/(exp(nu/Tr)-1) d\nu
// We use piecewise linear relationship between J and nnu2
Real NRRadiation::ConvertBBJNNu2(Real &bb_j, Real &nu_f) {
  Real jnu4 = bb_j/(nu_f*nu_f*nu_f*nu_f);
  Real log_jnu4=std::log10(jnu4);
  Real log_nnu2=1.0;
  if (log_jnu4 < -3.99793) {
    log_nnu2 = 0.9971*log_jnu4+0.7464;
  } else if (log_jnu4 >= -3.99793 && log_jnu4 < 2.73043) {
    Real jsq = log_jnu4*log_jnu4;
    log_nnu2 = 0.0003484*jsq*jsq+0.002327*jsq*log_jnu4-0.02137*jsq
               +0.8037*log_jnu4+0.3249;
  } else {
    log_nnu2 = 0.7506*log_jnu4+0.3773;
  }
  Real nnu2 = std::pow(10.0,log_nnu2)*nu_f*nu_f*nu_f;
  return nnu2;
}

// inverse conversion from nnu2 to J.
// The conversion from J to nnu2, and then from nnu2 to J
// Should be exact
Real NRRadiation::InverseConvertBBJNNu2(Real &nnu2, Real &nu_f) {
  Real nnu2_nu3=nnu2/(nu_f*nu_f*nu_f);
  Real log_nnu2=std::log(nnu2_nu3);
  Real log_jnu4=1.0;
  if (log_nnu2 > 0.746) {
    log_jnu4=(log_nnu2-0.746)/0.7644;
  } else {
    log_jnu4=(log_nnu2-0.746)/0.9846;
  }
  Real bbj = std::exp(log_jnu4)*nu_f*nu_f*nu_f*nu_f;

  return bbj;
}

// convert J to nnu2 assuming Wien profile spectrum

void NRRadiation::ConvertBBJWien(Real &bb_j, Real &nu_f, Real &tgas,
                              Real &nuj, Real &jonusq) {
  Real nu_t = nu_f/tgas;

  // return \int \nu Jd\nu
  Real j_coef = 6.0+6.0*nu_t+3.0*nu_t*nu_t+nu_t*nu_t*nu_t;
  Real nj_coef = 24.0 + 24.0*nu_t + 12.0*nu_t*nu_t
               + 4.0*nu_t*nu_t*nu_t + nu_t*nu_t*nu_t*nu_t;
  nuj = bb_j*tgas*nj_coef/j_coef;


  // return \int (J/\nu)^2 d\nu
  Real jonu2_coef = 0.25*(3.0+6.0*nu_t+6.0*nu_t*nu_t
                   +4.0*nu_t*nu_t*nu_t+2.0*nu_t*nu_t*nu_t*nu_t);
  Real bbj_ratio = bb_j/(tgas*tgas*j_coef);
  jonusq = bbj_ratio*bbj_ratio*tgas*jonu2_coef;
  return;
}

void NRRadiation::ConvertBBJWien2(Real &bb_j, Real &nu_f, Real &tgas,
                                  Real &nnu2, Real &n_nuf) {
  Real nu_t = nu_f/tgas;
  // return \int n\nu^2d\nu
  Real n_coef = 2.0+2.0*nu_t+nu_t*nu_t;
  Real j_coef = 6.0+6.0*nu_t+3.0*nu_t*nu_t+nu_t*nu_t*nu_t;
  nnu2 = bb_j*PI_FOUR_POWER*n_coef/(15.0*tgas*j_coef);
  // return n(nu_f) = exp(-nu_f/T)/\lambda
  Real tgas_four = tgas*tgas*tgas*tgas;
  n_nuf = bb_j*PI_FOUR_POWER/(15.0*tgas_four*j_coef);
  return;
}

// convert nnu2 to J assuming Wien profile spectrum

Real NRRadiation::InverseConvertBBJNNu2Wien(Real &nnu2, Real &nu_f, Real &tgas) {
  Real nu_t = nu_f/tgas;
  Real n_coef = 2.0+2.0*nu_t+nu_t*nu_t;
  Real j_coef = 6.0+6.0*nu_t+3.0*nu_t*nu_t+nu_t*nu_t*nu_t;

  Real bb_j = nnu2*tgas*15.0*ONE_PI_FOUR_POWER*j_coef/n_coef;

  return bb_j;
}
