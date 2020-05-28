//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file awa_z4c.cpp
//  \brief implementation of functions in the Z4c class for initializing Apples with Apples tests
//  See  https://arxiv.org/abs/gr-qc/0305023
//       https://arxiv.org/abs/1111.2177

// C++ standard headers
#include <cmath> // pow, rand, sin, sqrt
#include <ctime>
#include <random>
#include <iostream>
#include <fstream>

// Random number in [-1,1]
std::default_random_engine generator{137};
std::uniform_real_distribution<double> distribution(-1.,1.);
#define RANDOMNUMBER (distribution(generator))

// Sin wave for various wave tests
#define  SINWAVE(a,dx,dy,x,y) ( (a)*        std::sin(2.*M_PI*((x)/(dx) - (y)/(dy))) )
#define DSINWAVE(a,dx,dy,x,y) (-(a)*2.*M_PI*std::cos(2.*M_PI*((x)/(dx) - (y)/(dy))) )

// Gaussian profile for simple 3D gauge wave test
#define  GAUSSIAN(a,d,x,y,z) (           (a)*std::exp(-(SQR(x)+SQR(y)+SQR(z))/SQR(d)) )

// Gaussian profile for 1d linear test
#define   GAUSSIAN1(a,w,x) (               (a)*std::exp(-SQR(x)/(2 * SQR(w))) )
#define  DGAUSSIAN1(a,w,x) (  (a * x / SQR(W))*std::exp(-SQR(x)/(2 * SQR(w))) )


// Athena++ headers
#include "z4c.hpp"
#include "z4c_macro.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"

// External libraries
#ifdef GSL
#include <gsl/gsl_sf_bessel.h>   // Bessel functions
#endif

//----------------------------------------------------------------------------------------
// \!fn void Z4c::ADMRobustStability(AthenaArray<Real> & u)
// \brief Initialize ADM vars for robust stability test

// Note the amplitude of the noise (~1e-10) should be also rescaled by
// the square of the rho parameter

void Z4c::ADMRobustStability(AthenaArray<Real> & u_adm) {
  ADM_vars adm;
  SetADMAliases(u_adm, adm);

  Real const amp = opt.AwA_amplitude;

  // Flat spacetime
  ADMMinkowski(u_adm);

  GLOOP2(k,j) {
    // g_ab
    for(int a = 0; a < NDIM; ++a)
      for(int b = a; b < NDIM; ++b) {
        GLOOP1(i) {
          adm.g_dd(a,b,k,j,i) += RANDOMNUMBER * amp;
        }
      }
    // K_ab
    for(int a = 0; a < NDIM; ++a)
      for(int b = a; b < NDIM; ++b) {
        GLOOP1(i) {
          adm.K_dd(a,b,k,j,i) += RANDOMNUMBER * amp;
        }
      }
  }
}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::GaugeRobStab(AthenaArray<Real> & u, bool shifted)
// \brief Initialize lapse and shift for Robust Stability test

void Z4c::GaugeRobStab(AthenaArray<Real> & u) {
  Z4c_vars z4c;
  SetZ4cAliases(u, z4c);
  z4c.alpha.Fill(1.);
  z4c.beta_u.Fill(0.);

  Real const amp = opt.AwA_amplitude;

  GLOOP2(k,j) {
    GLOOP1(i) {
      // lapse
      z4c.alpha(k,j,i) += RANDOMNUMBER * amp;
      // shift
      z4c.beta_u(0,k,j,i) += RANDOMNUMBER * amp;
      z4c.beta_u(1,k,j,i) += RANDOMNUMBER * amp;
      z4c.beta_u(2,k,j,i) += RANDOMNUMBER * amp;
    }
  }
}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::LinearWave1(AthenaArray<Real> & u)
// \brief Initialize ADM vars for linear wave test in 1d

void Z4c::ADMLinearWave1(AthenaArray<Real> & u_adm) {
  ADM_vars adm;
  SetADMAliases(u_adm, adm);

  Real const amp = opt.AwA_amplitude;
  Real const d_x = opt.AwA_d_x;
  // y = 0 and d_y should therefore be irrelevant (set to non-zero to prevent nan)
  Real const d_y = 1;

  // Flat spacetime
  ADMMinkowski(u_adm);

  // For propagation along x ...
  GLOOP2(k,j) {
    // g_yy
    GLOOP1(i) {
      adm.g_dd(1,1,k,j,i) += SINWAVE(amp, d_x, d_y, mbi.x1(i), 0.);
    }
    // g_zz
    GLOOP1(i) {
      adm.g_dd(2,2,k,j,i) -= SINWAVE(amp, d_x, d_y, mbi.x1(i), 0.);
    }
    // K_yy
    GLOOP1(i) {
      adm.K_dd(1,1,k,j,i) += 0.5*DSINWAVE(amp, d_x, d_y, mbi.x1(i), 0.);
    }
    // K_zz
    GLOOP1(i) {
      adm.K_dd(2,2,k,j,i) -= 0.5*DSINWAVE(amp, d_x, d_y, mbi.x1(i), 0.);
    }
  }
}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::LinearWave1Gaussian(AthenaArray<Real> & u)
// \brief Initialize ADM vars for Gaussian linear wave test in 1d

void Z4c::ADMLinearWave1Gaussian(AthenaArray<Real> & u_adm) {
  ADM_vars adm;
  SetADMAliases(u_adm, adm);

  Real const amp = opt.AwA_amplitude;
  Real const w = opt.AwA_Gaussian_w;

  // Flat spacetime
  ADMMinkowski(u_adm);

  // For propagation along x ...
  GLOOP2(k,j) {
    // g_yy
    GLOOP1(i) {
      adm.g_dd(1,1,k,j,i) += GAUSSIAN1(amp, w, mbi.x1(i));
    }
    // g_zz
    GLOOP1(i) {
      adm.g_dd(2,2,k,j,i) -= GAUSSIAN1(amp, w, mbi.x1(i));
    }
    // K_yy
    GLOOP1(i) {
      adm.K_dd(1,1,k,j,i) += 0.5*GAUSSIAN1(amp, w, mbi.x1(i));
    }
    // K_zz
    GLOOP1(i) {
      adm.K_dd(2,2,k,j,i) -= 0.5*GAUSSIAN1(amp, w, mbi.x1(i));
    }
  }
}
//----------------------------------------------------------------------------------------
// \!fn void Z4c::LinearWave2(AthenaArray<Real> & u)
// \brief Initialize ADM vars for linear wave test in 2d

// BD:
// Note below SQRT2 factors added to have normalization coincide with testbed
// papers.

void Z4c::ADMLinearWave2(AthenaArray<Real> & u_adm) {
  ADM_vars adm;
  SetADMAliases(u_adm, adm);

  Real const amp = opt.AwA_amplitude;
  Real const d_x = opt.AwA_d_x;
  Real const d_y = opt.AwA_d_y;

  // Flat spacetime
  ADMMinkowski(u_adm);

  GLOOP2(k,j) {
    // g_xx, g_yy
    for(int a = 0; a < NDIM-1; ++a) {
      GLOOP1(i) {
        adm.g_dd(a,a,k,j,i) += 0.5*SINWAVE(amp, d_x, d_y,
                                           mbi.x1(i), mbi.x2(j));
      }
    }

    // g_xy, g_zz
    GLOOP1(i) {
      adm.g_dd(0,1,k,j,i)  = 0.5*SINWAVE(amp, d_x, d_y,
                                         mbi.x1(i), mbi.x2(j));
      adm.g_dd(2,2,k,j,i) -=     SINWAVE(amp, d_x, d_y,
                                         mbi.x1(i), mbi.x2(j));
    }

    // K_xx, K_xy, K_yy
    for(int a = 0; a < NDIM-1; ++a) {
      for(int b = a; b < NDIM-1; ++b) {
        GLOOP1(i) {
          adm.K_dd(a,b,k,j,i) = 0.25 * SQRT2 * DSINWAVE(amp, d_x, d_y,
            mbi.x1(i), mbi.x2(j));
        }
      }
    }
    // K_zz
    GLOOP1(i) {
      adm.K_dd(2,2,k,j,i) = -0.5 * SQRT2 * DSINWAVE(amp, d_x, d_y,
        mbi.x1(i), mbi.x2(j));
    }
  }
}


//----------------------------------------------------------------------------------------
// \!fn void Z4c::GaugeSimpleGaugeWave(AthenaArray<Real> & u)
// \brief Initializes lapse for simple 3D gauge wave test

void Z4c::GaugeSimpleGaugeWave(AthenaArray<Real> & u) {
  Z4c_vars z4c;
  SetZ4cAliases(u, z4c);
  z4c.alpha.Fill(1.);
  z4c.beta_u.Fill(0.);

  Real const amp = opt.AwA_amplitude;
  Real const d_x = opt.AwA_d_x;

  GLOOP2(k,j) {
      GLOOP1(i) {
        // lapse
        // z4c.alpha(k,j,i) += GAUSSIAN(amp, d_x,
        //                              mbi.x1(i), mbi.x2(j), mbi.x3(k));
        z4c.alpha(k,j,i) += amp * pow(sin(2. * M_PI * mbi.x1(i)) , 6) *
                                  pow(sin(2. * M_PI * mbi.x2(j)) , 6) *
                                  pow(sin(2. * M_PI * mbi.x3(k)) , 6);
      }
  }
}


//----------------------------------------------------------------------------------------
// \!fn void Z4c::ADMGaugeWave1(AthenaArray<Real> & u)
// \brief Initialize ADM vars for 1D gauge wave test

void Z4c::ADMGaugeWave1(AthenaArray<Real> & u_adm) {
  ADM_vars adm;
  SetADMAliases(u_adm, adm);

  Real const amp = opt.AwA_amplitude;
  Real const d_x = opt.AwA_d_x;
  Real const d_y = opt.AwA_d_y;

  // Flat spacetime
  ADMMinkowski(u_adm);

  // Propagation along x
  GLOOP2(k,j) {
     GLOOP1(i) {
        // g_xx
        adm.g_dd(0,0,k,j,i) -= SINWAVE(amp, d_x, d_y, mbi.x1(i), 0.);
        // K_xx
        adm.K_dd(0,0,k,j,i) =
                   0.5*DSINWAVE(amp, d_x, d_y, mbi.x1(i), 0.)/
        std::sqrt(1.0 - SINWAVE(amp, d_x, d_y, mbi.x1(i), 0.));
      }
  }
}


//----------------------------------------------------------------------------------------
// \!fn void Z4c::GaugeGaugeWave1(AthenaArray<Real> & u)
// \brief Initialize lapse for 1D gauge wave test

void Z4c::GaugeGaugeWave1(AthenaArray<Real> & u) {
  Z4c_vars z4c;
  SetZ4cAliases(u, z4c);
  z4c.alpha.Fill(1.);
  z4c.beta_u.Fill(0.);

  Real const amp = opt.AwA_amplitude;
  Real const d_x = opt.AwA_d_x;
  Real const d_y = opt.AwA_d_y;

  GLOOP2(k,j) {
     GLOOP1(i) {
       // lapse
       z4c.alpha(k,j,i) = std::sqrt(1.0 - SINWAVE(amp, d_x, d_y,
                                                  mbi.x1(i), 0.));
     }
  }
}


//----------------------------------------------------------------------------------------
// \!fn void Z4c::ADMGaugeWave1_shifted(AthenaArray<Real> & u)
// \brief Initialize ADM vars for shifted 1D gauge wave test

void Z4c::ADMGaugeWave1_shifted(AthenaArray<Real> & u_adm) {
  ADM_vars adm;
  SetADMAliases(u_adm, adm);

  Real const amp = opt.AwA_amplitude;
  Real const d_x = opt.AwA_d_x;
  Real const d_y = opt.AwA_d_y;

  // Flat spacetime
  ADMMinkowski(u_adm);

  // Propagation along x ...
  GLOOP2(k,j) {
     GLOOP1(i) {
        // g_xx
        adm.g_dd(0,0,k,j,i) += SINWAVE(amp, d_x, d_y, mbi.x1(i), 0.);
        // K_xx
        adm.K_dd(0,0,k,j,i) =
                   0.5*DSINWAVE(amp, d_x, d_y, mbi.x1(i), 0.)/
        std::sqrt(1.0 + SINWAVE(amp, d_x, d_y, mbi.x1(i), 0.));
    }
  }
}


//----------------------------------------------------------------------------------------
// \!fn void Z4c::GaugeGaugeWave1_shifted(AthenaArray<Real> & u)
// \brief Initialize lapse and shift for shifted 1D gauge wave test

void Z4c::GaugeGaugeWave1_shifted(AthenaArray<Real> & u) {
  Z4c_vars z4c;
  SetZ4cAliases(u, z4c);
  z4c.alpha.Fill(1.);
  z4c.beta_u.Fill(0.);

  Real const amp = opt.AwA_amplitude;
  Real const d_x = opt.AwA_d_x;
  Real const d_y = opt.AwA_d_y;

  GLOOP2(k,j) {
     GLOOP1(i) {
        // lapse
        z4c.alpha(k,j,i) = 1.0/(std::sqrt(1.0 + SINWAVE(amp, d_x, d_y,
                                                        mbi.x1(i),0.)));
        // shift
        z4c.beta_u(0,k,j,i) = - SINWAVE(amp, d_x, d_y, mbi.x1(i),0.)/
                          (1. + SINWAVE(amp, d_x, d_y, mbi.x1(i),0.));
      }
  }
}


//----------------------------------------------------------------------------------------
// \!fn void Z4c::ADMGaugeWave2(AthenaArray<Real> & u, bool shifted)
// \brief Initialize ADM vars for 2D gauge wave test

// BD:
// Note below SQRT2 factors added to have normalization coincide with testbed
// papers.

void Z4c::ADMGaugeWave2(AthenaArray<Real> & u_adm) {
  ADM_vars adm;
  SetADMAliases(u_adm, adm);

  Real const amp = opt.AwA_amplitude;
  Real const d_x = opt.AwA_d_x;
  Real const d_y = opt.AwA_d_y;

  // Flat spacetime
  ADMMinkowski(u_adm);

  GLOOP2(k,j) {
    // g_xx, g_yy
    for(int a = 0; a < NDIM-1; ++a) {
      GLOOP1(i) {
        adm.g_dd(a,a,k,j,i) -= 0.5*SINWAVE(amp, d_y, d_y,
                                           mbi.x1(i), mbi.x2(j));
      }
    }

    // g_xy
    GLOOP1(i) {
      adm.g_dd(0,1,k,j,i) = 0.5*SINWAVE(amp, d_y, d_y, mbi.x1(i), mbi.x2(j));
    }

    // K_xx, K_yy
    for(int a = 0; a < NDIM-1; ++a) {
      GLOOP1(i) {
        adm.K_dd(a,a,k,j,i) =
                  0.25 * SQRT2 * DSINWAVE(amp, d_y, d_y, mbi.x1(i), mbi.x2(j))/
          std::sqrt(1.0-SINWAVE(amp, d_y, d_y, mbi.x1(i), mbi.x2(j)));
      }
    }

    // K_xy, K_zz
    GLOOP1(i) {
      adm.K_dd(0,1,k,j,i) = - adm.K_dd(0,0,k,j,i);
      adm.K_dd(2,2,k,j,i) = -0.5 * SQRT2 * DSINWAVE(amp, d_y, d_y,
        mbi.x1(i), mbi.x2(j));
    }
  }
}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::GaugeGaugeWave2(AthenaArray<Real> & u)
// \brief Initialize lapse for 2D gauge wave test

void Z4c::GaugeGaugeWave2(AthenaArray<Real> & u) {
  Z4c_vars z4c;
  SetZ4cAliases(u, z4c);
  z4c.alpha.Fill(1.);
  z4c.beta_u.Fill(0.);

  Real const amp = opt.AwA_amplitude;
  Real const d_x = opt.AwA_d_x;
  Real const d_y = opt.AwA_d_y;

  GLOOP2(k,j) {
     GLOOP1(i) {
        // lapse
        z4c.alpha(k,j,i) = std::sqrt(1.0 - SINWAVE(amp, d_x, d_y,
                                                   mbi.x1(i), mbi.x2(j)));
      }
  }
}

#ifdef GSL
//----------------------------------------------------------------------------------------
// \!fn void Z4c::GaugePolarisedGowdy(AthenaArray<Real> & u)
// \brief Seed lapse for polarised Gowdy test.

void Z4c::GaugePolarisedGowdy(AthenaArray<Real> & u) {
  Z4c_vars z4c;
  SetZ4cAliases(u, z4c);
  z4c.alpha.Fill(1.);
  z4c.beta_u.Fill(0.);

  // compute info for \Lambda--------------------------------------------------
  Real t = opt.AwA_polarised_Gowdy_t0;

  Real J0 = gsl_sf_bessel_J0(2. * PI);
  Real J1 = gsl_sf_bessel_J1(2. * PI);
  Real sqr_J0 = SQR(J0);
  Real sqr_J1 = SQR(J1);

  Real J0_t = gsl_sf_bessel_J0(2. * PI * t);
  Real J1_t = gsl_sf_bessel_J1(2. * PI * t);
  Real sqr_J0_t = SQR(J0_t);
  Real sqr_J1_t = SQR(J1_t);

  Real dt_J0_t = -2. * PI * J1_t;
  Real dt_J1_t = 2. * PI * J0_t - J1_t / t;

  Real sqr_pi = SQR(PI);
  Real sqr_t = SQR(t);

  // Real L = -2. * PI * t * J0_t * J1_t * sqr_cos_x;
  // L += 2. * sqr_pi * sqr_t * (sqr_J0_t + sqr_J1_t);
  // L -= 1. / 2. * (4. * sqr_pi * (sqr_J0 + sqr_J1) - 2. * PI * J0 * J1);

  Real pow_t_m_1_4 = std::pow(t, -1./4.);
  //-----------


  GLOOP3(k,j,i) {
    Real cos_x = std::cos(2. * PI * mbi.x1(i));
    Real sqr_cos_x = SQR(cos_x);

    Real L = -2. * PI * t * J0_t * J1_t * sqr_cos_x;
    L += 2. * sqr_pi * sqr_t * (sqr_J0_t + sqr_J1_t);
    L -= 1. / 2. * (4. * sqr_pi * (sqr_J0 + sqr_J1) - 2. * PI * J0 * J1);

    z4c.alpha(k,j,i) = pow_t_m_1_4 * std::exp(L / 4.);
  }

}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::ADMPolarisedGowdy(AthenaArray<Real> & u)
// \brief Seed ADM variables for polarised Gowdy test.

void Z4c::ADMPolarisedGowdy(AthenaArray<Real> & u_adm) {
  ADM_vars adm;
  SetADMAliases(u_adm, adm);

  // Flat spacetime
  ADMMinkowski(u_adm);

  // compute P, \Lambda factors and derivatives--------------------------------
  Real sign_K = 1.;   // +1 expanding, -1 collapsing

  Real t = opt.AwA_polarised_Gowdy_t0;

  Real J0 = gsl_sf_bessel_J0(2. * PI);
  Real J1 = gsl_sf_bessel_J1(2. * PI);
  Real sqr_J0 = SQR(J0);
  Real sqr_J1 = SQR(J1);

  Real J0_t = gsl_sf_bessel_J0(2. * PI * t);
  Real J1_t = gsl_sf_bessel_J1(2. * PI * t);
  Real sqr_J0_t = SQR(J0_t);
  Real sqr_J1_t = SQR(J1_t);

  // Real dt_J0_t = -2. * PI * J1_t;
  // Real dt_J1_t = 2. * PI * J0_t - J1_t / t;

  Real sqr_pi = SQR(PI);
  Real sqr_t = SQR(t);

  Real pow_t_m_1_4 = std::pow(t, -1./4.);
  Real pow_t_p_1_4 = std::pow(t, 1./4.);
  Real pow_t_m_1_2 = std::pow(t, -1./2.);
  //---------------------------------------------------------------------------

  GLOOP3(k,j,i) {
    Real cos_x = std::cos(2. * PI * mbi.x1(i));
    Real sqr_cos_x = SQR(cos_x);

    Real sin_x = std::sin(2. * PI * mbi.x1(i));
    Real sqr_sin_x = SQR(sin_x);


    Real P = J0_t * cos_x;
    Real L = -2. * PI * t * J0_t * J1_t * sqr_cos_x;
    L += 2. * sqr_pi * sqr_t * (sqr_J0_t + sqr_J1_t);
    L -= 1. / 2. * (4. * sqr_pi * (sqr_J0 + sqr_J1) - 2. * PI * J0 * J1 );

    // Real dt_P = dt_J0_t * cos_x;
    // Real dt_L = - 2. * PI * sqr_cos_x * (t * J1_t * dt_J0_t +
    //                                     J0_t * (J1_t + t * dt_J1_t));
    // dt_L += 4. * sqr_pi * t * (sqr_J0_t + t * J0_t * dt_J0_t +
    //                           J1_t * (J1_t + t * dt_J1_t));

    Real dt_P = - 2 * PI * J1_t * cos_x;
    Real dt_L = 4 * sqr_pi * t * (sqr_J1_t * sqr_cos_x + sqr_J0_t * sqr_sin_x);

    // g_xx
    adm.g_dd(0,0,k,j,i) = pow_t_m_1_2 * std::exp(L / 2.);
    // g_yy
    adm.g_dd(1,1,k,j,i) = t * std::exp(P);
    // g_zz
    adm.g_dd(2,2,k,j,i) = t * std::exp(-P);

    // K_xx
    adm.K_dd(0,0,k,j,i) = sign_K / 4. * pow_t_m_1_4 * std::exp(L / 4.) *
      (1. / t - dt_L);

    // K_yy
    adm.K_dd(1,1,k,j,i) = sign_K / 2. * pow_t_p_1_4 * std::exp(-L / 4.) *
      std::exp(P) * (-1. - t * dt_P);

    // K_zz
    adm.K_dd(2,2,k,j,i) = sign_K / 2. * pow_t_p_1_4 * std::exp(-L / 4.) *
      std::exp(-P) * (-1. + t * dt_P);

    //psi4 [though this is inferred]
    adm.psi4(k,j,i) = std::pow(adm.g_dd(0,0,k,j,i) * adm.g_dd(1,1,k,j,i) *
      adm.g_dd(2,2,k,j,i), -1./3.);
  }

}

#endif //GSL