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
std::default_random_engine generator;
std::uniform_real_distribution<double> distribution(-1.,1.);
#define RANDOMNUMBER (distribution(generator))

// Sin wave for various wave tests
#define  SINWAVE(a,dx,dy,x,y) ( (a)*        std::sin(2.*M_PI*((x)/(dx) - (y)/(dy))) )
#define DSINWAVE(a,dx,dy,x,y) (-(a)*2.*M_PI*std::cos(2.*M_PI*((x)/(dx) - (y)/(dy))) )

// Gaussian profile for simple 3D gauge wave test
#define  GAUSSIAN(a,d,x,y,z) (           (a)*std::exp(-(SQR(x)+SQR(y)+SQR(z))/SQR(d)) )

// Athena++ headers
#include "z4c.hpp"
#include "z4c_macro.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"

//----------------------------------------------------------------------------------------
// \!fn void Z4c::ADMRobustStability(AthenaArray<Real> & u)
// \brief Initialize ADM vars for robust stability test

// Note the amplitude of the noise (~1e-10) should be also rescaled by
// the square of the rho parameter

void Z4c::ADMRobustStability(AthenaArray<Real> & u_adm)
{
  ADM_vars adm;
  SetADMAliases(u_adm, adm);

  // Flat spacetime
  ADMMinkowski(u_adm);

  GLOOP2(k,j) {
    // g_ab
    for(int a = 0; a < NDIM; ++a)
      for(int b = a; b < NDIM; ++b) {
        GLOOP1(i) {
          adm.g_dd(a,b,k,j,i) += RANDOMNUMBER*opt.AwA_amplitude;
        }
      }
    // K_ab
    for(int a = 0; a < NDIM; ++a)
      for(int b = a; b < NDIM; ++b) {
        GLOOP1(i) {
          adm.K_dd(a,b,k,j,i) += RANDOMNUMBER*opt.AwA_amplitude;
        }
      }
  }
}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::GaugeRobStab(AthenaArray<Real> & u, bool shifted)
// \brief Initialize lapse and shift for Robust Stability test

void Z4c::GaugeRobStab(AthenaArray<Real> & u)
{
  Z4c_vars z4c;
  SetZ4cAliases(u, z4c);
  z4c.alpha.Fill(1.);
  z4c.beta_u.Fill(0.);

  MeshBlock * pmb = pmy_block;
  Coordinates * pco = pmb->pcoord;

  GLOOP2(k,j) {
    GLOOP1(i) {
      // lapse
      z4c.alpha(k,j,i) += RANDOMNUMBER*opt.AwA_amplitude;
      // shift
      z4c.beta_u(0,k,j,i) += RANDOMNUMBER*opt.AwA_amplitude;
      z4c.beta_u(1,k,j,i) += RANDOMNUMBER*opt.AwA_amplitude;
      z4c.beta_u(2,k,j,i) += RANDOMNUMBER*opt.AwA_amplitude;
    }
  }
}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::LinearWave1(AthenaArray<Real> & u)
// \brief Initialize ADM vars for linear wave test in 1d

void Z4c::ADMLinearWave1(AthenaArray<Real> & u_adm)
{
  ADM_vars adm;
  SetADMAliases(u_adm, adm);

  MeshBlock * pmb = pmy_block;
  Coordinates * pco = pmb->pcoord;
  
  // Flat spacetime
  ADMMinkowski(u_adm);
  
  // For propagation along x ...
  GLOOP2(k,j) {    
    // g_yy
    GLOOP1(i) {
      adm.g_dd(1,1,k,j,i) += SINWAVE(opt.AwA_amplitude,opt.AwA_d_x,opt.AwA_d_y,pco->x1v(i),0.);
    }
    // g_zz
    GLOOP1(i) {
      adm.g_dd(2,2,k,j,i) -= SINWAVE(opt.AwA_amplitude,opt.AwA_d_x,opt.AwA_d_y,pco->x1v(i),0.);
    }
    // K_yy
    GLOOP1(i) {
      adm.K_dd(1,1,k,j,i) += 0.5*DSINWAVE(opt.AwA_amplitude,opt.AwA_d_x,opt.AwA_d_y,pco->x1v(i),0.);
    }
    // K_zz
    GLOOP1(i) {
      adm.K_dd(2,2,k,j,i) -= 0.5*DSINWAVE(opt.AwA_amplitude,opt.AwA_d_x,opt.AwA_d_y,pco->x1v(i),0.);
    }
  }
}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::LinearWave2(AthenaArray<Real> & u)
// \brief Initialize ADM vars for linear wave test in 2d

// Note we use the same macro as for the 1d,
// so there's a factor of sqrt(2) in the normalization different from literature

void Z4c::ADMLinearWave2(AthenaArray<Real> & u_adm)
{
  ADM_vars adm;
  SetADMAliases(u_adm, adm);

  MeshBlock * pmb = pmy_block;
  Coordinates * pco = pmb->pcoord;

  // Flat spacetime
  ADMMinkowski(u_adm);

  GLOOP2(k,j) {
    // g_xx, g_yy
    for(int a = 0; a < NDIM-1; ++a) {
      GLOOP1(i) {
        adm.g_dd(a,a,k,j,i) += 0.5*SINWAVE(opt.AwA_amplitude,opt.AwA_d_x,opt.AwA_d_y,pco->x1v(i),pco->x2v(j));
      }
    }

    // g_xy, g_zz
    GLOOP1(i) {
      adm.g_dd(0,1,k,j,i)  = 0.5*SINWAVE(opt.AwA_amplitude,opt.AwA_d_x,opt.AwA_d_y,pco->x1v(i),pco->x2v(j));
      adm.g_dd(2,2,k,j,i) -=     SINWAVE(opt.AwA_amplitude,opt.AwA_d_x,opt.AwA_d_y,pco->x1v(i),pco->x2v(j));
    }

    // K_xx, K_xy, K_yy
    for(int a = 0; a < NDIM-1; ++a) {
      for(int b = a; b < NDIM-1; ++b) {
        GLOOP1(i) {
          adm.K_dd(a,b,k,j,i) = 0.25*DSINWAVE(opt.AwA_amplitude,opt.AwA_d_x,opt.AwA_d_y,pco->x1v(i),pco->x2v(j));
        }
      }
    }
    // K_zz
    GLOOP1(i) {
      adm.K_dd(2,2,k,j,i) = -0.5*DSINWAVE(opt.AwA_amplitude,opt.AwA_d_x,opt.AwA_d_y,pco->x1v(i),pco->x2v(j));
    }
  }
}


//----------------------------------------------------------------------------------------
// \!fn void Z4c::GaugeSimpleGaugeWave(AthenaArray<Real> & u)
// \brief Initializes lapse for simple 3D gauge wave test

void Z4c::GaugeSimpleGaugeWave(AthenaArray<Real> & u)
{
  Z4c_vars z4c;
  SetZ4cAliases(u, z4c);
  z4c.alpha.Fill(1.);
  z4c.beta_u.Fill(0.);

  MeshBlock * pmb = pmy_block;
  Coordinates * pco = pmb->pcoord;

  GLOOP2(k,j) {
      GLOOP1(i) {
        // lapse
        //z4c.alpha(k,j,i) += GAUSSIAN(opt.AwA_amplitude,opt.AwA_d_x,pco->x1v(i),pco->x2v(j),pco->x3v(k));
        z4c.alpha(k,j,i) += opt.AwA_amplitude*pow(sin(2.*M_PI*pco->x1v(i)),6)*pow(sin(2.*M_PI*pco->x2v(j)),6)*pow(sin(2.*M_PI*pco->x3v(k)),6);
      }
  }
}


//----------------------------------------------------------------------------------------
// \!fn void Z4c::ADMGaugeWave1(AthenaArray<Real> & u)
// \brief Initialize ADM vars for 1D gauge wave test

void Z4c::ADMGaugeWave1(AthenaArray<Real> & u_adm)
{
  ADM_vars adm;
  SetADMAliases(u_adm, adm);

  MeshBlock * pmb = pmy_block;
  Coordinates * pco = pmb->pcoord;

  // Flat spacetime
  ADMMinkowski(u_adm);

  // Propagation along x
  GLOOP2(k,j) {
     GLOOP1(i) {
        // g_xx
        adm.g_dd(0,0,k,j,i) -= SINWAVE(opt.AwA_amplitude,opt.AwA_d_x,opt.AwA_d_y,pco->x1v(i),0.);
        // K_xx
        adm.K_dd(0,0,k,j,i) =
                   0.5*DSINWAVE(opt.AwA_amplitude,opt.AwA_d_x,opt.AwA_d_y,pco->x1v(i),0.)/
        std::sqrt(1.0 - SINWAVE(opt.AwA_amplitude,opt.AwA_d_x,opt.AwA_d_y,pco->x1v(i),0.));
      }
  }
}


//----------------------------------------------------------------------------------------
// \!fn void Z4c::GaugeGaugeWave1(AthenaArray<Real> & u)
// \brief Initialize lapse for 1D gauge wave test

void Z4c::GaugeGaugeWave1(AthenaArray<Real> & u)
{
  Z4c_vars z4c;
  SetZ4cAliases(u, z4c);
  z4c.alpha.Fill(1.);
  z4c.beta_u.Fill(0.);

  MeshBlock * pmb = pmy_block;
  Coordinates * pco = pmb->pcoord;

  GLOOP2(k,j) {
     GLOOP1(i) {
       // lapse
       z4c.alpha(k,j,i) = std::sqrt(1.0 - SINWAVE(opt.AwA_amplitude,opt.AwA_d_x,opt.AwA_d_y,pco->x1v(i),0.));
     }
  }
}


//----------------------------------------------------------------------------------------
// \!fn void Z4c::ADMGaugeWave1_shifted(AthenaArray<Real> & u)
// \brief Initialize ADM vars for shifted 1D gauge wave test

void Z4c::ADMGaugeWave1_shifted(AthenaArray<Real> & u_adm)
{
  ADM_vars adm;
  SetADMAliases(u_adm, adm);

  MeshBlock * pmb = pmy_block;
  Coordinates * pco = pmb->pcoord;

  // Flat spacetime
  ADMMinkowski(u_adm);

  // Propagation along x ...
  GLOOP2(k,j) {
     GLOOP1(i) {
        // g_xx
        adm.g_dd(0,0,k,j,i) += SINWAVE(opt.AwA_amplitude,opt.AwA_d_x,opt.AwA_d_y,pco->x1v(i),0.);
        // K_xx
        adm.K_dd(0,0,k,j,i) =
                   0.5*DSINWAVE(opt.AwA_amplitude,opt.AwA_d_x,opt.AwA_d_y,pco->x1v(i),0.)/
        std::sqrt(1.0 + SINWAVE(opt.AwA_amplitude,opt.AwA_d_x,opt.AwA_d_y,pco->x1v(i),0.));
    }
  }
}


//----------------------------------------------------------------------------------------
// \!fn void Z4c::GaugeGaugeWave1_shifted(AthenaArray<Real> & u)
// \brief Initialize lapse and shift for shifted 1D gauge wave test

void Z4c::GaugeGaugeWave1_shifted(AthenaArray<Real> & u)
{
  Z4c_vars z4c;
  SetZ4cAliases(u, z4c);
  z4c.alpha.Fill(1.);
  z4c.beta_u.Fill(0.);

  MeshBlock * pmb = pmy_block;
  Coordinates * pco = pmb->pcoord;

  GLOOP2(k,j) {
     GLOOP1(i) {
        // lapse
        z4c.alpha(k,j,i) = 1.0/(std::sqrt(1.0 + SINWAVE(opt.AwA_amplitude,opt.AwA_d_x,opt.AwA_d_y,pco->x1v(i),0.)));
        // shift
        z4c.beta_u(0,k,j,i) = - SINWAVE(opt.AwA_amplitude,opt.AwA_d_x,opt.AwA_d_y,pco->x1v(i),0.)/
                          (1. + SINWAVE(opt.AwA_amplitude,opt.AwA_d_x,opt.AwA_d_y,pco->x1v(i),0.));
      }
  }
}


//----------------------------------------------------------------------------------------
// \!fn void Z4c::ADMGaugeWave2(AthenaArray<Real> & u, bool shifted)
// \brief Initialize ADM vars for 2D gauge wave test

// Note we use the same macro as for the 1d,
// so there's a factor of sqrt(2) in the normalization different from literature

void Z4c::ADMGaugeWave2(AthenaArray<Real> & u_adm)
{
  ADM_vars adm;
  SetADMAliases(u_adm, adm);

  MeshBlock * pmb = pmy_block;
  Coordinates * pco = pmb->pcoord;

  // Flat spacetime
  ADMMinkowski(u_adm);

  GLOOP2(k,j) {
    // g_xx, g_yy
    for(int a = 0; a < NDIM-1; ++a) {
      GLOOP1(i) {
        adm.g_dd(a,a,k,j,i) -= 0.5*SINWAVE(opt.AwA_amplitude,opt.AwA_d_x,opt.AwA_d_y,pco->x1v(i),pco->x2v(j));
      }
    }

    // g_xy
    GLOOP1(i) {
      adm.g_dd(0,1,k,j,i) = 0.5*SINWAVE(opt.AwA_amplitude,opt.AwA_d_x,opt.AwA_d_y,pco->x1v(i),pco->x2v(j));
    }

    // K_xx, K_yy
    for(int a = 0; a < NDIM-1; ++a) {
      GLOOP1(i) {
        adm.K_dd(a,a,k,j,i) =
                  0.25*DSINWAVE(opt.AwA_amplitude,opt.AwA_d_x,opt.AwA_d_y,pco->x1v(i),pco->x2v(j))/
          std::sqrt(1.0-SINWAVE(opt.AwA_amplitude,opt.AwA_d_x,opt.AwA_d_y,pco->x1v(i),pco->x2v(j)));
      }
    }

    // K_xy, K_zz
    GLOOP1(i) {
      adm.K_dd(0,1,k,j,i) = - adm.K_dd(0,0,k,j,i);
      adm.K_dd(2,2,k,j,i) = -0.5*DSINWAVE(opt.AwA_amplitude,opt.AwA_d_x,opt.AwA_d_y,pco->x1v(i),pco->x2v(j));
    }    
  }
}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::GaugeGaugeWave2(AthenaArray<Real> & u)
// \brief Initialize lapse for 2D gauge wave test

void Z4c::GaugeGaugeWave2(AthenaArray<Real> & u)
{
  Z4c_vars z4c;
  SetZ4cAliases(u, z4c);
  z4c.alpha.Fill(1.);
  z4c.beta_u.Fill(0.);

  MeshBlock * pmb = pmy_block;
  Coordinates * pco = pmb->pcoord;

  GLOOP2(k,j) {
     GLOOP1(i) {
        // lapse
        z4c.alpha(k,j,i) = std::sqrt(1.0 - SINWAVE(opt.AwA_amplitude,opt.AwA_d_x,opt.AwA_d_y,pco->x1v(i),pco->x2v(j)));
      }
  }
}

