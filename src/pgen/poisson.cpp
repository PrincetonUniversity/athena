//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file poisson.cpp
//  \brief Problem generator to test Poisson's solver
//

// C++ headers
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <ctime>

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../gravity/gravity.hpp"
#include "../fft/athena_fft.hpp"
#include "../mesh/mesh.hpp"

#ifdef OPENMP_PARALLEL
#include "omp.h"
#endif

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  Real x0=0.0, y0=0.0, z0=0.0;

  if(SELF_GRAVITY_ENABLED){
    pgrav->gconst = 1.0;
    pgrav->four_pi_gconst = 4.0*PI*pgrav->gconst;
    pgrav->grav_mean_rho = 0.0;
    // Assigning Plummer density profile
    int iprob = pin->GetOrAddInteger("problem","iprob",1);
    switch(iprob){
      case 1: { // sine
        if(FFT_ENABLED){
          for (int k=ks; k<=ke; k++) {
          for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            Real x = pcoord->x1v(i);
            Real y = pcoord->x2v(j);
            Real z = pcoord->x3v(k);
            Real dx = pcoord->dx1v(i);
            Real dy = pcoord->dx2v(j);
            Real dz = pcoord->dx3v(k);
            Real den = std::sin(x*pfft->dkx/dx)
                      *std::sin(y*pfft->dky/dy)
                      *std::sin(z*pfft->dkz/dz); // dkx=2*PI/Nx
            phydro->u(IDN,k,j,i) = den;
            phydro->u(IM1,k,j,i) = den;
            phydro->u(IM2,k,j,i) = SQR(den);
            phydro->u(IM3,k,j,i) = den;
          }}} // for-loop
        } // fft
        break;
      } // case 1
      case 2: { // plummer density
        Real M = pin->GetOrAddReal("problem","M",1.0);
        Real a0 = pin->GetOrAddReal("problem","a0",1.0);
        Real da = 3.0*M/(4*PI*a0*a0*a0);
        for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          Real r,r2;
          if (COORDINATE_SYSTEM == "cartesian") {
            Real x = pcoord->x1v(i);
            Real y = pcoord->x2v(j);
            Real z = pcoord->x3v(k);
            r2 = sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
          }
          Real den = da*std::pow((1.0+r2/SQR(a0)),-2.5);
          phydro->u(IDN,k,j,i) = den;
          phydro->u(IM1,k,j,i) = den;
          phydro->u(IM2,k,j,i) = SQR(den);
          phydro->u(IM3,k,j,i) = den;
        }}} //for-loop
      } // case 2
    } // switch
    pgrav->Solver(phydro->u);
  } // self-gravity
}

//========================================================================================
//! \fn void MeshBlock::UserWorkInLoop(void)
//  \brief Function called once every time step for user-defined work.
//========================================================================================

void MeshBlock::UserWorkInLoop(void)
{
  // do nothing
  return;
}
//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{
  Hydro *phydro = pblock->phydro;
  Coordinates *pcoord = pblock->pcoord;
  Gravity *pgrav = pblock->pgrav;
  AthenaFFT *pfft = pblock->pfft;
  Real x0=0.0, y0=0.0, z0=0.0;
  int is=pblock->is, ie=pblock->ie;
  int js=pblock->js, je=pblock->je;
  int ks=pblock->ks, ke=pblock->ke;
  int cnt = (ke-ks+1)*(je-js+1)*(ie-is+1);
  Real phia;

  if(SELF_GRAVITY_ENABLED){
    Real err1=0.0,err2=0.0;
    int iprob = pin->GetOrAddInteger("problem","iprob",1);
    switch(iprob){
      case 1: { // sine
        if(FFT_ENABLED){
          for (int k=ks; k<=ke; k++) {
          for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            Real x = pcoord->x1v(i);
            Real y = pcoord->x2v(j);
            Real z = pcoord->x3v(k);
            Real dx = pcoord->dx1v(i);
            Real dy = pcoord->dx2v(j);
            Real dz = pcoord->dx3v(k);
            Real phi0 =-pgrav->four_pi_gconst
                       /(SQR(pfft->dkx/dx)+SQR(pfft->dky/dy)+SQR(pfft->dkz/dz));
            phia = phi0*phydro->u(IM1,k,j,i);
            err1 += std::abs(pgrav->phi(k,j,i) - phia);
            err2 += std::abs(phydro->u(IM1,k,j,i) - pgrav->phi(k,j,i));
          }}} // for-loop
        } // fft
        break;
      } // case 1
      case 2: { // plummer density
        Real M = pin->GetOrAddReal("problem","M",1.0);
        Real a0 = pin->GetOrAddReal("problem","a0",1.0);
        for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          Real r,r2;
          if (COORDINATE_SYSTEM == "cartesian") {
            Real x = pcoord->x1v(i);
            Real y = pcoord->x2v(j);
            Real z = pcoord->x3v(k);
            r2 = sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
          }
          phia = - pgrav->gconst*M/sqrt(r2+SQR(a0));
          err1 += std::abs(pgrav->phi(k,j,i) - phia);
          err2 += std::abs(phydro->u(IM1,k,j,i) - pgrav->phi(k,j,i));
        }}}
      }
    }

    err1 = err1/cnt;
    err2 = err2/cnt;

    std::cout << "=====================================================" << std::endl;
    std::cout << "L1 : " << err1 <<" L2: " << err2 << std::endl;
    std::cout << "=====================================================" << std::endl;
  }

  return;
}
