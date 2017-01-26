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
  Real x, y, z;
  Real r, r2;
  Real den, phia;
  RegionSize &mesh_size = pmy_mesh->mesh_size;

  Real x1size = mesh_size.x1max - mesh_size.x1min;
  Real x2size = mesh_size.x2max - mesh_size.x2min;
  Real x3size = mesh_size.x3max - mesh_size.x3min;

  Real gconst = 1.0;
  Real four_pi_gconst = 4.0*PI*gconst;
  Real grav_mean_rho = 0.0;
  // Assigning Plummer density profile
  int iprob = pin->GetOrAddInteger("problem","iprob",1);
  int nlim = pin->GetInteger("time","nlim");
  for (int k=ks; k<=ke; ++k) {
    std::cout << k << std::endl;
  for (int j=js; j<=je; ++j) {
  for (int i=is; i<=ie; ++i) {
    x = pcoord->x1v(i);
    y = pcoord->x2v(j);
    z = pcoord->x3v(k);
    r2 = sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
    if(iprob == 1){
      den = std::sin(2*PI*x/x1size)
                *std::sin(2*PI*y/x2size)
                *std::sin(2*PI*z/x3size); // dkx=2*PI/Nx
      phia =-den*four_pi_gconst
                 /(SQR(2*PI/x1size)+SQR(2*PI/x2size)+SQR(2*PI/x3size));
    } else if (iprob == 2){
      Real M = pin->GetOrAddReal("problem","M",1.0);
      Real a0 = pin->GetOrAddReal("problem","a0",1.0);
      Real da = 3.0*M/(4*PI*a0*a0*a0);
      den = da*std::pow((1.0+r2/SQR(a0)),-2.5);
      phia = - gconst*M/sqrt(r2+SQR(a0));
    } else if (iprob == 3){
      Real a0 = pin->GetOrAddReal("problem","a0",1.0);
      den = (4.0*SQR(a0)*r2-6.0*a0)*exp(-a0*r2);
      phia = four_pi_gconst*exp(-a0*r2);
    }

    if(nlim > 0){
      phydro->u(IDN,k,j,i) = den;
      phydro->u(IM1,k,j,i) = 0.0;
      phydro->u(IM2,k,j,i) = 0.0;
      phydro->u(IM3,k,j,i) = 0.0;
    } else {
      phydro->u(IDN,k,j,i) = den;
      phydro->u(IM1,k,j,i) = den;
      phydro->u(IM2,k,j,i) = phia;
      phydro->u(IM3,k,j,i) = 0.0;
    }
  }}}

  if(SELF_GRAVITY_ENABLED){
    pgrav->gconst = gconst;
    pgrav->four_pi_gconst = four_pi_gconst;
    pgrav->grav_mean_rho = grav_mean_rho;
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
  Real x0=0.0, y0=0.0, z0=0.0;
  int is=pblock->is, ie=pblock->ie;
  int js=pblock->js, je=pblock->je;
  int ks=pblock->ks, ke=pblock->ke;
  int cnt = (ke-ks+1)*(je-js+1)*(ie-is+1);

  int nlim = pin->GetInteger("time","nlim");

  if(SELF_GRAVITY_ENABLED && nlim == 0){
    Gravity *pgrav = pblock->pgrav;
    Real err1=0.0,err2=0.0;
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie; ++i) {
      err1 += std::abs(pgrav->phi(k,j,i) - phydro->u(IM2,k,j,i));
      err2 += pgrav->phi(k,j,i)/phydro->u(IM2,k,j,i);
    }}} // for-loop

    err1 = err1/cnt;
    err2 = err2/cnt;

    std::cout << "=====================================================" << std::endl;
    std::cout << "L1 : " << err1 <<" L2: " << err2 << std::endl;
    std::cout << "=====================================================" << std::endl;
  }

  return;
}
