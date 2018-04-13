//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file kh.cpp
//  \brief Problem generator for KH instability.
//
// Sets up two different problems:
//   - iprob=1: slip surface with random perturbations
//   - iprob=2: tanh profile at interface, with single-mode perturbation
//========================================================================================

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../utils/utils.hpp"

#include <algorithm>  // min, max
#include <cmath>


Real vflow;
int RefinementCondition(MeshBlock *pmb);

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  if(adaptive==true)
    EnrollUserRefinementCondition(RefinementCondition);
  vflow = pin->GetReal("problem","vflow");

  return;
}


//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the Kelvin-Helmholz test
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  int64_t iseed = -1 - gid;
  Real gm1 = peos->GetGamma() - 1.0;

  // Read problem parameters
  int iprob = pin->GetInteger("problem","iprob");
  Real drat = pin->GetReal("problem","drat");
  Real amp = pin->GetReal("problem","amp");

//--- iprob=1.  Two uniform streams moving at +/- vflow, random perturbations

  if (iprob == 1) {
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      phydro->u(IDN,k,j,i) = 1.0;
      phydro->u(IM1,k,j,i) = vflow + amp*(ran2(&iseed) - 0.5);
      phydro->u(IM2,k,j,i) = amp*(ran2(&iseed) - 0.5);
      phydro->u(IM3,k,j,i) = 0.0;
      if (fabs(pcoord->x2v(j)) < 0.25) {
        phydro->u(IDN,k,j,i) = drat;
        phydro->u(IM1,k,j,i) = -drat*(vflow + amp*(ran2(&iseed) - 0.5));
        phydro->u(IM2,k,j,i) = drat*amp*(ran2(&iseed) - 0.5);
      }
      // Pressure scaled to give a sound speed of 1 with gamma=1.4
      if (NON_BAROTROPIC_EOS) {
        phydro->u(IEN,k,j,i) = 2.5/gm1 + 0.5*(SQR(phydro->u(IM1,k,j,i)) +
          SQR(phydro->u(IM2,k,j,i)))/phydro->u(IDN,k,j,i);
      }
    }}}
  }

//--- iprob=2. Two uniform density flows with single mode pert., based on Ryu&Jones.

  if (iprob == 2) {
    Real a = 0.02;
    Real sigma = 0.2;
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      phydro->u(IDN,k,j,i) = 1.0;
      phydro->u(IM1,k,j,i) = vflow*tanh((pcoord->x2v(j))/a);
      phydro->u(IM2,k,j,i) = amp*sin(2.0*PI*pcoord->x1v(i))
        *exp(-(SQR(pcoord->x2v(j)))/SQR(sigma));
      phydro->u(IM3,k,j,i) = 0.0;
      if (NON_BAROTROPIC_EOS) {
        phydro->u(IEN,k,j,i) = 1.0/gm1 + 0.5*(SQR(phydro->u(IM1,k,j,i)) +
          SQR(phydro->u(IM2,k,j,i)))/phydro->u(IDN,k,j,i);
      }
    }}}
  }

//--- iprob=3.  Test in SR paper, based on iprob=2

  if (iprob == 3) {
    Real a = 0.01;
    Real sigma = 0.1;
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      phydro->u(IDN,k,j,i) = 0.505 + 0.495*tanh((fabs(pcoord->x2v(j))-0.5)/a);
      phydro->u(IM1,k,j,i) = vflow*tanh((fabs(pcoord->x2v(j))-0.5)/a);
      phydro->u(IM2,k,j,i) = amp*vflow*sin(2.0*PI*pcoord->x1v(i))
               *exp(-((fabs(pcoord->x2v(j))-0.5)*(fabs(pcoord->x2v(j))-0.5))/(sigma*sigma));
      if (pcoord->x2v(j) < 0.0) phydro->u(IM2,k,j,i) *= -1.0;
      phydro->u(IM1,k,j,i) *= phydro->u(IDN,k,j,i);
      phydro->u(IM2,k,j,i) *= phydro->u(IDN,k,j,i);
      phydro->u(IM3,k,j,i) = 0.0;
      if (NON_BAROTROPIC_EOS) {
        phydro->u(IEN,k,j,i) = 1.0/gm1 + 0.5*(SQR(phydro->u(IM1,k,j,i)) +
          SQR(phydro->u(IM2,k,j,i)))/phydro->u(IDN,k,j,i);
      }
    }}}
  }

//--- iprob=4.  Two uniform streams moving at +/- vflow, single mode
  if (iprob == 4) {
    Real a = 0.02;
    Real sigma = 0.2;
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      phydro->u(IDN,k,j,i) = (1.5-0.5*tanh((pcoord->x2v(j))/a));
      phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i)*vflow*tanh((pcoord->x2v(j))/a);
      phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i)*amp*cos(2.0*PI*pcoord->x1v(i))
                            *exp(-(SQR(pcoord->x2v(j)))/SQR(sigma));
      phydro->u(IM3,k,j,i) = 0.0;
      if (NON_BAROTROPIC_EOS) {
        phydro->u(IEN,k,j,i) = 2.5/gm1 + 0.5*(SQR(phydro->u(IM1,k,j,i)) +
          SQR(phydro->u(IM2,k,j,i)))/phydro->u(IDN,k,j,i);
      }
    }}}
  }

  // initialize interface B, same for all iprob
  if (MAGNETIC_FIELDS_ENABLED) {
    Real b0 = pin->GetReal("problem","b0");
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie+1; i++) {
      pfield->b.x1f(k,j,i) = b0;
    }}}
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je+1; j++) {
    for (int i=is; i<=ie; i++) {
      pfield->b.x2f(k,j,i) = 0.0;
    }}}
    for (int k=ks; k<=ke+1; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      pfield->b.x3f(k,j,i) = 0.0;
    }}}
    if (NON_BAROTROPIC_EOS) {
      for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        phydro->u(IEN,k,j,i) += 0.5*b0*b0;
      }}}
    }
  }

  return;
}


// refinement condition: velocity gradient
int RefinementCondition(MeshBlock *pmb)
{
  AthenaArray<Real> &w = pmb->phydro->w;
  Real vgmax=0.0;
  for(int k=pmb->ks; k<=pmb->ke; k++) {
    for(int j=pmb->js; j<=pmb->je; j++) {
      for(int i=pmb->is; i<=pmb->ie; i++) {
        Real vgy=std::fabs(w(IVY,k,j,i+1)-w(IVY,k,j,i-1))*0.5;
        Real vgx=std::fabs(w(IVX,k,j+1,i)-w(IVX,k,j-1,i))*0.5;
        if(vgy > vgmax) vgmax=vgy;
        if(vgx > vgmax) vgmax=vgx;
      }
    }
  }
  if(vgmax > 0.01) return 1;
  return -1;
}
