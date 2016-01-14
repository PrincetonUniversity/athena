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
//! \file kh.cpp
//  \brief Problem generator for KH instability. 
//
// Sets up two different problems:
//   - iprob=1: slip surface with random perturbations
//   - iprob=2: tanh profile at interface, with single-mode perturbation
//======================================================================================

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../hydro/eos/eos.hpp"
#include "../coordinates/coordinates.hpp"
#include "../utils/utils.hpp"

//======================================================================================
//! \fn void Mesh::InitUserMeshProperties(ParameterInput *pin)
//  \brief Init the Mesh properties
//======================================================================================

void Mesh::InitUserMeshProperties(ParameterInput *pin)
{
  return;
}


//======================================================================================
//! \fn void Mesh::TerminateUserMeshProperties(void)
//  \brief Clean up the Mesh properties
//======================================================================================

void Mesh::TerminateUserMeshProperties(void)
{
  // nothing to do
  return;
}


//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the Kelvin-Helmholz test
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  long int iseed = -1;
  Real gm1 = phydro->peos->GetGamma() - 1.0;

  // Read problem parameters
  int iprob = pin->GetInteger("problem","iprob");
  Real vflow = pin->GetReal("problem","vflow");
  Real drat = pin->GetReal("problem","drat");
  Real amp = pin->GetReal("problem","amp");

// iprob=1.  Two uniform streams moving at +/- vflow, random perturbations

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

// iprob=2. Two uniform density flows with single mode perturbation, based on Ryu&Jones.

  if (iprob == 2) {
    Real a = 0.05;
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

// iprob=3.  Test in SR paper, based on iprob=2

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


//======================================================================================
//! \fn void MeshBlock::UserWorkInLoop(void)
//  \brief User-defined work function for every time step
//======================================================================================

void MeshBlock::UserWorkInLoop(void)
{
  // nothing to do
  return;
}

