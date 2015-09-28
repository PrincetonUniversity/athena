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

void Mesh::ProblemGenerator(Hydro *pfl, Field *pfd, ParameterInput *pin)
{
  MeshBlock *pmb = pfl->pmy_block;
  Coordinates *pco = pmb->pcoord;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  long int iseed = -1;
  Real gm1 = pfl->pf_eos->GetGamma() - 1.0;

  // Read problem parameters
  int iprob = pin->GetInteger("problem","iprob");
  Real vflow = pin->GetReal("problem","vflow");
  Real drat = pin->GetReal("problem","drat");
  Real amp = pin->GetReal("problem","amp");
  Real b0  = pin->GetReal("problem","b0");

// iprob=1.  Two uniform streams moving at +/- vflow, random perturbations

  if (iprob == 1) {
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      pfl->u(IDN,k,j,i) = 1.0;
      pfl->u(IM1,k,j,i) = vflow + amp*(ran2(&iseed) - 0.5);
      pfl->u(IM2,k,j,i) = amp*(ran2(&iseed) - 0.5);
      pfl->u(IM3,k,j,i) = 0.0;
      if (fabs(pco->x2v(j)) < 0.25) {
        pfl->u(IDN,k,j,i) = drat;
        pfl->u(IM1,k,j,i) = -drat*(vflow + amp*(ran2(&iseed) - 0.5));
        pfl->u(IM2,k,j,i) = drat*amp*(ran2(&iseed) - 0.5);
      }
      // Pressure scaled to give a sound speed of 1 with gamma=1.4 
      if (NON_BAROTROPIC_EOS) {
        pfl->u(IEN,k,j,i) = 2.5/gm1 + 0.5*(SQR(pfl->u(IM1,k,j,i)) +
          SQR(pfl->u(IM2,k,j,i)))/pfl->u(IDN,k,j,i);
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
      pfl->u(IDN,k,j,i) = 1.0;
      pfl->u(IM1,k,j,i) = vflow*tanh((pco->x2v(j))/a);
      pfl->u(IM2,k,j,i) = amp*sin(2.0*PI*pco->x1v(i))
        *exp(-(SQR(pco->x2v(j)))/SQR(sigma));
      pfl->u(IM3,k,j,i) = 0.0;
      if (NON_BAROTROPIC_EOS) {
        pfl->u(IEN,k,j,i) = 1.0/gm1 + 0.5*(SQR(pfl->u(IM1,k,j,i)) +
          SQR(pfl->u(IM2,k,j,i)))/pfl->u(IDN,k,j,i);
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
      pfl->u(IDN,k,j,i) = 0.505 + 0.495*tanh((fabs(pco->x2v(j))-0.5)/a);
      pfl->u(IM1,k,j,i) = vflow*tanh((fabs(pco->x2v(j))-0.5)/a);
      pfl->u(IM2,k,j,i) = amp*vflow*sin(2.0*PI*pco->x1v(i))
               *exp(-((fabs(pco->x2v(j))-0.5)*(fabs(pco->x2v(j))-0.5))/(sigma*sigma));
      if (pco->x2v(j) < 0.0) pfl->u(IM2,k,j,i) *= -1.0;
      pfl->u(IM1,k,j,i) *= pfl->u(IDN,k,j,i);
      pfl->u(IM2,k,j,i) *= pfl->u(IDN,k,j,i);
      pfl->u(IM3,k,j,i) = 0.0;
      if (NON_BAROTROPIC_EOS) {
        pfl->u(IEN,k,j,i) = 1.0/gm1 + 0.5*(SQR(pfl->u(IM1,k,j,i)) +
          SQR(pfl->u(IM2,k,j,i)))/pfl->u(IDN,k,j,i);
      }
    }}}
  }

  // initialize interface B, same for all iprob
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie+1; i++) {
      pfd->b.x1f(k,j,i) = b0;
    }}}
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je+1; j++) {
    for (int i=is; i<=ie; i++) {
      pfd->b.x2f(k,j,i) = 0.0;
    }}}
    for (int k=ks; k<=ke+1; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      pfd->b.x3f(k,j,i) = 0.0;
    }}}
    if (NON_BAROTROPIC_EOS) {
      for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        pfl->u(IEN,k,j,i) += 0.5*b0*b0;
      }}}
    }
  }

  return;
}
