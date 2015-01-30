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

// Primary header
#include "../mesh.hpp"

// C++ headers
#include "float.h"    // DBL_EPSILON

// Athena headers
#include "../athena.hpp"           // enums, Real
#include "../athena_arrays.hpp"    // AthenaArray
#include "../parameter_input.hpp"  // ParameterInput
#include "../fluid/fluid.hpp"      // Fluid
#include "../fluid/eos/eos.hpp"    // GetGamma
#include "../field/field.hpp"      // magnetic field

static double ran2(long int *idum);  // internal random number generator

//======================================================================================
//! \file kh.cpp
//  \brief Problem generator for KH instability. 
//
// Sets up two different problems:
//   - iprob=1: slip surface with random perturbations
//   - iprob=2: tanh profile at interface, with single-mode perturbation
//======================================================================================

void Mesh::ProblemGenerator(Fluid *pfl, Field *pfd, ParameterInput *pin)
{
  MeshBlock *pmb = pfl->pmy_block;
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
      if (fabs(pmb->x2v(j)) < 0.25) {
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
      pfl->u(IM1,k,j,i) = vflow*tanh((pmb->x2v(j))/a);
      pfl->u(IM2,k,j,i) = amp*sin(2.0*PI*pmb->x1v(i))
        *exp(-(SQR(pmb->x2v(j)))/SQR(sigma));
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
      pfl->u(IDN,k,j,i) = 0.505 + 0.495*tanh((fabs(pmb->x2v(j))-0.5)/a);
      pfl->u(IM1,k,j,i) = vflow*tanh((fabs(pmb->x2v(j))-0.5)/a);
      pfl->u(IM2,k,j,i) = amp*vflow*sin(2.0*PI*pmb->x1v(i))
               *exp(-((fabs(pmb->x2v(j))-0.5)*(fabs(pmb->x2v(j))-0.5))/(sigma*sigma));
      if (pmb->x2v(j) < 0.0) pfl->u(IM2,k,j,i) *= -1.0;
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

//--------------------------------------------------------------------------------------
//! \fn double ran2(long int *idum)
//  \brief  Extracted from the Numerical Recipes in C (version 2) code. Modified
//   to use doubles instead of floats. -- T. A. Gardiner -- Aug. 12, 2003
//
// Long period (> 2 x 10^{18}) random number generator of L'Ecuyer with Bays-Durham
// shuffle and added safeguards.  Returns a uniform random deviate between 0.0 and 1.0
// (exclusive of the endpoint values).  Call with idum = a negative integer to
// initialize; thereafter, do not alter idum between successive deviates in a sequence.
//  RNMX should appriximate the largest floating point value that is less than 1. 

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-DBL_EPSILON)

double ran2(long int *idum){
  int j;
  long int k;
  static long int idum2=123456789;
  static long int iy=0;
  static long int iv[NTAB];
  double temp;

  if (*idum <= 0) { // Initialize
    if (-(*idum) < 1) *idum=1; // Be sure to prevent idum = 0
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { // Load the shuffle table (after 8 warm-ups)
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;                 // Start here when not initializing
  *idum=IA1*(*idum-k*IQ1)-k*IR1; // Compute idum=(IA1*idum) % IM1 without
  if (*idum < 0) *idum += IM1;   // overflows by Schrage's method
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; // Compute idum2=(IA2*idum) % IM2 likewise
  if (idum2 < 0) idum2 += IM2;
  j=(int)(iy/NDIV);              // Will be in the range 0...NTAB-1
  iy=iv[j]-idum2;                // Here idum is shuffled, idum and idum2
  iv[j] = *idum;                 // are combined to generate output
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; // No endpoint values
  else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef RNMX
