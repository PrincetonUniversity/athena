//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file ran2.cpp 

// C++ headers
#include "float.h"  // DBL_EPSILON

// Athena headers
#include "../athena.hpp"

//----------------------------------------------------------------------------------------
//! \fn double ran2(long int *idum)
//  \brief  Extracted from the Numerical Recipes in C (version 2) code. Modified
//   to use doubles instead of floats. -- T. A. Gardiner -- Aug. 12, 2003
//
// Long period (> 2 x 10^{18}) random number generator of L'Ecuyer with Bays-Durham
// shuffle and added safeguards.  Returns a uniform random deviate between 0.0 and 1.0
// (exclusive of the endpoint values).  Call with idum = a negative integer to
// initialize; thereafter, do not alter idum between successive deviates in a sequence.
// RNMX should appriximate the largest floating point value that is less than 1. 

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
