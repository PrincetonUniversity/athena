//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================

// Athena headers
#include "../athena.hpp"


//======================================================================================
//! \file lubksb.cpp
//======================================================================================



//--------------------------------------------------------------------------------------
//! \fn  void  Lubksb_nr(int n, AthenaArray<Real> &a, AthenaArray<int> &indx,
//      AthenaArray<Real> &b)
// \brief Backward substitution
// a is the input matrix done with LU decomposition, n is the matrix size
// indx id the history of row permutation
// b is the vector on the right (AX=b), and is returned with the solution
void Lubksb_nr(int n, AthenaArray<Real> &a, AthenaArray<int> &indx,
               AthenaArray<Real> &b) {
  int ii=-1,ip;
  Real sum;
  // Solve L*y=b
  for (int i=0; i<n; ++i) {
    ip = indx(i);
    sum=b(ip);
    b(ip)=b(i);
    if (ii>=0) {
      for (int j=ii; j<=i-1; ++j) sum -= a(i,j)*b(j);
    } else if (sum) {
      ii=i;
    }
    b(i)=sum;
  }
  // Solve U*x=y
  for (int i=n-1; i>=0; --i) {
    sum=b(i);
    for (int j=i+1; j<n; ++j) {
      Real& ap = a(i,j);
      Real& bp = b(j);
      sum -= ap * bp;
    }
    b(i)=sum/a(i,i);
  }
}
