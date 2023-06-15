//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================

// Athena headers
#include "../athena.hpp"


//======================================================================================
//! \file matrixmult.cpp
//======================================================================================

//--------------------------------------------------------------------------------------
//! \fn  void MatrixMult(int m, int n, AthenaArray<Real> &a,
//    AthenaArray<Real> &b, AthenaArray<Real> &c)
// \brief perform a(m,n)*b(n)=c(m)

void MatrixMult(int m, int n, AthenaArray<Real> &a,
                AthenaArray<Real> &b, AthenaArray<Real> &c) {
  for (int i=0; i<m; ++i) {
    c(i) = 0.0;
    for (int j=0; j<n; ++j) {
      Real& ap = a(i,j);
      Real& bp = b(j);
      c(i) += ap * bp;
    }
  }
}
