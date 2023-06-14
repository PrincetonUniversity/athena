//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================

// Athena headers
#include "../athena.hpp"
#include "./utils.hpp"


//======================================================================================
//! \file inverse_matrix.cpp
//======================================================================================


//--------------------------------------------------------------------------------------
//! \fn  void  InverseMatrix(int n, AthenaArray<Real> &a, AthenaArray<Real> &b)
// \brief Inverse matrix solver
// a: input matrix; n: matrix size, b: return matrix
// Note: the input matrix will be DESTROYED

void InverseMatrix(int n, AthenaArray<Real> &a, AthenaArray<Real> &b) {
  AthenaArray<int> indx;
  AthenaArray<Real> col;
  Real d;

  indx.NewAthenaArray(n);
  col.NewAthenaArray(n);
  Ludcmp_nr(n,a,indx,&d);

  for (int j=0; j<n; ++j) {
    for (int i=0; i<n; ++i) {
      col(i)=0.0;
    }
    col(j) = 1.0;
    Lubksb_nr(n, a, indx, col);
    for (int i=0; i<n; i++) {
      b(i,j) = col(i);
    }
  }
  indx.DeleteAthenaArray();
  col.DeleteAthenaArray();
  return;
}
