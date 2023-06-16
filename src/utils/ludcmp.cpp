//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================

// C headers

// C++ headers
#include <sstream>
#include <stdexcept>

// Athena headers
#include "../athena.hpp"

//======================================================================================
//! \file ludcmp.cpp
//======================================================================================


//--------------------------------------------------------------------------------------
//! \fn  void  Ludcmp_nr(int n, AthenaArray<Real> &a, AthenaArray<int> &indx,
//          Real *d);
// \brief LU decomposition from numerical recipes
// Using Crout's method with partial pivoting
//a is the input matrix, and is returned with LU decomposition readily made,
// is the matrix size, indx records the history of row permutation,
// whereas d =1(-1) for even(odd) number of permutations.

void Ludcmp_nr(int n, AthenaArray<Real> &a, AthenaArray<int> &indx,
               Real *d) {
  int imax;
  Real big,dum,sum,temp;
  AthenaArray<Real> rowscale;  // the implicit scaling of each row
  std::stringstream msg;
  rowscale.NewAthenaArray(n);
  *d=1.0;  // No row interchanges yet

  // Loop over rows to get the implicit scaling information
  for (int i=0; i<n; ++i) {
    big=0.0;
    for (int j=0; j<n; ++j)
      if ((temp=std::abs(a(i,j))) > big) big=temp;
    if (big == 0.0) {
      msg << "### [LUdecomp]:Input matrix is singular!"
      << "### FATAL ERROR in function [LUdecomp_nr]" << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    rowscale(i)=1.0/big;  // Save the scaling
  }

  for (int j=0; j<n; ++j) { // Loop over columns of Crout's method
    // Calculate the upper block
    for (int i=0; i<j; ++i) {
      sum=a(i,j);
      for (int k=0; k<i; ++k) {
        Real& aik = a(i,k);
        Real& akj = a(k,j);
        sum -= aik * akj;
      }
      a(i,j)=sum;
    }
    // Calculate the lower block (first step)
    big=0.0;
    for (int i=j; i<n; ++i) {
      sum=a(i,j);
      for (int k=0; k<j; ++k) {
        Real& aik = a(i,k);
        Real& akj = a(k,j);
        sum -= aik * akj;
      }
      a(i,j)=sum;
      // search for the largest pivot element
      if ( (dum=rowscale(i)*std::abs(sum)) >= big) {
        big=dum;
        imax=i;
      }
    }
    // row interchange
    if (j != imax) {
      for (int k=0; k<n; ++k) {
        dum=a(imax,k);
        a(imax,k)=a(j,k);
        a(j,k)=dum;
      }
      *d = -(*d);
      rowscale(imax)=rowscale(j);
    }
    indx(j)=imax; // record row interchange history
    // Calculate the lower block (second step)
    if (a(j,j) == 0.0) a(j,j)=TINY_NUMBER;
    dum=1.0/(a(j,j));
    for (int i=j+1; i<n; ++i)
      a(i,j) *= dum;
  }
  // free the memory
  rowscale.DeleteAthenaArray();
}
