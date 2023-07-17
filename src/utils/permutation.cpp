//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================

// Athena headers
#include "../athena.hpp"


//======================================================================================
//! \file permutation.cpp
//======================================================================================

//--------------------------------------------------------------------------------------
//! \fn int Permutation(int i, int j, int k, int np, AthenaArray<int> &pl)
// \brief permutate the array element

int Permutation(int i, int j, int k, int np, AthenaArray<int> &pl) {
  int ip = -1;
  for (int l=0; l<np; l++) {
    // check each permutation in the table
    for (int m=0; m<3; m++)
      if (i == pl(l,m)) {
        for (int n=0; n<3; n++) {
          if (n != m) {
            if (j == pl(l,n)) {
              for (int o=0; o<3; o++) {
                if ((o != m) && (o != n)) {
                  if (k == pl(l,o)) {
                    ip = l; {
                    }
                  }
                }
              }
            }
          }
        }
      }
  }
  return ip;
}
