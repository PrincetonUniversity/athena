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
#include "../integrators.hpp"

// C++ headers
#include <algorithm>  // max(), min()

// Athena headers
#include "../../../athena.hpp"         // enums, macros, Real
#include "../../../athena_arrays.hpp"  // AthenaArray
#include "../../fluid.hpp"             // Fluid
#include "../../eos/eos.hpp"           // GetGamma

//======================================================================================
//! \file ct.cpp
//  \brief


// 1-D update

  for (k=ks; k<=ke; ++k) {
  for (j=js; j<=je+1; ++j) {
    for (i=is; i<=ie; ++i) {
      b2i(k,j,i) += q1*(e3(k,j,i+1) - e3(k,j,i));
    }
  }}
  for (k=ks; k<=ke+1; ++k) {
  for (j=js; j<=je; ++j) {
    for (i=is; i<=ie; ++i) {
      b3i(k,j,i) -= q1*(e2(k,j,i+1) - e2(k,j,i));
    }
  }}

// 2-D update

  for (k=ks; k<=ke; ++k) {
  for (j=js; j<=je; ++j) {
    for (i=is; i<=ie+1; ++i) {
      b1i(k,j,i) -= q2*(e3(k,j+1,i) - e3(k,j,i));
    }
  }}
  for (k=ks; k<=ke+1; ++k) {
  for (j=js; j<=je; ++j) {
    for (i=is; i<=ie; ++i) {
      b3i(k,j,i) += q2*(e1(k,j+1,i) - e1(k,j,i));
    }
  }}

// 3-D update

  for (k=ks; k<=ke; ++k) {
  for (j=js; j<=je; ++j) {
    for (i=is; i<=ie+1; ++i) {
      b1i(k,j,i) += q3*(e2(k+1,j,i) - e2(k,j,i)) -
    }
  }}
  for (k=ks; k<=ke; ++k) {
  for (j=js; j<=je+1; ++j) {
    for (i=is; i<=ie; ++i) {
      b2i(k,j,i) -= q3*(e1(k+1,j,i) - e1(k,j,i));
    }
  }}


/*
  for (k=ks; k<=ke; ++k) {
    for (j=js; j<=je; ++j) {
      for (i=is; i<=ie; ++i) {
        b1i(k,j,i) += q3*(e2(k+1,j,i) - e2(k,j,i)) -
                      q2*(e3(k,j+1,i) - e3(k,j,i));
        b2i(k,j,i) += q1*(e3(k,j,i+1) - e3(k,j,i)) -
                      q3*(e1(k+1,j,i) - e1(k,j,i));
        b3i(k,j,i) += q2*(e1(k,j+1,i) - e1(k,j,i)) -
                      q1*(e2(k,j,i+1) - e2(k,j,i));
      }
    }
  }
*/
