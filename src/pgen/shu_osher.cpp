//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 *
 * This program is free software: you can redistribute and/or modify it under the terms
 * of the GNU General Public License (GPL) as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of GNU GPL in the file LICENSE included in
 * the code distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

// Primary header
#include "../fluid.hpp"

// C++ headers
#include <cmath>  // sin()

// Athena headers
#include "../athena.hpp"           // enums, macros, Real
#include "../athena_arrays.hpp"    // AthenaArray
#include "../mesh.hpp"             // Block
#include "../parameter_input.hpp"  // ParameterInput

//======================================================================================
/*! \file shu_osher.cpp
 *  \brief Problem generator for Shu-Osher shocktube test, involving
 *   interaction of a Mach 3 shock with a sine wave density distribution.  
 *
 * REFERENCE: C.W. Shu & S. Osher, "Efficient implementation of essentially
 *   non-oscillatory shock-capturing schemes, II", JCP, 83, 32 (1998)	     
 *====================================================================================*/

void Fluid::InitProblem(ParameterInput *pin)
{
  Block *pb = pparent_block;

  int is = pb->is; int js = pb->js; int ks = pb->ks;
  int ie = pb->ie; int je = pb->je; int ke = pb->ke;

// Read parameters from input file

  gamma_ = pin->GetReal("fluid","gamma");

// setup dependent variables

  Real dl = 3.857143;
  Real pl = 10.33333;
  Real ul = 2.629369;
  Real vl = 0.0;
  Real wl = 0.0;

  Real gm1 = (GetGamma()) - 1.0;

  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
#pragma simd
    for (int i=is; i<=ie; ++i) {

      if (pb->x1v(i) < -0.8) {
        u(IDN,k,j,i) = dl;
        u(IM1,k,j,i) = ul*dl;
        u(IM2,k,j,i) = vl*dl;
        u(IM3,k,j,i) = wl*dl;
        u(IEN,k,j,i) = pl/gm1 + 0.5*dl*(ul*ul + vl*vl + wl*wl);
      }
      else {
        u(IDN,k,j,i) = 1.0 + 0.2*sin(5.0*PI*(pb->x1v(i)));
        u(IM1,k,j,i) = 0.0;
        u(IM2,k,j,i) = 0.0;
        u(IM3,k,j,i) = 0.0;
        u(IEN,k,j,i) = 1.0/gm1;
      }
    }
  }}

  return;
}
