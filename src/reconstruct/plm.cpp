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

#include <iostream>
#include <string>
#include <math.h>
#include <algorithm>
#include <stdio.h>

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh.hpp"
#include "../fluid.hpp"

//======================================================================================
/*! \file plm.cpp
 *  \brief  piecewise linear reconstruction
 *====================================================================================*/

void Fluid::PiecewiseLinear(const int k, const int j, const int il, const int iu,
  const int dir, AthenaArray<Real> &w, AthenaArray<Real> &wl, AthenaArray<Real> &wr)
{
  int offset;
  int n1 = w.GetDim1();
  int n2 = w.GetDim2();
  int n3 = w.GetDim3();
  if (dir == 1) offset = 1;
  if (dir == 2) offset = n1;
  if (dir == 3) offset = n1*n2;

  for (int n=0; n<NVAR; ++n){
  int zero = n1*(j + n2*(k + n3*n));
#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& wli = wl(n,i);
    Real& wri = wr(n,i);
    Real& wim1 = w(i+zero-offset);
    Real& wi   = w(i+zero);
    Real& wip1 = w(i+zero+offset);

    Real dwl = wi - wim1;
    Real dwr = wip1 - wi;
    Real dw2 = dwl*dwr;

// Apply monotonicity constraints to differences in primitive vars

    Real dwm = dw2/(dwl + dwr + TINY_NUMBER);
    if (dw2 <= 0.0) dwm  = 0.0;
    
// Compute L/R values

    wri = wi + dwm;
    wli = wi - dwm;
  }}

  return;
}
