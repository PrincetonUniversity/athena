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
//! \file copy_average_hydro.cpp
//  \brief Copies or averages hydro variables, for us in RKX time integrators
//======================================================================================

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

// this class header
#include "hydro.hpp"

//--------------------------------------------------------------------------------------
// \!fn 
// \brief returns A = (1-f)B + fC, where A,B and C are equal-size AthenaArrays

void Hydro::CopyOrAverageHydro(AthenaArray<Real> &a, AthenaArray<Real> &b,
  AthenaArray<Real> &c, Real factor)
{
  if (factor == 0.0) {
    a = b;
  } else {
    int size = a.GetSize();
    for (int i=0; i<size; ++i) {
      a(i) = b(i) + factor*(c(i) - b(i));
    }
  }
}
