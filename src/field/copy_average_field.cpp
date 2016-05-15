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
//! \file copy_average_field.cpp
//  \brief Copies or averages face-centered field, for use in RKX time integrators
//======================================================================================

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

// this class header
#include "field.hpp"

//--------------------------------------------------------------------------------------
// \!fn 
// \brief returns A = (1-f)B + fC, where A,B and C are equal-size AthenaArrays

void Field::CopyOrAverageField(FaceField &a, FaceField &b, FaceField &c, Real factor)
{
  if (factor == 0.0) {
    a.x1f = b.x1f;
    a.x2f = b.x2f;
    a.x3f = b.x3f;
  } else {
    int size = a.x1f.GetSize();
    for (int i=0; i<size; ++i) {
      a.x1f(i) = b.x1f(i) + factor*(c.x1f(i) - b.x1f(i));
      a.x2f(i) = b.x2f(i) + factor*(c.x2f(i) - b.x2f(i));
      a.x3f(i) = b.x3f(i) + factor*(c.x3f(i) - b.x3f(i));
    }
  }
}
