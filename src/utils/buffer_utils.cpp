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

//======================================================================================
//! \file buffer_utils.cpp 
//  \brief namespace containing buffer utilities.
//======================================================================================

#include "../athena.hpp"
#include "../athena_arrays.hpp"

#include "buffer_utils.hpp"

namespace BufferUtility
{

//--------------------------------------------------------------------------------------
//! \fn int Pack4DData(AthenaArray<Real> &src, Real *buf, int sn, int en,
//                     int si, int ei, int sj, int ej, int sk, int ek)
//  \brief pack a 4D AthenaArray into a one-dimensional buffer
int Pack4DData(AthenaArray<Real> &src, Real *buf, int sn, int en,
               int si, int ei, int sj, int ej, int sk, int ek)
{
  int p=0;
  for (int n=sn; n<=en; ++n) {
    for (int k=sk; k<=ek; k++) {
      for (int j=sj; j<=ej; j++) {
#pragma simd
        for (int i=si; i<=ei; i++)
            buf[p++]=src(n,k,j,i);
      }
    }
  }
  return p;
}


//--------------------------------------------------------------------------------------
//! \fn void Unpack4DData(Real *buf, AthenaArray<Real> &dst, int sn, int en,
//                                    int si, int ei, int sj, int ej, int sk, int ek)
//  \brief unpack a one-dimensional buffer into a 4D AthenaArray
void Unpack4DData(Real *buf, AthenaArray<Real> &dst, int sn, int en,
                              int si, int ei, int sj, int ej, int sk, int ek)
{
  int p=0;
  for (int n=sn; n<=en; ++n) {
    for (int k=sk; k<=ek; ++k) {
      for (int j=sj; j<=ej; ++j) {
#pragma simd
        for (int i=si; i<=ei; ++i)
          dst(n,k,j,i) = buf[p++];
      }
    }
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn int Pack3DData(AthenaArray<Real> &src, Real *buf,
//                     int si, int ei, int sj, int ej, int sk, int ek)
//  \brief pack a 3D AthenaArray into a one-dimensional buffer
int Pack3DData(AthenaArray<Real> &src, Real *buf,
               int si, int ei, int sj, int ej, int sk, int ek)
{
  int p=0;
  for (int k=sk; k<=ek; k++) {
    for (int j=sj; j<=ej; j++) {
#pragma simd
      for (int i=si; i<=ei; i++)
          buf[p++]=src(k, j, i);
    }
  }
  return p;
}


//--------------------------------------------------------------------------------------
//! \fn void Unpack3DData(Real *buf, AthenaArray<Real> &dst,
//                        int si, int ei, int sj, int ej, int sk, int ek)
//  \brief unpack a one-dimensional buffer into a 3D AthenaArray
void Unpack3DData(Real *buf, AthenaArray<Real> &dst,
                  int si, int ei, int sj, int ej, int sk, int ek)
{
  int p=0;
  for (int k=sk; k<=ek; ++k) {
    for (int j=sj; j<=ej; ++j) {
#pragma simd
      for (int i=si; i<=ei; ++i)
        dst(k,j,i) = buf[p++];
    }
  }
  return;
}

}
