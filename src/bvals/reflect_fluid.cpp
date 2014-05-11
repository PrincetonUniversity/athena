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

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh.hpp"
#include "../fluid.hpp"

//======================================================================================
/*! \file reflect_fluid.cpp
 *  \brief implements reflecting BCs in each dimension for conserved fluid variables
 *====================================================================================*/
//--------------------------------------------------------------------------------------
//! \fn void ReflectInnerX1(Fluid *pf)
//  \brief  REFLECTING boundary conditions conserved vars, inner x1 boundary (ix1_bc=1)

void ReflectInnerX1(Fluid *pf, AthenaArray<Real> &a)
{
  Block *pb = pf->pmy_block;
  int is = pb->is;
  int js = pb->js, je = pb->je;
  int ks = pb->ks, ke = pb->ke;

  AthenaArray<Real> la = a.ShallowCopy();

  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for (int n=0; n<(NVAR); ++n) {

      if (n==(IM1)) {
#pragma simd
        for (int i=1; i<=(NGHOST); ++i) {
          la(IM1,k,j,is-i) = -la(IM1,k,j,(is+i-1));  // reflect 1-mom
        }

      } else {
#pragma simd
        for (int i=1; i<=(NGHOST); ++i) {
          la(n,k,j,is-i) = la(n,k,j,(is+i-1));
        }
      }

    }
  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void ReflectOuterX1(Fluid *pf)
//  \brief  REFLECTING boundary conditions conserved vars, outer x1 boundary (ox1_bc=1)

void ReflectOuterX1(Fluid *pf, AthenaArray<Real> &a)
{
  Block *pb = pf->pmy_block;
  int ie = pb->ie;
  int js = pb->js, je = pb->je;
  int ks = pb->ks, ke = pb->ke;

  AthenaArray<Real> la = a.ShallowCopy();

  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for (int n=0; n<(NVAR); ++n) {

      if (n==(IM1)) {
#pragma simd
        for (int i=1; i<=(NGHOST); ++i) {
          la(IM1,k,j,ie+i) = -la(IM1,k,j,(ie-i+1));  // reflect 1-mom
        }

      } else {
#pragma simd
        for (int i=1; i<=(NGHOST); ++i) {
          la(n,k,j,ie+i) = la(n,k,j,(ie-i+1));
        }
      }

    }
  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void ReflectInnerX2(Fluid *pf)
//  \brief  REFLECTING boundary conditions conserved vars, inner x2 boundary (ix2_bc=1)

void ReflectInnerX2(Fluid *pf, AthenaArray<Real> &a)
{
  Block *pb = pf->pmy_block;
  int is = pb->is, ie = pb->ie;
  int js = pb->js;
  int ks = pb->ks, ke = pb->ke;

  AthenaArray<Real> la = a.ShallowCopy();

  for (int k=ks; k<=ke; ++k) {
  for (int j=1; j<=(NGHOST); ++j) {
    for (int n=0; n<(NVAR); ++n) {

      if (n==(IM2)) {
#pragma simd
        for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
          la(IM2,k,js-j,i) = -la(IM2,k,js+j-1,i);  // reflect 2-mom
        }

      } else {
#pragma simd
        for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
          la(n,k,js-j,i) = la(n,k,js+j-1,i);
        }
      }

    }
  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void ReflectOuterX2(Fluid *pf)
//  \brief  REFLECTING boundary conditions conserved vars, outer x2 boundary (ox2_bc=1)

void ReflectOuterX2(Fluid *pf, AthenaArray<Real> &a)
{
  Block *pb = pf->pmy_block;
  int is = pb->is, ie = pb->ie;
  int je = pb->je;
  int ks = pb->ks, ke = pb->ke;

  AthenaArray<Real> la = a.ShallowCopy();

  for (int k=ks; k<=ke; ++k) {
  for (int j=1; j<=(NGHOST); ++j) {
    for (int n=0; n<(NVAR); ++n) {

      if (n==(IM2)) {
#pragma simd
        for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
          la(IM2,k,je+j,i) = -la(IM2,k,je-j+1,i);  // reflect 2-mom
        }

      } else {
#pragma simd
        for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
          la(n,k,je+j,i) = la(n,k,je-j+1,i);
        }
      }

    }
  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void ReflectInnerX3(Fluid *pf)
//  \brief  REFLECTING boundary conditions conserved vars, inner x3 boundary (ix3_bc=1)

void ReflectInnerX3(Fluid *pf, AthenaArray<Real> &a)
{
  Block *pb = pf->pmy_block;
  int is = pb->is, ie = pb->ie;
  int js = pb->js, je = pb->je;
  int ks = pb->ks;

  AthenaArray<Real> la = a.ShallowCopy();

  for (int k=1; k<=(NGHOST); ++k) {
  for (int j=js-(NGHOST); j<=je+(NGHOST); ++j) {
    for (int n=0; n<(NVAR); ++n) {

      if (n==(IM3)) {
#pragma simd
        for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
          la(IM3,ks-k,j,i) = -la(IM3,ks+k-1,j,i);  // reflect 3-mom
        }

      } else {
#pragma simd
        for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
          la(n,ks-k,j,i) = la(n,ks+k-1,j,i);
        }
      }

    }
  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void ReflectOuterX3(Fluid *pf)
//  \brief  REFLECTING boundary conditions conserved vars, outer x3 boundary (ox3_bc=1)

void ReflectOuterX3(Fluid *pf, AthenaArray<Real> &a)
{
  Block *pb = pf->pmy_block;
  int is = pb->is, ie = pb->ie;
  int js = pb->js, je = pb->je;
  int ke = pb->ke;

  AthenaArray<Real> la = a.ShallowCopy();

  for (int k=1; k<=(NGHOST); ++k) {
  for (int j=js-(NGHOST); j<=je+(NGHOST); ++j) {
    for (int n=0; n<(NVAR); ++n) {

      if (n==(IM3)) {
#pragma simd
        for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
          la(IM3,ke+k,j,i) = -la(IM3,ke-k+1,j,i);  // reflect 3-mom
        }

      } else {
#pragma simd
        for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
          la(n,ke+k,j,i) = la(n,ke-k+1,j,i);
        }
      }

    }
  }}

  return;
}
