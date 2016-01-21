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
#include "coordinates.hpp"

// Athena++ headers
#include "../athena.hpp"          // macros, Real
#include "../athena_arrays.hpp"   // AthenaArray
#include "../parameter_input.hpp" // ParameterInput
#include "../mesh.hpp"            // MeshBlock

//======================================================================================
//! \file cartesian.cpp
//  \brief implements Coordinates class functions for Cartesian (x-y-z) coordinates
//======================================================================================

//--------------------------------------------------------------------------------------
// Coordinates constructor

Coordinates::Coordinates(MeshBlock *pmb, ParameterInput *pin, int flag)
{
  pmy_block = pmb;
  cflag=flag;
  int is, ie, js, je, ks, ke, ng;
  if(cflag==0) {
    is = pmb->is; js = pmb->js; ks = pmb->ks;
    ie = pmb->ie; je = pmb->je; ke = pmb->ke;
    ng=NGHOST;
  }
  else {
    is = pmb->cis; js = pmb->cjs; ks = pmb->cks;
    ie = pmb->cie; je = pmb->cje; ke = pmb->cke;
    ng=pmb->cnghost;
  }

  // Set face centered positions and distances
  AllocateAndSetBasicCoordinates();

  // initialize volume-averaged positions and spacing
  // x1-direction: x1v = dx/2
  for (int i=is-ng; i<=ie+ng; ++i) {
    x1v(i) = 0.5*(x1f(i+1) + x1f(i));
  }
  for (int i=is-ng; i<=ie+ng-1; ++i) {
    dx1v(i) = x1v(i+1) - x1v(i);
  }

  // x2-direction: x2v = dy/2
  if (pmb->block_size.nx2 == 1) {
    x2v(js) = 0.5*(x2f(js+1) + x2f(js));
    dx2v(js) = dx2f(js);
  } else {
    for (int j=js-ng; j<=je+ng; ++j) {
      x2v(j) = 0.5*(x2f(j+1) + x2f(j));
    }
    for (int j=js-ng; j<=je+ng-1; ++j) {
      dx2v(j) = x2v(j+1) - x2v(j);
    }
  }

  // x3-direction: x3v = dz/2
  if (pmb->block_size.nx3 == 1) {
    x3v(ks) = 0.5*(x3f(ks+1) + x3f(ks));
    dx3v(ks) = dx3f(ks);
  } else {
    for (int k=ks-ng; k<=ke+ng; ++k) {
      x3v(k) = 0.5*(x3f(k+1) + x3f(k));
    }
    for (int k=ks-ng; k<=ke+ng-1; ++k) {
      dx3v(k) = x3v(k+1) - x3v(k);
    }
  }

  if((pmb->pmy_mesh->multilevel==true) && MAGNETIC_FIELDS_ENABLED) {
    for (int i=is-ng; i<=ie+ng; ++i)
      x1s2(i) = x1s3(i) = x1v(i);
    if (pmb->block_size.nx2 == 1) {
      x2s1(js) = x2s3(js) = x2v(js);
    }
    else {
      for (int j=js-ng; j<=je+ng; ++j)
        x2s1(j) = x2s3(j) = x2v(j);
    }
    if (pmb->block_size.nx3 == 1) {
      x3s1(ks) = x3s2(ks) = x3v(ks);
    }
    else {
      for (int k=ks-ng; k<=ke+ng; ++k)
        x3s1(k) = x3s2(k) = x3v(k);
    }
  }
}

// destructor
Coordinates::~Coordinates()
{
  DeleteBasicCoordinates();
}

//--------------------------------------------------------------------------------------
// Edge Length functions: returns physical length at cell edges
// Edge1(i,j,k) located at (i,j-1/2,k-1/2), i.e. (x1v(i), x2f(j), x3f(k))

void Coordinates::Edge1Length(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = dx1f(i);
  }
  return;
}

// Edge2(i,j,k) located at (i-1/2,j,k-1/2), i.e. (x1f(i), x2v(j), x3f(k))

void Coordinates::Edge2Length(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = dx2f(j);
  }
  return;
}

// Edge3(i,j,k) located at (i-1/2,j-1/2,k), i.e. (x1f(i), x2f(j), x3v(k))

void Coordinates::Edge3Length(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = dx3f(k);
  }
  return;
}


// GetEdge?Length functions: return one edge length at i
Real Coordinates::GetEdge2Length(const int k, const int j, const int i)
{
  return dx2f(j);
}

Real Coordinates::GetEdge3Length(const int k, const int j, const int i)
{
  return dx3f(k);
}

//--------------------------------------------------------------------------------------
// Cell-center Width functions: returns physical width at cell-center

Real Coordinates::CenterWidth1(const int k, const int j, const int i)
{
  return dx1f(i);
}

Real Coordinates::CenterWidth2(const int k, const int j, const int i)
{
  return dx2f(j);
}

Real Coordinates::CenterWidth3(const int k, const int j, const int i)
{
  return dx3f(k);
}

//--------------------------------------------------------------------------------------
// Face Area functions

void Coordinates::Face1Area(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    // area1 = dy dz
    Real& area_i = area(i);
    area_i = dx2f(j)*dx3f(k);
  }
  return;
}

void Coordinates::Face2Area(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    // area2 = dx dz
    Real& area_i = area(i);
    area_i = dx1f(i)*dx3f(k);
  }
  return;
}

void Coordinates::Face3Area(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    // area3 = dx dy
    Real& area_i = area(i);
    area_i = dx1f(i)*dx2f(j);
  }
  return;
}


// GetFace1Area returns only one Face1Area at i
Real Coordinates::GetFace1Area(const int k, const int j, const int i)
{
  return dx2f(j)*dx3f(k);
}



//--------------------------------------------------------------------------------------
// Cell Volume function

void Coordinates::CellVolume(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &vol)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    // volume = dx dy dz
    Real& vol_i = vol(i);
    vol_i = dx1f(i)*dx2f(j)*dx3f(k);
  }
  return;
}

// GetCellVolume returns only one CellVolume at i
Real Coordinates::GetCellVolume(const int k, const int j, const int i)
{
  return dx1f(i)*dx2f(j)*dx3f(k);
}

//--------------------------------------------------------------------------------------
// Coordinate (Geometric) source term functions

void Coordinates::CoordSrcTerms(const Real dt, const AthenaArray<Real> *flux,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u)
{
  return;
}


//-------------------------------------------------------------------------------------
// Calculate Divv 
void Coordinates::Divv(const AthenaArray<Real> &prim, AthenaArray<Real> &divv)
{

  int is = pmy_block->is; int js = pmy_block->js; int ks = pmy_block->ks;
  int ie = pmy_block->ie; int je = pmy_block->je; int ke = pmy_block->ke;
  int il = is-1; int iu = ie+1;
  int jl, ju, kl, ku;
  Real area_p1, area, vol;
  Real vel_p1, vel;

  if(pmy_block->block_size.nx2 == 1) // 1D
    jl=js, ju=je, kl=ks, ku=ke;
  else if(pmy_block->block_size.nx3 == 1) // 2D
    jl=js-1, ju=je+1, kl=ks, ku=ke;
  else // 3D
    jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;

  for (int k=kl; k<=ku; ++k){
    for (int j=jl; j<=ju; ++j){
#pragma simd
      for (int i=il; i<=iu; ++i){
        area_p1 = dx2f(j)*dx3f(k);
        area = area_p1; 
        vel_p1 = 0.5*(prim(IM1,k,j,i+1)+prim(IM1,k,j,i));
        vel = 0.5*(prim(IM1,k,j,i)+prim(IM1,k,j,i-1));
        divv(k,j,i) = area_p1*vel_p1 - area*vel;
      }
      if (pmy_block->block_size.nx2 > 1) {
        for (int i=il; i<=iu; ++i){
          area_p1 = dx1f(i)*dx3f(k); 
          area = area_p1;
          vel_p1 = 0.5*(prim(IM2,k,j+1,i)+prim(IM2,k,j,i));
          vel = 0.5*(prim(IM2,k,j,i)+prim(IM2,k,j-1,i));
          divv(k,j,i) += area_p1*vel_p1 - area*vel;
        }
      }
      if (pmy_block->block_size.nx3 > 1) {
        for (int i=il; i<=iu; ++i){
          area_p1 = dx1f(i)*dx2f(j); 
          area = area_p1;
          vel_p1 = 0.5*(prim(IM3,k+1,j,i)+prim(IM3,k,j,i));
          vel = 0.5*(prim(IM3,k,j,i)+prim(IM3,k-1,j,i));
          divv(k,j,i) += area_p1*vel_p1 - area*vel;
        }
      }
      for (int i=il; i<=iu; ++i){
        vol = dx1f(i)*dx2f(j)*dx3f(k); 
        divv(k,j,i) = divv(k,j,i)/vol;
      }
    }
  }

  return;
}


// v_{x1;x1}  covariant derivative at x1 interface
void Coordinates::FaceXdx(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = (prim(IM1,k,j,i)-prim(IM1,k,j,i-1))/dx1v(i-1);
  }
  return;
}

// v_{x2;x1}+v_{x1;x2}  covariant derivative at x1 interface
void Coordinates::FaceXdy(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = (prim(IM2,k,j,i)-prim(IM2,k,j,i-1))/dx1v(i-1)
             +0.5*(prim(IM1,k,j+1,i)+prim(IM1,k,j+1,i-1)-prim(IM1,k,j-1,i)-prim(IM1,k,j-1,i-1))
              /(dx2v(j-1)+dx2v(j));
  }
  return;
}

// v_{x3;x1}+v_{x1;x3}  covariant derivative at x1 interface
void Coordinates::FaceXdz(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = (prim(IM3,k,j,i)-prim(IM3,k,j,i-1))/dx1v(i-1)
             +0.5*(prim(IM1,k+1,j,i)+prim(IM1,k+1,j,i-1)-prim(IM1,k-1,j,i)-prim(IM1,k-1,j,i-1))
              /(dx3v(k-1)+dx3v(k));
  }
  return;
}

// v_{x1;x2}+v_{x2;x1}  covariant derivative at x2 interface
void Coordinates::FaceYdx(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = (prim(IM1,k,j,i)-prim(IM1,k,j-1,i))/dx2v(j-1)
             +0.5*(prim(IM2,k,j,i+1)+prim(IM2,k,j-1,i+1)-prim(IM2,k,j,i-1)-prim(IM2,k,j-1,i-1))
              /(dx1v(i-1)+dx1v(i));
  }
  return;
}

// v_{x2;x2}  covariant derivative at x2 interface
void Coordinates::FaceYdy(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = (prim(IM2,k,j,i)-prim(IM2,k,j-1,i))/dx2v(j-1);
  }
  return;
}

// v_{x3;x2}+v_{x2;x3}  covariant derivative at x2 interface
void Coordinates::FaceYdz(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = (prim(IM3,k,j,i)-prim(IM3,k,j-1,i))/dx2v(j-1)
             +0.5*(prim(IM2,k+1,j,i)+prim(IM2,k+1,j-1,i)-prim(IM2,k-1,j,i)-prim(IM2,k-1,j-1,i))
              /(dx3v(k-1)+dx3v(k));
  }
  return;
}

// v_{x1;x3}+v_{x3;x1}  covariant derivative at x3 interface
void Coordinates::FaceZdx(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = (prim(IM1,k,j,i)-prim(IM1,k-1,j,i))/dx3v(k-1)
             +0.5*(prim(IM3,k,j,i+1)+prim(IM3,k-1,j,i+1)-prim(IM3,k,j,i-1)-prim(IM3,k-1,j,i-1))
              /(dx1v(i-1)+dx1v(i));
  }
  return;
}

// v_{x2;x3}+v_{x3;x2}  covariant derivative at x3 interface
void Coordinates::FaceZdy(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = (prim(IM2,k,j,i)-prim(IM2,k-1,j,i))/dx3v(k-1)
             +0.5*(prim(IM3,k,j+1,i)+prim(IM3,k-1,j+1,i)-prim(IM3,k,j-1,i)-prim(IM3,k-1,j-1,i))
              /(dx2v(j-1)+dx2v(j));
  } 
  return;
}
    

// v_{x3;x3}  covariant derivative at x3 interface
void Coordinates::FaceZdz(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{   
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = (prim(IM3,k,j,i)-prim(IM3,k-1,j,i))/dx3v(k-1);
  }
  return;
}

//--------------------------------------------------------------------------------------
// Viscous (Geometric) source term functions

void Coordinates::VisSrcTermsX1(const int k, const int j, const Real dt,
  const AthenaArray<Real> &flx,
  const AthenaArray<Real> &prim, AthenaArray<Real> &u)
{
  return;
}

void Coordinates::VisSrcTermsX2(const int k, const int j, const Real dt,
  const AthenaArray<Real> &flx,  const AthenaArray<Real> &flx_p1,
  const AthenaArray<Real> &prim, AthenaArray<Real> &u)
{
  return;
}

void Coordinates::VisSrcTermsX3(const int k, const int j, const Real dt,
  const AthenaArray<Real> &flx,  const AthenaArray<Real> &flx_p1,
  const AthenaArray<Real> &prim, AthenaArray<Real> &u)
{
  return;
}

