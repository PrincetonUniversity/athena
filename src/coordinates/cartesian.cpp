//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cartesian.cpp
//  \brief implements Coordinates class functions for Cartesian (x-y-z) coordinates

// Athena++ headers
#include "coordinates.hpp"
#include "../athena.hpp"          // macros, Real
#include "../athena_arrays.hpp"   // AthenaArray
#include "../parameter_input.hpp" // ParameterInput
#include "../mesh/mesh.hpp"            // MeshBlock

//----------------------------------------------------------------------------------------
// Cartesian coordinates constructor

Cartesian::Cartesian(MeshBlock *pmb, ParameterInput *pin, int flag)
  : Coordinates(pmb, pin, flag)
{
  pmy_block = pmb;
  cflag=flag;
  int is, ie, js, je, ks, ke, ng;
  if(cflag==0) {
    is = pmb->is; js = pmb->js; ks = pmb->ks;
    ie = pmb->ie; je = pmb->je; ke = pmb->ke;
    ng=NGHOST;
  } else {
    is = pmb->cis; js = pmb->cjs; ks = pmb->cks;
    ie = pmb->cie; je = pmb->cje; ke = pmb->cke;
    ng=pmb->cnghost;
  }
  Mesh *pm=pmy_block->pmy_mesh;
  RegionSize& mesh_size  = pmy_block->pmy_mesh->mesh_size;
  RegionSize& block_size = pmy_block->block_size;

  // allocate arrays for volume-centered coordinates and positions of cells
  int ncells1 = (ie-is+1) + 2*ng;
  int ncells2 = 1, ncells3 = 1;
  if (block_size.nx2 > 1) ncells2 = (je-js+1) + 2*ng;
  if (block_size.nx3 > 1) ncells3 = (ke-ks+1) + 2*ng;
  dx1v.NewAthenaArray(ncells1);
  dx2v.NewAthenaArray(ncells2);
  dx3v.NewAthenaArray(ncells3);
  x1v.NewAthenaArray(ncells1);
  x2v.NewAthenaArray(ncells2);
  x3v.NewAthenaArray(ncells3);
  
  // allocate arrays for area weighted positions for AMR/SMR MHD
  if((pm->multilevel==true) && MAGNETIC_FIELDS_ENABLED) {
    x1s2.NewAthenaArray(ncells1);
    x1s3.NewAthenaArray(ncells1);
    x2s1.NewAthenaArray(ncells2);
    x2s3.NewAthenaArray(ncells2);
    x3s1.NewAthenaArray(ncells3);
    x3s2.NewAthenaArray(ncells3);
  }

  // initialize volume-averaged coordinates and spacing
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

  // initialize area-averaged coordinates used with MHD AMR
  if((pmb->pmy_mesh->multilevel==true) && MAGNETIC_FIELDS_ENABLED) {
    for (int i=is-ng; i<=ie+ng; ++i) {
      x1s2(i) = x1s3(i) = x1v(i);
    }
    if (pmb->block_size.nx2 == 1) {
      x2s1(js) = x2s3(js) = x2v(js);
    } else {
      for (int j=js-ng; j<=je+ng; ++j) {
        x2s1(j) = x2s3(j) = x2v(j);
      }
    }
    if (pmb->block_size.nx3 == 1) {
      x3s1(ks) = x3s2(ks) = x3v(ks);
    } else {
      for (int k=ks-ng; k<=ke+ng; ++k) {
        x3s1(k) = x3s2(k) = x3v(k);
      }
    }
  }
}

//----------------------------------------------------------------------------------------
// Cell Volume function: compute volume of cell as vector

void Cartesian::CellVolume(const int k, const int j, const int il, const int iu,
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

//----------------------------------------------------------------------------------------
// GetCellVolume: returns cell volume at (i,j,k)

Real Cartesian::GetCellVolume(const int k, const int j, const int i)
{
  return dx1f(i)*dx2f(j)*dx3f(k);
}

//----------------------------------------------------------------------------------------
// Coordinate (Geometric) source term function

void Cartesian::CoordSrcTerms(const Real dt, const AthenaArray<Real> *flux,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u)
{
  return;
}
