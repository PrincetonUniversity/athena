//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file spherical_grid.cpp
//  \brief implements SphericalGrid and SphericalPatch

#include <cmath>
#include <list>

#include "../coordinates/coordinates.hpp"
#include "spherical_grid.hpp"

#define SQ(X) ((X)*(X))

SphericalGrid::SphericalGrid(int nlev, Real rad_): GeodesicGrid(nlev) {
  rad.NewAthenaArray(GeodesicGrid::NumVertices());
  rad.Fill(rad_);
}

void SphericalGrid::Position(int ic, Real * x, Real * y, Real * z) const {
  Real theta, phi;
  GeodesicGrid::PositionPolar(ic, &theta, &phi);
  *x = rad(ic)*std::sin(theta)*std::cos(phi);
  *y = rad(ic)*std::sin(theta)*std::sin(phi);
  *z = rad(ic)*std::cos(theta);
}

Real SphericalGrid::ComputeWeight(int ic) const {
  return rad(ic)*rad(ic)*GeodesicGrid::ComputeWeight(ic);
}

Real SphericalGrid::ArcLength(int ic1, int ic2) const {
  Real x1, y1, z1;
  Real x2, y2, z2;
  Position(ic1, &x1, &y1, &z1);
  Position(ic2, &x2, &y2, &z2);
  return std::sqrt(SQ(x1 - x2) + SQ(y1 - y2) + SQ(z1 - z2));
}

SphericalPatch::SphericalPatch(SphericalGrid const * psphere, MeshBlock const * pblock, collocation_t coll):
    coll(coll), psphere(psphere), pblock(pblock), pinterp(nullptr) {
  MeshBlock const * pmb = pblock;
  Coordinates const * pmc = pblock->pcoord;

  Real xmin, xmax, ymin, ymax, zmin, zmax;
  Real origin[3];
  Real delta[3];
  int size[3];

  switch (coll) {
    case cell: 
      xmin = pmc->x1v(pmb->is);
      xmax = pmc->x1v(pmb->ie);
      ymin = pmc->x2v(pmb->js);
      ymax = pmc->x2v(pmb->je);
      zmin = pmc->x3v(pmb->ks);
      zmax = pmc->x3v(pmb->ke);
      origin[0] = pmc->x1v(0);
      origin[1] = pmc->x2v(0);
      origin[2] = pmc->x3v(0);
      size[0] = pmb->block_size.nx1 + 2*(NGHOST);
      size[1] = pmb->block_size.nx2 + 2*(NGHOST);
      size[2] = pmb->block_size.nx3 + 2*(NGHOST);
      break;
    case vertex:
      xmin = pmc->x1f(pmb->is);
      xmax = pmc->x1f(pmb->ie+1);
      ymin = pmc->x2f(pmb->js);
      ymax = pmc->x2f(pmb->je+1);
      zmin = pmc->x3f(pmb->ks);
      zmax = pmc->x3f(pmb->ke+1);
      origin[0] = pmc->x1f(0);
      origin[1] = pmc->x2f(0);
      origin[2] = pmc->x3f(0);
      size[0] = pmb->block_size.nx1 + 2*(NGHOST) + 1;
      size[1] = pmb->block_size.nx2 + 2*(NGHOST) + 1;
      size[2] = pmb->block_size.nx3 + 2*(NGHOST) + 1;
      break;
  }
  printf("xmin = %.16f,xmax = %.16f,ymin = %.16f,ymax = %.16f,zmin = %.16f,zmax = %.16f\n",xmin,xmax,ymin,ymax,zmin,zmax);
  delta[0] = pmc->dx1v(0);
  delta[1] = pmc->dx2v(0);
  delta[2] = pmc->dx3v(0);

  // Loop over all points to find those belonging to this spherical patch
  int const np = psphere->NumVertices();
  map.reserve(np);
  for (int ic = 0; ic < np; ++ic) {
    Real x, y, z;
    psphere->Position(ic, &x, &y, &z);
    if (x >= xmin && x < xmax && y >= ymin && y < ymax && z >= zmin && z < zmax) {
      map.push_back(ic);
    }
  }
  map.shrink_to_fit();
  n = map.size();

  pinterp = new LagrangeInterpND<2*NGHOST-1, 3> *[n];
  for (int i = 0; i < n; ++i) {
    Real coord[3];
    psphere->Position(map[i], &coord[0], &coord[1], &coord[2]);
    pinterp[i] = new LagrangeInterpND<2*NGHOST-1, 3>(origin, delta, size, coord);
  }
}

SphericalPatch::~SphericalPatch() {
  for (int i = 0; i < n; ++i) {
    delete pinterp[i];
  }
  delete[] pinterp;
}

void SphericalPatch::interpToSpherical(Real const * src, Real * dst) const {
  for (int i = 0; i < n; ++i) {
    dst[i] = pinterp[i]->eval(src);
  }
}

void SphericalPatch::InterpToSpherical(AthenaArray<Real> const & src, AthenaArray<Real> * dst) const {
  assert (src.GetDim4() == dst->GetDim2());
  assert (dst->GetDim1() == n);
  AthenaArray<Real> mySrc, myDst;
  int const nvars = src.GetDim4();
  for (int iv = 0; iv < nvars; ++iv) {
    mySrc.InitWithShallowSlice(const_cast<AthenaArray<Real>&>(src), iv, 1);
    myDst.InitWithShallowSlice(*dst, iv, 1);
    interpToSpherical(mySrc.data(), myDst.data());
  }
}

void SphericalPatch::mergeData(Real const * src, Real * dst) const {
  for (int i = 0; i < n; ++i) {
    dst[map[i]] = src[i];
  }
}

void SphericalPatch::MergeData(AthenaArray<Real> const & src, AthenaArray<Real> * dst) const {
  assert (src.GetDim2() == dst->GetDim2());
  assert (dst->GetDim1() == psphere->NumVertices());
  assert (src.GetDim1() == n);
  AthenaArray<Real> mySrc, myDst;
  int const nvars = src.GetDim2();
  for (int iv = 0; iv < nvars; ++iv) {
    mySrc.InitWithShallowSlice(const_cast<AthenaArray<Real>&>(src), iv, 1);
    myDst.InitWithShallowSlice(*dst, iv, 1);
    mergeData(mySrc.data(), myDst.data());
  }
}
