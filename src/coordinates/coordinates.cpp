//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file coordinates.cpp
//  \brief implements functions for Coordinates abstract base class

// Athena++ headers
#include "coordinates.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"

//----------------------------------------------------------------------------------------
// Coordinates constructor: sets coordinates and coordinate spacing of cell FACES

Coordinates::Coordinates(MeshBlock *pmb, ParameterInput *pin, bool flag)
{
  pmy_block = pmb;
  coarse_flag=flag;
  int is, ie, js, je, ks, ke, ng;
  if(coarse_flag==true) {
    is = pmb->cis; js = pmb->cjs; ks = pmb->cks;
    ie = pmb->cie; je = pmb->cje; ke = pmb->cke;
    ng=pmb->cnghost;
  } else {
    is = pmb->is; js = pmb->js; ks = pmb->ks;
    ie = pmb->ie; je = pmb->je; ke = pmb->ke;
    ng=NGHOST;
  }
  Mesh *pm=pmy_block->pmy_mesh;
  RegionSize& mesh_size  = pmy_block->pmy_mesh->mesh_size;
  RegionSize& block_size = pmy_block->block_size;

  // allocate arrays for face-centered coordinates and coordinate spacing
  int ncells1 = (ie-is+1) + 2*ng;
  int ncells2 = 1, ncells3 = 1;
  if (block_size.nx2 > 1) ncells2 = (je-js+1) + 2*ng;
  if (block_size.nx3 > 1) ncells3 = (ke-ks+1) + 2*ng;

  // note extra cell for face-positions
  dx1f.NewAthenaArray(ncells1);
  dx2f.NewAthenaArray(ncells2);
  dx3f.NewAthenaArray(ncells3);
  x1f.NewAthenaArray((ncells1+1));
  x2f.NewAthenaArray((ncells2+1));
  x3f.NewAthenaArray((ncells3+1));

  long long nrootmesh, noffset;
  long int &lx1=pmy_block->loc.lx1;
  long int &lx2=pmy_block->loc.lx2;
  long int &lx3=pmy_block->loc.lx3;
  int &ll=pmy_block->loc.level;

//--- X1-DIRECTION: initialize coordinates and spacing of cell FACES (x1f,dx1f)

  nrootmesh=mesh_size.nx1*(1L<<(ll-pm->root_level));

  if(pm->use_meshgen_fn_[X1DIR]==true) { // use nonuniform or user-defined meshgen fn
    for (int i=is-ng; i<=ie+ng+1; ++i) {
      // if there are too many levels, this won't work or be precise enough
      if (coarse_flag == false) {
        noffset = i-is + (long long)lx1*block_size.nx1;
      } else {
        noffset = (i-is)*2 + (long long)lx1*block_size.nx1;
      }
      Real rx=(Real)noffset/(Real)nrootmesh;
      x1f(i)=pm->MeshGenerator_[X1DIR](rx,mesh_size);
    }
    x1f(is) = block_size.x1min;
    x1f(ie+1) = block_size.x1max;
    for(int i=is-ng; i<=ie+ng; ++i) {
      dx1f(i)=x1f(i+1)-x1f(i);
    }

    // check that coordinate spacing is reasonable
    Real rmax=1.0, rmin=1.0;
    for(int i=is; i<=ie; i++) {
      rmax=std::max(dx1f(i+1)/dx1f(i),rmax);
      rmin=std::min(dx1f(i+1)/dx1f(i),rmin);
    }
    if(rmax > 1.1 || rmin  < 1.0/1.1) {
       std::cout << "### Warning in Coordinates constructor" << std::endl
         << "Neighboring cell sizes differ by more than 10% in the x1 direction."
         << std::endl;
    }

  } else {  // uniform grid
    Real dx=(block_size.x1max-block_size.x1min)/(ie-is+1);
    for(int i=is-ng; i<=ie+ng; ++i) {
      dx1f(i)=dx;
    }
    x1f(is-ng)=block_size.x1min-ng*dx;
    for(int i=is-ng+1;i<=ie+ng+1;i++) {
      x1f(i)=x1f(i-1)+dx;
    }
    x1f(is) = block_size.x1min;
    x1f(ie+1) = block_size.x1max;
  }

  // correct cell face coordinates in ghost zones for reflecting boundary condition
  if (pmy_block->block_bcs[INNER_X1] == REFLECTING_BNDRY) {
    for (int i=1; i<=ng; ++i) {
      dx1f(is-i) = dx1f(is+i-1);
       x1f(is-i) =  x1f(is-i+1) - dx1f(is-i);
    }
  }
  if (pmy_block->block_bcs[OUTER_X1] == REFLECTING_BNDRY) {
    for (int i=1; i<=ng; ++i) {
      dx1f(ie+i  ) = dx1f(ie-i+1);
       x1f(ie+i+1) =  x1f(ie+i) + dx1f(ie+i);
    }
  }

//--- X2-DIRECTION: initialize coordinates and spacing of cell FACES (x2f,dx2f)

  if(ncells2 > 1) {

    nrootmesh=mesh_size.nx2*(1L<<(ll-pm->root_level));

    if(pm->use_meshgen_fn_[X2DIR]==true) { // use nonuniform or user-defined meshgen fn
      for (int j=js-ng; j<=je+ng+1; ++j) {
        // if there are too many levels, this won't work or be precise enough
        if (coarse_flag == false) {
          noffset = j-js + (long long)lx2*block_size.nx2;
        } else {
          noffset = (j-js)*2 + (long long)lx2*block_size.nx2;
        }
        Real rx=(Real)noffset/(Real)nrootmesh;
        x2f(j)=pm->MeshGenerator_[X2DIR](rx,mesh_size);
      }
      x2f(js) = block_size.x2min;
      x2f(je+1) = block_size.x2max;
      for(int j=js-ng; j<=je+ng; ++j) {
        dx2f(j)=x2f(j+1)-x2f(j);
      }

      // check that coordinate spacing is reasonable
      Real rmax=1.0, rmin=1.0;
      for(int j=pmy_block->js; j<=pmy_block->je; j++) {
        rmax=std::max(dx2f(j+1)/dx2f(j),rmax);
        rmin=std::min(dx2f(j+1)/dx2f(j),rmin);
      }
      if(rmax > 1.1 || rmin  < 1.0/1.1) {
         std::cout << "### Warning in Coordinates constructor" << std::endl
           << "Neighboring cell sizes differ by more than 10% in the x2 direction."
           << std::endl;
      }

    } else {  // uniform grid
      Real dx=(block_size.x2max-block_size.x2min)/(je-js+1);
      for(int j=js-ng; j<=je+ng; ++j) {
        dx2f(j)=dx;
      }
      x2f(js-ng)=block_size.x2min-ng*dx;
      for(int j=js-ng+1;j<=je+ng+1;j++) {
        x2f(j)=x2f(j-1)+dx;
      }
      x2f(js) = block_size.x2min;
      x2f(je+1) = block_size.x2max;
    }

    // correct cell face coordinates in ghost zones for reflecting boundary condition
    if (pmy_block->block_bcs[INNER_X2] == REFLECTING_BNDRY
     || pmy_block->block_bcs[INNER_X2] == POLAR_BNDRY) { // also polar boundary
      for (int j=1; j<=ng; ++j) {
        dx2f(js-j) = dx2f(js+j-1);
         x2f(js-j) =  x2f(js-j+1) - dx2f(js-j);
      }
    }
    if (pmy_block->block_bcs[OUTER_X2] == REFLECTING_BNDRY
     || pmy_block->block_bcs[OUTER_X2] == POLAR_BNDRY) { // also polar boundary
      for (int j=1; j<=ng; ++j) {
        dx2f(je+j  ) = dx2f(je-j+1);
         x2f(je+j+1) =  x2f(je+j) + dx2f(je+j);
      }
    }

  // 1D problem
  } else {
    dx2f(js  ) = block_size.x2max-block_size.x2min;
    x2f (js  ) = block_size.x2min;
    x2f (je+1) = block_size.x2max;
  }

//--- X3-DIRECTION: initialize coordinates and spacing of cell FACES (x3f,dx3f)

  if(ncells3 > 1) {

    nrootmesh=mesh_size.nx3*(1L<<(ll-pm->root_level));

    if(pm->use_meshgen_fn_[X3DIR]==true) {  // use nonuniform or user-defined meshgen fn
      for (int k=ks-ng; k<=ke+ng+1; ++k) {
        // if there are too many levels, this won't work or be precise enough
        if (coarse_flag == false) {
          noffset = k-ks + (long long)lx3*block_size.nx3;
        } else {
          noffset = (k-ks)*2 + (long long)lx3*block_size.nx3;
        }
        Real rx=(Real)noffset/(Real)nrootmesh;
        x3f(k)=pm->MeshGenerator_[X3DIR](rx,mesh_size);
      }
      x3f(ks) = block_size.x3min;
      x3f(ke+1) = block_size.x3max;
      for(int k=ks-ng; k<=ke+ng; ++k) {
        dx3f(k)=x3f(k+1)-x3f(k);
      }

      // check that coordinate spacing is reasonable
      Real rmax=1.0, rmin=1.0;
      for(int k=pmy_block->ks; k<=pmy_block->ke; k++) {
        rmax=std::max(dx3f(k+1)/dx3f(k),rmax);
        rmin=std::min(dx3f(k+1)/dx3f(k),rmin);
      }
      if(rmax > 1.1 || rmin  < 1.0/1.1) {
         std::cout << "### Warning in Coordinates constructor" << std::endl
           << "Neighboring cell sizes differ by more than 10% in the x3 direction."
           << std::endl;
      }

    } else { // uniform grid
      Real dx=(block_size.x3max-block_size.x3min)/(ke-ks+1);
      for(int k=ks-ng; k<=ke+ng; ++k) {
        dx3f(k)=dx;
      }
      x3f(ks-ng)=block_size.x3min-ng*dx;
      for(int k=ks-ng+1;k<=ke+ng+1;k++) {
        x3f(k)=x3f(k-1)+dx;
      }
      x3f(ks) = block_size.x3min;
      x3f(ke+1) = block_size.x3max;
    }

    // correct cell face coordinates in ghost zones for reflecting boundary condition
    if (pmy_block->block_bcs[INNER_X3] == REFLECTING_BNDRY) {
      for (int k=1; k<=ng; ++k) {
        dx3f(ks-k) = dx3f(ks+k-1);
         x3f(ks-k) =  x3f(ks-k+1) - dx3f(ks-k);
      }
    }
    if (pmy_block->block_bcs[OUTER_X3] == REFLECTING_BNDRY) {
      for (int k=1; k<=ng; ++k) {
        dx3f(ke+k  ) = dx3f(ke-k+1);
         x3f(ke+k+1) =  x3f(ke+k) + dx3f(ke+k);
      }
    }

  // 1D or 2D problem
  } else {
    dx3f(ks) = block_size.x3max-block_size.x3min;
    x3f(ks  ) = block_size.x3min;
    x3f(ke+1) = block_size.x3max;
  }

}

// destructor

Coordinates::~Coordinates()
{
  dx1f.DeleteAthenaArray();
  dx2f.DeleteAthenaArray();
  dx3f.DeleteAthenaArray();
  x1f.DeleteAthenaArray();
  x2f.DeleteAthenaArray();
  x3f.DeleteAthenaArray();
}

//----------------------------------------------------------------------------------------
// EdgeXLength functions: compute physical length at cell edge-X as vector
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


//----------------------------------------------------------------------------------------
// GetEdgeXLength functions: return length of edge-X at (i,j,k)

Real Coordinates::GetEdge1Length(const int k, const int j, const int i)
{
  return dx1f(i);
}

Real Coordinates::GetEdge2Length(const int k, const int j, const int i)
{
  return dx2f(j);
}

Real Coordinates::GetEdge3Length(const int k, const int j, const int i)
{
  return dx3f(k);
}

//----------------------------------------------------------------------------------------
// CenterWidthX functions: return physical width in X-dir at (i,j,k) cell-center

void Coordinates::CenterWidth1(const int k, const int j, const int il, const int iu,
                               AthenaArray<Real> &dx1)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    dx1(i) = dx1f(i);
  }
  return;
}

void Coordinates::CenterWidth2(const int k, const int j, const int il, const int iu,
                               AthenaArray<Real> &dx2)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    dx2(i) = dx2f(j);
  }
  return;
}

void Coordinates::CenterWidth3(const int k, const int j, const int il, const int iu,
                               AthenaArray<Real> &dx3)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    dx3(i) = dx3f(k);
  }
  return;
}

//----------------------------------------------------------------------------------------
// FaceXArea functions: compute area of face with normal in X-dir as vector

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

//----------------------------------------------------------------------------------------
// GetFaceXArea functions: return area of face with normal in X-dir at (i,j,k)

Real Coordinates::GetFace1Area(const int k, const int j, const int i)
{
  return dx2f(j)*dx3f(k);
}

Real Coordinates::GetFace2Area(const int k, const int j, const int i)
{
  return dx1f(i)*dx3f(k);
}

Real Coordinates::GetFace3Area(const int k, const int j, const int i)
{
  return dx1f(i)*dx2f(j);
}

//----------------------------------------------------------------------------------------
// Cell Volume function: compute volume of cell as vector

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

//----------------------------------------------------------------------------------------
// GetCellVolume: returns cell volume at (i,j,k)

Real Coordinates::GetCellVolume(const int k, const int j, const int i)
{
  return dx1f(i)*dx2f(j)*dx3f(k);
}

//----------------------------------------------------------------------------------------
// Coordinate (Geometric) source term function

void Coordinates::CoordSrcTerms(const Real dt, const AthenaArray<Real> *flux,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u)
{
  return;
}

//----------------------------------------------------------------------------------------
// Function for determining if index corresponds to a polar boundary
// Inputs:
//   j: x2-index
// Outputs:
//   returned value: true if face indexed with j is on a pole; false otherwise

bool Coordinates::IsPole(int j)
{
  if (pmy_block->block_bcs[INNER_X2] == POLAR_BNDRY and j == pmy_block->js) {
    return true;
  }
  if (pmy_block->block_bcs[OUTER_X2] == POLAR_BNDRY and j == pmy_block->je+1) {
    return true;
  }
  return false;
}
