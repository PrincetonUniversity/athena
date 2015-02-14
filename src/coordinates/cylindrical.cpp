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

//Primary header
#include "coordinates.hpp"

// C headers
#include <math.h>   // pow function

// Athena++ headers
#include "../athena.hpp"          // macros, Real
#include "../athena_arrays.hpp"   // AthenaArray
#include "../parameter_input.hpp" // ParameterInput
#include "../mesh.hpp"            // MeshBlock
#include "../fluid/fluid.hpp"     // Fluid
#include "../fluid/eos/eos.hpp"   // SoundSpeed()

//======================================================================================
//! \file cylindrical.cpp
//  \brief implements Coordinates class functions for cylindrical (r-phi-z) coordinates
//======================================================================================

//--------------------------------------------------------------------------------------
// Coordinates constructor

Coordinates::Coordinates(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block = pmb;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  // initialize volume-averaged positions and spacing
  // x1-direction: x1v = (\int r dV / \int dV) = d(r^3/3)d(r^2/2)
  for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
    pmb->x1v(i) = (2.0/3.0)*(pow(pmb->x1f(i+1),3) - pow(pmb->x1f(i),3))
                     /(pow(pmb->x1f(i+1),2) - pow(pmb->x1f(i),2));
  }
  for (int i=is-(NGHOST); i<=ie+(NGHOST)-1; ++i) {
    pmb->dx1v(i) = pmb->x1v(i+1) - pmb->x1v(i);
  }

  // x2-direction: x2v = (\int phi dV / \int dV) = dphi/2
  if (pmb->block_size.nx2 == 1) {
    pmb->x2v(js) = 0.5*(pmb->x2f(js+1) + pmb->x2f(js));
    pmb->dx2v(js) = pmb->dx2f(js);
  } else {
    for (int j=js-(NGHOST); j<=je+(NGHOST); ++j) {
      pmb->x2v(j) = 0.5*(pmb->x2f(j+1) + pmb->x2f(j));
    }
    for (int j=js-(NGHOST); j<=je+(NGHOST)-1; ++j) {
      pmb->dx2v(j) = pmb->x2v(j+1) - pmb->x2v(j);
    }
  }

  // x3-direction: x3v = (\int z dV / \int dV) = dz/2
  if (pmb->block_size.nx3 == 1) {
    pmb->x3v(ks) = 0.5*(pmb->x3f(ks+1) + pmb->x3f(ks));
    pmb->dx3v(ks) = pmb->dx3f(ks);
  } else {
    for (int k=ks-(NGHOST); k<=ke+(NGHOST); ++k) {
      pmb->x3v(k) = 0.5*(pmb->x3f(k+1) + pmb->x3f(k));
    }
    for (int k=ks-(NGHOST); k<=ke+(NGHOST)-1; ++k) {
      pmb->dx3v(k) = pmb->x3v(k+1) - pmb->x3v(k);
    }
  }

  // Allocate memory for scratch arrays used in integrator, and internal scratch arrays
  int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
  coord_area3_i_.NewAthenaArray(ncells1);
  coord_vol_i_.NewAthenaArray(ncells1);
  coord_src1_i_.NewAthenaArray(ncells1);
  coord_src2_i_.NewAthenaArray(ncells1);

  // Compute and store constant coefficients needed for face-areas, cell-volumes, etc.
  // This helps improve performance.
#pragma simd
  for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i){
    Real rm = pmb->x1f(i  );
    Real rp = pmb->x1f(i+1);
    // dV = 0.5*(R_{i+1}^2 - R_{i}^2)
    coord_area3_i_(i)= 0.5*(rp*rp - rm*rm);
    // dV = 0.5*(R_{i+1}^2 - R_{i}^2)
    coord_vol_i_(i) = coord_area3_i_(i);
    // (A1^{+} - A1^{-})/dV
    coord_src1_i_(i) = pmb->dx1f(i)/coord_vol_i_(i);
    // (dR/2)/(R_c dV)
    coord_src2_i_(i) = pmb->dx1f(i)/((rm + rp)*coord_vol_i_(i));
  }

}

// destructor

Coordinates::~Coordinates()
{
  coord_area3_i_.DeleteAthenaArray();
  coord_vol_i_.DeleteAthenaArray();
  coord_src1_i_.DeleteAthenaArray();
  coord_src2_i_.DeleteAthenaArray();
}

//--------------------------------------------------------------------------------------
// Edge Length functions: returns physical length at cell edges
// Edge1(i,j,k) located at (i,j-1/2,k-1/2), i.e. (x1v(i), x2f(j), x3f(k))

void Coordinates::Edge1Length(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = pmy_block->dx1f(i);
  }
  return;
}

// Edge2(i,j,k) located at (i-1/2,j,k-1/2), i.e. (x1f(i), x2v(j), x3f(k))

void Coordinates::Edge2Length(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = pmy_block->x1f(i)*pmy_block->dx2f(j);
  }
  return;
}

// Edge3(i,j,k) located at (i-1/2,j-1/2,k), i.e. (x1f(i), x2f(j), x3v(k))

void Coordinates::Edge3Length(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = pmy_block->dx3f(k);
  }
  return;
}

//--------------------------------------------------------------------------------------
// Cell-center Width functions: returns physical width at cell-center

Real Coordinates::CenterWidth1(const int k, const int j, const int i)
{
  return (pmy_block->dx1f(i));
}

Real Coordinates::CenterWidth2(const int k, const int j, const int i)
{
  return (pmy_block->x1v(i)*pmy_block->dx2f(j));
}

Real Coordinates::CenterWidth3(const int k, const int j, const int i)
{
  return (pmy_block->dx3f(k));
}


//--------------------------------------------------------------------------------------
// Face Area functions

void Coordinates::Face1Area(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    // area1 = r dphi dz 
    Real& area_i = area(i);
    area_i = (pmy_block->x1f(i)*pmy_block->dx2f(j))*(pmy_block->dx3f(k));
  }
  return;
}

void Coordinates::Face2Area(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    // area2 = dr dz
    Real& area_i = area(i);
    area_i = (pmy_block->dx1f(i))*(pmy_block->dx3f(k));
  }
  return;
}

void Coordinates::Face3Area(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    // area3 = dr r dphi = d(r^2/2) dphi
    Real& area_i = area(i);
    area_i = coord_area3_i_(i)*(pmy_block->dx2f(j));
  }
  return;
}

//--------------------------------------------------------------------------------------
// Cell Volume function

void Coordinates::CellVolume(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &vol)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    // volume = dr dz r dphi = d(r^2/2) dphi dz
    Real& vol_i = vol(i);
    vol_i = coord_vol_i_(i)*(pmy_block->dx2f(j))*(pmy_block->dx3f(k));
  }
  return;
}

//--------------------------------------------------------------------------------------
// Coordinate (Geometric) source term functions

void Coordinates::CoordSrcTermsX1(const int k, const int j, const Real dt,
  const AthenaArray<Real> &flx,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u)
{
  Real iso_cs = pmy_block->pfluid->pf_eos->GetIsoSoundSpeed();

#pragma simd
  for (int i=(pmy_block->is); i<=(pmy_block->ie); ++i) {
    // src_1 = <M_{phi phi}><1/r>
    Real m_pp = prim(IDN,k,j,i)*prim(IM2,k,j,i)*prim(IM2,k,j,i);
    if (NON_BAROTROPIC_EOS) {
       m_pp += prim(IEN,k,j,i);
    } else {
       m_pp += (iso_cs*iso_cs)*prim(IDN,k,j,i);
    }
    if (MAGNETIC_FIELDS_ENABLED) {
       m_pp += 0.5*( SQR(bcc(IB1,k,j,i)) - SQR(bcc(IB2,k,j,i)) + SQR(bcc(IB3,k,j,i)) );
    }
    u(IM1,k,j,i) += dt*coord_src1_i_(i)*m_pp;

    // src_2 = -< M_{phi r} ><1/r>
    Real& x_i   = pmy_block->x1f(i);
    Real& x_ip1 = pmy_block->x1f(i+1);
    u(IM2,k,j,i) -= dt*coord_src2_i_(i)*(x_i*flx(IM2,i) + x_ip1*flx(IM2,i+1));
  }

  return;
}

void Coordinates::CoordSrcTermsX2(const int k, const int j, const Real dt,
  const AthenaArray<Real> &flx,  const AthenaArray<Real> &flx_m1,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u)
{
  return;
}

void Coordinates::CoordSrcTermsX3(const int k, const int j, const Real dt,
  const AthenaArray<Real> &flx,  const AthenaArray<Real> &flx_m1,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u)
{
  return;
}
