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

// C headers
#include <math.h>  // pow, trig functions

// Athena++ headers
#include "../athena.hpp"          // macros, Real
#include "../athena_arrays.hpp"   // AthenaArray
#include "../parameter_input.hpp" // ParameterInput
#include "../mesh.hpp"            // MeshBlock
#include "../fluid/fluid.hpp"     // Fluid
#include "../fluid/eos/eos.hpp"   // SoundSpeed()

//======================================================================================
//! \file spherical_polar.cpp
//  \brief implements Coordinates class functions for spherical polar (r-theta-phi)
//======================================================================================

//--------------------------------------------------------------------------------------
// Coordinates constructor

Coordinates::Coordinates(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block = pmb;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  // initialize volume-averaged positions and spacing
  // x1-direction: x1v = (\int r dV / \int dV) = d(r^4/4)/d(r^3/3)
  for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
    pmb->x1v(i) = 0.75*(pow(pmb->x1f(i+1),4) - pow(pmb->x1f(i),4))
                     /(pow(pmb->x1f(i+1),3) - pow(pmb->x1f(i),3));
  }
  for (int i=is-(NGHOST); i<=ie+(NGHOST)-1; ++i) {
    pmb->dx1v(i) = pmb->x1v(i+1) - pmb->x1v(i);
  }

  // x2-direction: x2v = (\int sin[theta] theta dV / \int dV) =
  //   d(sin[theta] - theta cos[theta])/d(-cos[theta])
  if (pmb->block_size.nx2 == 1) {
    pmb->x2v(js) = 0.5*(pmb->x2f(js+1) + pmb->x2f(js));
    pmb->dx2v(js) = pmb->dx2f(js);
  } else {
    for (int j=js-(NGHOST); j<=je+(NGHOST); ++j) {
      pmb->x2v(j) = ((sin(pmb->x2f(j+1)) - pmb->x2f(j+1)*cos(pmb->x2f(j+1))) -
                    (sin(pmb->x2f(j  )) - pmb->x2f(j  )*cos(pmb->x2f(j  ))))/
                         (cos(pmb->x2f(j  )) - cos(pmb->x2f(j+1)));
    }
    for (int j=js-(NGHOST); j<=je+(NGHOST)-1; ++j) {
      pmb->dx2v(j) = pmb->x2v(j+1) - pmb->x2v(j);
    }
  }

  // x3-direction: x3v = (\int phi dV / \int dV) = dphi/2
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
  coord_area1_i_.NewAthenaArray(ncells1);
  coord_area2_i_.NewAthenaArray(ncells1);
  coord_area3_i_.NewAthenaArray(ncells1);
  coord_vol_i_.NewAthenaArray(ncells1);
  coord_src1_i_.NewAthenaArray(ncells1);
  coord_src2_i_.NewAthenaArray(ncells1);

  int ncells2 = 1;
  if (pmb->block_size.nx2 > 1) ncells2 = pmb->block_size.nx2 + 2*(NGHOST);
  coord_area1_j_.NewAthenaArray(ncells2);
  coord_area2_j_.NewAthenaArray(ncells2);
  coord_vol_j_.NewAthenaArray(ncells2);
  coord_src1_j_.NewAthenaArray(ncells2);
  coord_src2_j_.NewAthenaArray(ncells2);

  // Compute and store constant coefficients needed for face-areas, cell-volumes, etc.
  // This helps improve performance.
#pragma simd
  for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i){
    Real rm = pmb->x1f(i  );
    Real rp = pmb->x1f(i+1);
    // R^2
    coord_area1_i_(i) = rm*rm;
    // 0.5*(R_{i+1}^2 - R_{i}^2)
    coord_area2_i_(i) = 0.5*(rp*rp - rm*rm);
    // 0.5*(R_{i+1}^2 - R_{i}^2)
    coord_area3_i_(i) = coord_area2_i_(i);
    // dV = (R_{i+1}^3 - R_{i}^3)/3
    coord_vol_i_(i) = (1.0/3.0)*(rp*rp*rp - rm*rm*rm);
    // (A1^{+} - A1^{-})/dV
    coord_src1_i_(i) = coord_area2_i_(i)/coord_vol_i_(i);
    // (dR/2)/(R_c dV)
    coord_src2_i_(i) = pmb->dx1f(i)/((rm + rp)*coord_vol_i_(i));
  }
  if (pmb->block_size.nx2 > 1) {
#pragma simd
    for (int j=js-(NGHOST); j<=je+(NGHOST); ++j){
      Real sm = sin(pmb->x2f(j  ));
      Real sp = sin(pmb->x2f(j+1));
      Real cm = cos(pmb->x2f(j  ));
      Real cp = cos(pmb->x2f(j+1));
      // d(sin theta) = d(-cos theta)
      coord_area1_j_(j) = cm - cp;
      // sin theta
      coord_area2_j_(j) = sm;
      // d(sin theta) = d(-cos theta)
      coord_vol_j_(j) = coord_area1_j_(j);
      // (A2^{+} - A2^{-})/dV
      coord_src1_j_(j) = (sp - sm)/coord_vol_j_(j);
      // (dS/2)/(S_c dV)
      coord_src2_j_(j) = (sp - sm)/((sm + sp)*coord_vol_j_(j));
    }
  } else {
    Real sm = sin(pmb->x2f(js  ));
    Real sp = sin(pmb->x2f(js+1));
    Real cm = cos(pmb->x2f(js  ));
    Real cp = cos(pmb->x2f(js+1));
    coord_area1_j_(js) = cm - cp;
    coord_area2_j_(js) = sm;
    coord_vol_j_(js) = coord_area1_j_(js);
    coord_src1_j_(js) = (sp - sm)/coord_vol_j_(js);
    coord_src2_j_(js) = (sp - sm)/((sm + sp)*coord_vol_j_(js));
  }

}

// destructor

Coordinates::~Coordinates()
{
  coord_area1_i_.DeleteAthenaArray();
  coord_area2_i_.DeleteAthenaArray();
  coord_area3_i_.DeleteAthenaArray();
  coord_vol_i_.DeleteAthenaArray();
  coord_src1_i_.DeleteAthenaArray();
  coord_src2_i_.DeleteAthenaArray();

  coord_area1_j_.DeleteAthenaArray();
  coord_area2_j_.DeleteAthenaArray();
  coord_vol_j_.DeleteAthenaArray();
  coord_src1_j_.DeleteAthenaArray();
  coord_src2_j_.DeleteAthenaArray();
}

//--------------------------------------------------------------------------------------
// Edge Length functions: returns physical length at cell edges
// Edge1(i,j,k) located at (i,j-1/2,k-1/2), i.e. (x1v(i), x2f(j), x3f(k))

void Coordinates::Edge1Length(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    // length1 = dr
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
    // length2 = r d(theta)
    len(i) = (pmy_block->x1f(i)*pmy_block->dx2f(j));
  }
  return;
}

// Edge3(i,j,k) located at (i-1/2,j-1/2,k), i.e. (x1f(i), x2f(j), x3v(k))

void Coordinates::Edge3Length(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    // length3 = r sin(theta) d(phi)
    len(i) = (pmy_block->x1f(i)*coord_area2_j_(j)*pmy_block->dx3f(k));
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
  return (pmy_block->x1v(i)*sin(pmy_block->x2v(j))*pmy_block->dx3f(k));
}

//--------------------------------------------------------------------------------------
// Face Area functions

void Coordinates::Face1Area(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    // area1 = r^2 sin[theta] dtheta dphi = r^2 d(-cos[theta]) dphi
    area(i) = coord_area1_i_(i)*coord_area1_j_(j)*(pmy_block->dx3f(k)); 
  }
  return;
}

void Coordinates::Face2Area(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    // area2 = dr r sin[theta] dphi = d(r^2/2) sin[theta] dphi
    area(i) = coord_area2_i_(i)*coord_area2_j_(j)*(pmy_block->dx3f(k));
  }
  return;
}

void Coordinates::Face3Area(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    // area3 = dr r dtheta = d(r^2/2) dtheta
    area(i) = coord_area3_i_(i)*(pmy_block->dx2f(j));
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
    // volume = r^2 sin(theta) dr dtheta dphi = d(r^3/3) d(-cos theta) dphi
    vol(i) = coord_vol_i_(i)*coord_vol_j_(j)*(pmy_block->dx3f(k));
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
    // src_1 = < M_{theta theta} + M_{phi phi} ><1/r>
    Real m_ii = prim(IDN,k,j,i)*(SQR(prim(IM2,k,j,i)) + SQR(prim(IM3,k,j,i)));
    if (NON_BAROTROPIC_EOS) {
       m_ii += 2.0*prim(IEN,k,j,i);
    } else {
       m_ii += 2.0*(iso_cs*iso_cs)*prim(IDN,k,j,i);
    }
    if (MAGNETIC_FIELDS_ENABLED) {
       m_ii += SQR(bcc(IB1,k,j,i));
    }
    u(IM1,k,j,i) += dt*coord_src1_i_(i)*m_ii;

    // src_2 = -< M_{theta r} ><1/r> 
    u(IM2,k,j,i) -= dt*coord_src2_i_(i)*
      (coord_area1_i_(i)*flx(IM2,i) + coord_area1_i_(i+1)*flx(IM2,i+1));

    // src_3 = -< M_{phi r} ><1/r> 
    u(IM3,k,j,i) -= dt*coord_src2_i_(i)*
      (coord_area1_i_(i)*flx(IM3,i) + coord_area1_i_(i+1)*flx(IM3,i+1));
  }

  return;
}

void Coordinates::CoordSrcTermsX2(const int k, const int j, const Real dt,
  const AthenaArray<Real> &flx,  const AthenaArray<Real> &flx_p1,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u)
{
  Real iso_cs = pmy_block->pfluid->pf_eos->GetIsoSoundSpeed();

#pragma simd
  for (int i=(pmy_block->is); i<=(pmy_block->ie); ++i) {
    // src_2 = < M_{phi phi} ><cot theta/r>
    Real m_pp = prim(IDN,k,j,i)*SQR(prim(IM3,k,j,i));
    if (NON_BAROTROPIC_EOS) {
       m_pp += prim(IEN,k,j,i);
    } else {
       m_pp += (iso_cs*iso_cs)*prim(IDN,k,j,i);
    }
    if (MAGNETIC_FIELDS_ENABLED) {
       m_pp += 0.5*( SQR(bcc(IB1,k,j,i)) + SQR(bcc(IB2,k,j,i)) - SQR(bcc(IB3,k,j,i)) );
    }
    u(IM2,k,j,i) += dt*coord_src1_j_(j)*m_pp;

    // src_3 = -< M_{phi theta} ><cot theta/r> 
    u(IM3,k,j,i) -= dt*coord_src2_j_(j)*
      (coord_area2_j_(j)*flx(IM3,i) + coord_area2_j_(j+1)*flx_p1(IM3,i));
  }

  return;
}

void Coordinates::CoordSrcTermsX3(const int k, const int j, const Real dt,
  const AthenaArray<Real> &flx,  const AthenaArray<Real> &flx_p1,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u)
{
  return;
}
