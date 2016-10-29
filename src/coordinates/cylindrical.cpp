//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cylindrical.cpp
//  \brief implements Coordinates class functions for cylindrical (r-phi-z) coordinates

// C/C++ headers
#include <math.h>   // pow function

// Athena++ headers
#include "coordinates.hpp"
#include "../athena.hpp"           // macros, Real
#include "../athena_arrays.hpp"    // AthenaArray
#include "../parameter_input.hpp"  // ParameterInput
#include "../mesh/mesh.hpp"        // MeshBlock
#include "../eos/eos.hpp"          // GetIsoSoundSpeed()

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
  // x1-direction: x1v = (\int r dV / \int dV) = d(r^3/3)d(r^2/2)
  for (int i=is-ng; i<=ie+ng; ++i) {
    x1v(i) = (2.0/3.0)*(pow(x1f(i+1),3) - pow(x1f(i),3))
                      /(pow(x1f(i+1),2) - pow(x1f(i),2));
  }
  for (int i=is-ng; i<=ie+ng-1; ++i) {
    dx1v(i) = x1v(i+1) - x1v(i);
  }

  // x2-direction: x2v = (\int phi dV / \int dV) = dphi/2
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

  // x3-direction: x3v = (\int z dV / \int dV) = dz/2
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

  // Allocate memory for scratch arrays used in integrator, and internal scratch arrays
  if(cflag==0) {
    int ncells1 = pmb->block_size.nx1 + 2*ng;
    coord_area3_i_.NewAthenaArray(ncells1);
    coord_vol_i_.NewAthenaArray(ncells1);
    coord_src1_i_.NewAthenaArray(ncells1);
    coord_src2_i_.NewAthenaArray(ncells1);
    phy_src1_i_.NewAthenaArray(ncells1);
    phy_src2_i_.NewAthenaArray(ncells1);

    // Compute and store constant coefficients needed for face-areas, cell-volumes, etc.
    // This helps improve performance.
#pragma simd
    for (int i=is-ng; i<=ie+ng; ++i){
      Real rm = x1f(i  );
      Real rp = x1f(i+1);
      // dV = 0.5*(R_{i+1}^2 - R_{i}^2)
      coord_area3_i_(i)= 0.5*(rp*rp - rm*rm);
      // dV = 0.5*(R_{i+1}^2 - R_{i}^2)
      coord_vol_i_(i) = coord_area3_i_(i);
      // (A1^{+} - A1^{-})/dV
      coord_src1_i_(i) = dx1f(i)/coord_vol_i_(i);
      // (dR/2)/(R_c dV)
      coord_src2_i_(i) = dx1f(i)/((rm + rp)*coord_vol_i_(i));
      // Rf_{i}/R_{i}/Rf_{i}^2
      phy_src1_i_(i) = 1.0/(x1v(i)*x1f(i));
    }
#pragma simd
    for (int i=is-ng; i<=ie+(ng-1); ++i){
       // Rf_{i+1}/R_{i}/Rf_{i+1}^2 
      phy_src2_i_(i) = 1.0/(x1v(i)*x1f(i+1));
    }
  }
}

// destructor

Coordinates::~Coordinates()
{
  DeleteBasicCoordinates();

  if(cflag==0) {
    coord_area3_i_.DeleteAthenaArray();
    coord_vol_i_.DeleteAthenaArray();
    coord_src1_i_.DeleteAthenaArray();
    coord_src2_i_.DeleteAthenaArray();
    phy_src1_i_.DeleteAthenaArray();
    phy_src2_i_.DeleteAthenaArray();
  }
}

//----------------------------------------------------------------------------------------
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
    len(i) = x1f(i)*dx2f(j);
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
  return x1f(i)*dx2f(j);
}

Real Coordinates::GetEdge3Length(const int k, const int j, const int i)
{
  return dx3f(k);
}


//----------------------------------------------------------------------------------------
// Cell-center Width functions: returns physical width at cell-center

Real Coordinates::CenterWidth1(const int k, const int j, const int i)
{
  return dx1f(i);
}

Real Coordinates::CenterWidth2(const int k, const int j, const int i)
{
  return x1v(i)*dx2f(j);
}

Real Coordinates::CenterWidth3(const int k, const int j, const int i)
{
  return dx3f(k);
}


//----------------------------------------------------------------------------------------
// Face Area functions

void Coordinates::Face1Area(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    // area1 = r dphi dz 
    area(i) = x1f(i)*dx2f(j)*dx3f(k);
  }
  return;
}

void Coordinates::Face2Area(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    // area2 = dr dz
    area(i) = dx1f(i)*dx3f(k);
  }
  return;
}

void Coordinates::Face3Area(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    // area3 = dr r dphi = d(r^2/2) dphi
    area(i) = coord_area3_i_(i)*dx2f(j);
  }
  return;
}

// GetFace1Area returns only one Face1Area at i
Real Coordinates::GetFace1Area(const int k, const int j, const int i)
{
  return x1f(i)*dx2f(j)*dx3f(k);
}

//----------------------------------------------------------------------------------------
// Cell Volume function

void Coordinates::CellVolume(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &vol)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    // volume = dr dz r dphi = d(r^2/2) dphi dz
    vol(i) = coord_vol_i_(i)*dx2f(j)*dx3f(k);
  }
  return;
}

// GetCellVolume returns only one CellVolume at i
Real Coordinates::GetCellVolume(const int k, const int j, const int i)
{
  return coord_vol_i_(i)*dx2f(j)*dx3f(k);
}

//----------------------------------------------------------------------------------------
// Coordinate (Geometric) source term functions

void Coordinates::CoordSrcTerms(const Real dt, const AthenaArray<Real> *flux,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u)
{
  Real iso_cs = pmy_block->peos->GetIsoSoundSpeed();

  for (int k=pmy_block->ks; k<=pmy_block->ke; ++k) {
#pragma omp parallel for schedule(static)
    for (int j=pmy_block->js; j<=pmy_block->je; ++j) {
#pragma simd
      for (int i=pmy_block->is; i<=pmy_block->ie; ++i) {
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
        Real& x_i   = x1f(i);
        Real& x_ip1 = x1f(i+1);
        u(IM2,k,j,i) -= dt*coord_src2_i_(i)*(x_i*flux[X1DIR](IM2,k,j,i)
                                           + x_ip1*flux[X1DIR](IM2,k,j,i+1));
      }
    }
  }

  return;
}

//---------------------------------------------------------------------------------------
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
        area_p1 = x1f(i+1)*dx2f(j)*dx3f(k); 
        area = x1f(i)*dx2f(j)*dx3f(k);
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
          area_p1 = coord_area3_i_(i)*(dx2f(j)); 
          area = area_p1;
          vel_p1 = 0.5*(prim(IM3,k+1,j,i)+prim(IM3,k,j,i));
          vel = 0.5*(prim(IM3,k,j,i)+prim(IM3,k-1,j,i));
          divv(k,j,i) += area_p1*vel_p1 - area*vel;
        }
      }
      for (int i=il; i<=iu; ++i){
        vol = coord_vol_i_(i)*dx2f(j)*dx3f(k); 
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
    len(i) = x1f(i)*(prim(IM2,k,j,i)/x1v(i)
                    -prim(IM2,k,j,i-1)/x1v(i-1))/dx1v(i-1)
             +0.5*(prim(IM1,k,j+1,i)+prim(IM1,k,j+1,i-1)-prim(IM1,k,j-1,i)-prim(IM1,k,j-1,i-1))
              /x1f(i)/(dx2v(j-1)+dx2v(j));
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
    len(i) = (prim(IM1,k,j,i)-prim(IM1,k,j-1,i))/x1v(i)/dx2v(j-1)
             +x1v(i)*0.5*((prim(IM2,k,j,i+1)+prim(IM2,k,j-1,i+1))/x1v(i+1)
                         -(prim(IM2,k,j,i-1)+prim(IM2,k,j-1,i-1))/x1v(i-1))
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
    len(i) = (prim(IM2,k,j,i)-prim(IM2,k,j-1,i))/x1v(i)/dx2v(j-1) 
              +0.5*(prim(IM1,k,j,i)+prim(IM1,k,j-1,i))/x1v(i);
  }
  return;
}

// v_{x3;x2}+v_{x2;x3}  covariant derivative at x2 interface
void Coordinates::FaceYdz(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = (prim(IM3,k,j,i)-prim(IM3,k,j-1,i))/x1v(i)/dx2v(j-1)
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
             /x1v(i)/(dx2v(j-1)+dx2v(j));
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

//----------------------------------------------------------------------------------------
// Viscous (Geometric) source term functions

void Coordinates::VisSrcTermsX1(const int k, const int j, const Real dt,
  const AthenaArray<Real> &flx,
  const AthenaArray<Real> &prim, AthenaArray<Real> &u)
{
  Real iso_cs = pmy_block->peos->GetIsoSoundSpeed();

#pragma simd
  for (int i=(pmy_block->is); i<=(pmy_block->ie); ++i) {
    // src_2 = < M_{phi r} ><1/r>
    Real& x_i   = x1f(i);
    Real& x_ip1 = x1f(i+1);
    u(IM2,k,j,i) += dt*coord_src2_i_(i)*(x_i*flx(IM2,i) + x_ip1*flx(IM2,i+1));
  }

  return;
}

void Coordinates::VisSrcTermsX2(const int k, const int j, const Real dt,
  const AthenaArray<Real> &flx,  const AthenaArray<Real> &flx_p1,
  const AthenaArray<Real> &prim, AthenaArray<Real> &u)
{
  #pragma simd
  for (int i=(pmy_block->is); i<=(pmy_block->ie); ++i) {
    // src_1 = -<M_{phi phi}><1/r>
    u(IM1,k,j,i) -= dt*coord_src1_i_(i)*0.5*(flx(IM2,i) + flx_p1(IM2,i));
  }

  return;
}

void Coordinates::VisSrcTermsX3(const int k, const int j, const Real dt,
  const AthenaArray<Real> &flx,  const AthenaArray<Real> &flx_p1,
  const AthenaArray<Real> &prim, AthenaArray<Real> &u)
{
  return;
}

