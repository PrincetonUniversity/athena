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
 * You should have received a copy of GNU GPL in the file LICENSE included in the code
 * distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

//Primary header
#include "coordinates.hpp"

// C headers
#include <math.h>   // pow function

// C++ headers
#include <stdexcept> // runtime_error
#include <sstream>   // sstream

// Athena++ headers
#include "../athena.hpp"          // macros, Real
#include "../athena_arrays.hpp"   // AthenaArray
#include "../parameter_input.hpp" // ParameterInput
#include "../mesh.hpp"            // MeshBlock
#include "../fluid/fluid.hpp"     // Fluid
#include "../fluid/eos/eos.hpp"   // SoundSpeed()

//======================================================================================
//! \file cylindrical-rphi.cpp
//  \brief implements functions in class Coordinates for 2D (r-phi) cylindrical coords
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

// x3-direction: These coordinates only work in 2D!

  if (pmb->block_size.nx3 == 1) {
    pmb->x3v(ks) = 0.5*(pmb->x3f(ks+1) + pmb->x3f(ks));
    pmb->dx3v(ks) = pmb->dx3f(ks);
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in Coordinates constructor" << std::endl
        << "Cylindrical r-phi coordinates restricted to 2D, but nx3="
        << pmb->pmy_domain->pmy_mesh->mesh_size.nx3 << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

// Allocate memory for scratch arrays used in integrator, and internal scratch arrays
// Allocate only those local scratch arrays needed for cylindrical coordinates

  int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
  face_area.NewAthenaArray(ATHENA_MAX_NUM_THREADS,ncells1);
  cell_volume.NewAthenaArray(ATHENA_MAX_NUM_THREADS,ncells1);
  volume_i_.NewAthenaArray(ncells1);
  src_terms_i_.NewAthenaArray(ncells1);

// compute constant factors used to compute face-areas and cell volumes and store in
// local scratch arrays.  This helps improve performance.

#pragma simd
  for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i){
    volume_i_(i)    = 0.5*(pmb->x1f(i+1)*pmb->x1f(i+1) - pmb->x1f(i)*pmb->x1f(i));
    src_terms_i_(i) = pmb->dx1f(i)/volume_i_(i);
  }

}

// destructor

Coordinates::~Coordinates()
{
  face_area.DeleteAthenaArray();
  cell_volume.DeleteAthenaArray();
  volume_i_.DeleteAthenaArray();
  src_terms_i_.DeleteAthenaArray();
}

//--------------------------------------------------------------------------------------
// \!fn void Coordinates::Area1Face(const int k,const int j, const int il, const int iu,
//        AthenaArray<Real> &area)
// \brief functions to compute area of cell faces in each direction

void Coordinates::Area1Face(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> *parea)
{
// area1 = r dphi 
#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& area_i = (*parea)(i);
    area_i = pmy_block->x1f(i)*pmy_block->dx2f(j);
  }
  return;
}

void Coordinates::Area2Face(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> *parea)
{
// area2 = dr
#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& area_i = (*parea)(i);
    area_i = pmy_block->dx1f(i);
  }
  return;
}

void Coordinates::Area3Face(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> *parea)
{
// 2D only!  area3 = 1.0
#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& area_i = (*parea)(i);
    area_i = 1.0;
  }
  return;
}

//--------------------------------------------------------------------------------------
// \!fn void Coordinates::CellVolume(const int k,const int j,const int il, const int iu,
//        AthenaArray<Real> &vol)
// \brief function to compute cell volume

void Coordinates::CellVolume(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> *pvol)
{
// volume = dr r dphi = d(r^2/2) dphi
#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& vol_i = (*pvol)(i);
    vol_i = volume_i_(i)*(pmy_block->dx2f(j));
  }
  return;
}

//--------------------------------------------------------------------------------------
// \!fn void Coordinates::CoordinateSourceTerms(const int k, const int j,
//        AthenaArray<Real> &prim, AthenaArray<Real> &src)
// \brief function to compute coordinate source term

void Coordinates::CoordinateSourceTerms(Real dt, AthenaArray<Real> &prim,
  AthenaArray<Real> &cons)
{
  Real src[NVAR],dummy_arg[NVAR];

// src_1 = <M_{phi phi}><1/r> = M_{phi phi} dr/d(r^2/2)
// src_2 = -< M_{phi r} ><1/r>  = -(<M_{pr}>) dr/d(r^2/2)

  for (int k=(pmy_block->ks); k<=(pmy_block->ke); ++k) {
  for (int j=(pmy_block->js); j<=(pmy_block->je); ++j) {
#pragma simd
    for (int i=(pmy_block->is); i<=(pmy_block->ie); ++i) {
      Real m_pp = prim(IDN,k,j,i)*prim(IM2,k,j,i)*prim(IM2,k,j,i);
      if (NON_BAROTROPIC_EOS) {
         m_pp += prim(IEN,k,j,i);
      } else {
         Real iso_cs = pmy_block->pfluid->pf_eos->SoundSpeed(dummy_arg);
         m_pp += (iso_cs*iso_cs)*prim(IDN,k,j,i);
      }
      src[IM1] = src_terms_i_(i)*m_pp;

      Real m_pr = prim(IDN,k,j,i)*prim(IM2,k,j,i)*prim(IM1,k,j,i);
      src[IM2] = (-1.0)*src_terms_i_(i)*m_pr;

      Real& uim1 = cons(IM1,k,j,i);
      Real& uim2 = cons(IM2,k,j,i);
      uim1 += dt*src[IM1];
      uim2 += dt*src[IM2];
    }
  }}

  return;
}
