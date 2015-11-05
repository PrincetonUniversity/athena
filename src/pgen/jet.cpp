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
//! \file jet.cpp
//  \brief Sets up a jet introduced through L-x1 boundary (left edge)
//======================================================================================

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../parameter_input.hpp"
#include "../mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../hydro/eos/eos.hpp"
#include "../coordinates/coordinates.hpp"

// BCs on L-x1 (left edge) of grid with jet inflow conditions
void jet_hydro_iib(MeshBlock *pmb, AthenaArray<Real> &a,
                   int is, int ie, int js, int je, int ks, int ke);
void jet_field_iib(MeshBlock *pmb, InterfaceField &a,
                   int is, int ie, int js, int je, int ks, int ke);

// Make radius of jet and jet variables global so they can be accessed by BC functions
static Real r_jet,d_jet,p_jet,vx_jet,vy_jet,vz_jet,bx_jet,by_jet,bz_jet;
static Real gm1,x2_0,x3_0;

void Mesh::ProblemGenerator(Hydro *phyd, Field *pfld, ParameterInput *pin)
{
  MeshBlock *pmb = phyd->pmy_block;

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  gm1 = phyd->pf_eos->GetGamma() - 1.0;

// read parameters from input file

  Real d_amb  = pin->GetReal("problem", "d");
  Real p_amb  = pin->GetReal("problem", "p");
  Real vx_amb = pin->GetReal("problem", "vx");
  Real vy_amb = pin->GetReal("problem", "vy");
  Real vz_amb = pin->GetReal("problem", "vz");
  Real bx_amb,by_amb,bz_amb;
  if (MAGNETIC_FIELDS_ENABLED) {
    bx_amb = pin->GetReal("problem", "bx");
    by_amb = pin->GetReal("problem", "by");
    bz_amb = pin->GetReal("problem", "bz");
  }

  d_jet  = pin->GetReal("problem", "djet");
  p_jet  = pin->GetReal("problem", "pjet");
  vx_jet = pin->GetReal("problem", "vxjet");
  vy_jet = pin->GetReal("problem", "vyjet");
  vz_jet = pin->GetReal("problem", "vzjet");
  if (MAGNETIC_FIELDS_ENABLED) {
    bx_jet = pin->GetReal("problem", "bxjet");
    by_jet = pin->GetReal("problem", "byjet");
    bz_jet = pin->GetReal("problem", "bzjet");
  }
  r_jet = pin->GetReal("problem", "rjet");
   
  x2_0 = 0.5*(pmb->pmy_mesh->mesh_size.x2max + pmb->pmy_mesh->mesh_size.x2min);
  x3_0 = 0.5*(pmb->pmy_mesh->mesh_size.x3max + pmb->pmy_mesh->mesh_size.x3min);

// initialize conserved variables
   
  for(int k=ks; k<=ke; ++k){
  for(int j=js; j<=je; ++j){
  for(int i=is; i<=ie; ++i){
    phyd->u(IDN,k,j,i) = d_amb;
    phyd->u(IM1,k,j,i) = d_amb*vx_amb;
    phyd->u(IM2,k,j,i) = d_amb*vy_amb;
    phyd->u(IM3,k,j,i) = d_amb*vz_amb;
    if (NON_BAROTROPIC_EOS) {
      phyd->u(IEN,k,j,i) = p_amb/gm1 + 0.5*d_amb*(SQR(vx_amb)+SQR(vy_amb)+SQR(vz_amb));
    }
  }}}

// initialize interface B

  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie+1; ++i) {
      pfld->b.x1f(k,j,i) = bx_amb;
    }}}
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
    for (int i=is; i<=ie; ++i) {
      pfld->b.x2f(k,j,i) = by_amb;
    }}}
    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie; ++i) {
      pfld->b.x3f(k,j,i) = bz_amb;
    }}}
    if (NON_BAROTROPIC_EOS) {
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        phyd->u(IEN,k,j,i) += 0.5*(SQR(bx_amb) + SQR(by_amb) + SQR(bz_amb));
      }}}
    }
  }

// Enroll boundary value function pointers
  pmb->pbval->EnrollHydroBoundaryFunction(inner_x1, jet_hydro_iib);

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void jet_hydro_iib()
//  \brief Sets boundary condition for hydro on left X boundary (iib) for jet problem

void jet_hydro_iib(MeshBlock *pmb, AthenaArray<Real> &a,
                   int is, int ie, int js, int je, int ks, int ke)
{
  for(int k=ks; k<=ke; ++k){
  for(int j=js; j<=je; ++j){
    for(int i=1; i<=(NGHOST); ++i){
      Real rad = sqrt(SQR(pmb->pcoord->x2v(j) - x2_0) + SQR(pmb->pcoord->x3v(k) - x3_0));
      if(rad <= r_jet){
        a(IDN,k,j,is-i) = d_jet;
        a(IM1,k,j,is-i) = d_jet*vx_jet;
        a(IM2,k,j,is-i) = d_jet*vy_jet;
        a(IM3,k,j,is-i) = d_jet*vz_jet;
        a(IEN,k,j,is-i) = p_jet/gm1 + 0.5*d_jet*(SQR(vx_jet)+SQR(vy_jet)+SQR(vz_jet));
        if (MAGNETIC_FIELDS_ENABLED){
          a(IEN,k,j,is-i) += 0.5*(SQR(bx_jet)+SQR(by_jet)+SQR(bz_jet));
        }
      } else{
        a(IDN,k,j,is-i) = a(IDN,k,j,is);
        a(IM1,k,j,is-i) = a(IM1,k,j,is);
        a(IM2,k,j,is-i) = a(IM2,k,j,is);
        a(IM3,k,j,is-i) = a(IM3,k,j,is);
        a(IEN,k,j,is-i) = a(IEN,k,j,is);
      }
    }
  }}
}

//--------------------------------------------------------------------------------------
//! \fn void jet_field_iib()
//  \brief Sets boundary condition for B field on left X boundary (iib) for jet problem

void jet_field_iib(MeshBlock *pmb, InterfaceField &a,
                   int is, int ie, int js, int je, int ks, int ke)
{
  for(int k=ks; k<=ke; ++k){
  for(int j=js; j<=je; ++j){
#pragma simd
    for(int i=1; i<=(NGHOST); ++i){
      Real rad = sqrt(SQR(pmb->pcoord->x2v(j) - x2_0) + SQR(pmb->pcoord->x3v(k) - x3_0));
      if(rad <= r_jet){
        a.x1f(k,j,is-i) = bx_jet;
      } else{
        a.x1f(k,j,is-i) = a.x1f(k,j,is);
      }
    }
  }}

  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je+1; ++j) {
#pragma simd
    for (int i=1; i<=(NGHOST); ++i) {
      Real rad = sqrt(SQR(pmb->pcoord->x2v(j) - x2_0) + SQR(pmb->pcoord->x3v(k) - x3_0));
      if(rad <= r_jet){
        a.x2f(k,j,is-i) = by_jet;
      } else{
        a.x2f(k,j,is-i) = a.x2f(k,j,is);
      }
    }
  }}

  for (int k=ks; k<=ke+1; ++k) {
  for (int j=js; j<=je; ++j) {
#pragma simd
    for (int i=1; i<=(NGHOST); ++i) {
      Real rad = sqrt(SQR(pmb->pcoord->x2v(j) - x2_0) + SQR(pmb->pcoord->x3v(k) - x3_0));
      if(rad <= r_jet){
        a.x3f(k,j,is-i) = bz_jet;
      } else{
        a.x3f(k,j,is-i) = a.x3f(k,j,is);
      }
    }
  }}

}
