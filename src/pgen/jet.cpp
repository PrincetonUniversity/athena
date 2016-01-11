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
void JetInnerX1(MeshBlock *pmb, AthenaArray<Real> &a, FaceField &b,
                   int is, int ie, int js, int je, int ks, int ke);

// Make radius of jet and jet variables global so they can be accessed by BC functions
static Real r_amb,d_amb,p_amb,vx_amb,vy_amb,vz_amb,bx_amb,by_amb,bz_amb;
static Real r_jet,d_jet,p_jet,vx_jet,vy_jet,vz_jet,bx_jet,by_jet,bz_jet;
static Real gm1,x2_0,x3_0;


//--------------------------------------------------------------------------------------
//! \fn void Mesh::InitUserMeshProperties(ParameterInput *pin)
//  \brief Init the Mesh properties
void Mesh::InitUserMeshProperties(ParameterInput *pin)
{
  // initialize global variables
  d_amb  = pin->GetReal("problem", "d");
  p_amb  = pin->GetReal("problem", "p");
  vx_amb = pin->GetReal("problem", "vx");
  vy_amb = pin->GetReal("problem", "vy");
  vz_amb = pin->GetReal("problem", "vz");
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
  x2_0 = 0.5*(mesh_size.x2max + mesh_size.x2min);
  x3_0 = 0.5*(mesh_size.x3max + mesh_size.x3min);
   
  // enroll boundary value function pointers
  EnrollUserBoundaryFunction(INNER_X1, JetInnerX1);
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void Mesh::TerminateUserMeshProperties(void)
//  \brief Clean up the Mesh properties
void Mesh::TerminateUserMeshProperties(void)
{
  // nothing to do
  return;
}



//--------------------------------------------------------------------------------------
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the Jet problem
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  gm1 = phydro->peos->GetGamma() - 1.0;

// initialize conserved variables
   
  for(int k=ks; k<=ke; ++k){
  for(int j=js; j<=je; ++j){
  for(int i=is; i<=ie; ++i){
    phydro->u(IDN,k,j,i) = d_amb;
    phydro->u(IM1,k,j,i) = d_amb*vx_amb;
    phydro->u(IM2,k,j,i) = d_amb*vy_amb;
    phydro->u(IM3,k,j,i) = d_amb*vz_amb;
    if (NON_BAROTROPIC_EOS) {
      phydro->u(IEN,k,j,i) = p_amb/gm1 + 0.5*d_amb*(SQR(vx_amb)+SQR(vy_amb)+SQR(vz_amb));
    }
  }}}

// initialize interface B

  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie+1; ++i) {
      pfield->b.x1f(k,j,i) = bx_amb;
    }}}
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
    for (int i=is; i<=ie; ++i) {
      pfield->b.x2f(k,j,i) = by_amb;
    }}}
    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie; ++i) {
      pfield->b.x3f(k,j,i) = bz_amb;
    }}}
    if (NON_BAROTROPIC_EOS) {
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        phydro->u(IEN,k,j,i) += 0.5*(SQR(bx_amb) + SQR(by_amb) + SQR(bz_amb));
      }}}
    }
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void MeshBlock::UserWorkInLoop(void)
//  \brief User-defined work function for every time step
void MeshBlock::UserWorkInLoop(void)
{
  // nothing to do
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void JetInnerX1()
//  \brief Sets boundary condition on left X boundary (iib) for jet problem

void JetInnerX1(MeshBlock *pmb, AthenaArray<Real> &a, FaceField &b,
                   int is, int ie, int js, int je, int ks, int ke)
{
  // set primitive variables in inlet ghost zones
  for(int k=ks; k<=ke; ++k){
  for(int j=js; j<=je; ++j){
    for(int i=1; i<=(NGHOST); ++i){
      Real rad = sqrt(SQR(pmb->pcoord->x2v(j)-x2_0) + SQR(pmb->pcoord->x3v(k)-x3_0));
      if(rad <= r_jet){
        a(IDN,k,j,is-i) = d_jet;
        a(IVX,k,j,is-i) = vx_jet;
        a(IVY,k,j,is-i) = vy_jet;
        a(IVZ,k,j,is-i) = vz_jet;
        a(IEN,k,j,is-i) = p_jet;
      } else{
        a(IDN,k,j,is-i) = a(IDN,k,j,is);
        a(IVX,k,j,is-i) = a(IVX,k,j,is);
        a(IVY,k,j,is-i) = a(IVY,k,j,is);
        a(IVZ,k,j,is-i) = a(IVZ,k,j,is);
        a(IEN,k,j,is-i) = a(IEN,k,j,is);
      }
    }
  }}

  // set magnetic field in inlet ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je; ++j){
#pragma simd
      for(int i=1; i<=(NGHOST); ++i){
        Real rad = sqrt(SQR(pmb->pcoord->x2v(j)-x2_0) + SQR(pmb->pcoord->x3v(k)-x3_0));
        if(rad <= r_jet){
          b.x1f(k,j,is-i) = bx_jet;
        } else{
          b.x1f(k,j,is-i) = b.x1f(k,j,is);
        }
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        Real rad = sqrt(SQR(pmb->pcoord->x2v(j)-x2_0) + SQR(pmb->pcoord->x3v(k)-x3_0));
        if(rad <= r_jet){
          b.x2f(k,j,is-i) = by_jet;
        } else{
          b.x2f(k,j,is-i) = b.x2f(k,j,is);
        }
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        Real rad = sqrt(SQR(pmb->pcoord->x2v(j)-x2_0) + SQR(pmb->pcoord->x3v(k)-x3_0));
        if(rad <= r_jet){
          b.x3f(k,j,is-i) = bz_jet;
        } else{
          b.x3f(k,j,is-i) = b.x3f(k,j,is);
        }
      }
    }}
  }

}
