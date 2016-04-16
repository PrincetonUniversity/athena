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
//! \file cpaw.c
//  \brief Circularly polarized Alfven wave (CPAW) for 1D/2D/3D problems
//
// In 1D, the problem is setup along one of the three coordinate axes (specified by
// setting [ang_2,ang_3] = 0.0 or PI/2 in the input file).  In 2D/3D this routine
// automatically sets the wavevector along the domain diagonal.
//
// Can be used for [standing/traveling] waves [(problem/v_par=1.0)/(problem/v_par=0.0)]
//
// REFERENCE: G. Toth,  "The div(B)=0 constraint in shock capturing MHD codes", JCP,
//   161, 605 (2000)
//======================================================================================

// C/C++ headers
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../hydro/eos/eos.hpp"
#include "../coordinates/coordinates.hpp"

#if !MAGNETIC_FIELDS_ENABLED
#error "This problem generator requires magnetic fields"
#endif

// Parameters which define initial solution -- made global so that they can be shared
// with functions A1,2,3 which compute vector potentials
static Real b_par, b_perp, den, v_perp, v_par;
static Real ang_2, ang_3; // Rotation angles about the y and z' axis
static Real fac, sin_a2, cos_a2, sin_a3, cos_a3;
static Real lambda, k_par; // Wavelength, 2*PI/wavelength

// functions to compute vector potential to initialize the solution
static Real A1(const Real x1, const Real x2, const Real x3);
static Real A2(const Real x1, const Real x2, const Real x3);
static Real A3(const Real x1, const Real x2, const Real x3);


//======================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//======================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin)
{
// Initialize magnetic field parameters
// For wavevector along coordinate axes, set desired values of ang_2/ang_3.
//    For example, for 1D problem use ang_2 = ang_3 = 0.0
//    For wavevector along grid diagonal, do not input values for ang_2/ang_3.
// Code below will automatically calculate these imposing periodicity and exactly one
// wavelength along each grid direction
  b_par = pin->GetReal("problem","b_par");
  b_perp = pin->GetReal("problem","b_perp");
  v_par = pin->GetReal("problem","v_par");
  ang_2 = pin->GetOrAddReal("problem","ang_2",-999.9);
  ang_3 = pin->GetOrAddReal("problem","ang_3",-999.9);
  Real dir = pin->GetOrAddReal("problem","dir",1); // right(1)/left(2) polarization
  
  Real x1size = mesh_size.x1max - mesh_size.x1min;
  Real x2size = mesh_size.x2max - mesh_size.x2min;
  Real x3size = mesh_size.x3max - mesh_size.x3min;

// User should never input -999.9 in angles
  if (ang_3 == -999.9) ang_3 = atan(x1size/x2size);
  sin_a3 = sin(ang_3);
  cos_a3 = cos(ang_3);
  
  if (ang_2 == -999.9) ang_2 = atan(0.5*(x1size*cos_a3 + x2size*sin_a3)/x3size);
  sin_a2 = sin(ang_2); 
  cos_a2 = cos(ang_2);
  
  Real x1 = x1size*cos_a2*cos_a3;
  Real x2 = x2size*cos_a2*sin_a3;
  Real x3 = x3size*sin_a2;

// For lambda choose the smaller of the 3
  lambda = x1;
  if (mesh_size.nx2 > 1 && ang_3 != 0.0) lambda = std::min(lambda,x2);
  if (mesh_size.nx3 > 1 && ang_2 != 0.0) lambda = std::min(lambda,x3);

// Initialize k_parallel
  k_par = 2.0*(PI)/lambda;
  v_perp = b_perp/sqrt(den);

  if (dir == 1) // right polarization
    fac = 1.0;
  else          // left polarization
    fac = -1.0;
  return;
}

//======================================================================================
//! \fn ProblemGenerator
//  \brief circularly polarized Alfven wave problem generator for 1D/2D/3D problems.
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  Real gm1 = (phydro->peos->GetGamma() - 1.0);
  Real pres = pin->GetReal("problem","pres");
  den = 1.0;

// Use the vector potential to initialize the interface magnetic fields

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie+1; i++) {
      pfield->b.x1f(k,j,i) =
        (A3(pcoord->x1f(i),pcoord->x2f(j+1),pcoord->x3v(k  )) -
         A3(pcoord->x1f(i),pcoord->x2f(j  ),pcoord->x3v(k  )))/pcoord->dx2f(j) -
        (A2(pcoord->x1f(i),pcoord->x2v(j  ),pcoord->x3f(k+1)) -
         A2(pcoord->x1f(i),pcoord->x2v(j  ),pcoord->x3f(k  )))/pcoord->dx3f(k);
    }
  }}

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je+1; j++) {
    for (int i=is; i<=ie; i++) {
      pfield->b.x2f(k,j,i) =
        (A1(pcoord->x1v(i  ),pcoord->x2f(j),pcoord->x3f(k+1)) -
         A1(pcoord->x1v(i  ),pcoord->x2f(j),pcoord->x3f(k  )))/pcoord->dx3f(k) -
        (A3(pcoord->x1f(i+1),pcoord->x2f(j),pcoord->x3v(k  )) -
         A3(pcoord->x1f(i  ),pcoord->x2f(j),pcoord->x3v(k  )))/pcoord->dx1f(i);
    }
  }}

  for (int k=ks; k<=ke+1; k++) {
  for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      pfield->b.x3f(k,j,i) =
        (A2(pcoord->x1f(i+1),pcoord->x2v(j  ),pcoord->x3f(k)) -
         A2(pcoord->x1f(i  ),pcoord->x2v(j  ),pcoord->x3f(k)))/pcoord->dx1f(i) -
        (A1(pcoord->x1v(i  ),pcoord->x2f(j+1),pcoord->x3f(k)) -
         A1(pcoord->x1v(i  ),pcoord->x2f(j  ),pcoord->x3f(k)))/pcoord->dx2f(j);
    }
  }}

/* Now initialize rest of the cell centered quantities */

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      Real x = cos_a2*(pcoord->x1v(i)*cos_a3 + pcoord->x2v(j)*sin_a3) + pcoord->x3v(k)*sin_a2;
      Real sn = sin(k_par*x);
      Real cs = fac*cos(k_par*x);

      phydro->u(IDN,k,j,i) = den;

      Real mx = den*v_par;
      Real my = -fac*den*v_perp*sn;
      Real mz = -den*v_perp*cs;

      phydro->u(IM1,k,j,i) = mx*cos_a2*cos_a3 - my*sin_a3 - mz*sin_a2*cos_a3;
      phydro->u(IM2,k,j,i) = mx*cos_a2*sin_a3 + my*cos_a3 - mz*sin_a2*sin_a3;
      phydro->u(IM3,k,j,i) = mx*sin_a2                    + mz*cos_a2;

      if (NON_BAROTROPIC_EOS) {
        phydro->u(IEN,k,j,i) = pres/gm1 +
          0.5*(SQR(0.5*(pfield->b.x1f(k,j,i) + pfield->b.x1f(k,j,i+1))) +
               SQR(0.5*(pfield->b.x2f(k,j,i) + pfield->b.x2f(k,j+1,i))) +
               SQR(0.5*(pfield->b.x3f(k,j,i) + pfield->b.x3f(k+1,j,i)))) + (0.5/den)*
          (SQR(phydro->u(IM1,k,j,i)) + SQR(phydro->u(IM2,k,j,i)) + SQR(phydro->u(IM3,k,j,i)));
      }
    }
  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn static Real A1(const Real x1,const Real x2,const Real x3)
//  \brief A1: 1-component of vector potential, using a gauge such that Ax = 0, and Ay,
//  Az are functions of x and y alone.

static Real A1(const Real x1, const Real x2, const Real x3)
{
  Real x =  x1*cos_a2*cos_a3 + x2*cos_a2*sin_a3 + x3*sin_a2;
  Real y = -x1*sin_a3        + x2*cos_a3;
  Real Ay = fac*(b_perp/k_par)*sin(k_par*(x));
  Real Az = (b_perp/k_par)*cos(k_par*(x)) + b_par*y;

  return -Ay*sin_a3 - Az*sin_a2*cos_a3;
}

//--------------------------------------------------------------------------------------
//! \fn static Real A2(const Real x1,const Real x2,const Real x3)
//  \brief A2: 2-component of vector potential

static Real A2(const Real x1, const Real x2, const Real x3)
{
  Real x =  x1*cos_a2*cos_a3 + x2*cos_a2*sin_a3 + x3*sin_a2;
  Real y = -x1*sin_a3        + x2*cos_a3;
  Real Ay = fac*(b_perp/k_par)*sin(k_par*(x));
  Real Az = (b_perp/k_par)*cos(k_par*(x)) + b_par*y;

  return Ay*cos_a3 - Az*sin_a2*sin_a3;
}

//--------------------------------------------------------------------------------------
//! \fn static Real A3(const Real x1,const Real x2,const Real x3)
//  \brief A3: 3-component of vector potential

static Real A3(const Real x1, const Real x2, const Real x3)
{
  Real x =  x1*cos_a2*cos_a3 + x2*cos_a2*sin_a3 + x3*sin_a2;
  Real y = -x1*sin_a3        + x2*cos_a3;
  Real Az = (b_perp/k_par)*cos(k_par*(x)) + b_par*y;

  return Az*cos_a2;
}
