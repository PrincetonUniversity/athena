//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file linear_wave.c
//  \brief Linear wave problem generator for 1D/2D/3D problems.
//
// In 1D, the problem is setup along one of the three coordinate axes (specified by
// setting [ang_2,ang_3] = 0.0 or PI/2 in the input file).  In 2D/3D this routine
// automatically sets the wavevector along the domain diagonal.
//========================================================================================

// C++ headers
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <algorithm>  // min, max

// Athena++ headers
#include "../globals.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../eos/eos.hpp"
#include "../coordinates/coordinates.hpp"
#include "../gravity/gravity.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

#ifdef OPENMP_PARALLEL
#include "omp.h"
#endif

// with functions A1,2,3 which compute vector potentials
static Real ang_2, ang_3; // Rotation angles about the y and z' axis
static Real sin_a2, cos_a2, sin_a3, cos_a3;
static Real amp, njeans, lambda, kwave; // amplitude, Wavelength, 2*PI/wavelength
static Real cs2;
static Real ev[NWAVE], rem[NWAVE][NWAVE], lem[NWAVE][NWAVE];


//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  Real x0=0.0, y0=0.0, z0=0.0;
  Real x, y, z;
  Real dx, dy, dz;
  Real xl, sinkx, coskx;
  RegionSize &mesh_size = pmy_mesh->mesh_size;

  Real x1size = mesh_size.x1max - mesh_size.x1min;
  Real x2size = mesh_size.x2max - mesh_size.x2min;
  Real x3size = mesh_size.x3max - mesh_size.x3min;
 
  amp = pin->GetReal("problem","amp");
  njeans = pin->GetReal("problem","njeans");
  ang_2 = pin->GetOrAddReal("problem","ang_2",-999.9);
  ang_3 = pin->GetOrAddReal("problem","ang_3",-999.9);

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

  Real d0 = 1.0, p0 = 1.0;
  Real u0 = 0.0, v0 = 0.0, w0 = 0.0;
  Real va = 0.0, b0 = 0.0;

  if (NON_BAROTROPIC_EOS) {
    Real gam = pin->GetReal("hydro","gamma");
    cs2 = gam*p0/d0;
  } else {
    Real iso_cs = pin->GetReal("hydro","iso_sound_speed");
    cs2 = SQR(iso_cs);
  }
  Real gconst = cs2*PI*njeans*njeans/(d0*lambda*lambda);
  Real four_pi_G = 4.0*PI*gconst;
  Real grav_mean_rho = d0;

  kwave = 2.0*PI/lambda;
  Real omega2 = SQR(kwave)*cs2*(1.0 - SQR(njeans));
  Real omega = sqrt(fabs(omega2));

  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
  for (int i=is; i<=ie; ++i) {
    Real x = cos_a2*(pcoord->x1v(i)*cos_a3 + pcoord->x2v(j)*sin_a3) + pcoord->x3v(k)*sin_a2;
    sinkx = sin(x*kwave);
    coskx = cos(x*kwave);

    phydro->u(IDN,k,j,i) = d0*(1.0+amp*sinkx);

    Real m = (omega2 < 0) ? d0*(omega/kwave)*amp*coskx:0.0;

    phydro->u(IM1,k,j,i) = m*cos_a3*cos_a2;
    phydro->u(IM2,k,j,i) = m*sin_a3*cos_a2;
    phydro->u(IM3,k,j,i) = m*sin_a2;
  }}}

  std::cout << "four_pi_G " << four_pi_G << std::endl;
  std::cout << "lambda " << lambda << std::endl;
  std::cout << "period " << (2*PI/omega) << std::endl;
  if(SELF_GRAVITY_ENABLED){
    pgrav->gconst = gconst;
    pgrav->four_pi_G = four_pi_G;
    pgrav->grav_mean_rho = grav_mean_rho;
  } // self-gravity
}

void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{
  if (!pin->GetOrAddBoolean("problem","compute_error",false)) return;
  // Initialize errors to zero
  Real l1_err[NHYDRO+NFIELD],max_err[NHYDRO+NFIELD];
  for (int i=0; i<(NHYDRO+NFIELD); ++i) {
    l1_err[i]=0.0;
    max_err[i]=0.0;
  }

  Hydro *phydro = pblock->phydro;
  Coordinates *pcoord = pblock->pcoord;
  Real sinkx, coskx;
  int is=pblock->is, ie=pblock->ie;
  int js=pblock->js, je=pblock->je;
  int ks=pblock->ks, ke=pblock->ke;

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
    Real x = cos_a2*(pcoord->x1v(i)*cos_a3 + pcoord->x2v(j)*sin_a3) + pcoord->x3v(k)*sin_a2;
    sinkx = sin(x*kwave);
    coskx = cos(x*kwave);

    if (omega2 > 0){
      l1_err[IDN] += fabs( d0*(1.0+amp*sinkx) - phydro->u(IDN,k,j,i));
      max_err[IDN] = std::max(fabs(d0*(1.0+amp*sinkx) - phydro->u(IDN,k,j,i)),max_err[IDN]);

      l1_err[IM1] += fabs(phydro->u(IM1,k,j,i));
      l1_err[IM2] += fabs(phydro->u(IM2,k,j,i));
      l1_err[IM3] += fabs(phydro->u(IM3,k,j,i));
      max_err[IM1] = std::max(fabs(phydro->u(IM1,k,j,i)),max_err[IM1]);
      max_err[IM2] = std::max(fabs(phydro->u(IM2,k,j,i)),max_err[IM2]);
      max_err[IM3] = std::max(fabs(phydro->u(IM3,k,j,i)),max_err[IM3]);
    } else {
    }

} 
