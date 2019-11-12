//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file wave_test.cpp
//  \brief Initial conditions for the wave equation; rightward propagating Gaussian

#include <cassert> // assert
#include <cmath> // abs, exp, sin, fmod
#include <iostream>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "../mesh/mesh_refinement.hpp"
#include "../wave/wave.hpp"

using namespace std;

int RefinementCondition(MeshBlock *pmb);

// standard deviation of packet
Real sigma = 1. / 16.;
Real phys_x1min = -1.0;
Real phys_x1max = 1.0;
Real amr_sigma_mul = 1;

bool allow_restrict = true;

Real grPhysToFund(Real rho){
  return (rho - phys_x1min) / (phys_x1max - phys_x1min);
}

Real grFundToPhys(Real nu){
  return (phys_x1max - phys_x1min) * nu + phys_x1min;
}

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  if(adaptive==true)
    EnrollUserRefinementCondition(RefinementCondition);
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Initialize the problem.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  pwave->use_Sommerfeld = pin->GetOrAddInteger("wave", "use_Sommerfeld",
                                               pwave->use_Sommerfeld);
  sigma = pin->GetOrAddReal("problem", "sigma", sigma);
  // for amr
  amr_sigma_mul = pin->GetOrAddReal("problem", "amr_sigma_mul", amr_sigma_mul);
  phys_x1min = pin->GetOrAddReal("mesh", "x1min", phys_x1min);
  phys_x1max = pin->GetOrAddReal("mesh", "x1max", phys_x1max);
  allow_restrict = pin->GetOrAddBoolean("problem", "allow_restrict", allow_restrict);
  //-

  for(int k = ks; k <= ke; ++k)
    for(int j = js; j <= je; ++j)
      for(int i = is; i <= ie; ++i) {
        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);


        Real x_2 = SQR(x);
        Real sigma_2 = SQR(sigma);
        Real Gaussian = exp(-x_2 / (2. * sigma_2));

        pwave->u(0,k,j,i) = Gaussian;
        pwave->u(1,k,j,i) = x * Gaussian / sigma_2;

        pwave->exact(k,j,i) = pwave->u(0,k,j,i);
        pwave->error(k,j,i) = 0.0;
      }

  return;
}

void MeshBlock::WaveUserWorkInLoop() {
  Real max_err = 0;
  Real fun_err = 0;

  Real t = pmy_mesh->time + pmy_mesh->dt;

  Real c = 1.0;
  Real sigma_2 = SQR(sigma);

  // for propagation
  Real del = (phys_x1max - phys_x1min);
  Real gc_x0 = phys_x1min + del / 2;
  // fractional distance of interval travelled [from left edge] (wrapped)
  Real gc_xfr = std::fmod((gc_x0 - phys_x1min) + t * c, del) / del;
  // propagated and wrapped physical coordinates
  Real gc_xc = grFundToPhys(gc_xfr);

  for(int k = ks; k <= ke; ++k)
    for(int j = js; j <= je; ++j)
      for(int i = is; i <= ie; ++i) {
        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);

        Real Gaussian = exp(-SQR(x - gc_xc) / (2. * sigma_2));

        pwave->exact(k,j,i) = Gaussian;
        pwave->error(k,j,i) = pwave->u(0,k,j,i) - pwave->exact(k,j,i);

        if (std::abs(pwave->error(k,j,i)) > max_err){
          max_err = std::abs(pwave->error(k,j,i));
          fun_err = pwave->u(0,k,j,i);
        }
      }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn
//  \brief refinement condition: simple time-dependent test

int RefinementCondition(MeshBlock *pmb){

  // physical parameters
  Real c = pmb->pwave->c;
  Real t = pmb->pmy_mesh->time;
  Real del = (phys_x1max - phys_x1min);

  Real r_width = sigma * amr_sigma_mul;

  // Gaussian refinement params
  Real gc_x0 = phys_x1min + del / 2;
  Real gl_x0 = phys_x1min + del / 2 - r_width;
  Real gr_x0 = phys_x1min + del / 2 + r_width;

  // fractional distance of interval travelled [from left edge] (wrapped)
  Real gc_xfr = std::fmod((gc_x0 - phys_x1min) + t * c, del) / del;

  Real gl_xfr = std::fmod((gl_x0 - phys_x1min) + t * c, del) / del;
  Real gr_xfr = std::fmod((gr_x0 - phys_x1min) + t * c, del) / del;

  // propagated and wrapped physical coordinates
  Real gc_xc = grFundToPhys(gc_xfr);
  Real gl_xc = grFundToPhys(gl_xfr);
  Real gr_xc = grFundToPhys(gr_xfr);


  // current block (physical) geometry
  Real x1min = pmb->block_size.x1min;
  Real x1max = pmb->block_size.x1max;

  // if left or right edge within current box then refine
  if ((x1min <= gl_xc) && (gl_xc <= x1max)) {
    return 1;
  }

  if ((x1min <= gc_xc) && (gc_xc <= x1max)) {
    return 1;
  }

  if ((x1min <= gr_xc) && (gr_xc <= x1max)) {
    return 1;
  }

  // otherwise derefine
  if (allow_restrict)
    return -1;
  return 0;
}
