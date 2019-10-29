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

  pwave->use_Sommerfeld = pin->GetOrAddBoolean("wave", "use_Sommerfeld",
                                               pwave->use_Sommerfeld);
  sigma = pin->GetOrAddReal("problem", "sigma", sigma);

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

  for(int k = ks; k <= ke; ++k)
    for(int j = js; j <= je; ++j)
      for(int i = is; i <= ie; ++i) {
        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);

        Real t = pmy_mesh->time + pmy_mesh->dt;
        Real c = 1.;

        Real sigma_2 = SQR(sigma);
        Real Gaussian = exp(-SQR(t - x) / (2. * sigma_2));

        pwave->exact(k,j,i) = Gaussian;
        pwave->error(k,j,i) = pwave->u(0,k,j,i) - pwave->exact(k,j,i);

        if (std::abs(pwave->error(k,j,i)) > max_err){
          max_err = std::abs(pwave->error(k,j,i));
          fun_err = pwave->u(0,k,j,i);
        }
      }

  if (pwave->use_Sommerfeld)
    printf("MB::UWIL: (max_err, fun_err)=(%1.7f, %1.7f)\n", max_err, fun_err);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn
//  \brief refinement condition: simple time-dependent test

int RefinementCondition(MeshBlock *pmb){
  // don't do anything
  return 0;
}
