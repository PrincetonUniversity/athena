//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file wave_test.cpp
//  \brief Initial conditions for the wave equation

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

  for(int k = ks; k <= ke; ++k)
    for(int j = js; j <= je; ++j)
#pragma omp simd
      for(int i = is; i <= ie; ++i) {
        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);
        Real c = pwave->c;

        Real cos_x = cos(PI * x);
        Real cos_y = cos(PI * y);
        Real cos_z = cos(PI * z);

        pwave->u(0,k,j,i) = cos_x * cos_y * cos_z;
        pwave->u(1,k,j,i) = 0.;

        pwave->exact(k,j,i) = pwave->u(0,k,j,i);
        pwave->error(k,j,i) = 0.0;
      }

  return;
}

void MeshBlock::WaveUserWorkInLoop() {
  Real max_err = 0;
  Real fun_err = 0;

  for(int k = ks; k <= ke; ++k){
    for(int j = js; j <= je; ++j){
#pragma omp simd
      for(int i = is; i <= ie; ++i) {
        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);
        Real t = pmy_mesh->time + pmy_mesh->dt;
        Real c = pwave->c;

        Real cos_s3t = cos(sqrt(3) * c * PI * t);
        Real cos_x = cos(PI * x);
        Real cos_y = cos(PI * y);
        Real cos_z = cos(PI * z);


        pwave->exact(k,j,i) = cos_x * cos_y * cos_z * cos_s3t;
        pwave->error(k,j,i) = pwave->u(0,k,j,i) - pwave->exact(k,j,i);

        if (std::abs(pwave->error(k,j,i)) > max_err){
          max_err = std::abs(pwave->error(k,j,i));
          fun_err = pwave->u(0,k,j,i);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn
//  \brief refinement condition: simple time-dependent test

int RefinementCondition(MeshBlock *pmb){
  // don't do anything
  return 0;
}
