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

namespace {

  Real linear(Real x) {
    return x;
  }

  Real linear_diff(Real x) {
    return 1;
  }

  Real bump(Real x) {
    if(abs(x) < 1.) {
      return exp(-1./(1. - SQR(x)));
    }
    else {
      return 0.;
    }
  }

  Real bump_diff(Real x) {
    if(abs(x) < 1.) {
      return -2.*x*bump(x)/SQR(-1. + SQR(x));
    }
    else {
      return 0.;
    }
  }

  Real unit_Gaussian(Real x, Real sigma) {
    return exp(-SQR(x) / (2. * SQR(sigma))) / (SQRT2 * SQRT_PI * sigma);
  }

  Real unit_Gaussian_diff(Real x, Real sigma) {
    return - x * exp(-SQR(x) / (2. * SQR(sigma))) / (SQRT2 * SQRT_PI * POW3(sigma));
  }

  typedef Real (*unary_function)(Real);

  unary_function prof = NULL;
  unary_function prof_diff = NULL;

  int direction = 0;

} // namespace


int RefinementCondition(MeshBlock *pmb);

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  if(adaptive==true)
    EnrollUserRefinementCondition(RefinementCondition);
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Initialize the problem.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin){
  direction = pin->GetOrAddInteger("problem", "direction", 1);
  if(abs(direction) > 1) {
    cerr << "Invalid direction: " << direction << endl;
    cerr << "Valid values are: -1, 0, and 1" << endl;
    cerr << flush;
    abort();
  }

  string profile = pin->GetOrAddString("problem", "profile", "linear");
  if(profile == "bump") {
    prof = bump;
    prof_diff = bump_diff;
  }
  else {
    prof = linear;
    prof_diff = linear_diff;
  }
  cout<<'\n'<<direction;
  for(int k = ks; k <= ke; ++k)
    for(int j = js; j <= je; ++j)
      for(int i = is; i <= ie; ++i) {
        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);
        Real c = pwave->c;

        //Sinusoidal initial profile
        Real sin_x = sin(M_PI*x);
        Real cos_x = cos(M_PI*x);

        // pwave->u(0,k,j,i) = prof(sin_x);
        // pwave->u(1,k,j,i) = -direction*M_PI*c*cos_x*prof_diff(sin_x);

        pwave->u(0,k,j,i) = prof(x);
        pwave->u(1,k,j,i) = -direction*prof_diff(sin_x);

        pwave->exact(k,j,i) = pwave->u(0,k,j,i);
        pwave->error(k,j,i) = 0.0;
      }

  return;
}

void MeshBlock::UserWorkInLoop(){
  Real max_err = 0;
  Real fun_err = 0;

  for(int k = ks; k <= ke; ++k)
      for(int j = js; j <= je; ++j)
      for(int i = is; i <= ie; ++i) {
        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);
        Real t = pmy_mesh->time + pmy_mesh->dt;
        Real c = pwave->c;

        Real xL, xR;

        switch(direction) {
        case -1:
          // left travelling clean profile
          xL = sin(M_PI * (x + c * t));

          pwave->exact(k,j,i) = prof(xL);
          break;
        case 0:

          // xL = sin(M_PI*(x + c*t));
          // xR = sin(M_PI*(x - c*t));

          xL = x + c*t;
          xR = x - c*t;

          // Average of left/right travelling
          pwave->exact(k,j,i) = 0.5*(prof(xL) + prof(xR));
          break;
        case 1:
          // right travelling clean profile
          xR = sin(M_PI*(x - c*t));

          pwave->exact(k,j,i) = prof(xR);
          break;
        default:
          assert(false); // you shouldn't be here
          abort();
        }
        pwave->error(k,j,i) = pwave->u(0,k,j,i) - pwave->exact(k,j,i);

        if (std::abs(pwave->error(k,j,i)) > max_err){
          max_err = std::abs(pwave->error(k,j,i));
          fun_err = pwave->u(0,k,j,i);
        }
      }

  printf("MB::UWIL: (max_err, fun_err)=(%1.7f, %1.7f)\n", max_err, fun_err);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn
//  \brief refinement condition: simple time-dependent test

int RefinementCondition(MeshBlock *pmb){
  // cout<<"wave_amr_2d RefinementCondition\n";

  // current time and maximum time
  Real t = pmb->pmy_mesh->time;
  Real tmax = pmb->pmy_mesh->tlim;

  // As a test we refine mesh globally between two times
  // if(t > 0.1){
  //   if(t <= 0.3)
  //     return 1;
  // }

  // if(t<=0.5)
  //   return -1;

  // current block (physical) geometry
  Real x1min = pmb->block_size.x1min;
  // Real x2min = pmb->block_size.x2min;
  Real x1max = pmb->block_size.x1max;
  // Real x2max = pmb->block_size.x2max;

  // box geometry to refine on
  Real br_x1min = -1.0;
  Real br_x1max = 1.0;

  // Real br_x2min = 0.2;
  // Real br_x2max = 0.3;

  // Slide refinement box in positive direction as time increases
  Real tiniref = 0.1; // time we start refining
  Real delta_x1 = 0.8; // distance to slide by
  // Real delta_x2 = -0.8;

  // velocity of box
  Real vx1 = delta_x1 / (tmax - tiniref);
  // Real vx2 = delta_x2 / (tmax - tiniref);

  br_x1min = br_x1min + vx1 * (t - tiniref);
  br_x1max = br_x1max + vx1 * (t - tiniref);

  // br_x2min = br_x2min + vx2 * (t - tiniref);
  // br_x2max = br_x2max + vx2 * (t - tiniref);

  // refine a box if t\in[0.1,0.3] otherwise derefine
  if(t > 0.0){
    if(t <= 0.9){
      // check x1 interval contains refinement box
      if(((x1min < br_x1min) && (br_x1min < x1max))
         || ((br_x1max < x1max) && (x1min < br_x1max))){
        cout << "Refinement Proceeding...\n";
        return 1;

        // // now for x2
        // if(((x2min < br_x2min) && (br_x2min < x2max))
        //    || ((br_x2max < x2max) && (x2min < br_x2max))){
        //   cout << "Refinement Proceeding...\n";
        //   return 1;
        // }
      }
      // derefine outside of box
      return -1;
    } else {
      // at later times simply derefine globally
      return -1;
    }
  }

  // otherwise don't do anything
  return 0;
}
