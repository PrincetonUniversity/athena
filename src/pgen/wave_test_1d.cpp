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
#include "../hydro/hydro.hpp"

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

  typedef Real (*unary_function)(Real);

  unary_function prof = NULL;
  unary_function prof_diff = NULL;

  int direction = 0;

} // namespace


//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Sets the initial conditions.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin){
  printf("->wave_test_1d MeshBlock::ProblemGenerator\n");

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

        pwave->u(0,k,j,i) = prof(sin_x);
        pwave->u(1,k,j,i) = -direction*M_PI*c*cos_x*prof_diff(sin_x);

        pwave->exact(0,k,j,i) = pwave->u(0,k,j,i);
        pwave->error(0,k,j,i) = 0.0;
      }

  printf("<-wave_test_1d MeshBlock::ProblemGenerator\n");
  return;
}
void MeshBlock::UserWorkInLoop(){
  printf("->wave_test_1d MeshBlock::UserWorkInLoop\n");
  // test parameter extraction
  printf("c: %f\n", pwave->c);

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

          pwave->exact(0,k,j,i) = prof(xL);
          break;
        case 0:

          xL = sin(M_PI*(x + c*t));
          xR = sin(M_PI*(x - c*t));

          // Average of left/right travelling
          pwave->exact(0,k,j,i) = 0.5*(prof(xL) + prof(xR));
          break;
        case 1:
          // right travelling clean profile
          xR = sin(M_PI*(x - c*t));

          pwave->exact(0,k,j,i) = prof(xR);
          break;
        default:
          assert(false); // you shouldn't be here
          abort();
        }
        pwave->error(0,k,j,i) = pwave->u(0,k,j,i) - pwave->exact(0,k,j,i);

        if (std::abs(pwave->error(0,k,j,i)) > max_err){
          max_err = std::abs(pwave->error(0,k,j,i));
          fun_err = pwave->u(0,k,j,i);
        }
      }

  printf("max_err=%1.7f\n", max_err);
  printf("fun_err=%1.7f\n", fun_err);

  printf("<-wave_test_1d MeshBlock::UserWorkInLoop\n");
  return;
}
