//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file resist.cpp
//  \brief Problem generator for resistivy diffusion of B-field.
//========================================================================================

// C headers

// C++ headers
#include <cmath>      // sqrt()
#include <cstring>    // strcmp()
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

#if !MAGNETIC_FIELDS_ENABLED
#error "This problem generator requires magnetic fields"
#endif

namespace {
Real x0, t0;
} // namespace

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  //Real gm1 = peos->GetGamma() - 1.0;

  // Read initial conditions, diffusion coefficients (if needed)
  Real etaO = pin->GetOrAddReal("problem","eta_ohm",0.03);
  Real amp = pin->GetOrAddReal("problem","amp",1e-3);
  t0 = pin->GetOrAddReal("problem","t0",0.5);
  x0 = pin->GetOrAddReal("problem","x0",0.0);
  int iprob = pin->GetOrAddInteger("problem","iprob",0);

  Real eta = etaO;   // set to 0.03 for debug only

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        phydro->u(IDN,k,j,i) = 1.0;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
      }
    }
  }

  // initialize interface B
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    if (iprob == 0) { // initialize B along y-axis
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie+1; i++) {
            pfield->b.x1f(k,j,i) = 0.0;
          }
        }
      }
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
          for (int i=is; i<=ie; i++) {
            Real x1 = pcoord->x1v(i);
            pfield->b.x2f(k,j,i) = amp/std::sqrt(4.0*PI*eta*t0)
                                   *std::exp(-SQR(x1-x0)/(4.0*eta*t0));
          }
        }
      }
      for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
    }

    if (iprob == 1) { // initialize B diagonally
      Real Lx = 2.0;
      Real Ly = 1.0;
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie+1; i++) {
            Real cost = Ly/std::sqrt(SQR(Lx)+SQR(Ly));
            Real sint = Lx/std::sqrt(SQR(Lx)+SQR(Ly));
            Real x1 = pcoord->x1f(i)*cost+pcoord->x2f(j)*sint;
            Real bprim = amp/std::sqrt(4.0*PI*eta*t0)*std::exp(-SQR(x1-x0)/(4.0*eta*t0));
            pfield->b.x1f(k,j,i) = -bprim*sint;
          }
        }
      }
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
          for (int i=is; i<=ie; i++) {
            Real cost = Ly/std::sqrt(SQR(Lx)+SQR(Ly));
            Real sint = Lx/std::sqrt(SQR(Lx)+SQR(Ly));
            Real x1 = pcoord->x1f(i)*cost+pcoord->x2f(j)*sint;
            Real bprim = amp/std::sqrt(4.0*PI*eta*t0)*std::exp(-SQR(x1-x0)/(4.0*eta*t0));
            pfield->b.x2f(k,j,i) = bprim*cost;
          }
        }
      }
      for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
    }
  }
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    if (x0 != 1.0) x0 = 1.0;
    if (iprob == 0) { // initialize B along y-axis
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie+1; i++) {
            Real rad=pcoord->x1v(i);
            Real phi=pcoord->x2v(j);
            Real z=pcoord->x3v(k);
            Real x1=rad*std::cos(phi);
            Real x2=rad*std::sin(phi);
            Real x3=z;
            Real bprim = amp/std::sqrt(4.0*PI*eta*t0)*std::exp(-SQR(x1-x0)/(4.0*eta*t0));
            pfield->b.x1f(k,j,i) = bprim*std::sin(phi);
          }
        }
      }
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
          for (int i=is; i<=ie; i++) {
            Real rad=pcoord->x1v(i);
            Real phi=pcoord->x2v(j);
            Real z=pcoord->x3v(k);
            Real x1=rad*std::cos(phi);
            Real x2=rad*std::sin(phi);
            Real x3=z;
            Real bprim = amp/std::sqrt(4.0*PI*eta*t0)*std::exp(-SQR(x1-x0)/(4.0*eta*t0));
            pfield->b.x2f(k,j,i) = bprim*std::cos(phi);
          }
        }
      }
      for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
    } else {
      std::cout << "only iprob = 0 allowed for cylindrical coord" << std::endl;
    }
  }
  return;
}
