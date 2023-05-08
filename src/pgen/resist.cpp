//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file resist.cpp
//! \brief Problem generator for resistive diffusion of B-field.
//!
//! - iprob = 0 - Resistive Diffusion of 1-D Gaussian
//! - iprob = 1 - Resistive Diffusion of 2-D Gaussian
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
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

#if !MAGNETIC_FIELDS_ENABLED
#error "This problem generator requires magnetic fields"
#endif

Real threshold;

int RefinementCondition(MeshBlock *pmb);

void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (adaptive) {
    EnrollUserRefinementCondition(RefinementCondition);
    threshold = pin->GetReal("problem","thr");
  }
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real x1, x2;
  Real r, phi;
  Real v1 = 0., v2 = 0., v3 = 0.;
  Real d0 = 1.;
  int iprob = pin->GetOrAddInteger("problem", "iprob", 0);
  Real amp = pin->GetOrAddReal("problem", "amp", 1.e-6);
  Real t0 = pin->GetOrAddReal("problem", "t0", 0.5);
  Real eta_ohm = pin->GetOrAddReal("problem", "eta_ohm", 0.25);
  Real x10 = pin->GetOrAddReal("problem", "x10", 0.);
  Real x20 = pin->GetOrAddReal("problem", "x20", 0.);

  // set hydro variables
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        phydro->u(IDN,k,j,i) = d0;
        phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i)*v1;
        phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i)*v2;
        phydro->u(IM3,k,j,i) = phydro->u(IDN,k,j,i)*v3;
      }
    }
  }

  // set face-centered B
  if (iprob == 0) {
    if (std::strcmp(COORDINATE_SYSTEM, "cartesian") != 0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in resist.cpp ProblemGenerator" << std::endl
          << "1-D diffusion test only compatible with cartesian coord" << std::endl;
      ATHENA_ERROR(msg);
    }
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie+1; i++) {
          pfield->b.x1f(k,j,i) = 0.;
        }
      }
    }
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie; i++) {
          Real x1 = pcoord->x1v(i);
          pfield->b.x2f(k,j,i) = amp/std::pow(std::sqrt(4.*PI*eta_ohm*t0),1.0)
                                  * std::exp(-(std::pow(x1-x10,2.))
                                             / (4.*eta_ohm*t0));
        }
      }
    }
    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          pfield->b.x3f(k,j,i) = 0.;
        }
      }
    }
  } else if (iprob == 1) { // 2-d diffusion
    if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie+1; i++) {
            pfield->b.x1f(k,j,i) = 0.;
          }
        }
      }
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x2f(k,j,i) = 0.;
          }
        }
      }
      for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            x1 = pcoord->x1v(i);
            x2 = pcoord->x2v(j);
            pfield->b.x3f(k,j,i) = amp/std::pow(std::sqrt(4.*PI*eta_ohm*t0),2.0)
                                   * std::exp(-(std::pow(x1-x10,2.)
                                                + std::pow(x2-x20,2.))
                                              / (4.*eta_ohm*t0));
          }
        }
      }
    } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie+1; i++) {
            pfield->b.x1f(k,j,i) = 0.;
          }
        }
      }
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x2f(k,j,i) = 0.;
          }
        }
      }
      for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            r = pcoord->x1v(i);
            phi = pcoord->x2v(j);
            x1 = r*std::cos(phi);
            x2 = r*std::sin(phi);
            pfield->b.x3f(k,j,i) = amp/std::pow(std::sqrt(4.*PI*eta_ohm*t0),2.0)
                                   * std::exp(-(std::pow(x1-x10,2.)
                                                + std::pow(x2-x20,2.))
                                              / (4.*eta_ohm*t0));
          }
        }
      }
    } else {
      std::stringstream msg;
      msg << "### FATAL ERROR in resist.cpp ProblemGenerator" << std::endl
          << "2-d diffusion test only compatible with cartesian or"
          << std::endl << "cylindrical coord" << std::endl;
      ATHENA_ERROR(msg);
    }
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in resist.cpp ProblemGenerator" << std::endl
        << "iprob must be set to 0 or 1" << std::endl;
    ATHENA_ERROR(msg);
  }

  return;
}

// refinement condition: check the maximum b-velocity gradient
int RefinementCondition(MeshBlock *pmb) {
  int f2 = pmb->pmy_mesh->f2, f3 = pmb->pmy_mesh->f3;
  AthenaArray<Real> &bcc = pmb->pfield->bcc;
  Real maxeps = 0.;
  if (f2) {
    int k = pmb->ks;
    for (int j=pmb->js-1; j<=pmb->je+1; j++) {
      for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
        Real eps = std::sqrt(SQR(0.5*(bcc(IB3,k,j,i+1) - bcc(IB3,k,j,i-1)))
                             + SQR(0.5*(bcc(IB3,k,j+1,i) - bcc(IB3,k,j-1,i))));
        maxeps = std::max(maxeps, eps);
      }
    }
  } else {
    int k = pmb->ks;
    int j = pmb->js;
    for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
      Real eps = std::sqrt(SQR(0.5*(bcc(IB2,k,j,i+1) - bcc(IB2,k,j,i-1))));
      maxeps = std::max(maxeps, eps);
    }
  }

  if (maxeps > threshold) return 1;
  if (maxeps < 0.25*threshold) return -1;
  return 0;
}
