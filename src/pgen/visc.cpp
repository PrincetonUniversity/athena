//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file visc.cpp
//! iprob = 0 - Viscous Diffusion of 1-D Gaussian
//! iprob = 1 - Viscous Diffusion of 2-D Gaussian
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
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

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
  Real p0 = 1.;
  int iprob = pin->GetOrAddInteger("problem", "iprob", 0);
  Real amp = pin->GetOrAddReal("problem", "amp", 1.e-6);
  Real t0 = pin->GetOrAddReal("problem", "t0", 0.5);
  Real nu_iso = pin->GetOrAddReal("problem", "nu_iso", 0.25);
  Real x10 = pin->GetOrAddReal("problem", "x10", 0.);
  Real x20 = pin->GetOrAddReal("problem", "x20", 0.);
  Real gamma = peos->GetGamma();
  Real gm1 = gamma - 1.0;

  if (iprob == 0) { // 1-D diffusion
    if (std::strcmp(COORDINATE_SYSTEM, "cartesian") != 0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in visc.cpp ProblemGenerator" << std::endl
          << "1-D diffusion test only compatible with cartesian coord" << std::endl;
      ATHENA_ERROR(msg);
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          x1=pcoord->x1v(i);
          phydro->u(IDN,k,j,i) = d0;
          phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i)*v1;
          phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i)
                                 * amp/std::pow(std::sqrt(4.*PI*nu_iso*t0),1.0)
                                 * std::exp(-(std::pow(x1-x10,2.))
                                            / (4.*nu_iso*t0));
          phydro->u(IM3,k,j,i) = phydro->u(IDN,k,j,i)*v3;
          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) =
                p0/gm1 + 0.5*(SQR(phydro->u(IM1,k,j,i)) +
                              SQR(phydro->u(IM2,k,j,i)) +
                              SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
          }
        }
      }
    }
  } else if (iprob == 1) { // 2-D diffusion
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
            x1 = pcoord->x1v(i);
            x2 = pcoord->x2v(j);
            phydro->u(IDN,k,j,i) = d0;
            phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i)*v1;
            phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i)*v2;
            phydro->u(IM3,k,j,i) = phydro->u(IDN,k,j,i)
                                   * amp/std::pow(std::sqrt(4.*PI*nu_iso*t0),2.0)
                                   * std::exp(-(std::pow(x1-x10,2.)
                                                + std::pow(x2-x20,2.))
                                              / (4.*nu_iso*t0));
            if (NON_BAROTROPIC_EOS) {
              phydro->u(IEN,k,j,i) =
                  p0/gm1 + 0.5*(SQR(phydro->u(IM1,k,j,i)) +
                                SQR(phydro->u(IM2,k,j,i)) +
                                SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
            }
          } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
            r = pcoord->x1v(i);
            phi = pcoord->x2v(j);
            x1 = r*std::cos(phi);
            x2 = r*std::sin(phi);
            phydro->u(IDN,k,j,i) = d0;
            phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i)*v1;
            phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i)*v2;
            phydro->u(IM3,k,j,i) = phydro->u(IDN,k,j,i)
                                   * amp/std::pow(std::sqrt(4.*PI*nu_iso*t0),2.0)
                                   * std::exp(-(std::pow(x1-x10,2.)
                                                + std::pow(x2-x20,2.))
                                              / (4.*nu_iso*t0));
            if (NON_BAROTROPIC_EOS) {
              phydro->u(IEN,k,j,i) =
                  p0/gm1 + 0.5*(SQR(phydro->u(IM1,k,j,i)) +
                                SQR(phydro->u(IM2,k,j,i)) +
                                SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
            }
          } else {
            std::stringstream msg;
            msg << "### FATAL ERROR in visc.cpp ProblemGenerator" << std::endl
                << "2-d diffusion test only compatible with cartesian or"
                << std::endl << "cylindrical coord" << std::endl;
            ATHENA_ERROR(msg);
          }
        }
      }
    }
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in visc.cpp ProblemGenerator" << std::endl
        << "iprob must be set to 0 or 1" << std::endl;
    ATHENA_ERROR(msg);
  }
  return;
}

// refinement condition: check the maximum z-velocity gradient
int RefinementCondition(MeshBlock *pmb) {
  int f2 = pmb->pmy_mesh->f2, f3 = pmb->pmy_mesh->f3;
  AthenaArray<Real> &w = pmb->phydro->w;
  Real maxeps = 0.;
  if (f2) {
    int k = pmb->ks;
    for (int j=pmb->js-1; j<=pmb->je+1; j++) {
      for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
        Real eps = std::sqrt(SQR(0.5*(w(IVZ,k,j,i+1) - w(IVZ,k,j,i-1)))
                             + SQR(0.5*(w(IVZ,k,j+1,i) - w(IVZ,k,j-1,i))));
        maxeps = std::max(maxeps, eps);
      }
    }
  } else {
    int k = pmb->ks;
    int j = pmb->js;
    for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
      Real eps = std::sqrt(SQR(0.5*(w(IVY,k,j,i+1) - w(IVY,k,j,i-1))));
      maxeps = std::max(maxeps, eps);
    }
  }

  if (maxeps > threshold) return 1;
  if (maxeps < 0.25*threshold) return -1;
  return 0;
}
