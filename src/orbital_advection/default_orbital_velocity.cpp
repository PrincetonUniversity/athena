//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file default_orbital_velocity.cpp
//! \brief define default orbital velocity functions

// C/C++ headers
#include <algorithm>  // max(), min()
#include <cfloat>     // FLT_MAX, FLT_MIN
#include <iostream>   // cout, endl
#include <sstream>    //
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

// this class header
#include "orbital_advection.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif


// Default Orbital Velocity
Real CartOrbitalVelocity(OrbitalAdvection *porb, Real x_, Real y_, Real z_) {
  return -porb->qshear*porb->Omega0*x_;
}

Real CartOrbitalVelocity_x(OrbitalAdvection *porb, Real x_, Real y_, Real z_) {
  return -porb->qshear*porb->Omega0;
}

Real CylOrbitalVelocity2D(OrbitalAdvection *porb, Real x_, Real y_, Real z_) {
  return std::sqrt(porb->gm/x_)-porb->Omega0*x_;
}

Real CylOrbitalVelocity2D_r(OrbitalAdvection *porb, Real x_, Real y_, Real z_) {
  return -0.5*std::sqrt(porb->gm/x_)/x_-porb->Omega0;
}

Real CylOrbitalVelocity3D(OrbitalAdvection *porb, Real x_, Real y_, Real z_) {
  Real d_ = std::sqrt(SQR(x_)+SQR(z_));
  return std::sqrt(porb->gm/d_)*x_/d_-porb->Omega0*x_;
}

Real CylOrbitalVelocity3D_r(OrbitalAdvection *porb, Real x_, Real y_, Real z_) {
  Real d_ = std::sqrt(SQR(x_)+SQR(z_));
  return std::sqrt(porb->gm/d_)/d_*(1.0-1.5*SQR(x_/d_))-porb->Omega0;
}

Real CylOrbitalVelocity3D_z(OrbitalAdvection *porb, Real x_, Real y_, Real z_) {
  Real d_ = std::sqrt(SQR(x_)+SQR(z_));
  return -1.5*std::sqrt(porb->gm/d_)*x_*z_/SQR(d_);
}

Real SphOrbitalVelocity(OrbitalAdvection *porb, Real x_, Real y_, Real z_) {
  Real sin_theta = std::sin(y_);
  if (sin_theta == 0.0)
    return 0.0;
  else
    return std::sqrt(porb->gm/(x_*sin_theta))-porb->Omega0*x_*sin_theta;
}

Real SphOrbitalVelocity_r(OrbitalAdvection *porb, Real x_, Real y_, Real z_) {
  Real sin_theta = std::sin(y_);
  if (sin_theta == 0.0)
    return 0.0;
  else
    return -0.5*std::sqrt(porb->gm/(sin_theta*x_*x_*x_))-porb->Omega0*sin_theta;
}

Real SphOrbitalVelocity_t(OrbitalAdvection *porb, Real x_, Real y_, Real z_) {
  Real sin_theta = std::sin(y_);
  if (sin_theta == 0.0)
    return 0.0;
  else
    return -0.5*std::sqrt(porb->gm/(sin_theta*x_))*std::cos(y_)/sin_theta
         -porb->Omega0*x_*std::cos(y_);
}

Real ZeroOrbitalVelocity(OrbitalAdvection *porb, Real x_, Real y_, Real z_) {
  return 0.0;
}
