//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file particle_gravity.cpp
//! \brief implements the members of the ParticleGravity class.

// C++ standard libraries
#include <cstring>  // strcmp()
#include <sstream>  // stringstream

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "particle_gravity.hpp"
#include "particle_mesh.hpp"
#include "particles.hpp"

//--------------------------------------------------------------------------------------
//! \fn ParticleGravity::ParticleGravity(Particles *ppar)
//! \brief constructs a new ParticleGravity instance.

ParticleGravity::ParticleGravity(Particles *ppar) :
  igx(ppar->igx), igy(ppar->igy), igz(ppar->igz) {
  // Remember my parent Particles instance.
  pmy_par = ppar;
  pmy_pm = ppar->ppm;

  // Remember the coordinates.
  pcoord = ppar->pmy_block->pcoord;

  // Remember the dimensions of my meshblock.
  MeshBlock *pmb(ppar->pmy_block);
  ncells1 = pmb->ncells1;
  ncells2 = pmb->ncells2;
  ncells3 = pmb->ncells3;
  active1 = ncells1 > 1;
  active2 = ncells2 > 1;
  active3 = ncells3 > 1;

  // Allocate space for gravitational force.
  gforce.NewAthenaArray(3, ncells3, ncells2, ncells1);
}

//--------------------------------------------------------------------------------------
//! \fn ParticleGravity::~ParticleGravity(Particles *ppar)
//! \brief constructs a new ParticleGravity instance.

ParticleGravity::~ParticleGravity() {
  // Deallocate space for gravitational force.
  gforce.DeleteAthenaArray();
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleGravity::InterpolateGravitationalForce(Real dt)
//! \brief interpolartes the gravitational force onto each particle.
void ParticleGravity::InterpolateGravitationalForce() {
  // Interpolate the gravitational force onto each particle.
  pmy_pm->InterpolateMeshToParticles(gforce, 0, pmy_par->work, igx, 3);
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleGravity::ExertGravitationalForce(Real dt)
//! \brief exerts the gravitational force on each particle.
void ParticleGravity::ExertGravitationalForce(Real dt) {
  // Add the force.
  for (int k = 0; k < pmy_par->npar; ++k) {
    pmy_par->vpx(k) += dt * pmy_par->work(igx,k);
    pmy_par->vpy(k) += dt * pmy_par->work(igy,k);
    pmy_par->vpz(k) += dt * pmy_par->work(igz,k);
  }
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleGravity::FindGravitationalForce(const AthenaArray<Real>& phi)
//! \brief computes the gravitational force from the potential phi.
//!
//! \todo (ccyang)
//! - currently only works with uniform cartesian.

void ParticleGravity::FindGravitationalForce(const AthenaArray<Real>& phi) {
  // Sanity check.
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") != 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [ParticleGravity::FindGravitationalForce]"
        << std::endl
        << "not implemented for non-Cartesian coordinates. " << std::endl;
    ATHENA_ERROR(msg);
    return;
  }
  if ((active1 && pcoord->dx1v(0) != pcoord->dx1v(1)) ||
      (active2 && pcoord->dx2v(0) != pcoord->dx2v(1)) ||
      (active3 && pcoord->dx3v(0) != pcoord->dx3v(1))) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [ParticleGravity::FindGravitationalForce]"
        << std::endl
        << "not implemented for non-uniform grid. " << std::endl;
    ATHENA_ERROR(msg);
    return;
  }

  // Set up the loop dimensions.
  const int is = active1 ? 1 : 0;
  const int ie = active1 ? ncells1 - 2 : 0;
  const int js = active2 ? 1 : 0;
  const int je = active2 ? ncells2 - 2 : 0;
  const int ks = active3 ? 1 : 0;
  const int ke = active3 ? ncells3 - 2 : 0;

  // Get the grid spacing.
  Real a1(0.5 / pcoord->dx1v(0)), a2(0.5 / pcoord->dx2v(0)), a3(0.5 / pcoord->dx3v(0));

  // Compute the force.
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        gforce(0,k,j,i) = active1 ? a1 * (phi(k,j,i-1) - phi(k,j,i+1)) : 0.0;
        gforce(1,k,j,i) = active2 ? a2 * (phi(k,j-1,i) - phi(k,j+1,i)) : 0.0;
        gforce(2,k,j,i) = active3 ? a3 * (phi(k-1,j,i) - phi(k+1,j,i)) : 0.0;
      }
}
