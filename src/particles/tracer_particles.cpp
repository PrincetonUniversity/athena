//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file tracer_particles.cpp
//! \brief implements functions in the TracerParticles class

// C++ headers
#include <algorithm>  // min()

// Athena++ headers
#include "../athena.hpp"
#include "../coordinates/coordinates.hpp"
#include "../gravity/gravity.hpp"
#include "../hydro/hydro.hpp"
#include "particles.hpp"

//--------------------------------------------------------------------------------------
//! \fn TracerParticles::TracerParticles(MeshBlock *pmb, ParameterInput *pin)
//! \brief constructs a TracerParticles instance.

TracerParticles::TracerParticles(MeshBlock *pmb, ParameterInput *pin,
  ParticleParameters *pp)
  : Particles(pmb, pin, pp) {
  // Add working array at particles for gas velocity/particle momentum change.
  iwx = AddWorkingArray();
  iwy = AddWorkingArray();
  iwz = AddWorkingArray();

  // Define mass.
  mass = pin->GetOrAddReal(input_block_name, "mass", 1.0);

  // allocate memory
  Particles::AllocateMemory();

  // Assign shorthands (need to do this for every constructor of a derived class)
  AssignShorthands();
}

//--------------------------------------------------------------------------------------
//! \fn TracerParticles::~TracerParticles()
//! \brief destroys a TracerParticles instance.

TracerParticles::~TracerParticles() {
  // nothing to do
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void TracerParticles::SetOneParticleMass(Real new_mass)
//! \brief sets the mass of each particle.

void TracerParticles::SetOneParticleMass(Real new_mass) {
  pinput->SetReal(input_block_name, "mass", mass = new_mass);
}

//--------------------------------------------------------------------------------------
//! \fn void TracerParticles::AssignShorthands()
//! \brief assigns shorthands by shallow coping slices of the data.

void TracerParticles::AssignShorthands() {
  Particles::AssignShorthands();
  wx.InitWithShallowSlice(work, 2, iwx, 1);
  wy.InitWithShallowSlice(work, 2, iwy, 1);
  wz.InitWithShallowSlice(work, 2, iwz, 1);
}

//--------------------------------------------------------------------------------------
//! \fn void TracerParticles::SourceTerms()
//! \brief adds acceleration to particles.

void TracerParticles::SourceTerms(Real t, Real dt, const AthenaArray<Real>& meshsrc) {
  ppm->InterpolateMeshToParticles(meshsrc, IVX, work, iwx, 3);

  // Transform the gas velocity into Cartesian.
  const Coordinates *pc = pmy_block->pcoord;

  // Tracer particles
  for (int k = 0; k < npar; ++k) {
    Real tmpx = vpx(k), tmpy = vpy(k), tmpz = vpz(k);
    vpx(k) = wx(k);
    vpy(k) = wy(k);
    vpz(k) = wz(k);
    vpx0(k) = tmpx; vpy0(k) = tmpy; vpz0(k) = tmpz;
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void TracerParticles::UserSourceTerms(Real t, Real dt,
//!                                         const AthenaArray<Real>& meshsrc)
//! \brief adds additional source terms to particles, overloaded by the user.

void __attribute__((weak)) TracerParticles::UserSourceTerms(
    Real t, Real dt, const AthenaArray<Real>& meshsrc) {
}

//--------------------------------------------------------------------------------------
//! \fn void TracerParticles::ReactToMeshAux(
//!              Real t, Real dt, const AthenaArray<Real>& meshsrc)
//! \brief Reacts to meshaux before boundary communications.

void TracerParticles::ReactToMeshAux(Real t, Real dt, const AthenaArray<Real>& meshsrc) {
  // Nothing to do for tracers
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void TracerParticles::DepositToMesh(Real t, Real dt,
//!              const AthenaArray<Real>& meshsrc, AthenaArray<Real>& meshdst);
//! \brief Deposits meshaux to Mesh.

void TracerParticles::DepositToMesh(
         Real t, Real dt, const AthenaArray<Real>& meshsrc, AthenaArray<Real>& meshdst) {
  // Nothing to do for tracers
  return;
}
