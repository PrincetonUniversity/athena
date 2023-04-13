#ifndef PARTICLES_PARTICLE_GRAVITY_HPP_
#define PARTICLES_PARTICLE_GRAVITY_HPP_
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file particle_gravity.hpp
//! \brief provides the header for the ParticleGravity class.

// Forward definition
class ParticleMesh;
class Particles;

//--------------------------------------------------------------------------------------
//! \class ParticleGravity
//! \brief defines the class for managing the gravity on particles.

class ParticleGravity {
 public:
  // Constructor
  explicit ParticleGravity(Particles *ppar);

  // Destructor
  ~ParticleGravity();

  // Instance methods
  void InterpolateGravitationalForce();
  void ExertGravitationalForce(Real dt);
  void FindGravitationalForce(const AthenaArray<Real>& phi);

 private:
  // Instance variables
  int igx, igy, igz;  // indices to working arrays

  // Attributes
  AthenaArray<Real> gforce;        //!> gravitational force
  Coordinates *pcoord;             //!> pointer to the coordinates
  Particles *pmy_par;              //!> pointer to parent Particles instance
  ParticleMesh *pmy_pm;            //!> pointer to my ParticleMesh instance
  bool active1, active2, active3;  // whether or not a direction is active
  int ncells1, ncells2, ncells3;   // block dimensions
};

#endif  // PARTICLES_PARTICLE_GRAVITY_HPP_
