#ifndef PARTICLES_PARTICLE_BUFFER_HPP_
#define PARTICLES_PARTICLE_BUFFER_HPP_
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
//======================================================================================
//! \file particle_buffer.hpp
//! \brief defines ParticleBuffer class for communication of particles.
//======================================================================================

// Athena++ headers
#include "../athena.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

//--------------------------------------------------------------------------------------
//! \class ParticleBuffer
//! \brief defines the class for managing buffers for transporting particles.

class ParticleBuffer {
friend class Particles;

 public:
  // Constructors and destructor
  ParticleBuffer();
  ParticleBuffer(int nparmax0, int nint, int nreal);
  ~ParticleBuffer();

  // Instance method
  void Reallocate(int new_nparmax, int nint, int nreal);

 protected:
  // Instance variables
  int* ibuf;   //!> ptr to integer buffer
  Real* rbuf;   //!> ptr to real buffer
  int nparmax;  //!> maximum number of particles
  int npar;     //!> actual number of particles in the buffer
#ifdef MPI_PARALLEL
  MPI_Request reqn, reqi, reqr;  //!> MPI request handles
  int flagn;               //!> Flag indicating if the incoming number is known
  int flagi, flagr;        //!> Flags indicating if the respective buffer is filled
  int tag;                 //!> MPI tag (allowing for from tag to tag + 2)
#endif
};
#endif  // PARTICLES_PARTICLE_BUFFER_HPP_
