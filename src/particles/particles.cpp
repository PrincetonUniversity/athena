//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file particles.cpp
//! \brief implements functions in particle classes

// C++ Standard Libraries
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "particles.hpp"

// Class variable initialization
bool Particles::initialized = false;
int Particles::num_particles = 0;
int Particles::num_particles_grav = 0;
int Particles::num_particles_output = 0;
ParameterInput* Particles::pinput = NULL;
std::vector<int> Particles::idmax;
#ifdef MPI_PARALLEL
MPI_Comm Particles::my_comm = MPI_COMM_NULL;
#endif

//--------------------------------------------------------------------------------------
//! \fn Particles::Initialize(Mesh *pm, ParameterInput *pin)
//! \brief initializes the class.

void Particles::Initialize(Mesh *pm, ParameterInput *pin) {
  if (initialized) return;

  InputBlock *pib = pin->pfirst_block;
  // pm->particle_params.reserve(1);
  // loop over input block names.  Find those that start with "particle", read parameters,
  // and construct singly linked list of ParticleTypes.
  while (pib != nullptr) {
    if (pib->block_name.compare(0, 8, "particle") == 0) {
      ParticleParameters pp;  // define temporary ParticleParameters struct

      // extract integer number of particle block.  Save name and number
      std::string parn = pib->block_name.substr(8); // 7 because counting starts at 0!
      pp.block_number = atoi(parn.c_str());
      pp.block_name.assign(pib->block_name);

      // set particle type = [tracer, star, dust, none]
      pp.partype = pin->GetString(pp.block_name, "type");
      if (pp.partype.compare("none") != 0) { // skip input block if the type is none
        if ((pp.partype.compare("dust") == 0) ||
            (pp.partype.compare("tracer") == 0) ||
            (pp.partype.compare("star") == 0)) {
          pp.ipar = num_particles++;
          idmax.push_back(0); // initialize idmax with 0
          pp.table_output = pin->GetOrAddBoolean(pp.block_name,"output",false);
          pp.gravity = pin->GetOrAddBoolean(pp.block_name,"gravity",false);
          if (pp.table_output) num_particles_output++;
          if (pp.gravity) num_particles_grav++;
          pm->particle_params.push_back(pp);
        } else { // unsupported particle type
          std::stringstream msg;
          msg << "### FATAL ERROR in Particle Initializer" << std::endl
              << "Unrecognized particle type = '" << pp.partype
              << "' in particle block '" << pp.block_name << "'" << std::endl;
          ATHENA_ERROR(msg);
        }
      }
    }
    pib = pib->pnext;  // move to next input block name
  }

  if (num_particles > 0) {
    pm->particle = true;

    if (SELF_GRAVITY_ENABLED && (num_particles_grav > 0))
      pm->particle_gravity = true;

    if (Globals::my_rank == 0) {
      std::cout << "Particles are initalized with "
                << "N containers = " << num_particles << std::endl;
      for (ParticleParameters ppnew : pm->particle_params)
        std::cout << "  ipar = " << ppnew.ipar
                  << "  partype = " << ppnew.partype
                  << "  output = " << ppnew.table_output
                  << "  block_name = " << ppnew.block_name
                  << std::endl;
    }
    // Remember the pointer to input parameters.
    pinput = pin;

#ifdef MPI_PARALLEL
    // Get my MPI communicator.
    MPI_Comm_dup(MPI_COMM_WORLD, &my_comm);
#endif
  }

  initialized = true;
}

//--------------------------------------------------------------------------------------
//! \fn Particles::PostInitialize(Mesh *pm, ParameterInput *pin)
//! \brief preprocesses the class after problem generator and before the main loop.

void Particles::PostInitialize(Mesh *pm, ParameterInput *pin) {
  // Set particle IDs.
  for (int ipar = 0; ipar < Particles::num_particles; ++ipar)
    ProcessNewParticles(pm, ipar);

  // Set position indices.
  for (int b = 0; b < pm->nblocal; ++b)
    for (Particles *ppar : pm->my_blocks(b)->ppar)
      ppar->SetPositionIndices();

  // Print particle csv
  for (int b = 0; b < pm->nblocal; ++b)
    for (Particles *ppar : pm->my_blocks(b)->ppar)
      if (ppar->parhstout_) ppar->OutputParticles(true);
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::FindDensityOnMesh(Mesh *pm, bool include_momentum)
//! \brief finds particle mesh densities for all particle containers.
//!
//! If include_momentum is true, the momentum density field is also computed.

void Particles::FindDensityOnMesh(Mesh *pm, bool include_momentum) {
  // Assign particle properties to mesh and send boundary.
  int nblocks(pm->nblocal);

}

//--------------------------------------------------------------------------------------
//! \fn void Particles::GetHistoryOutputNames(std::string output_names[])
//! \brief gets the names of the history output variables in history_output_names[].

void Particles::GetHistoryOutputNames(std::string output_names[], int ipar) {
  std::string head = "p";
  head.append(std::to_string(ipar));
  output_names[0] = head + "-n";
  output_names[1] = head + "-v1";
  output_names[2] = head + "-v2";
  output_names[3] = head + "-v3";
  output_names[4] = head + "-v1sq";
  output_names[5] = head + "-v2sq";
  output_names[6] = head + "-v3sq";
  output_names[7] = head + "-m";
}

//--------------------------------------------------------------------------------------
//! \fn int Particles::GetTotalNumber(Mesh *pm)
//! \brief returns total number of particles (from all processes).
//! \todo This should separately count different types of particles
std::int64_t Particles::GetTotalNumber(Mesh *pm) {
  std::int64_t npartot(0);
  for (int b = 0; b < pm->nblocal; ++b)
    for (Particles *ppar : pm->my_blocks(b)->ppar)
      npartot += ppar->npar;
#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, &npartot, 1, MPI_LONG, MPI_SUM, my_comm);
#endif
  return npartot;
}

//--------------------------------------------------------------------------------------
//! \fn Particles::Particles(MeshBlock *pmb, ParameterInput *pin)
//! \brief constructs a Particles instance.

Particles::Particles(MeshBlock *pmb, ParameterInput *pin, ParticleParameters *pp) :
  input_block_name(pp->block_name), partype(pp->partype),
  nint(0), nreal(0), naux(0), nwork(0),
  ipid(-1), ixp(-1), iyp(-1), izp(-1), ivpx(-1), ivpy(-1), ivpz(-1),
  ixp0(-1), iyp0(-1), izp0(-1), ivpx0(-1), ivpy0(-1), ivpz0(-1),
  ixi1(-1), ixi2(-1), ixi3(-1), imom1(-1), imom2(-1), imom3(-1), imass(-1), ish(-1),
  igx(-1), igy(-1), igz(-1),
  npar(0), nparmax(1),
  my_ipar_(pp->ipar), isgravity_(false), parhstout_(false), mass(1.0) {
  // Add particle ID.
  ipid = AddIntProperty();
  intfieldname.push_back("pid");


  // Actual memory allocation and shorthand assignment will be done in the derived class
  // Initialization of ParticleBuffer, ParticleGravity
  // has moved to the derived class
}

//--------------------------------------------------------------------------------------
//! \fn Particles::~Particles()
//! \brief destroys a Particles instance.

Particles::~Particles() {
  // Delete integer properties.
  intprop.DeleteAthenaArray();
  intfieldname.clear();

  // Delete real properties.
  realprop.DeleteAthenaArray();
  realfieldname.clear();

  // Delete auxiliary properties.
  if (naux > 0) {
    auxprop.DeleteAthenaArray();
    auxfieldname.clear();
  }

  // Delete working arrays.
  if (nwork > 0) work.DeleteAthenaArray();

  // Clear links to neighbors.
  ClearNeighbors();

  // Delete mesh auxiliaries.
  delete ppm;

  // Delete particle gravity.
  if (isgravity_) delete ppgrav;
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::AllocateMemory()
//! \brief memory allocation will be done at the end of derived class initialization
void Particles::AllocateMemory() {
  // Initiate ParticleBuffer class.
  nint_buf = nint;
  nreal_buf = nreal + naux;

  // Allocate mesh auxiliaries.
  ppm = new ParticleMesh(this, pmy_block);
  imom1 = ppm->imom1;
  imom2 = ppm->imom2;
  imom3 = ppm->imom3;

  // Allocate particle gravity
  if (isgravity_) {
    // Add working arrays for gravity forces
    igx = AddWorkingArray();
    igy = AddWorkingArray();
    igz = AddWorkingArray();
    // Activate particle gravity.
    ppgrav = new ParticleGravity(this);
  }

  // Allocate integer properties.
  intprop.NewAthenaArray(nint,nparmax);

  // Allocate integer properties.
  realprop.NewAthenaArray(nreal,nparmax);

  // Allocate auxiliary properties.
  if (naux > 0) auxprop.NewAthenaArray(naux,nparmax);

  // Allocate working arrays.
  if (nwork > 0) work.NewAthenaArray(nwork,nparmax);
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::FindHistoryOutput(Real data_sum[], int pos)
//! \brief finds the data sums of history output from particles in my process and assign
//!   them to data_sum beginning at index pos.

void Particles::AddHistoryOutput(Real data_sum[], int pos) {
  const int NSUM = NHISTORY - 1;

}

//--------------------------------------------------------------------------------------
//! \fn void Particles::CheckInMeshBlock()
//! \brief check whether given position is within the meshblock assuming Cartesian

bool Particles::CheckInMeshBlock(Real x1, Real x2, Real x3) {
  RegionSize& bsize = pmy_block->block_size;
  if ((x1>=bsize.x1min) && (x1<bsize.x1max) &&
      (x2>=bsize.x2min) && (x2<bsize.x2max) &&
      (x3>=bsize.x3min) && (x3<bsize.x3max)) {
    return true;
  } else {
    return false;
  }
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::AddOneParticle()
//! \brief add one particle if position is within the mesh block

void Particles::AddOneParticle(Real x1, Real x2, Real x3,
  Real v1, Real v2, Real v3) {
  if (CheckInMeshBlock(x1,x2,x3)) {
    if (npar == nparmax) UpdateCapacity(npar*2);
    pid(npar) = -1;
    xp(npar) = x1;
    yp(npar) = x2;
    zp(npar) = x3;
    vpx(npar) = v1;
    vpy(npar) = v2;
    vpz(npar) = v3;
    npar++;
  }
}

//--------------------------------------------------------------------------------------
//! \fn AthenaArray<Real> Particles::GetVelocityField()
//! \brief returns the particle velocity on the mesh.
//!
//! \note
//!   Precondition:
//!   The particle properties on mesh must be assigned using the class method
//!   Particles::FindDensityOnMesh().

AthenaArray<Real> Particles::GetVelocityField() const {
  AthenaArray<Real> vel(3, ppm->nx3_, ppm->nx2_, ppm->nx1_);
  for (int k = ppm->ks; k <= ppm->ke; ++k)
    for (int j = ppm->js; j <= ppm->je; ++j)
      for (int i = ppm->is; i <= ppm->ie; ++i) {
        Real rho;
        if (imass == -1) {
          rho = ppm->weight(k,j,i);
        } else {
          rho = ppm->density(k,j,i);
        }
        rho = (rho > 0.0) ? rho : 1.0;

        vel(0,k,j,i) = ppm->meshaux(imom1,k,j,i) / rho;
        vel(1,k,j,i) = ppm->meshaux(imom2,k,j,i) / rho;
        vel(2,k,j,i) = ppm->meshaux(imom3,k,j,i) / rho;
      }
  return vel;
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::Integrate(int step)
//! \brief updates all particle positions and velocities from t to t + dt.

void Particles::Integrate(int stage) {
  Real t = 0, dt = 0;


}

//--------------------------------------------------------------------------------------
//! \fn void Particles::RemoveOneParticle(int k)
//! \brief removes particle k in the block.

void Particles::RemoveOneParticle(int k) {
  if (0 <= k && k < npar && --npar != k) {
    xi1(k) = xi1(npar);
    xi2(k) = xi2(npar);
    xi3(k) = xi3(npar);
    for (int j = 0; j < nint; ++j)
      intprop(j,k) = intprop(j,npar);
    for (int j = 0; j < nreal; ++j)
      realprop(j,k) = realprop(j,npar);
    for (int j = 0; j < naux; ++j)
      auxprop(j,k) = auxprop(j,npar);
  }
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::SetPositionIndices()
//! \brief updates position indices of particles.

void Particles::SetPositionIndices() {
  GetPositionIndices(npar, xp, yp, zp, xi1, xi2, xi3);
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::ProcessNewParticles()
//! \brief searches for and books new particles.

void Particles::ProcessNewParticles(Mesh *pmesh, int ipar) {
  // Count new particles.
  const int nbtotal(pmesh->nbtotal), nblocks(pmesh->nblocal);
  std::vector<int> nnewpar(nbtotal, 0);
  for (int b = 0; b < nblocks; ++b) {
    const MeshBlock *pmb(pmesh->my_blocks(b));
    nnewpar[pmb->gid] = pmb->ppar[ipar]->CountNewParticles();
  }
#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, &nnewpar[0], nbtotal, MPI_INT, MPI_MAX, my_comm);
#endif

  // Make the counts cumulative.
  for (int i = 1; i < nbtotal; ++i)
    nnewpar[i] += nnewpar[i-1];

  // Set particle IDs.
  for (int b = 0; b < nblocks; ++b) {
    const MeshBlock *pmb(pmesh->my_blocks(b));
    int newid_start = idmax[ipar] + (pmb->gid > 0 ? nnewpar[pmb->gid - 1] : 0);
    pmb->ppar[ipar]->SetNewParticleID(newid_start);
  }
  idmax[ipar] += nnewpar[nbtotal - 1];
}

//--------------------------------------------------------------------------------------
//! \fn int Particles::CountNewParticles()
//! \brief counts new particles in the block.

int Particles::CountNewParticles() const {
  int n = 0;
  for (int i = 0; i < npar; ++i)
    if (pid(i) <= 0) ++n;
  return n;
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::EulerStep(Real t, Real dt, const AthenaArray<Real>& meshsrc)
//! \brief evolves the particle positions and velocities by one Euler step.

void Particles::EulerStep(Real t, Real dt, const AthenaArray<Real>& meshsrc) {
  // Update positions.
  for (int k = 0; k < npar; ++k) {
    //! \todo (ccyang):
    //! - This is a temporary hack.
    Real tmpx = xp(k), tmpy = yp(k), tmpz = zp(k);
    xp(k) = xp0(k) + dt * vpx(k);
    yp(k) = yp0(k) + dt * vpy(k);
    zp(k) = zp0(k) + dt * vpz(k);
    xp0(k) = tmpx;
    yp0(k) = tmpy;
    zp0(k) = tmpz;
  }

  // Integrate the source terms (e.g., acceleration).
  SourceTerms(t, dt, meshsrc);
  UserSourceTerms(t, dt, meshsrc);
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::GetPositionIndices(int npar,
//!                                        const AthenaArray<Real>& xp,
//!                                        const AthenaArray<Real>& yp,
//!                                        const AthenaArray<Real>& zp,
//!                                        AthenaArray<Real>& xi1,
//!                                        AthenaArray<Real>& xi2,
//!                                        AthenaArray<Real>& xi3)
//! \brief finds the position indices of each particle with respect to the local grid.

void Particles::GetPositionIndices(int npar,
                                   const AthenaArray<Real>& xp,
                                   const AthenaArray<Real>& yp,
                                   const AthenaArray<Real>& zp,
                                   AthenaArray<Real>& xi1,
                                   AthenaArray<Real>& xi2,
                                   AthenaArray<Real>& xi3) {

}

//--------------------------------------------------------------------------------------
//! \fn void Particles::SetNewParticleID(int id0)
//! \brief searches for new particles and assigns ID, beginning at id + 1.

void Particles::SetNewParticleID(int id) {
  for (int i = 0; i < npar; ++i)
    if (pid(i) <= 0) pid(i) = ++id;
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::SaveStatus()
//! \brief saves the current positions and velocities for later use.

void Particles::SaveStatus() {
  for (int k = 0; k < npar; ++k) {
    // Save current positions.
    xp0(k) = xp(k);
    yp0(k) = yp(k);
    zp0(k) = zp(k);

    // Save current velocities.
    vpx0(k) = vpx(k);
    vpy0(k) = vpy(k);
    vpz0(k) = vpz(k);
  }
}

//--------------------------------------------------------------------------------------
//! \fn int Particles::AddIntProperty()
//! \brief adds one integer property to the particles and returns the index.

int Particles::AddIntProperty() {
  return nint++;
}

//--------------------------------------------------------------------------------------
//! \fn int Particles::AddRealProperty()
//! \brief adds one real property to the particles and returns the index.

int Particles::AddRealProperty() {
  return nreal++;
}

//--------------------------------------------------------------------------------------
//! \fn int Particles::AddAuxProperty()
//! \brief adds one auxiliary property to the particles and returns the index.

int Particles::AddAuxProperty() {
  return naux++;
}

//--------------------------------------------------------------------------------------
//! \fn int Particles::AddWorkingArray()
//! \brief adds one working array to the particles and returns the index.

int Particles::AddWorkingArray() {
  return nwork++;
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::AssignShorthands()
//! \brief assigns shorthands by shallow copying slices of the data.

void Particles::AssignShorthands() {
  pid.InitWithShallowSlice(intprop, 2, ipid, 1);

  xp.InitWithShallowSlice(realprop, 2, ixp, 1);
  yp.InitWithShallowSlice(realprop, 2, iyp, 1);
  zp.InitWithShallowSlice(realprop, 2, izp, 1);
  vpx.InitWithShallowSlice(realprop, 2, ivpx, 1);
  vpy.InitWithShallowSlice(realprop, 2, ivpy, 1);
  vpz.InitWithShallowSlice(realprop, 2, ivpz, 1);

  xp0.InitWithShallowSlice(auxprop, 2, ixp0, 1);
  yp0.InitWithShallowSlice(auxprop, 2, iyp0, 1);
  zp0.InitWithShallowSlice(auxprop, 2, izp0, 1);
  vpx0.InitWithShallowSlice(auxprop, 2, ivpx0, 1);
  vpy0.InitWithShallowSlice(auxprop, 2, ivpy0, 1);
  vpz0.InitWithShallowSlice(auxprop, 2, ivpz0, 1);

  xi1.InitWithShallowSlice(work, 2, ixi1, 1);
  xi2.InitWithShallowSlice(work, 2, ixi2, 1);
  xi3.InitWithShallowSlice(work, 2, ixi3, 1);
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::UpdateCapacity(int new_nparmax)
//! \brief changes the capacity of particle arrays while preserving existing data.

void Particles::UpdateCapacity(int new_nparmax) {

}

//--------------------------------------------------------------------------------------
//! \fn Real Particles::NewBlockTimeStep();
//! \brief returns the time step required by particles in the block.

Real Particles::NewBlockTimeStep() {
  Coordinates *pc = pmy_block->pcoord;

  // Find the maximum coordinate speed.
  Real dt_inv2_max = 0.0;


  // Return the time step constrained by the coordinate speed.
  return dt_inv2_max > 0.0 ? cfl_par / std::sqrt(dt_inv2_max)
                           : std::numeric_limits<Real>::max();
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::FindLocalDensityOnMesh(Mesh *pm, bool include_momentum)
//! \brief finds the number density of particles on the mesh.
//!
//!   If include_momentum is true, the momentum density field is also computed,
//!   assuming mass of each particle is unity.
//! \note
//!   Postcondition: ppm->weight becomes the density in each cell, and
//!   if include_momentum is true, ppm->meshaux(imom1:imom3,:,:,:)
//!   becomes the momentum density.

void Particles::FindLocalDensityOnMesh(bool include_momentum) {
  Coordinates *pc(pmy_block->pcoord);


}
//--------------------------------------------------------------------------------------
//! \fn void Particles::ConvertToDensity(bool include_momentum)
//! \brief finds the number density of particles on the mesh.
//!
//!   If include_momentum is true, the momentum density field is also computed,
//!   assuming mass of each particle is unity.
//! \note
//!   Postcondition: ppm->weight becomes the density in each cell, and
//!   if include_momentum is true, ppm->meshaux(imom1:imom3,:,:,:)
//!   becomes the momentum density.

void Particles::ConvertToDensity(bool include_momentum) {
  Coordinates *pc(pmy_block->pcoord);
  // Convert to densities.
  int is = active1_ ? ppm->is-NGHOST: ppm->is, ie = active1_ ? ppm->ie+NGHOST: ppm->ie;
  int js = active2_ ? ppm->js-NGHOST: ppm->js, je = active2_ ? ppm->je+NGHOST: ppm->je;
  int ks = active3_ ? ppm->ks-NGHOST: ppm->ks, ke = active3_ ? ppm->ke+NGHOST: ppm->ke;
  if (include_momentum) {
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          Real vol(pc->GetCellVolume(k,j,i));
          Real rhop(mass/vol); // mass = 1.0 if imass != -1
          ppm->weight(k,j,i) *= rhop;
          ppm->meshaux(ppm->imom1,k,j,i) *= rhop;
          ppm->meshaux(ppm->imom2,k,j,i) *= rhop;
          ppm->meshaux(ppm->imom3,k,j,i) *= rhop;
          if (ppm->imass != -1) ppm->density(k,j,i) *= rhop;
        }
  } else {
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          Real vol(pc->GetCellVolume(k,j,i));
          Real rhop(mass/vol); // mass = 1.0 if imass != -1
          ppm->weight(k,j,i) *= rhop;
          if (ppm->imass != -1) ppm->density(k,j,i) *= rhop;
        }
  }
}

//--------------------------------------------------------------------------------------
//! \fn std::size_t Particles::GetSizeInBytes()
//! \brief returns the data size in bytes in the meshblock.

std::size_t Particles::GetSizeInBytes() {
  std::size_t size = sizeof(npar);
  if (npar > 0) size += npar * (nint * sizeof(int) + nreal * sizeof(Real));
  return size;
}

//--------------------------------------------------------------------------------------
//! \fn Particles::UnpackParticlesForRestart()
//! \brief reads the particle data from the restart file.

void Particles::UnpackParticlesForRestart(char *mbdata, std::size_t &os) {
  // Read number of particles.
  std::memcpy(&npar, &(mbdata[os]), sizeof(npar));
  os += sizeof(npar);
  if (nparmax < npar)
    UpdateCapacity(npar);

  if (npar > 0) {
    // Read integer properties.
    std::size_t size = npar * sizeof(int);
    for (int k = 0; k < nint; ++k) {
      std::memcpy(&(intprop(k,0)), &(mbdata[os]), size);
      os += size;
    }

    // Read real properties.
    size = npar * sizeof(Real);
    for (int k = 0; k < nreal; ++k) {
      std::memcpy(&(realprop(k,0)), &(mbdata[os]), size);
      os += size;
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn Particles::PackParticlesForRestart()
//! \brief pack the particle data for restart dump.

void Particles::PackParticlesForRestart(char *&pdata) {
  // Write number of particles.
  std::memcpy(pdata, &npar, sizeof(npar));
  pdata += sizeof(npar);

  if (npar > 0) {
    // Write integer properties.
    std::size_t size = npar * sizeof(int);
    for (int k = 0; k < nint; ++k) {
      std::memcpy(pdata, &(intprop(k,0)), size);
      pdata += size;
    }
    // Write real properties.
    size = npar * sizeof(Real);
    for (int k = 0; k < nreal; ++k) {
      std::memcpy(pdata, &(realprop(k,0)), size);
      pdata += size;
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn Particles::FormattedTableOutput()
//! \brief outputs the particle data in tabulated format.

void Particles::FormattedTableOutput(Mesh *pm, OutputParameters op) {
  std::stringstream fname, msg;
  std::ofstream os;

  // Loop over Particle containers
  for (int ipar = 0; ipar < num_particles ; ++ipar) {
    if (pm->particle_params[ipar].table_output) {
      // Loop over MeshBlocks
      for (int b = 0; b < pm->nblocal; ++b) {
        const MeshBlock *pmb(pm->my_blocks(b));
        const Particles *ppar(pmb->ppar[ipar]);

        // Create the filename.
        fname << op.file_basename
              << ".block" << pmb->gid << '.' << op.file_id
              << '.' << std::setw(5) << std::right << std::setfill('0') << op.file_number
              << '.' << "par" << ipar << ".tab";

        // Open the file for write.
        os.open(fname.str().data());
        if (!os.is_open()) {
          msg << "### FATAL ERROR in function [Particles::FormattedTableOutput]"
              << std::endl << "Output file '" << fname.str() << "' could not be opened"
              << std::endl;
          ATHENA_ERROR(msg);
        }

        // Write the time.
        os << std::scientific << std::showpoint << std::setprecision(18);
        os << "# Athena++ particle data at time = " << pm->time << std::endl;

        // Write header.
        os << "# ";
        for (int ip = 0; ip < ppar->nint; ++ip)
          os << ppar->intfieldname[ip] << "  ";
        for (int ip = 0; ip < ppar->nreal; ++ip)
          os << ppar->realfieldname[ip] << "  ";
        os << std::endl;

        // Write the particle data in the meshblock.
        for (int k = 0; k < ppar->npar; ++k) {
          for (int ip = 0; ip < ppar->nint; ++ip)
            os << ppar->intprop(ip,k) << "  ";
          for (int ip = 0; ip < ppar->nreal; ++ip)
            os << ppar->realprop(ip,k) << "  ";
          os << std::endl;
        }

        // Close the file and get the next meshblock.
        os.close();
        fname.str("");
      }
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn Particles::OutputParticles()
//! \brief outputs the particle data in tabulated format.
void Particles::OutputParticles(bool header) {
  std::stringstream fname, msg;
  std::ofstream os;
  std::string file_basename = pinput->GetString("job","problem_id");

  for (int k = 0; k < npar; ++k) {
    // Create the filename.
    fname << file_basename << ".par" << pid(k) << ".csv";

    // Open the file for write.
    if (header)
      os.open(fname.str().data(), std::ofstream::out);
    else
      os.open(fname.str().data(), std::ofstream::app);

    if (!os.is_open()) {
      msg << "### FATAL ERROR in function [Particles::OutputParticles]"
          << std::endl << "Output file '" << fname.str() << "' could not be opened"
          << std::endl;
      ATHENA_ERROR(msg);
    }

    OutputOneParticle(os, k, header);

    // Close the file
    os.close();
    // clear filename
    fname.str("");
  }
}

//--------------------------------------------------------------------------------------
//! \fn Particles::OutputParticles()
//! \brief outputs the particle data in tabulated format.
void Particles::OutputParticles(bool header, int kid) {
  std::stringstream fname, msg;
  std::ofstream os;
  std::string file_basename = pinput->GetString("job","problem_id");

  for (int k = 0; k < npar; ++k) {
    if (pid(k) != kid) continue;

    // Create the filename.
    fname << file_basename << ".pid" << pid(k) << ".par" << my_ipar_ << ".csv";

    // Open the file for write.
    if (header)
      os.open(fname.str().data(), std::ofstream::out);
    else
      os.open(fname.str().data(), std::ofstream::app);

    if (!os.is_open()) {
      msg << "### FATAL ERROR in function [Particles::OutputParticles]"
          << std::endl << "Output file '" << fname.str() << "' could not be opened"
          << std::endl;
      ATHENA_ERROR(msg);
    }

    OutputOneParticle(os, k, header);

    // Close the file
    os.close();
    // clear filename
    fname.str("");
  }
}

//--------------------------------------------------------------------------------------
//! \fn Particles::OutputParticle()
//! \brief outputs the particle data in tabulated format.
void Particles::OutputOneParticle(std::ostream &os, int k, bool header) {
  if (header) {
    os << "time,dt";
    for (int ip = 0; ip < nint; ++ip)
      os << "," << intfieldname[ip];
    for (int ip = 0; ip < nreal; ++ip)
      os << "," << realfieldname[ip];
    for (int ip = 0; ip < naux; ++ip)
      os << "," << auxfieldname[ip];
    os << std::endl;
  }

  // Write the time.
  os << std::scientific << std::showpoint << std::setprecision(18);
  os << pmy_mesh->time << "," << pmy_mesh->dt;

  // Write the particle data in the meshblock.
  for (int ip = 0; ip < nint; ++ip)
    os << "," << intprop(ip,k);
  for (int ip = 0; ip < nreal; ++ip)
    os << "," << realprop(ip,k);
  for (int ip = 0; ip < naux; ++ip)
    os << "," << auxprop(ip,k);
  os << std::endl;
}

//--------------------------------------------------------------------------------------
//! \fn Particles::ToggleParHstOutFlag()
//! \brief turn on individual particle history outputs
void Particles::ToggleParHstOutFlag() {
  if (npar < 100) {
    parhstout_ = true;
  } else {
    std::cout << "Warning [Particles]: npar = " << npar << " is too large to output"
      << "all individual particles' history automatically."
      << " Particle history output is turned off." << std::endl;
    parhstout_ = false;
  }
}

//--------------------------------------------------------------------------------------
//! \fn Particles::DepositPMtoMesh()
//! \brief deposit PM momentum to hydro vars
//!
//! this has to be tested
void Particles::DepositPMtoMesh(int stage) {
  // Deposit ParticleMesh meshaux to MeshBlock.
  Hydro *phydro = pmy_block->phydro;
  Real t = pmy_mesh->time, dt = pmy_mesh->dt;

  switch (stage) {
  case 1:
    dt = 0.5 * dt;
    break;

  case 2:
    t += 0.5 * dt;
    break;
  }

  DepositToMesh(t, dt, phydro->w, phydro->u);
}
