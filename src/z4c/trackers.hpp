#ifndef TRACKERS_HPP
#define TRACKERS_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file trackers.hpp
//  \brief definitions for the Tracker class

#include <string>

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../athena_tensor.hpp"

#define NDIM (3)

// Forward declaration
class Mesh;
class MeshBlock;
class MeshBlockTree;
class ParameterInput;
class Coordinates;
class ParameterInput;

//! \class Tracker
//! \brief Evolve the positions of the compact objects
//! This class performs the global reduction
class Tracker {
  public:
    //! Tracker position and velocity for each body
    struct Position_vars {
      Real pos[NDIM];
      Real betap[NDIM];
      };
    Position_vars pos_body[NPUNCT];
    int npunct;
    int root_lev;
    Real L_grid;
  public:
    //! Creates the Tracker object
    Tracker(Mesh * pmesh, ParameterInput * pin);
    //! Destructor (will close output file)
    ~Tracker();
    //! Call different initializations
    void Initialize(Mesh * pmesh, ParameterInput * pin);
    //! Initialize for one puncture case
    void InitializeOnepuncture(Mesh * pmesh, ParameterInput * pin);
    //! Initialize for two punctures case
    void InitializeTwopuncture(Mesh * pmesh, ParameterInput * pin);
    //! Reduces the data from all meshblocks and ranks
    void ReduceTracker();
    //! Call different integrators
    void EvolveTracker();
    //! Euler integrator
    void EvolveTrackerIntegrateEuler();
    //! Write data to file
    void WriteTracker(int iter, Real time) const;
  private:
    int root;
    Real times_in_block[NPUNCT];
    bool ioproc;
    std::string ofname;
    Mesh const * pmesh;
    FILE * pofile;
};

//! \class TrackerLocal
//! \brief Interpolate beta in current meshblock
//! This class performs interpolation of beta in current meshblock and check if block contains position
class TrackerLocal {
  MeshBlock * pmy_block;
  public:
    //! Store beta at previous timestep and presence in block
    struct Betap_vars {
      Real betap[NDIM];
      int  inblock;
     };
     Betap_vars betap[NPUNCT];
  public:
    //! Creates the TrackerLocal object
    TrackerLocal(MeshBlock * pmb, ParameterInput * pin);
     ~TrackerLocal();
    //! Check body's presence in block
    bool InBlock(int body);
    //! Interpolate and store beta at previous timestep
    void StoreBetaPrev(Betap_vars betap[NPUNCT], AthenaArray<Real> & u, int body);
};

#endif

