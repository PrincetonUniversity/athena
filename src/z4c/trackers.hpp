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
class ParameterInput;

//! \class Tracker
//! \brief Evolves the position of the compact objects
//! This class performs the global reduction
class Tracker {
  public:
    //! Monopole term
    struct Position_vars {
      Real pos[NDIM];
      Real betap[NDIM];       // position body one
      };
    Position_vars pos_body[NPUNCT];
    int npunct;
  public:
    //! Creates the WaveExtract object
    Tracker(Mesh * pmesh, ParameterInput * pin);
    //! Destructor (will close output file)
    ~Tracker();
    void Initialize(Mesh * pmesh, ParameterInput * pin);
    void InitializeOnepuncture(Mesh * pmesh, ParameterInput * pin);
    void InitializeTwopuncture(Mesh * pmesh, ParameterInput * pin);
    void ReduceTracker();
    //! Reduces the data from all of the SphericalPatches
    void EvolveTracker();
    //void EvolveTrackerIntegrate(Position_vars pos_gen[2]);
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
//! \brief Evolve tracker positions in corresponding meshblocks
//! This class performs the reduction on each TODO
class TrackerLocal {
  MeshBlock * pmy_block;
  public:
    struct Betap_vars {
      Real betap[NDIM];
      int  inblock;
     };
     Betap_vars betap[NPUNCT]; 
  public:
    //! Creates the TrackerLocal object
    TrackerLocal(MeshBlock * pmb, ParameterInput * pin);
    ~TrackerLocal();
    //! Evolve Tracker in correspondin meshblocks
    //void EvolveTracker(Position_vars pos_gen[2]);
    //void EvolveTrackerIntegrate(Position_vars pos_gen[2]); 
    //void InitializeTracker(ParameterInput * pin, Position_vars pos_gen[2], int body);
    bool InBlock(int body);
    void StoreBetaPrev(Betap_vars betap[2], AthenaArray<Real> & u, int body);
  private:
    //AthenaArray<Real> data;
    //AthenaArray<Real> weight;
    //Real rad;
};

#endif

