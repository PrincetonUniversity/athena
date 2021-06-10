#ifndef PUNCTURE_TRACKER_HPP
#define PUNCTURE_TRACKER_HPP

//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file puncture_tracker.hpp
//  \brief definitions for the PunctureTracker class

#include <string>

#include "../athena.hpp"
#include "z4c_macro.hpp"

// Forward declaration
class Mesh;
class ParameterInput;

//! \class PunctureTracker
//! \brief Tracks a single puncture
class PunctureTracker {
  public:
    //! Initialize a tracker
    PunctureTracker(Mesh * pmesh, ParameterInput * pin, int n);
    //! Destructor (will close output file)
    ~PunctureTracker();
    //! Interpolate the shift vector to the puncture position
    void InterpolateShift(MeshBlock * pmb, AthenaArray<Real> & u);
    //! Update and broadcast the puncture position
    void EvolveTracker();
    //! Write data to file
    void WriteTracker(int iter, Real time) const;
    //! Get position
    inline Real GetPos(int a) {
      return pos[a];
    }
    // These need to access the internals of PunctureTracker for checkpoint / recovery
    friend class Mesh;
    friend class RestartOutput;
  private:
    bool owns_puncture;
    Real pos[NDIM];
    Real betap[NDIM];
  private:
    Mesh const * pmesh;
  private:
    std::string ofname;
    FILE * pofile;
};

#endif
