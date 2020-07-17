#ifndef WAVE_EXTRACT_HPP
#define WAVE_EXTRACT_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file wave_extract.hpp
//  \brief definitions for the WaveExtract class

#include <string>

#include "../athena.hpp"
#include "../athena_arrays.hpp"

// Forward declaration
class Mesh;
class MeshBlock;
class SphericalGrid;
class SphericalPatch;
class ParameterInput;

//! \class WaveExtract
//! \brief Extracts the l m  components of the wave on a unit sphere
//! This class performs the global reduction
class WaveExtract {
  public:
    //! Creates the WaveExtract object
    WaveExtract(Mesh * pmesh, ParameterInput * pin);
    //! Destructor (will close output file)
    ~WaveExtract();
    //! Reduces the data from all of the SphericalPatches
    void ReduceMultipole();
    //! Write data to file
    void Write(int iter, Real time) const;
  public:
    //!  Array of lm modes
    AthenaArray<Real> psi;    
    //! SphericalGrid for wave extraction
    SphericalGrid * psphere;
  private:
    int root;
    int lmax;
    bool ioproc;
    std::string ofname;
    Mesh const * pmesh;
    FILE * pofile;
};

//! \class WaveExtractLocal
//! \brief Extracts the l m components of the wave on a unit sphere
//! This class performs the reduction on each SphericalPatch
class WaveExtractLocal {
  public:
    //! Creates the WaveExtractLocal object
    WaveExtractLocal(SphericalGrid * psphere, MeshBlock * pmb, ParameterInput * pin);
    ~WaveExtractLocal();
    //Calculates factorial - move to a general utility?
    Real fac(Real n);
    // Calculates spin weight -2 spherical harmonics real and imaginary parts
    void swsh(Real * ylmR, Real * ylmI, int l, int m, Real theta, Real phi);
    //! Computes the l m modes of the given grid function
    void Decompose_multipole(AthenaArray<Real> const & u_R,AthenaArray<Real> const & u_I);
  public:
    //! lm projections
    Real psilmR;
    Real psilmI;
    AthenaArray<Real> psi;    
//! Patch of the spherical grid on which we are working
    SphericalPatch * ppatch;
  private:
    AthenaArray<Real> datareal;
    AthenaArray<Real> dataim;
    AthenaArray<Real> weight;
    Real rad;
    int lmax;
};

#endif
