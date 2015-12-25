#ifndef VISCOSITY_HPP
#define VISCOSITY_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file viscosity.hpp 
//  \brief defines class Viscosity 
//  Contains data and functions that implement Viscosity terms
//======================================================================================

// Athena headers
#include "../../athena.hpp"         // Real
#include "../../athena_arrays.hpp"  // AthenaArray

// Declarations
class Hydro;
class ParameterInput;

class Viscosity {
public:
  Viscosity(Hydro *pf, ParameterInput *pin);
  ~Viscosity();

  void ViscosityTerms(const Real dt,
    const AthenaArray<Real> &prim, AthenaArray<Real> &cons); //update cons quantities due to viscosity
  Real VisDt(Real len, int k, int j, int i);
  Real nuiso1(int n, int k, int j, int i);
  Real cnuiso2(int k, int j, int i); 
private:
  Hydro *pmy_hydro_;  // ptr to Hydro containing this Viscosity 
  AthenaArray<Real> area,area_p1,vol; //
  AthenaArray<Real> visflx_, jvisflx_j_, kvisflx_k_; // stress tensor at the cell face 
  AthenaArray<Real> dx_,dy_,dz_; // velocity derivatives for x,y,z directions at each face 
  AthenaArray<Real> divv_; // divv
  Real nuiso_; 
};
#endif
