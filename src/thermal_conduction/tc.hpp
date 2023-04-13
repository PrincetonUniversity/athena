#ifndef TC_HPP
#define TC_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file radiation.hpp
//  \brief definitions for Radiation class
//======================================================================================

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/cc/bvals_cc.hpp"
#include <string>


class MeshBlock;
class ParameterInput;
class TCIntegrator;

//! \class NewThermalConduction
//  \brief Thermal Conduction data and functions


// Array indices for  moments
enum {TCE=0, TCF1=1, TCF2=2, TCF3=3, TCRHO=4, TCT=5};

class ThermalConduction {
  friend class TCIntegrator;
  friend class BoundaryValues;
public:
  ThermalConduction(MeshBlock *pmb, ParameterInput *pin);
 // ~NewThermalConduction();
    
  AthenaArray<Real> u_tc, u_tc1, u_tc2; //thermal conduction flux
  AthenaArray<Real> coarse_tc_;

  //  three components of conduction coefficients
  AthenaArray<Real> kappa; 
  AthenaArray<Real> b_angle;

  AthenaArray<Real> flux[3]; // store transport flux, also need for refinement
  int refinement_idx{-1};


  Real vmax; // the maximum velocity (effective speed of light)
  Real min_kappa;


  MeshBlock* pmy_block;    // ptr to the parent hydro

  CellCenteredBoundaryVariable tc_bvar;

  TCIntegrator *ptcintegrator;
  
  
  //Function in problem generators to update opacity
  void EnrollOpacityFunction(TCOpacityFunc MyOpacityFunction);


  void Initialize(MeshBlock *pmb, AthenaArray<Real> &prim, 
                                  AthenaArray<Real> &u_tc);

  // The function pointer for the diffusion coefficient
  TCOpacityFunc UpdateOpacity;
 

private:

  

};

#endif // TC_HPP
