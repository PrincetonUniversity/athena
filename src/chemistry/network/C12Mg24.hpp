#ifndef CHEMISTRY_NETWORK_C12MG24_HPP_
#define CHEMISTRY_NETWORK_C12MG24_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file H2.hpp
//! \brief definitions for a very simple chemical network to test nuclear
//! reations. Created by Goni Halevi.

//c++ headers
#include <string> //std::string

// Athena++ classes headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "network.hpp"

//! \class ChemNetwork
//!  \brief Chemical Network for nuclear reaction test
class ChemNetwork : public NetworkWrapper {
  friend class MeshBlock;
 public:
  ChemNetwork(MeshBlock *pmb, ParameterInput *pin);
  ~ChemNetwork();

  //a list of species name, used in output
  static const std::string species_names[NSCALARS];

  void InitializeNextStep(const int k, const int j, const int i);

  void RHS(const Real t, const Real y[NSCALARS], const Real ED,
           Real ydot[NSCALARS]);

  Real Edot(const Real t, const Real y[NSCALARS], const Real ED);
 private:
  PassiveScalars *pmy_spec_;
  MeshBlock *pmy_mb_;

  std::string species_names_all_[NSCALARS];//all species
  //index of species
  static const int iC12_;
  static const int iMg24_;
  Real Q12_;
  Real mn_;
  Real k_;
  Real unit_density_in_cgs;
  Real unit_temp_in_K;
  Real unit_edot_in_cgs;
  Real density_;
  Real mu_;
  Real rate_C12;
  void OutputRates(FILE *pf) const;
};

#endif // CHEMISTRY_NETWORK_C12MG24_HPP_
