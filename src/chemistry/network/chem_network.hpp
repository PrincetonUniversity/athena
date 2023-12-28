#ifndef CHEMISTRY_NETWORK_CHEM_NETWORK_HPP_
#define CHEMISTRY_NETWORK_CHEM_NETWORK_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file chem_network.hpp
//! \brief definitions for empty ChemNetwork when chemistry is not included.

// C headers

// C++ headers
#include <array>
#include <string>

// Athena++ classes headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "network.hpp"

class ChemNetwork : public NetworkWrapper {
  friend class MeshBlock;
 public:
  ChemNetwork(MeshBlock *pmb, ParameterInput *pin) {}
  ~ChemNetwork() {}

  // species names required by output
  //static const std::string species_names[NSPECIES];
  static const std::array<std::string, NSPECIES> species_names;

  void InitializeNextStep(const int k, const int j, const int i) {}
  void RHS(const Real t, const Real *y, const Real ED, Real *ydot) {}
  Real Edot(const Real t, const Real *y, const Real ED) {return 0.0;}

 private:
  // cr indexing required by problem generatorx chem_*.cpp, not used.
  static constexpr int index_cr_ = -1;
};

#endif // CHEMISTRY_NETWORK_CHEM_NETWORK_HPP_
