#ifndef CHEMISTRY_NETWORK_G14SOD_HPP_
#define CHEMISTRY_NETWORK_G14SOD_HPP_
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file G14Sod.hpp
//  \brief implementation of functions in Grassi 2014 Fig. 21
//======================================================================================

// C headers

// C++ headers
#include <array>
#include <string>

// Athena++ classes headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "network.hpp"

class ChemNetwork : public NetworkWrapper {
  // OutputProperties in problem generator called by Mesh::UserWorkAfterLoop.
  friend class Mesh;
  // It would be convenient to know the species names in
  // initialization of chemical species in problem
  friend class MeshBlock;
 public:
  static constexpr int ngs_     = 1;

  ChemNetwork(MeshBlock *pmb, ParameterInput *pin);
  ~ChemNetwork();

  // a list of species name, used in output
  static const std::array<std::string, NSCALARS> species_names;

  // Set the rates of chemical reactions, eg. through density.
  // k, j, i are the corresponding index of the grid
  void InitializeNextStep(const int k, const int j, const int i);

  // RHS: right-hand-side of ODE. dy/dt = ydot(t, y). Here y are the abundance
  // of species. details see CVODE package documentation.
  // all input/output variables are in code units
  void RHS(const Real t, const Real *y, const Real ED, Real *ydot);

  // energy equation dE/dt, all input/output variables are in code units
  // (ED is the energy density)
  Real Edot(const Real t, const Real *y, const Real ED);
  void GetGhostSpecies(const Real *y, Real *yghost);
  void UpdateRates(const Real *y, const Real E);
  void OutputRates(FILE *pf) const;

 private:
  PassiveScalars *pmy_spec_;
  MeshBlock *pmy_mb_;

  // constants
  static constexpr int n_2body_ = 20;
  static const std::array<std::string, ngs_> ghost_species_names_;
  Real mu_; // mean molecular weight, assume to be constant
  Real muH_; // convert mass density to number density of H neuclei
  Real gamma_; // adiabatic index
  std::array<std::string, NSCALARS+ngs_> species_names_all_;

  Real nH_;   // hydrogen nuclei number denisty, updated at InitializeNextStep

  // index of species
  static const int iH_;
  static const int iHplus_;
  static const int iHe_;
  static const int iHeplus_;
  static const int iHe2plus_;
  static const int iHmin_;
  static const int iH2_;
  static const int iH2plus_;
  static const int ige_;
  static const int igr_;
  // -------------------chemical network---------------------
  // 2body reactions
  static const int in2body1_[n_2body_];
  static const int in2body2_[n_2body_];
  static const int out2body1_[n_2body_];
  static const int out2body2_[n_2body_];
  static const Real stoich_in2body1[n_2body_];
  static const Real stoich_in2body2[n_2body_];
  static const Real stoich_out2body1[n_2body_];
  static const Real stoich_out2body2[n_2body_];
  std::array<Real, n_2body_> k2body_;
  // TODO(KGF): cpplint false positive VLA detection?
  // Real k2body_[n_2body_]; // rates for 2 body reacrtions.
};

#endif // CHEMISTRY_NETWORK_G14SOD_HPP_
