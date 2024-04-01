//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file H2.cpp
//! \brief implementation of functions in class ChemNetwork, using the simple
//! network for H2 formation and destruction.


// this class header
#include "H2.hpp"

// C headers

// C++ header
#include <iostream>   // endl
#include <limits>     // inf
#include <sstream>    // stringstream

// Athena++ header
#include "../../defs.hpp"
#include "../../eos/eos.hpp"
#include "../../hydro/hydro.hpp"
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"
#include "../../scalars/scalars.hpp"
#include "../../units/units.hpp"
#include "../utils/chemistry_utils.hpp"
#include "../utils/thermo.hpp"
#include "network.hpp"

// constants
const Real ChemNetwork::kgr_ = 3e-17;

// species names
const std::array<std::string, NSPECIES> ChemNetwork::species_names = {"H", "H2"};

const int ChemNetwork::iH_ =
  ChemistryUtility::FindStrIndex(species_names.data(), NSPECIES, "H");
const int ChemNetwork::iH2_ =
  ChemistryUtility::FindStrIndex(species_names.data(), NSPECIES, "H2");

// flag for Cv
static bool is_const_Cv;


//----------------------------------------------------------------------------------------
//! \brief ChemNetwork constructor

ChemNetwork::ChemNetwork(MeshBlock *pmb, ParameterInput *pin) {
  // number of species and a list of name of species
  pmy_spec_ = pmb->pscalars;
  pmy_mb_ = pmb;

  // set the parameters from input file
  xi_cr_ = pin->GetOrAddReal("chemistry", "xi_cr", 2e-16);
  kcr_ = xi_cr_ * 3.;
  // set Cv: constant or H2 abundance dependent
  is_const_Cv = pin->GetOrAddBoolean("problem", "is_const_Cv", true);
}

//----------------------------------------------------------------------------------------
//! \brief ChemNetwork destructor

ChemNetwork::~ChemNetwork() {
}

//----------------------------------------------------------------------------------------
//! \fn void ChemNetwork::InitializeNextStep(const int k, const int j, const int i)
//! \brief Set the rates of chemical reactions, eg. through density and radiation field.
//!
//! k, j, i are the corresponding index of the grid

void ChemNetwork::InitializeNextStep(const int k, const int j, const int i) {
  Real rho, rho_floor;
  // density
  rho = pmy_mb_->phydro->w(IDN, k, j, i);
  // apply density floor
  rho_floor = pmy_mb_->peos->GetDensityFloor();
  rho = (rho > rho_floor) ?  rho : rho_floor;
  // hydrogen atom number density
  nH_ =  rho;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ChemNetwork::RHS(const Real t, const Real *y, const Real ED,
//!                       Real *ydot)
//! \brief RHS: right-hand-side of ODE.
//!
//! dy/dt = ydot(t, y). Here y are the abundance
//! of species. details see CVODE package documentation.
//! all input/output variables are in code units

void ChemNetwork::RHS(const Real t, const Real *y, const Real ED, Real *ydot) {
  const Real rate_cr = kcr_ * y[iH2_];
  const Real rate_gr = kgr_ * nH_ * y[iH_];
  ydot[iH2_] = rate_gr - rate_cr;
  ydot[iH_] = -2*rate_gr + 2*rate_cr;
  for (int i=0; i<NSPECIES; i++) {
    // return in code units
    ydot[i] *= pmy_mb_->pmy_mesh->punit->code_time_cgs;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Real ChemNetwork::Edot(const Real t, const Real *y, const Real ED)
//! \brief energy equation dED/dt
//!
//! all input/output variables are in code units (ED is the energy density)

Real ChemNetwork::Edot(const Real t, const Real *y, const Real ED) {
  // isothermal
  if (!NON_BAROTROPIC_EOS) {
    return 0;
  }
  const Real x_He = 0.1;
  const Real x_e = 0.;
  Real x_H2;
  if (is_const_Cv) {
    x_H2 = 0;
  } else {
    x_H2 = y[iH2_];
  }
  const Real T_floor = 1.; // temperature floor for cooling
  // ernergy per hydrogen atom
  const Real E_ergs = ED * pmy_mb_->pmy_mesh->punit->code_energydensity_cgs / nH_;
  Real T = E_ergs / Thermo::CvCold(x_H2, x_He, x_e);
  if (T < T_floor) {
    return 0;
  }
  Real dEdt = - Thermo::alpha_GD_ * nH_ * std::sqrt(T) * T;
  // return in code units
  Real dEDdt = (dEdt * nH_ / pmy_mb_->pmy_mesh->punit->code_energydensity_cgs)
                * pmy_mb_->pmy_mesh->punit->code_time_cgs;
  return dEDdt;
}
