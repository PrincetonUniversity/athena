//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file G14Sod.cpp
//  \brief implementation of functions in Grassi 2014 Fig. 21
//======================================================================================

// this class header
#include "G14Sod.hpp"

// C header

// C++ header
#include <array>
#include <cmath>
#include <cstdio>
#include <iostream>  // endl
#include <limits>    // inf
#include <sstream>   // stringstream

// Athena++ header
#include "../../coordinates/coordinates.hpp"
#include "../../defs.hpp"
#include "../../eos/eos.hpp"
#include "../../hydro/hydro.hpp"
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"
#include "../../scalars/scalars.hpp"
#include "../utils/chemistry_utils.hpp"
#include "../utils/thermo.hpp"
#include "network.hpp"

const std::array<std::string, NSPECIES> ChemNetwork::species_names =
{"H","H+","He","He+","He2+","H-","H2","H2+"};

const std::array<std::string, ChemNetwork::ngs_> ChemNetwork::ghost_species_names_ =
{"e"};

// index of species
const int ChemNetwork::iH_ =
  ChemistryUtility::FindStrIndex(species_names.data(), NSPECIES, "H");
const int ChemNetwork::iHplus_ =
  ChemistryUtility::FindStrIndex(species_names.data(), NSPECIES, "H+");
const int ChemNetwork::iHe_ =
  ChemistryUtility::FindStrIndex(species_names.data(), NSPECIES, "He");
const int ChemNetwork::iHeplus_ =
  ChemistryUtility::FindStrIndex(species_names.data(), NSPECIES, "He+");
const int ChemNetwork::iHe2plus_ =
  ChemistryUtility::FindStrIndex(species_names.data(), NSPECIES, "He2+");
const int ChemNetwork::iHmin_ =
  ChemistryUtility::FindStrIndex(species_names.data(), NSPECIES, "H-");
const int ChemNetwork::iH2_ =
  ChemistryUtility::FindStrIndex(species_names.data(), NSPECIES, "H2");
const int ChemNetwork::iH2plus_ =
  ChemistryUtility::FindStrIndex(species_names.data(), NSPECIES, "H2+");

const int ChemNetwork::ige_ =
      ChemistryUtility::FindStrIndex(ghost_species_names_.data(), ngs_, "e") + NSPECIES;
const int ChemNetwork::igr_ = ige_;

//-------------------chemical network---------------------
//  (1) H    + e   -> H+   + 2e
//  (2) H+   + e   -> H    (+ γ)
//  (3) He   + e   -> He+  + 2e
//  (4) He+  + e   -> He   (+ γ)
//  (5) He+  + e   -> He2+ + 2e
//  (6) He2+ + e   -> He+  (+ γ)
//  (7) H    + e   -> H−   (+ γ)
//  (8) H−   + H   -> H2   + e
//  (9) H    + H+  -> H2+  (+ γ)
// (10) H    + H2+ -> H2   + H+
// (11) H+   + H2  -> H2+  + H
// (12) H2   + e   -> 2 H  + e
// (13) H−   + e   -> H    + 2 e
// (14) H−   + H   -> 2 H  + e
// (15) H−   + H+  -> 2H   (+ γ)
// (16) H−   + H+  -> H2+  + e
// (17) H2+  + e   -> 2H   (+ γ)
// (18) H2+  + H-  -> H    + H2
// (19) 2H   + H   -> H2   + H
// (20) H2   + H   -> 2H   + H

const int ChemNetwork::out2body1_[n_2body_] =
         {iHplus_,iH_,iHeplus_,iHe_,iHe2plus_,
          iHeplus_,iHmin_,iH2_,iH2plus_,iH2_,
          iH2plus_,iH2_,iH2plus_,iH_,iH_,
          iH2plus_,iH_,iH_,iH2_,iH_};

const int ChemNetwork::in2body1_[n_2body_] =
        {iH_,iHplus_,iHe_,iHeplus_,iHeplus_,
         iHe2plus_,iH_,iHmin_,iH_,iH_,
          iHplus_,iH2_,iHmin_,iHmin_,iHmin_,
          iHmin_,iH2plus_,iH2plus_,iH2_,iH2_};

// Note: output to ghost species doesn't matter. The abundances of ghost species
// are updated using the other species at every timestep
const int ChemNetwork::in2body2_[n_2body_] =
        {ige_,ige_,ige_,ige_,ige_,
         ige_,ige_,iH_,iHplus_,iH2plus_,
         iH2_,ige_,ige_,iH_,iHplus_,
        iHplus_,ige_,iHmin_,iH_,iH_};

const int ChemNetwork::out2body2_[n_2body_]
           {ige_,igr_,ige_,igr_,ige_,
            igr_,igr_,ige_,igr_,iHplus_,
            iH_,ige_,ige_,ige_,igr_,
            ige_,igr_,iH2_,iH_,iH_};

// Note output of delta (beta-alpha) : stoichiometric coefficients
const Real ChemNetwork::stoich_in2body1[n_2body_] =
                    {1, 1, 1, 1, 1,
                     1, 1, 1, 1, 1,
                     1, 1, 1, 1, 1,
                     1, 1, 1, 2, 1};

const Real ChemNetwork::stoich_in2body2[n_2body_] =
                    {1, 1, 1, 1, 1,
                     1, 1, 1, 1, 1,
                     1, 1, 1, 1, 1,
                     1, 1, 1, 1};

const Real ChemNetwork::stoich_out2body1[n_2body_] =
                    {1, 1, 1, 1, 1,
                     1, 1, 1, 1, 1,
                     1, 2, 1, 2, 2,
                     1, 2, 1, 1, 2};

const Real ChemNetwork::stoich_out2body2[n_2body_] =
                    {2, 0, 2, 0, 2,
                     0, 0, 1, 0, 1,
                     1, 1, 2, 1, 0,
                     1, 0, 1, 1, 1};

ChemNetwork::ChemNetwork(MeshBlock *pmb, ParameterInput *pin) {
  // number of species and a list of name of species
  pmy_spec_ = pmb->pscalars;
  pmy_mb_   = pmb;

  // constants
  mu_ = pin->GetReal("problem", "mu");
  muH_ = pin->GetReal("problem", "muH");
  gamma_ = pin->GetReal("hydro", "gamma");

  // initialize rates to zero
  for (int i=0; i<n_2body_; i++) {
    k2body_[i] = 0;
  }

  // copy species to a full list of species names
  for (int i=0; i<NSPECIES; i++) {
    species_names_all_[i] = species_names[i];
  }
  for (int i=NSPECIES; i<NSPECIES+ngs_; i++) {
    species_names_all_[i] = ghost_species_names_[i-NSPECIES];
  }
}

ChemNetwork::~ChemNetwork() {}

void ChemNetwork::RHS(const Real t, const Real *y, const Real ED, Real *ydot) {
  // function of evolution of the abundance of the element
  Real rate;
  // store previous y includeing ghost species
  Real *yprev = new Real[NSPECIES+ngs_];
  Real *yprev0 = new Real[NSPECIES+ngs_];
  Real *ydotg = new Real[NSPECIES+ngs_];
  Real E_ergs = ED * pmy_mb_->pmy_mesh->punit->code_energydensity_cgs / nH_;
  // Real Xe = y[iHplus_] + y[iHeplus_] + 2.*y[iHe2plus_] + y[iH2plus_] - y[iHmin_];

  // temperature, definition follows G14, T = p/rho/kb*mu, mu = 1.25
  Real T = E_ergs/Thermo::kb_*(gamma_ - 1.0)*mu_/muH_;

  // copy y to yprev and set ghost species
  GetGhostSpecies(y, yprev);

  for(int i=0; i<NSPECIES+ngs_; i++) {
    ydotg[i] = 0.0;
    if (yprev[i] < 0) {
      yprev0[i] = 0;
    } else {
      yprev0[i] = yprev[i];
    }
    // throw error if nan, or inf, or large negative value occurs
    if ( std::isnan(yprev[i]) || std::isinf(yprev[i]) ) {
      printf("RHS: ");
      for (int j=0; j<NSPECIES+ngs_; j++) {
        printf("%s: %.2e  ", species_names_all_[j].c_str(), yprev[j]);
      }
      printf("\n");
      OutputRates(stdout);
      printf("\n");
      printf("nH_ = %.2e, E_ergs = %.2e, T= %.2e \n", nH_,E_ergs, T);
      std::stringstream msg;
      msg << "ChemNetwork (G14): RHS(yprev): nan or inf" << std::endl;
      ATHENA_ERROR(msg);
    }
  }

  // update the rate
  UpdateRates(yprev0, E_ergs);

  // 2body reactions
  for (int i=0; i<n_2body_; i++) {
    rate =  k2body_[i] * yprev0[in2body1_[i]] * yprev0[in2body2_[i]];
    if (yprev0[in2body1_[i]] < 0 && yprev0[in2body2_[i]] < 0) {
      rate *= -1.;
    }
    if (DEBUG) {
      if (species_names_all_[in2body1_[i]] == "H2"
          || species_names_all_[in2body2_[i]]== "H2") {
        printf("Dec: k2body_[%.2i]*(X_%s)*(X_%s) = %.2e* %.2e * %.2e = %.2e \n",
               i, species_names_all_[in2body1_[i]].c_str(),
               species_names_all_[in2body2_[i]].c_str(),
               k2body_[i], yprev0[in2body1_[i]], yprev0[in2body2_[i]],rate);
      }
      if (species_names_all_[out2body1_[i]] == "H2"
          || species_names_all_[out2body2_[i]]== "H2") {
        printf("Inc: k2body_[%.2i]*(X_%s)*(X_%s) = %.2e* %.2e * %.2e = %.2e \n",
               i, species_names_all_[in2body1_[i]].c_str(),
               species_names_all_[in2body2_[i]].c_str(),
               k2body_[i] ,yprev0[in2body1_[i]], yprev0[in2body2_[i]],rate);
      }
    }
    ydotg[in2body1_[i]]  -= stoich_in2body1[i]*rate;
    ydotg[in2body2_[i]]  -= stoich_in2body2[i]*rate;
    ydotg[out2body1_[i]] += stoich_out2body1[i]*rate;
    ydotg[out2body2_[i]] += stoich_out2body2[i]*rate;
  }
  // set ydot to return
  for (int i=0; i<NSPECIES; i++) {
    // return in code units
    ydot[i] = ydotg[i] * pmy_mb_->pmy_mesh->punit->code_time_cgs * nH_;
  }
  if (DEBUG) {
    OutputRates(stdout);
    printf("nH_ = %2.e \n",nH_);
    printf("one loop finish, ydot, ydotg yprev = \n");
    for (int j=0; j<NSPECIES+ngs_; j++) {
      printf("%s: %.2e , %.2e, %.2e ", species_names_all_[j].c_str(),
             ydot[j],ydotg[j],yprev[j]);
      printf("\n");
    }
  }
  delete[] yprev;
  delete[] yprev0;
  delete[] ydotg;
  return;
}

void ChemNetwork::GetGhostSpecies(const Real *y, Real *yghost) {
  // copy the aboundances in y to yghost
  for (int i=0; i<NSPECIES; i++) {
    yghost[i] = y[i];
  }
  // set the ghost species
  yghost[ige_] = y[iHplus_] + y[iHeplus_] + 2.*y[iHe2plus_] + y[iH2plus_] - y[iHmin_];
  return;
}

void ChemNetwork::UpdateRates(const Real *y, const Real E) {
  // temperature, definition follows G14, T = p/rho/kb*mu, mu = 1.25
  Real T = E/Thermo::kb_*(gamma_ - 1.0)*mu_/muH_;

  // const Real logT    = std::log10(T);
  const Real lnTe    = std::log(T* 8.6163e-5);
  const Real Te      = T* 8.6163e-5;

  // (1) H    + e   -> H+   + 2e
  k2body_[0]  =  std::exp(-3.271396786e1 +
                    ( 1.35365560e1  + (-5.73932875e0 + (1.56315498e0   +
                    (-2.877056e-1   + (3.48255977e-2 + (-2.63197617e-3 +
                    ( 1.11954395e-4 + (-2.03914985e-6)
                    *lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)
                    *lnTe);
  //  (2) H+   + e   -> H    (+ γ)
  if (T<=5500) {
    k2body_[1]  =  3.92e-13*std::pow(Te,-0.6353);
  } else {
    k2body_[1]  =  std::exp(-2.861303380689232e1 +
                (-7.241125657826851e-1 + (-2.026044731984691e-2 + (-2.380861877349834e-3 +
                (-3.212605213188796e-4 + (-1.421502914054107e-5 + ( 4.989108920299513e-6 +
                ( 5.755614137575758e-7 + (-1.856767039775261e-8 + (-3.071135243196595e-9 )
                *lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)
                *lnTe)*lnTe);
  }
  //  (3) He   + e   -> He+  + 2e [unit in cgs]
  k2body_[2]  = std::exp(-4.409864886e1 +
                   ( 2.391596563e1 + (-1.07532302e1  + (+3.05803875e0 +
                   (-5.6851189e-1  + ( 6.79539123e-2 + (-5.0090561e-3 +
                   ( 2.06723616e-4 + (-3.64916141e-6)
                    *lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe);

  //  (4) He+  + e   -> He   (+ γ)  k_4 = k_4a + k_4b unit in cgs
  k2body_[3]  = 3.92e-13*std::pow(Te,-0.6353);
  if (T>9280) {
    k2body_[3] += 1.54e-9*std::pow(Te,-1.5)
                  *(1.0+0.3/std::exp(8.099328789667/Te))/std::exp(40.49664394833662/Te);
  }

  //  (5) He+  + e   -> He2+ + 2e  unit in cgs
  k2body_[4]  = std::exp(-6.8710409e1  +
                  ( 4.393347632635e1   + (-1.848066993568e1  + (+4.701626486759e0  +
                  (-7.692466334492e-1  + (8.113042097303e-2  + (-5.324020628287e-3 +
                  (+1.975705312221e-4  + (-3.165581065665e-6)
                  *lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe);

  //  (6) He2+ + e   -> He+  (+ γ) unit in cgs
  k2body_[5]  = 3.36e-10*std::pow(T,-0.5)*std::pow(T/1e3,-0.2)
                *std::pow(1+std::pow(T/1e6,0.7),-1);

  //  (7) H    + e   -> H−   (+ γ)
  k2body_[6]  = 6.77e-15*std::pow(Te,0.8779);

  //  (8) H−   + H   -> H2   + e unit in cgs
  if (T<=1160) {
    k2body_[7]  = 1.43e-9;
  } else {
    k2body_[7]  = std::exp(-2.00691389e1  +
             ( 2.289800603272916e-1  + ( 3.599837721023835e-2  + (-4.555120027032095e-3 +
             (-3.105115447124016e-4  + ( 1.073294010367247e-4  + (-8.36671960467864e-6  +
             ( 2.238306228891639e-7)
             *lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe);
  }

  //  (9) H    + H+  -> H2+  (+ γ)
  if (T<=6700) {
    k2body_[8]  = 1.85e-23*std::pow(T,1.8);
  } else {
    k2body_[8]  = 5.81e-16*std::pow((T/56200), -0.6657*std::log10(T/56200));
  }
  // (10) H    + H2+ -> H2   + H+ unit in cgs
  k2body_[9]  = 6e-10;

  // (11) H+   + H2  -> H2+  + H
  if (T<=3480) {
    k2body_[10] = 0;
  } else {
    k2body_[10] = std::exp(-2.424914687731536e1  +
              ( 3.400824447095291e0   + (-3.898003964650152e0  + ( 2.045587822403071e0  +
              (-5.416182856220388e-1  + ( 8.41077503763412e-2  + (-7.879026154483455e-3 +
              ( 4.138398421504563e-4  + (-9.36345888928611e-6)
              *lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe);
  }

  // (12) H2   + e   -> 2 H  + e
  k2body_[11] = 5.6e-11*std::sqrt(T)*std::exp(-102124.0/T);

  // (13) H−   + e   -> H    + 2 e
  k2body_[12] = std::exp(-1.801849334273e1  +
                   ( 2.360852208681e0     + (-2.827443061704e-1 + ( 1.623316639567-2  +
                   (-3.365012031362999e-2 + ( 1.178329782711e-2 + (-1.656194699504e-3 +
                   ( 1.068275202678e-4    + (-2.631285809207e-6)
                   *lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe);

  // (14) H−   + H   -> 2 H  + e
  if (T<=1160) {
    k2body_[13]  = 2.56e-9*std::pow(Te,1.78186);
  } else {
    k2body_[13]  =  std::exp(-2.037260896203573e1 +
                    (+1.139449335841631e0 + (-1.421013521554148e-1 + ( 8.46445538663e-3 +
                    (-1.4327641212992e-3  + ( 2.012250284791e-4    + (+8.66396324309e-5 +
                    (-2.585009680264e-5   + (+2.4555011970392e-6   + (-8.06838246118e-8 )
                    *lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)
                    *lnTe)*lnTe);
  }
  // (15) H−   + H+  -> 2H  (+ γ)
  k2body_[14] = 6.5e-9/std::sqrt(Te);

  // (16) H−   + H+  -> H2+  + e
  k2body_[15] = 1e-8*std::pow(T,-0.4);

  // (17) H2+  + e   -> 2H  (+ γ)
  if (T<=617) {
    k2body_[16] = 1e-8;
  } else {
    k2body_[16] = 1.32e-6*std::pow(T,-0.76);
  }

  // (18) H2+  + H-   -> H    + H2
  k2body_[17]  = 5e-7*std::sqrt(1.0e2/T);

  // (19) 3H → H2 + H
  if (T<=300) {
    k2body_[18] = 1.3e-32*std::pow(T/300,-0.38);
  } else {
    k2body_[18] = 1.3e-32*std::pow(T/300,-1.0);
  }

  // (20) H2   + H   -> 3 H
  k2body_[19] = 1.067e-10*std::pow(Te,2.012)*std::exp(-4.463/Te)
                *std::pow(1. + 0.2472*Te,-3.512);
  return;
}


void ChemNetwork::InitializeNextStep(const int k, const int j, const int i) {
  // density
  nH_ = pmy_mb_->phydro->w(IDN, k, j, i);
  return;
}


Real ChemNetwork::Edot(const Real t, const Real *y, const Real ED) {
  // function of evolution of energy
  // return dEdt;
  Real LH2; // Define H2 Cooling Term
  Real E_ergs = ED * pmy_mb_->pmy_mesh->punit->code_energydensity_cgs / nH_;

  // temperature, definition follows G14, T = p/rho/kb*mu, mu = 1.25
  Real T = E_ergs/Thermo::kb_*(gamma_ - 1.0)*mu_/muH_;

  // Real dEdt = 0.;
  Real *yprev = new Real[NSPECIES+ngs_];
  // copy y to yprev and set ghost species
  GetGhostSpecies(y, yprev);
  // correct negative abundance to zero
  for (int i=0; i<NSPECIES+ngs_; i++) {
    if (yprev[i] < 0) {
      yprev[i] = 0;
    }
  }

  // H2 CE(ro-vibrational) Cooling
  LH2  = Thermo::CoolingH2(yprev[iH2_], nH_*yprev[iH_],  nH_*yprev[iH2_],
          nH_*yprev[iHe_],  nH_*yprev[iHplus_], nH_*yprev[ige_],
          T);

  // return in code units
  Real dEDdt = - LH2* nH_ / pmy_mb_->pmy_mesh->punit->code_energydensity_cgs
               * pmy_mb_->pmy_mesh->punit->code_time_cgs;
  if (DEBUG) {
    printf("T=%.4e, dEDdt=%.2e, E=%.2e, E_ergs=%.2e, nH=%.2e\n",
           T, dEDdt, ED, E_ergs, nH_);
    for (int i=0; i<NSPECIES+ngs_; i++) {
      printf("%s: %.2e  ", species_names_all_[i].c_str(), yprev[i]);
    }
    printf("\n");
    printf("=============================\n");
  }
  delete[] yprev;
  return dEDdt;
}

void ChemNetwork::OutputRates(FILE *pf) const {
  // output the reactions and base rates
  for (int i=0; i<n_2body_; i++) {
    fprintf(pf, "%4s  + %4s -> %4s  + %4s,     k2body = %.2e\n",
     species_names_all_[in2body1_[i]].c_str(),
     species_names_all_[in2body2_[i]].c_str(),
     species_names_all_[out2body1_[i]].c_str(),
     species_names_all_[out2body2_[i]].c_str(),
     k2body_[i]);
  }
  return;
}
