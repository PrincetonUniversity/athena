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
//! \file gow17.cpp
//  \brief implementation of functions in class ChemNetwork, using the GOW17
//  network, see paper by Gong, Ostriker, Wolfire 2017 
//======================================================================================

// this class header
#include "gow17.hpp"

//athena++ header
#include "network.hpp"
#include "../../scalars/scalars.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../parameter_input.hpp"       //ParameterInput
#include "../../mesh/mesh.hpp"
#include "../../hydro/hydro.hpp"
#include "../../radiation/radiation.hpp"
#include "../../radiation/integrators/rad_integrators.hpp"
#include "../utils/chemistry_utils.hpp"
#include "../utils/thermo.hpp"
#include "../../defs.hpp"
#include "../../eos/eos.hpp"

//c++ header
#include <sstream>    // stringstream
#include <iostream>   // endl
#include <math.h> //a^x = pow(a,x)
#include <stdio.h> //FILE, fprintf()
#include <limits> //inf

#ifdef DEBUG
static bool output_flag = true;
#endif

//constants
const Real ChemNetwork::temp_coll_ = 7.0e2;
//small number
const Real ChemNetwork::small_ = 1e-50;

//species names
const std::string ChemNetwork::species_names[NSCALARS] = 
{"He+", "OHx", "CHx", "CO", "C+", "HCO+", "H2", "H+", "H3+", "H2+", "O+", "Si+"};

//below are ghost species. The aboundances of ghost species are
// recalculated in RHS everytime by other species.
const std::string ChemNetwork::ghost_species_names_[ngs_] = 
{"*Si", "*C", "*O", "*He", "*e", "*H"};

//index of species
const int ChemNetwork::iHeplus_ =
	ChemistryUtility::FindStrIndex(species_names, NSCALARS, "He+");
const int ChemNetwork::iOHx_ = 
	ChemistryUtility::FindStrIndex(species_names, NSCALARS, "OHx");
const int ChemNetwork::iCHx_ = 
	ChemistryUtility::FindStrIndex(species_names, NSCALARS, "CHx");
const int ChemNetwork::iCO_ =
	ChemistryUtility::FindStrIndex(species_names, NSCALARS, "CO");
const int ChemNetwork::iCplus_ =
  ChemistryUtility::FindStrIndex(species_names, NSCALARS, "C+");
const int ChemNetwork::iHCOplus_ =
  ChemistryUtility::FindStrIndex(species_names, NSCALARS, "HCO+");
const int ChemNetwork::iH2_ =
  ChemistryUtility::FindStrIndex(species_names, NSCALARS, "H2");
const int ChemNetwork::iHplus_ =
  ChemistryUtility::FindStrIndex(species_names, NSCALARS, "H+");
const int ChemNetwork::iH3plus_ =
  ChemistryUtility::FindStrIndex(species_names, NSCALARS, "H3+");
const int ChemNetwork::iH2plus_ =
  ChemistryUtility::FindStrIndex(species_names, NSCALARS, "H2+");
const int ChemNetwork::iOplus_ =
  ChemistryUtility::FindStrIndex(species_names, NSCALARS, "O+");
const int ChemNetwork::iSiplus_ =
  ChemistryUtility::FindStrIndex(species_names, NSCALARS, "Si+");
//index of ghost species
const int ChemNetwork::igSi_ =
  ChemistryUtility::FindStrIndex(ghost_species_names_, ngs_, "*Si") + NSCALARS;
const int ChemNetwork::igC_ =
  ChemistryUtility::FindStrIndex(ghost_species_names_, ngs_, "*C") + NSCALARS;
const int ChemNetwork::igO_ =
  ChemistryUtility::FindStrIndex(ghost_species_names_, ngs_, "*O") + NSCALARS;
const int ChemNetwork::igHe_ =
  ChemistryUtility::FindStrIndex(ghost_species_names_, ngs_, "*He") + NSCALARS;
const int ChemNetwork::ige_ =
  ChemistryUtility::FindStrIndex(ghost_species_names_, ngs_, "*e") + NSCALARS;
const int ChemNetwork::igH_ =
  ChemistryUtility::FindStrIndex(ghost_species_names_, ngs_, "*H") + NSCALARS;


//-------------------chemical network---------------------
//cosmic ray chemistry network
// (0) cr + H2 -> H2+ + *e
// (1) cr + *He -> He+ + *e 
// (2) cr + *H  -> H+ + *e 
// -----added as Clark + Glover 2015---- 
// (3) cr + *C -> C+ + *e     --including direct and cr induce photo reactions 
// (4) crphoto + CO -> *O + *C       
// (5) cr + CO -> HCO+ + e  --schematic for cr + CO -> CO+ + e
// -----Si, CR induced photo ionization, experimenting----
// (6) cr + Si -> Si+ + e, UMIST12
const int ChemNetwork::icr_H2_ = 0;
const int ChemNetwork::icr_He_ = 1;
const int ChemNetwork::icr_H_ = 2;
const int ChemNetwork::incr_[n_cr_] = 
												 {iH2_, igHe_, igH_, 
													igC_, iCO_, iCO_,
													igSi_};
const int ChemNetwork::outcr_[n_cr_] =
												 {iH2plus_, iHeplus_, iHplus_, 
													iCplus_, igO_, iHCOplus_,
													iSiplus_};
const Real ChemNetwork::kcr_base_[n_cr_] = 
												 {2.0, 1.1, 1.0, 
													560., 90., 6.52,
													8400.}; 

//2 body reactions
//NOTE: photons from recombination are ignored
// Reactions are, in order.
//  -- are equations of special rate treatment in Glover, Federrath+ 2010:
// (0) H3+ + *C -> CH + H2         --Vissapragada2016 new rates
// (1) H3+ + *O -> OH + H2        
// (2) H3+ + CO -> HCO+ + H2
// (3) He+ + H2 -> H+ + *He + *H    --fit to Schauer1989
// (4) He+ + CO -> C+ + *O + *He   
// (5) C+ + H2 -> CH + *H         -- schematic reaction for C+ + H2 -> CH2+
// (6) C+ + OH -> HCO+             -- Schematic equation for C+ + OH -> CO+ + H.
// Use rates in KIDA website.
// (7) CH + *O -> CO + *H
// (8) OH + *C -> CO + *H          --exp(0.108/T)
// (9) He+ + *e -> *He             --(17) Case B
// (10) H3+ + *e -> H2 + *H
// (11) C+ + *e -> *C              -- Include RR and DR, Badnell2003, 2006.
// (12) HCO+ + *e -> CO + *H
// ----added in GO2012--------
// (13) H2+ + H2 -> H3+ + *H       --(54) exp(-T/46600)
// (14) H+ + *e -> *H              --(12) Case B
// ---collisional dissociation, only important at high temperature T>1e3---
// (15) H2 + *H -> 3 *H            --(9) Density dependent. See Glover+MacLow2007
// (16) H2 + H2 -> H2 + 2 *H       --(10) Density dependent. See Glover+MacLow2007
// (17) *H + *e -> H+ + 2 *e       --(11) Relates to Te
// ----added for H3+ destruction in addition to (10)----
// (18) H3+ + *e -> *3H            --(111)
// ----added He+ destruction in addtion to (3), from UMIST12----
// (19) He+ + H2 -> H2+ + *He
// ----added CH reaction to match for abundances of CH---
// (20) CH + *H -> H2 + *C         
// ----added to match the Meudon code ---
// (21) OH + *O -> *O + *O + *H
// ---branching of C+ + H2 ------
// (22) C+ + H2 + *e -> *C + *H + *H
// ---Si , rate from UMIST12---
// (23) Si+ + *e -> *Si
// --- H2O+ + e reaction ---
// (24) H3+ + *O + *e -> H2 + *O + *H
// --- OH destruction with He+
// (25) He+ + OH -> O+ + *He + *H 
// --- H2+ charge exchange with H ---
// (26) H2+ + *H -> H+ + H2 
//  --- O+ reactions ---
// (27) H+ + *O -> O+ + *H -- exp(-232/T)
// (28) O+ + *H -> H+ + *O 
// (29) O+ + H2 -> OH + *H     -- branching of H2O+
// (30) O+ + H2 -> *O + *H + *H  -- branching of H2O+

const int ChemNetwork::i2body_H2_H = 15;
const int ChemNetwork::i2body_H2_H2 = 16;
const int ChemNetwork::i2body_H_e = 17;
const int ChemNetwork::in2body1_[n_2body_] = 
          {iH3plus_, iH3plus_, iH3plus_, iHeplus_, iHeplus_,    
           iCplus_, iCplus_, iCHx_, iOHx_, iHeplus_,
           iH3plus_, iCplus_, iHCOplus_, iH2plus_, iHplus_,
           iH2_, iH2_, igH_, iH3plus_, iHeplus_, 
           iCHx_, iOHx_, iCplus_, iSiplus_, iH3plus_,
           iHeplus_, iH2plus_, iHplus_, iOplus_, iOplus_,
           iOplus_};
const int ChemNetwork::in2body2_[n_2body_] = 
          {igC_, igO_, iCO_, iH2_, iCO_,   
           iH2_, iOHx_, igO_, igC_, ige_,   
           ige_, ige_, ige_, iH2_, ige_,
           igH_, iH2_, ige_, ige_, iH2_, 
           igH_, igO_, iH2_, ige_, igO_,
           iOHx_, igH_, igO_, igH_, iH2_,
           iH2_};
//Note: output to ghost species doesn't matter. The abundances of ghost species
// are updated using the other species at every timestep
const int ChemNetwork::out2body1_[n_2body_] = 
          {iCHx_, iOHx_, iHCOplus_, iHplus_, iCplus_,   
           iCHx_, iHCOplus_, iCO_, iCO_, igHe_,   
           iH2_, igC_, iCO_, iH3plus_, igH_,
           igH_, iH2_, iHplus_, igH_, iH2plus_, 
           iH2_, igO_, igC_, igSi_, iH2_,
           iOplus_, iHplus_, iOplus_, iHplus_, iOHx_,
           igO_};
const int ChemNetwork::out2body2_[n_2body_] = 
          {iH2_, iH2_, iH2_, igHe_, igO_,   
           igH_, igH_, igH_, igH_, igH_,   
           igH_, igH_, igH_, igH_, igH_,
           igH_, igH_, ige_, igH_, igHe_, 
           igC_, igH_, igH_, igH_, igO_,
           igHe_, iH2_, igH_, igO_, igH_,
           igH_};
const Real ChemNetwork::k2Texp_[n_2body_] = 
 {0.0, -0.190, 0.0, 0.0, 0.0, 
  -1.3, 0.0, 0.0, -0.339, -0.5, 
  -0.52, 0.0, -0.64, 0.042, 0.0,
  0.0, 0.0, 0.0, -0.52, 0.0,
  0.26, 0.0, -1.3, -0.62, -0.190,
  0.0, 0.0, 0.0, 0.0, 0.0,
  0.0};
const Real ChemNetwork::k2body_base_[n_2body_] = 
                {1.00, 1.99e-9, 1.7e-9, 1.26e-13, 1.6e-9, 
                 3.3e-13 * 0.7, 1.00, 7.0e-11, 7.95e-10, 1.0e-11, 
                 4.54e-7, 1.00, 1.06e-5, 1.76e-9, 2.753e-14,
                 1.00, 1.00, 1.00, 8.46e-7, 7.20e-15, 
                 2.81e-11, 3.5e-11, 3.3e-13 * 0.3, 1.46e-10, 1.99e-9,
                 1.00, 6.4e-10, 1.00, 1.00, 1.6e-9,
                 1.6e-9};
//rates for H3+ + C forming CH+ and CH2+
const Real ChemNetwork::A_kCHx_ = 1.04e-9;
const Real ChemNetwork::n_kCHx_ = 2.31e-3;
const Real ChemNetwork::c_kCHx_[4] = {3.4e-8, 6.97e-9, 1.31e-7, 1.51e-4};
const Real ChemNetwork::Ti_kCHx_[4] = {7.62, 1.38, 2.66e1, 8.11e3};

// photo reactions.
// Reaction rates in Drain 1978 field units.
// Reactions are, in order:
// (0) h nu + *C -> C+ + *e
// (1) h nu + CH -> *C + *H
// (2) h nu + CO -> *C + *O            --self-shielding and shielding by H2
// (3) h nu + OH -> *O + *H
// ----added in GO2012--------
// (4) h nu + H2 -> *H + *H            --self- and dust shielding
// ----Si, from UMIST12
// (5) h nu + *Si -> Si+
const int ChemNetwork::iph_C_ = 0;
const int ChemNetwork::iph_CHx_ = 1;
const int ChemNetwork::iph_CO_ = 2;
const int ChemNetwork::iph_OHx_ = 3;
const int ChemNetwork::iph_H2_ = 4;
const int ChemNetwork::iph_Si_ = 5;
const int ChemNetwork::inph_[n_ph_] = {
              igC_, iCHx_, iCO_,
              iOHx_, iH2_, igSi_};
const int ChemNetwork::outph1_[n_ph_] = {
              iCplus_, igC_, igC_,
              igO_, igH_, iSiplus_};
const Real ChemNetwork::kph_base_[n_ph_] = {3.5e-10, 9.1e-10, 2.4e-10, 
																			      3.8e-10, 5.7e-11, 4.5e-9}; 
const Real ChemNetwork::kph_avfac_[n_ph_] = {3.76, 2.12, 3.88,
	                                           2.66, 4.18, 2.61};

// Grain assisted recombination of H, H2, C+ and H+
// (0) *H + *H + gr -> H2 + gr
// (1) H+ + *e + gr -> *H + gr
// (2) C+ + *e + gr -> *C + gr
// (3) He+ + *e + gr -> *He + gr
// ------Si, from WD2001-----
// (4) Si+ + *e + gr -> *Si + gr
const int ChemNetwork::igr_H_ = 0;
const int ChemNetwork::ingr_[n_gr_] = {igH_, iHplus_, iCplus_, iHeplus_, 
                                       iSiplus_};
const int ChemNetwork::outgr_[n_gr_] = {iH2_, igH_, igC_, igHe_, 
                                        igSi_};
const Real ChemNetwork::cHp_[7] = {12.25, 8.074e-6, 1.378, 5.087e2, 1.586e-2,
 											             0.4723, 1.102e-5}; 
const Real ChemNetwork::cCp_[7] = {45.58, 6.089e-3, 1.128, 4.331e2, 4.845e-2,
                                   0.8120, 1.333e-4};
const Real ChemNetwork::cHep_[7] = {5.572, 3.185e-7, 1.512, 5.115e3, 3.903e-7,
                                    0.4956, 5.494e-7};
const Real ChemNetwork::cSip_[7] = {2.166, 5.678e-8, 1.874, 4.375e4, 1.635e-6,
                                    0.8964, 7.538e-5};
//-----------------end of chemical network---------------------


ChemNetwork::ChemNetwork(MeshBlock *pmb, ParameterInput *pin) {
	//number of species and a list of name of species
  pmy_spec_ = pmb->pscalars;
	pmy_mb_ = pmb;

	//set the parameters from input file
	zdg_ = pin->GetOrAddReal("chemistry", "Zdg", 1.);//dust and gas metallicity
	xHe_ = pin->GetOrAddReal("chemistry", "xHe", 0.1);//He aboundance per H
	//metal abundance at Z=1
	xC_std_ = pin->GetOrAddReal("chemistry", "xC", 1.6e-4); 
	xO_std_ = pin->GetOrAddReal("chemistry", "xO", 3.2e-4);
	xSi_std_ = pin->GetOrAddReal("chemistry", "xSi", 1.7e-6);
	//cosmic ray ionization rate per H
	cr_rate0_ = pin->GetOrAddReal("chemistry", "CR", 2e-16);
  //units of density and radiation
	unit_density_in_nH_ = pin->GetReal("chemistry", "unit_density_in_nH");
	unit_length_in_cm_ = pin->GetReal("chemistry", "unit_length_in_cm");
	unit_vel_in_cms_ = pin->GetReal("chemistry", "unit_vel_in_cms");
  unit_time_in_s_ = unit_length_in_cm_/unit_vel_in_cms_;
  unit_E_in_cgs_ = 1.67e-24 * (
                            1 + xHe_*4 + xC_std_*12 + xO_std_*16 + xSi_std_*28
                    ) * unit_density_in_nH_ * unit_vel_in_cms_ * unit_vel_in_cms_;
	unit_radiation_in_draine1987_ = pin->GetReal(
                                "chemistry", "unit_radiation_in_draine1987");
  //check whether number of frequencies equal to the input file specification
  const int nfreq = pin->GetOrAddInteger("radiation", "n_frequency", 1);
  std::stringstream msg; //error message
  if (nfreq != n_freq_) {
    msg << "### FATAL ERROR in ChemNetwork constructor" << std::endl
      << "number of frequencies in radiation: " << nfreq 
      << " not equal to that in chemistry: " << n_freq_  << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  if (NON_BAROTROPIC_EOS) {
    temperature_ = 0.;
  } else {
    //isothermal
    temperature_ = pin->GetReal("chemistry", "temperature");
  }
  //CR shielding
  is_cr_shielding_ = pin->GetOrAddInteger("chemistry", "is_cr_shielding", 0);
  //H2 rovibrational cooling
  is_H2_rovib_cooling_ = pin->GetOrAddInteger("chemistry", "isH2RVcooling", 1);
	//temperature above or below which heating and cooling is turned off
	Real inf = std::numeric_limits<Real>::infinity();
	temp_max_heat_ = pin->GetOrAddReal("chemistry", "temp_max_heat", inf);
	temp_min_cool_ = pin->GetOrAddReal("chemistry", "temp_min_cool", 1.);
	//minimum temperature for reaction rates, also applied to energy equation
	temp_min_rates_ = pin->GetOrAddReal("chemistry", "temp_min_rates", 1.);
  //cap temperature when calculating rates 
  //do not apply for collisional dissociation reactions and energy equation
	temp_max_rates_ = pin->GetOrAddReal("chemistry", "temp_max_rates", inf);
	//CO cooling parameters
	//Maximum CO cooling length. default 100pc.
	Leff_CO_max_ = pin->GetOrAddReal("chemistry", "Leff_CO_max", 3.0e20);
	gradv_ = 0.;
	
  //atomic abundance
  xC_ = zdg_ * xC_std_;
  xO_ = zdg_ * xO_std_;
  xSi_ = zdg_ * xSi_std_;

  //initialize rates to zero
  for (int i=0; i<n_cr_; i++) {
    kcr_[i] = 0;
  }
  for (int i=0; i<n_2body_; i++) {
    k2body_[i] = 0;
  }
  for (int i=0; i<n_ph_; i++) {
    kph_[i] = 0;
  }
  for (int i=0; i<n_gr_; i++) {
    kgr_[i] = 0;
  }
  //copy species to a full list of species names
  for (int i=0; i<NSCALARS; i++) {
    species_names_all_[i] = species_names[i];
  }
  for (int i=NSCALARS; i<NSCALARS+ngs_; i++) {
    species_names_all_[i] = ghost_species_names_[i-NSCALARS];
  }

}

ChemNetwork::~ChemNetwork() {}

void ChemNetwork::RHS(const Real t, const Real y[NSCALARS], const Real ED, 
                      Real ydot[NSCALARS]) {
	Real rate;
	//store previous y includeing ghost species
	Real yprev[NSCALARS+ngs_];
	Real yprev0[NSCALARS+ngs_];//correct negative abundance
	Real ydotg[NSCALARS+ngs_];
  Real E_ergs = ED * unit_E_in_cgs_ / nH_; //ernergy per hydrogen atom

	for(int i=0; i<NSCALARS+ngs_; i++) {
		ydotg[i] = 0.0;
  }

  // copy y to yprev and set ghost species
  GetGhostSpecies(y, yprev);
  //correct negative abundance to zero, used in rate update
  for (int i=0; i<NSCALARS+ngs_; i++) {
    if (yprev[i] < 0) {
      yprev0[i] = 0;
    } else {
      yprev0[i] = yprev[i];
    }
    //throw error if nan, or inf, or large negative value occurs
    if ( isnan(yprev[i]) || isinf(yprev[i]) ) {
      printf("RHS: ");
      for (int j=0; j<NSCALARS+ngs_; j++) {
        printf("%s: %.2e  ", species_names_all_[j].c_str(), yprev[j]);
      }
      printf("\n");
      OutputRates(stdout);
      printf("rad_ = ");
      for (int ifreq=0; ifreq < n_freq_; ++ifreq) {
        printf("%.2e  ", rad_[ifreq]);
      }
      printf("\n");
      printf("nH_ = %.2e\n", nH_);
      std::stringstream msg;
      msg << "ChemNetwork (gow17): RHS(yprev): nan or inf" << std::endl;
      ATHENA_ERROR(msg);
    }
  }

  UpdateRates(yprev0, E_ergs);

#ifdef DEBUG
  if (output_flag) {
    FILE *pf = fopen("chem_network.dat", "w");
    OutputRates(pf);
    fclose(pf);
    output_flag = false;
  }
#endif

  //cosmic ray reactions
  for (int i=0; i<n_cr_; i++) {
    rate = kcr_[i] * yprev[incr_[i]];
    ydotg[incr_[i]] -= rate;
    ydotg[outcr_[i]] += rate;
  }

  //2body reactions
  for (int i=0; i<n_2body_; i++) {
    rate =  k2body_[i] * yprev[in2body1_[i]] * yprev[in2body2_[i]];
    if (yprev[in2body1_[i]] < 0 && yprev[in2body2_[i]] < 0) {
      rate *= -1.;
    }
    ydotg[in2body1_[i]] -= rate;
    ydotg[in2body2_[i]] -= rate;
    ydotg[out2body1_[i]] += rate;
    ydotg[out2body2_[i]] += rate;
  }

  //photo reactions
  for (int i=0; i<n_ph_; i++) {
    rate = kph_[i] * yprev[inph_[i]];
    ydotg[inph_[i]] -= rate;
    ydotg[outph1_[i]] += rate;
  }

  //grain assisted reactions
  for (int i=0; i<n_gr_; i++) {
    rate = kgr_[i] * yprev[ingr_[i]];
    ydotg[ingr_[i]] -= rate;
    ydotg[outgr_[i]] += rate;
  }

	//set ydot to return
	for (int i=0; i<NSCALARS; i++) {
    //return in code units
		ydot[i] = ydotg[i] * unit_time_in_s_;
	}

  //throw error if nan, or inf, or large value occurs
  for (int i=0; i<NSCALARS; i++) {
    if ( isnan(ydot[i]) || isinf(ydot[i]) ) {
      printf("ydot: ");
      for (int j=0; j<NSCALARS; j++) {
        printf("%s: %.2e  ", species_names[j].c_str(), ydot[j]);
      }
      printf("abundances: ");
      for (int j=0; j<NSCALARS+ngs_; j++) {
        printf("%s: %.2e  ", species_names_all_[j].c_str(), yprev[j]);
      }
      printf("\n");
      OutputRates(stdout);
      printf("rad_ = ");
      for (int ifreq=0; ifreq < n_freq_; ++ifreq) {
        printf("%.2e  ", rad_[ifreq]);
      }
      printf("\n");
      printf("nH_ = %.2e\n", nH_);
      printf("ED = %.2e\n", ED);
      printf("E_ergs = %.2e\n", E_ergs);
      printf("unit_E_in_cgs_ = %.2e\n", unit_E_in_cgs_);
      printf("T = %.2e\n", E_ergs/Thermo::CvCold(yprev0[iH2_], xHe_, yprev0[ige_]));
      std::stringstream msg;
      msg << "ChemNetwork (gow17): RHS(ydot): nan or inf" << std::endl;
      ATHENA_ERROR(msg);
    }
  }

  return;
}


void ChemNetwork::InitializeNextStep(const int k, const int j, const int i) {
  Real rad_sum, temp, NH, rho, rho_floor;
  int nang = pmy_mb_->prad->nang;
  //density
  rho = pmy_mb_->phydro->w(IDN, k, j, i);
  //apply density floor
  rho_floor = pmy_mb_->peos->GetDensityFloor();
  rho = (rho > rho_floor) ?  rho : rho_floor;
  //hydrogen atom number density
  nH_ =  rho * unit_density_in_nH_;
  //average radiation field of all angles
  //TODO: put floor on radiation variables?
  for (int ifreq=0; ifreq < n_freq_; ++ifreq) {
    rad_sum = 0;
    //radiation
    for (int iang=0; iang < nang; ++iang) {
      rad_sum += pmy_mb_->prad->ir(k, j, i, ifreq * nang + iang);
    }
    if (ifreq == index_cr_) {
      rad_[index_cr_] = rad_sum / float(nang);
    } else {
      rad_[ifreq] = rad_sum * unit_radiation_in_draine1987_ / float(nang) ;
    }
#ifdef DEBUG
    if (isnan(rad_[ifreq])) {
      printf("InitializeNextStep: ");
      printf("ifreq=%d, nang=%d, rad_sum=%.2e\n", ifreq, nang, rad_sum);
      OutputRates(stdout);
    }
#endif
  }
  //CO cooling paramters
  //TODO: for six-ray, this should be in the right units
  SetGrad_v(k, j, i);
	return;
}

void ChemNetwork::GetGhostSpecies(const Real *y, Real yghost[NSCALARS+ngs_]) {
	//copy the aboundances in y to yghost
	for (int i=0; i<NSCALARS; i++) {
		yghost[i] = y[i];
	}
	//set the ghost species
 	yghost[igC_] = xC_ - yghost[iHCOplus_] -  yghost[iCHx_] 
                     - yghost[iCO_] - yghost[iCplus_]; 
	yghost[igO_] = xO_ - yghost[iHCOplus_] -  yghost[iOHx_] 
                     - yghost[iCO_] - yghost[iOplus_]; 
	yghost[igHe_] = xHe_ - yghost[iHeplus_]; 
	yghost[igSi_] = xSi_ - yghost[iSiplus_]; 
  yghost[ige_] = yghost[iHeplus_] + yghost[iCplus_] + yghost[iHCOplus_]
                     + yghost[iH3plus_] + yghost[iH2plus_] + yghost[iHplus_]
                     + yghost[iOplus_] + yghost[iSiplus_]; 
	yghost[igH_] = 1.0 - (yghost[iOHx_] + yghost[iCHx_] + yghost[iHCOplus_]
                     + 3.0*yghost[iH3plus_] + 2.0*yghost[iH2plus_] + yghost[iHplus_]
										 + 2.0*yghost[iH2_]);
	return;
}

Real ChemNetwork::CII_rec_rate_(const Real temp) {
  Real A, B, T0, T1, C, T2, BN, term1, term2, alpharr, alphadr;
  A = 2.995e-9;
  B = 0.7849;
  T0 =  6.670e-3;
  T1 = 1.943e6;
  C = 0.1597;
  T2 = 4.955e4;
  BN = B + C * exp(-T2/temp);
  term1 = sqrt(temp/T0);
  term2 = sqrt(temp/T1);
  alpharr = A / ( term1*pow(1.0+term1, 1.0-BN) * pow(1.0+term2, 1.0+BN) );
  alphadr = pow( temp, -3.0/2.0 ) * ( 6.346e-9 * exp(-1.217e1/temp) +
        9.793e-09 * exp(-7.38e1/temp) + 1.634e-06 * exp(-1.523e+04/temp) );
  return (alpharr+alphadr);
}

void ChemNetwork::UpdateRates(const Real y[NSCALARS+ngs_], const Real E) {
  Real T, Tcoll;
  //constant or evolve temperature
  if (NON_BAROTROPIC_EOS) {
    T = E / Thermo::CvCold(y[iH2_], xHe_, y[ige_]);
  } else {
    //isohermal EOS
    T = temperature_;
  }
  Tcoll = T;
	//cap T above some minimum temperature
	if (T < temp_min_rates_) {
		T = temp_min_rates_;
    Tcoll = T;
	} else if (T > temp_max_rates_) {
    //do not put upper limit on the temperature for collisional ionization and
    //dissociation rates Tcoll
    T = temp_max_rates_;
  }
  const Real logT = log10(T);
	const Real logT4coll = log10(Tcoll/1.0e4);
	const Real lnTecoll = log(Tcoll * 8.6173e-5);
  Real ncr, n2ncr;
	Real psi; //H+ grain recombination parameter
  Real kcr_H_fac;//ratio of total rate to primary rate
  Real psi_gr_fac_;
  const Real kida_fac = ( 0.62 + 45.41 / sqrt(T) ) * nH_;
  Real t1_CHx, t2_CHx;
	//cosmic ray reactions
	for (int i=0; i<n_cr_; i++) {
		kcr_[i] = kcr_base_[i] * rad_[index_cr_];
	}
  //cosmic ray induced photo-reactions, proportional to x(H2)
  //(0) cr + H2 -> H2+ + *e
  //(1) cr + *He -> He+ + *e 
  //(2) cr + *H  -> H+ + *e 
  //(3) cr + *C -> C+ + *e     --including direct and cr induce photo reactions 
  //(4) crphoto + CO -> *O + *C 
  //(5) cr + CO -> HCO+ + e  --schematic for cr + CO -> CO+ + e
  //(6) cr + Si -> Si+ + e, UMIST12 
  kcr_H_fac = 1.15 * 2*y[iH2_] + 1.5 * y[igH_];
  kcr_[0] *= kcr_H_fac;
  kcr_[2] *= kcr_H_fac;
  kcr_[3] *= (2*y[iH2_] + 3.85/kcr_base_[3]);
  kcr_[4] *= 2*y[iH2_];
  kcr_[6] *= 2*y[iH2_];
	//2 body reactions
	for (int i=0; i<n_2body_; i++){
		k2body_[i] = k2body_base_[i] * pow(T, k2Texp_[i]) * nH_;
	}
	//Special treatment of rates for some equations
  /*(0) H3+ + *C -> CH + H2         --Vissapragada2016 new rates*/
  t1_CHx = A_kCHx_ * pow( 300./T, n_kCHx_);
  t2_CHx = c_kCHx_[0] * exp(-Ti_kCHx_[0]/T) + c_kCHx_[1] * exp(-Ti_kCHx_[1]/T)
           + c_kCHx_[2]*exp(-Ti_kCHx_[2]/T) + c_kCHx_[3] *exp(-Ti_kCHx_[3]/T);
  k2body_[0] *= t1_CHx + pow(T, -1.5) * t2_CHx;
	/*(3) He+ + H2 -> H+ + *He + *H   --fit to Schauer1989 */
	k2body_[3] *= exp(-22.5/T);
  //(5) C+ + H2 -> CH + *H         -- schematic reaction for C+ + H2 -> CH2+
  k2body_[5] *= exp(-23./T);
  // ---branching of C+ + H2 ------
  //(22) C+ + H2 + *e -> *C + *H + *H
  k2body_[22] *= exp(-23./T);
  // (6) C+ + OH -> HCO+             -- Schematic equation for C+ + OH -> CO+ + H.
  // Use rates in KIDA website.
  k2body_[6] = 9.15e-10 * kida_fac;
  //(8) OH + *C -> CO + *H          --exp(0.108/T)
  k2body_[8] *= exp(0.108/T);
	//(9) He+ + *e -> *He             --(17) Case B 
	k2body_[9] *= 11.19 + (-1.676 + (-0.2852 + 0.04433*logT) * logT )* logT;
  // (11) C+ + *e -> *C              -- Include RR and DR, Badnell2003, 2006. 
  k2body_[11] = CII_rec_rate_(T) * nH_;
  // (13) H2+ + H2 -> H3+ + *H       --(54) exp(-T/46600) 
  k2body_[13] *= exp(- T/46600.);
	// (14) H+ + *e -> *H              --(12) Case B 
	k2body_[14] *= pow( 315614.0 / T, 1.5) 
									 * pow(  1.0 + pow( 115188.0 / T, 0.407) , -2.242 );
  //--- H2O+ + e branching--
  //(1) H3+ + *O -> OH + H2        
  //(24) H3+ + *O + *e -> H2 + *O + *H     
	Real h2oplus_ratio, fac_H2Oplus_H2, fac_H2Oplus_e;
	if (y[ige_] < small_) {
		h2oplus_ratio = 1.0e10;
	} else {
		h2oplus_ratio = 6e-10 * y[iH2_] / ( 5.3e-6 / sqrt(T) * y[ige_] );
	}
  fac_H2Oplus_H2 = h2oplus_ratio / (h2oplus_ratio + 1.);
  fac_H2Oplus_e = 1. / (h2oplus_ratio + 1.);
  k2body_[1] *= fac_H2Oplus_H2;
  k2body_[24] *= fac_H2Oplus_e;
  // (25) He+ + OH -> *H + *He + *O(O+)
  k2body_[25] = 1.35e-9 * kida_fac;
  //  --- O+ reactions ---
  //  (27) H+ + *O -> O+ + *H -- exp(-227/T)
  //  (28) O+ + *H -> H+ + *O 
  //  (29) O+ + H2 -> OH + *H     -- branching of H2O+
  //  (30) O+ + H2 -> *O + *H + *H  -- branching of H2O+ */
  k2body_[27] *= ( 1.1e-11 * pow(T, 0.517) + 4.0e-10 * pow(T, 6.69e-3) )*exp(-227./T);
  k2body_[28] *= 4.99e-11* pow(T, 0.405) + 7.5e-10 * pow(T, -0.458);
  k2body_[29] *= fac_H2Oplus_H2;
  k2body_[30] *= fac_H2Oplus_e;


  //Collisional dissociation, k>~1.0e-30 at T>~5e2.
  Real k9l, k9h, k10l, k10h, ncrH, ncrH2, div_ncr;
  if (Tcoll > temp_coll_) {
    //(15) H2 + *H -> 3 *H   
    //(16) H2 + H2 -> H2 + 2 *H
    // --(9) Density dependent. See Glover+MacLow2007
  	k9l = 6.67e-12 * sqrt(Tcoll) * exp(-(1. + 63590./Tcoll)); 
    k9h = 3.52e-9 * exp(-43900.0 / Tcoll);
    k10l = 5.996e-30 * pow(Tcoll, 4.1881) / pow((1.0 + 6.761e-6 * Tcoll), 5.6881)  
            * exp(-54657.4 / Tcoll);
    k10h = 1.3e-9 * exp(-53300.0 / Tcoll); 
    ncrH = pow(10, (3.0 - 0.416 * logT4coll - 0.327 * logT4coll*logT4coll));
    ncrH2 = pow(10, (4.845 - 1.3 * logT4coll + 1.62 * logT4coll*logT4coll));
		div_ncr = y[igH_]/ncrH + y[iH2_]/ncrH2;
		if (div_ncr < small_) {
			ncr = 1./ small_;
		} else {
			ncr = 1. / div_ncr;
		}
    n2ncr = nH_ / ncr;
    k2body_[15] = pow(10, log10(k9h) *  n2ncr/(1. + n2ncr) 
                         + log10(k9l) / (1. + n2ncr)) * nH_;
    k2body_[16] = pow(10, log10(k10h) *  n2ncr/(1. + n2ncr) 
                         + log10(k10l) / (1. + n2ncr)) * nH_;
    // (17) *H + *e -> H+ + 2 *e       --(11) Relates to Te 
    k2body_[17] *= exp( -3.271396786e1 + 
                      (1.35365560e1 + (- 5.73932875 + (1.56315498 
                    + (- 2.877056e-1 + (3.48255977e-2 + (- 2.63197617e-3
                    + (1.11954395e-4 + (-2.03914985e-6)
        *lnTecoll)*lnTecoll)*lnTecoll)*lnTecoll)*lnTecoll)*lnTecoll)*lnTecoll)
                      *lnTecoll);
  } else {
    k2body_[15] = 0.;
    k2body_[16] = 0.;
    k2body_[17] = 0.;
  }
  
	//photo reactions
	for (int i=0; i<n_ph_; i++) {
    kph_[i] = kph_base_[i] * rad_[i];
	}

	// Grain assisted recombination of H and H2
	//	 (0) *H + *H + gr -> H2 + gr , from Draine book chapter 31.2 page 346,
	//	 Jura 1975
	kgr_[0] = 3.0e-17 * nH_ * zdg_;
	//	 (1) H+ + *e + gr -> *H + gr
  //	 (2) C+ + *e + gr -> *C + gr
  //   (3) He+ + *e + gr -> *He + gr
  //   (4) Si+ + *e + gr -> *Si + gr
  //   , rate dependent on e aboundance. 
	if (y[ige_] > small_ ) {
    //set lower limit to radiation field in calculating kgr_ to avoid nan
    //values.
    Real GPE_limit = 1.0e-10;
    Real GPE0 = rad_[index_gpe_];
    if (GPE0 < GPE_limit) {
      GPE0 = GPE_limit;
    }
		psi_gr_fac_ = 1.7 * GPE0 * sqrt(T) / nH_; 
		psi = psi_gr_fac_ / y[ige_];
		kgr_[1] = 1.0e-14 * cHp_[0] / 
								 (
									 1.0 + cHp_[1]*pow(psi, cHp_[2]) * 
										 (1.0 + cHp_[3] * pow(T, cHp_[4])
																	 *pow( psi, -cHp_[5]-cHp_[6]*log(T) ) 
										 ) 
									) * nH_ * zdg_;
		kgr_[2] = 1.0e-14 * cCp_[0] / 
								 (
									 1.0 + cCp_[1]*pow(psi, cCp_[2]) * 
										 (1.0 + cCp_[3] * pow(T, cCp_[4])
																	 *pow( psi, -cCp_[5]-cCp_[6]*log(T) ) 
										 ) 
									) * nH_ * zdg_;
		kgr_[3] = 1.0e-14 * cHep_[0] / 
								 (
									 1.0 + cHep_[1]*pow(psi, cHep_[2]) * 
										 (1.0 + cHep_[3] * pow(T, cHep_[4])
																	 *pow( psi, -cHep_[5]-cHep_[6]*log(T) ) 
										 ) 
									) * nH_ * zdg_;
		kgr_[4] = 1.0e-14 * cSip_[0] / 
								 (
									 1.0 + cSip_[1]*pow(psi, cSip_[2]) * 
										 (1.0 + cSip_[3] * pow(T, cSip_[4])
																	 *pow( psi, -cSip_[5]-cSip_[6]*log(T) ) 
										 ) 
									) * nH_ * zdg_;
	} else {
		for (int i=1; i<5; i++) {
			kgr_[i] = 0.;
		}
	}

  return;
}

Real ChemNetwork::Edot(const Real t, const Real y[NSCALARS], const Real ED) {
  Real E_ergs = ED * unit_E_in_cgs_ / nH_; //ernergy per hydrogen atom
  //isothermal
  if (!NON_BAROTROPIC_EOS) {
    return 0;
  }

  Real T = 0.;
  Real dEdt = 0.;
	Real yprev[NSCALARS+ngs_];
  // copy y to yprev and set ghost species
  GetGhostSpecies(y, yprev);
  //correct negative abundance to zero
  for (int i=0; i<NSCALARS+ngs_; i++) {
    if (yprev[i] < 0) {
      yprev[i] = 0;
    }
  }

  //temperature
  T = E_ergs / Thermo::CvCold(yprev[iH2_], xHe_, yprev[ige_]);
  //apply temperature floor, incase of very small or negative energy
	if (T < temp_min_rates_) {
		T = temp_min_rates_;
  }

  //--------------------------heating-----------------------------
  //cosmic ray ionization of H, He, and H2 
  //NOTE: because these depends on rates, make sure ChemInit is called before.
  //NOTE: the kcr_[i] assume the order of equastions are not changed
	//cut-off heating at high temperature
	Real GCR, GPE, GH2gr, dot_xH2_photo, GH2pump, GH2diss;
	if (T > temp_max_heat_) {
		GCR = 0.;
		GPE = 0;
		GH2gr = 0;
		GH2pump = 0.;
		GH2diss = 0.;
	} else {
		GCR = Thermo::HeatingCr(yprev[ige_],  nH_,
				yprev[igH_],  yprev[igHe_],  yprev[iH2_],
				kcr_[icr_H_],  kcr_[icr_He_],  kcr_[icr_H2_]);
		//photo electric effect on dust
		GPE = Thermo::HeatingPE(rad_[index_gpe_], zdg_, T, nH_*yprev[ige_]);
		//H2 formation on dust grains
		GH2gr = Thermo::HeatingH2gr(yprev[igH_],  yprev[iH2_],  nH_,
				T,  kgr_[igr_H_]);
		//H2 UV pumping
		dot_xH2_photo = kph_[iph_H2_] * yprev[iH2_];
		GH2pump = Thermo::HeatingH2pump(yprev[igH_],  yprev[iH2_],  nH_,
				T,  dot_xH2_photo);
		//H2 photo dissiociation.
		GH2diss = Thermo::HeatingH2diss(dot_xH2_photo);
	}
  //--------------------------cooling-----------------------------
	//cut-off cooling at low temperature
	Real LCII, LCI, LOI, LHotGas, LCOR, LH2, LDust, LRec, LH2diss, LHIion;
	Real vth, nCO, grad_small_;
  Real NCOeff, gradeff;
	if (T < temp_min_cool_) {
		LCII = 0.;
		LCI = 0;
		LOI = 0.;
		LHotGas = 0;
		LCOR = 0;
		LH2 = 0;
		LDust = 0;
		LRec = 0;
		LH2diss = 0;
		LHIion = 0;
	} else {
		// C+ fine structure line 
		LCII = Thermo::CoolingCII(yprev[iCplus_],  nH_*yprev[igH_],  nH_*yprev[iH2_],
				nH_*yprev[ige_],  T);
		// CI fine structure line 
		LCI = Thermo:: CoolingCI(yprev[igC_],  nH_*yprev[igH_],  nH_*yprev[iH2_],
				nH_*yprev[ige_],  T);
		// OI fine structure line 
		LOI = Thermo:: CoolingOI(yprev[igO_],  nH_*yprev[igH_],  nH_*yprev[iH2_],
				nH_*yprev[ige_],  T);
		// cooling of hot gas: radiative cooling, free-free.
		LHotGas = Thermo::CoolingHotGas(nH_,  T, zdg_);
		// CO rotational lines 
		// Calculate effective CO column density
		vth = sqrt(2. * Thermo::kb_ * T / ChemistryUtility::mCO);
		nCO = nH_ * yprev[iCO_];
		grad_small_ = vth/Leff_CO_max_;
    gradeff = std::max(gradv_, grad_small_);
    NCOeff = nCO / gradeff;
		LCOR = Thermo::CoolingCOR(yprev[iCO_], nH_*yprev[igH_],  nH_*yprev[iH2_],
				nH_*yprev[ige_],  T,  NCOeff);
		// H2 vibration and rotation lines 
    if (is_H2_rovib_cooling_ != 0) {
      LH2 = Thermo::CoolingH2(yprev[iH2_], nH_*yprev[igH_],  nH_*yprev[iH2_],
          nH_*yprev[igHe_],  nH_*yprev[iHplus_], nH_*yprev[ige_],
          T);
    } else {
      LH2 = 0.;
    }
		// dust thermo emission 
		LDust = Thermo::CoolingDustTd(zdg_,  nH_, T, 10.);
		// reconbination of e on PAHs 
		LRec = Thermo::CoolingRec(zdg_,  T,  nH_*yprev[ige_], rad_[index_gpe_]);
		// collisional dissociation of H2 
		LH2diss = Thermo::CoolingH2diss(yprev[igH_],  yprev[iH2_], k2body_[i2body_H2_H],
				k2body_[i2body_H2_H2]);
		// collisional ionization of HI 
		LHIion = Thermo::CoolingHIion(yprev[igH_],  yprev[ige_],
				k2body_[i2body_H_e]);
	}
  dEdt = (GCR + GPE + GH2gr + GH2pump + GH2diss)
            - (LCII + LCI + LOI + LHotGas + LCOR 
                + LH2 + LDust + LRec + LH2diss + LHIion);
	if ( isnan(dEdt) || isinf(dEdt) ) {
    if ( isnan(LCOR) || isinf(LCOR) ) {
      printf("NCOeff=%.2e, gradeff=%.2e, gradv_=%.2e, vth=%.2e, nH_=%.2e, nCO=%.2e\n",
          NCOeff, gradeff, gradv_, vth, nH_, nCO);
    }
		printf("GCR=%.2e, GPE=%.2e, GH2gr=%.2e, GH2pump=%.2e GH2diss=%.2e\n",
				GCR , GPE , GH2gr , GH2pump , GH2diss);
		printf("LCII=%.2e, LCI=%.2e, LOI=%.2e, LHotGas=%.2e, LCOR=%.2e\n",
				LCII , LCI , LOI , LHotGas , LCOR);
		printf("LH2=%.2e, LDust=%.2e, LRec=%.2e, LH2diss=%.2e, LHIion=%.2e\n",
				LH2 , LDust , LRec , LH2diss , LHIion);
		printf("T=%.2e, dEdt=%.2e, E=%.2e, Cv=%.2e, nH=%.2e\n", T, dEdt, E_ergs,
				Thermo::CvCold(yprev[iH2_], xHe_, yprev[ige_]), nH_);
		for (int i=0; i<NSCALARS+ngs_; i++) {
			printf("%s: %.2e  ", species_names_all_[i].c_str(), yprev[i]);
		}
		printf("\n");
    std::stringstream msg;
    msg << "ChemNetwork (gow17): dEdt: nan or inf number" << std::endl;
    ATHENA_ERROR(msg);
	}
  //return in code units
  Real dEDdt = dEdt * nH_ / unit_E_in_cgs_ * unit_time_in_s_;
  return dEDdt;
}

void ChemNetwork::OutputRates(FILE *pf) const {
  //output the reactions and base rates
	for (int i=0; i<n_cr_; i++) {
		fprintf(pf, "cr  + %4s -> %4s,     kcr = %.2e\n", 
		 species_names_all_[incr_[i]].c_str(), species_names_all_[outcr_[i]].c_str(),
		 kcr_[i]);
	}
	for (int i=0; i<n_2body_; i++) {
		fprintf(pf, "%4s  + %4s -> %4s  + %4s,     k2body = %.2e\n", 
		 species_names_all_[in2body1_[i]].c_str(),
		 species_names_all_[in2body2_[i]].c_str(),
		 species_names_all_[out2body1_[i]].c_str(),
		 species_names_all_[out2body2_[i]].c_str(),
		 k2body_[i]);
	}
	for (int i=0; i<n_ph_; i++) {
		fprintf(pf, "h nu  + %4s -> %4s,     kph = %.2e\n", 
		 species_names_all_[inph_[i]].c_str(), species_names_all_[outph1_[i]].c_str(),
		 kph_[i]);
	}
	for (int i=0; i<n_gr_; i++) {
		fprintf(pf, "gr  + %4s -> %4s,       kgr = %.2e\n", 
		 species_names_all_[ingr_[i]].c_str(), species_names_all_[outgr_[i]].c_str(),
		 kgr_[i]);
	}
  return;
}

Real ChemNetwork::GetStddev(Real arr[], const int len) {
  Real sum=0, avg=0, sum_sq=0, avg_sq=0;
  for (int i=0; i<len; i++) {
    sum += arr[i];
    sum_sq += arr[i] * arr[i];
  }
  avg = sum/Real(len);
  avg_sq = sum_sq/Real(len);
  return sqrt(avg_sq - avg*avg);
}

void ChemNetwork::SetGrad_v(const int k, const int j, const int i) {
  AthenaArray<Real> &w = pmy_mb_->phydro->w;
  Real dvdx, dvdy, dvdz, dvdr_avg, di1, di2;
  Real dx1, dx2, dy1, dy2, dz1, dz2;
  Real dndx, dndy, dndz, gradn;
  //velocity gradient, same as LVG approximation in RADMC-3D when calculating
  //CO line emission.
  //vx
  di1 = w(IVX, k, j, i+1) - w(IVX, k, j, i);
  dx1 = ( pmy_mb_->pcoord->dx1f(i+1)+pmy_mb_->pcoord->dx1f(i) )/2.;
  di2 = w(IVX, k, j, i) - w(IVX, k, j, i-1);
  dx2 = ( pmy_mb_->pcoord->dx1f(i)+pmy_mb_->pcoord->dx1f(i-1) )/2.;
  dvdx = (di1/dx1 + di2/dx2)/2.;
  //vy
  di1 = w(IVY, k, j+1, i) - w(IVY, k, j, i);
  dy1 = ( pmy_mb_->pcoord->dx2f(j+1)+pmy_mb_->pcoord->dx2f(j) )/2.;
  di2 = w(IVY, k, j, i) - w(IVY, k, j-1, i);
  dy2 = ( pmy_mb_->pcoord->dx2f(j)+pmy_mb_->pcoord->dx2f(j-1) )/2.;
  dvdy = (di1/dy1 + di2/dy2)/2.;
  //vz
  di1 = w(IVZ, k+1, j, i) - w(IVZ, k, j, i);
  dz1 = ( pmy_mb_->pcoord->dx3f(k+1)+pmy_mb_->pcoord->dx3f(k) )/2.;
  di2 = w(IVZ, k, j, i) - w(IVZ, k-1, j, i);
  dz2 = ( pmy_mb_->pcoord->dx3f(k)+pmy_mb_->pcoord->dx3f(k-1) )/2.;
  dvdz = (di1/dz1 + di2/dz2)/2.;
  dvdr_avg = ( fabs(dvdx) + fabs(dvdy) + fabs(dvdz) ) / 3.;
  //asign gradv_, in cgs.
  gradv_ = dvdr_avg * unit_vel_in_cms_ / unit_length_in_cm_;
  return;
}
