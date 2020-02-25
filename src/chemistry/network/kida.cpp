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
//! \file kida.cpp
//  \brief implementation of functions in class ChemNetwork, using the simple
//  network for kida style network files.
//======================================================================================

// this class header
#include "kida.hpp"

//athena++ header
#include "network.hpp"
#include "../../scalars/scalars.hpp"
#include "../../parameter_input.hpp"       //ParameterInput
#include "../../mesh/mesh.hpp"
#include "../../hydro/hydro.hpp"
#include "../utils/chemistry_utils.hpp"
#include "../../utils/string_utils.hpp"
#include "../../defs.hpp"
#include "../../eos/eos.hpp"

//c++ header
#include <sstream>    // stringstream
#include <iostream>   // endl
#include <limits>    //inf
#include <fstream>   //file()
#include <stdio.h>    // c style file
#include <algorithm>    // std::find()

ChemNetwork::ChemNetwork(MeshBlock *pmb, ParameterInput *pin) {
	//number of species and a list of name of species
  pmy_spec_ = pmb->pscalars;
	pmy_mb_ = pmb;

	//set the parameters from input file
	xi_cr_ = pin->GetOrAddReal("chemistry", "xi_cr", 2e-16);
	zdg_ = pin->GetOrAddReal("chemistry", "Zdg", 1.);//dust and gas metallicity
  //units
	unit_density_in_nH_ = pin->GetReal("chemistry", "unit_density_in_nH");
	unit_length_in_cm_ = pin->GetReal("chemistry", "unit_length_in_cm");
	unit_vel_in_cms_ = pin->GetReal("chemistry", "unit_vel_in_cms");
  unit_time_in_s_ = unit_length_in_cm_/unit_vel_in_cms_;
  unit_E_in_cgs_ = 1.67e-24 * 1.4 * unit_density_in_nH_
                           * unit_vel_in_cms_ * unit_vel_in_cms_;
  //folder of the network
  network_dir_ = pin->GetString("chemistry", "network_dir");

  //read in the species
  std::string species_file_name = network_dir_ + "/species.dat";
  std::ifstream species_file(species_file_name);
  if (!species_file) {
    std::stringstream msg; //error message
    msg << "### FATAL ERROR in ChemNetwork constructor [ChemNetwork]" << std::endl
      << "Cannot open file" << species_file_name << std::endl;
    ATHENA_ERROR(msg);
  }
  std::string line;
  int nline = 0;
  while (getline(species_file, line)) {
    //trim white spaces
    StringUtils::trim(line);
    //skip blank lines and comments
    if(line.empty() || (line.find("!") == 0)) {
        continue;
    }
    KidaSpecies si(line, nline);
    species_.push_back(si);
    species_names[nline] = si.name;
    ispec_map_[si.name] = nline;
    nline++;
  }
  if (nline != NSCALARS) {
    std::stringstream msg; //error message
    msg << "### FATAL ERROR in ChemNetwork constructor [ChemNetwork]" << std::endl
      << "number of species in species.dat does not match the number of scalars (" 
      << NSCALARS << ")" << std::endl;
    ATHENA_ERROR(msg);
  }

  //read in the reactions
  std::string reactions_file_name = network_dir_ + "/reactions.dat";
  std::ifstream reactions_file(reactions_file_name);
  std::vector<int> rids; //array of id for reactions
  if (!reactions_file) {
    std::stringstream msg; //error message
    msg << "### FATAL ERROR in ChemNetwork constructor [ChemNetwork]" << std::endl
      << "Cannot open file" << reactions_file_name << std::endl;
    ATHENA_ERROR(msg);
  }
  nline = 0;
  while (getline(reactions_file, line)) {
    //trim white spaces
    StringUtils::trim(line);
    //skip blank lines and comments
    if(line.empty() || (line.find("!") == 0)) {
        continue;
    }
    KidaReactions ri(line);
    if (std::find(rids.begin(), rids.end(), ri.id_) == rids.end()) {
      rids.push_back(ri.id_);
      reactions_.push_back(ri);
    } else {
      std::stringstream msg; //error message
      msg << "### FATAL ERROR in ChemNetwork constructor [ChemNetwork]" << std::endl
        << "reaction ID ( " << ri.id_ << ") not unique" << std::endl;
      ATHENA_ERROR(msg);
    }
    nline++;
  }
  nr_ = nline; //number of reactions

  //initialize coefficients of reactions
  InitializeReactions();

#ifdef DEBUG
  PrintProperties();
#endif
}

ChemNetwork::~ChemNetwork() {}

void ChemNetwork::InitializeReactions() {
  KidaReactions *pr = NULL;
  //error message
  bool error=false;
  for (int ir=0; ir<nr_; ir++) {
    pr = &reactions_[ir];
    //---------------- 1 - direct cosmic-ray ionization --------------
    if (pr->itype_ == 1) {
      n_cr_ = 0;
      //check format of reaction 
      std::string in_spec; //input species
      if (pr->reactants_.size() == 2 && pr->products_.size() == 2
          && (pr->reactants_[0] == "CR" || pr->reactants_[1] == "CR") ) {
        if (pr->reactants_[0] == "CR") {
          in_spec = pr->reactants_[1];
        } else {
          in_spec = pr->reactants_[0];
        }
      } else {
        std::stringstream msg; 
        msg << "### FATAL ERROR in ChemNetwork InitializeReactions() [ChemNetwork]"
            << std::endl << "Wrong format in CR reaction ID=" << pr->id_ << std::endl;
      }
      //set indexing arrays
      if (pr->formula_ == 1) {
        incr_.push_back(ispec_map_[in_spec]);
        outcr1_.push_back(ispec_map_[ pr->products_[0]]);
        outcr2_.push_back(ispec_map_[ pr->products_[1]]);
        kcr_base_.push_back(pr->alpha_);
        kcr_.push_back(0.);
        n_cr_++;
      } else {
        error = true;
      }

    //-------------------- 9 - grain assisted reaction ----------------
    } else if (pr->itype_ == 9) {
      n_gr_ = 0;
      //check format
      if (pr->reactants_.size() != 2 || pr->products_.size() != 1) {
        std::stringstream msg; 
        msg << "### FATAL ERROR in ChemNetwork InitializeReactions() [ChemNetwork]"
            << std::endl << "Wrong format in gr reaction ID=" << pr->id_ << std::endl;
        ATHENA_ERROR(msg);
      }
      if (pr->reactants_[1] != "e" && pr->reactants_[1] != "H") {
        std::stringstream msg; 
        msg << "### FATAL ERROR in ChemNetwork InitializeReactions() [ChemNetwork]"
            << std::endl << "second reactant must be H or e." << std::endl;
        ATHENA_ERROR(msg);
      }
      //set indexing arrays
      if (pr->formula_ == 7) {
        ingr1_.push_back(ispec_map_[pr->reactants_[0]]);
        ingr2_.push_back(ispec_map_[pr->reactants_[1]]);
        outgr_.push_back(ispec_map_[pr->products_[0]]);
        idmap_gr_[pr->id_] = n_gr_;
        kgr_.push_back(0.);
        n_gr_++;
      } else{
        error = true;
      }
    } else{
      error = true;
    }

    if (error) {
      std::stringstream msg; 
      msg << "### FATAL ERROR in ChemNetwork InitializeReactions() [ChemNetwork]"
          << std::endl
          << "reaction ID=" << pr->id_ << ", itype=" << pr->itype_ 
          << " and forumla=" << pr->formula_ << " undefined." << std::endl;
      ATHENA_ERROR(msg);
    }
  }
  return;
}

void ChemNetwork::UpdateRates(const Real y[NSCALARS], const Real E) {
	//cosmic ray reactions
	for (int i=0; i<n_cr_; i++) {
		kcr_[i] = kcr_base_[i] * xi_cr_;
	}

	// Grain assisted reactions
	//(2) H + H + gr -> H2 + gr , from Draine book chapter 31.2 page 346, Jura 1975
	kgr_[idmap_gr_[2]] = 3.0e-17 * nH_ * zdg_;
  return;
}


void ChemNetwork::PrintProperties() const {
  //print each species.
  for (int i=0; i<NSCALARS; i++) {
    std::cout << "species: " << i << std::endl;
    std::cout << "name=" << species_[i].name << ", index=" << species_[i].index 
      << ", charge=" << species_[i].charge_ << std::endl;
    std::cout << "atom_count_ = ";
    for (int j=0; j < species_[i].natom_; j++) {
      std::cout << species_[i].atom_count_[j] << " ";
    }
    std::cout << std::endl;
  }

  //print each reactions.
  std::cout << "number of reactions: " << nr_ << std::endl;
  for (int i=0; i<reactions_.size(); i++) {
    std::cout << "reaction ID=" << reactions_[i].id_ << ": ";
    for (int j=0; j<reactions_[i].reactants_.size()-1; j++) {
      std::cout << reactions_[i].reactants_[j] << "+";
    }
    std::cout << reactions_[i].reactants_[reactions_[i].reactants_.size()-1]
              << " -> ";

    for (int j=0; j<reactions_[i].products_.size()-1; j++) {
      std::cout << reactions_[i].products_[j] << "+";
    }
    std::cout << reactions_[i].products_[reactions_[i].products_.size()-1]
              << std::endl;

    std::cout << "alpha=" << reactions_[i].alpha_ << "," 
              << "beta=" << reactions_[i]. beta_ << "," 
              << "gamma=" << reactions_[i].gamma_ << "," 
              << "itype=" << reactions_[i].itype_ << ","
              << "forumla=" << reactions_[i].formula_ << std::endl;
  }

  //print reaction coefficients
  std::cout << "CR reations:" << std::endl;
  for (int i=0; i<n_cr_; i++) {
    std::cout<< species_names[incr_[i]] << " + CR -> "
      << species_names[outcr1_[i]] << " + " << species_names[outcr2_[i]] << ", "
      << "kcr_base_=" << kcr_base_[i] << std::endl;
  }
  std::cout << "gr reations:" << std::endl;
  for (int i=0; i<n_gr_; i++) {
    std::cout<< species_names[ingr1_[i]] << " + " << species_names[ingr2_[i]] 
      <<" -> " << species_names[outgr_[i]] << std::endl;
  }
  std::cout << "idmap_gr_: " << std::endl;
  for (std::map<int,int>::const_iterator it=idmap_gr_.begin();
       it!=idmap_gr_.end(); it++) {
    std::cout << it->first << " => " << it->second << std::endl;
  }

  return;
}

void ChemNetwork::OutputRates(FILE *pf) const {
  //output the reactions and base rates
	for (int i=0; i<n_cr_; i++) {
		fprintf(pf, "%4s + CR -> %4s + %4s,     kcr = %.2e\n", 
		 species_names[incr_[i]].c_str(), species_names[outcr1_[i]].c_str(),
     species_names[outcr2_[i]].c_str(), kcr_[i]);
	}
	//for (int i=0; i<n_2body_; i++) {
	//	fprintf(pf, "%4s  + %4s -> %4s  + %4s,     k2body = %.2e\n", 
	//	 species_names[in2body1_[i]].c_str(),
	//	 species_names[in2body2_[i]].c_str(),
	//	 species_names[out2body1_[i]].c_str(),
	//	 species_names[out2body2_[i]].c_str(),
	//	 k2body_[i]);
	//}
	//for (int i=0; i<n_ph_; i++) {
	//	fprintf(pf, "%4s + h nu -> %4s,     kph = %.2e\n", 
	//	 species_names[inph_[i]].c_str(), species_names[outph1_[i]].c_str(),
	//	 kph_[i]);
	//}
	for (int i=0; i<n_gr_; i++) {
		fprintf(pf, "%4s + %4s (+ gr) -> %4s (+ gr),       kgr = %.2e\n", 
		 species_names[ingr1_[i]].c_str(), species_names[ingr2_[i]].c_str(),
     species_names[outgr_[i]].c_str(), kgr_[i]);
	}
  return;
}

void ChemNetwork::InitializeNextStep(const int k, const int j, const int i) {
  Real rho, rho_floor;
  //density
  rho = pmy_mb_->phydro->w(IDN, k, j, i);
  //apply density floor
  rho_floor = pmy_mb_->peos->GetDensityFloor();
  rho = (rho > rho_floor) ?  rho : rho_floor;
  //hydrogen atom number density
  nH_ =  rho * unit_density_in_nH_;
  return;
}

void ChemNetwork::RHS(const Real t, const Real y[NSCALARS], const Real ED,
                      Real ydot[NSCALARS]){
	Real rate = 0;
  Real E_ergs = ED * unit_E_in_cgs_ / nH_; //ernergy per hydrogen atom
	//store previous y includeing negative abundance correction
	Real yprev0[NSCALARS];//correct negative abundance, only for UpdateRates()
	Real ydotg[NSCALARS];

	for(int i=0; i<NSCALARS; i++) {
		ydotg[i] = 0.0;
  }

  //correct negative abundance to zero, used in rate update
  for (int i=0; i<NSCALARS; i++) {
    if (y[i] < 0) {
      yprev0[i] = 0;
    } else {
      yprev0[i] = y[i];
    }
    //throw error if nan, or inf, or large negative value occurs
    if ( isnan(y[i]) || isinf(y[i]) ) {
      printf("RHS: ");
      for (int j=0; j<NSCALARS; j++) {
        printf("%s: %.2e  ", species_names[j].c_str(), y[j]);
      }
      printf("\n");
      OutputRates(stdout);
      //printf("rad_ = ");
      //for (int ifreq=0; ifreq < n_freq_; ++ifreq) {
      //  printf("%.2e  ", rad_[ifreq]);
      //}
      printf("\n");
      printf("nH_ = %.2e\n", nH_);
      std::stringstream msg;
      msg << "ChemNetwork (kida): RHS(y): nan or inf" << std::endl;
      ATHENA_ERROR(msg);
    }
  }

  UpdateRates(yprev0, E_ergs);

  //cosmic ray reactions
  for (int i=0; i<n_cr_; i++) {
    rate = kcr_[i] * y[incr_[i]];
    ydotg[incr_[i]] -= rate;
    ydotg[outcr1_[i]] += rate;
    ydotg[outcr2_[i]] += rate;
  }

  //grain assisted reactions
  for (int i=0; i<n_gr_; i++) {
    rate = kgr_[i] * y[ingr1_[i]];
    ydotg[ingr1_[i]] -= rate;
    ydotg[ingr2_[i]] -= rate;
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
      for (int j=0; j<NSCALARS; j++) {
        printf("%s: %.2e  ", species_names[j].c_str(), y[j]);
      }
      printf("\n");
      OutputRates(stdout);
      //printf("rad_ = ");
      //for (int ifreq=0; ifreq < n_freq_; ++ifreq) {
      //  printf("%.2e  ", rad_[ifreq]);
      //}
      printf("\n");
      printf("nH_ = %.2e\n", nH_);
      printf("ED = %.2e\n", ED);
      printf("E_ergs = %.2e\n", E_ergs);
      printf("unit_E_in_cgs_ = %.2e\n", unit_E_in_cgs_);
      //printf("T = %.2e\n", E_ergs/Thermo::CvCold(yprev0[iH2_], xHe_, yprev0[ige_]));
      std::stringstream msg;
      msg << "ChemNetwork (kida): RHS(ydot): nan or inf" << std::endl;
      ATHENA_ERROR(msg);
    }
  }

  return;
}

Real ChemNetwork::Edot(const Real t, const Real y[NSCALARS], const Real ED){
  const Real dEDdt = 0;
  return dEDdt;
}
