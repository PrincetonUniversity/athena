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
#include "../../radiation/radiation.hpp"
#include "../utils/chemistry_utils.hpp"
#include "../../utils/string_utils.hpp"
#include "../utils/thermo.hpp"
#include "../../defs.hpp"
#include "../../eos/eos.hpp"

//c++ header
#include <sstream>    // stringstream
#include <iostream>   // endl
#include <limits>    //inf
#include <fstream>   //file()
#include <stdio.h>    // c style file
#include <algorithm>    // std::find()

#ifdef DEBUG
static bool output_flag = true;
#endif

ChemNetwork::ChemNetwork(MeshBlock *pmb, ParameterInput *pin) :  
  pmy_spec_(pmb->pscalars), pmy_mb_(pmb), n_cr_(0), n_crp_(0), n_ph_(0), n_2body_(0),
  n_gr_(0), n_sr_(0), n_freq_(0), index_gpe_(0), index_cr_(0) {

	//set the parameters from input file
	zdg_ = pin->GetOrAddReal("chemistry", "Zdg", 1.);//dust and gas metallicity
  //units
	unit_density_in_nH_ = pin->GetReal("chemistry", "unit_density_in_nH");
	unit_length_in_cm_ = pin->GetReal("chemistry", "unit_length_in_cm");
	unit_vel_in_cms_ = pin->GetReal("chemistry", "unit_vel_in_cms");
	unit_radiation_in_draine1987_ = pin->GetReal(
                                "chemistry", "unit_radiation_in_draine1987");
  unit_time_in_s_ = unit_length_in_cm_/unit_vel_in_cms_;
  unit_E_in_cgs_ = 1.67e-24 * 1.4 * unit_density_in_nH_
                           * unit_vel_in_cms_ * unit_vel_in_cms_;
  //temperature
  if (NON_BAROTROPIC_EOS) {
    temperature_ = 0.;
  } else {
    //isothermal
    temperature_ = pin->GetReal("chemistry", "temperature");
  }
	//minimum temperature for reaction rates, also applied to energy equation
	temp_min_rates_ = pin->GetOrAddReal("chemistry", "temp_min_rates", 1.);
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
    KidaReaction ri(line);
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

  //radiation related variables
  const int nfreq = pin->GetOrAddInteger("radiation", "n_frequency", 1);
  n_freq_ = n_ph_ + 2;
  std::stringstream msg;
  //check whether number of frequencies equal to the input file specification
  if (nfreq != n_freq_) {
    msg << "### FATAL ERROR in ChemNetwork constructor" << std::endl
      << "number of frequencies in radiation: " << nfreq 
      << " not equal to that in chemistry: " << n_freq_  << std::endl;
    ATHENA_ERROR(msg);
  }
  index_gpe_ = n_ph_;
  index_cr_ = n_ph_ + 1;
  rad_.NewAthenaArray(n_freq_);

#ifdef DEBUG
  PrintProperties();
#endif
}

ChemNetwork::~ChemNetwork() {}

void ChemNetwork::InitializeReactions() {
  KidaReaction *pr = NULL;
  //error message
  bool error=false;
  for (int ir=0; ir<nr_; ir++) {
    CheckReaction(reactions_[ir]);
    pr = &reactions_[ir];
    //---------------- 1 - direct cosmic-ray ionization --------------
    if (pr->itype_ == 1) {
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
      } else if (pr->formula_ == 7) {
        incr_.push_back(ispec_map_[in_spec]);
        outcr1_.push_back(ispec_map_[ pr->products_[0]]);
        outcr2_.push_back(ispec_map_[ pr->products_[1]]);
        idmap_cr_[pr->id_] = n_cr_;
        kcr_base_.push_back(0.);
        kcr_.push_back(0.);
        n_cr_++;
      } else {
        error = true;
      }

    //---------------- 2 - cosmic-ray induced photo ionization --------
    } else if (pr->itype_ == 2) {
      //check format of reaction 
      std::string in_spec; //input species
      if (pr->reactants_.size() == 2 && pr->products_.size() == 2
          && (pr->reactants_[0] == "CRP" || pr->reactants_[1] == "CRP") ) {
        if (pr->reactants_[0] == "CRP") {
          in_spec = pr->reactants_[1];
        } else {
          in_spec = pr->reactants_[0];
        }
      } else {
        std::stringstream msg; 
        msg << "### FATAL ERROR in ChemNetwork InitializeReactions() [ChemNetwork]"
           << std::endl << "Wrong format in CRP reaction ID=" << pr->id_ << std::endl;
      }
      //set indexing arrays
      if (pr->formula_ == 1) {
        incrp_.push_back(ispec_map_[in_spec]);
        outcrp1_.push_back(ispec_map_[ pr->products_[0]]);
        outcrp2_.push_back(ispec_map_[ pr->products_[1]]);
        kcrp_base_.push_back(pr->alpha_);
        kcrp_.push_back(0.);
        n_crp_++;
      } else if (pr->formula_ == 7) {
        incrp_.push_back(ispec_map_[in_spec]);
        outcrp1_.push_back(ispec_map_[ pr->products_[0]]);
        outcrp2_.push_back(ispec_map_[ pr->products_[1]]);
        idmap_crp_[pr->id_] = n_crp_;
        kcrp_base_.push_back(0.);
        kcrp_.push_back(0.);
        n_crp_++;
      } else {
        error = true;
      }

    //---------------- 3 - FUV ionization/dissociation ----------------
    } else if (pr->itype_ == 3) {
      //check format of reaction 
      std::string in_spec; //input species
      if (pr->reactants_.size() == 2 && pr->products_.size() == 2
          && (pr->reactants_[0] == "Photon" || pr->reactants_[1] == "Photon") ) {
        if (pr->reactants_[0] == "Photon") {
          in_spec = pr->reactants_[1];
        } else {
          in_spec = pr->reactants_[0];
        }
      } else {
        std::stringstream msg; 
        msg << "### FATAL ERROR in ChemNetwork InitializeReactions() [ChemNetwork]"
           << std::endl << "Wrong format in FUV reaction ID=" << pr->id_ << std::endl;
      }
      //set indexing arrays
      if (pr->formula_ == 2) {
        inph_.push_back(ispec_map_[in_spec]);
        outph1_.push_back(ispec_map_[ pr->products_[0]]);
        outph2_.push_back(ispec_map_[ pr->products_[1]]);
        kph_base_.push_back(pr->alpha_);
        kph_avfac_.push_back(pr->gamma_);
        smap_ph_[in_spec] = n_ph_;
        kph_.push_back(0.);
        n_ph_++;
      } else {
        error = true;
      }
    //---------------- 4-8 - 2body reaction ---------------------------
    } else if (pr->itype_ == 4 || pr->itype_ == 5 || pr->itype_ == 6 
               || pr->itype_ == 7 || pr->itype_ == 8) {
      //check format
      if (pr->reactants_.size() != 2 || 
          (pr->products_.size() != 1 && pr->products_.size() != 2
           && pr->products_.size() != 3)) {
        std::stringstream msg; 
        msg << "### FATAL ERROR in ChemNetwork InitializeReactions() [ChemNetwork]"
            << std::endl << "Wrong format in 2body reaction ID=" << pr->id_
            << std::endl;
        ATHENA_ERROR(msg);
      }
      //set indexing arrays
      if (pr->formula_ == 3 || pr->formula_ == 4 || pr->formula_ == 5) {
        in2body1_.push_back(ispec_map_[ pr->reactants_[0]]);
        in2body2_.push_back(ispec_map_[ pr->reactants_[1]]);
        out2body1_.push_back(ispec_map_[ pr->products_[0]]);
        if (pr->products_.size() >= 2) {
          out2body2_.push_back(ispec_map_[ pr->products_[1]]);
        } else {
          out2body2_.push_back(-1);
        }
        if (pr->products_.size() == 3) {
          out2body3_.push_back(ispec_map_[ pr->products_[2]]);
        } else {
          out2body3_.push_back(-1);
        }
        frml_2body_.push_back(pr->formula_);
        a2body_.push_back(pr->alpha_);
        b2body_.push_back(pr->beta_);
        c2body_.push_back(pr->gamma_);
        k2body_.push_back(0.);
        n_2body_++;
      } else if (pr->formula_ == 7) {
        in2body1_.push_back(ispec_map_[ pr->reactants_[0]]);
        in2body2_.push_back(ispec_map_[ pr->reactants_[1]]);
        out2body1_.push_back(ispec_map_[ pr->products_[0]]);
        if (pr->products_.size() >= 2) {
          out2body2_.push_back(ispec_map_[ pr->products_[1]]);
        } else {
          out2body2_.push_back(-1);
        }
        if (pr->products_.size() == 3) {
          out2body3_.push_back(ispec_map_[ pr->products_[2]]);
        } else {
          out2body3_.push_back(-1);
        }
        frml_2body_.push_back(pr->formula_);
        a2body_.push_back(0.);
        b2body_.push_back(0.);
        c2body_.push_back(0.);
        k2body_.push_back(0.);
        idmap_2body_[pr->id_] = n_2body_;
        n_2body_++;
      } else {
        error = true;
      }

    //-------------------- 9 - grain assisted reaction ----------------
    } else if (pr->itype_ == 9) {
      //check format
      if (pr->reactants_.size() != 2 || pr->products_.size() != 1) {
        std::stringstream msg; 
        msg << "### FATAL ERROR in ChemNetwork InitializeReactions() [ChemNetwork]"
            << std::endl << "Wrong format in gr reaction ID=" << pr->id_ << std::endl;
        ATHENA_ERROR(msg);
      }
      if (pr->reactants_[1] != "e-" && pr->reactants_[1] != "H") {
        std::stringstream msg; 
        msg << "### FATAL ERROR in ChemNetwork InitializeReactions() [ChemNetwork]"
            << std::endl << "second reactant must be H or e-." << std::endl;
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

    //------------------ 10 - special reactions -----------------------
    } else if (pr->itype_ == 10) {
      if (pr->reactants_.size() > n_insr_ || pr->products_.size() > n_outsr_) {
        std::stringstream msg; 
        msg << "### FATAL ERROR in ChemNetwork InitializeReactions() [ChemNetwork]"
            << std::endl << "Wrong format in special reaction ID=" << pr->id_ 
            << std::endl;
      }
      if (pr->formula_ == 7) {
        idmap_sr_[pr->id_] = n_sr_;
        ksr_.push_back(0.);
        pr_sr_.push_back(pr);
        n_sr_++;
      } else {
        error = true;
      }

    //------------------ type not recogonized -------------------------
    } else{
      error = true;
    }

    //special reaction coefficients
    if (n_sr_ > 0) {
      insr_.NewAthenaArray(n_sr_, n_insr_);
      outsr_.NewAthenaArray(n_sr_, n_outsr_);
    }
    for (int i=0; i<n_sr_; i++) {
      for (int jin=0; jin<n_insr_; jin++) {
        if (jin < pr_sr_[i]->reactants_.size()) {
          insr_(i, jin) = ispec_map_[pr_sr_[i]->reactants_[jin]];
        } else {
          insr_(i, jin) = -1;
        }
      }
      for (int jout=0; jout<n_outsr_; jout++) {
        if (jout < pr_sr_[i]->products_.size()) {
          outsr_(i, jout) = ispec_map_[pr_sr_[i]->products_[jout]];
        } else {
          outsr_(i, jout) = -1;
        }
      }
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
  const Real y_H2 = y[ispec_map_["H2"]];
  const Real y_H = y[ispec_map_["H"]];
  const Real y_e = y[ispec_map_["e-"]];
  Real T;
  if (NON_BAROTROPIC_EOS) {
    T = E / Thermo::CvCold(y_H2, xHe_, y_e);
  } else {
    //isohermal EOS
    T = temperature_;
  }
	//cosmic ray reactions
	for (int i=0; i<n_cr_; i++) {
		kcr_[i] = kcr_base_[i] * rad_(index_cr_);
	}

	//cosmic ray induced photo reactions
	for (int i=0; i<n_crp_; i++) {
		kcrp_[i] = kcrp_base_[i] * rad_(index_cr_) * 2*y_H2;
	}

	//FUV reactions
	for (int i=0; i<n_ph_; i++) {
    kph_[i] = kph_base_[i] * rad_(i);
	}

  //2body reactions
	for (int i=0; i<n_2body_; i++) {
    if (frml_2body_[i] == 3) {
      k2body_[i] = a2body_[i]*pow(T/300., b2body_[i])*exp(-c2body_[i]/T) * nH_;
    } else if (frml_2body_[i] == 4) {
      k2body_[i] = a2body_[i]*b2body_[i]*( 0.62 
                                          + 0.4767*c2body_[i]*sqrt(300./T) ) * nH_;
    } else if (frml_2body_[i] == 5) {
      k2body_[i] = a2body_[i]*b2body_[i]*( 1 + 0.0967*c2body_[i]*sqrt(300./T) 
                                            + 28.501*c2body_[i]*c2body_[i]/T ) * nH_;
    }
  }

  //special rates and grain assisted reactions
  UpdateRatesSpecial(y, E);
  return;
}

void ChemNetwork::CheckReaction(KidaReaction reaction) {
  int atom_count_in[KidaSpecies::natom_];
  int atom_count_out[KidaSpecies::natom_];
  int charge_in = 0;
  int charge_out = 0;
  for (int ia=0; ia<KidaSpecies::natom_; ia++) {
    atom_count_in[ia] = 0;
    atom_count_out[ia] = 0;
  }
  for (int i=0; i<reaction.reactants_.size(); i++) {
    if (reaction.reactants_[i] == "CR" || reaction.reactants_[i] == "CRP" 
        || reaction.reactants_[i] == "Photon") {
      continue;
    }
    for (int ia=0; ia<KidaSpecies::natom_; ia++) {
      atom_count_in[ia] +=
        species_[ispec_map_[reaction.reactants_[i]]].atom_count_[ia];
    }
    charge_in += species_[ispec_map_[reaction.reactants_[i]]].charge_;
  }

  for (int i=0; i<reaction.products_.size(); i++) {
    for (int ia=0; ia<KidaSpecies::natom_; ia++) {
      atom_count_out[ia] +=
        species_[ispec_map_[reaction.products_[i]]].atom_count_[ia];
    }
    charge_out += species_[ispec_map_[reaction.products_[i]]].charge_;
  }

  if (charge_in != charge_out) {
    reaction.Print();
    std::stringstream msg; 
    msg << "### FATAL ERROR in ChemNetwork CheckReaction() [ChemNetwork] :"
        << "charge not conserved." << std::endl;
    ATHENA_ERROR(msg);
  }

  for (int ia=0; ia<KidaSpecies::natom_; ia++) {
    if (atom_count_in[ia] != atom_count_out[ia]) {
      reaction.Print();
      std::stringstream msg; 
      msg << "### FATAL ERROR in ChemNetwork CheckReaction() [ChemNetwork] :"
          << "atoms not conserved." << std::endl;
      ATHENA_ERROR(msg);
    }
  }
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
    reactions_[i].Print();
    std::cout << "alpha=" << reactions_[i].alpha_ << "," 
              << "beta=" << reactions_[i]. beta_ << "," 
              << "gamma=" << reactions_[i].gamma_ << "," 
              << "itype=" << reactions_[i].itype_ << ","
              << "forumla=" << reactions_[i].formula_ << std::endl;
  }

  //print reaction coefficients
  //cosmic-ray reactions
  std::cout << "CR reations:" << std::endl;
  for (int i=0; i<n_cr_; i++) {
    std::cout<< species_names[incr_[i]] << " + CR -> "
      << species_names[outcr1_[i]] << " + " << species_names[outcr2_[i]] << ", "
      << "kcr_base_=" << kcr_base_[i] << std::endl;
  }
  std::cout << "idmap_cr_: " << std::endl;
  for (std::map<int,int>::const_iterator it=idmap_cr_.begin();
       it!=idmap_cr_.end(); it++) {
    std::cout << it->first << " => " << it->second << std::endl;
  }

  //cosmic-ray induced photo reactions
  std::cout << "CRP reations:" << std::endl;
  for (int i=0; i<n_crp_; i++) {
    std::cout<< species_names[incrp_[i]] << " + CRP -> "
      << species_names[outcrp1_[i]] << " + " << species_names[outcrp2_[i]] << ", "
      << "kcrp_base_=" << kcrp_base_[i] << std::endl;
  }
  std::cout << "idmap_crp_: " << std::endl;
  for (std::map<int,int>::const_iterator it=idmap_crp_.begin();
       it!=idmap_crp_.end(); it++) {
    std::cout << it->first << " => " << it->second << std::endl;
  }

  //FUV reactions
  std::cout << "FUV photo- ionization/dissociation:" << std::endl;
  for (int i=0; i<n_ph_; i++) {
    std::cout<< species_names[inph_[i]] << " + Photon -> "
      << species_names[outph1_[i]] << " + " << species_names[outph2_[i]] << ", "
      << "kph_base_=" << kph_base_[i] << ", kph_avfac_=" << kph_avfac_[i]
      << std::endl;
  }
  std::cout << "smap_ph_: " << std::endl;
  for (std::map<std::string,int>::const_iterator it=smap_ph_.begin();
       it!=smap_ph_.end(); it++) {
    std::cout << it->first << " => " << it->second << std::endl;
  }

  //2body reactions
  std::cout << "2body reactions:" << std::endl;
  for (int i=0; i<n_2body_; i++) {
    std::cout<< species_names[in2body1_[i]] << " + "
      << species_names[in2body2_[i]]<< " -> " << species_names[out2body1_[i]];
    if (out2body2_[i] >= 0) {
      std::cout<< " + " << species_names[out2body2_[i]];
    }
    if (out2body3_[i] >= 0) {
      std::cout<< " + " << species_names[out2body3_[i]];
    }
    std::cout<< ", " << "alpha=" << a2body_[i] << ", beta=" << b2body_[i]
      << ", gamma=" << c2body_[i] << std::endl;
  }
  std::cout << "idmap_2body_: " << std::endl;
  for (std::map<int,int>::const_iterator it=idmap_2body_.begin();
       it!=idmap_2body_.end(); it++) {
    std::cout << it->first << " => " << it->second << std::endl;
  }

  //grain assisted reactions
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

  //special reactions
  std::cout << "special reations:" << std::endl;
  for (int i=0; i<n_sr_; i++) {
    for (int jin=0; jin<n_insr_; jin++) {
      if (insr_(i, jin) >= 0) {
        std::cout << species_names[insr_(i, jin)];
        if (jin < n_insr_-1 && insr_(i, jin+1) >= 0) {
          std::cout << " + ";
        }
      }
    }
    std::cout << " -> ";
    for (int jout=0; jout<n_outsr_; jout++) {
      if (outsr_(i, jout) >= 0) {
        std::cout << species_names[outsr_(i, jout)];
        if (jout < n_outsr_-1 && outsr_(i, jout+1) >= 0) {
          std::cout << " + ";
        }
      }
    }
    std::cout << std::endl;
  }
  std::cout << "idmap_sr_: " << std::endl;
  for (std::map<int,int>::const_iterator it=idmap_sr_.begin();
       it!=idmap_sr_.end(); it++) {
    std::cout << it->first << " => " << it->second << std::endl;
  }

  return;
}

void ChemNetwork::OutputRates(FILE *pf) const {
  //output the reactions and rates
	for (int i=0; i<n_cr_; i++) {
		fprintf(pf, "%4s + CR -> %4s + %4s,     kcr = %.2e\n", 
		 species_names[incr_[i]].c_str(), species_names[outcr1_[i]].c_str(),
     species_names[outcr2_[i]].c_str(), kcr_[i]);
	}
	for (int i=0; i<n_crp_; i++) {
		fprintf(pf, "%4s + CRP -> %4s + %4s,     kcrp = %.2e\n", 
		 species_names[incrp_[i]].c_str(), species_names[outcrp1_[i]].c_str(),
     species_names[outcrp2_[i]].c_str(), kcrp_[i]);
	}
	for (int i=0; i<n_ph_; i++) {
		fprintf(pf, "%4s + Photon -> %4s + %4s,     kph = %.2e\n", 
		 species_names[inph_[i]].c_str(), species_names[outph1_[i]].c_str(),
     species_names[outph2_[i]].c_str(), kph_[i]);
	}
	for (int i=0; i<n_2body_; i++) {
    fprintf(pf, "%4s + %4s -> %4s",
        species_names[in2body1_[i]].c_str(),
        species_names[in2body2_[i]].c_str(),
        species_names[out2body1_[i]].c_str());
    if (out2body2_[i] >= 0) {
      fprintf(pf, " + %4s", species_names[out2body2_[i]].c_str());
    }
    if (out2body3_[i] >= 0) {
      fprintf(pf, " + %4s", species_names[out2body3_[i]].c_str());
    }
    fprintf(pf,   ",     k2body = %.2e\n", k2body_[i]);
	}
	for (int i=0; i<n_gr_; i++) {
		fprintf(pf, "%4s + %4s (+ gr) -> %4s (+ gr),       kgr = %.2e\n", 
		 species_names[ingr1_[i]].c_str(), species_names[ingr2_[i]].c_str(),
     species_names[outgr_[i]].c_str(), kgr_[i]);
	}
  for (int i=0; i<n_sr_; i++) {
    for (int jin=0; jin<n_insr_; jin++) {
      if (insr_(i, jin) >= 0) {
        fprintf(pf, "%4s", species_names[insr_(i, jin)].c_str());
        if (jin < n_insr_-1 && insr_(i, jin+1) >= 0) {
          fprintf(pf, " + ");
        }
      }
    }
    fprintf(pf, " -> ");
    for (int jout=0; jout<n_outsr_; jout++) {
      if (outsr_(i, jout) >= 0) {
        fprintf(pf, "%4s", species_names[outsr_(i, jout)].c_str());
        if (jout < n_outsr_-1 && outsr_(i, jout+1) >= 0) {
          fprintf(pf, " + ");
        }
      }
    }
		fprintf(pf, ",       ksr = %.2e\n", ksr_[i]);
  }
  return;
}

void ChemNetwork::InitializeNextStep(const int k, const int j, const int i) {
  Real rho, rho_floor, rad_sum;
  const int nang = pmy_mb_->prad->nang;
  //density
  rho = pmy_mb_->phydro->w(IDN, k, j, i);
  //apply density floor
  rho_floor = pmy_mb_->peos->GetDensityFloor();
  rho = (rho > rho_floor) ?  rho : rho_floor;
  //hydrogen atom number density
  nH_ =  rho * unit_density_in_nH_;
  //average radiation field of all angles
  for (int ifreq=0; ifreq < n_freq_; ++ifreq) {
    rad_sum = 0;
    //radiation
    for (int iang=0; iang < nang; ++iang) {
      rad_sum += pmy_mb_->prad->ir(k, j, i, ifreq * nang + iang);
    }
    if (ifreq == index_cr_) {
      rad_(index_cr_) = rad_sum / float(nang);
    } else {
      rad_(ifreq) = rad_sum * unit_radiation_in_draine1987_ / float(nang) ;
    }
#ifdef DEBUG
    if (isnan(rad_(ifreq))) {
      printf("InitializeNextStep: ");
      printf("ifreq=%d, nang=%d, rad_sum=%.2e\n", ifreq, nang, rad_sum);
      OutputRates(stdout);
    }
#endif
  }
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
      printf("rad_ = ");
      for (int ifreq=0; ifreq < n_freq_; ++ifreq) {
        printf("%.2e  ", rad_(ifreq));
      }
      printf("\n");
      printf("nH_ = %.2e\n", nH_);
      std::stringstream msg;
      msg << "ChemNetwork (kida): RHS(y): nan or inf" << std::endl;
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
    rate = kcr_[i] * y[incr_[i]];
    ydotg[incr_[i]] -= rate;
    ydotg[outcr1_[i]] += rate;
    ydotg[outcr2_[i]] += rate;
  }

  //cosmic ray induced photo reactions
  for (int i=0; i<n_crp_; i++) {
    rate = kcrp_[i] * y[incrp_[i]];
    ydotg[incrp_[i]] -= rate;
    ydotg[outcrp1_[i]] += rate;
    ydotg[outcrp2_[i]] += rate;
  }

  //cosmic ray induced photo reactions
  for (int i=0; i<n_ph_; i++) {
    rate = kph_[i] * y[inph_[i]];
    ydotg[inph_[i]] -= rate;
    ydotg[outph1_[i]] += rate;
    ydotg[outph2_[i]] += rate;
  }

  //2body reactions
  for (int i=0; i<n_2body_; i++) {
    rate =  k2body_[i] * y[in2body1_[i]] * y[in2body2_[i]];
    if (y[in2body1_[i]] < 0 && y[in2body2_[i]] < 0) {
      rate *= -1.;
    }
    ydotg[in2body1_[i]] -= rate;
    ydotg[in2body2_[i]] -= rate;
    ydotg[out2body1_[i]] += rate;
    if (out2body2_[i] >= 0) {
      ydotg[out2body2_[i]] += rate;
    }
    if (out2body3_[i] >= 0) {
      ydotg[out2body3_[i]] += rate;
    }
  }

  //grain assisted reactions
  for (int i=0; i<n_gr_; i++) {
    rate = kgr_[i] * y[ingr1_[i]];
    ydotg[ingr1_[i]] -= rate;
    ydotg[ingr2_[i]] -= rate;
    ydotg[outgr_[i]] += rate;
  }

  //special reactions
  for (int i=0; i<n_sr_; i++) {
    rate = ksr_[i];
    for (int jin=0; jin<n_insr_; jin++) {
      if (insr_(i, jin) >= 0) {
        ydotg[insr_(i, jin)] -= rate;
      }
    }
    for (int jout=0; jout<n_outsr_; jout++) {
      if (outsr_(i, jout) >= 0) {
        ydotg[outsr_(i, jout)] += rate;
      }
    }
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
      printf("rad_ = ");
      for (int ifreq=0; ifreq < n_freq_; ++ifreq) {
        printf("%.2e  ", rad_(ifreq));
      }
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


//default: no special rates
void __attribute__((weak)) ChemNetwork::UpdateRatesSpecial(const Real y[NSCALARS],
                                                           const Real E) {
  // do nothing
  return;
}
