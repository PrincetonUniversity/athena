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
//! \file G14ShockTube.cpp
//  \brief implementation of functions in Grassi 2014 Fig. 21
//======================================================================================


// this class header
#include "G14Sod.hpp"

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
#include <sstream>   // stringstream
#include <iostream>  // endl
#include <math.h>    // a^x = pow(a,x)
#include <stdio.h>   // FILE, fprintf()
#include <limits>    // inf


//species namesspecies_names
const std::string ChemNetwork::species_names[NSPECIES] = 
{"H","H+","He","He+","He2+","H-","H2","H2+"};

const std::string ChemNetwork::ghost_species_names_[ngs_] = 
{"e"};


//index of species
const int ChemNetwork::iH_ =
  ChemistryUtility::FindStrIndex(species_names, NSPECIES, "H");
const int ChemNetwork::iHplus_ =
  ChemistryUtility::FindStrIndex(species_names, NSPECIES, "H+");
const int ChemNetwork::iHe_ =
  ChemistryUtility::FindStrIndex(species_names, NSPECIES, "He");
const int ChemNetwork::iHeplus_ =
  ChemistryUtility::FindStrIndex(species_names, NSPECIES, "He+");
const int ChemNetwork::iHe2plus_ =
  ChemistryUtility::FindStrIndex(species_names, NSPECIES, "He2+");
const int ChemNetwork::iHmin_ =
  ChemistryUtility::FindStrIndex(species_names, NSPECIES, "H-");
const int ChemNetwork::iH2_ =
  ChemistryUtility::FindStrIndex(species_names, NSPECIES, "H2");
const int ChemNetwork::iH2plus_ =
  ChemistryUtility::FindStrIndex(species_names, NSPECIES, "H2+");


//index of ghost species _> change the netwrok definition?
const int ChemNetwork::ige_ =
      ChemistryUtility::FindStrIndex(ghost_species_names_, ngs_, "e") + NSPECIES;
const int  ChemNetwork::igr_ = ige_;

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

//Note: output to ghost species doesn't matter. The abundances of ghost species
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
	//number of species and a list of name of species
  pmy_spec_ = pmb->pscalars;
	pmy_mb_   = pmb;

  //units of density and radiation
	unit_density_in_nH_ = pin->GetReal("chemistry", "unit_density_in_nH");
	unit_length_in_cm_ = pin->GetReal("chemistry", "unit_length_in_cm");
	unit_vel_in_cms_ = pin->GetReal("chemistry", "unit_vel_in_cms");
	unit_time_in_s_ = unit_length_in_cm_/unit_vel_in_cms_;
  // Copy from G14 setup doc
	unit_E_in_cgs_ = 1.67e-24  * mu* unit_density_in_nH_ 
                    * unit_vel_in_cms_ * unit_vel_in_cms_;

  if (NON_BAROTROPIC_EOS) {
    temperature_ = 0.;
  } else {
    //isothermal
    temperature_ = pin->GetReal("chemistry", "temperature");
  }
    //temperature above or below which heating and cooling is turned off
    Real inf = std::numeric_limits<Real>::infinity();
    temp_max_heat_ = pin->GetOrAddReal("chemistry", "temp_max_heat", inf);
    temp_min_cool_ = pin->GetOrAddReal("chemistry", "temp_min_cool", 1.);
    //minimum temperature for reaction rates, also applied to energy equation
    temp_min_rates_ = pin->GetOrAddReal("chemistry", "temp_min_rates", 1.);

  //initialize rates to zero
  for (int i=0; i<n_2body_; i++) {
    k2body_[i] = 0;
  }

  //copy species to a full list of species names
  for (int i=0; i<NSPECIES; i++) {
    species_names_all_[i] = species_names[i];
  }
  for (int i=NSPECIES; i<NSPECIES+ngs_; i++) {
    species_names_all_[i] = ghost_species_names_[i-NSPECIES];
  }

}
ChemNetwork::~ChemNetwork() {}

//-----------------end of chemical network---------------------
void ChemNetwork::RHS(const Real t, const Real y[NSPECIES], const Real ED, 
                      Real ydot[NSPECIES]) {
  //function of evolution of the abundance of the element
  Real rate;
  //store previous y includeing ghost species
  Real yprev[NSPECIES+ngs_];
  Real yprev0[NSPECIES+ngs_]; //correct negative abundance
  Real ydotg[NSPECIES+ngs_];
  Real E_ergs = ED * unit_E_in_cgs_ / rho; //ernergy per hydrogen atom
  Real Xe = y[iHplus_] + y[iHeplus_] + 2.*y[iHe2plus_] + y[iH2plus_] - y[iHmin_];
  // tgas = p / rho / kboltzmann * mu * pmass
  Real T = E_ergs/Thermo::kb_*(gamma-1);
  // copy y to yprev and set ghost species
  GetGhostSpecies(y, yprev);

  for(int i=0; i<NSPECIES+ngs_; i++) {
	ydotg[i] = 0.0;
	if (yprev[i] < 0) {
      yprev0[i] = 0;
    } else {
      yprev0[i] = yprev[i];
    }
    //throw error if nan, or inf, or large negative value occurs
    if ( isnan(yprev[i]) || isinf(yprev[i]) ) {
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

  //2body reactions
  for (int i=0; i<n_2body_; i++) {
    rate =  k2body_[i] * yprev0[in2body1_[i]] * yprev0[in2body2_[i]];
    if (yprev0[in2body1_[i]] < 0 && yprev0[in2body2_[i]] < 0) {
      rate *= -1.;
    }
#ifdef DEBUG
    if (indi == indj == indk == 1){
    if (species_names_all_[in2body1_[i]] == "H2" || species_names_all_[in2body2_[i]]== "H2"){
    	printf("Dec: k2body_[%.2i]*(X_%s)*(X_%s) = %.2e* %.2e * %.2e = %.2e \n",
    		i, species_names_all_[in2body1_[i]].c_str(), species_names_all_[in2body2_[i]].c_str() ,k2body_[i] ,yprev0[in2body1_[i]], yprev0[in2body2_[i]],rate);
    }
    if (species_names_all_[out2body1_[i]] == "H2" || species_names_all_[out2body2_[i]]== "H2"){
//        printf("%s + %s , k2body_[%.2i]\n",species_names_all_[out2body1_[i]].c_str(), species_names_all_[out2body2_[i]].c_str(),i);
    	printf("Inc: k2body_[%.2i]*(X_%s)*(X_%s) = %.2e* %.2e * %.2e = %.2e \n",
    		i, species_names_all_[in2body1_[i]].c_str(), species_names_all_[in2body2_[i]].c_str() ,k2body_[i] ,yprev0[in2body1_[i]], yprev0[in2body2_[i]],rate);
    }
	}

#endif
    ydotg[in2body1_[i]]  -= stoich_in2body1[i]*rate;
    ydotg[in2body2_[i]]  -= stoich_in2body2[i]*rate;
    ydotg[out2body1_[i]] += stoich_out2body1[i]*rate;
    ydotg[out2body2_[i]] += stoich_out2body2[i]*rate;
  }
   //set ydot to return
   for (int i=0; i<NSPECIES; i++) {
     //return in code units
     ydot[i] = ydotg[i] * unit_time_in_s_ * nH_;
   }
#ifdef DEBUG
  OutputRates(stdout);
  printf("nH_ = %2.e \n",nH_);
  printf("one loop finish, ydot, ydotg yprev = \n");
  for (int j=0; j<NSPECIES+ngs_; j++) {
    printf("%s: %.2e , %.2e, %.2e ", species_names_all_[j].c_str(), ydot[j],ydotg[j],yprev[j]);
    printf("\n");
  	}
#endif
  return;
}

void ChemNetwork::GetGhostSpecies(const Real *y, Real yghost[NSPECIES+ngs_]) {
  //copy the aboundances in y to yghost
  for (int i=0; i<NSPECIES; i++) {
    yghost[i] = y[i];
  }
  //set the ghost species
  yghost[ige_] = y[iHplus_] + y[iHeplus_] + 2.*y[iHe2plus_] + y[iH2plus_] - y[iHmin_];

return;
}

void ChemNetwork::UpdateRates(const Real y[NSPECIES+ngs_], const Real E) {
  // The following rate coff. Ziegler is in the SI unit, converting it back to cgs.
  const Real rate_in_cgs_ = 1e6;
  Real T = E/Thermo::kb_*(gamma - 1.0);

  const Real logT    = log10(T);
  const Real lnTe    = log(T* 8.6163e-5);
  const Real Te      = T* 8.6163e-5;

  //(1) H    + e   -> H+   + 2e
  k2body_[0]  =  exp(-3.271396786e1 + 
                    ( 1.35365560e1  + (-5.73932875e0 + (1.56315498e0   +
                    (-2.877056e-1   + (3.48255977e-2 + (-2.63197617e-3 +
                    ( 1.11954395e-4 + (-2.03914985e-6)
                    *lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)
                    *lnTe);
  //  (2) H+   + e   -> H    (+ γ)
  if (T<=5500){
    k2body_[1]  =  3.92e-13*pow(Te,-0.6353);
  }else{
  	k2body_[1]  =  exp(-2.861303380689232e1 +
  		              (-7.241125657826851e-1 + (-2.026044731984691e-2 + (-2.380861877349834e-3 +
  		              (-3.212605213188796e-4 + (-1.421502914054107e-5 + ( 4.989108920299513e-6 +
  		              ( 5.755614137575758e-7 + (-1.856767039775261e-8 + (-3.071135243196595e-9 )
  		              *lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)
  		              *lnTe)*lnTe);
  }
  //  (3) He   + e   -> He+  + 2e [unit in cgs]
  k2body_[2]  = exp(-4.409864886e1 +
  	               ( 2.391596563e1 + (-1.07532302e1  + (+3.05803875e0 +
  	               (-5.6851189e-1  + ( 6.79539123e-2 + (-5.0090561e-3 +
  	               ( 2.06723616e-4 + (-3.64916141e-6)
  	               	*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe);

  //  (4) He+  + e   -> He   (+ γ)  k_4 = k_4a + k_4b unit in cgs
  k2body_[3]  = 3.92e-13*pow(Te,-0.6353);
  if (T>9280){

  	k2body_[3]  += 1.54e-9*pow(Te,-1.5)*(1.0+0.3/exp(8.099328789667/Te))/exp(40.49664394833662/Te);

  }

  //  (5) He+  + e   -> He2+ + 2e  unit in cgs
  k2body_[4]  = exp(-6.8710409e1  +
                  ( 4.393347632635e1   + (-1.848066993568e1  + (+4.701626486759e0  +
                  (-7.692466334492e-1  + (8.113042097303e-2  + (-5.324020628287e-3 +
                  (+1.975705312221e-4  + (-3.165581065665e-6)
                  *lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe);

  //  (6) He2+ + e   -> He+  (+ γ) unit in cgs
  k2body_[5]  = 3.36e-10*pow(T,-0.5)*pow(T/1e3,-0.2)*pow(1+pow(T/1e6,0.7),-1);

  //  (7) H    + e   -> H−   (+ γ)
  k2body_[6]  = 6.77e-15*pow(Te,0.8779);

  //  (8) H−   + H   -> H2   + e unit in cgs
  if (T<=1160){
  	k2body_[7]  = 1.43e-9;
  }else{
  	k2body_[7]  = exp(-2.00691389e1  +
                    ( 2.289800603272916e-1  + ( 3.599837721023835e-2  + (-4.555120027032095e-3 +
                    (-3.105115447124016e-4  + ( 1.073294010367247e-4  + (-8.36671960467864e-6  +
                    ( 2.238306228891639e-7)
                    *lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe);
  }

  //  (9) H    + H+  -> H2+  (+ γ)
  if (T<=6700){
  	k2body_[8]  = 1.85e-23*pow(T,1.8);
  }else{
  	k2body_[8]  = 5.81e-16*pow((T/56200), -0.6657*log10(T/56200));
  }
  // (10) H    + H2+ -> H2   + H+ unit in cgs
  k2body_[9]  = 6e-10;

  // (11) H+   + H2  -> H2+  + H
  if (T<=3480){
  	k2body_[10] = 0;
  }else{
  	k2body_[10] = exp(-2.424914687731536e1  +
                    ( 3.400824447095291e0   + (-3.898003964650152e0  + ( 2.045587822403071e0  +
                    (-5.416182856220388e-1  + ( 8.41077503763412e-2  + (-7.879026154483455e-3 +
                    ( 4.138398421504563e-4  + (-9.36345888928611e-6)
                    *lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe);
  }
  
  // (12) H2   + e   -> 2 H  + e
  k2body_[11] = 5.6e-11*sqrt(T)*exp(-102124.0/T);

  // (13) H−   + e   -> H    + 2 e
  k2body_[12] = exp(-1.801849334273e1  +
                   ( 2.360852208681e0     + (-2.827443061704e-1 + ( 1.623316639567-2  +
                   (-3.365012031362999e-2 + ( 1.178329782711e-2 + (-1.656194699504e-3 +
                   ( 1.068275202678e-4    + (-2.631285809207e-6)
                   *lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe);

  // (14) H−   + H   -> 2 H  + e
  if (T<=1160){
  	k2body_[13]  = 2.56e-9*pow(Te,1.78186);
  }else{
  	k2body_[13]  =  exp(-2.037260896203573e1 +
  		              (+1.139449335841631e0 + (-1.421013521554148e-1 + ( 8.46445538663e-3 +
  		              (-1.4327641212992e-3  + ( 2.012250284791e-4    + (+8.66396324309e-5 + 
  		              (-2.585009680264e-5   + (+2.4555011970392e-6   + (-8.06838246118e-8 )
  		              *lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)
  		              *lnTe)*lnTe);

  }
  // (15) H−   + H+  -> 2H  (+ γ)
  k2body_[14] = 6.5e-9/sqrt(Te);

  // (16) H−   + H+  -> H2+  + e
  k2body_[15] = 1e-8*pow(T,-0.4);

  // (17) H2+  + e   -> 2H  (+ γ)
  if (T<=617){
    k2body_[16] = 1e-8;
  }else{
  	k2body_[16] = 1.32e-6*pow(T,-0.76);
  }

  // (18) H2+  + H-   -> H    + H2
  k2body_[17]  = 5e-7*sqrt(1.0e2/T);

  // (19) 3H → H2 + H
  if (T<=300){
    k2body_[18] = 1.3e-32*pow(T/300,-0.38);
  }else{
    k2body_[18] = 1.3e-32*pow(T/300,-1.0);
  }

  // (20) H2   + H   -> 3 H
  k2body_[19] = 1.067e-10*pow(Te,2.012)*exp(-4.463/Te)*pow(1. + 0.2472*Te,-3.512);

  return;
}

void ChemNetwork::InitializeNextStep(const int k, const int j, const int i) {
  //density
  rho = pmy_mb_->phydro->w(IDN, k, j, i);
  //hydrogen atom number density
  nH_ =  rho * unit_density_in_nH_/mu;
	return;
}

Real ChemNetwork::Edot(const Real t, const Real y[NSPECIES], const Real ED){
  //function of evolution of energy
  //return dEdt;
  Real LH2; // Define H2 Cooling Term
  Real E_ergs = ED * unit_E_in_cgs_ / rho; //ernergy per hydrogen atom

  // Define Temperature
  Real T = E_ergs/Thermo::kb_*(gamma - 1.0);
  Real dEdt = 0.;
  Real yprev[NSPECIES+ngs_];
  // copy y to yprev and set ghost species
  GetGhostSpecies(y, yprev);
  //correct negative abundance to zero
  for (int i=0; i<NSPECIES+ngs_; i++) {
    if (yprev[i] < 0) {
      yprev[i] = 0;
    }
  }

  // H2 CE(ro-vibrational) Cooling
  LH2  = Thermo::CoolingH2(yprev[iH2_], nH_*yprev[iH_],  nH_*yprev[iH2_],
          nH_*yprev[iHe_],  nH_*yprev[iHplus_], nH_*yprev[ige_],
          T);

  //return in code units
  Real dEDdt = - LH2* rho/ unit_E_in_cgs_ * unit_time_in_s_;

#ifdef DEBUG
  	printf("Cooling = %.2e \n", Cooling );
  	printf("T=%.4e, dEDdt=%.2e, E=%.2e, dEergsdt=%.2e, E_ergs=%.2e, nH=%.2e\n",
            T, dEDdt, ED, dEdt, E_ergs,nH_);
			for (int i=0; i<NSPECIES+ngs_; i++) {
			printf("%s: %.2e  ", species_names_all_[i].c_str(), yprev[i]);
		}
		printf("\n");
		printf("=============================\n");
#endif

  return dEDdt;

}

void ChemNetwork::OutputRates(FILE *pf) const {
  //output the reactions and base rates
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