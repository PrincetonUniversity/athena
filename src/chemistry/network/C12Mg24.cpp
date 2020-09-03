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
//! \file H2.cpp
//  \brief implementation of functions in class ChemNetwork, using the simple
//  network for H2 formation and destruction.
//======================================================================================

// this class header
#include "C12Mg24.hpp"

//athena++ header
#include "network.hpp"
#include "../../scalars/scalars.hpp"
#include "../../parameter_input.hpp"       //ParameterInput
#include "../../mesh/mesh.hpp"
#include "../../hydro/hydro.hpp"
#include "../utils/chemistry_utils.hpp"
#include "../../defs.hpp"
#include "../../eos/eos.hpp"

//c++ header
#include <sstream>    // stringstream
#include <iostream>   // endl
#include <limits>    //inf

#ifdef DEBUG
static bool output_rates = true;
#endif

//constants
// const Real gm1  = 1.666666666666667 - 1.0;   //not a good way to do this

//species names
const std::string ChemNetwork::species_names[NSCALARS] = 
{"C12", "Mg24"};

const int ChemNetwork::iC12_ =
  ChemistryUtility::FindStrIndex(species_names, NSCALARS, "C12");
const int ChemNetwork::iMg24_ =
  ChemistryUtility::FindStrIndex(species_names, NSCALARS, "Mg24");

ChemNetwork::ChemNetwork(MeshBlock *pmb, ParameterInput *pin) {
	//number of species and a list of name of species
  pmy_spec_ = pmb->pscalars;
	pmy_mb_ = pmb;
	//set the parameters from input file
  mn_ = 1.674920e-24;
  Q12_ = 2.2323e-5;
  k_ = 1.380658e-16;
  //units
  unit_density_in_cgs = pin->GetReal("chemistry","unit_density_in_cgs");
  unit_temp_in_K = pin->GetReal("chemistry","unit_temp_in_K");
  unit_edot_in_cgs = pin->GetReal("chemistry","unit_edot_in_cgs");
//  const Real gm1 = pmy_mb_->peos->GetGamma() - 1;
}

ChemNetwork::~ChemNetwork() {}

void ChemNetwork::InitializeNextStep(const int k, const int j, const int i) {
  Real rho, rho_floor;
  //density
  rho = pmy_mb_->phydro->w(IDN, k, j, i);
  //apply density floor
  rho_floor = pmy_mb_->peos->GetDensityFloor();
  rho = (rho > rho_floor) ?  rho : rho_floor;
  // density_ = unit_density_in_cgs * rho;
  density_ = rho;
  return;
}

void ChemNetwork::RHS(const Real t, const Real y[NSCALARS], const Real ED,
                      Real ydot[NSCALARS]) {
  //correct negative rates
  Real y_corr[NSCALARS];
  for (int i=0; i<NSCALARS; i++) {
    if (y[i] < 0) {
      y_corr[i] = 0;
    } else {
      y_corr[i] = y[i];
    }
    //throw error if nan, or inf, or large negative value occurs
//    if ( isnan(y[i]) || isinf(y[i]) ) {
//      printf("RHS: ");
//      for (int j=0; j<NSCALARS; j++) {
//        printf("%s: %.2e  ", species_names[j].c_str(), y[j]);
//      }
//      printf("\n");
//      OutputRates(stdout);
//      printf("Y_C12 = %.2e\n", y[0]);
//      std::stringstream msg;
//      msg << "ChemNetwork (C12Mg24): RHS(yprev): nan or inf" << std::endl;
//     // ATHENA_ERROR(msg);
//      }
    }                                                                                                 
  Real mu_ = (y_corr[0]*12 + y_corr[1]*24);
  const Real gm1 = pmy_mb_->peos->GetGamma() - 1;
  Real temp = ED*gm1*mu_*mn_/(density_*k_);
  Real T9 = temp*1e-9;
//  Real T9 = 5.0;
  Real TA9 = T9/(1+0.0396*T9);
//  const Real rate_C12 = -(12*mn_/Q12_)*unit_edot_in_cgs*pow(y[0],2)*density_*pow(T9,29);
  const Real rate_C12 = -(12*mn_/Q12_)*3.96e43*density_*pow(y_corr[0],2)*pow(TA9,5./6.)*exp(-84.165*pow(TA9,-1./3.)
    - 2.12*1e-3*pow(T9,3));
  const Real rate_Mg24 = -rate_C12;
  Real tnuc = ED/(rate_C12*Q12_/(12*mn_));
if ((-1*tnuc<1e-3) or (-1*rate_C12>0.1)) {
  printf("tnuc, ED, density, mu, T(1e9 K), rate_C12 =  %.2e, %.2e, %.2e, %.2e, %.2e, %.2e\n", tnuc, ED, density_, mu_, T9, rate_C12);
}// if (T9 > 10){
//     printf("HIGH T || ED, density, mu, T(1e9 K), rate_C12 =  %.2e, %.2e, %.2e, %.2e, %.2e\n", ED, density_, mu_, T9, rate_C12);
// }
// if (T9 < 1e-2){
//     printf("LOW T || ED, density, mu, T(1e9 K), rate_C12 =  %.2e, %.2e, %.2e, %.2e, %.2e\n", ED, density_, mu_, T9, rate_C12);
// }
//  else {printf("GOAL T || ED, density, mu, T(1e9 K), rate_C12 =  %.2e, %.2e, %.2e, %.2e, %.2e\n", ED, density_, mu_, T9, rate_C12);}
  ydot[0] = rate_C12;
  ydot[1] = rate_Mg24;

#ifdef DEBUG
  if (output_rates) {
    FILE *pf = fopen("network.dat", "w");
    OutputRates(pf);
    fclose(pf);
    output_rates = false;
  }
#endif

  return;
}

Real ChemNetwork::Edot(const Real t, const Real y[NSCALARS], const Real ED) {
  Real mu_ = (y[0]*12 + y[1]*24)*mn_;
  const Real gm1 = pmy_mb_->peos->GetGamma() - 1;
  Real temp = ED*gm1*mu_/(density_*k_);
  Real T9 = temp*1e-9;
  // if (temp<0){
  //   temp = 1e-5;
  // }
//  const Real dEDdt = unit_edot_in_cgs*pow(y[0],2)*density_*pow(T9,29);
  const Real dEDdt = 0.0;
  return dEDdt;
}

void ChemNetwork::OutputRates(FILE *pf) const {
  //output the reactions and base rates
  fprintf(pf,"rate_C12 = %4s",rate_C12);
  return;
}
