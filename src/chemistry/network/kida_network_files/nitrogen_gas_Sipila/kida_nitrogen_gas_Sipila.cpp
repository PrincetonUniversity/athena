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
//! \file kida_gow17.cpp
//  \brief implementation of special rates in gow17 network in kida format
//======================================================================================

// this class header
#include "../../kida.hpp"

//athena++ header
#include "../../network.hpp"
#include "../../../utils/chemistry_utils.hpp"
#include "../../../utils/thermo.hpp"
#include "../../../../defs.hpp"

//c++ header
#include <iostream>   // endl
#include <sstream>    // stringstream

static bool check_index = true;
Real CII_rec_rate(const Real temp);

void ChemNetwork::UpdateRatesSpecial(const Real y[NSCALARS], const Real E) {
  //constant or evolve temperature
  const Real y_H2 = y[ispec_map_["H2"]];
  const Real y_H = y[ispec_map_["H"]];
  const Real y_e = y[ispec_map_["e-"]];
  const Real small = 1e-50;
  Real T;
  if (NON_BAROTROPIC_EOS) {
    T = E / Thermo::CvCold(y_H2, xHe_, y_e);
  } else {
    //isohermal EOS
    T = temperature_;
  }
	//cap T above some minimum temperature
	if (T < temp_min_rates_) {
		T = temp_min_rates_;
	} 

  //2body reactions
  const int ns_2body = 1;
  const int indices_2body[ns_2body] = {17};
  //(17) H2 + N+ -> H + NH+      -- depends on o/p, Dislaire et al. (2012)
  const Real x_oH2 = o2pH2_ / (1. + o2pH2_);
  const Real ko = 4.2e-10 * pow(T/300., -0.17) * exp(-44.5/T);
  const Real kp = 8.35e-10 * exp(-168.5/T);
  k2body_(id7map_(17)) = ( x_oH2*ko + (1.-x_oH2)*kp ) * nH_;

  //sanity check
  if (check_index) {
    std::stringstream msg; //error message
    msg << "### FATAL ERROR in ChemNetwork UpdateRatesSpecial [ChemNetwork]: "
      << "Specital rate reaction type does not match" << std::endl;
    for (int i=0; i<ns_2body; i++) {
      if (id7type_(indices_2body[i]) != ReactionType::twobody) {
        ATHENA_ERROR(msg);
      }
    }
    check_index = false;
  }
  return;
}
