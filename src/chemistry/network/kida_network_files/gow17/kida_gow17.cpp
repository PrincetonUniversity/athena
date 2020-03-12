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
  const Real y_O = y[ispec_map_["O"]];
  const Real y_H3plus = y[ispec_map_["H3+"]];
  const Real y_Oplus = y[ispec_map_["O+"]];
  const Real y_Cplus = y[ispec_map_["C+"]];
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
  const Real logT = log10(T);
	const Real logT4 = log10(T/1.0e4);
	const Real lnTe = log(T * 8.6173e-5);
  const Real kcr_H_fac = 2.3*y_H2 + 1.5 * y_H; //ratio of total to primary rate
  //grain assisted reactions
  const Real GPE0 = rad_(index_gpe_);
	const Real psi_gr_fac = 1.7 * GPE0 * sqrt(T) / nH_; 
  const Real psi = psi_gr_fac / y_e;
  const Real cHp[7] = {12.25, 8.074e-6, 1.378, 5.087e2, 1.586e-2, 0.4723, 1.102e-5}; 
  const Real cCp[7] = {45.58, 6.089e-3, 1.128, 4.331e2, 4.845e-2, 0.8120, 1.333e-4};
  const Real cHep[7] = {5.572, 3.185e-7, 1.512, 5.115e3, 3.903e-7, 0.4956, 5.494e-7};
  const Real cSip[7] = {2.166, 5.678e-8, 1.874, 4.375e4, 1.635e-6, 0.8964, 7.538e-5};

	//cosmic ray reactions
  const int ns_cr = 2;
  const int indices_cr[ns_cr] = {1, 2};
  //(1) H + CR -> H+ + e-
  kcr_(idmap_(1)) = kcr_H_fac * rad_(index_cr_); 
  //(2) H2 + CR -> H2+ + e-
  kcr_(idmap_(2)) = 2. * kcr_H_fac * rad_(index_cr_); 

	// Grain assisted reactions
  const int ns_gr = 5;
  const int indices_gr[ns_gr] = {15, 16, 17, 18, 19};
	//(15) H + H + gr -> H2 + gr , from Draine book chapter 31.2 page 346, Jura 1975
	kgr_(idmap_(15)) = 3.0e-17 * nH_ * zdg_;
	//(16) H+ + e- + gr -> H + gr
  kgr_(idmap_(16)) = 1.0e-14 * cHp[0] / 
               (
                 1.0 + cHp[1]*pow(psi, cHp[2]) * 
                   (1.0 + cHp[3] * pow(T, cHp[4])
                                 *pow( psi, -cHp[5]-cHp[6]*log(T) ) 
                   ) 
                ) * nH_ * zdg_;
	//(17) C+ + e- + gr -> C + gr
  kgr_(idmap_(17)) = 1.0e-14 * cCp[0] / 
               (
                 1.0 + cCp[1]*pow(psi, cCp[2]) * 
                   (1.0 + cCp[3] * pow(T, cCp[4])
                                 *pow( psi, -cCp[5]-cCp[6]*log(T) ) 
                   ) 
                ) * nH_ * zdg_;
	//(18) He+ + e- + gr -> He + gr
  kgr_(idmap_(18)) = 1.0e-14 * cHep[0] / 
               (
                 1.0 + cHep[1]*pow(psi, cHep[2]) * 
                   (1.0 + cHep[3] * pow(T, cHep[4])
                                 *pow( psi, -cHep[5]-cHep[6]*log(T) ) 
                   ) 
                ) * nH_ * zdg_;
	//(19) Si+ + e- + gr -> Si + gr
  kgr_(idmap_(19)) = 1.0e-14 * cSip[0] / 
               (
                 1.0 + cSip[1]*pow(psi, cSip[2]) * 
                   (1.0 + cSip[3] * pow(T, cSip[4])
                                 *pow( psi, -cSip[5]-cSip[6]*log(T) ) 
                   ) 
                ) * nH_ * zdg_;

  //2body reactions
  const int ns_2body = 9;
  const int indices_2body[ns_2body] = {33, 36, 38, 40, 41, 42, 43, 49, 50};
  //(33) He+ + e- -> He  -- Case B
	k2body_(idmap_(33)) = 1e-11*pow(T, -0.5)*(11.19 + 
                             (-1.676 + (-0.2852 + 0.04433*logT) * logT )* logT) * nH_;
  //(36) C+ + e- -> C    -- Include RR and DR, Badnell2003, 2006.
  k2body_(idmap_(36)) = CII_rec_rate(T) * nH_;
  //(38) H2 + H2+ -> H + H3+
  k2body_(idmap_(38)) = 1.76e-9 * pow(T, 0.042) * exp(- T/46600.) * nH_; 
  //(40) H+ + e- -> H  -- Case B
	k2body_(idmap_(40)) = 2.753e-14 * pow( 315614.0 / T, 1.5) 
									 * pow(  1.0 + pow( 115188.0 / T, 0.407) , -2.242 ) * nH_;
  //Collisional dissociation, k>~1.0e-30 at T>~5e2.
  //(41) H2 + H -> H + H + H
  //(42) H2 + H2 -> H + H + H2
  //(43) H + e- -> H+ + e- + e-
  Real k9l, k9h, k10l, k10h, ncrH, ncrH2, div_ncr, ncr, n2ncr;
  const Real temp_coll = 7.0e2;
  if (T > temp_coll) {
    // Density dependent. See Glover+MacLow2007
  	k9l = 6.67e-12 * sqrt(T) * exp(-(1. + 63590./T)); 
    k9h = 3.52e-9 * exp(-43900.0 / T);
    k10l = 5.996e-30 * pow(T, 4.1881) / pow((1.0 + 6.761e-6 * T), 5.6881)  
            * exp(-54657.4 / T);
    k10h = 1.3e-9 * exp(-53300.0 / T); 
    ncrH = pow(10, (3.0 - 0.416 * logT4 - 0.327 * logT4*logT4));
    ncrH2 = pow(10, (4.845 - 1.3 * logT4 + 1.62 * logT4*logT4));
		div_ncr = y_H/ncrH + y_H2/ncrH2;
		if (div_ncr < small) {
			ncr = 1./ small;
		} else {
			ncr = 1. / div_ncr;
		}
    n2ncr = nH_ / ncr;
    k2body_(idmap_(41)) = pow(10, log10(k9h) *  n2ncr/(1. + n2ncr) 
                         + log10(k9l) / (1. + n2ncr)) * nH_;
    k2body_(idmap_(42)) = pow(10, log10(k10h) *  n2ncr/(1. + n2ncr) 
                         + log10(k10l) / (1. + n2ncr)) * nH_;
    k2body_(idmap_(43)) = exp( -3.271396786e1 + 
                      (1.35365560e1 + (- 5.73932875 + (1.56315498 
                    + (- 2.877056e-1 + (3.48255977e-2 + (- 2.63197617e-3
                    + (1.11954395e-4 + (-2.03914985e-6)
                       *lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe) *lnTe) * nH_;
  } else {
    k2body_(idmap_(41)) = 0.;
    k2body_(idmap_(42)) = 0.;
    k2body_(idmap_(43)) = 0.;
  }
  //(49) O + H+ -> H + O+
  //(50) H + O+ -> O + H+
  k2body_(idmap_(49)) = ( 1.1e-11 * pow(T, 0.517) + 4.0e-10 * pow(T, 6.69e-3) 
                             ) * exp(-227./T) * nH_;
  k2body_(idmap_(50)) = (4.99e-11* pow(T, 0.405) + 7.5e-10 * pow(T, -0.458) 
                             )* nH_;
  
  //special reactions
  const int ns_sr = 8;
  const int indices_sr[ns_sr] = {5, 20, 21, 22, 23, 24, 28, 29};
  //(5) CO + H + CR -> HCO+ + e-, CO + CR -> CO+ + e- and CO+ + H/H2 -> HCO+
  ksr_(idmap_(5)) = 6.52 * rad_(index_cr_) * y[ispec_map_["CO"]];
  //(20) C + H3+ e- -> H2 + CHx       --Vissapragada2016 new rates*/
  const Real n_kCHx = 2.31e-3;
  const Real A_kCHx = 1.04e-9;
  const Real c_kCHx[4] = {3.4e-8, 6.97e-9, 1.31e-7, 1.51e-4};
  const Real Ti_kCHx[4] = {7.62, 1.38, 2.66e1, 8.11e3};
  const Real t1_CHx = A_kCHx * pow( 300./T, n_kCHx);
  const Real t2_CHx = c_kCHx[0] * exp(-Ti_kCHx[0]/T) + c_kCHx[1] * exp(-Ti_kCHx[1]/T)
           + c_kCHx[2]*exp(-Ti_kCHx[2]/T) + c_kCHx[3] *exp(-Ti_kCHx[3]/T);
  ksr_(idmap_(20)) = (t1_CHx + pow(T, -1.5) * t2_CHx) * nH_ 
                          * y[ispec_map_["C"]] * y_H3plus;
  //--- H2O+ + e branching--
  //(21) O + H3+ e- -> H2 + OHx
  //(22) O + H3+ e- -> H2 + O + H
  //(23) H2 + O+ + e- -> H + OHx
  //(24) H2 + O+ + e- -> H + H + O
	Real h2oplus_ratio, fac_H2Oplus_H2, fac_H2Oplus_e;
	if (y_e < small) {
		h2oplus_ratio = 1.0e10;
	} else {
		h2oplus_ratio = 6e-10 * y_H2 / ( 5.3e-6 / sqrt(T) * y_e );
	}
  fac_H2Oplus_H2 = h2oplus_ratio / (h2oplus_ratio + 1.);
  fac_H2Oplus_e = 1. / (h2oplus_ratio + 1.);
  ksr_(idmap_(21)) = 1.99e-9 * pow(T, -0.190) * fac_H2Oplus_H2 * nH_
                         * y_O * y_H3plus;
  ksr_(idmap_(22)) = 1.99e-9 * pow(T, -0.190) * fac_H2Oplus_e * nH_
                         * y_O * y_H3plus;
  ksr_(idmap_(23)) = 1.6e-9 * fac_H2Oplus_H2 * nH_
                         * y_H2 * y_Oplus;
  ksr_(idmap_(24)) = 1.6e-9 * fac_H2Oplus_e * nH_
                         * y_H2 * y_Oplus;
  //(28) H2 + C+ + e- -> CHx + H        -- schematic reaction for C+ + H2 -> CH2+
  //(29) H2 + C+ + e- -> C + H + H      -- schematic reaction for C+ + H2 -> CH2+
  ksr_(idmap_(28)) = 2.31e-13 * pow(T, -1.3) * exp(-23./T) * nH_
                         * y_H2 * y_Cplus;
  ksr_(idmap_(29)) = 0.99e-13 * pow(T, -1.3) * exp(-23./T) * nH_
                         * y_H2 * y_Cplus;
  //sanity check
  if (check_index) {
    std::stringstream msg; //error message
    msg << "### FATAL ERROR in ChemNetwork UpdateRatesSpecial [ChemNetwork]: "
      << "Specital rate reaction type does not match" << std::endl;
    for (int i=0; i<ns_cr; i++) {
      if (idtype_(indices_cr[i]) != ReactionType::cr) {
        ATHENA_ERROR(msg);
      }
    }
    for (int i=0; i<ns_gr; i++) {
      if (idtype_(indices_gr[i]) != ReactionType::grain) {
        ATHENA_ERROR(msg);
      }
    }
    for (int i=0; i<ns_2body; i++) {
      if (idtype_(indices_2body[i]) != ReactionType::twobody) {
        ATHENA_ERROR(msg);
      }
    }
    for (int i=0; i<ns_sr; i++) {
      if (idtype_(indices_sr[i]) != ReactionType::special) {
        ATHENA_ERROR(msg);
      }
    }
    check_index = false;
  }
  return;
}

Real CII_rec_rate(const Real temp) {
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
