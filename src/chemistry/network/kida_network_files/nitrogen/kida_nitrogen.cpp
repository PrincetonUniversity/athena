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
  const Real logT = log10(T);
	const Real logT4 = log10(T/1.0e4);
	const Real lnTe = log(T * 8.6173e-5);
  const Real sqrtT = sqrt(T);
  const Real kcr_H_fac = 2.3*y_H2 + 1.5 * y_H; //ratio of total to primary rate

	//cosmic ray reactions
  const int ns_cr = 2;
  const int indices_cr[ns_cr] = {1, 2};
  //(1) H + CR -> H+ + e-
  kcr_(id7map_(1)) = kcr_H_fac * rad_(index_cr_); 
  //(2) H2 + CR -> H2+ + e-
  kcr_(id7map_(2)) = 2. * kcr_H_fac * rad_(index_cr_); 

	// Grain assisted reactions
  const int ns_gr = 1;
  const int indices_gr[ns_gr] = {16};
	//(16) H + H + gr -> H2 + gr , from Draine book chapter 31.2 page 346, Jura 1975
	kgr_(id7map_(16)) = 3.0e-17 * nH_ * Z_d_ / (a_d_/1e-5);

  //2body reactions
  const int ns_2body = 14;
  const int indices_2body[ns_2body] = {3, 4, 5, 6, 7,
                                       8, 9, 10, 11, 12,
                                       13, 14, 15, 17};
  //(3) H+ + e- -> H  -- Case B
	k2body_(id7map_(3)) = 2.753e-14 * pow( 315614.0 / T, 1.5) 
									 * pow(  1.0 + pow( 115188.0 / T, 0.407) , -2.242 );
  //(4) O + H+ -> H + O+
  //(10) H + O+ -> O + H+
  k2body_(id7map_(4)) = ( 1.1e-11 * pow(T, 0.517) + 4.0e-10 * pow(T, 6.69e-3) 
                             ) * exp(-227./T);
  k2body_(id7map_(10)) = (4.99e-11* pow(T, 0.405) + 7.5e-10 * pow(T, -0.458) 
                             );

  //(5) H2 + H2+ -> H + H3+
  k2body_(id7map_(5)) = 1.76e-9 * pow(T, 0.042) * exp(- T/46600.); 

  //(6) C + H3+ -> H2 + CH+       --Vissapragada2016 new rates
  //(7) C + H3+ -> H + CH2+       --Vissapragada2016 new rates
  //valid between 1-1e4 K
  const Real A_kCHp = 6.93e-10;
  const Real n_kCHp = -8.34e-2;
  const Real c_kCHp[4] = {-7.24e-9, -9.07e-10, 7.48e-8, 9.93e-5};
  const Real Ti_kCHp[4] = {8.01, 1.92, 4.19e1, 8.08e3};
  const Real A_kCH2p =  3.35e-10;
  const Real n_kCH2p =  1.89e-1;
  const Real c_kCH2p[4] = {2.73e-8, 5.58e-9, 7.46e-8, -1.92e-4};
  const Real Ti_kCH2p[4] = {6.49, 1.30, 1.90e1, 1.62e4};
  const Real t1_CHp = A_kCHp * pow( 300./T, n_kCHp);
  const Real t2_CHp = c_kCHp[0] * exp(-Ti_kCHp[0]/T) + c_kCHp[1] * exp(-Ti_kCHp[1]/T)
           + c_kCHp[2]*exp(-Ti_kCHp[2]/T) + c_kCHp[3] *exp(-Ti_kCHp[3]/T);
  const Real t1_CH2p = A_kCH2p * pow( 300./T, n_kCH2p);
  const Real t2_CH2p = c_kCH2p[0] * exp(-Ti_kCH2p[0]/T)
           + c_kCH2p[1] * exp(-Ti_kCH2p[1]/T)
           + c_kCH2p[2]*exp(-Ti_kCH2p[2]/T) + c_kCH2p[3] *exp(-Ti_kCH2p[3]/T);
  k2body_(id7map_(6)) = (t1_CHp + pow(T, -1.5) * t2_CHp);
  k2body_(id7map_(7)) = (t1_CH2p + pow(T, -1.5) * t2_CH2p);
  //(8) He+ + e- -> He  -- Case B
	k2body_(id7map_(8)) = 1e-11*pow(T, -0.5)*(11.19 + 
                             (-1.676 + (-0.2852 + 0.04433*logT) * logT )* logT);
  //(9) C+ + e- -> C    -- Include RR and DR, Badnell2003, 2006.
  k2body_(id7map_(9)) = CII_rec_rate(T);

  //(11)O + H3+ -> H2 + OH+       --de Routte et al. (2016) new rates
  //(12)O + H3+ -> H + H2O+       --de Routte et al. (2016) new rates
  const Real a_OHp[3] = {5.1142e-10, 2.6568e-11, 9.8503e-15};
  const Real b_OHp[3] = {1.6747e-2, -9.9613e-5, 1.1006e-6};
  const Real a0_H2Op = 4.2253e-10;
  const Real b_H2Op[3] = {3.4977e-3, -1.4126e-4, 6.3584e-5};
  k2body_(id7map_(11)) = (a_OHp[0] + a_OHp[1]*sqrtT + a_OHp[2]*T) / (
        pow(T,1./6.) + b_OHp[0]*sqrtT + b_OHp[1]*T + b_OHp[2]*T*sqrtT);
  k2body_(id7map_(12)) = a0_H2Op / (
        pow(T,1./6.) + b_H2Op[0]*sqrtT + b_H2Op[1]*T + b_H2Op[2]*T*sqrtT);

  //Collisional dissociation, k>~1.0e-30 at T>~5e2.
  //(13) H + e- -> H+ + e- + e-
  //(14) H2 + H -> H + H + H
  //(15) H2 + H2 -> H + H + H2
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
    k2body_(id7map_(13)) = exp( -3.271396786e1 + 
                      (1.35365560e1 + (- 5.73932875 + (1.56315498 
                    + (- 2.877056e-1 + (3.48255977e-2 + (- 2.63197617e-3
                    + (1.11954395e-4 + (-2.03914985e-6)
                       *lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe) *lnTe);

    k2body_(id7map_(14)) = pow(10, log10(k9h) *  n2ncr/(1. + n2ncr) 
                         + log10(k9l) / (1. + n2ncr));
    k2body_(id7map_(15)) = pow(10, log10(k10h) *  n2ncr/(1. + n2ncr) 
                         + log10(k10l) / (1. + n2ncr));
  } else {
    k2body_(id7map_(13)) = 0.;
    k2body_(id7map_(14)) = 0.;
    k2body_(id7map_(15)) = 0.;
  }
  //(17) H2 + N+ -> H + NH+      -- depends on o/p, Dislaire et al. (2012)
  const Real x_oH2 = o2pH2_ / (1. + o2pH2_);
  const Real ko = 4.2e-10 * pow(T/300., -0.17) * exp(-44.5/T);
  const Real kp = 8.35e-10 * exp(-168.5/T);
  k2body_(id7map_(17)) = ( x_oH2*ko + (1.-x_oH2)*kp );

  //sanity check
  if (check_index) {
    std::stringstream msg; //error message
    msg << "### FATAL ERROR in ChemNetwork UpdateRatesSpecial [ChemNetwork]: "
      << "Specital rate reaction type does not match" << std::endl;
    for (int i=0; i<ns_cr; i++) {
      if (id7type_(indices_cr[i]) != ReactionType::cr) {
        ATHENA_ERROR(msg);
      }
    }
    for (int i=0; i<ns_gr; i++) {
      if (id7type_(indices_gr[i]) != ReactionType::grain_implicit) {
        ATHENA_ERROR(msg);
      }
    }
    for (int i=0; i<ns_2body; i++) {
      if (id7type_(indices_2body[i]) != ReactionType::twobody) {
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
