#ifndef SHIELDING_HPP
#define SHIELDING_HPP

//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file shielding.hpp
//  \brief definition of shielding class.
//======================================================================================
#include "../../athena.hpp"
#include <math.h> //a^x = pow(a,x)

//! \class Shielding
//  \brief Shielding functions.
class Shielding {
  public:
		//van Dishoeck & Black's shielding function.
		// NCO, NH2: column densith of CO and H2 in cm^-2.
    Shielding();
		static Real fShield_CO_vDB(const Real NCO, const Real NH2);
    //CO sheilding from Visser+2009, Table 5
		static Real fShield_CO_V09(const Real NCO, const Real NH2);
    //H2 self shielding from Draine+Bertoldi1996
		static Real fShield_H2(const Real NH2, const Real bH2);
    //CI self shielding. 
		static Real fShield_C(const Real NC, const Real NH2);
    //CI shielding of CO
    static Real fShield_CO_C(const Real NC);
  private:
		//CO column density for DB table
		static const int len_NCO_DB_ = 8;
		static const Real logNCOvDB_[len_NCO_DB_];
		// H2 column densities for DB table
		static const int len_NH2_DB_ = 6;
		static const Real logNH2vDB_[len_NH2_DB_];
		// Tabulated shielding factors 
		static const Real ThetavDB_[len_NH2_DB_][len_NCO_DB_];

    //Visser+ 2009 Table 5, b(CO)=0.3, Tex(CO)=5K
    static const int len_NCO_V09_ = 47;
    static const Real logNCOV09_[len_NCO_V09_];
    static const int len_NH2_V09_ = 42;
    static const Real logNH2V09_[len_NH2_V09_];
		static const Real ThetaV09_[len_NH2_V09_][len_NCO_V09_];

		// Find index of linear interpretation, return the first index i for i,
		// i+1.
		static int LinearInterpIndex(const int len, const Real xarr[], 
                                 const Real x){
      int i = 0;
      if ( x < xarr[0]) {
        return 0;
      } else if ( x > xarr[len-1]) {
        return len-2;
      } else {
        for (i=0; x>xarr[i]; i++) {}
        return i-1;
      }
    }

		static Real LinearInterp(const Real x0, const Real x1,
                               const Real y0, const Real y1,
                               const Real x){
      return y0 + ( (y1-y0)/(x1-x0) ) * (x-x0);
    }
};

#endif //SHIELDING_HPP
