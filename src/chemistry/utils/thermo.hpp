#ifndef CHEMISTRY_UTILS_THERMO_HPP_
#define CHEMISTRY_UTILS_THERMO_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file thermo.hpp
//! \brief definitions for heating and cooling processes

// Athena++ classes headers
#include "../../athena.hpp"

//! \class Thermo
//! \brief Heating and cooling functions
class Thermo {
  friend class ChemNetwork;
  friend class ChemRadIntegrator;
 public:
    Thermo();
    static Real HeatingCr(const Real xe, const Real nH,
                          const Real xHI, const Real xH2, const Real crir_prim);
    static Real HeatingPE(const Real G, const Real Zd, const Real T,
                            const Real ne);
    static Real HeatingPE_W03(const Real G, const Real Z_PAH, const Real T,
                              const Real ne, const Real phi_PAH);
    static Real HeatingH2gr(const Real xHI, const Real xH2, const Real nH,
                       const Real T, const Real kgr, const Real k_xH2_photo);
    static Real HeatingH2pump(const Real xHI, const Real xH2, const Real nH,
                         const Real T, const Real k_xH2_photo);
    static Real HeatingH2diss(const Real k_xH2_photo, const Real xH2);
    static Real CoolingCII(const Real xCII, const Real nHI, const Real nH2,
                      const Real ne, const Real T);
    static Real CoolingCI(const Real xCI, const Real nHI, const Real nH2,
                     const Real ne, const Real T);
    static Real CoolingOI(const Real xOI, const Real nHI, const Real nH2,
                     const Real ne, const Real T);
    static Real CoolingLya(const Real xHI, const Real ne, const Real T);
    static Real CoolingCOR(const Real xCO, const Real nHI, const Real nH2,
                      const Real ne, const Real T, const Real NCOeff);
    static Real CoolingH2(const Real xH2, const Real nHI, const Real nH2,
                     const Real nHe, const Real nHplus, const Real ne,
                     const Real T);
    static Real CoolingDust(const Real Zd, const Real nH, const Real Tg,
                            const Real GISRF);
    static Real CoolingDustTd(const Real Zd, const Real nH, const Real Tg,
                              const Real Td);
    static Real CoolingRec(const Real Zd, const Real T, const Real ne,
                             const Real G);
    static Real CoolingRec_W03(const Real Z_PAH, const Real T, const Real ne,
                               const Real G, const Real phi_PAH);
    static Real CoolingH2diss(const Real xHI, const Real xH2,
                         const Real k_H2_H, const Real k_H2_H2);
    static Real CoolingHIion(const Real xHI, const Real xe, const Real k_H_e);
    static Real CoolingHotGas(const Real nH, const Real T, const Real Zg);
    static Real CvCold(const Real xH2, const Real xHe_total, const Real xe);

 private:
    static const Real eV_; //eV in erg
    static const Real kb_; //boltzmann constant in erg/K
    static const Real ca_; //speed of light * radiation constant, or
                            // stephan-bolzmann constant*4
    static const Real TCMB_; //CMB temperature
    static const Real o2p_;//ratio of ortho to para H2
    static const Real fo_; //ortho H2 fraction
    static const Real fp_; //para H2 fraction
    static const Real sigmaPE_; //dust cross-section for 8-13.6eV photons in cm2
    static const Real sigmaISRF_; //dust cross-section for ISRF in cm2
    static const Real sigmad10_; //dust cross-section for IR at T=10K
    static const Real alpha_GD_;
    static Real q10CII_(const Real nHI, const Real nH2, const Real ne,
                   const Real T);
    //-----C+ atomic data------
    static const Real A10CII_;
    static const Real E10CII_;
    static const Real g0CII_;
    static const Real g1CII_;
    //-----HI atomic data------
    static const Real A10HI_;
    static const Real E10HI_;
    static const Real g0HI_;
    static const Real g1HI_;
    //-----CI atomic data------
    static const Real g0CI_;
    static const Real g1CI_;
    static const Real g2CI_;
    static const Real A10CI_;
    static const Real A20CI_;
    static const Real A21CI_;
    static const Real E10CI_;
    static const Real E20CI_;
    static const Real E21CI_;
    //-----OI atomic data------
    static const Real g0OI_;
    static const Real g1OI_;
    static const Real g2OI_;
    static const Real A10OI_;
    static const Real A20OI_;
    static const Real A21OI_;
    static const Real E10OI_;
    static const Real E20OI_;
    static const Real E21OI_;
    //-----CO cooling table data, from Omukai+2010-----
    static const int lenTCO_ = 11;
    static const int lenNeffCO_ = 11;
    static const Real TCO_[lenTCO_];
    static const Real NeffCO_[lenNeffCO_];
    static const Real L0CO_[lenTCO_];
    static const Real LLTECO_[lenNeffCO_*lenTCO_];
    static const Real nhalfCO_[lenNeffCO_*lenTCO_];
    static const Real alphaCO_[lenNeffCO_*lenTCO_];
    //-----radiative cooling from Schure 2009 -------
    static const int len_rad_cool_ = 110;
    static const Real log_Trad_[len_rad_cool_];
    static const Real log_gamma_H_He_[len_rad_cool_];
    static const Real log_gamma_Z_[len_rad_cool_];
    //-----PE heating coefficients from WD2001 Table 2, second last line ---
    static const Real CPE_[7];
    //-----Collisional cooling included in PE heating, WD2001 Table3, second
    // last line---
    static const Real DPE_[5];
    //------dust cooling table from Despotic--------
    // dedt_dust = L_CMB - L_dust, thermo - PsiGD = 0
    // Ignored ISRF heating.
    // Tabulate PsiGD as a function of Tg and nH.
    static const int lenTg_ = 10;
    static const int lennH_ = 15;
    static const Real logTg_[lenTg_];
    static const Real lognH_[lennH_];
    static const Real logps_[lennH_ * lenTg_];
    static Real Cooling2Level_(const Real q01, const Real q10,
                          const Real A10, const Real E10,
                          const Real xs);
    static Real Cooling3Level_(const Real q01, const Real q10,
                          const Real q02, const Real q20,
                          const Real q12, const Real q21,
                          const Real A10, const Real A20,
                          const Real A21, const Real E10,
                          const Real E20, const Real E21,
                          const Real xs);
};
#endif //CHEMISTRY_UTILS_THERMO_HPP_
