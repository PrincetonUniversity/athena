#ifndef THERMO_H_
#define THERMO_H_

//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file thermo.hpp
//  \brief definitions for heating and cooling processes
//======================================================================================

// Athena++ classes headers
#include "../../athena.hpp"

class Thermo {
  friend class ChemNetwork;
  friend class RadIntegrator;
	public:
		Thermo();
		//----heating and cooling rate per H----
		// Heating by cosmic ray ionization of H, He, and H2.
		// Arguments:
		// xs = ns/nH, abundances of species s.
		// ks: cosmic ray ioniztion rate per s particle.
		// nH: number density of H atom.
		// Return: 
		// cosmic ray heating in erg H^-1 s^-1
		static Real HeatingCr(const Real xe, const Real nH,
										 const Real xHI, const Real xHe, const Real xH2,
		   							 const Real kHI, const Real kHe, const Real kH2);
		// Heating by photo electric effect on dust, including collsional cooling. 
    // WD2001 Table 2 and 3. Use the diffuse ISM: Rv=3.1, ISRF, bc=4.
		// Arguments:
		// G: UV radiation scaled by solar neighbourhood value, including
    // extinction.  G=G0 * exp(-NH*sigmaPE_) at one line of sight.
		// Zd: dust abundance scaled by solar neighbourhood value.
    // T: temperature in K
    // ne: electron number density in cm^-3.
		// Return:
		// Photo electric by dust heating rate in erg H^-1 s^-1.
		static Real HeatingPE(const Real G, const Real Zd, const Real T,
                            const Real ne);
		// Heating by photo electric effect on dust, including collsional cooling. 
    // Wolfire et al. 2003 Equation (19)
		// Arguments:
		// G: UV radiation scaled by solar neighbourhood value, including
    // extinction.  G=G0 * exp(-NH*sigmaPE_) at one line of sight.
		// Z_PAH: dust (PAH) abundance scaled by solar neighbourhood value.
    // T: temperature in K
    // ne: electron number density in cm^-3.
		// Return:
		// Photo electric by dust heating rate in erg H^-1 s^-1.
		static Real HeatingPE_W03(const Real G, const Real Z_PAH, const Real T,
                              const Real ne, const Real phi_PAH);
    // Heating by H2 formation on dust grains.
    // From Hollenbach + McKee 1979
    // (0) *H + *H + gr -> H2 + gr
    // Arguments:
    // xi = ni/nH
    // T: temperature in K
    // kgr: grain reaction rates, kgr_ in NL99p. 
    // Return:
    // Heating rate by H2 formation on dust grains in erg H^-1 s^-1.
    static Real HeatingH2gr(const Real xHI, const Real xH2, const Real nH,
                       const Real T, const Real kgr);
    // Heating by H2 UV pumping.
    // From Hollenbach + McKee 1979
    // Arguments:
    // xi = ni/nH
    // T: temperature in K
    // dot_xH2_photo = dxH2/dt by photo dissociation of H2 by UV light.
    // Calculated in RHS in NL99p.
    // Return:
    // Heating rate by H2 UV pumping in erg H^-1 s^-1.
    static Real HeatingH2pump(const Real xHI, const Real xH2, const Real nH,
                         const Real T, const Real dot_xH2_photo);
    // Heating by H2 photo dissiociation.
    // From Black + Dalgarno 1977, 0.4eV per reaction.
    // Arguments:
    // dot_xH2_photo = dxH2/dt by photo dissociation of H2 by UV light.
    // Calculated in RHS in NL99p.
    // Return:
    // Heating rate by H2 photo dissiociation in erg H^-1 s^-1.
    static Real HeatingH2diss(const Real dot_xH2_photo);
		// Cooling by C+ fine structure line.
		// Collisional species: HI, H2, e.
		// Arguments:
		// xCII = nC+/nH.
		// ni: number density of species i, in cm^-3.
		// T: temperature in K
		// Return:
		// Cooling rate for C+ fine structure line in erg H^-1 s^-1
		static Real CoolingCII(const Real xCII, const Real nHI, const Real nH2,
										  const Real ne, const Real T);
		// Cooling by CI fine structure line.
		// Collisional species: HI, H2
    // Note: ignored H+.
    // Arguments:
		// xCI = nCI/nH.
		// ni: number density of species i, in cm^-3.
		// T: temperature in K
		// Return:
		// Cooling rate for C fine structure line in erg H^-1 s^-1
		static Real CoolingCI(const Real xCI, const Real nHI, const Real nH2,
                     const Real ne, const Real T);
		// Cooling by OI fine structure line.
		// Collisional species: HI, H2, e
    // Note: the cooling rate is very insensitive to ne, for xe <~0.5.
    // At xe >~0.5 region, the O+ and other cooling will start to be important
    // anyway.
		// Arguments:
		// xOI = nOI/nH.
		// ni: number density of species i, in cm^-3.
		// T: temperature in K
		// Return:
		// Cooling rate for OI fine structure line in erg H^-1 s^-1
		static Real CoolingOI(const Real xOI, const Real nHI, const Real nH2,
										 const Real ne, const Real T);
		// Cooling by collisional exicited lyman alphya line.
		// Collisional species: e
		// Arguments:
		// xHI = nHI/nH.
		// ni: number density of species i, in cm^-3.
		// T: temperature in K
		// Return:
		// Cooling rate for Lyman alpha line in erg H^-1 s^-1
    static Real CoolingLya(const Real xHI, const Real ne, const Real T);
    // Cooling by CO rotational lines.
    // Collision species: HI, H2, e
    // Note: from Omukai+2010
    // Auguments:
    // xCO = nCO/nH
    // ni: number density of species i, in cm^-3
    // T: temperature in K
    // NCOeff: effective column density of CO using LVG approximation. See
    // notes.
    // Return:
    // Cooling rate for CO rotational lines in erg H^-1 s^-1
    static Real CoolingCOR(const Real xCO, const Real nHI, const Real nH2, 
                      const Real ne, const Real T, const Real NCOeff);
    // Cooling by H2 vibration and rotation lines.
    // Collision species: HI, H2, He, H+, e
    // Note: Using Glover + Abel 2008 fitting formulas in Table 8, assuming
    // that ortho to para ration of H2 is 3:1.
    // Auguments:
    // xH2 = nH2 / nH
    // ni: number density of species i, in cm^-3
    // T: temperature in K
    // Return:
    // Cooling rate for H2 vibrational and rotational lines in erg H^-1 s^-1///
    static Real CoolingH2(const Real xH2, const Real nHI, const Real nH2,
                     const Real nHe, const Real nHplus, const Real ne,
                     const Real T);
    // Cooling by dust thermo emission.
    // Tabulated from Despotic.
    // dedt_dust = L_CMB - G_dust, thermo - PsiGD = 0
    // dedt_dust = L_CMB + L_dust, ISRF - PsiGD = 0
    // maxmium rate of above.
    // Ignored ISRF and IR heating. Assume Zd = 1.
    // Auguments:
    // Zd: dust metalicity compared to solar neighbourhood.
    // nH: hydrogen number density. Here implicitly assume all in H2, which
    // determines alpha_gd in despotic (see eq B8 in despotic paper).
    // Tg: gas temperature
    // GISRF: strength of ISRF = chi * exp(-sigma_{d, ISRF} * NH)
    // Return:
    // Cooling rate for dust in erg H^-1 s^-1 
    static Real CoolingDust(const Real Zd, const Real nH, const Real Tg,
                            const Real GISRF);
    // Cooling by dust thermo emission, assume a constant dust temperature
    // Auguments:
    // Zd: dust metalicity compared to solar neighbourhood.
    // nH: hydrogen number density. Here implicitly assume all in H2, which
    // determines alpha_gd in despotic (see eq B8 in despotic paper).
    // Tg: gas temperature
    // Td: dust temperature
    // Cooling rate for dust in erg H^-1 s^-1 
    static Real CoolingDustTd(const Real Zd, const Real nH, const Real Tg,
                              const Real Td);
    // Cooling by reconbination of e on PAHs.
    // From WD2001 Eq(45).
    // Arguments:
    // Zd: dust metalicity compared to solar neighbourhood.
    // T: temperature in K.
    // ne: number denstiy of electrons.
		// G: UV radiation scaled by solar neighbourhood value, including
    // extinction. G=G0 * exp(-NH*sigmaPE_) at one line of sight.
    // Return:
    // Cooling rate for  recombination of e on PAHs in erg H^-1 s^-1 
    static Real CoolingRec(const Real Zd, const Real T, const Real ne, 
                             const Real G);
    // Cooling by reconbination of e on PAHs.
    // Wolfire et al. (2003) eq. (21)
    // Arguments:
    // Z_PAH: dust (PAH) metalicity compared to solar neighbourhood.
    // T: temperature in K.
    // ne: number denstiy of electrons.
		// G: UV radiation scaled by solar neighbourhood value, including
    // extinction. G=G0 * exp(-NH*sigmaPE_) at one line of sight.
    // Return:
    // Cooling rate for  recombination of e on PAHs in erg H^-1 s^-1 
    static Real CoolingRec_W03(const Real Z_PAH, const Real T, const Real ne, 
                               const Real G, const Real phi_PAH);
    // Cooling by collisional dissociation of H2
    //  H2 + *H -> 3 *H 
    //  H2 + H2 -> H2 + 2 *H
    // reaction heat: 4.48 eV from Krome Paper
    // Arguments:
    // xi = ni/nH
    // k_H2_H, k_H2_H2: reaction rate cooefficients (k2body_ in NL99p)
    // Return:
    // Cooling rate for H2 collisional dissociation in erg H^-1 s^-1.
    static Real CoolingH2diss(const Real xHI, const Real xH2,
                         const Real k_H2_H, const Real k_H2_H2);
    // Cooling by collisional ionization of HI
    //  *H + *e -> H+ + 2 *e
    // reaction heat: 13.6 eV from Krome Paper
    // Arguments:
    // xi = ni/nH
    // k_H_e: reaction rate cooefficients (k2body_ in NL99p)
    // Return:
    // Cooling rate for collisional ionization of HI in erg H^-1 s^-1.
    static Real CoolingHIion(const Real xHI, const Real xe, const Real k_H_e);
		// Cooling by Radiative recombination in hot gas, use CIE
		// T = 10^3.8-10^8. For T>10^8, use free-free emmission.
		// Arguments:
		// nH: hydrogen number density in cm^-3
		// T: temperature in Kelvin
		// Zg: gas metalicity relative to solar
		// Return: Cooling rate for hot gas in erg H^-1 s^-1
		static Real CoolingHotGas(const Real nH, const Real T, const Real Zg);
    // specific heat, assume that H2 rotational and vibrational levels not
    // excited.
    // xH2, xe = nH2 or ne / nH
    // xHe_total = xHeI + xHeII = 0.1 for solar value.
    // Return: specific heat per H atom.
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
		//-----C+ atomic data------
		// Return Collisional rate for C+ atom in s^-1.
		// Collisional species: HI, H2, e.
		// T: gas temeperature.
		// ni: number density of species i, in cm^-3.
		static Real q10CII_(const Real nHI, const Real nH2, const Real ne,
									 const Real T);
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
		//line cooling rate per H for 2 level atom. Ignore radiative excitation
		// and de-excitation, and assume optically thin.
		// Arguments:
		// q01 = \sum (nc * k_{s, 01}). Collisional excitation rate per second 
		// for all the collider species sumed up together from level 0 to 1.
		// q10 = \sum (nc * k_{s, 10}). Collisional de-excitation rate per second,
		// similar to q01.
		// A10: Enstein A coefficent for spontanious emission from level 1 to 0,
		// in sec^-1.
		// E10: energy difference of levels E1 - E0, in erg.
		// xs = ns/nH, abundances of species s.
		// Return:
		// Line cooling rate in erg H^-1 s^-1.
		static Real Cooling2Level_(const Real q01, const Real q10,
											    const Real A10, const Real E10,
													const Real xs);
		//line cooling rate per H for 3 level atom. Ignore radiative excitation
		// and de-excitation, and assume optically thin.
		// Arguments:
		// qij = \sum (nc * k_{s, ij}). Collisional excitation rate per second 
		// for all the collider species sumed up together from level i to j.
		// Aij: Enstein A coefficent for spontanious emission from level i to j,
		// in sec^-1, i > j.
		// Eij: energy difference of levels Ei - Ej, in erg.
		// xs = ns/nH, abundances of species s.
		// Return:
		// Total line cooling rate in erg H^-1 s^-1.
		static Real Cooling3Level_(const Real q01, const Real q10,
													const Real q02, const Real q20,
													const Real q12, const Real q21,
													const Real A10, const Real A20,
													const Real A21, const Real E10,
													const Real E20, const Real E21,
													const Real xs);
};
#endif //THERMO_H_
