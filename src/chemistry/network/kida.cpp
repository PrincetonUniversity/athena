//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file kida.cpp
//! \brief implementation of functions in class ChemNetwork, using the simple
//! network for kida style network files.

// this class header
#include "kida.hpp"

// C header

// C++ header
#include <algorithm>    // std::find()
#include <cmath>       //M_PI
#include <cstdio>
#include <fstream>   //file()
#include <iostream>   // endl
#include <iterator>     // std::distance()
#include <limits>    //inf
#include <sstream>    // stringstream

// Athena++ header
#include "../../chem_rad/chem_rad.hpp"
#include "../../defs.hpp"
#include "../../eos/eos.hpp"
#include "../../globals.hpp"
#include "../../hydro/hydro.hpp"
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"
#include "../../scalars/scalars.hpp"
#include "../../units/units.hpp"
#include "../../utils/string_utils.hpp"
#include "../utils/chemistry_utils.hpp"
#include "../utils/thermo.hpp"
#include "network.hpp"

static bool output_rates = true;
static bool output_thermo = true;  // only used if DEBUG

//----------------------------------------------------------------------------------------
//! \brief ChemNetwork constructor
ChemNetwork::ChemNetwork(MeshBlock *pmb, ParameterInput *pin) :
    pmy_spec_(pmb->pscalars), pmy_mb_(pmb), flag_T_rates_(true), id7max_(0), n_cr_(0),
    icr_H_(-1), icr_H2_(-1), icr_He_(-1), n_crp_(0), n_ph_(0), iph_H2_(-1),
    n_2body_(0), i2body_H2_H_(-1), i2body_H2_H2_(-1), i2body_H_e_(-1),
    n_2bodytr_(0), n_gr_(0), igr_H_(-1), n_sr_(0), n_gc_(0),
    n_freq_(0), index_gpe_(0), index_cr_(0), gradv_(0.) {
  // set the parameters from input file
  Z_g_ = pin->GetOrAddReal("chemistry", "Z_g", 1.); // gas metallicity
  // PAH recombination efficiency
  phi_PAH_ = pin->GetOrAddReal("chemistry", "phi_PAH", 0.4);
  // size of the dust grain in cm, default 0.1 micron
  a_d_ = pin->GetOrAddReal("chemistry", "a_d", 1e-5);
  // density of grain in cgs, default 2 g/cm3
  rho_d_ = pin->GetOrAddReal("chemistry", "rho_d", 2.);
  // mass of the dust grain in g, assuming density of 2 g/cm3
  m_d_ = rho_d_ * 4.*M_PI * a_d_*a_d_*a_d_ / 3.;
  // fixing grain abundance for implicity grain assisted reactions
  is_fixed_Zd_ = pin->GetOrAddBoolean("chemistry", "is_fixed_Zd", false);
  // fixing PAH abundance for PE heating and Rec cooling
  is_fixed_PAH_ = pin->GetOrAddBoolean("chemistry", "is_fixed_PAH", false);
  if (is_fixed_PAH_) {
    Z_PAH_ = pin->GetReal("chemistry", "Z_PAH");
  } else {
    Z_PAH_ = 0.;
  }
  // grain abundance, read from initialization
  const Real sinit = pin->GetOrAddReal("problem", "s_init", 0.);
  const Real xg0 = pin->GetOrAddReal("problem", "s_init_g0", sinit);
  const Real xgp = pin->GetOrAddReal("problem", "s_init_g+", sinit);
  const Real xgm = pin->GetOrAddReal("problem", "s_init_g-", sinit);
  x_d_ = xg0 + xgp + xgm; // total dust abundance
  // setting dust metallicity for implicit grain assisted reactions
  if (is_fixed_Zd_) { // fixing Z_d
    Z_d_ = pin->GetReal("chemistry", "Z_d");
  } else { // setting Z_d according to the grain abundance
    // relative abundance of all dust at Zd=1
    const Real xdZ1 = 0.013 * 1.4 * 1.67e-24 / m_d_;
    Z_d_ = x_d_ / xdZ1;
  }
  if (Globals::my_rank == 0) {
    std::cout << "Z_d=" << Z_d_ << ", Z_PAH=" << Z_PAH_ << std::endl;
  }
  o2pH2_ = pin->GetOrAddReal("chemistry", "o2pH2", 3.); // ortho to para H2 ratio
  Yi_ = pin->GetOrAddReal("chemistry", "Yi", 1e-3); // ortho to para H2 ratio

  // temperature
  if (NON_BAROTROPIC_EOS) {
    // check adiabatic index, this is for calling CvCold later
    const Real eps = 0.01;
    const Real gm = pin->GetReal("hydro", "gamma");
    const Real diff = std::abs(gm - 5./3.) / (5./3.);
    if (diff > eps) {
      std::stringstream msg1; // error message
      msg1 << "### FATAL ERROR in ChemNetwork constructor" << std::endl
           << "kida network with energy equation: adiabatic index must be 5/3."
           << std::endl;
      ATHENA_ERROR(msg1);
    }
    temperature_ = 0.;
  } else {
    // isothermal, temperature is calculated from a fixed mean molecular weight
    const Real mu_iso = pin->GetReal("chemistry", "mu_iso");
    const Real cs = pin->GetReal("hydro", "iso_sound_speed");
    temperature_ = cs * cs * pmy_mb_->pmy_mesh->punit->code_temperature_mu_cgs * mu_iso;
    std::cout << "isothermal temperature = " << temperature_ << " K" << std::endl;
  }
  // whether to cap temperature if the reaction is outside of the temperature range
  // only for 2 body reactions
  is_Tcap_2body_ = pin->GetOrAddBoolean("chemistry", "is_Tcap_2body", false);
  // minimum temperature for reaction rates, also applied to energy equation
  temp_min_rates_ = pin->GetOrAddReal("chemistry", "temp_min_rates", 1.);
  // minimum temperature below which cooling is turned off
  temp_min_cool_ = pin->GetOrAddReal("chemistry", "temp_min_cool", 1.);
  // cooling for neutral medium is capped at this temperature
  temp_max_cool_nm_ = pin->GetOrAddReal("chemistry", "temp_max_cool_nm", 1.0e9);
  // dust temperature for dust thermo cooling of the gas at high densities
  temp_dust_thermo_ = pin->GetOrAddReal("chemistry", "temp_dust_thermo", 10.);
  // folder of the network
  network_dir_ = pin->GetString("chemistry", "network_dir");
  // CO cooling parameters
  // Maximum CO cooling length in cm. default 100pc.
  Leff_CO_max_ = pin->GetOrAddReal("chemistry", "Leff_CO_max", 3.0e20);

  // read in the species
  std::string species_file_name = network_dir_ + "/species.dat";
  std::ifstream species_file(species_file_name);
  if (!species_file) {
    std::stringstream msg; // error message
    msg << "### FATAL ERROR in ChemNetwork constructor [ChemNetwork]" << std::endl
        << "Cannot open file" << species_file_name << std::endl;
    ATHENA_ERROR(msg);
  }
  std::string line;
  int nline = 0;
  int nices = 0;
  while (getline(species_file, line)) {
    // trim white spaces
    StringUtils::trim(line);
    // skip blank lines and comments
    if(line.empty() || (line.find("!") == 0)) {
      continue;
    }
    KidaSpecies si(line, nline);
    // set mass for dust grains
    if (si.name.find("g") == 0) {
      si.SetMass(m_d_);
    }
    species_names[nline] = si.name;
    ispec_map_[si.name] = nline;
    // cound ice species
    if (si.name.find("s") == 0) {
      nices++;
    }
    species_.push_back(si);
    nline++;
  }
  if (nline != NSPECIES) {
    std::stringstream msg; // error message
    msg << "### FATAL ERROR in ChemNetwork constructor [ChemNetwork]" << std::endl
        << "number of species in species.dat does not match the number of scalars ("
        << NSPECIES << ")" << std::endl;
    ATHENA_ERROR(msg);
  }
  // set ice species index
  nices_ = nices;
  id_ices_.NewAthenaArray(nices);
  nices = 0;
  for (int i=0; i<NSPECIES; i++) {
    if (species_[i].name.find("s") == 0) {
      id_ices_(nices) = i;
      nices++;
    }
  }

  // read in the reactions
  std::string reactions_file_name = network_dir_ + "/reactions.dat";
  std::ifstream reactions_file(reactions_file_name);
  std::vector<int> rids; // array of id for reactions
  if (!reactions_file) {
    std::stringstream msg; // error message
    msg << "### FATAL ERROR in ChemNetwork constructor [ChemNetwork]" << std::endl
        << "Cannot open file" << reactions_file_name << std::endl;
    ATHENA_ERROR(msg);
  }
  nline = 0;
  while (getline(reactions_file, line)) {
    // trim white spaces
    StringUtils::trim(line);
    // skip blank lines and comments
    if(line.empty() || (line.find("!") == 0)) {
      continue;
    }
    KidaReaction ri(line);
    if (ri.id_ <= 0) {
      std::stringstream msg; // error message
      msg << "### FATAL ERROR in ChemNetwork constructor [ChemNetwork]" << std::endl
          << "reaction ID ( " << ri.id_ << ") not positive" << std::endl;
      ATHENA_ERROR(msg);
    }
    auto itr = std::find(rids.rbegin(), rids.rend(), ri.id_);
    const int ifind = rids.rend() - itr - 1;
    if (ifind < 0) {
      rids.push_back(ri.id_);
      if (ri.id_ > id7max_ && ri.formula_ == 7) {
        id7max_ = ri.id_;
      }
      reactions_.push_back(ri);
    } else {
      if (reactions_[ifind].reactants_ == ri.reactants_
          && reactions_[ifind].products_ == ri.products_
          && reactions_[ifind].itype_ == ri.itype_
          && reactions_[ifind].Tmax_ < ri.Tmin_
          && reactions_[ifind].formula_ != 7 && ri.formula_ != 7
          && ri.itype_ >= 4 && ri.itype_<= 8 ) {
        if (ifind != nline - 1) {
          std::stringstream msg; // error message
          msg << "### FATAL ERROR in ChemNetwork constructor [ChemNetwork]"
              << std::endl << "reactions ID ( " << ri.id_ << " ) of different"
              << " temperature ranges are not arranged next to each other."
              << std::endl;
          ATHENA_ERROR(msg);
        }
        rids.push_back(ri.id_);
        reactions_.push_back(ri);
        const int nprev = std::count(id_2bodytr_.begin(), id_2bodytr_.end(), ri.id_);
        if (nprev == 0) {
          id_2bodytr_.push_back(ri.id_);
          n_2bodytr_++;
        } else if (nprev < n_range_ - 1) {
          id_2bodytr_.push_back(ri.id_);
        } else {
          std::stringstream msg; // error message
          msg << "### FATAL ERROR in ChemNetwork constructor [ChemNetwork]"
              << std::endl << "reaction ID ( " << ri.id_ << " )"
              << "too many temperature ranges" << std::endl;
          ATHENA_ERROR(msg);
        }
      } else {
        std::stringstream msg; // error message
        msg << "### FATAL ERROR in ChemNetwork constructor [ChemNetwork]"
            << std::endl << "reaction ID ( " << ri.id_ << ") not unique" << std::endl;
        ATHENA_ERROR(msg);
      }
    }
    nline++;
  }
  nr_ = nline; // number of reactions
  if (id7max_ > 0) {
    id7map_.NewAthenaArray(id7max_+1);
    id7type_.NewAthenaArray(id7max_+1);
    for (int i=0; i<id7max_+1; i++) {
      id7map_(i) = -1;
      id7type_(i) = ReactionType::none;
    }
  }

  // initialize coefficients of reactions
  InitializeReactions();

  // radiation related variables
  const int nfreq = pin->GetOrAddInteger("chem_radiation", "n_frequency", 1);
  n_freq_ = n_ph_ + 2;
  std::stringstream msg;
  // check whether number of frequencies equal to the input file specification
  if (nfreq != n_freq_) {
    msg << "### FATAL ERROR in ChemNetwork constructor" << std::endl
        << "number of frequencies in radiation: " << nfreq
        << " not equal to that in chemistry: " << n_freq_  << std::endl;
    ATHENA_ERROR(msg);
  }
  index_gpe_ = n_ph_;
  index_cr_ = n_ph_ + 1;
  rad_.NewAthenaArray(n_freq_);

  if (DEBUG) {
    PrintProperties();
  }
}

//----------------------------------------------------------------------------------------
//! \brief ChemNetwork destructor
ChemNetwork::~ChemNetwork() {}

//----------------------------------------------------------------------------------------
//! \fn void ChemNetwork::InitializeNextStep(const int k, const int j, const int i)
//! \brief Set the rates of chemical reactions, eg. through density and radiation field.
//!
//! k, j, i are the corresponding index of the grid
void ChemNetwork::InitializeNextStep(const int k, const int j, const int i) {
  Real rho, rho_floor, rad_sum;
  const int nang = pmy_mb_->pchemrad->nang;
  // density
  rho = pmy_mb_->phydro->w(IDN, k, j, i);
  // apply density floor
  rho_floor = pmy_mb_->peos->GetDensityFloor();
  rho = (rho > rho_floor) ?  rho : rho_floor;
  // hydrogen atom number density
  nH_ =  rho;
  // average radiation field of all angles
  for (int ifreq=0; ifreq < n_freq_; ++ifreq) {
    rad_sum = 0;
    // radiation
    for (int iang=0; iang < nang; ++iang) {
      rad_sum += pmy_mb_->pchemrad->ir(k, j, i, ifreq * nang + iang);
    }
    if (ifreq == index_cr_) {
      rad_(index_cr_) = rad_sum / static_cast<float>(nang);
    } else {
      rad_(ifreq) = rad_sum / static_cast<float>(nang);
    }
    if (DEBUG) {
      if (std::isnan(rad_(ifreq))) {
        printf("InitializeNextStep: ");
        printf("ifreq=%d, nang=%d, rad_sum=%.2e\n", ifreq, nang, rad_sum);
        OutputRates(stdout);
      }
    }
  }
  // CO cooling paramters
  SetGrad_v(k, j, i);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ChemNetwork::RHS(const Real t, const Real *y, const Real ED,
//!                       Real *ydot)
//! \brief RHS: right-hand-side of ODE.
//!
//! dy/dt = ydot(t, y). Here y are the abundance
//! of species. details see CVODE package documentation.
//! all input/output variables are in code units

void ChemNetwork::RHS(const Real t, const Real *y, const Real ED,
                      Real *ydot) {
  Real rate = 0;
  // energy per hydrogen atom
  Real E_ergs = ED * pmy_mb_->pmy_mesh->punit->code_energydensity_cgs / nH_;
  // store previous y including negative abundance correction
  Real *y0 = new Real[NSPECIES]; // correct negative abundance, only for UpdateRates()
  Real *ydotg = new Real[NSPECIES];

  for(int i=0; i<NSPECIES; i++) {
    ydotg[i] = 0.0;
  }

  // correct negative abundance to zero, used in rate update
  for (int i=0; i<NSPECIES; i++) {
    if (y[i] < 0) {
      y0[i] = 0;
    } else {
      y0[i] = y[i];
    }
    // throw error if nan, or inf, or large negative value occurs
    if ( std::isnan(y[i]) || std::isinf(y[i]) ) {
      printf("RHS: ");
      for (int j=0; j<NSPECIES; j++) {
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

  UpdateRates(y0, E_ergs);

  if (output_rates) {
    FILE *pf = fopen("chem_network_init.dat", "w");
    OutputRates(pf);
    fclose(pf);
    output_rates = false;
  }

  // cosmic ray reactions
  for (int i=0; i<n_cr_; i++) {
    rate = kcr_(i) * y[incr_(i)];
    ydotg[incr_(i)] -= rate;
    ydotg[outcr1_(i)] += rate;
    ydotg[outcr2_(i)] += rate;
    if (outcr3_(i) >= 0) {
      ydotg[outcr3_(i)] += rate;
    }
  }

  // cosmic ray induced photo reactions
  for (int i=0; i<n_crp_; i++) {
    rate = kcrp_(i) * y[incrp_(i)];
    ydotg[incrp_(i)] -= rate;
    ydotg[outcrp1_(i)] += rate;
    ydotg[outcrp2_(i)] += rate;
  }

  // FUV photo-dissociation and photo-ionisation
  for (int i=0; i<n_ph_; i++) {
    rate = kph_(i) * y[inph_(i)];
    ydotg[inph_(i)] -= rate;
    ydotg[outph1_(i)] += rate;
    ydotg[outph2_(i)] += rate;
  }

  // 2body reactions
  for (int i=0; i<n_2body_; i++) {
    rate =  k2body_(i) * y[in2body_(i, 0)] * y[in2body_(i, 1)] * nH_;
    if (y[in2body_(i, 0)] < 0 && y[in2body_(i, 1)] < 0) {
      rate *= -1.;
    }
    for (int jin=0; jin<n_in2body_; jin++) {
      if (in2body_(i, jin) >= 0) {
        ydotg[in2body_(i, jin)] -= rate;
      }
    }
    for (int jout=0; jout<n_out2body_; jout++) {
      if (out2body_(i, jout) >= 0) {
        ydotg[out2body_(i, jout)] += rate;
      }
    }
  }

  // 2bodytr reactions
  for (int i=0; i<n_2bodytr_; i++) {
    rate =  k2bodytr_(i) * y[in2bodytr1_(i)] * y[in2bodytr2_(i)] * nH_;
    if (y[in2bodytr1_(i)] < 0 && y[in2bodytr2_(i)] < 0) {
      rate *= -1.;
    }
    ydotg[in2bodytr1_(i)] -= rate;
    ydotg[in2bodytr2_(i)] -= rate;
    ydotg[out2bodytr1_(i)] += rate;
    if (out2bodytr2_(i) >= 0) {
      ydotg[out2bodytr2_(i)] += rate;
    }
    if (out2bodytr3_(i) >= 0) {
      ydotg[out2bodytr3_(i)] += rate;
    }
  }

  // grain assisted reactions
  for (int i=0; i<n_gr_; i++) {
    rate = kgr_(i) * y[ingr1_(i)];
    ydotg[ingr1_(i)] -= rate;
    if (ingr2_(i) >= 0) {
      ydotg[ingr2_(i)] -= rate;
    }
    ydotg[outgr_(i)] += rate;
  }

  // special reactions
  for (int i=0; i<n_sr_; i++) {
    rate = ksr_(i);
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

  // grain collision reactions
  for (int i=0; i<n_gc_; i++) {
    rate =  kgc_(i) * y[ingc_(i, 0)] * y[ingc_(i, 1)] * nH_;
    if (y[ingc_(i, 0)] < 0 && y[ingc_(i, 1)] < 0) {
      rate *= -1.;
    }
    for (int jin=0; jin<n_ingc_; jin++) {
      if (ingc_(i, jin) >= 0) {
        ydotg[ingc_(i, jin)] -= rate;
      }
    }
    for (int jout=0; jout<n_outgc_; jout++) {
      if (outgc_(i, jout) >= 0) {
        ydotg[outgc_(i, jout)] += rate;
      }
    }
  }

  // set ydot to return
  for (int i=0; i<NSPECIES; i++) {
    // return in code units
    ydot[i] = ydotg[i] * pmy_mb_->pmy_mesh->punit->code_time_cgs;
  }

  // throw error if nan, or inf, or large value occurs
  for (int i=0; i<NSPECIES; i++) {
    if ( std::isnan(ydot[i]) || std::isinf(ydot[i]) ) {
      printf("ydot: ");
      for (int j=0; j<NSPECIES; j++) {
        printf("%s: %.2e  ", species_names[j].c_str(), ydot[j]);
      }
      printf("abundances: ");
      for (int j=0; j<NSPECIES; j++) {
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
      Real y_e;
      if (ispec_map_.find("e-") != ispec_map_.end()) {
        y_e = y0[ispec_map_["e-"]];
      } else {
        y_e = 0.;
      }
      printf("T = %.2e\n", E_ergs/Thermo::CvCold(y0[ispec_map_["H2"]], xHe_, y_e));
      std::stringstream msg;
      msg << "ChemNetwork (kida): RHS(ydot): nan or inf" << std::endl;
      ATHENA_ERROR(msg);
    }
  }
  delete[] y0;
  delete[] ydotg;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Real ChemNetwork::Edot(const Real t, const Real *y, const Real ED)
//! \brief energy equation dED/dt
//!
//! all input/output variables are in code units (ED is the energy density)

Real ChemNetwork::Edot(const Real t, const Real *y, const Real ED) {
  // ernergy per hydrogen atom
  Real E_ergs = ED * pmy_mb_->pmy_mesh->punit->code_energydensity_cgs / nH_;
  // isothermal
  if (!NON_BAROTROPIC_EOS) {
    return 0;
  }
  Real T = 0.;
  Real dEdt = 0.;
  Real y0[NSPECIES];
  // correct negative abundance to zero
  for (int i=0; i<NSPECIES; i++) {
    if (y[i] < 0) {
      y0[i] = 0;
    } else {
      y0[i] = y[i];
    }
  }
  // abundances
  const Real y_H2 = y0[ispec_map_["H2"]];
  const Real y_H = y0[ispec_map_["H"]];
  Real y_e, y_He, y_Hplus;
  // Real kcr_H;
  // Real kcr_H2, kcr_He;
  Real kgr_H;
  Real kph_H2, y_PAH;
  // explicit PAH abundance for heating and cooling
  if (!is_fixed_PAH_) {
    y_PAH = y0[ispec_map_["PAH0"]] + y0[ispec_map_["PAH+"]] + y0[ispec_map_["PAH-"]];
    Z_PAH_ = y_PAH/6e-7;
  }
  if (ispec_map_.find("e-") != ispec_map_.end()) {
    y_e = y0[ispec_map_["e-"]];
  } else {
    y_e = 0.;
  }
  if (ispec_map_.find("He") != ispec_map_.end()) {
    y_He = y0[ispec_map_["He"]];
  } else {
    y_He = 0.;
  }
  if (ispec_map_.find("H+") != ispec_map_.end()) {
    y_Hplus = y0[ispec_map_["H+"]];
  } else {
    y_Hplus = 0.;
  }
  // if (icr_H_ >= 0) {
  //   kcr_H = kcr_(icr_H_);
  // } else {
  //   kcr_H = 0.;
  // }
  // if (icr_H2_ >= 0) {
  //   kcr_H2 = kcr_(icr_H2_);
  // } else {
  //   kcr_H2 = 0.;
  // }
  // if (icr_He_ >= 0) {
  //   kcr_He = kcr_(icr_He_);
  // } else {
  //   kcr_He = 0.;
  // }
  if (igr_H_ >= 0) {
    kgr_H = kgr_(igr_H_);
  } else {
    kgr_H = 0.;
  }
  if (iph_H2_ >= 0) {
    kph_H2 = kph_(iph_H2_);
  } else {
    kph_H2 = 0.;
  }
  // temperature
  T = E_ergs / Thermo::CvCold(y_H2, xHe_, y_e);
  // apply temperature floor, incase of very small or negative energy
  if (T < temp_min_rates_) {
    T = temp_min_rates_;
  }


  // --------------------------heating-----------------------------
  Real GCR, GPE, GH2gr, k_xH2_photo, GH2pump, GH2diss;
  // CR heating
  GCR = Thermo::HeatingCr(y_e,  nH_, y_H, y_H2, rad_(index_cr_));
  // photo electric effect on dust
  if (is_fixed_PAH_) {
    GPE = Thermo::HeatingPE(rad_(index_gpe_), Z_PAH_, T, nH_*y_e);
  } else {
    GPE = Thermo::HeatingPE_W03(rad_(index_gpe_), Z_PAH_, T, nH_*y_e, phi_PAH_);
  }
  // H2 formation on dust grains
  k_xH2_photo = kph_H2;
  GH2gr = Thermo::HeatingH2gr(y_H,  y_H2, nH_, T, kgr_H, k_xH2_photo);
  // H2 UV pumping
  GH2pump = Thermo::HeatingH2pump(y_H,  y_H2, nH_, T, k_xH2_photo);
  // H2 photo dissiociation.
  GH2diss = Thermo::HeatingH2diss(k_xH2_photo, y_H2);

  // --------------------------cooling-----------------------------
  Real LCII, LCI, LOI, LHotGas, LCOR, LH2, LDust, LRec, LH2diss, LHIion;
  Real vth, nCO, grad_small;
  Real NCOeff, gradeff;
  Real k2body_H2_H, k2body_H2_H2, k2body_H_e;
  Real Tcool_nm;
  if (i2body_H2_H_ >= 0) {
    k2body_H2_H = k2body_(i2body_H2_H_);
  } else {
    k2body_H2_H = 0.;
  }
  if (T > temp_max_cool_nm_) {
    Tcool_nm = temp_max_cool_nm_;
  } else {
    Tcool_nm = T;
  }
  if (i2body_H2_H2_ >= 0) {
    k2body_H2_H2 = k2body_(i2body_H2_H2_);
  } else {
    k2body_H2_H2 = 0.;
  }
  if (i2body_H_e_ >= 0) {
    k2body_H_e = k2body_(i2body_H_e_);
  } else {
    k2body_H_e = 0.;
  }
  if (T < temp_min_cool_) {
    LCII = 0.;
    LCI = 0;
    LOI = 0.;
    LHotGas = 0;
    LCOR = 0;
    LH2 = 0;
    LDust = 0;
    LRec = 0;
    LH2diss = 0;
    LHIion = 0;
  } else {
    // C+ fine structure line
    if (ispec_map_.find("C+") != ispec_map_.end()) {
      LCII = Thermo::CoolingCII(y0[ispec_map_["C+"]],
                                nH_*y_H,  nH_*y_H2, nH_*y_e, Tcool_nm);
    } else {
      LCII = 0.;
    }
    // CI fine structure line
    if (ispec_map_.find("C") != ispec_map_.end()) {
      LCI = Thermo::CoolingCI(y0[ispec_map_["C"]], nH_*y_H, nH_*y_H2, nH_*y_e, Tcool_nm);
    } else {
      LCI = 0.;
    }
    // OI fine structure line
    if (ispec_map_.find("O") != ispec_map_.end()) {
      LOI = Thermo::CoolingOI(y0[ispec_map_["O"]], nH_*y_H, nH_*y_H2, nH_*y_e, Tcool_nm);
    } else {
      LOI = 0.;
    }
    // cooling of hot gas: radiative cooling, free-free.
    LHotGas = Thermo::CoolingLya(y_H, nH_*y_e, T); // Thermo::CoolingHotGas(nH_, T, Z_g_);
    // CO rotational lines
    if (ispec_map_.find("CO") != ispec_map_.end()) {
      // Calculate effective CO column density
      Real y_CO = y0[ispec_map_["CO"]];
      vth = std::sqrt(2. * Constants::k_boltzmann_cgs * Tcool_nm / ChemistryUtility::mCO);
      nCO = nH_ * y_CO;
      grad_small = vth/Leff_CO_max_;
      gradeff = std::max(gradv_, grad_small);
      NCOeff = nCO / gradeff;
      LCOR = Thermo::CoolingCOR(y_CO, nH_*y_H,  nH_*y_H2, nH_*y_e, Tcool_nm, NCOeff);
    } else {
      LCOR = 0.;
    }
    // H2 vibration and rotation lines
    LH2 = Thermo::CoolingH2(y_H2, nH_*y_H, nH_*y_H2, nH_*y_He,
                            nH_*y_Hplus, nH_*y_e, Tcool_nm);
    // dust thermo emission
    LDust = 0.; // Thermo::CoolingDustTd(Z_d_, nH_, Tcool_nm, temp_dust_thermo_);
    // reconbination of e on PAHs
    if (is_fixed_PAH_) {
      LRec = Thermo::CoolingRec(Z_PAH_, Tcool_nm, nH_*y_e, rad_(index_gpe_));
    } else {
      LRec = Thermo::CoolingRec_W03(Z_PAH_, Tcool_nm, nH_*y_e, rad_(index_gpe_),
                                    phi_PAH_);
    }
    // collisional dissociation of H2
    LH2diss = Thermo::CoolingH2diss(y_H, y_H2, k2body_H2_H, k2body_H2_H2);
    // collisional ionization of HI
    LHIion = Thermo::CoolingHIion(y_H,  y_e, k2body_H_e);
  }

  dEdt = (GCR + GPE + GH2gr + GH2pump + GH2diss)
         - (LCII + LCI + LOI + LHotGas + LCOR
            + LH2 + LDust + LRec + LH2diss + LHIion);
  // return in code units
  Real dEDdt = dEdt * nH_ / pmy_mb_->pmy_mesh->punit->code_energydensity_cgs
               * pmy_mb_->pmy_mesh->punit->code_time_cgs;
  if ( std::isnan(dEdt) || std::isinf(dEdt) ) {
    if ( std::isnan(LCOR) || std::isinf(LCOR) ) {
      printf("NCOeff=%.2e, gradeff=%.2e, gradv_=%.2e, vth=%.2e, nH_=%.2e, nCO=%.2e\n",
             NCOeff, gradeff, gradv_, vth, nH_, nCO);
    }
    printf("GCR=%.2e, GPE=%.2e, GH2gr=%.2e, GH2pump=%.2e GH2diss=%.2e\n",
           GCR , GPE , GH2gr , GH2pump , GH2diss);
    printf("LCII=%.2e, LCI=%.2e, LOI=%.2e, LHotGas=%.2e, LCOR=%.2e\n",
           LCII , LCI , LOI , LHotGas , LCOR);
    printf("LH2=%.2e, LDust=%.2e, LRec=%.2e, LH2diss=%.2e, LHIion=%.2e\n",
           LH2 , LDust , LRec , LH2diss , LHIion);
    printf("T=%.2e, dEdt=%.2e, E=%.2e, dEergsdt=%.2e, E_ergs=%.2e, Cv=%.2e, nH=%.2e\n",
           T, dEDdt, ED, dEdt, E_ergs, Thermo::CvCold(y_H2, xHe_, y_e), nH_);
    for (int i=0; i<NSPECIES; i++) {
      printf("%s: %.2e  ", species_names[i].c_str(), y0[i]);
    }
    printf("\n");
    std::stringstream msg;
    msg << "ChemNetwork (kida): dEdt: nan or inf number" << std::endl;
    ATHENA_ERROR(msg);
  }
  if (DEBUG) {
    if (output_thermo) {
      printf("NCOeff=%.2e, gradeff=%.2e, gradv_=%.2e, vth=%.2e, nH_=%.2e, nCO=%.2e\n",
             NCOeff, gradeff, gradv_, vth, nH_, nCO);
      printf("GCR=%.2e, GPE=%.2e, GH2gr=%.2e, GH2pump=%.2e GH2diss=%.2e\n",
             GCR , GPE , GH2gr , GH2pump , GH2diss);
      printf("LCII=%.2e, LCI=%.2e, LOI=%.2e, LHotGas=%.2e, LCOR=%.2e\n",
             LCII , LCI , LOI , LHotGas , LCOR);
      printf("LH2=%.2e, LDust=%.2e, LRec=%.2e, LH2diss=%.2e, LHIion=%.2e\n",
             LH2 , LDust , LRec , LH2diss , LHIion);
      printf("T=%.2e, dEdt=%.2e, E=%.2e, dEergsdt=%.2e, E_ergs=%.2e, Cv=%.2e, nH=%.2e\n",
             T, dEDdt, ED, dEdt, E_ergs, Thermo::CvCold(y_H2, xHe_, y_e), nH_);
      for (int i=0; i<NSPECIES; i++) {
        printf("%s: %.2e  ", species_names[i].c_str(), y0[i]);
      }
      printf("\n");
      output_thermo = false;
    }
  }
  return dEDdt;
}


//----------------------------------------------------------------------------------------
//! \fn void ChemNetwork::Jacobian_isothermal(const Real t, const Real *y,
//!           const Real *ydot, AthenaArray<Real> &jac)
//! \brief Jacobian for isothermal EOS

void ChemNetwork::Jacobian_isothermal(const Real t, const Real *y,
                                      const Real *ydot,
                                      AthenaArray<Real> &jac) {
  Real rate = 0;
  Real y_ices, Nl; // for desorption reactions
  const Real eps = 1e-50;
  const int i_H2 = ispec_map_["H2"];
  const Real xi = rad_(index_cr_);
  // initialize
  for (int i=0; i<NSPECIES; i++) {
    for (int j=0; j<NSPECIES; j++) {
      jac(i,j) = 0.;
    }
  }

  // cosmic ray reactions
  for (int i=0; i<n_cr_; i++) {
    rate = kcr_(i);
    jac(incr_(i),incr_(i)) -= rate;
    jac(outcr1_(i),incr_(i)) += rate;
    jac(outcr2_(i),incr_(i)) += rate;
    if (outcr3_(i) >= 0) {
      jac(outcr3_(i),incr_(i)) += rate;
    }
  }

  // cosmic ray induced photo reactions
  for (int i=0; i<n_crp_; i++) {
    rate = kcrp_(i);
    jac(incrp_(i),incrp_(i)) -= rate;
    jac(outcrp1_(i),incrp_(i)) += rate;
    jac(outcrp2_(i),incrp_(i)) += rate;
    rate = y[incrp_(i)] * kcrp_base_(i) * xi * 2.;
    jac(incrp_(i),i_H2) -= rate;
    jac(outcrp1_(i),i_H2) += rate;
    jac(outcrp2_(i),i_H2) += rate;
  }

  // FUV photo-dissociation and photo-ionisation
  for (int i=0; i<n_ph_; i++) {
    rate = kph_(i);
    jac(inph_(i),inph_(i)) -= rate;
    jac(outph1_(i),inph_(i)) += rate;
    jac(outph2_(i),inph_(i)) += rate;
  }

  // 2body reactions
  for (int i=0; i<n_2body_; i++) {
    // df/dy1
    rate =  k2body_(i) * y[in2body_(i, 1)] * nH_;
    if (y[in2body_(i, 0)] < 0 && y[in2body_(i, 1)] < 0) {
      rate *= -1.;
    }
    for (int jin=0; jin<n_in2body_; jin++) {
      if (in2body_(i, jin) >= 0) {
        jac(in2body_(i, jin),in2body_(i, 0)) -= rate;
      }
    }
    for (int jout=0; jout<n_out2body_; jout++) {
      if (out2body_(i, jout) >= 0) {
        jac(out2body_(i, jout),in2body_(i, 0)) += rate;
      }
    }
    // df/dy2
    rate =  k2body_(i) * y[in2body_(i, 0)] * nH_;
    if (y[in2body_(i, 0)] < 0 && y[in2body_(i, 1)] < 0) {
      rate *= -1.;
    }
    for (int jin=0; jin<n_in2body_; jin++) {
      if (in2body_(i, jin) >= 0) {
        jac(in2body_(i, jin),in2body_(i, 1)) -= rate;
      }
    }
    for (int jout=0; jout<n_out2body_; jout++) {
      if (out2body_(i, jout) >= 0) {
        jac(out2body_(i, jout),in2body_(i, 1)) += rate;
      }
    }
  }

  // 2bodytr reactions
  for (int i=0; i<n_2bodytr_; i++) {
    // df/dy1
    rate =  k2bodytr_(i) * y[in2bodytr2_(i)] * nH_;
    if (y[in2bodytr1_(i)] < 0 && y[in2bodytr2_(i)] < 0) {
      rate *= -1.;
    }
    jac(in2bodytr1_(i),in2bodytr1_(i)) -= rate;
    jac(in2bodytr2_(i),in2bodytr1_(i)) -= rate;
    jac(out2bodytr1_(i),in2bodytr1_(i)) += rate;
    if (out2bodytr2_(i) >= 0) {
      jac(out2bodytr2_(i),in2bodytr1_(i)) += rate;
    }
    if (out2bodytr3_(i) >= 0) {
      jac(out2bodytr3_(i),in2bodytr1_(i)) += rate;
    }
    // df/dy2
    rate =  k2bodytr_(i) * y[in2bodytr1_(i)] * nH_;
    if (y[in2bodytr1_(i)] < 0 && y[in2bodytr2_(i)] < 0) {
      rate *= -1.;
    }
    jac(in2bodytr1_(i),in2bodytr2_(i)) -= rate;
    jac(in2bodytr2_(i),in2bodytr2_(i)) -= rate;
    jac(out2bodytr1_(i),in2bodytr2_(i)) += rate;
    if (out2bodytr2_(i) >= 0) {
      jac(out2bodytr2_(i),in2bodytr2_(i)) += rate;
    }
    if (out2bodytr3_(i) >= 0) {
      jac(out2bodytr3_(i),in2bodytr2_(i)) += rate;
    }
  }

  // grain assisted reactions
  // desorption reactions: dependence on ice thickness
  y_ices = 0; // total ice abundance
  for (int i=0; i<nices_; i++) {
    y_ices += y[id_ices_(i)];
  }
  if (x_d_ < eps) { // control for very small dust abundance
    Nl = 0.;
  } else {
    Nl = y_ices / (6.0e15*M_PI*a_d_*a_d_*x_d_); // number of layers
  }
  for (int i=0; i<n_gr_; i++) {
    rate = kgr_(i);
    jac(ingr1_(i),ingr1_(i)) -= rate;
    if (ingr2_(i) >= 0) {
      jac(ingr2_(i),ingr1_(i)) -= rate;
    }
    jac(outgr_(i),ingr1_(i)) += rate;
    if (frml_gr_(i) == 10 && Nl > 1.) {// desorption
      rate = - kgr_(i)*y[ingr1_(i)]/y_ices;
      for (int j=0; j<nices_; j++) {
        jac(ingr1_(i),id_ices_(j)) -= rate;
        if (ingr2_(i) >= 0) {
          jac(ingr2_(i),id_ices_(j)) -= rate;
        }
        jac(outgr_(i),id_ices_(j)) += rate;
      }
    }
  }

  // grain collision reactions
  for (int i=0; i<n_gc_; i++) {
    // df/dy1
    rate =  kgc_(i) * y[ingc_(i, 1)] * nH_;
    if (y[ingc_(i, 0)] < 0 && y[ingc_(i, 1)] < 0) {
      rate *= -1.;
    }
    for (int jin=0; jin<n_ingc_; jin++) {
      if (ingc_(i, jin) >= 0) {
        jac(ingc_(i, jin),ingc_(i, 0)) -= rate;
      }
    }
    for (int jout=0; jout<n_outgc_; jout++) {
      if (outgc_(i, jout) >= 0) {
        jac(outgc_(i, jout),ingc_(i, 0)) += rate;
      }
    }
    // df/dy2
    rate =  kgc_(i) * y[ingc_(i, 0)] * nH_;
    if (y[ingc_(i, 0)] < 0 && y[ingc_(i, 1)] < 0) {
      rate *= -1.;
    }
    for (int jin=0; jin<n_ingc_; jin++) {
      if (ingc_(i, jin) >= 0) {
        jac(ingc_(i, jin),ingc_(i, 1)) -= rate;
      }
    }
    for (int jout=0; jout<n_outgc_; jout++) {
      if (outgc_(i, jout) >= 0) {
        jac(outgc_(i, jout),ingc_(i, 1)) += rate;
      }
    }
  }

  // special reactions with frml=7
  UpdateJacobianSpecial(y, 0., jac);

  // set unit for jacobian
  for (int i=0; i<NSPECIES; i++) {
    for (int j=0; j<NSPECIES; j++) {
      jac(i,j) *= pmy_mb_->pmy_mesh->punit->code_time_cgs;
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ChemNetwork::InitializeReactions()
//! \brief set up coefficients of chemical reactions

void ChemNetwork::InitializeReactions() {
  KidaReaction *pr = NULL;
  // error message
  bool error=false;
  ReactionType rtype;
  // count reactions
  for (int ir=0; ir<nr_; ir++) {
    CheckReaction(reactions_[ir]);
    pr = &reactions_[ir];
    rtype = SortReaction(pr);
    switch(rtype) {
      case ReactionType::cr: n_cr_++; break;
      case ReactionType::crp: n_crp_++; break;
      case ReactionType::photo: n_ph_++; break;
      case ReactionType::twobody: n_2body_++; break;
      case ReactionType::twobodytr: break;
      case ReactionType::grain_implicit: n_gr_++; break;
      case ReactionType::special: n_sr_++; break;
      case ReactionType::grain_collision: n_gc_++; break;
      default: std::stringstream msg;
        msg << "### FATAL ERROR in ChemNetwork InitializeReactions()"
            << " [ChemNetwork]: reaction type not recognized." << std::endl;
        ATHENA_ERROR(msg);
        break;
    }
  }
  // create arrays
  if (n_cr_ > 0) {
    incr_.NewAthenaArray(n_cr_);
    outcr1_.NewAthenaArray(n_cr_);
    outcr2_.NewAthenaArray(n_cr_);
    outcr3_.NewAthenaArray(n_cr_);
    kcr_base_.NewAthenaArray(n_cr_);
    kcr_.NewAthenaArray(n_cr_);
  }
  if (n_crp_ > 0) {
    incrp_.NewAthenaArray(n_crp_);
    outcrp1_.NewAthenaArray(n_crp_);
    outcrp2_.NewAthenaArray(n_crp_);
    kcrp_base_.NewAthenaArray(n_crp_);
    kcrp_.NewAthenaArray(n_crp_);
  }
  if (n_ph_ > 0) {
    inph_.NewAthenaArray(n_ph_);
    outph1_.NewAthenaArray(n_ph_);
    outph2_.NewAthenaArray(n_ph_);
    kph_base_.NewAthenaArray(n_ph_);
    kph_avfac_.NewAthenaArray(n_ph_);
    kph_.NewAthenaArray(n_ph_);
  }
  if (n_2body_ > 0) {
    in2body_.NewAthenaArray(n_2body_, n_in2body_);
    out2body_.NewAthenaArray(n_2body_, n_out2body_);
    frml_2body_.NewAthenaArray(n_2body_);
    a2body_.NewAthenaArray(n_2body_);
    b2body_.NewAthenaArray(n_2body_);
    c2body_.NewAthenaArray(n_2body_);
    Tmin_2body_.NewAthenaArray(n_2body_);
    Tmax_2body_.NewAthenaArray(n_2body_);
    k2body_.NewAthenaArray(n_2body_);
  }
  if (n_2bodytr_ > 0) {
    in2bodytr1_.NewAthenaArray(n_2bodytr_);
    in2bodytr2_.NewAthenaArray(n_2bodytr_);
    out2bodytr1_.NewAthenaArray(n_2bodytr_);
    out2bodytr2_.NewAthenaArray(n_2bodytr_);
    out2bodytr3_.NewAthenaArray(n_2bodytr_);
    nr_2bodytr_.NewAthenaArray(n_2bodytr_);
    frml_2bodytr_.NewAthenaArray(n_2bodytr_, n_range_);
    a2bodytr_.NewAthenaArray(n_2bodytr_, n_range_);
    b2bodytr_.NewAthenaArray(n_2bodytr_, n_range_);
    c2bodytr_.NewAthenaArray(n_2bodytr_, n_range_);
    Tmin_2bodytr_.NewAthenaArray(n_2bodytr_, n_range_);
    Tmax_2bodytr_.NewAthenaArray(n_2bodytr_, n_range_);
    k2bodytr_.NewAthenaArray(n_2bodytr_);
  }
  if (n_gr_ > 0) {
    ingr1_.NewAthenaArray(n_gr_);
    ingr2_.NewAthenaArray(n_gr_);
    outgr_.NewAthenaArray(n_gr_);
    frml_gr_.NewAthenaArray(n_gr_);
    kgr_.NewAthenaArray(n_gr_);
    TDgr_.NewAthenaArray(n_gr_);
    nu0gr_.NewAthenaArray(n_gr_);
  }
  if (n_sr_ > 0) {
    insr_.NewAthenaArray(n_sr_, n_insr_);
    outsr_.NewAthenaArray(n_sr_, n_outsr_);
    ksr_.NewAthenaArray(n_sr_);
  }
  if (n_gc_ > 0) {
    ingc_.NewAthenaArray(n_gc_, n_ingc_);
    outgc_.NewAthenaArray(n_gc_, n_outgc_);
    nu_gc_.NewAthenaArray(n_gc_);
    r1_gc_.NewAthenaArray(n_gc_);
    t1_gc_.NewAthenaArray(n_gc_);
    kgc_.NewAthenaArray(n_gc_);
  }

  int icr=0, icrp=0, iph=0, i2body=0, i2bodytr=0, igr=0, isr=0, igc=0;
  std::vector<int> idtr;
  for (int ir=0; ir<nr_; ir++) {
    pr = &reactions_[ir];
    rtype = SortReaction(pr);
    // ---------------- cr - direct cosmic-ray ionization --------------
    if (rtype == ReactionType::cr) {
      std::string in_spec; // input species
      if (pr->reactants_[0] == "CR") {
        in_spec = pr->reactants_[1];
      } else {
        in_spec = pr->reactants_[0];
      }
      if (in_spec == "H") {
        icr_H_ = icr;
      } else if (in_spec == "H2") {
        icr_H2_ = icr;
      } else if (in_spec == "He") {
        icr_He_ = icr;
      }
      if (pr->formula_ == 1) {
        incr_(icr) = ispec_map_[in_spec];
        outcr1_(icr) = ispec_map_[ pr->products_[0]];
        outcr2_(icr) = ispec_map_[ pr->products_[1]];
        if (pr->products_.size() == 3) {
          outcr3_(icr) = ispec_map_[ pr->products_[2]];
        } else {
          outcr3_(icr) = -1;
        }
        kcr_base_(icr) = pr->alpha_;
        kcr_(icr) = 0.;
        icr++;
      } else if (pr->formula_ == 7) {
        incr_(icr) = ispec_map_[in_spec];
        outcr1_(icr) = ispec_map_[ pr->products_[0]];
        outcr2_(icr) = ispec_map_[ pr->products_[1]];
        if (pr->products_.size() == 3) {
          outcr3_(icr) = ispec_map_[ pr->products_[2]];
        } else {
          outcr3_(icr) = -1;
        }
        id7map_(pr->id_) = icr;
        id7type_(pr->id_) = ReactionType::cr;
        kcr_base_(icr) = 0.;
        kcr_(icr) = 0.;
        icr++;
      } else {
        error = true;
      }

      // ---------------- crp - cosmic-ray induced photo ionization --------
    } else if (rtype == ReactionType::crp) {
      std::string in_spec; // input species
      if (pr->reactants_[0] == "CRP") {
        in_spec = pr->reactants_[1];
      } else {
        in_spec = pr->reactants_[0];
      }
      if (pr->formula_ == 1) {
        incrp_(icrp) = ispec_map_[in_spec];
        outcrp1_(icrp) = ispec_map_[ pr->products_[0]];
        outcrp2_(icrp) = ispec_map_[ pr->products_[1]];
        kcrp_base_(icrp) = pr->alpha_;
        kcrp_(icrp) = 0.;
        icrp++;
      } else if (pr->formula_ == 7) {
        incrp_(icrp) = ispec_map_[in_spec];
        outcrp1_(icrp) = ispec_map_[ pr->products_[0]];
        outcrp2_(icrp) = ispec_map_[ pr->products_[1]];
        id7map_(pr->id_) = icrp;
        id7type_(pr->id_) = ReactionType::crp;
        kcrp_base_(icrp) = 0.;
        kcrp_(icrp) = 0.;
        icrp++;
      } else {
        error = true;
      }

      // ---------------- photo - FUV ionization/dissociation ----------------
    } else if (rtype == ReactionType::photo) {
      std::string in_spec; // input species
      if (pr->reactants_[0] == "Photon") {
        in_spec = pr->reactants_[1];
      } else {
        in_spec = pr->reactants_[0];
      }
      if (in_spec == "H2" && pr->products_[0] == "H" && pr->products_[1] == "H") {
        iph_H2_ = iph;
      }
      if (pr->formula_ == 2) {
        inph_(iph) = ispec_map_[in_spec];
        outph1_(iph) = ispec_map_[ pr->products_[0]];
        outph2_(iph) = ispec_map_[ pr->products_[1]];
        kph_base_(iph) = pr->alpha_;
        kph_avfac_(iph) = pr->gamma_;
        smap_ph_[in_spec] = iph;
        kph_(iph) = 0.;
        iph++;
      } else {
        error = true;
      }

      // ---------------- twobody - 2body reaction ---------------------------
    } else if (rtype == ReactionType::twobody) {
      if (pr->reactants_[0] == "H2" && pr->reactants_[1] == "H") {
        i2body_H2_H_ = i2body;
      }
      if (pr->reactants_[0] == "H2" && pr->reactants_[1] == "H2") {
        i2body_H2_H2_ = i2body;
      }
      if (pr->reactants_[0] == "H" && pr->reactants_[1] == "e-") {
        i2body_H_e_ = i2body;
      }
      if (pr->formula_ == 3 || pr->formula_ == 4 || pr->formula_ == 5) {
        for (int jin=0; jin<n_in2body_; jin++) {
          if (jin < static_cast<int>(pr->reactants_.size())) {
            in2body_(i2body, jin) = ispec_map_[pr->reactants_[jin]];
          } else {
            in2body_(i2body, jin) = -1;
          }
        }
        for (int jout=0; jout<n_out2body_; jout++) {
          if (jout < static_cast<int>(pr->products_.size())) {
            out2body_(i2body, jout) = ispec_map_[pr->products_[jout]];
          } else {
            out2body_(i2body, jout) = -1;
          }
        }
        frml_2body_(i2body) = pr->formula_;
        a2body_(i2body) = pr->alpha_;
        b2body_(i2body) = pr->beta_;
        c2body_(i2body) = pr->gamma_;
        Tmin_2body_(i2body) = pr->Tmin_;
        Tmax_2body_(i2body) = pr->Tmax_;
        k2body_(i2body) = 0.;
        i2body++;
      } else if (pr->formula_ == 7) {
        for (int jin=0; jin<n_in2body_; jin++) {
          if (jin < static_cast<int>(pr->reactants_.size())) {
            in2body_(i2body, jin) = ispec_map_[pr->reactants_[jin]];
          } else {
            in2body_(i2body, jin) = -1;
          }
        }
        for (int jout=0; jout<n_out2body_; jout++) {
          if (jout < static_cast<int>(pr->products_.size())) {
            out2body_(i2body, jout) = ispec_map_[pr->products_[jout]];
          } else {
            out2body_(i2body, jout) = -1;
          }
        }
        frml_2body_(i2body) = pr->formula_;
        a2body_(i2body) = 0.;
        b2body_(i2body) = 0.;
        c2body_(i2body) = 0.;
        Tmin_2body_(i2body) = pr->Tmin_;
        Tmax_2body_(i2body) = pr->Tmax_;
        k2body_(i2body) = 0.;
        id7map_(pr->id_) = i2body;
        id7type_(pr->id_) = ReactionType::twobody;
        i2body++;
      } else {
        error = true;
      }
      // ------- twobodytr - 2body reaction with temperature range ----------
    } else if (rtype == ReactionType::twobodytr) {
      const int nprev = std::count(idtr.begin(), idtr.end(), pr->id_);
      if (nprev == 0) {
        in2bodytr1_(i2bodytr) = ispec_map_[ pr->reactants_[0]];
        in2bodytr2_(i2bodytr) = ispec_map_[ pr->reactants_[1]];
        out2bodytr1_(i2bodytr) = ispec_map_[ pr->products_[0]];
        if (pr->products_.size() >= 2) {
          out2bodytr2_(i2bodytr) = ispec_map_[ pr->products_[1]];
        } else {
          out2bodytr2_(i2bodytr) = -1;
        }
        if (pr->products_.size() >= 3) {
          out2bodytr3_(i2bodytr) = ispec_map_[ pr->products_[2]];
        } else {
          out2bodytr3_(i2bodytr) = -1;
        }
        k2bodytr_(i2bodytr) = 0.;
        i2bodytr++;
      }
      nr_2bodytr_(i2bodytr-1) = nprev + 1;
      frml_2bodytr_(i2bodytr-1, nprev) = pr->formula_;
      a2bodytr_(i2bodytr-1, nprev) = pr->alpha_;
      b2bodytr_(i2bodytr-1, nprev) = pr->beta_;
      c2bodytr_(i2bodytr-1, nprev) = pr->gamma_;
      Tmin_2bodytr_(i2bodytr-1, nprev) = pr->Tmin_;
      Tmax_2bodytr_(i2bodytr-1, nprev) = pr->Tmax_;
      idtr.push_back(pr->id_);

      // ----------- grain_implicit - implicit grain assisted reaction -------
    } else if (rtype == ReactionType::grain_implicit) {
      if (pr->reactants_[0] == "H" && pr->reactants_[1] == "H") {
        igr_H_ = igr;
      }
      if (pr->formula_ == 7) { // special rates
        ingr1_(igr) = ispec_map_[pr->reactants_[0]];
        if (pr->reactants_.size() == 2) {
          ingr2_(igr) = ispec_map_[pr->reactants_[1]];
        } else {
          ingr2_(igr) = -1;
        }
        outgr_(igr) = ispec_map_[pr->products_[0]];
        id7map_(pr->id_) = igr;
        id7type_(pr->id_) = ReactionType::grain_implicit;
        frml_gr_(igr) = pr->formula_;
        kgr_(igr) = 0.;
        igr++;
      } else if (pr->formula_ == 10 && pr->reactants_.size() == 1) { // desorption
        const Real mi = species_[ispec_map_[pr->reactants_[0]]].mass_;
        ingr1_(igr) = ispec_map_[pr->reactants_[0]];
        ingr2_(igr) = -1;
        outgr_(igr) = ispec_map_[pr->products_[0]];
        frml_gr_(igr) = pr->formula_;
        kgr_(igr) = 0.;
        TDgr_(igr) = pr->gamma_;
        nu0gr_(igr) = std::sqrt( 3.0e15*pr->gamma_*Constants::k_boltzmann_cgs
                                 /(M_PI*M_PI*mi) );
        igr++;
      } else {
        error = true;
      }

      // ------------------ special - special reactions -----------------------
    } else if (rtype == ReactionType::special) {
      if (pr->formula_ == 7) {
        id7map_(pr->id_) = isr;
        id7type_(pr->id_) = ReactionType::special;
        ksr_(isr) = 0.;
        for (int jin=0; jin<n_insr_; jin++) {
          if (jin < static_cast<int>(pr->reactants_.size())) {
            insr_(isr, jin) = ispec_map_[pr->reactants_[jin]];
          } else {
            insr_(isr, jin) = -1;
          }
        }
        for (int jout=0; jout<n_outsr_; jout++) {
          if (jout < static_cast<int>(pr->products_.size())) {
            outsr_(isr, jout) = ispec_map_[pr->products_[jout]];
          } else {
            outsr_(isr, jout) = -1;
          }
        }
        isr++;
      } else {
        error = true;
      }

      // ------------- grain_collision - grain collisional reactions -----------
    } else if (rtype == ReactionType::grain_collision) {
      if (pr->formula_ == 8) { // electron and ion
        const Real br = pr->alpha_;
        const Real se = pr->beta_;
        const Real ag = pr->gamma_;
        const Real mi = species_[ispec_map_[pr->reactants_[1]]].mass_;
        const Real q_charge = species_[ispec_map_[pr->reactants_[1]]].charge_;
        const Real qi = q_charge * Constants::echarge_cgs;
        for (int jin=0; jin<n_ingc_; jin++) {
          if (jin < static_cast<int>(pr->reactants_.size())) {
            ingc_(igc, jin) = ispec_map_[pr->reactants_[jin]];
          } else {
            ingc_(igc, jin) = -1;
          }
        }
        for (int jout=0; jout<n_outgc_; jout++) {
          if (jout < static_cast<int>(pr->products_.size())) {
            outgc_(igc, jout) = ispec_map_[pr->products_[jout]];
          } else {
            outgc_(igc, jout) = -1;
          }
        }
        kgc_(igc) = 0.;
        if (q_charge != 1 && q_charge != -1) {
          std::stringstream msg;
          msg << "### FATAL ERROR in ChemNetwork InitializeReactions() [ChemNetwork]"
              << std::endl
              << "grain collsion reaction: charge of electron/ion != +-1, "
              << "reaction ID=" << std::endl;
          ATHENA_ERROR(msg);
        }
        nu_gc_(igc) = species_[ispec_map_[pr->reactants_[0]]].charge_ / q_charge;
        r1_gc_(igc) = br*se* M_PI *ag*ag
                      * std::sqrt(8.*Constants::k_boltzmann_cgs/(M_PI*mi));
        t1_gc_(igc) = ag * Constants::k_boltzmann_cgs / (qi*qi);
        igc++;
      } else if (pr->formula_ == 9) { // neutral freeze-out
        const Real ag = pr->gamma_;
        const Real mi = species_[ispec_map_[pr->reactants_[1]]].mass_;
        for (int jin=0; jin<n_ingc_; jin++) {
          if (jin < static_cast<int>(pr->reactants_.size())) {
            ingc_(igc, jin) = ispec_map_[pr->reactants_[jin]];
          } else {
            ingc_(igc, jin) = -1;
          }
        }
        for (int jout=0; jout<n_outgc_; jout++) {
          if (jout < static_cast<int>(pr->products_.size())) {
            outgc_(igc, jout) = ispec_map_[pr->products_[jout]];
          } else {
            outgc_(igc, jout) = -1;
          }
        }
        kgc_(igc) = 0.;
        nu_gc_(igc) = 9; // flag for freeze-out reaction
        r1_gc_(igc) = M_PI*ag*ag*std::sqrt( 8.*Constants::k_boltzmann_cgs/(M_PI*mi) );
        t1_gc_(igc) = 0.;
        igc++;
      } else {
        error = true;
      }

      //------------------ formula not recogonized -------------------------
    } else {
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

  // sanity check
  if (icr != n_cr_ || icrp != n_crp_ || iph != n_ph_ || i2body != n_2body_
      || i2bodytr != n_2bodytr_ || igr != n_gr_ || isr != n_sr_ || igc != n_gc_) {
    std::stringstream msg;
    msg << "### FATAL ERROR in ChemNetwork InitializeReactions() [ChemNetwork]"
        << ": counts of reactions does not match." << std::endl;
    ATHENA_ERROR(msg);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ChemNetwork::UpdateRates(const Real *y, const Real E)
//! \brief update the rates for chemical reactions.
//!
//! This is called at the beginning of the RHS. E is in the unit of ergs.

void ChemNetwork::UpdateRates(const Real *y, const Real E) {
  const Real y_H2 = y[ispec_map_["H2"]];
  Real y_e;
  if (ispec_map_.find("e-") != ispec_map_.end()) {
    y_e = y[ispec_map_["e-"]];
  } else {
    y_e = 0.;
  }

  Real T, Tcap;
  if (NON_BAROTROPIC_EOS) {
    T = E / Thermo::CvCold(y_H2, xHe_, y_e);
  } else {
    // isothermal EOS
    T = temperature_;
  }
  // cap T above some minimum temperature
  if (T < temp_min_rates_) {
    T = temp_min_rates_;
  }

  // cosmic ray reactions
  for (int i=0; i<n_cr_; i++) {
    kcr_(i) = kcr_base_(i) * rad_(index_cr_);
  }

  // cosmic ray induced photo reactions
  for (int i=0; i<n_crp_; i++) {
    kcrp_(i) = kcrp_base_(i) * rad_(index_cr_) * 2.*y_H2;
  }

  // FUV reactions
  for (int i=0; i<n_ph_; i++) {
    kph_(i) = kph_base_(i) * rad_(i);
  }

  // grain collision reactions
  if (flag_T_rates_) {
    for (int i=0; i<n_gc_; i++) {
      if (nu_gc_(i) == 0) { // polarisation
        kgc_(i) = r1_gc_(i) * std::sqrt(T) * (1. + std::sqrt( M_PI/(2*t1_gc_(i)*T) ) );
      } else if (nu_gc_(i) == -1) { // Coulomb focusing
        kgc_(i) = r1_gc_(i) * std::sqrt(T) * (1. + 1./(t1_gc_(i)*T) )
                  * (1. + std::sqrt( 2./(2. + t1_gc_(i)*T) ) );
      } else if (nu_gc_(i) == 9) { // freeze-out
        kgc_(i) = r1_gc_(i) * std::sqrt(T);
      } else {
        std::stringstream msg;
        msg << "### fatal error in chemnetwork UpdateRates() [chemnetwork]: "
            << "grain collsion reaction type not implemented."
            << std::endl;
        ATHENA_ERROR(msg);
      }
    }
  }

  // freeze-out reactions
  Real kth, kcr, kcrp, kFUV, y_ices, Nl, fl;
  const Real eps = 1e-50;
  y_ices = 0; // total ice abundance
  for (int i=0; i<nices_; i++) {
    y_ices += y[id_ices_(i)];
  }
  // control for very small dust abundance
  if (x_d_ < eps) {
    fl = 0.;
  } else {
    Nl = y_ices / (6.0e15*M_PI*a_d_*a_d_*x_d_); // number of layers
    if (Nl <= 1.) {
      fl = 1.;
    } else {
      fl = 1./Nl;
    }
  }
  for (int i=0; i<n_gr_; i++) {
    if (frml_gr_(i) == 10) {
      kth = nu0gr_(i) * std::exp( -TDgr_(i)/T );
      kcr = 0.108 * rad_(index_cr_) * nu0gr_(i) * std::exp( -TDgr_(i)/70. );
      kcrp = 1.037e5 * rad_(index_cr_) * Yi_;
      kFUV = 3.23e-11 * rad_(index_gpe_);
      kgr_(i) = (kth + kcr + kcrp + kFUV) * fl;
    }
  }

  // 2body reactions
  if (flag_T_rates_) {
    if (is_Tcap_2body_) {
      for (int i=0; i<n_2body_; i++) {
        if (T < Tmin_2body_(i)) {
          Tcap = Tmin_2body_(i);
        } else if (T > Tmax_2body_(i)) {
          Tcap = Tmax_2body_(i);
        } else {
          Tcap = T;
        }
        if (frml_2body_(i) == 3) {
          k2body_(i) = a2body_(i)*std::pow(Tcap/300.,
                                           b2body_(i))*std::exp(-c2body_(i)/Tcap);
        } else if (frml_2body_(i) == 4) {
          k2body_(i) = a2body_(i)*b2body_(i)*( 0.62
                                               + 0.4767*c2body_(i)*std::sqrt(300./Tcap) );
        } else if (frml_2body_(i) == 5) {
          k2body_(i) = a2body_(i)*b2body_(i)*( 1 + 0.0967*c2body_(i)*std::sqrt(300./Tcap)
                                               + 28.501*c2body_(i)*c2body_(i)/Tcap );
        }
      }
    } else {
      for (int i=0; i<n_2body_; i++) {
        if (frml_2body_(i) == 3) {
          k2body_(i) = a2body_(i)*std::pow(T/300., b2body_(i))*std::exp(-c2body_(i)/T);
        } else if (frml_2body_(i) == 4) {
          k2body_(i) = a2body_(i)*b2body_(i)*( 0.62
                                               + 0.4767*c2body_(i)*std::sqrt(300./T) );
        } else if (frml_2body_(i) == 5) {
          k2body_(i) = a2body_(i)*b2body_(i)*( 1 + 0.0967*c2body_(i)*std::sqrt(300./T)
                                               + 28.501*c2body_(i)*c2body_(i)/T );
        }
      }
    }
  }

  // 2bodytr reactions
  if (flag_T_rates_) {
    if (is_Tcap_2body_) {
      for (int i=0; i<n_2bodytr_; i++) {
        int nr = nr_2bodytr_(i);
        int irange1 = 0;
        int irange2 = 0;
        Real rate1 = 0.;
        Real rate2 = 0.;
        if ( T < Tmin_2bodytr_(i,0) ) {
          Tcap = Tmin_2bodytr_(i,0);
        } else if ( T > Tmax_2bodytr_(i,nr-1) ) {
          Tcap = Tmax_2bodytr_(i,nr-1);
        } else {
          Tcap = T;
        }
        // select which temperature range to use
        if ( Tcap <= Tmax_2bodytr_(i,0) ) {
          irange1 = 0;
          irange2 = 0;
        } else if ( Tcap <= Tmin_2bodytr_(i,1) ) {
          irange1 = 0;
          irange2 = 1;
        } else if ( Tcap <= Tmax_2bodytr_(i,1) ) {
          irange1 = 1;
          irange2 = 1;
        } else {
          if (nr == 2) {
            irange1 = 1;
            irange2 = 1;
          } else if (nr == 3) {
            if ( Tcap <= Tmin_2bodytr_(i,2) ) {
              irange1 = 1;
              irange2 = 2;
            } else {
              irange1 = 2;
              irange2 = 2;
            }
          } else {
            std::stringstream msg;
            msg << "### fatal error in chemnetwork UpdateRates() [chemnetwork]: "
                << "2bodytr reaction with more than 3 temperature ranges not implemented."
                << std::endl;
            ATHENA_ERROR(msg);
          }
        }
        // calculate rates
        if (frml_2bodytr_(i,irange1) == 3) {
          rate1 = a2bodytr_(i,irange1)*std::pow(Tcap/300., b2bodytr_(i,irange1))
                  *std::exp(-c2bodytr_(i,irange1)/Tcap);
        } else if (frml_2bodytr_(i,irange1) == 4) {
          rate1 = a2bodytr_(i,irange1)*b2bodytr_(i,irange1)
                  *( 0.62 + 0.4767*c2bodytr_(i,irange1)*std::sqrt(300./Tcap) );
        } else if (frml_2bodytr_(i,irange1) == 5) {
          rate1 = a2bodytr_(i,irange1)*b2bodytr_(i,irange1)*(
              1 + 0.0967*c2bodytr_(i,irange1)*std::sqrt(300./Tcap)
              + 28.501*c2bodytr_(i,irange1)*c2bodytr_(i,irange1)/Tcap );
        }
        if (irange1 == irange2) {
          rate2 = rate1;
        } else {
          if (frml_2bodytr_(i,irange2) == 3) {
            rate2 = a2bodytr_(i,irange2)*std::pow(Tcap/300., b2bodytr_(i,irange2))
                    *std::exp(-c2bodytr_(i,irange2)/Tcap);
          } else if (frml_2bodytr_(i,irange2) == 4) {
            rate2 = a2bodytr_(i,irange2)*b2bodytr_(i,irange2)
                    *( 0.62 + 0.4767*c2bodytr_(i,irange2)*std::sqrt(300./Tcap) );
          } else if (frml_2bodytr_(i,irange2) == 5) {
            rate2 = a2bodytr_(i,irange2)*b2bodytr_(i,irange2)*(
                1 + 0.0967*c2bodytr_(i,irange2)*std::sqrt(300./Tcap)
                + 28.501*c2bodytr_(i,irange2)*c2bodytr_(i,irange2)/Tcap );
          }
        }
        // assign reaction rate
        k2bodytr_(i) = (rate1 + rate2) * 0.5;
      }
    } else {
      for (int i=0; i<n_2bodytr_; i++) {
        int nr = nr_2bodytr_(i);
        int irange1 = 0;
        int irange2 = 0;
        Real rate1 = 0.;
        Real rate2 = 0.;
        // select which temperature range to use
        if ( T <= Tmax_2bodytr_(i,0) ) {
          irange1 = 0;
          irange2 = 0;
        } else if ( T <= Tmin_2bodytr_(i,1) ) {
          irange1 = 0;
          irange2 = 1;
        } else if ( T <= Tmax_2bodytr_(i,1) ) {
          irange1 = 1;
          irange2 = 1;
        } else {
          if (nr == 2) {
            irange1 = 1;
            irange2 = 1;
          } else if (nr == 3) {
            if ( T <= Tmin_2bodytr_(i,2) ) {
              irange1 = 1;
              irange2 = 2;
            } else {
              irange1 = 2;
              irange2 = 2;
            }
          } else {
            std::stringstream msg;
            msg << "### fatal error in chemnetwork UpdateRates() [chemnetwork]: "
                << "2bodytr reaction with more than 3 temperature ranges not implemented."
                << std::endl;
            ATHENA_ERROR(msg);
          }
        }
        // calculate rates
        // TODO(Gong): duplicate code; modularize in a function?
        if (frml_2bodytr_(i,irange1) == 3) {
          rate1 = a2bodytr_(i,irange1)*std::pow(T/300., b2bodytr_(i,irange1))
                  *std::exp(-c2bodytr_(i,irange1)/T);
        } else if (frml_2bodytr_(i,irange1) == 4) {
          rate1 = a2bodytr_(i,irange1)*b2bodytr_(i,irange1)*
                  ( 0.62 + 0.4767*c2bodytr_(i,irange1)*std::sqrt(300./T) );
        } else if (frml_2bodytr_(i,irange1) == 5) {
          rate1 = a2bodytr_(i,irange1)*b2bodytr_(i,irange1)*(
              1 + 0.0967*c2bodytr_(i,irange1)*std::sqrt(300./T)
              + 28.501*c2bodytr_(i,irange1)*c2bodytr_(i,irange1)/T );
        }
        if (irange1 == irange2) {
          rate2 = rate1;
        } else {
          if (frml_2bodytr_(i,irange2) == 3) {
            rate2 = a2bodytr_(i,irange2)*std::pow(T/300., b2bodytr_(i,irange2))
                    *std::exp(-c2bodytr_(i,irange2)/T);
          } else if (frml_2bodytr_(i,irange2) == 4) {
            rate2 = a2bodytr_(i,irange2)*b2bodytr_(i,irange2)
                    *( 0.62 + 0.4767*c2bodytr_(i,irange2)*std::sqrt(300./T) );
          } else if (frml_2bodytr_(i,irange2) == 5) {
            rate2 = a2bodytr_(i,irange2)*b2bodytr_(i,irange2)*(
                1 + 0.0967*c2bodytr_(i,irange2)*std::sqrt(300./T)
                + 28.501*c2bodytr_(i,irange2)*c2bodytr_(i,irange2)/T );
          }
        }
        // assign reaction rate
        k2bodytr_(i) = (rate1 + rate2) * 0.5;
      }
    }
  }

  // special rates and grain assisted reactions
  UpdateRatesSpecial(y, E);

  // isothermal case: temperature dependent rates only calculated once
  if (!NON_BAROTROPIC_EOS) {
    flag_T_rates_ = false;
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn ReactionType ChemNetwork::SortReaction(KidaReaction* pr) const
//! \brief sort the type of the reaction, check format

ReactionType ChemNetwork::SortReaction(KidaReaction* pr) const {
  //---------------- 1 - direct cosmic-ray ionization --------------
  if (pr->itype_ == 1) {
    //check format of reaction
    if (pr->reactants_.size() == 2
        && (pr->products_.size() == 2 || pr->products_.size() == 3)
        && (pr->reactants_[0] == "CR" || pr->reactants_[1] == "CR") ) {
    } else {
      std::stringstream msg;
      msg << "### fatal error in chemnetwork sortreaction() [chemnetwork]"
          << std::endl << "wrong format in cr reaction id=" << pr->id_ << std::endl;
      ATHENA_ERROR(msg);
    }
    return ReactionType::cr;

    //---------------- 2 - cosmic-ray induced photo ionization --------
  } else if (pr->itype_ == 2) {
    //check format of reaction
    if (pr->reactants_.size() == 2 && pr->products_.size() == 2
        && (pr->reactants_[0] == "CRP" || pr->reactants_[1] == "CRP") ) {
    } else {
      std::stringstream msg;
      msg << "### FATAL ERROR in ChemNetwork SortReaction() [ChemNetwork]"
          << std::endl << "Wrong format in CRP reaction ID=" << pr->id_ << std::endl;
      ATHENA_ERROR(msg);
    }
    return ReactionType::crp;

    //---------------- 3 - FUV ionization/dissociation ----------------
  } else if (pr->itype_ == 3) {
    //check format of reaction
    if (pr->reactants_.size() == 2 && pr->products_.size() == 2
        && (pr->reactants_[0] == "Photon" || pr->reactants_[1] == "Photon") ) {
    } else {
      std::stringstream msg;
      msg << "### FATAL ERROR in ChemNetwork SortReaction() [ChemNetwork]"
          << std::endl << "Wrong format in FUV reaction ID=" << pr->id_ << std::endl;
      ATHENA_ERROR(msg);
    }
    return ReactionType::photo;

    //---------------- 4-8 - 2body reaction ---------------------------
  } else if (pr->itype_ >= 4 && pr->itype_ <= 8) {
    //check format
    if (std::find(id_2bodytr_.begin(), id_2bodytr_.end(), pr->id_)
        == id_2bodytr_.end()) {
      if (pr->reactants_.size() != 2 ||
          (pr->products_.size() != 1 && pr->products_.size() != 2
           && pr->products_.size() != 3 && pr->products_.size() != 4)) {
        std::stringstream msg;
        msg << "### FATAL ERROR in ChemNetwork SortReaction() [ChemNetwork]"
            << std::endl << "Wrong format in 2body reaction ID=" << pr->id_
            << std::endl;
        ATHENA_ERROR(msg);
      }
      return ReactionType::twobody;
    } else { //2 body reaction with temperature range
      if (pr->reactants_.size() != 2 ||
          (pr->products_.size() != 1 && pr->products_.size() != 2
           && pr->products_.size() != 3)) {
        std::stringstream msg;
        msg << "### FATAL ERROR in ChemNetwork SortReaction() [ChemNetwork]"
            << std::endl << "Wrong format in 2bodytr reaction ID=" << pr->id_
            << std::endl;
        ATHENA_ERROR(msg);
      }
      return ReactionType::twobodytr;
    }

    //-------------------- 9 - grain assisted reaction ----------------
  } else if (pr->itype_ == 9) {
    //check format
    if ( (pr->reactants_.size() != 1 && pr->reactants_.size() != 2)
         || pr->products_.size() != 1 ) {
      std::stringstream msg;
      msg << "### FATAL ERROR in ChemNetwork SortReaction() [ChemNetwork]"
          << std::endl << "Wrong format in implicit grain reaction ID="
          << pr->id_ << std::endl;
      ATHENA_ERROR(msg);
    }
    return ReactionType::grain_implicit;

    //------------------ 10 - special reactions -----------------------
  } else if (pr->itype_ == 10) {
    if (static_cast<int>(pr->reactants_.size()) > n_insr_
        || static_cast<int>(pr->products_.size()) > n_outsr_) {
      std::stringstream msg;
      msg << "### FATAL ERROR in ChemNetwork SortReaction() [ChemNetwork]"
          << std::endl << "Wrong format in special reaction ID=" << pr->id_
          << std::endl;
      ATHENA_ERROR(msg);
    }
    return ReactionType::special;

    //-------11 - grain collision reactions with electron/ion  -----------
  } else if (pr->itype_ == 11) {
    if (static_cast<int>(pr->reactants_.size()) != 2
        || (pr->products_.size() != 1 && pr->products_.size() != 2
            && pr->products_.size() != 3)) {
      std::stringstream msg;
      msg << "### FATAL ERROR in ChemNetwork SortReaction() [ChemNetwork]"
          << std::endl << "Wrong format in grain collsion reaction ID="
          << pr->id_ << std::endl;
      ATHENA_ERROR(msg);
    }
    if (pr->reactants_[0].find("g") != 0 && pr->reactants_[0].find("PAH") != 0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in ChemNetwork SortReaction() [ChemNetwork]"
          << std::endl << "grain collision reactions must be with g or PAH."
          << std::endl;
      ATHENA_ERROR(msg);
    }
    return ReactionType::grain_collision;


    //------------------ type not recogonized -------------------------
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in ChemNetwork SortReaction() [ChemNetwork]"
        << std::endl << "Wrong format in special reaction ID=" << pr->id_
        << std::endl;
    ATHENA_ERROR(msg);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void ChemNetwork::CheckReaction(KidaReaction reaction)
//! \brief check the conservation of atoms and charge for chemical reaction

void ChemNetwork::CheckReaction(KidaReaction reaction) {
  int atom_count_in[KidaSpecies::natom_];
  int atom_count_out[KidaSpecies::natom_];
  int charge_in = 0;
  int charge_out = 0;
  for (int ia=0; ia<KidaSpecies::natom_; ia++) {
    atom_count_in[ia] = 0;
    atom_count_out[ia] = 0;
  }
  for (std::size_t i=0; i<reaction.reactants_.size(); i++) {
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

  for (std::size_t i=0; i<reaction.products_.size(); i++) {
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

//----------------------------------------------------------------------------------------
//! \fn void ChemNetwork::PrintProperties() const
//! \brief print out reactions and rates, for debug

void ChemNetwork::PrintProperties() const {
  // print each species.
  for (int i=0; i<NSPECIES; i++) {
    std::cout << "species: " << i << std::endl;
    std::cout << "name=" << species_[i].name << ", index=" << species_[i].index
              << ", charge=" << species_[i].charge_ << std::endl;
    std::cout << "atom_count_ = ";
    for (int j=0; j < species_[i].natom_; j++) {
      std::cout << species_[i].atom_count_[j] << " ";
    }
    std::cout << std::endl;
  }

  // print each reactions.
  std::cout << "number of reactions: " << nr_ << std::endl;
  for (std::uint64_t i=0; i<reactions_.size(); i++) {
    reactions_[i].Print();
    std::cout << "alpha=" << reactions_[i].alpha_ << ","
              << "beta=" << reactions_[i]. beta_ << ","
              << "gamma=" << reactions_[i].gamma_ << ","
              << "Tmin=" << reactions_[i].Tmin_ << ","
              << "Tmax=" << reactions_[i].Tmax_ << ","
              << "itype=" << reactions_[i].itype_ << ","
              << "forumla=" << reactions_[i].formula_ << std::endl;
  }

  // print reaction coefficients
  // cosmic-ray reactions
  std::cout << "CR reations:" << std::endl;
  for (int i=0; i<n_cr_; i++) {
    std::cout<< species_names[incr_(i)] << " + CR -> "
             << species_names[outcr1_(i)] << " + " << species_names[outcr2_(i)];
    if (outcr3_(i) >= 0) {
      std::cout<< " + " << species_names[outcr3_(i)];
    }
    std::cout << ", kcr_base_=" << kcr_base_(i) << std::endl;
  }

  // cosmic-ray induced photo reactions
  std::cout << "CRP reations:" << std::endl;
  for (int i=0; i<n_crp_; i++) {
    std::cout<< species_names[incrp_(i)] << " + CRP -> "
             << species_names[outcrp1_(i)] << " + " << species_names[outcrp2_(i)] << ", "
             << "kcrp_base_=" << kcrp_base_(i) << std::endl;
  }

  // FUV reactions
  std::cout << "FUV photo- ionization/dissociation:" << std::endl;
  for (int i=0; i<n_ph_; i++) {
    std::cout<< species_names[inph_(i)] << " + Photon -> "
             << species_names[outph1_(i)] << " + " << species_names[outph2_(i)] << ", "
             << "kph_base_=" << kph_base_(i) << ", kph_avfac_=" << kph_avfac_(i)
             << std::endl;
  }
  std::cout << "smap_ph_: " << std::endl;
  for (std::map<std::string,int>::const_iterator it=smap_ph_.begin();
       it!=smap_ph_.end(); it++) {
    std::cout << it->first << " => " << it->second << std::endl;
  }

  // 2body reactions
  std::cout << "2body reactions:" << std::endl;
  for (int i=0; i<n_2body_; i++) {
    for (int jin=0; jin<n_in2body_; jin++) {
      if (in2body_(i, jin) >= 0) {
        std::cout << species_names[in2body_(i, jin)];
        if (jin < n_in2body_-1 && in2body_(i, jin+1) >= 0) {
          std::cout << " + ";
        }
      }
    }
    std::cout << " -> ";
    for (int jout=0; jout<n_out2body_; jout++) {
      if (out2body_(i, jout) >= 0) {
        std::cout << species_names[out2body_(i, jout)];
        if (jout < n_out2body_-1 && out2body_(i, jout+1) >= 0) {
          std::cout << " + ";
        }
      }
    }
    std::cout<< ", " << "alpha=" << a2body_(i) << ", beta=" << b2body_(i)
             << ", gamma=" << c2body_(i) << ", Trange=[" << Tmin_2body_(i) << ","
             << Tmax_2body_(i) << "]" << std::endl;
  }

  // 2body reactions with temperature ranges
  std::cout << "2bodytr reactions:" << std::endl;
  for (int i=0; i<n_2bodytr_; i++) {
    std::cout<< species_names[in2bodytr1_(i)] << " + "
             << species_names[in2bodytr2_(i)]<< " -> " << species_names[out2bodytr1_(i)];
    if (out2bodytr2_(i) >= 0) {
      std::cout<< " + " << species_names[out2bodytr2_(i)];
    }
    if (out2bodytr3_(i) >= 0) {
      std::cout<< " + " << species_names[out2bodytr3_(i)];
    }
    std::cout << "    ,nr_2bodytr_=" << nr_2bodytr_(i) << std::endl;
    for (int j=0; j<nr_2bodytr_(i); j++) {
      std::cout<< "alpha=" << a2bodytr_(i, j) << ", beta="
               << b2bodytr_(i, j) << ", gamma=" << c2bodytr_(i, j)
               << ", Trange=[" << Tmin_2bodytr_(i, j) << "," << Tmax_2bodytr_(i, j)
               << "], formula=" << frml_2bodytr_(i, j) << std::endl;
    }
  }

  // grain assisted reactions
  std::cout << "gr reations:" << std::endl;
  for (int i=0; i<n_gr_; i++) {
    std::cout << species_names[ingr1_(i)];
    if (ingr2_(i) >= 0) {
      std::cout << " + " << species_names[ingr2_(i)];
    }
    std::cout <<" -> " << species_names[outgr_(i)] << std::endl;
  }

  // special reactions
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

  // grain collision reactions
  std::cout << "grain collsion reactions:" << std::endl;
  for (int i=0; i<n_gc_; i++) {
    for (int jin=0; jin<n_ingc_; jin++) {
      if (ingc_(i, jin) >= 0) {
        std::cout << species_names[ingc_(i, jin)];
        if (jin < n_ingc_-1 && ingc_(i, jin+1) >= 0) {
          std::cout << " + ";
        }
      }
    }
    std::cout << " -> ";
    for (int jout=0; jout<n_outgc_; jout++) {
      if (outgc_(i, jout) >= 0) {
        std::cout << species_names[outgc_(i, jout)];
        if (jout < n_outgc_-1 && outgc_(i, jout+1) >= 0) {
          std::cout << " + ";
        }
      }
    }
    std::cout << ", r1_gc_=" << r1_gc_(i) << ", t1_gc_=" << t1_gc_(i) << std::endl;
  }

  for (int i=0; i<id7max_+1; i++) {
    if (id7map_(i) >= 0) {
      std::cout << i << " => " << id7map_(i) << ", ";
      switch (id7type_(i)) {
        case ReactionType::cr: std::cout << "cr"; break;
        case ReactionType::crp: std::cout << "crp"; break;
        case ReactionType::twobody: std::cout << "2body"; break;
        case ReactionType::grain_implicit: std::cout << "gr"; break;
        case ReactionType::special: std::cout << "sr"; break;
        default: std::stringstream msg;
          msg << "### FATAL ERROR in ChemNetwork PrintPropeties() "
              << "[ChemNetwork]: reaction type not recognized for special rate."
              << std::endl;
          ATHENA_ERROR(msg);
          break;
      }
      std::cout << std::endl;
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ChemNetwork::OutputRates(FILE *pf) const
//! \brief output the reactions and rates

void ChemNetwork::OutputRates(FILE *pf) const {
  for (int i=0; i<n_cr_; i++) {
    fprintf(pf, "%4s + CR -> %4s + %4s",
            species_names[incr_(i)].c_str(), species_names[outcr1_(i)].c_str(),
            species_names[outcr2_(i)].c_str());
    if (outcr3_(i) >= 0) {
      fprintf(pf, " + %4s", species_names[outcr3_(i)].c_str());
    }
    fprintf(pf,   ",     kcr = %.2e\n", kcr_(i));
  }
  for (int i=0; i<n_crp_; i++) {
    fprintf(pf, "%4s + CRP -> %4s + %4s,     kcrp = %.2e\n",
            species_names[incrp_(i)].c_str(), species_names[outcrp1_(i)].c_str(),
            species_names[outcrp2_(i)].c_str(), kcrp_(i));
  }
  for (int i=0; i<n_ph_; i++) {
    fprintf(pf, "%4s + Photon -> %4s + %4s,     kph = %.2e\n",
            species_names[inph_(i)].c_str(), species_names[outph1_(i)].c_str(),
            species_names[outph2_(i)].c_str(), kph_(i));
  }
  for (int i=0; i<n_2body_; i++) {
    for (int jin=0; jin<n_in2body_; jin++) {
      if (in2body_(i, jin) >= 0) {
        fprintf(pf, "%4s", species_names[in2body_(i, jin)].c_str());
        if (jin < n_in2body_-1 && in2body_(i, jin+1) >= 0) {
          fprintf(pf, " + ");
        }
      }
    }
    fprintf(pf, " -> ");
    for (int jout=0; jout<n_out2body_; jout++) {
      if (out2body_(i, jout) >= 0) {
        fprintf(pf, "%4s", species_names[out2body_(i, jout)].c_str());
        if (jout < n_out2body_-1 && out2body_(i, jout+1) >= 0) {
          fprintf(pf, " + ");
        }
      }
    }
    fprintf(pf,   ",     nk2body = %.2e\n", k2body_(i)*nH_);
  }
  for (int i=0; i<n_2bodytr_; i++) {
    fprintf(pf, "%4s + %4s -> %4s",
            species_names[in2bodytr1_(i)].c_str(),
            species_names[in2bodytr2_(i)].c_str(),
            species_names[out2bodytr1_(i)].c_str());
    if (out2bodytr2_(i) >= 0) {
      fprintf(pf, " + %4s", species_names[out2bodytr2_(i)].c_str());
    }
    if (out2bodytr3_(i) >= 0) {
      fprintf(pf, " + %4s", species_names[out2bodytr3_(i)].c_str());
    }
    fprintf(pf,   ",     nk2bodytr = %.2e\n", k2bodytr_(i)*nH_);
  }
  for (int i=0; i<n_gr_; i++) {
    if (ingr2_(i) >= 0) {
      fprintf(pf, "%4s + %4s -> %4s,       kgr = %.2e\n",
              species_names[ingr1_(i)].c_str(), species_names[ingr2_(i)].c_str(),
              species_names[outgr_(i)].c_str(), kgr_(i));
    } else {
      fprintf(pf, "%4s -> %4s,       kgr = %.2e\n",
              species_names[ingr1_(i)].c_str(),
              species_names[outgr_(i)].c_str(), kgr_(i));
    }
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
    fprintf(pf, ",       ksr = %.2e\n", ksr_(i));
  }
  for (int i=0; i<n_gc_; i++) {
    for (int jin=0; jin<n_ingc_; jin++) {
      if (ingc_(i, jin) >= 0) {
        fprintf(pf, "%4s", species_names[ingc_(i, jin)].c_str());
        if (jin < n_ingc_-1 && ingc_(i, jin+1) >= 0) {
          fprintf(pf, " + ");
        }
      }
    }
    fprintf(pf, " -> ");
    for (int jout=0; jout<n_outgc_; jout++) {
      if (outgc_(i, jout) >= 0) {
        fprintf(pf, "%4s", species_names[outgc_(i, jout)].c_str());
        if (jout < n_outgc_-1 && outgc_(i, jout+1) >= 0) {
          fprintf(pf, " + ");
        }
      }
    }
    fprintf(pf, ",       nkgc = %.2e\n", kgc_(i)*nH_);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ChemNetwork::OutputJacobian(FILE *pf, const AthenaArray<Real> &jac) const
//! \brief output jacobian coefficients, for debug

void ChemNetwork::OutputJacobian(FILE *pf, const AthenaArray<Real> &jac) const {
  const int dim = jac.GetDim1();
  for (int i=0; i<dim; i++) {
    for (int j=0; j<dim; j++) {
      fprintf(pf, "%12.4e ", jac(i,j));
    }
    fprintf(pf, "\n");
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ChemNetwork::Jacobian_isothermal_numerical(const Real t,
//!    const Real *y, const Real *ydot, AthenaArray<Real> &jac)
//! \brief calculate Jacobian with numerical differentiation
void ChemNetwork::Jacobian_isothermal_numerical(
    const Real t, const Real *y, const Real *ydot, AthenaArray<Real> &jac) {
  const Real dy = 1e-3;
  Real *y1 = new Real[NSPECIES];
  Real *y2 = new Real[NSPECIES];
  Real *ydot1 = new Real[NSPECIES];
  Real *ydot2 = new Real[NSPECIES];
  Real *eps = new Real[NSPECIES];
  for (int i=0; i<NSPECIES; i++) {
    eps[i] = dy;
  }
  for (int i=0; i<NSPECIES; i++) {
    for (int j=0; j<NSPECIES; j++) {
      for (int k=0; k<NSPECIES; k++) {
        y1[k] = y[k];
        y2[k] = y[k];
      }
      y1[j] -= eps[j];
      y2[j] += eps[j];
      RHS(t, y1, 0., ydot1);
      RHS(t, y2, 0., ydot2);
      jac(i, j) = (ydot2[i] - ydot1[i])/(2.*eps[j]);
    }
  }
  delete[] y1;
  delete[] y2;
  delete[] ydot1;
  delete[] ydot2;
  delete[] eps;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ChemNetwork::SetGrad_v(const int k, const int j, const int i)
//! \brief set gradients of v and nH for CO cooling

void ChemNetwork::SetGrad_v(const int k, const int j, const int i) {
  AthenaArray<Real> &w = pmy_mb_->phydro->w;
  Real dvdx, dvdy, dvdz, dvdr_avg, di1, di2;
  Real dx1, dx2, dy1, dy2, dz1, dz2;
  // Real dndx, dndy, dndz, gradn;

  // velocity gradient, same as LVG approximation in RADMC-3D when calculating
  // CO line emission.
  // vx
  di1 = w(IVX, k, j, i+1) - w(IVX, k, j, i);
  dx1 = ( pmy_mb_->pcoord->dx1f(i+1)+pmy_mb_->pcoord->dx1f(i) )/2.;
  di2 = w(IVX, k, j, i) - w(IVX, k, j, i-1);
  dx2 = ( pmy_mb_->pcoord->dx1f(i)+pmy_mb_->pcoord->dx1f(i-1) )/2.;
  dvdx = (di1/dx1 + di2/dx2)/2.;
  // vy
  di1 = w(IVY, k, j+1, i) - w(IVY, k, j, i);
  dy1 = ( pmy_mb_->pcoord->dx2f(j+1)+pmy_mb_->pcoord->dx2f(j) )/2.;
  di2 = w(IVY, k, j, i) - w(IVY, k, j-1, i);
  dy2 = ( pmy_mb_->pcoord->dx2f(j)+pmy_mb_->pcoord->dx2f(j-1) )/2.;
  dvdy = (di1/dy1 + di2/dy2)/2.;
  // vz
  di1 = w(IVZ, k+1, j, i) - w(IVZ, k, j, i);
  dz1 = ( pmy_mb_->pcoord->dx3f(k+1)+pmy_mb_->pcoord->dx3f(k) )/2.;
  di2 = w(IVZ, k, j, i) - w(IVZ, k-1, j, i);
  dz2 = ( pmy_mb_->pcoord->dx3f(k)+pmy_mb_->pcoord->dx3f(k-1) )/2.;
  dvdz = (di1/dz1 + di2/dz2)/2.;
  dvdr_avg = ( std::abs(dvdx) + std::abs(dvdy) + std::abs(dvdz) ) / 3.;
  // asign gradv_, in cgs.
  gradv_ = dvdr_avg * pmy_mb_->pmy_mesh->punit->code_velocity_cgs
           / pmy_mb_->pmy_mesh->punit->code_length_cgs;
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void ChemNetwork::UpdateRatesSpecial(const Real *y, const Real E)
//! \brief update reaction rates for special reactions with formula = 7
//!
//! default: no special rates
void __attribute__((weak)) ChemNetwork::UpdateRatesSpecial(const Real *y, const Real E) {
  // do nothing
  return;
}

//----------------------------------------------------------------------------------------
//! \fn ChemNetwork::UpdateJacobianSpecial( const Real *y, const Real E,
//!       AthenaArray<Real> &jac)
//! \brief update jacobian coefficients for special reactions with formula = 7
//!
//! default: no special rates, and therefore no special terms for Jacobian
void __attribute__((weak)) ChemNetwork::UpdateJacobianSpecial(
    const Real *y, const Real E, AthenaArray<Real> &jac) {
  // do nothing
  return;
}
