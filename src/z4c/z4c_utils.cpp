//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file z4c_utils.cpp
//  \brief support utilities for z4c such as debug etc

#include <cassert> // assert
#include <iostream>

// Athena++ headers
#include "z4c.hpp"

//----------------------------------------------------------------------------------------
// \!fn void Z4c::is_finite_adm()
// \brief Aggregate 'is_finite' conditional
bool Z4c::is_finite_adm() {
  bool finite = true;
  finite &= adm.psi4.is_finite();
  finite &= adm.g_dd.is_finite();
  finite &= adm.K_dd.is_finite();
  return finite;
}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::is_finite_con()
// \brief Aggregate 'is_finite' conditional
bool Z4c::is_finite_con() {
  bool finite = true;
  finite &= con.C.is_finite();
  finite &= con.H.is_finite();
  finite &= con.M.is_finite();
  finite &= con.Z.is_finite();
  finite &= con.M_d.is_finite();
  return finite;
}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::is_finite_mat()
// \brief Aggregate 'is_finite' conditional
bool Z4c::is_finite_mat() {
  bool finite = true;
  finite &= mat.rho.is_finite();
  finite &= mat.S_d.is_finite();
  finite &= mat.S_dd.is_finite();
  return finite;
}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::is_finite_Z4c()
// \brief Aggregate 'is_finite' conditional
bool Z4c::is_finite_z4c() {
  bool finite = true;
  finite &= z4c.chi.is_finite();
  finite &= z4c.Khat.is_finite();
  finite &= z4c.Theta.is_finite();
  finite &= z4c.alpha.is_finite();
  finite &= z4c.Gam_u.is_finite();
  finite &= z4c.beta_u.is_finite();
  finite &= z4c.g_dd.is_finite();
  finite &= z4c.A_dd.is_finite();
  return finite;
}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::assert_is_finite_adm()
// \brief Aggregate 'is_finite' conditional as assertion
void Z4c::assert_is_finite_adm() {
  if (not adm.psi4.is_finite()) {
    coutBoldRed("adm.psi4 not finite, terminating...\n");
    std::exit(EXIT_FAILURE);
  }
  if (not adm.g_dd.is_finite()) {
    coutBoldRed("adm.g_dd not finite, terminating...\n");
    std::exit(EXIT_FAILURE);
  }
  if (not adm.K_dd.is_finite()) {
    coutBoldRed("adm.K_dd not finite, terminating...\n");
    std::exit(EXIT_FAILURE);
  }
}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::assert_is_finite_con()
// \brief Aggregate 'is_finite' conditional as assertion
void Z4c::assert_is_finite_con() {
  if (not con.C.is_finite()) {
    coutBoldRed("con.C not finite, terminating...\n");
    std::exit(EXIT_FAILURE);
  }
  if (not con.H.is_finite()) {
    coutBoldRed("con.H not finite, terminating...\n");
    std::exit(EXIT_FAILURE);
  }
  if (not con.M.is_finite()) {
    coutBoldRed("con.M not finite, terminating...\n");
    std::exit(EXIT_FAILURE);
  }
  if (not con.Z.is_finite()) {
    coutBoldRed("con.Z not finite, terminating...\n");
    std::exit(EXIT_FAILURE);
  }
  if (not con.M_d.is_finite()) {
    coutBoldRed("con.M_d not finite, terminating...\n");
    std::exit(EXIT_FAILURE);
  }
}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::assert_is_finite_mat()
// \brief Aggregate 'is_finite' conditional as assertion
void Z4c::assert_is_finite_mat() {
  if (not mat.rho.is_finite()) {
    coutBoldRed("mat.rho not finite, terminating...\n");
    std::exit(EXIT_FAILURE);
  }
  if (not mat.S_d.is_finite()) {
    coutBoldRed("mat.S_d not finite, terminating...\n");
    std::exit(EXIT_FAILURE);
  }
  if (not mat.S_dd.is_finite()) {
    coutBoldRed("mat.S_dd not finite, terminating...\n");
    std::exit(EXIT_FAILURE);
  }
}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::assert_is_finite_z4c()
// \brief Aggregate 'is_finite' conditional as assertion
void Z4c::assert_is_finite_z4c() {
  if (not z4c.chi.is_finite()) {
    coutBoldRed("z4c.chi not finite, terminating...\n");
    std::exit(EXIT_FAILURE);
  }
  if (not z4c.Khat.is_finite()) {
    coutBoldRed("z4c.Khat not finite, terminating...\n");
    std::exit(EXIT_FAILURE);
  }
  if (not z4c.Theta.is_finite()) {
    coutBoldRed("z4c.Theta not finite, terminating...\n");
    std::exit(EXIT_FAILURE);
  }
  if (not z4c.alpha.is_finite()) {
    coutBoldRed("z4c.alpha not finite, terminating...\n");
    std::exit(EXIT_FAILURE);
  }
  if (not z4c.Gam_u.is_finite()) {
    coutBoldRed("z4c.Gam_u not finite, terminating...\n");
    std::exit(EXIT_FAILURE);
  }
  if (not z4c.beta_u.is_finite()) {
    coutBoldRed("z4c.beta_u not finite, terminating...\n");
    std::exit(EXIT_FAILURE);
  }
  if (not z4c.g_dd.is_finite()) {
    coutBoldRed("z4c.g_dd not finite, terminating...\n");
    std::exit(EXIT_FAILURE);
  }
  if (not z4c.A_dd.is_finite()) {
    coutBoldRed("z4c.A_dd not finite, terminating...\n");
    std::exit(EXIT_FAILURE);
  }
}
