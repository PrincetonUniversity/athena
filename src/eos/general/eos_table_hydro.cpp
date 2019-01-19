//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file eos_table_hydro.cpp
//  \brief implements functions in class EquationOfState for an EOS lookup table
//======================================================================================

// C headers

// C++ headers
#include <cfloat>  // FLT_MIN
#include <cmath>   // sqrt()
#include <fstream>
#include <iostream> // ifstream
#include <sstream>
#include <stdexcept> // std::invalid_argument
#include <string>

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../field/field.hpp"
#include "../../parameter_input.hpp"
#include "../../utils/interp_table.hpp"

// this class header
#include "../eos.hpp"

//#define EOSDEBUG0
//#define EOSDEBUG1

void EosTestLoop(EquationOfState *peos);

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::SimplePres(Real rho, Real egas)
//  \brief Return interpolated gas pressure
Real EquationOfState::SimplePres(Real rho, Real egas) {
  return ptable->GetEosData(0, egas, rho) * egas;
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::SimpleEgas(Real rho, Real pres)
//  \brief Return interpolated internal energy density
Real EquationOfState::SimpleEgas(Real rho, Real pres) {
  return ptable->GetEosData(1, pres, rho) * pres;
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::SimpleAsq(Real rho, Real pres)
//  \brief Return interpolated adiabatic sound speed squared
Real EquationOfState::SimpleAsq(Real rho, Real pres) {
  return ptable->GetEosData(2, pres, rho) * pres / rho;
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::RiemannAsq(Real rho, Real hint)
//  \brief Return interpolated adiabatic sound speed squared for use in
//         Riemann solver.
Real EquationOfState::RiemannAsq(Real rho, Real hint) {
  return std::pow(static_cast<Real>(10.0),
                  ptable->table.interpolate(3, std::log10(hint*ptable->hUnit),
                                            std::log10(rho*ptable->rhoUnit))) * hint;
}
//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::PrepEOS(ParameterInput *pin)
//  \brief Read data and initialize interpolated table.
void EquationOfState::PrepEOS(ParameterInput *pin) {
  // Debugging/testing code
#ifdef EOSDEBUG1
  std::cout << "Eos table: " << ptable->nVar << ", " << ptable->nRho << ", "
            << ptable->nEgas << std::endl;
  std::cout << "logRhoMin, logRhoMax: " << ptable->logRhoMin << ", " << ptable->logRhoMax
            << std::endl;
  std::cout << "logEgasMin, logEgasMax: " << ptable->logEgasMin << ", "
            << ptable->logEgasMax << std::endl;
  std::cout << "Ratios: ";
  for (int i=0; i<ptable->nVar; i++) std::cout << ptable->EosRatios(i) << ", ";
  std::cout << std::endl;
  for (int i=0; i<ptable->nVar; i++) {
    std::cout << "var = " << i << std::endl;
    for (int j=0; j<ptable->nRho; j++) {
      for (int k=0; k<ptable->nEgas; k++) {
        std::cout << std::pow((Real) 10, ptable->table.data(i,j,k)) << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
#endif

  // Debugging/testing code
#ifdef EOSDEBUG0
  std::cout << "prepEOS: " << ptable->nVar << ", " << ptable->nEgas
            <<  ", " << ptable->nRho << std::endl;
  std::cout << "eUnit, rhoUnit, hUnit: " << ptable->eUnit << ", " << ptable->rhoUnit
            << ", " << ptable->hUnit << '\n';
  std::cout << "p(1e-7,1e-7)= " << ptable->GetEosData(0, 1e-7, 1e-7) << std::endl;
  EosTestLoop(this);
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::CleanEOS()
//  \brief Clear memory/objects used to store EOS data
void EquationOfState::CleanEOS() {
}

//----------------------------------------------------------------------------------------
//! \fn void EosTestRhoEgas(Real rho, Real egas, AthenaArray<Real> &data)
//  \brief Debugging/testing function
void EosTestRhoEgas(EquationOfState *peos, Real rho, Real egas, AthenaArray<Real> &data) {
  Real idn = 1./rho;
  Real ien = 1./egas;
  data(0) = peos->ptable->GetEosData(0, egas, rho) * egas; // pressure
  data(1) = (data(0) + egas) * idn; // specific enthalpy
  data(2) = peos->ptable->GetEosData(2, data(0), rho) * data(0) * idn; // Asq
  // PrimToCons error
  data(3) = (egas-peos->ptable->GetEosData(1, data(0), rho) * data(0)) * ien;
  data(4) = 1.0 - peos->RiemannAsq(rho, data(1)) / data(2); // Asq error
}

//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::EosTestLoop()
//  \brief Debugging/testing function
void EosTestLoop(EquationOfState *peos) {
  Real rho = 0;
  Real egas;
  EosTable *ptable = peos->ptable;
  AthenaArray<Real> data;
  data.NewAthenaArray(5);
  std::cout << "logRhoMin, logRhoMax: " << ptable->logRhoMin << ", " << ptable->logRhoMax
            << std::endl;
  std::cout << "logEgasMin, logEgasMax: " << ptable->logEgasMin << ", "
            << ptable->logEgasMax << std::endl;
  std::cout << "EosRatios = ";
  for (int i=0; i < ptable->nVar; ++i) std::cout << ptable->EosRatios(i) << ", ";
  std::cout << std::endl;
  std::cout << "Input rho (g/cc):";
  std::cin >> rho;
  std::cout << "Input egas (erg/cc):";
  std::cin >> egas;
  while (rho >= 0) {
    EosTestRhoEgas(peos, rho, egas, data);
    std::cout << rho << ", " << egas << std::endl;
    std::cout << "P(d, e)    , h(d, e)    , ASq(d, P)  ,PErr    , ASqErr\n";
    for (int i=0; i<5; i++) std::cout << data(i) << ", ";
    std::cout << std::endl;
    std::cout << "P, e, Asq, Asq: " << peos->SimplePres(rho, egas) << ", "
              << peos->SimpleEgas(rho, data(0)) << ", "
              << peos->SimpleAsq(rho, data(0)) << ", "
              << peos->RiemannAsq(rho, data(1)) << "\n\n";
    std::cout << "Input rho (g/cc):";
    std::cin >> rho;
    std::cout << "Input egas (erg/cc):";
    std::cin >> egas;
  }
  data.DeleteAthenaArray();
}
