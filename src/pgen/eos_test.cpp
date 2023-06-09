//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file eos_test.cpp
//! \brief Problem generator to test/verify equation of state.
//!
//! Problem generator designed to test, debug, and verify EOS tables and functions
//
//========================================================================================

// C headers

// C++ headers
#include <cmath>      // sqrt()
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../utils/interp_table.hpp"

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Problem Generator for eos tests
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  EosTable *ptable = pmy_mesh->peos_table;
  bool print_table = pin->GetOrAddBoolean("problem","print_table",false);
  bool exp_table = pin->GetOrAddBoolean("problem","exponentiate_table",false);
  bool eos_loop = pin->GetOrAddBoolean("problem","eos_loop",false);

  // Print EOS info
  std::cout << "Equation of state (EOS) diagnostics:" << '\n';
  if (GENERAL_EOS) {
    std::cout << "General EOS enabled with " << EQUATION_OF_STATE << '\n';
    if (EOS_TABLE_ENABLED) {
      // Print table info
      std::cout << "Using table '" << pin->GetString("hydro", "eos_file_name") << "'.\n";
      std::cout << "Shape (nVar, nEgas, nRho): " << ptable->nVar << ", " << ptable->nEgas
                <<  ", " << ptable->nRho << '\n';
      std::cout << "logEgasMin, logEgasMax: " << ptable->logEgasMin << ", "
                << ptable->logEgasMax << '\n';
      std::cout << "logRhoMin, logRhoMax: " << ptable->logRhoMin << ", "
                << ptable->logRhoMax << '\n';
      std::cout << "eUnit, rhoUnit, hUnit: " << ptable->eUnit << ", " << ptable->rhoUnit
                << ", " << ptable->hUnit << '\n';
      std::cout << "EosRatios: ";
      for (int i=0; i < ptable->nVar; ++i) std::cout << ptable->EosRatios(i) << ", ";
      std::cout << '\n';
    }
  } else if (NON_BAROTROPIC_EOS) {
    std::cout << "Adiabatic/ideal EOS enabled." << '\n';
  } else {
    std::cout << "Isothermal EOS enabled." << '\n';
  }
  std::cout << std::endl;

  if (print_table) {
    if (EOS_TABLE_ENABLED) {
      // Print EOS table
      for (int i=0; i<ptable->nVar; i++) {
        std::cout << "var = " << i << std::endl;
        if (exp_table) {
          for (int j=0; j<ptable->nEgas; j++) {
            for (int k=0; k<ptable->nRho; k++) {
              std::cout << std::pow((Real) 10, ptable->table.data(i,j,k)) << " ";
            }
            std::cout << '\n';
          }
          std::cout << '\n';
        } else {
          for (int j=0; j<ptable->nEgas; j++) {
            for (int k=0; k<ptable->nRho; k++) {
              std::cout << ptable->table.data(i,j,k) << " ";
            }
            std::cout << '\n';
          }
          std::cout << '\n';
        }
      }
    } else {
      std::cout << "Warning: print table specified, but no table to print." << '\n';
    }
    std::cout << std::endl;
  }

  if (eos_loop) {
    if (!NON_BAROTROPIC_EOS) {
      std::stringstream msg;
      msg << "### FATAL ERROR in Problem Generator" << std::endl
          << "Isothermal EOS incompatible with eos_loop."<< std::endl;
      ATHENA_ERROR(msg);
    }
    Real rho=0.0, egas=0.0;
    AthenaArray<Real> u, w, zeros;
    u.NewAthenaArray(NHYDRO,1,1,1);
    w.NewAthenaArray(NHYDRO,1,1,1);
    zeros.NewAthenaArray(NHYDRO,1,1,1);
    Real w2[(NHYDRO)];
    for (int i=0; i<NHYDRO; i++) {
      u(i)=0.0;
      zeros(i)=0.0;
      w2[i]=0.0;
    }
    std::cout << "Input fluid parameters and retrieve EOS parameters." << '\n'
              << "Non-positive inputs will exit loop." << '\n';
    std::cout << "Input density (mass/volume): ";
    std::cin >> rho;
    std::cout << "Input internal energy (energy/volume): ";
    std::cin >> egas;
    while (rho > 0 && std::isfinite(rho) && egas >0 && std::isfinite(egas)) {
      Real p, asq, perr;  //, h;
      FaceField f;
      std::cout << "Density, internal energy: " << rho << ", " << egas << '\n';
      u(IDN,0,0,0) = rho;
      u(IEN,0,0,0) = egas;
      peos->ConservedToPrimitive(u, zeros, f, w, zeros, pcoord, 0, 0, 0, 0, 0, 0);
      peos->PrimitiveToConserved(w, zeros, u, pcoord, 0, 0, 0, 0, 0, 0);
      p = w(IPR,0,0,0);
      w2[IPR]=p;
      w2[IDN]=rho;
      asq = SQR(peos->SoundSpeed(w2));
      //h = p + egas;
      perr = 1.0 - u(IEN)/egas;
      std::cout << "P(d, e)    , ASq(d, P)  , PErr\n";
      std::cout << p << ", " << asq  << ", " << perr << '\n' << std::endl;
      std::cout << "Input density (mass/volume): ";
      std::cin >> rho;
      std::cout << "Input internal energy (energy/volume): ";
      std::cin >> egas;
    }
    std::cout << std::endl;
  }

  // Set all fluid parameters to zero
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        phydro->u(IDN,k,j,i) = 0.0;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS) phydro->u(IEN,k,j,i) = 0.0;
      }
    }
  }
  return;
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
}
