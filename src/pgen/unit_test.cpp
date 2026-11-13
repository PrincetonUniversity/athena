//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file unit_test.cpp
//! \brief Problem generator demonstrating the units module: a uniform medium with a
//!        dense sphere, set up from physical values in the <problem> block.

// C headers

// C++ headers
#include <cmath>      // sqrt()
#include <iostream>   // cout, endl
#include <string>
#include <tuple>      // get()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../units/units.hpp"

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Print the unit system and example conversions of <problem> values
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (Globals::my_rank != 0) return;

  punit->PrintBasisUnits();
  punit->PrintCodeUnits();
  punit->PrintConstantsInCodeUnits();

  // Values in c.g.s. are converted with code_X_cgs
  Real pres = pin->GetReal("problem", "pamb");
  Real en = pin->GetReal("problem", "energy");
  std::cout << "pressure = " << pres << " dyne/cm^2 = "
            << pres/punit->code_pressure_cgs << " code pressure" << std::endl;
  std::cout << "energy = " << en << " erg = "
            << en/punit->code_energy_cgs << " code energy" << std::endl;

  // Values in the basis units are converted with std::get<0>(basis_X)
  Real mass = pin->GetReal("problem", "mass_val");
  Real rad = pin->GetReal("problem", "radius");
  Real vel = pin->GetReal("problem", "vel");
  Real ndin = pin->GetReal("problem", "ndin");
  std::cout << "mass = " << mass << " " << std::get<1>(punit->basis_mass) << " = "
            << mass/std::get<0>(punit->basis_mass) << " code mass" << std::endl;
  std::cout << "length = " << rad << " " << std::get<1>(punit->basis_length) << " = "
            << rad/std::get<0>(punit->basis_length) << " code length" << std::endl;
  std::cout << "velocity = " << vel << " " << std::get<1>(punit->basis_velocity)
            << " = " << vel/std::get<0>(punit->basis_velocity) << " code velocity"
            << std::endl;
  std::cout << "ndensity = " << ndin << " " << std::get<1>(punit->basis_ndensity)
            << " = " << ndin/std::get<0>(punit->basis_ndensity) << " code density"
            << std::endl;
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Uniform medium with a dense sphere of radius rin, moving with velocity
//!        (vel, vel/2, vel/3)
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Units *punit = pmy_mesh->punit;

  Real pres = pin->GetReal("problem", "pamb")/punit->code_pressure_cgs;
  Real da   = pin->GetReal("problem", "ndamb")/std::get<0>(punit->basis_ndensity);
  Real din  = pin->GetReal("problem", "ndin")/std::get<0>(punit->basis_ndensity);
  Real rin  = pin->GetReal("problem", "radius")/std::get<0>(punit->basis_length);
  Real vel  = pin->GetReal("problem", "vel")/std::get<0>(punit->basis_velocity);
  Real v1 = vel, v2 = vel/2.0, v3 = vel/3.0;

  Real gm1 = peos->GetGamma() - 1.0;

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);
        Real rad = std::sqrt(SQR(x) + SQR(y) + SQR(z));
        Real den = (rad < rin) ? din : da;
        phydro->u(IDN,k,j,i) = den;
        phydro->u(IM1,k,j,i) = den*v1;
        phydro->u(IM2,k,j,i) = den*v2;
        phydro->u(IM3,k,j,i) = den*v3;
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = pres/gm1 + 0.5*den*(SQR(v1) + SQR(v2) + SQR(v3));
        }
      }
    }
  }
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//! \brief Print the yt units_override dictionary for loading the HDF5 outputs
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  if (Globals::my_rank == 0) punit->PrintYtUnitsOverride();
}
