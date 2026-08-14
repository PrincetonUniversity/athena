//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file unit_test.cpp
//! \brief Problem file to demonstrate unit module
//!

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdio>    // fopen(), fprintf(), freopen()
#include <cstring>   // strcmp()
#include <sstream>
#include <stdexcept>
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Get units pointer
  Units *punit = Mesh::punit;
  // Print units info
  punit->PrintBasisUnits();
  punit->PrintCodeUnits();
  punit->PrintConstantsInCodeUnits();

  std::string lunit = std::get<1>(punit->basis_length);
  std::string tunit = std::get<1>(punit->basis_time);
  std::string munit = std::get<1>(punit->basis_mass);
  std::string vunit = std::get<1>(punit->basis_velocity);
  std::string nunit = std::get<1>(punit->basis_ndensity);

  // Convert time and mesh values from input file
  Real tlim0  = pin->GetReal("time", "tlim");
  Real dt0    = pin->GetReal("output1", "dt");
  Real x1max0 = pin->GetReal("mesh", "x1max");
  punit->ConvertInputFile(pin);
  Real tlim1  = pin->GetReal("time", "tlim");
  Real dt1    = pin->GetReal("output1", "dt");
  Real x1max1 = pin->GetReal("mesh", "x1max");

  std::cout << "=== Values from input file converted ===" << std::endl;
  std::cout << "From the <time> block:" << std::endl;
  std::cout << "  tlim before conversion  = " << tlim0 << " " << tunit << std::endl;
  std::cout << "  tlim after conversion   = " << tlim1 << " code time" << std::endl;
  std::cout << std::endl;
  std::cout << "From the <output1> block:" << std::endl;
  std::cout << "  dt before conversion    = " << dt0 << " " << tunit << std::endl;
  std::cout << "  dt after conversion     = " << dt1 << " code time" << std::endl;
  std::cout << std::endl;
  std::cout << "From the <mesh> block:" << std::endl;
  std::cout << "  x1max before conversion = " << x1max0 << " " << lunit << std::endl;
  std::cout << "  x1max after conversion  = " << x1max1 << " code length" << std::endl;
  std::cout << std::endl;


  // Pressure
  Real pres = pin->GetReal("problem", "pamb");
  Real pres_code = pres/punit->code_pressure_cgs;
  std::cout << "Pressure value in cgs units  = " << pres << " dyne" << std::endl;
  std::cout << "Pressure value in code units = " << pres_code
            << " code pressure" << std::endl;
  std::cout << std::endl;

  // Energy
  Real en = pin->GetReal("problem", "energy");
  Real en_code = en/punit->code_energy_cgs;
  std::cout << "Energy value in cgs units  = " << en << " erg" << std::endl;
  std::cout << "Energy value in code units = " << en_code << " code energy" << std::endl;
  std::cout << std::endl;

  // Mass
  Real mass = pin->GetReal("problem", "mass_val");
  Real mass_code = mass/std::get<0>(punit->basis_mass);
  std::cout << "Mass value in basis units = " << mass << " " << munit << std::endl;
  std::cout << "Mass value in code units  = " << mass_code << " code mass" << std::endl;
  std::cout << std::endl;

  // Length
  Real rad = pin->GetReal("problem", "radius");
  Real rad_code = rad/std::get<0>(punit->basis_length);
  std::cout << "Length value in basis units = " << rad << " " << lunit << std::endl;
  std::cout << "Length value in code units  = " << rad_code
            << " code length" << std::endl;
  std::cout << std::endl;

  // Velocity
  Real vel = pin->GetReal("problem", "vel");
  Real vel_code = vel/std::get<0>(punit->basis_velocity);
  std::cout << "Velocity value in basis units = " << vel << " " << vunit << std::endl;
  std::cout << "Velocity value in code units  = " << vel_code
            << " code velocity" << std::endl;
  std::cout << std::endl;

  // ndensity
  Real ndin = pin->GetReal("problem", "ndin");
  Real ndin_code = ndin/std::get<0>(punit->basis_ndensity);
  std::cout << "ndensity value in basis units = " << ndin << " " << nunit << std::endl;
  std::cout << "ndensity value in code units  = " << ndin_code
            << " code ndensity" << std::endl;
  std::cout << std::endl;

  std::cout << "Gravitational constant in code units = "
            << punit->Gconst_code << std::endl;
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Get units pointer
  Units *punit = MeshBlock::pmy_mesh->punit;

  Real pres = pin->GetReal("problem", "pamb")/punit->code_pressure_cgs;
  Real da   = pin->GetReal("problem", "ndamb")/std::get<0>(punit->basis_ndensity);
  Real din  = pin->GetReal("problem", "ndin")/std::get<0>(punit->basis_ndensity);
  Real rin  = pin->GetReal("problem", "radius")/std::get<0>(punit->basis_length);
  Real vel  = pin->GetReal("problem", "vel")/std::get<0>(punit->basis_velocity);

  Real gamma = peos->GetGamma();
  Real gm1 = gamma - 1.0;

  // setup uniform ambient medium with spherical over-pressured region
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);
        Real rad = std::sqrt(SQR(x) + SQR(y) + SQR(z));
        Real den = da;
        if (rad < rin) {
          den = din;
        }
        phydro->u(IDN,k,j,i) = den;
        phydro->u(IM1,k,j,i) = vel;
        phydro->u(IM2,k,j,i) = vel/2;
        phydro->u(IM3,k,j,i) = vel/3;
        phydro->u(IEN,k,j,i) = pres/gm1;
      }
    }
  }
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//! \brief Print yt units override
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  Units *punit = Mesh::punit;
  // Prints a Python dictionary to use units_override
  //  when loading hdf5 data into yt.
  punit->PrintytUnitsOverride();
}
