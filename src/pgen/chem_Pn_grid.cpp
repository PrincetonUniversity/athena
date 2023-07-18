//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file chem_Pn_grid.cpp
//! \brief problem generator, for producing a grid of P-n curves.
//!
//! x - nH, hydrogen nucleus number density in cm-3
//! y - chi, FUV radiation intensity relative to solar neighborhood
//! z - CRIR, primary CRIR per hydrogen nucleus
//! Gas and dust metallicities needs to be specified by the input file.
//======================================================================================

// C headers

// C++ headers
#include <algorithm>  // std::find()
#include <cstdio>     // snprintf
#include <fstream>    // ifstream
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // std::runtime_error()
#include <string>     // c_str()
#include <vector>     // vector container

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../chem_rad/chem_rad.hpp"
#include "../chem_rad/integrators/rad_integrators.hpp"
#include "../chemistry/utils/thermo.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"


//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief initialize problem of a grid of radiation and density conditions
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // dimensions of meshblock
  const int Nx = ie - is + 1;
  const int Ny = je - js + 1;
  const int Nz = ke - ks + 1;
  // read density, FUV radiation, and CRIR grids
  const Real nH_min = pin->GetReal("problem", "nH_min"); // minimum density
  const Real nH_max = pin->GetReal("problem", "nH_max"); // maximum density
  // initial abundance and sound speed (which sets the initial temperature)
  const Real r_init = pin->GetReal("problem", "r_init");
  const Real iso_cs = pin->GetReal("hydro", "iso_sound_speed");
  const Real gm1  = peos->GetGamma() - 1.0;
  Real pres = 0.;
  // factors between the logspace grids cells
  const Real fac_nH = (std::log10(nH_max) - std::log10(nH_min) ) / (Nx-1.);
  // arrays for storing chi and crir
  AthenaArray<Real> chi;
  AthenaArray<Real> cr;
  chi.NewAthenaArray(Nz, Ny);
  cr.NewAthenaArray(Nz, Ny);

  // read chi and crir from text files
  std::string dir_input = pin->GetString("problem", "dir_input");
  const Real Zdg = pin->GetReal("chemistry", "Zdg");
  char dir_z[20];
  std::snprintf(dir_z, sizeof(dir_z), "Z%.1f/", Zdg);
  // std::string dir_z = std::to_string(dir_z_buf);
  std::string fn_chi = dir_input + dir_z + "chi.txt";
  std::string fn_cr = dir_input + dir_z + "cr.txt";
  std::cout << "filename chi: " << fn_chi << std::endl;
  std::cout << "filename cr: " << fn_cr << std::endl;
  // open files
  std::ifstream infile_chi(fn_chi);
  if (!infile_chi.is_open()) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function MeshBlock::ProblemGenerator"
      << std::endl << "Cannot open file: " << fn_chi << std::endl;
    ATHENA_ERROR(msg);
  }
  std::ifstream infile_cr(fn_cr);
  if (!infile_cr.is_open()) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function MeshBlock::ProblemGenerator"
      << std::endl << "Cannot open file: " << fn_cr << std::endl;
    ATHENA_ERROR(msg);
  }
  // read data into array
  for (int k=0; k<Nz; ++k) {
    for (int j=0; j<Ny; ++j) {
      infile_chi >> chi(k, j);
      infile_cr >> cr(k, j);
    }
  }
  // close files
  infile_chi.close();
  infile_cr.close();


  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        // density
        phydro->u(IDN, k, j, i) = nH_min * std::pow(10, (i-is)*fac_nH);
        // energy
        if (NON_BAROTROPIC_EOS) {
          pres = phydro->u(IDN, k, j, i) * SQR(iso_cs);
          phydro->u(IEN, k, j, i) = pres/gm1;
        }
      }
    }
  }

  // intialize radiation field
  if (CHEMRADIATION_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          for (int ifreq=0; ifreq < pchemrad->nfreq; ++ifreq) {
            for (int iang=0; iang < pchemrad->nang; ++iang) {
              pchemrad->ir(k, j, i, ifreq * pchemrad->nang + iang)
                = chi(k-ks, j-js);
            }
          }
          if (CHEMISTRY_ENABLED) {
            for (int iang=0; iang < pchemrad->nang; ++iang) {
              // cr rate
              pchemrad->ir(k, j, i, pscalars->chemnet.index_cr_ * pchemrad->nang + iang)
                = cr(k-ks, j-js);
            }
          }
        }
      }
    }
    // calculate the average radiation field for output of the initial condition
    pchemrad->pchemradintegrator->CopyToOutput();
  }

  // intialize chemical species
  if (NSPECIES > 0) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          for (int ispec=0; ispec < NSPECIES; ++ispec) {
            pscalars->s(ispec, k, j, i) = r_init * phydro->u(IDN, k, j, i);
            if (CHEMISTRY_ENABLED) {
              Real s_ispec = pin->GetOrAddReal("problem",
                  "r_init_"+pscalars->chemnet.species_names[ispec], -1);
              if (s_ispec >= 0.) {
                pscalars->s(ispec, k, j, i) = s_ispec * phydro->u(IDN, k, j, i);
              }
            }
          }
        }
      }
    }
  }

  return;
}
