//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file from_array.cpp
//  \brief Problem generator for initializing with preexisting array

// C++ headers
#include <string>     // c_str(), string

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "../athena.hpp"                   // Real
#include "../athena_arrays.hpp"            // AthenaArray
#include "../parameter_input.hpp"          // ParameterInput
#include "../field/field.hpp"              // Field
#include "../hydro/hydro.hpp"              // Hydro
#include "../inputs/hdf5_reader.hpp"       // HDF5ReadRealArray()

//----------------------------------------------------------------------------------------
// Function for setting initial conditions
// Inputs:
//   pin: parameters
// Outputs: (none)
// Notes:
//   uses input parameters to determine which file contains array of conserved values
//   dataset must be 5-dimensional array with the following sizes:
//     NHYDRO
//     total number of MeshBlocks
//     MeshBlock/nx3
//     MeshBlock/nx2
//     MeshBlock/nx1

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  // Determine locations of initial values
  std::string input_filename = pin->GetString("problem", "input_filename");
  std::string dataset_cons = pin->GetString("problem", "dataset_cons");
  int index_dens = pin->GetInteger("problem", "index_dens");
  int index_mom1 = pin->GetInteger("problem", "index_mom1");
  int index_mom2 = pin->GetInteger("problem", "index_mom2");
  int index_mom3 = pin->GetInteger("problem", "index_mom3");
  int index_etot = pin->GetInteger("problem", "index_etot");

  // Prepare array selections
  int start_file[5] = {0};
  start_file[1] = gid;
  int start_mem[5] = {0};
  start_mem[2] = ks;
  start_mem[3] = js;
  start_mem[4] = is;
  int start_indices[5];
  start_indices[IDN] = index_dens;
  start_indices[IM1] = index_mom1;
  start_indices[IM2] = index_mom2;
  start_indices[IM3] = index_mom3;
  start_indices[IEN] = index_etot;
  int count[5];
  count[0] = 1;
  count[1] = 1;
  count[2] = block_size.nx3;
  count[3] = block_size.nx2;
  count[4] = block_size.nx1;

  // Read conserved values from file
  for (int n = 0; n < NHYDRO; ++n) {
    start_file[0] = start_indices[n];
    start_mem[1] = n;
    HDF5ReadRealArray(input_filename.c_str(), dataset_cons.c_str(), 5, start_file, count,
        start_mem, count, phydro->u);
  }
  return;
}
