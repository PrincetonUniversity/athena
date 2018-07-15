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

  // Determine location of initial values
  std::string input_filename = pin->GetString("problem", "input_filename");
  std::string dataset_name = pin->GetString("problem", "dataset_name");

  // Prepare array selections
  int start_file[5] = {0};
  start_file[1] = gid;
  int count_file[5];
  count_file[0] = NHYDRO;
  count_file[1] = 1;
  count_file[2] = block_size.nx3;
  count_file[3] = block_size.nx2;
  count_file[4] = block_size.nx1;
  int start_mem[5] = {0};
  start_mem[2] = ks;
  start_mem[3] = js;
  start_mem[4] = is;
  int count_mem[5];
  count_mem[0] = 1;
  count_mem[1] = NHYDRO;
  count_mem[2] = block_size.nx3;
  count_mem[3] = block_size.nx2;
  count_mem[4] = block_size.nx1;

  // Read conserved values from file
  HDF5ReadRealArray(input_filename.c_str(), dataset_name.c_str(), 5, start_file,
      count_file, start_mem, count_mem, phydro->u);
  return;
}
