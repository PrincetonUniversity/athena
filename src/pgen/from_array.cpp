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

  // Prepare index bounds
  int il = is - NGHOST;
  int iu = ie + NGHOST;
  int jl = js;
  int ju = je;
  if (block_size.nx2 > 1) {
    jl -= (NGHOST);
    ju += (NGHOST);
  }
  int kl = ks;
  int ku = ke;
  if (block_size.nx3 > 1) {
    kl -= (NGHOST);
    ku += (NGHOST);
  }

  // Prepare scratch arrays for conserved values
  AthenaArray<Real> cons;
  cons.NewAthenaArray(NHYDRO, 1, block_size.nx3, block_size.nx2, block_size.nx1);

  // Read conserved values from file
  std::string input_filename = pin->GetString("problem", "input_filename");
  std::string dataset_name = pin->GetString("problem", "dataset_name");
  int start_file[5] = {0};
  start_file[1] = gid;
  int start_mem[5] = {0};
  int count[5];
  count[0] = NHYDRO;
  count[1] = 1;
  count[2] = block_size.nx3;
  count[3] = block_size.nx2;
  count[4] = block_size.nx1;
  HDF5ReadRealArray(input_filename.c_str(), dataset_name.c_str(), 5, start_file,
      start_mem, count, cons);

  // Set conserved values
  for (int n = 0; n <= NHYDRO; ++n) {
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie; ++i) {
          phydro->u(n, k, j, i) = cons(n, 0, k-NGHOST, j-NGHOST, i-NGHOST);
        }
      }
    }
  }
  return;
}
