//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file from_array.cpp
//! \brief Problem generator for initializing with preexisting array from HDF5 input

// C headers

// C++ headers
#include <algorithm>  // max()
#include <string>     // c_str(), string

// Athena++ headers
#include "../athena.hpp"              // Real
#include "../athena_arrays.hpp"       // AthenaArray
#include "../field/field.hpp"         // Field
#include "../globals.hpp"             // Globals
#include "../hydro/hydro.hpp"         // Hydro
#include "../inputs/hdf5_reader.hpp"  // HDF5ReadRealArray()
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"     // ParameterInput

//----------------------------------------------------------------------------------------
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Function for setting initial conditions
//!
//! Inputs:
//! - pin: parameters
//! Outputs: (none)
//! Notes:
//! - uses input parameters to determine which file contains array of conserved values
//!   dataset must be 5-dimensional array with the following sizes:
//!   - NHYDRO
//!   - total number of MeshBlocks
//!   - MeshBlock/nx3
//!   - MeshBlock/nx2
//!   - MeshBlock/nx1

void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (SELF_GRAVITY_ENABLED) {
    Real four_pi_G = pin->GetReal("problem","four_pi_G");
    Real eps = pin->GetOrAddReal("problem","grav_eps", 0.0);
    SetFourPiG(four_pi_G);
    SetGravityThreshold(eps);
  }
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Determine locations of initial values
  std::string input_filename = pin->GetString("problem", "input_filename");
  std::string dataset_cons = pin->GetString("problem", "dataset_cons");
  int index_dens = pin->GetInteger("problem", "index_dens");
  int index_mom1 = pin->GetInteger("problem", "index_mom1");
  int index_mom2 = pin->GetInteger("problem", "index_mom2");
  int index_mom3 = pin->GetInteger("problem", "index_mom3");
  int index_etot = pin->GetInteger("problem", "index_etot");
  std::string dataset_b1 = pin->GetString("problem", "dataset_b1");
  std::string dataset_b2 = pin->GetString("problem", "dataset_b2");
  std::string dataset_b3 = pin->GetString("problem", "dataset_b3");

  // Set conserved array selections
  int start_cons_file[5];
  start_cons_file[1] = gid;
  start_cons_file[2] = 0;
  start_cons_file[3] = 0;
  start_cons_file[4] = 0;
  int start_cons_indices[5];
  start_cons_indices[IDN] = index_dens;
  start_cons_indices[IM1] = index_mom1;
  start_cons_indices[IM2] = index_mom2;
  start_cons_indices[IM3] = index_mom3;
  start_cons_indices[IEN] = index_etot;
  int count_cons_file[5];
  count_cons_file[0] = 1;
  count_cons_file[1] = 1;
  count_cons_file[2] = block_size.nx3;
  count_cons_file[3] = block_size.nx2;
  count_cons_file[4] = block_size.nx1;
  int start_cons_mem[4];
  start_cons_mem[1] = ks;
  start_cons_mem[2] = js;
  start_cons_mem[3] = is;
  int count_cons_mem[4];
  count_cons_mem[0] = 1;
  count_cons_mem[1] = block_size.nx3;
  count_cons_mem[2] = block_size.nx2;
  count_cons_mem[3] = block_size.nx1;

  // Set conserved values from file
  for (int n = 0; n < NHYDRO; ++n) {
    start_cons_file[0] = start_cons_indices[n];
    start_cons_mem[0] = n;
    HDF5ReadRealArray(input_filename.c_str(), dataset_cons.c_str(), 5, start_cons_file,
                      count_cons_file, 4, start_cons_mem,
                      count_cons_mem, phydro->u, true);
  }

  // Set field array selections
  int start_field_file[4];
  start_field_file[0] = gid;
  start_field_file[1] = 0;
  start_field_file[2] = 0;
  start_field_file[3] = 0;
  int count_field_file[4];
  count_field_file[0] = 1;
  int start_field_mem[3];
  start_field_mem[0] = ks;
  start_field_mem[1] = js;
  start_field_mem[2] = is;
  int count_field_mem[3];

  // Set magnetic field values from file
  if (MAGNETIC_FIELDS_ENABLED) {
    // Set B1
    count_field_file[1] = block_size.nx3;
    count_field_file[2] = block_size.nx2;
    count_field_file[3] = block_size.nx1 + 1;
    count_field_mem[0] = block_size.nx3;
    count_field_mem[1] = block_size.nx2;
    count_field_mem[2] = block_size.nx1 + 1;
    HDF5ReadRealArray(input_filename.c_str(), dataset_b1.c_str(), 4, start_field_file,
                      count_field_file, 3, start_field_mem,
                      count_field_mem, pfield->b.x1f, true);

    // Set B2
    count_field_file[1] = block_size.nx3;
    count_field_file[2] = block_size.nx2 + 1;
    count_field_file[3] = block_size.nx1;
    count_field_mem[0] = block_size.nx3;
    count_field_mem[1] = block_size.nx2 + 1;
    count_field_mem[2] = block_size.nx1;
    HDF5ReadRealArray(input_filename.c_str(), dataset_b2.c_str(), 4, start_field_file,
                      count_field_file, 3, start_field_mem,
                      count_field_mem, pfield->b.x2f, true);

    // Set B3
    count_field_file[1] = block_size.nx3 + 1;
    count_field_file[2] = block_size.nx2;
    count_field_file[3] = block_size.nx1;
    count_field_mem[0] = block_size.nx3 + 1;
    count_field_mem[1] = block_size.nx2;
    count_field_mem[2] = block_size.nx1;
    HDF5ReadRealArray(input_filename.c_str(), dataset_b3.c_str(), 4, start_field_file,
                      count_field_file, 3, start_field_mem,
                      count_field_mem, pfield->b.x3f, true);
  }

  // Make no-op collective reads if using MPI and ranks have unequal numbers of blocks
#ifdef MPI_PARALLEL
  {
    int num_blocks_this_rank = pmy_mesh->nblist[Globals::my_rank];
    if (lid == num_blocks_this_rank - 1) {
      int block_shortage_this_rank = 0;
      for (int rank = 0; rank < Globals::nranks; ++rank) {
        block_shortage_this_rank =
            std::max(block_shortage_this_rank,
                     pmy_mesh->nblist[rank] - num_blocks_this_rank);
      }
      for (int block = 0; block < block_shortage_this_rank; ++block) {
        for (int n = 0; n < NHYDRO; ++n) {
          start_cons_file[0] = start_cons_indices[n];
          start_cons_mem[0] = n;
          HDF5ReadRealArray(input_filename.c_str(), dataset_cons.c_str(), 5,
                            start_cons_file, count_cons_file, 4,
                            start_cons_mem, count_cons_mem,
                            phydro->u, true, true);
        }
        if (MAGNETIC_FIELDS_ENABLED) {
          count_field_file[1] = block_size.nx3;
          count_field_file[2] = block_size.nx2;
          count_field_file[3] = block_size.nx1 + 1;
          count_field_mem[0] = block_size.nx3;
          count_field_mem[1] = block_size.nx2;
          count_field_mem[2] = block_size.nx1 + 1;
          HDF5ReadRealArray(input_filename.c_str(), dataset_b1.c_str(), 4,
                            start_field_file, count_field_file, 3,
                            start_field_mem, count_field_mem,
                            pfield->b.x1f, true, true);
          count_field_file[1] = block_size.nx3;
          count_field_file[2] = block_size.nx2 + 1;
          count_field_file[3] = block_size.nx1;
          count_field_mem[0] = block_size.nx3;
          count_field_mem[1] = block_size.nx2 + 1;
          count_field_mem[2] = block_size.nx1;
          HDF5ReadRealArray(input_filename.c_str(), dataset_b2.c_str(), 4,
                            start_field_file, count_field_file, 3,
                            start_field_mem, count_field_mem,
                            pfield->b.x2f, true, true);
          count_field_file[1] = block_size.nx3 + 1;
          count_field_file[2] = block_size.nx2;
          count_field_file[3] = block_size.nx1;
          count_field_mem[0] = block_size.nx3 + 1;
          count_field_mem[1] = block_size.nx2;
          count_field_mem[2] = block_size.nx1;
          HDF5ReadRealArray(input_filename.c_str(), dataset_b3.c_str(), 4,
                            start_field_file, count_field_file, 3,
                            start_field_mem, count_field_mem,
                            pfield->b.x3f, true, true);
        }
      }
    }
  }
#endif
  return;
}
