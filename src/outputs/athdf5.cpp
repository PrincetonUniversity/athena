//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================

// HDF5 output

// Only proceed if HDF5 output enabled
#include "../athena.hpp"  // enums, macros, LogicalLocation
#ifdef HDF5OUTPUT

// Primary header
#include "outputs.hpp"

// C++ headers
#include <cstdio>     // sprintf()
#include <cstring>    // strlen(), strncpy()
#include <fstream>    // ofstream
#include <iomanip>    // setfill(), setw()
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // string

// External library headers
#include <hdf5.h>  // H5[F|P|S|T]_*, H5[A|D|F|P|S|T]*(), hid_t
#ifdef MPI_PARALLEL
#include <mpi.h>   // MPI_COMM_WORLD, MPI_INFO_NULL
#endif

// Athena++ headers
#include "../athena_arrays.hpp"            // AthenaArray
#include "../globals.hpp"                  // Globals
#include "../mesh.hpp"                     // Mesh, MeshBlock, RegionSize
#include "../parameter_input.hpp"          // ParameterInput
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../field/field.hpp"              // Field
#include "../hydro/hydro.hpp"              // Hydro

//--------------------------------------------------------------------------------------

// Mesh-level initialization of output
// Inputs:
//   pmesh: pointer to Mesh
//   pin: pointer to inputs (unused)
//   walltime_limit: flag indicating file written because walltime ran out (unused)
// Outputs: (none)
// Notes:
//   filename: <problem_id>.out<n>.<ddddd>.athdf
//     n: arbitrary-length integer corresponding to output block
//     ddddd: 0-padded length-5 integer incremented each output
//   opens file and writes file-level attributes
//   creates datasets in file
//   prepares dataspaces for both file and memory
void ATHDF5Output::Initialize(Mesh *pmesh, ParameterInput *pin,
    bool walltime_limit=false)
{
  // Determine number and sizes of blocks
  num_blocks_global = pmesh->nbtotal;
  num_blocks_local = pmesh->nblist[Globals::my_rank];
  int first_block = pmesh->nslist[Globals::my_rank];
  MeshBlock *pblock = pmesh->pblock;
  nx1 = pblock->block_size.nx1;
  nx2 = pblock->block_size.nx2;
  nx3 = pblock->block_size.nx3;
  is = pblock->is;
  ie = pblock->ie;
  js = pblock->js;
  je = pblock->je;
  ks = pblock->ks;
  ke = pblock->ke;
  root_level = pmesh->root_level;

  // Determine what variables to output
  std::string variable = output_params.variable;
  if (variable.compare("prim") == 0 or variable.compare("cons") == 0)
  {
    num_datasets = 1;
    if (MAGNETIC_FIELDS_ENABLED)
      num_datasets = 2;
    num_variables = new int[num_datasets];
    num_variables[0] = NHYDRO;
    num_total_variables = NHYDRO;
    if (MAGNETIC_FIELDS_ENABLED)
    {
      num_variables[1] = 3;
      num_total_variables += 3;
    }
    dataset_names = new char[num_datasets][max_name_length+1];
    variable_names = new char[num_total_variables][max_name_length+1];
    if (variable.compare("prim") == 0)
    {
      std::strncpy(dataset_names[0], "prim", max_name_length+1);
      for (int n = 0; n < NHYDRO; ++n)
        switch (n)
        {
          case IDN:
            std::strncpy(variable_names[n], "rho", max_name_length+1);
            break;
          case IEN:
            std::strncpy(variable_names[n], "pgas", max_name_length+1);
            break;
          case IM1:
            std::strncpy(variable_names[n], "vel1", max_name_length+1);
            break;
          case IM2:
            std::strncpy(variable_names[n], "vel2", max_name_length+1);
            break;
          case IM3:
            std::strncpy(variable_names[n], "vel3", max_name_length+1);
            break;
        }
    }
    else
    {
      std::strncpy(dataset_names[0], "cons", max_name_length+1);
      for (int n = 0; n < NHYDRO; ++n)
        switch (n)
        {
          case IDN:
            std::strncpy(variable_names[n], "dens", max_name_length+1);
            break;
          case IEN:
            std::strncpy(variable_names[n], "Etot", max_name_length+1);
            break;
          case IM1:
            std::strncpy(variable_names[n], "mom1", max_name_length+1);
            break;
          case IM2:
            std::strncpy(variable_names[n], "mom2", max_name_length+1);
            break;
          case IM3:
            std::strncpy(variable_names[n], "mom3", max_name_length+1);
            break;
        }
    }
    if (MAGNETIC_FIELDS_ENABLED)
    {
      std::strncpy(dataset_names[1], "B", max_name_length+1);
      std::strncpy(variable_names[NHYDRO], "B1", max_name_length+1);
      std::strncpy(variable_names[NHYDRO+1], "B2", max_name_length+1);
      std::strncpy(variable_names[NHYDRO+2], "B3", max_name_length+1);
    }
  }
  else if (variable.compare("d") == 0)
  {
    num_datasets = 1;
    num_variables = new int[num_datasets];
    num_variables[0] = 1;
    num_total_variables = 1;
    dataset_names = new char[num_datasets][max_name_length+1];
    variable_names = new char[num_total_variables][max_name_length+1];
    std::strncpy(dataset_names[0], "d", max_name_length+1);
    std::strncpy(variable_names[0], "rho", max_name_length+1);
  }
  else if (NON_BAROTROPIC_EOS and variable.compare("p") == 0)
  {
    num_datasets = 1;
    num_variables = new int[num_datasets];
    num_variables[0] = 1;
    num_total_variables = 1;
    dataset_names = new char[num_datasets][max_name_length+1];
    variable_names = new char[num_total_variables][max_name_length+1];
    std::strncpy(dataset_names[0], "p", max_name_length+1);
    std::strncpy(variable_names[0], "pgas", max_name_length+1);
  }
  else if (variable.compare("v") == 0)
  {
    num_datasets = 1;
    num_variables = new int[num_datasets];
    num_variables[0] = 3;
    num_total_variables = 3;
    dataset_names = new char[num_datasets][max_name_length+1];
    variable_names = new char[num_total_variables][max_name_length+1];
    std::strncpy(dataset_names[0], "v", max_name_length+1);
    std::strncpy(variable_names[0], "vel1", max_name_length+1);
    std::strncpy(variable_names[1], "vel2", max_name_length+1);
    std::strncpy(variable_names[2], "vel3", max_name_length+1);
  }
  else if (variable.compare("D") == 0)
  {
    num_datasets = 1;
    num_variables = new int[num_datasets];
    num_variables[0] = 1;
    num_total_variables = 1;
    dataset_names = new char[num_datasets][max_name_length+1];
    variable_names = new char[num_total_variables][max_name_length+1];
    std::strncpy(dataset_names[0], "D", max_name_length+1);
    std::strncpy(variable_names[0], "dens", max_name_length+1);
  }
  else if (NON_BAROTROPIC_EOS and variable.compare("E") == 0)
  {
    num_datasets = 1;
    num_variables = new int[num_datasets];
    num_variables[0] = 1;
    num_total_variables = 1;
    dataset_names = new char[num_datasets][max_name_length+1];
    variable_names = new char[num_total_variables][max_name_length+1];
    std::strncpy(dataset_names[0], "E", max_name_length+1);
    std::strncpy(variable_names[0], "Etot", max_name_length+1);
  }
  else if (variable.compare("m") == 0)
  {
    num_datasets = 1;
    num_variables = new int[num_datasets];
    num_variables[0] = 3;
    num_total_variables = 3;
    dataset_names = new char[num_datasets][max_name_length+1];
    variable_names = new char[num_total_variables][max_name_length+1];
    std::strncpy(dataset_names[0], "m", max_name_length+1);
    std::strncpy(variable_names[0], "mom1", max_name_length+1);
    std::strncpy(variable_names[1], "mom2", max_name_length+1);
    std::strncpy(variable_names[2], "mom3", max_name_length+1);
  }
  else if (MAGNETIC_FIELDS_ENABLED and variable.compare("b") == 0)
  {
    num_datasets = 1;
    num_variables = new int[num_datasets];
    num_variables[0] = 3;
    num_total_variables = 3;
    dataset_names = new char[num_datasets][max_name_length+1];
    variable_names = new char[num_total_variables][max_name_length+1];
    std::strncpy(dataset_names[0], "B", max_name_length+1);
    std::strncpy(variable_names[0], "B1", max_name_length+1);
    std::strncpy(variable_names[1], "B2", max_name_length+1);
    std::strncpy(variable_names[2], "B3", max_name_length+1);
  }
  else if (variable.compare("ifov") == 0)
  {
    if (NIFOV <= 0)
    {
      std::stringstream message;
      message << "### FATAL ERROR in athdf5 initialization\n"
              << "No variables to output\n";
      throw std::runtime_error(message.str().c_str());
    }
    int max_ifov_digits = max_name_length - std::strlen("ifov");
    int max_ifov = 1;
    for (int n = 0; n < max_ifov_digits; ++n)
      max_ifov *= 10;
    if (max_ifov_digits <= 0)
      max_ifov = 0;
    if (NIFOV > max_ifov)
    {
      std::stringstream message;
      message << "### FATAL ERROR in athdf5 initialization\n"
              << "Can only support " << max_ifov << " ifov outputs\n";
      throw std::runtime_error(message.str().c_str());
    }
    num_datasets = 1;
    num_variables = new int[num_datasets];
    num_variables[0] = NIFOV;
    num_total_variables = NIFOV;
    dataset_names = new char[num_datasets][max_name_length+1];
    variable_names = new char[num_total_variables][max_name_length+1];
    std::strncpy(dataset_names[0], "ifov", max_name_length+1);
    for (int n = 0; n < NIFOV; ++n)
      std::sprintf(variable_names[n], "ifov%d", n);
  }
  else
  {
    std::stringstream message;
    message << "### FATAL ERROR in athdf5 initialization\n"
            << "Output variable " << variable << " unrecognized\n";
    throw std::runtime_error(message.str().c_str());
  }

  // Make sure C-strings are null-terminated
  for (int n = 0; n < num_datasets; ++n)
    dataset_names[n][max_name_length] = '\0';
  for (int n = 0; n < num_total_variables; ++n)
    variable_names[n][max_name_length] = '\0';

  // Define output filename
  filename = std::string(output_params.file_basename);
  filename.append(".");
  filename.append(output_params.file_id);
  filename.append(".");
  std::stringstream file_number;
  file_number << std::setw(5) << std::setfill('0') << output_params.file_number;
  filename.append(file_number.str());
  filename.append(".athdf");

  // Create new file
  #ifdef MPI_PARALLEL
  {
    hid_t property_list_file = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(property_list_file, MPI_COMM_WORLD, MPI_INFO_NULL);
    file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, property_list_file);
    H5Pclose(property_list_file);
  }
  #else
    file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  #endif
  if (file < 0)
  {
    std::stringstream message;
    message << "### FATAL ERROR in athdf5 initialization\n"
            << "Could not open " << filename << "\n";
    throw std::runtime_error(message.str().c_str());
  }

  // Prepare datatypes and dataspaces for writing attributes
  hid_t string_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type, max_name_length+1);
  hid_t dataspace_scalar = H5Screate(H5S_SCALAR);
  dims_count[0] = 3;
  hid_t dataspace_triple = H5Screate_simple(1, dims_count, NULL);
  dims_count[0] = num_datasets;
  hid_t dataspace_dataset_list = H5Screate_simple(1, dims_count, NULL);
  dims_count[0] = num_total_variables;
  hid_t dataspace_variable_list = H5Screate_simple(1, dims_count, NULL);

  // Write cycle number
  int num_cycles = pmesh->ncycle;
  hid_t attribute = H5Acreate2(file, "NumCycles", H5T_STD_I32BE, dataspace_scalar,
      H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute, H5T_NATIVE_INT, &num_cycles);
  H5Aclose(attribute);

  // Write simulation time
  double time = pmesh->time;
  attribute = H5Acreate2(file, "Time", H5T_IEEE_F64BE, dataspace_scalar, H5P_DEFAULT,
      H5P_DEFAULT);
  H5Awrite(attribute, H5T_NATIVE_DOUBLE, &time);
  H5Aclose(attribute);

  // Write coordinate system
  if (std::strlen(COORDINATE_SYSTEM) > max_name_length)
  {
    std::stringstream message;
    message << "### FATAL ERROR in athdf5 initialization\n"
            << "Coordinate name too long\n";
    throw std::runtime_error(message.str().c_str());
  }
  attribute = H5Acreate2(file, "Coordinates", string_type, dataspace_scalar,
      H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute, string_type, COORDINATE_SYSTEM);
  H5Aclose(attribute);

  // Write extent of grid in x1-direction
  double coord_range[3];
  coord_range[0] = pmesh->mesh_size.x1min;
  coord_range[1] = pmesh->mesh_size.x1max;
  coord_range[2] = pmesh->mesh_size.x1rat;
  attribute = H5Acreate2(file, "RootGridX1", H5T_IEEE_F64BE, dataspace_triple,
      H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute, H5T_NATIVE_DOUBLE, coord_range);
  H5Aclose(attribute);

  // Write extent of grid in x2-direction
  coord_range[0] = pmesh->mesh_size.x2min;
  coord_range[1] = pmesh->mesh_size.x2max;
  coord_range[2] = pmesh->mesh_size.x2rat;
  attribute = H5Acreate2(file, "RootGridX2", H5T_IEEE_F64BE, dataspace_triple,
      H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute, H5T_NATIVE_DOUBLE, coord_range);
  H5Aclose(attribute);

  // Write extent of grid in x3-direction
  coord_range[0] = pmesh->mesh_size.x3min;
  coord_range[1] = pmesh->mesh_size.x3max;
  coord_range[2] = pmesh->mesh_size.x3rat;
  attribute = H5Acreate2(file, "RootGridX3", H5T_IEEE_F64BE, dataspace_triple,
      H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute, H5T_NATIVE_DOUBLE, coord_range);
  H5Aclose(attribute);

  // Write root grid size
  int root_grid_size[3];
  root_grid_size[0] = pmesh->mesh_size.nx1;
  root_grid_size[1] = pmesh->mesh_size.nx2;
  root_grid_size[2] = pmesh->mesh_size.nx3;
  attribute = H5Acreate2(file, "RootGridSize", H5T_STD_I32BE, dataspace_triple,
      H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute, H5T_NATIVE_INT, root_grid_size);
  H5Aclose(attribute);

  // Write number of MeshBlocks
  attribute = H5Acreate2(file, "NumMeshBlocks", H5T_STD_I32BE, dataspace_scalar,
      H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute, H5T_NATIVE_INT, &num_blocks_global);
  H5Aclose(attribute);

  // Write MeshBlock size
  int meshblock_size[3];
  meshblock_size[0] = nx1;
  meshblock_size[1] = nx2;
  meshblock_size[2] = nx3;
  attribute = H5Acreate2(file, "MeshBlockSize", H5T_STD_I32BE, dataspace_triple,
      H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute, H5T_NATIVE_INT, meshblock_size);
  H5Aclose(attribute);

  // Write maximum refinement level
  int max_level = pmesh->current_level - pmesh->root_level;
  attribute = H5Acreate2(file, "MaxLevel", H5T_STD_I32BE, dataspace_scalar, H5P_DEFAULT,
      H5P_DEFAULT);
  H5Awrite(attribute, H5T_NATIVE_INT, &max_level);
  H5Aclose(attribute);

  // Write number of output cell-centered variables
  attribute = H5Acreate2(file, "NumVariables", H5T_STD_I32BE, dataspace_dataset_list,
      H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute, H5T_NATIVE_INT, num_variables);
  H5Aclose(attribute);

  // Write names of datasets in same order
  attribute = H5Acreate2(file, "DatasetNames", string_type, dataspace_dataset_list,
      H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute, string_type, dataset_names);
  H5Aclose(attribute);

  // Write array of variable names
  attribute = H5Acreate2(file, "VariableNames", string_type, dataspace_variable_list,
      H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute, string_type, variable_names);
  H5Aclose(attribute);

  // Close attribute dataspaces
  H5Sclose(dataspace_scalar);
  H5Sclose(dataspace_triple);
  H5Sclose(dataspace_dataset_list);
  H5Sclose(dataspace_variable_list);

  // Prepare global (full) dataspaces for writing datasets
  dims_count[0] = num_blocks_global;
  filespace_blocks = H5Screate_simple(1, dims_count, NULL);
  dims_count[1] = 3;
  filespace_blocks_3 = H5Screate_simple(2, dims_count, NULL);
  dims_count[1] = nx1+1;
  filespace_blocks_nx1 = H5Screate_simple(2, dims_count, NULL);
  dims_count[1] = nx2+1;
  filespace_blocks_nx2 = H5Screate_simple(2, dims_count, NULL);
  dims_count[1] = nx3+1;
  filespace_blocks_nx3 = H5Screate_simple(2, dims_count, NULL);
  dims_count[2] = nx3;
  dims_count[3] = nx2;
  dims_count[4] = nx1;
  filespaces_blocks_vars_nx3_nx2_nx1 = new hid_t[num_datasets];
  for (int n = 0; n < num_datasets; ++n)
  {
    dims_count[1] = num_variables[n];
    filespaces_blocks_vars_nx3_nx2_nx1[n] = H5Screate_simple(5, dims_count, NULL);
  }

  // Create datasets
  dataset_levels = H5Dcreate(file, "Levels", H5T_STD_I32BE, filespace_blocks,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dataset_locations = H5Dcreate(file, "LogicalLocations", H5T_STD_I64BE,
      filespace_blocks_3, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dataset_x1f = H5Dcreate(file, "x1f", H5T_IEEE_F32BE, filespace_blocks_nx1,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dataset_x2f = H5Dcreate(file, "x2f", H5T_IEEE_F32BE, filespace_blocks_nx2,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dataset_x3f = H5Dcreate(file, "x3f", H5T_IEEE_F32BE, filespace_blocks_nx3,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  datasets_celldata = new hid_t[num_datasets];
  for (int n = 0; n < num_datasets; ++n)
    datasets_celldata[n] = H5Dcreate(file, dataset_names[n], H5T_IEEE_F32BE,
        filespaces_blocks_vars_nx3_nx2_nx1[n], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Prepare local (hyperslabbed) dataspaces for writing datasets to file
  dims_start[0] = first_block;
  dims_count[0] = num_blocks_local;
  H5Sselect_hyperslab(filespace_blocks, H5S_SELECT_SET, dims_start, NULL, dims_count,
      NULL);
  dims_start[1] = 0;
  dims_count[1] = 3;
  H5Sselect_hyperslab(filespace_blocks_3, H5S_SELECT_SET, dims_start, NULL, dims_count,
      NULL);
  dims_count[1] = nx1+1;
  H5Sselect_hyperslab(filespace_blocks_nx1, H5S_SELECT_SET, dims_start, NULL,
      dims_count, NULL);
  dims_count[1] = nx2+1;
  H5Sselect_hyperslab(filespace_blocks_nx2, H5S_SELECT_SET, dims_start, NULL,
      dims_count, NULL);
  dims_count[1] = nx3+1;
  H5Sselect_hyperslab(filespace_blocks_nx3, H5S_SELECT_SET, dims_start, NULL,
      dims_count, NULL);
  dims_start[2] = 0;
  dims_start[3] = 0;
  dims_start[4] = 0;
  dims_count[2] = nx3;
  dims_count[3] = nx2;
  dims_count[4] = nx1;
  for (int n = 0; n < num_datasets; ++n)
  {
    dims_count[1] = num_variables[n];
    H5Sselect_hyperslab(filespaces_blocks_vars_nx3_nx2_nx1[n], H5S_SELECT_SET,
        dims_start, NULL, dims_count, NULL);
  }

  // Allocate contiguous buffers for data in memory
  levels_mesh = new int[num_blocks_local];
  locations_mesh = new long int[num_blocks_local * 3];
  x1f_mesh = new float[num_blocks_local * (nx1+1)];
  x2f_mesh = new float[num_blocks_local * (nx2+1)];
  x3f_mesh = new float[num_blocks_local * (nx3+1)];
  data_buffers = new float *[num_datasets];
  for (int n = 0; n < num_datasets; ++n)
    data_buffers[n] = new float[num_blocks_local * num_variables[n] * nx3 * nx2 * nx1];
  data_arrays = new AthenaArray<Real>[num_datasets];

  // Prepare dataspaces for describing data in memory
  dims_count[0] = num_blocks_local;
  memspace_blocks = H5Screate_simple(1, dims_count, NULL);
  dims_count[1] = 3;
  memspace_blocks_3 = H5Screate_simple(2, dims_count, NULL);
  dims_count[1] = nx1+1;
  memspace_blocks_nx1 = H5Screate_simple(2, dims_count, NULL);
  dims_count[1] = nx2+1;
  memspace_blocks_nx2 = H5Screate_simple(2, dims_count, NULL);
  dims_count[1] = nx3+1;
  memspace_blocks_nx3 = H5Screate_simple(2, dims_count, NULL);
  memspaces_blocks_vars_nx3_nx2_nx1 = new hid_t[num_datasets];
  for (int n = 0; n < num_datasets; ++n)
  {
    dims_count[1] = num_variables[n];
    dims_count[2] = nx3;
    dims_count[3] = nx2;
    dims_count[4] = nx1;
    memspaces_blocks_vars_nx3_nx2_nx1[n] = H5Screate_simple(5, dims_count, NULL);
  }

  // Use collective MPI calls if applicable
  property_list = H5Pcreate(H5P_DATASET_XFER);
  #ifdef MPI_PARALLEL
    H5Pset_dxpl_mpio(property_list, H5FD_MPIO_COLLECTIVE);
  #endif
  return;
}

//--------------------------------------------------------------------------------------

// Mesh-level finalization of output
// Inputs:
//   pin: pointer to inputs
// Outputs: (none)
// Notes:
//   has single process write .athdf.xdmf file
void ATHDF5Output::Finalize(ParameterInput *pin)
{
  // Close property list
  H5Pclose(property_list);

  // Close dataspaces for describing memory
  H5Sclose(memspace_blocks);
  H5Sclose(memspace_blocks_3);
  H5Sclose(memspace_blocks_nx1);
  H5Sclose(memspace_blocks_nx2);
  H5Sclose(memspace_blocks_nx3);
  for (int n = 0; n < num_datasets; ++n)
    H5Sclose(memspaces_blocks_vars_nx3_nx2_nx1[n]);
  delete[] memspaces_blocks_vars_nx3_nx2_nx1;

  // Close dataspaces for describing file
  H5Sclose(filespace_blocks);
  H5Sclose(filespace_blocks_3);
  H5Sclose(filespace_blocks_nx1);
  H5Sclose(filespace_blocks_nx2);
  H5Sclose(filespace_blocks_nx3);
  for (int n = 0; n < num_datasets; ++n)
    H5Sclose(filespaces_blocks_vars_nx3_nx2_nx1[n]);
  delete[] filespaces_blocks_vars_nx3_nx2_nx1;

  // Close datasets
  H5Dclose(dataset_levels);
  H5Dclose(dataset_locations);
  H5Dclose(dataset_x1f);
  H5Dclose(dataset_x2f);
  H5Dclose(dataset_x3f);
  for (int n = 0; n < num_datasets; ++n)
    H5Dclose(datasets_celldata[n]);
  delete[] datasets_celldata;

  // Close .athdf file
  H5Fclose(file);

  // Write .athdf.xdmf file
  int write_xdmf = pin->GetOrAddInteger(output_params.block_name, "xdmf", 1);
  if(Globals::my_rank == 0 && write_xdmf != 0)
    MakeXDMF();

  // Delete data storage
  delete[] num_variables;
  delete[] dataset_names;
  delete[] variable_names;
  delete[] levels_mesh;
  delete[] locations_mesh;
  delete[] x1f_mesh;
  delete[] x2f_mesh;
  delete[] x3f_mesh;
  for (int n = 0; n < num_datasets; ++n)
    delete[] data_buffers[n];
  delete[] data_buffers;
  delete[] data_arrays;

  // Reset parameters for next time file is written
  output_params.file_number++;
  output_params.next_time += output_params.dt;
  pin->SetInteger(output_params.block_name, "file_number", output_params.file_number);
  pin->SetReal(output_params.block_name, "next_time", output_params.next_time);
  return;
}

//--------------------------------------------------------------------------------------

// Block-level preparation of data to be output
// Inputs:
//   pout_data: pointer to output data (unused)
//   pblock: pointer to MeshBlock
// Outputs: (none)
// Notes:
//   operates on all blocks when called with first block
//   does nothing when called for all subsequent blocks
//   all data must be loaded to allow WriteOutputFile() to write all data the first time
//       it is called
void ATHDF5Output::LoadOutputData(OutputData *pout_data, MeshBlock *pblock)
{
  // Do nothing except for first block on Mesh
  if (pblock->lid != 0)
    return;

  // Go through all blocks
  MeshBlock *pblock_current = pblock;
  for (int lid = 0; lid < num_blocks_local; ++lid, pblock_current=pblock_current->next)
  {
    // Determine where variables to be output are stored
    std::string variable = output_params.variable;
    if (variable.compare("prim") == 0)
    {
      data_arrays[0].InitWithShallowSlice(pblock_current->phydro->w, 4, 0, NHYDRO);
      if (MAGNETIC_FIELDS_ENABLED)
        data_arrays[1].InitWithShallowSlice(pblock_current->pfield->bcc, 4, 0, 3);
    }
    else if (variable.compare("cons") == 0)
    {
      data_arrays[0].InitWithShallowSlice(pblock_current->phydro->u, 4, 0, NHYDRO);
      if (MAGNETIC_FIELDS_ENABLED)
        data_arrays[1].InitWithShallowSlice(pblock_current->pfield->bcc, 4, 0, 3);
    }
    else if (variable.compare("d") == 0)
      data_arrays[0].InitWithShallowSlice(pblock_current->phydro->w, 4, IDN, 1);
    else if (variable.compare("p") == 0)
      data_arrays[0].InitWithShallowSlice(pblock_current->phydro->w, 4, IEN, 1);
    else if (variable.compare("v") == 0)
      data_arrays[0].InitWithShallowSlice(pblock_current->phydro->w, 4, IM1, 3);
    else if (variable.compare("D") == 0)
      data_arrays[0].InitWithShallowSlice(pblock_current->phydro->u, 4, IDN, 1);
    else if (variable.compare("E") == 0)
      data_arrays[0].InitWithShallowSlice(pblock_current->phydro->u, 4, IEN, 1);
    else if (variable.compare("m") == 0)
      data_arrays[0].InitWithShallowSlice(pblock_current->phydro->u, 4, IM1, 3);
    else if (variable.compare("b") == 0)
      data_arrays[0].InitWithShallowSlice(pblock_current->pfield->bcc, 4, 0, 3);
    else
      data_arrays[0].InitWithShallowSlice(pblock_current->phydro->ifov, 4, 0, NIFOV);

    // Load location information
    levels_mesh[lid] = pblock_current->loc.level - root_level;
    locations_mesh[lid*3 + 0] = pblock_current->loc.lx1;
    locations_mesh[lid*3 + 1] = pblock_current->loc.lx2;
    locations_mesh[lid*3 + 2] = pblock_current->loc.lx3;

    // Load coordinates (implicitly converting to float)
    for (int i=is, index=0; i <= ie+1; ++i, ++index)
      x1f_mesh[lid*(nx1+1) + index] = pblock_current->pcoord->x1f(i);
    for (int j=js, index=0; j <= je+1; ++j, ++index)
      x2f_mesh[lid*(nx2+1) + index] = pblock_current->pcoord->x2f(j);
    for (int k=ks, index=0; k <= ke+1; ++k, ++index)
      x3f_mesh[lid*(nx3+1) + index] = pblock_current->pcoord->x3f(k);

    // Load celldata (implicitly converting to float)
    for (int n = 0; n < num_datasets; ++n)
    {
      int index = 0;
      for (int v = 0; v < num_variables[n]; ++v)
        for (int k = ks; k <= ke; ++k)
          for (int j = js; j <= je; ++j)
            for (int i = is; i <= ie; ++i)
              data_buffers[n][lid*num_variables[n]*nx3*nx2*nx1 + index++]
                  = data_arrays[n](v,k,j,i);
      data_arrays[n].DeleteAthenaArray();
    }
  }
  return;
}

//--------------------------------------------------------------------------------------

// Block-level writing of data to be output
// Inputs:
//   pout_data: pointer to output data (unused)
//   pblock: pointer to MeshBlock
// Outputs: (none)
// Notes:
//   writes the following data to file for all N blocks on this Mesh:
//     refinement levels (N, int)
//     logical locations (N x 3, long int)
//     x1f (N x (nx1+1), float)
//     x2f (N x (nx2+1), float)
//     x3f (N x (nx3+1), float)
//     each desired dataset (N x nvar x nx3 x nx2 x nx1, float)
//   operates on all blocks when called with first block
//   does nothing when called for all subsequent blocks
//   assumes LoadOutputData() does the same
void ATHDF5Output::WriteOutputFile(OutputData *pout_data, MeshBlock *pblock)
{
  // Do nothing except for first block on Mesh
  if (pblock->lid != 0)
    return;

  // Write refinement level and logical location
  H5Dwrite(dataset_levels, H5T_NATIVE_INT, memspace_blocks, filespace_blocks,
      property_list, levels_mesh);
  H5Dwrite(dataset_locations, H5T_NATIVE_LONG, memspace_blocks_3, filespace_blocks_3,
      property_list, locations_mesh);

  // Write coordinates
  H5Dwrite(dataset_x1f, H5T_NATIVE_FLOAT, memspace_blocks_nx1, filespace_blocks_nx1,
      property_list, x1f_mesh);
  H5Dwrite(dataset_x2f, H5T_NATIVE_FLOAT, memspace_blocks_nx2, filespace_blocks_nx2,
      property_list, x2f_mesh);
  H5Dwrite(dataset_x3f, H5T_NATIVE_FLOAT, memspace_blocks_nx3, filespace_blocks_nx3,
      property_list, x3f_mesh);

  // Write cell data
  for (int n = 0; n < num_datasets; ++n)
    H5Dwrite(datasets_celldata[n], H5T_NATIVE_FLOAT,
        memspaces_blocks_vars_nx3_nx2_nx1[n], filespaces_blocks_vars_nx3_nx2_nx1[n],
        property_list, data_buffers[n]);
  return;
}

//--------------------------------------------------------------------------------------

// Function for writing auxiliary XDMF metadata file
// Inputs: (none)
// Outputs: (none)
// Notes:
//   writes .athdf.xdmf file for describing .athdf file
//   should only be called by single process
//   file size scales proportional to total number of MeshBlocks
//   for many small MeshBlocks, this can take most of the output writing time
void ATHDF5Output::MakeXDMF()
{
  // Open auxiliary file for writing
  std::string filename_aux(filename);
  filename_aux.append(".xdmf");
  std::ofstream xdmf(filename_aux.c_str());

  // Write header
  xdmf << "<?xml version=\"1.0\" ?>\n";
  xdmf << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
  xdmf << "<Xdmf Version=\"2.0\">\n";
  xdmf << "<Domain>\n";
  xdmf << "<Grid Name=\"Mesh\" GridType=\"Collection\">\n";

  // Go through all MeshBlocks
  for (int n_block = 0; n_block < num_blocks_global; ++n_block)
  {
    // Begin block
    xdmf << "  <Grid Name=\"MeshBlock" << n_block << "\" GridType=\"Uniform\">\n";

    // Write topology
    if (nx3 > 1)
      xdmf << "    <Topology TopologyType=\"3DRectMesh\" NumberOfElements=\"" << nx3+1
          << " " << nx2+1 << " " << nx1+1 << "\"/>\n";
    else
      xdmf << "    <Topology TopologyType=\"2DRectMesh\" NumberOfElements=\"" << nx2+1
          << " " << nx1+1 << "\"/>\n";

    // Write geometry
    if (nx3 > 1)
      xdmf << "    <Geometry GeometryType=\"VXVYVZ\">\n";
    else
      xdmf << "    <Geometry GeometryType=\"VXVY\">\n";
    xdmf << "      <DataItem ItemType=\"HyperSlab\" Dimensions=\"" << nx1+1 << "\">\n";
    xdmf << "        <DataItem Dimensions=\"3 2\" NumberType=\"Int\"> " << n_block
        << " 0 1 1 1 " << nx1+1 << " </DataItem>\n";
    xdmf << "        <DataItem Dimensions=\"" << num_blocks_global << " " << nx1+1
        << "\" Format=\"HDF\"> " << filename << ":/x1f </DataItem>\n";
    xdmf << "      </DataItem>\n";
    xdmf << "      <DataItem ItemType=\"HyperSlab\" Dimensions=\"" << nx2+1 << "\">\n";
    xdmf << "        <DataItem Dimensions=\"3 2\" NumberType=\"Int\"> " << n_block
        << " 0 1 1 1 " << nx2+1 << " </DataItem>\n";
    xdmf << "        <DataItem Dimensions=\"" << num_blocks_global << " " << nx2+1
        << "\" Format=\"HDF\"> " << filename << ":/x2f </DataItem>\n";
    xdmf << "      </DataItem>\n";
    if (nx3 > 1)
    {
      xdmf << "      <DataItem ItemType=\"HyperSlab\" Dimensions=\"" << nx3+1
          << "\">\n";
      xdmf << "        <DataItem Dimensions=\"3 2\" NumberType=\"Int\"> " << n_block
          << " 0 1 1 1 " << nx3+1 << " </DataItem>\n";
      xdmf << "        <DataItem Dimensions=\"" << num_blocks_global << " " << nx3+1
          << "\" Format=\"HDF\"> " << filename << ":/x3f </DataItem>\n";
      xdmf << "      </DataItem>\n";
    }
    xdmf << "    </Geometry>\n";

    // Write description of cell-centered data
    int n_quantity = 0;
    for (int n_dataset = 0; n_dataset < num_datasets; ++n_dataset)
      for (int n_variable = 0; n_variable < num_variables[n_dataset]; ++n_variable)
      {
        xdmf << "    <Attribute Name=\"" << variable_names[n_quantity]
            << "\" Center=\"Cell\">\n";
        if (nx3 > 1)
        {
          xdmf << "      <DataItem ItemType=\"HyperSlab\" Dimensions=\"" << nx3 << " "
              << nx2 << " " << nx1 << "\">\n";
          xdmf << "        <DataItem Dimensions=\"3 5\" NumberType=\"Int\"> " << n_block
              << " " << n_quantity << " 0 0 0 1 1 1 1 1 1 1 " << nx3 << " " << nx2
              << " " << nx1 << " </DataItem>\n";
        }
        else
        {
          xdmf << "      <DataItem ItemType=\"HyperSlab\" Dimensions=\"" << nx2 << " "
              << nx1 << "\">\n";
          xdmf << "        <DataItem Dimensions=\"3 5\" NumberType=\"Int\"> " << n_block
              << " " << n_quantity << " 0 0 0 1 1 1 1 1 1 1 1 " << nx2 << " " << nx1
              << " </DataItem>\n";
        }
        xdmf << "        <DataItem Dimensions=\"" << num_blocks_global << " "
            << num_total_variables << " " << nx3 << " " << nx2 << " " << nx1
            << "\" Format=\"HDF\"> " << filename << ":/" << dataset_names[n_dataset]
            << " </DataItem>\n";
        xdmf << "      </DataItem>\n";
        xdmf << "    </Attribute>\n";
        ++n_quantity;
      }

    // End block
    xdmf << "  </Grid>\n";
  }

  // Complete header elements
  xdmf << "</Grid>\n";
  xdmf << "</Domain>\n";
  xdmf << "</Xdmf>";

  // Close file
  xdmf.close();
  return;
}

#endif  // HDF5OUTPUT
