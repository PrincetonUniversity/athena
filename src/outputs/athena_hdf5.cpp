//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file athena_hdf5.cpp
//! \brief hdf5 outputs

// C headers

// C++ headers
#include <cstdio>     // snprintf()
#include <cstring>    // strlen(), strncpy()
#include <fstream>    // ofstream
#include <iomanip>    // setfill(), setw()
#include <iostream>   // cout
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // string

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "outputs.hpp"

// Only proceed if HDF5 output enabled
#ifdef HDF5OUTPUT

// External library headers
#include <hdf5.h>  // H5[F|P|S|T]_*, H5[A|D|F|P|S|T]*(), hid_t
#ifdef MPI_PARALLEL
#include <mpi.h>   // MPI_COMM_WORLD, MPI_INFO_NULL
#endif

// type alias that allows HDF5 output to be written in either floats or doubles
#if H5_DOUBLE_PRECISION_ENABLED
using H5Real = double;
#if SINGLE_PRECISION_ENABLED
#error "Cannot create HDF5 output at higher precision than internal representation"
#endif
#define H5T_NATIVE_REAL H5T_NATIVE_DOUBLE

#else
using H5Real = float;
#define H5T_NATIVE_REAL H5T_NATIVE_FLOAT
#endif


//----------------------------------------------------------------------------------------
//! \fn void ATHDF5Output:::WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag)
//! \brief Cycles over all MeshBlocks and writes OutputData in the Athena++ HDF5 format,
//!        one file per output using parallel IO.

void ATHDF5Output::WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag) {
  // HDF5 structures
  hid_t file;                                  // file to be written to
  hsize_t dims_start[5], dims_count[5];        // array sizes
  hid_t dataset_levels;                        // datasets to be written
  hid_t dataset_locations;
  hid_t dataset_x1f, dataset_x2f, dataset_x3f;
  hid_t dataset_x1v, dataset_x2v, dataset_x3v;
  hid_t *datasets_celldata;
  hid_t filespace_blocks;                      // local dataspaces for file
  hid_t filespace_blocks_3;
  hid_t filespace_blocks_nx1,  filespace_blocks_nx2,  filespace_blocks_nx3;
  hid_t filespace_blocks_nx1v, filespace_blocks_nx2v, filespace_blocks_nx3v;
  hid_t *filespaces_vars_blocks_nx3_nx2_nx1;
  hid_t memspace_blocks;                       // local dataspaces for memory
  hid_t memspace_blocks_3;
  hid_t memspace_blocks_nx1,  memspace_blocks_nx2,  memspace_blocks_nx3;
  hid_t memspace_blocks_nx1v, memspace_blocks_nx2v, memspace_blocks_nx3v;
  hid_t *memspaces_vars_blocks_nx3_nx2_nx1;
  hid_t property_list;                         // properties for writing

  int num_blocks_local;                        // number of MeshBlocks on this Mesh
  int *levels_mesh;                            // array of refinement levels on Mesh
  std::int64_t *locations_mesh;                // array of logical locations on Mesh
  H5Real *x1f_mesh;                            // array of x1 values on Mesh
  H5Real *x2f_mesh;                            // array of x2 values on Mesh
  H5Real *x3f_mesh;                            // array of x3 values on Mesh
  H5Real *x1v_mesh;                            // array of x1 values on Mesh
  H5Real *x2v_mesh;                            // array of x2 values on Mesh
  H5Real *x3v_mesh;                            // array of x3 values on Mesh
  H5Real **data_buffers;                       // array of data buffers

  MeshBlock *pmb = pm->my_blocks(0);
  OutputData* pod;
  int max_blocks_global = pm->nbtotal;
  int max_blocks_local = pm->nblist[Globals::my_rank];
  int first_block = pm->nslist[Globals::my_rank];
  bool *active_flags = new bool[max_blocks_local];
  std::string variable = output_params.variable;

  for (int i=0; i<max_blocks_local; i++)
    active_flags[i] = true;

  // shooting a blank just for getting the variable names
  out_is = pmb->is; out_ie = pmb->ie;
  out_js = pmb->js; out_je = pmb->je;
  out_ks = pmb->ks; out_ke = pmb->ke;
  if (output_params.include_ghost_zones) {
    out_is -= NGHOST; out_ie += NGHOST;
    if (out_js != out_je) {out_js -= NGHOST; out_je += NGHOST;}
    if (out_ks != out_ke) {out_ks -= NGHOST; out_ke += NGHOST;}
  }
  LoadOutputData(pmb);
  // set num_datasets and num_variables

  // NEW_OUTPUT_TYPES:
  // this must be expanded when new variables are introduced

  if (variable.compare("prim") == 0 || variable.compare("cons") == 0) {
    int n_dataset = 0;
    // cell-centered data packed into 1 dataset: (prim, phi, s0, etc.)
    num_datasets = 1;
    // face-centered data packed into 1 dataset:
    if (MAGNETIC_FIELDS_ENABLED)
      num_datasets += 1;
    num_variables = new int[num_datasets];

    // n_dataset = 0: all cell-centered AthenaArray variable data of the same size
    // Hydro conserved variables:
    num_variables[n_dataset] = NHYDRO;
    if (output_params.cartesian_vector)
      num_variables[n_dataset] += 3;
    // Graviatational potential:
    if (SELF_GRAVITY_ENABLED)
      num_variables[n_dataset] += 1;
    // Passive scalars:
    if (NSCALARS > 0)
      num_variables[n_dataset] += NSCALARS;

    // n_dataset = 1: face-centered FaceField variable data
    n_dataset++;
    // Longitudinal, face-centered magnetic field components:
    if (MAGNETIC_FIELDS_ENABLED) {
      num_variables[n_dataset] = 3;
      if (output_params.cartesian_vector)
        num_variables[n_dataset] += 3;
    }
  } else {
    num_datasets = 1;
    num_variables = new int[num_datasets];
    num_variables[0] = num_vars_;
  }
  dataset_names = new char[num_datasets][max_name_length+1];
  variable_names = new char[num_vars_][max_name_length+1];

  // set dataset names
  int n_dataset_names = 0;
  if (variable.compare("prim") == 0 || variable.compare("cons") == 0) {
    if (variable.compare("prim") == 0)
      std::strncpy(dataset_names[n_dataset_names++], "prim", max_name_length+1);
    else
      std::strncpy(dataset_names[n_dataset_names++], "cons", max_name_length+1);
    if (MAGNETIC_FIELDS_ENABLED)
      std::strncpy(dataset_names[n_dataset_names++], "B", max_name_length+1);
  } else { // single data
    if (variable.compare(0,1,"B") == 0 && MAGNETIC_FIELDS_ENABLED)
      std::strncpy(dataset_names[n_dataset_names++], "B", max_name_length+1);
    else if (variable.compare(0,3,"uov") == 0
             || variable.compare(0,12,"user_out_var") == 0)
      std::strncpy(dataset_names[n_dataset_names++], "user_out_var", max_name_length+1);
    else
      std::strncpy(dataset_names[n_dataset_names++], "hydro", max_name_length+1);
  }

  // set variable names, loop over outputdata
  int n_variable = 0;
  pod = pfirst_data_;
  while (pod != nullptr) {
    if (pod->type == "VECTORS") {
      for (int i=1; i<=3; i++) {
        char sn[3];
        std::snprintf(sn, sizeof(sn), "%d", i);
        std::string vname = pod->name + sn;
        std::strncpy(variable_names[n_variable++], vname.c_str(), max_name_length+1);
      }
    } else {
      std::strncpy(variable_names[n_variable++], pod->name.c_str(), max_name_length+1);
    }
    pod = pod->pnext;
  }

  // Make sure C-strings are null-terminated
  for (int n = 0; n < num_datasets; ++n)
    dataset_names[n][max_name_length] = '\0';
  for (int n = 0; n < num_vars_; ++n)
    variable_names[n][max_name_length] = '\0';

  ClearOutputData();

  // count the number of active blocks if slicing
  if (output_params.output_slicex1 || output_params.output_slicex2
      || output_params.output_slicex3) {
    int nb = 0, nba = 0;
    for (int b=0; b<pm->nblocal; ++b) {
      pmb = pm->my_blocks(b);
      if (output_params.output_slicex1) {
        if (pmb->block_size.x1min >  output_params.x1_slice
            || pmb->block_size.x1max <= output_params.x1_slice)
          active_flags[nb] = false;
      }
      if (output_params.output_slicex2) {
        if (pmb->block_size.x2min >  output_params.x2_slice
            || pmb->block_size.x2max <= output_params.x2_slice)
          active_flags[nb] = false;
      }
      if (output_params.output_slicex3) {
        if (pmb->block_size.x3min >  output_params.x3_slice
            || pmb->block_size.x3max <= output_params.x3_slice)
          active_flags[nb] = false;
      }
      if (active_flags[nb]) nba++;
      nb++;
    }
#ifdef MPI_PARALLEL
    int *n_active = new int[Globals::nranks];
    MPI_Allgather(&nba, 1, MPI_INT, n_active, 1, MPI_INT, MPI_COMM_WORLD);
    num_blocks_local = n_active[Globals::my_rank];
    first_block = 0;
    for (int n=0; n<Globals::my_rank; n++)
      first_block += n_active[n];
    num_blocks_global = 0;
    for (int n=0; n<Globals::nranks; n++)
      num_blocks_global += n_active[n];
    delete [] n_active;
#else
    num_blocks_global=num_blocks_local=nba;
    first_block = 0;
#endif
  } else {
    num_blocks_global = max_blocks_global;
    num_blocks_local = max_blocks_local;
  }

  pm->my_blocks(0);
  // set output size
  nx1 = pmb->block_size.nx1;
  nx2 = pmb->block_size.nx2;
  nx3 = pmb->block_size.nx3;
  if (output_params.include_ghost_zones) {
    nx1 += 2*NGHOST;
    if (nx2 > 1) nx2 += 2*NGHOST;
    if (nx3 > 1) nx3 += 2*NGHOST;
  }
  if (output_params.output_slicex1) nx1=1;
  if (output_params.output_slicex2) nx2=1;
  if (output_params.output_slicex3) nx3=1;
  if (output_params.output_sumx1) nx1=1;
  if (output_params.output_sumx2) nx2=1;
  if (output_params.output_sumx3) nx3=1;

  // Allocate contiguous buffers for data in memory
  levels_mesh = new int[num_blocks_local];
  locations_mesh = new std::int64_t[num_blocks_local * 3];
  x1f_mesh = new H5Real[num_blocks_local * (nx1+1)];
  x2f_mesh = new H5Real[num_blocks_local * (nx2+1)];
  x3f_mesh = new H5Real[num_blocks_local * (nx3+1)];
  x1v_mesh = new H5Real[num_blocks_local * nx1];
  x2v_mesh = new H5Real[num_blocks_local * nx2];
  x3v_mesh = new H5Real[num_blocks_local * nx3];
  data_buffers = new H5Real *[num_datasets];
  for (int n = 0; n < num_datasets; ++n)
    data_buffers[n] = new H5Real[num_variables[n]*num_blocks_local*nx3*nx2*nx1];

  int nb = 0, nba = 0;
  for (int b=0; b<pm->nblocal; ++b) {
    pmb = pm->my_blocks(b);
    // Load the output data
    if (active_flags[nb]) {
      // set the default size because TransformOutputData will override it
      out_is = pmb->is; out_ie = pmb->ie;
      out_js = pmb->js; out_je = pmb->je;
      out_ks = pmb->ks; out_ke = pmb->ke;
      if (output_params.include_ghost_zones) {
        out_is -= NGHOST; out_ie += NGHOST;
        if (out_js != out_je) {out_js -= NGHOST; out_je += NGHOST;}
        if (out_ks != out_ke) {out_ks -= NGHOST; out_ke += NGHOST;}
      }
      LoadOutputData(pmb);
      TransformOutputData(pmb);
      if (output_params.output_sumx1) {
        out_ie = out_is;
      }
      if (output_params.output_sumx2) {
        out_je = out_js;
      }
      if (output_params.output_sumx3) {
        out_ke = out_ks;
      }

      // Load location information
      levels_mesh[nba] = pmb->loc.level - pm->root_level;
      locations_mesh[nba*3 + 0] = pmb->loc.lx1;
      locations_mesh[nba*3 + 1] = pmb->loc.lx2;
      locations_mesh[nba*3 + 2] = pmb->loc.lx3;

      // Load coordinates
      if (output_params.output_slicex1) {
        x1f_mesh[nba*(nx1+1)] =
            static_cast<H5Real>(pmb->pcoord->x1f(output_params.islice));
        x1f_mesh[nba*(nx1+1)+1] = static_cast<H5Real>(
            pmb->pcoord->x1f(output_params.islice+1));
        x1v_mesh[nba*nx1] = static_cast<H5Real>(pmb->pcoord->x1v(output_params.islice));
      } else if (output_params.output_sumx1) {
        x1f_mesh[nba*(nx1+1)] = pmb->pcoord->x1f(pmb->is);
        x1f_mesh[nba*(nx1+1)+1] = pmb->pcoord->x1f(pmb->ie+1);
        if (pmb->block_size.nx1 % 2 == 0) {
          x1v_mesh[nba*nx1] = pmb->pcoord->x1f((pmb->is + pmb->ie + 1) / 2);
        } else {
          x1v_mesh[nba*nx1] = pmb->pcoord->x1v((pmb->is + pmb->ie) / 2);
        }
      } else {
        for (int i=out_is, index=0; i <= out_ie+1; ++i, ++index)
          x1f_mesh[nba*(nx1+1) + index] = static_cast<H5Real>(pmb->pcoord->x1f(i));
        for (int i=out_is, index=0; i <= out_ie; ++i, ++index)
          x1v_mesh[nba*nx1 + index] = static_cast<H5Real>(pmb->pcoord->x1v(i));
      }
      if (output_params.output_slicex2) {
        x2f_mesh[nba*(nx2+1)]
            = static_cast<H5Real>(pmb->pcoord->x2f(output_params.jslice));
        x2f_mesh[nba*(nx2+1)+1]
            = static_cast<H5Real>(pmb->pcoord->x2f(output_params.jslice+1));
        x2v_mesh[nba*nx2]
            = static_cast<H5Real>(pmb->pcoord->x2v(output_params.jslice));
      } else if (output_params.output_sumx2) {
        x2f_mesh[nba*(nx2+1)] = pmb->pcoord->x2f(pmb->js);
        x2f_mesh[nba*(nx2+1)+1] = pmb->pcoord->x2f(pmb->je+1);
        if (pmb->block_size.nx2 % 2 == 0) {
          x2v_mesh[nba*nx2] = pmb->pcoord->x2f((pmb->js + pmb->je + 1) / 2);
        } else {
          x2v_mesh[nba*nx2] = pmb->pcoord->x2v((pmb->js + pmb->je) / 2);
        }
      } else {
        for (int j=out_js, index=0; j <= out_je+1; ++j, ++index)
          x2f_mesh[nba*(nx2+1) + index] = static_cast<H5Real>(pmb->pcoord->x2f(j));
        for (int j=out_js, index=0; j <= out_je; ++j, ++index)
          x2v_mesh[nba*nx2 + index] = static_cast<H5Real>(pmb->pcoord->x2v(j));
      }
      if (output_params.output_slicex3) {
        x3f_mesh[nba*(nx3+1)]
            = static_cast<H5Real>(pmb->pcoord->x3f(output_params.kslice));
        x3f_mesh[nba*(nx3+1)+1]
            = static_cast<H5Real>(pmb->pcoord->x3f(output_params.kslice+1));
        x3v_mesh[nba*nx3]
            = static_cast<H5Real>(pmb->pcoord->x3v(output_params.kslice));
      } else if (output_params.output_sumx3) {
        x3f_mesh[nba*(nx3+1)] = pmb->pcoord->x3f(pmb->ks);
        x3f_mesh[nba*(nx3+1)+1] = pmb->pcoord->x3f(pmb->ke+1);
        if (pmb->block_size.nx3 % 2 == 0) {
          x3v_mesh[nba*nx3] = pmb->pcoord->x3f((pmb->ks + pmb->ke + 1) / 2);
        } else {
          x3v_mesh[nba*nx3] = pmb->pcoord->x3v((pmb->ks + pmb->ke) / 2);
        }
      } else {
        for (int k=out_ks, index=0; k <= out_ke+1; ++k, ++index)
          x3f_mesh[nba*(nx3+1) + index] = static_cast<H5Real>(pmb->pcoord->x3f(k));
        for (int k=out_ks, index=0; k <= out_ke; ++k, ++index)
          x3v_mesh[nba*nx3 + index] = static_cast<H5Real>(pmb->pcoord->x3v(k));
      }

      // store the data into the data_buffers
      if (variable.compare("prim") == 0 || variable.compare("cons") == 0) {
        int n_dataset = 0;
        int ndv = 0;
        pod = pfirst_data_;
        while (pod != nullptr) {
          if (pod->name == "Bcc") {
            n_dataset++;
            ndv = 0;
          }
          int nv=1;
          if (pod->type == "VECTORS") nv=3;
          for (int v=0; v < nv; v++, ndv++) {
            int index = 0;
            for (int k = out_ks; k <= out_ke; k++) {
              for (int j = out_js; j <= out_je; j++) {
                for (int i = out_is; i <= out_ie; i++, index++)
                  data_buffers[n_dataset][(ndv*num_blocks_local+nba)*nx3*nx2*nx1+index]
                      = pod->data(v,k,j,i);
              }
            }
          }
          pod = pod->pnext;
        }
      } else {
        int ndv = 0;
        pod = pfirst_data_;
        while (pod != nullptr) {
          int nv=1;
          if (pod->type == "VECTORS") nv=3;
          for (int v=0; v < nv; v++, ndv++) {
            int index = 0;
            for (int k = out_ks; k <= out_ke; k++) {
              for (int j = out_js; j <= out_je; j++) {
                for (int i = out_is; i <= out_ie; i++, index++)
                  data_buffers[0][(ndv*num_blocks_local+nba)*nx3*nx2*nx1+index]
                      = pod->data(v,k,j,i);
              }
            }
          }
          pod = pod->pnext;
        }
      }
      nba++;
      ClearOutputData();  // required when LoadOutputData() is used.
    }
    nb++;
  }

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
  hid_t property_list_file = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(property_list_file, MPI_COMM_WORLD, MPI_INFO_NULL);
  file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, property_list_file);
  H5Pclose(property_list_file);
#else
  file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
#endif
  if (file < 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in athdf5 initialization\n"
        << "Could not open " << filename << std::endl;
    ATHENA_ERROR(msg);
  }

  // Prepare datatypes and dataspaces for writing attributes
  hid_t string_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type, max_name_length+1);
  hid_t dataspace_scalar = H5Screate(H5S_SCALAR);
  dims_count[0] = 3;
  hid_t dataspace_triple = H5Screate_simple(1, dims_count, NULL);
  dims_count[0] = num_datasets;
  hid_t dataspace_dataset_list = H5Screate_simple(1, dims_count, NULL);
  dims_count[0] = num_vars_;
  hid_t dataspace_variable_list = H5Screate_simple(1, dims_count, NULL);

  // Write cycle number
  int num_cycles = pm->ncycle;
  hid_t attribute = H5Acreate2(file, "NumCycles", H5T_STD_I32BE, dataspace_scalar,
                               H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute, H5T_NATIVE_INT, &num_cycles);
  H5Aclose(attribute);

  // Write simulation time
  double time = pm->time;
  attribute = H5Acreate2(file, "Time", H5T_NATIVE_REAL, dataspace_scalar, H5P_DEFAULT,
                         H5P_DEFAULT);
  H5Awrite(attribute, H5T_NATIVE_DOUBLE, &time);
  H5Aclose(attribute);
  code_time = static_cast<H5Real>(time); // output time for xdmf

  // Write coordinate system
  if (std::strlen(COORDINATE_SYSTEM) > max_name_length) {
    std::stringstream msg;
    msg << "### FATAL ERROR in athdf5 initialization\n"
        << "Coordinate name too long\n";
    ATHENA_ERROR(msg);
  }
  attribute = H5Acreate2(file, "Coordinates", string_type, dataspace_scalar,
                         H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute, string_type, COORDINATE_SYSTEM);
  H5Aclose(attribute);

  // Write extent of grid in x1-direction
  double coord_range[3];
  coord_range[0] = pm->mesh_size.x1min;
  coord_range[1] = pm->mesh_size.x1max;
  coord_range[2] = pm->mesh_size.x1rat;
  attribute = H5Acreate2(file, "RootGridX1", H5T_NATIVE_REAL, dataspace_triple,
                         H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute, H5T_NATIVE_DOUBLE, coord_range);
  H5Aclose(attribute);

  // Write extent of grid in x2-direction
  coord_range[0] = pm->mesh_size.x2min;
  coord_range[1] = pm->mesh_size.x2max;
  coord_range[2] = pm->mesh_size.x2rat;
  attribute = H5Acreate2(file, "RootGridX2", H5T_NATIVE_REAL, dataspace_triple,
                         H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute, H5T_NATIVE_DOUBLE, coord_range);
  H5Aclose(attribute);

  // Write extent of grid in x3-direction
  coord_range[0] = pm->mesh_size.x3min;
  coord_range[1] = pm->mesh_size.x3max;
  coord_range[2] = pm->mesh_size.x3rat;
  attribute = H5Acreate2(file, "RootGridX3", H5T_NATIVE_REAL, dataspace_triple,
                         H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute, H5T_NATIVE_DOUBLE, coord_range);
  H5Aclose(attribute);

  // Write root grid size
  int root_grid_size[3];
  root_grid_size[0] = pm->mesh_size.nx1;
  root_grid_size[1] = pm->mesh_size.nx2;
  root_grid_size[2] = pm->mesh_size.nx3;
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
  int max_level = pm->current_level - pm->root_level;
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
  dims_count[1] = nx1;
  filespace_blocks_nx1v = H5Screate_simple(2, dims_count, NULL);
  dims_count[1] = nx2;
  filespace_blocks_nx2v = H5Screate_simple(2, dims_count, NULL);
  dims_count[1] = nx3;
  filespace_blocks_nx3v = H5Screate_simple(2, dims_count, NULL);
  dims_count[1] = num_blocks_global;
  dims_count[2] = nx3;
  dims_count[3] = nx2;
  dims_count[4] = nx1;
  filespaces_vars_blocks_nx3_nx2_nx1 = new hid_t[num_datasets];
  for (int n = 0; n < num_datasets; ++n) {
    dims_count[0] = num_variables[n];
    filespaces_vars_blocks_nx3_nx2_nx1[n] = H5Screate_simple(5, dims_count, NULL);
  }

  // Create datasets
  dataset_levels = H5Dcreate(file, "Levels", H5T_STD_I32BE, filespace_blocks,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dataset_locations = H5Dcreate(file, "LogicalLocations", H5T_STD_I64BE,
                                filespace_blocks_3,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dataset_x1f = H5Dcreate(file, "x1f", H5T_NATIVE_REAL, filespace_blocks_nx1,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dataset_x2f = H5Dcreate(file, "x2f", H5T_NATIVE_REAL, filespace_blocks_nx2,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dataset_x3f = H5Dcreate(file, "x3f", H5T_NATIVE_REAL, filespace_blocks_nx3,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dataset_x1v = H5Dcreate(file, "x1v", H5T_NATIVE_REAL, filespace_blocks_nx1v,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dataset_x2v = H5Dcreate(file, "x2v", H5T_NATIVE_REAL, filespace_blocks_nx2v,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dataset_x3v = H5Dcreate(file, "x3v", H5T_NATIVE_REAL, filespace_blocks_nx3v,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  datasets_celldata = new hid_t[num_datasets];
  for (int n = 0; n < num_datasets; ++n)
    datasets_celldata[n] = H5Dcreate(file, dataset_names[n], H5T_NATIVE_REAL,
                                     filespaces_vars_blocks_nx3_nx2_nx1[n],
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

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
  dims_count[1] = nx1;
  H5Sselect_hyperslab(filespace_blocks_nx1v, H5S_SELECT_SET, dims_start, NULL,
                      dims_count, NULL);
  dims_count[1] = nx2;
  H5Sselect_hyperslab(filespace_blocks_nx2v, H5S_SELECT_SET, dims_start, NULL,
                      dims_count, NULL);
  dims_count[1] = nx3;
  H5Sselect_hyperslab(filespace_blocks_nx3v, H5S_SELECT_SET, dims_start, NULL,
                      dims_count, NULL);
  dims_start[0] = 0;
  dims_start[1] = first_block;
  dims_start[2] = 0;
  dims_start[3] = 0;
  dims_start[4] = 0;
  dims_count[1] = num_blocks_local;
  dims_count[2] = nx3;
  dims_count[3] = nx2;
  dims_count[4] = nx1;
  for (int n = 0; n < num_datasets; ++n) {
    dims_count[0] = num_variables[n];
    H5Sselect_hyperslab(filespaces_vars_blocks_nx3_nx2_nx1[n], H5S_SELECT_SET,
                        dims_start, NULL, dims_count, NULL);
  }

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
  dims_count[1] = nx1;
  memspace_blocks_nx1v = H5Screate_simple(2, dims_count, NULL);
  dims_count[1] = nx2;
  memspace_blocks_nx2v = H5Screate_simple(2, dims_count, NULL);
  dims_count[1] = nx3;
  memspace_blocks_nx3v = H5Screate_simple(2, dims_count, NULL);
  dims_count[1] = num_blocks_local;
  dims_count[2] = nx3;
  dims_count[3] = nx2;
  dims_count[4] = nx1;
  memspaces_vars_blocks_nx3_nx2_nx1 = new hid_t[num_datasets];
  for (int n = 0; n < num_datasets; ++n) {
    dims_count[0] = num_variables[n];
    memspaces_vars_blocks_nx3_nx2_nx1[n] = H5Screate_simple(5, dims_count, NULL);
  }

  // Use collective MPI calls if applicable
  property_list = H5Pcreate(H5P_DATASET_XFER);
#ifdef MPI_PARALLEL
  H5Pset_dxpl_mpio(property_list, H5FD_MPIO_COLLECTIVE);
#endif

  // dump all the data
  // Write refinement level and logical location
  H5Dwrite(dataset_levels, H5T_NATIVE_INT, memspace_blocks, filespace_blocks,
           property_list, levels_mesh);
  H5Dwrite(dataset_locations, H5T_NATIVE_LONG, memspace_blocks_3, filespace_blocks_3,
           property_list, locations_mesh);

  // Write coordinates
  H5Dwrite(dataset_x1f, H5T_NATIVE_REAL, memspace_blocks_nx1, filespace_blocks_nx1,
           property_list, x1f_mesh);
  H5Dwrite(dataset_x2f, H5T_NATIVE_REAL, memspace_blocks_nx2, filespace_blocks_nx2,
           property_list, x2f_mesh);
  H5Dwrite(dataset_x3f, H5T_NATIVE_REAL, memspace_blocks_nx3, filespace_blocks_nx3,
           property_list, x3f_mesh);
  H5Dwrite(dataset_x1v, H5T_NATIVE_REAL, memspace_blocks_nx1v, filespace_blocks_nx1v,
           property_list, x1v_mesh);
  H5Dwrite(dataset_x2v, H5T_NATIVE_REAL, memspace_blocks_nx2v, filespace_blocks_nx2v,
           property_list, x2v_mesh);
  H5Dwrite(dataset_x3v, H5T_NATIVE_REAL, memspace_blocks_nx3v, filespace_blocks_nx3v,
           property_list, x3v_mesh);

  // Write cell data
  for (int n = 0; n < num_datasets; ++n)
    H5Dwrite(datasets_celldata[n], H5T_NATIVE_REAL,
             memspaces_vars_blocks_nx3_nx2_nx1[n], filespaces_vars_blocks_nx3_nx2_nx1[n],
             property_list, data_buffers[n]);


  // Close property list
  H5Pclose(property_list);

  // Close dataspaces for describing memory
  H5Sclose(memspace_blocks);
  H5Sclose(memspace_blocks_3);
  H5Sclose(memspace_blocks_nx1);
  H5Sclose(memspace_blocks_nx2);
  H5Sclose(memspace_blocks_nx3);
  H5Sclose(memspace_blocks_nx1v);
  H5Sclose(memspace_blocks_nx2v);
  H5Sclose(memspace_blocks_nx3v);
  for (int n = 0; n < num_datasets; ++n)
    H5Sclose(memspaces_vars_blocks_nx3_nx2_nx1[n]);
  delete[] memspaces_vars_blocks_nx3_nx2_nx1;

  // Close dataspaces for describing file
  H5Sclose(filespace_blocks);
  H5Sclose(filespace_blocks_3);
  H5Sclose(filespace_blocks_nx1);
  H5Sclose(filespace_blocks_nx2);
  H5Sclose(filespace_blocks_nx3);
  H5Sclose(filespace_blocks_nx1v);
  H5Sclose(filespace_blocks_nx2v);
  H5Sclose(filespace_blocks_nx3v);
  for (int n = 0; n < num_datasets; ++n)
    H5Sclose(filespaces_vars_blocks_nx3_nx2_nx1[n]);
  delete[] filespaces_vars_blocks_nx3_nx2_nx1;

  // Close datasets
  H5Dclose(dataset_levels);
  H5Dclose(dataset_locations);
  H5Dclose(dataset_x1f);
  H5Dclose(dataset_x2f);
  H5Dclose(dataset_x3f);
  H5Dclose(dataset_x1v);
  H5Dclose(dataset_x2v);
  H5Dclose(dataset_x3v);
  for (int n = 0; n < num_datasets; ++n)
    H5Dclose(datasets_celldata[n]);
  delete[] datasets_celldata;

  // Close .athdf file
  H5Fclose(file);

  // Write .athdf.xdmf file
  int write_xdmf = pin->GetOrAddInteger(output_params.block_name, "xdmf", 1);
  if (Globals::my_rank == 0 && write_xdmf != 0)
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
  delete[] x1v_mesh;
  delete[] x2v_mesh;
  delete[] x3v_mesh;
  for (int n = 0; n < num_datasets; ++n)
    delete[] data_buffers[n];
  delete[] data_buffers;
  delete[] active_flags;

  // Reset parameters for next time file is written
  output_params.file_number++;
  output_params.next_time += output_params.dt;
  pin->SetInteger(output_params.block_name, "file_number", output_params.file_number);
  pin->SetReal(output_params.block_name, "next_time", output_params.next_time);
}

//----------------------------------------------------------------------------------------
//! \fn void ATHDF5Output::MakeXDMF()
//! \brief Function for writing auxiliary XDMF metadata file
//!
//! Inputs: (none)
//! Outputs: (none)
//! Notes:
//!   writes .athdf.xdmf file for describing .athdf file
//!   should only be called by single process
//!   file size scales proportional to total number of MeshBlocks
//!   for many small MeshBlocks, this can take most of the output writing time

void ATHDF5Output::MakeXDMF() {
  std::string filename_aux(filename);
  filename_aux.append(".xdmf");
  std::ofstream xdmf(filename_aux.c_str());

  // Write header
  xdmf << "<?xml version=\"1.0\" ?>\n";
  xdmf << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
  xdmf << "<Xdmf Version=\"2.0\">\n";
  xdmf << "<Information Name=\"TimeVaryingMetaData\" Value=\"True\"/>\n";
  xdmf << "<Domain>\n";
  xdmf << "<Grid Name=\"Mesh\" GridType=\"Collection\">\n";
  xdmf << " <Time Value=\"" << code_time << "\"/>\n";

  // Go through all MeshBlocks
  for (int n_block = 0; n_block < num_blocks_global; ++n_block) {
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
    if (nx3 > 1) {
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
    for (int n_dataset = 0; n_dataset < num_datasets; ++n_dataset) {
      for (int n_variable = 0; n_variable < num_variables[n_dataset]; ++n_variable) {
        xdmf << "    <Attribute Name=\"" << variable_names[n_quantity++]
             << "\" Center=\"Cell\">\n";
        if (nx3 > 1) {
          xdmf << "      <DataItem ItemType=\"HyperSlab\" Dimensions=\"" << nx3 << " "
               << nx2 << " " << nx1 << "\">\n";
          xdmf << "        <DataItem Dimensions=\"3 5\" NumberType=\"Int\"> "
               << n_variable << " " << n_block << " 0 0 0 1 1 1 1 1 1 1 " << nx3 << " "
               << nx2 << " " << nx1 << " </DataItem>\n";
        } else {
          xdmf << "      <DataItem ItemType=\"HyperSlab\" Dimensions=\"" << nx2 << " "
               << nx1 << "\">\n";
          xdmf << "        <DataItem Dimensions=\"3 5\" NumberType=\"Int\"> "
               << n_variable << " " << n_block << " 0 0 0 1 1 1 1 1 1 1 1 " << nx2 << " "
               << nx1 << " </DataItem>\n";
        }
        xdmf << "        <DataItem Dimensions=\"" << num_variables[n_dataset] << " "
             << num_blocks_global << " " << nx3 << " " << nx2 << " " << nx1
             << "\" Format=\"HDF\"> " << filename << ":/" << dataset_names[n_dataset]
             << " </DataItem>\n";
        xdmf << "      </DataItem>\n";
        xdmf << "    </Attribute>\n";
      }
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
