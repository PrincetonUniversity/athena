//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hdf5_reader.cpp
//! \brief Implements HDF5 reader functions

// C headers

// C++ headers
#include <iostream>   // cout
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // string

// Athena++ headers
#include "../athena.hpp"         // Real
#include "../athena_arrays.hpp"  // AthenaArray
#include "../defs.hpp"           // SINGLE_PRECISION_ENABLED
#include "hdf5_reader.hpp"

// Only proceed if HDF5 enabled
#ifdef HDF5OUTPUT

// External library headers
#include <hdf5.h>  // H5[F|P|S|T]_*, H5[D|F|P|S]*(), hid_t
#ifdef MPI_PARALLEL
#include <mpi.h>  // MPI_COMM_WORLD, MPI_INFO_NULL
#endif

// Determine floating-point precision (in memory, not file)
#if SINGLE_PRECISION_ENABLED
#define H5T_REAL H5T_NATIVE_FLOAT
#else
#define H5T_REAL H5T_NATIVE_DOUBLE
#endif

//----------------------------------------------------------------------------------------
//! \fn void HDF5ReadArray(const char *filename, const char *dataset_name, int rank_file,
//!     const int *start_file, const int *count_file, int rank_mem, const int *start_mem,
//!     const int *count_mem, AthenaArray<Real> &array, bool collective=false,
//!     bool noop=false)
//! \brief Read a single dataset from an HDF5 file into a pre-allocated array.

void HDF5ReadRealArray(const char *filename, const char *dataset_name, int rank_file,
                       const int *start_file, const int *count_file, int rank_mem,
                       const int *start_mem, const int *count_mem,
                       AthenaArray<Real> &array,
                       bool collective, bool noop) {
  // Check that user is not trying to exceed limits of HDF5 array or AthenaArray
  // dimensionality
  if (rank_file > MAX_RANK_FILE) {
    std::stringstream msg;
    msg << "### FATAL ERROR\nAttempting to read HDF5 array of ndim= " << rank_file
        << "\nExceeding MAX_RANK_FILE=" << MAX_RANK_FILE << std::endl;
    ATHENA_ERROR(msg);
  }
  if (rank_mem > MAX_RANK_MEM) {
    std::stringstream msg;
    msg << "### FATAL ERROR\nAttempting to read HDF5 array of ndim= " << rank_mem
        << "\nExceeding MAX_RANK_MEM=" << MAX_RANK_MEM << std::endl;
    ATHENA_ERROR(msg);
  }

  // Cast selection arrays to appropriate types
  hsize_t start_file_hid[MAX_RANK_FILE];
  hsize_t count_file_hid[MAX_RANK_FILE];
  for (int n = 0; n < rank_file; ++n) {
    start_file_hid[n] = start_file[n];
    count_file_hid[n] = count_file[n];
  }
  hsize_t start_mem_hid[MAX_RANK_MEM];
  hsize_t count_mem_hid[MAX_RANK_MEM];
  for (int n = 0; n < rank_mem; ++n) {
    start_mem_hid[n] = start_mem[n];
    count_mem_hid[n] = count_mem[n];
  }

  // Determine AthenaArray dimensions
  hsize_t dims_mem_base[5];
  dims_mem_base[0] = array.GetDim5();
  dims_mem_base[1] = array.GetDim4();
  dims_mem_base[2] = array.GetDim3();
  dims_mem_base[3] = array.GetDim2();
  dims_mem_base[4] = array.GetDim1();
  hsize_t *dims_mem = dims_mem_base + 5 - rank_mem;

  // Open data file
  hid_t property_list_file = H5Pcreate(H5P_FILE_ACCESS);
#ifdef MPI_PARALLEL
  {
    if (collective) {
      H5Pset_fapl_mpio(property_list_file, MPI_COMM_WORLD, MPI_INFO_NULL);
    }
  }
#endif
  hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, property_list_file);
  H5Pclose(property_list_file);
  if (file < 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR\nCould not open " << filename << std::endl;
    ATHENA_ERROR(msg);
  }
  hid_t property_list_transfer = H5Pcreate(H5P_DATASET_XFER);
#ifdef MPI_PARALLEL
  {
    if (collective) {
      H5Pset_dxpl_mpio(property_list_transfer, H5FD_MPIO_COLLECTIVE);
    }
  }
#endif

  // Read dataset into array
  hid_t dataset = H5Dopen(file, dataset_name, H5P_DEFAULT);
  hid_t dataspace_file = H5Dget_space(dataset);
  if (noop) {
    H5Sselect_none(dataspace_file);
  }
  H5Sselect_hyperslab(dataspace_file, H5S_SELECT_SET, start_file_hid, NULL,
                      count_file_hid, NULL);
  hid_t dataspace_mem = H5Screate_simple(rank_mem, dims_mem, NULL);
  if (noop) {
    H5Sselect_none(dataspace_mem);
  }
  H5Sselect_hyperslab(dataspace_mem, H5S_SELECT_SET, start_mem_hid, NULL, count_mem_hid,
                      NULL);
  H5Dread(dataset, H5T_REAL, dataspace_mem, dataspace_file, property_list_transfer,
          array.data());
  H5Dclose(dataset);
  H5Sclose(dataspace_file);
  H5Sclose(dataspace_mem);

  // Close data file
  H5Pclose(property_list_transfer);
  H5Fclose(file);
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void HDF5TableLoader(const char *filename, InterpTable2D* ptable, const int nvar,
//!                    const char **var_names, char *x2lim_name, char *x1lim_name) {
//! \brief Reads datasets from an HDF5 file into a InterpTable2D.

void HDF5TableLoader(const char *filename, InterpTable2D* ptable, const int nvar,
                     const char **var_names, const char *x2lim_name,
                     const char *x1lim_name) {
  hsize_t dims[2];
  int tmp[2];
  int count_file[2];
  hid_t dataset, dspace;
  hid_t property_list_file = H5Pcreate(H5P_FILE_ACCESS);
  hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, property_list_file);
  for (int i = 0; i < nvar; ++i) {
    dataset = H5Dopen(file, var_names[i], H5P_DEFAULT);
    dspace = H5Dget_space(dataset);
    int ndims = H5Sget_simple_extent_ndims(dspace);
    if (ndims != 2) {
      std::stringstream msg;
      msg << "### FATAL ERROR in HDF5TableLoader" << std::endl
          << "Rank of data field '" << var_names[i] << "' in file '" << filename
          << "' must be 2. Rank is " << ndims << "." << std::endl;
      ATHENA_ERROR(msg);
    }
    H5Sget_simple_extent_dims(dspace, dims, NULL);
    tmp[0] = static_cast<int>(dims[0]);
    tmp[1] = static_cast<int>(dims[1]);
    if (i == 0) {
      count_file[0] = tmp[0];
      count_file[1] = tmp[1];
    } else if (count_file[0]!=tmp[0] || count_file[1]!=tmp[1]) {
      std::stringstream msg;
      msg << "### FATAL ERROR in HDF5TableLoader" << std::endl
          << "Inconsistent data field shape in file '" << filename << "'." << std::endl;
      ATHENA_ERROR(msg);
    }
  }
  H5Fclose(file);
  ptable->SetSize(nvar, count_file[0], count_file[1]);
  int start_file[2];
  start_file[0] = 0;
  start_file[1] = 0;
  int start_mem[3];
  start_mem[1] = 0;
  start_mem[2] = 0;
  int count_mem[3];
  count_mem[0] = 1;
  count_mem[1] = count_file[0];
  count_mem[2] = count_file[1];
  for (int i = 0; i < nvar; ++i) {
    start_mem[0] = i;
    HDF5ReadRealArray(filename, var_names[i], 2, start_file, count_file,
                      3, start_mem, count_mem, ptable->data);
  }
  if (x2lim_name) {
    AthenaArray<Real> lim;
    lim.NewAthenaArray(2);
    int zero[] = {0};
    int two[] = {2};
    HDF5ReadRealArray(filename, x2lim_name, 1, zero, two, 1, zero, two, lim);
    ptable->SetX2lim(lim(0), lim(1));
  }
  if (x1lim_name) {
    AthenaArray<Real> lim;
    lim.NewAthenaArray(2);
    int zero[] = {0};
    int two[] = {2};
    HDF5ReadRealArray(filename, x1lim_name, 1, zero, two, 1, zero, two, lim);
    ptable->SetX1lim(lim(0), lim(1));
  }
  return;
}
#endif  // HDF5OUTPUT
