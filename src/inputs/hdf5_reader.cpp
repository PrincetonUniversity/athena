//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hdf5_reader.cpp
//  \brief Implements HDF5 reader functions

// C++ headers
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // string

// Athena++ headers
#include "hdf5_reader.hpp"
#include "../athena.hpp"         // Real
#include "../athena_arrays.hpp"  // AthenaArray
#include "../defs.hpp"           // SINGLE_PRECISION_ENABLED

// Only proceed if HDF5 enabled
#ifdef HDF5OUTPUT

// External library headers
#include <hdf5.h>  // H5[F|P|S|T]_*, H5[D|F|P|S]*(), hid_t
#ifdef MPI_PARALLEL
  #include <mpi.h>  // MPI_COMM_WORLD, MPI_INFO_NULL
#endif

// Determine floating point precision (in memory, not file)
#if SINGLE_PRECISION_ENABLED
  #define H5T_REAL H5T_NATIVE_FLOAT
#else
  #define H5T_REAL H5T_NATIVE_DOUBLE
#endif

//----------------------------------------------------------------------------------------
//! \fn void HDF5ReadArray(const char *filename, const char *dataset_name, int rank_file,
//      const int *start_file, const int *count_file, int rank_mem, const int *start_mem,
//      const int *count_mem, AthenaArray<Real> &array, bool collective=false,
//      bool noop=false)
//  \brief Read a single dataset from an HDF5 file into a pre-allocated array.

void HDF5ReadRealArray(const char *filename, const char *dataset_name, int rank_file,
    const int *start_file, const int *count_file, int rank_mem, const int *start_mem,
    const int *count_mem, AthenaArray<Real> &array, bool collective, bool noop) {

  // Check that user is not trying to exceed limits of HDF5 array or AthenaArray
  // dimensionality
  if (rank_file > MAX_RANK_FILE) {
    std::stringstream message;
    message << "### FATAL ERROR\nAttempting to read HDF5 array of ndim= " << rank_file
            << "\nExceeding MAX_RANK_FILE=" << MAX_RANK_FILE << std::endl;
    throw std::runtime_error(message.str().c_str());
  }
  if (rank_mem > MAX_RANK_MEM) {
    std::stringstream message;
    message << "### FATAL ERROR\nAttempting to read HDF5 array of ndim= " << rank_mem
            << "\nExceeding MAX_RANK_MEM=" << MAX_RANK_MEM << std::endl;
    throw std::runtime_error(message.str().c_str());
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
    std::stringstream message;
    message << "### FATAL ERROR\nCould not open " << filename << "\n";
    throw std::runtime_error(message.str().c_str());
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

#endif  // HDF5OUTPUT
