//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//
//! \file hdf5_reader.cpp
//  \brief prototypes of functions and class definitions for utils/*.cpp files

// C++ headers
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // string

// Athena++ headers
#include "../athena.hpp"         // Real
#include "../athena_arrays.hpp"  // AthenaArray
#include "../defs.hpp"           // SINGLE_PRECISION_ENABLED

// Only proceed if HDF5 enabled
#ifdef HDF5OUTPUT

// External library headers
#include <hdf5.h>  // H5[F|P|T]_*, H5[D|F|P|S]*(), hid_t
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
//! \fn void HDF5ReadArray(const char *filename, const char *dataset_name, int rank,
//      const hsize_t *dims, AthenaArray<Real> &array)
//  \brief Read a single dataset from an HDF5 file into a pre-allocated array.

void HDF5ReadRealArray(const char *filename, const char *dataset_name, int rank,
    const hsize_t *dims, AthenaArray<Real> &array) {

  // Open data file
  hid_t property_list_file = H5Pcreate(H5P_FILE_ACCESS);
  #ifdef MPI_PARALLEL
    H5Pset_fapl_mpio(property_list_file, MPI_COMM_WORLD, MPI_INFO_NULL);
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
    H5Pset_dxpl_mpio(property_list_transfer, H5FD_MPIO_COLLECTIVE);
  #endif

  // Read dataset into array
  hid_t dataset = H5Dopen(file, dataset_name, H5P_DEFAULT);
  hid_t dataspace = H5Screate_simple(rank, dims, NULL);
  H5Dread(dataset, H5T_REAL, dataspace, dataspace, property_list_transfer, array.data());
  H5Dclose(dataset);
  H5Sclose(dataspace);

  // Close data file
  H5Pclose(property_list_transfer);
  H5Fclose(file);
  return;
}

#endif  // HDF5OUTPUT
