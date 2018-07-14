#ifndef INPUTS_HDF5_READER_HPP_
#define INPUTS_HDF5_READER_HPP_

//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hdf5_reader.hpp
//  \brief Declares HDF5 reader functions

// Declarations
void HDF5ReadRealArray(const char *filename, const char *dataset_name, int rank,
    const int *dims, const int *offset, AthenaArray<Real> &array);

#endif  // INPUTS_HDF5_READER_HPP_
