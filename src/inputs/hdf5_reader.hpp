#ifndef INPUTS_HDF5_READER_HPP_
#define INPUTS_HDF5_READER_HPP_

//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hdf5_reader.hpp
//! \brief Declares HDF5 reader functions

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"             // Real
#include "../athena_arrays.hpp"      // AthenaArray
#include "../utils/interp_table.hpp" // InterpTable2D

// upper limit on the dimensionality of the resulting AthenaArray
#define MAX_RANK_MEM 5
// upper limit on the dimensionality of the external dataset
#define MAX_RANK_FILE 5

// Declarations
void HDF5ReadRealArray(const char *filename, const char *dataset_name, int rank_file,
                       const int *start_file, const int *count_file, int rank_mem,
                       const int *start_mem, const int *count_mem,
                       AthenaArray<Real> &array,
                       bool collective=false, bool noop=false);

void HDF5TableLoader(const char *filename, InterpTable2D* ptable, const int nvar,
                     const char **var_names,
                     const char *x2lim_name=NULL,
                     const char *x1lim_name=NULL);

#endif  // INPUTS_HDF5_READER_HPP_
