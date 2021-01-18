#ifndef INPUTS_ASCII_TABLE_READER_HPP_
#define INPUTS_ASCII_TABLE_READER_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file ascii_table_reader.hpp
//! \brief Declares ASCII table reader functions

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"             // Real
#include "../athena_arrays.hpp"      // AthenaArray
#include "../utils/interp_table.hpp" // InterpTable2D

//!Declarations
void ASCIITableLoader(const char *filename, InterpTable2D &table,
                      AthenaArray<Real>* pratios=nullptr);

#endif // INPUTS_ASCII_TABLE_READER_HPP_
