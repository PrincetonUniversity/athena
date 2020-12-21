//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file ascii_table_reader.cpp
//! \brief Implements ASCII table reader functions

// C headers

// C++ headers
#include <fstream>
#include <iostream>   // ifstream
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // string

// Athena++ headers
#include "../athena.hpp"             // Real
#include "../athena_arrays.hpp"      // AthenaArray
#include "../utils/interp_table.hpp" // InterpTable2D
#include "ascii_table_reader.hpp"

//----------------------------------------------------------------------------------------
//! \fn void ASCIITableLoader(const char *filename, InterpTable2D* ptable,
//!                           AthenaArray<Real>* pratios)
//! \brief Load a table stored in ASCII form and initialize an interpolated table.
//!        Fastest index corresponds to column
void ASCIITableLoader(const char *filename, InterpTable2D &table,
                      AthenaArray<Real>* pratios) {
  std::ifstream file(filename, std::ios::in);
  std::string line;

  while (std::getline(file, line) && (line[0] == '#')) continue; // skip comments
  int nvar, nx2, nx1;
  std::stringstream stream;
  stream.str(line);
  // get data shape
  stream >> nvar;
  stream >> nx2;
  stream >> nx1;
  if (nvar<1 || nx2<2 || nx1<2) {
    std::stringstream msg;
    msg << "### FATAL ERROR in ASCIITableLoader" << std::endl
        << "Invalid shape: (" << nvar << ", " << nx2 << ", " << nx1 << ")" << std::endl;
    ATHENA_ERROR(msg);
  }
  table.SetSize(nvar, nx2, nx1);

  // Read and store x2lim
  Real min_, max_;
  while (std::getline(file, line) && (line[0] == '#')) continue; // skip comments
  stream.str(line);
  stream.clear();
  stream >> min_;
  stream >> max_;
  if (min_>=max_) {
    std::stringstream msg;
    msg << "### FATAL ERROR in ASCIITableLoader" << std::endl
        << "x2min>=x2max." << std::endl;
    ATHENA_ERROR(msg);
  }
  table.SetX2lim(min_, max_);

  // Read and store x1lim
  while (std::getline(file, line) && (line[0] == '#')) continue; // skip comments
  stream.str(line);
  stream.clear();
  stream >> min_;
  stream >> max_;
  if (min_>=max_) {
    std::stringstream msg;
    msg << "### FATAL ERROR in ASCIITableLoader" << std::endl
        << "x1min>=x1max." << std::endl;
    ATHENA_ERROR(msg);
  }
  table.SetX1lim(min_, max_);

  // read ratios for each (#=nvar) x2
  if (pratios != nullptr) {
    while (std::getline(file, line) && (line[0] == '#')) continue;
    stream.str(line);
    stream.clear();
    pratios->NewAthenaArray(nvar);
    for (int i = 0; i < nvar; ++i) {
      stream >> (*pratios)(i);
    }
  }

  // read table data
  for (int row = 0; row < nx2 * nvar; ++row) {
    while (std::getline(file, line) && (line[0] == '#')) continue;
    std::stringstream lstream(line);
    for (int col = 0; col < nx1; ++col) {
      lstream >> table.data(row, col);
    }
  }
  return;
}
