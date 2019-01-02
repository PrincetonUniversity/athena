#ifndef INPUTS_ASCII_READER_HPP_
#define INPUTS_ASCII_READER_HPP_

//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file ascii_table_reader.cpp
//  \brief Implements ASCII table reader functions

// C++ headers
#include <iostream>   // ifstream
#include <fstream>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // string

// Athena++ headers
#include "../athena.hpp"             // Real
#include "../athena_arrays.hpp"      // AthenaArray
#include "../utils/interp_table.hpp" // InterpTable2D

//----------------------------------------------------------------------------------------
//! \fn void ASCIITableLoader(const char *filename, InterpTable2D* ptable,
//                            AthenaArray<Real>* pratios)
//  \brief Load a table stored in ASCII form and initialize an interpolated table.
//         Fastest index corresponds to column
void ASCIITableLoader(const char *filename, InterpTable2D* ptable,
                      AthenaArray<Real>* pratios) {
  std::ifstream file(filename, std::ios::in);
  std::string line;

  while (std::getline(file, line) && (line[0] == '#')); // skip comments
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
    throw std::runtime_error(msg.str().c_str());
  }
  ptable->SetSize(nvar, nx2, nx1);

  // Read and store x2lim
  Real min_, max_;
  while (std::getline(file, line) && (line[0] == '#')); // skip comments
  stream.str(line);
  stream.clear();
  stream >> min_;
  stream >> max_;
  if (min_>=max_) {
    std::stringstream msg;
    msg << "### FATAL ERROR in ASCIITableLoader" << std::endl
        << "x2min>=x2max." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  ptable->SetX2lim(min_, max_);

  // Read and store x1lim
  while (std::getline(file, line) && (line[0] == '#')); // skip comments
  stream.str(line);
  stream.clear();
  stream >> min_;
  stream >> max_;
  if (min_>=max_) {
    std::stringstream msg;
    msg << "### FATAL ERROR in ASCIITableLoader" << std::endl
        << "x1min>=x1max." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  ptable->SetX1lim(min_, max_);

  // read ratios for each (#=nvar) x2
  if (pratios) {
    while (std::getline(file, line) && (line[0] == '#'));
    stream.str(line);
    stream.clear();
    pratios->NewAthenaArray(nvar);
    for (int i = 0; i < nvar; ++i) {
      stream >> (*pratios)(i);
    }
  }

  // read table data
  for (int row = 0; row < nx2 * nvar; ++row) {
    while (std::getline(file, line) && (line[0] == '#'));
    std::stringstream stream(line);
    for (int col = 0; col < nx1; ++col) {
      stream >> ptable->data(row, col);
    }
  }
  return;
}

#endif
