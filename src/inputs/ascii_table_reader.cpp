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

// Declarations
void ASCIITableLoader(const char *filename, InterpTable2D* ptable, AthenaArray<Real>* pratios) {
  std::ifstream file(filename, std::ios::in);
  std::string line;

  std::getline(file, line);
  while ((line[0] = '#') && std::getline(file, line));
  int nvar, nx2, nx1;
  std::stringstream stream;
  stream.str(line);
  stream >> nvar;
  stream >> nx2;
  stream >> nx1;
  ptable->SetSize(nvar, nx2, nx1);

  std::getline(file, line);
  while ((line[0] = '#') && std::getline(file, line));
  Real min_, max_;
  stream.str(line);
  stream >> min_;
  stream >> max_;
  ptable->SetX2lim(min_, max_);

  std::getline(file, line);
  while ((line[0] = '#') && std::getline(file, line));
  stream.str(line);
  stream >> min_;
  stream >> max_;
  ptable->SetX1lim(min_, max_);

  if (pratios) {
    std::getline(file, line);
    while ((line[0] = '#') && std::getline(file, line));
    stream.str(line);
    pratios->NewAthenaArray(nvar);
    for (int i = 0; i < nvar; ++i) {
      stream >> (*pratios)(i);
    }
  }

  for (int row = 0; row < nx2; ++row) {
    std::getline(file, line);
    while ((line[0] = '#') && std::getline(file, line));
    std::stringstream stream(line);
    for (int col = 0; col < nx1; ++col) {
      stream >> ptable->data(row, col);
    }
  }
  return;
}

#endif
