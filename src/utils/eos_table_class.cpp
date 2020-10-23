//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file eos_table_class.cpp
//! \brief Implements class EosTable for an EOS lookup table
//========================================================================================

// C headers

// C++ headers
#include <cmath>   // sqrt()
#include <fstream>
#include <iostream> // ifstream
#include <sstream>
#include <stdexcept> // std::invalid_argument
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../field/field.hpp"
#include "../inputs/ascii_table_reader.hpp"
#include "../inputs/hdf5_reader.hpp"
#include "../parameter_input.hpp"
#include "interp_table.hpp"

// Order of datafields for HDF5 EOS tables
const char *var_names[] = {"p/e(e/rho,rho)", "e/p(p/rho,rho)", "asq*rho/p(p/rho,rho)"};

//----------------------------------------------------------------------------------------
//! \fn void ReadBinaryTable(std::string fn, EosTable *peos_table)
//! \brief Read data from binary EOS table and initialize interpolated table.

void ReadBinaryTable(std::string fn, EosTable *peos_table) {
  std::ifstream eos_file(fn.c_str(), std::ios::binary);
  if (eos_file.is_open()) {
    eos_file.seekg(0, std::ios::beg);
    eos_file.read(reinterpret_cast<char*>(&peos_table->nVar), sizeof(peos_table->nVar));
    eos_file.read(reinterpret_cast<char*>(&peos_table->nEgas), sizeof(peos_table->nEgas));
    eos_file.read(reinterpret_cast<char*>(&peos_table->nRho), sizeof(peos_table->nRho));
    eos_file.read(reinterpret_cast<char*>(&peos_table->logEgasMin),
                  sizeof(peos_table->logEgasMin));
    eos_file.read(reinterpret_cast<char*>(&peos_table->logEgasMax),
                  sizeof(peos_table->logEgasMax));
    eos_file.read(reinterpret_cast<char*>(&peos_table->logRhoMin),
                  sizeof(peos_table->logRhoMin));
    eos_file.read(reinterpret_cast<char*>(&peos_table->logRhoMax),
                  sizeof(peos_table->logRhoMax));
    peos_table->EosRatios.NewAthenaArray(peos_table->nVar);
    eos_file.read(reinterpret_cast<char*>(peos_table->EosRatios.data()),
                  peos_table->nVar * sizeof(peos_table->logRhoMin));
    peos_table->table.SetSize(peos_table->nVar, peos_table->nEgas, peos_table->nRho);
    peos_table->table.SetX1lim(peos_table->logRhoMin, peos_table->logRhoMax);
    peos_table->table.SetX2lim(peos_table->logEgasMin, peos_table->logEgasMax);
    eos_file.read(reinterpret_cast<char*>(peos_table->table.data.data()),
                  peos_table->nVar * peos_table->nRho * peos_table->nEgas
                  * sizeof(peos_table->logRhoMin));
    eos_file.close();
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in EosTable::EosTable, ReadBinaryTable" << std::endl
        << "Unable to open eos table: " << fn << std::endl;
    ATHENA_ERROR(msg);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void ReadBinaryTable(std::string fn, EosTable *peos_table, ParameterInput *pin)
//! \brief Read data from HDF5 EOS table and initialize interpolated table.

void ReadHDF5Table(std::string fn, EosTable *peos_table, ParameterInput *pin) {
#ifndef HDF5OUTPUT
  {
    std::stringstream msg;
    msg << "### FATAL ERROR in EosTable::EosTable, ReadHDF5Table" << std::endl
        << "HDF5 EOS table specified, but HDF5 flag is not enabled."  << std::endl;
    ATHENA_ERROR(msg);
  }
#endif
  bool read_ratios = pin->GetOrAddBoolean("hydro", "eos_read_ratios", true);
  std::string dens_lim_field =
      pin->GetOrAddString("hydro", "EOS_dens_lim_field", "LogDensLim");
  std::string espec_lim_field =
      pin->GetOrAddString("hydro", "EOS_espec_lim_field", "LogEspecLim");
  HDF5TableLoader(fn.c_str(), &peos_table->table, 3, var_names,
                  espec_lim_field.c_str(), dens_lim_field.c_str());
  peos_table->table.GetSize(peos_table->nVar, peos_table->nEgas, peos_table->nRho);
  peos_table->table.GetX2lim(peos_table->logEgasMin, peos_table->logEgasMax);
  peos_table->table.GetX1lim(peos_table->logRhoMin, peos_table->logRhoMax);
  peos_table->EosRatios.NewAthenaArray(peos_table->nVar);
  if (read_ratios) {
    std::string ratio_field=pin->GetOrAddString("hydro", "EOS_ratio_field", "ratios");
    int zero[] = {0};
    int pnVar[] = {peos_table->nVar};
    HDF5ReadRealArray(fn.c_str(), ratio_field.c_str(), 1, zero, pnVar,
                      1, zero, pnVar, peos_table->EosRatios);
    if (peos_table->EosRatios(0) <= 0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in EosTable::EosTable, ReadHDF5Table" << std::endl
          << "Invalid ratio. " << fn.c_str() << ", " << ratio_field << ", "
          << peos_table->EosRatios(0) << std::endl;
      ATHENA_ERROR(msg);
    }
  } else {
    for (int i=0; i<peos_table->nVar; ++i) peos_table->EosRatios(i) = 1.0;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void ReadAsciiTable(std::string fn, EosTable *peos_table, ParameterInput *pin)
//! \brief Read data from HDF5 EOS table and initialize interpolated table.

void ReadAsciiTable(std::string fn, EosTable *peos_table, ParameterInput *pin) {
  bool read_ratios = pin->GetOrAddBoolean("hydro", "eos_read_ratios", true);
  AthenaArray<Real> *pratios = nullptr;
  if (read_ratios) pratios = &peos_table->EosRatios;
  // If read_ratios then EosRatios.NewAthenaArray is called in ASCIITableLoader
  ASCIITableLoader(fn.c_str(), peos_table->table, pratios);
  peos_table->table.GetSize(peos_table->nVar, peos_table->nEgas, peos_table->nRho);
  peos_table->table.GetX2lim(peos_table->logEgasMin, peos_table->logEgasMax);
  peos_table->table.GetX1lim(peos_table->logRhoMin, peos_table->logRhoMax);
  if (!read_ratios) {
    peos_table->EosRatios.NewAthenaArray(peos_table->nVar);
    for (int i=0; i<peos_table->nVar; ++i) peos_table->EosRatios(i) = 1.0;
  }
}

// ctor
EosTable::EosTable(ParameterInput *pin) :
    table(), logRhoMin(), logRhoMax(),
    rhoUnit(pin->GetOrAddReal("hydro", "eos_rho_unit", 1.0)),
    eUnit(pin->GetOrAddReal("hydro", "eos_egas_unit", 1.0)),
    hUnit(eUnit/rhoUnit) {
  std::string eos_fn, eos_file_type;
  eos_fn = pin->GetString("hydro", "eos_file_name");
  eos_file_type = pin->GetOrAddString("hydro", "eos_file_type", "auto");

  if (eos_file_type.compare("auto") == 0) {
    std::string ext = eos_fn.substr(eos_fn.find_last_of(".") + 1);
    if (ext.compare("data")*ext.compare("bin") == 0) {
      eos_file_type.assign("binary");
    } else if (ext.compare("hdf5") == 0) {
      eos_file_type.assign("hdf5");
    } else if (ext.compare("tab")*ext.compare("txt")*ext.compare("ascii") == 0) {
      eos_file_type.assign("ascii");
    }
  }
  if (eos_file_type.compare("binary") == 0) { //Raw binary
    ReadBinaryTable(eos_fn, this);
  } else if (eos_file_type.compare("hdf5") == 0) { // HDF5 table
    ReadHDF5Table(eos_fn, this, pin);
  } else if (eos_file_type.compare("ascii") == 0) { // ASCII/text table
    ReadAsciiTable(eos_fn, this, pin);
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in EosTable::EosTable" << std::endl
        << "EOS table of type '" << eos_file_type << "' not recognized."  << std::endl
        << "Options are 'ascii', 'binary', and 'hdf5'." << std::endl;
    ATHENA_ERROR(msg);
  }
}
