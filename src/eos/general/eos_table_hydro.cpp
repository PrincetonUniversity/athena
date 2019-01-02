//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file eos_table.cpp
//  \brief implements functions in class EquationOfState for an EOS lookup table
//======================================================================================

// C++ headers
#include <cmath>   // sqrt()
#include <cfloat>  // FLT_MIN
#include <iostream> // ifstream
#include <fstream>
#include <sstream>
#include <stdexcept> // std::invalid_argument

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../parameter_input.hpp"
#include "../../field/field.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../inputs/hdf5_reader.hpp"
#include "../../inputs/ascii_table_reader.hpp"

// this class header
#include "../eos.hpp"

#if EOS_TABLE_ENABLED
//#define EOSDEBUG0
//#define EOSDEBUG1

const char *var_names[] = {"p/e(e/rho,rho)", "e/p(p/rho,rho)", "asq*rho/p(p/rho,rho)", "asq*rho/h(h/rho,rho)"};

Real EquationOfState::GetEosData(int kOut, Real var, Real rho) {
  Real x1 = std::log10(rho * ts.rhoUnit);
  Real x2 = std::log10(var * ts.EosRatios(kOut) * ts.eUnit) - x1;
  return std::pow((Real)10, ptable_->interpolate(kOut, x2, x1));
}

Real EquationOfState::SimplePres(Real rho, Real egas) {
  return GetEosData(0, egas, rho) * egas;
}

Real EquationOfState::SimpleEgas(Real rho, Real pres) {
  return GetEosData(1, pres, rho) * pres;
}

Real EquationOfState::SimpleAsq(Real rho, Real pres) {
  return GetEosData(2, pres, rho) * pres / rho;
}

Real EquationOfState::RiemannAsq(Real rho, Real hint) {
  return std::pow((Real)10, ptable_->interpolate(3, std::log10(hint*ts.hUnit), std::log10(rho*ts.rhoUnit))) * hint;
}

void ReadBinaryTable(std::string fn, InterpTable2D *ptable, TableSize &ts) {
  std::ifstream eos_table(fn.c_str(), std::ios::binary);
  Real egasOverPres;
  if (eos_table.is_open())
  {
    eos_table.seekg(0, std::ios::beg);
    eos_table.read((char*)&ts.nEgas, sizeof(ts.nEgas));
    eos_table.read((char*)&ts.logEgasMin, sizeof(ts.logEgasMin));
    eos_table.read((char*)&ts.logEgasMax, sizeof(ts.logEgasMax));
    eos_table.read((char*)&ts.nRho, sizeof(ts.nRho));
    eos_table.read((char*)&ts.logRhoMin, sizeof(ts.logRhoMin));
    eos_table.read((char*)&ts.logRhoMax, sizeof(ts.logRhoMax));
    eos_table.read((char*)&ts.nVar, sizeof(ts.nVar));
    ts.EosRatios.NewAthenaArray(ts.nVar);
    //eos_table.read((char*)&egasOverPres, sizeof(egasOverPres));
    eos_table.read((char*)ts.EosRatios.data(), ts.nVar * sizeof(ts.logRhoMin));
    ptable->SetSize(ts.nVar, ts.nEgas, ts.nRho);
    ptable->SetX1lim(ts.logRhoMin, ts.logRhoMax);
    ptable->SetX2lim(ts.logEgasMin, ts.logEgasMax);
    eos_table.read((char*)ptable->data.data(), ts.nVar * ts.nRho * ts.nEgas * sizeof(ts.logRhoMin));
    eos_table.close();
  }
  else throw std::invalid_argument("Unable to open eos table: " + fn);
}

void EquationOfState::PrepEOS(ParameterInput *pin) {
  std::string EOS_fn, EOS_file_type, dens_lim_field, espec_lim_field;
  EOS_fn = pin->GetString("hydro", "EOS_file_name");
  EOS_file_type = pin->GetString("hydro", "EOS_file_type");
  bool read_ratios = pin->GetOrAddBoolean("hydro", "EOS_read_ratios", true);
  dens_lim_field = pin->GetOrAddString("hydro", "EOS_dens_lim_field", "LogDensLim");
  espec_lim_field = pin->GetOrAddString("hydro", "EOS_espec_lim_field", "LogEspecLim");
  ptable_ = new InterpTable2D();

  if (EOS_file_type.compare("binary") == 0) {
    ReadBinaryTable(EOS_fn, ptable_, ts);
  } else if (EOS_file_type.compare("hdf5") == 0) {
#ifndef HDF5OUTPUT
    std::stringstream msg;
    msg << "### FATAL ERROR in EquationOfState::PrepEOS" << std::endl
        << "HDF5 EOS table specified, but HDF5 flag is not enabled."  << std::endl;
    throw std::runtime_error(msg.str().c_str());
#endif
    HDF5TableLoader(EOS_fn.c_str(), ptable_, 4, var_names, espec_lim_field.c_str(), dens_lim_field.c_str());
    ptable_->GetSize(ts.nVar, ts.nEgas, ts.nRho);
    ptable_->GetX2lim(ts.logEgasMin, ts.logEgasMax);
    ptable_->GetX1lim(ts.logRhoMin, ts.logRhoMax);
    ts.EosRatios.NewAthenaArray(ts.nVar);
    if (read_ratios) {
      std::string ratio_field = pin->GetOrAddString("hydro", "EOS_ratio_field", "ratios");
      int zero[] = {0};
      int pnVar[] = {ts.nVar};
      HDF5ReadRealArray(EOS_fn.c_str(), ratio_field.c_str(), 1, zero, pnVar, 1, zero, pnVar, ts.EosRatios);
      if (ts.EosRatios(0) <= 0){
        std::stringstream msg;
        msg << "### FATAL ERROR in EquationOfState::PrepEOS" << std::endl
            << "Invalid ratio. " << EOS_fn.c_str() << ", " << ratio_field << ", " << ts.EosRatios(0) << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      //std::cout << "Ratios: " << ratio_array(0) << ", " << ratio_array(1) << ", " << ratio_array(2) << '\n';
    } else {
      for (int i=0; i<ts.nVar; ++i) ts.EosRatios(i) = 1.0;
    }
  } else if (EOS_file_type.compare("ascii") == 0) {
    AthenaArray<Real> *pratios = NULL;
    if (read_ratios) pratios = &ts.EosRatios;
    ASCIITableLoader(EOS_fn.c_str(), ptable_, pratios);
    ptable_->GetSize(ts.nVar, ts.nEgas, ts.nRho);
    ptable_->GetX2lim(ts.logEgasMin, ts.logEgasMax);
    ptable_->GetX1lim(ts.logRhoMin, ts.logRhoMax);
    if (!read_ratios) for (int i=0; i<ts.nVar; ++i) ts.EosRatios(i) = 1.0;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in EquationOfState::PrepEOS" << std::endl
        << "EOS table of type '" << EOS_file_type << "' not recognized."  << std::endl
        << "Options are 'ascii', 'binary', and 'hdf5'." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }


  ts.rhoUnit = pin->GetOrAddReal("hydro", "EosRhoUnit", 1.0);
  ts.eUnit = pin->GetOrAddReal("hydro", "EosEgasUnit", 1.0);
  ts.hUnit = ts.eUnit/ts.rhoUnit;

#ifdef EOSDEBUG1
  std::cout << "Eos table: " << ts.nVar << ", " << ts.nRho << ", " << ts.nEgas << "\n";
  std::cout << "logRhoMin, logRhoMax: " << ts.logRhoMin << ", " << ts.logRhoMax << "\n";
  std::cout << "logEgasMin, logEgasMax: " << ts.logEgasMin << ", " << ts.logEgasMax << "\n";
  std::cout << "Ratios: ";
  for (int i=0;i<ts.nVar;i++) std::cout << ts.EosRatios(i) << ", ";
  std::cout << "\n";
  for (int i=0;i<ts.nVar;i++) {
    std::cout << "var = " << i << "\n";
    for (int j=0;j<ts.nRho;j++){
      for (int k=0;k<ts.nEgas;k++){
        std::cout << std::pow((Real) 10, ptable_->data(i,j,k)) << " ";
      }
      std::cout << "\n";
    }
    std::cout << "\n";
  }
#endif

#ifdef EOSDEBUG0
  std::cout << "prepEOS: " << EOS_fn << ", " << ts.nVar << ", " << ts.nEgas <<  ", " << ts.nRho << "\n";
  EosTestLoop();
#endif
}

void EquationOfState::CleanEOS()
{
  ts.EosRatios.DeleteAthenaArray();
  ptable_->~InterpTable2D();
}

void EquationOfState::EosTestRhoEgas(Real rho, Real egas, AthenaArray<Real> &data)
{
  Real idn = 1./rho;
  Real ien = 1./egas;
  data(0) = GetEosData(0, egas, rho) * egas;
  data(1) = (data(0) + egas) * idn;
  data(2) = GetEosData(2, data(0), rho) * data(0) * idn;
  data(3) = (egas -  GetEosData(1, data(0), rho) * data(0)) * ien;
  data(4) = 1.0 - RiemannAsq(rho, data(1)) / data(2);
}

void EquationOfState::EosTestLoop()
{
  Real rho = 0;
  Real egas;
  AthenaArray<Real> data;
  data.NewAthenaArray(5);
  std::cout << "logRhoMin, logRhoMax: " << ts.logRhoMin << ", " << ts.logRhoMax << "\n";
  std::cout << "logEgasMin, logEgasMax: " << ts.logEgasMin << ", " << ts.logEgasMax << "\n";
  std::cout << "EosRatios = ";
  for (int i=0; i<ts.nVar; ++i) std::cout << ts.EosRatios(i) << ", ";
  std::cout << "\n";
  std::cout << "Input rho (g/cc):";
  std::cin >> rho;
  std::cout << "Input egas (erg/cc):";
  std::cin >> egas;
  while (rho >= 0){
    Real Lrho = std::log10(rho*ts.rhoUnit);
    std::cout << "log10[p/e(e/d,d)] = " << ptable_->interpolate(0, std::log10(egas*ts.eUnit)-Lrho, Lrho) << "\n";
    EosTestRhoEgas(rho, egas, data);
    std::cout << "log10[e/p(p/d,d))] = " << ptable_->interpolate(1, std::log10(data(0)*ts.eUnit)-Lrho, Lrho) << "\n";
    std::cout << rho << ", " << egas << "\n";
    std::cout << "P(d, e)    , h(d, e)    , ASq(d, P)  ,PErr    , ASqErr\n";
    for (int i=0;i<5;i++) std::cout << data(i) << ", ";
    std::cout << "\n";
    std::cout << "P, e, Asq, Asq: " << SimplePres(rho, egas) << ", "
                                    << SimpleEgas(rho, data(0)) << ", "
                                    << SimpleAsq(rho, data(0)) << ", "
                                    << RiemannAsq(rho, data(1)) << "\n\n";
    std::cout << "Input rho (g/cc):";
    std::cin >> rho;
    std::cout << "Input egas (erg/cc):";
    std::cin >> egas;
  }
  data.DeleteAthenaArray();
}

#endif
