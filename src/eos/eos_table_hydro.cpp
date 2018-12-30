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
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../field/field.hpp"
#include "../coordinates/coordinates.hpp"
#include "../inputs/hdf5_reader.hpp"
#include "../inputs/ascii_table_reader.hpp"

// this class header
#include "eos.hpp"

#if EOS_TABLE_ENABLED
//#define EOSDEBUG0
#define EOSDEBUG1

const char *var_names[] = {"e/p(p/rho,rho)", "asq*rho/p(p/rho,rho)", "T*rho/p(p/rho,rho)",
                           "p/e(e/rho,rho)", "asq*rho/e(e/rho,rho)", "T*rho/e(e/rho,rho)",
                           "e/h(h/rho,rho)", "asq*rho/h(h/rho,rho)", "T*rho/h(h/rho,rho)"};

Real EquationOfState::SimplePres(Real rho, Real egas) {
  return GetEosData(rho, egas, ts.axisEgas, ts.iPres) * egas;
}

Real EquationOfState::SimpleEgas(Real rho, Real pres) {
  return GetEosData(rho, pres, ts.axisPres, ts.iPres) * pres;
}

Real EquationOfState::SimpleAsq(Real rho, Real pres) {
  return GetEosData(rho, pres, ts.axisPres, ts.iASq) * pres / rho;
}

Real EquationOfState::RiemannAsq(Real rho, Real hint) {
  return GetEosData(rho, hint, ts.axisHint, ts.iASq) * hint;
}

void ReadBinaryTable(std::string fn, InterpTable2D *ptable, TableSize &ts) {
  std::ifstream eos_table(fn.c_str(), std::ios::binary);
  if (eos_table.is_open())
  {
    eos_table.seekg(0, std::ios::beg);
    eos_table.read((char*)&ts.nRho, sizeof(ts.nRho));
    eos_table.read((char*)&ts.logRhoMin, sizeof(ts.logRhoMin));
    eos_table.read((char*)&ts.logRhoMax, sizeof(ts.logRhoMax));
    eos_table.read((char*)&ts.nEgas, sizeof(ts.nEgas));
    eos_table.read((char*)&ts.logEgasMin, sizeof(ts.logEgasMin));
    eos_table.read((char*)&ts.logEgasMax, sizeof(ts.logEgasMax));
    eos_table.read((char*)&ts.egasOverPres, sizeof(ts.egasOverPres));
    eos_table.read((char*)&ts.nVar, sizeof(ts.nVar));
    ptable->SetSize(ts.nVar, ts.nEgas, ts.nRho);
    ptable->SetX1lim(ts.logRhoMin, ts.logRhoMax);
    ptable->SetX2lim(ts.logEgasMin, ts.logEgasMax);
    eos_table.read((char*)ptable->data.data(), ts.nVar * ts.nRho * ts.nEgas * sizeof(ts.logRhoMin));
    eos_table.close();
  }
  else throw std::invalid_argument("Unable to open eos table: " + fn);

  ts.EosRatios[0] = 1;
  ts.EosRatios[1] = ts.egasOverPres;
  ts.EosRatios[2] = ts.egasOverPres / (1. + ts.egasOverPres);
}

void EquationOfState::PrepEOS(ParameterInput *pin) {
  std::string EOS_fn, EOS_file_type, dens_lim_field, espec_lim_field;
  EOS_fn = pin->GetString("hydro", "EOS_file_name");
  EOS_file_type = pin->GetString("hydro", "EOS_file_type");
  bool read_ratios = pin->GetOrAddBoolean("hydro", "EOS_read_ratios", true);
  dens_lim_field = pin->GetOrAddString("hydro", "EOS_dens_lim_field", "LogDensLim");
  espec_lim_field = pin->GetOrAddString("hydro", "EOS_espec_lim_field", "LogEspecLim");
  ptable_ = new InterpTable2D();
  ts.egasOverPres = 1;
  ts.EosRatios[0] = 1;
  ts.EosRatios[1] = 1;
  ts.EosRatios[2] = 1;
  ts.iPres = 0;
  ts.iASq = 1;
  ts.iTemp = 2;
  ts.iOffset = 3;
  ts.axisEgas = 0;
  ts.axisPres = 1;
  ts.axisHint = 2;

  if (EOS_file_type.compare("binary") == 0) {
    ReadBinaryTable(EOS_fn, ptable_, ts);
  } else if (EOS_file_type.compare("hdf5") == 0) {
#ifndef HDF5OUTPUT
    std::stringstream msg;
    msg << "### FATAL ERROR in EquationOfState::PrepEOS" << std::endl
        << "HDF5 EOS table specified, but HDF5 flag is not enabled."  << std::endl;
    throw std::runtime_error(msg.str().c_str());
#endif
    HDF5TableLoader(EOS_fn.c_str(), ptable_, 9, var_names, dens_lim_field.c_str(), espec_lim_field.c_str());
    ptable_->GetSize(ts.nVar, ts.nEgas, ts.nRho);
    ptable_->GetX2lim(ts.logEgasMin, ts.logEgasMax);
    ptable_->GetX1lim(ts.logRhoMin, ts.logRhoMax);
    if (read_ratios) {
      AthenaArray<Real> ratio_array;
      ratio_array.NewAthenaArray(3);
      std::string ratio_field = pin->GetOrAddString("hydro", "EOS_ratio_field", "ratios");
      int zero[] = {0};
      const int one = 1;
      int three[] = {3};
      HDF5ReadRealArray(EOS_fn.c_str(), ratio_field.c_str(), one, zero, three, one, zero, three, ratio_array);
      ts.EosRatios[0] = ratio_array(0);
      ts.EosRatios[1] = ratio_array(1);
      ts.EosRatios[2] = ratio_array(2);
      ts.egasOverPres = ts.EosRatios[1];
      //std::cout << "Ratios: " << ratio_array(0) << ", " << ratio_array(1) << ", " << ratio_array(2) << '\n';
      ratio_array.DeleteAthenaArray();
    }
  } else if (EOS_file_type.compare("ascii") == 0) {
    AthenaArray<Real> *pratios = NULL;
    if (read_ratios) pratios = new AthenaArray<Real>;
    ASCIITableLoader(EOS_fn.c_str(), ptable_, pratios);
    if (read_ratios) {
      ts.EosRatios[0] = (*pratios)(0);
      ts.EosRatios[1] = (*pratios)(1);
      ts.EosRatios[2] = (*pratios)(2);
      ts.egasOverPres = ts.EosRatios[1];
      pratios->DeleteAthenaArray();
      delete pratios;
    }
    ptable_->GetSize(ts.nVar, ts.nEgas, ts.nRho);
    ptable_->GetX2lim(ts.logEgasMin, ts.logEgasMax);
    ptable_->GetX1lim(ts.logRhoMin, ts.logRhoMax);
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in EquationOfState::PrepEOS" << std::endl
        << "EOS table of type '" << EOS_file_type << "' not recognized."  << std::endl
        << "Options are 'ascii', 'binary', and 'hdf5'." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  ts.rhoNorm = (ts.nRho - 1.) / (ts.logRhoMax - ts.logRhoMin);
  ts.eNorm = (ts.nEgas - 1.) / (ts.logEgasMax - ts.logEgasMin);
  ts.rhoUnit = pin->GetOrAddReal("hydro", "EosRhoUnit", 1.);
  ts.eUnit = pin->GetOrAddReal("hydro", "EosEgasUnit", 1.);

#ifdef EOSDEBUG1
  std::cout << "Eos table: " << ts.nVar << ", " << ts.nRho << ", " << ts.nEgas << "\n";
  std::cout << "logRhoMin, logRhoMax: " << ts.logRhoMin << ", " << ts.logRhoMax << "\n";
  std::cout << "logEgasMin, logEgasMax: " << ts.logEgasMin << ", " << ts.logEgasMax << "\n";
  for (int i=0;i<ts.nVar;i++){
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
  std::cout << "prepEOS: " << EOS_fn << ", " << ts.nVar << ", " << ts.nRho << ", " << ts.nEgas << "\n";
  EosTestLoop();
#endif
}

void EquationOfState::CleanEOS()
{
  ptable_->~InterpTable2D();
}

Real EquationOfState::GetEosData(Real rho, Real var, int axis, int kOut)
{
  Real rhoIndex;
  Real varIndex;
  rhoIndex = log10(rho * ts.rhoUnit);
  varIndex = log10(var * ts.EosRatios[axis] * ts.eUnit) - rhoIndex;
  rhoIndex = (rhoIndex - ts.logRhoMin) * ts.rhoNorm;
  return std::pow((Real)10, ptable_->interpolate(kOut + ts.iOffset * axis, rhoIndex, varIndex));
}

void EquationOfState::EosTestRhoEgas(Real rho, Real egas, AthenaArray<Real> &data)
{
  Real idn = 1./rho;
  Real ien = 1./egas;
  data(0) = GetEosData(rho, egas, ts.axisEgas, ts.iPres) * egas;
  data(1) = data(0) + egas;
  data(2) = GetEosData(rho, egas, ts.axisEgas, ts.iASq) * egas * idn;
  data(3) = GetEosData(rho, egas, ts.axisEgas, ts.iTemp) * egas * idn;
  data(4) = (egas -  GetEosData(rho, data(0), ts.axisPres, ts.iPres) * data(0)) * ien;
  data(5) = (data(2) - GetEosData(rho, data(0), ts.axisPres, ts.iASq) * data(0) * idn) / data(2);
  data(6) = (egas - GetEosData(rho, data(1), ts.axisHint, ts.iPres) * data(1)) * ien;
  data(7) = (data(2) - GetEosData(rho, data(1), ts.axisHint, ts.iASq) * data(1) * idn) / data(2);
}

void EquationOfState::EosTestLoop()
{
  Real rho = 0;
  Real egas;
  AthenaArray<Real> data;
  data.NewAthenaArray(8);
  std::cout << "logRhoMin, logRhoMax: " << ts.logRhoMin << ", " << ts.logRhoMax << "\n";
  std::cout << "logEgasMin, logEgasMax: " << ts.logEgasMin << ", " << ts.logEgasMax << "\n";
  std::cout << "egasOverPres = " << ts.egasOverPres << "\n";
  std::cout << "Input rho (g/cc):";
  std::cin >> rho;
  std::cout << "Input egas (erg/cc):";
  std::cin >> egas;
  while (rho >= 0){
    std::cout << "log10(e/p) = " << GetEosData(rho, egas, ts.axisEgas, ts.iPres) << "\n";
    EosTestRhoEgas(rho, egas, data);
    std::cout << rho << ", " << egas << "\n";
    std::cout << "P(d, e)    , h(d, e)    , ASq(d, e)  , T(d, e) , PErr    , ASqErr  , hErr    , AsqErr\n";
    for (int i=0;i<8;i++){
      std::cout << data(i) << ", ";
    }
    std::cout << "\n";
    std::cout << "P, e, Asq, Asq: " << SimplePres(rho, egas) << ", "
                                    << SimpleEgas(rho, data(0)) << ", "
                                    << SimpleAsq(rho, egas) << ", "
                                    << RiemannAsq(rho, data(1)) << "\n";
    std::cout << "Input rho (g/cc):";
    std::cin >> rho;
    std::cout << "Input egas (erg/cc):";
    std::cin >> egas;
  }
  data.DeleteAthenaArray();
}

#endif
