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
#include <stdexcept> // std::invalid_argument

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../field/field.hpp"
#include "../coordinates/coordinates.hpp"

// this class header
#include "eos.hpp"

#if EOS_TABLE_ENABLED
//#define EOSDEBUG0
//#define EOSDEBUG1

Real EquationOfState::SimplePres(Real rho, Real egas) {
  return GetEosData(rho, egas, axisEgas, iPresEOS) * egas;
}

Real EquationOfState::SimpleEgas(Real rho, Real pres) {
  return GetEosData(rho, pres, axisPres, iPresEOS) * pres;
}

Real EquationOfState::SimpleAsq(Real rho, Real pres) {
  return GetEosData(rho, pres, axisPres, iASqEOS) * pres / rho;
}

Real EquationOfState::RiemannAsq(Real rho, Real hint) {
  return GetEosData(rho, hint, axisHint, iASqEOS) * hint;
}

void EquationOfState::PrepEOS(ParameterInput *pin) {
  std::string EosFn;
  EosFn = pin->GetString("hydro", "EosFn");
  std::ifstream eos_table(EosFn.c_str(), std::ios::binary);
  if (eos_table.is_open())
  {
    eos_table.seekg(0, std::ios::beg);
    eos_table.read((char*)&nRho_, sizeof(nRho_));
    eos_table.read((char*)&logRhoMin_, sizeof(logRhoMin_));
    eos_table.read((char*)&logRhoMax_, sizeof(logRhoMax_));
    eos_table.read((char*)&nEgas_, sizeof(nEgas_));
    eos_table.read((char*)&logEgasMin_, sizeof(logEgasMin_));
    eos_table.read((char*)&logEgasMax_, sizeof(logEgasMax_));
    eos_table.read((char*)&egasOverPres_, sizeof(egasOverPres_));
    eos_table.read((char*)&nVar_, sizeof(nVar_));
    std::cout << "prepEOS: " << EosFn << ", " << nVar_ << ", " << nRho_ << ", " << nEgas_ << "\n";
    eos_data_.NewAthenaArray(nVar_, nRho_, nEgas_);
    eos_table.read((char*)eos_data_.data(), nVar_ * nRho_ * nEgas_ * sizeof(logRhoMin_));
    eos_table.close();
  }
  else throw std::invalid_argument("Unable to open eos table: " + EosFn);

  rhoNorm_ = (nRho_ - 1.) / (logRhoMax_ - logRhoMin_);
  eNorm_ = (nEgas_ - 1.) / (logEgasMax_ - logEgasMin_);

  iPresEOS = 0;
  iASqEOS = 1;
  iTempEOS = 2;
  iOffsetEOS = 3;
  axisEgas = 0;
  axisPres = 1;
  axisHint = 2;
  EosRatios_[0] = 1;
  EosRatios_[1] = egasOverPres_;
  EosRatios_[2] = egasOverPres_ / (1. + egasOverPres_);

  rhoUnit_ = pin->GetOrAddReal("hydro", "EosRhoUnit", 1.);
  eUnit_ = pin->GetOrAddReal("hydro", "EosEgasUnit", 1.);

#ifdef EOSDEBUG1
  std::cout << "Eos table:\n";
  for (int i=0;i<nVar_;i++){
    std::cout << "var = " << i << "\n";
    for (int j=0;j<nRho_;j++){
      for (int k=0;k<nEgas_;k++){
        std::cout << std::pow((Real) 10, eos_data_(i,j,k)) << " ";
      }
      std::cout << "\n";
    }
    std::cout << "\n";
  }
#endif

#ifdef EOSDEBUG0
  EosTestLoop();
#endif
}

void EquationOfState::CleanEOS()
{
  eos_data_.DeleteAthenaArray();
}

/* bilinear interpolation */
/* 0 <= x <= nx-1, 0 <= y <= ny-1 */
Real EquationOfState::BilinearInterp(Real x, Real y, int var)
{
  Real xrl, yrl, out;
  int xil = (int) x; /* index lower */
  int yil = (int) y;
  int nx = nRho_;
  int ny = nEgas_;
  if (xil < 0)
  {
    xil = 0;
#ifdef EOSDEBUG0
    printf("WARNING!\nX value off table, extrapolating.\n");
#endif
  }
  else if (xil >= nx - 1)
  {
    xil = nx - 2;
#ifdef EOSDEBUG0
    printf("WARNING!\nX value off table, extrapolating.\n");
#endif
  }
  xrl = 1 + xil - x;  /* residual lower */

  if (yil < 0)
  {
    yil = 0;
#ifdef EOSDEBUG0
    printf("WARNING!\nY value off table, extrapolating.\n");
#endif
  }
  else if (yil >= ny - 1)
  {
    yil = ny - 2;
#ifdef EOSDEBUG0
    printf("WARNING!\nY value off table, extrapolating.\n");
#endif
  }
  yrl = 1 + yil - y;  /* residual lower */

#ifdef EOSDEBUG1
  std::cout << "xil, yil = " << xil <<  ", " << yil << "\n";
  std::cout << "xrl, yrl = " << xrl <<  ", " << yrl << "\n";
  std::cout << "Interp: var = " << var << "\n";
#endif
  out =   xrl  *  yrl  *eos_data_(var, xil , yil )
      +   xrl  *(1-yrl)*eos_data_(var, xil ,yil+1)
      + (1-xrl)*  yrl  *eos_data_(var,xil+1, yil )
      + (1-xrl)*(1-yrl)*eos_data_(var,xil+1,yil+1);
#ifdef EOSDEBUG1
  std::cout << eos_data_(var, xil ,yil) << ", " << eos_data_(var, xil ,yil+1) << ", "
            << eos_data_(var,xil+1,yil) << ", " << eos_data_(var,xil+1,yil+1) << "\n";
#endif
  //std::cout << "\n";
  return out;
}

void EquationOfState::GetEosIndices(Real rho, Real var, int axis, Real &rhoIndex, Real &varIndex)
{
  rhoIndex = log10(rho * rhoUnit_);
  varIndex = (log10(var * EosRatios_[axis] * eUnit_) - rhoIndex - logEgasMin_) * eNorm_;
  rhoIndex = (rhoIndex - logRhoMin_) * rhoNorm_;
#ifdef EOSDEBUG1
  std::cout << rho << ", " << var * EosRatios_[axis] << ", (" << rhoIndex << ", " << varIndex << ")" << ", Axis: " << axis << "\n";
#endif
}

Real EquationOfState::GetEosData(Real rho, Real var, int axis, int kOut)
{
  Real rhoIndex;
  Real varIndex;
  GetEosIndices(rho, var, axis, rhoIndex, varIndex);
  return std::pow((Real)10, BilinearInterp(rhoIndex, varIndex, kOut + iOffsetEOS * axis));
}

void EquationOfState::EosTestRhoEgas(Real rho, Real egas, AthenaArray<Real> &data)
{
  Real idn = 1./rho;
  Real ien = 1./egas;
  data(0) = GetEosData(rho, egas, axisEgas, iPresEOS) * egas;
  data(1) = data(0) + egas;
  data(2) = GetEosData(rho, egas, axisEgas, iASqEOS) * egas * idn;
  data(3) = GetEosData(rho, egas, axisEgas, iTempEOS) * egas * idn;
  data(4) = (egas -  GetEosData(rho, data(0), axisPres, iPresEOS) * data(0)) * ien;
  data(5) = (data(2) - GetEosData(rho, data(0), axisPres, iASqEOS) * data(0) * idn) / data(2);
  data(6) = (egas - GetEosData(rho, data(1), axisHint, iPresEOS) * data(1)) * ien;
  data(7) = (data(2) - GetEosData(rho, data(1), axisHint, iASqEOS) * data(1) * idn) / data(2);
}

void EquationOfState::EosTestLoop()
{
  Real rho = 0;
  Real egas, iRho, iEgas;
  AthenaArray<Real> data;
  data.NewAthenaArray(8);
  std::cout << "logRhoMin_, logRhoMax_: " << logRhoMin_ << ", " << logRhoMax_ << "\n";
  std::cout << "logEgasMin_, logEgasMax_: " << logEgasMin_ << ", " << logEgasMax_ << "\n";
  std::cout << "egasOverPres_ = " << egasOverPres_ << "\n";
  std::cout << "Input rho (g/cc):";
  std::cin >> rho;
  std::cout << "Input egas (erg/cc):";
  std::cin >> egas;
  while (rho >= 0){
    GetEosIndices(rho, egas, axisEgas, iRho, iEgas);
    std::cout << "log10(e/p) = " << GetEosData(rho, egas, axisEgas, iPresEOS) << ", " << std::pow((Real) 10, eos_data_(0, int(iRho + .5), int(iEgas + .5))) << "\n";
    EosTestRhoEgas(rho, egas, data);
    std::cout << rho << ", " << egas << ", (" << iRho << ", " << iEgas << "): \n";
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
