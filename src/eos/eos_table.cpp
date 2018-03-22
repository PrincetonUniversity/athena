//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
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
//#include "../hydro.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
//#include "../mesh.hpp"
#include "../parameter_input.hpp"
#include "../field/field.hpp"
#include "../coordinates/coordinates.hpp"

// this class header
#include "eos.hpp"

#if EOS_TABLE_ENABLED
//#define EOSDEBUG0
//#define EOSDEBUG1

std::string DefaultEOS()
{
  std::string fn("eos_table.data");
  return fn;
}

Real SimplePres(Real rho, Real egas, EquationOfState* peos)
{
  return peos->GetEosData(rho, egas, peos->axisEgas, peos->iPresEOS) * egas;
}

Real SimpleEgas(Real rho, Real pres, EquationOfState* peos)
{
  return peos->GetEosData(rho, pres, peos->axisPres, peos->iPresEOS) * pres;
}

Real SimpleAsq(Real rho, Real pres, EquationOfState* peos)
{
  return peos->GetEosData(rho, pres, peos->axisPres, peos->iASqEOS) * pres / rho;
}

Real RoeAsq(Real rho, Real hint, EquationOfState* peos)
{
  return peos->GetEosData(rho, hint, peos->axisHint, peos->iASqEOS) * hint / rho;
}

void EquationOfState::EnrollEosTable(EosFn_t GenEosTableFilename)
{
  GetEosFn = GenEosTableFilename;
}

void EquationOfState::PrepEOS(ParameterInput *pin)
{
  std::string EosFn;
  GetEosFn = ( GetEosFn == NULL ) ? DefaultEOS : GetEosFn;
  EosFn = pin->GetOrAddString("hydro", "EosFn", (*GetEosFn)());
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

  EnrollSimplePres(&SimplePres);
  EnrollSimpleEgas(&SimpleEgas);
  EnrollSimpleAsq(&SimpleAsq);
  EnrollAsqFromHint(&RoeAsq);

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
    //printf("WARNING!\nY value off table, extrapolating.\n");
    std::cout << "WARNING!\nY value off table, extrapolating. " << y <<"\n";
    if ( y != y ) throw std::runtime_error("Meh");
#endif
  }
  else if (yil >= ny - 1)
  {
    yil = ny - 2;
#ifdef EOSDEBUG0
    //printf("WARNING!\nY value off table, extrapolating.\n");
    std::cout << "WARNING!\nY value off table, extrapolating. " << y <<"\n";
    if ( y != y ) throw std::runtime_error("Meh");
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
    std::cout << "P, e, Asq, Asq: " << (*SimplePres_)(rho, egas, this) << ", "
                                    << (*SimpleEgas_)(rho, data(0), this) << ", "
                                    << (*AsqFromPres_)(rho, egas, this) << ", "
                                    << (*AsqFromHint_)(rho, data(1), this) << "\n";
    std::cout << "Input rho (g/cc):";
    std::cin >> rho;
    std::cout << "Input egas (erg/cc):";
    std::cin >> egas;
  }
  data.DeleteAthenaArray();
  //throw std::runtime_error("Meh");
}

#endif
