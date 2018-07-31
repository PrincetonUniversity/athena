#ifndef FFT_TURBULENCE_HPP_
#define FFT_TURBULENCE_HPP_

//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file turbulence.hpp
//  \brief defines Turbulence class

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "athena_fft.hpp"

class Mesh;
class MeshBlock;
class ParameterInput;
class Coordinates;
class FFTBlock;
class FFTDriver;

//! \class TurbulenceDriver
//  \brief Turbulence Driver

class TurbulenceDriver : public FFTDriver{
public:
  TurbulenceDriver(Mesh *pm, ParameterInput *pin);
  ~TurbulenceDriver();
  void Driving(void);
  void Generate(void);
  void PowerSpectrum(AthenaFFTComplex *amp);
  void Perturb(Real dt);
  int64_t GetKcomp(int idx, int disp, int Nx);
private:
  int64_t rseed;
  int nlow,nhigh;
  Real dtdrive,tdrive;
  Real expo,dedt,dvol;
  AthenaArray<Real> *vel;
};

#endif // FFT_TURBULENCE_HPP_
