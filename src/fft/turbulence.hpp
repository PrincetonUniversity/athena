#ifndef FFT_TURBULENCE_HPP_
#define FFT_TURBULENCE_HPP_

//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file turbulence.hpp
//  \brief defines Turbulence class

// C headers

// C++ headers

// Athena++ headers
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
  void Driving();
  void Generate();
  void PowerSpectrum(AthenaFFTComplex *amp);
  void Perturb(Real dt);
  void OUProcess(Real dt);
  void Project(AthenaFFTComplex **fv, Real f_shear);
  void Project(AthenaFFTComplex **fv, AthenaFFTComplex **fv_sh, AthenaFFTComplex **fv_co);
  std::int64_t GetKcomp(int idx, int disp, int Nx);
 private:
  std::int64_t rseed;
  int nlow,nhigh;
  Real tdrive,dtdrive,tcorr,f_shear;
  Real expo,dedt,dvol;
  AthenaArray<Real> *vel;
  AthenaFFTComplex **fv_, **fv_new_;
  AthenaFFTComplex **fv_sh_, **fv_co_;
  bool initialized_ = false;
};

#endif // FFT_TURBULENCE_HPP_
