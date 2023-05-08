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
#include <random>     // mt19937, normal_distribution, uniform_real_distribution

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
  void PowerSpectrum(std::complex<Real> *amp);
  void Perturb(Real dt);
  void OUProcess(Real dt);
  void Project(std::complex<Real> **fv, Real f_shear);
  void Project(std::complex<Real> **fv, std::complex<Real> **fv_sh,
               std::complex<Real> **fv_co);
  std::int64_t GetKcomp(int idx, int disp, int Nx);
 private:
  std::int64_t rseed;
  int nlow, nhigh;
  Real tdrive, dtdrive, tcorr, f_shear;
  Real expo, dedt, dvol;
  AthenaArray<Real> vel[3];
  std::complex<Real> **fv_, **fv_new_;
  std::complex<Real> **fv_sh_, **fv_co_;
  bool initialized_ = false;
  bool global_ps_ = false;
  std::mt19937_64 rng_generator;
};

#endif // FFT_TURBULENCE_HPP_
