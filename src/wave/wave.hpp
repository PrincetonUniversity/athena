#ifndef WAVE_HPP
#define WAVE_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file wave.hpp
//  \brief definitions for the Wave class

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../finite_differencing.hpp"
#include "../bvals/cc/bvals_cc.hpp"

class MeshBlock;
class ParameterInput;

//! \class Wave
//  \brief Wave data and functions

class Wave {
public:
  Wave(MeshBlock *pmb, ParameterInput *pin);
  ~Wave();

  // data
  MeshBlock *pmy_block;         // ptr to MeshBlock containing this Wave
  AthenaArray<Real> u;          // solution of the wave equation
  AthenaArray<Real> u1, u2;     // auxiliary arrays at intermediate steps
  AthenaArray<Real> rhs;        // wave equation rhs

  AthenaArray<Real> exact;      // exact solution of of the wave equation
  AthenaArray<Real> error;      // error with respect to the exact solution

  Real c;                       // light speed

  // control whether radiative condition is applied for outflow  or
  // extrapolate_outflow BC;
  // 0: not applied
  // 1,2,3: applied in respective dimensions
  int use_Sommerfeld = 0;

  // boundary and grid data
  CellCenteredBoundaryVariable ubvar;
  AthenaArray<Real> empty_flux[3];

  // storage for SMR/AMR
  // TODO(KGF): remove trailing underscore or revert to private:
  AthenaArray<Real> coarse_u_;
  int refinement_idx{-1};

  // functions
  Real NewBlockTimeStep(void);  // compute new timestep on a MeshBlock
  void WaveRHS(AthenaArray<Real> &u);
  void WaveBoundaryRHS(AthenaArray<Real> &u);
  void AddWaveRHS(const Real wght, AthenaArray<Real> &u_out);

  static const int NWAVE_CPT = 2;      // num. of wave equation field components

private:
  AthenaArray<Real> dt1_,dt2_,dt3_;    // scratch arrays used in NewTimeStep

private:
  void WaveSommerfeld_1d_L_(AthenaArray<Real> & u,
                            int const is, int const ie,
                            int const js, int const je,
                            int const ks, int const ke);
  void WaveSommerfeld_1d_R_(AthenaArray<Real> & u,
                            int const is, int const ie,
                            int const js, int const je,
                            int const ks, int const ke);
  void WaveSommerfeld_2d_(AthenaArray<Real> & u,
                          int const is, int const ie,
                          int const js, int const je,
                          int const ks, int const ke);
  void WaveSommerfeld_3d_(AthenaArray<Real> & u,
                          int const is, int const ie,
                          int const js, int const je,
                          int const ks, int const ke);
private:
  struct {
    typedef FDCenteredStencil<2, NGHOST> stencil;

    int stride[3];
    Real idx[3];

    inline Real Ds(int dir, Real & u) {
      Real * pu = &u;
      return 0.5 * idx[dir] * (pu[stride[dir]] - pu[-stride[dir]]);
    }
    inline Real Dxx(int dir, Real & u) {
      Real * pu = &u - stencil::offset*stride[dir];

      Real out(0.);
      for(int n1 = 0; n1 < stencil::nghost; ++n1) {
        int const n2  = stencil::width - n1 - 1;
        Real const c1 = stencil::coeff[n1] * pu[n1*stride[dir]];
        Real const c2 = stencil::coeff[n2] * pu[n2*stride[dir]];
        out += (c1 + c2);
      }
      out += stencil::coeff[stencil::nghost] * pu[stencil::nghost*stride[dir]];

      return out*SQR(idx[dir]);
    }
  } FD;

};
#endif // WAVE_HPP
