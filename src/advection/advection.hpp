#ifndef ADVECTION_HPP
#define ADVECTION_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file advection.hpp
//  \brief definitions for the Advection class

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../finite_differencing.hpp"
#include "../bvals/cc/bvals_cc.hpp"

class MeshBlock;
class ParameterInput;

//! \class Advection
//  \brief Advection data and functions

class Advection {
public:
  Advection(MeshBlock *pmb, ParameterInput *pin);
  ~Advection();

  // data
  MeshBlock *pmy_block;         // ptr to MeshBlock containing this Advection
  AthenaArray<Real> u;          // solution of the advection equation
  AthenaArray<Real> u1, u2;     // auxiliary arrays at intermediate steps
  AthenaArray<Real> rhs;        // advection equation rhs

  AthenaArray<Real> exact;      // exact solution of of the advection equation
  AthenaArray<Real> error;      // error with respect to the exact solution

  // propagation velocity components
  Real cx1;
  Real cx2;
  Real cx3;

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
  void AdvectionRHS(AthenaArray<Real> &u);
  void AdvectionBoundaryRHS(AthenaArray<Real> &u);
  void AddAdvectionRHS(const Real wght, AthenaArray<Real> &u_out);

private:
  AthenaArray<Real> dt1_,dt2_,dt3_;    // scratch arrays used in NewTimeStep

private:
  struct {
    // 1st derivative stecil
    typedef FDCenteredStencil<1, NGHOST> s1;

    typedef FDCenteredStencil<1, NGHOST> stencil;

    // left-biased derivative
    typedef FDLeftBiasedStencil<
      FDBiasedChoice<1, NGHOST-1>::degree,
      FDBiasedChoice<1, NGHOST-1>::nghost,
      FDBiasedChoice<1, NGHOST-1>::lopsize
      > sl;
    // right-biased derivative
    typedef FDRightBiasedStencil<
      FDBiasedChoice<1, NGHOST-1>::degree,
      FDBiasedChoice<1, NGHOST-1>::nghost,
      FDBiasedChoice<1, NGHOST-1>::lopsize
      > sr;

    int stride[3];
    Real idx[3];

    // 1st derivative (high order centered)
    inline Real Dx(int dir, Real & u) {
      Real * pu = &u - s1::offset*stride[dir];

      Real out(0.);
      for(int n1 = 0; n1 < s1::nghost; ++n1) {
        int const n2  = s1::width - n1 - 1;
        Real const c1 = s1::coeff[n1] * pu[n1*stride[dir]];
        Real const c2 = s1::coeff[n2] * pu[n2*stride[dir]];
        out += (c1 + c2);
      }
      out += s1::coeff[s1::nghost] * pu[s1::nghost*stride[dir]];

      return out * idx[dir];
    }
    // 1st derivative 2nd order centered
    inline Real Ds(int dir, Real & u) {
      Real * pu = &u;
      return 0.5 * idx[dir] * (pu[stride[dir]] - pu[-stride[dir]]);
    }

    // Advective derivative
    // The advective derivative is for an equation in the form
    //    d_t u = -vx d_x u
    // vx>0 => solution propagates rightward
    //
    // vx>0 => we use left-sided stencil
    // vx<0 => we use right-sided stencil
    //
    // internally both are computed
    inline Real Da_x(int dir, Real & vx, Real & u) {
      Real * pu = &u;

      Real dl(0.);
      for(int n = 0; n < sl::width; ++n) {
        dl += sl::coeff[n] * pu[(n - sl::offset)*stride[dir]];
      }

      Real dr(0.);
      for(int n = sr::width-1; n >= 0; --n) {
        dr += sr::coeff[n] * pu[(n - sr::offset)*stride[dir]];
      }

      return ((vx < 0) ? (vx * dr) : (vx * dl)) * idx[dir];
    }


  } FD;

};
#endif // ADVECTION_HPP
