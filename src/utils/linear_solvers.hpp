//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file linear_solvers.hpp
//  \brief Functionality for solution of linear systems

#ifndef LINEAR_SOLVERS_HPP
#define LINEAR_SOLVERS_HPP

// C headers

// C++ headers
#include<vector>     // std::vector
#include<algorithm>  // std::reverse

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

class Linear_banded_mask {
  public:
    Linear_banded_mask(
      // size along each dimension
      const std::vector<unsigned int> & dims_N,
      // coefficient arrays
      const std::vector<std::vector<Real>> & coeffL,
      const std::vector<Real> & coeffC,
      const std::vector<std::vector<Real>> & coeffR,
      // Stencil base-point pos relative to diag. ix
      const std::vector<std::vector<unsigned int>> & offset,
      // Range of idxs to apply over on target vec.
      // [relative to start / end of vec.]
      const std::vector<std::vector<int>> & rix,
      // {Left, center, right} stencil widths
      const std::vector<std::vector<unsigned int>> & width
    );
    virtual ~Linear_banded_mask();

    void Mask(
      const unsigned int axis,
      const AthenaArray<Real> & y,
      AthenaArray<Real> & x
    );

  private:
    unsigned int * N;
    unsigned int dim;
    unsigned int * strides;

    std::vector<std::vector<Real>> coeffL, coeffR;
    std::vector<Real> coeffC;

    std::vector<std::vector<unsigned int>> offset;
    std::vector<std::vector<int>> rix;
    std::vector<std::vector<unsigned int>> width;

    void imp_Mask(
      const unsigned int axis,
      const unsigned int idx_offset,
      const AthenaArray<Real> & y,
      AthenaArray<Real> & x
    );

};

class Linear_banded_direct {
  // Direct, non-destructive implementations of banded solvers.
  // Scratch memory is automatically assigned
  public:
    static void Tridiagonal(
      const AthenaArray<Real> &T_a,
      const AthenaArray<Real> &T_b,
      const AthenaArray<Real> &T_c,
      const AthenaArray<Real> &y,
      AthenaArray<Real> &x,
      const unsigned int sz
    );

    static void Pentadiagonal(
      const AthenaArray<Real> &P_a,
      const AthenaArray<Real> &P_b,
      const AthenaArray<Real> &P_c,
      const AthenaArray<Real> &P_d,
      const AthenaArray<Real> &P_e,
      const AthenaArray<Real> &y,
      AthenaArray<Real> &x,
      const unsigned int sz
    );

};

class solver_Tridiagonal {
  public:
    solver_Tridiagonal(const std::vector<unsigned int> & dims_N);
    virtual ~solver_Tridiagonal();

    void Prepare(
      const unsigned int axis,
      const std::vector<AthenaArray<Real> *> & bands
    );
    void Solve(
      const unsigned int axis,
      const AthenaArray<Real> & y,
      AthenaArray<Real> & x
    );

  private:
    // for precomputed factors
    AthenaArray<Real> * _c_star, * _m, * _a_m, * _z;

    void factor_computation(
      const unsigned int dim,
      const AthenaArray<Real> &T_a,
      const AthenaArray<Real> &T_b,
      const AthenaArray<Real> &T_c
    );

    void imp_Solve(
      const unsigned int axis,
      const unsigned int idx_offset,
      const AthenaArray<Real> & y,
      AthenaArray<Real> & x
    );

    unsigned int * N;
    unsigned int dim;
    unsigned int * strides;

};

class solver_Pentadiagonal {
  public:
    solver_Pentadiagonal(const std::vector<unsigned int> & dims_N);
    virtual ~solver_Pentadiagonal();

    void Prepare(
      const unsigned int axis,
      const std::vector<AthenaArray<Real> *> & bands
    );
    void Solve(
      const unsigned int axis,
      const AthenaArray<Real> & y,
      AthenaArray<Real> & x
    );

  private:
    // for precomputed factors
    AthenaArray<Real> * _al, *_be, *_ga_mu, *_e_mu, *_i_mu, *_z;

    void factor_computation(
      const unsigned int dim,
      const AthenaArray<Real> &P_a,
      const AthenaArray<Real> &P_b,
      const AthenaArray<Real> &P_c,
      const AthenaArray<Real> &P_d,
      const AthenaArray<Real> &P_e
    );

    void imp_Solve(
      const unsigned int axis,
      const unsigned int idx_offset,
      const AthenaArray<Real> & y,
      AthenaArray<Real> & x
    );

    unsigned int * N;
    unsigned int dim;
    unsigned int * strides;
};

class Linear_solver {
  // General interface for solvers that precompute factors
  public:
    Linear_solver(
      unsigned int bandwidth,
      const std::vector<unsigned int> & dims_N
    );

    void Prepare(
      const unsigned int axis,
      const std::vector<AthenaArray<Real> *> & bands
    );

    void Solve(
      const unsigned int axis,
      const AthenaArray<Real> & y,
      AthenaArray<Real> & x
    );

    virtual ~Linear_solver();

  private:
    unsigned int dim;
    unsigned int bandwidth;
    solver_Tridiagonal * solver_Tri;
    solver_Pentadiagonal * solver_Pent;
};

class Linear_solver_utils {
  // some utilities for interacting with arrays
  public:
    Linear_solver_utils();
    virtual ~Linear_solver_utils();

    static void BandedAthenaArrayToVector(
      const unsigned int bandwidth,
      const AthenaArray<Real> & arr,
      std::vector<AthenaArray<Real> *> & bands
    );

    static void DeleteBandedVector(
      std::vector<AthenaArray<Real> *> & bands
    );

};

#endif // LINEAR_SOLVERS_HPP