//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file linear_solvers.hpp
//  \brief Functionality for solution of linear systems

#include "linear_solvers.hpp"

// ============================================================================
// Solver for x in Tx=y where T is tridiagonal with structure:
//
//   b_0     c_0      0      ...
//   a_0     b_1     c_1      0      ...
//    0      a_1     b_2     c_2      0      ...
//    .       0       .       .       .       .
//            .       0    a_{n-3} b_{n-2} c_{n-2}
//                    .       0    a_{n-2} b_{n-1}
// Therefore for P = [nxn]:
//   a = [n-1]
//   b = [n]
//   c = [n-1]
//
// This uses the Thomas algorithm.
void Linear_banded_direct::Tridiagonal(
  const AthenaArray<Real> &T_a,
  const AthenaArray<Real> &T_b,
  const AthenaArray<Real> &T_c,
  const AthenaArray<Real> &y,
  AthenaArray<Real> &x,
  const unsigned int sz)
{

  // Create temporary arrays
  AthenaArray<Real> c_star;
  AthenaArray<Real> d_star;
  c_star.NewAthenaArray(sz-1);
  d_star.NewAthenaArray(sz);

  // Forward sweep ------------------------------------------------------------
  c_star(0) = T_c(0) / T_b(0);
  d_star(0) = y(0) / T_b(0);

#pragma omp simd
  for(int ix=1; ix<sz-1; ++ix){
    const Real m{1. / (T_b(ix) - T_a(ix-1) * c_star(ix-1))};
    c_star(ix) = T_c(ix) * m;
    d_star(ix) = (y(ix) - T_a(ix-1) * d_star(ix-1)) * m;
  }

  const Real m{1. / (T_b(sz-1) - T_a(sz-2) * c_star(sz-2))};
  d_star(sz-1) = (y(sz-1) - T_a(sz-2) * d_star(sz-2)) * m;

  // backward substitution ----------------------------------------------------
  x(sz-1) = d_star(sz-1);
#pragma omp simd
  for (int ix=sz-2; ix >= 0; --ix) {
    x(ix) = d_star(ix) - c_star(ix) * x(ix+1);
  }

  // clean up temporary arrays
  c_star.DeleteAthenaArray();
  d_star.DeleteAthenaArray();

  return;
}
// ============================================================================

// ============================================================================
// Solver for x in Px=y where P is pentadiagonal with structure:
//
//   d_0     a_0     b_0      0      ...
//   c_0     d_1     a_1     b_1      0      ...
//   e_0     c_1     d_2     a_2     b_2      0      ...
//    0      e_1     c_2     d_3     a_3     b_3      0     ...
//    .       0       .       .       .       .       .      .
//            .       0    e_{n-5} c_{n-4} d_{n-3} a_{n-3} b_{n-3}
//                    .       0    e_{n-4} c_{n-3} d_{n-1} a_{n-2}
//                                    0    e_{n-3} c_{n-2} d_{n-1}
//
// Therefore for P = [nxn]:
//   a = [n-1]
//   b = [n-2]
//   c = [n-1]
//   d = [n]
//   e = [n-2]
//
// Particular algorithm is PTRANS-I of the reference below.
//
// Inputs are preserved.
//
// Note
// ----
// No checks are performed as to whether P is invertible.
//
// Reference
// ---------
//   S. S. Askar and A. A. Karawia
//   ‘On Solving Pentadiagonal Linear Systems via Transformations’,
//   Research Article, Mathematical Problems in Engineering
//   (Hindawi, 31 March 2015), https://doi.org/10.1155/2015/232456.
void Linear_banded_direct::Pentadiagonal(
  const AthenaArray<Real> &P_a,
  const AthenaArray<Real> &P_b,
  const AthenaArray<Real> &P_c,
  const AthenaArray<Real> &P_d,
  const AthenaArray<Real> &P_e,
  const AthenaArray<Real> &y,
  AthenaArray<Real> &x,
  const unsigned int sz)
{
  // Create temporary arrays
  AthenaArray<Real> al, be, ga, mu, z;
  al.NewAthenaArray(sz-1);
  be.NewAthenaArray(sz-2);
  ga.NewAthenaArray(sz-1);
  mu.NewAthenaArray(sz);
  z.NewAthenaArray(sz);

  // step 3
  mu(0) = P_d(0);
  al(0) = P_a(0) / mu(0);
  be(0) = P_b(0) / mu(0);
  z(0) = y(0) / mu(0);

  // step 4
  ga(0) = P_c(0);
  mu(1) = P_d(1) - al(0) * ga(0);
  al(1) = (P_a(1) - be(0) * ga(0)) / mu(1);
  be(1) = P_b(1) / mu(1);
  z(1) = (y(1) - z(0) * ga(0)) / mu(1);

  // step 5
  for(int ix=2; ix < sz - 2; ++ix){
    ga(ix-1) = P_c(ix-1) - al(ix-2) * P_e(ix-2);
    mu(ix) = P_d(ix) - be(ix-2) * P_e(ix-2) - al(ix-1) * ga(ix-1);
    al(ix) = (P_a(ix) - be(ix-1) * ga(ix-1)) / mu(ix);
    be(ix) = P_b(ix) / mu(ix);
    z(ix) = (y(ix) - z(ix-2) * P_e(ix-2) - z(ix-1) * ga(ix-1)) / mu(ix);
  }

  ga(sz-3) = P_c(sz-3) - al(sz-4) * P_e(sz-4);
  mu(sz-2) = P_d(sz-2) - be(sz-4) * P_e(sz-4) - al(sz-3) * ga(sz-3);
  al(sz-2) = (P_a(sz-2) - be(sz-3) * ga(sz-3)) / mu(sz-2);

  ga(sz-2) = P_c(sz-2) - al(sz-3) * P_e(sz-3);
  mu(sz-1) = P_d(sz-1) - be(sz-3) * P_e(sz-3) - al(sz-2) * ga(-1+sz-1);

  // z(sz-3) -> z(sz-4) in first appearance cf. algo
  z(sz-2) = (y(sz-2) - z(sz-4) * P_e(sz-4) - z(sz-3) * ga(sz-3)) / mu(sz-2);
  // z(sz-2) -> z(sz-3) in first appearance cf. algo
  z(sz-1) = (y(sz-1) - z(sz-3) * P_e(sz-3) - z(sz-2) * ga(sz-2)) / mu(sz-1);

  // step 6
  x(sz-1) = z(sz-1);
  x(sz-2) = z(sz-2) - al(sz-2) * x(sz-1);

  for(int ix=sz-3; ix >= 0; --ix) {
    x(ix) = z(ix) - al(ix) * x(ix+1) - be(ix) * x(ix+2);
  }
  // clean up temporary arrays
  al.DeleteAthenaArray();
  be.DeleteAthenaArray();
  ga.DeleteAthenaArray();
  mu.DeleteAthenaArray();
  z.DeleteAthenaArray();

  return;
}
// ============================================================================


// ============================================================================
// Tridiagonal solver; uses factor precomputation
solver_Tridiagonal::solver_Tridiagonal(
  const std::vector<unsigned int> & idims_N)
{

  dim = idims_N.size();

  // allocate the scratch arrays
  _c_star = new AthenaArray<Real> [dim];
  _m = new AthenaArray<Real> [dim];
  _a_m = new AthenaArray<Real> [dim];
  _z = new AthenaArray<Real> [dim];

  N = new unsigned int [dim];
  strides = new unsigned int[dim];

  for(int ix=0; ix<dim; ++ix) {
    N[ix] = idims_N[ix];

    _c_star[ix].NewAthenaArray(N[ix]-1);
    _m[ix].NewAthenaArray(N[ix]);
    _a_m[ix].NewAthenaArray(N[ix]-1);
    _z[ix].NewAthenaArray(N[ix]);
  }

  strides[0] = 1;
  for(int ix=1; ix<dim; ++ix)
    strides[ix] = N[ix-1] * strides[ix-1];
}

// precomputation of factors (based on banded input)
void solver_Tridiagonal::factor_computation(
  const unsigned int axis,
  const AthenaArray<Real> & T_a,
  const AthenaArray<Real> & T_b,
  const AthenaArray<Real> & T_c)
{
  const unsigned int sz = N[axis];

  AthenaArray<Real> & c_star {_c_star[axis]};
  AthenaArray<Real> & m {_m[axis]};
  AthenaArray<Real> & a_m {_a_m[axis]};

  c_star(0) = T_c(0) / T_b(0);
  m(0) = 1. / T_b(0);
  for(int ix=1; ix<sz-1; ++ix){
    m(ix) = 1. / (T_b(ix) - T_a(ix-1) * c_star(ix-1));
    c_star(ix) = T_c(ix) * m(ix);
    a_m(ix-1) = T_a(ix-1) * m(ix);
  }

  m(sz-1) = 1. / (T_b(sz-1) - T_a(sz-2) * c_star(sz-2));
  a_m(sz-2) = m(sz-1) * T_a(sz-2);
}

// prepare the solver
void solver_Tridiagonal::Prepare(
  const unsigned int dim,
  const std::vector<AthenaArray<Real> *> & bands)
{
  const AthenaArray<Real> * T_a {bands[0]};
  const AthenaArray<Real> * T_b {bands[1]};
  const AthenaArray<Real> * T_c {bands[2]};

  factor_computation(dim, *T_a, *T_b, *T_c);
  return;
}

// implementation details of the solver
void solver_Tridiagonal::imp_Solve(
  const unsigned int axis,
  const unsigned int idx_offset,
  const AthenaArray<Real> & y,
  AthenaArray<Real> & x)
{
  // Scratch and precomputed factors
  const AthenaArray<Real> & c_star {_c_star[axis]};
  const AthenaArray<Real> & m {_m[axis]};
  const AthenaArray<Real> & a_m {_a_m[axis]};
  AthenaArray<Real> & z {_z[axis]};

  // Deal with different strides along different dimensions
  Real * px = &(x(0));
  const Real * py = &(const_cast<AthenaArray<Real>&>(y)(0));

  const unsigned int sz = N[axis];

  // Forward sweep
  z(0) = py[0 + idx_offset] * m(0);

// #pragma omp simd
  for(int ix=1; ix<sz; ++ix){
    z(ix) = py[(ix) * strides[axis] + idx_offset] * m(ix) \
      - a_m(ix-1) * z(ix-1);
  }

  // Backward sweep
  px[(sz-1) * strides[axis] + idx_offset] = z(sz-1);
// #pragma omp simd
  for(int ix=sz-2; ix >= 0; --ix) {
    px[(ix) * strides[axis] + idx_offset] = z(ix) - c_star(ix) * px[(ix+1) \
      * strides[axis] + idx_offset];
  }
}

// wrapper to apply solver over various dimensions based on precomputed factors
void solver_Tridiagonal::Solve(
  const unsigned int axis,
  const AthenaArray<Real> & y,
  AthenaArray<Real> & x)
{
  if (dim==3) {
    const unsigned cJ[] = {2, 2, 1};
    const unsigned cI[] = {1, 0, 0};
    unsigned int J = cJ[axis];
    unsigned int I = cI[axis];
    for(int jx=0; jx<N[cJ[axis]]; ++jx) {
      for(int ix=0; ix<N[cI[axis]]; ++ix) {
        imp_Solve(axis, ix * strides[cI[axis]] + jx * strides[cJ[axis]], y, x);
      }
    }
  } else if (dim==2) {
    const unsigned int cI[] = {1, 0}; // 1-axis

    for(int ix=0; ix<N[cI[axis]]; ++ix) {
      imp_Solve(axis, ix * strides[cI[axis]], y, x);
    }
  } else {
    imp_Solve(axis, 0, y, x);
  }
}

solver_Tridiagonal::~solver_Tridiagonal()
{
  for(int ix=0; ix<dim; ++ix) {
    _c_star[ix].DeleteAthenaArray();
    _m[ix].DeleteAthenaArray();
    _a_m[ix].DeleteAthenaArray();
    _z[ix].DeleteAthenaArray();
  }

  delete[] _c_star;
  delete[] _m;
  delete[] _a_m;
  delete[] _z;

  delete[] N;
  delete[] strides;
}
// ============================================================================

// ============================================================================
// Pentadiagonal solver; uses factor precomputation
solver_Pentadiagonal::solver_Pentadiagonal(
  const std::vector<unsigned int> & idims_N)
{

  dim = idims_N.size();

  // allocate the scratch arrays
  _al = new AthenaArray<Real> [dim];
  _be = new AthenaArray<Real> [dim];
  _ga_mu = new AthenaArray<Real> [dim];
  _e_mu = new AthenaArray<Real> [dim];
  _i_mu = new AthenaArray<Real> [dim];
  _z = new AthenaArray<Real> [dim];

  N = new unsigned int [dim];
  strides = new unsigned int[dim];

  for(int ix=0; ix<dim; ++ix) {
    N[ix] = idims_N[ix];

    _al[ix].NewAthenaArray(N[ix]-1);
    _be[ix].NewAthenaArray(N[ix]-2);
    _ga_mu[ix].NewAthenaArray(N[ix]-1);
    _e_mu[ix].NewAthenaArray(N[ix]-2);
    _i_mu[ix].NewAthenaArray(N[ix]);
    _z[ix].NewAthenaArray(N[ix]);
  }

  strides[0] = 1;
  for(int ix=1; ix<dim; ++ix)
    strides[ix] = N[ix-1] * strides[ix-1];
}

// precomputation of factors (based on banded input)
void solver_Pentadiagonal::factor_computation(
  const unsigned int axis,
  const AthenaArray<Real> & P_a,
  const AthenaArray<Real> & P_b,
  const AthenaArray<Real> & P_c,
  const AthenaArray<Real> & P_d,
  const AthenaArray<Real> & P_e)
{
  const unsigned int sz = N[axis];

  // Create temporary arrays
  AthenaArray<Real> ga, mu, z;
  ga.NewAthenaArray(sz-1);
  mu.NewAthenaArray(sz);

  AthenaArray<Real> & al {_al[axis]};
  AthenaArray<Real> & be {_be[axis]};
  AthenaArray<Real> & ga_mu {_ga_mu[axis]};
  AthenaArray<Real> & e_mu {_e_mu[axis]};
  AthenaArray<Real> & i_mu {_i_mu[axis]};

  // step 3
  mu(0) = P_d(0);
  al(0) = P_a(0) / mu(0);
  be(0) = P_b(0) / mu(0);

  // step 4
  ga(0) = P_c(0);
  mu(1) = P_d(1) - al(0) * ga(0);
  al(1) = (P_a(1) - be(0) * ga(0)) / mu(1);
  be(1) = P_b(1) / mu(1);

  // step 5
  for(int ix=2; ix < sz - 2; ++ix){
    ga(ix-1) = P_c(ix-1) - al(ix-2) * P_e(ix-2);
    mu(ix) = P_d(ix) - be(ix-2) * P_e(ix-2) - al(ix-1) * ga(ix-1);
    al(ix) = (P_a(ix) - be(ix-1) * ga(ix-1)) / mu(ix);
    be(ix) = P_b(ix) / mu(ix);
  }

  ga(sz-3) = P_c(sz-3) - al(sz-4) * P_e(sz-4);
  mu(sz-2) = P_d(sz-2) - be(sz-4) * P_e(sz-4) - al(sz-3) * ga(sz-3);
  al(sz-2) = (P_a(sz-2) - be(sz-3) * ga(sz-3)) / mu(sz-2);

  ga(sz-2) = P_c(sz-2) - al(sz-3) * P_e(sz-3);
  mu(sz-1) = P_d(sz-1) - be(sz-3) * P_e(sz-3) - al(sz-2) * ga(-1+sz-1);

  // seed factors
  for(int ix=0; ix<sz-2; ++ix) {
    i_mu(ix) = 1. / mu(ix);
    ga_mu(ix) = ga(ix) / mu(ix+1);
    e_mu(ix) = P_e(ix) / mu(ix+2);
  }
  i_mu(sz-2) = 1. / mu(sz-2);
  i_mu(sz-1) = 1. / mu(sz-1);
  ga_mu(sz-2) = ga(sz-2) / mu(sz-1);

  // clean up temporary arrays
  ga.DeleteAthenaArray();
  mu.DeleteAthenaArray();
}

// prepare the solver
void solver_Pentadiagonal::Prepare(
  const unsigned int dim,
  const std::vector<AthenaArray<Real> *> & bands)
{
  // reorder for algo.:
  // e<-a, c<-b, d<-c, a<-d, b<-e
  const AthenaArray<Real> * P_a {bands[3]};
  const AthenaArray<Real> * P_b {bands[4]};
  const AthenaArray<Real> * P_c {bands[1]};
  const AthenaArray<Real> * P_d {bands[2]};
  const AthenaArray<Real> * P_e {bands[0]};

  factor_computation(dim, *P_a, *P_b, *P_c, *P_d, *P_e);
  return;
}

// implementation details of the solver
void solver_Pentadiagonal::imp_Solve(
  const unsigned int axis,
  const unsigned int idx_offset,
  const AthenaArray<Real> & y,
  AthenaArray<Real> & x)
{
  // Scratch and precomputed factors
  const AthenaArray<Real> & al {_al[axis]};
  const AthenaArray<Real> & be {_be[axis]};
  const AthenaArray<Real> & ga_mu {_ga_mu[axis]};
  const AthenaArray<Real> & e_mu {_e_mu[axis]};
  const AthenaArray<Real> & i_mu {_i_mu[axis]};

  AthenaArray<Real> & z {_z[axis]};

  // Deal with different strides along different dimensions
  Real * px = &(x(0));
  const Real * py = &(const_cast<AthenaArray<Real>&>(y)(0));

  const unsigned int sz = N[axis];

  // step 3
  z(0) = py[0 + idx_offset] * i_mu(0);

  // step 4
  z(1) = py[1 * strides[axis] + idx_offset] * i_mu(1) - z(0) * ga_mu(0);

  // step 5
// #pragma omp simd
  for(int ix=2; ix < sz - 2; ++ix){
    z(ix) = py[ix * strides[axis] + idx_offset] * i_mu(ix) \
      - z(ix-2) * e_mu(ix-2) - z(ix-1) * ga_mu(ix-1);
  }

  // z[sz-3] -> z[sz-4] in first appearance cf. algo
  z(sz-2) = (py[(sz-2) * strides[axis] + idx_offset] * i_mu(sz-2)
             -z(sz-4) * e_mu(sz-4)
             -z(sz-3) * ga_mu(sz-3));
  // z[sz-2] -> z[sz-3] in first appearance cf. algo
  z(sz-1) = (py[(sz-1) * strides[axis] + idx_offset] * i_mu(sz-1)
             -z(sz-3) * e_mu(sz-3)
             -z(sz-2) * ga_mu(sz-2));

  // step 6
  px[(sz-1) * strides[axis] + idx_offset] = z(sz-1);
  px[(sz-2) * strides[axis] + idx_offset] = z(sz-2) - al(sz-2) \
    * px[(sz-1) * strides[axis] + idx_offset];

// #pragma omp simd
  for(int ix=sz-3; ix >= 0; --ix) {
    px[ix * strides[axis] + idx_offset] = z(ix) - al(ix) \
      * px[(ix+1) * strides[axis] + idx_offset] \
      - be(ix) * px[(ix+2) * strides[axis] + idx_offset];
  }
}

// wrapper to apply solver over various dimensions based on precomputed factors
void solver_Pentadiagonal::Solve(
  const unsigned int axis,
  const AthenaArray<Real> & y,
  AthenaArray<Real> & x)
{
  if (dim==3) {
    const unsigned cJ[] = {2, 2, 1};
    const unsigned cI[] = {1, 0, 0};
    unsigned int J = cJ[axis];
    unsigned int I = cI[axis];
    for(int jx=0; jx<N[cJ[axis]]; ++jx) {
      for(int ix=0; ix<N[cI[axis]]; ++ix) {
        imp_Solve(axis, ix * strides[cI[axis]] + jx * strides[cJ[axis]], y, x);
      }
    }
  } else if (dim==2) {
    const unsigned int cI[] = {1, 0}; // 1-axis

    for(int ix=0; ix<N[cI[axis]]; ++ix) {
      imp_Solve(axis, ix * strides[cI[axis]], y, x);
    }
  } else {
    imp_Solve(axis, 0, y, x);
  }
}

solver_Pentadiagonal::~solver_Pentadiagonal()
{
  for(int ix=0; ix<dim; ++ix) {
    _al[ix].DeleteAthenaArray();
    _be[ix].DeleteAthenaArray();
    _ga_mu[ix].DeleteAthenaArray();
    _e_mu[ix].DeleteAthenaArray();
    _i_mu[ix].DeleteAthenaArray();
    _z[ix].DeleteAthenaArray();
  }

  delete[] _al;
  delete[] _be;
  delete[] _ga_mu;
  delete[] _e_mu;
  delete[] _i_mu;
  delete[] _z;

  delete[] N;
  delete[] strides;
}
// ============================================================================

// ============================================================================
// Mask with a matrix of banded structure
Linear_banded_mask::Linear_banded_mask(
  const std::vector<unsigned int> & idims_N,
  const std::vector<std::vector<Real>> & icoeffL,
  const std::vector<Real> & icoeffC,
  const std::vector<std::vector<Real>> & icoeffR,
  const std::vector<std::vector<unsigned int>> & ioffset,
  const std::vector<std::vector<int>> & irix,
  const std::vector<std::vector<unsigned int>> & iwidth)
{

  dim = idims_N.size();

  N = new unsigned int [dim];
  strides = new unsigned int[dim];

  coeffL = icoeffL;
  coeffC = icoeffC;
  coeffR = icoeffR;

  offset = ioffset;
  rix = irix;
  width = iwidth;

  for(int ix=0; ix<dim; ++ix) {
    N[ix] = idims_N[ix];
  }

  strides[0] = 1;
  for(int ix=1; ix<dim; ++ix)
    strides[ix] = N[ix-1] * strides[ix-1];
}

void Linear_banded_mask::imp_Mask(
      const unsigned int axis,
      const unsigned int idx_offset,
      const AthenaArray<Real> & y,
      AthenaArray<Real> & x)
{
  // left edge --------------------------------------------
  const unsigned int ixL_l = rix[0][0];
  const unsigned int ixL_u = rix[0][1];

  for(int ix=ixL_l; ix<ixL_u; ++ix) {
    const unsigned int lix = ix-ixL_l;

    for(int six=0; six<width[0][lix]; ++six) {
      const unsigned int ix_Y = ix+six-offset[0][lix];

      x(ix * strides[axis] + idx_offset) += coeffL[lix][six] \
        * y(ix_Y * strides[axis] + idx_offset);
    }

  }
  // ------------------------------------------------------

  // center -----------------------------------------------
  const unsigned int ix_l = rix[1][0];
  const unsigned int ix_u = N[axis] + rix[1][1];

  for(int ix=ix_l; ix<ix_u; ++ix) {
    for(int six=0; six<width[1][0]; ++six) {
      const unsigned int ix_Y = ix+six-offset[1][0];

      x(ix * strides[axis] + idx_offset) += coeffC[six] \
        * y(ix_Y * strides[axis] + idx_offset);
    }
  }
  // ------------------------------------------------------

  // right edge -------------------------------------------
  const unsigned int ixR_l = rix[2][0] + N[axis];
  const unsigned int ixR_u = rix[2][1] + N[axis];

  for(int ix=ixR_l; ix<ixR_u; ++ix) {
    const unsigned int lix = ix-ixR_l;

    for(int six=0; six<width[2][lix]; ++six) {
      const unsigned int ix_Y = ix+six-offset[2][lix];

      x(ix * strides[axis] + idx_offset) += coeffR[lix][six] \
        * y(ix_Y * strides[axis] + idx_offset);
    }
  }
  // ------------------------------------------------------
}

// wrapper to apply mask over various dimensions
void Linear_banded_mask::Mask(
  const unsigned int axis,
  const AthenaArray<Real> & y,
  AthenaArray<Real> & x)
{

  if (dim==3) {
    const unsigned cJ[] = {2, 2, 1};
    const unsigned cI[] = {1, 0, 0};
    unsigned int J = cJ[axis];
    unsigned int I = cI[axis];
    for(int jx=0; jx<N[cJ[axis]]; ++jx) {
      for(int ix=0; ix<N[cI[axis]]; ++ix) {
        imp_Mask(axis, ix * strides[cI[axis]] + jx * strides[cJ[axis]], y, x);
      }
    }
  } else if (dim==2) {
    const unsigned int cI[] = {1, 0}; // 1-axis

    for(int ix=0; ix<N[cI[axis]]; ++ix) {
      imp_Mask(axis, ix * strides[cI[axis]], y, x);
    }
  } else {
    imp_Mask(axis, 0, y, x);
  }

}

Linear_banded_mask::~Linear_banded_mask()
{
  delete[] N;
  delete[] strides;
}
// ============================================================================


// ============================================================================
// General interface for 3, 5 banded solvers; uses factor precomputation
//
// Summary of conventions:
//
// 3-bands:
//  a: subdiagonal [n-1], b: diagonal [n], c: supdiagonal [n-1]
//
// 5-bands:
//  a: subsubdiagonal [n-2]
//  b: diagonal [n-1]
//  c: diagonal [n]
//  d: diagonal [n-1]
//  e: supdiagonal [n-2]
//  Note - conventions here are different to the direct analogue.
Linear_solver::Linear_solver(
  unsigned int ibandwidth,
  const std::vector<unsigned int> & idims_N)
{
  bandwidth = ibandwidth;
  dim = idims_N.size();

  if (bandwidth==3) {
    solver_Tri = new solver_Tridiagonal(idims_N);
  } else if (bandwidth==5) {
    solver_Pent = new solver_Pentadiagonal(idims_N);
  }

}

void Linear_solver::Prepare(
  const unsigned int axis,
  const std::vector<AthenaArray<Real> *> & bands)
{
  if (bandwidth == 3) {
    solver_Tri->Prepare(axis, bands);
  } else if (bandwidth == 5) {
    solver_Pent->Prepare(axis, bands);
  }
}

void Linear_solver::Solve(
  const unsigned int axis,
  const AthenaArray<Real> & y,
  AthenaArray<Real> & x)
{
  if (bandwidth == 3) {
    solver_Tri->Solve(axis, y, x);
  } else if (bandwidth == 5) {
    solver_Pent->Solve(axis, y, x);
  }
}

// destructor
Linear_solver::~Linear_solver(){
  if (bandwidth == 3) {
    delete solver_Tri;
  } else if (bandwidth == 5) {
    delete solver_Pent;
  }
}
// ============================================================================

// ============================================================================
void Linear_solver_utils::BandedAthenaArrayToVector(
  const unsigned int bandwidth,
  const AthenaArray<Real> & arr,
  std::vector<AthenaArray<Real> *> & bands)
{
  unsigned int N = arr.GetDim1();

  if (bandwidth>=3) {
    AthenaArray<Real> * M_sp, * M_d, * M_su;

    M_sp = new AthenaArray<Real>;
    M_d = new AthenaArray<Real>;
    M_su = new AthenaArray<Real>;

    M_sp->NewAthenaArray(N-1);
    M_d->NewAthenaArray(N);
    M_su->NewAthenaArray(N-1);

    for(int ix=0; ix<N; ++ix) {
      (*M_d)(ix) = arr(ix, ix);
    }

    for(int ix=0; ix<N-1; ++ix) {
      (*M_sp)(ix) = arr(ix, ix+1);
      (*M_su)(ix) = arr(ix+1, ix);
    }
    bands.push_back(M_su);
    bands.push_back(M_d);
    bands.push_back(M_sp);
  }

  if (bandwidth==5) {
    AthenaArray<Real> * M_spsp, * M_susu;

    M_spsp = new AthenaArray<Real>;
    M_susu = new AthenaArray<Real>;

    M_spsp->NewAthenaArray(N-2);
    M_susu->NewAthenaArray(N-2);

    for(int ix=0; ix<N-2; ++ix) {
      (*M_spsp)(ix) = arr(ix, ix+2);
      (*M_susu)(ix) = arr(ix+2, ix);
    }

    bands.insert(bands.begin(), M_susu);
    bands.push_back(M_spsp);
  }

  return;
}

void Linear_solver_utils::DeleteBandedVector(
    std::vector<AthenaArray<Real> *> & bands
){
  for(int ix=0; ix<bands.size(); ++ix) {
    delete bands[ix];
  }
}
// ============================================================================
