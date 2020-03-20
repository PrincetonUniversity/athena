#ifndef LAGRANGE_INTERP_HPP_
#define LAGRANGE_INTERP_HPP_

//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file lagrange.hpp
//  \brief Lagrange interpolation

#include <cassert>
#include <cmath>
#include <cstring>
#include <new>

#include "athena.hpp"

// If this is uncommented always use symmetric operators
#ifndef LOCALINTERP_SYMMETRIC
#define LOCALINTERP_SYMMETRIC
#endif

// Use fuzzy comparison for whether a coord is "half-way" between known knodes
// Comparison is done based on relative error and below factor of machine eps
#ifndef LOCALINTERP_MID_FUZZY
#define LOCALINTERP_MID_FUZZY
#define LOCALINTERP_MID_FUZZY_FAC 10
#endif


template<int order>
class LagrangeInterp1D {
  public:
    LagrangeInterp1D(
        //! [in] Grid origin
        Real const origin,
        //! [in] Grid spacing
        Real const delta,
        //! [in] Number of grid points
        int siz,
        //! [in] Interpolation point
        Real const coord):
      m_origin(origin),
      m_delta(delta),
      m_delta_inv(1.0/m_delta),
      m_siz(siz),
      m_coord(coord),
      m_mid_flag(false),
      m_npoint(order + 1),
      m_out_of_bounds(false),
      point(m_point),
      npoint(m_npoint),
      out_of_bounds(m_out_of_bounds) {
      // Check if we are in the middle between two grid points
#ifdef LOCALINTERP_SYMMETRIC
      if(0 == order % 2) {
        int idx = std::lrint(std::floor((m_coord - m_origin)*m_delta_inv));

        // original conditional read:
        // m_coord - (idx*m_delta + m_origin) == 0.5 * m_delta
        Real cond_l = m_coord - (idx*m_delta + m_origin);
        Real cond_r = 0.5 * m_delta;

        bool cond;

#ifdef LOCALINTERP_MID_FUZZY
        Real eps = std::numeric_limits<Real>::epsilon();
        cond = std::abs(1 - cond_l / cond_r) < eps * LOCALINTERP_MID_FUZZY_FAC;
#else
        cond = cond_l == cond_r
#endif

        if(cond) {
          m_mid_flag = true;
          ++m_npoint;
        }
        else {
          m_mid_flag = false;
        }

      }
#endif
      // First point (from the left) of the interpolation stencil
      m_point = std::lrint(std::floor((m_coord - m_origin)*m_delta_inv
          - 0.5*(order - 1)));
#ifdef LOCALINTERP_SYMMETRIC
      if(m_mid_flag) {
          m_point -= 1;
      }
#endif

      // Shift the interpolation stencil if out of the grid
      int shift = m_point;
      if(shift < 0) {
        m_point -= shift;
        m_out_of_bounds = true;
      }
      shift = m_point + order - (m_siz - 1) + m_mid_flag;
      if(shift > 0) {
        m_point -= shift;
        m_out_of_bounds = true;
      }

      Real xp[order+2];
      for(int i = 0; i <= order; ++i) {
        xp[i] = (m_point + i) * m_delta + m_origin;
      }
#ifdef LOCALINTERP_SYMMETRIC
      if(0 == order % 2 && m_mid_flag) {
        xp[order+1] = (m_point + order + 1) * m_delta + m_origin;
      }
#endif
      m_calc_coeff_lr(xp, &m_coeff_lr[0]);
#ifdef LOCALINTERP_SYMMETRIC
      if(0 == order % 2 && m_mid_flag) {
        m_calc_coeff_rl(&xp[1], &m_coeff_rl[0]);
      }
      else {
        std::memcpy(m_coeff_rl, m_coeff_lr, sizeof(m_coeff_lr));
      }
#endif
    }

    //! Evaluates the interpolator
    template<typename T>
    T eval(
        //! [in] must be offset so that vals[0] = vals["point"]
        T const * const vals,
        //! [in] stride used to access vals
        int const stride
        ) const {

      // BD: debug
      coutBoldYellow("LagrangeInterp1D.eval:\n");
      coutBoldYellow("  (origin, delta, coord) = ");
      printf("(%1.5f, %1.5f, %1.5f)\n", m_origin, m_delta, m_coord);
      coutBoldYellow("  (m_siz, m_npoint, order) = ");
      printf("(%d, %d, %d)\n", m_siz, m_npoint, order);
      coutBoldYellow("  m_out_of_bounds = ");
      printf("%d\n", m_out_of_bounds);
      coutBoldYellow("  m_mid_flag = ");
      printf("%d\n", m_mid_flag);
      coutBoldYellow("  LOCALINTERP_SYMMETRIC = ");
#ifdef LOCALINTERP_SYMMETRIC
      printf("%d\n", true);
#else
      printf("%d\n", false);
#endif

      //-

      T out_lr = 0;
      for(int i = 0; i <= order; ++i) {
        out_lr += static_cast<T>(m_coeff_lr[i]) * vals[i*stride];
      }
#ifdef LOCALINTERP_SYMMETRIC
      T out_rl = 0;
      int shift = static_cast<int>(0 == order % 2 && m_mid_flag);
      for(int i = order; i >= 0; --i) {
        out_rl += static_cast<T>(m_coeff_rl[i]) * vals[(i + shift)*stride];
      }
      return T(0.5)*(out_lr + out_rl);
#else
      return out_lr;
#endif
    }
  private:
    // Compute the Lagrange interpolation coefficients on a given stencil
    void m_calc_coeff_lr(
        Real const * const xp,
        Real * const coeff
        ) const {
#define TYPECASE(I0, TEST, I1, OP)                                            \
      for(int j = 0; j <= order; ++j) {                                       \
        Real num = 1.0;                                                       \
        Real den = 1.0;                                                       \
        for(int i = I0; i TEST I1; OP i) {                                    \
          if(i == j) {                                                        \
            continue;                                                         \
          }                                                                   \
          num = num * (m_coord - xp[i]);                                      \
          den = den * (xp[j] - xp[i]);                                        \
        }                                                                     \
        coeff[j] = num/den;                                                   \
      }
      TYPECASE(0, <=, order, ++)
    }
    void m_calc_coeff_rl(
        Real const * const xp,
        Real * const coeff
        ) const {
      TYPECASE(order, >=, 0, --)
#undef TYPECASE
    }
  private:
    Real m_origin;
    Real m_delta;
    Real m_delta_inv;
    int m_siz;

    Real m_coord;

    // If true we have an asymmetric stencil, but the interpolation point is
    // exactly in the middle between two grid points. In this case we need to
    // average the results obtained from the interpolation on two separate
    // stencils.
    bool m_mid_flag;
    // First point (going from left to right) of the stencil
    int m_point;
    // Number of points needed for the interpolation
    int m_npoint;
    // The stencil was shifted to avoid going out of bounds
    bool m_out_of_bounds;

    // Interpolation coefficients for interpolation from left to right
    Real m_coeff_lr[order+1];
#ifdef LOCALINTERP_SYMMETRIC
    // Interpolation coefficients for interpolation from right to left
    Real m_coeff_rl[order+1];
#endif
  public:
    //! Index of the first point of the interpolation stencil
    int const & point;
    //! Number of points needed for the interpolation
    int const & npoint;
    //! The stencil
    bool const & out_of_bounds;
};

template<int ndim, int D>
class NextStencil {
  public:
    enum { value = D - 1 };
};

template<int ndim>
class NextStencil<ndim, 0> {
  public:
    enum { value = 0 };
};

// Multi-dimensional Lagrange interpolation
//
// This class assumes that
// . the grid is uniformly spaced in each direction, different grid spacing in
//   different direction is allowed;
// . the fastest running index in the data is that corrsponding to the first dimension
// . the stride used to access the data along the first dimension is 1
template<int order, int ndim>
class LagrangeInterpND {
  public:
    LagrangeInterpND(
        //! [in] Grid origin
        Real const origin[ndim],
        //! [in] Grid spacing
        Real const delta[ndim],
        //! [in] Number of grid points
        int const siz[ndim],
        //! [in] Interpolation point
        Real const coord[ndim]):
        m_out_of_bounds(false),
        out_of_bounds(m_out_of_bounds) {
      for(int d = 0; d < ndim; ++d) {
        m_origin[d] = origin[d];
        m_delta[d] = delta[d];
        m_siz[d] = siz[d];
        m_coord[d] = coord[d];
        mp_interp[d] = new (&m_interp_scratch[d][0])
          LagrangeInterp1D<order>(m_origin[d], m_delta[d],
              m_siz[d], m_coord[d]);
        m_out_of_bounds = m_out_of_bounds || mp_interp[d]->out_of_bounds;
      }
    }

    template<typename T>
    T eval(
        //! [in] Grid function to interpolate
        T const * const gf
        ) const {
      T vals[ndim][order+2];
      int pos[ndim];
      m_fill_stencil<T, ndim-1>(gf, pos, vals);
      return mp_interp[ndim-1]->eval(vals[ndim-1], 1);
    }
  private:
    // Recursively fill the stencil used for the interpolation
    template<typename T, int D>
    void m_fill_stencil(
        T const * const gf,
        int pos[ndim],
        T vals[ndim][order+2]
        ) const {
      assert(D >= 0 && D < ndim);
      if(D == 0) {
        int gidx = mp_interp[0]->point;
        int stride = 1;
        for(int d = 1; d < ndim; ++d) {
          stride *= m_siz[d-1];
          gidx += stride * (mp_interp[d]->point + pos[d]);
        }
        std::memcpy(&vals[0][0], &gf[gidx], mp_interp[0]->npoint*sizeof(T));
      }
      else {
        for(pos[D] = 0; pos[D] < mp_interp[D]->npoint; ++pos[D]) {
          m_fill_stencil<T, NextStencil<ndim, D>::value>(gf, pos, vals);
          vals[D][pos[D]] = mp_interp[D-1]->eval(vals[D-1], 1);
        }
      }
    }
  private:
    Real m_origin[ndim];
    Real m_delta[ndim];
    int m_siz[ndim];

    Real m_coord[ndim];
    bool m_out_of_bounds;

    // 1D interpolators (will be placed on the stack)
    LagrangeInterp1D<order> * mp_interp[ndim];
    // Scratch space used for placement new
    char m_interp_scratch[ndim][sizeof(LagrangeInterp1D<order>)];
  public:
    bool const & out_of_bounds;
};


#endif
