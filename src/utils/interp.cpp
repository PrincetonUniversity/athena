//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file interp.cpp
//! \brief interpolation function implementations

//this class header
#include "interp.hpp"

//c header
#include <math.h>

namespace Interpolation {
  //--------------------------------------------------------------------------------------
  //! \fn Real LP1D(const int len, const Real *xarr, const Real *data, const Real x)
  //! \brief ID arrray linear interpolation
  Real LP1D(const int len, const Real *xarr, const Real *data,
      const Real x) {
    const int ix = LinearInterpIndex(len, xarr, x);
    return LP1Di(xarr, data, ix, x);
  }

  //--------------------------------------------------------------------------------------
  //! \fn Real LP2D(const int lenx, const Real *xarr, const int leny,
  //!               const Real *yarr, const Real *data, const Real x, const Real y)
  //! \brief 2D array bi-linear interpolation
  Real LP2D(const int lenx, const Real *xarr,
      const int leny, const Real *yarr,
      const Real *data, const Real x, const Real y) {
    const int ix = LinearInterpIndex(lenx, xarr, x);
    const int iy = LinearInterpIndex(leny, yarr, y);
    return LP2Di(xarr, yarr, lenx, ix, iy, data, x, y);
  }

  //--------------------------------------------------------------------------------------
  //! \fn Real LP1Di(const Real *xarr, const Real *data, const int ix, const Real x)
  //! \brief Interpolation with index provided.
  //!
  //!  ix, iy: return from LinearInterpIndex
  Real LP1Di(const Real *xarr, const Real *data, const int ix,
      const Real x) {
    return LinearInterp(xarr[ix], xarr[ix+1], data[ix], data[ix+1], x);
  }

  //--------------------------------------------------------------------------------------
  //! \fn Real LP2Di(const Real *xarr, const Real *yarr, const int lenx,
  //!          const int ix, const int iy, const Real *data, const Real x, const Real y)
  //! \brief 2D array bi-linear interpolation with index provided
  Real LP2Di(const Real *xarr, const Real *yarr,
      const int lenx, const int ix, const int iy,
      const Real *data, const Real x, const Real y) {
    Real fl1, fl2;
    const Real x0 = xarr[ix];
    const Real x1 = xarr[ix+1];
    fl1 = LinearInterp(x0, x1, data[iy*lenx + ix], data[iy*lenx + ix+1], x);
    fl2 = LinearInterp(x0, x1, data[(iy+1)*lenx + ix], data[(iy+1)*lenx + ix+1], x);
    return LinearInterp(yarr[iy], yarr[iy+1], fl1, fl2, y);
  }
} // namespace Interpolation
