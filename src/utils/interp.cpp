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
#include "interp.hpp"
#include <math.h>
//======================================================================================
//! \file interp.cpp 
//  \brief interpolation function implementations
//======================================================================================

namespace Interpolation {
	Real LP1D(const int len, const Real *xarr, const Real *data,
			const Real x) {
		const int ix = LinearInterpIndex(len, xarr, x);
		return LP1Di(xarr, data, ix, x);

	}

	Real LP2D(const int lenx, const Real *xarr, 
			const int leny, const Real *yarr,
			const Real *data, const Real x, const Real y) {
		const int ix = LinearInterpIndex(lenx, xarr, x);
		const int iy = LinearInterpIndex(leny, yarr, y);
		return LP2Di(xarr, yarr, lenx, ix, iy, data, x, y);
	}

	Real LP1Di(const Real *xarr, const Real *data, const int ix,
			const Real x) {
		return LinearInterp(xarr[ix], xarr[ix+1], data[ix], data[ix+1], x);
	}

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
}
