#ifndef INTERP_HPP
#define INTERP_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file Interp.hpp
//  \brief prototypes of linear interpolation functions
//======================================================================================

#include "../athena.hpp"

namespace Interpolation {
	// ID arrray linear interpolation
	Real LP1D(const int len, const Real *xarr, const Real *data,
			const Real x);
	// 2D array bi-linear interpolation
	Real LP2D(const int lenx, const Real *xarr, 
			const int leny, const Real *yarr,
			const Real *data, const Real x, const Real y);
	// Interpolation with index provided.
	// ix, iy: return from LinearInterpIndex 
	Real LP1Di(const Real *xarr, const Real *data, const int ix,
			const Real x);
	// 2D array bi-linear interpolation
	Real LP2Di(const Real *xarr, const Real *yarr,
			const int lenx, const int ix, const int iy,
			const Real *data, const Real x, const Real y);
	inline int LinearInterpIndex(const int len, const Real xarr[], 
			const Real x){
		int i = 0;
		if ( x < xarr[0]) {
			return 0;
		} else if ( x > xarr[len-1]) {
			return len-2;
		} else {
			for (i=0; x>xarr[i]; i++) {}
			return i-1;
		}
	}

	inline Real LinearInterp(const Real x0, const Real x1,
			const Real y0, const Real y1,
			const Real x){
		return y0 + ( (y1-y0)/(x1-x0) ) * (x-x0);
	}
}

#endif //INTERP_HPP
