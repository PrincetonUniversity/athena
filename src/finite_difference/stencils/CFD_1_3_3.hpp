#ifndef CFD_1_3_3_HPP_
#define CFD_1_3_3_HPP_

#include "../compact_finite_differencing.hpp"

// ============================================================================
// Unfiltered: gamma[eta] = 1 (eta<=pi)
// Parameters A: M_L, M_R; B: M_L, M_R; p
// C: 1, 1; 1, 1; 3
// L: 0, 3; 0, 1; 3
// R: 3, 0; 1, 0; 3

// LHS:
// center -------------------------------------------------
template<>
const std::vector<Real> FDCompactStencil<
  1, 3, 3, 0, 0
>::LcoeffC = {
  0.25,
  1.,
  0.25
};

// left ---------------------------------------------------
template<>
const std::vector<std::vector<Real>> FDCompactStencil<
  1, 3, 3, 0, 0
>::LcoeffL = {
  {
    1.000000000000000000000000000000000000000000000000,
    3.00000000000000000000000000000000000000000000000
  }
};

// right --------------------------------------------------
template<>
const std::vector<std::vector<Real>> FDCompactStencil<
  1, 3, 3, 0, 0
>::LcoeffR = {
  {
    3.00000000000000000000000000000000000000000000000,
    1.000000000000000000000000000000000000000000000000
  }
};

// collected ----------------------------------------------
template<>
const std::vector<std::vector<unsigned int>> FDCompactStencil<
  1, 3, 3, 0, 0
>::Loffset = {
  {0}, // left
  {1}, // center
  {1}  // right
};

template<>
const std::vector<std::vector<unsigned int>> FDCompactStencil<
  1, 3, 3, 0, 0
>::Lwidth = {
  {2}, // left
  {3}, // center
  {2}  // right
};

template<>
const std::vector<std::vector<int>> FDCompactStencil<
  1, 3, 3, 0, 0
>::Lrix = {
  {0, 1}, {1, -1}, {-1, 0}  // left, center, right
};
// --------------------------------------------------------

// RHS:
// center -------------------------------------------------
template<>
const std::vector<Real> FDCompactStencil<
  1, 3, 3, 0, 0
>::RcoeffC = {
  -0.75,
  0.,
  0.75
};

// left ---------------------------------------------------
template<>
const std::vector<std::vector<Real>> FDCompactStencil<
  1, 3, 3, 0, 0
>::RcoeffL = {
  {
    -2.83333333333333333333333333333333333333333333333,
    1.500000000000000000000000000000000000000000000000,
    1.500000000000000000000000000000000000000000000000,
    -0.166666666666666666666666666666666666666666666667
  }
};

// right --------------------------------------------------
template<>
const std::vector<std::vector<Real>> FDCompactStencil<
  1, 3, 3, 0, 0
>::RcoeffR = {
  {
    0.166666666666666666666666666666666666666666666667,
    -1.500000000000000000000000000000000000000000000000,
    -1.500000000000000000000000000000000000000000000000,
    2.83333333333333333333333333333333333333333333333
  }
};

// collected ----------------------------------------------
template<>
const std::vector<std::vector<unsigned int>> FDCompactStencil<
  1, 3, 3, 0, 0
>::Roffset = {
  {0}, // left
  {1}, // center
  {3}  // right
};

template<>
const std::vector<std::vector<unsigned int>> FDCompactStencil<
  1, 3, 3, 0, 0
>::Rwidth = {
  {4}, // left
  {3}, // center
  {4}  // right
};

template<>
const std::vector<std::vector<int>> FDCompactStencil<
  1, 3, 3, 0, 0
>::Rrix = {
  {0, 1}, {1, -1}, {-1, 0}  // left, center, right
};
// --------------------------------------------------------

// ============================================================================
// Unfiltered: gamma[eta] = 1 (eta<=pi); bandwidth matched @ max p
// Parameters A: M_L, M_R; B: M_L, M_R; p
// C: 1, 1; 1, 1; 3
// L: 0, 1; 0, 1; 1
// R: 1, 0; 1, 0; 1

// LHS:
// center -------------------------------------------------
template<>
const std::vector<Real> FDCompactStencil<
  1, 3, 3, 0, -1
>::LcoeffC = {
  0.25,
  1.,
  0.25
};

// left ---------------------------------------------------
template<>
const std::vector<std::vector<Real>> FDCompactStencil<
  1, 3, 3, 0, -1
>::LcoeffL = {
  {
    1.0000000000000000000000000000000000000000000000000,
    1.0000000000000000000000000000000000000000000000000
  }
};

// right --------------------------------------------------
template<>
const std::vector<std::vector<Real>> FDCompactStencil<
  1, 3, 3, 0, -1
>::LcoeffR = {
  {
    1.0000000000000000000000000000000000000000000000000,
    1.0000000000000000000000000000000000000000000000000
  }
};

// collected ----------------------------------------------
template<>
const std::vector<std::vector<unsigned int>> FDCompactStencil<
  1, 3, 3, 0, -1
>::Loffset = {
  {0}, // left
  {1}, // center
  {1}  // right
};

template<>
const std::vector<std::vector<unsigned int>> FDCompactStencil<
  1, 3, 3, 0, -1
>::Lwidth = {
  {2}, // left
  {3}, // center
  {2}  // right
};

template<>
const std::vector<std::vector<int>> FDCompactStencil<
  1, 3, 3, 0, -1
>::Lrix = {
  {0, 1}, {1, -1}, {-1, 0}  // left, center, right
};
// --------------------------------------------------------

// RHS:
// center -------------------------------------------------
template<>
const std::vector<Real> FDCompactStencil<
  1, 3, 3, 0, -1
>::RcoeffC = {
  -0.75,
  0.,
  0.75
};

// left ---------------------------------------------------
template<>
const std::vector<std::vector<Real>> FDCompactStencil<
  1, 3, 3, 0, -1
>::RcoeffL = {
  {
    -2.000000000000000000000000000000000000000000000000,
    2.000000000000000000000000000000000000000000000000
  }
};

// right --------------------------------------------------
template<>
const std::vector<std::vector<Real>> FDCompactStencil<
  1, 3, 3, 0, -1
>::RcoeffR = {
  {
    -2.000000000000000000000000000000000000000000000000,
    2.000000000000000000000000000000000000000000000000
  }
};

// collected ----------------------------------------------
template<>
const std::vector<std::vector<unsigned int>> FDCompactStencil<
  1, 3, 3, 0, -1
>::Roffset = {
  {0}, // left
  {1}, // center
  {1}  // right
};

template<>
const std::vector<std::vector<unsigned int>> FDCompactStencil<
  1, 3, 3, 0, -1
>::Rwidth = {
  {2}, // left
  {3}, // center
  {2}  // right
};

template<>
const std::vector<std::vector<int>> FDCompactStencil<
  1, 3, 3, 0, -1
>::Rrix = {
  {0, 1}, {1, -1}, {-1, 0}  // left, center, right
};
// --------------------------------------------------------

// ============================================================================
// Filtered

#endif