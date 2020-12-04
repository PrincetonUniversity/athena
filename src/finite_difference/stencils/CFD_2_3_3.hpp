#ifndef CFD_2_3_3_HPP_
#define CFD_2_3_3_HPP_

#include "../compact_finite_differencing.hpp"

// ============================================================================
// Unfiltered: gamma[eta] = 1 (eta<=pi)
// Parameters A: M_L, M_R; B: M_L, M_R; p
// C: 1, 1; 1, 1; 3
// L: 0, 5; 0, 1; 3
// R: 5, 0; 1, 0; 3

// LHS:
// center -------------------------------------------------
template<>
const std::vector<Real> FDCompactStencil<
  2, 3, 3, 0, 0
>::LcoeffC = {
  0.10000000000000000000000000000000000000000000000000,
  1.0000000000000000000000000000000000000000000000000,
  0.10000000000000000000000000000000000000000000000000
};

// left ---------------------------------------------------
template<>
const std::vector<std::vector<Real>> FDCompactStencil<
  2, 3, 3, 0, 0
>::LcoeffL = {
  {
    1.00000000000000000000000000000000000000000000000,
    1.815478549935985000301441561195440444022829326
  }
};

// right --------------------------------------------------
template<>
const std::vector<std::vector<Real>> FDCompactStencil<
  2, 3, 3, 0, 0
>::LcoeffR = {
  {
    1.815478549935985000301441561195440444022829326,
    1.00000000000000000000000000000000000000000000000
  }
};

// collected ----------------------------------------------
template<>
const std::vector<std::vector<unsigned int>> FDCompactStencil<
  2, 3, 3, 0, 0
>::Loffset = {
  {0}, // left
  {1}, // center
  {1}  // right
};

template<>
const std::vector<std::vector<unsigned int>> FDCompactStencil<
  2, 3, 3, 0, 0
>::Lwidth = {
  {2}, // left
  {3}, // center
  {2}  // right
};

template<>
const std::vector<std::vector<int>> FDCompactStencil<
  2, 3, 3, 0, 0
>::Lrix = {
  {0, 1}, {1, -1}, {-1, 0}  // left, center, right
};
// --------------------------------------------------------

// RHS:
// center -------------------------------------------------
template<>
const std::vector<Real> FDCompactStencil<
  2, 3, 3, 0, 0
>::RcoeffC = {
  1.2000000000000000000000000000000000000000000000000,
  -2.400000000000000000000000000000000000000000000000,
  1.2000000000000000000000000000000000000000000000000
};

// left ---------------------------------------------------
template<>
const std::vector<std::vector<Real>> FDCompactStencil<
  2, 3, 3, 0, 0
>::RcoeffL = {
  {
    5.2628987916133208335845346343295337033523577717,
    -15.102681520753314583710135284827633888361869990947,
    17.228173816688004999899519479601519851992390225,
    -10.881941691741350832981651511938652815306699120,
    4.1755940583653408331826125527356131113219186703,
    -0.6820434541720012499748798699003799629980975562
  }
};

// right --------------------------------------------------
template<>
const std::vector<std::vector<Real>> FDCompactStencil<
  2, 3, 3, 0, 0
>::RcoeffR = {
  {
    -0.6820434541720012499748798699003799629980975562,
    4.1755940583653408331826125527356131113219186703,
    -10.881941691741350832981651511938652815306699120,
    17.228173816688004999899519479601519851992390225,
    -15.102681520753314583710135284827633888361869990947,
    5.2628987916133208335845346343295337033523577717
  }
};

// collected ----------------------------------------------
template<>
const std::vector<std::vector<unsigned int>> FDCompactStencil<
  2, 3, 3, 0, 0
>::Roffset = {
  {0}, // left
  {1}, // center
  {5}  // right
};

template<>
const std::vector<std::vector<unsigned int>> FDCompactStencil<
  2, 3, 3, 0, 0
>::Rwidth = {
  {6}, // left
  {3}, // center
  {6}  // right
};

template<>
const std::vector<std::vector<int>> FDCompactStencil<
  2, 3, 3, 0, 0
>::Rrix = {
  {0, 1}, {1, -1}, {-1, 0}  // left, center, right
};
// --------------------------------------------------------

// ============================================================================
// Unfiltered: gamma[eta] = 1 (eta<=pi)
// Parameters A: M_L, M_R; B: M_L, M_R; p
// C: 1, 1; 1, 1; 3
// L: 0, 1; 0, 1; 1
// R: 1, 0; 1, 0; 1

// LHS:
// center -------------------------------------------------
template<>
const std::vector<Real> FDCompactStencil<
  2, 3, 3, 0, -1
>::LcoeffC = {
  0.10000000000000000000000000000000000000000000000000,
  1.0000000000000000000000000000000000000000000000000,
  0.10000000000000000000000000000000000000000000000000
};

// left ---------------------------------------------------
template<>
const std::vector<std::vector<Real>> FDCompactStencil<
  2, 3, 3, 0, -1
>::LcoeffL = {
  {
    1.0000000000000000000000000000000000000000000000000,
    -1.0000000000000000000000000000000000000000000000000
  }
};

// right --------------------------------------------------
template<>
const std::vector<std::vector<Real>> FDCompactStencil<
  2, 3, 3, 0, -1
>::LcoeffR = {
  {
    -1.0000000000000000000000000000000000000000000000000,
    1.0000000000000000000000000000000000000000000000000
  }
};

// collected ----------------------------------------------
template<>
const std::vector<std::vector<unsigned int>> FDCompactStencil<
  2, 3, 3, 0, -1
>::Loffset = {
  {0}, // left
  {1}, // center
  {1}  // right
};

template<>
const std::vector<std::vector<unsigned int>> FDCompactStencil<
  2, 3, 3, 0, -1
>::Lwidth = {
  {2}, // left
  {3}, // center
  {2}  // right
};

template<>
const std::vector<std::vector<int>> FDCompactStencil<
  2, 3, 3, 0, -1
>::Lrix = {
  {0, 1}, {1, -1}, {-1, 0}  // left, center, right
};
// --------------------------------------------------------

// RHS:
// center -------------------------------------------------
template<>
const std::vector<Real> FDCompactStencil<
  2, 3, 3, 0, -1
>::RcoeffC = {
  1.2000000000000000000000000000000000000000000000000,
  -2.400000000000000000000000000000000000000000000000,
  1.2000000000000000000000000000000000000000000000000
};

// left ---------------------------------------------------
template<>
const std::vector<std::vector<Real>> FDCompactStencil<
  2, 3, 3, 0, -1
>::RcoeffL = {
  {
    0.,
    0.
  }
};

// right --------------------------------------------------
template<>
const std::vector<std::vector<Real>> FDCompactStencil<
  2, 3, 3, 0, -1
>::RcoeffR = {
  {
    0.,
    0.
  }
};

// collected ----------------------------------------------
template<>
const std::vector<std::vector<unsigned int>> FDCompactStencil<
  2, 3, 3, 0, -1
>::Roffset = {
  {0}, // left
  {1}, // center
  {1}  // right
};

template<>
const std::vector<std::vector<unsigned int>> FDCompactStencil<
  2, 3, 3, 0, -1
>::Rwidth = {
  {2}, // left
  {3}, // center
  {2}  // right
};

template<>
const std::vector<std::vector<int>> FDCompactStencil<
  2, 3, 3, 0, -1
>::Rrix = {
  {0, 1}, {1, -1}, {-1, 0}  // left, center, right
};
// --------------------------------------------------------

// ============================================================================
// Filtered

#endif