#ifndef COMPACT_FINITE_DIFFERENCING_HPP_
#define COMPACT_FINITE_DIFFERENCING_HPP_

// C headers

// C++ headers
#include<vector>     // std::vector
// #include<algorithm>  // std::reverse

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../utils/linear_solvers.hpp"

// LHS, central part of CFD
template<
  unsigned int order_,   // derivative order
  unsigned int bwL_,     // width of LHS stencil
  unsigned int bwR_,     // width of LHS stencil
  unsigned int bias_,    // relative to center idx
  int filter_            // specify filter [weight] type
>
class FDCompactStencil {
  public:
    enum {order = order_};
    enum {bwL = bwL_};
    enum {bwR = bwR_};
    enum {filter = filter_};

    // LHS: -----------------------------------------------
    // coefficient arrays
    static const std::vector<std::vector<Real>> LcoeffL;
    static const std::vector<Real> LcoeffC;
    static const std::vector<std::vector<Real>> LcoeffR;
    // Stencil base-point pos relative to diag. ix
    // i.e., where in the stencil does the diagonal appear?
    static const std::vector<std::vector<unsigned int>> Loffset;
    // Range of idxs to apply over on target vec.
    // [relative to start / end of vec.]
    // {Left, center, right} {start, end}
    static const std::vector<std::vector<int>> Lrix;
    // {Left, center, right} stencil widths
    static const std::vector<std::vector<unsigned int>> Lwidth;

    // RHS: -----------------------------------------------
    // see above for purpose
    static const std::vector<std::vector<Real>> RcoeffL;
    static const std::vector<Real> RcoeffC;
    static const std::vector<std::vector<Real>> RcoeffR;
    static const std::vector<std::vector<unsigned int>> Roffset;
    static const std::vector<std::vector<int>> Rrix;
    static const std::vector<std::vector<unsigned int>> Rwidth;

  public:
    void Prepare_Banded(
      unsigned int N,
      std::vector<AthenaArray<Real> *> & M_bands,
      Real mul_factor);
};

class FDCompact {
  public:
    FDCompact(
      const unsigned int order,
      const std::vector<unsigned int> & dims_N,
      const std::vector<Real> & dx
    );

    void Prepare_Banded(unsigned int axis);

    void diff(
      const unsigned int axis,
      const AthenaArray<Real> & y,
      AthenaArray<Real> & x
    );

    virtual ~FDCompact();

    Linear_solver * Linear_sol;
    Linear_banded_mask * Linear_bm;
  private:
    unsigned int order;
    unsigned int dim;
    unsigned int * N;
    Real * dx;
    bool prepared_bm;

    AthenaArray<Real> scratch;
};

#endif