#include "compact_finite_differencing.hpp"

// All the stencils
#include "stencils/CFD_1_3_3.hpp"
#include "stencils/CFD_1_5_5.hpp"
#include "stencils/CFD_2_3_3.hpp"
#include "stencils/CFD_2_5_5.hpp"

// ============================================================================
template<
  unsigned int order_,   // derivative order
  unsigned int bwL_,     // width of LHS stencil
  unsigned int bwR_,     // width of LHS stencil
  unsigned int bias_,    // relative to center idx
  int filter_            // using filtered?
>
void FDCompactStencil<order_, bwL_, bwR_, bias_, filter_>::Prepare_Banded(
  unsigned int N,
  std::vector<AthenaArray<Real> *> & M_bands,
  Real mul_factor)
{
  AthenaArray<Real> M;
  M.NewAthenaArray(N, N);

  // populate left edges ----------------------------------
  const unsigned int LixL_l = Lrix[0][0];
  const unsigned int LixL_u = Lrix[0][1];

  for(int ix=LixL_l; ix<LixL_u; ++ix) {
    const unsigned int lix = ix-LixL_l;

    for(int six=0; six<Lwidth[0][lix]; ++six) {
      const unsigned int ix_M = ix+six-Loffset[0][lix];
      M(ix, ix_M) = mul_factor * LcoeffL[lix][six];
    }
  }
  // ------------------------------------------------------

  // populate center --------------------------------------
  const unsigned int LixC_l = Lrix[1][0];
  const unsigned int LixC_u = N + Lrix[1][1];

  for(int ix=LixC_l; ix<LixC_u; ++ix) {
    for(int six=0; six<Lwidth[1][0]; ++six) {
      M(ix, ix+six-Loffset[1][0]) = mul_factor * LcoeffC[six];
    }
  }
  // ------------------------------------------------------

  // populate right edges ---------------------------------
  const unsigned int LixR_l = Lrix[2][0] + N;
  const unsigned int LixR_u = Lrix[2][1] + N;

  for(int ix=LixR_l; ix<LixR_u; ++ix) {
    const unsigned int lix = ix-LixR_l;
    for(int six=0; six<Lwidth[2][lix]; ++six) {
      const unsigned int ix_M = ix+six-Loffset[2][lix];
      M(ix, ix_M) = mul_factor * LcoeffR[lix][six];
    }
  }
  // ------------------------------------------------------

  // Extract bands ----------------------------------------
  Linear_solver_utils::BandedAthenaArrayToVector(NCFDWIDTH_L, M, M_bands);

  // Cleanup
  // Linear_solver_utils::DeleteBandedVector(M_bands);
  M.DeleteAthenaArray();

}

// ============================================================================
// Interface for using compact operators
FDCompact::FDCompact(
  const unsigned int iorder,
  const std::vector<unsigned int> & idims_N,
  const std::vector<Real> & idx)
{

  order = iorder;
  dim = idims_N.size();
  N = new unsigned int [dim];
  dx = new Real [dim];

  prepared_bm = false;

  for(int ix=0; ix<dim; ++ix) {
    N[ix] = idims_N[ix];
    dx[ix] = idx[ix];
  }

  Linear_sol = new Linear_solver(NCFDWIDTH_L, idims_N);

  if(dim==3) {
    scratch.NewAthenaArray(N[2], N[1], N[0]);
  } else if (dim==2) {
    scratch.NewAthenaArray(N[1], N[0]);
  } else {
    scratch.NewAthenaArray(N[0]);
  }


  // finalize
  for(int ix=0; ix<dim; ++ix) {
    Prepare_Banded(ix);
  }

}

// destructor
FDCompact::~FDCompact(){
  delete Linear_sol;

  if (prepared_bm)
    delete Linear_bm;

  delete[] N;
  delete[] dx;

  scratch.DeleteAthenaArray();
}


void FDCompact::Prepare_Banded(unsigned int axis)
{
  std::vector<AthenaArray<Real> *> M_bands;
  Real mul_factor = std::pow(dx[axis], order);

  if (order == 1) {
    FDCompactStencil<1, NCFDWIDTH_L, NCFDWIDTH_R, 0, CFDFILTER> FD;

    FD.Prepare_Banded(N[axis], M_bands, mul_factor);
    Linear_sol->Prepare(axis, M_bands);
    Linear_solver_utils::DeleteBandedVector(M_bands);

    if (!prepared_bm) {
      std::vector<unsigned int> dims_N;
      for(int ix=0; ix<dim; ++ix) {
        dims_N.push_back(N[ix]);
      }

      Linear_bm = new Linear_banded_mask(
        dims_N,
        FD.RcoeffL,
        FD.RcoeffC,
        FD.RcoeffR,
        FD.Roffset,
        FD.Rrix,
        FD.Rwidth);

      prepared_bm = true;
    }

  } else if (order == 2) {
    FDCompactStencil<2, NCFDWIDTH_L, NCFDWIDTH_R, 0, CFDFILTER> FD;

    FD.Prepare_Banded(N[axis], M_bands, mul_factor);
    Linear_sol->Prepare(axis, M_bands);
    Linear_solver_utils::DeleteBandedVector(M_bands);

    if (!prepared_bm) {
      std::vector<unsigned int> dims_N;
      for(int ix=0; ix<dim; ++ix) {
        dims_N.push_back(N[ix]);
      }

      Linear_bm = new Linear_banded_mask(
        dims_N,
        FD.RcoeffL,
        FD.RcoeffC,
        FD.RcoeffR,
        FD.Roffset,
        FD.Rrix,
        FD.Rwidth);

      prepared_bm = true;
    }
  }

}

void FDCompact::diff(
  const unsigned int axis,
  const AthenaArray<Real> & y,
  AthenaArray<Real> & x)
{
  // pre-fill
  scratch.Fill(0.);
  x.Fill(0.);

  // mask rhs
  Linear_bm->Mask(axis, y, scratch);
  // solve linear system
  Linear_sol->Solve(axis, scratch, x);
}