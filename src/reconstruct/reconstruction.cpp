//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file reconstruction.cpp
//! \brief

// C headers

// C++ headers
#include <algorithm>  // max()
#include <cmath>      // abs()
#include <cstring>    // strcmp()
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../globals.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "reconstruction.hpp"

namespace {
// TODO(felker): replace these hand-rolled linear algebra routines with a real library
constexpr Real lu_tol = 3e-16;
int DoolittleLUPDecompose(Real **a, int n, int *pivot);
void DoolittleLUPSolve(Real **lu, int *pivot, Real *b, int n, Real *x);
} // namespace

//!constructor

Reconstruction::Reconstruction(MeshBlock *pmb, ParameterInput *pin) :
    characteristic_projection{false}, uniform{true, true, true},
    curvilinear{false, false},
    // read fourth-order solver switches
    correct_ic{pin->GetOrAddBoolean("time", "correct_ic", false)},
    correct_err{pin->GetOrAddBoolean("time", "correct_err", false)}, pmy_block_{pmb}
{
  // Read and set type of spatial reconstruction
  // --------------------------------
  std::string input_recon = pin->GetOrAddString("time", "xorder", "2");

  if (input_recon == "1") {
    xorder = 1;
  } else if (input_recon == "2") {
    xorder = 2;
  } else if (input_recon == "2c") {
    xorder = 2;
    characteristic_projection = true;
  } else if (input_recon == "3") {
    // PPM approximates interfaces with 4th-order accurate stencils, but use xorder=3
    // to denote that the overall scheme is "between 2nd and 4th" order w/o flux terms
    xorder = 3;
  } else if (input_recon == "3c") {
    xorder = 3;
    characteristic_projection = true;
  } else if ((input_recon == "4") || (input_recon == "4c")) {
    // Full 4th-order scheme for hydro or MHD on uniform Cartesian grids
    xorder = 4;
    if (input_recon == "4c")
      characteristic_projection = true;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in Reconstruction constructor" << std::endl
        << "xorder=" << input_recon << " not valid choice for reconstruction"<< std::endl;
    ATHENA_ERROR(msg);
  }
  // Check for incompatible choices with broader solver configuration
  // --------------------------------
  if (GENERAL_EOS && characteristic_projection) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Reconstruction constructor" << std::endl
        << "General EOS does not support characteristic reconstruction."<< std::endl;
    ATHENA_ERROR(msg);
  }
  if (GENERAL_EOS && xorder > 2) {
    std::cout << "### Warning in Reconstruction constructor" << std::endl
              << "The combination of general EOS and xorder > 2 has documented bugs"
              << " in some cases." << std::endl;
  }

  // check for necessary number of ghost zones for PPM w/o fourth-order flux corrections
  if (xorder == 3) {
    int req_nghost = 3;
    if (NGHOST < req_nghost) {
      std::stringstream msg;
      msg << "### FATAL ERROR in Reconstruction constructor" << std::endl
          << "xorder=" << input_recon <<
          " (PPM) reconstruction selected, but nghost=" << NGHOST << std::endl
          << "Reconfigure with --nghost=XXX with XXX > " << req_nghost-1 << std::endl;
      ATHENA_ERROR(msg);
    }
  }

  // perform checks of fourth-order solver configuration restrictions:
  if (xorder == 4) {
    // Uniform, Cartesian mesh with square cells (dx1f=dx2f=dx3f)
    if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
      if (pmb->block_size.x1rat != 1.0 || pmb->block_size.x2rat != 1.0 ||
          pmb->block_size.x3rat != 1.0) {
        std::stringstream msg;
        msg << "### FATAL ERROR in Reconstruction constructor" << std::endl
            << "Selected time/xorder=" << input_recon << " flux calculations"
            << " require a uniform (x1rat=x2rat=x3rat=1.0), " << std::endl
            << "Carteisan mesh with square cells. Rerun with uniform cell spacing "
            << std::endl
            << "Current values are:" << std::endl
            << std::scientific
            << std::setprecision(std::numeric_limits<Real>::max_digits10 -1)
            << "x1rat= " << pmb->block_size.x1rat << std::endl
            << "x2rat= " << pmb->block_size.x2rat << std::endl
            << "x3rat= " << pmb->block_size.x3rat << std::endl;
        ATHENA_ERROR(msg);
      }
      Real& dx_i   = pmb->pcoord->dx1f(pmb->is);
      Real& dx_j   = pmb->pcoord->dx2f(pmb->js);
      Real& dx_k   = pmb->pcoord->dx3f(pmb->ks);
      // Note, probably want to make the following condition less strict (signal warning
      // for small differences due to floating-point issues) but upgrade to error for
      // large deviations from a square mesh. Currently signals a warning for each
      // MeshBlock with non-square cells.
      if ((pmb->block_size.nx2 > 1 && dx_i != dx_j) ||
          (pmb->block_size.nx3 > 1 && dx_j != dx_k)) {
        // It is possible for small floating-point differences to arise despite equal
        // analytic values for grid spacings in the coordinates.cpp calculation of:
        // Real dx=(block_size.x1max-block_size.x1min)/(ie-is+1);
        // due to the 3x rounding operations in numerator, e.g.
        // float(float(x1max) - float((x1min))
        // if mesh/x1max != mesh/x2max, etc. and/or if an asymmetric MeshBlock
        // decomposition is used
        if (Globals::my_rank == 0) {
          // std::stringstream msg;
          std::cout
              << "### Warning in Reconstruction constructor" << std::endl
              << "Selected time/xorder=" << input_recon << " flux calculations"
              << " require a uniform, Carteisan mesh with" << std::endl
              << "square cells (dx1f=dx2f=dx3f). "
              << "Change mesh limits and/or number of cells for equal spacings\n"
              << "Current values are:" << std::endl
              << std::scientific
              << std::setprecision(std::numeric_limits<Real>::max_digits10 - 1)
              << "dx1f=" << dx_i << std::endl
              << "dx2f=" << dx_j << std::endl
              << "dx3f=" << dx_k << std::endl;
          // ATHENA_ERROR(msg);
        }
      }
      if (pmb->pmy_mesh->multilevel) {
        std::stringstream msg;
        msg << "### FATAL ERROR in Reconstruction constructor" << std::endl
            << "Selected time/xorder=" << input_recon << " flux calculations"
            << " currently does not support SMR/AMR " << std::endl;
        ATHENA_ERROR(msg);
      }
    } else {
      std::stringstream msg;
      msg << "### FATAL ERROR in Reconstruction constructor" << std::endl
          << "Specified COORDINATE_SYSTEM=" << COORDINATE_SYSTEM << "\n"
          << "is incompatible with selected time/xorder=" << input_recon << std::endl
          << "Reconfigure with Cartesian coordinates " << std::endl;
      ATHENA_ERROR(msg);
    }

    if (pmb->pmy_mesh->shear_periodic) {
      std::stringstream msg;
      msg << "### FATAL ERROR in Reconstruction constructor" << std::endl
          << "Selected time/xorder=" << input_recon << " flux calculations"
          << "currently does not support shearing box boundary conditions " << std::endl;
      ATHENA_ERROR(msg);
      return;
    }

    // check for necessary number of ghost zones for PPM w/ fourth-order flux corrections
    int req_nghost = 4;
    // until new algorithm for face-averaged Field->bf to cell-averaged Hydro->bcc
    // conversion is added, NGHOST>=6
    if (MAGNETIC_FIELDS_ENABLED)
      req_nghost += 2;
    if (NGHOST < req_nghost) {
      std::stringstream msg;
      msg << "### FATAL ERROR in Reconstruction constructor" << std::endl
          << "time/xorder=" << input_recon
          << " reconstruction selected, but nghost=" << NGHOST << std::endl
          << "Reconfigure with --nghost=XXX with XXX > " << req_nghost-1 << std::endl;
      ATHENA_ERROR(msg);
    }
  }

  // switch to secondary PLM and PPM variants for nonuniform and/or curvilinear meshes
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    // cylindrical: x1=r requires special treatment. x2=phi, x3=z do not.
    curvilinear[X1DIR] = true;
  }
  if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    // spherical_polar: x1=r and x2=theta require special treatment. x3=phi does not
    curvilinear[X1DIR] = true;
    curvilinear[X2DIR] = true;
  }
  // for all coordinate systems, nonuniform geometric spacing or user-defined
  // MeshGenerator ---> use nonuniform reconstruction weights and limiter terms
  if (pmb->block_size.x1rat != 1.0)
    uniform[X1DIR] = false;
  if (pmb->block_size.x2rat != 1.0)
    uniform[X2DIR] = false;
  if (pmb->block_size.x3rat != 1.0)
    uniform[X3DIR] = false;

  // Uniform mesh with --coord=cartesian or GR: Minkowski, Schwarzschild, Kerr-Schild,
  // GR-User will use the uniform Cartesian limiter and reconstruction weights
  // TODO(c-white): use modified version of curvilinear PPM reconstruction weights and
  // limiter formulations for Schwarzschild, Kerr metrics instead of Cartesian-like wghts

  // Allocate memory for scratch arrays used in PLM and PPM
  int nc1 = pmb->ncells1;
  scr01_i_.NewAthenaArray(nc1);
  scr02_i_.NewAthenaArray(nc1);

  scr1_ni_.NewAthenaArray(std::max(NWAVE, NSCALARS), nc1);
  scr2_ni_.NewAthenaArray(std::max(NWAVE, NSCALARS), nc1);
  scr3_ni_.NewAthenaArray(std::max(NWAVE, NSCALARS), nc1);
  scr4_ni_.NewAthenaArray(std::max(NWAVE, NSCALARS), nc1);

  if ((xorder == 3) || (xorder == 4)) {
    Coordinates *pco = pmb->pcoord;
    scr03_i_.NewAthenaArray(nc1);
    scr04_i_.NewAthenaArray(nc1);
    scr05_i_.NewAthenaArray(nc1);
    scr06_i_.NewAthenaArray(nc1);
    scr07_i_.NewAthenaArray(nc1);
    scr08_i_.NewAthenaArray(nc1);
    scr09_i_.NewAthenaArray(nc1);
    scr10_i_.NewAthenaArray(nc1);
    scr11_i_.NewAthenaArray(nc1);
    scr12_i_.NewAthenaArray(nc1);
    scr13_i_.NewAthenaArray(nc1);
    scr14_i_.NewAthenaArray(nc1);

    scr5_ni_.NewAthenaArray(std::max(NWAVE, NSCALARS), nc1);
    scr6_ni_.NewAthenaArray(std::max(NWAVE, NSCALARS), nc1);
    scr7_ni_.NewAthenaArray(std::max(NWAVE, NSCALARS), nc1);
    scr8_ni_.NewAthenaArray(std::max(NWAVE, NSCALARS), nc1);

    // Precompute PPM coefficients in x1-direction ---------------------------------------
    c1i.NewAthenaArray(nc1);
    c2i.NewAthenaArray(nc1);
    c3i.NewAthenaArray(nc1);
    c4i.NewAthenaArray(nc1);
    c5i.NewAthenaArray(nc1);
    c6i.NewAthenaArray(nc1);
    hplus_ratio_i.NewAthenaArray(nc1);
    hminus_ratio_i.NewAthenaArray(nc1);

    // Greedily allocate tiny 4x4 matrix + 4x1 vectors (RHS, solution, and permutation
    // indices) in case PPMx1 and/or PPMx2 require them for computing the curvilinear
    // coorddinate reconstruction weights. Same data structures are reused at each spatial
    // index (i or j) and for both PPMx1 and PPMx2 weight calculations:
    constexpr int kNrows = 4;       // = [i-i_L, i+i_R] stencil of reconstruction
    constexpr int kNcols = 4;       // = [0, p-1], p=order of reconstruction
    // system in Mignone equation 21
    Real **beta = new Real*[kNrows];
    for (int i=0; i<kNrows; ++i) {
      beta[i] = new Real[kNcols];
    }

    Real w_sol[kNrows], b_rhs[kNrows];
    int permute[kNrows];
    int m_coord = 0;    // power of simple Jacobian term;  0 ---> Cartesian coordinates
    if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
      m_coord = 1;
    } else { // if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
      m_coord = 2;
    }
    if (curvilinear[X1DIR]) {
      for (int i=(pmb->is)-NGHOST+2; i<=(pmb->ie)+NGHOST-1; ++i) {
        // nonuniform Beta matrix entries reach uppermost face pmb->ncells1+1
        if (!uniform[X1DIR]) {  // nonuniform spacing in curvilinear coordinate
          for (int row=0; row<kNrows; row++) {
            b_rhs[row] = std::pow(pco->x1f(i), row);
            for (int col=0; col<kNcols; col++) {
              // Mignone equation 23: simplification of entries in \beta matrix when
              // Jacobian is a simple power of the coordinate, e.g. radial coordinate:
              Real coeff = (m_coord + 1.0)/(col + m_coord + 1.0);
              int pow_num = col + m_coord + 1;
              int pow_denom = m_coord + 1;
              int s = row - kNrows/2;  // rescale index from -2 to +1
              // initializing transpose:
              beta[col][row] =
                  coeff*(
                      std::pow(pco->x1f(i+s+1), pow_num)
                      - std::pow(pco->x1f(i+s), pow_num)) / (
                          std::pow(pco->x1f(i+s+1), pow_denom)
                          - std::pow(pco->x1f(i+s), pow_denom));
              // TODO(kfelker): add NaN detection to matrix entries. They are challenging
              // to compute accurately w/ finite-precision arithmetic.
              // Both the numerator and denominator consist of differences of nearby x1f
              // values (that get closer the as resolution increases) raised to the same
              // power, where the numerator has pow_num=3:6 and pow_denom=3 in
              // spherical-polar coordinates, e.g.
            }
          }

          if (DoolittleLUPDecompose(beta, 4, permute)) {
            // LUP decomposition succeeded; solve the linear system w/ forward/back sub:
            DoolittleLUPSolve(beta, permute, b_rhs, 4, w_sol);
          } else {
            std::stringstream msg;
            msg << "### FATAL ERROR in Reconstruction constructor" << std::endl
                << "Failure in LUP decomposition step for computing PPMx1 curvilinear\n"
                << "reconstruction weights at cell index i=" << i << " for nx1="
                << pmb->block_size.nx1 << std::endl;
            std::cout << msg.str();
            ATHENA_ERROR(msg);
          }
          // TODO(felker): probably need to change the following condition to ensure that
          // the ghost cells have negative radius. Add conditional for r < 0
          if (i < pmb->is) {
            // cant simply take abs(r) for negative radius ghost cells when we have to
            // solve a matrix-vector system of equations to get the coefficients, since
            // that introduces a singular matrix. Highest indexed ghost cell (with upper
            // face at r=0) appear to have same radial face locations (reversed) as the
            // first real cell. instead, reverse the coefficients:
            c1i(i) = w_sol[3];
            c2i(i) = w_sol[2];
            c3i(i) = w_sol[1];
            c4i(i) = w_sol[0];
          } else {
            c1i(i) = w_sol[0];
            c2i(i) = w_sol[1];
            c3i(i) = w_sol[2];
            c4i(i) = w_sol[3];
          }
        } else { // if (uniform[X1DIR]) {
          // Mignone section 2.2: conservative reconstruction from volume averages
          Real io = std::abs(i - pmy_block_->is);  // il=is-1 ---> io = 2
          // Notes:
          // - io (i offset) must be floating-point, not integer type. io^4 and io^8 terms
          // in below lines quickly cause overflows of 32-bit and 64-bit integer limtis in
          // the intermediate calculations of RHS expressions
          // - io=1 should correspond to "is" (first real cell face)
          // - take absolute value to handle lower x1 ghost zone cell faces properly
          if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
            Real delta = 120.0*SQR(SQR(io)) - 360.0*SQR(io) + 96.0;
            // Mignone equation B.9:
            // w_im1
            c1i(i) = -(2.0*io - 3.0)*(5.0*SQR(io)*io + 8.0*SQR(io) - 3.0*io - 4.0)/delta;
            // w_i
            c2i(i) = (2.0*io - 1.0)*(35.0*SQR(io)*io + 24.0*SQR(io)
                                     - 93.0*io - 60.0)/delta;
            // w_ip1
            c3i(i) = (2.0*io + 1.0)*(35.0*SQR(io)*io - 24.0*SQR(io)
                                     - 93.0*io + 60.0)/delta;
            // w_ip2
            c4i(i) = -(2.0*io + 3.0)*(5.0*SQR(io)*io - 8.0*SQR(io) - 3.0*io + 4.0)/delta;
          } else { // if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
            // Mignone equation B.14:
            Real delta = 36.0*(15.0*SQR(SQR(SQR(io))) - 85.0*SQR(SQR(io))*SQR(io)
                               + 150.0*SQR(SQR(io)) - 60.0*SQR(io) + 16);
            c1i(i) = -(3.0*SQR(io) - 9.0*io + 7.0)*(
                15.0*SQR(SQR(io))*SQR(io) + 48.0*SQR(SQR(io))*io
                + 23.0*SQR(SQR(io)) - 48.0*SQR(io)*io - 30.0*SQR(io)
                + 16.0*io + 12.0)/delta;
            c2i(i) = (3.0*SQR(io) - 3.0*io + 1.0)*(
                105.0*SQR(SQR(io))*SQR(io) + 144.0*SQR(SQR(io))*io
                - 487.0*SQR(SQR(io)) - 720.0*SQR(io)*io +510.0*SQR(io)
                + 1008.0*io + 372.0)/delta;
            c3i(i) = (3.0*SQR(io) + 3*io + 1.0)*(
                105.0*SQR(SQR(io))*SQR(io) - 144.0*SQR(SQR(io))*io
                - 487.0*SQR(SQR(io)) + 720.0*SQR(io)*io +510.0*SQR(io)
                - 1008.0*io + 372.0)/delta;
            c4i(i) = -(3.0*SQR(io) + 9.0*io + 7.0)*(
                15.0*SQR(SQR(io))*SQR(io) - 48.0*SQR(SQR(io))*io
                + 23.0*SQR(SQR(io)) + 48.0*SQR(io)*io - 30.0*SQR(io)
                - 16.0*io + 12.0)/delta;
          } // end "spherical_polar"
          // TODO(felker): add check for normalization condition, Mignone eq 22.
          // (typical deviations of sum(wghts)!=1.0 are around machine precision)
        }   // end "uniform[X1DIR]"
      }  // end loop over i

      // Compute curvilinear geometric factors for limiter (Mignone eq 48): radial
      // direction in cylindrical and spherical-polar coordinates. Same formulas
      // for nonuniform and uniform radial mesh spacings.
      for (int i=(pmb->is)-1; i<=(pmb->ie)+1; ++i) {
        Real h_plus, h_minus;
        Real& dx_i   = pco->dx1f(i);
        Real& xv_i   = pco->x1v(i);

        // radius may beomce negative in the lower x1 ghost cells:
        xv_i = 0.5*(pco->x1f(i+1) + pco->x1f(i));
        xv_i = std::abs(xv_i);

        if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          // cylindrical radial coordinate
          h_plus = 3.0 + dx_i/(2.0*xv_i);
          h_minus = 3.0 - dx_i/(2.0*xv_i);
        } else { // if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
          // spherical radial coordinate
          h_plus = 3.0 + (2.0*dx_i*(10.0*xv_i + dx_i))/(20.0*SQR(xv_i) + SQR(dx_i));
          h_minus = 3.0 + (2.0*dx_i*(-10.0*xv_i + dx_i))/(20.0*SQR(xv_i) + SQR(dx_i));
        }
        hplus_ratio_i(i) = (h_plus + 1.0)/(h_minus - 1.0);
        hminus_ratio_i(i) = (h_minus + 1.0)/(h_plus - 1.0);
      }
    } else { // Cartesian-like x2 coordinate
      // zero-curvature PPM limiter does not depend on mesh uniformity:
      for (int i=(pmb->is)-1; i<=(pmb->ie)+1; ++i) {
        // h_plus = 3.0;
        // h_minus = 3.0;
        // Ratios are = 2 for Cartesian coords, as in the original PPM limiter's
        // overshoot conditions
        hplus_ratio_i(i) = 2.0;
        hminus_ratio_i(i) = 2.0;
      }
      // 4th order reconstruction weights along Cartesian-like x1 w/ uniform spacing
      if (uniform[X1DIR]) {
#pragma omp simd
        for (int i=(pmb->is)-NGHOST; i<=(pmb->ie)+NGHOST; ++i) {
          // reducing general formula in ppm.cpp corresonds to Mignone eq B.4 weights:
          // (-1/12, 7/12, 7/12, -1/12)
          c1i(i) = 0.5;
          c2i(i) = 0.5;
          c3i(i) = 0.5;
          c4i(i) = 0.5;
          c5i(i) = 1.0/6.0;
          c6i(i) = -1.0/6.0;
        }
      } else { // coeffcients along Cartesian-like x1 with nonuniform mesh spacing
#pragma omp simd
        for (int i=(pmb->is)-NGHOST+1; i<=(pmb->ie)+NGHOST-1; ++i) {
          Real& dx_im1 = pco->dx1f(i-1);
          Real& dx_i   = pco->dx1f(i  );
          Real& dx_ip1 = pco->dx1f(i+1);
          Real qe = dx_i/(dx_im1 + dx_i + dx_ip1);       // Outermost coeff in CW eq 1.7
          c1i(i) = qe*(2.0*dx_im1+dx_i)/(dx_ip1 + dx_i); // First term in CW eq 1.7
          c2i(i) = qe*(2.0*dx_ip1+dx_i)/(dx_im1 + dx_i); // Second term in CW eq 1.7
          if (i > (pmb->is)-NGHOST+1) {  // c3-c6 are not computed in first iteration
            Real& dx_im2 = pco->dx1f(i-2);
            Real qa = dx_im2 + dx_im1 + dx_i + dx_ip1;
            Real qb = dx_im1/(dx_im1 + dx_i);
            Real qc = (dx_im2 + dx_im1)/(2.0*dx_im1 + dx_i);
            Real qd = (dx_ip1 + dx_i)/(2.0*dx_i + dx_im1);
            qb = qb + 2.0*dx_i*qb/qa*(qc-qd);
            c3i(i) = 1.0 - qb;
            c4i(i) = qb;
            c5i(i) = dx_i/qa*qd;
            c6i(i) = -dx_im1/qa*qc;
          }
        }
      }
    } // end !curvilinear[X1DIR]

    // Precompute PPM coefficients in x2-direction ---------------------------------------
    if (pmb->block_size.nx2 > 1) {
      int nc2 = pmb->ncells2;
      c1j.NewAthenaArray(nc2);
      c2j.NewAthenaArray(nc2);
      c3j.NewAthenaArray(nc2);
      c4j.NewAthenaArray(nc2);
      c5j.NewAthenaArray(nc2);
      c6j.NewAthenaArray(nc2);
      hplus_ratio_j.NewAthenaArray(nc2);
      hminus_ratio_j.NewAthenaArray(nc2);

      if (curvilinear[X2DIR]) {
        // meridional/polar/theta direction in spherical-polar coordinates
        std::function<int(int)> factorial;  // recursive lambda function
                                            // (C++14 would make this easier)
        factorial = [&factorial] (int k) -> int { return (k <=1) ? 1 : factorial(k-1)*k;};
        for (int j=(pmb->js)-NGHOST+2; j<=(pmb->je)+NGHOST-1; ++j) {
          // need to numerically invert system of equations given by Mignone
          // equations 21 and 27 with specific nonuniform or uniform meridional mesh
          // spacings in order to precompute reconstruction weights
          for (int row=0; row<kNrows; row++) {
            b_rhs[row] = std::pow(pco->x2f(j), row);
            int s = row - kNrows/2;  // rescale index from -2 to +1
            Real coeff = 1.0/(std::cos(pco->x2f(j+s)) - std::cos(pco->x2f(j+s+1)));
            for (int col=0; col<kNcols; col++) {
              beta[col][row] = 0.0;
              for (int k=0; k<=col; k++) {
                // Mignone equation 27:
                // transpose:
                beta[col][row] += factorial(col)/(factorial(col - k))*(
                    std::pow(pco->x2f(j+s), col-k)
                    *std::cos(pco->x2f(j+s) + 0.5*PI*k)
                    - std::pow(pco->x2f(j+s+1), col-k)
                    *std::cos(pco->x2f(j+s+1) + 0.5*PI*k));
              }
              beta[col][row] *= coeff;
              // TODO(felker): add NaN detection for beta matrix entries
            }
          }

          if (DoolittleLUPDecompose(beta, 4, permute)) {
            // LUP decomposition succeeded; solve the linear system w/ forward/back sub:
            DoolittleLUPSolve(beta, permute, b_rhs, 4, w_sol);
          } else {
            std::stringstream msg;
            msg << "### FATAL ERROR in Reconstruction constructor" << std::endl
                << "Failure in LUP decomposition step for computing PPMx2 curvilinear\n"
                << "reconstruction weights at cell index j=" << j << " for nx2="
                << pmb->block_size.nx2 << std::endl;
            std::cout << msg.str();
            ATHENA_ERROR(msg);
          }
          // TODO(felker): possibly transform weights at theta <0 and/or theta >pi
          c1j(j) = w_sol[0];
          c2j(j) = w_sol[1];
          c3j(j) = w_sol[2];
          c4j(j) = w_sol[3];
        }
        // Compute curvilinear geometric factors for limiter (Mignone eq 48)
        for (int j=(pmb->js)-1; j<=(pmb->je)+1; ++j) {
          // corrections to PPMx2 only for spherical-polar coordinates
          // if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
          Real h_plus, h_minus;
          // note: x2 may become negative at lower boundary ghost cells and may exceedc
          // pi at upper boundary ghost cells
          // TODO(felker): may need to take abs() or change signs of terms in ghost zone
          Real& dx_j   = pco->dx2f(j);
          Real& xf_j   = pco->x2f(j);
          Real& xf_jp1   = pco->x2f(j+1);
          Real dmu = std::cos(xf_j) - std::cos(xf_jp1);
          Real dmu_tilde = std::sin(xf_j) - std::sin(xf_jp1);
          h_plus = (dx_j*(dmu_tilde + dx_j*std::cos(xf_jp1)))/(
              dx_j*(std::sin(xf_j) + std::sin(xf_jp1)) - 2.0*dmu);
          h_minus = -(dx_j*(dmu_tilde + dx_j*std::cos(xf_j)))/(
              dx_j*(std::sin(xf_j) + std::sin(xf_jp1)) - 2.0*dmu);
          hplus_ratio_j(j) = (h_plus + 1.0)/(h_minus - 1.0);
          hminus_ratio_j(j) = (h_minus + 1.0)/(h_plus - 1.0);
        }
      } else { // Cartesian-like x2 coordinate
        // zero-curvature PPM limiter does not depend on mesh uniformity:
        for (int j=(pmb->js)-1; j<=(pmb->je)+1; ++j) {
          // h_plus = 3.0;
          // h_minus = 3.0;
          // Ratios are = 2 for Cartesian coords, as in the original PPM limiter's
          // overshoot conditions
          hplus_ratio_j(j) = 2.0;
          hminus_ratio_j(j) = 2.0;
        }
        // 4th order reconstruction weights along Cartesian-like x2 w/ uniform spacing
        if (uniform[X2DIR]) {
#pragma omp simd
          for (int j=(pmb->js)-NGHOST; j<=(pmb->je)+NGHOST; ++j) {
            c1j(j) = 0.5;
            c2j(j) = 0.5;
            c3j(j) = 0.5;
            c4j(j) = 0.5;
            c5j(j) = 1.0/6.0;
            c6j(j) = -1.0/6.0;
          }
        } else { // coeffcients along Cartesian-like x2 with nonuniform mesh spacing
#pragma omp simd
          for (int j=(pmb->js)-NGHOST+2; j<=(pmb->je)+NGHOST-1; ++j) {
            Real& dx_jm1 = pco->dx2f(j-1);
            Real& dx_j   = pco->dx2f(j  );
            Real& dx_jp1 = pco->dx2f(j+1);
            Real qe = dx_j/(dx_jm1 + dx_j + dx_jp1);       // Outermost coeff in CW eq 1.7
            c1j(j) = qe*(2.0*dx_jm1 + dx_j)/(dx_jp1 + dx_j); // First term in CW eq 1.7
            c2j(j) = qe*(2.0*dx_jp1 + dx_j)/(dx_jm1 + dx_j); // Second term in CW eq 1.7

            if (j > (pmb->js)-NGHOST+1) {  // c3-c6 are not computed in first iteration
              Real& dx_jm2 = pco->dx2f(j-2);
              Real qa = dx_jm2 + dx_jm1 + dx_j + dx_jp1;
              Real qb = dx_jm1/(dx_jm1 + dx_j);
              Real qc = (dx_jm2 + dx_jm1)/(2.0*dx_jm1 + dx_j);
              Real qd = (dx_jp1 + dx_j)/(2.0*dx_j + dx_jm1);
              qb = qb + 2.0*dx_j*qb/qa*(qc-qd);
              c3j(j) = 1.0 - qb;
              c4j(j) = qb;
              c5j(j) = dx_j/qa*qd;
              c6j(j) = -dx_jm1/qa*qc;
            }
          }
        } // end nonuniform Cartesian-like
      } // end !curvilinear[X2DIR]
    } // end 2D or 3D

    // Precompute PPM coefficients in x3-direction
    if (pmb->block_size.nx3 > 1) {
      int nc3 = pmb->ncells3;
      c1k.NewAthenaArray(nc3);
      c2k.NewAthenaArray(nc3);
      c3k.NewAthenaArray(nc3);
      c4k.NewAthenaArray(nc3);
      c5k.NewAthenaArray(nc3);
      c6k.NewAthenaArray(nc3);
      hplus_ratio_k.NewAthenaArray(nc3);
      hminus_ratio_k.NewAthenaArray(nc3);

      // reconstruction coeffiencients in x3, Cartesian-like coordinate:
      if (uniform[X3DIR]) { // uniform spacing
#pragma omp simd
        for (int k=(pmb->ks)-NGHOST; k<=(pmb->ke)+NGHOST; ++k) {
          c1k(k) = 0.5;
          c2k(k) = 0.5;
          c3k(k) = 0.5;
          c4k(k) = 0.5;
          c5k(k) = 1.0/6.0;
          c6k(k) = -1.0/6.0;
        }

      } else { // nonuniform spacing
#pragma omp simd
        for (int k=(pmb->ks)-NGHOST+2; k<=(pmb->ke)+NGHOST-1; ++k) {
          Real& dx_km1 = pco->dx3f(k-1);
          Real& dx_k   = pco->dx3f(k  );
          Real& dx_kp1 = pco->dx3f(k+1);
          Real qe = dx_k/(dx_km1 + dx_k + dx_kp1);       // Outermost coeff in CW eq 1.7
          c1k(k) = qe*(2.0*dx_km1+dx_k)/(dx_kp1 + dx_k); // First term in CW eq 1.7
          c2k(k) = qe*(2.0*dx_kp1+dx_k)/(dx_km1 + dx_k); // Second term in CW eq 1.7

          if (k > (pmb->ks)-NGHOST+1) {  // c3-c6 are not computed in first iteration
            Real& dx_km2 = pco->dx3f(k-2);
            Real qa = dx_km2 + dx_km1 + dx_k + dx_kp1;
            Real qb = dx_km1/(dx_km1 + dx_k);
            Real qc = (dx_km2 + dx_km1)/(2.0*dx_km1 + dx_k);
            Real qd = (dx_kp1 + dx_k)/(2.0*dx_k + dx_km1);
            qb = qb + 2.0*dx_k*qb/qa*(qc-qd);
            c3k(k) = 1.0 - qb;
            c4k(k) = qb;
            c5k(k) = dx_k/qa*qd;
            c6k(k) = -dx_km1/qa*qc;
          }
        }
        // Compute geometric factors for x3 limiter (Mignone eq 48)
        // (no curvilinear corrections in x3)
        for (int k=(pmb->ks)-1; k<=(pmb->ke)+1; ++k) {
          // h_plus = 3.0;
          // h_minus = 3.0;
          // Ratios are both = 2 for Cartesian and all curviliniear coords
          hplus_ratio_k(k) = 2.0;
          hminus_ratio_k(k) = 2.0;
        }
      }
    }
    for (int i=0; i<kNrows; ++i) {
      delete[] beta[i];
    }
    delete[] beta;
  } // end "if PPM or full 4th order spatial integrator"
}


namespace {

//----------------------------------------------------------------------------------------
//! \fn void DoolittleLUPDecompose(Real **a, int n, int *pivot)
//! \brief perform LU decomposition with partial (row) pivoting using Doolittle's
//! algorithm. Partial pivoting is required for stability.
//!
//! Let D be a diagonal matrix, L be a unit lower triangular matrix (main diagonal is all
//! 1's), and U be a unit upper triangular matrix
//! Crout = (LD)U  ---> unit upper triangular U and L'=LD non-unit lower triangular
//! Doolittle = L(DU) ---> unit lower triangular L and U'=DU non-unit upper triangular
//!
//! INPUT:
//!   - a: square nxn matrix A of real numbers. Must be a mutable pointer-to-pointer/rows.
//!   - n: number of rows and columns in "a", rows in pivot output
//!
//! Also expects "const Real lu_tol >=0" file-scope variable to be defined = criterion
//! for detecting degenerate input "a" (or nearly-degenerate).
//!
//! OUTPUT:
//!   - a: modified in-place to contain both lower- and upper-triangular matrices L, U
//!     as A <- L + U (the 1's on the diagonal of L are not stored) in the decomposition
//!     PA=LU. See NR pg 50; even though they claim to use Crout, they are probably
//!     use Doolittle. They assume unit diagonal in Lx=Pb forward substitution.
//!   - pivot: nx1 int vector that is a sparse representation of the nxn permutation
//!     matrix P. For each row/vector entry, the value = the column # of the nonzero
//!     pivot element
//!
//! RETURN:
//!  - failure=0: routine detected that "a" matrix was nearly-singular
//!  - success=1: LUP decomposition completed
//!
//! Both "a", "pivot" can then be passed with RHS vector "b" to DoolittleLUPSolve in
//! order to solve Ax=b system of linear equations
//!
//! REFERENCES:
//!   - Numerical Recipes, 3rd ed. (NR) section 2.3 "LU Decomposition & its Applications"

int DoolittleLUPDecompose(Real **a, int n, int *pivot) {
  constexpr int failure = 0, success = 1;
  // initialize unit permutation matrix P=I
  for (int i=0; i<n; i++)
    pivot[i] = i;
  // In our sparse representation, could let pivot be (n+1)x1 and init. pivot[n]=n,
  // increment for each pivot so that it equals n+s upon return, s=# permutations
  // Useful for determinant calculation.

  // loop over rows of input matrix:
  for (int i=0; i<n; i++) {
    Real a_max = 0.0, a_abs = 0.0;
    int i_max = i;
    // search for largest pivot element, located at row i_max
    for (int k=i; k<n; k++) {
      a_abs = std::abs(a[k][i]);
      if (a_abs > a_max) { // found larger pivot element
        a_max = a_abs;
        i_max = k;
      }
    }

    // if the pivot element is near zero, the matrix is likely singular
    if (a_max < lu_tol) {  // 0.0) { // see NR comment in ludcmp.h
      // do not divide by 0
      std::cout << std::scientific
                << std::setprecision(std::numeric_limits<Real>::max_digits10 -1)
                << "DoolittleLUPDecompose detects singular matrix with\n"
                << "pivot element=" << a_max << " < tolerance=" << lu_tol << std::endl;
      return failure;
    }

    if (i != i_max) {  // need to pivot rows:
      // pivoting "pivot" vector
      int row_idx = pivot[i];
      pivot[i] = pivot[i_max];
      pivot[i_max] = row_idx;

      // pivoting rows of A
      Real *pivot_ptr = a[i];
      a[i] = a[i_max];
      a[i_max] = pivot_ptr;
    }

    // these lines are the only difference from Crout's in-place approach w/ pivoting
    for (int j=i+1; j<n; j++) { // loop over rows; NR has the same approach as here
      // fill lower triangular matrix L elements at column "i":
      a[j][i] /= a[i][i];
      // (Crout finds upper triangular matrix U elemens at row "i" in this step)
      for (int k=i+1; k<n; k++) // update remaining submatrix
        a[j][k] -= a[j][i]*a[i][k];
    }
  }
  // in-place LU factorization with partial pivoting is complete
  return success;
}


//----------------------------------------------------------------------------------------
//! \fn void DoolittleLUPSolve(Real **lu, int *pivot, Real *b, int n, Real *x)
//! \brief after DoolittleLUPDecompose() function has transformed input the LHS of Ax=b
//! system to partially-row pivoted, LUP decomposed equivalent PAx=LUx=Pb, solve for x
//!
//! INPUT:
//!  - lu: square nxn matrix of real numbers containing output "a" of successful
//!    DoolittleLUPDecompose() function call. See notes in that function for details.
//!  - pivot: nx1 vector of integers produced by DoolittleLUPDecompose()
//!  - b: RHS column vector of real numbers in original Ax=b system of linear equations
//!
//! OUTPUT:
//!  - x: nx1 column vector of real numbers containing solution in original Ax=b system

void DoolittleLUPSolve(Real **lu, int *pivot, Real *b, int n, Real *x) {
  // forward substitution, Ly=Pb (L must be a UNIT lower-triangular matrix)
  for (int i=0; i<n; i++) {
    // initialize the solution to the RHS values, repeating permutation from LU decomp.
    x[i] = b[pivot[i]];
    for (int j=0; j<i; j++)
      x[i] -= lu[i][j]*x[j];
  }

  // back substitution, Ux=y (U is a NOT a unit upper-triangular matrix)
  for (int i=(n-1); i>=0; i--) {
    for (int j=(i+1); j<n; j++) {
      x[i] -= x[j]*lu[i][j];
    }
    x[i] /= lu[i][i];
  }
  return;
}
} // namespace
