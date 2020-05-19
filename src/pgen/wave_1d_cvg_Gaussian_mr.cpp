//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file wave_test.cpp
//  \brief Initial conditions for the wave equation

#include <cassert> // assert
#include <cmath> // abs, exp, sin, fmod
#include <iostream>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "../mesh/mesh_refinement.hpp"
#include "../wave/wave.hpp"

// using namespace std;

namespace {

  // standard deviation of packet
  Real sigma = 1. / 16.;
  Real phys_x1min = -1.0;
  Real phys_x1max = 1.0;
  Real amr_sigma_mul = 1;

  bool allow_restrict = true;

  // StackOverflow: 4633177
  Real wrapMax (Real x, Real max) {
    // wrap x -> [0, max)
    // integer math: (max + x % max) % max
    return std::fmod(max + std::fmod(x, max), max);
  }

  Real wrapMinMax (Real x, Real min, Real max) {
    // wrap x -> [min, max)
    return min + wrapMax(x - min, max - min);
  }

} // namespace

int RefinementCondition(MeshBlock *pmb);

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  if(adaptive==true)
    EnrollUserRefinementCondition(RefinementCondition);
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Initialize the problem.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  pwave->use_Sommerfeld = false;
  sigma = pin->GetOrAddReal("wave", "sigma", sigma);
  // for amr
  amr_sigma_mul = pin->GetOrAddReal("wave", "amr_sigma_mul", amr_sigma_mul);
  phys_x1min = pin->GetOrAddReal("mesh", "x1min", phys_x1min);
  phys_x1max = pin->GetOrAddReal("mesh", "x1max", phys_x1max);
  allow_restrict = pin->GetOrAddBoolean("wave", "allow_restrict", allow_restrict);
  //-

  int il = pwave->mbi.il, iu = pwave->mbi.iu;
  int kl = pwave->mbi.kl, ku = pwave->mbi.ku;
  int jl = pwave->mbi.jl, ju = pwave->mbi.ju;

  Real c = pwave->c;

  for(int k = kl; k <= ku; ++k)
    for(int j = jl; j <= ju; ++j)
      for(int i = il; i <= iu; ++i) {

        Real x = pwave->mbi.x1(i);
        Real y = pwave->mbi.x2(j);
        Real z = pwave->mbi.x3(k);

        Real x_2 = SQR(x);
        Real sigma_2 = SQR(sigma);
        Real Gaussian = exp(-x_2 / (2. * sigma_2));

        //Real cos_x = cos(PI / 10. * x);

        pwave->u(0,k,j,i) = Gaussian;// + cos_x;
        pwave->u(1,k,j,i) = x * Gaussian / sigma_2;// - PI / 10. * sin(PI / 10. * x);

        pwave->exact(k,j,i) = pwave->u(0,k,j,i);
        pwave->error(k,j,i) = 0.0;
      }


  return;
}

void MeshBlock::WaveUserWorkInLoop() {
  Real max_err = 0;
  Real fun_max = 0;

  Real t = pmy_mesh->time + pmy_mesh->dt;
  bool const debug_inspect_error = pwave->debug_inspect_error;
  Real const debug_abort_threshold = pwave->debug_abort_threshold;

  Real c = pwave->c;
  Real sigma_2 = SQR(sigma);


  int il = 0, iu = pwave->mbi.nn1 - 1;
  int jl = 0, ju = pwave->mbi.nn2 - 1;
  int kl = 0, ku = pwave->mbi.nn3 - 1;

  for(int k = kl; k <= ku; ++k)
    for(int j = jl; j <= ju; ++j)
      for(int i = il; i <= iu; ++i) {

        Real x = pwave->mbi.x1(i);
        Real y = pwave->mbi.x2(j);
        Real z = pwave->mbi.x3(k);

        // wrap argument to domain
        Real arg = wrapMinMax(x - t * c, phys_x1min, phys_x1max);
        Real Gaussian = exp(-SQR(arg) / (2. * sigma_2));

        pwave->exact(k,j,i) = Gaussian;
        // pwave->exact(k,j,i) += cos(PI / 10. * arg);

        pwave->error(k,j,i) = pwave->u(0,k,j,i) - pwave->exact(k,j,i);

        if (std::abs(pwave->error(k,j,i)) > max_err){
          max_err = std::abs(pwave->error(k,j,i));
          fun_max = pwave->u(0,k,j,i);
        }
      }

  if (debug_inspect_error) {
    printf(">>>\n");
    coutBoldRed("MB::UWIL gid = ");
    printf("%d\n", gid);
    printf("(max_err, fun_max, t)=(%1.18f, %1.18f, %1.18f)\n",
          max_err, fun_max, t);

    if (max_err > debug_abort_threshold) {
      printf("pwave->u:\n");
      pwave->u.print_all("%1.5f");

      printf("pwave->exact:\n");
      pwave->exact.print_all("%1.5f");

      printf("pwave->error:\n");
      pwave->error.print_all("%1.5f");

      Q();
    }

    printf("<<<\n");
  }
  return;
}

// //----------------------------------------------------------------------------------------
// //! \fn
// //  \brief debug for vertex dev
// void MeshBlock::DebugWaveMeshBlock(AthenaArray<Real> &u_wave,
//                                    int il, int iu,
//                                    int jl, int ju,
//                                    int kl, int ku,
//                                    bool is_additive,
//                                    bool is_coarse) {
//   Mesh *pm = pmy_mesh;

//   // compute next time-step (ensure not past limit)
//   // not initially available from mesh
//   Real dt = new_block_dt_;
//   if (pm->time < pm->tlim && (pm->tlim - pm->time) < dt)
//     dt = pm->tlim - pm->time;

//   Real c = pwave->c;
//   Real t = pm->time + dt;

//   // u_wave.ZeroClear(); // kill all data
//   // u_wave.print_all();

//   // solution on fundamental grid of block
//   for(int k = kl; k <= ku; ++k)
//     for(int j = jl; j <= ju; ++j)
//       for(int i = il; i <= iu; ++i) {

//         Real x, y, z;

//         if (is_coarse) {
//           x = pwave->mbi.cx1(i);
//           y = pwave->mbi.cx2(j);
//           z = pwave->mbi.cx3(k);
//         } else {
//           x = pwave->mbi.x1(i);
//           y = pwave->mbi.x2(j);
//           z = pwave->mbi.x3(k);
//         }

//         Real cos_2x = cos(2.*PI*x);
//         Real cos_4x = cos(4.*PI*x);
//         Real cos_4ct = cos(4.*PI*c*t);
//         Real sin_2ct = sin(2.*PI*c*t);

//         Real dt_cos_4ct = -4.*PI*c*sin(4.*PI*c*t);
//         Real dt_sin_2ct = 2.*PI*c*cos(2.*PI*c*t);

//         pwave->exact(k,j,i) = (2. + 2. * cos_4ct * cos_4x -
//                                cos_2x * sin_2ct / ( c * PI )) / 4.;

//         Real u_0 = pwave->exact(k,j,i);
//         Real u_1 = (2. * dt_cos_4ct * cos_4x -
//                     cos_2x * dt_sin_2ct / ( c * PI )) / 4.;


//         if (DBG_VC_CONSISTENCY) {
//           // BD: Debug- vertex consistency
//           u_0 = 1;
//           u_1 = 1;
//           //-
//         }

//         // unpack buffer additive logic emu
//         if (is_additive) {
//           u_wave(0, k, j, i) += u_0;
//           u_wave(1, k, j, i) += u_1;
//         } else {
//           u_wave(0, k, j, i) = u_0;
//           u_wave(1, k, j, i) = u_1;
//         }

//       }
// }

//----------------------------------------------------------------------------------------
//! \fn
//  \brief refinement condition: simple time-dependent test

int RefinementCondition(MeshBlock *pmb){
  // physical parameters
  Real c = pmb->pwave->c;
  Real t = pmb->pmy_mesh->time;
  Real del = (phys_x1max - phys_x1min);

  Real r_width = sigma * amr_sigma_mul;

  // Gaussian refinement params
  Real gc_x0 = phys_x1min + del / 2;
  Real gl_x0 = phys_x1min + del / 2 - r_width;
  Real gr_x0 = phys_x1min + del / 2 + r_width;

  // // fractional distance of interval travelled [from left edge] (wrapped)
  // Real gc_xfr = std::fmod((gc_x0 - phys_x1min) + t * c, del) / del;

  // Real gl_xfr = std::fmod((gl_x0 - phys_x1min) + t * c, del) / del;
  // Real gr_xfr = std::fmod((gr_x0 - phys_x1min) + t * c, del) / del;

  // // propagated and wrapped physical coordinates
  // Real gc_xc = grFundToPhys(gc_xfr);
  // Real gl_xc = grFundToPhys(gl_xfr);
  // Real gr_xc = grFundToPhys(gr_xfr);

  // propagated and wrapped physical coordinates
  Real gl_xc = wrapMinMax(gl_x0 + t * c, phys_x1min, phys_x1max);
  Real gc_xc = wrapMinMax(gc_x0 + t * c, phys_x1min, phys_x1max);
  Real gr_xc = wrapMinMax(gr_x0 + t * c, phys_x1min, phys_x1max);

  // current block (physical) geometry
  Real x1min = pmb->block_size.x1min;
  Real x1max = pmb->block_size.x1max;

  // if left or right edge within current box then refine
  if ((x1min <= gl_xc) && (gl_xc <= x1max)) {
    return 1;
  }

  if ((x1min <= gc_xc) && (gc_xc <= x1max)) {
    return 1;
  }

  if ((x1min <= gr_xc) && (gr_xc <= x1max)) {
    return 1;
  }

  // otherwise derefine
  if (allow_restrict) {
    // printf("r @ [%1.2f, %1.2f] ~ (%1.2f, %1.2f, %1.2f)\n",
    //        x1min, x1max, gl_xc, gc_xc, gr_xc);
    return -1;
  }
  return 0;
}
