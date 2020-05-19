//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file wave_1d_cvg_Dirichlet.cpp
//  \brief Problem generator for one-dimensional Dirichlet problem wave equation

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

using namespace std;

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

  int il = pwave->mbi.il, iu = pwave->mbi.iu;
  int kl = pwave->mbi.kl, ku = pwave->mbi.ku;
  int jl = pwave->mbi.jl, ju = pwave->mbi.ju;

  Real c = pwave->c;

  int const num_ix = pwave->A_.GetDim1();
  Real const Lx1 = pwave->Lx1_;

  for(int k = kl; k <= ku; ++k)
    for(int j = jl; j <= ju; ++j)
      for(int i = il; i <= iu; ++i) {

        Real x = pwave->mbi.x1(i);
        Real y = pwave->mbi.x2(j);
        Real z = pwave->mbi.x3(k);

        pwave->u(0,k,j,i) = 0.;
        pwave->u(1,k,j,i) = 0.;

        for(int ix_s=1; ix_s<=num_ix; ++ix_s) {
          Real const beta_ix_s = PI * ix_s / Lx1;
          Real const sin_beta_x = sin(beta_ix_s * x);
          Real const A = pwave->A_(ix_s-1);
          Real const B = pwave->B_(ix_s-1);

          pwave->u(0,k,j,i) += A * sin_beta_x;
          pwave->u(1,k,j,i) += B * c * beta_ix_s * sin_beta_x;
        }


        if (DBG_VC_CONSISTENCY) {
          // BD: Debug- vertex consistency
          pwave->u(0,k,j,i) = 1;
          pwave->u(1,k,j,i) = 1;
          //-
        }

        pwave->exact(k,j,i) = pwave->u(0,k,j,i);
        pwave->error(k,j,i) = 0.0;
      }


  return;
}

void MeshBlock::WaveUserWorkInLoop() {
  Real max_err = 0;
  Real fun_max = 0;

  int il = pwave->mbi.il, iu = pwave->mbi.iu;
  int kl = pwave->mbi.kl, ku = pwave->mbi.ku;
  int jl = pwave->mbi.jl, ju = pwave->mbi.ju;

  // // test full block
  // int il = 0, iu = pwave->mbi.nn1 - 1;
  // int jl = 0, ju = pwave->mbi.nn2 - 1;
  // int kl = 0, ku = pwave->mbi.nn3 - 1;

  Real const c = pwave->c;
  Real const t = pmy_mesh->time + pmy_mesh->dt;
  Real const debug_abort_threshold = pwave->debug_abort_threshold;

  int const num_ix = pwave->A_.GetDim1();
  Real const Lx1 = pwave->Lx1_;

  for(int k = kl; k <= ku; ++k)
    for(int j = jl; j <= ju; ++j)
      for(int i = il; i <= iu; ++i) {

        Real const x = pwave->mbi.x1(i);
        Real const y = pwave->mbi.x2(j);
        Real const z = pwave->mbi.x3(k);

        pwave->exact(k,j,i) = 0.;
        if(false)
        if (i<pwave->mbi.il+NGHOST+1 || i>pwave->mbi.iu-NGHOST-1) {
          pwave->u(0,k,j,i) = 0.;
          pwave->u(1,k,j,i) = 0.;
        }
        for(int ix_s=1; ix_s<=num_ix; ++ix_s) {
          Real const beta_ix_s = PI * ix_s / Lx1;
          Real const cos_beta_ct = cos(beta_ix_s * c * t);
          Real const sin_beta_ct = sin(beta_ix_s * c * t);
          Real const sin_beta_x = sin(beta_ix_s * x);

          Real const A = pwave->A_(ix_s-1);
          Real const B = pwave->B_(ix_s-1);
          Real sol = (A * cos_beta_ct +
                      B * sin_beta_ct) * sin_beta_x;
          pwave->exact(k,j,i) += sol;

          // BD: debug
          if(false)
          if (i<pwave->mbi.il+NGHOST+1 || i>pwave->mbi.iu-NGHOST-1) {
            Real dt_cos_beta_ct = -beta_ix_s * c * sin_beta_ct;
            Real dt_sin_beta_ct = beta_ix_s * c * cos_beta_ct;
            Real dt_sol = (A * dt_cos_beta_ct +
                           B * dt_sin_beta_ct) * sin_beta_x;
            pwave->u(0,k,j,i) += sol;
            pwave->u(1,k,j,i) += dt_sol;
          }
        }

        pwave->error(k,j,i) = pwave->u(0,k,j,i) - pwave->exact(k,j,i);

        if (std::abs(pwave->error(k,j,i)) > max_err){
          max_err = std::abs(pwave->error(k,j,i));
          fun_max = pwave->u(0,k,j,i);
        }
      }

  if (max_err > debug_abort_threshold) {
    printf(">>>\n");
    coutBoldRed("MB::UWIL gid = ");
    printf("%d\n", gid);
    printf("(max_err, fun_max, t)=(%1.18f, %1.18f, %1.18f)\n",
          max_err, fun_max, t);

    if (max_err > 0.1) {
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

//----------------------------------------------------------------------------------------
//! \fn
//  \brief refinement condition: simple time-dependent test

int RefinementCondition(MeshBlock *pmb){
  // don't do anything
  return 0;
}

//----------------------------------------------------------------------------------------
//! \fn
//  \brief debug for vertex dev
void MeshBlock::DebugWaveMeshBlock(AthenaArray<Real> &u_wave,
                                   int il, int iu,
                                   int jl, int ju,
                                   int kl, int ku,
                                   bool is_additive,
                                   bool is_coarse) {
  Mesh *pm = pmy_mesh;

  // compute next time-step (ensure not past limit)
  // not initially available from mesh
  Real dt = new_block_dt_;
  if (pm->time < pm->tlim && (pm->tlim - pm->time) < dt)
    dt = pm->tlim - pm->time;

  Real const c = pwave->c;
  Real const t = pm->time + dt;

  int const num_ix = pwave->A_.GetDim1();
  Real const Lx1 = pwave->Lx1_;

  // solution on fundamental grid of block
  for(int k = kl; k <= ku; ++k)
    for(int j = jl; j <= ju; ++j)
      for(int i = il; i <= iu; ++i) {

        Real x, y, z;

        if (is_coarse) {
          x = pwave->mbi.cx1(i);
          y = pwave->mbi.cx2(j);
          z = pwave->mbi.cx3(k);
        } else {
          x = pwave->mbi.x1(i);
          y = pwave->mbi.x2(j);
          z = pwave->mbi.x3(k);
        }

        Real u_0 = 0.;
        Real u_1 = 0.;

        for(int ix_s=1; ix_s<=num_ix; ++ix_s) {
          Real const beta_ix_s = PI * ix_s / Lx1;
          Real const cos_beta_ct = cos(beta_ix_s * c * t);
          Real const sin_beta_ct = sin(beta_ix_s * c * t);

          Real const sin_beta_x = sin(beta_ix_s * x);

          Real const dt_cos_beta_ct = -beta_ix_s * c * sin(beta_ix_s * c * t);
          Real const dt_sin_beta_ct = beta_ix_s * c * cos(beta_ix_s * c * t);

          Real const A = pwave->A_(ix_s-1);
          Real const B = pwave->B_(ix_s-1);

          u_0 += (A * cos_beta_ct +
                  B * sin_beta_ct) * sin_beta_x;

          u_1 += (A * dt_cos_beta_ct +
                  B * dt_sin_beta_ct) * sin_beta_x;

        }

        pwave->exact(k,j,i) = u_0;

        if (DBG_VC_CONSISTENCY) {
          // BD: Debug- vertex consistency
          u_0 = 1;
          u_1 = 1;
          //-
        }

        // unpack buffer additive logic emu
        if (is_additive) {
          u_wave(0, k, j, i) += u_0;
          u_wave(1, k, j, i) += u_1;
        } else {
          u_wave(0, k, j, i) = u_0;
          u_wave(1, k, j, i) = u_1;
        }

      }
}
