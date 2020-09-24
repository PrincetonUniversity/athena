//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file awa_test.cpp
//  \brief Initial conditions for Apples with Apples Test

#include <cassert> // assert
#include <iostream>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "../z4c/z4c.hpp"
#include "../z4c/trackers.hpp"

// twopuncturesc: Stand-alone library ripped from Cactus
#include "TwoPunctures.h"

using namespace std;

int RefinementCondition(MeshBlock *pmb);
// QUESTION: is it better to setup two different problems instead of using ifdef?
static ini_data *data;

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) //, int res_flag)
{

    //if (!res_flag) {     
      string set_name = "problem";
      //printf("BeforeSetDefault\n");
      TwoPunctures_params_set_default();
      //printf("AfterSetDefault\n");
      TwoPunctures_params_set_Boolean((char *) "verbose",
                                  pin->GetOrAddBoolean(set_name, "verbose", 0));
      TwoPunctures_params_set_Real((char *) "par_b",
                                   pin->GetOrAddReal(set_name, "par_b", 1.));
      TwoPunctures_params_set_Real((char *) "par_m_plus",
                                   pin->GetOrAddReal(set_name, "par_m_plus", 1.));
      TwoPunctures_params_set_Real((char *) "par_m_minus",
                                   pin->GetOrAddReal(set_name, "par_m_minus", 1.));

      TwoPunctures_params_set_Real((char *) "target_M_plus",
                                   pin->GetOrAddReal(set_name, "target_M_plus", 1.));

      TwoPunctures_params_set_Real((char *) "target_M_minus",
                                   pin->GetOrAddReal(set_name, "target_M_minus", 1.));

      TwoPunctures_params_set_Real((char *) "par_P_plus1",
                                   pin->GetOrAddReal(set_name, "par_P_plus1", 0.));
      TwoPunctures_params_set_Real((char *) "par_P_plus2",
                                   pin->GetOrAddReal(set_name, "par_P_plus2", 0.5));
      TwoPunctures_params_set_Real((char *) "par_P_plus3",
                                   pin->GetOrAddReal(set_name, "par_P_plus3", 0.));


      TwoPunctures_params_set_Real((char *) "par_P_minus1",
                                   pin->GetOrAddReal(set_name, "par_P_minus1", 0.));
      TwoPunctures_params_set_Real((char *) "par_P_minus2",
                                   pin->GetOrAddReal(set_name, "par_P_minus2", 0.5));
      TwoPunctures_params_set_Real((char *) "par_P_minus3",
                                   pin->GetOrAddReal(set_name, "par_P_minus3", 0.));


      TwoPunctures_params_set_Real((char *) "par_S_plus1",
                                   pin->GetOrAddReal(set_name, "par_S_plus1", 0.));
      TwoPunctures_params_set_Real((char *) "par_S_plus2",
                                   pin->GetOrAddReal(set_name, "par_S_plus2", 0.));
      TwoPunctures_params_set_Real((char *) "par_S_plus3",
                                   pin->GetOrAddReal(set_name, "par_S_plus3", 0.));


      TwoPunctures_params_set_Real((char *) "par_S_minus1",
                                   pin->GetOrAddReal(set_name, "par_S_minus1", 0.));
      TwoPunctures_params_set_Real((char *) "par_S_minus2",
                                   pin->GetOrAddReal(set_name, "par_S_minus2", 0.));
      TwoPunctures_params_set_Real((char *) "par_S_minus3",
                                   pin->GetOrAddReal(set_name, "par_S_minus3", 0.));
      TwoPunctures_params_set_Real((char *) "center_offset1",
                                   pin->GetOrAddReal(set_name, "center_offset1", 0.));

      TwoPunctures_params_set_Real((char *) "center_offset2",
                                   pin->GetOrAddReal(set_name, "center_offset2", 0.));
      TwoPunctures_params_set_Real((char *) "center_offset3",
                                   pin->GetOrAddReal(set_name, "center_offset3", 0.));

      TwoPunctures_params_set_Boolean((char *) "give_bare_mass",
                                   pin->GetOrAddBoolean(set_name, "give_bare_mass", 1));

      TwoPunctures_params_set_Int((char *) "npoints_A",
                                   pin->GetOrAddInteger(set_name, "npoints_A", 30));
      TwoPunctures_params_set_Int((char *) "npoints_B",
                                   pin->GetOrAddInteger(set_name, "npoints_B", 30));
      TwoPunctures_params_set_Int((char *) "npoints_phi",
                                   pin->GetOrAddInteger(set_name, "npoints_phi", 16));


      TwoPunctures_params_set_Real((char *) "Newton_tol",
                                   pin->GetOrAddReal(set_name, "Newton_tol", 1.e-10));

      TwoPunctures_params_set_Int((char *) "Newton_maxit",
                                   pin->GetOrAddInteger(set_name, "Newton_maxit", 5));


      TwoPunctures_params_set_Real((char *) "TP_epsilon",
                                   pin->GetOrAddReal(set_name, "TP_epsilon", 0.));

      TwoPunctures_params_set_Real((char *) "TP_Tiny",
                                   pin->GetOrAddReal(set_name, "TP_Tiny", 0.));
      TwoPunctures_params_set_Real((char *) "TP_Extend_Radius",
                                   pin->GetOrAddReal(set_name, "TP_Extend_Radius", 0.));


      TwoPunctures_params_set_Real((char *) "adm_tol",
                                   pin->GetOrAddReal(set_name, "adm_tol", 1.e-10));


      TwoPunctures_params_set_Boolean((char *) "do_residuum_debug_output",
                                   pin->GetOrAddBoolean(set_name, "do_residuum_debug_output", 0));

      TwoPunctures_params_set_Boolean((char *) "solve_momentum_constraint",
                                   pin->GetOrAddBoolean(set_name, "solve_momentum_constraint", 0));

      TwoPunctures_params_set_Real((char *) "initial_lapse_psi_exponent",
                                   pin->GetOrAddReal(set_name, "initial_lapse_psi_exponent", -2.0));

      TwoPunctures_params_set_Boolean((char *) "swap_xz",
                                   pin->GetOrAddBoolean(set_name, "swap_xz", 0));
      //printf("AfterParSetting\n");
      data = TwoPunctures_make_initial_data();
      //printf("AfterMakeInitData\n");
    //}
    if(adaptive==true)
      EnrollUserRefinementCondition(RefinementCondition);

    return;
}

//void Mesh::UserWorkAfterLoop(ParameterInput *pin, int res_flag)
void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{
  //if (!res_flag)
  TwoPunctures_finalise(data);
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Sets the initial conditions.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{

#ifdef Z4C_ASSERT_FINITE
  // as a sanity check (these should be over-written)
  pz4c->adm.psi4.Fill(NAN);
  pz4c->adm.g_dd.Fill(NAN);
  pz4c->adm.K_dd.Fill(NAN);

  pz4c->z4c.chi.Fill(NAN);
  pz4c->z4c.Khat.Fill(NAN);
  pz4c->z4c.Theta.Fill(NAN);
  pz4c->z4c.alpha.Fill(NAN);
  pz4c->z4c.Gam_u.Fill(NAN);
  pz4c->z4c.beta_u.Fill(NAN);
  pz4c->z4c.g_dd.Fill(NAN);
  pz4c->z4c.A_dd.Fill(NAN);

  /*
  pz4c->con.C.Fill(NAN);
  pz4c->con.H.Fill(NAN);
  pz4c->con.M.Fill(NAN);
  pz4c->con.Z.Fill(NAN);
  pz4c->con.M_d.Fill(NAN);

  */
#endif //Z4C_ASSERT_FINITE
  //---------------------------------------------------------------------------

  // call the interpolation
  pz4c->ADMTwoPunctures(pin, pz4c->storage.adm, data);

  // collapse in both cases
  pz4c->GaugePreCollapsedLapse(pz4c->storage.adm, pz4c->storage.u);

  //std::cout << "Two punctures initialized." << std::endl;
  pz4c->ADMToZ4c(pz4c->storage.adm, pz4c->storage.u);

  // debug!
  /*
  pz4c->Z4cToADM(pz4c->storage.u, pz4c->storage.adm);

  pz4c->ADMConstraints(pz4c->storage.con, pz4c->storage.adm,
                       pz4c->storage.mat, pz4c->storage.u);

  */

#ifdef Z4C_ASSERT_FINITE
  pz4c->assert_is_finite_adm();
  pz4c->assert_is_finite_con();
  pz4c->assert_is_finite_mat();
  pz4c->assert_is_finite_z4c();
#endif //Z4C_ASSERT_FINITE

  return;
}

int RefinementCondition(MeshBlock *pmb)
{
#ifdef Z4C_TRACKER
  //Initial distance between one of the punctures and the edge of the full mesh, needed to
  //calculate the box-in-box grid structure
  Real L = pmb->pmy_mesh->pz4c_tracker->L_grid;
  int root_lev = pmb->pmy_mesh->pz4c_tracker->root_lev;
#ifdef DEBUG
  printf("Root lev = %d\n", root_lev);
#endif
  Real xv[24];
#ifdef DEBUG
  printf("Max x = %g\n", pmb->block_size.x1max);
  printf("Min x = %g\n", pmb->block_size.x1min);
#endif
  //Needed to calculate coordinates of vertices of a block with same center but
  //edge of 1/8th of the original size
  Real x1sum_sup = (5*pmb->block_size.x1max+3*pmb->block_size.x1min)/8.;
  Real x1sum_inf = (3*pmb->block_size.x1max+5*pmb->block_size.x1min)/8.;
  Real x2sum_sup = (5*pmb->block_size.x2max+3*pmb->block_size.x2min)/8.;
  Real x2sum_inf = (3*pmb->block_size.x2max+5*pmb->block_size.x2min)/8.;
  Real x3sum_sup = (5*pmb->block_size.x3max+3*pmb->block_size.x3min)/8.;
  Real x3sum_inf = (3*pmb->block_size.x3max+5*pmb->block_size.x3min)/8.;

  xv[0] = x1sum_sup;
  xv[1] = x2sum_sup;
  xv[2] = x3sum_sup;

  xv[3] = x1sum_sup;
  xv[4] = x2sum_sup;
  xv[5] = x3sum_inf;

  xv[6] = x1sum_sup;
  xv[7] = x2sum_inf;
  xv[8] = x3sum_sup;

  xv[9] = x1sum_sup;
  xv[10] = x2sum_inf;
  xv[11] = x3sum_inf;

  xv[12] = x1sum_inf;
  xv[13] = x2sum_sup;
  xv[14] = x3sum_sup;

  xv[15] = x1sum_inf;
  xv[16] = x2sum_sup;
  xv[17] = x3sum_inf;

  xv[18] = x1sum_inf;
  xv[19] = x2sum_inf;
  xv[20] = x3sum_sup;

  xv[21] = x1sum_inf;
  xv[22] = x2sum_inf;
  xv[23] = x3sum_inf;

  //Level of current block
  int level = pmb->loc.level-root_lev;
#ifdef DEBUG
  printf("\n<===================================================>\n");
  printf("L = %g\n",L);
  printf("lev = %d\n",level);
#endif
  // Min distance between the two punctures
  Real d = 1000000;
  for (int i_punct = 0; i_punct < NPUNCT; ++i_punct) {
    // Abs difference
    Real diff;
    // Max norm_inf
    Real dmin_punct = 1000000;
#ifdef DEBUG
    printf("==> Punc = %d\n", i_punct);
#endif
    for (int i_vert = 0; i_vert < 8; ++i_vert) {
      // Norm_inf
      Real norm_inf = -1;
      for (int i_diff = 0; i_diff < 3; ++ i_diff) {
        diff = std::abs(pmb->pmy_mesh->pz4c_tracker->pos_body[i_punct].pos[i_diff] - xv[i_vert*3+i_diff]);
#ifdef DEBUG
        printf("======> Coordpos = %g, coordblock = %g\n",pmb->pmy_mesh->pz4c_tracker->pos_body[i_punct].pos[i_diff], xv[i_vert*3+i_diff]);
#endif
        if (diff > norm_inf) {
          norm_inf = diff;
        }
#ifdef DEBUG
	printf("======> Dist = %g\n", diff);
#endif
      }
#ifdef DEBUG
      printf("====> Inf norm = %g\n", norm_inf);
#endif
      //Calculate minimum of the distances of the 8 vertices above
      if (dmin_punct > norm_inf) {
        dmin_punct = norm_inf;
      }
    }
#ifdef DEBUG
    printf("====> dmin_punct = %g\n", dmin_punct);
#endif
    //Calculate minimum of the distances between the n punctures
    if (d > dmin_punct) {
      d = dmin_punct;
    }
  }
#ifdef DEBUG
  printf("Min dist = %g\n", d);
#endif
  Real ratio = L/d;
  if (ratio < 1) return -1;
  //Calculate level that the block should be in, given a box-in-box theoretical structure of the grid
  Real th_level = std::floor(std::log2(ratio));
#ifdef DEBUG
  printf("Level = %d, th_level = %g\n", level, th_level);
  printf("<===================================================>\n");
#endif
  if (th_level > level) {
#ifdef DEBUG
    printf("Refine\n");
#endif
    return 1;
  } else if (th_level < level) {
#ifdef DEBUG
    printf("Derefine\n");
#endif
    return -1;
  } else
#ifdef DEBUG
    printf("Do nothing\n");
#endif
    return 0;

#else // Z4C_TRACKER
  return 0;
#endif // Z4C_TRACKER
}
