//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file awa_test.cpp
//  \brief Initial conditions for Apples with Apples Test

#include <cassert> // assert
#include <limits>
#include <iostream>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "../z4c/z4c.hpp"
#include "../z4c/puncture_tracker.hpp"

// twopuncturesc: Stand-alone library ripped from Cactus
#include "TwoPunctures.h"

using namespace std;

int RefinementCondition(MeshBlock *pmb);
// QUESTION: is it better to setup two different problems instead of using ifdef?
static ini_data *data = NULL;

static Real par_b;
static Real L_grid;

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin)
{
    par_b = pin->GetOrAddReal("problem", "par_b", 1.);
    if (!resume_flag) {
      string set_name = "problem";
      TwoPunctures_params_set_default();
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
      data = TwoPunctures_make_initial_data();
    }
    if(adaptive==true)
      EnrollUserRefinementCondition(RefinementCondition);

    return;
}

void Mesh::DeleteTemporaryUserMeshData() {
  if (NULL != data) {
    TwoPunctures_finalise(data);
    data = NULL;
  }
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Sets the initial conditions.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  assert(NULL != data);

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

  pz4c->con.C.Fill(NAN);
  pz4c->con.H.Fill(NAN);
  pz4c->con.M.Fill(NAN);
  pz4c->con.Z.Fill(NAN);
  pz4c->con.M_d.Fill(NAN);

  // call the interpolation
  pz4c->ADMTwoPunctures(pin, pz4c->storage.adm, data);

  // collapse in both cases
  pz4c->GaugePreCollapsedLapse(pz4c->storage.adm, pz4c->storage.u);

  //std::cout << "Two punctures initialized." << std::endl;
  pz4c->ADMToZ4c(pz4c->storage.adm, pz4c->storage.u);

  pz4c->Z4cToADM(pz4c->storage.u, pz4c->storage.adm);
  pz4c->ADMConstraints(pz4c->storage.con, pz4c->storage.adm,
                       pz4c->storage.mat, pz4c->storage.u);

  pz4c->assert_is_finite_adm();
  pz4c->assert_is_finite_con();
  pz4c->assert_is_finite_mat();
  pz4c->assert_is_finite_z4c();

  return;
}

int RefinementCondition(MeshBlock *pmb)
{
  if (pmb->pmy_mesh->pz4c_tracker.size() != 2) {
    return 0;
  }

  int root_lev = pmb->pmy_mesh->GetRootLevel();
  int level = pmb->loc.level - root_lev;

  // Box in box ---------------------------------------------------------------
#ifdef Z4C_REF_BOX_IN_BOX
  L_grid = (pmb->pmy_mesh->mesh_size.x1max - pmb->pmy_mesh->mesh_size.x1min)/2. - par_b; 
  //Initial distance between one of the punctures and the edge of the full mesh, needed to
  //calculate the box-in-box grid structure
  Real L = L_grid;
  Real xv[24];

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

  // Min distance between the two punctures
  Real d = std::numeric_limits<Real>::max();
  for (auto ptracker : pmb->pmy_mesh->pz4c_tracker) {
    // Abs difference
    Real diff;
    // Max norm_inf
    Real dmin_punct = std::numeric_limits<Real>::max();
    for (int i_vert = 0; i_vert < 8; ++i_vert) {
      // Norm_inf
      Real norm_inf = -1;
      for (int i_diff = 0; i_diff < 3; ++ i_diff) {
        diff = std::abs(ptracker->GetPos(i_diff) - xv[i_vert*3+i_diff]);
        if (diff > norm_inf) {
          norm_inf = diff;
        }
      }
      //Calculate minimum of the distances of the 8 vertices above
      if (dmin_punct > norm_inf) {
        dmin_punct = norm_inf;
      }
    }
    //Calculate minimum of the closest between the n punctures
    if (d > dmin_punct) {
      d = dmin_punct;
    }
  }
  Real ratio = L/d;

  if (ratio < 1) {
    return -1;
  }

  //Calculate level that the block should be in, given a box-in-box theoretical structure of the grid
  Real th_level = std::floor(std::log2(ratio));
  if (th_level > level) {
    return 1;
  } else if (th_level < level) {
    return -1;
  }
  return 0;

#endif // Z4C_REF_BOX_IN_BOX

#ifdef Z4C_REF_SPHERES
  for (int six=0; six<pmb->pz4c->opt.sphere_zone_number; ++six) {
    Real xyz_wz[3] = {0., 0., 0.};

    // use tracker if we can and if it is relevant
    int const pix = pmb->pz4c->opt.sphere_zone_puncture(six);
    if (pix != -1) {
      xyz_wz[0] = pmb->pmy_mesh->pz4c_tracker[pix].pos[0];
      xyz_wz[1] = pmb->pmy_mesh->pz4c_tracker[pix].pos[1];
      xyz_wz[2] = pmb->pmy_mesh->pz4c_tracker[pix].pos[2];
    } else {
      xyz_wz[0] = pmb->pz4c->opt.sphere_zone_center1(six);
      xyz_wz[1] = pmb->pz4c->opt.sphere_zone_center2(six);
      xyz_wz[2] = pmb->pz4c->opt.sphere_zone_center3(six);
    }

    int const lev_wz = pmb->pz4c->opt.sphere_zone_levels(six);
    Real const R_wz = pmb->pz4c->opt.sphere_zone_radii(six);

    if (lev_wz > 0) { // ensure currently iterated sphere actually has non-trivial level
      if (pmb->SphereIntersects(xyz_wz[0], xyz_wz[1], xyz_wz[2], R_wz)) {
        need_ref = need_ref or (level < lev_wz);
        satisfied_ref = satisfied_ref or (level == lev_wz);
      }
    }
  }

  if (need_ref) {
    return 1;
  } else if (satisfied_ref) {
    return 0;
  }
  // force de-refine if no condition satisfied
  return -1;
#endif // Z4C_REF_SPHERES

  return 0;

}
