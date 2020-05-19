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
// #include "../athena_tensor.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
// #include "../mesh/mesh_refinement.hpp"
#include "../z4c/z4c.hpp"

using namespace std;

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Sets the initial conditions.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  string puncture = pin->GetOrAddString("problem", "puncture", "one_puncture");

  pz4c->ADMOnePuncture(pz4c->storage.adm);
  pz4c->GaugePreCollapsedLapse(pz4c->storage.u);

  std::cout << "One puncture initialized." << std::endl;

  // Constructing Z4c vars from ADM ones
  pz4c->ADMToZ4c(pz4c->storage.adm, pz4c->storage.u);
  std::cout << "z4c vars generated." << std::endl;
  return;
}

// void MeshBlock::Z4cUserWorkInLoop() {
//   coutBoldRed("T\n");
//   pz4c->storage.u.print_dim();
// //  pz4c->con.H

// }

/*
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(22);
  return;
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {
        user_out_var(0,k,j,i) = pz4c->rhs.alpha(k,j,i);
        user_out_var(1,k,j,i) = pz4c->rhs.beta_u(0,k,j,i);
        user_out_var(2,k,j,i) = pz4c->rhs.beta_u(1,k,j,i);
        user_out_var(3,k,j,i) = pz4c->rhs.beta_u(2,k,j,i);
        user_out_var(4,k,j,i) = pz4c->rhs.A_dd(0,0,k,j,i); 
        user_out_var(5,k,j,i) = pz4c->rhs.A_dd(0,1,k,j,i);
        user_out_var(6,k,j,i) = pz4c->rhs.A_dd(0,2,k,j,i);
        user_out_var(7,k,j,i) = pz4c->rhs.A_dd(1,1,k,j,i);
        user_out_var(8,k,j,i) = pz4c->rhs.A_dd(1,2,k,j,i);
        user_out_var(9,k,j,i) = pz4c->rhs.A_dd(2,2,k,j,i);
        user_out_var(10,k,j,i) = pz4c->rhs.chi(k,j,i);
        user_out_var(11,k,j,i) = pz4c->rhs.g_dd(0,0,k,j,i);
        user_out_var(12,k,j,i) = pz4c->rhs.Gam_u(0,k,j,i);
        user_out_var(13,k,j,i) = pz4c->rhs.g_dd(0,1,k,j,i);
        user_out_var(14,k,j,i) = pz4c->rhs.g_dd(0,2,k,j,i);
        user_out_var(15,k,j,i) = pz4c->rhs.Gam_u(1,k,j,i);
        user_out_var(16,k,j,i) = pz4c->rhs.g_dd(1,1,k,j,i);
        user_out_var(17,k,j,i) = pz4c->rhs.g_dd(1,2,k,j,i);
        user_out_var(18,k,j,i) = pz4c->rhs.Gam_u(2,k,j,i);
        user_out_var(19,k,j,i) = pz4c->rhs.g_dd(2,2,k,j,i);
        user_out_var(20,k,j,i) = pz4c->rhs.Khat(k,j,i);
	user_out_var(21,k,j,i) = pz4c->rhs.Theta(k,j,i);
      }
    }
  }
}
*/
