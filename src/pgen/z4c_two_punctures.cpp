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

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  
      
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
    
    if(adaptive==true)
      EnrollUserRefinementCondition(RefinementCondition);
    
    return;
}

void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{ 
  TwoPunctures_finalise(data);
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Sets the initial conditions.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  // we shall default to one_puncture

  // call the interpolation
  pz4c->ADMTwoPunctures(pin, pz4c->storage.adm, data);

  // collapse in both cases
  pz4c->GaugePreCollapsedLapse(pz4c->storage.adm, pz4c->storage.u);

  //std::cout << "Two punctures initialized." << std::endl;
  pz4c->ADMToZ4c(pz4c->storage.adm, pz4c->storage.u);


  return;
}

int RefinementCondition(MeshBlock *pmb)
{ 
  for (int i_punct = 0; i_punct < NPUNCT; ++i_punct) {
    if (pmb->pz4c_tracker_loc->InBlock(i_punct))
      return 1;
    else
      return -1;
  }
  return 0;
}
