//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file time_integrator.cpp
//  \brief derived class for time integrator task list. Can create task lists for one
//  of many different time integrators (e.g. van Leer, RK2, RK3, etc.)

// C headers

// C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../bvals/bvals.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../field/field_diffusion/field_diffusion.hpp"
#include "../gravity/gravity.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/hydro_diffusion/hydro_diffusion.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../reconstruct/reconstruction.hpp"
#include "../scalars/scalars.hpp"
#include "task_list.hpp"

//----------------------------------------------------------------------------------------
//  TimeIntegratorTaskList constructor

TimeIntegratorTaskList::TimeIntegratorTaskList(ParameterInput *pin, Mesh *pm) {
  // First, define each time-integrator by setting weights for each step of the algorithm
  // and the CFL number stability limit when coupled to the single-stage spatial operator.
  // Currently, the explicit, multistage time-integrators must be expressed as 2S-type
  // algorithms as in Ketcheson (2010) Algorithm 3, which incudes 2N (Williamson) and 2R
  // (van der Houwen) popular 2-register low-storage RK methods. The 2S-type integrators
  // depend on a bidiagonally sparse Shu-Osher representation; at each stage l:
  //
  //    U^{l} = a_{l,l-2}*U^{l-2} + a_{l-1}*U^{l-1}
  //          + b_{l,l-2}*dt*Div(F_{l-2}) + b_{l,l-1}*dt*Div(F_{l-1}),
  //
  // where U^{l-1} and U^{l-2} are previous stages and a_{l,l-2}, a_{l,l-1}=(1-a_{l,l-2}),
  // and b_{l,l-2}, b_{l,l-1} are weights that are different for each stage and
  // integrator. Previous timestep U^{0} = U^n is given, and the integrator solves
  // for U^{l} for 1 <= l <= nstages.
  //
  // The 2x RHS evaluations of Div(F) and source terms per stage is avoided by adding
  // another weighted average / caching of these terms each stage. The API and framework
  // is extensible to three register 3S* methods, although none are currently implemented.

  // Notation: exclusively using "stage", equivalent in lit. to "substage" or "substep"
  // (infrequently "step"), to refer to the intermediate values of U^{l} between each
  // "timestep" = "cycle" in explicit, multistage methods. This is to disambiguate the
  // temporal integration from other iterative sequences; "Step" is often used for generic
  // sequences in code, e.g. main.cpp: "Step 1: MPI"
  //
  // main.cpp invokes the tasklist in a for () loop from stage=1 to stage=ptlist->nstages

  // TODO(felker): validate Field and Hydro diffusion with RK3, RK4, SSPRK(5,4)
  integrator = pin->GetOrAddString("time", "integrator", "vl2");

  if (integrator == "vl2") {
    // VL: second-order van Leer integrator (Stone & Gardiner, NewA 14, 139 2009)
    // Simple predictor-corrector scheme similar to MUSCL-Hancock
    // Expressed in 2S or 3S* algorithm form
    nstages = 2;
    cfl_limit = 1.0;
    // Modify VL2 stability limit in 2D, 3D
    if (pm->ndim == 2) cfl_limit = 0.5;
    if (pm->ndim == 3) cfl_limit = ONE_3RD;

    stage_wghts[0].delta = 1.0; // required for consistency
    stage_wghts[0].gamma_1 = 0.0;
    stage_wghts[0].gamma_2 = 1.0;
    stage_wghts[0].gamma_3 = 0.0;
    stage_wghts[0].beta = 0.5;

    stage_wghts[1].delta = 0.0;
    stage_wghts[1].gamma_1 = 0.0;
    stage_wghts[1].gamma_2 = 1.0;
    stage_wghts[1].gamma_3 = 0.0;
    stage_wghts[1].beta = 1.0;
  } else if (integrator == "rk1") {
    // RK1: first-order Runge-Kutta / the forward Euler (FE) method
    nstages = 1;
    cfl_limit = 1.0;
    stage_wghts[0].delta = 1.0;
    stage_wghts[0].gamma_1 = 0.0;
    stage_wghts[0].gamma_2 = 1.0;
    stage_wghts[0].gamma_3 = 0.0;
    stage_wghts[0].beta = 1.0;
  } else if (integrator == "rk2") {
    // Heun's method / SSPRK (2,2): Gottlieb (2009) equation 3.1
    // Optimal (in error bounds) explicit two-stage, second-order SSPRK
    nstages = 2;
    cfl_limit = 1.0;
    stage_wghts[0].delta = 1.0;
    stage_wghts[0].gamma_1 = 0.0;
    stage_wghts[0].gamma_2 = 1.0;
    stage_wghts[0].gamma_3 = 0.0;
    stage_wghts[0].beta = 1.0;

    stage_wghts[1].delta = 0.0;
    stage_wghts[1].gamma_1 = 0.5;
    stage_wghts[1].gamma_2 = 0.5;
    stage_wghts[1].gamma_3 = 0.0;
    stage_wghts[1].beta = 0.5;
  } else if (integrator == "rk3") {
    // SSPRK (3,3): Gottlieb (2009) equation 3.2
    // Optimal (in error bounds) explicit three-stage, third-order SSPRK
    nstages = 3;
    cfl_limit = 1.0;
    stage_wghts[0].delta = 1.0;
    stage_wghts[0].gamma_1 = 0.0;
    stage_wghts[0].gamma_2 = 1.0;
    stage_wghts[0].gamma_3 = 0.0;
    stage_wghts[0].beta = 1.0;

    stage_wghts[1].delta = 0.0;
    stage_wghts[1].gamma_1 = 0.25;
    stage_wghts[1].gamma_2 = 0.75;
    stage_wghts[1].gamma_3 = 0.0;
    stage_wghts[1].beta = 0.25;

    stage_wghts[2].delta = 0.0;
    stage_wghts[2].gamma_1 = TWO_3RD;
    stage_wghts[2].gamma_2 = ONE_3RD;
    stage_wghts[2].gamma_3 = 0.0;
    stage_wghts[2].beta = TWO_3RD;
    //} else if (integrator == "ssprk5_3") {
    //} else if (integrator == "ssprk10_4") {
  } else if (integrator == "rk4") {
    // RK4()4[2S] from Table 2 of Ketcheson (2010)
    // Non-SSP, explicit four-stage, fourth-order RK
    nstages = 4;
    // Stability properties are similar to classical RK4
    // Refer to Colella (2011) for constant advection with 4th order fluxes
    // linear stability analysis
    cfl_limit = 1.3925;
    stage_wghts[0].delta = 1.0;
    stage_wghts[0].gamma_1 = 0.0;
    stage_wghts[0].gamma_2 = 1.0;
    stage_wghts[0].gamma_3 = 0.0;
    stage_wghts[0].beta = 1.193743905974738;

    stage_wghts[1].delta = 0.217683334308543;
    stage_wghts[1].gamma_1 = 0.121098479554482;
    stage_wghts[1].gamma_2 = 0.721781678111411;
    stage_wghts[1].gamma_3 = 0.0;
    stage_wghts[1].beta = 0.099279895495783;

    stage_wghts[2].delta = 1.065841341361089;
    stage_wghts[2].gamma_1 = -3.843833699660025;
    stage_wghts[2].gamma_2 = 2.121209265338722;
    stage_wghts[2].gamma_3 = 0.0;
    stage_wghts[2].beta = 1.131678018054042;

    stage_wghts[3].delta = 0.0;
    stage_wghts[3].gamma_1 = 0.546370891121863;
    stage_wghts[3].gamma_2 = 0.198653035682705;
    stage_wghts[3].gamma_3 = 0.0;
    stage_wghts[3].beta = 0.310665766509336;
  } else if (integrator == "ssprk5_4") {
    // SSPRK (5,4): Gottlieb (2009) section 3.1; between eq 3.3 and 3.4
    // Optimal (in error bounds) explicit five-stage, fourth-order SSPRK
    // 3N method, but there is no 3S* formulation due to irregular sparsity
    // of Shu-Osher form matrix, alpha
    nstages = 5;
    cfl_limit = 1.3925;
    // u^(1)
    stage_wghts[0].delta = 1.0; // u1 = u^n
    stage_wghts[0].gamma_1 = 0.0;
    stage_wghts[0].gamma_2 = 1.0;
    stage_wghts[0].gamma_3 = 0.0;
    stage_wghts[0].beta = 0.391752226571890;

    // u^(2)
    stage_wghts[1].delta = 0.0; // u1 = u^n
    stage_wghts[1].gamma_1 = 0.555629506348765;
    stage_wghts[1].gamma_2 = 0.444370493651235;
    stage_wghts[1].gamma_3 = 0.0;
    stage_wghts[1].beta = 0.368410593050371;

    // u^(3)
    stage_wghts[2].delta = 0.517231671970585; // u1 <- (u^n + d*u^(2))
    stage_wghts[2].gamma_1 = 0.379898148511597;
    stage_wghts[2].gamma_2 = 0.0;
    stage_wghts[2].gamma_3 = 0.620101851488403; // u^(n) coeff =  u2
    stage_wghts[2].beta = 0.251891774271694;

    // u^(4)
    stage_wghts[3].delta = 0.096059710526147; // u1 <- (u^n + d*u^(2) + d'*u^(3))
    stage_wghts[3].gamma_1 = 0.821920045606868;
    stage_wghts[3].gamma_2 = 0.0;
    stage_wghts[3].gamma_3 = 0.178079954393132; // u^(n) coeff =  u2
    stage_wghts[3].beta = 0.544974750228521;

    // u^(n+1) partial expression
    stage_wghts[4].delta = 0.0;
    stage_wghts[4].gamma_1 = 0.386708617503268; // 1 ulp lower than Gottlieb u^(4) coeff
    stage_wghts[4].gamma_2 = 1.0; // u1 <- (u^n + d*u^(2) + d'*u^(3))
    stage_wghts[4].gamma_3 = 1.0; // partial sum from hardcoded extra stage=4
    stage_wghts[4].beta = 0.226007483236906; // F(u^(4)) coeff.
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in TimeIntegratorTaskList constructor" << std::endl
        << "integrator=" << integrator << " not valid time integrator" << std::endl;
    ATHENA_ERROR(msg);
  }

  // Set cfl_number based on user input and time integrator CFL limit
  Real cfl_number = pin->GetReal("time", "cfl_number");
  if (cfl_number > cfl_limit
      && pm->fluid_setup == FluidFormulation::evolve) {
    std::cout << "### Warning in TimeIntegratorTaskList constructor" << std::endl
              << "User CFL number " << cfl_number << " must be smaller than " << cfl_limit
              << " for integrator=" << integrator << " in " << pm->ndim
              << "D simulation" << std::endl << "Setting to limit" << std::endl;
    cfl_number = cfl_limit;
  }
  // Save to Mesh class
  pm->cfl_number = cfl_number;

  // Now assemble list of tasks for each stage of time integrator
  {using namespace HydroIntegratorTaskNames; // NOLINT (build/namespace)
    // calculate hydro/field diffusive fluxes
    if (!STS_ENABLED) {
      AddTask(DIFFUSE_HYD,NONE);
      if (MAGNETIC_FIELDS_ENABLED) {
        AddTask(DIFFUSE_FLD,NONE);
        // compute hydro fluxes, integrate hydro variables
        AddTask(CALC_HYDFLX,(DIFFUSE_HYD|DIFFUSE_FLD));
      } else { // Hydro
        AddTask(CALC_HYDFLX,DIFFUSE_HYD);
      }
      if (NSCALARS > 0) {
        AddTask(DIFFUSE_SCLR,NONE);
        AddTask(CALC_SCLRFLX,(CALC_HYDFLX|DIFFUSE_SCLR));
      }
    } else { // STS enabled:
      AddTask(CALC_HYDFLX,NONE);
      if (NSCALARS > 0)
        AddTask(CALC_SCLRFLX,CALC_HYDFLX);
    }
    if (pm->multilevel) { // SMR or AMR
      AddTask(SEND_HYDFLX,CALC_HYDFLX);
      AddTask(RECV_HYDFLX,CALC_HYDFLX);
      AddTask(INT_HYD,RECV_HYDFLX);
    } else {
      AddTask(INT_HYD, CALC_HYDFLX);
    }
    AddTask(SRCTERM_HYD,INT_HYD);
    AddTask(SEND_HYD,SRCTERM_HYD);
    AddTask(RECV_HYD,NONE);
    AddTask(SETB_HYD,(RECV_HYD|SRCTERM_HYD));
    if (SHEARING_BOX) { // Shearingbox BC for Hydro
      AddTask(SEND_HYDSH,SETB_HYD);
      AddTask(RECV_HYDSH,SETB_HYD);
    }

    if (NSCALARS > 0) {
      if (pm->multilevel) {
        AddTask(SEND_SCLRFLX,CALC_SCLRFLX);
        AddTask(RECV_SCLRFLX,CALC_SCLRFLX);
        AddTask(INT_SCLR,RECV_SCLRFLX);
      } else {
        AddTask(INT_SCLR,CALC_SCLRFLX);
      }
      // there is no SRCTERM_SCLR task
      AddTask(SEND_SCLR,INT_SCLR);
      AddTask(RECV_SCLR,NONE);
      AddTask(SETB_SCLR,(RECV_SCLR|INT_SCLR));
      // if (SHEARING_BOX) {
      //   AddTask(SEND_SCLRSH,SETB_SCLR);
      //   AddTask(RECV_SCLRSH,SETB_SCLR);
      // }
    }

    if (MAGNETIC_FIELDS_ENABLED) { // MHD
      // compute MHD fluxes, integrate field
      AddTask(CALC_FLDFLX,CALC_HYDFLX);
      AddTask(SEND_FLDFLX,CALC_FLDFLX);
      AddTask(RECV_FLDFLX,SEND_FLDFLX);
      if (SHEARING_BOX) {// Shearingbox BC for EMF
        AddTask(SEND_EMFSH,RECV_FLDFLX);
        AddTask(RECV_EMFSH,RECV_FLDFLX);
        AddTask(RMAP_EMFSH,RECV_EMFSH);
        AddTask(INT_FLD,RMAP_EMFSH);
      } else {
        AddTask(INT_FLD,RECV_FLDFLX);
      }

      AddTask(SEND_FLD,INT_FLD);
      AddTask(RECV_FLD,NONE);
      AddTask(SETB_FLD,(RECV_FLD|INT_FLD));
      if (SHEARING_BOX) { // Shearingbox BC for Bfield
        AddTask(SEND_FLDSH,SETB_FLD);
        AddTask(RECV_FLDSH,SETB_FLD);
      }

      // TODO(felker): these nested conditionals are horrible now. Add option to AddTask
      // for "wait for all previously added tasks"?

      // prolongate, compute new primitives
      if (pm->multilevel) { // SMR or AMR
        if (NSCALARS > 0) {
          AddTask(PROLONG,(SEND_HYD|SETB_HYD|SEND_FLD|SETB_FLD|SEND_SCLR|SETB_SCLR));
        } else {
          AddTask(PROLONG,(SEND_HYD|SETB_HYD|SEND_FLD|SETB_FLD));
        }
        AddTask(CONS2PRIM,PROLONG);
      } else {
        if (SHEARING_BOX) {
          if (NSCALARS > 0) {
            AddTask(CONS2PRIM,
                    (SETB_HYD|SETB_FLD|SETB_SCLR|RECV_HYDSH|RECV_FLDSH|RMAP_EMFSH));
          } else {
            AddTask(CONS2PRIM,(SETB_HYD|SETB_FLD|RECV_HYDSH|RECV_FLDSH|RMAP_EMFSH));
          }
        } else {
          if (NSCALARS > 0) {
            AddTask(CONS2PRIM,(SETB_HYD|SETB_FLD|SETB_SCLR));
          } else {
            AddTask(CONS2PRIM,(SETB_HYD|SETB_FLD));
          }
        }
      }
    } else {  // HYDRO
      // prolongate, compute new primitives
      if (pm->multilevel) { // SMR or AMR
        if (NSCALARS > 0) {
          AddTask(PROLONG,(SEND_HYD|SETB_HYD|SETB_SCLR|SEND_SCLR));
        } else {
          AddTask(PROLONG,(SEND_HYD|SETB_HYD));
        }
        AddTask(CONS2PRIM,PROLONG);
      } else {
        if (SHEARING_BOX) {
          if (NSCALARS > 0) {
            AddTask(CONS2PRIM,(SETB_HYD|RECV_HYDSH|SETB_SCLR));  // RECV_SCLRSH
          } else {
            AddTask(CONS2PRIM,(SETB_HYD|RECV_HYDSH));
          }
        } else {
          if (NSCALARS > 0) {
            AddTask(CONS2PRIM,(SETB_HYD|SETB_SCLR));
          } else {
            AddTask(CONS2PRIM,(SETB_HYD));
          }
        }
      }
    }

    // everything else
    AddTask(PHY_BVAL,CONS2PRIM);
    AddTask(USERWORK,PHY_BVAL);
    AddTask(NEW_DT,USERWORK);
    if (pm->adaptive) {
      AddTask(FLAG_AMR,USERWORK);
      AddTask(CLEAR_ALLBND,FLAG_AMR);
    } else {
      AddTask(CLEAR_ALLBND,NEW_DT);
    }
  } // end of using namespace block
}

//---------------------------------------------------------------------------------------
//  Sets id and dependency for "ntask" member of task_list_ array, then iterates value of
//  ntask.

void TimeIntegratorTaskList::AddTask(const TaskID& id, const TaskID& dep) {
  task_list_[ntasks].task_id = id;
  task_list_[ntasks].dependency = dep;
  // TODO(felker): change naming convention of either/both of TASK_NAME and TaskFunc
  // There are some issues with the current names:
  // 1) VERB_OBJECT is confusing with ObjectVerb(). E.g. seeing SEND_HYD in the task list
  // assembly would lead the user to believe the corresponding function is SendHydro(),
  // when it is actually HydroSend()--- Probaby change function names to active voice
  // VerbObject() since "HydroFluxCalculate()" doesn't sound quite right.

  // Note, there are exceptions to the "verb+object" convention in some TASK_NAMES and
  // TaskFunc, e.g. NEW_DT + NewBlockTimeStep() and AMR_FLAG + CheckRefinement(),
  // SRCTERM_HYD and HydroSourceTerms(), USERWORK, PHY_BVAL, PROLONG, CONS2PRIM,
  // ... Although, AMR_FLAG = "flag blocks for AMR" should be FLAG_AMR in VERB_OBJECT
  using namespace HydroIntegratorTaskNames; // NOLINT (build/namespace)
  if (id == CLEAR_ALLBND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ClearAllBoundary);
    task_list_[ntasks].lb_time = false;
  } else if (id == CALC_HYDFLX) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CalculateHydroFlux);
    task_list_[ntasks].lb_time = true;
  } else if (id == CALC_FLDFLX) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CalculateEMF);
    task_list_[ntasks].lb_time = true;
  } else if (id == SEND_HYDFLX) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendHydroFlux);
    task_list_[ntasks].lb_time = true;
  } else if (id == SEND_FLDFLX) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendEMF);
    task_list_[ntasks].lb_time = true;
  } else if (id == RECV_HYDFLX) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveAndCorrectHydroFlux);
    task_list_[ntasks].lb_time = false;
  } else if (id == RECV_FLDFLX) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveAndCorrectEMF);
    task_list_[ntasks].lb_time = false;
  } else if (id == INT_HYD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::IntegrateHydro);
    task_list_[ntasks].lb_time = true;
  } else if (id == INT_FLD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::IntegrateField);
    task_list_[ntasks].lb_time = true;
  } else if (id == SRCTERM_HYD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::AddSourceTermsHydro);
    task_list_[ntasks].lb_time = true;
  } else if (id == SEND_HYD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendHydro);
    task_list_[ntasks].lb_time = true;
  } else if (id == SEND_FLD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendField);
    task_list_[ntasks].lb_time = true;
  } else if (id == RECV_HYD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveHydro);
    task_list_[ntasks].lb_time = false;
  } else if (id == RECV_FLD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveField);
    task_list_[ntasks].lb_time = false;
  } else if (id == SETB_HYD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SetBoundariesHydro);
    task_list_[ntasks].lb_time = true;
  } else if (id == SETB_FLD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SetBoundariesField);
    task_list_[ntasks].lb_time = true;
  } else if (id == SEND_HYDSH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendHydroShear);
    task_list_[ntasks].lb_time = true;
  } else if (id == RECV_HYDSH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveHydroShear);
    task_list_[ntasks].lb_time = false;
  } else if (id == SEND_FLDSH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendFieldShear);
    task_list_[ntasks].lb_time = true;
  } else if (id == RECV_FLDSH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveFieldShear);
    task_list_[ntasks].lb_time = false;
  } else if (id == SEND_EMFSH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendEMFShear);
    task_list_[ntasks].lb_time = true;
  } else if (id == RECV_EMFSH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveEMFShear);
    task_list_[ntasks].lb_time = false;
  } else if (id == RMAP_EMFSH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::RemapEMFShear);
    task_list_[ntasks].lb_time = true;
  } else if (id == PROLONG) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::Prolongation);
    task_list_[ntasks].lb_time = true;
  } else if (id == CONS2PRIM) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::Primitives);
    task_list_[ntasks].lb_time = true;
  } else if (id == PHY_BVAL) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::PhysicalBoundary);
    task_list_[ntasks].lb_time = true;
  } else if (id == USERWORK) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::UserWork);
    task_list_[ntasks].lb_time = true;
  } else if (id == NEW_DT) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::NewBlockTimeStep);
    task_list_[ntasks].lb_time = true;
  } else if (id == FLAG_AMR) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CheckRefinement);
    task_list_[ntasks].lb_time = true;
  } else if (id == DIFFUSE_HYD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::DiffuseHydro);
    task_list_[ntasks].lb_time = true;
  } else if (id == DIFFUSE_FLD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::DiffuseField);
    task_list_[ntasks].lb_time = true;
  } else if (id == CALC_SCLRFLX) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CalculateScalarFlux);
    task_list_[ntasks].lb_time = true;
  } else if (id == SEND_SCLRFLX) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendScalarFlux);
    task_list_[ntasks].lb_time = true;
  } else if (id == RECV_SCLRFLX) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveScalarFlux);
    task_list_[ntasks].lb_time = false;
  } else if (id == INT_SCLR) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::IntegrateScalars);
    task_list_[ntasks].lb_time = true;
  } else if (id == SEND_SCLR) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendScalars);
    task_list_[ntasks].lb_time = true;
  } else if (id == RECV_SCLR) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveScalars);
    task_list_[ntasks].lb_time = false;
  } else if (id == SETB_SCLR) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SetBoundariesScalars);
    task_list_[ntasks].lb_time = true;
  } else if (id == DIFFUSE_SCLR) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::DiffuseScalars);
    task_list_[ntasks].lb_time = true;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in AddTask" << std::endl
        << "Invalid Task is specified" << std::endl;
    ATHENA_ERROR(msg);
  }
  ntasks++;
  return;
}


void TimeIntegratorTaskList::StartupTaskList(MeshBlock *pmb, int stage) {
  if (stage == 1) {
    // For each Meshblock, initialize time abscissae of each memory register pair (u,b)
    // at stage=0 to correspond to the beginning of the interval [t^n, t^{n+1}]
    pmb->stage_abscissae[0][0] = 0.0;
    pmb->stage_abscissae[0][1] = 0.0; // u1 advances to u1 = 0*u1 + 1.0*u in stage=1
    pmb->stage_abscissae[0][2] = 0.0; // u2 = u cached for all stages in 3S* methods

    // Given overall timestep dt, compute the time abscissae for all registers, stages
    for (int l=1; l<=nstages; l++) {
      // Update the dt abscissae of each memory register to values at end of this stage
      const IntegratorWeight w = stage_wghts[l-1];

      // u1 = u1 + delta*u
      pmb->stage_abscissae[l][1] = pmb->stage_abscissae[l-1][1]
                                   + w.delta*pmb->stage_abscissae[l-1][0];
      // u = gamma_1*u + gamma_2*u1 + gamma_3*u2 + beta*dt*F(u)
      pmb->stage_abscissae[l][0] = w.gamma_1*pmb->stage_abscissae[l-1][0]
                                   + w.gamma_2*pmb->stage_abscissae[l][1]
                                   + w.gamma_3*pmb->stage_abscissae[l-1][2]
                                   + w.beta*pmb->pmy_mesh->dt;
      // u2 = u^n
      pmb->stage_abscissae[l][2] = 0.0;
    }

    // Initialize storage registers
    Hydro *ph = pmb->phydro;
    ph->u1.ZeroClear();
    if (integrator == "ssprk5_4")
      ph->u2 = ph->u;

    if (MAGNETIC_FIELDS_ENABLED) { // MHD
      Field *pf = pmb->pfield;
      pf->b1.x1f.ZeroClear();
      pf->b1.x2f.ZeroClear();
      pf->b1.x3f.ZeroClear();
    }
    if (NSCALARS > 0) {
      PassiveScalars *ps = pmb->pscalars;
      ps->s1.ZeroClear();
      if (integrator == "ssprk5_4")
        ps->s2 = ps->s;
    }
  }

  if (SHEARING_BOX) {
    Real dt = (stage_wghts[(stage-1)].beta)*(pmb->pmy_mesh->dt);
    Real time = pmb->pmy_mesh->time+dt;
    pmb->pbval->ComputeShear(time);
  }
  pmb->pbval->StartReceiving(BoundaryCommSubset::all);

  return;
}

//----------------------------------------------------------------------------------------
// Functions to end MPI communication

TaskStatus TimeIntegratorTaskList::ClearAllBoundary(MeshBlock *pmb, int stage) {
  pmb->pbval->ClearBoundary(BoundaryCommSubset::all);
  return TaskStatus::success;
}

//----------------------------------------------------------------------------------------
// Functions to calculates fluxes

TaskStatus TimeIntegratorTaskList::CalculateHydroFlux(MeshBlock *pmb, int stage) {
  Hydro *phydro = pmb->phydro;
  Field *pfield = pmb->pfield;

  if (stage <= nstages) {
    if ((stage == 1) && (integrator == "vl2")) {
      phydro->CalculateFluxes(phydro->w,  pfield->b,  pfield->bcc, 1);
      return TaskStatus::next;
    } else {
      phydro->CalculateFluxes(phydro->w,  pfield->b,  pfield->bcc, pmb->precon->xorder);
      return TaskStatus::next;
    }
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::CalculateEMF(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    pmb->pfield->ComputeCornerE(pmb->phydro->w,  pmb->pfield->bcc);
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}

//----------------------------------------------------------------------------------------
// Functions to communicate fluxes between MeshBlocks for flux correction with AMR

TaskStatus TimeIntegratorTaskList::SendHydroFlux(MeshBlock *pmb, int stage) {
  pmb->phydro->hbvar.SendFluxCorrection();
  return TaskStatus::success;
}


TaskStatus TimeIntegratorTaskList::SendEMF(MeshBlock *pmb, int stage) {
  pmb->pfield->fbvar.SendFluxCorrection();
  return TaskStatus::success;
}

//----------------------------------------------------------------------------------------
// Functions to receive fluxes between MeshBlocks

TaskStatus TimeIntegratorTaskList::ReceiveAndCorrectHydroFlux(MeshBlock *pmb, int stage) {
  if (pmb->phydro->hbvar.ReceiveFluxCorrection()) {
    return TaskStatus::next;
  } else {
    return TaskStatus::fail;
  }
}

TaskStatus TimeIntegratorTaskList::ReceiveAndCorrectEMF(MeshBlock *pmb, int stage) {
  if (pmb->pfield->fbvar.ReceiveFluxCorrection()) {
    return TaskStatus::next;
  } else {
    return TaskStatus::fail;
  }
}

//----------------------------------------------------------------------------------------
// Functions to integrate conserved variables

TaskStatus TimeIntegratorTaskList::IntegrateHydro(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;
  Field *pf = pmb->pfield;

  if (pmb->pmy_mesh->fluid_setup != FluidFormulation::evolve) return TaskStatus::next;

  if (stage <= nstages) {
    // This time-integrator-specific averaging operation logic is identical to FieldInt
    Real ave_wghts[3];
    ave_wghts[0] = 1.0;
    ave_wghts[1] = stage_wghts[stage-1].delta;
    ave_wghts[2] = 0.0;
    pmb->WeightedAve(ph->u1, ph->u, ph->u2, ave_wghts);

    ave_wghts[0] = stage_wghts[stage-1].gamma_1;
    ave_wghts[1] = stage_wghts[stage-1].gamma_2;
    ave_wghts[2] = stage_wghts[stage-1].gamma_3;
    if (ave_wghts[0] == 0.0 && ave_wghts[1] == 1.0 && ave_wghts[2] == 0.0)
      ph->u.SwapAthenaArray(ph->u1);
    else
      pmb->WeightedAve(ph->u, ph->u1, ph->u2, ave_wghts);

    const Real wght = stage_wghts[stage-1].beta*pmb->pmy_mesh->dt;
    ph->AddFluxDivergence(wght, ph->u);
    // add coordinate (geometric) source terms
    pmb->pcoord->AddCoordTermsDivergence(wght, ph->flux, ph->w, pf->bcc, ph->u);

    // Hardcode an additional flux divergence weighted average for the penultimate
    // stage of SSPRK(5,4) since it cannot be expressed in a 3S* framework
    if (stage == 4 && integrator == "ssprk5_4") {
      // From Gottlieb (2009), u^(n+1) partial calculation
      ave_wghts[0] = -1.0; // -u^(n) coeff.
      ave_wghts[1] = 0.0;
      ave_wghts[2] = 0.0;
      const Real beta = 0.063692468666290; // F(u^(3)) coeff.
      const Real wght = beta*pmb->pmy_mesh->dt;
      // writing out to u2 register
      pmb->WeightedAve(ph->u2, ph->u1, ph->u2, ave_wghts);

       ph->AddFluxDivergence(wght, ph->u2);
      // add coordinate (geometric) source terms
      pmb->pcoord->AddCoordTermsDivergence(wght, ph->flux, ph->w, pf->bcc, ph->u2);
    }
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::IntegrateField(MeshBlock *pmb, int stage) {
  Field *pf = pmb->pfield;

  if (pmb->pmy_mesh->fluid_setup != FluidFormulation::evolve) return TaskStatus::next;

  if (stage <= nstages) {
    // This time-integrator-specific averaging operation logic is identical to HydroInt
    Real ave_wghts[3];
    ave_wghts[0] = 1.0;
    ave_wghts[1] = stage_wghts[stage-1].delta;
    ave_wghts[2] = 0.0;
    pmb->WeightedAve(pf->b1, pf->b, pf->b2, ave_wghts);

    ave_wghts[0] = stage_wghts[stage-1].gamma_1;
    ave_wghts[1] = stage_wghts[stage-1].gamma_2;
    ave_wghts[2] = stage_wghts[stage-1].gamma_3;
    if (ave_wghts[0] == 0.0 && ave_wghts[1] == 1.0 && ave_wghts[2] == 0.0) {
      pf->b.x1f.SwapAthenaArray(pf->b1.x1f);
      pf->b.x2f.SwapAthenaArray(pf->b1.x2f);
      pf->b.x3f.SwapAthenaArray(pf->b1.x3f);
    } else {
      pmb->WeightedAve(pf->b, pf->b1, pf->b2, ave_wghts);
    }

    pf->CT(stage_wghts[stage-1].beta*pmb->pmy_mesh->dt, pf->b);

    return TaskStatus::next;
  }

  return TaskStatus::fail;
}

//----------------------------------------------------------------------------------------
// Functions to add source terms

TaskStatus TimeIntegratorTaskList::AddSourceTermsHydro(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;
  Field *pf = pmb->pfield;

  // return if there are no source terms to be added
  if (!(ph->hsrc.hydro_sourceterms_defined)
      || pmb->pmy_mesh->fluid_setup != FluidFormulation::evolve) return TaskStatus::next;

  if (stage <= nstages) {
    // Time at beginning of stage for u()
    Real t_start_stage = pmb->pmy_mesh->time + pmb->stage_abscissae[stage-1][0];
    // Scaled coefficient for RHS update
    Real dt = (stage_wghts[(stage-1)].beta)*(pmb->pmy_mesh->dt);
    // Evaluate the time-dependent source terms at the time at the beginning of the stage
    ph->hsrc.AddHydroSourceTerms(t_start_stage, dt, ph->flux, ph->w, pf->bcc, ph->u);
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::next;
}

//----------------------------------------------------------------------------------------
// Functions to calculate hydro diffusion fluxes (stored in HydroDiffusion::visflx[],
// cndflx[], added at the end of Hydro::CalculateFluxes()

TaskStatus TimeIntegratorTaskList::DiffuseHydro(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;

  // return if there are no diffusion to be added
  if (!(ph->hdif.hydro_diffusion_defined)
      || pmb->pmy_mesh->fluid_setup != FluidFormulation::evolve) return TaskStatus::next;

  if (stage <= nstages) {
    ph->hdif.CalcDiffusionFlux(ph->w, ph->u, ph->flux);
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::next;
}

//----------------------------------------------------------------------------------------
// Functions to calculate diffusion EMF

TaskStatus TimeIntegratorTaskList::DiffuseField(MeshBlock *pmb, int stage) {
  Field *pf = pmb->pfield;

  // return if there are no diffusion to be added
  if (!(pf->fdif.field_diffusion_defined)) return TaskStatus::next;

  if (stage <= nstages) {
    // TODO(pdmullen): DiffuseField is also called in SuperTimeStepTaskLsit. It must skip
    // Hall effect (once implemented) diffusion process in STS and always calculate those
    // terms in the main integrator.
    pf->fdif.CalcDiffusionEMF(pf->b, pf->bcc, pf->e);
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::next;
}

//----------------------------------------------------------------------------------------
// Functions to communicate conserved variables between MeshBlocks

TaskStatus TimeIntegratorTaskList::SendHydro(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    // Swap Hydro quantity in BoundaryVariable interface back to conserved var formulation
    // (also needed in SetBoundariesHydro(), since the tasks are independent)
    pmb->phydro->hbvar.SwapHydroQuantity(pmb->phydro->u, HydroBoundaryQuantity::cons);
    pmb->phydro->hbvar.SendBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::success;
}


TaskStatus TimeIntegratorTaskList::SendField(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    pmb->pfield->fbvar.SendBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::success;
}

//----------------------------------------------------------------------------------------
// Functions to receive conserved variables between MeshBlocks

TaskStatus TimeIntegratorTaskList::ReceiveHydro(MeshBlock *pmb, int stage) {
  bool ret;
  if (stage <= nstages) {
    ret = pmb->phydro->hbvar.ReceiveBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  if (ret) {
    return TaskStatus::success;
  } else {
    return TaskStatus::fail;
  }
}


TaskStatus TimeIntegratorTaskList::ReceiveField(MeshBlock *pmb, int stage) {
  bool ret;
  if (stage <= nstages) {
    ret = pmb->pfield->fbvar.ReceiveBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  if (ret) {
    return TaskStatus::success;
  } else {
    return TaskStatus::fail;
  }
}


TaskStatus TimeIntegratorTaskList::SetBoundariesHydro(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    pmb->phydro->hbvar.SwapHydroQuantity(pmb->phydro->u, HydroBoundaryQuantity::cons);
    pmb->phydro->hbvar.SetBoundaries();
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::SetBoundariesField(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    pmb->pfield->fbvar.SetBoundaries();
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::SendHydroShear(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    pmb->phydro->hbvar.SendShearingBoxBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::success;
}


TaskStatus TimeIntegratorTaskList::ReceiveHydroShear(MeshBlock *pmb, int stage) {
  bool ret;
  ret = false;
  if (stage <= nstages) {
    ret = pmb->phydro->hbvar.ReceiveShearingBoxBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  if (ret) {
    return TaskStatus::success;
  } else {
    return TaskStatus::fail;
  }
}


TaskStatus TimeIntegratorTaskList::SendFieldShear(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    pmb->pfield->fbvar.SendShearingBoxBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::success;
}


TaskStatus TimeIntegratorTaskList::ReceiveFieldShear(MeshBlock *pmb, int stage) {
  bool ret;
  ret = false;
  if (stage <= nstages) {
    ret = pmb->pfield->fbvar.ReceiveShearingBoxBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  if (ret) {
    return TaskStatus::success;
  } else {
    return TaskStatus::fail;
  }
}


TaskStatus TimeIntegratorTaskList::SendEMFShear(MeshBlock *pmb, int stage) {
  pmb->pfield->fbvar.SendEMFShearingBoxBoundaryCorrection();
  return TaskStatus::success;
}


TaskStatus TimeIntegratorTaskList::ReceiveEMFShear(MeshBlock *pmb, int stage) {
  if (pmb->pfield->fbvar.ReceiveEMFShearingBoxBoundaryCorrection()) {
    return TaskStatus::next;
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::RemapEMFShear(MeshBlock *pmb, int stage) {
  pmb->pfield->fbvar.RemapEMFShearingBoxBoundary();
  return TaskStatus::success;
}

//--------------------------------------------------------------------------------------
// Functions for everything else

TaskStatus TimeIntegratorTaskList::Prolongation(MeshBlock *pmb, int stage) {
  BoundaryValues *pbval = pmb->pbval;

  if (stage <= nstages) {
    // Time at the end of stage for (u, b) register pair
    Real t_end_stage = pmb->pmy_mesh->time + pmb->stage_abscissae[stage][0];
    // Scaled coefficient for RHS time-advance within stage
    Real dt = (stage_wghts[(stage-1)].beta)*(pmb->pmy_mesh->dt);
    pbval->ProlongateBoundaries(t_end_stage, dt);
  } else {
    return TaskStatus::fail;
  }

  return TaskStatus::success;
}


TaskStatus TimeIntegratorTaskList::Primitives(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;
  Field *pf = pmb->pfield;
  PassiveScalars *ps = pmb->pscalars;
  BoundaryValues *pbval = pmb->pbval;

  int il = pmb->is, iu = pmb->ie, jl = pmb->js, ju = pmb->je, kl = pmb->ks, ku = pmb->ke;
  if (pbval->nblevel[1][1][0] != -1) il -= NGHOST;
  if (pbval->nblevel[1][1][2] != -1) iu += NGHOST;
  if (pbval->nblevel[1][0][1] != -1) jl -= NGHOST;
  if (pbval->nblevel[1][2][1] != -1) ju += NGHOST;
  if (pbval->nblevel[0][1][1] != -1) kl -= NGHOST;
  if (pbval->nblevel[2][1][1] != -1) ku += NGHOST;

  if (stage <= nstages) {
    // At beginning of this task, ph->w contains previous stage's W(U) output
    // and ph->w1 is used as a register to store the current stage's output.
    // For the second order integrators VL2 and RK2, the prim_old initial guess for the
    // Newton-Raphson solver in GR EOS uses the following abscissae:
    // stage=1: W at t^n and
    // stage=2: W at t^{n+1/2} (VL2) or t^{n+1} (RK2)
    pmb->peos->ConservedToPrimitive(ph->u, ph->w, pf->b,
                                    ph->w1, pf->bcc, pmb->pcoord,
                                    il, iu, jl, ju, kl, ku);
    if (NSCALARS > 0) {
      // r1/r_old for GR is currently unused:
      pmb->peos->PassiveScalarConservedToPrimitive(ps->s, ph->w1, // ph->u, (updated rho)
                                                   ps->r, ps->r,
                                                   pmb->pcoord, il, iu, jl, ju, kl, ku);
    }
    // fourth-order EOS:
    if (pmb->precon->xorder == 4) {
      // for hydro, shrink buffer by 1 on all sides
      if (pbval->nblevel[1][1][0] != -1) il += 1;
      if (pbval->nblevel[1][1][2] != -1) iu -= 1;
      if (pbval->nblevel[1][0][1] != -1) jl += 1;
      if (pbval->nblevel[1][2][1] != -1) ju -= 1;
      if (pbval->nblevel[0][1][1] != -1) kl += 1;
      if (pbval->nblevel[2][1][1] != -1) ku -= 1;
      // for MHD, shrink buffer by 3
      // TODO(felker): add MHD loop limit calculation for 4th order W(U)
      pmb->peos->ConservedToPrimitiveCellAverage(ph->u, ph->w, pf->b,
                                                 ph->w1, pf->bcc, pmb->pcoord,
                                                 il, iu, jl, ju, kl, ku);
      if (NSCALARS > 0) {
        pmb->peos->PassiveScalarConservedToPrimitiveCellAverage(
            ps->s, ps->r, ps->r, pmb->pcoord, il, iu, jl, ju, kl, ku);
      }
    }
    // swap AthenaArray data pointers so that w now contains the updated w_out
    ph->w.SwapAthenaArray(ph->w1);
    // r1/r_old for GR is currently unused:
    // ps->r.SwapAthenaArray(ps->r1);
  } else {
    return TaskStatus::fail;
  }

  return TaskStatus::success;
}


TaskStatus TimeIntegratorTaskList::PhysicalBoundary(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;
  PassiveScalars *ps = pmb->pscalars;
  BoundaryValues *pbval = pmb->pbval;

  if (stage <= nstages) {
    // Time at the end of stage for (u, b) register pair
    Real t_end_stage = pmb->pmy_mesh->time + pmb->stage_abscissae[stage][0];
    // Scaled coefficient for RHS time-advance within stage
    Real dt = (stage_wghts[(stage-1)].beta)*(pmb->pmy_mesh->dt);
    // Swap Hydro and (possibly) passive scalar quantities in BoundaryVariable interface
    // from conserved to primitive formulations:
    ph->hbvar.SwapHydroQuantity(ph->w, HydroBoundaryQuantity::prim);
    if (NSCALARS > 0)
      ps->sbvar.var_cc = &(ps->r);
    pbval->ApplyPhysicalBoundaries(t_end_stage, dt);
  } else {
    return TaskStatus::fail;
  }

  return TaskStatus::success;
}


TaskStatus TimeIntegratorTaskList::UserWork(MeshBlock *pmb, int stage) {
  if (stage != nstages) return TaskStatus::success; // only do on last stage

  pmb->UserWorkInLoop();
  return TaskStatus::success;
}


TaskStatus TimeIntegratorTaskList::NewBlockTimeStep(MeshBlock *pmb, int stage) {
  if (stage != nstages) return TaskStatus::success; // only do on last stage

  pmb->phydro->NewBlockTimeStep();
  return TaskStatus::success;
}


TaskStatus TimeIntegratorTaskList::CheckRefinement(MeshBlock *pmb, int stage) {
  if (stage != nstages) return TaskStatus::success; // only do on last stage

  pmb->pmr->CheckRefinementCondition();
  return TaskStatus::success;
}


TaskStatus TimeIntegratorTaskList::CalculateScalarFlux(MeshBlock *pmb, int stage) {
  PassiveScalars *ps = pmb->pscalars;
  if (stage <= nstages) {
    if ((stage == 1) && (integrator == "vl2")) {
      ps->CalculateFluxes(ps->r, 1);
      return TaskStatus::next;
    } else {
      ps->CalculateFluxes(ps->r, pmb->precon->xorder);
      return TaskStatus::next;
    }
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::SendScalarFlux(MeshBlock *pmb, int stage) {
  pmb->pscalars->sbvar.SendFluxCorrection();
  return TaskStatus::success;
}


TaskStatus TimeIntegratorTaskList::ReceiveScalarFlux(MeshBlock *pmb, int stage) {
  if (pmb->pscalars->sbvar.ReceiveFluxCorrection()) {
    return TaskStatus::next;
  } else {
    return TaskStatus::fail;
  }
}


TaskStatus TimeIntegratorTaskList::IntegrateScalars(MeshBlock *pmb, int stage) {
  PassiveScalars *ps = pmb->pscalars;
  if (stage <= nstages) {
    // This time-integrator-specific averaging operation logic is identical to
    // IntegrateHydro, IntegrateField
    Real ave_wghts[3];
    ave_wghts[0] = 1.0;
    ave_wghts[1] = stage_wghts[stage-1].delta;
    ave_wghts[2] = 0.0;
    pmb->WeightedAve(ps->s1, ps->s, ps->s2, ave_wghts);

    ave_wghts[0] = stage_wghts[stage-1].gamma_1;
    ave_wghts[1] = stage_wghts[stage-1].gamma_2;
    ave_wghts[2] = stage_wghts[stage-1].gamma_3;
    if (ave_wghts[0] == 0.0 && ave_wghts[1] == 1.0 && ave_wghts[2] == 0.0)
      ps->s.SwapAthenaArray(ps->s1);
    else
      pmb->WeightedAve(ps->s, ps->s1, ps->s2, ave_wghts);

    const Real wght = stage_wghts[stage-1].beta*pmb->pmy_mesh->dt;
    ps->AddFluxDivergence(wght, ps->s);

    // Hardcode an additional flux divergence weighted average for the penultimate
    // stage of SSPRK(5,4) since it cannot be expressed in a 3S* framework
    if (stage == 4 && integrator == "ssprk5_4") {
      // From Gottlieb (2009), u^(n+1) partial calculation
      ave_wghts[0] = -1.0; // -u^(n) coeff.
      ave_wghts[1] = 0.0;
      ave_wghts[2] = 0.0;
      const Real beta = 0.063692468666290; // F(u^(3)) coeff.
      const Real wght = beta*pmb->pmy_mesh->dt;
      // writing out to s2 register
      pmb->WeightedAve(ps->s2, ps->s1, ps->s2, ave_wghts);
      ps->AddFluxDivergence(beta, ps->s2);
    }
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::SendScalars(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    // Swap PassiveScalars quantity in BoundaryVariable interface back to conserved var
    // formulation (also needed in SetBoundariesScalars() since the tasks are independent)
    pmb->pscalars->sbvar.var_cc = &(pmb->pscalars->s);
    pmb->pscalars->sbvar.SendBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::success;
}


TaskStatus TimeIntegratorTaskList::ReceiveScalars(MeshBlock *pmb, int stage) {
  bool ret;
  if (stage <= nstages) {
    ret = pmb->pscalars->sbvar.ReceiveBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  if (ret) {
    return TaskStatus::success;
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::success;
}


TaskStatus TimeIntegratorTaskList::SetBoundariesScalars(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    // Set PassiveScalars quantity in BoundaryVariable interface to cons var formulation
    pmb->pscalars->sbvar.var_cc = &(pmb->pscalars->s);
    pmb->pscalars->sbvar.SetBoundaries();
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::DiffuseScalars(MeshBlock *pmb, int stage) {
  PassiveScalars *ps = pmb->pscalars;
  Hydro *ph = pmb->phydro;
  // return if there are no diffusion to be added
  if (!(ps->scalar_diffusion_defined))
    return TaskStatus::next;

  if (stage <= nstages) {
    // TODO(felker): adapted directly from HydroDiffusion::ClearFlux. Deduplicate
    ps->diffusion_flx[X1DIR].ZeroClear();
    ps->diffusion_flx[X2DIR].ZeroClear();
    ps->diffusion_flx[X3DIR].ZeroClear();

    // unlike HydroDiffusion, only 1x passive scalar diffusive process is allowed, so
    // there is no need for counterpart to wrapper fn HydroDiffusion::CalcDiffusionFlux
    ps->DiffusiveFluxIso(ps->r, ph->w, ps->diffusion_flx);
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::next;
}
