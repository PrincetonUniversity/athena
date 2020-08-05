//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file z4c_task_list.cpp
//  \brief time integrator for the z4c system (based on time_integrator.cpp)

// C/C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ classes headers
#include "../athena.hpp"
#include "../bvals/bvals.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../z4c/z4c.hpp"
//WGC wext
#ifdef Z4C_WEXT
#include "../z4c/wave_extract.hpp"
#endif
//WGC end
#ifdef Z4C_TRACKER
#include "../z4c/trackers.hpp"
#endif // Z4C_TRACKER
#include "../parameter_input.hpp"
#include "task_list.hpp"

// BD TODO: Significant code duplication with time_integrator, leave decoupled

//----------------------------------------------------------------------------------------
//  Z4cIntegratorTaskList constructor

Z4cIntegratorTaskList::Z4cIntegratorTaskList(ParameterInput *pin, Mesh *pm){
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
  integrator = pin->GetOrAddString("time", "integrator", "vl2");

  if (integrator == "vl2") {
    // VL: second-order van Leer integrator (Stone & Gardiner, NewA 14, 139 2009)
    // Simple predictor-corrector scheme similar to MUSCL-Hancock
    // Expressed in 2S or 3S* algorithm form
    nstages = 2;
    cfl_limit = 1.0;
    // Modify VL2 stability limit in 2D, 3D
    if (pm->ndim == 2) cfl_limit = 0.5;
    if (pm->ndim == 3) cfl_limit = 0.5;

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
    cfl_limit = 1.0;  // c_eff = c/nstages = 1/2 (Gottlieb (2009), pg 271)
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
    cfl_limit = 1.0;  // c_eff = c/nstages = 1/3 (Gottlieb (2009), pg 271)
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
    // Stability properties are similar to classical (non-SSP) RK4 (but ~2x L2 principal
    // error norm). Refer to Colella (2011) for linear stability analysis of constant
    // coeff. advection of classical RK4 + 4th or 1st order (limiter engaged) fluxes
    cfl_limit = 1.3925; // Colella (2011) eq 101; 1st order flux is most severe constraint
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
    // of Shu-Osher form matrix, alpha.
    nstages = 5;
    // Because it is an SSP method, we can use the SSP coefficient c=1.508 to to trivially
    // relate the CFL constraint to the RK1 CFL=1 (for first-order fluxes). There is no
    // need to perform stability analysis from scratch (unlike e.g. the linear stability
    // analysis for classical/non-SSP RK4 in Colella (2011)) However, PLM and PPM w/o the
    // limiter engaged are unconditionally unstable under RK1 integration, so the SSP
    // guarantees do not hold for the Athena++ spatial discretizations.
    cfl_limit = 1.508;         //  (effective SSP coeff = 0.302) Gottlieb (2009) pg 272
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

  // Now assemble list of tasks for each stage of z4c integrator
  {using namespace Z4cIntegratorTaskNames;
    AddTask(CALC_Z4CRHS, NONE);                // CalculateZ4cRHS
    AddTask(INT_Z4C, CALC_Z4CRHS);             // IntegrateZ4c

    AddTask(SEND_Z4C, INT_Z4C);                // SendZ4c
    AddTask(RECV_Z4C, NONE);                   // ReceiveZ4c

    AddTask(SETB_Z4C, (RECV_Z4C|INT_Z4C));     // SetBoundariesZ4c
    if (pm->multilevel) { // SMR or AMR
      AddTask(PROLONG, (SEND_Z4C|SETB_Z4C));   // Prolongation
      AddTask(PHY_BVAL, PROLONG);              // PhysicalBoundary
    } else {
      AddTask(PHY_BVAL, SETB_Z4C);             // PhysicalBoundary
    }

    AddTask(ALG_CONSTR, PHY_BVAL);             // EnforceAlgConstr
    AddTask(Z4C_TO_ADM, ALG_CONSTR);           // Z4cToADM
    AddTask(ADM_CONSTR, Z4C_TO_ADM);           // ADM_Constraints
//WGC wext
    AddTask(Z4C_WEYL, Z4C_TO_ADM);           // Calc Psi4
    AddTask(WAVE_EXTR, Z4C_WEYL);           // Project Psi4 multipoles
//WGC end
    AddTask(USERWORK, ADM_CONSTR);             // UserWork

    AddTask(NEW_DT, USERWORK);                 // NewBlockTimeStep
    if (pm->adaptive) {
      AddTask(FLAG_AMR, USERWORK);             // CheckRefinement
      AddTask(CLEAR_ALLBND, FLAG_AMR);         // ClearAllBoundary
    } else {
      AddTask(CLEAR_ALLBND, NEW_DT);           // ClearAllBoundary
    }

#ifdef Z4C_ASSERT_FINITE
    AddTask(ASSERT_FIN, CLEAR_ALLBND);         // AssertFinite
#endif // Z4C_ASSERT_FINITE

  } // end of using namespace block

}

//---------------------------------------------------------------------------------------
//  Sets id and dependency for "ntask" member of task_list_ array, then iterates value of
//  ntask.

void Z4cIntegratorTaskList::AddTask(const TaskID& id, const TaskID& dep) {
    task_list_[ntasks].task_id = id;
    task_list_[ntasks].dependency = dep;

    using namespace Z4cIntegratorTaskNames; // NOLINT (build/namespace)

    if (id == CLEAR_ALLBND) {
      task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&Z4cIntegratorTaskList::ClearAllBoundary);
      task_list_[ntasks].lb_time = false;
    } else if (id == CALC_Z4CRHS) {
      task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&Z4cIntegratorTaskList::CalculateZ4cRHS);
      task_list_[ntasks].lb_time = true;
    } else if (id == INT_Z4C) {
      task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&Z4cIntegratorTaskList::IntegrateZ4c);
      task_list_[ntasks].lb_time = true;
    } else if (id == SEND_Z4C) {
      task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&Z4cIntegratorTaskList::SendZ4c);
      task_list_[ntasks].lb_time = true;
    } else if (id == RECV_Z4C) {
      task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&Z4cIntegratorTaskList::ReceiveZ4c);
      task_list_[ntasks].lb_time = false;
    } else if (id == SETB_Z4C) {
      task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&Z4cIntegratorTaskList::SetBoundariesZ4c);
      task_list_[ntasks].lb_time = true;
    } else if (id == PROLONG) {
      task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&Z4cIntegratorTaskList::Prolongation);
      task_list_[ntasks].lb_time = true;
    } else if (id == PHY_BVAL) {
      task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&Z4cIntegratorTaskList::PhysicalBoundary);
      task_list_[ntasks].lb_time = true;
    } else if (id == ALG_CONSTR) {
      task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&Z4cIntegratorTaskList::EnforceAlgConstr);
      task_list_[ntasks].lb_time = true;
    } else if (id == Z4C_TO_ADM) {
      task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&Z4cIntegratorTaskList::Z4cToADM);
      task_list_[ntasks].lb_time = true;
    } else if (id == USERWORK) {
      task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&Z4cIntegratorTaskList::UserWork);
      task_list_[ntasks].lb_time = true;
    } else if (id == ADM_CONSTR) {
      task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&Z4cIntegratorTaskList::ADM_Constraints);
      task_list_[ntasks].lb_time = true;
//WGC wext
    } else if (id == Z4C_WEYL) {
      task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&Z4cIntegratorTaskList::Z4c_Weyl);
      task_list_[ntasks].lb_time = true;
    } else if (id == WAVE_EXTR) {
      task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&Z4cIntegratorTaskList::WaveExtract);
      task_list_[ntasks].lb_time = true;
//WGC end
    } else if (id == NEW_DT) {
      task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&Z4cIntegratorTaskList::NewBlockTimeStep);
      task_list_[ntasks].lb_time = true;
    } else if (id == FLAG_AMR) {
      task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&Z4cIntegratorTaskList::CheckRefinement);
      task_list_[ntasks].lb_time = true;
    }
#ifdef Z4C_ASSERT_FINITE
    else if (id == ASSERT_FIN) {
      task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&Z4cIntegratorTaskList::AssertFinite);
      task_list_[ntasks].lb_time = false;
    }
#endif // Z4C_ASSERT_FINITE
    else {
      std::stringstream msg;
      msg << "### FATAL ERROR in AddTask" << std::endl
          << "Invalid Task is specified" << std::endl;
      ATHENA_ERROR(msg);
    }

    ntasks++;
    return;
}


void Z4cIntegratorTaskList::StartupTaskList(MeshBlock *pmb, int stage) {
  // BD: debug-

  // application of Sommerfeld boundary conditions
  pmb->pz4c->Z4cBoundaryRHS(pmb->pz4c->storage.u,
                            pmb->pz4c->storage.mat,
                            pmb->pz4c->storage.rhs);

  BoundaryValues *pbval = pmb->pbval;

  Real t_end_stage = pmb->pmy_mesh->time + pmb->stage_abscissae[stage][0];
  // Scaled coefficient for RHS time-advance within stage
  Real dt = (stage_wghts[(stage-1)].beta)*(pmb->pmy_mesh->dt);

  if (PREFER_VC) {
    pbval->ApplyPhysicalVertexCenteredBoundaries(t_end_stage, dt);
  } else {
    pbval->ApplyPhysicalBoundaries(t_end_stage, dt);
  }
  //--------

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

    // Auxiliar var u1 needs to be initialized to 0 at the beginning of each cycle
    // Change to emulate PassiveScalars logic
    pmb->pz4c->storage.u1.ZeroClear();
    if (integrator == "ssprk5_4")
      pmb->pz4c->storage.u2 = pmb->pz4c->storage.u2;
  }

  pmb->pbval->StartReceiving(BoundaryCommSubset::all);
  return;
}

//----------------------------------------------------------------------------------------
// Functions to end MPI communication

TaskStatus Z4cIntegratorTaskList::ClearAllBoundary(MeshBlock *pmb, int stage) {
  pmb->pbval->ClearBoundary(BoundaryCommSubset::all);
  return TaskStatus::success;
}
//----------------------------------------------------------------------------------------
// Functions to calculate the RHS

TaskStatus Z4cIntegratorTaskList::CalculateZ4cRHS(MeshBlock *pmb, int stage) {

#ifdef Z4C_TRACKER
  // Tracker: interpolate beta at puncture position before evolution
  if (stage==1) {
    for (int i_punc = 0; i_punc<NPUNCT; i_punc++) {
      pmb->pz4c_tracker_loc->StoreBetaPrev(pmb->pz4c_tracker_loc->betap, pmb->pz4c->storage.u, i_punc);
    }
  }
#endif // Z4C_TRACKER

  if (stage <= nstages) {
    pmb->pz4c->Z4cRHS(pmb->pz4c->storage.u,
                      pmb->pz4c->storage.mat,
                      pmb->pz4c->storage.rhs);

    // application of Sommerfeld boundary conditions
    pmb->pz4c->Z4cBoundaryRHS(pmb->pz4c->storage.u,
                              pmb->pz4c->storage.mat,
                              pmb->pz4c->storage.rhs);
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}

//----------------------------------------------------------------------------------------
// Functions to integrate variables
TaskStatus Z4cIntegratorTaskList::IntegrateZ4c(MeshBlock *pmb, int stage) {
  Z4c *pz4c = pmb->pz4c;

  if (stage <= nstages) {
    // This time-integrator-specific averaging operation logic is identical
    // to IntegrateField
    Real ave_wghts[3];
    ave_wghts[0] = 1.0;
    ave_wghts[1] = stage_wghts[stage-1].delta;
    ave_wghts[2] = 0.0;
    pz4c->WeightedAve(pz4c->storage.u1, pz4c->storage.u,
                      pz4c->storage.u2, ave_wghts);

    ave_wghts[0] = stage_wghts[stage-1].gamma_1;
    ave_wghts[1] = stage_wghts[stage-1].gamma_2;
    ave_wghts[2] = stage_wghts[stage-1].gamma_3;

    pz4c->WeightedAve(pz4c->storage.u, pz4c->storage.u1,
                      pz4c->storage.u2, ave_wghts);
    pz4c->AddZ4cRHS(pz4c->storage.rhs, stage_wghts[stage-1].beta,
                    pz4c->storage.u);

    return TaskStatus::next;
  }
  return TaskStatus::fail;
}
//----------------------------------------------------------------------------------------
// Functions to communicate conserved variables between MeshBlocks

TaskStatus Z4cIntegratorTaskList::SendZ4c(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    pmb->pz4c->ubvar.SendBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::success;
}

//----------------------------------------------------------------------------------------
// Functions to receive conserved variables between MeshBlocks

TaskStatus Z4cIntegratorTaskList::ReceiveZ4c(MeshBlock *pmb, int stage) {
  bool ret;
  if (stage <= nstages) {
    ret = pmb->pz4c->ubvar.ReceiveBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  if (ret) {
    return TaskStatus::success;
  } else {
    return TaskStatus::fail;
  }
}

TaskStatus Z4cIntegratorTaskList::SetBoundariesZ4c(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    pmb->pz4c->ubvar.SetBoundaries();
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}
//--------------------------------------------------------------------------------------
// Functions for everything else
TaskStatus Z4cIntegratorTaskList::Prolongation(MeshBlock *pmb, int stage) {
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

TaskStatus Z4cIntegratorTaskList::PhysicalBoundary(MeshBlock *pmb, int stage) {
  BoundaryValues *pbval = pmb->pbval;

  if (stage <= nstages) {
    // Time at the end of stage for (u, b) register pair
    Real t_end_stage = pmb->pmy_mesh->time + pmb->stage_abscissae[stage][0];
    // Scaled coefficient for RHS time-advance within stage
    Real dt = (stage_wghts[(stage-1)].beta)*(pmb->pmy_mesh->dt);

    if (PREFER_VC) {
      pbval->ApplyPhysicalVertexCenteredBoundaries(t_end_stage, dt);
    } else {
      pbval->ApplyPhysicalBoundaries(t_end_stage, dt);
    }

  } else {
    return TaskStatus::fail;
  }

  return TaskStatus::success;
}

TaskStatus Z4cIntegratorTaskList::UserWork(MeshBlock *pmb, int stage) {
  if (stage != nstages) return TaskStatus::success; // only do on last stage

  pmb->Z4cUserWorkInLoop();
  return TaskStatus::success;
}

TaskStatus Z4cIntegratorTaskList::EnforceAlgConstr(MeshBlock *pmb, int stage) {
  if (stage != nstages) return TaskStatus::success; // only do on last stage

  pmb->pz4c->AlgConstr(pmb->pz4c->storage.u);
  return TaskStatus::success;
}

TaskStatus Z4cIntegratorTaskList::Z4cToADM(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    pmb->pz4c->Z4cToADM(pmb->pz4c->storage.u, pmb->pz4c->storage.adm);
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}
//WGC wext
TaskStatus Z4cIntegratorTaskList::Z4c_Weyl(MeshBlock *pmb, int stage) {
 if (stage <= nstages) { 
#ifdef Z4C_WEXT
 pmb->pz4c->Z4cWeyl(pmb->pz4c->storage.adm, pmb->pz4c->storage.mat,
    pmb->pz4c->storage.weyl);
#endif
   return TaskStatus::success;
 }
 else {
 return TaskStatus::fail;
 }
}  


TaskStatus Z4cIntegratorTaskList::WaveExtract(MeshBlock *pmb, int stage) {
 if (stage != nstages) return TaskStatus::success;  
#ifdef Z4C_WEXT
 AthenaArray<Real> u_R;
 AthenaArray<Real> u_I;
 u_R.InitWithShallowSlice(pmb->pz4c->storage.weyl, Z4c::I_WEY_rpsi4, 1);
 u_I.InitWithShallowSlice(pmb->pz4c->storage.weyl, Z4c::I_WEY_ipsi4, 1);
 for(int n = 0; n<NRAD;++n){
   pmb->pwave_extr_loc[n]->Decompose_multipole(u_R,u_I);
 }
#endif
  return TaskStatus::success;
  }


//WGC end


TaskStatus Z4cIntegratorTaskList::ADM_Constraints(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    pmb->pz4c->ADMConstraints(pmb->pz4c->storage.con, pmb->pz4c->storage.adm,
                              pmb->pz4c->storage.mat, pmb->pz4c->storage.u);
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}

TaskStatus Z4cIntegratorTaskList::NewBlockTimeStep(MeshBlock *pmb, int stage) {
  if (stage != nstages) return TaskStatus::success; // only do on last stage

  pmb->pz4c->NewBlockTimeStep();
  return TaskStatus::success;
}

TaskStatus Z4cIntegratorTaskList::CheckRefinement(MeshBlock *pmb, int stage) {
  if (stage != nstages) return TaskStatus::success; // only do on last stage

  pmb->pmr->CheckRefinementCondition();
  return TaskStatus::success;
}

#ifdef Z4C_ASSERT_FINITE
TaskStatus Z4cIntegratorTaskList::AssertFinite(MeshBlock *pmb, int stage) {
  if (stage != nstages) return TaskStatus::success; // only do on last stage

  pmb->pz4c->assert_is_finite_adm();
  pmb->pz4c->assert_is_finite_con();
  pmb->pz4c->assert_is_finite_mat();
  pmb->pz4c->assert_is_finite_z4c();

  return TaskStatus::success;
}
#endif // Z4C_ASSERT_FINITE
