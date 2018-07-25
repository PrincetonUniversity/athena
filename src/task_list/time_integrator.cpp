//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file time_integrator.cpp
//  \brief derived class for time integrator task list.  Can create task lists for one
//  of many different time integrators (e.g. van Leer, RK2, RK3, etc.)

// C/C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ classes headers
#include "../athena.hpp"
#include "../bvals/bvals.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../field/field_diffusion/field_diffusion.hpp"
#include "../gravity/gravity.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../hydro/hydro_diffusion/hydro_diffusion.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../reconstruct/reconstruction.hpp"
#include "task_list.hpp"
//----------------------------------------------------------------------------------------
//  TimeIntegratorTaskList constructor

TimeIntegratorTaskList::TimeIntegratorTaskList(ParameterInput *pin, Mesh *pm)
  : TaskList(pm) {
  // First, define each time-integrator by setting weights for each step of the algorithm
  // and the CFL number stability limit when coupled to the single-stage spatial operator.
  // Currently, the time-integrators must be expressed as 2S-type algorithms as in
  // Ketchenson (2010) Algorithm 3, which incudes 2N (Williamson) and 2R (van der Houwen)
  // popular 2-register low-storage RK methods. The 2S-type integrators depend on a
  // bidiagonally sparse Shu-Osher representation; at each stage l:
  //
  //    U^{l} = a_{l,l-2}*U^{l-2} + a_{l-1}*U^{l-1}
  //          + b_{l,l-2}*dt*Div(F_{l-2}) + b_{l,l-1}*dt*Div(F_{l-1}),
  //
  // where U^{l-1} and U^{l-2} are previous stages and a_{l,l-2}, a_{l,l-1}=(1-a_{l,l-2}),
  // and b_{l,l-2}, b_{l,l-1} are weights that are different for each stage and integrator
  //
  // The 2x RHS evaluations of Div(F) and source terms per stage is avoided by adding
  // another weighted average / caching of these terms each stage. The API and framework
  // is extensible to three register 3S* methods.

  // Notation: exclusively using "stage", equivalent in lit. to "substage" or "substep"
  // (sometimes "step"), to refer to intermediate values between "timesteps" = "cycles"
  // "Step" is often used for generic sequences in code, e.g. main.cpp: "Step 1: MPI"

  integrator = pin->GetOrAddString("time","integrator","vl2");
  int dim = 1;
  if (pm->mesh_size.nx2 > 1) dim = 2;
  if (pm->mesh_size.nx3 > 1) dim = 3;

  if (integrator == "vl2") {
    // VL: second-order van Leer integrator (Stone & Gardiner, NewA 14, 139 2009)
    // Simple predictor-corrector scheme similar to MUSCL-Hancock
    // Expressed in 2S or 3S* algorithm form
    nstages = 2;
    cfl_limit = 1.0;
    // Modify VL2 stability limit in 2D, 3D
    if (dim == 2) cfl_limit = 0.5;
    if (dim == 3) cfl_limit = ONE_3RD;

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
    // RK4()4[2S] from Table 2 of Ketchenson (2010)
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
    // SSPRK (5,4): Gottlieb (2009) section 3.1
    // Optimal (in error bounds) explicit five-stage, fourth-order SSPRK
    // 3N method, but there is no 3S* formulation due to irregular sparsity
    // of Shu-Osher form matrix, alpha
    nstages = 5;
    cfl_limit = 1.3925;
    stage_wghts[0].delta = 1.0;
    stage_wghts[0].gamma_1 = 0.0;
    stage_wghts[0].gamma_2 = 1.0;
    stage_wghts[0].gamma_3 = 0.0;
    stage_wghts[0].beta = 0.391752226571890;

    stage_wghts[1].delta = 0.0; // u1 = u^n
    stage_wghts[1].gamma_1 = 0.555629506348765;
    stage_wghts[1].gamma_2 = 0.444370493651235;
    stage_wghts[1].gamma_3 = 0.0;
    stage_wghts[1].beta = 0.368410593050371;

    stage_wghts[2].delta = 0.0;
    stage_wghts[2].gamma_1 = 0.379898148511597;
    stage_wghts[2].gamma_2 = 0.0;
    stage_wghts[2].gamma_3 = 0.620101851488403; // u2 = u^n
    stage_wghts[2].beta = 0.251891774271694;

    stage_wghts[3].delta = 0.0;
    stage_wghts[3].gamma_1 = TWO_3RD;
    stage_wghts[3].gamma_2 = ONE_3RD;
    stage_wghts[3].gamma_3 = 0.178079954393132; // u2 = u^n
    stage_wghts[3].beta = 0.544974750228521;

    stage_wghts[4].delta = 0.0;
    stage_wghts[4].gamma_1 = 0.386708617503268; // from Gottlieb (2009), u^(4) coeff.
    stage_wghts[4].gamma_2 = ONE_3RD;
    stage_wghts[4].gamma_3 = 0.0;
    stage_wghts[4].beta = 0.226007483236906; // from Gottlieb (2009), F(u^(4)) coeff.
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in CreateTimeIntegrator" << std::endl
        << "integrator=" << integrator << " not valid time integrator" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // Set cfl_number based on user input and time integrator CFL limit
  Real cfl_number = pin->GetReal("time","cfl_number");
  if (cfl_number > cfl_limit) {
    std::cout << "### Warning in CreateTimeIntegrator" << std::endl
        << "User CFL number " << cfl_number << " must be smaller than " << cfl_limit
        << " for integrator=" << integrator << " in "
        << dim << "D simulation" << std::endl << "Setting to limit" << std::endl;
    cfl_number = cfl_limit;
  }
  // Save to Mesh class
  pm->cfl_number = cfl_number;

  // Now assemble list of tasks for each stage of time integrator
  {using namespace HydroIntegratorTaskNames; // NOLINT (build/namespace)
    AddTimeIntegratorTask(STARTUP_INT,NONE);
    AddTimeIntegratorTask(START_ALLRECV,STARTUP_INT);
    // calculate hydro/field diffusive fluxes
    AddTimeIntegratorTask(DIFFUSE_HYD,START_ALLRECV);
    if (MAGNETIC_FIELDS_ENABLED)
      AddTimeIntegratorTask(DIFFUSE_FLD,START_ALLRECV);
    // compute hydro fluxes, integrate hydro variables
    if (MAGNETIC_FIELDS_ENABLED)
      AddTimeIntegratorTask(CALC_HYDFLX,(START_ALLRECV|DIFFUSE_HYD|DIFFUSE_FLD));
    else
      AddTimeIntegratorTask(CALC_HYDFLX,(START_ALLRECV|DIFFUSE_HYD));
    if (pm->multilevel==true) { // SMR or AMR
      AddTimeIntegratorTask(SEND_HYDFLX,CALC_HYDFLX);
      AddTimeIntegratorTask(RECV_HYDFLX,CALC_HYDFLX);
      AddTimeIntegratorTask(INT_HYD,RECV_HYDFLX);
    } else {
      AddTimeIntegratorTask(INT_HYD, CALC_HYDFLX);
    }
    AddTimeIntegratorTask(SRCTERM_HYD,INT_HYD);
    AddTimeIntegratorTask(SEND_HYD,SRCTERM_HYD);
    AddTimeIntegratorTask(RECV_HYD,START_ALLRECV);
    if (SHEARING_BOX) { // Shearingbox BC for Hydro
      AddTimeIntegratorTask(SEND_HYDSH,RECV_HYD);
      AddTimeIntegratorTask(RECV_HYDSH,RECV_HYD);
    }

    // compute MHD fluxes, integrate field
    if (MAGNETIC_FIELDS_ENABLED) { // MHD
      AddTimeIntegratorTask(CALC_FLDFLX,CALC_HYDFLX);
      AddTimeIntegratorTask(SEND_FLDFLX,CALC_FLDFLX);
      AddTimeIntegratorTask(RECV_FLDFLX,SEND_FLDFLX);
      if (SHEARING_BOX) {// Shearingbox BC for EMF
        AddTimeIntegratorTask(SEND_EMFSH,RECV_FLDFLX);
        AddTimeIntegratorTask(RECV_EMFSH,RECV_FLDFLX);
        AddTimeIntegratorTask(RMAP_EMFSH,RECV_EMFSH);
        AddTimeIntegratorTask(INT_FLD,RMAP_EMFSH);
      } else {
        AddTimeIntegratorTask(INT_FLD,RECV_FLDFLX);
      }

      AddTimeIntegratorTask(SEND_FLD,INT_FLD);
      AddTimeIntegratorTask(RECV_FLD,START_ALLRECV);
      if (SHEARING_BOX) { // Shearingbox BC for Bfield
        AddTimeIntegratorTask(SEND_FLDSH,RECV_FLD);
        AddTimeIntegratorTask(RECV_FLDSH,RECV_FLD);
      }
    }

    // prolongate, compute new primitives
    if (MAGNETIC_FIELDS_ENABLED) { // MHD
      if (pm->multilevel==true) { // SMR or AMR
        AddTimeIntegratorTask(PROLONG,(SEND_HYD|RECV_HYD|SEND_FLD|RECV_FLD));
        AddTimeIntegratorTask(CON2PRIM,PROLONG);
      } else {
        if (SHEARING_BOX) {
          AddTimeIntegratorTask(CON2PRIM,(INT_HYD|RECV_HYD|INT_FLD|RECV_FLD|
                                         RECV_HYDSH|RECV_FLDSH|RMAP_EMFSH));
        } else {
          AddTimeIntegratorTask(CON2PRIM,(INT_HYD|RECV_HYD|INT_FLD|RECV_FLD));
        }
      }
    } else {  // HYDRO
      if (pm->multilevel==true) { // SMR or AMR
        AddTimeIntegratorTask(PROLONG,(SEND_HYD|RECV_HYD));
        AddTimeIntegratorTask(CON2PRIM,PROLONG);
      } else {
        if (SHEARING_BOX) {
          AddTimeIntegratorTask(CON2PRIM,(INT_HYD|RECV_HYD|RECV_HYDSH));
        } else {
          AddTimeIntegratorTask(CON2PRIM,(INT_HYD|RECV_HYD));
        }
      }
    }

    // everything else
    AddTimeIntegratorTask(PHY_BVAL,CON2PRIM);
//    if (SELF_GRAVITY_ENABLED == 1) {
//      AddTimeIntegratorTask(CORR_GFLX,PHY_BVAL);
//      AddTimeIntegratorTask(USERWORK,CORR_GFLX);
//    } else {
    AddTimeIntegratorTask(USERWORK,PHY_BVAL);
//    }
    AddTimeIntegratorTask(NEW_DT,USERWORK);
    if (pm->adaptive==true) {
      AddTimeIntegratorTask(AMR_FLAG,USERWORK);
      AddTimeIntegratorTask(CLEAR_ALLBND,AMR_FLAG);
    } else {
      AddTimeIntegratorTask(CLEAR_ALLBND,NEW_DT);
    }

  } // end of using namespace block
}

//---------------------------------------------------------------------------------------
//  Sets id and dependency for "ntask" member of task_list_ array, then iterates value of
//  ntask.

void TimeIntegratorTaskList::AddTimeIntegratorTask(uint64_t id, uint64_t dep) {
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;

  using namespace HydroIntegratorTaskNames; // NOLINT (build/namespace)
  switch((id)) {
    case (START_ALLRECV):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::StartAllReceive);
      break;
    case (CLEAR_ALLBND):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ClearAllBoundary);
      break;

    case (CALC_HYDFLX):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CalculateFluxes);
      break;
    case (CALC_FLDFLX):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CalculateEMF);
      break;

    case (SEND_HYDFLX):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::FluxCorrectSend);
      break;
    case (SEND_FLDFLX):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::EMFCorrectSend);
      break;

    case (RECV_HYDFLX):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::FluxCorrectReceive);
      break;
    case (RECV_FLDFLX):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::EMFCorrectReceive);
      break;

    case (INT_HYD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::HydroIntegrate);
      break;
    case (INT_FLD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::FieldIntegrate);
      break;

    case (SRCTERM_HYD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::HydroSourceTerms);
      break;

    case (SEND_HYD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::HydroSend);
      break;
    case (SEND_FLD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::FieldSend);
      break;

    case (RECV_HYD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::HydroReceive);
      break;
    case (RECV_FLD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::FieldReceive);
      break;
    case (SEND_HYDSH):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::HydroShearSend);
      break;
    case (RECV_HYDSH):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::HydroShearReceive);
      break;
    case (SEND_FLDSH):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::FieldShearSend);
      break;
    case (RECV_FLDSH):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::FieldShearReceive);
      break;
    case (SEND_EMFSH):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::EMFShearSend);
      break;
    case (RECV_EMFSH):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::EMFShearReceive);
      break;
    case (RMAP_EMFSH):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::EMFShearRemap);
      break;

    case (PROLONG):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::Prolongation);
      break;
    case (CON2PRIM):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::Primitives);
      break;
    case (PHY_BVAL):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::PhysicalBoundary);
      break;
    case (USERWORK):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::UserWork);
      break;
    case (NEW_DT):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::NewBlockTimeStep);
      break;
    case (AMR_FLAG):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CheckRefinement);
      break;
    case (CORR_GFLX):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::GravFluxCorrection);
      break;
    case (STARTUP_INT):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::StartupIntegrator);
      break;
    case (DIFFUSE_HYD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::HydroDiffusion);
      break;
    case (DIFFUSE_FLD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::FieldDiffusion);
      break;

    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in AddTimeIntegratorTask" << std::endl
          << "Invalid Task "<< id << " is specified" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  ntasks++;
  return;
}

//----------------------------------------------------------------------------------------
// Functions to start/end MPI communication

enum TaskStatus TimeIntegratorTaskList::StartAllReceive(MeshBlock *pmb, int stage) {
  Real dt = (stage_wghts[(stage-1)].beta)*(pmb->pmy_mesh->dt);
  Real time = pmb->pmy_mesh->time+dt;
  pmb->pbval->StartReceivingAll(time);
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::ClearAllBoundary(MeshBlock *pmb, int stage) {
  pmb->pbval->ClearBoundaryAll();
  return TASK_SUCCESS;
}

//----------------------------------------------------------------------------------------
// Functions to calculates fluxes

enum TaskStatus TimeIntegratorTaskList::CalculateFluxes(MeshBlock *pmb, int stage) {
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;

  if (stage <= nstages) {
    if ((stage == 1) && (integrator == "vl2")) {
      phydro->CalculateFluxes(phydro->w,  pfield->b,  pfield->bcc, 1);
      return TASK_NEXT;
    } else {
      phydro->CalculateFluxes(phydro->w,  pfield->b,  pfield->bcc, pmb->precon->xorder);
      return TASK_NEXT;
    }
  }
  return TASK_FAIL;
}

enum TaskStatus TimeIntegratorTaskList::CalculateEMF(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    pmb->pfield->ComputeCornerE(pmb->phydro->w,  pmb->pfield->bcc);
    return TASK_NEXT;
  }
  return TASK_FAIL;
}

//----------------------------------------------------------------------------------------
// Functions to communicate fluxes between MeshBlocks for flux correction with AMR

enum TaskStatus TimeIntegratorTaskList::FluxCorrectSend(MeshBlock *pmb, int stage) {
  pmb->pbval->SendFluxCorrection(FLUX_HYDRO);
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::EMFCorrectSend(MeshBlock *pmb, int stage) {
  pmb->pbval->SendEMFCorrection();
  return TASK_SUCCESS;
}

//----------------------------------------------------------------------------------------
// Functions to receive fluxes between MeshBlocks

enum TaskStatus TimeIntegratorTaskList::FluxCorrectReceive(MeshBlock *pmb, int stage) {
  if (pmb->pbval->ReceiveFluxCorrection(FLUX_HYDRO) == true) {
    return TASK_NEXT;
  } else {
    return TASK_FAIL;
  }
}

enum TaskStatus TimeIntegratorTaskList::EMFCorrectReceive(MeshBlock *pmb, int stage) {
  if (pmb->pbval->ReceiveEMFCorrection() == true) {
    return TASK_NEXT;
  } else {
    return TASK_FAIL;
  }
}

//----------------------------------------------------------------------------------------
// Functions to integrate conserved variables

enum TaskStatus TimeIntegratorTaskList::HydroIntegrate(MeshBlock *pmb, int stage) {
  Hydro *ph=pmb->phydro;
  Field *pf=pmb->pfield;
  if (stage <= nstages) {
    // This time-integrator-specific averaging operation logic is identical to FieldInt
    Real ave_wghts[3];
    ave_wghts[0] = 1.0;
    ave_wghts[1] = stage_wghts[stage-1].delta;
    ave_wghts[2] = 0.0;
    ph->WeightedAveU(ph->u1,ph->u,ph->u2,ave_wghts);

    ave_wghts[0] = stage_wghts[stage-1].gamma_1;
    ave_wghts[1] = stage_wghts[stage-1].gamma_2;
    ave_wghts[2] = stage_wghts[stage-1].gamma_3;
    ph->WeightedAveU(ph->u,ph->u1,ph->u2,ave_wghts);
    ph->AddFluxDivergenceToAverage(ph->w,pf->bcc,stage_wghts[stage-1].beta,ph->u);

    return TASK_NEXT;
  }

  return TASK_FAIL;
}

enum TaskStatus TimeIntegratorTaskList::FieldIntegrate(MeshBlock *pmb, int stage) {
  Field *pf=pmb->pfield;

  if (stage <= nstages) {
    // This time-integrator-specific averaging operation logic is identical to HydroInt
    Real ave_wghts[3];
    ave_wghts[0] = 1.0;
    ave_wghts[1] = stage_wghts[stage-1].delta;
    ave_wghts[2] = 0.0;
    pf->WeightedAveB(pf->b1,pf->b,pf->b2,ave_wghts);

    ave_wghts[0] = stage_wghts[stage-1].gamma_1;
    ave_wghts[1] = stage_wghts[stage-1].gamma_2;
    ave_wghts[2] = stage_wghts[stage-1].gamma_3;
    pf->WeightedAveB(pf->b,pf->b1,pf->b2,ave_wghts);
    pf->CT(stage_wghts[stage-1].beta, pf->b);

    return TASK_NEXT;
  }

  return TASK_FAIL;
}

//----------------------------------------------------------------------------------------
// Functions to add source terms

enum TaskStatus TimeIntegratorTaskList::HydroSourceTerms(MeshBlock *pmb, int stage) {
  Hydro *ph=pmb->phydro;
  Field *pf=pmb->pfield;

  // return if there are no source terms to be added
  if (ph->psrc->hydro_sourceterms_defined == false) return TASK_NEXT;

  if (stage <= nstages) {
    // Time at beginning of stage for u()
    Real t_start_stage = pmb->pmy_mesh->time + pmb->stage_abscissae[stage-1][0];
    // Scaled coefficient for RHS update
    Real dt = (stage_wghts[(stage-1)].beta)*(pmb->pmy_mesh->dt);
    ph->psrc->AddHydroSourceTerms(t_start_stage, dt, ph->flux, ph->w, pf->bcc, ph->u);
  } else {
    // Evaluate the source terms at the beginning of the
    return TASK_FAIL;
  }
  return TASK_NEXT;
}

//----------------------------------------------------------------------------------------
// Functions to calculate hydro diffusion fluxes

enum TaskStatus TimeIntegratorTaskList::HydroDiffusion(MeshBlock *pmb, int stage) {
  Hydro *ph=pmb->phydro;

// return if there are no diffusion to be added
  if (ph->phdif->hydro_diffusion_defined == false) return TASK_NEXT;

  // *** this must be changed for the RK3 integrator
  if(stage <= nstages) {
    ph->phdif->CalcHydroDiffusionFlux(ph->w, ph->u, ph->flux);
  } else {
    return TASK_FAIL;
  }
  return TASK_NEXT;
}

//----------------------------------------------------------------------------------------
// Functions to calculate diffusion EMF

enum TaskStatus TimeIntegratorTaskList::FieldDiffusion(MeshBlock *pmb, int stage) {
  Field *pf=pmb->pfield;

// return if there are no diffusion to be added
  if (pf->pfdif->field_diffusion_defined == false) return TASK_NEXT;

  // *** this must be changed for the RK3 integrator
  if(stage <= nstages) {
    pf->pfdif->CalcFieldDiffusionEMF(pf->b,pf->bcc,pf->e);
  } else {
    return TASK_FAIL;
  }
  return TASK_NEXT;
}

//----------------------------------------------------------------------------------------
// Functions to communicate conserved variables between MeshBlocks

enum TaskStatus TimeIntegratorTaskList::HydroSend(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    pmb->pbval->SendCellCenteredBoundaryBuffers(pmb->phydro->u, HYDRO_CONS);
  } else {
    return TASK_FAIL;
  }
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::FieldSend(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    pmb->pbval->SendFieldBoundaryBuffers(pmb->pfield->b);
  } else {
    return TASK_FAIL;
  }
  return TASK_SUCCESS;
}

//----------------------------------------------------------------------------------------
// Functions to receive conserved variables between MeshBlocks

enum TaskStatus TimeIntegratorTaskList::HydroReceive(MeshBlock *pmb, int stage) {
  bool ret;
  if (stage <= nstages) {
    ret=pmb->pbval->ReceiveCellCenteredBoundaryBuffers(pmb->phydro->u, HYDRO_CONS);
  } else {
    return TASK_FAIL;
  }

  if (ret==true) {
    return TASK_SUCCESS;
  } else {
    return TASK_FAIL;
  }
}

enum TaskStatus TimeIntegratorTaskList::FieldReceive(MeshBlock *pmb, int stage) {
  bool ret;
  if (stage <= nstages) {
    ret=pmb->pbval->ReceiveFieldBoundaryBuffers(pmb->pfield->b);
  } else {
    return TASK_FAIL;
  }

  if (ret==true) {
    return TASK_SUCCESS;
  } else {
    return TASK_FAIL;
  }
}

enum TaskStatus TimeIntegratorTaskList::HydroShearSend(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    pmb->pbval->SendHydroShearingboxBoundaryBuffers(pmb->phydro->u, true);
  } else {
    return TASK_FAIL;
  }
  return TASK_SUCCESS;
}
enum TaskStatus TimeIntegratorTaskList::HydroShearReceive(MeshBlock *pmb, int stage) {
  bool ret;
  if (stage <= nstages) {
    ret=pmb->pbval->ReceiveHydroShearingboxBoundaryBuffers(pmb->phydro->u);
  } else {
    return TASK_FAIL;
  }

  if (ret==true) {
    return TASK_SUCCESS;
  } else {
    return TASK_FAIL;
  }
}
enum TaskStatus TimeIntegratorTaskList::FieldShearSend(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    pmb->pbval->SendFieldShearingboxBoundaryBuffers(pmb->pfield->b, true);
  } else {
    return TASK_FAIL;
  }
  return TASK_SUCCESS;
}
enum TaskStatus TimeIntegratorTaskList::FieldShearReceive(MeshBlock *pmb, int stage) {
  bool ret;
  if (stage <= nstages) {
    ret=pmb->pbval->ReceiveFieldShearingboxBoundaryBuffers(pmb->pfield->b);
  } else {
    return TASK_FAIL;
  }
  if (ret==true) {
    return TASK_SUCCESS;
  } else {
    return TASK_FAIL;
  }
}
enum TaskStatus TimeIntegratorTaskList::EMFShearSend(MeshBlock *pmb, int stage) {
  pmb->pbval->SendEMFShearingboxBoundaryCorrection();
  return TASK_SUCCESS;
}
enum TaskStatus TimeIntegratorTaskList::EMFShearReceive(MeshBlock *pmb, int stage) {
  if (pmb->pbval->ReceiveEMFShearingboxBoundaryCorrection() == true) {
    return TASK_NEXT;
  } else {
    return TASK_FAIL;
  }
}
enum TaskStatus TimeIntegratorTaskList::EMFShearRemap(MeshBlock *pmb, int stage) {
  pmb->pbval->RemapEMFShearingboxBoundary();
  return TASK_SUCCESS;
}

//--------------------------------------------------------------------------------------
// Functions for everything else

enum TaskStatus TimeIntegratorTaskList::Prolongation(MeshBlock *pmb, int stage) {
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;

  if (stage <= nstages) {
    // Time at the end of stage for (u, b) register pair
    Real t_end_stage = pmb->pmy_mesh->time + pmb->stage_abscissae[stage][0];
    // Scaled coefficient for RHS time-advance within stage
    Real dt = (stage_wghts[(stage-1)].beta)*(pmb->pmy_mesh->dt);
    pbval->ProlongateBoundaries(phydro->w,  phydro->u,  pfield->b,  pfield->bcc,
                                t_end_stage, dt);
  } else {
    return TASK_FAIL;
  }

  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::Primitives(MeshBlock *pmb, int stage) {
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;
  int il=pmb->is, iu=pmb->ie, jl=pmb->js, ju=pmb->je, kl=pmb->ks, ku=pmb->ke;
  if (pbval->nblevel[1][1][0] != -1) il-=NGHOST;
  if (pbval->nblevel[1][1][2] != -1) iu+=NGHOST;
  if (pbval->nblevel[1][0][1] != -1) jl-=NGHOST;
  if (pbval->nblevel[1][2][1] != -1) ju+=NGHOST;
  if (pbval->nblevel[0][1][1] != -1) kl-=NGHOST;
  if (pbval->nblevel[2][1][1] != -1) ku+=NGHOST;

  if (stage <= nstages) {
    // At beginning of this task, phydro->w contains previous stage's W(U) output
    // and phydro->w1 is used as a register to store the current stage's output.
    // For the second order integrators VL2 and RK2, the prim_old initial guess for the
    // Newton-Raphson solver in GR EOS uses the following abscissae:
    // stage=1: W at t^n and
    // stage=2: W at t^{n+1/2} (VL2) or t^{n+1} (RK2)
    pmb->peos->ConservedToPrimitive(phydro->u, phydro->w, pfield->b,
                                    phydro->w1, pfield->bcc, pmb->pcoord,
                                    il, iu, jl, ju, kl, ku);
    // swap AthenaArray data pointers so that w now contains the updated w_out
    phydro->w.SwapAthenaArray(phydro->w1);
  } else {
    return TASK_FAIL;
  }

  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::PhysicalBoundary(MeshBlock *pmb, int stage) {
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;

  if (stage <= nstages) {
    // Time at the end of stage for (u, b) register pair
    Real t_end_stage = pmb->pmy_mesh->time + pmb->stage_abscissae[stage][0];
    // Scaled coefficient for RHS time-advance within stage
    Real dt = (stage_wghts[(stage-1)].beta)*(pmb->pmy_mesh->dt);
    pbval->ApplyPhysicalBoundaries(phydro->w,  phydro->u,  pfield->b,  pfield->bcc,
                                   t_end_stage, dt);
  } else {
    return TASK_FAIL;
  }

  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::UserWork(MeshBlock *pmb, int stage) {
  if (stage != nstages) return TASK_SUCCESS; // only do on last stage

  pmb->UserWorkInLoop();
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::NewBlockTimeStep(MeshBlock *pmb, int stage) {
  if (stage != nstages) return TASK_SUCCESS; // only do on last stage

  pmb->phydro->NewBlockTimeStep();
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::CheckRefinement(MeshBlock *pmb, int stage) {
  if (stage != nstages) return TASK_SUCCESS; // only do on last stage

  pmb->pmr->CheckRefinementCondition();
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::GravFluxCorrection(MeshBlock *pmb, int stage) {
  if (stage != nstages) return TASK_SUCCESS; // only do on last stage

  pmb->phydro->CorrectGravityFlux();
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::StartupIntegrator(MeshBlock *pmb, int stage) {
  // Initialize time-integrator only on first stage
  if (stage != 1) {
    return TASK_SUCCESS;
  } else {
    // For each Meshblock, initialize time abscissae of each memory register pair (u,b)
    // at stage=0 to correspond to the beginning of the interval [t^n, t^{n+1}]
    pmb->stage_abscissae[0][0] = 0.0;
    pmb->stage_abscissae[0][1] = 0.0; // u1 advances to u1 = 0*u1 + 1.0*u in stage=1
    pmb->stage_abscissae[0][2] = 0.0; // u2 = u cached for all stages in 3S* methods

    // Given overall timestep dt, precompute the time abscissae for all registers, stages
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
    Hydro *ph=pmb->phydro;
    // Cache U^n in third memory register, u2, via deep copy
    // (if using a 3S* time-integrator)
    // ph->u2 = ph->u;

    if (MAGNETIC_FIELDS_ENABLED) { // MHD
      Field *pf=pmb->pfield;
      // Cache face-averaged B^n in third memory register, b2, via AthenaArray deep copy
      // (if using a 3S* time-integrator)
      // pf->b2.x1f = pf->b.x1f;
      // pf->b2.x2f = pf->b.x2f;
      // pf->b2.x3f = pf->b.x3f;

      // 2nd set of registers, including b1, need to be initialized to 0 each cycle
      Real ave_wghts[3];
      ave_wghts[0] = 0.0;
      ave_wghts[1] = 0.0;
      ave_wghts[2] = 0.0;
      pf->WeightedAveB(pf->b1,pf->b,pf->b,ave_wghts);
    }
    // 2nd set of registers, including u1, need to be initialized to 0 each cycle
    Real ave_wghts[3];
    ave_wghts[0] = 0.0;
    ave_wghts[1] = 0.0;
    ave_wghts[2] = 0.0;
    ph->WeightedAveU(ph->u1,ph->u,ph->u,ave_wghts);
    return TASK_SUCCESS;
  }
}
