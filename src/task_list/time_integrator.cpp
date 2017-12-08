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
#include "task_list.hpp"
#include "../athena.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../bvals/bvals.hpp"
#include "../gravity/gravity.hpp"
#include "../eos/eos.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"

//----------------------------------------------------------------------------------------
//  TimeIntegratorTaskList constructor

TimeIntegratorTaskList::TimeIntegratorTaskList(ParameterInput *pin, Mesh *pm)
  : TaskList(pm)
{
  // First, set weights for each step of time-integration algorithm.  Each step is
  //    U^{2} = a*U^0 + b*U^1 + c*dt*Div(F), where U^0 and U^1 are previous steps
  // a,b=(1-a),and c are weights that are different for each step and each integrator
  // These are stored as: time_int_wght1 = a, time_int_wght2 = b, time_int_wght3 = c

  integrator = pin->GetOrAddString("time","integrator","vl2");
  int dim = 1;
  if (pm->mesh_size.nx2 > 1) dim = 2;
  if (pm->mesh_size.nx3 > 1) dim = 3;

  if (integrator == "vl2") {
    // VL: second-order van Leer integrator (Stone & Gardiner, NewA 14, 139 2009)
    // Simple predictor-corrector scheme similar to MUSCL-Hancock
    // Expressed in 2S or 3S* algorithm form
    nsub_steps = 2;
    cfl_limit = 1.0;
    // Modify VL2 stability limit in 2D, 3D
    if (dim == 2) cfl_limit = 0.5;
    if (dim == 3) cfl_limit = 1.0/3.0;

    step_wghts[0].delta = 1.0; // required for consistency
    step_wghts[0].gamma_1 = 0.0;
    step_wghts[0].gamma_2 = 1.0;
    step_wghts[0].gamma_3 = 0.0;
    step_wghts[0].beta = 0.5;

    step_wghts[1].delta = 0.0;
    step_wghts[1].gamma_1 = 0.0;
    step_wghts[1].gamma_2 = 1.0;
    step_wghts[1].gamma_3 = 0.0;
    step_wghts[1].beta = 1.0;
  } else if (integrator == "rk2") {
    // Heun's method / SSPRK (2,2): Gottlieb (2009) equation 3.1
    // Optimal (in error bounds) explicit two-stage, second-order SSPRK
    nsub_steps = 2;
    cfl_limit = 1.0;
    step_wghts[0].delta = 1.0;
    step_wghts[0].gamma_1 = 0.0;
    step_wghts[0].gamma_2 = 1.0;
    step_wghts[0].gamma_3 = 0.0;
    step_wghts[0].beta = 1.0;

    step_wghts[1].delta = 0.0;
    step_wghts[1].gamma_1 = 0.5;
    step_wghts[1].gamma_2 = 0.5;
    step_wghts[1].gamma_3 = 0.0;
    step_wghts[1].beta = 0.5;
  } else if (integrator == "rk3") {
    // SSPRK (3,3): Gottlieb (2009) equation 3.2
    // Optimal (in error bounds) explicit three-stage, third-order SSPRK
    nsub_steps = 3;
    cfl_limit = 1.0;
    step_wghts[0].delta = 1.0;
    step_wghts[0].gamma_1 = 0.0;
    step_wghts[0].gamma_2 = 1.0;
    step_wghts[0].gamma_3 = 0.0;
    step_wghts[0].beta = 1.0;

    step_wghts[1].delta = 0.0;
    step_wghts[1].gamma_1 = 0.25;
    step_wghts[1].gamma_2 = 0.75;
    step_wghts[1].gamma_3 = 0.0;
    step_wghts[1].beta = 0.25;

    step_wghts[2].delta = 0.0;
    step_wghts[2].gamma_1 = TWO_3RD;
    step_wghts[2].gamma_2 = ONE_3RD;
    step_wghts[2].gamma_3 = 0.0;
    step_wghts[2].beta = TWO_3RD;
    //} else if (integrator == "ssprk5_3") {
    //} else if (integrator == "ssprk10_4") {
  } else if (integrator == "rk4") {
    // RK4()4[2S] from Table 2 of Ketchenson (2010)
    // Non-SSP, explicit four-stage, fourth-order RK
    nsub_steps = 4;
    // Stability properties are similar to classical RK4
    // Refer to Colella (2011) for constant advection with 4th order fluxes
    // linear stability analysis
    cfl_limit = 1.3925;
    step_wghts[0].delta = 1.0;
    step_wghts[0].gamma_1 = 0.0;
    step_wghts[0].gamma_2 = 1.0;
    step_wghts[0].gamma_3 = 0.0;
    step_wghts[0].beta = 1.193743905974738;

    step_wghts[1].delta = 0.217683334308543;
    step_wghts[1].gamma_1 = 0.121098479554482;
    step_wghts[1].gamma_2 = 0.721781678111411;
    step_wghts[1].gamma_3 = 0.0;
    step_wghts[1].beta = 0.099279895495783;

    step_wghts[2].delta = 1.065841341361089;
    step_wghts[2].gamma_1 = -3.843833699660025;
    step_wghts[2].gamma_2 = 2.121209265338722;
    step_wghts[2].gamma_3 = 0.0;
    step_wghts[2].beta = 1.131678018054042;

    step_wghts[3].delta = 0.0;
    step_wghts[3].gamma_1 = 0.546370891121863;
    step_wghts[3].gamma_2 = 0.198653035682705;
    step_wghts[3].gamma_3 = 0.0;
    step_wghts[3].beta = 0.310665766509336;
  } else if (integrator == "ssprk5_4") {
    // SSPRK (5,4): Gottlieb (2009) section 3.1
    // Optimal (in error bounds) explicit five-stage, fourth-order SSPRK
    // 3N method, but there is no 3S* formulation due to irregular sparsity
    // of Shu-Osher form matrix, alpha
    nsub_steps = 5;
    cfl_limit = 1.3925;
    step_wghts[0].delta = 1.0;
    step_wghts[0].gamma_1 = 0.0;
    step_wghts[0].gamma_2 = 1.0;
    step_wghts[0].gamma_3 = 0.0;
    step_wghts[0].beta = 0.391752226571890;

    step_wghts[1].delta = 0.0; // u1 = u^n
    step_wghts[1].gamma_1 = 0.555629506348765;
    step_wghts[1].gamma_2 = 0.444370493651235;
    step_wghts[1].gamma_3 = 0.0;
    step_wghts[1].beta = 0.368410593050371;

    step_wghts[2].delta = 0.0;
    step_wghts[2].gamma_1 = 0.379898148511597;
    step_wghts[2].gamma_2 = 0.0;
    step_wghts[2].gamma_3 = 0.620101851488403; // u2 = u^n
    step_wghts[2].beta = 0.251891774271694;

    step_wghts[3].delta = 0.0;
    step_wghts[3].gamma_1 = TWO_3RD;
    step_wghts[3].gamma_2 = ONE_3RD;
    step_wghts[3].gamma_3 = 0.178079954393132; // u2 = u^n
    step_wghts[3].beta = 0.544974750228521;

    step_wghts[4].delta = 0.0;
    step_wghts[4].gamma_1 = 0.386708617503268; // from Gottlieb (2009), u^(4) coeff.
    step_wghts[4].gamma_2 = ONE_3RD;
    step_wghts[4].gamma_3 = 0.0;
    step_wghts[4].beta = 0.226007483236906; // from Gottlieb (2009), F(u^(4)) coeff.
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in CreateTimeIntegrator" << std::endl
        << "integrator=" << integrator << " not valid time integrator" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // Set cfl_number based on user input and time integrator CFL limit
  Real cfl_number = pin->GetReal("time","cfl_number");
  if(cfl_number > cfl_limit) {
    std::cout << "### Warning in CreateTimeIntegrator" << std::endl
        << "User CFL number " << cfl_number << " must be smaller than " << cfl_limit
        << " for integrator=" << integrator << " in "
        << dim << "D simulation" << std::endl << "Setting to limit" << std::endl;
    cfl_number = cfl_limit;
  }
  // Save to Mesh class
  pm->cfl_number = cfl_number;

  // Now assemble list of tasks for each step of time integrator
  {using namespace HydroIntegratorTaskNames;
    AddTimeIntegratorTask(STARTUP_INT,NONE);
    AddTimeIntegratorTask(START_ALLRECV,STARTUP_INT);

    // compute hydro fluxes, integrate hydro variables
    AddTimeIntegratorTask(CALC_HYDFLX,START_ALLRECV);
    if(pm->multilevel==true) { // SMR or AMR
      AddTimeIntegratorTask(SEND_HYDFLX,CALC_HYDFLX);
      AddTimeIntegratorTask(RECV_HYDFLX,CALC_HYDFLX);
      AddTimeIntegratorTask(INT_HYD, RECV_HYDFLX);
    } else {
      AddTimeIntegratorTask(INT_HYD, CALC_HYDFLX);
    }
    AddTimeIntegratorTask(UPDATE_DT,INT_HYD);
    AddTimeIntegratorTask(SRCTERM_HYD,UPDATE_DT);
    AddTimeIntegratorTask(SEND_HYD,SRCTERM_HYD);
    AddTimeIntegratorTask(RECV_HYD,START_ALLRECV);

    // compute MHD fluxes, integrate field
    if (MAGNETIC_FIELDS_ENABLED) { // MHD
      AddTimeIntegratorTask(CALC_FLDFLX,CALC_HYDFLX);
      AddTimeIntegratorTask(SEND_FLDFLX,CALC_FLDFLX);
      AddTimeIntegratorTask(RECV_FLDFLX,SEND_FLDFLX);
      AddTimeIntegratorTask(INT_FLD, RECV_FLDFLX);
      AddTimeIntegratorTask(SEND_FLD,INT_FLD);
      AddTimeIntegratorTask(RECV_FLD,START_ALLRECV);
    }

    // prolongate, compute new primitives
    if (MAGNETIC_FIELDS_ENABLED) { // MHD
      if(pm->multilevel==true) { // SMR or AMR
        AddTimeIntegratorTask(PROLONG, (SEND_HYD|RECV_HYD|SEND_FLD|RECV_FLD));
        AddTimeIntegratorTask(CON2PRIM,PROLONG);
      } else {
        AddTimeIntegratorTask(CON2PRIM,(INT_HYD|RECV_HYD|INT_FLD|RECV_FLD));
      }
    } else {  // HYDRO
      if(pm->multilevel==true) { // SMR or AMR
        AddTimeIntegratorTask(PROLONG,(SEND_HYD|RECV_HYD));
        AddTimeIntegratorTask(CON2PRIM,PROLONG);
      } else {
        AddTimeIntegratorTask(CON2PRIM,(INT_HYD|RECV_HYD));
      }
    }

    // everything else
    AddTimeIntegratorTask(PHY_BVAL,CON2PRIM);
//    if (SELF_GRAVITY_ENABLED == 1){
//      AddTimeIntegratorTask(CORR_GFLX,PHY_BVAL);
//      AddTimeIntegratorTask(USERWORK,CORR_GFLX);
//    } else {
    AddTimeIntegratorTask(USERWORK,PHY_BVAL);
//    }
    AddTimeIntegratorTask(NEW_DT,USERWORK);
    if(pm->adaptive==true) {
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

void TimeIntegratorTaskList::AddTimeIntegratorTask(uint64_t id, uint64_t dep)
{
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;

  using namespace HydroIntegratorTaskNames;
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

    case (UPDATE_DT):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::UpdateTimeStep);
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
//! \fn
//  \brief

//----------------------------------------------------------------------------------------
// Functions to start/end MPI communication

enum TaskStatus TimeIntegratorTaskList::StartAllReceive(MeshBlock *pmb, int step)
{
  pmb->pbval->StartReceivingAll();
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::ClearAllBoundary(MeshBlock *pmb, int step)
{
  pmb->pbval->ClearBoundaryAll();
  return TASK_SUCCESS;
}

//----------------------------------------------------------------------------------------
// Functions to calculates fluxes

enum TaskStatus TimeIntegratorTaskList::CalculateFluxes(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;

  if((step == 1) && (integrator == "vl2")) {
    phydro->CalculateFluxes(phydro->w,  pfield->b,  pfield->bcc, true);
    return TASK_NEXT;
  }

  else if (step != 1 && step <= nsub_steps) {
    phydro->CalculateFluxes(phydro->w,  pfield->b,  pfield->bcc, false);
    return TASK_NEXT;
  }

  return TASK_FAIL;
}

enum TaskStatus TimeIntegratorTaskList::CalculateEMF(MeshBlock *pmb, int step)
{
  if (step <= nsub_steps) {
    pmb->pfield->ComputeCornerE(pmb->phydro->w,  pmb->pfield->bcc);
    return TASK_NEXT;
  }

  return TASK_FAIL;
}

//----------------------------------------------------------------------------------------
// Functions to communicate fluxes between MeshBlocks for flux correction step with AMR

enum TaskStatus TimeIntegratorTaskList::FluxCorrectSend(MeshBlock *pmb, int step)
{
  pmb->pbval->SendFluxCorrection(FLUX_HYDRO);
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::EMFCorrectSend(MeshBlock *pmb, int step)
{
  pmb->pbval->SendEMFCorrection();
  return TASK_SUCCESS;
}

//----------------------------------------------------------------------------------------
// Functions to receive fluxes between MeshBlocks

enum TaskStatus TimeIntegratorTaskList::FluxCorrectReceive(MeshBlock *pmb, int step)
{
  if(pmb->pbval->ReceiveFluxCorrection(FLUX_HYDRO) == true) {
    return TASK_NEXT;
  } else {
    return TASK_FAIL;
  }
}

enum TaskStatus TimeIntegratorTaskList::EMFCorrectReceive(MeshBlock *pmb, int step)
{
  if(pmb->pbval->ReceiveEMFCorrection() == true) {
    return TASK_NEXT;
  } else {
    return TASK_FAIL;
  }
}

//----------------------------------------------------------------------------------------
// Functions to integrate conserved variables

enum TaskStatus TimeIntegratorTaskList::HydroIntegrate(MeshBlock *pmb, int step)
{
  Hydro *ph=pmb->phydro;
  Field *pf=pmb->pfield;
  if (step <= nsub_steps) {
    Real ave_wghts[3];
    ave_wghts[0] = 1.0;
    ave_wghts[1] = step_wghts[step-1].delta;
    ave_wghts[2] = 0.0;
    ph->WeightedAveU(ph->u1,ph->u,ph->u2,ave_wghts);

    ave_wghts[0] = step_wghts[step-1].gamma_1;
    ave_wghts[1] = step_wghts[step-1].gamma_2;
    ave_wghts[2] = step_wghts[step-1].gamma_3;
    ph->WeightedAveU(ph->u,ph->u1,ph->u2,ave_wghts);
    ph->AddFluxDivergenceToAverage(ph->w,pf->bcc,step_wghts[step-1].beta,ph->u);

    return TASK_NEXT;
  }

  return TASK_FAIL;
}

enum TaskStatus TimeIntegratorTaskList::FieldIntegrate(MeshBlock *pmb, int step)
{
  Field *pf=pmb->pfield;

  if (step <= nsub_steps) {
    // This time-integrator-specific averaging is redundant with Hydro
    Real ave_wghts[3];
    ave_wghts[0] = 1.0;
    ave_wghts[1] = step_wghts[step-1].delta;
    ave_wghts[2] = 0.0;
    pf->WeightedAveB(pf->b1,pf->b,pf->b2,ave_wghts);

    ave_wghts[0] = step_wghts[step-1].gamma_1;
    ave_wghts[1] = step_wghts[step-1].gamma_2;
    ave_wghts[2] = step_wghts[step-1].gamma_3;
    pf->WeightedAveB(pf->b,pf->b1,pf->b2,ave_wghts);

    pf->CT(step_wghts[step-1].beta, pf->b);

    return TASK_NEXT;
  }

  return TASK_FAIL;
}

//----------------------------------------------------------------------------------------
// Functions to add source terms

enum TaskStatus TimeIntegratorTaskList::HydroSourceTerms(MeshBlock *pmb, int step)
{
  Hydro *ph=pmb->phydro;
  Field *pf=pmb->pfield;

  // return if there are no source terms to be added
  if (ph->psrc->hydro_sourceterms_defined == false) return TASK_NEXT;

  if (step <= nsub_steps) {
    Real time=pmb->pmy_mesh->time;
    Real dt = step_dt[0];
    ph->psrc->AddHydroSourceTerms(time,dt,ph->flux,ph->w,pf->bcc,ph->u);
  } else {
    // Evaluate the source terms at the beginning of the
    return TASK_FAIL;
  }
  return TASK_NEXT;
}

//----------------------------------------------------------------------------------------
// Functions to communicate conserved variables between MeshBlocks

enum TaskStatus TimeIntegratorTaskList::HydroSend(MeshBlock *pmb, int step)
{
  if (step <= nsub_steps) {
    pmb->pbval->SendCellCenteredBoundaryBuffers(pmb->phydro->u, HYDRO_CONS);
  }
  else {
    return TASK_FAIL;
  }
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::FieldSend(MeshBlock *pmb, int step)
{
  if (step <= nsub_steps) {
    pmb->pbval->SendFieldBoundaryBuffers(pmb->pfield->b);
  }
  else {
    return TASK_FAIL;
  }
  return TASK_SUCCESS;
}

//----------------------------------------------------------------------------------------
// Functions to receive conserved variables between MeshBlocks

enum TaskStatus TimeIntegratorTaskList::HydroReceive(MeshBlock *pmb, int step)
{
  bool ret;
  if (step <= nsub_steps) {
    ret=pmb->pbval->ReceiveCellCenteredBoundaryBuffers(pmb->phydro->u, HYDRO_CONS);
  }
  else {
    return TASK_FAIL;
  }

  if(ret==true) {
    return TASK_SUCCESS;
  } else {
    return TASK_FAIL;
  }
}

enum TaskStatus TimeIntegratorTaskList::FieldReceive(MeshBlock *pmb, int step)
{
  bool ret;
  if (step <= nsub_steps) {
    ret=pmb->pbval->ReceiveFieldBoundaryBuffers(pmb->pfield->b);
  }
  else {
    return TASK_FAIL;
  }

  if(ret==true) {
    return TASK_SUCCESS;
  } else {
    return TASK_FAIL;
  }
}

//----------------------------------------------------------------------------------------
// Functions for everything else

enum TaskStatus TimeIntegratorTaskList::Prolongation(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;
  Real dt;

  if (step <= nsub_steps) {
    dt = step_dt[0];
    pbval->ProlongateBoundaries(phydro->w,  phydro->u,  pfield->b,  pfield->bcc,
                                pmb->pmy_mesh->time+dt, dt);
  } else {
    return TASK_FAIL;
  }

  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::Primitives(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  if(pbval->nblevel[1][1][0]!=-1) is-=NGHOST;
  if(pbval->nblevel[1][1][2]!=-1) ie+=NGHOST;
  if(pbval->nblevel[1][0][1]!=-1) js-=NGHOST;
  if(pbval->nblevel[1][2][1]!=-1) je+=NGHOST;
  if(pbval->nblevel[0][1][1]!=-1) ks-=NGHOST;
  if(pbval->nblevel[2][1][1]!=-1) ke+=NGHOST;

  if (step <= nsub_steps) {
    // Incompatible with GR right now due to hardcoded w_old=w1 usage
    // Cache w from previous substep in w1
    phydro->w1 = phydro->w;
    pmb->peos->ConservedToPrimitive(phydro->u, phydro->w1, pfield->b,
                                    phydro->w, pfield->bcc, pmb->pcoord,
                                    is, ie, js, je, ks, ke);
  }
  else {
    return TASK_FAIL;
  }
  // if(step == 1) {
  //   pmb->peos->ConservedToPrimitive(phydro->u1, phydro->w, pfield->b1,
  //                                   phydro->w1, pfield->bcc1, pmb->pcoord,
  //                                   is, ie, js, je, ks, ke);
  // } else if(step == 2) {
  //   pmb->peos->ConservedToPrimitive(phydro->u, phydro->w1, pfield->b,
  //                                   phydro->w, pfield->bcc, pmb->pcoord,
  //                                   is, ie, js, je, ks, ke);
  // } else {
  //   return TASK_FAIL;
  // }
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::PhysicalBoundary(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;
  Real dt;

  if (step <= nsub_steps) {
    dt = step_dt[0];
    pbval->ApplyPhysicalBoundaries(phydro->w,  phydro->u,  pfield->b,  pfield->bcc,
                                   pmb->pmy_mesh->time+dt, dt);
  }
  else {
    return TASK_FAIL;
  }

  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::UserWork(MeshBlock *pmb, int step)
{
  if (step != nsub_steps) return TASK_SUCCESS; // only do on last sub-step

  pmb->UserWorkInLoop();
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::NewBlockTimeStep(MeshBlock *pmb, int step)
{
  if (step != nsub_steps) return TASK_SUCCESS; // only do on last sub-step

  pmb->phydro->NewBlockTimeStep();
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::CheckRefinement(MeshBlock *pmb, int step)
{
  if (step != nsub_steps) return TASK_SUCCESS; // only do on last sub-step

  pmb->pmr->CheckRefinementCondition();
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::GravFluxCorrection(MeshBlock *pmb, int step)
{
  if (step != nsub_steps) return TASK_SUCCESS; // only do on last sub-step

  pmb->phydro->CorrectGravityFlux();
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::StartupIntegrator(MeshBlock *pmb, int step)
{
  // Initialize registers only on first sub-step
  if (step != 1) {
    return TASK_SUCCESS;
  }
  else {
    // if (nsub_steps <= 3) return TASK_SUCCESS; // not necessary for third-order or lower
    Hydro *ph=pmb->phydro;
    // Cache U^n in third memory register, u2, via deep copy
    ph->u2 = ph->u;

    if (MAGNETIC_FIELDS_ENABLED) { // MHD
      Field *pf=pmb->pfield;
      // Cache face-averaged B^n in third memory register, b2, via deep copy
      pf->b2.x1f = pf->b.x1f;
      pf->b2.x2f = pf->b.x2f;
      pf->b2.x3f = pf->b.x3f;
      // Cache cell-averaged B^n in third memory register, bcc2, via deep copy
      // However, bcc is not computed until end of substep, in W(U)
      //      pf->bcc2 = pf->bcc;

      // 2nd registers, including u1, need to be initialized to 0
      pf->b1.x1f = pf->b.x1f;
      pf->b1.x2f = pf->b.x2f;
      pf->b1.x3f = pf->b.x3f;
      Real ave_wghts[3];
      ave_wghts[0] = 0.0;
      ave_wghts[1] = 0.0;
      ave_wghts[2] = 0.0;
      pf->WeightedAveB(pf->b1,pf->b,pf->b,ave_wghts);
    }
    // 2nd registers, including u1, need to be initialized to 0
    ph->u1 = ph->u;
    Real ave_wghts[3];
    ave_wghts[0] = 0.0;
    ave_wghts[1] = 0.0;
    ave_wghts[2] = 0.0;
    ph->WeightedAveU(ph->u1,ph->u,ph->u,ave_wghts);
    return TASK_SUCCESS;
  }
}


enum TaskStatus TimeIntegratorTaskList::UpdateTimeStep(MeshBlock *pmb, int step){
  // Occurs after HydroIntegrate(), but before HydroSourceTerms() and FieldIntegrate()
  if (step != 1) {
    // Update the dt abscissae of each memory register to values at end of this substep
    Real dt, dt1, dt2;
    const IntegratorWeight w = step_wghts[step-1];
    // u1 = u1 + delta*u
    dt1 = step_dt[1] + w.delta*step_dt[0];
    // u = gamma_1*u + gamma_2*u1 + gamma_3*u2 + beta*dt*F(u)
    dt = w.gamma_1*step_dt[0] +
        w.gamma_2*dt1 +
        w.gamma_3*step_dt[2] +
        w.beta*pmb->pmy_mesh->dt;
    // u2 = u^n
    dt2 = 0.0;

    step_dt[0]= dt;
    step_dt[1]= dt1;
    step_dt[2]= dt2;
    return TASK_SUCCESS;
  }
  else {
    // Initialize the dt abscissae of each memory register
    step_dt[0]= 0.0;
    step_dt[1]= 0.0;
    step_dt[2]= 0.0;
    return TASK_SUCCESS;
  }
}
