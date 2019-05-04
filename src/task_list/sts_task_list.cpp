//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file sts_task_list.cpp
//  \brief RKL1 Super-Time-Stepping Algorithm
//
// REFERENCE:
// Meyer, C. D., Balsara, D. S., & Aslam, T. D. 2014, J. Comput. Phys.,
//    257A, 594-626

// C headers

// C++ headers
#include <cstring>    // strcmp()
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
//  SuperTimeStepTaskList constructor

SuperTimeStepTaskList::SuperTimeStepTaskList(
    ParameterInput *pin, Mesh *pm, TimeIntegratorTaskList *ptlist) : ptlist_(ptlist) {
  // STS Incompatiblities
  if (MAGNETIC_FIELDS_ENABLED
      && !(pm->pblock->pfield->fdif.field_diffusion_defined)
      && !(pm->pblock->phydro->hdif.hydro_diffusion_defined)) {
    std::stringstream msg;
    msg << "### FATAL ERROR in SuperTimeStepTaskList" << std::endl
        << "Super-time-stepping requires setting parameters for "
        << "diffusive processes in input file." << std::endl;
    ATHENA_ERROR(msg);
  }
  // TODO(pdmullen): time-dep BC's require knowing the time within
  //                 an RKL1 operator-split STS, what is the time?
  if (SHEARING_BOX) {
    std::stringstream msg;
    msg << "### FATAL ERROR in SuperTimeStepTaskList" << std::endl
        << "Super-time-stepping is not yet compatible "
        << "with shearing box BC's." << std::endl;
    ATHENA_ERROR(msg);
  }
  // TODO(pdmullen): how should source terms be handled inside
  //                 operator-split RKL1 STS?
  if (pm->pblock->phydro->hsrc.hydro_sourceterms_defined) {
    std::stringstream msg;
    msg << "### FATAL ERROR in SuperTimeStepTaskList" << std::endl
        << "Super-time-stepping is not yet compatible "
        << "with source terms." << std::endl;
    ATHENA_ERROR(msg);
  }
  // TODO(pdmullen): fix non-Cartesian compatibility; this requires the
  //                 handling of coordinate source terms.
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") != 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in SuperTimeStepTaskList" << std::endl
        << "Super-time-stepping is not yet compatibile "
        << "with non-Cartesian coordinates." << std::endl;
    ATHENA_ERROR(msg);
  }
  // TODO(pdmullen): add mesh-refinement functionality
  if (pm->multilevel) {
    std::stringstream msg;
    msg << "### FATAL ERROR in SuperTimeStepTaskList" << std::endl
        << "Super-time-stepping is not yet compatibile "
        << "with mesh refinement." << std::endl;
    ATHENA_ERROR(msg);
  }

  // Now assemble list of tasks for each stage of SuperTimeStep integrator
  {using namespace HydroIntegratorTaskNames; // NOLINT (build/namespace)
    // calculate hydro/field diffusive fluxes
    AddTask(DIFFUSE_HYD,NONE);
    if (NSCALARS > 0)
      AddTask(DIFFUSE_SCLR,NONE);
    // compute hydro fluxes, integrate hydro variables
    if (MAGNETIC_FIELDS_ENABLED) {
      AddTask(DIFFUSE_FLD,NONE);
      AddTask(CALC_HYDFLX,(DIFFUSE_HYD|DIFFUSE_FLD));
    } else {
      AddTask(CALC_HYDFLX,DIFFUSE_HYD);
    }
    AddTask(INT_HYD, CALC_HYDFLX);
    AddTask(SEND_HYD,INT_HYD);
    AddTask(RECV_HYD,NONE);
    AddTask(SETB_HYD,(RECV_HYD|INT_HYD));

    // compute passive scalar fluxes and integrate the density of each species:
    if (NSCALARS > 0) {
      AddTask(CALC_SCLRFLX,(CALC_HYDFLX|DIFFUSE_SCLR));
      // TODO(felker): uncomment after S/AMR prohibitiion for STS is removed
      // if (pm->multilevel) {
      //   AddTask(SEND_SCLRFLX,CALC_SCLRFLX);
      //   AddTask(RECV_SCLRFLX,CALC_SCLRFLX);
      //   AddTask(INT_SCLR,RECV_SCLRFLX);
      // } else {
      AddTask(INT_SCLR, CALC_SCLRFLX);
        //}
      AddTask(SEND_SCLR,INT_SCLR);
      AddTask(RECV_SCLR,NONE);
      AddTask(SETB_SCLR,RECV_SCLR);
    }

    // compute MHD fluxes, integrate field
    if (MAGNETIC_FIELDS_ENABLED) { // MHD
      AddTask(CALC_FLDFLX,CALC_HYDFLX);
      AddTask(SEND_FLDFLX,CALC_FLDFLX);
      AddTask(RECV_FLDFLX,SEND_FLDFLX);
      AddTask(INT_FLD,RECV_FLDFLX);
      AddTask(SEND_FLD,INT_FLD);
      AddTask(RECV_FLD,NONE);
      AddTask(SETB_FLD,(RECV_FLD|INT_FLD));
      // compute new primitives
      if (NSCALARS > 0) {
        AddTask(CONS2PRIM,(SETB_HYD|SETB_FLD|SETB_SCLR));
      } else {
        AddTask(CONS2PRIM,(SETB_HYD|SETB_FLD));
      }
    } else {  // Hydro
      if (NSCALARS > 0) {
        AddTask(CONS2PRIM,(SETB_HYD|SETB_SCLR));
      } else {
        AddTask(CONS2PRIM,SETB_HYD);
      }
    }

    // everything else
    AddTask(PHY_BVAL,CONS2PRIM);
    AddTask(CLEAR_ALLBND,PHY_BVAL);
  } // end of using namespace block
}

//---------------------------------------------------------------------------------------
//  Sets id and dependency for "ntask" member of task_list_ array, then iterates value of
//  ntask.

void SuperTimeStepTaskList::AddTask(std::uint64_t id, std::uint64_t dep) {
  task_list_[ntasks].task_id = id;
  task_list_[ntasks].dependency = dep;

  using namespace HydroIntegratorTaskNames; // NOLINT (build/namespace)
  switch (id) {
    case (CLEAR_ALLBND):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&TimeIntegratorTaskList::ClearAllBoundary);
      task_list_[ntasks].lb_time = false;
      break;

    case (CALC_HYDFLX):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&SuperTimeStepTaskList::CalculateHydroFlux_STS);
      task_list_[ntasks].lb_time = true;
      break;
    case (CALC_FLDFLX):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&SuperTimeStepTaskList::CalculateEMF_STS);
      task_list_[ntasks].lb_time = true;
      break;

    case (SEND_FLDFLX):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&TimeIntegratorTaskList::SendEMF);
      task_list_[ntasks].lb_time = true;
      break;

    case (RECV_FLDFLX):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&TimeIntegratorTaskList::ReceiveAndCorrectEMF);
      task_list_[ntasks].lb_time = false;
      break;

    case (INT_HYD):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&SuperTimeStepTaskList::IntegrateHydro_STS);
      task_list_[ntasks].lb_time = true;
      break;
    case (INT_FLD):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&SuperTimeStepTaskList::IntegrateField_STS);
      task_list_[ntasks].lb_time = true;
      break;

    case (SEND_HYD):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&TimeIntegratorTaskList::SendHydro);
      task_list_[ntasks].lb_time = true;
      break;
    case (SEND_FLD):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&TimeIntegratorTaskList::SendField);
      task_list_[ntasks].lb_time = true;
      break;

    case (RECV_HYD):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&TimeIntegratorTaskList::ReceiveHydro);
      task_list_[ntasks].lb_time = false;
      break;
    case (RECV_FLD):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&TimeIntegratorTaskList::ReceiveField);
      task_list_[ntasks].lb_time = false;
      break;

    case (SETB_HYD):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&TimeIntegratorTaskList::SetBoundariesHydro);
      task_list_[ntasks].lb_time = true;
      break;
    case (SETB_FLD):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&TimeIntegratorTaskList::SetBoundariesField);
      task_list_[ntasks].lb_time = true;
      break;

    case (CALC_SCLRFLX):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&SuperTimeStepTaskList::CalculateScalarFlux_STS);
      task_list_[ntasks].lb_time = true;
      break;
    case (INT_SCLR):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&SuperTimeStepTaskList::IntegrateScalars_STS);
      task_list_[ntasks].lb_time = true;
      break;
    case (SEND_SCLR):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&TimeIntegratorTaskList::SendScalars);
      task_list_[ntasks].lb_time = true;
      break;
    case (RECV_SCLR):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&TimeIntegratorTaskList::ReceiveScalars);
      task_list_[ntasks].lb_time = false;
      break;
    case (SETB_SCLR):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&TimeIntegratorTaskList::SetBoundariesScalars);
      task_list_[ntasks].lb_time = true;
      break;

    case (CONS2PRIM):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&TimeIntegratorTaskList::Primitives);
      task_list_[ntasks].lb_time = true;
      break;
    case (PHY_BVAL):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&SuperTimeStepTaskList::PhysicalBoundary_STS);
      task_list_[ntasks].lb_time = true;
      break;

    case (DIFFUSE_HYD):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&TimeIntegratorTaskList::DiffuseHydro);
      task_list_[ntasks].lb_time = true;
      break;
    case (DIFFUSE_FLD):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&TimeIntegratorTaskList::DiffuseField);
      task_list_[ntasks].lb_time = true;
      break;
    case (DIFFUSE_SCLR):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&TimeIntegratorTaskList::DiffuseScalars);
      task_list_[ntasks].lb_time = true;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in AddTask" << std::endl
          << "Invalid Task "<< id << " is specified" << std::endl;
      ATHENA_ERROR(msg);
  }
  ntasks++;
  return;
}


void SuperTimeStepTaskList::StartupTaskList(MeshBlock *pmb, int stage) {
#pragma omp single
  {
    // Set RKL1 params
    pmb->pmy_mesh->muj = (2.*stage-1.)/stage;
    pmb->pmy_mesh->nuj = (1.-stage)/stage;
    pmb->pmy_mesh->muj_tilde = pmb->pmy_mesh->muj*2./(std::pow(nstages,2.)+nstages);
  }

  // Clear flux arrays from previous stage
  pmb->phydro->hdif.ClearFlux(pmb->phydro->flux);
  if (MAGNETIC_FIELDS_ENABLED)
    pmb->pfield->fdif.ClearEMF(pmb->pfield->e);
  if (NSCALARS > 0) {
    PassiveScalars *ps = pmb->pscalars;
    ps->diffusion_flx[X1DIR].ZeroClear();
    ps->diffusion_flx[X2DIR].ZeroClear();
    ps->diffusion_flx[X3DIR].ZeroClear();
  }

  // Real time = pmb->pmy_mesh->time;
  pmb->pbval->StartReceiving(BoundaryCommSubset::all);

  return;
}

//----------------------------------------------------------------------------------------
// Functions to calculates fluxes

TaskStatus SuperTimeStepTaskList::CalculateHydroFlux_STS(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;
  if (stage <= nstages) {
    ph->CalculateFluxes_STS();
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}


TaskStatus SuperTimeStepTaskList::CalculateScalarFlux_STS(MeshBlock *pmb, int stage) {
  PassiveScalars *ps = pmb->pscalars;
  if (stage <= nstages) {
    std::cout << "In STS CALC_SCLRFLX, stage= " << stage << std::endl;
    ps->CalculateFluxes_STS();
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}


TaskStatus SuperTimeStepTaskList::CalculateEMF_STS(MeshBlock *pmb, int stage) {
  Field *pf = pmb->pfield;
  if (stage <= nstages) {
    pf->ComputeCornerE_STS();
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}

//----------------------------------------------------------------------------------------
// Functions to integrate conserved variables

TaskStatus SuperTimeStepTaskList::IntegrateHydro_STS(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;
  Field *pf = pmb->pfield;

  // set registers
  ph->u2.SwapAthenaArray(ph->u1);
  ph->u1.SwapAthenaArray(ph->u);
  if (MAGNETIC_FIELDS_ENABLED) {
    pf->b2.x1f.SwapAthenaArray(pf->b1.x1f);
    pf->b2.x2f.SwapAthenaArray(pf->b1.x2f);
    pf->b2.x3f.SwapAthenaArray(pf->b1.x3f);

    pf->b1.x1f.SwapAthenaArray(pf->b.x1f);
    pf->b1.x2f.SwapAthenaArray(pf->b.x2f);
    pf->b1.x3f.SwapAthenaArray(pf->b.x3f);
  }

  // update u
  if (stage <= nstages) {
    Real ave_wghts[3];
    ave_wghts[0] = 0.0;
    ave_wghts[1] = pmb->pmy_mesh->muj;
    ave_wghts[2] = pmb->pmy_mesh->nuj;
    pmb->WeightedAve(ph->u, ph->u1, ph->u2, ave_wghts);
    Real wght = pmb->pmy_mesh->muj_tilde;
    ph->AddFluxDivergence(wght, ph->u);
    // TODO(pdmullen): check this after disabling ATHENA_ERROR for src terms
    pmb->pcoord->AddCoordTermsDivergence(wght, ph->flux, ph->w, pf->bcc, ph->u);
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}


TaskStatus SuperTimeStepTaskList::IntegrateScalars_STS(MeshBlock *pmb, int stage) {
  PassiveScalars *ps = pmb->pscalars;
  // set registers
  ps->s2.SwapAthenaArray(ps->s1);
  ps->s1.SwapAthenaArray(ps->s);

  // update u
  if (stage <= nstages) {
    Real ave_wghts[3];
    ave_wghts[0] = 0.0;
    ave_wghts[1] = pmb->pmy_mesh->muj;
    ave_wghts[2] = pmb->pmy_mesh->nuj;
    pmb->WeightedAve(ps->s, ps->s1, ps->s2, ave_wghts);
    Real wght = pmb->pmy_mesh->muj_tilde;
    ps->AddFluxDivergence(wght, ps->s);

    return TaskStatus::next;
  }
  return TaskStatus::fail;
}


TaskStatus SuperTimeStepTaskList::IntegrateField_STS(MeshBlock *pmb, int stage) {
  Field *pf = pmb->pfield;

  if (stage <= nstages) {
    Real ave_wghts[3];
    ave_wghts[0] = 0.0;
    ave_wghts[1] = pmb->pmy_mesh->muj;
    ave_wghts[2] = pmb->pmy_mesh->nuj;
    pmb->WeightedAve(pf->b, pf->b1, pf->b2, ave_wghts);
    pf->CT(pmb->pmy_mesh->muj_tilde, pf->b);
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}

TaskStatus SuperTimeStepTaskList::PhysicalBoundary_STS(MeshBlock *pmb, int stage) {
  BoundaryValues *pbval=pmb->pbval;
  if (stage <= nstages) {
    // TODO(pdmullen): for time-dependent BC's, what is the time inside of an
    //                 operator-split RKL1 STS? For now, disable time-dep BCs.
    // Real t_end_stage = pmb->pmy_mesh->time;
    // Real dt = pmb->pmy_mesh->dt;
    pmb->phydro->hbvar.SwapHydroQuantity(pmb->phydro->w, HydroBoundaryQuantity::prim);
    pbval->ApplyPhysicalBoundaries(pmb->pmy_mesh->time, pmb->pmy_mesh->dt);
  } else {
    return TaskStatus::fail;
  }

  return TaskStatus::success;
}
