//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file time_integrator.hpp
//  \brief implements task list for any time integrator (e.g. van Leer, RK2, RK3, etc.)
//======================================================================================

// C/C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ classes headers
#include "../athena.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../bvals/bvals.hpp"
#include "../eos/eos.hpp"
#include "../hydro/srcterms/srcterms.hpp"

// this class header
#include "task_list.hpp"

//--------------------------------------------------------------------------------------
//! \fn
//  \brief

namespace TaskFunctions {
enum TaskStatus CalculateFluxes(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;

  if(step == 1) {
//    phydro->u1 = phydro->u;
    phydro->CopyOrAverageHydro(phydro->u1, phydro->u, phydro->u, 0.0);

    phydro->CalculateFluxes(pmb, phydro->u1, phydro->w, pfield->b,
                                           pfield->bcc, 1);
  } else if(step == 2) {
    phydro->CalculateFluxes(pmb, phydro->u, phydro->w1, pfield->b1,
                                           pfield->bcc1, 2);
  } else {
    return TASK_FAIL;
  }

  return TASK_NEXT;
}

enum TaskStatus FluxCorrectSend(MeshBlock *pmb, int step)
{
  pmb->pbval->SendFluxCorrection();
  return TASK_SUCCESS;
}

enum TaskStatus FluxCorrectReceive(MeshBlock *pmb, int step)
{
  if(pmb->pbval->ReceiveFluxCorrection() == true) {
    return TASK_NEXT;
  } else {
    return TASK_FAIL;
  }
}

enum TaskStatus CalculateEMF(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  if(step == 1) {
    pfield->ComputeCornerE(pmb, phydro->w, pfield->bcc);
  } else if(step == 2) {
    pfield->ComputeCornerE(pmb, phydro->w1, pfield->bcc1);
  } else {
    return TASK_FAIL;
  }
  return TASK_NEXT;
}

enum TaskStatus EMFCorrectSend(MeshBlock *pmb, int step)
{
  pmb->pbval->SendEMFCorrection();
  return TASK_SUCCESS;
}

enum TaskStatus EMFCorrectReceive(MeshBlock *pmb, int step)
{
  if(pmb->pbval->ReceiveEMFCorrection() == true) {
    return TASK_NEXT;
  } else {
    return TASK_FAIL;
  }
}

enum TaskStatus HydroIntegrate(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;

  if(step == 1) {
    phydro->FluxDivergence(pmb, phydro->u1, phydro->w, pfield->b,
                                          pfield->bcc, 1);
    // source terms
    Real dt = 0.5*(pmb->pmy_mesh->dt);
    phydro->psrc->AddHydroSourceTerms(dt,pmb->phydro->flux,phydro->w,pfield->bcc,phydro->u);
  } else if(step == 2) {
    phydro->FluxDivergence(pmb, phydro->u, phydro->w1, pfield->b1,
                                          pfield->bcc1, 2);
    // source terms
    Real dt = (pmb->pmy_mesh->dt);
    phydro->psrc->AddHydroSourceTerms(dt,pmb->phydro->flux,phydro->w,pfield->bcc,phydro->u);
  } else {
    return TASK_FAIL;
  }

  return TASK_NEXT;
}

enum TaskStatus HydroSend(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  BoundaryValues *pbval=pmb->pbval;
  if(step == 1) {
    pbval->SendHydroBoundaryBuffers(phydro->u1, 1);
  } else if(step == 2) {
    pbval->SendHydroBoundaryBuffers(phydro->u, 0);
  } else {
    return TASK_FAIL;
  }
  return TASK_SUCCESS;
}

enum TaskStatus HydroReceive(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  BoundaryValues *pbval=pmb->pbval;
  bool ret;
  if(step == 1) {
    ret=pbval->ReceiveHydroBoundaryBuffers(phydro->u1,1);
  } else if(step == 2) {
    ret=pbval->ReceiveHydroBoundaryBuffers(phydro->u,0);
  } else {
    return TASK_FAIL;
  }
  if(ret==true) {
    return TASK_SUCCESS;
  } else {
    return TASK_FAIL;
  }
}

enum TaskStatus FieldIntegrate(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  if(step == 1) {
//    pfield->b1.x1f = pfield->b.x1f;
//    pfield->b1.x2f = pfield->b.x2f;
//    pfield->b1.x3f = pfield->b.x3f;
    pfield->CopyOrAverageField(pfield->b1, pfield->b, pfield->b, 0.0);

    pfield->CT(pmb, pfield->b1, phydro->w, pfield->bcc, 1);
  } else if(step == 2) {
    pfield->CT(pmb, pfield->b, phydro->w1, pfield->bcc1, 2);
  } else {
    return TASK_FAIL;
  }
  return TASK_NEXT;
}

enum TaskStatus FieldSend(MeshBlock *pmb, int step)
{
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;
  if(step == 1) {
    pbval->SendFieldBoundaryBuffers(pfield->b1,1);
  } else if(step == 2) {
    pbval->SendFieldBoundaryBuffers(pfield->b,0);
  } else {
    return TASK_FAIL;
  }
  return TASK_SUCCESS;
}

enum TaskStatus FieldReceive(MeshBlock *pmb, int step)
{
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;
  bool ret;
  if(step == 1) {
    ret=pbval->ReceiveFieldBoundaryBuffers(pfield->b1,1);
  } else if(step == 2) {
    ret=pbval->ReceiveFieldBoundaryBuffers(pfield->b,0);
  } else {
    return TASK_FAIL;
  }
  if(ret==true) {
    return TASK_SUCCESS;
  } else {
    return TASK_FAIL;
  }
}

enum TaskStatus Prolongation(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;
  if(step == 1) {
    pbval->ProlongateBoundaries(phydro->w1, phydro->u1, pfield->b1, pfield->bcc1);
  } else if(step == 2) {
    pbval->ProlongateBoundaries(phydro->w,  phydro->u,  pfield->b,  pfield->bcc);
  } else {
    return TASK_FAIL;
  }
  return TASK_SUCCESS;
}

enum TaskStatus Primitives(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  if(pmb->nblevel[1][1][0]!=-1) is-=NGHOST;
  if(pmb->nblevel[1][1][2]!=-1) ie+=NGHOST;
  if(pmb->nblevel[1][0][1]!=-1) js-=NGHOST;
  if(pmb->nblevel[1][2][1]!=-1) je+=NGHOST;
  if(pmb->nblevel[0][1][1]!=-1) ks-=NGHOST;
  if(pmb->nblevel[2][1][1]!=-1) ke+=NGHOST;

  if(step == 1) {
    pmb->peos->ConservedToPrimitive(phydro->u1, phydro->w, pfield->b1,
                                    phydro->w1, pfield->bcc1, pmb->pcoord,
                                    is, ie, js, je, ks, ke);
  } else if(step == 2) {
    pmb->peos->ConservedToPrimitive(phydro->u, phydro->w1, pfield->b,
                                    phydro->w, pfield->bcc, pmb->pcoord,
                                    is, ie, js, je, ks, ke);
  } else {
    return TASK_FAIL;
  }
  return TASK_SUCCESS;
}

enum TaskStatus PhysicalBoundary(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;
  if(step == 1) {
    pbval->ApplyPhysicalBoundaries(phydro->w1, phydro->u1, pfield->b1, pfield->bcc1);
  } else if(step == 2) {
    pbval->ApplyPhysicalBoundaries(phydro->w,  phydro->u,  pfield->b,  pfield->bcc);
  } else {
    return TASK_FAIL;
  }
  return TASK_SUCCESS;
}

enum TaskStatus UserWork(MeshBlock *pmb, int step)
{
  if (step != 2) return TASK_SUCCESS;

  pmb->UserWorkInLoop();
  return TASK_SUCCESS;
}

enum TaskStatus NewBlockTimeStep(MeshBlock *pmb, int step)
{
  if (step != 2) return TASK_SUCCESS;

  pmb->phydro->NewBlockTimeStep(pmb);
  return TASK_SUCCESS;
}

enum TaskStatus CheckRefinement(MeshBlock *pmb, int step)
{
  if (step != 2) return TASK_SUCCESS;

  pmb->pmr->CheckRefinementCondition();
  return TASK_SUCCESS;
}

} // namespace TaskFunctions

//--------------------------------------------------------------------------------------
//! \fn
//  \brief

void TaskList::CreateTimeIntegrator(Mesh *pm)
{
  if (MAGNETIC_FIELDS_ENABLED) { // MHD
    AddTask(CALC_FLX,NONE);
    if(pm->multilevel==true) { // SMR or AMR
      AddTask(FLX_SEND,CALC_FLX);
      AddTask(FLX_RECV,CALC_FLX);
      AddTask(HYD_INT, FLX_RECV);
      AddTask(HYD_SEND,HYD_INT);
      AddTask(CALC_EMF,CALC_FLX);
      AddTask(EMF_SEND,CALC_EMF);
      AddTask(EMF_RECV,EMF_SEND);
      AddTask(FLD_INT, EMF_RECV);
    } else {
      AddTask(HYD_INT, CALC_FLX);
      AddTask(HYD_SEND,HYD_INT);
      AddTask(CALC_EMF,CALC_FLX);
      AddTask(EMF_SEND,CALC_EMF);
      AddTask(EMF_RECV,EMF_SEND);
      AddTask(FLD_INT, EMF_RECV);
    }
    AddTask(FLD_SEND,FLD_INT);
    AddTask(HYD_RECV,NONE);
    AddTask(FLD_RECV,NONE);
    if(pm->multilevel==true) { // SMR or AMR
      AddTask(PROLONG, (HYD_SEND|HYD_RECV|FLD_SEND|FLD_RECV));
      AddTask(CON2PRIM,PROLONG);
    } else {
      AddTask(CON2PRIM,(HYD_INT|HYD_RECV|FLD_INT|FLD_RECV));
    }
    AddTask(PHY_BVAL,CON2PRIM);

  } else {  // HYDRO
    AddTask(CALC_FLX,NONE);
    if(pm->multilevel==true) { // SMR or AMR
      AddTask(FLX_SEND,CALC_FLX);
      AddTask(FLX_RECV,CALC_FLX);
      AddTask(HYD_INT, FLX_RECV);
    }
    else {
      AddTask(HYD_INT,CALC_FLX);
    }
    AddTask(HYD_SEND,HYD_INT);
    AddTask(HYD_RECV,NONE);
    if(pm->multilevel==true) { // SMR or AMR
      AddTask(PROLONG,(HYD_SEND|HYD_RECV));
      AddTask(CON2PRIM,PROLONG);
    } else {
      AddTask(CON2PRIM,(HYD_INT|HYD_RECV));
    }
    AddTask(PHY_BVAL,CON2PRIM);
  }

  AddTask(USERWORK,PHY_BVAL);
  AddTask(NEW_DT,USERWORK);
  if(pm->adaptive==true) {
    AddTask(AMR_FLAG,USERWORK);
  }
}

//--------------------------------------------------------------------------------------//! \fn
//  \brief

void TaskList::AddTask(unsigned long int id,unsigned long int dep)
{
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;

  switch((id))
  {
  case (CALC_FLX):
    task_list_[ntasks].TaskFunc=TaskFunctions::CalculateFluxes;
    break;

  case (FLX_SEND):
    task_list_[ntasks].TaskFunc=TaskFunctions::FluxCorrectSend;
    break;

  case (FLX_RECV):
    task_list_[ntasks].TaskFunc=TaskFunctions::FluxCorrectReceive;
    break;

  case (CALC_EMF):
    task_list_[ntasks].TaskFunc=TaskFunctions::CalculateEMF;
    break;

  case (EMF_SEND):
    task_list_[ntasks].TaskFunc=TaskFunctions::EMFCorrectSend;
    break;

  case (EMF_RECV):
    task_list_[ntasks].TaskFunc=TaskFunctions::EMFCorrectReceive;
    break;

  case (HYD_INT):
    task_list_[ntasks].TaskFunc=TaskFunctions::HydroIntegrate;
    break;

  case (HYD_SEND):
    task_list_[ntasks].TaskFunc=TaskFunctions::HydroSend;
    break;

  case (HYD_RECV):
    task_list_[ntasks].TaskFunc=TaskFunctions::HydroReceive;
    break;

  case (FLD_INT):
    task_list_[ntasks].TaskFunc=TaskFunctions::FieldIntegrate;
    break;

  case (FLD_SEND):
    task_list_[ntasks].TaskFunc=TaskFunctions::FieldSend;
    break;

  case (FLD_RECV):
    task_list_[ntasks].TaskFunc=TaskFunctions::FieldReceive;
    break;

  case (PROLONG):
    task_list_[ntasks].TaskFunc=TaskFunctions::Prolongation;
    break;

  case (PHY_BVAL):
    task_list_[ntasks].TaskFunc=TaskFunctions::PhysicalBoundary;
    break;

  case (CON2PRIM):
    task_list_[ntasks].TaskFunc=TaskFunctions::Primitives;
    break;

  case (USERWORK):
    task_list_[ntasks].TaskFunc=TaskFunctions::UserWork;
    break;

  case (NEW_DT):
    task_list_[ntasks].TaskFunc=TaskFunctions::NewBlockTimeStep;
    break;

  case (AMR_FLAG):
    task_list_[ntasks].TaskFunc=TaskFunctions::CheckRefinement;
    break;

  default:
    std::stringstream msg;
    msg << "### FATAL ERROR in AddTask" << std::endl
        << "Invalid Task "<< id << " is specified" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  ntasks++;
  return;
}

