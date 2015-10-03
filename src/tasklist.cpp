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
//! \file tasklist.hpp
//  \brief task functions
//======================================================================================

// C/C++ headers
#include <iostream>   // cout, endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ classes headers
#include "athena.hpp"
#include "mesh.hpp"
#include "hydro/hydro.hpp"
#include "field/field.hpp"
#include "bvals/bvals.hpp"
#include "hydro/eos/eos.hpp"
#include "hydro/integrators/hydro_integrator.hpp"
#include "field/integrators/field_integrator.hpp"

// this class header
#include "tasklist.hpp"

namespace taskfunc {

enum task_status HydroIntegrate(MeshBlock *pmb, int task_flag)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  if(task_flag==0) {
    phydro->pf_integrator->OneStep(pmb, phydro->u, phydro->w1, pfield->b1,
                                   pfield->bcc1, 2);
  }
  else if(task_flag==1) {
    phydro->u1 = phydro->u;
    phydro->pf_integrator->OneStep(pmb, phydro->u1, phydro->w, pfield->b,
                                   pfield->bcc, 1);
  }

  return task_do_next;
}

enum task_status CalculateEMF(MeshBlock *pmb, int task_flag)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  if(task_flag==0)
    pfield->pint->ComputeCornerE(pmb, phydro->w1, pfield->bcc1);
  else if(task_flag==1)
    pfield->pint->ComputeCornerE(pmb, phydro->w, pfield->bcc);
  return task_do_next;
}

enum task_status FieldIntegrate(MeshBlock *pmb, int task_flag)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  if(task_flag==0) {
    pfield->pint->CT(pmb, pfield->b, phydro->w1, pfield->bcc1, 2);
  }
  else if(task_flag==1) {
    pfield->b1.x1f = pfield->b.x1f;
    pfield->b1.x2f = pfield->b.x2f;
    pfield->b1.x3f = pfield->b.x3f;
    pfield->pint->CT(pmb, pfield->b1, phydro->w, pfield->bcc, 1);
  }
  return task_do_next;
}

enum task_status HydroSend(MeshBlock *pmb, int task_flag)
{
  Hydro *phydro=pmb->phydro;
  BoundaryValues *pbval=pmb->pbval;
  if(task_flag==0)
    pbval->SendHydroBoundaryBuffers(phydro->u,0);
  else if(task_flag==1)
    pbval->SendHydroBoundaryBuffers(phydro->u1,1);
  return task_success;
}

enum task_status HydroReceive(MeshBlock *pmb, int task_flag)
{
  Hydro *phydro=pmb->phydro;
  BoundaryValues *pbval=pmb->pbval;
  bool ret;
  if(task_flag==0)
    ret=pbval->ReceiveHydroBoundaryBuffers(phydro->u,0);
  else if(task_flag==1)
    ret=pbval->ReceiveHydroBoundaryBuffers(phydro->u1,1);
  if(ret==true)
    return task_success;
  return task_failure;
}

enum task_status FluxCorrectionSend(MeshBlock *pmb, int task_flag)
{
  pmb->pbval->SendFluxCorrection(task_flag);
  return task_success;
}

enum task_status FluxCorrectionReceive(MeshBlock *pmb, int task_flag)
{
  Hydro *phydro=pmb->phydro;
  BoundaryValues *pbval=pmb->pbval;
  bool ret;
  if(task_flag==0)
    ret=pbval->ReceiveFluxCorrection(phydro->u,0);
  else if(task_flag==1)
    ret=pbval->ReceiveFluxCorrection(phydro->u1,1);
  if(ret==true) return task_do_next;
  return task_failure;
}

enum task_status HydroProlongation(MeshBlock *pmb, int task_flag)
{
  Hydro *phydro=pmb->phydro;
  BoundaryValues *pbval=pmb->pbval;
  if(task_flag==0)
    pbval->ProlongateHydroBoundaries(phydro->u);
  else if(task_flag==1)
    pbval->ProlongateHydroBoundaries(phydro->u1);
  return task_success;
}

enum task_status HydroPhysicalBoundary(MeshBlock *pmb, int task_flag)
{
  Hydro *phydro=pmb->phydro;
  BoundaryValues *pbval=pmb->pbval;
  if(task_flag==0)
    pbval->HydroPhysicalBoundaries(phydro->u);
  else if(task_flag==1)
    pbval->HydroPhysicalBoundaries(phydro->u1);
  return task_success;
}

enum task_status FieldSend(MeshBlock *pmb, int task_flag)
{
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;
  if(task_flag==0)
    pbval->SendFieldBoundaryBuffers(pfield->b,0);
  else if(task_flag==1)
    pbval->SendFieldBoundaryBuffers(pfield->b1,1);
  return task_success;
}

enum task_status FieldReceive(MeshBlock *pmb, int task_flag)
{
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;
  bool ret;
  if(task_flag==0)
    ret=pbval->ReceiveFieldBoundaryBuffers(pfield->b,0);
  else if(task_flag==1)
    ret=pbval->ReceiveFieldBoundaryBuffers(pfield->b1,1);
  if(ret==true)
    return task_success;
  return task_failure;
}

enum task_status EMFCorrectionSend(MeshBlock *pmb, int task_flag)
{
  pmb->pbval->SendEMFCorrection(task_flag);
  return task_success;
}

enum task_status EMFCorrectionReceive(MeshBlock *pmb, int task_flag)
{
  BoundaryValues *pbval=pmb->pbval;
  if(pbval->ReceiveEMFCorrection(task_flag)==true) return task_do_next;
  else return task_failure;
}

enum task_status FieldProlongation(MeshBlock *pmb, int task_flag)
{
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;
  if(task_flag==0)
    pbval->ProlongateFieldBoundaries(pfield->b);
  else if(task_flag==1)
    pbval->ProlongateFieldBoundaries(pfield->b1);
  return task_success;
}

enum task_status FieldPhysicalBoundary(MeshBlock *pmb, int task_flag)
{
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;
  if(task_flag==0)
    pbval->FieldPhysicalBoundaries(pfield->b);
  else if(task_flag==1)
    pbval->FieldPhysicalBoundaries(pfield->b1);
  return task_success;
}

enum task_status Primitives(MeshBlock *pmb, int task_flag)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  if(task_flag==0)
    phydro->pf_eos->ConservedToPrimitive(phydro->u, phydro->w1, pfield->b,
                                         phydro->w, pfield->bcc);
  else if(task_flag==1)
    phydro->pf_eos->ConservedToPrimitive(phydro->u1, phydro->w, pfield->b1,
                                         phydro->w1, pfield->bcc1);

  return task_success;
}

enum task_status NewBlockTimeStep(MeshBlock *pmb, int task_flag)
{
  pmb->phydro->NewBlockTimeStep(pmb);
  return task_success;
}

} // namespace task


void TaskList::AddTask(enum task_code t, unsigned long int dependency)
{
  std::stringstream msg;
  task_list_[ntasks_].task_id=t;
  task_list_[ntasks_].dependency=dependency;
  switch(t)
  {
  case hydro_integrate_1:
    task_list_[ntasks_].TaskFunc=taskfunc::HydroIntegrate;
    task_list_[ntasks_].task_flag=1;
    break;

  case calculate_emf_1:
    task_list_[ntasks_].TaskFunc=taskfunc::CalculateEMF;
    task_list_[ntasks_].task_flag=1;
    break;

  case field_integrate_1:
    task_list_[ntasks_].TaskFunc=taskfunc::FieldIntegrate;
    task_list_[ntasks_].task_flag=1;
    break;

  case hydro_send_1:
    task_list_[ntasks_].TaskFunc=taskfunc::HydroSend;
    task_list_[ntasks_].task_flag=1;
    break;

  case hydro_recv_1:
    task_list_[ntasks_].TaskFunc=taskfunc::HydroReceive;
    task_list_[ntasks_].task_flag=1;
    break;

  case flux_correct_send_1:
    task_list_[ntasks_].TaskFunc=taskfunc::FluxCorrectionSend;
    task_list_[ntasks_].task_flag=1;
    break;

  case flux_correct_recv_1:
    task_list_[ntasks_].TaskFunc=taskfunc::FluxCorrectionReceive;
    task_list_[ntasks_].task_flag=1;
    break;

  case hydro_prolong_1:
    task_list_[ntasks_].TaskFunc=taskfunc::HydroProlongation;
    task_list_[ntasks_].task_flag=1;
    break;

  case hydro_boundary_1:
    task_list_[ntasks_].TaskFunc=taskfunc::HydroPhysicalBoundary;
    task_list_[ntasks_].task_flag=1;
    break;

  case field_send_1:
    task_list_[ntasks_].TaskFunc=taskfunc::FieldSend;
    task_list_[ntasks_].task_flag=1;
    break;

  case field_recv_1:
    task_list_[ntasks_].TaskFunc=taskfunc::FieldReceive;
    task_list_[ntasks_].task_flag=1;
    break;

  case emf_correct_send_1:
    task_list_[ntasks_].TaskFunc=taskfunc::EMFCorrectionSend;
    task_list_[ntasks_].task_flag=1;
    break;

  case emf_correct_recv_1:
    task_list_[ntasks_].TaskFunc=taskfunc::EMFCorrectionReceive;
    task_list_[ntasks_].task_flag=1;
    break;

  case field_prolong_1:
    task_list_[ntasks_].TaskFunc=taskfunc::FieldProlongation;
    task_list_[ntasks_].task_flag=1;
    break;

  case field_boundary_1:
    task_list_[ntasks_].TaskFunc=taskfunc::FieldPhysicalBoundary;
    task_list_[ntasks_].task_flag=1;
    break;

  case primitives_1:
    task_list_[ntasks_].TaskFunc=taskfunc::Primitives;
    task_list_[ntasks_].task_flag=1;
    break;

  case hydro_integrate_0:
    task_list_[ntasks_].TaskFunc=taskfunc::HydroIntegrate;
    task_list_[ntasks_].task_flag=0;
    break;

  case calculate_emf_0:
    task_list_[ntasks_].TaskFunc=taskfunc::CalculateEMF;
    task_list_[ntasks_].task_flag=0;
    break;

  case field_integrate_0:
    task_list_[ntasks_].TaskFunc=taskfunc::FieldIntegrate;
    task_list_[ntasks_].task_flag=0;
    break;

  case hydro_send_0:
    task_list_[ntasks_].TaskFunc=taskfunc::HydroSend;
    task_list_[ntasks_].task_flag=0;
    break;

  case hydro_recv_0:
    task_list_[ntasks_].TaskFunc=taskfunc::HydroReceive;
    task_list_[ntasks_].task_flag=0;
    break;

  case flux_correct_send_0:
    task_list_[ntasks_].TaskFunc=taskfunc::FluxCorrectionSend;
    task_list_[ntasks_].task_flag=0;
    break;

  case flux_correct_recv_0:
    task_list_[ntasks_].TaskFunc=taskfunc::FluxCorrectionReceive;
    task_list_[ntasks_].task_flag=0;
    break;

  case hydro_prolong_0:
    task_list_[ntasks_].TaskFunc=taskfunc::HydroProlongation;
    task_list_[ntasks_].task_flag=0;
    break;

  case hydro_boundary_0:
    task_list_[ntasks_].TaskFunc=taskfunc::HydroPhysicalBoundary;
    task_list_[ntasks_].task_flag=0;
    break;

  case field_send_0:
    task_list_[ntasks_].TaskFunc=taskfunc::FieldSend;
    task_list_[ntasks_].task_flag=0;
    break;

  case field_recv_0:
    task_list_[ntasks_].TaskFunc=taskfunc::FieldReceive;
    task_list_[ntasks_].task_flag=0;
    break;

  case emf_correct_send_0:
    task_list_[ntasks_].TaskFunc=taskfunc::EMFCorrectionSend;
    task_list_[ntasks_].task_flag=0;
    break;

  case emf_correct_recv_0:
    task_list_[ntasks_].TaskFunc=taskfunc::EMFCorrectionReceive;
    task_list_[ntasks_].task_flag=0;
    break;

  case field_prolong_0:
    task_list_[ntasks_].TaskFunc=taskfunc::FieldProlongation;
    task_list_[ntasks_].task_flag=0;
    break;

  case field_boundary_0:
    task_list_[ntasks_].TaskFunc=taskfunc::FieldPhysicalBoundary;
    task_list_[ntasks_].task_flag=0;
    break;

  case primitives_0:
    task_list_[ntasks_].TaskFunc=taskfunc::Primitives;
    task_list_[ntasks_].task_flag=0;
    break;

  case new_blocktimestep:
    task_list_[ntasks_].TaskFunc=taskfunc::NewBlockTimeStep;
    task_list_[ntasks_].task_flag=0;
    break;

  default:
    msg << "### FATAL ERROR in AddTask" << std::endl
        << "Invalid Task "<< t << " is specified" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  
  }
  ntasks_++;
  return;
}

/*
void TaskList::CreateTaskList()
{
  if (MAGNETIC_FIELDS_ENABLED) { // MHD
    task_list.AddTask(hydro_integrate_1, none); // predict hydro
    if(pmesh->multilevel==true) { // SMR or AMR
      task_list.AddTask(flux_correct_send_1, hydro_integrate_1);
      task_list.AddTask(flux_correct_recv_1, hydro_integrate_1);
      task_list.AddTask(hydro_send_1, flux_correct_recv_1); // send block boundaries
      task_list.AddTask(calculate_emf_1, hydro_integrate_1); // calculate EMF
      task_list.AddTask(emf_correct_send_1, calculate_emf_1);
      task_list.AddTask(emf_correct_recv_1, emf_correct_send_1);
      task_list.AddTask(field_integrate_1, emf_correct_recv_1); // predict b-field
    }
    else {
      task_list.AddTask(hydro_send_1, hydro_integrate_1); // send block boundaries
      task_list.AddTask(calculate_emf_1, hydro_integrate_1); // calculate EMF
      task_list.AddTask(field_integrate_1, calculate_emf_1); // predict b-field
    }
    task_list.AddTask(field_send_1, field_integrate_1); // send block boundaries
    task_list.AddTask(hydro_recv_1, none); // receive block boundaries
    task_list.AddTask(hydro_boundary_1, hydro_recv_1 | hydro_integrate_1); // physical boundaries
    if(pmesh->multilevel==true) // SMR or AMR
      task_list.AddTask(hydro_prolong_1, hydro_boundary_1);
    task_list.AddTask(field_recv_1, none); // receive block boundaries
    task_list.AddTask(field_boundary_1, field_recv_1 | field_integrate_1); // physical boundaries
    if(pmesh->multilevel==true) {// SMR or AMR
      task_list.AddTask(field_prolong_1, field_boundary_1);
      task_list.AddTask(primitives_1, hydro_prolong_1 | field_prolong_1);
    }
    else 
      task_list.AddTask(primitives_1, hydro_boundary_1 | field_boundary_1);

    task_list.AddTask(hydro_integrate_0, primitives_1); // predict hydro
    if(pmesh->multilevel==true) { // SMR or AMR
      task_list.AddTask(flux_correct_send_0, hydro_integrate_0);
      task_list.AddTask(flux_correct_recv_0, hydro_integrate_0);
      task_list.AddTask(hydro_send_0, flux_correct_recv_0); // send block boundaries
      task_list.AddTask(calculate_emf_0, hydro_integrate_0); // calculate EMF
      task_list.AddTask(emf_correct_send_0, calculate_emf_0);
      task_list.AddTask(emf_correct_recv_0, emf_correct_send_0);
      task_list.AddTask(field_integrate_0, emf_correct_recv_0); // predict b-field
    }
    else {
      task_list.AddTask(hydro_send_0, hydro_integrate_0); // send block boundaries
      task_list.AddTask(calculate_emf_0, hydro_integrate_0); // calculate EMF
      task_list.AddTask(field_integrate_0, calculate_emf_0); // predict b-field
    }
    task_list.AddTask(field_send_0, field_integrate_0); // send block boundaries
    task_list.AddTask(hydro_recv_0, primitives_1); // receive block boundaries
    task_list.AddTask(hydro_boundary_0, hydro_recv_0 | hydro_integrate_0); // physical boundaries
    if(pmesh->multilevel==true) // SMR or AMR
      task_list.AddTask(hydro_prolong_0, hydro_boundary_0);
    task_list.AddTask(field_recv_0, primitives_1); // receive block boundaries
    task_list.AddTask(field_boundary_0, field_recv_0 | field_integrate_0); // physical boundaries
    if(pmesh->multilevel==true) {// SMR or AMR
      task_list.AddTask(field_prolong_0, field_boundary_0);
      task_list.AddTask(primitives_0, hydro_prolong_0 | field_prolong_0);
    }
    else 
      task_list.AddTask(primitives_0, hydro_boundary_0 | field_boundary_0);
  }
  else { // hydro
    task_list.AddTask(hydro_integrate_1, none); // predict hydro
    if(pmesh->multilevel==true) { // SMR or AMR
      task_list.AddTask(flux_correct_send_1, hydro_integrate_1);
      task_list.AddTask(flux_correct_recv_1, hydro_integrate_1);
      task_list.AddTask(hydro_send_1, flux_correct_recv_1); // send block boundaries
    }
    else
      task_list.AddTask(hydro_send_1, hydro_integrate_1); // send block boundaries
    task_list.AddTask(hydro_recv_1, none); // receive block boundaries
    task_list.AddTask(hydro_boundary_1, hydro_recv_1 | hydro_integrate_1); // physical boundaries
    if(pmesh->multilevel==true) { // SMR or AMR
      task_list.AddTask(hydro_prolong_1, hydro_boundary_1);
      task_list.AddTask(primitives_1, hydro_prolong_1);
    }
    else
      task_list.AddTask(primitives_1, hydro_boundary_1);
    task_list.AddTask(hydro_integrate_0, primitives_1); // predict correct
    if(pmesh->multilevel==true) { // SMR or AMR
      task_list.AddTask(flux_correct_send_0, hydro_integrate_0);
      task_list.AddTask(flux_correct_recv_0, hydro_integrate_0);
      task_list.AddTask(hydro_send_0, flux_correct_recv_0); // send block boundaries
    }
    else
      task_list.AddTask(hydro_send_0, hydro_integrate_0); // send block boundaries
    task_list.AddTask(hydro_recv_0, primitives_1); // receive block boundaries
    task_list.AddTask(hydro_boundary_0, hydro_recv_0 | hydro_integrate_0); // physical boundaries
    if(pmesh->multilevel==true) { // SMR or AMR
      task_list.AddTask(hydro_prolong_0, hydro_boundary_0);
      task_list.AddTask(primitives_0, hydro_prolong_0);
    }
    else
      task_list.AddTask(primitives_0, hydro_boundary_0);
  }
  task_list.AddTask(new_blocktimestep,primitives_0);

  pmesh->SetTaskList(task_list);
}
*/
