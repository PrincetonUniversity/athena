//======================================================================================
//! \file tasklist.hpp
//  \brief task functions
//======================================================================================
#include <iostream>   // cout, endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

#include "athena.hpp"
#include "tasklist.hpp"
#include "mesh.hpp"
#include "fluid/fluid.hpp"
#include "field/field.hpp"
#include "bvals/bvals.hpp"
#include "fluid/eos/eos.hpp"
#include "fluid/integrators/fluid_integrator.hpp"
#include "field/integrators/field_integrator.hpp"

namespace taskfunc {

enum task_status FluidIntegrate(MeshBlock *pmb, int task_arg)
{
  Fluid *pfluid=pmb->pfluid;
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;
  if(task_arg==0) {
    pfluid->pf_integrator->OneStep(pmb, pfluid->u, pfluid->w1, pfield->b1,
                                   pfield->bcc1, 2);
  }
  else if(task_arg==1) {
    pfluid->u1 = pfluid->u;
    pfluid->pf_integrator->OneStep(pmb, pfluid->u1, pfluid->w, pfield->b,
                                   pfield->bcc, 1);
  }

  return task_donext;
}

enum task_status FieldIntegrate(MeshBlock *pmb, int task_arg)
{
  Fluid *pfluid=pmb->pfluid;
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;
  if(task_arg==0) {
    pfield->pint->CT(pmb, pfield->b, pfluid->w1, pfield->bcc1, 2);
  }
  else if(task_arg==1) {
    pfield->b1.x1f = pfield->b.x1f;
    pfield->b1.x2f = pfield->b.x2f;
    pfield->b1.x3f = pfield->b.x3f;
    pfield->pint->CT(pmb, pfield->b1, pfluid->w, pfield->bcc, 1);
  }
  return task_donext;
}

enum task_status FluidSend(MeshBlock *pmb, int task_arg)
{
  Fluid *pfluid=pmb->pfluid;
  BoundaryValues *pbval=pmb->pbval;
  if(task_arg==0)
    pbval->SendFluidBoundaryBuffers(pfluid->u,0);
  else if(task_arg==1)
    pbval->SendFluidBoundaryBuffers(pfluid->u1,1);
  return task_success;
}

enum task_status FluidReceive(MeshBlock *pmb, int task_arg)
{
  Fluid *pfluid=pmb->pfluid;
  BoundaryValues *pbval=pmb->pbval;
  bool ret;
  if(task_arg==0)
    ret=pbval->ReceiveFluidBoundaryBuffers(pfluid->u,0);
  else if(task_arg==1)
    ret=pbval->ReceiveFluidBoundaryBuffers(pfluid->u1,1);
  if(ret==true)
    return task_success;
  return task_failure;
}

enum task_status FluidPhysicalBoundary(MeshBlock *pmb, int task_arg)
{
  if(task_arg==0)
    FluidPhysicalBoundaries(pfluid->u);
  else if(task_arg==1)
    FluidPhysicalBoundaries(pfluid->u1);
  return task_success;
}

enum task_status FieldPhysicalBoundary(MeshBlock *pmb, int task_arg)
{
  if(task_arg==0)
    FieldPhysicalBoundaries(pfluid->b);
  else if(task_arg==1)
    FieldPhysicalBoundaries(pfluid->b1);
  return task_success;
}

enum task_status FluxCorrectionSend(MeshBlock *pmb, int task_arg)
{
  return task_success;
}

enum task_status FluxCorrectionRecv(MeshBlock *pmb, int task_arg)
{
  return task_success;
}

enum task_status FluidProlongation(MeshBlock *pmb, int task_arg)
{
  return task_success;
}

enum task_status FieldSend(MeshBlock *pmb, int task_arg)
{
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;
  if(task_arg==0)
    pbval->SendFieldBoundaryBuffers(pfield->b,0);
  else if(task_arg==1)
    pbval->SendFieldBoundaryBuffers(pfield->b1,1);
  return task_success;
}

enum task_status FieldReceive(MeshBlock *pmb, int task_arg)
{
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;
  bool ret;
  if(task_arg==0)
    ret=pbval->ReceiveFieldBoundaryBuffers(pfield->b,0);
  else if(task_arg==1)
    ret=pbval->ReceiveFluidBoundaryBuffers(pfield->b1,1);
  if(ret==true)
    return task_success;
  return task_failure;
}

enum task_status EMFCorrectionSend(MeshBlock *pmb, int task_arg)
{
  return task_success;
}

enum task_status EMFCorrectionRecv(MeshBlock *pmb, int task_arg)
{
  return task_success;
}

enum task_status FieldProlongation(MeshBlock *pmb, int task_arg)
{
  return task_success;
}

enum task_status Primitives(MeshBlock *pmb, int task_arg)
{
  Fluid *pfluid=pmb->pfluid;
  Field *pfield=pmb->pfield;
  if(task_arg==0)
    pfluid->pf_eos->ConservedToPrimitive(pfluid->u, pfluid->w1, pfield->b,
                                         pfluid->w, pfield->bcc);
  else if(task_arg==1)
    pfluid->pf_eos->ConservedToPrimitive(pfluid->u1, pfluid->w, pfield->b1,
                                         pfluid->w1, pfield->bcc1);

  return task_success;
}

enum task_status NewBlockTimeStep(MeshBlock *pmb, int task_arg)
{
  pmb->pfluid->NewBlockTimeStep(pmb);
  return true;
}

} // namespace task


void TaskList::AddTask(enum task t, unsigned long int dependence)
{
  std::stringstream msg;
  task[ntask].taskid=t;
  task[ntask].depend=dependence;
  switch(t)
  {
  case fluid_integrate_1:
    task[ntask].TaskFunc=taskfunc::FluidIntegrate;
    task[ntask].task_arg=1;
    break;

  case field_integrate_1:
    task[ntask].TaskFunc=taskfunc::FieldIntegrate;
    task[ntask].task_arg=1;
    break;

  case fluid_send_1:
    task[ntask].TaskFunc=taskfunc::FluidSend;
    task[ntask].task_arg=1;
    break;

  case fluid_recv_1:
    task[ntask].TaskFunc=taskfunc::FluidRecv;
    task[ntask].task_arg=1;
    break;

  case flux_correction_send_1:
    task[ntask].TaskFunc=taskfunc::FluxCorrectionSend;
    task[ntask].task_arg=1;
    break;

  case flux_correction_recv_1:
    task[ntask].TaskFunc=taskfunc::FluxCorrectionRecv;
    task[ntask].task_arg=1;
    break;

  case fluid_prolongation_1:
    task[ntask].TaskFunc=taskfunc::FluidProlongation;
    task[ntask].task_arg=1;
    break;

  case fluid_boundary_1:
    task[ntask].TaskFunc=taskfunc::FluidPhysicalBoundary;
    task[ntask].task_arg=1;
    break;

  case field_send_1:
    task[ntask].TaskFunc=taskfunc::FieldSend;
    task[ntask].task_arg=1;
    break;

  case field_recv_1:
    task[ntask].TaskFunc=taskfunc::FieldRecv;
    task[ntask].task_arg=1;
    break;

  case emf_correction_send_1:
    task[ntask].TaskFunc=taskfunc::EMFCorrectionSend;
    task[ntask].task_arg=1;
    break;

  case emf_correction_recv_1:
    task[ntask].TaskFunc=taskfunc::EMFCorrectionRecv;
    task[ntask].task_arg=1;
    break;

  case field_prolongation_1:
    task[ntask].TaskFunc=taskfunc::FieldProlongation;
    task[ntask].task_arg=1;
    break;

  case field_boundary_1:
    task[ntask].TaskFunc=taskfunc::FieldPhysicalBoundary;
    task[ntask].task_arg=1;
    break;

  case primitives_1:
    task[ntask].TaskFunc=taskfunc::Primitives;
    task[ntask].task_arg=1;
    break;

  case fluid_integrate_0:
    task[ntask].TaskFunc=taskfunc::FluidIntegrate;
    task[ntask].task_arg=0;
    break;

  case field_integrate_0:
    task[ntask].TaskFunc=taskfunc::FieldIntegrate;
    task[ntask].task_arg=0;
    break;

  case fluid_send_0:
    task[ntask].TaskFunc=taskfunc::FluidSend;
    task[ntask].task_arg=0;
    break;

  case fluid_recv_0:
    task[ntask].TaskFunc=taskfunc::FluidRecv;
    task[ntask].task_arg=0;
    break;

  case flux_correction_send_0:
    task[ntask].TaskFunc=taskfunc::FluxCorrectionSend;
    task[ntask].task_arg=0;
    break;

  case flux_correction_recv_0:
    task[ntask].TaskFunc=taskfunc::FluxCorrectionRecv;
    task[ntask].task_arg=0;
    break;

  case fluid_prolongation_0:
    task[ntask].TaskFunc=taskfunc::FluidProlongation;
    task[ntask].task_arg=0;
    break;

  case fluid_boundary_0:
    task[ntask].TaskFunc=taskfunc::FluidPhysicalBoundary;
    task[ntask].task_arg=0;
    break;

  case field_send_0:
    task[ntask].TaskFunc=taskfunc::FieldSend;
    task[ntask].task_arg=0;
    break;

  case field_recv_0:
    task[ntask].TaskFunc=taskfunc::FieldRecv;
    task[ntask].task_arg=0;
    break;

  case emf_correction_send_0:
    task[ntask].TaskFunc=taskfunc::EMFCorrectionSend;
    task[ntask].task_arg=0;
    break;

  case emf_correction_recv_0:
    task[ntask].TaskFunc=taskfunc::EMFCorrectionRecv;
    task[ntask].task_arg=0;
    break;

  case field_prolongation_0:
    task[ntask].TaskFunc=taskfunc::FieldProlongation;
    task[ntask].task_arg=0;
    break;

  case field_boundary_0:
    task[ntask].TaskFunc=taskfunc::FieldPhysicalBoundary;
    task[ntask].task_arg=0;
    break;

  case primitives_0:
    task[ntask].TaskFunc=taskfunc::Primitives;
    task[ntask].task_arg=0;
    break;

  case new_blocktimestep:
    task[ntask].TaskFunc=taskfunc::NewBlockTimeStep;
    task[ntask].task_arg=0;
    break;

  default:
    msg << "### FATAL ERROR in AddTask" << std::endl
        << "Invalid Task "<< t << " is specified" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  
  }
  ntask++;
  return;
}

