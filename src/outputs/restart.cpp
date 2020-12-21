//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file restart.cpp
//! \brief writes restart files

// C headers

// C++ headers
#include <cstdio>    // snprintf()
#include <cstring>   // memcpy()
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"
#include "outputs.hpp"


//----------------------------------------------------------------------------------------
//! \fn void RestartOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag)
//! \brief Cycles over all MeshBlocks and writes data to a single restart file.

void RestartOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin, bool force_write) {
  IOWrapper resfile;
  IOWrapperSizeT listsize, headeroffset, datasize;

  // create single output filename:"file_basename"+"."+XXXXX+".rst",
  // where XXXXX = 5-digit file_number
  std::string fname;
  char number[6];
  std::snprintf(number, sizeof(number), "%05d", output_params.file_number);

  fname.assign(output_params.file_basename);
  fname.append(".");
  // add file number to name, unless write is forced by terminate signal, in which case
  // replace number in the name by the string "final".  This keeps the restart file
  // numbers consistent with output.dt when a job is restarted many times.
  if (!force_write)
    fname.append(number);
  else
    fname.append("final");
  fname.append(".rst");

  // increment counters now so values for *next* dump are stored in restart file
  if (!force_write) {
    output_params.file_number++;
    output_params.next_time += output_params.dt;
    pin->SetInteger(output_params.block_name, "file_number", output_params.file_number);
    pin->SetReal(output_params.block_name, "next_time", output_params.next_time);
  }
  resfile.Open(fname.c_str(), IOWrapper::FileMode::write);

  // prepare the input parameters
  std::stringstream ost;
  pin->ParameterDump(ost);
  std::string sbuf = ost.str();

  // calculate the header size
  IOWrapperSizeT udsize = 0;
  for (int n=0; n<pm->nint_user_mesh_data_; n++)
    udsize += pm->iuser_mesh_data[n].GetSizeInBytes();
  for (int n=0; n<pm->nreal_user_mesh_data_; n++)
    udsize += pm->ruser_mesh_data[n].GetSizeInBytes();

  headeroffset = sbuf.size()*sizeof(char) + 3*sizeof(int)+sizeof(RegionSize)
                 + 2*sizeof(Real)+sizeof(IOWrapperSizeT)+udsize;
  // the size of an element of the ID and cost list
  listsize = sizeof(LogicalLocation)+sizeof(double);
  // the size of each MeshBlock
  datasize = pm->my_blocks(0)->GetBlockSizeInBytes();
  int nbtotal = pm->nbtotal;
  int myns = pm->nslist[Globals::my_rank];
  int mynb = pm->nblist[Globals::my_rank];

  // write the header; this part is serial
  if (Globals::my_rank == 0) {
    // output the input parameters
    resfile.Write(sbuf.c_str(),sizeof(char),sbuf.size());

    // output Mesh information
    resfile.Write(&(pm->nbtotal), sizeof(int), 1);
    resfile.Write(&(pm->root_level), sizeof(int), 1);
    resfile.Write(&(pm->mesh_size), sizeof(RegionSize), 1);
    resfile.Write(&(pm->time), sizeof(Real), 1);
    resfile.Write(&(pm->dt), sizeof(Real), 1);
    resfile.Write(&(pm->ncycle), sizeof(int), 1);
    resfile.Write(&(datasize), sizeof(IOWrapperSizeT), 1);

    // collect and write user Mesh data
    if (udsize != 0) {
      char *ud = new char[udsize];
      IOWrapperSizeT udoffset = 0;
      for (int n=0; n<pm->nint_user_mesh_data_; n++) {
        std::memcpy(&(ud[udoffset]), pm->iuser_mesh_data[n].data(),
                    pm->iuser_mesh_data[n].GetSizeInBytes());
        udoffset += pm->iuser_mesh_data[n].GetSizeInBytes();
      }
      for (int n=0; n<pm->nreal_user_mesh_data_; n++) {
        std::memcpy(&(ud[udoffset]), pm->ruser_mesh_data[n].data(),
                    pm->ruser_mesh_data[n].GetSizeInBytes());
        udoffset += pm->ruser_mesh_data[n].GetSizeInBytes();
      }
      resfile.Write(ud, 1, udsize);
      delete [] ud;
    }
  }

  // allocate memory for the ID list and the data
  char *idlist = new char[listsize*mynb];
  char *data = new char[mynb*datasize];

  // Loop over MeshBlocks and pack the meta data
  int os=0;
  for (int b=0; b<pm->nblocal; ++b) {
    MeshBlock *pmb = pm->my_blocks(b);
    std::memcpy(&(idlist[os]), &(pmb->loc), sizeof(LogicalLocation));
    os += sizeof(LogicalLocation);
    std::memcpy(&(idlist[os]), &(pmb->cost_), sizeof(double));
    os += sizeof(double);
  }

  // write the ID list collectively
  IOWrapperSizeT myoffset = headeroffset + listsize*myns;
  resfile.Write_at_all(idlist, listsize, mynb, myoffset);

  // deallocate the idlist array
  delete [] idlist;

  // Loop over MeshBlocks and pack the data
  for (int b=0; b<pm->nblocal; ++b) {
    MeshBlock *pmb = pm->my_blocks(b);
    char *pdata = &(data[pmb->lid*datasize]);

    // NEW_OUTPUT_TYPES: add output of additional physics to restarts here also update
    // MeshBlock::GetBlockSizeInBytes accordingly and MeshBlock constructor for restarts.

    // Hydro conserved variables:
    std::memcpy(pdata, pmb->phydro->u.data(), pmb->phydro->u.GetSizeInBytes());
    pdata += pmb->phydro->u.GetSizeInBytes();

    // Hydro primitive variables (at current and previous step):
    if (GENERAL_RELATIVITY) {
      std::memcpy(pdata, pmb->phydro->w.data(), pmb->phydro->w.GetSizeInBytes());
      pdata += pmb->phydro->w.GetSizeInBytes();
      std::memcpy(pdata, pmb->phydro->w1.data(), pmb->phydro->w1.GetSizeInBytes());
      pdata += pmb->phydro->w1.GetSizeInBytes();
    }

    // Longitudinal, face-centered magnetic field components:
    if (MAGNETIC_FIELDS_ENABLED) {
      std::memcpy(pdata, pmb->pfield->b.x1f.data(), pmb->pfield->b.x1f.GetSizeInBytes());
      pdata += pmb->pfield->b.x1f.GetSizeInBytes();
      std::memcpy(pdata, pmb->pfield->b.x2f.data(), pmb->pfield->b.x2f.GetSizeInBytes());
      pdata += pmb->pfield->b.x2f.GetSizeInBytes();
      std::memcpy(pdata, pmb->pfield->b.x3f.data(), pmb->pfield->b.x3f.GetSizeInBytes());
      pdata += pmb->pfield->b.x3f.GetSizeInBytes();
    }

    // (conserved variable) Passive scalars:
    if (NSCALARS > 0) {
      AthenaArray<Real> &s = pmb->pscalars->s;
      std::memcpy(pdata, s.data(), s.GetSizeInBytes());
      pdata += s.GetSizeInBytes();
    }
    // (primitive variable) density-normalized passive scalar concentrations
    // if ???
    // for (int n=0; n<NSCALARS; n++) {
    //   AthenaArray<Real> &r = pmb->pscalars->r;
    //   std::memcpy(pdata, r.data(), r.GetSizeInBytes());
    //   pdata += r.GetSizeInBytes();
    // }

    // User MeshBlock data:
    // integer data:
    for (int n=0; n<pmb->nint_user_meshblock_data_; n++) {
      std::memcpy(pdata, pmb->iuser_meshblock_data[n].data(),
                  pmb->iuser_meshblock_data[n].GetSizeInBytes());
      pdata += pmb->iuser_meshblock_data[n].GetSizeInBytes();
    }
    // floating-point data:
    for (int n=0; n<pmb->nreal_user_meshblock_data_; n++) {
      std::memcpy(pdata, pmb->ruser_meshblock_data[n].data(),
                  pmb->ruser_meshblock_data[n].GetSizeInBytes());
      pdata += pmb->ruser_meshblock_data[n].GetSizeInBytes();
    }
  }

  // now write restart data in parallel
  myoffset = headeroffset + listsize*nbtotal + datasize*myns;
  resfile.Write_at_all(data, datasize, mynb, myoffset);
  resfile.Close();
  delete [] data;
}
