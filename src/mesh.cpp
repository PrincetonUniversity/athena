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
//! \file mesh.cpp
//  \brief implementation of functions in classes Mesh, and MeshBlock
//======================================================================================

// C/C++ headers
#include <cfloat>     // FLT_MAX
#include <cmath>      // std::abs(), pow()
#include <iostream>
#include <sstream>
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <algorithm>  // sort
#include <iomanip>
#include <stdlib.h>
#include <string.h>  // memcpy

// Athena++ classes headers
#include "athena.hpp"
#include "globals.hpp"
#include "athena_arrays.hpp"
#include "coordinates/coordinates.hpp"
#include "hydro/hydro.hpp" 
#include "field/field.hpp"
#include "bvals/bvals.hpp"
#include "hydro/eos/eos.hpp"
#include "hydro/integrators/hydro_integrator.hpp" 
#include "field/integrators/field_integrator.hpp"
#include "parameter_input.hpp"
#include "meshblocktree.hpp"
#include "outputs/wrapper.hpp"
#include "task_list.hpp"
#include "mesh_refinement/mesh_refinement.hpp"
#include "utils/buffer_utils.hpp"

// this class header
#include "mesh.hpp"

// MPI/OpenMP header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//--------------------------------------------------------------------------------------
// Mesh constructor, builds mesh at start of calculation using parameters in input file

Mesh::Mesh(ParameterInput *pin, int test_flag)
{
  std::stringstream msg;
  RegionSize block_size;
  MeshBlockTree *neibt;
  MeshBlock *pfirst;
  enum BoundaryFlag block_bcs[6];
  int nbmax, dim;

// mesh test
  if(test_flag>0) Globals::nranks=test_flag;

// read time and cycle limits from input file

  start_time = pin->GetOrAddReal("time","start_time",0.0);
  tlim       = pin->GetReal("time","tlim");
  cfl_number = pin->GetReal("time","cfl_number");
  time = start_time;
  dt   = (FLT_MAX*0.4);

  nlim = pin->GetOrAddInteger("time","nlim",-1);
  ncycle = 0;

// read number of OpenMP threads for mesh

  num_mesh_threads_ = pin->GetOrAddInteger("mesh","num_threads",1);
  if (num_mesh_threads_ < 1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Number of OpenMP threads must be >= 1, but num_threads=" 
        << num_mesh_threads_ << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

// read number of grid cells in root level of mesh from input file.  

  mesh_size.nx1 = pin->GetInteger("mesh","nx1");
  if (mesh_size.nx1 < 4) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "In mesh block in input file nx1 must be >= 4, but nx1=" 
        << mesh_size.nx1 << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  mesh_size.nx2 = pin->GetInteger("mesh","nx2");
  if (mesh_size.nx2 < 1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "In mesh block in input file nx2 must be >= 1, but nx2=" 
        << mesh_size.nx2 << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  mesh_size.nx3 = pin->GetInteger("mesh","nx3");
  if (mesh_size.nx3 < 1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "In mesh block in input file nx3 must be >= 1, but nx3=" 
        << mesh_size.nx3 << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (mesh_size.nx2 == 1 && mesh_size.nx3 > 1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "In mesh block in input file: nx2=1, nx3=" << mesh_size.nx3 
        << ", 2D problems in x1-x3 plane not supported" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

// check cfl_number
  if(cfl_number > 1.0 && mesh_size.nx2==1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "The CFL number must be smaller than 1.0 in 1D simulation" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if(cfl_number > 0.5 && mesh_size.nx2 > 1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "The CFL number must be smaller than 0.5 in 2D/3D simulation" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  dim=1;
  if(mesh_size.nx2>1) dim=2;
  if(mesh_size.nx3>1) dim=3;

// read physical size of mesh (root level) from input file.  

  mesh_size.x1min = pin->GetReal("mesh","x1min");
  mesh_size.x2min = pin->GetReal("mesh","x2min");
  mesh_size.x3min = pin->GetReal("mesh","x3min");

  mesh_size.x1max = pin->GetReal("mesh","x1max");
  mesh_size.x2max = pin->GetReal("mesh","x2max");
  mesh_size.x3max = pin->GetReal("mesh","x3max");

  if (mesh_size.x1max <= mesh_size.x1min) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Input x1max must be larger than x1min: x1min=" << mesh_size.x1min 
        << " x1max=" << mesh_size.x1max << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (mesh_size.x2max <= mesh_size.x2min) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Input x2max must be larger than x2min: x2min=" << mesh_size.x2min 
        << " x2max=" << mesh_size.x2max << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (mesh_size.x3max <= mesh_size.x3min) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Input x3max must be larger than x3min: x3min=" << mesh_size.x3min 
        << " x3max=" << mesh_size.x3max << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

// read ratios of grid cell size in each direction

  block_size.x1rat = mesh_size.x1rat = pin->GetOrAddReal("mesh","x1rat",1.0);
  block_size.x2rat = mesh_size.x2rat = pin->GetOrAddReal("mesh","x2rat",1.0);
  block_size.x3rat = mesh_size.x3rat = pin->GetOrAddReal("mesh","x3rat",1.0);

  if (std::abs(mesh_size.x1rat - 1.0) > 0.1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Ratio of cell sizes must be 0.9 <= x1rat <= 1.1, x1rat=" 
        << mesh_size.x1rat << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (std::abs(mesh_size.x2rat - 1.0) > 0.1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Ratio of cell sizes must be 0.9 <= x2rat <= 1.1, x2rat=" 
        << mesh_size.x2rat << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (std::abs(mesh_size.x3rat - 1.0) > 0.1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Ratio of cell sizes must be 0.9 <= x3rat <= 1.1, x3rat=" 
        << mesh_size.x3rat << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // read BC flags for each of the 6 boundaries in turn.
  mesh_bcs[INNER_X1] = GetBoundaryFlag(pin->GetOrAddString("mesh","ix1_bc","none"));
  mesh_bcs[OUTER_X1] = GetBoundaryFlag(pin->GetOrAddString("mesh","ox1_bc","none"));
  mesh_bcs[INNER_X2] = GetBoundaryFlag(pin->GetOrAddString("mesh","ix2_bc","none"));
  mesh_bcs[OUTER_X2] = GetBoundaryFlag(pin->GetOrAddString("mesh","ox2_bc","none"));
  mesh_bcs[INNER_X3] = GetBoundaryFlag(pin->GetOrAddString("mesh","ix3_bc","none"));
  mesh_bcs[OUTER_X3] = GetBoundaryFlag(pin->GetOrAddString("mesh","ox3_bc","none"));

// read MeshBlock parameters
  block_size.nx1 = pin->GetOrAddInteger("meshblock","nx1",mesh_size.nx1);
  if(dim>=2)
    block_size.nx2 = pin->GetOrAddInteger("meshblock","nx2",mesh_size.nx2);
  else
    block_size.nx2=mesh_size.nx2;
  if(dim==3)
    block_size.nx3 = pin->GetOrAddInteger("meshblock","nx3",mesh_size.nx3);
  else
    block_size.nx3=mesh_size.nx3;

// check consistency of the block and mesh
  if(mesh_size.nx1%block_size.nx1 != 0
  || mesh_size.nx2%block_size.nx2 != 0
  || mesh_size.nx3%block_size.nx3 != 0) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "the mesh must be evenly divisible by the meshblock" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if(block_size.nx1 <4 || (block_size.nx2<4 && dim>=2)
     || (block_size.nx3<4 && dim==3)) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "block_size must be larger than or equal to 4 meshes." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

// calculate the number of the blocks
  nrbx1=mesh_size.nx1/block_size.nx1;
  nrbx2=mesh_size.nx2/block_size.nx2;
  nrbx3=mesh_size.nx3/block_size.nx3;
  nbmax=(nrbx1>nrbx2)?nrbx1:nrbx2;
  nbmax=(nbmax>nrbx3)?nbmax:nrbx3;

  if(Globals::my_rank==0)
    std::cout << "RootGrid = " << nrbx1 << " x " << nrbx2
              << " x " << nrbx3 << std::endl;

// calculate the logical root level and maximum level
  for(root_level=0;(1<<root_level)<nbmax;root_level++);
  current_level=root_level;

// create the root grid
  tree.CreateRootGrid(nrbx1,nrbx2,nrbx3,root_level);

// SMR / AMR: create finer grids here
  multilevel=false;
  adaptive=false;
  if(pin->GetOrAddString("mesh","refinement","static")=="adaptive")
    adaptive=true, multilevel=true;
  if(adaptive==true) {
    max_level = pin->GetOrAddInteger("mesh","numlevel",1)+root_level-1;
    if(max_level > 63) {
      msg << "### FATAL ERROR in Mesh constructor" << std::endl
          << "The number of the refinement level must be smaller than "
          << 63-root_level+1 << "." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
  }
  else
    max_level = 63;

  //initialize user-enrollable functions
  if(mesh_size.x1rat!=1.0)
    user_meshgen_[x1dir]=true;
  else
    user_meshgen_[x1dir]=false;
  if(mesh_size.x2rat!=1.0)
    user_meshgen_[x2dir]=true;
  else
    user_meshgen_[x2dir]=false;
  if(mesh_size.x3rat!=1.0)
    user_meshgen_[x3dir]=true;
  else
    user_meshgen_[x3dir]=false;
  MeshGenerator_[x1dir]=DefaultMeshGeneratorX1;
  MeshGenerator_[x2dir]=DefaultMeshGeneratorX2;
  MeshGenerator_[x3dir]=DefaultMeshGeneratorX3;
  for(int dir=0; dir<6; dir++)
    BoundaryFunction_[dir]=NULL;
  AMRFlag_=NULL;
  UserSourceTerm_=NULL;
  InitUserMeshProperties(pin);

  InputBlock *pib = pin->pfirst_block;
  while (pib != NULL) {
    if (pib->block_name.compare(0,10,"refinement") == 0) {
      RegionSize ref_size;
      ref_size.x1min=pin->GetReal(pib->block_name,"x1min");
      ref_size.x1max=pin->GetReal(pib->block_name,"x1max");
      if(dim>=2) {
        ref_size.x2min=pin->GetReal(pib->block_name,"x2min");
        ref_size.x2max=pin->GetReal(pib->block_name,"x2max");
      }
      else {
        ref_size.x2min=mesh_size.x2min;
        ref_size.x2max=mesh_size.x2max;
      }
      if(dim>=3) {
        ref_size.x3min=pin->GetReal(pib->block_name,"x3min");
        ref_size.x3max=pin->GetReal(pib->block_name,"x3max");
      }
      else {
        ref_size.x3min=mesh_size.x3min;
        ref_size.x3max=mesh_size.x3max;
      }
      int ref_lev=pin->GetReal(pib->block_name,"level");
      int lrlev=ref_lev+root_level;
      if(lrlev>current_level) current_level=lrlev;
      if(lrlev!=root_level)
        multilevel=true;
      // range check
      if(ref_lev<1) {
        msg << "### FATAL ERROR in Mesh constructor" << std::endl
            << "Refinement level must be larger than 0 (root level = 0)" << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      if(lrlev > max_level) {
        msg << "### FATAL ERROR in Mesh constructor" << std::endl
            << "Refinement level exceeds the maximum level (specify maxlevel in <mesh> if adaptive)."
            << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      if(ref_size.x1min > ref_size.x1max || ref_size.x2min > ref_size.x2max
      || ref_size.x3min > ref_size.x3max)  {
        msg << "### FATAL ERROR in Mesh constructor" << std::endl
            << "Invalid refinement region is specified."<<  std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      if(ref_size.x1min < mesh_size.x1min || ref_size.x1max > mesh_size.x1max
      || ref_size.x2min < mesh_size.x2min || ref_size.x2max > mesh_size.x2max
      || ref_size.x3min < mesh_size.x3min || ref_size.x3max > mesh_size.x3max) {
        msg << "### FATAL ERROR in Mesh constructor" << std::endl
            << "Refinement region must be smaller than the whole mesh." << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      // find the logical range in the ref_level
      // note: if this is too slow, this should be replaced with bi-section search.
      long int lx1min=0, lx1max=0, lx2min=0, lx2max=0, lx3min=0, lx3max=0;
      long int lxmax=nrbx1*(1L<<ref_lev);
      for(lx1min=0;lx1min<lxmax;lx1min++) {
        if(MeshGenerator_[x1dir]((Real)(lx1min+1)/lxmax,mesh_size)>ref_size.x1min)
          break;
      }
      for(lx1max=lx1min;lx1max<lxmax;lx1max++) {
        if(MeshGenerator_[x1dir]((Real)(lx1max+1)/lxmax,mesh_size)>=ref_size.x1max)
          break;
      }
      if(lx1min%2==1) lx1min--;
      if(lx1max%2==0) lx1max++;
      if(dim>=2) { // 2D or 3D
        lxmax=nrbx2*(1L<<ref_lev);
        for(lx2min=0;lx2min<lxmax;lx2min++) {
          if(MeshGenerator_[x2dir]((Real)(lx2min+1)/lxmax,mesh_size)>ref_size.x2min)
            break;
        }
        for(lx2max=lx2min;lx2max<lxmax;lx2max++) {
          if(MeshGenerator_[x2dir]((Real)(lx2max+1)/lxmax,mesh_size)>=ref_size.x2max)
            break;
        }
        if(lx2min%2==1) lx2min--;
        if(lx2max%2==0) lx2max++;
      }
      if(dim==3) { // 3D
        lxmax=nrbx3*(1L<<ref_lev);
        for(lx3min=0;lx3min<lxmax;lx3min++) {
          if(MeshGenerator_[x3dir]((Real)(lx3min+1)/lxmax,mesh_size)>ref_size.x3min)
            break;
        }
        for(lx3max=lx3min;lx3max<lxmax;lx3max++) {
          if(MeshGenerator_[x3dir]((Real)(lx3max+1)/lxmax,mesh_size)>=ref_size.x3max)
            break;
        }
        if(lx3min%2==1) lx3min--;
        if(lx3max%2==0) lx3max++;
      }
      // create the finest level
      std::cout << "refinenment: logical level = " << lrlev << ", lx1min = "
                << lx1min << ", lx1max = " << lx1max << ", lx2min = " << lx2min
                << ", lx2max = " << lx2max << ", lx3min = " << lx3min << ", lx3max = "
                << lx3max << std::endl;
      if(dim==1) {
        for(long int i=lx1min; i<lx1max; i+=2) {
          LogicalLocation nloc;
          nloc.level=lrlev, nloc.lx1=i, nloc.lx2=0, nloc.lx3=0;
          int nnew;
          tree.AddMeshBlock(tree,nloc,dim,mesh_bcs,nrbx1,nrbx2,nrbx3,root_level,nnew);
        }
      }
      if(dim==2) {
        for(long int j=lx2min; j<lx2max; j+=2) {
          for(long int i=lx1min; i<lx1max; i+=2) {
            LogicalLocation nloc;
            nloc.level=lrlev, nloc.lx1=i, nloc.lx2=j, nloc.lx3=0;
            int nnew;
            tree.AddMeshBlock(tree,nloc,dim,mesh_bcs,nrbx1,nrbx2,nrbx3,root_level,nnew);
          }
        }
      }
      if(dim==3) {
        for(long int k=lx3min; k<lx3max; k+=2) {
          for(long int j=lx2min; j<lx2max; j+=2) {
            for(long int i=lx1min; i<lx1max; i+=2) {
              LogicalLocation nloc;
              nloc.level=lrlev, nloc.lx1=i, nloc.lx2=j, nloc.lx3=k;
              int nnew;
              tree.AddMeshBlock(tree,nloc,dim,mesh_bcs,nrbx1,nrbx2,nrbx3,root_level,nnew);
            }
          }
        }
      }
    }
    pib=pib->pnext;
  }

  if(multilevel==true) {
    if(block_size.nx1%2==1 || (block_size.nx2%2==1 && block_size.nx2>1)
                           || (block_size.nx3%2==1 && block_size.nx3>1)) {
      msg << "### FATAL ERROR in Mesh constructor" << std::endl
      << "The size of MeshBlock must be divisible by 2 in order to use SMR or AMR."
      << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
  }

  face_only=true;
  if (MAGNETIC_FIELDS_ENABLED || multilevel==true || VISCOSITY)
    face_only=false;

  maxneighbor_=BufferID(dim, multilevel, face_only);

  // initial mesh hierarchy construction is completed here

  tree.CountMeshBlock(nbtotal);
  loclist=new LogicalLocation[nbtotal];
  tree.GetMeshBlockList(loclist,NULL,nbtotal);

// check if there are sufficient blocks
#ifdef MPI_PARALLEL
  if(nbtotal < Globals::nranks) {
    if(test_flag==0) {
      msg << "### FATAL ERROR in Mesh constructor" << std::endl
          << "Too few blocks: nbtotal (" << nbtotal << ") < nranks ("<< Globals::nranks
          << ")" << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    else { // test
      std::cout << "### Warning in Mesh constructor" << std::endl
          << "Too few blocks: nbtotal (" << nbtotal << ") < nranks ("<< Globals::nranks
          << ")" << std::endl;
    }
  }
#endif

  ranklist=new int[nbtotal];
  nslist=new int[Globals::nranks];
  nblist=new int[Globals::nranks];
  costlist=new Real[nbtotal];
  if(adaptive==true) { // allocate arrays for AMR
    nref = new int [Globals::nranks];
    nderef = new int [Globals::nranks];
    rdisp = new int [Globals::nranks];
    ddisp = new int [Globals::nranks];
    bnref = new int [Globals::nranks];
    bnderef = new int [Globals::nranks];
    brdisp = new int [Globals::nranks];
    bddisp = new int [Globals::nranks];
  }

  for(int i=0;i<nbtotal;i++)
    costlist[i]=1.0; // the simplest estimate; all the blocks are equal

  LoadBalancing(costlist, ranklist, nslist, nblist, nbtotal);

  // Mesh test only; do not create meshes
  if(test_flag>0) {
    if(Globals::my_rank==0)
      MeshTest(dim);
    return;
  }

  // create MeshBlock list for this process
  int nbs=nslist[Globals::my_rank];
  int nbe=nbs+nblist[Globals::my_rank]-1;
// create MeshBlock list for this process
  for(int i=nbs;i<=nbe;i++) {
    SetBlockSizeAndBoundaries(loclist[i], block_size, block_bcs);
    // create a block and add into the link list
    if(i==nbs) {
      pblock = new MeshBlock(i, i-nbs, loclist[i], block_size, block_bcs, this, pin);
      pfirst = pblock;
    }
    else {
      pblock->next = new MeshBlock(i, i-nbs, loclist[i], block_size, block_bcs, this, pin);
      pblock->next->prev = pblock;
      pblock = pblock->next;
    }

    pblock->SearchAndSetNeighbors(tree, ranklist, nslist);
  }
  pblock=pfirst;

// create new Task List, requires mesh to already be constructed
  ptlist = new TaskList(this);

}


//--------------------------------------------------------------------------------------
// Mesh constructor for restarting. Load the restarting file

Mesh::Mesh(ParameterInput *pin, IOWrapper& resfile, int test_flag)
{
  std::stringstream msg;
  RegionSize block_size;
  MeshBlock *pfirst;
  int i, j, nerr, dim;
  IOWrapperSize_t *offset;

// mesh test
  if(test_flag>0) Globals::nranks=test_flag;

// read time and cycle limits from input file

  start_time = pin->GetOrAddReal("time","start_time",0.0);
  tlim       = pin->GetReal("time","tlim");
  cfl_number = pin->GetReal("time","cfl_number");
  nlim = pin->GetOrAddInteger("time","nlim",-1);

// read number of OpenMP threads for mesh
  num_mesh_threads_ = pin->GetOrAddInteger("mesh","num_threads",1);
  if (num_mesh_threads_ < 1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Number of OpenMP threads must be >= 1, but num_threads=" 
        << num_mesh_threads_ << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // read from the restarting file (everyone)
  // the file is already open and the pointer is set to after <par_end>
  nerr=0;
  if(resfile.Read(&nbtotal, sizeof(int), 1)!=1) nerr++;
  if(resfile.Read(&root_level, sizeof(int), 1)!=1) nerr++;
  current_level=root_level;
  if(resfile.Read(&mesh_size, sizeof(RegionSize), 1)!=1) nerr++;
  if(resfile.Read(mesh_bcs, sizeof(enum BoundaryFlag), 6)!=6) nerr++;
  if(resfile.Read(&time, sizeof(Real), 1)!=1) nerr++;
  if(resfile.Read(&dt, sizeof(Real), 1)!=1) nerr++;
  if(resfile.Read(&ncycle, sizeof(int), 1)!=1) nerr++;
  if(nerr>0) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "The restarting file is broken." << std::endl;
    resfile.Close();
    throw std::runtime_error(msg.str().c_str());
  }

  max_level = pin->GetOrAddInteger("mesh","maxlevel",1)+root_level-1;

  dim=1;
  if(mesh_size.nx2>1) dim=2;
  if(mesh_size.nx3>1) dim=3;

// check cfl_number
  if(cfl_number > 1.0 && mesh_size.nx2==1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "The CFL number must be smaller than 1.0 in 1D simulation" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if(cfl_number > 0.5 && mesh_size.nx2 > 1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "The CFL number must be smaller than 0.5 in 2D/3D simulation" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  //initialize
  loclist=new LogicalLocation[nbtotal];
  offset=new IOWrapperSize_t[nbtotal];
  costlist=new Real[nbtotal];
  ranklist=new int[nbtotal];
  nslist=new int[Globals::nranks];
  nblist=new int[Globals::nranks];

  int nx1 = pin->GetOrAddReal("meshblock","nx1",mesh_size.nx1);
  int nx2 = pin->GetOrAddReal("meshblock","nx2",mesh_size.nx2);
  int nx3 = pin->GetOrAddReal("meshblock","nx3",mesh_size.nx3);

// calculate the number of the blocks
  nrbx1=mesh_size.nx1/nx1;
  nrbx2=mesh_size.nx2/nx2;
  nrbx3=mesh_size.nx3/nx3;

  // read the id list (serial, because we need the costs for load balancing)
  // ... perhaps I should pack them.
  multilevel=false;
  nerr=0;
  for(int i=0;i<nbtotal;i++) {
    int bgid;
    if(resfile.Read(&bgid,sizeof(int),1)!=1) nerr++;
    if(resfile.Read(&(loclist[i]),sizeof(LogicalLocation),1)!=1) nerr++;
    if(loclist[i].level!=root_level) multilevel=true;
    if(loclist[i].level>current_level) current_level=loclist[i].level;
    if(resfile.Read(&(costlist[i]),sizeof(Real),1)!=1) nerr++;
    if(resfile.Read(&(offset[i]),sizeof(IOWrapperSize_t),1)!=1) nerr++;
  }
  if(nerr>0) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "The restarting file is broken." << std::endl;
    resfile.Close();
    throw std::runtime_error(msg.str().c_str());
  }

  adaptive=false;
  if(pin->GetOrAddString("mesh","refinement","static")=="adaptive")
    adaptive=true, multilevel=true;

  //initialize user-enrollable functions
  if(mesh_size.x1rat!=1.0)
    user_meshgen_[x1dir]=true;
  else
    user_meshgen_[x1dir]=false;
  if(mesh_size.x2rat!=1.0)
    user_meshgen_[x2dir]=true;
  else
    user_meshgen_[x2dir]=false;
  if(mesh_size.x3rat!=1.0)
    user_meshgen_[x3dir]=true;
  else
    user_meshgen_[x3dir]=false;
  MeshGenerator_[x1dir]=DefaultMeshGeneratorX1;
  MeshGenerator_[x2dir]=DefaultMeshGeneratorX2;
  MeshGenerator_[x3dir]=DefaultMeshGeneratorX3;
  for(int dir=0; dir<6; dir++)
    BoundaryFunction_[dir]=NULL;
  AMRFlag_=NULL;
  UserSourceTerm_=NULL;
  InitUserMeshProperties(pin);

  face_only=true;
  if (MAGNETIC_FIELDS_ENABLED || multilevel==true || VISCOSITY)
    face_only=false;

  maxneighbor_=BufferID(dim, multilevel, face_only);

  // rebuild the Block Tree
  for(int i=0;i<nbtotal;i++)
    tree.AddMeshBlockWithoutRefine(loclist[i],nrbx1,nrbx2,nrbx3,root_level);
  int nnb;
  // check the tree structure, and assign GID
  tree.GetMeshBlockList(loclist, NULL, nnb);
  if(nnb!=nbtotal) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Tree reconstruction failed. The total numbers of the blocks do not match. ("
        << nbtotal << " != " << nnb << ")" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

#ifdef MPI_PARALLEL
  if(nbtotal < Globals::nranks) {
    if(test_flag==0) {
      msg << "### FATAL ERROR in Mesh constructor" << std::endl
          << "Too few blocks: nbtotal (" << nbtotal << ") < nranks ("<< Globals::nranks
          << ")" << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    else { // test
      std::cout << "### Warning in Mesh constructor" << std::endl
          << "Too few blocks: nbtotal (" << nbtotal << ") < nranks ("<< Globals::nranks
          << ")" << std::endl;
      return;
    }
  }
#endif

  if(adaptive==true) { // allocate arrays for AMR
    nref = new int [Globals::nranks];
    nderef = new int [Globals::nranks];
    rdisp = new int [Globals::nranks];
    ddisp = new int [Globals::nranks];
    bnref = new int [Globals::nranks];
    bnderef = new int [Globals::nranks];
    brdisp = new int [Globals::nranks];
    bddisp = new int [Globals::nranks];
  }

  LoadBalancing(costlist, ranklist, nslist, nblist, nbtotal);

  // Mesh test only; do not create meshes
  if(test_flag>0) {
    if(Globals::my_rank==0)
      MeshTest(dim);
    delete [] offset;
    return;
  }

  // load MeshBlocks (parallel)
  int nbs=nslist[Globals::my_rank];
  int nbe=nbs+nblist[Globals::my_rank]-1;  
  for(i=nbs;i<=nbe;i++) {
    // create a block and add into the link list
    if(i==nbs) {
      pblock = new MeshBlock(i, i-nbs, this, pin, loclist[i], resfile, offset[i],
                             costlist[i], ranklist, nslist);
      pfirst = pblock;
    }
    else {
      pblock->next = new MeshBlock(i, i-nbs, this, pin, loclist[i], resfile,
                                   offset[i], costlist[i], ranklist, nslist);
      pblock->next->prev = pblock;
      pblock = pblock->next;
    }
    pblock->SearchAndSetNeighbors(tree, ranklist, nslist);
  }
  pblock=pfirst;

// create new Task List
  ptlist = new TaskList(this);

// clean up
  delete [] offset;
}


// destructor

Mesh::~Mesh()
{
  TerminateUserMeshProperties();
  while(pblock->prev != NULL) // should not be true
    delete pblock->prev;
  while(pblock->next != NULL)
    delete pblock->next;
  delete pblock;
  delete ptlist;
  delete [] nslist;
  delete [] nblist;
  delete [] ranklist;
  delete [] costlist;
  delete [] loclist;
  if(adaptive==true) { // allocate arrays for AMR
    delete [] nref;
    delete [] nderef;
    delete [] rdisp;
    delete [] ddisp;
    delete [] bnref;
    delete [] bnderef;
    delete [] brdisp;
    delete [] bddisp;
  }

}


//--------------------------------------------------------------------------------------
//! \fn void Mesh::MeshTest(int dim)
//  \brief print the mesh structure information

void Mesh::MeshTest(int dim)
{
  int i, j, nbt=0;
  long int lx1, lx2, lx3;
  int ll;
  Real mycost=0, mincost=FLT_MAX, maxcost=0.0, totalcost=0.0;
  int *nb=new int [max_level-root_level+1];
  FILE *fp;
  if(dim>=2) {
    if ((fp = fopen("meshtest.dat","wb")) == NULL) {
      std::cout << "### ERROR in function Mesh::MeshTest" << std::endl
                << "Cannot open meshtest.dat" << std::endl;
      return;
    }
  }

  std::cout << "Logical level of the physical root grid = "<< root_level << std::endl;
  std::cout << "Logical level of maximum refinement = "<< current_level << std::endl;
  std::cout << "List of MeshBlocks" << std::endl;
  for(i=root_level;i<=max_level;i++) {
    Real dx=1.0/(Real)(1L<<i);
    nb[i-root_level]=0;
    for(j=0;j<nbtotal;j++) {
      if(loclist[j].level==i) {
        long int &lx1=loclist[j].lx1;
        long int &lx2=loclist[j].lx2;
        long int &lx3=loclist[j].lx3;
        int &ll=loclist[j].level;
        std::cout << "MeshBlock " << j << ", lx1 = "
                  << loclist[j].lx1 << ", lx2 = " << lx2 <<", lx3 = " << lx3
                  << ", logical level = " << ll << ", physical level = "
                  << ll-root_level << ", cost = " << costlist[j]
                  << ", rank = " << ranklist[j] << std::endl;
        mincost=std::min(mincost,costlist[i]);
        maxcost=std::max(maxcost,costlist[i]);
        totalcost+=costlist[i];
        nb[i-root_level]++;
        if(dim==2) {
          fprintf(fp, "#MeshBlock %d at %ld %ld %ld %d\n", j, lx1, lx2, lx3, ll);
          fprintf(fp, "%g %g %d %d\n", lx1*dx, lx2*dx, ll, ranklist[j]);
          fprintf(fp, "%g %g %d %d\n", lx1*dx+dx, lx2*dx, ll, ranklist[j]);
          fprintf(fp, "%g %g %d %d\n", lx1*dx+dx, lx2*dx+dx, ll, ranklist[j]);
          fprintf(fp, "%g %g %d %d\n", lx1*dx, lx2*dx+dx, ll, ranklist[j]);
          fprintf(fp, "%g %g %d %d\n\n\n", lx1*dx, lx2*dx, ll, ranklist[j]);
        }
        if(dim==3) {
          fprintf(fp, "#MeshBlock %d at %ld %ld %ld %d\n", j, lx1, lx2, lx3, ll);
          fprintf(fp, "%g %g %g %d %d\n", lx1*dx, lx2*dx, lx3*dx, ll, ranklist[j]);
          fprintf(fp, "%g %g %g %d %d\n", lx1*dx+dx, lx2*dx, lx3*dx, ll, ranklist[j]);
          fprintf(fp, "%g %g %g %d %d\n", lx1*dx+dx, lx2*dx+dx, lx3*dx, ll, ranklist[j]);
          fprintf(fp, "%g %g %g %d %d\n", lx1*dx, lx2*dx+dx, lx3*dx, ll, ranklist[j]);
          fprintf(fp, "%g %g %g %d %d\n", lx1*dx, lx2*dx, lx3*dx, ll, ranklist[j]);
          fprintf(fp, "%g %g %g %d %d\n", lx1*dx, lx2*dx, lx3*dx+dx, ll, ranklist[j]);
          fprintf(fp, "%g %g %g %d %d\n", lx1*dx+dx, lx2*dx, lx3*dx+dx, ll, ranklist[j]);
          fprintf(fp, "%g %g %g %d %d\n", lx1*dx+dx, lx2*dx, lx3*dx, ll, ranklist[j]);
          fprintf(fp, "%g %g %g %d %d\n", lx1*dx+dx, lx2*dx, lx3*dx+dx, ll, ranklist[j]);
          fprintf(fp, "%g %g %g %d %d\n", lx1*dx+dx, lx2*dx+dx, lx3*dx+dx, ll, ranklist[j]);
          fprintf(fp, "%g %g %g %d %d\n", lx1*dx+dx, lx2*dx+dx, lx3*dx, ll, ranklist[j]);
          fprintf(fp, "%g %g %g %d %d\n", lx1*dx+dx, lx2*dx+dx, lx3*dx+dx, ll, ranklist[j]);
          fprintf(fp, "%g %g %g %d %d\n", lx1*dx, lx2*dx+dx, lx3*dx+dx, ll, ranklist[j]);
          fprintf(fp, "%g %g %g %d %d\n", lx1*dx, lx2*dx+dx, lx3*dx, ll, ranklist[j]);
          fprintf(fp, "%g %g %g %d %d\n", lx1*dx, lx2*dx+dx, lx3*dx+dx, ll, ranklist[j]);
          fprintf(fp, "%g %g %g %d %d\n", lx1*dx, lx2*dx, lx3*dx+dx, ll, ranklist[j]);
          fprintf(fp, "%g %g %g %d %d\n\n\n", lx1*dx, lx2*dx, lx3*dx, ll, ranklist[j]);
        }
      }
    }
  }
  if(dim>=2) fclose(fp);

  std::cout << std::endl;

  for(i=root_level;i<=max_level;i++) {
    if(nb[i-root_level]!=0)
      std::cout << "Level " << i-root_level << " (logical level " << i << ") : "
        << nb[i-root_level] << " MeshBlocks" << std::endl;
  }

  std::cout << "Total : " << nbtotal << " MeshBlocks" << std::endl << std::endl;
  std::cout << "Load Balance :" << std::endl;
  std::cout << "Minimum cost = " << mincost << ", Maximum cost = " << maxcost
            << ", Average cost = " << totalcost/nbtotal << std::endl;
  j=0;
  nbt=0;
  for(i=0;i<nbtotal;i++) {
    if(ranklist[i]==j) {
      mycost+=costlist[i];
      nbt++;
    }
    else if(ranklist[i]!=j) {
      std::cout << "Rank " << j << ": " << nbt <<" MeshBlocks, cost = " << mycost << std::endl;
      mycost=costlist[i];
      nbt=1;
      j++;
    }
  }
  std::cout << "Rank " << j << ": " << nbt <<" MeshBlocks, cost = " << mycost << std::endl;

  delete [] nb;
  return;
}

//--------------------------------------------------------------------------------------
// MeshBlock constructor: builds 1D vectors of cell positions and spacings, and
// constructs coordinate, boundary condition, hydro and field objects.

MeshBlock::MeshBlock(int igid, int ilid, LogicalLocation iloc, RegionSize input_block,
                     enum BoundaryFlag *input_bcs, Mesh *pm, ParameterInput *pin)
{
  std::stringstream msg;
  int root_level;
  pmy_mesh = pm;
  root_level = pm->root_level;
  block_size = input_block;
  for(int i=0; i<6; i++) block_bcs[i] = input_bcs[i];
  prev=NULL;
  next=NULL;
  gid=igid;
  lid=ilid;
  loc=iloc;
  cost=1.0;

// initialize grid indices

  is = NGHOST;
  ie = is + block_size.nx1 - 1;

  if (block_size.nx2 > 1) {
    js = NGHOST;
    je = js + block_size.nx2 - 1;
  } else {
    js = je = 0;
  }

  if (block_size.nx3 > 1) {
    ks = NGHOST;
    ke = ks + block_size.nx3 - 1;
  } else {
    ks = ke = 0;
  }

  if(pm->multilevel==true) {
    cnghost=(NGHOST+1)/2+1;
    cis=cnghost; cie=cis+block_size.nx1/2-1;
    cjs=cje=cks=cke=0;
    if(block_size.nx2>1) // 2D or 3D
      cjs=cnghost, cje=cjs+block_size.nx2/2-1;
    if(block_size.nx3>1) // 3D
      cks=cnghost, cke=cks+block_size.nx3/2-1;
  }

  if(pm->adaptive==false) { // too noisy for AMR
    std::cout << "MeshBlock " << gid << ", rank = " << Globals::my_rank << ", lx1 = "
              << loc.lx1 << ", lx2 = " << loc.lx2 <<", lx3 = " << loc.lx3
              << ", level = " << loc.level << std::endl;
    std::cout << "is=" << is << " ie=" << ie << " x1min=" << block_size.x1min
              << " x1max=" << block_size.x1max << std::endl;
    std::cout << "js=" << js << " je=" << je << " x2min=" << block_size.x2min
              << " x2max=" << block_size.x2max << std::endl;
    std::cout << "ks=" << ks << " ke=" << ke << " x3min=" << block_size.x3min
              << " x3max=" << block_size.x3max << std::endl;
  }

// construct Coordinates and Hydro objects stored in MeshBlock class.  Note that the
// initial conditions for the hydro are set in problem generator called from main, not
// in the Hydro constructor
 
  pcoord = new Coordinates(this, pin);
  if(pm->multilevel==true)
    pmr = new MeshRefinement(this, pin);
  phydro = new Hydro(this, pin);
  if (MAGNETIC_FIELDS_ENABLED)
    pfield = new Field(this, pin);
  pbval  = new BoundaryValues(this, pin);

  return;
}

//--------------------------------------------------------------------------------------
// MeshBlock constructor for restarting

MeshBlock::MeshBlock(int igid, int ilid, Mesh *pm, ParameterInput *pin,
                     LogicalLocation iloc, IOWrapper& resfile, IOWrapperSize_t offset,
                     Real icost, int *ranklist, int *nslist)
{
  std::stringstream msg;
  pmy_mesh = pm;
  prev=NULL;
  next=NULL;
  gid=igid;
  lid=ilid;
  loc=iloc;
  cost=icost;
  int nerr=0;
//  task=NULL;

  // seek the file
  resfile.Seek(offset);
  // load block structure
  if(resfile.Read(&block_size, sizeof(RegionSize), 1)!=1) nerr++;
  if(resfile.Read(block_bcs, sizeof(int), 6)!=6) nerr++;

  if(nerr>0) {
    msg << "### FATAL ERROR in MeshBlock constructor" << std::endl
        << "The restarting file is broken." << std::endl;
    resfile.Close();
    throw std::runtime_error(msg.str().c_str());
  }

// initialize grid indices

  is = NGHOST;
  ie = is + block_size.nx1 - 1;

  if (block_size.nx2 > 1) {
    js = NGHOST;
    je = js + block_size.nx2 - 1;
  } else {
    js = je = 0;
  }

  if (block_size.nx3 > 1) {
    ks = NGHOST;
    ke = ks + block_size.nx3 - 1;
  } else {
    ks = ke = 0;
  }

  if(pm->multilevel==true) {
    cnghost=(NGHOST+1)/2+1;
    cis=cnghost; cie=cis+block_size.nx1/2-1;
    cjs=cje=cks=cke=0;
    if(block_size.nx2>1) // 2D or 3D
      cjs=cnghost, cje=cjs+block_size.nx2/2-1;
    if(block_size.nx3>1) // 3D
      cks=cnghost, cke=cks+block_size.nx3/2-1;
  }

  std::cout << "MeshBlock " << gid << ", rank = " << Globals::my_rank << ", lx1 = "
            << loc.lx1 << ", lx2 = " << loc.lx2 <<", lx3 = " << loc.lx3
            << ", level = " << loc.level << std::endl;
  std::cout << "is=" << is << " ie=" << ie << " x1min=" << block_size.x1min
            << " x1max=" << block_size.x1max << std::endl;
  std::cout << "js=" << js << " je=" << je << " x2min=" << block_size.x2min
            << " x2max=" << block_size.x2max << std::endl;
  std::cout << "ks=" << ks << " ke=" << ke << " x3min=" << block_size.x3min
            << " x3max=" << block_size.x3max << std::endl;

  // create coordinates, hydro, field, and boundary conditions
  pcoord = new Coordinates(this, pin);
  if(pm->multilevel==true)
    pmr = new MeshRefinement(this, pin);
  phydro = new Hydro(this, pin);
  if (MAGNETIC_FIELDS_ENABLED)
    pfield = new Field(this, pin);
  pbval  = new BoundaryValues(this, pin);

  // load hydro and field data
  nerr=0;
  if(resfile.Read(phydro->u.GetArrayPointer(),sizeof(Real),
                         phydro->u.GetSize())!=phydro->u.GetSize()) nerr++;
  if (GENERAL_RELATIVITY) {
    if(resfile.Read(phydro->w.GetArrayPointer(),sizeof(Real),
                           phydro->w.GetSize())!=phydro->w.GetSize()) nerr++;
    if(resfile.Read(phydro->w1.GetArrayPointer(),sizeof(Real),
                           phydro->w1.GetSize())!=phydro->w1.GetSize()) nerr++;
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    if(resfile.Read(pfield->b.x1f.GetArrayPointer(),sizeof(Real),
               pfield->b.x1f.GetSize())!=pfield->b.x1f.GetSize()) nerr++;
    if(resfile.Read(pfield->b.x2f.GetArrayPointer(),sizeof(Real),
               pfield->b.x2f.GetSize())!=pfield->b.x2f.GetSize()) nerr++;
    if(resfile.Read(pfield->b.x3f.GetArrayPointer(),sizeof(Real),
               pfield->b.x3f.GetSize())!=pfield->b.x3f.GetSize()) nerr++;
  }
  if(nerr>0) {
    msg << "### FATAL ERROR in MeshBlock constructor" << std::endl
        << "The restarting file is broken." << std::endl;
    resfile.Close();
    throw std::runtime_error(msg.str().c_str());
  }
  return;
}

// destructor

MeshBlock::~MeshBlock()
{
  if(prev!=NULL) prev->next=next;
  if(next!=NULL) next->prev=prev;

  delete pcoord;
  delete phydro;
  if (MAGNETIC_FIELDS_ENABLED)
    delete pfield;
  delete pbval;
}


//--------------------------------------------------------------------------------------
// \!fn void Mesh::NewTimeStep(void)
// \brief function that loops over all MeshBlocks and find new timestep
//        this assumes that phydro->NewBlockTimeStep is already called

void Mesh::NewTimeStep(void)
{
  MeshBlock *pmb = pblock;
  Real min_dt=pmb->new_block_dt;
  pmb=pmb->next;
  while (pmb != NULL)  {
    min_dt=std::min(min_dt,pmb->new_block_dt);
    pmb=pmb->next;
  }
#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE,&min_dt,1,MPI_ATHENA_REAL,MPI_MIN,MPI_COMM_WORLD);
#endif
  // set it
  dt=std::min(min_dt*cfl_number,2.0*dt);
  if (time < tlim && tlim-time < dt)  // timestep would take us past desired endpoint
    dt = tlim-time;
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollUserBoundaryFunction(enum BoundaryFace dir, BValHydro_t my_bc)
//  \brief Enroll a user-defined boundary function

void Mesh::EnrollUserBoundaryFunction(enum BoundaryFace dir, BValFunc_t my_bc)
{
  std::stringstream msg;
  if(dir<0 || dir>5) {
    msg << "### FATAL ERROR in EnrollBoundaryCondition function" << std::endl
        << "dirName = " << dir << " not valid" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if(mesh_bcs[dir]!=USER_BNDRY) {
    msg << "### FATAL ERROR in EnrollUserBoundaryFunction" << std::endl
        << "The boundary condition flag must be set to the string 'user' in the "
        << " <mesh> block in the input file to use user-enrolled BCs" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  BoundaryFunction_[dir]=my_bc;
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollUserRefinementCondition(AMRFlag_t amrflag)
//  \brief Enroll a user-defined function for checking refinement criteria
void Mesh::EnrollUserRefinementCondition(AMRFlag_t amrflag)
{
  if(adaptive==true)
    AMRFlag_=amrflag;
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollUserMeshGenerator(enum direction, MeshGenFunc_t my_mg)
//  \brief Enroll a user-defined function for Mesh generation
void Mesh::EnrollUserMeshGenerator(enum direction dir, MeshGenFunc_t my_mg)
{
  std::stringstream msg;
  if(dir<0 || dir>3) {
    msg << "### FATAL ERROR in EnrollUserMeshGenerator function" << std::endl
        << "dirName = " << dir << " not valid" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  user_meshgen_[dir]=true;
  MeshGenerator_[dir]=my_mg;
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollUserSourceTermFunction(SrcTermFunc_t my_func)
//  \brief Enroll a user-defined source function

void Mesh::EnrollUserSourceTermFunction(SrcTermFunc_t my_func)
{
  UserSourceTerm_ = my_func;
  return;
}


//--------------------------------------------------------------------------------------
// \!fn void Mesh::Initialize(int res_flag, ParameterInput *pin)
// \brief  initialization before the main loop

void Mesh::Initialize(int res_flag, ParameterInput *pin)
{
  MeshBlock *pmb;
  Hydro *phydro;
  Field *pfield;
  BoundaryValues *pbval;
  std::stringstream msg;
  int inb=nbtotal;

  bool iflag=true;
  do {
    if(res_flag==0) {
      pmb = pblock;
      while (pmb != NULL)  {
        pmb->ProblemGenerator(pin);
        pmb->pbval->CheckBoundary();
        pmb=pmb->next;
      }
    }

    pmb = pblock;
    while (pmb != NULL)  {
      pmb->pbval->Initialize();
      pmb->pbval->StartReceivingForInit();
      pmb=pmb->next;
    }

    pmb = pblock;
    while (pmb != NULL)  {
      phydro=pmb->phydro;
      pfield=pmb->pfield;
      pbval=pmb->pbval;
      pbval->SendHydroBoundaryBuffers(phydro->u,0);
      if (MAGNETIC_FIELDS_ENABLED)
        pbval->SendFieldBoundaryBuffers(pfield->b,0);
      pmb=pmb->next;
    }

    pmb = pblock;
    while (pmb != NULL)  {
      phydro=pmb->phydro;
      pfield=pmb->pfield;
      pbval=pmb->pbval;
      pbval->ReceiveHydroBoundaryBuffersWithWait(phydro->u ,0);
      if (MAGNETIC_FIELDS_ENABLED)
        pbval->ReceiveFieldBoundaryBuffersWithWait(pfield->b ,0);
      pmb->pbval->ClearBoundaryForInit();
      if(multilevel==true)
        pbval->ProlongateBoundaries(phydro->w, phydro->u, pfield->b, pfield->bcc);

      int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
      if(pmb->nblevel[1][1][0]!=-1) is-=NGHOST;
      if(pmb->nblevel[1][1][2]!=-1) ie+=NGHOST;
      if(pmb->block_size.nx2 > 1) {
        if(pmb->nblevel[1][0][1]!=-1) js-=NGHOST;
        if(pmb->nblevel[1][2][1]!=-1) je+=NGHOST;
      }
      if(pmb->block_size.nx3 > 1) {
        if(pmb->nblevel[0][1][1]!=-1) ks-=NGHOST;
        if(pmb->nblevel[2][1][1]!=-1) ke+=NGHOST;
      }
      phydro->peos->ConservedToPrimitive(phydro->u, phydro->w1, pfield->b, 
                                         phydro->w, pfield->bcc, pmb->pcoord,
                                         is, ie, js, je, ks, ke);
      pbval->ApplyPhysicalBoundaries(phydro->w, phydro->u, pfield->b, pfield->bcc);
      pmb=pmb->next;
    }
    if((res_flag==0) && (adaptive==true)) {
      iflag=false;
      int onb=nbtotal;
      pmb = pblock;
      while (pmb != NULL)  {
        pmb->pmr->CheckRefinementCondition();
        pmb=pmb->next;
      }
      AdaptiveMeshRefinement(pin);
      if(nbtotal==onb) iflag=true;
      else if(nbtotal < onb && Globals::my_rank==0) {
         std::cout << "### Warning in Mesh::Initialize" << std::endl
         << "The number of MeshBlocks decreased during AMR grid initialization." << std::endl
         << "Possibly the refinement criteria have a problem." << std::endl;
      }
      if(nbtotal > 2*inb && Globals::my_rank==0) {
         std::cout << "### Warning in Mesh::Initialize" << std::endl
         << "The number of MeshBlocks increased more than twice during initialization."<< std::endl
         << "More computing power than you expected may be required." << std::endl;
      }
    }
  } while(iflag==false);

  if(res_flag==0 || res_flag==2) {
    pmb = pblock;
    while (pmb != NULL)  {
      pmb->phydro->NewBlockTimeStep(pmb);
      pmb=pmb->next;
    }
    NewTimeStep();
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn int64_t Mesh::GetTotalCells(void)
//  \brief return the total number of cells for performance counting

int64_t Mesh::GetTotalCells(void)
{
  return (int64_t)nbtotal*pblock->block_size.nx1*pblock->block_size.nx2*pblock->block_size.nx3;
}

//--------------------------------------------------------------------------------------
//! \fn long int MeshBlock::GetBlockSizeInBytes(void)
//  \brief Calculate the block data size required for restarting.

size_t MeshBlock::GetBlockSizeInBytes(void)
{
  size_t size;

  size =sizeof(RegionSize)+sizeof(int)*6;
  size+=sizeof(Real)*phydro->u.GetSize();
  if (GENERAL_RELATIVITY) {
    size+=sizeof(Real)*phydro->w.GetSize();
    size+=sizeof(Real)*phydro->w1.GetSize();
  }
  if (MAGNETIC_FIELDS_ENABLED)
    size+=sizeof(Real)*(pfield->b.x1f.GetSize()+pfield->b.x2f.GetSize()
                       +pfield->b.x3f.GetSize());
  // please add the size counter here when new physics is introduced

  return size;
}

//--------------------------------------------------------------------------------------
//! \fn void Mesh::UpdateOneStep(void)
//  \brief process the task list and advance one time step

void Mesh::UpdateOneStep(void)
{
  MeshBlock *pmb = pblock;
  int nb=nblist[Globals::my_rank];

  // initialize
  while (pmb != NULL)  {
    pmb->first_task=0;
    pmb->num_tasks_todo=ptlist->ntasks;
    for(int i=0; i<4; ++i) pmb->finished_tasks[i]=0; // encodes which tasks are done
    pmb->pbval->StartReceivingAll();
    pmb=pmb->next;
  }

  // main loop
  while(nb>0) {
    pmb = pblock;
    while (pmb != NULL)  {
      if(ptlist->DoOneTask(pmb)==TL_COMPLETE) // task list completed
        nb--;
      pmb=pmb->next;
    }
  }

  pmb = pblock;
  while (pmb != NULL)  {
    pmb->pbval->ClearBoundaryAll();
    pmb=pmb->next;
  }
  return;
}

//--------------------------------------------------------------------------------------
//! \fn MeshBlock* Mesh::FindMeshBlock(int tgid)
//  \brief return the MeshBlock whose gid is tgid

MeshBlock* Mesh::FindMeshBlock(int tgid)
{
  MeshBlock *pbl=pblock;
  while(pbl!=NULL)
  {
    if(pbl->gid==tgid)
      break;
    pbl=pbl->next;
  }
  return pbl;
}

//--------------------------------------------------------------------------------------
// \!fn void NeighborBlock::SetNeighbor(int irank, int ilevel, int igid, int ilid,
//                          int iox1, int iox2, int iox3, enum neighbor_type itype,
//                          int ibid, int itargetid, int ifi1=0, int ifi2=0,
//                          bool ipolar=false)
// \brief Set neighbor information

void NeighborBlock::SetNeighbor(int irank, int ilevel, int igid, int ilid,
  int iox1, int iox2, int iox3, enum neighbor_type itype, int ibid, int itargetid,
  bool ipolar, int ifi1=0, int ifi2=0)
{
  rank=irank; level=ilevel; gid=igid; lid=ilid; ox1=iox1; ox2=iox2; ox3=iox3;
  type=itype; bufid=ibid; targetid=itargetid; polar=ipolar; fi1=ifi1; fi2=ifi2;
  if(type==neighbor_face) {
    if(ox1==-1)      fid=INNER_X1;
    else if(ox1==1)  fid=OUTER_X1;
    else if(ox2==-1) fid=INNER_X2;
    else if(ox2==1)  fid=OUTER_X2;
    else if(ox3==-1) fid=INNER_X3;
    else if(ox3==1)  fid=OUTER_X3;
  }
  if(type==neighbor_edge) {
    if(ox3==0)      eid=(edgeid)(   ((ox1+1)>>1) | ((ox2+1)&2));
    else if(ox2==0) eid=(edgeid)(4+(((ox1+1)>>1) | ((ox3+1)&2)));
    else if(ox1==0) eid=(edgeid)(8+(((ox2+1)>>1) | ((ox3+1)&2)));
  }
  return;
}

//--------------------------------------------------------------------------------------
// \!fn void MeshBlock::SearchAndSetNeighbors(MeshBlockTree &tree, int *ranklist, int *nslist)
// \brief Search and set all the neighbor blocks

void MeshBlock::SearchAndSetNeighbors(MeshBlockTree &tree, int *ranklist, int *nslist)
{
  MeshBlockTree* neibt;
  int myox1, myox2=0, myox3=0, myfx1, myfx2, myfx3;
  myfx1=(int)(loc.lx1&1L);
  myfx2=(int)(loc.lx2&1L);
  myfx3=(int)(loc.lx3&1L);
  myox1=((int)(loc.lx1&1L))*2-1;
  if(block_size.nx2>1) myox2=((int)(loc.lx2&1L))*2-1;
  if(block_size.nx3>1) myox3=((int)(loc.lx3&1L))*2-1;
  long int nrbx1=pmy_mesh->nrbx1, nrbx2=pmy_mesh->nrbx2, nrbx3=pmy_mesh->nrbx3;

  int nf1=1, nf2=1;
  if(pmy_mesh->multilevel==true) {
    if(block_size.nx2>1) nf1=2;
    if(block_size.nx3>1) nf2=2;
  }
  int bufid=0;
  nneighbor=0;
  for(int k=0; k<=2; k++) {
    for(int j=0; j<=2; j++) {
      for(int i=0; i<=2; i++)
        nblevel[k][j][i]=-1;
    }
  }
  nblevel[1][1][1]=loc.level;

  // x1 face
  for(int n=-1; n<=1; n+=2) {
    neibt=tree.FindNeighbor(loc,n,0,0,block_bcs,nrbx1,nrbx2,nrbx3,pmy_mesh->root_level);
    if(neibt==NULL) { bufid+=nf1*nf2; continue;}
    if(neibt->flag==false) { // neighbor at finer level
      int fface=1-(n+1)/2; // 0 for OUTER_X1, 1 for INNER_X1
      nblevel[1][1][n+1]=neibt->loc.level+1;
      for(int f2=0;f2<nf2;f2++) {
        for(int f1=0;f1<nf1;f1++) {
          MeshBlockTree* nf=neibt->GetLeaf(fface,f1,f2);
          int fid = nf->gid;
          int nlevel=nf->loc.level;
          int tbid=FindBufferID(-n,0,0,0,0,pmy_mesh->maxneighbor_);
          neighbor[nneighbor].SetNeighbor(ranklist[fid], nlevel, fid,
              fid-nslist[ranklist[fid]], n, 0, 0, neighbor_face, bufid, tbid, false, f1,
              f2);
          bufid++; nneighbor++;
        }
      }
    }
    else { // neighbor at same or coarser level
      int nlevel=neibt->loc.level;
      int nid=neibt->gid;
      nblevel[1][1][n+1]=nlevel;
      int tbid;
      if(nlevel==loc.level) { // neighbor at same level
        tbid=FindBufferID(-n,0,0,0,0,pmy_mesh->maxneighbor_);
      }
      else { // neighbor at coarser level
        tbid=FindBufferID(-n,0,0,myfx2,myfx3,pmy_mesh->maxneighbor_);
      }
      neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
          nid-nslist[ranklist[nid]], n, 0, 0, neighbor_face, bufid, tbid, false);
      bufid+=nf1*nf2; nneighbor++;
    }
  }
  if(block_size.nx2==1) return;

  // x2 face
  for(int n=-1; n<=1; n+=2) {
    neibt=tree.FindNeighbor(loc,0,n,0,block_bcs,nrbx1,nrbx2,nrbx3,pmy_mesh->root_level);
    if(neibt==NULL) { bufid+=nf1*nf2; continue;}
    if(neibt->flag==false) { // neighbor at finer level
      int fface=1-(n+1)/2; // 0 for OUTER_X2, 1 for INNER_X2
      nblevel[1][n+1][1]=neibt->loc.level+1;
      for(int f2=0;f2<nf2;f2++) {
        for(int f1=0;f1<nf1;f1++) {
          MeshBlockTree* nf=neibt->GetLeaf(f1,fface,f2);
          int fid = nf->gid;
          int nlevel=nf->loc.level;
          int tbid=FindBufferID(0,-n,0,0,0,pmy_mesh->maxneighbor_);
          neighbor[nneighbor].SetNeighbor(ranklist[fid], nlevel, fid,
              fid-nslist[ranklist[fid]], 0, n, 0, neighbor_face, bufid, tbid, false, f1,
              f2);
          bufid++; nneighbor++;
        }
      }
    }
    else { // neighbor at same or coarser level
      int nlevel=neibt->loc.level;
      int nid=neibt->gid;
      nblevel[1][n+1][1]=nlevel;
      int tbid;
      bool polar=false;
      if(nlevel==loc.level) { // neighbor at same level
        long int num_x2 = nrbx2<<(loc.level-pmy_mesh->root_level);
        if ((loc.lx2+n<0 and block_bcs[INNER_X2]==POLAR_BNDRY)
            or (loc.lx2+n>=num_x2 and block_bcs[OUTER_X2]==POLAR_BNDRY)) {
          polar = true; // neighbor is across top or bottom pole
        }
        tbid=FindBufferID(0,polar?n:-n,0,0,0,pmy_mesh->maxneighbor_);
      }
      else { // neighbor at coarser level
        tbid=FindBufferID(0,-n,0,myfx1,myfx3,pmy_mesh->maxneighbor_);
      }
      neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
          nid-nslist[ranklist[nid]], 0, n, 0, neighbor_face, bufid, tbid, polar);
      bufid+=nf1*nf2; nneighbor++;
    }
  }

  // x3 face
  if(block_size.nx3>1) {
    for(int n=-1; n<=1; n+=2) {
      neibt=tree.FindNeighbor(loc,0,0,n,block_bcs,nrbx1,nrbx2,nrbx3,pmy_mesh->root_level);
      if(neibt==NULL) { bufid+=nf1*nf2; continue;}
      if(neibt->flag==false) { // neighbor at finer level
        int fface=1-(n+1)/2; // 0 for OUTER_X3, 1 for INNER_X3
        nblevel[n+1][1][1]=neibt->loc.level+1;
        for(int f2=0;f2<nf2;f2++) {
          for(int f1=0;f1<nf1;f1++) {
            MeshBlockTree* nf=neibt->GetLeaf(f1,f2,fface);
            int fid = nf->gid;
            int nlevel=nf->loc.level;
            int tbid=FindBufferID(0,0,-n,0,0,pmy_mesh->maxneighbor_);
            neighbor[nneighbor].SetNeighbor(ranklist[fid], nlevel, fid,
                fid-nslist[ranklist[fid]], 0, 0, n, neighbor_face, bufid, tbid, false,
                f1, f2);
            bufid++; nneighbor++;
          }
        }
      }
      else { // neighbor at same or coarser level
        int nlevel=neibt->loc.level;
        int nid=neibt->gid;
        nblevel[n+1][1][1]=nlevel;
        int tbid;
        if(nlevel==loc.level) { // neighbor at same level
          tbid=FindBufferID(0,0,-n,0,0,pmy_mesh->maxneighbor_);
        }
        else { // neighbor at coarser level
          tbid=FindBufferID(0,0,-n,myfx1,myfx2,pmy_mesh->maxneighbor_);
        }
        neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
            nid-nslist[ranklist[nid]], 0, 0, n, neighbor_face, bufid, tbid, false);
        bufid+=nf1*nf2; nneighbor++;
      }
    }
  }
  if(pmy_mesh->face_only==true) return;

  // x1x2 edge
  for(int m=-1; m<=1; m+=2) {
    for(int n=-1; n<=1; n+=2) {
      neibt=tree.FindNeighbor(loc,n,m,0,block_bcs,nrbx1,nrbx2,nrbx3,pmy_mesh->root_level);
      if(neibt==NULL) { bufid+=nf2; continue;}
      if(neibt->flag==false) { // neighbor at finer level
        int ff1=1-(n+1)/2; // 0 for OUTER_X1, 1 for INNER_X1
        int ff2=1-(m+1)/2; // 0 for OUTER_X2, 1 for INNER_X2
        nblevel[1][m+1][n+1]=neibt->loc.level+1;
        for(int f1=0;f1<nf2;f1++) {
          MeshBlockTree* nf=neibt->GetLeaf(ff1,ff2,f1);
          int fid = nf->gid;
          int nlevel=nf->loc.level;
          int tbid=FindBufferID(-n,-m,0,0,0,pmy_mesh->maxneighbor_);
          neighbor[nneighbor].SetNeighbor(ranklist[fid], nlevel, fid,
              fid-nslist[ranklist[fid]], n, m, 0, neighbor_edge, bufid, tbid, false, f1,
              0);
          bufid++; nneighbor++;
        }
      }
      else { // neighbor at same or coarser level
        int nlevel=neibt->loc.level;
        int nid=neibt->gid;
        nblevel[1][m+1][n+1]=nlevel;
        int tbid;
        bool polar=false;
        if(nlevel==loc.level) { // neighbor at same level
          long int num_x2 = nrbx2<<(loc.level-pmy_mesh->root_level);
          if ((loc.lx2+m<0 and block_bcs[INNER_X2]==POLAR_BNDRY)
              or (loc.lx2+m>=num_x2 and block_bcs[OUTER_X2]==POLAR_BNDRY)) {
            polar = true; // neighbor is across top or bottom pole
          }
          tbid=FindBufferID(-n,polar?m:-m,0,0,0,pmy_mesh->maxneighbor_);
        }
        else { // neighbor at coarser level
          tbid=FindBufferID(-n,polar?m:-m,0,myfx3,0,pmy_mesh->maxneighbor_);
        }
        if(nlevel>=loc.level || (myox1==n && myox2==m)) {
          neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
              nid-nslist[ranklist[nid]], n, m, 0, neighbor_edge, bufid, tbid, polar);
          nneighbor++;
        }
        bufid+=nf2;
      }
    }
  }
  if(block_size.nx3==1) return;

  // x1x3 edge
  for(int m=-1; m<=1; m+=2) {
    for(int n=-1; n<=1; n+=2) {
      neibt=tree.FindNeighbor(loc,n,0,m,block_bcs,nrbx1,nrbx2,nrbx3,pmy_mesh->root_level);
      if(neibt==NULL) { bufid+=nf1; continue;}
      if(neibt->flag==false) { // neighbor at finer level
        int ff1=1-(n+1)/2; // 0 for OUTER_X1, 1 for INNER_X1
        int ff2=1-(m+1)/2; // 0 for OUTER_X3, 1 for INNER_X3
        nblevel[m+1][1][n+1]=neibt->loc.level+1;
        for(int f1=0;f1<nf1;f1++) {
          MeshBlockTree* nf=neibt->GetLeaf(ff1,f1,ff2);
          int fid = nf->gid;
          int nlevel=nf->loc.level;
          int tbid=FindBufferID(-n,0,-m,0,0,pmy_mesh->maxneighbor_);
          neighbor[nneighbor].SetNeighbor(ranklist[fid], nlevel, fid,
              fid-nslist[ranklist[fid]], n, 0, m, neighbor_edge, bufid, tbid, false, f1,
              0);
          bufid++; nneighbor++;
        }
      }
      else { // neighbor at same or coarser level
        int nlevel=neibt->loc.level;
        int nid=neibt->gid;
        nblevel[m+1][1][n+1]=nlevel;
        int tbid;
        if(nlevel==loc.level) { // neighbor at same level
          tbid=FindBufferID(-n,0,-m,0,0,pmy_mesh->maxneighbor_);
        }
        else { // neighbor at coarser level
          tbid=FindBufferID(-n,0,-m,myfx2,0,pmy_mesh->maxneighbor_);
        }
        if(nlevel>=loc.level || (myox1==n && myox3==m)) {
          neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
              nid-nslist[ranklist[nid]], n, 0, m, neighbor_edge, bufid, tbid, false);
          nneighbor++;
        }
        bufid+=nf1;
      }
    }
  }

  // x2x3 edge
  for(int m=-1; m<=1; m+=2) {
    for(int n=-1; n<=1; n+=2) {
      neibt=tree.FindNeighbor(loc,0,n,m,block_bcs,nrbx1,nrbx2,nrbx3,pmy_mesh->root_level);
      if(neibt==NULL) { bufid+=nf1; continue;}
      if(neibt->flag==false) { // neighbor at finer level
        int ff1=1-(n+1)/2; // 0 for OUTER_X2, 1 for INNER_X2
        int ff2=1-(m+1)/2; // 0 for OUTER_X3, 1 for INNER_X3
        nblevel[m+1][n+1][1]=neibt->loc.level+1;
        for(int f1=0;f1<nf1;f1++) {
          MeshBlockTree* nf=neibt->GetLeaf(f1,ff1,ff2);
          int fid = nf->gid;
          int nlevel=nf->loc.level;
          int tbid=FindBufferID(0,-n,-m,0,0,pmy_mesh->maxneighbor_);
          neighbor[nneighbor].SetNeighbor(ranklist[fid], nlevel, fid,
              fid-nslist[ranklist[fid]], 0, n, m, neighbor_edge, bufid, tbid, false, f1,
              0);
          bufid++; nneighbor++;
        }
      }
      else { // neighbor at same or coarser level
        int nlevel=neibt->loc.level;
        int nid=neibt->gid;
        nblevel[m+1][n+1][1]=nlevel;
        int tbid;
        bool polar=false;
        if(nlevel==loc.level) { // neighbor at same level
          long int num_x2 = nrbx2<<(loc.level-pmy_mesh->root_level);
          if ((loc.lx2+n<0 and block_bcs[INNER_X2]==POLAR_BNDRY)
              or (loc.lx2+n>=num_x2 and block_bcs[OUTER_X2]==POLAR_BNDRY)) {
            polar = true; // neighbor is across top or bottom pole
          }
          tbid=FindBufferID(0,polar?n:-n,-m,0,0,pmy_mesh->maxneighbor_);
        }
        else { // neighbor at coarser level
          tbid=FindBufferID(0,-n,-m,myfx1,0,pmy_mesh->maxneighbor_);
        }
        if(nlevel>=loc.level || (myox2==n && myox3==m)) {
          neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
              nid-nslist[ranklist[nid]], 0, n, m, neighbor_edge, bufid, tbid, polar);
          nneighbor++;
        }
        bufid+=nf1;
      }
    }
  }

  // corners
  for(int l=-1; l<=1; l+=2) {
    for(int m=-1; m<=1; m+=2) {
      for(int n=-1; n<=1; n+=2) {
        neibt=tree.FindNeighbor(loc,n,m,l,block_bcs,nrbx1,nrbx2,nrbx3,pmy_mesh->root_level);
        if(neibt==NULL) { bufid++; continue;}
        bool polar=false;
        long int num_x2 = nrbx2<<(loc.level-pmy_mesh->root_level);
        if ((loc.lx2+m<0 and block_bcs[INNER_X2]==POLAR_BNDRY)
            or (loc.lx2+m>=num_x2 and block_bcs[OUTER_X2]==POLAR_BNDRY)) {
          polar = true; // neighbor is across top or bottom pole
        }
        if(neibt->flag==false) { // neighbor at finer level
          int ff1=1-(n+1)/2; // 0 for OUTER_X1, 1 for INNER_X1
          int ff2=1-(m+1)/2; // 0 for OUTER_X2, 1 for INNER_X2
          int ff3=1-(l+1)/2; // 0 for OUTER_X3, 1 for INNER_X3
          neibt=neibt->GetLeaf(ff1,ff2,ff3);
        }
        int nlevel=neibt->loc.level;
        nblevel[l+1][m+1][n+1]=nlevel;
        if(nlevel>=loc.level || (myox1==n && myox2==m && myox3==l)) {
          int nid=neibt->gid;
          int tbid=FindBufferID(-n,polar?m:-m,-l,0,0,pmy_mesh->maxneighbor_);
          neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
              nid-nslist[ranklist[nid]], n, m, l, neighbor_corner, bufid, tbid, polar);
          nneighbor++;
        }
        bufid++;
      }
    }
  }
  return;
}

//--------------------------------------------------------------------------------------
// \!fn void Mesh::TestConservation(void)
// \brief Calculate and print the total of conservative variables

void Mesh::TestConservation(void)
{
  MeshBlock *pmb = pblock;
  Real tcons[NHYDRO];
  for(int n=0;n<NHYDRO;n++) tcons[n]=0.0;
  while(pmb!=NULL) {
    pmb->IntegrateConservative(tcons);
    pmb=pmb->next;
  }

#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE,tcons,NHYDRO,MPI_ATHENA_REAL,MPI_SUM,MPI_COMM_WORLD);
#endif

  if(Globals::my_rank==0) {
    std::cout << "Total Conservative : " ;
    for(int n=0;n<NHYDRO;n++)
      std::cout << tcons[n] << " ";
    std::cout << std::endl;
  }

  return;
}

//--------------------------------------------------------------------------------------
// \!fn void MeshBlock::IntegrateConservative(Real *tcons)
// \brief Calculate and print the total of conservative variables

void MeshBlock::IntegrateConservative(Real *tcons)
{
  for(int n=0;n<NHYDRO;n++) {
    for(int k=ks;k<=ke;k++) {
      for(int j=js;j<=je;j++) {
        for(int i=is;i<=ie;i++)
          tcons[n]+=phydro->u(n,k,j,i)*pcoord->GetCellVolume(k,j,i);
      }
    }
  }
  return;
}


//--------------------------------------------------------------------------------------
// \!fn void Mesh::LoadBalancing(Real *clist, int *rlist, int *slist, int *nlist, int nb)
// \brief Calculate distribution of MeshBlocks based on the cost list
void Mesh::LoadBalancing(Real *clist, int *rlist, int *slist, int *nlist, int nb)
{
  std::stringstream msg;
  Real totalcost=0, maxcost=0.0, mincost=(FLT_MAX);

  for(int i=0; i<nb; i++) {
    totalcost+=clist[i];
    mincost=std::min(mincost,clist[i]);
    maxcost=std::max(maxcost,clist[i]);
  }
  int j=(Globals::nranks)-1;
  Real targetcost=totalcost/Globals::nranks;
  Real mycost=0.0;
  // create rank list from the end: the master node should have less load
  for(int i=nb-1;i>=0;i--) {
    if(targetcost==0.0) {
      msg << "### FATAL ERROR in LoadBalancing" << std::endl
          << "There is at least one process which has no MeshBlock" << std::endl
          << "Decrease the number of processes or use smaller MeshBlocks." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    mycost+=clist[i];
    rlist[i]=j;
    if(mycost >= targetcost && j>0) {
      j--;
      totalcost-=mycost;
      mycost=0.0;
      targetcost=totalcost/(j+1);
    }
  }
  slist[0]=0;
  j=0;
  for(int i=1;i<nb;i++) { // make the list of nbstart and nblocks
    if(rlist[i]!=rlist[i-1]) {
      nlist[j]=i-nslist[j];
      slist[++j]=i;
    }
  }
  nlist[j]=nb-slist[j];

#ifdef MPI_PARALLEL
  if(nb % Globals::nranks != 0 && adaptive == false
  && maxcost == mincost && Globals::my_rank==0) {
    std::cout << "### Warning in LoadBalancing" << std::endl
              << "The number of MeshBlocks cannot be divided evenly. "
              << "This will cause a poor load balance." << std::endl;
  }
#endif
  return;
}


//--------------------------------------------------------------------------------------
// \!fn void Mesh::SetBlockSizeAndBoundaries(LogicalLocation loc,
//                 RegionSize &block_size, enum BundaryFlag *block_bcs)
// \brief Set the physical part of a block_size structure and block boundary conditions
void Mesh::SetBlockSizeAndBoundaries(LogicalLocation loc, RegionSize &block_size,
                                     enum BoundaryFlag *block_bcs)
{
  long int &lx1=loc.lx1;
  long int &lx2=loc.lx2;
  long int &lx3=loc.lx3;
  int &ll=loc.level;
  // calculate physical block size, x1
  if(lx1==0) {
    block_size.x1min=mesh_size.x1min;
    block_bcs[INNER_X1]=mesh_bcs[INNER_X1];
  }
  else {
    Real rx=(Real)lx1/(Real)(nrbx1<<(ll-root_level));
    block_size.x1min=MeshGenerator_[x1dir](rx,mesh_size);
    block_bcs[INNER_X1]=BLOCK_BNDRY;
  }
  if(lx1==(nrbx1<<(ll-root_level))-1) {
    block_size.x1max=mesh_size.x1max;
    block_bcs[OUTER_X1]=mesh_bcs[OUTER_X1];
  }
  else {
    Real rx=(Real)(lx1+1)/(Real)(nrbx1<<(ll-root_level));
    block_size.x1max=MeshGenerator_[x1dir](rx,mesh_size);
    block_bcs[OUTER_X1]=BLOCK_BNDRY;
  }

  // calculate physical block size, x2
  if(mesh_size.nx2 == 1) {
    block_size.x2min=mesh_size.x2min;
    block_size.x2max=mesh_size.x2max;
    block_bcs[INNER_X2]=mesh_bcs[INNER_X2];
    block_bcs[OUTER_X2]=mesh_bcs[OUTER_X2];
  }
  else {
    if(lx2==0) {
      block_size.x2min=mesh_size.x2min;
      block_bcs[INNER_X2]=mesh_bcs[INNER_X2];
    }
    else {
      Real rx=(Real)lx2/(Real)(nrbx2<<(ll-root_level));
      block_size.x2min=MeshGenerator_[x2dir](rx,mesh_size);
      block_bcs[INNER_X2]=BLOCK_BNDRY;
    }
    if(lx2==(nrbx2<<(ll-root_level))-1) {
      block_size.x2max=mesh_size.x2max;
      block_bcs[OUTER_X2]=mesh_bcs[OUTER_X2];
    }
    else {
      Real rx=(Real)(lx2+1)/(Real)(nrbx2<<(ll-root_level));
      block_size.x2max=MeshGenerator_[x2dir](rx,mesh_size);
      block_bcs[OUTER_X2]=BLOCK_BNDRY;
    }
  }

  // calculate physical block size, x3
  if(mesh_size.nx3 == 1) {
    block_size.x3min=mesh_size.x3min;
    block_size.x3max=mesh_size.x3max;
    block_bcs[INNER_X3]=mesh_bcs[INNER_X3];
    block_bcs[OUTER_X3]=mesh_bcs[OUTER_X3];
  }
  else {
    if(lx3==0) {
      block_size.x3min=mesh_size.x3min;
      block_bcs[INNER_X3]=mesh_bcs[INNER_X3];
    }
    else {
      Real rx=(Real)lx3/(Real)(nrbx3<<(ll-root_level));
      block_size.x3min=MeshGenerator_[x3dir](rx,mesh_size);
      block_bcs[INNER_X3]=BLOCK_BNDRY;
    }
    if(lx3==(nrbx3<<(ll-root_level))-1) {
      block_size.x3max=mesh_size.x3max;
      block_bcs[OUTER_X3]=mesh_bcs[OUTER_X3];
    }
    else {
      Real rx=(Real)(lx3+1)/(Real)(nrbx3<<(ll-root_level));
      block_size.x3max=MeshGenerator_[x3dir](rx,mesh_size);
      block_bcs[OUTER_X3]=BLOCK_BNDRY;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------
// \!fn void Mesh::AdaptiveMeshRefinement(ParameterInput *pin)
// \brief Main function for adaptive mesh refinement
void Mesh::AdaptiveMeshRefinement(ParameterInput *pin)
{
  MeshBlock *pmb;
#ifdef MPI_PARALLEL
  MPI_Request areq[3];
  // start sharing the cost list
  MPI_Iallgatherv(MPI_IN_PLACE, nblist[Globals::my_rank], MPI_INT,
                  costlist, nblist, nslist, MPI_INT, MPI_COMM_WORLD, &areq[2]);
#endif

  // collect refinement flags from all the meshblocks
  // count the number of the blocks to be (de)refined
  nref[Globals::my_rank]=0;
  nderef[Globals::my_rank]=0;
  pmb=pblock;
  while(pmb!=NULL) {
    if(pmb->pmr->refine_flag_== 1) nref[Globals::my_rank]++;
    if(pmb->pmr->refine_flag_==-1) nderef[Globals::my_rank]++;
    pmb=pmb->next;
  }
#ifdef MPI_PARALLEL
  // if this does not work due to a version issue, replace these with blocking AllGather
  MPI_Iallgather(MPI_IN_PLACE, 1, MPI_INT, nref,   1, MPI_INT, MPI_COMM_WORLD, &areq[0]);
  MPI_Iallgather(MPI_IN_PLACE, 1, MPI_INT, nderef, 1, MPI_INT, MPI_COMM_WORLD, &areq[1]);
  MPI_Waitall(2, areq, MPI_STATUSES_IGNORE);
#endif

  // count the number of the blocks to be (de)refined and displacement
  int tnref=0, tnderef=0;
  for(int n=0; n<Globals::nranks; n++) {
    tnref  += nref[n];
    tnderef+= nderef[n];
  }
  if(tnref==0 && tnderef==0) // nothing to do
    return;

  int rd=0, dd=0;
  for(int n=0; n<Globals::nranks; n++) {
    bnref[n]   = nref[n]*sizeof(LogicalLocation);
    bnderef[n] = nderef[n]*sizeof(LogicalLocation);
    rdisp[n] = rd;
    ddisp[n] = dd;
    brdisp[n] = rd*sizeof(LogicalLocation);
    bddisp[n] = dd*sizeof(LogicalLocation);
    rd+=nref[n];
    dd+=nderef[n];
  }

  // allocate memory for the location arrays
  int nlbl=2, dim=1;
  if(mesh_size.nx2 > 1) nlbl=4, dim=2;
  if(mesh_size.nx3 > 1) nlbl=8, dim=3;
  LogicalLocation *lref, *lderef, *clderef;
  if(tnref!=0)
    lref = new LogicalLocation[tnref];
  if(tnderef>nlbl) {
    lderef = new LogicalLocation[tnderef];
    clderef = new LogicalLocation[tnderef/nlbl];
  }

  // collect the locations and costs
  int iref = rdisp[Globals::my_rank], ideref = ddisp[Globals::my_rank];
  pmb=pblock;
  while(pmb!=NULL) {
    if(pmb->pmr->refine_flag_== 1)
      lref[iref++]=pmb->loc;
    if(pmb->pmr->refine_flag_==-1 && tnderef>nlbl)
      lderef[ideref++]=pmb->loc;
    pmb=pmb->next;
  }
#ifdef MPI_PARALLEL
  if(tnref>0 && tnderef>nlbl) {
    MPI_Iallgatherv(MPI_IN_PLACE, bnref[Globals::my_rank],   MPI_BYTE,
                    lref,   bnref,   brdisp, MPI_BYTE, MPI_COMM_WORLD, &areq[0]);
    MPI_Iallgatherv(MPI_IN_PLACE, bnderef[Globals::my_rank], MPI_BYTE,
                    lderef, bnderef, bddisp, MPI_BYTE, MPI_COMM_WORLD, &areq[1]);
    MPI_Waitall(2, areq, MPI_STATUSES_IGNORE);
  }
  else if(tnref>0) {
    MPI_Allgatherv(MPI_IN_PLACE, bnref[Globals::my_rank],   MPI_BYTE,
                    lref,   bnref,   brdisp, MPI_BYTE, MPI_COMM_WORLD);
  }
  else if(tnderef>nlbl) {
    MPI_Allgatherv(MPI_IN_PLACE, bnderef[Globals::my_rank], MPI_BYTE,
                    lderef, bnderef, bddisp, MPI_BYTE, MPI_COMM_WORLD);
  }
#endif

  // calculate the list of the newly derefined blocks
  int ctnd=0;
  if(tnderef>nlbl) {
    int lk=0, lj=0;
    if(mesh_size.nx2 > 1) lj=1;
    if(mesh_size.nx3 > 1) lk=1;
    for(int n=0; n<tnderef; n++) {
      if((lderef[n].lx1&1L)==0 && (lderef[n].lx2&1L)==0 && (lderef[n].lx3&1L)==0) {
        int r=n, rr=0;
        for(long int k=0;k<=lk;k++) {
          for(long int j=0;j<=lj;j++) {
            for(long int i=0;i<=1;i++) {
              if((lderef[n].lx1+i)==lderef[r].lx1
              && (lderef[n].lx2+j)==lderef[r].lx2
              && (lderef[n].lx3+k)==lderef[r].lx3
              &&  lderef[n].level ==lderef[r].level)
                rr++;
              r++;
            }
          }
        }
        if(rr==nlbl) {
          clderef[ctnd].lx1  =(lderef[n].lx1>>1);
          clderef[ctnd].lx2  =(lderef[n].lx2>>1);
          clderef[ctnd].lx3  =(lderef[n].lx3>>1);
          clderef[ctnd].level=lderef[n].level-1;
          ctnd++;
        }
      }
    }
  }
  // sort the lists by level
  if(ctnd>1)
    std::sort(clderef, &(clderef[ctnd-1]), LogicalLocation::Greater);

/*  if(Globals::my_rank==0) {
    for(int n=0; n<tnref; n++)
      std::cout << "Refine flag:" << lref[n].lx1 << " " << lref[n].lx2 << " " << lref[n].lx3 << " " << lref[n].level << std::endl;
    for(int n=0; n<ctnd; n++)
      std::cout << "Derefine flag:" << clderef[n].lx1 << " " << clderef[n].lx2 << " " << clderef[n].lx3 << " " << clderef[n].level << std::endl;
  }*/

  if(tnderef>nlbl)
    delete [] lderef;

  // Now the lists of the blocks to be refined and derefined are completed
  // Start tree manipulation
  // Step 1. perform refinement
  int nnew=0, ndel=0, ntot=0;
  for(int n=0; n<tnref; n++) {
    MeshBlockTree *bt=tree.FindMeshBlock(lref[n]);
    bt->Refine(tree, dim, mesh_bcs, nrbx1, nrbx2, nrbx3, root_level, nnew);
  }
  if(tnref!=0)
    delete [] lref;
  // Step 2. perform derefinement
  for(int n=0; n<ctnd; n++) {
    MeshBlockTree *bt=tree.FindMeshBlock(clderef[n]);
    bt->Derefine(tree, dim, mesh_bcs, nrbx1, nrbx2, nrbx3, root_level, ndel);
  }
  if(tnderef>nlbl)
    delete [] clderef;
  ntot=nbtotal+nnew-ndel;
  if(nnew==0 && ndel==0)
    return; // nothing to do
  // Tree manipulation completed

  // Block exchange
  // Step 1. construct new lists
  LogicalLocation *newloc = new LogicalLocation[ntot];
  int *newrank = new int[ntot];
  Real *newcost = new Real[ntot];
  int *newtoold = new int[ntot];
  int *oldtonew = new int[nbtotal];
  tree.GetMeshBlockList(newloc,newtoold,nbtotal);
  // create a list mapping the previous gid to the current one
  oldtonew[0]=0;
  int k=1;
  for(int n=1; n<ntot; n++) {
    if(newtoold[n]==newtoold[n-1]+1) { // normal
      oldtonew[k++]=n;
    }
    else if(newtoold[n]==newtoold[n-1]+nlbl) { // derefined
      for(int j=0; j<nlbl-1; j++)
        oldtonew[k++]=n-1;
      oldtonew[k++]=n;
    }
  }

#ifdef MPI_PARALLEL
  // receive the old cost list
  MPI_Wait(&areq[2], MPI_STATUS_IGNORE);
#endif
  for(int n=0; n<ntot; n++) {
    int on=newtoold[n];
    if(newloc[n].level>=loclist[on].level) // same or refined
      newcost[n]=costlist[on];
    else {
      Real acost=0.0;
      for(int l=0; l<nlbl; l++)
        acost+=costlist[on+l];
      newcost[n]=acost/nlbl;
    }
  }

  // store old nbstart and nbend
  int onbs=nslist[Globals::my_rank];
  int onbe=onbs+nblist[Globals::my_rank]-1;

  // Step 2. Calculate new load balance
  LoadBalancing(newcost, newrank, nslist, nblist, ntot);

  int nbs=nslist[Globals::my_rank];
  int nbe=nbs+nblist[Globals::my_rank]-1;

  int f2, f3;
  int &bnx1=pblock->block_size.nx1;
  int &bnx2=pblock->block_size.nx2;
  int &bnx3=pblock->block_size.nx3;
  if(mesh_size.nx2>1) f2=1;
  else f2=0;
  if(mesh_size.nx3>1) f3=1;
  else f3=0;

#ifdef MPI_PARALLEL
  // Step 3. count the number of the blocks to be sent / received
  int nsend=0, nrecv=0;
  for(int n=nbs; n<=nbe; n++) { 
    int on=newtoold[n];
    if(loclist[on].level > newloc[n].level) { // f2c
      for(int k=0; k<nlbl; k++) {
        if(ranklist[on+k]!=Globals::my_rank)
          nrecv++;
      }
    }
    else {
      if(ranklist[on]!=Globals::my_rank)
        nrecv++;
    }
  }
  for(int n=onbs; n<=onbe; n++) { 
    int nn=oldtonew[n];
    if(loclist[n].level < newloc[nn].level) { // c2f
      for(int k=0; k<nlbl; k++) {
        if(newrank[nn+k]!=Globals::my_rank)
          nsend++;
      }
    }
    else {
      if(newrank[nn]!=Globals::my_rank)
        nsend++;
    }
  }

  // Step 4. calculate buffer sizes
  Real **sendbuf, **recvbuf;
  int bssame=bnx1*bnx2*bnx3*NHYDRO;
  int bsf2c=(bnx1/2)*((bnx2+1)/2)*((bnx3+1)/2)*NHYDRO;
  int bsc2f=(bnx1/2+2)*((bnx2+1)/2+2*f2)*((bnx3+1)/2+2*f3)*NHYDRO;
  if(MAGNETIC_FIELDS_ENABLED) {
    bssame+=(bnx1+1)*bnx2*bnx3+bnx1*(bnx2+f2)*bnx3+bnx1*bnx2*(bnx3+f3);
    bsf2c+=((bnx1/2)+1)*((bnx2+1)/2)*((bnx3+1)/2)
          +(bnx1/2)*(((bnx2+1)/2)+f2)*((bnx3+1)/2)
          +(bnx1/2)*((bnx2+1)/2)*(((bnx3+1)/2)+f3);
    bsc2f+=((bnx1/2)+1+2)*((bnx2+1)/2+2*f2)*((bnx3+1)/2+2*f3)
          +(bnx1/2+2)*(((bnx2+1)/2)+f2+2*f2)*((bnx3+1)/2+2*f3)
          +(bnx1/2+2)*((bnx2+1)/2+2*f2)*(((bnx3+1)/2)+f3+2*f3);
  }

  MPI_Request *req_send, *req_recv;
  // Step 5. allocate and start receiving buffers
  if(nrecv!=0) {
    recvbuf = new Real*[nrecv];
    req_recv = new MPI_Request[nrecv];
    int k=0;
    for(int n=nbs; n<=nbe; n++) { 
      int on=newtoold[n];
      LogicalLocation &oloc=loclist[on];
      LogicalLocation &nloc=newloc[n];
      if(oloc.level>nloc.level) { // f2c
        for(int l=0; l<nlbl; l++) {
          if(ranklist[on+l]==Globals::my_rank) continue;
          LogicalLocation &lloc=loclist[on+l];
          int ox1=lloc.lx1&1L, ox2=lloc.lx2&1L, ox3=lloc.lx3&1L;
          recvbuf[k] = new Real[bsf2c];
          int tag=CreateAMRMPITag(n-nbs, ox1, ox2, ox3);
          MPI_Irecv(recvbuf[k], bsf2c, MPI_ATHENA_REAL, ranklist[on+l],
                    tag, MPI_COMM_WORLD, &(req_recv[k]));
          k++;
        }
      }
      else { // same or c2f
        if(ranklist[on]==Globals::my_rank) continue;
        int size;
        if(oloc.level == nloc.level) size=bssame;
        else size=bsc2f;
        recvbuf[k] = new Real[size];
        int tag=CreateAMRMPITag(n-nbs, 0, 0, 0);
        MPI_Irecv(recvbuf[k], size, MPI_ATHENA_REAL, ranklist[on],
                  tag, MPI_COMM_WORLD, &(req_recv[k]));
        k++;
      }
    }
  }
  // Step 6. allocate, pack and start sending buffers
  if(nsend!=0) {
    sendbuf = new Real*[nsend];
    req_send = new MPI_Request[nsend];
    int k=0;
    for(int n=onbs; n<=onbe; n++) { 
      int nn=oldtonew[n];
      LogicalLocation &oloc=loclist[n];
      LogicalLocation &nloc=newloc[nn];
      MeshBlock* pb=FindMeshBlock(n);
      if(nloc.level==oloc.level) { // same
        if(newrank[nn]==Globals::my_rank) continue;
        sendbuf[k] = new Real[bssame];
        // pack
        int p=0;
        BufferUtility::Pack4DData(pb->phydro->u, sendbuf[k], 0, NHYDRO-1,
                       pb->is, pb->ie, pb->js, pb->je, pb->ks, pb->ke, p);
        if(MAGNETIC_FIELDS_ENABLED) {
          BufferUtility::Pack3DData(pb->pfield->b.x1f, sendbuf[k],
                         pb->is, pb->ie+1, pb->js, pb->je, pb->ks, pb->ke, p);
          BufferUtility::Pack3DData(pb->pfield->b.x2f, sendbuf[k],
                         pb->is, pb->ie+1, pb->js, pb->je, pb->ks, pb->ke, p);
          BufferUtility::Pack3DData(pb->pfield->b.x3f, sendbuf[k],
                         pb->is, pb->ie+1, pb->js, pb->je, pb->ks, pb->ke, p);
        }
        int tag=CreateAMRMPITag(nn-nslist[newrank[nn]], 0, 0, 0);
        MPI_Isend(sendbuf[k], bssame, MPI_ATHENA_REAL, newrank[nn],
                  tag, MPI_COMM_WORLD, &(req_send[k]));
        k++;
      }
      else if(nloc.level>oloc.level) { // c2f
        for(int l=0; l<nlbl; l++) {
          if(newrank[nn+l]==Globals::my_rank) continue;
          LogicalLocation &lloc=newloc[nn+l];
          int ox1=lloc.lx1&1L, ox2=lloc.lx2&1L, ox3=lloc.lx3&1L;
          sendbuf[k] = new Real[bsc2f];
          // pack
          int is, ie, js, je, ks, ke;
          if(ox1==0) is=pb->is-1,                       ie=pb->is+pb->block_size.nx1/2;
          else       is=pb->is+pb->block_size.nx1/2-1,  ie=pb->ie+1;
          if(ox2==0) js=pb->js-f2,                      je=pb->js+pb->block_size.nx2/2;
          if(ox2==0) js=pb->js-f2,                      je=pb->js+pb->block_size.nx2/2;
          else       js=pb->js+pb->block_size.nx2/2-f2, je=pb->je+f2;
          if(ox3==0) ks=pb->ks-f3,                      ke=pb->ks+pb->block_size.nx3/2;
          else       ks=pb->ks+pb->block_size.nx3/2-f3, ke=pb->ke+f3;
          int p=0;
          BufferUtility::Pack4DData(pb->phydro->u, sendbuf[k], 0, NHYDRO-1,
                                    is, ie, js, je, ks, ke, p);
          if(MAGNETIC_FIELDS_ENABLED) {
            BufferUtility::Pack3DData(pb->pfield->b.x1f, sendbuf[k],
                                      is, ie+1, js, je, ks, ke, p);
            BufferUtility::Pack3DData(pb->pfield->b.x2f, sendbuf[k],
                                      is, ie, js, je+f2, ks, ke, p);
            BufferUtility::Pack3DData(pb->pfield->b.x3f, sendbuf[k],
                                      is, ie, js, je, ks, ke+f3, p);
          }
          int tag=CreateAMRMPITag(nn+l-nslist[newrank[nn+l]], 0, 0, 0);
          MPI_Isend(sendbuf[k], bsc2f, MPI_ATHENA_REAL, newrank[nn+l],
                    tag, MPI_COMM_WORLD, &(req_send[k]));
          k++;
        }
      }
      else { // f2c
        if(newrank[nn]==Globals::my_rank) continue;
        int ox1=oloc.lx1&1L, ox2=oloc.lx2&1L, ox3=oloc.lx3&1L;
        sendbuf[k] = new Real[bsf2c];
        // restrict and pack
        MeshRefinement *pmr=pb->pmr;
        pmr->RestrictCellCenteredValues(pb->phydro->u, pmr->coarse_cons_,
             0, NHYDRO-1, pb->cis, pb->cie, pb->cjs, pb->cje, pb->cks, pb->cke);
        int p=0;
        BufferUtility::Pack4DData(pmr->coarse_cons_, sendbuf[k], 0, NHYDRO-1,
                       pb->cis, pb->cie, pb->cjs, pb->cje, pb->cks, pb->cke, p);
        if(MAGNETIC_FIELDS_ENABLED) {
          pmr->RestrictFieldX1(pb->pfield->b.x1f, pmr->coarse_b_.x1f,
                               pb->cis, pb->cie+1, pb->cjs, pb->cje, pb->cks, pb->cke);
          BufferUtility::Pack3DData(pmr->coarse_b_.x1f, sendbuf[k],
                         pb->cis, pb->cie+1, pb->cjs, pb->cje, pb->cks, pb->cke, p);
          pmr->RestrictFieldX2(pb->pfield->b.x2f, pmr->coarse_b_.x2f,
                               pb->cis, pb->cie, pb->cjs, pb->cje+f2, pb->cks, pb->cke);
          BufferUtility::Pack3DData(pmr->coarse_b_.x2f, sendbuf[k],
                         pb->cis, pb->cie, pb->cjs, pb->cje+f2, pb->cks, pb->cke, p);
          pmr->RestrictFieldX3(pb->pfield->b.x3f, pmr->coarse_b_.x3f,
                               pb->cis, pb->cie, pb->cjs, pb->cje, pb->cks, pb->cke+f3);
          BufferUtility::Pack3DData(pmr->coarse_b_.x3f, sendbuf[k],
                         pb->cis, pb->cie, pb->cjs, pb->cje, pb->cks, pb->cke+f3, p);
        }
        int tag=CreateAMRMPITag(nn-nslist[newrank[nn]], ox1, ox2, ox3);
        MPI_Isend(sendbuf[k], bsf2c, MPI_ATHENA_REAL, newrank[nn],
                  tag, MPI_COMM_WORLD, &(req_send[k]));
        k++;
      }
    }
  }
#endif

  // Step 7. construct a new MeshBlock list
  // move the data within the node
  MeshBlock *newlist=NULL;
  RegionSize block_size=pblock->block_size;
  enum BoundaryFlag block_bcs[6];

  for(int n=nbs; n<=nbe; n++) {
    int on=newtoold[n];
    if((ranklist[on]==Globals::my_rank) && (loclist[on].level == newloc[n].level)) {
      // on the same node and same level -> just move it
      MeshBlock* pob=FindMeshBlock(on);
      if(pob->prev==NULL) pblock=pob->next;
      else pob->prev->next=pob->next;
      if(pob->next!=NULL) pob->next->prev=pob->prev;
      pob->next=NULL;
      if(n==nbs) { // first
        pob->prev=NULL;
        newlist=pmb=pob;
      }
      else {
        pmb->next=pob;
        pob->prev=pmb;
        pmb=pmb->next;
      }
      pmb->gid=n; pmb->lid=n-nbs;
    }
    else {
      // on a different level or node - create a new block
      SetBlockSizeAndBoundaries(newloc[n], block_size, block_bcs);
      if(n==nbs) { // first
        newlist = new MeshBlock(n, n-nbs, newloc[n], block_size, block_bcs, this, pin);
        pmb=newlist;
      }
      else {
        pmb->next = new MeshBlock(n, n-nbs, newloc[n], block_size, block_bcs, this, pin);
        pmb->next->prev=pmb;
        pmb=pmb->next;
      }
      // fill the conservative variables
      if((loclist[on].level>newloc[n].level)) { // fine to coarse
        for(int ll=0; ll<nlbl; ll++) {
          if(ranklist[on+ll]!=Globals::my_rank) continue;
          // on the same node - restriction
          MeshBlock* pob=FindMeshBlock(on+ll);
          MeshRefinement *pmr=pob->pmr;
          pmr->RestrictCellCenteredValues(pob->phydro->u, pmr->coarse_cons_,
               0, NHYDRO-1, pob->cis, pob->cie, pob->cjs, pob->cje, pob->cks, pob->cke);
          int is=pmb->is+(loclist[on+ll].lx1&1L)*pmb->block_size.nx1/2;
          int js=pmb->js+(loclist[on+ll].lx2&1L)*pmb->block_size.nx2/2;
          int ks=pmb->ks+(loclist[on+ll].lx3&1L)*pmb->block_size.nx3/2;
          AthenaArray<Real> &src=pmr->coarse_cons_;
          AthenaArray<Real> &dst=pmb->phydro->u;
          for(int nn=0; nn<NHYDRO; nn++) {
            for(int k=ks, fk=pob->cks; fk<=pob->cke; k++, fk++) {
              for(int j=js, fj=pob->cjs; fj<=pob->cje; j++, fj++) {
                for(int i=is, fi=pob->cis; fi<=pob->cie; i++, fi++)
                  dst(nn, k, j, i)=src(nn, fk, fj, fi);
          }}}
          if(MAGNETIC_FIELDS_ENABLED) {
            pmr->RestrictFieldX1(pob->pfield->b.x1f, pmr->coarse_b_.x1f,
                         pob->cis, pob->cie+1, pob->cjs, pob->cje, pob->cks, pob->cke);
            pmr->RestrictFieldX2(pob->pfield->b.x2f, pmr->coarse_b_.x2f,
                         pob->cis, pob->cie, pob->cjs, pob->cje+f2, pob->cks, pob->cke);
            pmr->RestrictFieldX3(pob->pfield->b.x3f, pmr->coarse_b_.x3f,
                         pob->cis, pob->cie, pob->cjs, pob->cje, pob->cks, pob->cke+f3);
            FaceField &src=pmr->coarse_b_;
            FaceField &dst=pmb->pfield->b;
            for(int k=ks, fk=pob->cks; fk<=pob->cke; k++, fk++) {
              for(int j=js, fj=pob->cjs; fj<=pob->cje; j++, fj++) {
                for(int i=is, fi=pob->cis; fi<=pob->cie+1; i++, fi++)
                  dst.x1f(k, j, i)=src.x1f(fk, fj, fi);
            }}
            for(int k=ks, fk=pob->cks; fk<=pob->cke; k++, fk++) {
              for(int j=js, fj=pob->cjs; fj<=pob->cje+f2; j++, fj++) {
                for(int i=is, fi=pob->cis; fi<=pob->cie; i++, fi++)
                  dst.x2f(k, j, i)=src.x2f(fk, fj, fi);
            }}
            if(pmb->block_size.nx2==1) {
              int ie=is+block_size.nx1/2-1;
              for(int i=is; i<=ie; i++)
                dst.x2f(pmb->ks, pmb->js+1, i)=dst.x2f(pmb->ks, pmb->js, i);
            }
            for(int k=ks, fk=pob->cks; fk<=pob->cke+f3; k++, fk++) {
              for(int j=js, fj=pob->cjs; fj<=pob->cje; j++, fj++) {
                for(int i=is, fi=pob->cis; fi<=pob->cie; i++, fi++)
                  dst.x3f(k, j, i)=src.x3f(fk, fj, fi);
            }}
            if(pmb->block_size.nx3==1) {
              int ie=is+block_size.nx1/2-1, je=js+block_size.nx2/2-1;
              for(int j=js; j<=je; j++) {
                for(int i=is; i<=ie; i++)
                  dst.x3f(pmb->ks+1, j, i)=dst.x3f(pmb->ks, j, i);
              }
            }
          }
        }
      }
      else if((loclist[on].level < newloc[n].level) && (ranklist[on]==Globals::my_rank)) {
        // coarse to fine on the same node - prolongation
        if(ranklist[on]!=Globals::my_rank) continue;
        MeshBlock* pob=FindMeshBlock(on);
        MeshRefinement *pmr=pmb->pmr;
        int is=pob->cis-1, ie=pob->cie+1, js=pob->cjs-f2,
            je=pob->cje+f2, ks=pob->cks-f3, ke=pob->cke+f3;
        int cis=(newloc[n].lx1&1L)*pob->block_size.nx1/2+pob->is-1;
        int cjs=(newloc[n].lx2&1L)*pob->block_size.nx2/2+pob->js-f2;
        int cks=(newloc[n].lx3&1L)*pob->block_size.nx3/2+pob->ks-f3;
        AthenaArray<Real> &src=pob->phydro->u;
        AthenaArray<Real> &dst=pmr->coarse_cons_;
        // fill the coarse buffer
        for(int nn=0; nn<NHYDRO; nn++) {
          for(int k=ks, ck=cks; k<=ke; k++, ck++) {
            for(int j=js, cj=cjs; j<=je; j++, cj++) {
              for(int i=is, ci=cis; i<=ie; i++, ci++)
                dst(nn, k, j, i)=src(nn, ck, cj, ci);
        }}}
        pmr->ProlongateCellCenteredValues(dst, pmb->phydro->u, 0, NHYDRO-1,
                                          is, ie, js, je, ks, ke);
        if(MAGNETIC_FIELDS_ENABLED) {
          FaceField &src=pob->pfield->b;
          FaceField &dst=pmr->coarse_b_;
          for(int k=ks, ck=cks; k<=ke; k++, ck++) {
            for(int j=js, cj=cjs; j<=je; j++, cj++) {
              for(int i=is, ci=cis; i<=ie+1; i++, ci++)
                dst.x1f(k, j, i)=src.x1f(ck, cj, ci);
          }}
          for(int k=ks, ck=cks; k<=ke; k++, ck++) {
            for(int j=js, cj=cjs; j<=je+f2; j++, cj++) {
              for(int i=is, ci=cis; i<=ie; i++, ci++)
                dst.x2f(k, j, i)=src.x2f(ck, cj, ci);
          }}
          for(int k=ks, ck=cks; k<=ke+f3; k++, ck++) {
            for(int j=js, cj=cjs; j<=je; j++, cj++) {
              for(int i=is, ci=cis; i<=ie; i++, ci++)
                dst.x3f(k, j, i)=src.x3f(ck, cj, ci);
          }}
          pmr->ProlongateSharedFieldX1(dst.x1f, pmb->pfield->b.x1f,
                                       pob->is, ie+1, js, je, ks, ke);
          pmr->ProlongateSharedFieldX2(dst.x2f, pmb->pfield->b.x2f,
                                       is, ie, js, je+f2, ks, ke);
          pmr->ProlongateSharedFieldX3(dst.x3f, pmb->pfield->b.x3f,
                                       is, ie, js, je, ks, ke+f3);
          pmr->ProlongateInternalField(pmb->pfield->b, is, ie, js, je, ks, ke);
        }
      }
    }
  }

  // discard remaining MeshBlocks
  // they could be reused, but for the moment, just throw them away for simplicity
  if(pblock!=NULL) {
    while(pblock->next != NULL)
      delete pblock->next;
    delete pblock;
  }

  // Replace the MeshBlock list
  pblock=newlist;

  // Step 8. Receive the data and load into MeshBlocks
  // This is a test: try MPI_Waitall later.
#ifdef MPI_PARALLEL
  if(nrecv!=0) {
    int k=0;
    for(int n=nbs; n<=nbe; n++) { 
      int on=newtoold[n];
      LogicalLocation &oloc=loclist[on];
      LogicalLocation &nloc=newloc[n];
      MeshBlock *pb=FindMeshBlock(n);
      if(oloc.level==nloc.level) { // same
        if(ranklist[on]==Globals::my_rank) continue;
        MPI_Wait(&(req_recv[k]), MPI_STATUS_IGNORE);
        int p=0;
        BufferUtility::Unpack4DData(recvbuf[k], pb->phydro->u, 0, NHYDRO-1,
                       pb->is, pb->ie, pb->js, pb->je, pb->ks, pb->ke, p);
        if(MAGNETIC_FIELDS_ENABLED) {
          FaceField &dst=pb->pfield->b;
          BufferUtility::Unpack3DData(recvbuf[k], dst.x1f,
                         pb->is, pb->ie+1, pb->js, pb->je, pb->ks, pb->ke, p);
          BufferUtility::Unpack3DData(recvbuf[k], dst.x2f,
                         pb->is, pb->ie+1, pb->js, pb->je, pb->ks, pb->ke, p);
          BufferUtility::Unpack3DData(recvbuf[k], dst.x3f,
                         pb->is, pb->ie+1, pb->js, pb->je, pb->ks, pb->ke, p);
          if(pb->block_size.nx2==1) {
            for(int i=pb->is; i<=pb->ie; i++)
              dst.x2f(pb->ks, pb->js+1, i)=dst.x2f(pb->ks, pb->js, i);
          }
          if(pb->block_size.nx3==1) {
            for(int j=pb->js; j<=pb->je; j++) {
              for(int i=pb->is; i<=pb->ie; i++)
                dst.x3f(pb->ks+1, j, i)=dst.x3f(pb->ks, j, i);
            }
          }
        }
        k++;
      }
      else if(oloc.level>nloc.level) { // f2c
        for(int l=0; l<nlbl; l++) {
          if(ranklist[on+l]==Globals::my_rank) continue;
          LogicalLocation &lloc=loclist[on+l];
          int ox1=lloc.lx1&1L, ox2=lloc.lx2&1L, ox3=lloc.lx3&1L;
          int p=0, is, ie, js, je, ks, ke;
          if(ox1==0) is=pb->is,                      ie=pb->is+pb->block_size.nx1/2-1;
          else       is=pb->is+pb->block_size.nx1/2, ie=pb->ie;
          if(ox2==0) js=pb->js,                      je=pb->js+pb->block_size.nx2/2-f2;
          else       js=pb->js+pb->block_size.nx2/2, je=pb->je;
          if(ox3==0) ks=pb->ks,                      ke=pb->ks+pb->block_size.nx3/2-f3;
          else       ks=pb->ks+pb->block_size.nx3/2, ke=pb->ke;
          MPI_Wait(&(req_recv[k]), MPI_STATUS_IGNORE);
          BufferUtility::Unpack4DData(recvbuf[k], pb->phydro->u, 0, NHYDRO-1,
                         is, ie, js, je, ks, ke, p);
          if(MAGNETIC_FIELDS_ENABLED) {
            FaceField &dst=pb->pfield->b;
            BufferUtility::Unpack3DData(recvbuf[k], dst.x1f,
                           is, ie+1, js, je, ks, ke, p);
            BufferUtility::Unpack3DData(recvbuf[k], dst.x2f,
                           is, ie+1, js, je, ks, ke, p);
            BufferUtility::Unpack3DData(recvbuf[k], dst.x3f,
                           is, ie+1, js, je, ks, ke, p);
            if(pb->block_size.nx2==1) {
              for(int i=is; i<=ie; i++)
                dst.x2f(pb->ks, pb->js+1, i)=dst.x2f(pb->ks, pb->js, i);
            }
            if(pb->block_size.nx3==1) {
              for(int j=js; j<=je; j++) {
                for(int i=is; i<=ie; i++)
                  dst.x3f(pb->ks+1, j, i)=dst.x3f(pb->ks, j, i);
              }
            }
          }
          k++;
        }
      }
      else { // c2f
        if(ranklist[on]==Globals::my_rank) continue;
        MeshRefinement *pmr=pb->pmr;
        int p=0;
        int is=pb->cis-1, ie=pb->cie+1, js=pb->cjs-f2,
            je=pb->cje+f2, ks=pb->cks-f3, ke=pb->cke+f3;
        MPI_Wait(&(req_recv[k]), MPI_STATUS_IGNORE);
        BufferUtility::Unpack4DData(recvbuf[k], pmr->coarse_cons_,
                                    0, NHYDRO-1, is, ie, js, je, ks, ke, p);
        pmr->ProlongateCellCenteredValues(pmr->coarse_cons_, pb->phydro->u, 0, NHYDRO-1,
                                          is, ie, js, je, ks, ke);
        if(MAGNETIC_FIELDS_ENABLED) {
          BufferUtility::Unpack3DData(recvbuf[k], pmr->coarse_b_.x1f,
                                      is, ie+1, js, je, ks, ke, p);
          BufferUtility::Unpack3DData(recvbuf[k], pmr->coarse_b_.x2f,
                                      is, ie, js, je+f2, ks, ke, p);
          BufferUtility::Unpack3DData(recvbuf[k], pmr->coarse_b_.x3f,
                                      is, ie, js, je, ks, ke+f3, p);
          pmr->ProlongateSharedFieldX1(pmr->coarse_b_.x1f, pb->pfield->b.x1f,
                                       is, ie+1, js, je, ks, ke);
          pmr->ProlongateSharedFieldX2(pmr->coarse_b_.x2f, pb->pfield->b.x2f,
                                       is, ie, js, je+f2, ks, ke);
          pmr->ProlongateSharedFieldX3(pmr->coarse_b_.x3f, pb->pfield->b.x3f,
                                       is, ie, js, je, ks, ke+f3);
          pmr->ProlongateInternalField(pb->pfield->b, is, ie, js, je, ks, ke);
        }
        k++;
      }
    }
  }
#endif

  // deallocate arrays
  delete [] loclist;
  delete [] ranklist;
  delete [] costlist;
  delete [] newtoold;
  delete [] oldtonew;
#ifdef MPI_PARALLEL
  if(nsend!=0) {
    MPI_Waitall(nsend, req_send, MPI_STATUSES_IGNORE);
    for(int n=0;n<nsend;n++)
      delete [] sendbuf[n];
    delete [] sendbuf;
    delete [] req_send;
  }
  if(nrecv!=0) {
    for(int n=0;n<nrecv;n++)
      delete [] recvbuf[n];
    delete [] recvbuf;
    delete [] req_recv;
  }
#endif

  // update the lists
  loclist = newloc;
  ranklist = newrank;
  costlist = newcost;

  // re-initialize the MeshBlocks
  pmb=pblock;
  while(pmb!=NULL) {
    pmb->SearchAndSetNeighbors(tree, ranklist, nslist);
    pmb=pmb->next;
  }
  Initialize(2, pin);

  return;
}

