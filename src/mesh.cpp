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

// Primary header
#include "mesh.hpp"

// C++ headers
#include <cfloat>     // FLT_MAX
#include <cmath>      // std::abs(), pow()
#include <iostream>   // cout, endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <algorithm>  // sort, find
#include <iomanip>

#include <stdlib.h>
#include <string.h>  // memcpy

// Athena headers
#include "athena.hpp"                   // enums, macros, Real
#include "athena_arrays.hpp"            // AthenaArray
#include "coordinates/coordinates.hpp"  // Coordinates
#include "fluid/fluid.hpp"              // Fluid
#include "field/field.hpp"              // Field
#include "bvals/bvals.hpp"              // BoundaryValues
#include "fluid/eos/eos.hpp"              // FluidEqnOfState
#include "fluid/integrators/fluid_integrator.hpp"  // FluidIntegrator
#include "field/integrators/field_integrator.hpp"  // FieldIntegrator
#include "parameter_input.hpp"          // ParameterInput
#include "blockuid.hpp"        // BlockUID
#include "wrapio.hpp"
#include "tasklist.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif


//======================================================================================
//! \file mesh.cpp
//  \brief implementation of functions in classes Mesh, and MeshBlock
//======================================================================================

//--------------------------------------------------------------------------------------
// Mesh constructor, builds mesh at start of calculation using parameters in input file

Mesh::Mesh(ParameterInput *pin, int test_flag)
{
  std::stringstream msg;
  RegionSize block_size;
  BlockTree *neibt;
  MeshBlock *pfirst;
  int block_bcs[6];
  Real totalcost, maxcost, mincost, mycost, targetcost;
  int nbmax, dim;
  long int lx1, lx2, lx3;
  int ll, i, j;

// mesh test
  if(test_flag>0) nproc=test_flag;

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

// read BC flags for each of the 6 boundaries in turn.  Error tests performed in
// BoundaryValues constructor

  mesh_bcs[inner_x1] = pin->GetOrAddInteger("mesh","ix1_bc",0);
  mesh_bcs[outer_x1] = pin->GetOrAddInteger("mesh","ox1_bc",0);
  mesh_bcs[inner_x2] = pin->GetOrAddInteger("mesh","ix2_bc",0);
  mesh_bcs[outer_x2] = pin->GetOrAddInteger("mesh","ox2_bc",0);
  mesh_bcs[inner_x3] = pin->GetOrAddInteger("mesh","ix3_bc",0);
  mesh_bcs[outer_x3] = pin->GetOrAddInteger("mesh","ox3_bc",0);


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

  if(myrank==0)
    std::cout << "RootGrid = " << nrbx1 << " x " << nrbx2
              << " x " << nrbx3 << std::endl;

// calculate the logical root level and maximum level
  for(root_level=0;(1<<root_level)<nbmax;root_level++);

// create Block UID list
  tree.CreateRootGrid(nrbx1,nrbx2,nrbx3,root_level);

// SMR / AMR: create finer grids here
  multilevel=false;
  adaptive=false;
  if(pin->GetOrAddString("mesh","refinement","static")=="adaptive")
    adaptive=true, multilevel=true;
  if(adaptive==true)
    max_level = pin->GetOrAddInteger("mesh","maxlevel",1)+root_level-1;
  else
    max_level = IDLENGTH*20;

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
        if(MeshGeneratorX1((Real)(lx1min+1)/lxmax,mesh_size)>ref_size.x1min)
          break;
      }
      for(lx1max=lx1min;lx1max<lxmax;lx1max++) {
        if(MeshGeneratorX1((Real)(lx1max+1)/lxmax,mesh_size)>=ref_size.x1max)
          break;
      }
      if(lx1min%2==1) lx1min--;
      if(lx1max%2==0) lx1max++;
      if(dim>=2) { // 2D or 3D
        lxmax=nrbx2*(1L<<ref_lev);
        for(lx2min=0;lx2min<lxmax;lx2min++) {
          if(MeshGeneratorX2((Real)(lx2min+1)/lxmax,mesh_size)>ref_size.x2min)
            break;
        }
        for(lx2max=lx2min;lx2max<lxmax;lx2max++) {
          if(MeshGeneratorX2((Real)(lx2max+1)/lxmax,mesh_size)>=ref_size.x2max)
            break;
        }
        if(lx2min%2==1) lx2min--;
        if(lx2max%2==0) lx2max++;
      }
      if(dim==3) { // 3D
        lxmax=nrbx3*(1L<<ref_lev);
        for(lx3min=0;lx3min<lxmax;lx3min++) {
          if(MeshGeneratorX3((Real)(lx3min+1)/lxmax,mesh_size)>ref_size.x3min)
            break;
        }
        for(lx3max=lx3min;lx3max<lxmax;lx3max++) {
          if(MeshGeneratorX3((Real)(lx3max+1)/lxmax,mesh_size)>=ref_size.x3max)
            break;
        }
        if(lx3min%2==1) lx3min--;
        if(lx3max%2==0) lx3max++;
      }
      // create the finest level
      std::cout << "refinenment: logical level = " << lrlev << ", lx1min = "
                << lx1min << ", lx1max = " << lx1max;
      if(dim==1) {
        for(long int i=lx1min; i<lx1max; i+=2) {
          BlockUID tid;
          tid.CreateUIDfromLocation(i,0,0,lrlev);
          tree.AddMeshBlock(tree,tid,dim,mesh_bcs,nrbx1,nrbx2,nrbx3,root_level);
        }
      }
      if(dim==2) {
        std::cout << ", lx2min = " << lx2min << ", lx2max = " << lx2max;
        for(long int j=lx2min; j<lx2max; j+=2) {
          for(long int i=lx1min; i<lx1max; i+=2) {
            BlockUID tid;
            tid.CreateUIDfromLocation(i,j,0,lrlev);
            tree.AddMeshBlock(tree,tid,dim,mesh_bcs,nrbx1,nrbx2,nrbx3,root_level);
          }
        }
      }
      if(dim==3) {
        std::cout << ", lx3min = " << lx3min << ", lx3max = " << lx3max;
        for(long int k=lx3min; k<lx3max; k+=2) {
          for(long int j=lx2min; j<lx2max; j+=2) {
            for(long int i=lx1min; i<lx1max; i+=2) {
              BlockUID tid;
              tid.CreateUIDfromLocation(i,j,k,lrlev);
              tree.AddMeshBlock(tree,tid,dim,mesh_bcs,nrbx1,nrbx2,nrbx3,root_level);
            }
          }
        }
      }
      std::cout << std::endl;
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
  if (MAGNETIC_FIELDS_ENABLED || multilevel==true)
    face_only=false;

  maxneighbor_=BufferID(dim, multilevel, face_only);

  // initial mesh hierarchy construction is completed here

  tree.AssignGID(nbtotal); // count blocks at the same time
  buid=new BlockUID[nbtotal];
  tree.GetIDList(buid,nbtotal);

  ranklist=new int[nbtotal];
  nslist=new int[nproc];
  nblist=new int[nproc];
  costlist=new Real[nbtotal];
  maxcost=0.0;
  mincost=(FLT_MAX);
  totalcost=0.0;
  for(i=0;i<nbtotal;i++) {
    costlist[i]=1.0; // the simplest estimate; all the blocks are equal
    totalcost+=costlist[i];
    mincost=std::min(mincost,costlist[i]);
    maxcost=std::max(maxcost,costlist[i]);
  }
  j=nproc-1;
  targetcost=totalcost/nproc;
  mycost=0.0;
  // create rank list from the end: the master node should have less load
  for(i=nbtotal-1;i>=0;i--) {
    mycost+=costlist[i];
    ranklist[i]=j;
    if(mycost >= targetcost && j>0) {
      j--;
      totalcost-=mycost;
      mycost=0.0;
      targetcost=totalcost/(j+1);
    }
  }
  nslist[0]=0;
  j=0;
  for(i=1;i<nbtotal;i++) { // make the list of nbstart
    if(ranklist[i]!=ranklist[i-1]) {
      nblist[j]=i-nslist[j];
      nslist[++j]=i;
    }
  }
  nblist[j]=nbtotal-nslist[j];

  // store my nbstart and nbend
  nbstart=nslist[myrank];
  if(myrank+1==nproc)
    nbend=nbtotal-1;
  else 
    nbend=nslist[myrank+1]-1;

// check if there are sufficient blocks
#ifdef MPI_PARALLEL
  if(nbtotal < nproc) {
    if(test_flag==0) {
      msg << "### FATAL ERROR in Mesh constructor" << std::endl
          << "Too few blocks: nbtotal (" << nbtotal << ") < nproc ("<< nproc
          << ")" << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    else { // test
      std::cout << "### Warning in Mesh constructor" << std::endl
          << "Too few blocks: nbtotal (" << nbtotal << ") < nproc ("<< nproc
          << ")" << std::endl;
    }
  }
  if(nbtotal % nproc != 0 && adaptive == false && maxcost == mincost && myrank==0) {
    std::cout << "### Warning in Mesh constructor" << std::endl
              << "The number of MeshBlocks cannot be divided evenly. "
              << "This will cause a poor load balance." << std::endl;
  }
#endif

  // Mesh test only; do not create meshes
  if(test_flag>0) {
    if(myrank==0)
      MeshTest(dim);
    return;
  }

// create MeshBlock list for this process
  for(i=nbstart;i<=nbend;i++) {
    buid[i].GetLocation(lx1,lx2,lx3,ll);

    // calculate physical block size, x1
    if(lx1==0) {
      block_size.x1min=mesh_size.x1min;
      block_bcs[inner_x1]=mesh_bcs[inner_x1];
    }
    else {
      Real rx=(Real)lx1/(Real)(nrbx1<<(ll-root_level));
      block_size.x1min=MeshGeneratorX1(rx,mesh_size);
      block_bcs[inner_x1]=-1;
    }
    if(lx1==(nrbx1<<(ll-root_level))-1) {
      block_size.x1max=mesh_size.x1max;
      block_bcs[outer_x1]=mesh_bcs[outer_x1];
    }
    else {
      Real rx=(Real)(lx1+1)/(Real)(nrbx1<<(ll-root_level));
      block_size.x1max=MeshGeneratorX1(rx,mesh_size);
      block_bcs[outer_x1]=-1;
    }

    // calculate physical block size, x2
    if(dim==1) {
      block_size.x2min=mesh_size.x2min;
      block_size.x2max=mesh_size.x2max;
      block_bcs[inner_x2]=mesh_bcs[inner_x2];
      block_bcs[outer_x2]=mesh_bcs[outer_x2];
    }
    else {
      if(lx2==0) {
        block_size.x2min=mesh_size.x2min;
        block_bcs[inner_x2]=mesh_bcs[inner_x2];
      }
      else {
        Real rx=(Real)lx2/(Real)(nrbx2<<(ll-root_level));
        block_size.x2min=MeshGeneratorX2(rx,mesh_size);
        block_bcs[inner_x2]=-1;
      }
      if(lx2==(nrbx2<<(ll-root_level))-1) {
        block_size.x2max=mesh_size.x2max;
        block_bcs[outer_x2]=mesh_bcs[outer_x2];
      }
      else {
        Real rx=(Real)(lx2+1)/(Real)(nrbx2<<(ll-root_level));
        block_size.x2max=MeshGeneratorX2(rx,mesh_size);
        block_bcs[outer_x2]=-1;
      }
    }

    // calculate physical block size, x3
    if(dim<=2) {
      block_size.x3min=mesh_size.x3min;
      block_size.x3max=mesh_size.x3max;
      block_bcs[inner_x3]=mesh_bcs[inner_x3];
      block_bcs[outer_x3]=mesh_bcs[outer_x3];
    }
    else {
      if(lx3==0) {
        block_size.x3min=mesh_size.x3min;
        block_bcs[inner_x3]=mesh_bcs[inner_x3];
      }
      else {
        Real rx=(Real)lx3/(Real)(nrbx3<<(ll-root_level));
        block_size.x3min=MeshGeneratorX3(rx,mesh_size);
        block_bcs[inner_x3]=-1;
      }
      if(lx3==(nrbx3<<(ll-root_level))-1) {
        block_size.x3max=mesh_size.x3max;
        block_bcs[outer_x3]=mesh_bcs[outer_x3];
      }
      else {
        Real rx=(Real)(lx3+1)/(Real)(nrbx3<<(ll-root_level));
        block_size.x3max=MeshGeneratorX3(rx,mesh_size);
        block_bcs[outer_x3]=-1;
      }
    }

    // create a block and add into the link list
    if(i==nbstart) {
      pblock = new MeshBlock(i, i-nbstart, buid[i], block_size, block_bcs, this, pin);
      pfirst = pblock;
    }
    else {
      pblock->next = new MeshBlock(i, i-nbstart, buid[i], block_size, block_bcs, this, pin);
      pblock->next->prev = pblock;
      pblock = pblock->next;
    }

    pblock->SearchAndSetNeighbors(tree, ranklist, nslist);
  }
  pblock=pfirst;
}


//--------------------------------------------------------------------------------------
// Mesh constructor for restarting. Load the restarting file

Mesh::Mesh(ParameterInput *pin, WrapIO& resfile, int test_flag)
{
  std::stringstream msg;
  RegionSize block_size;
  MeshBlock *pfirst;
  int idl, i, j, ll, nerr, dim;
  long int lx1, lx2, lx3;
  WrapIOSize_t *offset;
  Real totalcost, targetcost, maxcost, mincost, mycost;
  ID_t *rawid;

// mesh test
  if(test_flag>0) nproc=test_flag;

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
  if(resfile.Read(&idl, sizeof(int), 1)!=1) nerr++;
  if(resfile.Read(&root_level, sizeof(int), 1)!=1) nerr++;
  if(resfile.Read(&mesh_size, sizeof(RegionSize), 1)!=1) nerr++;
  if(resfile.Read(mesh_bcs, sizeof(int), 6)!=6) nerr++;
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

  if(idl>IDLENGTH) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "IDLENGTH in the restarting files is larger than the current configuration"
        << std::endl << "Please reconfigure the code accordingly." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

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
  buid=new BlockUID[nbtotal];
  offset=new WrapIOSize_t[nbtotal];
  costlist=new Real[nbtotal];
  ranklist=new int[nbtotal];
  nslist=new int[nproc];
  nblist=new int[nproc];
  rawid=new ID_t[IDLENGTH];
  for(int i=0;i<IDLENGTH;i++) rawid[i]=0;

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
  totalcost=0.0;
  nerr=0;
  maxcost=0.0;
  mincost=(FLT_MAX);
  for(int i=0;i<nbtotal;i++) {
    int bgid,level;
    if(resfile.Read(&bgid,sizeof(int),1)!=1) nerr++;
    if(resfile.Read(&level,sizeof(int),1)!=1) nerr++;
    if(level!=root_level) multilevel=true;
    if(resfile.Read(rawid,sizeof(ID_t),idl)!=idl) nerr++;
    if(resfile.Read(&(costlist[i]),sizeof(Real),1)!=1) nerr++;
    if(resfile.Read(&(offset[i]),sizeof(WrapIOSize_t),1)!=1) nerr++;
    buid[i].SetUID(rawid,level);
    totalcost+=costlist[i];
    mincost=std::min(mincost,costlist[i]);
    maxcost=std::max(maxcost,costlist[i]);
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

  delete [] rawid; // no longer needed

  face_only=true;
  if (MAGNETIC_FIELDS_ENABLED || multilevel==true)
    face_only=false;

  maxneighbor_=BufferID(dim, multilevel, face_only);

  // rebuild the Block Tree
  for(int i=0;i<nbtotal;i++)
    tree.AddMeshBlockWithoutRefine(tree,buid[i],dim,mesh_bcs,nrbx1,nrbx2,nrbx3,root_level);
  int nnb;
  tree.AssignGID(nnb); // count blocks at the same time
  if(nnb!=nbtotal) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Tree reconstruction failed. The total numbers of the blocks do not match. ("
        << nbtotal << " != " << nnb << ")" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

#ifdef MPI_PARALLEL
  if(nbtotal < nproc) {
    if(test_flag==0) {
      msg << "### FATAL ERROR in Mesh constructor" << std::endl
          << "Too few blocks: nbtotal (" << nbtotal << ") < nproc ("<< nproc
          << ")" << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    else { // test
      std::cout << "### Warning in Mesh constructor" << std::endl
          << "Too few blocks: nbtotal (" << nbtotal << ") < nproc ("<< nproc
          << ")" << std::endl;
      return;
    }
  }
  if(nbtotal % nproc != 0 && adaptive == false && maxcost == mincost && myrank==0) {
    std::cout << "### Warning in Mesh constructor" << std::endl
              << "The number of MeshBlocks cannot be divided evenly. "
              << "This will cause a poor load balance." << std::endl;
  }
#endif

  // divide the list evenly and distribute among the processes
  // note: ordering should be maintained, although it might not be optimal.
  j=nproc-1;
  targetcost=totalcost/nproc;
  mycost=0.0;
  // create rank list from the end: the master node should have less load
  for(i=nbtotal-1;i>=0;i--) {
    mycost+=costlist[i];
    ranklist[i]=j;
    if(mycost >= targetcost && j>0) {
      j--;
      totalcost-=mycost;
      mycost=0.0;
      targetcost=totalcost/(j+1);
    }
  }

  if(nbtotal < nproc && test_flag==0) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Too few blocks: nbtotal (" << nbtotal << ") < nproc ("<< nproc
        << ")" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if(nbtotal % nproc != 0 && adaptive == false && maxcost == mincost && myrank==0) {
    std::cout << "### Warning in Mesh constructor" << std::endl
              << "The number of MeshBlocks cannot be divided evenly. "
              << "This will cause a poor load balance." << std::endl;
  }

  nslist[0]=0;
  j=0;
  for(i=1;i<nbtotal;i++) { // make the list of nbstart
    if(ranklist[i]!=ranklist[i-1]) {
      nblist[j]=i-nslist[j];
      nslist[++j]=i;
    }
  }
  nblist[j]=nbtotal-nslist[j];
  // store my nbstart and nbend
  nbstart=nslist[myrank];
  if(myrank+1==nproc)
    nbend=nbtotal-1;
  else 
    nbend=nslist[myrank+1]-1;

  // Mesh test only; do not create meshes
  if(test_flag>0) {
    if(myrank==0)
      MeshTest(dim);
    delete [] offset;
    return;
  }

  // load MeshBlocks (parallel)
  for(i=nbstart;i<=nbend;i++) {
    // create a block and add into the link list
    if(i==nbstart) {
      pblock = new MeshBlock(i, i-nbstart, this, pin, buid, resfile, offset[i],
                             costlist[i], ranklist, nslist);
      pfirst = pblock;
    }
    else {
      pblock->next = new MeshBlock(i, i-nbstart, this, pin, buid, resfile,
                                   offset[i], costlist[i], ranklist, nslist);
      pblock->next->prev = pblock;
      pblock = pblock->next;
    }
    pblock->SearchAndSetNeighbors(tree, ranklist, nslist);
  }
  pblock=pfirst;

// clean up
  delete [] offset;
}


// destructor

Mesh::~Mesh()
{
  while(pblock->prev != NULL) // should not be true
    delete pblock->prev;
  while(pblock->next != NULL)
    delete pblock->next;
  delete pblock;
  delete [] nslist;
  delete [] nblist;
  delete [] ranklist;
  delete [] costlist;
  delete [] buid;
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
  std::cout << "Logical level of maximum refinement = "<< max_level << std::endl;
  std::cout << "List of MeshBlocks" << std::endl;
  for(i=root_level;i<=max_level;i++) {
    Real dx=1.0/(Real)(1L<<i);
    nb[i-root_level]=0;
    for(j=0;j<nbtotal;j++) {
      if(buid[j].GetLevel()==i) {
        buid[j].GetLocation(lx1,lx2,lx3,ll);
        std::cout << "MeshBlock " << j << ", lx1 = "
                  << lx1 << ", lx2 = " << lx2 <<", lx3 = " << lx3
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
// constructs coordinate, boundary condition, fluid and field objects.

MeshBlock::MeshBlock(int igid, int ilid, BlockUID iuid, RegionSize input_block,
                     int *input_bcs, Mesh *pm, ParameterInput *pin)
{
  std::stringstream msg;
  long int lx1, lx2, lx3;
  int ll, root_level;
  RegionSize& mesh_size  = pm->mesh_size;
  long long nrootmesh, noffset;
  pmy_mesh = pm;
  root_level = pm->root_level;
  block_size = input_block;
  for(int i=0; i<6; i++) block_bcs[i] = input_bcs[i];
  prev=NULL;
  next=NULL;
  gid=igid;
  lid=ilid;
  uid=iuid;
  cost=1.0;
  task=NULL;

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

  uid.GetLocation(lx1,lx2,lx3,ll);
  std::cout << "MeshBlock " << gid << ", rank = " << myrank << ", lx1 = "
            << lx1 << ", lx2 = " << lx2 <<", lx3 = " << lx3 << ", level = "
            << ll << std::endl;
  std::cout << "is=" << is << " ie=" << ie << " x1min=" << block_size.x1min
            << " x1max=" << block_size.x1max << std::endl;
  std::cout << "js=" << js << " je=" << je << " x2min=" << block_size.x2min
            << " x2max=" << block_size.x2max << std::endl;
  std::cout << "ks=" << ks << " ke=" << ke << " x3min=" << block_size.x3min
            << " x3max=" << block_size.x3max << std::endl;

// allocate arrays for sizes and positions of cells

  int ncells1 = block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (block_size.nx2 > 1) ncells2 = block_size.nx2 + 2*(NGHOST);
  if (block_size.nx3 > 1) ncells3 = block_size.nx3 + 2*(NGHOST);

// cell sizes
  dx1f.NewAthenaArray(ncells1);
  dx2f.NewAthenaArray(ncells2);
  dx3f.NewAthenaArray(ncells3);
  dx1v.NewAthenaArray(ncells1);
  dx2v.NewAthenaArray(ncells2);
  dx3v.NewAthenaArray(ncells3);

// cell positions. Note the extra element for cell face positions
  x1f.NewAthenaArray((ncells1+1));
  x2f.NewAthenaArray((ncells2+1));
  x3f.NewAthenaArray((ncells3+1));
  x1v.NewAthenaArray(ncells1);
  x2v.NewAthenaArray(ncells2);
  x3v.NewAthenaArray(ncells3);

// X1-DIRECTION: initialize sizes and positions of cell FACES (dx1f,x1f)

  nrootmesh=mesh_size.nx1*(1L<<(ll-root_level));
  if(block_size.x1rat == 1.0) { // uniform
    Real dx=(block_size.x1max-block_size.x1min)/block_size.nx1;
    for(int i=is-NGHOST; i<=ie+NGHOST; ++i)
      dx1f(i)=dx;
    x1f(is-NGHOST)=block_size.x1min-NGHOST*dx;
    for(int i=is-NGHOST+1;i<=ie+NGHOST+1;i++)
      x1f(i)=x1f(i-1)+dx;
    x1f(is) = block_size.x1min;
    x1f(ie+1) = block_size.x1max;
  }
  else {
    for (int i=is-NGHOST; i<=ie+NGHOST+1; ++i) {
       // if there are too many levels, this won't work or be precise enough
      noffset=(i-is)+(long long)lx1*block_size.nx1;
      Real rx=(Real)noffset/(Real)nrootmesh;
      x1f(i)=pm->MeshGeneratorX1(rx,mesh_size);
    }
    x1f(is) = block_size.x1min;
    x1f(ie+1) = block_size.x1max;
    for(int i=is-NGHOST; i<=ie+NGHOST; ++i)
      dx1f(i)=x1f(i+1)-x1f(i);
  }

// correct cell face positions in ghost zones for reflecting boundary condition
  if (block_bcs[inner_x1] == 1) {
    for (int i=1; i<=(NGHOST); ++i) {
      dx1f(is-i) = dx1f(is+i-1);
       x1f(is-i) =  x1f(is-i+1) - dx1f(is-i);
    }
  }
  if (block_bcs[outer_x1] == 1) {
    for (int i=1; i<=(NGHOST); ++i) {
      dx1f(ie+i  ) = dx1f(ie-i+1);
       x1f(ie+i+1) =  x1f(ie+i) + dx1f(ie+i);
    }
  }

// X2-DIRECTION: initialize spacing and positions of cell FACES (dx2f,x2f)

  if(block_size.nx2 > 1) {
    nrootmesh=mesh_size.nx2*(1L<<(ll-root_level));
    if(block_size.x2rat == 1.0) { // uniform
      Real dx=(block_size.x2max-block_size.x2min)/block_size.nx2;
      for(int j=js-NGHOST; j<=je+NGHOST; ++j)
        dx2f(j)=dx;
      x2f(js-NGHOST)=block_size.x2min-NGHOST*dx;
      for(int j=js-NGHOST+1;j<=je+NGHOST+1;j++)
        x2f(j)=x2f(j-1)+dx;
      x2f(js) = block_size.x2min;
      x2f(je+1) = block_size.x2max;
    }
    else {
      for (int j=js-NGHOST; j<=je+NGHOST+1; ++j) {
         // if there are too many levels, this won't work or be precise enough
        noffset=(j-js)+(long long)lx2*block_size.nx2;
        Real rx=(Real)noffset/(Real)nrootmesh;
        x2f(j)=pm->MeshGeneratorX2(rx,mesh_size);
      }
      x2f(js) = block_size.x2min;
      x2f(je+1) = block_size.x2max;
      for(int j=js-NGHOST; j<=je+NGHOST; ++j)
        dx2f(j)=x2f(j+1)-x2f(j);
    }

  // correct cell face positions in ghost zones for reflecting boundary condition
    if (block_bcs[inner_x2] == 1) {
      for (int j=1; j<=(NGHOST); ++j) {
        dx2f(js-j) = dx2f(js+j-1);
         x2f(js-j) =  x2f(js-j+1) - dx2f(js-j);
      }
    }
    if (block_bcs[outer_x2] == 1) {
      for (int j=1; j<=(NGHOST); ++j) {
        dx2f(je+j  ) = dx2f(je-j+1);
         x2f(je+j+1) =  x2f(je+j) + dx2f(je+j);
      }
    }
  }
  else {
    dx2f(js) = block_size.x2max-block_size.x2min;
    x2f(js  ) = block_size.x2min;
    x2f(je+1) = block_size.x2max;
  }


// X3-DIRECTION: initialize spacing and positions of cell FACES (dx3f,x3f)


  if(block_size.nx3 > 1) {
    nrootmesh=mesh_size.nx3*(1L<<(ll-root_level));
    if(block_size.x3rat == 1.0) { // uniform
      Real dx=(block_size.x3max-block_size.x3min)/block_size.nx3;
      for(int k=ks-NGHOST; k<=ke+NGHOST; ++k)
        dx3f(k)=dx;
      x3f(ks-NGHOST)=block_size.x3min-NGHOST*dx;
      for(int k=ks-NGHOST+1;k<=ke+NGHOST+1;k++)
        x3f(k)=x3f(k-1)+dx;
      x3f(ks) = block_size.x3min;
      x3f(ke+1) = block_size.x3max;
    }
    else {
      for (int k=ks-NGHOST; k<=ke+NGHOST+1; ++k) {
         // if there are too many levels, this won't work or be precise enough
        noffset=(k-ks)+(long long)lx3*block_size.nx3;
        Real rx=(Real)noffset/(Real)nrootmesh;
        x3f(k)=pm->MeshGeneratorX3(rx,mesh_size);
      }
      x3f(ks) = block_size.x3min;
      x3f(ke+1) = block_size.x3max;
      for(int k=ks-NGHOST; k<=ke+NGHOST; ++k)
        dx3f(k)=x3f(k+1)-x3f(k);
    }

  // correct cell face positions in ghost zones for reflecting boundary condition
    if (block_bcs[inner_x3] == 1) {
      for (int k=1; k<=(NGHOST); ++k) {
        dx3f(ks-k) = dx3f(ks+k-1);
         x3f(ks-k) =  x3f(ks-k+1) - dx3f(ks-k);
      }
    }
    if (block_bcs[outer_x3] == 1) {
      for (int k=1; k<=(NGHOST); ++k) {
        dx3f(ke+k  ) = dx3f(ke-k+1);
         x3f(ke+k+1) =  x3f(ke+k) + dx3f(ke+k);
      }
    }
  }
  else {
    dx3f(ks) = block_size.x3max-block_size.x3min;
    x3f(ks  ) = block_size.x3min;
    x3f(ke+1) = block_size.x3max;
  }

  if(pm->multilevel==true) // set coarse coordinates for SMR/AMR
    SetCoarserCoordinates();

// construct Coordinates and Fluid objects stored in MeshBlock class.  Note that the
// initial conditions for the fluid are set in problem generator called from main, not
// in the Fluid constructor
 
  pcoord = new Coordinates(this, pin);
  pfluid = new Fluid(this, pin);
  pfield = new Field(this, pin);
  pbval  = new BoundaryValues(this, pin);

  return;
}

//--------------------------------------------------------------------------------------
// MeshBlock constructor for restarting

MeshBlock::MeshBlock(int igid, int ilid, Mesh *pm, ParameterInput *pin, BlockUID *list,
           WrapIO& resfile, WrapIOSize_t offset, Real icost, int *ranklist, int *nslist)
{
  std::stringstream msg;
  pmy_mesh = pm;
  prev=NULL;
  next=NULL;
  gid=igid;
  lid=ilid;
  uid=list[gid];
  cost=icost;
  int nerr=0;
  task=NULL;

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

  long int lx1, lx2, lx3;
  int ll;
  uid.GetLocation(lx1,lx2,lx3,ll);
  std::cout << "MeshBlock " << gid << ", rank = " << myrank << ", lx1 = "
            << lx1 << ", lx2 = " << lx2 <<", lx3 = " << lx3 << ", level = "
            << ll << std::endl;
  std::cout << "is=" << is << " ie=" << ie << " x1min=" << block_size.x1min
            << " x1max=" << block_size.x1max << std::endl;
  std::cout << "js=" << js << " je=" << je << " x2min=" << block_size.x2min
            << " x2max=" << block_size.x2max << std::endl;
  std::cout << "ks=" << ks << " ke=" << ke << " x3min=" << block_size.x3min
            << " x3max=" << block_size.x3max << std::endl;

// allocate arrays for sizes and positions of cells

  int ncells1 = block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (block_size.nx2 > 1) ncells2 = block_size.nx2 + 2*(NGHOST);
  if (block_size.nx3 > 1) ncells3 = block_size.nx3 + 2*(NGHOST);

// cell sizes
  dx1f.NewAthenaArray(ncells1);
  dx2f.NewAthenaArray(ncells2);
  dx3f.NewAthenaArray(ncells3);
  dx1v.NewAthenaArray(ncells1);
  dx2v.NewAthenaArray(ncells2);
  dx3v.NewAthenaArray(ncells3);

// cell positions. Note the extra element for cell face positions
  x1f.NewAthenaArray((ncells1+1));
  x2f.NewAthenaArray((ncells2+1));
  x3f.NewAthenaArray((ncells3+1));
  x1v.NewAthenaArray(ncells1);
  x2v.NewAthenaArray(ncells2);
  x3v.NewAthenaArray(ncells3);

  //load x1f, x2f, x3f
  nerr=0;
  if(resfile.Read(x1f.GetArrayPointer(),sizeof(Real),x1f.GetDim1())
     !=x1f.GetDim1()) nerr++;
  if(resfile.Read(x2f.GetArrayPointer(),sizeof(Real),x2f.GetDim1())
     !=x2f.GetDim1()) nerr++;
  if(resfile.Read(x3f.GetArrayPointer(),sizeof(Real),x3f.GetDim1())
     !=x3f.GetDim1()) nerr++;
  if(resfile.Read(dx1f.GetArrayPointer(),sizeof(Real),dx1f.GetDim1())
     !=dx1f.GetDim1()) nerr++;
  if(resfile.Read(dx2f.GetArrayPointer(),sizeof(Real),dx2f.GetDim1())
     !=dx2f.GetDim1()) nerr++;
  if(resfile.Read(dx3f.GetArrayPointer(),sizeof(Real),dx3f.GetDim1())
     !=dx3f.GetDim1()) nerr++;
  if(nerr>0) {
    msg << "### FATAL ERROR in MeshBlock constructor" << std::endl
        << "The restarting file is broken." << std::endl;
    resfile.Close();
    throw std::runtime_error(msg.str().c_str());
  }

  if(pmy_mesh->multilevel==true) // set coarse coordinates for SMR/AMR
    SetCoarserCoordinates();

  // create coordinates, fluid, field, and boundary conditions
  pcoord = new Coordinates(this, pin);
  pfluid = new Fluid(this, pin);
  pfield = new Field(this, pin);
  pbval  = new BoundaryValues(this, pin);

  // load fluid and field data
  nerr=0;
  if(resfile.Read(pfluid->u.GetArrayPointer(),sizeof(Real),
                         pfluid->u.GetSize())!=pfluid->u.GetSize()) nerr++;
  if (GENERAL_RELATIVITY) {
    if(resfile.Read(pfluid->w.GetArrayPointer(),sizeof(Real),
                           pfluid->w.GetSize())!=pfluid->w.GetSize()) nerr++;
    if(resfile.Read(pfluid->w1.GetArrayPointer(),sizeof(Real),
                           pfluid->w1.GetSize())!=pfluid->w1.GetSize()) nerr++;
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
  dx1f.DeleteAthenaArray();  
  dx2f.DeleteAthenaArray();  
  dx3f.DeleteAthenaArray();  
  dx1v.DeleteAthenaArray();  
  dx2v.DeleteAthenaArray();  
  dx3v.DeleteAthenaArray();  
  x1f.DeleteAthenaArray();
  x2f.DeleteAthenaArray();
  x3f.DeleteAthenaArray();
  x1v.DeleteAthenaArray();
  x2v.DeleteAthenaArray();
  x3v.DeleteAthenaArray();
  if(pmy_mesh->multilevel==true) {
    coarse_dx1f.DeleteAthenaArray();  
    coarse_dx2f.DeleteAthenaArray();  
    coarse_dx3f.DeleteAthenaArray();  
    coarse_dx1v.DeleteAthenaArray();  
    coarse_dx2v.DeleteAthenaArray();  
    coarse_dx3v.DeleteAthenaArray();  
    coarse_x1f.DeleteAthenaArray();
    coarse_x2f.DeleteAthenaArray();
    coarse_x3f.DeleteAthenaArray();
    coarse_x1v.DeleteAthenaArray();
    coarse_x2v.DeleteAthenaArray();
    coarse_x3v.DeleteAthenaArray();
    coarse_data.DeleteAthenaArray();
  }

  if(prev!=NULL) prev->next=next;
  if(next!=NULL) next->prev=prev;

  delete pcoord;
  delete pfluid;
  delete pfield;
  delete pbval;
  delete [] task;
}



//--------------------------------------------------------------------------------------
// \!fn void MeshBlock::SetCoarserCoordinates(void)
// \brief set 1-level coarser coordinates for SMR/AMR
void MeshBlock::SetCoarserCoordinates(void) {
  Mesh *pm=pmy_mesh;
  long long nrootmesh, noffset;
  long int lx1, lx2, lx3;
  int ll;
  uid.GetLocation(lx1,lx2,lx3,ll);
  nrootmesh=pm->mesh_size.nx1*(1L<<(uid.GetLevel()-pm->root_level));
  cnghost=(NGHOST+1)/2+1;
  
  int ncc1=block_size.nx1/2+2*cnghost;
  cis=cnghost; cie=cis+block_size.nx1/2-1;
  int ncc2=1, ncc3=1;
  cjs=cje=cks=cke=0;
  if(block_size.nx2>1) { // 2D or 3D
    ncc2=block_size.nx2/2+2*cnghost;
    cjs=cnghost; cje=cjs+block_size.nx2/2-1;
  }
  if(block_size.nx3>1) { // 3D
    ncc3=block_size.nx3/2+2*cnghost;
    cks=cnghost; cke=cks+block_size.nx3/2-1;
  }
  // cell sizes
  coarse_dx1f.NewAthenaArray(ncc1);
  coarse_dx2f.NewAthenaArray(ncc2);
  coarse_dx3f.NewAthenaArray(ncc3);
  coarse_dx1v.NewAthenaArray(ncc1);
  coarse_dx2v.NewAthenaArray(ncc2);
  coarse_dx3v.NewAthenaArray(ncc3);
  // cell positions. Note the extra element for cell face positions
  coarse_x1f.NewAthenaArray((ncc1+1));
  coarse_x2f.NewAthenaArray((ncc2+1));
  coarse_x3f.NewAthenaArray((ncc3+1));
  coarse_x1v.NewAthenaArray(ncc1);
  coarse_x2v.NewAthenaArray(ncc2);
  coarse_x3v.NewAthenaArray(ncc3);
  coarse_data.NewAthenaArray(NFLUID, ncc3, ncc2, ncc1);

  // x1
  for(int i=cis; i<=cie; i++) { // active region
    int ifl=(i-cis)*2+is;
    coarse_dx1f(i)=dx1f(ifl)+dx1f(ifl+1);
    coarse_x1f(i)=x1f(ifl);
  }
  coarse_x1f(cie+1)=x1f(ie+1);
  // left ghost zone
  if(block_bcs[inner_x1]==1) { // reflecting
    for(int i=1; i<=cnghost; i++) {
      coarse_dx1f(cis-i) = coarse_dx1f(cis+i-1);
      coarse_x1f(cis-i) =  coarse_x1f(cis-i+1) - coarse_dx1f(cis-i);
    }
  }
  else {
    if(block_size.x1rat==1.0) { //uniform
      for(int i=1; i<=cnghost; i++) {
        coarse_dx1f(cis-i) = coarse_dx1f(cis-i+1);
        coarse_x1f(cis-i) =  coarse_x1f(cis-i+1) - coarse_dx1f(cis-i);
      }
    }
    else {
      for(int i=1; i<=cnghost; i++) { 
        noffset=(long long)lx1*block_size.nx1-i*2;
        Real rx=(Real)noffset/(Real)nrootmesh;
        coarse_x1f(cis-i)=pm->MeshGeneratorX1(rx,pm->mesh_size);
        coarse_dx1f(cis-i) = coarse_x1f(cis-i+1)-coarse_x1f(cis-i);
      }
    }
  }
  // right ghost zone
  if(block_bcs[outer_x1]==1) { // reflecting
    for(int i=1; i<=cnghost; i++) {
      coarse_dx1f(cie+i) = coarse_dx1f(cie-i+1);
      coarse_x1f(cie+i+1) =  coarse_x1f(cie+i) + coarse_dx1f(cie+i);
    }
  }
  else {
    if(block_size.x1rat==1.0) { //uniform
      for(int i=1; i<=cnghost; i++) {
        coarse_dx1f(cie+i) = coarse_dx1f(cie+i-1);
        coarse_x1f(cie+i+1) =  coarse_x1f(cie+i) + coarse_dx1f(cie+i);
      }
    }
    else {
      for(int i=1; i<=cnghost; i++) { 
        noffset=(long long)(lx1+1L)*block_size.nx1+i*2;
        Real rx=(Real)noffset/(Real)nrootmesh;
        coarse_x1f(cie+i+1)=pm->MeshGeneratorX1(rx,pm->mesh_size);
        coarse_dx1f(cie+i) = coarse_x1f(cie+i+1)-coarse_x1f(cie+i);
      }
    }
  }

  // x2
  for(int j=cjs; j<=cje; j++) { // active region
    int jfl=(j-cjs)*2+js;
    coarse_dx2f(j)=dx2f(jfl)+dx2f(jfl+1);
    coarse_x2f(j)=x2f(jfl);
  }
  coarse_x2f(cje+1)=x2f(je+1);
  // left ghost zone
  if(block_bcs[inner_x2]==1) { // reflecting
    for(int j=1; j<=cnghost; j++) {
      coarse_dx2f(cjs-j) = coarse_dx2f(cjs+j-1);
      coarse_x2f(cjs-j) =  coarse_x2f(cjs-j+1) - coarse_dx2f(cjs-j);
    }
  }
  else {
    if(block_size.x2rat==1.0) { //uniform
      for(int j=1; j<=cnghost; j++) {
        coarse_dx2f(cjs-j) = coarse_dx2f(cjs-j+1);
        coarse_x2f(cjs-j) =  coarse_x2f(cjs-j+1) - coarse_dx2f(cjs-j);
      }
    }
    else {
      for(int j=1; j<=cnghost; j++) { 
        noffset=(long long)lx2*block_size.nx2-j*2;
        Real rx=(Real)noffset/(Real)nrootmesh;
        coarse_x2f(cjs-j)=pm->MeshGeneratorX2(rx,pm->mesh_size);
        coarse_dx2f(cjs-j) = coarse_x2f(cjs-j+1)-coarse_x2f(cjs-j);
      }
    }
  }
  // right ghost zone
  if(block_bcs[outer_x2]==1) { // reflecting
    for(int j=1; j<=cnghost; j++) {
      coarse_dx2f(cje+j) = coarse_dx2f(cje-j+1);
      coarse_x2f(cje+j+1) =  coarse_x2f(cje+j) + coarse_dx2f(cje+j);
    }
  }
  else {
    if(block_size.x2rat==1.0) { //uniform
      for(int j=1; j<=cnghost; j++) {
        coarse_dx2f(cje+j) = coarse_dx2f(cje+j-1);
        coarse_x2f(cje+j+1) =  coarse_x2f(cje+j) + coarse_dx2f(cje+j);
      }
    }
    else {
      for(int j=1; j<=cnghost; j++) { 
        noffset=(long long)(lx2+1L)*block_size.nx2+j*2;
        Real rx=(Real)noffset/(Real)nrootmesh;
        coarse_x2f(cje+j+1)=pm->MeshGeneratorX2(rx,pm->mesh_size);
        coarse_dx2f(cje+j) = coarse_x2f(cje+j+1)-coarse_x2f(cje+j);
      }
    }
  }

  // x3
  for(int k=cks; k<=cke; k++) { // active region
    int kfl=(k-cks)*2+ks;
    coarse_dx3f(k)=dx3f(kfl)+dx3f(kfl+1);
    coarse_x3f(k)=x3f(kfl);
  }
  coarse_x3f(cke+1)=x3f(ke+1);
  // left ghost zone
  if(block_bcs[inner_x3]==1) { // reflecting
    for(int k=1; k<=cnghost; k++) {
      coarse_dx3f(cks-k) = coarse_dx3f(cks+k-1);
      coarse_x3f(cks-k) =  coarse_x3f(cks-k+1) - coarse_dx3f(cks-k);
    }
  }
  else {
    if(block_size.x3rat==1.0) { //uniform
      for(int k=1; k<=cnghost; k++) {
        coarse_dx3f(cks-k) = coarse_dx3f(cks-k+1);
        coarse_x3f(cks-k) =  coarse_x3f(cks-k+1) - coarse_dx3f(cks-k);
      }
    }
    else {
      for(int k=1; k<=cnghost; k++) { 
        noffset=(long long)lx3*block_size.nx3-k*2;
        Real rx=(Real)noffset/(Real)nrootmesh;
        coarse_x3f(cks-k)=pm->MeshGeneratorX3(rx,pm->mesh_size);
        coarse_dx3f(cks-k) = coarse_x3f(cks-k+1)-coarse_x3f(cks-k);
      }
    }
  }
  // right ghost zone
  if(block_bcs[outer_x3]==1) { // reflecting
    for(int k=1; k<=cnghost; k++) {
      coarse_dx3f(cke+k) = coarse_dx3f(cke-k+1);
      coarse_x3f(cke+k+1) =  coarse_x3f(cke+k) + coarse_dx3f(cke+k);
    }
  }
  else {
    if(block_size.x3rat==1.0) { //uniform
      for(int k=1; k<=cnghost; k++) {
        coarse_dx3f(cke+k) = coarse_dx3f(cke+k-1);
        coarse_x3f(cke+k+1) =  coarse_x3f(cke+k) + coarse_dx3f(cke+k);
      }
    }
    else {
      for(int k=1; k<=cnghost; k++) { 
        noffset=(long long)(lx3+1L)*block_size.nx3+k*2;
        Real rx=(Real)noffset/(Real)nrootmesh;
        coarse_x3f(cke+k+1)=pm->MeshGeneratorX3(rx,pm->mesh_size);
        coarse_dx3f(cke+k) = coarse_x3f(cke+k+1)-coarse_x3f(cke+k);
      }
    }
  }
}

//--------------------------------------------------------------------------------------
// \!fn void Mesh::NewTimeStep(void)
// \brief function that loops over all MeshBlocks and find new timestep
//        this assumes that pfluid->NewBlockTimeStep is already called
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
// \!fn void Mesh::Initialize(int res_flag, ParameterInput *pin)
// \brief  initialization before the main loop

void Mesh::Initialize(int res_flag, ParameterInput *pin)
{
  MeshBlock *pmb;
  Fluid *pfluid;
  Field *pfield;
  BoundaryValues *pbval;

  if(res_flag==0) {
    pmb = pblock;
    while (pmb != NULL)  {
      pfluid=pmb->pfluid;
      pfield=pmb->pfield;
      pbval=pmb->pbval;
      ProblemGenerator(pfluid,pfield,pin);
      pbval->CheckBoundary();
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
    pfluid=pmb->pfluid;
    pfield=pmb->pfield;
    pbval=pmb->pbval;
    pbval->SendFluidBoundaryBuffers(pfluid->u,0);
    if (MAGNETIC_FIELDS_ENABLED)
      pbval->SendFieldBoundaryBuffers(pfield->b,0);
    pmb=pmb->next;
  }

  pmb = pblock;
  while (pmb != NULL)  {
    pfluid=pmb->pfluid;
    pfield=pmb->pfield;
    pbval=pmb->pbval;
    pbval->ReceiveFluidBoundaryBuffersWithWait(pfluid->u ,0);
    if (MAGNETIC_FIELDS_ENABLED)
      pbval->ReceiveFieldBoundaryBuffersWithWait(pfield->b ,0);
    pmb->pbval->ClearBoundaryForInit();
    pbval->FluidPhysicalBoundaries(pfluid->u);
    if(multilevel==true)
      pbval->ProlongateFluidBoundaries(pfluid->u);
    if (MAGNETIC_FIELDS_ENABLED) {
      pbval->FieldPhysicalBoundaries(pfield->b);
//      if(multilevel==true)
//        pbval->ProlongateFieldBoundaries(pfield->b);
    }
    pfluid->pf_eos->ConservedToPrimitive(pfluid->u, pfluid->w1, pfield->b, 
                                         pfluid->w, pfield->bcc);
    pmb=pmb->next;
  }

  if(res_flag==0) {
    pmb = pblock;
    while (pmb != NULL)  {
      pmb->pfluid->NewBlockTimeStep(pmb);
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
  size+=sizeof(Real)*(x1f.GetSize()+x2f.GetSize()+x3f.GetSize());
  size+=sizeof(Real)*(dx1f.GetSize()+dx2f.GetSize()+dx3f.GetSize());
  size+=sizeof(Real)*pfluid->u.GetSize();
  if (GENERAL_RELATIVITY) {
    size+=sizeof(Real)*pfluid->w.GetSize();
    size+=sizeof(Real)*pfluid->w1.GetSize();
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
  int nb=nbend-nbstart+1;
  MeshBlock *pmb = pblock;
  // initialize
  while (pmb != NULL)  {
    pmb->firsttask=0;
    pmb->ntodo=pmb->ntask;
    pmb->task_flag=0L;
    pmb->pbval->StartReceivingAll();
    pmb=pmb->next;
  }

  // main loop
  while(nb>0) {
    pmb = pblock;
    while (pmb != NULL)  {
      if(pmb->DoOneTask()==tl_complete) // task list completed
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
//! \fn void Mesh::SetTaskList(TaskList& tl)
//  \brief set task list for all the meshblocks
void Mesh::SetTaskList(TaskList& tl)
{
  MeshBlock *pmb = pblock;
  while (pmb != NULL)  {
    pmb->SetTaskList(tl);
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
//! \fn void MeshBlock::SetTaskList(TaskList& tl)
//  \brief set task list for the meshblock
void MeshBlock::SetTaskList(TaskList& tl)
{
  if(task!=NULL)
    delete [] task;
  task=new Task[tl.ntask];
  ntask=tl.ntask;
  memcpy(task, tl.task, sizeof(Task)*ntask);
  return;
}


//--------------------------------------------------------------------------------------
//! \fn enum tlstatus MeshBlock::DoOneTask(void)
//  \brief process one task (if possible), return true if the list is completed
enum tasklist_status MeshBlock::DoOneTask(void) {
  int skip=0;
  enum task_status ret;
  std::stringstream msg;
  if(ntodo==0) return tl_nothing;
  for(int i=firsttask; i<ntask; i++) {
    Task &ti=task[i];
    if((ti.taskid & task_flag)==0L) { // this task is not done
      if (((ti.depend & task_flag) == ti.depend)) { // dependency clear
        ret=ti.TaskFunc(this,ti.task_arg);
//        std::cout << "Meshblock " << gid << " task "<< ti.taskid << " returns " << ret << std::endl;
        if(ret!=task_failure) { // success
          ntodo--;
          task_flag |= ti.taskid;
          if(skip==0)
            firsttask++;
          if(ntodo==0)
            return tl_complete;
          if(ret==task_donext) continue;
          return tl_running;
        }
      }
      skip++; // count the number of skipped tasks
    }
    else if(skip==0) // this task is done and at the top of the list
      firsttask++;
  }
  return tl_stuck; // there are still something to do but nothing can be done now
}

//--------------------------------------------------------------------------------------
// \!fn void NeighborBlock::SetNeighbor(int irank, int ilevel, int igid, int ilid,
//                          int iox1, int iox2, int iox3, enum neighbor_type itype,
//                          int ibid, int itargetid, int ifi1=0, int ifi2=0)
// \brief Set neighbor information
void NeighborBlock::SetNeighbor(int irank, int ilevel, int igid, int ilid,
  int iox1, int iox2, int iox3, enum neighbor_type itype, int ibid, int itargetid,
  int ifi1=0, int ifi2=0)
{
  rank=irank; level=ilevel; gid=igid; lid=ilid; ox1=iox1; ox2=iox2; ox3=iox3; type=itype;
  bufid=ibid; targetid=itargetid; fi1=ifi1; fi2=ifi2;
  if(type==neighbor_face) {
    if(ox1==-1)      fid=inner_x1;
    else if(ox1==1)  fid=outer_x1;
    else if(ox2==-1) fid=inner_x2;
    else if(ox2==1)  fid=outer_x2;
    else if(ox3==-1) fid=inner_x3;
    else if(ox3==1)  fid=outer_x3;
  }
  return;
}

//--------------------------------------------------------------------------------------
// \!fn void MeshBlock::SearchAndSetNeighbors(BlockTree &tree, int *ranklist, int *nslist)
// \brief Search and set all the neighbor blocks
void MeshBlock::SearchAndSetNeighbors(BlockTree &tree, int *ranklist, int *nslist)
{
  BlockTree* neibt;
  long int lx1, lx2, lx3;
  int ll, myox1, myox2, myox3, myfx1, myfx2, myfx3;
  uid.GetLocation(lx1,lx2,lx3,ll);
  myox1=((int)(lx1&1L))*2-1;
  myox2=((int)(lx2&1L))*2-1;
  myox3=((int)(lx3&1L))*2-1;
  myfx1=(int)(lx1&1L);
  myfx2=(int)(lx2&1L);
  myfx3=(int)(lx3&1L);
  long int nrbx1=pmy_mesh->nrbx1, nrbx2=pmy_mesh->nrbx2, nrbx3=pmy_mesh->nrbx3;

  int nf1=1, nf2=1;
  if(pmy_mesh->multilevel==true) {
    if(block_size.nx2>1) nf1=2;
    if(block_size.nx3>1) nf2=2;
  }
  int bufid=0;
  nneighbor=0;
  nblevel[1][1][1]=ll;
  // x1 face
  for(int n=-1; n<=1; n+=2) {
    neibt=tree.FindNeighbor(uid,n,0,0,block_bcs,nrbx1,nrbx2,nrbx3,pmy_mesh->root_level);
    if(neibt==NULL) { bufid+=nf1*nf2; nblevel[1][1][n+1]=-2; continue;}
    if(neibt->flag==false) { // finer
      int fface=1-(n+1)/2; // 0 for outer_x3, 1 for inner_x3
      nblevel[1][1][n+1]=neibt->uid.GetLevel()+1;
      for(int f2=0;f2<nf2;f2++) {
        for(int f1=0;f1<nf1;f1++) {
          BlockTree* nf=neibt->GetLeaf(fface,f1,f2);
          int fid = nf->gid;
          int nlevel=nf->uid.GetLevel();
          int tbid=FindBufferID(-n,0,0,0,0,pmy_mesh->maxneighbor_);
          neighbor[nneighbor].SetNeighbor(ranklist[fid], nlevel, fid,
            fid-nslist[ranklist[fid]], n, 0, 0, neighbor_face, bufid, tbid, f1, f2);
          bufid++; nneighbor++;
        }
      }
    }
    else {
      int nlevel=neibt->uid.GetLevel();
      int nid=neibt->gid;
      nblevel[1][1][n+1]=nlevel;
      int tbid;
      if(nlevel==ll) tbid=FindBufferID(-n,0,0,0,0,pmy_mesh->maxneighbor_);
      else tbid=FindBufferID(-n,0,0,myfx2,myfx3,pmy_mesh->maxneighbor_);
      neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
        nid-nslist[ranklist[nid]], n, 0, 0, neighbor_face, bufid, tbid);
      bufid+=nf1*nf2; nneighbor++;
    }
  }
  if(block_size.nx2==1) return;
  // x2 face
  for(int n=-1; n<=1; n+=2) {
    neibt=tree.FindNeighbor(uid,0,n,0,block_bcs,nrbx1,nrbx2,nrbx3,pmy_mesh->root_level);
    if(neibt==NULL) { bufid+=nf1*nf2; nblevel[1][n+1][1]=-2; continue;}
    if(neibt->flag==false) { // finer
      int fface=1-(n+1)/2; // 0 for outer_x3, 1 for inner_x3
      nblevel[1][n+1][1]=neibt->uid.GetLevel()+1;
      for(int f2=0;f2<nf2;f2++) {
        for(int f1=0;f1<nf1;f1++) {
          BlockTree* nf=neibt->GetLeaf(f1,fface,f2);
          int fid = nf->gid;
          int nlevel=nf->uid.GetLevel();
          int tbid=FindBufferID(0,-n,0,0,0,pmy_mesh->maxneighbor_);
          neighbor[nneighbor].SetNeighbor(ranklist[fid], nlevel, fid,
            fid-nslist[ranklist[fid]], 0, n, 0, neighbor_face, bufid, tbid, f1, f2);
          bufid++; nneighbor++;
        }
      }
    }
    else {
      int nlevel=neibt->uid.GetLevel();
      int nid=neibt->gid;
      nblevel[1][n+1][1]=nlevel;
      int tbid;
      if(nlevel==ll) tbid=FindBufferID(0,-n,0,0,0,pmy_mesh->maxneighbor_);
      else tbid=FindBufferID(0,-n,0,myfx1,myfx3,pmy_mesh->maxneighbor_);
      neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
        nid-nslist[ranklist[nid]], 0, n, 0, neighbor_face, bufid, tbid);
      bufid+=nf1*nf2; nneighbor++;
    }
  }
  if(block_size.nx3>1) {
    // x3 face
    for(int n=-1; n<=1; n+=2) {
      neibt=tree.FindNeighbor(uid,0,0,n,block_bcs,nrbx1,nrbx2,nrbx3,pmy_mesh->root_level);
      if(neibt==NULL) { bufid+=nf1*nf2; nblevel[n+1][1][1]=-2; continue;}
      if(neibt->flag==false) { // finer
        int fface=1-(n+1)/2; // 0 for outer_x3, 1 for inner_x3
        nblevel[n+1][1][1]=neibt->uid.GetLevel()+1;
        for(int f2=0;f2<nf2;f2++) {
          for(int f1=0;f1<nf1;f1++) {
            BlockTree* nf=neibt->GetLeaf(f1,f2,fface);
            int fid = nf->gid;
            int nlevel=nf->uid.GetLevel();
            int tbid=FindBufferID(0,0,-n,0,0,pmy_mesh->maxneighbor_);
            neighbor[nneighbor].SetNeighbor(ranklist[fid], nlevel, fid,
              fid-nslist[ranklist[fid]], 0, 0, n, neighbor_face, bufid, tbid, f1, f2);
            bufid++; nneighbor++;
          }
        }
      }
      else {
        int nlevel=neibt->uid.GetLevel();
        int nid=neibt->gid;
        nblevel[n+1][1][1]=nlevel;
        int tbid;
        if(nlevel==ll) tbid=FindBufferID(0,0,-n,0,0,pmy_mesh->maxneighbor_);
        else tbid=FindBufferID(0,0,-n,myfx1,myfx2,pmy_mesh->maxneighbor_);
        neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
          nid-nslist[ranklist[nid]], 0, 0, n, neighbor_face, bufid, tbid);
        bufid+=nf1*nf2; nneighbor++;
      }
    }
  }
  if(pmy_mesh->face_only==true) return;
  // edges
  // x1x2
  for(int m=-1; m<=1; m+=2) {
    for(int n=-1; n<=1; n+=2) {
      neibt=tree.FindNeighbor(uid,n,m,0,block_bcs,nrbx1,nrbx2,nrbx3,pmy_mesh->root_level);
      if(neibt==NULL) { bufid+=nf1; nblevel[1][m+1][n+1]=-2; continue;}
      if(neibt->flag==false) { // finer
        int ff1=1-(n+1)/2; // 0 for outer_x1, 1 for inner_x1
        int ff2=1-(m+1)/2; // 0 for outer_x2, 1 for inner_x2
        nblevel[1][m+1][n+1]=neibt->uid.GetLevel()+1;
        for(int f1=0;f1<nf1;f1++) {
          BlockTree* nf=neibt->GetLeaf(ff1,ff2,f1);
          int fid = nf->gid;
          int nlevel=nf->uid.GetLevel();
          int tbid=FindBufferID(-n,-m,0,0,0,pmy_mesh->maxneighbor_);
          neighbor[nneighbor].SetNeighbor(ranklist[fid], nlevel, fid,
            fid-nslist[ranklist[fid]], n, m, 0, neighbor_edge, bufid, tbid, f1, 0);
          bufid++; nneighbor++;
        }
      }
      else {
        int nlevel=neibt->uid.GetLevel();
        int nid=neibt->gid;
        nblevel[1][m+1][n+1]=nlevel;
        int tbid;
        if(nlevel==ll) tbid=FindBufferID(-n,-m,0,0,0,pmy_mesh->maxneighbor_);
        else tbid=FindBufferID(-n,-m,0,myfx3,0,pmy_mesh->maxneighbor_);
        if(nlevel>=ll || (myox1==n && myox2==m)) {
          neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
            nid-nslist[ranklist[nid]], n, m, 0, neighbor_edge, bufid, tbid);
          nneighbor++;
        }
        bufid+=nf1;
      }
    }
  }
  if(block_size.nx3==1) return;
  // x1x3
  for(int m=-1; m<=1; m+=2) {
    for(int n=-1; n<=1; n+=2) {
      neibt=tree.FindNeighbor(uid,n,0,m,block_bcs,nrbx1,nrbx2,nrbx3,pmy_mesh->root_level);
      if(neibt==NULL) { bufid+=nf1; nblevel[m+1][1][n+1]=-2; continue;}
      if(neibt->flag==false) { // finer
        int ff1=1-(n+1)/2; // 0 for outer_x1, 1 for inner_x1
        int ff2=1-(m+1)/2; // 0 for outer_x3, 1 for inner_x3
        nblevel[m+1][1][n+1]=neibt->uid.GetLevel()+1;
        for(int f1=0;f1<nf1;f1++) {
          BlockTree* nf=neibt->GetLeaf(ff1,f1,ff2);
          int fid = nf->gid;
          int nlevel=nf->uid.GetLevel();
          int tbid=FindBufferID(-n,0,-m,0,0,pmy_mesh->maxneighbor_);
          neighbor[nneighbor].SetNeighbor(ranklist[fid], nlevel, fid,
            fid-nslist[ranklist[fid]], n, 0, m, neighbor_edge, bufid, tbid, f1, 0);
          bufid++; nneighbor++;
        }
      }
      else {
        int nlevel=neibt->uid.GetLevel();
        int nid=neibt->gid;
        nblevel[m+1][1][n+1]=nlevel;
        int tbid;
        if(nlevel==ll) tbid=FindBufferID(-n,0,-m,0,0,pmy_mesh->maxneighbor_);
        else tbid=FindBufferID(-n,0,-m,myfx2,0,pmy_mesh->maxneighbor_);
        if(nlevel>=ll || (myox1==n && myox3==m)) {
          neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
            nid-nslist[ranklist[nid]], n, 0, m, neighbor_edge, bufid, tbid);
          nneighbor++;
        }
        bufid+=nf1;
      }
    }
  }
  // x2x3
  for(int m=-1; m<=1; m+=2) {
    for(int n=-1; n<=1; n+=2) {
      neibt=tree.FindNeighbor(uid,0,n,m,block_bcs,nrbx1,nrbx2,nrbx3,pmy_mesh->root_level);
      if(neibt==NULL) { bufid+=nf1; nblevel[m+1][n+1][1]=-2; continue;}
      if(neibt->flag==false) { // finer
        int ff1=1-(n+1)/2; // 0 for outer_x1, 1 for inner_x1
        int ff2=1-(m+1)/2; // 0 for outer_x3, 1 for inner_x3
        nblevel[m+1][n+1][1]=neibt->uid.GetLevel()+1;
        for(int f1=0;f1<nf1;f1++) {
          BlockTree* nf=neibt->GetLeaf(f1,ff1,ff2);
          int fid = nf->gid;
          int nlevel=nf->uid.GetLevel();
          int tbid=FindBufferID(0,-n,-m,0,0,pmy_mesh->maxneighbor_);
          neighbor[nneighbor].SetNeighbor(ranklist[fid], nlevel, fid,
            fid-nslist[ranklist[fid]], 0, n, m, neighbor_edge, bufid, tbid, f1, 0);
          bufid++; nneighbor++;
        }
      }
      else {
        int nlevel=neibt->uid.GetLevel();
        int nid=neibt->gid;
        nblevel[m+1][n+1][1]=nlevel;
        int tbid;
        if(nlevel==ll) tbid=FindBufferID(0,-n,-m,0,0,pmy_mesh->maxneighbor_);
        else tbid=FindBufferID(0,-n,-m,myfx1,0,pmy_mesh->maxneighbor_);
        if(nlevel>=ll || (myox2==n && myox3==m)) {
          neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
            nid-nslist[ranklist[nid]], 0, n, m, neighbor_edge, bufid, tbid);
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
        neibt=tree.FindNeighbor(uid,n,m,l,block_bcs,nrbx1,nrbx2,nrbx3,pmy_mesh->root_level);
        if(neibt==NULL) { bufid++; nblevel[l+1][m+1][n+1]=-2; continue;}
        if(neibt->flag==false) { // finer
          int ff1=1-(n+1)/2; // 0 for outer_x1, 1 for inner_x1
          int ff2=1-(m+1)/2; // 0 for outer_x2, 1 for inner_x2
          int ff3=1-(l+1)/2; // 0 for outer_x3, 1 for inner_x3
          neibt=neibt->GetLeaf(ff1,ff2,ff3);
        }
        int nlevel=neibt->uid.GetLevel();
        nblevel[l+1][m+1][n+1]=nlevel;
        if(nlevel>=ll || (myox1==n && myox2==m && myox3==l)) {
          int nid=neibt->gid;
          int tbid=FindBufferID(-n,-m,-l,0,0,pmy_mesh->maxneighbor_);
          neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
            nid-nslist[ranklist[nid]], n, m, l, neighbor_corner, bufid, tbid);
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
  Real tcons[NFLUID];
  for(int n=0;n<NFLUID;n++) tcons[n]=0.0;
  while(pmb!=NULL) {
    pmb->IntegrateConservative(tcons);
    pmb=pmb->next;
  }

#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE,tcons,NFLUID,MPI_ATHENA_REAL,MPI_SUM,MPI_COMM_WORLD);
#endif

  if(myrank==0) {
    std::cout << "Total Conservative : " ;
    for(int n=0;n<NFLUID;n++)
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
  for(int n=0;n<NFLUID;n++) {
    for(int k=ks;k<=ke;k++) {
      for(int j=js;j<=je;j++) {
        for(int i=is;i<=ie;i++)
          tcons[n]+=pfluid->u(n,k,j,i)*pcoord->GetCellVolume(k,j,i);
      }
    }
  }
  return;
}
