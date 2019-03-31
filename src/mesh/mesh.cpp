//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mesh.cpp
//  \brief implementation of functions in Mesh class

// C headers
// pre-C11: needed before including inttypes.h, else won't define int64_t for C++ code
// #define __STDC_FORMAT_MACROS

// C++ headers
#include <algorithm>  // std::sort()
#include <cinttypes>  // std::int64_t format macro PRId64 for fixed-width integer types
#include <cmath>      // std::abs(), std::pow()
#include <cstdlib>
#include <cstring>    // std::memcpy()
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../fft/athena_fft.hpp"
#include "../fft/turbulence.hpp"
#include "../field/field.hpp"
#include "../field/field_diffusion/field_diffusion.hpp"
#include "../globals.hpp"
#include "../gravity/fft_gravity.hpp"
#include "../gravity/gravity.hpp"
#include "../gravity/mg_gravity.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/hydro_diffusion/hydro_diffusion.hpp"
#include "../multigrid/multigrid.hpp"
#include "../outputs/io_wrapper.hpp"
#include "../parameter_input.hpp"
#include "../reconstruct/reconstruction.hpp"
#include "../utils/buffer_utils.hpp"
#include "mesh.hpp"
#include "mesh_refinement.hpp"
#include "meshblock_tree.hpp"

// MPI/OpenMP header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

//----------------------------------------------------------------------------------------
// Mesh constructor, builds mesh at start of calculation using parameters in input file

Mesh::Mesh(ParameterInput *pin, int mesh_test) {
  std::stringstream msg;
  RegionSize block_size;
  MeshBlock *pfirst{};
  BoundaryFlag block_bcs[6];
  std::int64_t nbmax;
  int dim;

  // mesh test
  if (mesh_test>0) Globals::nranks=mesh_test;

  // read time and cycle limits from input file
  start_time = pin->GetOrAddReal("time", "start_time", 0.0);
  tlim       = pin->GetReal("time", "tlim");
  cfl_number = pin->GetReal("time", "cfl_number");
  ncycle_out = pin->GetOrAddInteger("time", "ncycle_out", 1);
  time = start_time;
  Real real_max = std::numeric_limits<Real>::max();
  dt = dt_diff = (real_max);
  muj = 0.0;
  nuj = 0.0;
  muj_tilde = 0.0;
  nbnew=0; nbdel=0;

  four_pi_G_=0.0, grav_eps_=-1.0, grav_mean_rho_=-1.0;

  turb_flag = 0;

  nlim = pin->GetOrAddInteger("time", "nlim", -1);
  ncycle = 0;
  nint_user_mesh_data_=0;
  nreal_user_mesh_data_=0;
  nuser_history_output_=0;

  // read number of OpenMP threads for mesh
  num_mesh_threads_ = pin->GetOrAddInteger("mesh", "num_threads", 1);
  if (num_mesh_threads_ < 1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Number of OpenMP threads must be >= 1, but num_threads="
        << num_mesh_threads_ << std::endl;
    ATHENA_ERROR(msg);
  }

  // read number of grid cells in root level of mesh from input file.
  mesh_size.nx1 = pin->GetInteger("mesh","nx1");
  if (mesh_size.nx1 < 4) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "In mesh block in input file nx1 must be >= 4, but nx1="
        << mesh_size.nx1 << std::endl;
    ATHENA_ERROR(msg);
  }

  mesh_size.nx2 = pin->GetInteger("mesh","nx2");
  if (mesh_size.nx2 < 1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "In mesh block in input file nx2 must be >= 1, but nx2="
        << mesh_size.nx2 << std::endl;
    ATHENA_ERROR(msg);
  }

  mesh_size.nx3 = pin->GetInteger("mesh","nx3");
  if (mesh_size.nx3 < 1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "In mesh block in input file nx3 must be >= 1, but nx3="
        << mesh_size.nx3 << std::endl;
    ATHENA_ERROR(msg);
  }
  if (mesh_size.nx2 == 1 && mesh_size.nx3 > 1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "In mesh block in input file: nx2=1, nx3=" << mesh_size.nx3
        << ", 2D problems in x1-x3 plane not supported" << std::endl;
    ATHENA_ERROR(msg);
  }

  dim=1;
  if (mesh_size.nx2>1) dim=2;
  if (mesh_size.nx3>1) dim=3;

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
    ATHENA_ERROR(msg);
  }
  if (mesh_size.x2max <= mesh_size.x2min) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Input x2max must be larger than x2min: x2min=" << mesh_size.x2min
        << " x2max=" << mesh_size.x2max << std::endl;
    ATHENA_ERROR(msg);
  }
  if (mesh_size.x3max <= mesh_size.x3min) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Input x3max must be larger than x3min: x3min=" << mesh_size.x3min
        << " x3max=" << mesh_size.x3max << std::endl;
    ATHENA_ERROR(msg);
  }

  // read ratios of grid cell size in each direction
  block_size.x1rat = mesh_size.x1rat = pin->GetOrAddReal("mesh", "x1rat", 1.0);
  block_size.x2rat = mesh_size.x2rat = pin->GetOrAddReal("mesh", "x2rat", 1.0);
  block_size.x3rat = mesh_size.x3rat = pin->GetOrAddReal("mesh", "x3rat", 1.0);

  // read BC flags for each of the 6 boundaries in turn.
  mesh_bcs[BoundaryFace::inner_x1] =
      GetBoundaryFlag(pin->GetOrAddString("mesh", "ix1_bc", "none"));
  mesh_bcs[BoundaryFace::outer_x1] =
      GetBoundaryFlag(pin->GetOrAddString("mesh", "ox1_bc", "none"));
  mesh_bcs[BoundaryFace::inner_x2] =
      GetBoundaryFlag(pin->GetOrAddString("mesh", "ix2_bc", "none"));
  mesh_bcs[BoundaryFace::outer_x2] =
      GetBoundaryFlag(pin->GetOrAddString("mesh", "ox2_bc", "none"));
  mesh_bcs[BoundaryFace::inner_x3] =
      GetBoundaryFlag(pin->GetOrAddString("mesh", "ix3_bc", "none"));
  mesh_bcs[BoundaryFace::outer_x3] =
      GetBoundaryFlag(pin->GetOrAddString("mesh", "ox3_bc", "none"));

  // read MeshBlock parameters
  block_size.nx1 = pin->GetOrAddInteger("meshblock", "nx1", mesh_size.nx1);
  if (dim>=2)
    block_size.nx2 = pin->GetOrAddInteger("meshblock", "nx2", mesh_size.nx2);
  else
    block_size.nx2=mesh_size.nx2;
  if (dim==3)
    block_size.nx3 = pin->GetOrAddInteger("meshblock", "nx3", mesh_size.nx3);
  else
    block_size.nx3=mesh_size.nx3;

  // check consistency of the block and mesh
  if (mesh_size.nx1 % block_size.nx1 != 0
      || mesh_size.nx2 % block_size.nx2 != 0
      || mesh_size.nx3 % block_size.nx3 != 0) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "the Mesh must be evenly divisible by the MeshBlock" << std::endl;
    ATHENA_ERROR(msg);
  }
  if (block_size.nx1 < 4 || (block_size.nx2 < 4 && dim >= 2)
      || (block_size.nx3 < 4 && dim==3)) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "block_size must be larger than or equal to 4 cells." << std::endl;
    ATHENA_ERROR(msg);
  }

  // calculate the number of the blocks
  nrbx1 = mesh_size.nx1/block_size.nx1;
  nrbx2 = mesh_size.nx2/block_size.nx2;
  nrbx3 = mesh_size.nx3/block_size.nx3;
  nbmax = (nrbx1>nrbx2) ? nrbx1:nrbx2;
  nbmax = (nbmax>nrbx3) ? nbmax:nrbx3;

  // initialize user-enrollable functions
  if (mesh_size.x1rat != 1.0) {
    use_uniform_meshgen_fn_[X1DIR]=false;
    MeshGenerator_[X1DIR]=DefaultMeshGeneratorX1;
  } else {
    use_uniform_meshgen_fn_[X1DIR]=true;
    MeshGenerator_[X1DIR]=UniformMeshGeneratorX1;
  }
  if (mesh_size.x2rat != 1.0) {
    use_uniform_meshgen_fn_[X2DIR]=false;
    MeshGenerator_[X2DIR]=DefaultMeshGeneratorX2;
  } else {
    use_uniform_meshgen_fn_[X2DIR]=true;
    MeshGenerator_[X2DIR]=UniformMeshGeneratorX2;
  }
  if (mesh_size.x3rat != 1.0) {
    use_uniform_meshgen_fn_[X3DIR]=false;
    MeshGenerator_[X3DIR]=DefaultMeshGeneratorX3;
  } else {
    use_uniform_meshgen_fn_[X3DIR]=true;
    MeshGenerator_[X3DIR]=UniformMeshGeneratorX3;
  }

  for (int dir=0; dir<6; dir++)
    BoundaryFunction_[dir]=nullptr;
  AMRFlag_=nullptr;
  UserSourceTerm_=nullptr;
  UserTimeStep_=nullptr;
  ViscosityCoeff_=nullptr;
  ConductionCoeff_=nullptr;
  FieldDiffusivity_=nullptr;
  MGBoundaryFunction_[BoundaryFace::inner_x1]=MGPeriodicInnerX1;
  MGBoundaryFunction_[BoundaryFace::outer_x1]=MGPeriodicOuterX1;
  MGBoundaryFunction_[BoundaryFace::inner_x2]=MGPeriodicInnerX2;
  MGBoundaryFunction_[BoundaryFace::outer_x2]=MGPeriodicOuterX2;
  MGBoundaryFunction_[BoundaryFace::inner_x3]=MGPeriodicInnerX3;
  MGBoundaryFunction_[BoundaryFace::outer_x3]=MGPeriodicOuterX3;


  // calculate the logical root level and maximum level
  for (root_level=0; (1<<root_level)<nbmax; root_level++) {}
  current_level=root_level;

  // create the root grid
  tree.CreateRootGrid(nrbx1, nrbx2, nrbx3, root_level);

  // SMR / AMR: create finer grids here
  multilevel=false;
  adaptive=false;
  if (pin->GetOrAddString("mesh", "refinement", "none")=="adaptive")
    adaptive=true, multilevel=true;
  else if (pin->GetOrAddString("mesh", "refinement", "none")=="static")
    multilevel=true;
  if (adaptive==true) {
    max_level = pin->GetOrAddInteger("mesh", "numlevel", 1)+root_level-1;
    if (max_level > 63) {
      msg << "### FATAL ERROR in Mesh constructor" << std::endl
          << "The number of the refinement level must be smaller than "
          << 63-root_level+1 << "." << std::endl;
      ATHENA_ERROR(msg);
    }
  } else {
    max_level = 63;
  }

  if (EOS_TABLE_ENABLED) peos_table = new EosTable(pin);
  InitUserMeshData(pin);

  if (multilevel==true) {
    if (block_size.nx1 % 2==1 || (block_size.nx2 % 2==1 && block_size.nx2>1)
        || (block_size.nx3 % 2==1 && block_size.nx3>1)) {
      msg << "### FATAL ERROR in Mesh constructor" << std::endl
          << "The size of MeshBlock must be divisible by 2 in order to use SMR or AMR."
          << std::endl;
      ATHENA_ERROR(msg);
    }

    InputBlock *pib = pin->pfirst_block;
    while (pib != nullptr) {
      if (pib->block_name.compare(0, 10, "refinement") == 0) {
        RegionSize ref_size;
        ref_size.x1min = pin->GetReal(pib->block_name, "x1min");
        ref_size.x1max = pin->GetReal(pib->block_name, "x1max");
        if (dim>=2) {
          ref_size.x2min = pin->GetReal(pib->block_name, "x2min");
          ref_size.x2max = pin->GetReal(pib->block_name, "x2max");
        } else {
          ref_size.x2min=mesh_size.x2min;
          ref_size.x2max=mesh_size.x2max;
        }
        if (dim>=3) {
          ref_size.x3min = pin->GetReal(pib->block_name, "x3min");
          ref_size.x3max = pin->GetReal(pib->block_name, "x3max");
        } else {
          ref_size.x3min=mesh_size.x3min;
          ref_size.x3max=mesh_size.x3max;
        }
        int ref_lev = pin->GetInteger(pib->block_name, "level");
        int lrlev=ref_lev+root_level;
        if (lrlev>current_level) current_level=lrlev;
        // range check
        if (ref_lev<1) {
          msg << "### FATAL ERROR in Mesh constructor" << std::endl
              << "Refinement level must be larger than 0 (root level = 0)" << std::endl;
          ATHENA_ERROR(msg);
        }
        if (lrlev > max_level) {
          msg << "### FATAL ERROR in Mesh constructor" << std::endl
              << "Refinement level exceeds the maximum level (specify"
              << "maxlevel in <mesh> if adaptive)."
              << std::endl;
          ATHENA_ERROR(msg);
        }
        if (ref_size.x1min > ref_size.x1max || ref_size.x2min > ref_size.x2max
            || ref_size.x3min > ref_size.x3max)  {
          msg << "### FATAL ERROR in Mesh constructor" << std::endl
              << "Invalid refinement region is specified."<<  std::endl;
          ATHENA_ERROR(msg);
        }
        if (ref_size.x1min < mesh_size.x1min || ref_size.x1max > mesh_size.x1max
            || ref_size.x2min < mesh_size.x2min || ref_size.x2max > mesh_size.x2max
            || ref_size.x3min < mesh_size.x3min || ref_size.x3max > mesh_size.x3max) {
          msg << "### FATAL ERROR in Mesh constructor" << std::endl
              << "Refinement region must be smaller than the whole mesh." << std::endl;
          ATHENA_ERROR(msg);
        }
        // find the logical range in the ref_level
        // note: if this is too slow, this should be replaced with bi-section search.
        std::int64_t lx1min=0, lx1max=0, lx2min=0, lx2max=0, lx3min=0, lx3max=0;
        std::int64_t lxmax=nrbx1*(1LL<<ref_lev);
        for (lx1min=0; lx1min<lxmax; lx1min++) {
          Real rx=ComputeMeshGeneratorX(lx1min+1, lxmax, use_uniform_meshgen_fn_[X1DIR]);
          if (MeshGenerator_[X1DIR](rx, mesh_size) > ref_size.x1min)
            break;
        }
        for (lx1max=lx1min; lx1max<lxmax; lx1max++) {
          Real rx=ComputeMeshGeneratorX(lx1max+1, lxmax, use_uniform_meshgen_fn_[X1DIR]);
          if (MeshGenerator_[X1DIR](rx, mesh_size) >= ref_size.x1max)
            break;
        }
        if (lx1min % 2==1) lx1min--;
        if (lx1max % 2==0) lx1max++;
        if (dim>=2) { // 2D or 3D
          lxmax=nrbx2*(1LL<<ref_lev);
          for (lx2min=0; lx2min<lxmax; lx2min++) {
            Real rx=ComputeMeshGeneratorX(lx2min+1,lxmax,use_uniform_meshgen_fn_[X2DIR]);
            if (MeshGenerator_[X2DIR](rx, mesh_size) > ref_size.x2min)
              break;
          }
          for (lx2max=lx2min; lx2max<lxmax; lx2max++) {
            Real rx=ComputeMeshGeneratorX(lx2max+1,lxmax,use_uniform_meshgen_fn_[X2DIR]);
            if (MeshGenerator_[X2DIR](rx, mesh_size) >= ref_size.x2max)
              break;
          }
          if (lx2min % 2==1) lx2min--;
          if (lx2max % 2==0) lx2max++;
        }
        if (dim==3) { // 3D
          lxmax=nrbx3*(1LL<<ref_lev);
          for (lx3min=0; lx3min<lxmax; lx3min++) {
            Real rx=ComputeMeshGeneratorX(lx3min+1,lxmax,use_uniform_meshgen_fn_[X3DIR]);
            if (MeshGenerator_[X3DIR](rx, mesh_size) > ref_size.x3min)
              break;
          }
          for (lx3max=lx3min; lx3max<lxmax; lx3max++) {
            Real rx=ComputeMeshGeneratorX(lx3max+1,lxmax,use_uniform_meshgen_fn_[X3DIR]);
            if (MeshGenerator_[X3DIR](rx, mesh_size) >= ref_size.x3max)
              break;
          }
          if (lx3min % 2==1) lx3min--;
          if (lx3max % 2==0) lx3max++;
        }
        // create the finest level
        if (dim==1) {
          for (std::int64_t i=lx1min; i<lx1max; i+=2) {
            LogicalLocation nloc;
            nloc.level=lrlev, nloc.lx1=i, nloc.lx2=0, nloc.lx3=0;
            int nnew;
            tree.AddMeshBlock(tree, nloc, dim, mesh_bcs, nrbx1, nrbx2, nrbx3, root_level,
                              nnew);
          }
        }
        if (dim==2) {
          for (std::int64_t j=lx2min; j<lx2max; j+=2) {
            for (std::int64_t i=lx1min; i<lx1max; i+=2) {
              LogicalLocation nloc;
              nloc.level=lrlev, nloc.lx1=i, nloc.lx2=j, nloc.lx3=0;
              int nnew;
              tree.AddMeshBlock(tree, nloc, dim, mesh_bcs, nrbx1, nrbx2, nrbx3,
                                root_level, nnew);
            }
          }
        }
        if (dim==3) {
          for (std::int64_t k=lx3min; k<lx3max; k+=2) {
            for (std::int64_t j=lx2min; j<lx2max; j+=2) {
              for (std::int64_t i=lx1min; i<lx1max; i+=2) {
                LogicalLocation nloc;
                nloc.level=lrlev, nloc.lx1=i, nloc.lx2=j, nloc.lx3=k;
                int nnew;
                tree.AddMeshBlock(tree, nloc, dim, mesh_bcs, nrbx1, nrbx2, nrbx3,
                                  root_level, nnew);
              }
            }
          }
        }
      }
      pib = pib->pnext;
    }
  }

  // initial mesh hierarchy construction is completed here

  tree.CountMeshBlock(nbtotal);
  loclist=new LogicalLocation[nbtotal];
  tree.GetMeshBlockList(loclist,nullptr,nbtotal);

#ifdef MPI_PARALLEL
  // check if there are sufficient blocks
  if (nbtotal < Globals::nranks) {
    if (mesh_test==0) {
      msg << "### FATAL ERROR in Mesh constructor" << std::endl
          << "Too few mesh blocks: nbtotal ("<< nbtotal <<") < nranks ("
          << Globals::nranks << ")" << std::endl;
      ATHENA_ERROR(msg);
    } else { // test
      std::cout << "### Warning in Mesh constructor" << std::endl
                << "Too few mesh blocks: nbtotal ("<< nbtotal <<") < nranks ("
                << Globals::nranks << ")" << std::endl;
    }
  }
#endif

  ranklist=new int[nbtotal];
  nslist=new int[Globals::nranks];
  nblist=new int[Globals::nranks];
  costlist=new Real[nbtotal];
  if (adaptive==true) { // allocate arrays for AMR
    nref = new int[Globals::nranks];
    nderef = new int[Globals::nranks];
    rdisp = new int[Globals::nranks];
    ddisp = new int[Globals::nranks];
    bnref = new int[Globals::nranks];
    bnderef = new int[Globals::nranks];
    brdisp = new int[Globals::nranks];
    bddisp = new int[Globals::nranks];
  }

  // initialize cost array with the simplest estimate; all the blocks are equal
  for (int i=0; i<nbtotal; i++) costlist[i]=1.0;

  LoadBalance(costlist, ranklist, nslist, nblist, nbtotal);

  // Output some diagnostic information to terminal

  // Output MeshBlock list and quit (mesh test only); do not create meshes
  if (mesh_test>0) {
    if (Globals::my_rank==0) OutputMeshStructure(dim);
    return;
  }

  // set gravity flag
  gflag=0;
  if (SELF_GRAVITY_ENABLED) gflag=1;
  //  if (SELF_GRAVITY_ENABLED==2 && ...) // independent allocation
  //    gflag=2;

  // create MeshBlock list for this process
  int nbs=nslist[Globals::my_rank];
  int nbe=nbs+nblist[Globals::my_rank]-1;
  // create MeshBlock list for this process
  for (int i=nbs; i<=nbe; i++) {
    SetBlockSizeAndBoundaries(loclist[i], block_size, block_bcs);
    // create a block and add into the link list
    if (i==nbs) {
      pblock = new MeshBlock(i, i-nbs, loclist[i], block_size, block_bcs, this,
                             pin, gflag);
      pfirst = pblock;
    } else {
      pblock->next = new MeshBlock(i, i-nbs, loclist[i], block_size, block_bcs,
                                   this, pin, gflag);
      pblock->next->prev = pblock;
      pblock = pblock->next;
    }
    pblock->pbval->SearchAndSetNeighbors(tree, ranklist, nslist);
  }
  pblock = pfirst;

  if (SELF_GRAVITY_ENABLED==1)
    pfgrd = new FFTGravityDriver(this, pin);
  else if (SELF_GRAVITY_ENABLED==2)
    pmgrd = new MGGravityDriver(this, MGBoundaryFunction_, pin);

  if (turb_flag > 0)
    ptrbd = new TurbulenceDriver(this, pin);
}

//----------------------------------------------------------------------------------------
// Mesh constructor for restarts. Load the restart file

Mesh::Mesh(ParameterInput *pin, IOWrapper& resfile, int mesh_test) {
  std::stringstream msg;
  RegionSize block_size;
  BoundaryFlag block_bcs[6];
  MeshBlock *pfirst{};
  IOWrapperSizeT *offset{};
  IOWrapperSizeT datasize, listsize, headeroffset;

  // mesh test
  if (mesh_test>0) Globals::nranks=mesh_test;

  // read time and cycle limits from input file
  start_time = pin->GetOrAddReal("time","start_time",0.0);
  tlim       = pin->GetReal("time","tlim");
  ncycle_out = pin->GetOrAddInteger("time","ncycle_out",1);
  nlim = pin->GetOrAddInteger("time","nlim",-1);
  nint_user_mesh_data_=0;
  nreal_user_mesh_data_=0;
  nuser_history_output_=0;

  four_pi_G_=0.0, grav_eps_=-1.0, grav_mean_rho_=-1.0;

  turb_flag = 0;

  nbnew=0; nbdel=0;

  // read number of OpenMP threads for mesh
  num_mesh_threads_ = pin->GetOrAddInteger("mesh","num_threads",1);
  if (num_mesh_threads_ < 1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Number of OpenMP threads must be >= 1, but num_threads="
        << num_mesh_threads_ << std::endl;
    ATHENA_ERROR(msg);
  }

  // read BC flags for each of the 6 boundaries
  mesh_bcs[BoundaryFace::inner_x1] =
      GetBoundaryFlag(pin->GetOrAddString("mesh", "ix1_bc", "none"));
  mesh_bcs[BoundaryFace::outer_x1] =
      GetBoundaryFlag(pin->GetOrAddString("mesh", "ox1_bc", "none"));
  mesh_bcs[BoundaryFace::inner_x2] =
      GetBoundaryFlag(pin->GetOrAddString("mesh", "ix2_bc", "none"));
  mesh_bcs[BoundaryFace::outer_x2] =
      GetBoundaryFlag(pin->GetOrAddString("mesh", "ox2_bc", "none"));
  mesh_bcs[BoundaryFace::inner_x3] =
      GetBoundaryFlag(pin->GetOrAddString("mesh", "ix3_bc", "none"));
  mesh_bcs[BoundaryFace::outer_x3] =
      GetBoundaryFlag(pin->GetOrAddString("mesh", "ox3_bc", "none"));

  // get the end of the header
  headeroffset=resfile.GetPosition();
  // read the restart file
  // the file is already open and the pointer is set to after <par_end>
  IOWrapperSizeT headersize = sizeof(int)*3+sizeof(Real)*2
                              + sizeof(RegionSize)+sizeof(IOWrapperSizeT);
  char *headerdata = new char[headersize];
  if (Globals::my_rank==0) { // the master process reads the header data
    if (resfile.Read(headerdata, 1, headersize) != headersize) {
      msg << "### FATAL ERROR in Mesh constructor" << std::endl
          << "The restart file is broken." << std::endl;
      ATHENA_ERROR(msg);
    }
  }
#ifdef MPI_PARALLEL
  // then broadcast the header data
  MPI_Bcast(headerdata, headersize, MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
  IOWrapperSizeT hdos = 0;
  std::memcpy(&nbtotal, &(headerdata[hdos]), sizeof(int));
  hdos += sizeof(int);
  std::memcpy(&root_level, &(headerdata[hdos]), sizeof(int));
  hdos += sizeof(int);
  current_level=root_level;
  std::memcpy(&mesh_size, &(headerdata[hdos]), sizeof(RegionSize));
  hdos += sizeof(RegionSize);
  std::memcpy(&time, &(headerdata[hdos]), sizeof(Real));
  hdos += sizeof(Real);
  std::memcpy(&dt, &(headerdata[hdos]), sizeof(Real));
  hdos += sizeof(Real);
  std::memcpy(&ncycle, &(headerdata[hdos]), sizeof(int));
  hdos += sizeof(int);
  std::memcpy(&datasize, &(headerdata[hdos]), sizeof(IOWrapperSizeT));
  hdos += sizeof(IOWrapperSizeT);   // KGF: this updated value is never used

  delete [] headerdata;

  int dim = 1;
  if (mesh_size.nx2 > 1) dim=2;
  if (mesh_size.nx3 > 1) dim=3;

  // initialize
  loclist=new LogicalLocation[nbtotal];
  offset=new IOWrapperSizeT[nbtotal];
  costlist=new Real[nbtotal];
  ranklist=new int[nbtotal];
  nslist=new int[Globals::nranks];
  nblist=new int[Globals::nranks];

  block_size.nx1 = pin->GetOrAddInteger("meshblock","nx1",mesh_size.nx1);
  block_size.nx2 = pin->GetOrAddInteger("meshblock","nx2",mesh_size.nx2);
  block_size.nx3 = pin->GetOrAddInteger("meshblock","nx3",mesh_size.nx3);

  // calculate the number of the blocks
  nrbx1=mesh_size.nx1/block_size.nx1;
  nrbx2=mesh_size.nx2/block_size.nx2;
  nrbx3=mesh_size.nx3/block_size.nx3;

  // initialize user-enrollable functions
  if (mesh_size.x1rat != 1.0) {
    use_uniform_meshgen_fn_[X1DIR]=false;
    MeshGenerator_[X1DIR]=DefaultMeshGeneratorX1;
  } else {
    use_uniform_meshgen_fn_[X1DIR]=true;
    MeshGenerator_[X1DIR]=UniformMeshGeneratorX1;
  }
  if (mesh_size.x2rat != 1.0) {
    use_uniform_meshgen_fn_[X2DIR]=false;
    MeshGenerator_[X2DIR]=DefaultMeshGeneratorX2;
  } else {
    use_uniform_meshgen_fn_[X2DIR]=true;
    MeshGenerator_[X2DIR]=UniformMeshGeneratorX2;
  }
  if (mesh_size.x3rat != 1.0) {
    use_uniform_meshgen_fn_[X3DIR]=false;
    MeshGenerator_[X3DIR]=DefaultMeshGeneratorX3;
  } else {
    use_uniform_meshgen_fn_[X3DIR]=true;
    MeshGenerator_[X3DIR]=UniformMeshGeneratorX3;
  }

  for (int dir=0; dir<6; dir++)
    BoundaryFunction_[dir]=nullptr;
  AMRFlag_=nullptr;
  UserSourceTerm_=nullptr;
  UserTimeStep_=nullptr;
  ViscosityCoeff_=nullptr;
  ConductionCoeff_=nullptr;
  FieldDiffusivity_=nullptr;
  MGBoundaryFunction_[BoundaryFace::inner_x1]=MGPeriodicInnerX1;
  MGBoundaryFunction_[BoundaryFace::outer_x1]=MGPeriodicOuterX1;
  MGBoundaryFunction_[BoundaryFace::inner_x2]=MGPeriodicInnerX2;
  MGBoundaryFunction_[BoundaryFace::outer_x2]=MGPeriodicOuterX2;
  MGBoundaryFunction_[BoundaryFace::inner_x3]=MGPeriodicInnerX3;
  MGBoundaryFunction_[BoundaryFace::outer_x3]=MGPeriodicOuterX3;

  multilevel=false;
  adaptive=false;
  if (pin->GetOrAddString("mesh","refinement","none")=="adaptive")
    adaptive=true, multilevel=true;
  else if (pin->GetOrAddString("mesh","refinement","none")=="static")
    multilevel=true;
  if (adaptive==true) {
    max_level = pin->GetOrAddInteger("mesh","numlevel",1)+root_level-1;
    if (max_level > 63) {
      msg << "### FATAL ERROR in Mesh constructor" << std::endl
          << "The number of the refinement level must be smaller than "
          << 63-root_level+1 << "." << std::endl;
      ATHENA_ERROR(msg);
    }
  } else {
    max_level = 63;
  }

  if (EOS_TABLE_ENABLED) peos_table = new EosTable(pin);
  InitUserMeshData(pin);

  // read user Mesh data
  IOWrapperSizeT udsize = 0;
  for (int n=0; n<nint_user_mesh_data_; n++)
    udsize += iuser_mesh_data[n].GetSizeInBytes();
  for (int n=0; n<nreal_user_mesh_data_; n++)
    udsize += ruser_mesh_data[n].GetSizeInBytes();
  if (udsize != 0) {
    char *userdata = new char[udsize];
    if (Globals::my_rank==0) { // only the master process reads the ID list
      if (resfile.Read(userdata,1,udsize) != udsize) {
        msg << "### FATAL ERROR in Mesh constructor" << std::endl
            << "The restart file is broken." << std::endl;
        ATHENA_ERROR(msg);
      }
    }
#ifdef MPI_PARALLEL
    // then broadcast the ID list
    MPI_Bcast(userdata, udsize, MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

    IOWrapperSizeT udoffset=0;
    for (int n=0; n<nint_user_mesh_data_; n++) {
      std::memcpy(iuser_mesh_data[n].data(), &(userdata[udoffset]),
                  iuser_mesh_data[n].GetSizeInBytes());
      udoffset += iuser_mesh_data[n].GetSizeInBytes();
    }
    for (int n=0; n<nreal_user_mesh_data_; n++) {
      std::memcpy(ruser_mesh_data[n].data(), &(userdata[udoffset]),
                  ruser_mesh_data[n].GetSizeInBytes());
      udoffset += ruser_mesh_data[n].GetSizeInBytes();
    }
    delete [] userdata;
  }

  // read the ID list
  listsize=sizeof(LogicalLocation)+sizeof(Real);
  //allocate the idlist buffer
  char *idlist = new char[listsize*nbtotal];
  if (Globals::my_rank==0) { // only the master process reads the ID list
    if (resfile.Read(idlist,listsize,nbtotal) != static_cast<unsigned int>(nbtotal)) {
      msg << "### FATAL ERROR in Mesh constructor" << std::endl
          << "The restart file is broken." << std::endl;
      ATHENA_ERROR(msg);
    }
  }
#ifdef MPI_PARALLEL
  // then broadcast the ID list
  MPI_Bcast(idlist, listsize*nbtotal, MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

  int os=0;
  for (int i=0; i<nbtotal; i++) {
    std::memcpy(&(loclist[i]), &(idlist[os]), sizeof(LogicalLocation));
    os += sizeof(LogicalLocation);
    std::memcpy(&(costlist[i]), &(idlist[os]), sizeof(Real));
    os += sizeof(Real);
    if (loclist[i].level>current_level) current_level=loclist[i].level;
  }
  delete [] idlist;

  // calculate the header offset and seek
  headeroffset += headersize+udsize+listsize*nbtotal;
  if (Globals::my_rank != 0)
    resfile.Seek(headeroffset);

  // rebuild the Block Tree
  for (int i=0; i<nbtotal; i++)
    tree.AddMeshBlockWithoutRefine(loclist[i],nrbx1,nrbx2,nrbx3,root_level);
  int nnb;
  // check the tree structure, and assign GID
  tree.GetMeshBlockList(loclist, nullptr, nnb);
  if (nnb != nbtotal) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Tree reconstruction failed. The total numbers of the blocks do not match. ("
        << nbtotal << " != " << nnb << ")" << std::endl;
    ATHENA_ERROR(msg);
  }

#ifdef MPI_PARALLEL
  if (nbtotal < Globals::nranks) {
    if (mesh_test==0) {
      msg << "### FATAL ERROR in Mesh constructor" << std::endl
          << "Too few mesh blocks: nbtotal ("<< nbtotal <<") < nranks ("
          << Globals::nranks << ")" << std::endl;
      ATHENA_ERROR(msg);
    } else { // test
      std::cout << "### Warning in Mesh constructor" << std::endl
                << "Too few mesh blocks: nbtotal ("<< nbtotal <<") < nranks ("
                << Globals::nranks << ")" << std::endl;
      delete [] offset;
      return;
    }
  }
#endif

  if (adaptive==true) { // allocate arrays for AMR
    nref = new int[Globals::nranks];
    nderef = new int[Globals::nranks];
    rdisp = new int[Globals::nranks];
    ddisp = new int[Globals::nranks];
    bnref = new int[Globals::nranks];
    bnderef = new int[Globals::nranks];
    brdisp = new int[Globals::nranks];
    bddisp = new int[Globals::nranks];
  }

  LoadBalance(costlist, ranklist, nslist, nblist, nbtotal);

  // Output MeshBlock list and quit (mesh test only); do not create meshes
  if (mesh_test>0) {
    if (Globals::my_rank==0) OutputMeshStructure(dim);
    delete [] offset;
    return;
  }

  // set gravity flag
  gflag=0;
  if (SELF_GRAVITY_ENABLED) gflag=1;
  //  if (SELF_GRAVITY_ENABLED==2 && ...) // independent allocation
  //    gflag=2;

  // allocate data buffer
  int nb=nblist[Globals::my_rank];
  int nbs=nslist[Globals::my_rank];
  int nbe=nbs+nb-1;
  char *mbdata = new char[datasize*nb];
  // load MeshBlocks (parallel)
  if (resfile.Read_at_all(mbdata, datasize, nb, headeroffset+nbs*datasize) !=
      static_cast<unsigned int>(nb)) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "The restart file is broken or input parameters are inconsistent."
        << std::endl;
    ATHENA_ERROR(msg);
  }
  for (int i=nbs; i<=nbe; i++) {
    // Match fixed-width integer precision of IOWrapperSizeT datasize
    std::uint64_t buff_os = datasize * (i-nbs);
    SetBlockSizeAndBoundaries(loclist[i], block_size, block_bcs);
    // create a block and add into the link list
    if (i==nbs) {
      pblock = new MeshBlock(i, i-nbs, this, pin, loclist[i], block_size,
                             block_bcs, costlist[i], mbdata+buff_os, gflag);
      pfirst = pblock;
    } else {
      pblock->next = new MeshBlock(i, i-nbs, this, pin, loclist[i], block_size,
                                   block_bcs, costlist[i], mbdata+buff_os, gflag);
      pblock->next->prev = pblock;
      pblock = pblock->next;
    }
    pblock->pbval->SearchAndSetNeighbors(tree, ranklist, nslist);
  }
  pblock = pfirst;
  delete [] mbdata;
  // check consistency
  if (datasize != pblock->GetBlockSizeInBytes()) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "The restart file is broken or input parameters are inconsistent."
        << std::endl;
    ATHENA_ERROR(msg);
  }

  // clean up
  delete [] offset;

  if (SELF_GRAVITY_ENABLED==1)
    pfgrd = new FFTGravityDriver(this, pin);
  else if (SELF_GRAVITY_ENABLED==2)
    pmgrd = new MGGravityDriver(this, MGBoundaryFunction_, pin);

  if (turb_flag > 0)
    ptrbd = new TurbulenceDriver(this, pin);
}

//----------------------------------------------------------------------------------------
// destructor

Mesh::~Mesh() {
  while (pblock->prev != nullptr) // should not be true
    delete pblock->prev;
  while (pblock->next != nullptr)
    delete pblock->next;
  delete pblock;
  delete [] nslist;
  delete [] nblist;
  delete [] ranklist;
  delete [] costlist;
  delete [] loclist;
  if (SELF_GRAVITY_ENABLED==1) delete pfgrd;
  else if (SELF_GRAVITY_ENABLED==2) delete pmgrd;
  if (turb_flag > 0) delete ptrbd;
  if (adaptive==true) { // deallocate arrays for AMR
    delete [] nref;
    delete [] nderef;
    delete [] rdisp;
    delete [] ddisp;
    delete [] bnref;
    delete [] bnderef;
    delete [] brdisp;
    delete [] bddisp;
  }
  // delete user Mesh data
  for (int n=0; n<nreal_user_mesh_data_; n++)
    ruser_mesh_data[n].DeleteAthenaArray();
  if (nreal_user_mesh_data_>0) delete [] ruser_mesh_data;
  for (int n=0; n<nint_user_mesh_data_; n++)
    iuser_mesh_data[n].DeleteAthenaArray();
  if (nint_user_mesh_data_>0) delete [] iuser_mesh_data;
  if (EOS_TABLE_ENABLED) delete peos_table;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::OutputMeshStructure(int dim)
//  \brief print the mesh structure information

void Mesh::OutputMeshStructure(int dim) {
  RegionSize block_size;
  BoundaryFlag block_bcs[6];
  FILE *fp = nullptr;

  // open 'mesh_structure.dat' file
  if (dim>=2) {
    if ((fp = std::fopen("mesh_structure.dat","wb")) == nullptr) {
      std::cout << "### ERROR in function Mesh::OutputMeshStructure" << std::endl
                << "Cannot open mesh_structure.dat" << std::endl;
      return;
    }
  }

  // Write overall Mesh structure to stdout and file
  std::cout << std::endl;
  std::cout << "Root grid = " << nrbx1 << " x " << nrbx2 << " x " << nrbx3
            << " MeshBlocks" << std::endl;
  std::cout << "Total number of MeshBlocks = " << nbtotal << std::endl;
  std::cout << "Number of physical refinement levels = "
            << (current_level - root_level) << std::endl;
  std::cout << "Number of logical  refinement levels = " << current_level << std::endl;

  // compute/output number of blocks per level, and cost per level
  int *nb_per_plevel = new int[max_level];
  int *cost_per_plevel = new int[max_level];
  for (int i=0; i<=max_level; ++i) {
    nb_per_plevel[i]=0;
    cost_per_plevel[i]=0;
  }
  for (int i=0; i<nbtotal; i++) {
    nb_per_plevel[(loclist[i].level - root_level)]++;
    cost_per_plevel[(loclist[i].level - root_level)] += costlist[i];
  }
  for (int i=root_level; i<=max_level; i++) {
    if (nb_per_plevel[i-root_level] != 0) {
      std::cout << "  Physical level = " << i-root_level << " (logical level = " << i
                << "): " << nb_per_plevel[i-root_level] << " MeshBlocks, cost = "
                << cost_per_plevel[i-root_level] <<  std::endl;
    }
  }

  // compute/output number of blocks per rank, and cost per rank
  std::cout << "Number of parallel ranks = " << Globals::nranks << std::endl;
  int *nb_per_rank = new int[Globals::nranks];
  int *cost_per_rank = new int[Globals::nranks];
  for (int i=0; i<Globals::nranks; ++i) {
    nb_per_rank[i]=0;
    cost_per_rank[i]=0;
  }
  for (int i=0; i<nbtotal; i++) {
    nb_per_rank[ranklist[i]]++;
    cost_per_rank[ranklist[i]] += costlist[i];
  }
  for (int i=0; i<Globals::nranks; ++i) {
    std::cout << "  Rank = " << i << ": " << nb_per_rank[i] <<" MeshBlocks, cost = "
              << cost_per_rank[i] << std::endl;
  }

  // output relative size/locations of meshblock to file, for plotting
  Real real_max = std::numeric_limits<Real>::max();
  Real mincost=real_max, maxcost=0.0, totalcost=0.0;
  for (int i=root_level; i<=max_level; i++) {
    for (int j=0; j<nbtotal; j++) {
      if (loclist[j].level==i) {
        SetBlockSizeAndBoundaries(loclist[j], block_size, block_bcs);
        std::int64_t &lx1=loclist[j].lx1;
        std::int64_t &lx2=loclist[j].lx2;
        std::int64_t &lx3=loclist[j].lx3;
        int &ll=loclist[j].level;
        mincost=std::min(mincost,costlist[i]);
        maxcost=std::max(maxcost,costlist[i]);
        totalcost+=costlist[i];
        std::fprintf(fp,"#MeshBlock %d on rank=%d with cost=%g\n", j, ranklist[j],
                     costlist[j]);
        std::fprintf(
            fp, "#  Logical level %d, location = (%" PRId64 " %" PRId64 " %" PRId64")\n",
            ll, lx1, lx2, lx3);
        if (dim==2) {
          std::fprintf(fp, "%g %g\n", block_size.x1min, block_size.x2min);
          std::fprintf(fp, "%g %g\n", block_size.x1max, block_size.x2min);
          std::fprintf(fp, "%g %g\n", block_size.x1max, block_size.x2max);
          std::fprintf(fp, "%g %g\n", block_size.x1min, block_size.x2max);
          std::fprintf(fp, "%g %g\n", block_size.x1min, block_size.x2min);
          std::fprintf(fp, "\n\n");
        }
        if (dim==3) {
          std::fprintf(fp, "%g %g %g\n", block_size.x1min, block_size.x2min,
                       block_size.x3min);
          std::fprintf(fp, "%g %g %g\n", block_size.x1max, block_size.x2min,
                       block_size.x3min);
          std::fprintf(fp, "%g %g %g\n", block_size.x1max, block_size.x2max,
                       block_size.x3min);
          std::fprintf(fp, "%g %g %g\n", block_size.x1min, block_size.x2max,
                       block_size.x3min);
          std::fprintf(fp, "%g %g %g\n", block_size.x1min, block_size.x2min,
                       block_size.x3min);
          std::fprintf(fp, "%g %g %g\n", block_size.x1min, block_size.x2min,
                       block_size.x3max);
          std::fprintf(fp, "%g %g %g\n", block_size.x1max, block_size.x2min,
                       block_size.x3max);
          std::fprintf(fp, "%g %g %g\n", block_size.x1max, block_size.x2min,
                       block_size.x3min);
          std::fprintf(fp, "%g %g %g\n", block_size.x1max, block_size.x2min,
                       block_size.x3max);
          std::fprintf(fp, "%g %g %g\n", block_size.x1max, block_size.x2max,
                       block_size.x3max);
          std::fprintf(fp, "%g %g %g\n", block_size.x1max, block_size.x2max,
                       block_size.x3min);
          std::fprintf(fp, "%g %g %g\n", block_size.x1max, block_size.x2max,
                       block_size.x3max);
          std::fprintf(fp, "%g %g %g\n", block_size.x1min, block_size.x2max,
                       block_size.x3max);
          std::fprintf(fp, "%g %g %g\n", block_size.x1min, block_size.x2max,
                       block_size.x3min);
          std::fprintf(fp, "%g %g %g\n", block_size.x1min, block_size.x2max,
                       block_size.x3max);
          std::fprintf(fp, "%g %g %g\n", block_size.x1min, block_size.x2min,
                       block_size.x3max);
          std::fprintf(fp, "%g %g %g\n", block_size.x1min, block_size.x2min,
                       block_size.x3min);
          std::fprintf(fp, "\n\n");
        }
      }
    }
  }

  // close file, final outputs
  if (dim>=2) std::fclose(fp);
  std::cout << "Load Balancing:" << std::endl;
  std::cout << "  Minimum cost = " << mincost << ", Maximum cost = " << maxcost
            << ", Average cost = " << totalcost/nbtotal << std::endl << std::endl;
  std::cout << "See the 'mesh_structure.dat' file for a complete list"
            << " of MeshBlocks." << std::endl;
  std::cout << "Use 'python ../vis/python/plot_mesh.py' or gnuplot"
            << " to visualize mesh structure." << std::endl << std::endl;

  delete [] nb_per_plevel;
  delete [] cost_per_plevel;
  delete [] nb_per_rank;
  delete [] cost_per_rank;

  return;
}

//----------------------------------------------------------------------------------------
// \!fn void Mesh::NewTimeStep()
// \brief function that loops over all MeshBlocks and find new timestep
//        this assumes that phydro->NewBlockTimeStep is already called

void Mesh::NewTimeStep() {
  MeshBlock *pmb = pblock;

  dt_diff=dt=static_cast<Real>(2.0)*dt;

  while (pmb != nullptr)  {
    dt = std::min(dt,pmb->new_block_dt_);
    dt_diff  = std::min(dt_diff, pmb->new_block_dt_diff_);
    pmb = pmb->next;
  }

#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE,&dt,1,MPI_ATHENA_REAL,MPI_MIN,MPI_COMM_WORLD);
  if (STS_ENABLED)
    MPI_Allreduce(MPI_IN_PLACE,&dt_diff,1,MPI_ATHENA_REAL,MPI_MIN,MPI_COMM_WORLD);
#endif

  if (time < tlim && tlim-time < dt) // timestep would take us past desired endpoint
    dt = tlim-time;

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollUserBoundaryFunction(BoundaryFace dir, BValHydro my_bc)
//  \brief Enroll a user-defined boundary function

void Mesh::EnrollUserBoundaryFunction(BoundaryFace dir, BValFunc my_bc) {
  std::stringstream msg;
  if (dir<0 || dir>5) {
    msg << "### FATAL ERROR in EnrollBoundaryCondition function" << std::endl
        << "dirName = " << dir << " not valid" << std::endl;
    ATHENA_ERROR(msg);
  }
  if (mesh_bcs[dir] != BoundaryFlag::user) {
    msg << "### FATAL ERROR in EnrollUserBoundaryFunction" << std::endl
        << "The boundary condition flag must be set to the string 'user' in the "
        << " <mesh> block in the input file to use user-enrolled BCs" << std::endl;
    ATHENA_ERROR(msg);
  }
  BoundaryFunction_[static_cast<int>(dir)]=my_bc;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollUserMGBoundaryFunction(BoundaryFace dir
//                                              MGBoundaryFunc my_bc)
//  \brief Enroll a user-defined Multigrid boundary function

void Mesh::EnrollUserMGBoundaryFunction(BoundaryFace dir, MGBoundaryFunc my_bc) {
  std::stringstream msg;
  if (dir<0 || dir>5) {
    msg << "### FATAL ERROR in EnrollBoundaryCondition function" << std::endl
        << "dirName = " << dir << " not valid" << std::endl;
    ATHENA_ERROR(msg);
  }
  MGBoundaryFunction_[static_cast<int>(dir)]=my_bc;
  return;
}

// DEPRECATED(felker): provide trivial overloads for old-style BoundaryFace enum argument
void Mesh::EnrollUserBoundaryFunction(int dir, BValFunc my_bc) {
  EnrollUserBoundaryFunction(static_cast<BoundaryFace>(dir), my_bc);
  return;
}

void Mesh::EnrollUserMGBoundaryFunction(int dir, MGBoundaryFunc my_bc) {
  EnrollUserMGBoundaryFunction(static_cast<BoundaryFace>(dir), my_bc);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollUserRefinementCondition(AMRFlagFunc amrflag)
//  \brief Enroll a user-defined function for checking refinement criteria

void Mesh::EnrollUserRefinementCondition(AMRFlagFunc amrflag) {
  if (adaptive==true)
    AMRFlag_=amrflag;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollUserMeshGenerator(CoordinateDirection,MeshGenFunc my_mg)
//  \brief Enroll a user-defined function for Mesh generation

void Mesh::EnrollUserMeshGenerator(CoordinateDirection dir, MeshGenFunc my_mg) {
  std::stringstream msg;
  if (dir<0 || dir>=3) {
    msg << "### FATAL ERROR in EnrollUserMeshGenerator function" << std::endl
        << "dirName = " << dir << " not valid" << std::endl;
    ATHENA_ERROR(msg);
  }
  if (dir == X1DIR && mesh_size.x1rat > 0.0) {
    msg << "### FATAL ERROR in EnrollUserMeshGenerator function" << std::endl
        << "x1rat = " << mesh_size.x1rat <<
        " must be negative for user-defined mesh generator in X1DIR " << std::endl;
    ATHENA_ERROR(msg);
  }
  if (dir == X2DIR && mesh_size.x2rat > 0.0) {
    msg << "### FATAL ERROR in EnrollUserMeshGenerator function" << std::endl
        << "x2rat = " << mesh_size.x2rat <<
        " must be negative for user-defined mesh generator in X2DIR " << std::endl;
    ATHENA_ERROR(msg);
  }
  if (dir == X3DIR && mesh_size.x3rat > 0.0) {
    msg << "### FATAL ERROR in EnrollUserMeshGenerator function" << std::endl
        << "x3rat = " << mesh_size.x3rat <<
        " must be negative for user-defined mesh generator in X3DIR " << std::endl;
    ATHENA_ERROR(msg);
  }
  use_uniform_meshgen_fn_[dir]=false;
  MeshGenerator_[dir]=my_mg;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollUserExplicitSourceFunction(SrcTermFunc my_func)
//  \brief Enroll a user-defined source function

void Mesh::EnrollUserExplicitSourceFunction(SrcTermFunc my_func) {
  UserSourceTerm_ = my_func;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollUserTimeStepFunction(TimeStepFunc my_func)
//  \brief Enroll a user-defined time step function

void Mesh::EnrollUserTimeStepFunction(TimeStepFunc my_func) {
  UserTimeStep_ = my_func;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::AllocateUserHistoryOutput(int n)
//  \brief set the number of user-defined history outputs

void Mesh::AllocateUserHistoryOutput(int n) {
  nuser_history_output_ = n;
  user_history_output_names_ = new std::string[n];
  user_history_func_ = new HistoryOutputFunc[n];
  for (int i=0; i<n; i++) user_history_func_[i] = nullptr;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollUserHistoryOutput(int i, HistoryOutputFunc my_func,
//                                         const char *name)
//  \brief Enroll a user-defined history output function and set its name

void Mesh::EnrollUserHistoryOutput(int i, HistoryOutputFunc my_func, const char *name) {
  std::stringstream msg;
  if (i>=nuser_history_output_) {
    msg << "### FATAL ERROR in EnrollUserHistoryOutput function" << std::endl
        << "The number of the user-defined history output (" << i << ") "
        << "exceeds the declared number (" << nuser_history_output_ << ")." << std::endl;
    ATHENA_ERROR(msg);
  }
  user_history_output_names_[i] = name;
  user_history_func_[i] = my_func;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollUserMetric(MetricFunc my_func)
//  \brief Enroll a user-defined metric for arbitrary GR coordinates

void Mesh::EnrollUserMetric(MetricFunc my_func) {
  UserMetric_ = my_func;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollViscosityCoefficient(ViscosityCoeff my_func)
//  \brief Enroll a user-defined magnetic field diffusivity function

void Mesh::EnrollViscosityCoefficient(ViscosityCoeffFunc my_func) {
  ViscosityCoeff_ = my_func;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollConductionCoefficient(ConductionCoeff my_func)
//  \brief Enroll a user-defined thermal conduction function

void Mesh::EnrollConductionCoefficient(ConductionCoeffFunc my_func) {
  ConductionCoeff_ = my_func;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollFieldDiffusivity(FieldDiffusionCoeff my_func)
//  \brief Enroll a user-defined magnetic field diffusivity function

void Mesh::EnrollFieldDiffusivity(FieldDiffusionCoeffFunc my_func) {
  FieldDiffusivity_ = my_func;
  return;
}
//----------------------------------------------------------------------------------------
//! \fn void Mesh::AllocateRealUserMeshDataField(int n)
//  \brief Allocate Real AthenaArrays for user-defned data in Mesh

void Mesh::AllocateRealUserMeshDataField(int n) {
  if (nreal_user_mesh_data_ != 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Mesh::AllocateRealUserMeshDataField"
        << std::endl << "User Mesh data arrays are already allocated" << std::endl;
    ATHENA_ERROR(msg);
  }
  nreal_user_mesh_data_=n;
  ruser_mesh_data = new AthenaArray<Real>[n];
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::AllocateIntUserMeshDataField(int n)
//  \brief Allocate integer AthenaArrays for user-defned data in Mesh

void Mesh::AllocateIntUserMeshDataField(int n) {
  if (nint_user_mesh_data_ != 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Mesh::AllocateIntUserMeshDataField"
        << std::endl << "User Mesh data arrays are already allocated" << std::endl;
    ATHENA_ERROR(msg);
  }
  nint_user_mesh_data_=n;
  iuser_mesh_data = new AthenaArray<int>[n];
  return;
}


//----------------------------------------------------------------------------------------
// \!fn void Mesh::ApplyUserWorkBeforeOutput(ParameterInput *pin)
// \brief Apply MeshBlock::UserWorkBeforeOutput

void Mesh::ApplyUserWorkBeforeOutput(ParameterInput *pin) {
  MeshBlock *pmb = pblock;
  while (pmb != nullptr)  {
    pmb->UserWorkBeforeOutput(pin);
    pmb = pmb->next;
  }
}

//----------------------------------------------------------------------------------------
// \!fn void Mesh::Initialize(int res_flag, ParameterInput *pin)
// \brief  initialization before the main loop

void Mesh::Initialize(int res_flag, ParameterInput *pin) {
  bool iflag = true;
  int inb = nbtotal;
  int nthreads = GetNumMeshThreads();
  int nmb = GetNumMeshBlocksThisRank(Globals::my_rank);
  std::vector<MeshBlock*> pmb_array(nmb);

  do {
    // initialize a vector of MeshBlock pointers
    nmb = GetNumMeshBlocksThisRank(Globals::my_rank);
    if (static_cast<unsigned int>(nmb) != pmb_array.size()) pmb_array.resize(nmb);
    MeshBlock *pmbl = pblock;
    for (int i=0; i<nmb; ++i) {
      pmb_array[i] = pmbl;
      pmbl = pmbl->next;
    }

    if (res_flag==0) {
#pragma omp parallel for num_threads(nthreads)
      for (int i=0; i<nmb; ++i) {
        MeshBlock *pmb = pmb_array[i];
        pmb->ProblemGenerator(pin);
        pmb->pbval->CheckBoundary();
      }
    }

    // add initial perturbation for decaying or impulsive turbulence
    if (((turb_flag == 1) || (turb_flag == 2)) && (res_flag == 0))
      ptrbd->Driving();

    // Create send/recv MPI_Requests for all BoundaryData objects
#pragma omp parallel for num_threads(nthreads)
    for (int i=0; i<nmb; ++i) {
      MeshBlock *pmb = pmb_array[i];
      // BoundaryVariable objects evolved in main TimeIntegratorTaskList:
      pmb->pbval->SetupPersistentMPI();
      // other BoundaryVariable objects:
      if (SELF_GRAVITY_ENABLED == 1)
        pmb->pgrav->pgbval->SetupPersistentMPI();
    }

    // solve gravity for the first time
    if (SELF_GRAVITY_ENABLED == 1)
      pfgrd->Solve(1, 0);
    else if (SELF_GRAVITY_ENABLED == 2)
      pmgrd->Solve(1);

#pragma omp parallel num_threads(nthreads)
    {
      MeshBlock *pmb;
      Hydro *phydro;
      Field *pfield;
      BoundaryValues *pbval;

      // prepare to receive conserved variables
#pragma omp for private(pmb,pbval)
      for (int i=0; i<nmb; ++i) {
        pmb = pmb_array[i]; pbval = pmb->pbval;
        pbval->StartReceiving(BoundaryCommSubset::mesh_init);
      }

      // send conserved variables
#pragma omp for private(pmb,pbval)
      for (int i=0; i<nmb; ++i) {
        pmb = pmb_array[i]; pbval = pmb->pbval;
        pmb->phydro->phbval->SwapHydroQuantity(pmb->phydro->u,
                                               HydroBoundaryQuantity::cons);
        // KGF: (pmb->phydro->u, HydroBoundaryQuantity::cons); where u was bound to &dst
        pmb->phydro->phbval->SendBoundaryBuffers();
        if (MAGNETIC_FIELDS_ENABLED)
          // KGF: (pmb->pfield->b); where b was bound to &dst
          pmb->pfield->pfbval->SendBoundaryBuffers();
      }

      // wait to receive conserved variables
#pragma omp for private(pmb,pbval)
      for (int i=0; i<nmb; ++i) {
        pmb = pmb_array[i]; pbval = pmb->pbval;
        pmb->phydro->phbval->SwapHydroQuantity(pmb->phydro->u,
                                               HydroBoundaryQuantity::cons);
        // KGF: (pmb->phydro->u, HydroBoundaryQuantity::cons); where u was bound to &dst
        pmb->phydro->phbval->ReceiveAndSetBoundariesWithWait();
        if (MAGNETIC_FIELDS_ENABLED)
          // KGF: (pmb->pfield->b); where b was bound to &dst
          pmb->pfield->pfbval->ReceiveAndSetBoundariesWithWait();
        // KGF: disable shearing box bvals/ calls
        // send and receive shearingbox boundary conditions
        // if (SHEARING_BOX)
        //   pmb->phydro->phbval->
        //   SendHydroShearingboxBoundaryBuffersForInit(pmb->phydro->u, true);
        pbval->ClearBoundary(BoundaryCommSubset::mesh_init);
      }

      // With AMR/SMR GR send primitives to enable cons->prim before prolongation
      if (GENERAL_RELATIVITY && multilevel) {
        // prepare to receive primitives
#pragma omp for private(pmb,pbval)
        for (int i=0; i<nmb; ++i) {
          pmb = pmb_array[i]; pbval = pmb->pbval;
          pbval->StartReceiving(BoundaryCommSubset::gr_amr);
        }

        // KGF: the below 2x loops are the only places where "prim" is passed to calls
        // send primitives
#pragma omp for private(pmb,pbval)
        for (int i=0; i<nmb; ++i) {
          pmb = pmb_array[i]; pbval = pmb->pbval;
          pmb->phydro->phbval->SwapHydroQuantity(pmb->phydro->w,
                                                 HydroBoundaryQuantity::prim);
          // KGF: (pmb->phydro->w, HydroBoundaryQuantity::prim); where w was bound to &dst
          pmb->phydro->phbval->SendBoundaryBuffers();
        }

        // wait to receive AMR/SMR GR primitives
#pragma omp for private(pmb,pbval)
        for (int i=0; i<nmb; ++i) {
          pmb = pmb_array[i]; pbval = pmb->pbval;
          pmb->phydro->phbval->SwapHydroQuantity(pmb->phydro->w,
                                                 HydroBoundaryQuantity::prim);
          // KGF: (pmb->phydro->w, HydroBoundaryQuantity::prim); where w was bound to &dst
          pmb->phydro->phbval->ReceiveAndSetBoundariesWithWait();
          pbval->ClearBoundary(BoundaryCommSubset::gr_amr);
          // KGF: just to be sure?
          pmb->phydro->phbval->SwapHydroQuantity(pmb->phydro->u,
                                                 HydroBoundaryQuantity::cons);
        }
      }

      // begin fourth-order correction of midpoint initial condition:
      // --------------------------

      // correct IC on all MeshBlocks or none; switch cannot be toggled independently
      bool correct_ic = pmb_array[0]->precon->correct_ic;
      if (correct_ic == true) {
#pragma omp for private(pmb, phydro, pfield, pbval)
        for (int nb=0; nb<nmb; ++nb) {
          pmb = pmb_array[nb];
          phydro = pmb->phydro;
          pfield = pmb->pfield;
          pbval = pmb->pbval;

          // Assume cell-centered analytic value is computed at all real cells, and ghost
          // cells with the cell-centered U have been exchanged
          int il = pmb->is, iu = pmb->ie, jl = pmb->js, ju = pmb->je,
              kl = pmb->ks, ku = pmb->ke;

          // Laplacian of cell-averaged conserved variables
          AthenaArray<Real> delta_cons_;

          // Allocate memory for 4D Laplacian
          int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
          int ncells2 = 1, ncells3 = 1;
          if (pmb->block_size.nx2 > 1) ncells2 = pmb->block_size.nx2 + 2*(NGHOST);
          if (pmb->block_size.nx3 > 1) ncells3 = pmb->block_size.nx3 + 2*(NGHOST);
          int ncells4 = NHYDRO;
          int nl = 0;
          int nu = ncells4-1;
          delta_cons_.NewAthenaArray(ncells4, ncells3, ncells2, ncells1);

          // Compute and store Laplacian of cell-averaged conserved variables
          pmb->pcoord->Laplacian(phydro->u, delta_cons_, il, iu, jl, ju, kl, ku, nl, nu);
          // TODO(felker): assuming uniform mesh with dx1f=dx2f=dx3f, so this factors out
          // TODO(felker): also, this may need to be dx1v, since Laplacian is cell-center
          Real h = pmb->pcoord->dx1f(il);  // pco->dx1f(i); inside loop
          Real C = (h*h)/24.0;

          // Compute fourth-order approximation to cell-centered conserved variables
          for (int n=nl; n<=nu; ++n) {
            for (int k=kl; k<=ku; ++k) {
              for (int j=jl; j<=ju; ++j) {
                for (int i=il; i<=iu; ++i) {
                  // We do not actually need to store all cell-centered cons. variables,
                  // but the ConservedToPrimitivePointwise() implementation operates on 4D
                  phydro->u(n,k,j,i) = phydro->u(n,k,j,i) + C*delta_cons_(n,k,j,i);
                }
              }
            }
          }
          delta_cons_.DeleteAthenaArray();
        }

        // begin second exchange of ghost cells with corrected cell-averaged <U>
        // -----------------  (mostly copied from above)
        // prepare to receive conserved variables
#pragma omp for private(pmb,pbval)
        for (int i=0; i<nmb; ++i) {
          pmb = pmb_array[i]; pbval = pmb->pbval;
          // no need to re-SetupPersistentMPI() the MPI requests for boundary values
          pbval->StartReceiving(BoundaryCommSubset::mesh_init);
        }

#pragma omp for private(pmb,pbval)
        for (int i=0; i<nmb; ++i) {
          pmb = pmb_array[i]; pbval = pmb->pbval;
          pmb->phydro->phbval->SwapHydroQuantity(pmb->phydro->u,
                                                 HydroBoundaryQuantity::cons);
          // KGF: (pmb->phydro->u, HydroBoundaryQuantity::cons); where u was bound to &dst
          pmb->phydro->phbval->SendBoundaryBuffers();
          if (MAGNETIC_FIELDS_ENABLED)
            // KGF: (pmb->pfield->b); where b was bound to &dst
            pmb->pfield->pfbval->SendBoundaryBuffers();
        }

        // wait to receive conserved variables
#pragma omp for private(pmb,pbval)
        for (int i=0; i<nmb; ++i) {
          pmb = pmb_array[i]; pbval = pmb->pbval;
          pmb->phydro->phbval->SwapHydroQuantity(pmb->phydro->u,
                                                 HydroBoundaryQuantity::cons);
          // KGF: (pmb->phydro->u, HydroBoundaryQuantity::cons); where u was bound to &dst
          pmb->phydro->phbval->ReceiveAndSetBoundariesWithWait();
          if (MAGNETIC_FIELDS_ENABLED)
            // KGF: (pmb->pfield->b); where b was bound to &dst
            pmb->pfield->pfbval->ReceiveAndSetBoundariesWithWait();
          // KGF: disable shearing box bvals/ calls
          // send and receive shearingbox boundary conditions
          // if (SHEARING_BOX)
          //   pmb->phydro->phbval->
          //   SendHydroShearingboxBoundaryBuffersForInit(pmb->phydro->u, true);
          pbval->ClearBoundary(BoundaryCommSubset::mesh_init);
        }
        // -----------------  (verbatim copied from above)
        // end second exchange of ghost cells
      } // end if (correct_ic == true)
      // --------------------------
      // end fourth-order correction of midpoint initial condition

      // Now do prolongation, compute primitives, apply BCs
#pragma omp for private(pmb,pbval,phydro,pfield)
      for (int i=0; i<nmb; ++i) {
        pmb = pmb_array[i]; pbval = pmb->pbval, phydro = pmb->phydro, pfield = pmb->pfield;
        if (multilevel==true)
          pbval->ProlongateBoundaries(time, 0.0);

        int il = pmb->is, iu = pmb->ie, jl = pmb->js, ju = pmb->je, kl = pmb->ks, ku = pmb->ke;
        if (pbval->nblevel[1][1][0] != -1) il-=NGHOST;
        if (pbval->nblevel[1][1][2] != -1) iu+=NGHOST;
        if (pmb->block_size.nx2 > 1) {
          if (pbval->nblevel[1][0][1] != -1) jl-=NGHOST;
          if (pbval->nblevel[1][2][1] != -1) ju+=NGHOST;
        }
        if (pmb->block_size.nx3 > 1) {
          if (pbval->nblevel[0][1][1] != -1) kl-=NGHOST;
          if (pbval->nblevel[2][1][1] != -1) ku+=NGHOST;
        }
        pmb->peos->ConservedToPrimitive(phydro->u, phydro->w1, pfield->b,
                                        phydro->w, pfield->bcc, pmb->pcoord,
                                        il, iu, jl, ju, kl, ku);
        // --------------------------
        int order = pmb->precon->xorder;
        if (order == 4) {
          // fourth-order EOS:
          // for hydro, shrink buffer by 1 on all sides
          if (pbval->nblevel[1][1][0] != -1) il+=1;
          if (pbval->nblevel[1][1][2] != -1) iu-=1;
          if (pbval->nblevel[1][0][1] != -1) jl+=1;
          if (pbval->nblevel[1][2][1] != -1) ju-=1;
          if (pbval->nblevel[0][1][1] != -1) kl+=1;
          if (pbval->nblevel[2][1][1] != -1) ku-=1;
          // for MHD, shrink buffer by 3
          // TODO(felker): add MHD loop limit calculation for 4th order W(U)
          pmb->peos->ConservedToPrimitiveCellAverage(phydro->u, phydro->w1, pfield->b,
                                                     phydro->w, pfield->bcc, pmb->pcoord,
                                                     il, iu, jl, ju, kl, ku);
        }
        // --------------------------
        // end fourth-order EOS
        pmb->phydro->phbval->SwapHydroQuantity(pmb->phydro->w,
                                               HydroBoundaryQuantity::prim);
        pbval->ApplyPhysicalBoundaries(time, 0.0);
      }

      // Calc initial diffusion coefficients
#pragma omp for private(pmb,phydro,pfield)
      for (int i=0; i<nmb; ++i) {
        pmb = pmb_array[i]; phydro = pmb->phydro, pfield = pmb->pfield;
        if (phydro->phdif->hydro_diffusion_defined)
          phydro->phdif->SetHydroDiffusivity(phydro->w, pfield->bcc);
        if (MAGNETIC_FIELDS_ENABLED) {
          if (pfield->pfdif->field_diffusion_defined)
            pfield->pfdif->SetFieldDiffusivity(phydro->w, pfield->bcc);
        }
      }

      if ((res_flag == 0) && (adaptive == true)) {
#pragma omp for
        for (int i=0; i<nmb; ++i) {
          pmb_array[i]->pmr->CheckRefinementCondition();
        }
      }
    } // omp parallel

    if ((res_flag == 0) && (adaptive == true)) {
      iflag = false;
      int onb = nbtotal;
      AdaptiveMeshRefinement(pin);
      if (nbtotal == onb) {
        iflag=true;
      } else if (nbtotal < onb && Globals::my_rank == 0) {
        std::cout << "### Warning in Mesh::Initialize" << std::endl
                  << "The number of MeshBlocks decreased during AMR grid initialization."
                  << std::endl
                  << "Possibly the refinement criteria have a problem." << std::endl;
      }
      if (nbtotal > 2*inb && Globals::my_rank == 0) {
        std::cout
            << "### Warning in Mesh::Initialize" << std::endl
            << "The number of MeshBlocks increased more than twice during initialization."
            << std::endl
            << "More computing power than you expected may be required." << std::endl;
      }
    }
  } while (iflag == false);

  // calculate the first time step
#pragma omp parallel for num_threads(nthreads)
  for (int i=0; i<nmb; ++i) {
    pmb_array[i]->phydro->NewBlockTimeStep();
  }

  NewTimeStep();
  return;
}

//----------------------------------------------------------------------------------------
//! \fn MeshBlock* Mesh::FindMeshBlock(int tgid)
//  \brief return the MeshBlock whose gid is tgid

MeshBlock* Mesh::FindMeshBlock(int tgid) {
  MeshBlock *pbl = pblock;
  while (pbl != nullptr) {
    if (pbl->gid == tgid)
      break;
    pbl = pbl->next;
  }
  return pbl;
}

//----------------------------------------------------------------------------------------
// \!fn void Mesh::LoadBalance(Real *clist, int *rlist, int *slist, int *nlist, int nb)
// \brief Calculate distribution of MeshBlocks based on the cost list

void Mesh::LoadBalance(Real *clist, int *rlist, int *slist, int *nlist, int nb) {
  std::stringstream msg;
  Real real_max = std::numeric_limits<Real>::max();
  Real totalcost=0, maxcost=0.0, mincost=(real_max);

  for (int i=0; i<nb; i++) {
    totalcost+=clist[i];
    mincost=std::min(mincost,clist[i]);
    maxcost=std::max(maxcost,clist[i]);
  }
  int j=(Globals::nranks)-1;
  Real targetcost=totalcost/Globals::nranks;
  Real mycost=0.0;
  // create rank list from the end: the master node should have less load
  for (int i=nb-1; i>=0; i--) {
    if (targetcost == 0.0) {
      msg << "### FATAL ERROR in LoadBalance" << std::endl
          << "There is at least one process which has no MeshBlock" << std::endl
          << "Decrease the number of processes or use smaller MeshBlocks." << std::endl;
      ATHENA_ERROR(msg);
    }
    mycost+=clist[i];
    rlist[i]=j;
    if (mycost >= targetcost && j>0) {
      j--;
      totalcost-=mycost;
      mycost=0.0;
      targetcost=totalcost/(j+1);
    }
  }
  slist[0]=0;
  j=0;
  for (int i=1; i<nb; i++) { // make the list of nbstart and nblocks
    if (rlist[i] != rlist[i-1]) {
      nlist[j]=i-nslist[j];
      slist[++j]=i;
    }
  }
  nlist[j]=nb-slist[j];

#ifdef MPI_PARALLEL
  if (nb % (Globals::nranks * num_mesh_threads_) != 0 && adaptive == false
      && maxcost == mincost && Globals::my_rank == 0) {
    std::cout << "### Warning in LoadBalance" << std::endl
              << "The number of MeshBlocks cannot be divided evenly. "
              << "This will result in poor load balancing." << std::endl;
  }
#endif
  if ((Globals::nranks)*(num_mesh_threads_) > nb) {
    msg << "### FATAL ERROR in LoadBalance" << std::endl
        << "There are fewer MeshBlocks than OpenMP threads on each MPI rank" << std::endl
        << "Decrease the number of threads or use more MeshBlocks." << std::endl;
    ATHENA_ERROR(msg);
  }

  return;
}

//----------------------------------------------------------------------------------------
// \!fn void Mesh::SetBlockSizeAndBoundaries(LogicalLocation loc,
//                 RegionSize &block_size, BundaryFlag *block_bcs)
// \brief Set the physical part of a block_size structure and block boundary conditions

void Mesh::SetBlockSizeAndBoundaries(LogicalLocation loc, RegionSize &block_size,
                                     BoundaryFlag *block_bcs) {
  std::int64_t &lx1 = loc.lx1;
  std::int64_t &lx2 = loc.lx2;
  std::int64_t &lx3 = loc.lx3;
  int &ll = loc.level;
  std::int64_t nrbx_ll = nrbx1<<(ll-root_level);

  // calculate physical block size, x1
  if (lx1 == 0) {
    block_size.x1min=mesh_size.x1min;
    block_bcs[BoundaryFace::inner_x1]=mesh_bcs[BoundaryFace::inner_x1];
  } else {
    Real rx = ComputeMeshGeneratorX(lx1, nrbx_ll, use_uniform_meshgen_fn_[X1DIR]);
    block_size.x1min=MeshGenerator_[X1DIR](rx,mesh_size);
    block_bcs[BoundaryFace::inner_x1]=BoundaryFlag::block;
  }
  if (lx1 == nrbx_ll-1) {
    block_size.x1max=mesh_size.x1max;
    block_bcs[BoundaryFace::outer_x1]=mesh_bcs[BoundaryFace::outer_x1];
  } else {
    Real rx = ComputeMeshGeneratorX(lx1+1, nrbx_ll, use_uniform_meshgen_fn_[X1DIR]);
    block_size.x1max=MeshGenerator_[X1DIR](rx,mesh_size);
    block_bcs[BoundaryFace::outer_x1]=BoundaryFlag::block;
  }

  // calculate physical block size, x2
  if (mesh_size.nx2 == 1) {
    block_size.x2min=mesh_size.x2min;
    block_size.x2max=mesh_size.x2max;
    block_bcs[BoundaryFace::inner_x2]=mesh_bcs[BoundaryFace::inner_x2];
    block_bcs[BoundaryFace::outer_x2]=mesh_bcs[BoundaryFace::outer_x2];
  } else {
    nrbx_ll = nrbx2<<(ll-root_level);
    if (lx2 == 0) {
      block_size.x2min=mesh_size.x2min;
      block_bcs[BoundaryFace::inner_x2]=mesh_bcs[BoundaryFace::inner_x2];
    } else {
      Real rx = ComputeMeshGeneratorX(lx2, nrbx_ll, use_uniform_meshgen_fn_[X2DIR]);
      block_size.x2min=MeshGenerator_[X2DIR](rx,mesh_size);
      block_bcs[BoundaryFace::inner_x2]=BoundaryFlag::block;
    }
    if (lx2 == (nrbx_ll)-1) {
      block_size.x2max=mesh_size.x2max;
      block_bcs[BoundaryFace::outer_x2]=mesh_bcs[BoundaryFace::outer_x2];
    } else {
      Real rx = ComputeMeshGeneratorX(lx2+1, nrbx_ll, use_uniform_meshgen_fn_[X2DIR]);
      block_size.x2max=MeshGenerator_[X2DIR](rx,mesh_size);
      block_bcs[BoundaryFace::outer_x2]=BoundaryFlag::block;
    }
  }

  // calculate physical block size, x3
  if (mesh_size.nx3 == 1) {
    block_size.x3min=mesh_size.x3min;
    block_size.x3max=mesh_size.x3max;
    block_bcs[BoundaryFace::inner_x3]=mesh_bcs[BoundaryFace::inner_x3];
    block_bcs[BoundaryFace::outer_x3]=mesh_bcs[BoundaryFace::outer_x3];
  } else {
    nrbx_ll = nrbx3<<(ll-root_level);
    if (lx3 == 0) {
      block_size.x3min=mesh_size.x3min;
      block_bcs[BoundaryFace::inner_x3]=mesh_bcs[BoundaryFace::inner_x3];
    } else {
      Real rx = ComputeMeshGeneratorX(lx3, nrbx_ll, use_uniform_meshgen_fn_[X3DIR]);
      block_size.x3min=MeshGenerator_[X3DIR](rx,mesh_size);
      block_bcs[BoundaryFace::inner_x3]=BoundaryFlag::block;
    }
    if (lx3 == (nrbx_ll)-1) {
      block_size.x3max=mesh_size.x3max;
      block_bcs[BoundaryFace::outer_x3]=mesh_bcs[BoundaryFace::outer_x3];
    } else {
      Real rx = ComputeMeshGeneratorX(lx3+1, nrbx_ll, use_uniform_meshgen_fn_[X3DIR]);
      block_size.x3max = MeshGenerator_[X3DIR](rx,mesh_size);
      block_bcs[BoundaryFace::outer_x3]=BoundaryFlag::block;
    }
  }

  block_size.x1rat = mesh_size.x1rat;
  block_size.x2rat = mesh_size.x2rat;
  block_size.x3rat = mesh_size.x3rat;

  return;
}

//----------------------------------------------------------------------------------------
// \!fn void Mesh::AdaptiveMeshRefinement(ParameterInput *pin)
// \brief Main function for adaptive mesh refinement

void Mesh::AdaptiveMeshRefinement(ParameterInput *pin) {
  MeshBlock *pmb;
  // compute nleaf= number of leaf MeshBlocks per refined block
  int nleaf = 2, dim = 1;
  if (mesh_size.nx2 > 1) nleaf = 4, dim = 2;
  if (mesh_size.nx3 > 1) nleaf = 8, dim = 3;

  // collect refinement flags from all the meshblocks
  // count the number of the blocks to be (de)refined
  nref[Globals::my_rank] = 0;
  nderef[Globals::my_rank] = 0;
  pmb = pblock;
  while (pmb != nullptr) {
    if (pmb->pmr->refine_flag_ ==  1) nref[Globals::my_rank]++;
    if (pmb->pmr->refine_flag_ == -1) nderef[Globals::my_rank]++;
    pmb = pmb->next;
  }
#ifdef MPI_PARALLEL
  MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, nref,   1, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, nderef, 1, MPI_INT, MPI_COMM_WORLD);
#endif

  // count the number of the blocks to be (de)refined and displacement
  int tnref = 0, tnderef = 0;
  for (int n=0; n<Globals::nranks; n++) {
    tnref  += nref[n];
    tnderef += nderef[n];
  }
  if (tnref == 0 && tnderef < nleaf) // nothing to do
    return;

  int rd = 0, dd = 0;
  for (int n=0; n<Globals::nranks; n++) {
    rdisp[n] = rd;
    ddisp[n] = dd;
    // technically could overflow, since sizeof() operator returns
    // std::size_t = long unsigned int > int
    // on many platforms (LP64). However, these are used below in MPI calls for
    // integer arguments (recvcounts, displs). MPI does not support > 64-bit count ranges
    bnref[n] = static_cast<int>(nref[n]*sizeof(LogicalLocation));
    bnderef[n] = static_cast<int>(nderef[n]*sizeof(LogicalLocation));
    brdisp[n] = static_cast<int>(rd*sizeof(LogicalLocation));
    bddisp[n] = static_cast<int>(dd*sizeof(LogicalLocation));
    rd += nref[n];
    dd += nderef[n];
  }

  // allocate memory for the location arrays
  LogicalLocation *lref{}, *lderef{}, *clderef{};
  if (tnref > 0)
    lref = new LogicalLocation[tnref];
  if (tnderef >= nleaf) {
    lderef = new LogicalLocation[tnderef];
    clderef = new LogicalLocation[tnderef/nleaf];
  }

  // collect the locations and costs
  int iref = rdisp[Globals::my_rank], ideref = ddisp[Globals::my_rank];
  pmb = pblock;
  while (pmb != nullptr) {
    if (pmb->pmr->refine_flag_ ==  1)
      lref[iref++] = pmb->loc;
    if (pmb->pmr->refine_flag_ == -1 && tnderef >= nleaf)
      lderef[ideref++] = pmb->loc;
    pmb = pmb->next;
  }
#ifdef MPI_PARALLEL
  if (tnref > 0) {
    MPI_Allgatherv(MPI_IN_PLACE, bnref[Globals::my_rank],   MPI_BYTE,
                   lref,   bnref,   brdisp, MPI_BYTE, MPI_COMM_WORLD);
  }
  if (tnderef >= nleaf) {
    MPI_Allgatherv(MPI_IN_PLACE, bnderef[Globals::my_rank], MPI_BYTE,
                   lderef, bnderef, bddisp, MPI_BYTE, MPI_COMM_WORLD);
  }
#endif

  // calculate the list of the newly derefined blocks
  int ctnd = 0;
  if (tnderef >= nleaf) {
    int lk = 0, lj = 0;
    if (mesh_size.nx2 > 1) lj = 1;
    if (mesh_size.nx3 > 1) lk = 1;
    for (int n=0; n<tnderef; n++) {
      if ((lderef[n].lx1 & 1LL) == 0LL &&
          (lderef[n].lx2 & 1LL) == 0LL &&
          (lderef[n].lx3 & 1LL) == 0LL) {
        int r = n, rr = 0;
        for (std::int64_t k=0; k<=lk; k++) {
          for (std::int64_t j=0; j<=lj; j++) {
            for (std::int64_t i=0; i<=1; i++) {
              if (r < tnderef) {
                if ((lderef[n].lx1+i) == lderef[r].lx1
                    && (lderef[n].lx2+j) == lderef[r].lx2
                    && (lderef[n].lx3+k) == lderef[r].lx3
                    &&  lderef[n].level  == lderef[r].level)
                  rr++;
                r++;
              }
            }
          }
        }
        if (rr == nleaf) {
          clderef[ctnd].lx1   = (lderef[n].lx1>>1);
          clderef[ctnd].lx2   = (lderef[n].lx2>>1);
          clderef[ctnd].lx3   = (lderef[n].lx3>>1);
          clderef[ctnd].level = lderef[n].level-1;
          ctnd++;
        }
      }
    }
  }
  // sort the lists by level
  if (ctnd > 1)
    std::sort(clderef, &(clderef[ctnd-1]), LogicalLocation::Greater);

  if (tnderef >= nleaf)
    delete [] lderef;

  // Now the lists of the blocks to be refined and derefined are completed
  // Start tree manipulation
  // Step 1. perform refinement
  int nnew = 0, ndel = 0, ntot = 0;
  for (int n=0; n<tnref; n++) {
    MeshBlockTree *bt=tree.FindMeshBlock(lref[n]);
    bt->Refine(tree, dim, mesh_bcs, nrbx1, nrbx2, nrbx3, root_level, nnew);
  }
  if (tnref != 0)
    delete [] lref;

  // Step 2. perform derefinement
  for (int n=0; n<ctnd; n++) {
    MeshBlockTree *bt = tree.FindMeshBlock(clderef[n]);
    bt->Derefine(tree, dim, mesh_bcs, nrbx1, nrbx2, nrbx3, root_level, ndel);
  }
  if (tnderef >= nleaf)
    delete [] clderef;
  ntot = nbtotal + nnew - ndel;
  if (nnew == 0 && ndel == 0)
    return; // nothing to do
  // Tree manipulation completed
  nbnew += nnew; nbdel += ndel;

  // Block exchange
  // Step 1. construct new lists
  LogicalLocation *newloc = new LogicalLocation[ntot];
  int *newrank = new int[ntot];
  Real *newcost = new Real[ntot];
  int *newtoold = new int[ntot];
  int *oldtonew = new int[nbtotal];
  int nbtold = nbtotal;
  tree.GetMeshBlockList(newloc,newtoold,nbtotal);

  // create a list mapping the previous gid to the current one
  oldtonew[0] = 0;
  int mb_idx = 1;
  for (int n=1; n<ntot; n++) {
    if (newtoold[n] == newtoold[n-1] + 1) { // normal
      oldtonew[mb_idx++] = n;
    } else if (newtoold[n] == newtoold[n-1] + nleaf) { // derefined
      for (int j=0; j<nleaf-1; j++)
        oldtonew[mb_idx++] = n-1;
      oldtonew[mb_idx++] = n;
    }
  }
  // fill the last block
  for ( ; mb_idx<nbtold; mb_idx++)
    oldtonew[mb_idx] = ntot-1;

#ifdef MPI_PARALLEL
  // share the cost list
  MPI_Allgatherv(MPI_IN_PLACE, nblist[Globals::my_rank], MPI_ATHENA_REAL,
                 costlist, nblist, nslist, MPI_ATHENA_REAL, MPI_COMM_WORLD);
#endif

  current_level = 0;
  for (int n=0; n<ntot; n++) {
    // "on" = "old n" = "old gid" = "old global MeshBlock ID"
    int on = newtoold[n];
    if (newloc[n].level>current_level) // set the current max level
      current_level = newloc[n].level;
    if (newloc[n].level >= loclist[on].level) { // same or refined
      newcost[n] = costlist[on];
    } else {
      Real acost = 0.0;
      for (int l=0; l<nleaf; l++)
        acost += costlist[on+l];
      newcost[n] = acost/nleaf;
    }
  }
#ifdef MPI_PARALLEL
  // store old nbstart and nbend before load balancing in Step 2.
  int onbs = nslist[Globals::my_rank];
  int onbe = onbs + nblist[Globals::my_rank] - 1;
#endif
  // Step 2. Calculate new load balance
  LoadBalance(newcost, newrank, nslist, nblist, ntot);

  int nbs = nslist[Globals::my_rank];
  int nbe = nbs + nblist[Globals::my_rank] - 1;

  int f2, f3;
  int bnx1 = pblock->block_size.nx1;
  int bnx2 = pblock->block_size.nx2;
  int bnx3 = pblock->block_size.nx3;
  if (mesh_size.nx2 > 1) {
    f2 = 1;
  } else {
    f2 = 0;
  }
  if (mesh_size.nx3 > 1) {
    f3 = 1;
  } else {
    f3 = 0;
  }

#ifdef MPI_PARALLEL
  // Step 3. count the number of the blocks to be sent / received
  int nsend = 0, nrecv = 0;
  for (int n=nbs; n<=nbe; n++) {
    int on = newtoold[n];
    if (loclist[on].level > newloc[n].level) { // f2c
      for (int k=0; k<nleaf; k++) {
        if (ranklist[on+k] != Globals::my_rank)
          nrecv++;
      }
    } else {
      if (ranklist[on] != Globals::my_rank)
        nrecv++;
    }
  }
  for (int n=onbs; n<=onbe; n++) {
    int nn = oldtonew[n];
    if (loclist[n].level < newloc[nn].level) { // c2f
      for (int k=0; k<nleaf; k++) {
        if (newrank[nn+k] != Globals::my_rank)
          nsend++;
      }
    } else {
      if (newrank[nn] != Globals::my_rank)
        nsend++;
    }
  }

  // Step 4. calculate buffer sizes
  Real **sendbuf, **recvbuf;
  // use the first MeshBlock in the linked list of blocks belonging to this MPI rank as a
  // representative of all MeshBlocks for counting the "SMR/AMR-enrolled" quantities

  // int num_cc = pblock->pmr->pvars_cc_.size();
  int num_fc = pblock->pmr->pvars_fc_.size();

  // KGF: need to calculate weighted sum of nx4 overall num_cc (originally just NHYDRO)
  int nx4_tot = 0; // NHYDRO;
  for (auto cc_pair : pblock->pmr->pvars_cc_) {
    AthenaArray<Real> *var_cc = std::get<0>(cc_pair);
    nx4_tot += var_cc->GetDim4();
  }

  // cell-centered quantities enrolled in SMR/AMR
  int bssame = bnx1*bnx2*bnx3*nx4_tot;
  int bsf2c = (bnx1/2)*((bnx2 + 1)/2)*((bnx3 + 1)/2)*nx4_tot;
  int bsc2f = (bnx1/2 + 2)*((bnx2 + 1)/2 + 2*f2)*((bnx3 + 1)/2 + 2*f3)*nx4_tot;
  // face-centered quantities enrolled in SMR/AMR
  bssame += num_fc*((bnx1 + 1)*bnx2*bnx3 + bnx1*(bnx2 + f2)*bnx3 + bnx1*bnx2*(bnx3 + f3));
  bsf2c += num_fc*(((bnx1/2) + 1)*((bnx2 + 1)/2)*((bnx3 + 1)/2)
                   + (bnx1/2)*(((bnx2 + 1)/2) + f2)*((bnx3 + 1)/2)
                   + (bnx1/2)*((bnx2 + 1)/2)*(((bnx3 + 1)/2) + f3));
  bsc2f += num_fc*(((bnx1/2) + 1 + 2)*((bnx2 + 1)/2 + 2*f2)*((bnx3 + 1)/2 + 2*f3)
                   + (bnx1/2 + 2)*(((bnx2 + 1)/2) + f2 + 2*f2)*((bnx3 + 1)/2 + 2*f3)
                   + (bnx1/2 + 2)*((bnx2 + 1)/2 + 2*f2)*(((bnx3 + 1)/2) + f3 +2*f3));
  // add one more element to buffer size for storing the derefinement counter
  bssame++;

  MPI_Request *req_send, *req_recv;
  // Step 5. allocate and start receiving buffers
  if (nrecv != 0) {
    recvbuf = new Real*[nrecv];
    req_recv = new MPI_Request[nrecv];
    int rb_idx = 0;     // recv buffer index
    for (int n=nbs; n<=nbe; n++) {
      int on = newtoold[n];
      LogicalLocation &oloc = loclist[on];
      LogicalLocation &nloc = newloc[n];
      if (oloc.level > nloc.level) { // f2c
        for (int l=0; l<nleaf; l++) {
          if (ranklist[on+l] == Globals::my_rank) continue;
          LogicalLocation &lloc = loclist[on+l];
          int ox1 = lloc.lx1 & 1LL, ox2 = lloc.lx2 & 1LL, ox3 = lloc.lx3 & 1LL;
          recvbuf[rb_idx] = new Real[bsf2c];
          int tag = CreateAMRMPITag(n-nbs, ox1, ox2, ox3);
          MPI_Irecv(recvbuf[rb_idx], bsf2c, MPI_ATHENA_REAL, ranklist[on+l],
                    tag, MPI_COMM_WORLD, &(req_recv[rb_idx]));
          rb_idx++;
        }
      } else { // same or c2f
        if (ranklist[on] == Globals::my_rank) continue;
        int size;
        if (oloc.level == nloc.level) {
          size = bssame;
        } else {
          size = bsc2f;
        }
        recvbuf[rb_idx] = new Real[size];
        int tag = CreateAMRMPITag(n-nbs, 0, 0, 0);
        MPI_Irecv(recvbuf[rb_idx], size, MPI_ATHENA_REAL, ranklist[on],
                  tag, MPI_COMM_WORLD, &(req_recv[rb_idx]));
        rb_idx++;
      }
    }
  }
  // Step 6. allocate, pack and start sending buffers
  if (nsend != 0) {
    sendbuf = new Real*[nsend];
    req_send = new MPI_Request[nsend];
    int sb_idx = 0;      // send buffer index
    for (int n=onbs; n<=onbe; n++) {
      int nn = oldtonew[n];
      LogicalLocation &oloc = loclist[n];
      LogicalLocation &nloc = newloc[nn];
      MeshBlock* pb = FindMeshBlock(n);
      if (nloc.level == oloc.level) { // same level
        if (newrank[nn] == Globals::my_rank) continue;
        sendbuf[sb_idx] = new Real[bssame];

        // pack
        int p = 0;
        // KGF: CellCentered step 6, branch 1 (same2same: just pack+send)

        // for (auto cc_it = pb->pmr->pvars_cc_.begin();
        //      cc_it != pb->pmr->pvars_cc_.end(); ++cc_it) {
        //  AthenaArray<Real> *var_cc = std::get<0>(*cc_it);

        // KGF: C++11 range-based for loop
        for (auto cc_pair : pb->pmr->pvars_cc_) {
        // for (std::tuple<AthenaArray<Real>*, AthenaArray<Real>*> cc_pair
          //   : pb->pmr->pvars_cc_) {

          // KGF: No need to dereference an Iterator object in range-based for loop:
          AthenaArray<Real> *var_cc = std::get<0>(cc_pair);
          int nu = var_cc->GetDim4() - 1;
          BufferUtility::PackData(*var_cc, sendbuf[sb_idx], 0, nu,
                                  pb->is, pb->ie, pb->js, pb->je, pb->ks, pb->ke, p);
        }

        // KGF: FaceCentered step 6, branch 1 (same2same: just pack+send)
        // for (auto fc_it = pb->pmr->pvars_fc_.begin();
        //      fc_it != pb->pmr->pvars_fc_.end(); ++fc_it) {

        // KGF: C++11 range-based for loop
        for (auto fc_pair : pb->pmr->pvars_fc_) {
          FaceField *var_fc = std::get<0>(fc_pair);
          BufferUtility::PackData((*var_fc).x1f, sendbuf[sb_idx],
                                  pb->is, pb->ie+1, pb->js, pb->je, pb->ks, pb->ke, p);
          BufferUtility::PackData((*var_fc).x2f, sendbuf[sb_idx],
                                  pb->is, pb->ie, pb->js, pb->je+f2, pb->ks, pb->ke, p);
          BufferUtility::PackData((*var_fc).x3f, sendbuf[sb_idx],
                                  pb->is, pb->ie, pb->js, pb->je, pb->ks, pb->ke+f3, p);
        }

        // KGF: dangerous? casting from "Real *" to "int *"
        int *dcp = reinterpret_cast<int *>(&(sendbuf[sb_idx][p]));
        *dcp = pb->pmr->deref_count_;
        int tag = CreateAMRMPITag(nn-nslist[newrank[nn]], 0, 0, 0);
        MPI_Isend(sendbuf[sb_idx], bssame, MPI_ATHENA_REAL, newrank[nn],
                  tag, MPI_COMM_WORLD, &(req_send[sb_idx]));
        sb_idx++;
      } else if (nloc.level > oloc.level) { // c2f
        for (int l=0; l<nleaf; l++) {
          if (newrank[nn+l] == Globals::my_rank) continue;
          LogicalLocation &lloc = newloc[nn+l];
          int ox1 = lloc.lx1 & 1LL, ox2 = lloc.lx2 & 1LL, ox3 = lloc.lx3 & 1LL;
          sendbuf[sb_idx] = new Real[bsc2f];
          // pack
          int il, iu, jl, ju, kl, ku;
          if (ox1 == 0) il = pb->is - 1,               iu = pb->is + pb->block_size.nx1/2;
          else        il = pb->is + pb->block_size.nx1/2-1,  iu = pb->ie + 1;
          if (ox2 == 0) jl = pb->js - f2,              ju = pb->js + pb->block_size.nx2/2;
          else        jl = pb->js + pb->block_size.nx2/2-f2, ju = pb->je + f2;
          if (ox3 == 0) kl = pb->ks - f3,              ku = pb->ks + pb->block_size.nx3/2;
          else        kl = pb->ks + pb->block_size.nx3/2-f3, ku = pb->ke + f3;
          int p = 0;

          // KGF: CellCentered step 6, branch 2 (c2f: just pack+send)
          for (auto cc_pair : pb->pmr->pvars_cc_) {
            AthenaArray<Real> *var_cc = std::get<0>(cc_pair);
            int nu = var_cc->GetDim4() - 1;
            BufferUtility::PackData((*var_cc), sendbuf[sb_idx], 0, nu,
                                    il, iu, jl, ju, kl, ku, p);
          }

          // KGF: FaceCentered step 6, branch 2 (c2f: just pack+send)
          for (auto fc_pair : pb->pmr->pvars_fc_) {
            FaceField *var_fc = std::get<0>(fc_pair);
            BufferUtility::PackData((*var_fc).x1f, sendbuf[sb_idx],
                                    il, iu+1, jl, ju, kl, ku, p);
            BufferUtility::PackData((*var_fc).x2f, sendbuf[sb_idx],
                                    il, iu, jl, ju+f2, kl, ku, p);
            BufferUtility::PackData((*var_fc).x3f, sendbuf[sb_idx],
                                    il, iu, jl, ju, kl, ku+f3, p);
          }
          int tag = CreateAMRMPITag(nn+l-nslist[newrank[nn+l]], 0, 0, 0);
          MPI_Isend(sendbuf[sb_idx], bsc2f, MPI_ATHENA_REAL, newrank[nn+l],
                    tag, MPI_COMM_WORLD, &(req_send[sb_idx]));
          sb_idx++;
        }
      } else { // f2c
        if (newrank[nn] == Globals::my_rank) continue;
        int ox1 = oloc.lx1 & 1LL, ox2 = oloc.lx2 & 1LL, ox3 = oloc.lx3 & 1LL;
        sendbuf[sb_idx] = new Real[bsf2c];
        // restrict and pack
        MeshRefinement *pmr = pb->pmr;

        int p = 0;
        // KGF: CellCentered step 6, branch 3 (f2c: restrict, pack, send)
        for (auto cc_pair : pb->pmr->pvars_cc_) {
          AthenaArray<Real> *var_cc = std::get<0>(cc_pair);
          AthenaArray<Real> *coarse_cc = std::get<1>(cc_pair);
          int nu = var_cc->GetDim4() - 1;
          pmr->RestrictCellCenteredValues((*var_cc), (*coarse_cc),
                                          0, nu,
                                          pb->cis, pb->cie,
                                          pb->cjs, pb->cje,
                                          pb->cks, pb->cke);
          BufferUtility::PackData((*coarse_cc), sendbuf[sb_idx], 0, nu,
                                  pb->cis, pb->cie,
                                  pb->cjs, pb->cje,
                                  pb->cks, pb->cke, p);
        }
        // KGF: FaceCentered step 6, branch 3 (f2c: restrict, pack, send)
        for (auto fc_pair : pb->pmr->pvars_fc_) {
          FaceField *var_fc = std::get<0>(fc_pair);
          FaceField *coarse_fc = std::get<1>(fc_pair);
          pmr->RestrictFieldX1((*var_fc).x1f, (*coarse_fc).x1f,
                               pb->cis, pb->cie+1,
                               pb->cjs, pb->cje,
                               pb->cks, pb->cke);
          BufferUtility::PackData((*coarse_fc).x1f, sendbuf[sb_idx],
                                  pb->cis, pb->cie+1,
                                  pb->cjs, pb->cje,
                                  pb->cks, pb->cke, p);
          pmr->RestrictFieldX2((*var_fc).x2f, (*coarse_fc).x2f,
                               pb->cis, pb->cie,
                               pb->cjs, pb->cje+f2,
                               pb->cks, pb->cke);
          BufferUtility::PackData((*coarse_fc).x2f, sendbuf[sb_idx],
                                  pb->cis, pb->cie,
                                  pb->cjs, pb->cje+f2,
                                  pb->cks, pb->cke, p);
          pmr->RestrictFieldX3((*var_fc).x3f, (*coarse_fc).x3f,
                               pb->cis, pb->cie,
                               pb->cjs, pb->cje,
                               pb->cks, pb->cke+f3);
          BufferUtility::PackData((*coarse_fc).x3f, sendbuf[sb_idx],
                                  pb->cis, pb->cie,
                                  pb->cjs, pb->cje,
                                  pb->cks, pb->cke+f3, p);
        }

        int tag = CreateAMRMPITag(nn-nslist[newrank[nn]], ox1, ox2, ox3);
        MPI_Isend(sendbuf[sb_idx], bsf2c, MPI_ATHENA_REAL, newrank[nn],
                  tag, MPI_COMM_WORLD, &(req_send[sb_idx]));
        sb_idx++;
      }
    }
  }
#endif // MPI_PARALLEL

  // Step 7. construct a new MeshBlock list
  // move the data within the node
  MeshBlock *newlist = nullptr;
  RegionSize block_size = pblock->block_size;

  for (int n=nbs; n<=nbe; n++) {
    int on = newtoold[n];
    if ((ranklist[on] == Globals::my_rank) && (loclist[on].level == newloc[n].level)) {
      // on the same node and same level -> just move it
      MeshBlock* pob = FindMeshBlock(on);
      if (pob->prev == nullptr) {
        pblock = pob->next;
      } else {
        pob->prev->next = pob->next;
      }
      if (pob->next != nullptr) pob->next->prev = pob->prev;
      pob->next = nullptr;
      if (n == nbs) { // first
        pob->prev = nullptr;
        newlist = pob;
        pmb = newlist;
      } else {
        pmb->next = pob;
        pob->prev = pmb;
        pmb = pmb->next;
      }
      pmb->gid = n;
      pmb->lid = n - nbs;
    } else {
      // on a different level or node - create a new block
      BoundaryFlag block_bcs[6];
      block_size.nx1 = bnx1, block_size.nx2 = bnx2, block_size.nx3 = bnx3;
      SetBlockSizeAndBoundaries(newloc[n], block_size, block_bcs);
      if (n == nbs) { // first
        newlist = new MeshBlock(n, n-nbs, newloc[n], block_size, block_bcs, this,
                                pin, gflag, true);
        pmb = newlist;
      } else {
        pmb->next = new MeshBlock(n, n-nbs, newloc[n], block_size, block_bcs, this,
                                  pin, gflag, true);
        pmb->next->prev = pmb;
        pmb = pmb->next;
      }
      // fill the conservative variables
      if ((loclist[on].level > newloc[n].level)) { // fine to coarse
        for (int ll=0; ll<nleaf; ll++) {
          if (ranklist[on+ll] != Globals::my_rank) continue;
          // on the same node - restriction
          MeshBlock* pob = FindMeshBlock(on+ll);
          MeshRefinement *pmr = pob->pmr;
          int il = pmb->is + ((loclist[on+ll].lx1 & 1LL) == 1LL)*pmb->block_size.nx1/2;
          int jl = pmb->js + ((loclist[on+ll].lx2 & 1LL) == 1LL)*pmb->block_size.nx2/2;
          int kl = pmb->ks + ((loclist[on+ll].lx3 & 1LL) == 1LL)*pmb->block_size.nx3/2;

          // KGF: CellCentered step 7: f2c, same node (just restrict+copy, no pack/send)

          // KGF: absent a zip() feature for range-based for loops, manually advance the
          // iterator over "SMR/AMR-enrolled" cell-centered quantities on the new
          // MeshBlock in lock-step with pob
          auto pmb_cc_it = pmb->pmr->pvars_cc_.begin();
          // iterate MeshRefinement std::vectors on pob
          for (auto cc_pair : pmr->pvars_cc_) {
            AthenaArray<Real> *var_cc = std::get<0>(cc_pair);
            AthenaArray<Real> *coarse_cc = std::get<1>(cc_pair);
            int nu = var_cc->GetDim4() - 1;
            pmr->RestrictCellCenteredValues((*var_cc), (*coarse_cc),
                                            0, nu,
                                            pob->cis, pob->cie,
                                            pob->cjs, pob->cje,
                                            pob->cks, pob->cke);
            // KGF:
            // copy from old/original/other MeshBlock (pob) to newly created block (pmb)
            AthenaArray<Real> &src = (*coarse_cc);
            AthenaArray<Real> &dst = *(std::get<0>(*pmb_cc_it)); // pmb->phydro->u;
            for (int nv=0; nv<=nu; nv++) {
              for (int k=kl, fk=pob->cks; fk<=pob->cke; k++, fk++) {
                for (int j=jl, fj=pob->cjs; fj<=pob->cje; j++, fj++) {
                  for (int i=il, fi=pob->cis; fi<=pob->cie; i++, fi++)
                    dst(nv, k, j, i) = src(nv, fk, fj, fi);
                }
              }
            }
            pmb_cc_it++;
          } // end loop over var_cc

          // KGF: FaceCentered step 7: f2c, same node (just restrict+copy, no pack/send)
          auto pmb_fc_it = pmb->pmr->pvars_fc_.begin();
          for (auto fc_pair : pmr->pvars_fc_) {
            FaceField *var_fc = std::get<0>(fc_pair);
            FaceField *coarse_fc = std::get<1>(fc_pair);
            pmr->RestrictFieldX1((*var_fc).x1f, (*coarse_fc).x1f,
                                 pob->cis, pob->cie+1,
                                 pob->cjs, pob->cje,
                                 pob->cks, pob->cke);
            pmr->RestrictFieldX2((*var_fc).x2f, (*coarse_fc).x2f,
                                 pob->cis, pob->cie,
                                 pob->cjs, pob->cje+f2,
                                 pob->cks, pob->cke);
            pmr->RestrictFieldX3((*var_fc).x3f, (*coarse_fc).x3f,
                                 pob->cis, pob->cie,
                                 pob->cjs, pob->cje,
                                 pob->cks, pob->cke+f3);
            FaceField &src_b = (*coarse_fc);
            FaceField &dst_b = *(std::get<0>(*pmb_fc_it)); // pmb->pfield->b;
            for (int k=kl, fk=pob->cks; fk<=pob->cke; k++, fk++) {
              for (int j=jl, fj=pob->cjs; fj<=pob->cje; j++, fj++) {
                for (int i=il, fi=pob->cis; fi<=pob->cie+1; i++, fi++)
                  dst_b.x1f(k, j, i) = src_b.x1f(fk, fj, fi);
              }
            }
            for (int k=kl, fk=pob->cks; fk<=pob->cke; k++, fk++) {
              for (int j=jl, fj=pob->cjs; fj<=pob->cje+f2; j++, fj++) {
                for (int i=il, fi=pob->cis; fi<=pob->cie; i++, fi++)
                  dst_b.x2f(k, j, i) = src_b.x2f(fk, fj, fi);
              }
            }
            if (pmb->block_size.nx2 == 1) {
              int iu = il + block_size.nx1/2 - 1;
              for (int i=il; i<=iu; i++)
                dst_b.x2f(pmb->ks, pmb->js+1, i) = dst_b.x2f(pmb->ks, pmb->js, i);
            }
            for (int k=kl, fk=pob->cks; fk<=pob->cke+f3; k++, fk++) {
              for (int j=jl, fj=pob->cjs; fj<=pob->cje; j++, fj++) {
                for (int i=il, fi=pob->cis; fi<=pob->cie; i++, fi++)
                  dst_b.x3f(k, j, i) = src_b.x3f(fk, fj, fi);
              }
            }
            if (pmb->block_size.nx3 == 1) {
              int iu = il + block_size.nx1/2-1, ju = jl + block_size.nx2/2-1;
              if (pmb->block_size.nx2 == 1) ju = jl;
              for (int j=jl; j<=ju; j++) {
                for (int i=il; i<=iu; i++)
                  dst_b.x3f(pmb->ks+1, j, i) = dst_b.x3f(pmb->ks, j, i);
              }
            }
            pmb_fc_it++;
          } // end loop over fc
        }
      } else if ((loclist[on].level < newloc[n].level) &&
                 (ranklist[on] == Globals::my_rank)) {
        // coarse to fine on the same node - prolongation
        MeshBlock* pob = FindMeshBlock(on);
        MeshRefinement *pmr = pmb->pmr;
        int il = pob->cis - 1, iu = pob->cie + 1, jl = pob->cjs - f2,
            ju = pob->cje + f2, kl = pob->cks - f3, ku = pob->cke + f3;
        int cis = ((newloc[n].lx1 & 1LL) == 1LL)*pob->block_size.nx1/2 + pob->is - 1;
        int cjs = ((newloc[n].lx2 & 1LL) == 1LL)*pob->block_size.nx2/2 + pob->js - f2;
        int cks = ((newloc[n].lx3 & 1LL) == 1LL)*pob->block_size.nx3/2 + pob->ks - f3;

        // KGF: CellCentered step 7: c2f, same node (just copy+prolongate, no pack/send)
        auto pob_cc_it = pob->pmr->pvars_cc_.begin();
        // iterate MeshRefinement std::vectors on new pmb
        for (auto cc_pair : pmr->pvars_cc_) {
          AthenaArray<Real> *var_cc = std::get<0>(cc_pair);
          AthenaArray<Real> *coarse_cc = std::get<1>(cc_pair);
          int nu = var_cc->GetDim4() - 1;

          AthenaArray<Real> &src = *(std::get<0>(*pob_cc_it)); // pob->phydro->u;
          AthenaArray<Real> &dst = (*coarse_cc); // pmb->phydro->coarse_cons_;
          // fill the coarse buffer
          for (int nv=0; nv<=nu; nv++) {
            for (int k=kl, ck=cks; k<=ku; k++, ck++) {
              for (int j=jl, cj=cjs; j<=ju; j++, cj++) {
                for (int i=il, ci=cis; i<=iu; i++, ci++)
                  dst(nv, k, j, i) = src(nv, ck, cj, ci);
              }
            }
          }
          pmr->ProlongateCellCenteredValues(
              dst, (*var_cc), 0, nu,
              pob->cis, pob->cie, pob->cjs, pob->cje, pob->cks, pob->cke);
          pob_cc_it++;
        } // end loop over var_cc

        // KGF: FaceCentered step 7: c2f, same node (just copy+prolongate, no pack/send)
        auto pob_fc_it = pob->pmr->pvars_fc_.begin();
        // iterate MeshRefinement std::vectors on new pmb
        for (auto fc_pair : pmr->pvars_fc_) {
          FaceField *var_fc = std::get<0>(fc_pair);
          FaceField *coarse_fc = std::get<1>(fc_pair);

          FaceField &src_b = *(std::get<0>(*pob_fc_it));
          FaceField &dst_b = (*coarse_fc);
          for (int k=kl, ck=cks; k<=ku; k++, ck++) {
            for (int j=jl, cj=cjs; j<=ju; j++, cj++) {
              for (int i=il, ci=cis; i<=iu+1; i++, ci++)
                dst_b.x1f(k, j, i) = src_b.x1f(ck, cj, ci);
            }
          }
          for (int k=kl, ck=cks; k<=ku; k++, ck++) {
            for (int j=jl, cj=cjs; j<=ju+f2; j++, cj++) {
              for (int i=il, ci=cis; i<=iu; i++, ci++)
                dst_b.x2f(k, j, i) = src_b.x2f(ck, cj, ci);
            }
          }
          for (int k=kl, ck=cks; k<=ku+f3; k++, ck++) {
            for (int j=jl, cj=cjs; j<=ju; j++, cj++) {
              for (int i=il, ci=cis; i<=iu; i++, ci++)
                dst_b.x3f(k, j, i) = src_b.x3f(ck, cj, ci);
            }
          }
          pmr->ProlongateSharedFieldX1(
              dst_b.x1f, (*var_fc).x1f,
              pob->cis, pob->cie+1, pob->cjs, pob->cje, pob->cks, pob->cke);
          pmr->ProlongateSharedFieldX2(
              dst_b.x2f, (*var_fc).x2f,
              pob->cis, pob->cie, pob->cjs, pob->cje+f2, pob->cks, pob->cke);
          pmr->ProlongateSharedFieldX3(
              dst_b.x3f, (*var_fc).x3f,
              pob->cis, pob->cie, pob->cjs, pob->cje, pob->cks, pob->cke+f3);
          pmr->ProlongateInternalField(
              (*var_fc), pob->cis, pob->cie,
              pob->cjs, pob->cje, pob->cks, pob->cke);
          pob_fc_it++;
        } // end loop over var_fc

        // KGF: what are these loops?
      }
    }
  }

  // discard remaining MeshBlocks
  // they could be reused, but for the moment, just throw them away for simplicity
  if (pblock != nullptr) {
    while (pblock->next  !=  nullptr)
      delete pblock->next;
    delete pblock;
  }

  // Replace the MeshBlock list
  pblock = newlist;

  // Step 8. Receive the data and load into MeshBlocks
  // This is a test: try MPI_Waitall later.
#ifdef MPI_PARALLEL
  if (nrecv != 0) {
    int rb_idx = 0;     // recv buffer index
    for (int n=nbs; n<=nbe; n++) {
      int on = newtoold[n];
      LogicalLocation &oloc = loclist[on];
      LogicalLocation &nloc = newloc[n];
      MeshBlock *pb = FindMeshBlock(n);
      if (oloc.level == nloc.level) { // same
        if (ranklist[on] == Globals::my_rank) continue;
        MPI_Wait(&(req_recv[rb_idx]), MPI_STATUS_IGNORE);
        int p = 0;
        // KGF: CellCentered step 8 (receive and load), branch 1 (same2same: unpack)
        for (auto cc_pair : pb->pmr->pvars_cc_) {
          AthenaArray<Real> *var_cc = std::get<0>(cc_pair);
          int nu = var_cc->GetDim4() - 1;
          BufferUtility::UnpackData(recvbuf[rb_idx], (*var_cc), 0, nu,
                                    pb->is, pb->ie, pb->js, pb->je, pb->ks, pb->ke, p);
        }
        // KGF: FaceCentered step 8 (receive and load), branch 1 (same2same: unpack)
        for (auto fc_pair : pb->pmr->pvars_fc_) {
          FaceField *var_fc = std::get<0>(fc_pair);
          FaceField &dst_b = (*var_fc);
          BufferUtility::UnpackData(
              recvbuf[rb_idx], dst_b.x1f,
              pb->is, pb->ie+1, pb->js, pb->je, pb->ks, pb->ke, p);
          BufferUtility::UnpackData(
              recvbuf[rb_idx], dst_b.x2f,
              pb->is, pb->ie, pb->js, pb->je+f2, pb->ks, pb->ke, p);
          BufferUtility::UnpackData(
              recvbuf[rb_idx], dst_b.x3f,
              pb->is, pb->ie, pb->js, pb->je, pb->ks, pb->ke+f3, p);
          if (pb->block_size.nx2 == 1) {
            for (int i=pb->is; i<=pb->ie; i++)
              dst_b.x2f(pb->ks, pb->js+1, i) = dst_b.x2f(pb->ks, pb->js, i);
          }
          if (pb->block_size.nx3 == 1) {
            for (int j=pb->js; j<=pb->je; j++) {
              for (int i=pb->is; i<=pb->ie; i++)
                dst_b.x3f(pb->ks+1, j, i) = dst_b.x3f(pb->ks, j, i);
            }
          }
        } // end loop over var_fc
        // KGF: dangerous? casting from "Real *" to "int *"
        int *dcp = reinterpret_cast<int *>(&(recvbuf[rb_idx][p]));
        pb->pmr->deref_count_ = *dcp;
        rb_idx++;
      } else if (oloc.level > nloc.level) { // f2c
        for (int l=0; l<nleaf; l++) {
          if (ranklist[on+l] == Globals::my_rank) continue;
          LogicalLocation &lloc = loclist[on+l];
          int ox1 = lloc.lx1 & 1LL, ox2 = lloc.lx2 & 1LL, ox3 = lloc.lx3 & 1LL;
          int p = 0, il, iu, jl, ju, kl, ku;
          if (ox1 == 0) il = pb->is,              iu = pb->is + pb->block_size.nx1/2 - 1;
          else        il = pb->is + pb->block_size.nx1/2, iu = pb->ie;
          if (ox2 == 0) jl = pb->js,              ju = pb->js + pb->block_size.nx2/2 - f2;
          else        jl = pb->js + pb->block_size.nx2/2, ju = pb->je;
          if (ox3 == 0) kl = pb->ks,              ku = pb->ks + pb->block_size.nx3/2 - f3;
          else        kl = pb->ks + pb->block_size.nx3/2, ku = pb->ke;
          MPI_Wait(&(req_recv[rb_idx]), MPI_STATUS_IGNORE);

          // KGF: CellCentered step 8 (receive and load), branch 2 (f2c: unpack)
          for (auto cc_pair : pb->pmr->pvars_cc_) {
            AthenaArray<Real> *var_cc = std::get<0>(cc_pair);
            int nu = var_cc->GetDim4() - 1;
            BufferUtility::UnpackData(recvbuf[rb_idx], (*var_cc), 0, nu,
                                      il, iu, jl, ju, kl, ku, p);
          }
          // KGF: FaceCentered step 8 (receive and load), branch 2 (f2c: unpack)
          for (auto fc_pair : pb->pmr->pvars_fc_) {
            FaceField *var_fc = std::get<0>(fc_pair);
            FaceField &dst_b = (*var_fc);
            BufferUtility::UnpackData(recvbuf[rb_idx], dst_b.x1f,
                                        il, iu+1, jl, ju, kl, ku, p);
            BufferUtility::UnpackData(recvbuf[rb_idx], dst_b.x2f,
                                        il, iu, jl, ju+f2, kl, ku, p);
            BufferUtility::UnpackData(recvbuf[rb_idx], dst_b.x3f,
                                        il, iu, jl, ju, kl, ku+f3, p);
            if (pb->block_size.nx2 == 1) {
              for (int i=il; i<=iu; i++)
                dst_b.x2f(pb->ks, pb->js+1, i) = dst_b.x2f(pb->ks, pb->js, i);
            }
            if (pb->block_size.nx3 == 1) {
              for (int j=jl; j<=ju; j++) {
                for (int i=il; i<=iu; i++)
                  dst_b.x3f(pb->ks+1, j, i) = dst_b.x3f(pb->ks, j, i);
              }
            }
          } // end loop over var_fc
          rb_idx++;
        }
      } else { // c2f
        if (ranklist[on] == Globals::my_rank) continue;
        MeshRefinement *pmr = pb->pmr;
        int p = 0;
        int il = pb->cis - 1, iu = pb->cie+1, jl = pb->cjs - f2,
            ju = pb->cje + f2, kl = pb->cks - f3, ku = pb->cke + f3;
        MPI_Wait(&(req_recv[rb_idx]), MPI_STATUS_IGNORE);
        // KGF: CellCentered step 8 (receive and load), branch 2 (c2f: unpack+prolongate)
        for (auto cc_pair : pb->pmr->pvars_cc_) {
          AthenaArray<Real> *var_cc = std::get<0>(cc_pair);
          AthenaArray<Real> *coarse_cc = std::get<1>(cc_pair);
          int nu = var_cc->GetDim4() - 1;
          BufferUtility::UnpackData(recvbuf[rb_idx], (*coarse_cc),
                                    0, nu, il, iu, jl, ju, kl, ku, p);
          pmr->ProlongateCellCenteredValues(
              (*coarse_cc), (*var_cc), 0, nu,
              pb->cis, pb->cie, pb->cjs, pb->cje, pb->cks, pb->cke);
        }
        // KGF: FaceCentered step 8 (receive and load), branch 2 (c2f: unpack+prolongate)
        for (auto fc_pair : pb->pmr->pvars_fc_) {
          FaceField *var_fc = std::get<0>(fc_pair);
          FaceField *coarse_fc = std::get<1>(fc_pair);

          BufferUtility::UnpackData(recvbuf[rb_idx], (*coarse_fc).x1f,
                                      il, iu+1, jl, ju, kl, ku, p);
          BufferUtility::UnpackData(recvbuf[rb_idx], (*coarse_fc).x2f,
                                      il, iu, jl, ju+f2, kl, ku, p);
          BufferUtility::UnpackData(recvbuf[rb_idx], (*coarse_fc).x3f,
                                      il, iu, jl, ju, kl, ku+f3, p);
          pmr->ProlongateSharedFieldX1(
              (*coarse_fc).x1f, (*var_fc).x1f,
              pb->cis, pb->cie+1, pb->cjs, pb->cje, pb->cks, pb->cke);
          pmr->ProlongateSharedFieldX2(
              (*coarse_fc).x2f, (*var_fc).x2f,
              pb->cis, pb->cie, pb->cjs, pb->cje+f2, pb->cks, pb->cke);
          pmr->ProlongateSharedFieldX3(
              (*coarse_fc).x3f, (*var_fc).x3f,
              pb->cis, pb->cie, pb->cjs, pb->cje, pb->cks, pb->cke+f3);
          pmr->ProlongateInternalField(
              (*var_fc), pb->cis, pb->cie,
              pb->cjs, pb->cje, pb->cks, pb->cke);
        }
        rb_idx++;
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
  if (nsend != 0) {
    MPI_Waitall(nsend, req_send, MPI_STATUSES_IGNORE);
    for (int n=0; n<nsend; n++)
      delete [] sendbuf[n];
    delete [] sendbuf;
    delete [] req_send;
  }
  if (nrecv != 0) {
    for (int n=0; n<nrecv; n++)
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
  pmb = pblock;
  while (pmb != nullptr) {
    pmb->pbval->SearchAndSetNeighbors(tree, ranklist, nslist);
    pmb = pmb->next;
  }
  Initialize(2, pin);

  return;
}

//----------------------------------------------------------------------------------------
//! \fn int CreateAMRMPITag(int lid, int ox1, int ox2, int ox3)
//  \brief calculate an MPI tag for AMR block transfer
// tag = local id of destination (remaining bits) + ox1(1 bit) + ox2(1 bit) + ox3(1 bit)
//       + physics(5 bits)

// See comments on BoundaryBase::CreateBvalsMPITag()

int Mesh::CreateAMRMPITag(int lid, int ox1, int ox2, int ox3) {
  // KGF: former "AthenaTagMPI" AthenaTagMPI::amr=8 redefined to 0
  return (lid<<8) | (ox1<<7)| (ox2<<6) | (ox3<<5) | 0;
}
