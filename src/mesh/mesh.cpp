//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mesh.cpp
//! \brief implementation of functions in Mesh class

// C headers
// pre-C11: needed before including inttypes.h, else won't define int64_t for C++ code
// #define __STDC_FORMAT_MACROS

// C++ headers
#include <algorithm>
#include <cinttypes>  // format macro "PRId64" for fixed-width integer type std::int64_t
#include <cmath>      // std::abs(), std::pow()
#include <cstdint>    // std::int64_t fixed-wdith integer type alias
#include <cstdlib>
#include <cstring>    // std::memcpy()
#include <iomanip>    // std::setprecision()
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
#include "../orbital_advection/orbital_advection.hpp"
#include "../outputs/io_wrapper.hpp"
#include "../parameter_input.hpp"
#include "../reconstruct/reconstruction.hpp"
#include "../scalars/scalars.hpp"
#include "../utils/buffer_utils.hpp"
#include "mesh.hpp"
#include "mesh_refinement.hpp"
#include "meshblock_tree.hpp"

// MPI/OpenMP header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

//----------------------------------------------------------------------------------------
//! Mesh constructor, builds mesh at start of calculation using parameters in input file

Mesh::Mesh(ParameterInput *pin, int mesh_test) :
    // public members:
    // aggregate initialization of RegionSize struct:
    mesh_size{pin->GetReal("mesh", "x1min"), pin->GetReal("mesh", "x2min"),
              pin->GetReal("mesh", "x3min"), pin->GetReal("mesh", "x1max"),
              pin->GetReal("mesh", "x2max"), pin->GetReal("mesh", "x3max"),
              pin->GetOrAddReal("mesh", "x1rat", 1.0),
              pin->GetOrAddReal("mesh", "x2rat", 1.0),
              pin->GetOrAddReal("mesh", "x3rat", 1.0),
              pin->GetInteger("mesh", "nx1"), pin->GetInteger("mesh", "nx2"),
              pin->GetInteger("mesh", "nx3") },
    mesh_bcs{GetBoundaryFlag(pin->GetOrAddString("mesh", "ix1_bc", "none")),
             GetBoundaryFlag(pin->GetOrAddString("mesh", "ox1_bc", "none")),
             GetBoundaryFlag(pin->GetOrAddString("mesh", "ix2_bc", "none")),
             GetBoundaryFlag(pin->GetOrAddString("mesh", "ox2_bc", "none")),
             GetBoundaryFlag(pin->GetOrAddString("mesh", "ix3_bc", "none")),
             GetBoundaryFlag(pin->GetOrAddString("mesh", "ox3_bc", "none"))},
    f2(mesh_size.nx2 > 1 ? true : false), f3(mesh_size.nx3 > 1 ? true : false),
    ndim(f3 ? 3 : (f2 ? 2 : 1)),
    adaptive(pin->GetOrAddString("mesh", "refinement", "none") == "adaptive"
             ? true : false),
    multilevel((adaptive || pin->GetOrAddString("mesh", "refinement", "none") == "static")
               ? true : false),
    orbital_advection(pin->GetOrAddInteger("orbital_advection","OAorder",0)),
    shear_periodic(GetBoundaryFlag(pin->GetOrAddString("mesh", "ix1_bc", "none"))
                   == BoundaryFlag::shear_periodic ? true : false),
    fluid_setup(GetFluidFormulation(pin->GetOrAddString("hydro", "active", "true"))),
    start_time(pin->GetOrAddReal("time", "start_time", 0.0)), time(start_time),
    tlim(pin->GetReal("time", "tlim")), dt(std::numeric_limits<Real>::max()),
    dt_hyperbolic(dt), dt_parabolic(dt), dt_user(dt),
    cfl_number(pin->GetReal("time", "cfl_number")),
    nlim(pin->GetOrAddInteger("time", "nlim", -1)), ncycle(),
    ncycle_out(pin->GetOrAddInteger("time", "ncycle_out", 1)),
    dt_diagnostics(pin->GetOrAddInteger("time", "dt_diagnostics", -1)),
    sts_integrator(pin->GetOrAddString("time", "sts_integrator", "rkl2")),
    sts_max_dt_ratio(pin->GetOrAddReal("time", "sts_max_dt_ratio", -1.0)),
    sts_loc(TaskType::main_int),
    muj(), nuj(), muj_tilde(), gammaj_tilde(),
    nbnew(), nbdel(),
    step_since_lb(), gflag(), turb_flag(), amr_updated(multilevel),
    // private members:
    next_phys_id_(), num_mesh_threads_(pin->GetOrAddInteger("mesh", "num_threads", 1)),
    gids_(), gide_(),
    tree(this),
    use_uniform_meshgen_fn_{true, true, true},
    nreal_user_mesh_data_(), nint_user_mesh_data_(), nuser_history_output_(),
    four_pi_G_(), grav_eps_(-1.0),
    lb_flag_(true), lb_automatic_(), lb_manual_(),
    MeshGenerator_{UniformMeshGeneratorX1, UniformMeshGeneratorX2,
                   UniformMeshGeneratorX3},
    BoundaryFunction_{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    AMRFlag_{}, UserSourceTerm_{}, UserTimeStep_{}, ViscosityCoeff_{},
    ConductionCoeff_{}, FieldDiffusivity_{},
    OrbitalVelocity_{}, OrbitalVelocityDerivative_{nullptr, nullptr},
    MGGravityBoundaryFunction_{MGPeriodicInnerX1, MGPeriodicOuterX1, MGPeriodicInnerX2,
                               MGPeriodicOuterX2, MGPeriodicInnerX3, MGPeriodicOuterX3} {
  std::stringstream msg;
  RegionSize block_size;
  MeshBlock *pfirst{};
  BoundaryFlag block_bcs[6];
  std::int64_t nbmax;

  // mesh test
  if (mesh_test > 0) Globals::nranks = mesh_test;

#ifdef MPI_PARALLEL
  // reserve phys=0 for former TAG_AMR=8; now hard-coded in Mesh::CreateAMRMPITag()
  next_phys_id_  = 1;
  ReserveMeshBlockPhysIDs();
#endif

  // check number of OpenMP threads for mesh
  if (num_mesh_threads_ < 1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Number of OpenMP threads must be >= 1, but num_threads="
        << num_mesh_threads_ << std::endl;
    ATHENA_ERROR(msg);
  }

  // check number of grid cells in root level of mesh from input file.
  if (mesh_size.nx1 < 4) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "In mesh block in input file nx1 must be >= 4, but nx1="
        << mesh_size.nx1 << std::endl;
    ATHENA_ERROR(msg);
  }
  if (mesh_size.nx2 < 1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "In mesh block in input file nx2 must be >= 1, but nx2="
        << mesh_size.nx2 << std::endl;
    ATHENA_ERROR(msg);
  }
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

  // check physical size of mesh (root level) from input file.
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

  // check the consistency of the periodic boundaries
  if ( ((mesh_bcs[BoundaryFace::inner_x1] == BoundaryFlag::periodic
    &&   mesh_bcs[BoundaryFace::outer_x1] != BoundaryFlag::periodic)
    ||  (mesh_bcs[BoundaryFace::inner_x1] != BoundaryFlag::periodic
    &&   mesh_bcs[BoundaryFace::outer_x1] == BoundaryFlag::periodic))
    ||  (mesh_size.nx2 > 1
    && ((mesh_bcs[BoundaryFace::inner_x2] == BoundaryFlag::periodic
    &&   mesh_bcs[BoundaryFace::outer_x2] != BoundaryFlag::periodic)
    ||  (mesh_bcs[BoundaryFace::inner_x2] != BoundaryFlag::periodic
    &&   mesh_bcs[BoundaryFace::outer_x2] == BoundaryFlag::periodic)))
    ||  (mesh_size.nx3 > 1
    && ((mesh_bcs[BoundaryFace::inner_x3] == BoundaryFlag::periodic
    &&   mesh_bcs[BoundaryFace::outer_x3] != BoundaryFlag::periodic)
    ||  (mesh_bcs[BoundaryFace::inner_x3] != BoundaryFlag::periodic
    &&   mesh_bcs[BoundaryFace::outer_x3] == BoundaryFlag::periodic)))) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "When periodic boundaries are in use, both sides must be periodic."
        << std::endl;
    ATHENA_ERROR(msg);
  }
  if ( ((mesh_bcs[BoundaryFace::inner_x1] == BoundaryFlag::shear_periodic
    &&   mesh_bcs[BoundaryFace::outer_x1] != BoundaryFlag::shear_periodic)
    ||  (mesh_bcs[BoundaryFace::inner_x1] != BoundaryFlag::shear_periodic
    &&   mesh_bcs[BoundaryFace::outer_x1] == BoundaryFlag::shear_periodic))) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "When shear_periodic boundaries are in use, "
        << "both sides must be shear_periodic." << std::endl;
    ATHENA_ERROR(msg);
  }

  // read and set MeshBlock parameters
  block_size.x1rat = mesh_size.x1rat;
  block_size.x2rat = mesh_size.x2rat;
  block_size.x3rat = mesh_size.x3rat;
  block_size.nx1 = pin->GetOrAddInteger("meshblock", "nx1", mesh_size.nx1);
  if (f2)
    block_size.nx2 = pin->GetOrAddInteger("meshblock", "nx2", mesh_size.nx2);
  else
    block_size.nx2 = mesh_size.nx2;
  if (f3)
    block_size.nx3 = pin->GetOrAddInteger("meshblock", "nx3", mesh_size.nx3);
  else
    block_size.nx3 = mesh_size.nx3;

  // check consistency of the block and mesh
  if (mesh_size.nx1 % block_size.nx1 != 0
      || mesh_size.nx2 % block_size.nx2 != 0
      || mesh_size.nx3 % block_size.nx3 != 0) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "the Mesh must be evenly divisible by the MeshBlock" << std::endl;
    ATHENA_ERROR(msg);
  }
  if (multilevel) { // SMR/AMR
    // It is required to know xorder before pmb->precon is defined.
    // This assumes that the stencil size is equal to xorder
    int xorder = 1;
    std::string input_recon = pin->GetOrAddString("time", "xorder", "2");
    if ((input_recon == "2") || (input_recon == "2c")) {
      xorder = 2;
    } else if ((input_recon == "3") || (input_recon == "3c")) {
      xorder = 3;
    } else if ((input_recon == "4") || (input_recon == "4c")) {
      xorder = 4;
    }
    int slimit = std::max(2*(xorder-1), NGHOST);
    if (block_size.nx1 < slimit || (block_size.nx2 < slimit && f2)
        || (block_size.nx3 < slimit && f3)) {
      msg << "### FATAL ERROR in Mesh constructor" << std::endl
          << "block_size must be larger than or equal to " << slimit
          << " cells with SMR/AMR, when xorder = " << input_recon
          << " and NGHOST = " << NGHOST << ". "<<std::endl;
      ATHENA_ERROR(msg);
    }
  } else {
    if (block_size.nx1 < NGHOST || (block_size.nx2 < NGHOST && f2)
        || (block_size.nx3 < NGHOST && f3)) {
      msg << "### FATAL ERROR in Mesh constructor" << std::endl
          << "block_size must be larger than or equal to NGHOST = " << NGHOST
          << " cells with uniform grid." << std::endl;
      ATHENA_ERROR(msg);
    }
  }

  // calculate the number of the blocks
  nrbx1 = mesh_size.nx1/block_size.nx1;
  nrbx2 = mesh_size.nx2/block_size.nx2;
  nrbx3 = mesh_size.nx3/block_size.nx3;
  nbmax = (nrbx1 > nrbx2) ? nrbx1:nrbx2;
  nbmax = (nbmax > nrbx3) ? nbmax:nrbx3;

  // initialize user-enrollable functions
  if (mesh_size.x1rat != 1.0) {
    use_uniform_meshgen_fn_[X1DIR] = false;
    MeshGenerator_[X1DIR] = DefaultMeshGeneratorX1;
  }
  if (mesh_size.x2rat != 1.0) {
    use_uniform_meshgen_fn_[X2DIR] = false;
    MeshGenerator_[X2DIR] = DefaultMeshGeneratorX2;
  }
  if (mesh_size.x3rat != 1.0) {
    use_uniform_meshgen_fn_[X3DIR] = false;
    MeshGenerator_[X3DIR] = DefaultMeshGeneratorX3;
  }

  // calculate the logical root level and maximum level
  for (root_level=0; (1<<root_level) < nbmax; root_level++) {}
  current_level = root_level;

  tree.CreateRootGrid();

  // Load balancing flag and parameters
#ifdef MPI_PARALLEL
  if (pin->GetOrAddString("loadbalancing","balancer","default") == "automatic")
    lb_automatic_ = true;
  else if (pin->GetOrAddString("loadbalancing","balancer","default") == "manual")
    lb_manual_ = true;
  lb_tolerance_ = pin->GetOrAddReal("loadbalancing","tolerance",0.5);
  lb_interval_ = pin->GetOrAddReal("loadbalancing","interval",10);
#endif

  // SMR / AMR:
  if (adaptive) {
    max_level = pin->GetOrAddInteger("mesh", "numlevel", 1) + root_level - 1;
    if (max_level > 63) {
      msg << "### FATAL ERROR in Mesh constructor" << std::endl
          << "The number of the refinement level must be smaller than "
          << 63 - root_level + 1 << "." << std::endl;
      ATHENA_ERROR(msg);
    }
  } else {
    max_level = 63;
  }

  if (EOS_TABLE_ENABLED) peos_table = new EosTable(pin);
  InitUserMeshData(pin);

  if (multilevel) {
    if (block_size.nx1 % 2 == 1 || (block_size.nx2 % 2 == 1 && f2)
        || (block_size.nx3 % 2 == 1 && f3)) {
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
        if (f2) {
          ref_size.x2min = pin->GetReal(pib->block_name, "x2min");
          ref_size.x2max = pin->GetReal(pib->block_name, "x2max");
        } else {
          ref_size.x2min = mesh_size.x2min;
          ref_size.x2max = mesh_size.x2max;
        }
        if (ndim == 3) {
          ref_size.x3min = pin->GetReal(pib->block_name, "x3min");
          ref_size.x3max = pin->GetReal(pib->block_name, "x3max");
        } else {
          ref_size.x3min = mesh_size.x3min;
          ref_size.x3max = mesh_size.x3max;
        }
        int ref_lev = pin->GetInteger(pib->block_name, "level");
        int lrlev = ref_lev + root_level;
        if (lrlev > current_level) current_level = lrlev;
        // range check
        if (ref_lev < 1) {
          msg << "### FATAL ERROR in Mesh constructor" << std::endl
              << "Refinement level must be larger than 0 (root level = 0)" << std::endl;
          ATHENA_ERROR(msg);
        }
        if (lrlev > max_level) {
          msg << "### FATAL ERROR in Mesh constructor" << std::endl
              << "Refinement level exceeds the maximum level (specify "
              << "'numlevel' parameter in <mesh> input block if adaptive)."
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
        std::int64_t lx1min = 0, lx1max = 0, lx2min = 0, lx2max = 0,
                     lx3min = 0, lx3max = 0;
        std::int64_t lxmax = nrbx1*(1LL<<ref_lev);
        for (lx1min=0; lx1min<lxmax; lx1min++) {
          Real rx = ComputeMeshGeneratorX(lx1min+1, lxmax,
                                          use_uniform_meshgen_fn_[X1DIR]);
          if (MeshGenerator_[X1DIR](rx, mesh_size) > ref_size.x1min)
            break;
        }
        for (lx1max=lx1min; lx1max<lxmax; lx1max++) {
          Real rx = ComputeMeshGeneratorX(lx1max+1, lxmax,
                                          use_uniform_meshgen_fn_[X1DIR]);
          if (MeshGenerator_[X1DIR](rx, mesh_size) >= ref_size.x1max)
            break;
        }
        if (lx1min % 2 == 1) lx1min--;
        if (lx1max % 2 == 0) lx1max++;
        if (f2) { // 2D or 3D
          lxmax = nrbx2*(1LL << ref_lev);
          for (lx2min=0; lx2min<lxmax; lx2min++) {
            Real rx = ComputeMeshGeneratorX(lx2min+1, lxmax,
                                            use_uniform_meshgen_fn_[X2DIR]);
            if (MeshGenerator_[X2DIR](rx, mesh_size) > ref_size.x2min)
              break;
          }
          for (lx2max=lx2min; lx2max<lxmax; lx2max++) {
            Real rx = ComputeMeshGeneratorX(lx2max+1, lxmax,
                                            use_uniform_meshgen_fn_[X2DIR]);
            if (MeshGenerator_[X2DIR](rx, mesh_size) >= ref_size.x2max)
              break;
          }
          if (lx2min % 2 == 1) lx2min--;
          if (lx2max % 2 == 0) lx2max++;
        }
        if (ndim == 3) { // 3D
          lxmax = nrbx3*(1LL<<ref_lev);
          for (lx3min=0; lx3min<lxmax; lx3min++) {
            Real rx = ComputeMeshGeneratorX(lx3min+1, lxmax,
                                            use_uniform_meshgen_fn_[X3DIR]);
            if (MeshGenerator_[X3DIR](rx, mesh_size) > ref_size.x3min)
              break;
          }
          for (lx3max=lx3min; lx3max<lxmax; lx3max++) {
            Real rx = ComputeMeshGeneratorX(lx3max+1, lxmax,
                                            use_uniform_meshgen_fn_[X3DIR]);
            if (MeshGenerator_[X3DIR](rx, mesh_size) >= ref_size.x3max)
              break;
          }
          if (lx3min % 2 == 1) lx3min--;
          if (lx3max % 2 == 0) lx3max++;
        }
        // create the finest level
        if (ndim == 1) {
          for (std::int64_t i=lx1min; i<lx1max; i+=2) {
            LogicalLocation nloc;
            nloc.level=lrlev, nloc.lx1=i, nloc.lx2=0, nloc.lx3=0;
            int nnew;
            tree.AddMeshBlock(nloc, nnew);
          }
        }
        if (ndim == 2) {
          for (std::int64_t j=lx2min; j<lx2max; j+=2) {
            for (std::int64_t i=lx1min; i<lx1max; i+=2) {
              LogicalLocation nloc;
              nloc.level=lrlev, nloc.lx1=i, nloc.lx2=j, nloc.lx3=0;
              int nnew;
              tree.AddMeshBlock(nloc, nnew);
            }
          }
        }
        if (ndim == 3) {
          for (std::int64_t k=lx3min; k<lx3max; k+=2) {
            for (std::int64_t j=lx2min; j<lx2max; j+=2) {
              for (std::int64_t i=lx1min; i<lx1max; i+=2) {
                LogicalLocation nloc;
                nloc.level = lrlev, nloc.lx1 = i, nloc.lx2 = j, nloc.lx3 = k;
                int nnew;
                tree.AddMeshBlock(nloc, nnew);
              }
            }
          }
        }
      }
      pib = pib->pnext;
    }
  }

  if (!adaptive) max_level = current_level;

  // initial mesh hierarchy construction is completed here
  tree.CountMeshBlock(nbtotal);
  loclist = new LogicalLocation[nbtotal];
  tree.GetMeshBlockList(loclist, nullptr, nbtotal);

#ifdef MPI_PARALLEL
  // check if there are sufficient blocks
  if (nbtotal < Globals::nranks) {
    if (mesh_test == 0) {
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

  ranklist = new int[nbtotal];
  nslist = new int[Globals::nranks];
  nblist = new int[Globals::nranks];
  costlist = new double[nbtotal];
  if (adaptive) { // allocate arrays for AMR
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
  for (int i=0; i<nbtotal; i++) costlist[i] = 1.0;

  CalculateLoadBalance(costlist, ranklist, nslist, nblist, nbtotal);

  // Output some diagnostic information to terminal

  // Output MeshBlock list and quit (mesh test only); do not create meshes
  if (mesh_test > 0) {
    if (Globals::my_rank == 0) OutputMeshStructure(ndim);
    return;
  }


  if (SELF_GRAVITY_ENABLED == 1) {
    gflag = 1; // set gravity flag
    pfgrd = new FFTGravityDriver(this, pin);
  } else if (SELF_GRAVITY_ENABLED == 2) {
    // MGDriver must be initialzied before MeshBlocks
    pmgrd = new MGGravityDriver(this, pin);
  }
  //  if (SELF_GRAVITY_ENABLED == 2 && ...) // independent allocation
  //    gflag = 2;

  // create MeshBlock list for this process
  gids_ = nslist[Globals::my_rank];
  gide_ = gids_ + nblist[Globals::my_rank] - 1;
  nblocal = nblist[Globals::my_rank];
  my_blocks.NewAthenaArray(nblocal);
  // create MeshBlocks for this node
  for (int i=gids_; i<=gide_; i++) {
    SetBlockSizeAndBoundaries(loclist[i], block_size, block_bcs);
    my_blocks(i-gids_) = new MeshBlock(i, i-gids_, loclist[i], block_size, block_bcs,
                                       this, pin, gflag);
    my_blocks(i-gids_)->pbval->SearchAndSetNeighbors(tree, ranklist, nslist);
  }

  ResetLoadBalanceVariables();

  if (turb_flag > 0) // TurbulenceDriver depends on the MeshBlock ctor
    ptrbd = new TurbulenceDriver(this, pin);
}

//----------------------------------------------------------------------------------------
//! Mesh constructor for restarts. Load the restart file

Mesh::Mesh(ParameterInput *pin, IOWrapper& resfile, int mesh_test) :
    // public members:
    // aggregate initialization of RegionSize struct:
    // (will be overwritten by memcpy from restart file, in this case)
    mesh_size{pin->GetReal("mesh", "x1min"), pin->GetReal("mesh", "x2min"),
              pin->GetReal("mesh", "x3min"), pin->GetReal("mesh", "x1max"),
              pin->GetReal("mesh", "x2max"), pin->GetReal("mesh", "x3max"),
              pin->GetOrAddReal("mesh", "x1rat", 1.0),
              pin->GetOrAddReal("mesh", "x2rat", 1.0),
              pin->GetOrAddReal("mesh", "x3rat", 1.0),
              pin->GetInteger("mesh", "nx1"), pin->GetInteger("mesh", "nx2"),
              pin->GetInteger("mesh", "nx3") },
    mesh_bcs{GetBoundaryFlag(pin->GetOrAddString("mesh", "ix1_bc", "none")),
             GetBoundaryFlag(pin->GetOrAddString("mesh", "ox1_bc", "none")),
             GetBoundaryFlag(pin->GetOrAddString("mesh", "ix2_bc", "none")),
             GetBoundaryFlag(pin->GetOrAddString("mesh", "ox2_bc", "none")),
             GetBoundaryFlag(pin->GetOrAddString("mesh", "ix3_bc", "none")),
             GetBoundaryFlag(pin->GetOrAddString("mesh", "ox3_bc", "none"))},
    f2(mesh_size.nx2 > 1 ? true : false), f3(mesh_size.nx3 > 1 ? true : false),
    ndim(f3 ? 3 : (f2 ? 2 : 1)),
    adaptive(pin->GetOrAddString("mesh", "refinement", "none") == "adaptive"
             ? true : false),
    multilevel((adaptive || pin->GetOrAddString("mesh", "refinement", "none") == "static")
               ? true : false),
    orbital_advection(pin->GetOrAddInteger("orbital_advection","OAorder",0)),
    shear_periodic(GetBoundaryFlag(pin->GetOrAddString("mesh", "ix1_bc", "none"))
                   == BoundaryFlag::shear_periodic ? true : false),
    fluid_setup(GetFluidFormulation(pin->GetOrAddString("hydro", "active", "true"))),
    start_time(pin->GetOrAddReal("time", "start_time", 0.0)), time(start_time),
    tlim(pin->GetReal("time", "tlim")), dt(std::numeric_limits<Real>::max()),
    dt_hyperbolic(dt), dt_parabolic(dt), dt_user(dt),
    cfl_number(pin->GetReal("time", "cfl_number")),
    nlim(pin->GetOrAddInteger("time", "nlim", -1)), ncycle(),
    ncycle_out(pin->GetOrAddInteger("time", "ncycle_out", 1)),
    dt_diagnostics(pin->GetOrAddInteger("time", "dt_diagnostics", -1)),
    sts_integrator(pin->GetOrAddString("time", "sts_integrator", "rkl2")),
    sts_max_dt_ratio(pin->GetOrAddReal("time", "sts_max_dt_ratio", -1.0)),
    sts_loc(TaskType::main_int),
    muj(), nuj(), muj_tilde(), gammaj_tilde(),
    nbnew(), nbdel(),
    step_since_lb(), gflag(), turb_flag(), amr_updated(multilevel),
    // private members:
    next_phys_id_(), num_mesh_threads_(pin->GetOrAddInteger("mesh", "num_threads", 1)),
    gids_(), gide_(),
    tree(this),
    use_uniform_meshgen_fn_{true, true, true},
    nreal_user_mesh_data_(), nint_user_mesh_data_(), nuser_history_output_(),
    four_pi_G_(), grav_eps_(-1.0),
    lb_flag_(true), lb_automatic_(), lb_manual_(),
    MeshGenerator_{UniformMeshGeneratorX1, UniformMeshGeneratorX2,
                   UniformMeshGeneratorX3},
    BoundaryFunction_{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    AMRFlag_{}, UserSourceTerm_{}, UserTimeStep_{}, ViscosityCoeff_{},
    ConductionCoeff_{}, FieldDiffusivity_{},
    OrbitalVelocity_{}, OrbitalVelocityDerivative_{nullptr, nullptr},
    MGGravityBoundaryFunction_{MGPeriodicInnerX1, MGPeriodicOuterX1, MGPeriodicInnerX2,
                        MGPeriodicOuterX2, MGPeriodicInnerX3, MGPeriodicOuterX3} {
  std::stringstream msg;
  RegionSize block_size;
  BoundaryFlag block_bcs[6];
  MeshBlock *pfirst{};
  IOWrapperSizeT *offset{};
  IOWrapperSizeT datasize, listsize, headeroffset;

  // mesh test
  if (mesh_test > 0) Globals::nranks = mesh_test;

#ifdef MPI_PARALLEL
  // reserve phys=0 for former TAG_AMR=8; now hard-coded in Mesh::CreateAMRMPITag()
  next_phys_id_  = 1;
  ReserveMeshBlockPhysIDs();
#endif

  // check the number of OpenMP threads for mesh
  if (num_mesh_threads_ < 1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Number of OpenMP threads must be >= 1, but num_threads="
        << num_mesh_threads_ << std::endl;
    ATHENA_ERROR(msg);
  }

  // get the end of the header
  headeroffset = resfile.GetPosition();
  // read the restart file
  // the file is already open and the pointer is set to after <par_end>
  IOWrapperSizeT headersize = sizeof(int)*3+sizeof(Real)*2
                              + sizeof(RegionSize)+sizeof(IOWrapperSizeT);
  char *headerdata = new char[headersize];
  if (Globals::my_rank == 0) { // the master process reads the header data
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
  current_level = root_level;
  std::memcpy(&mesh_size, &(headerdata[hdos]), sizeof(RegionSize));
  hdos += sizeof(RegionSize);
  std::memcpy(&time, &(headerdata[hdos]), sizeof(Real));
  hdos += sizeof(Real);
  std::memcpy(&dt, &(headerdata[hdos]), sizeof(Real));
  hdos += sizeof(Real);
  std::memcpy(&ncycle, &(headerdata[hdos]), sizeof(int));
  hdos += sizeof(int);
  std::memcpy(&datasize, &(headerdata[hdos]), sizeof(IOWrapperSizeT));
  hdos += sizeof(IOWrapperSizeT);   // (this updated value is never used)

  delete [] headerdata;

  // initialize
  loclist = new LogicalLocation[nbtotal];
  offset = new IOWrapperSizeT[nbtotal];
  costlist = new double[nbtotal];
  ranklist = new int[nbtotal];
  nslist = new int[Globals::nranks];
  nblist = new int[Globals::nranks];

  block_size.nx1 = pin->GetOrAddInteger("meshblock", "nx1", mesh_size.nx1);
  block_size.nx2 = pin->GetOrAddInteger("meshblock", "nx2", mesh_size.nx2);
  block_size.nx3 = pin->GetOrAddInteger("meshblock", "nx3", mesh_size.nx3);

  // calculate the number of the blocks
  nrbx1 = mesh_size.nx1/block_size.nx1;
  nrbx2 = mesh_size.nx2/block_size.nx2;
  nrbx3 = mesh_size.nx3/block_size.nx3;

  // initialize user-enrollable functions
  if (mesh_size.x1rat != 1.0) {
    use_uniform_meshgen_fn_[X1DIR] = false;
    MeshGenerator_[X1DIR] = DefaultMeshGeneratorX1;
  }
  if (mesh_size.x2rat != 1.0) {
    use_uniform_meshgen_fn_[X2DIR] = false;
    MeshGenerator_[X2DIR] = DefaultMeshGeneratorX2;
  }
  if (mesh_size.x3rat != 1.0) {
    use_uniform_meshgen_fn_[X3DIR] = false;
    MeshGenerator_[X3DIR] = DefaultMeshGeneratorX3;
  }

  // Load balancing flag and parameters
#ifdef MPI_PARALLEL
  if (pin->GetOrAddString("loadbalancing", "balancer", "default") == "automatic")
    lb_automatic_ = true;
  else if (pin->GetOrAddString("loadbalancing", "balancer", "default") == "manual")
    lb_manual_ = true;
  lb_tolerance_ = pin->GetOrAddReal("loadbalancing", "tolerance", 0.5);
  lb_interval_ = pin->GetOrAddReal("loadbalancing", "interval", 10);
#endif

  // SMR / AMR
  if (adaptive) {
    max_level = pin->GetOrAddInteger("mesh", "numlevel", 1) + root_level - 1;
    if (max_level > 63) {
      msg << "### FATAL ERROR in Mesh constructor" << std::endl
          << "The number of the refinement level must be smaller than "
          << 63 - root_level + 1 << "." << std::endl;
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
    if (Globals::my_rank == 0) { // only the master process reads the ID list
      if (resfile.Read(userdata, 1, udsize) != udsize) {
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
  listsize = sizeof(LogicalLocation)+sizeof(Real);
  //allocate the idlist buffer
  char *idlist = new char[listsize*nbtotal];
  if (Globals::my_rank == 0) { // only the master process reads the ID list
    if (resfile.Read(idlist, listsize, nbtotal) != static_cast<unsigned int>(nbtotal)) {
      msg << "### FATAL ERROR in Mesh constructor" << std::endl
          << "The restart file is broken." << std::endl;
      ATHENA_ERROR(msg);
    }
  }
#ifdef MPI_PARALLEL
  // then broadcast the ID list
  MPI_Bcast(idlist, listsize*nbtotal, MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

  int os = 0;
  for (int i=0; i<nbtotal; i++) {
    std::memcpy(&(loclist[i]), &(idlist[os]), sizeof(LogicalLocation));
    os += sizeof(LogicalLocation);
    std::memcpy(&(costlist[i]), &(idlist[os]), sizeof(double));
    os += sizeof(double);
    if (loclist[i].level > current_level) current_level = loclist[i].level;
  }
  delete [] idlist;

  if (!adaptive) max_level = current_level;

  // calculate the header offset and seek
  headeroffset += headersize + udsize + listsize*nbtotal;
  if (Globals::my_rank != 0)
    resfile.Seek(headeroffset);

  // rebuild the Block Tree
  tree.CreateRootGrid();
  for (int i=0; i<nbtotal; i++)
    tree.AddMeshBlockWithoutRefine(loclist[i]);
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
    if (mesh_test == 0) {
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

  if (adaptive) { // allocate arrays for AMR
    nref = new int[Globals::nranks];
    nderef = new int[Globals::nranks];
    rdisp = new int[Globals::nranks];
    ddisp = new int[Globals::nranks];
    bnref = new int[Globals::nranks];
    bnderef = new int[Globals::nranks];
    brdisp = new int[Globals::nranks];
    bddisp = new int[Globals::nranks];
  }

  CalculateLoadBalance(costlist, ranklist, nslist, nblist, nbtotal);

  // Output MeshBlock list and quit (mesh test only); do not create meshes
  if (mesh_test > 0) {
    if (Globals::my_rank == 0) OutputMeshStructure(ndim);
    delete [] offset;
    return;
  }

  if (SELF_GRAVITY_ENABLED == 1) {
    gflag = 1; // set gravity flag
    pfgrd = new FFTGravityDriver(this, pin);
  } else if (SELF_GRAVITY_ENABLED == 2) {
    // MGDriver must be initialzied before MeshBlocks
    pmgrd = new MGGravityDriver(this, pin);
  }
  //  if (SELF_GRAVITY_ENABLED == 2 && ...) // independent allocation
  //    gflag=2;

  // allocate data buffer
  nblocal = nblist[Globals::my_rank];
  gids_ = nslist[Globals::my_rank];
  gide_ = gids_ + nblocal - 1;
  char *mbdata = new char[datasize*nblocal];
  my_blocks.NewAthenaArray(nblocal);
  // load MeshBlocks (parallel)
  if (resfile.Read_at_all(mbdata, datasize, nblocal, headeroffset+gids_*datasize) !=
      static_cast<unsigned int>(nblocal)) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "The restart file is broken or input parameters are inconsistent."
        << std::endl;
    ATHENA_ERROR(msg);
  }
  for (int i=gids_; i<=gide_; i++) {
    // Match fixed-width integer precision of IOWrapperSizeT datasize
    std::uint64_t buff_os = datasize * (i-gids_);
    SetBlockSizeAndBoundaries(loclist[i], block_size, block_bcs);
    my_blocks(i-gids_) = new MeshBlock(i, i-gids_, this, pin, loclist[i], block_size,
                                       block_bcs, costlist[i], mbdata+buff_os, gflag);
    my_blocks(i-gids_)->pbval->SearchAndSetNeighbors(tree, ranklist, nslist);
  }
  delete [] mbdata;
  // check consistency
  if (datasize != my_blocks(0)->GetBlockSizeInBytes()) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "The restart file is broken or input parameters are inconsistent."
        << std::endl;
    ATHENA_ERROR(msg);
  }

  ResetLoadBalanceVariables();

  // clean up
  delete [] offset;

  if (turb_flag > 0) // TurbulenceDriver depends on the MeshBlock ctor
    ptrbd = new TurbulenceDriver(this, pin);
}

//----------------------------------------------------------------------------------------
//! destructor

Mesh::~Mesh() {
  for (int b=0; b<nblocal; ++b)
    delete my_blocks(b);
  delete [] nslist;
  delete [] nblist;
  delete [] ranklist;
  delete [] costlist;
  delete [] loclist;
  if (SELF_GRAVITY_ENABLED == 1) delete pfgrd;
  else if (SELF_GRAVITY_ENABLED == 2) delete pmgrd;
  if (turb_flag > 0) delete ptrbd;
  if (adaptive) { // deallocate arrays for AMR
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
  if (nreal_user_mesh_data_>0) delete [] ruser_mesh_data;
  if (nuser_history_output_ > 0) {
    delete [] user_history_output_names_;
    delete [] user_history_func_;
    delete [] user_history_ops_;
  }
  if (nint_user_mesh_data_>0) delete [] iuser_mesh_data;
  if (EOS_TABLE_ENABLED) delete peos_table;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::OutputMeshStructure(int ndim)
//! \brief print the mesh structure information

void Mesh::OutputMeshStructure(int ndim) {
  RegionSize block_size;
  BoundaryFlag block_bcs[6];
  FILE *fp = nullptr;

  // open 'mesh_structure.dat' file
  if (f2) {
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
  int *nb_per_plevel = new int[max_level+1];
  int *cost_per_plevel = new int[max_level+1];
  for (int i=0; i<=max_level; ++i) {
    nb_per_plevel[i] = 0;
    cost_per_plevel[i] = 0;
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
    nb_per_rank[i] = 0;
    cost_per_rank[i] = 0;
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
  double real_max = std::numeric_limits<double>::max();
  double mincost = real_max, maxcost = 0.0, totalcost = 0.0;
  for (int i=root_level; i<=max_level; i++) {
    for (int j=0; j<nbtotal; j++) {
      if (loclist[j].level == i) {
        SetBlockSizeAndBoundaries(loclist[j], block_size, block_bcs);
        std::int64_t &lx1 = loclist[j].lx1;
        std::int64_t &lx2 = loclist[j].lx2;
        std::int64_t &lx3 = loclist[j].lx3;
        int &ll = loclist[j].level;
        mincost = std::min(mincost,costlist[i]);
        maxcost = std::max(maxcost,costlist[i]);
        totalcost += costlist[i];
        std::fprintf(fp,"#MeshBlock %d on rank=%d with cost=%g\n", j, ranklist[j],
                     costlist[j]);
        std::fprintf(
            fp, "#  Logical level %d, location = (%" PRId64 " %" PRId64 " %" PRId64")\n",
            ll, lx1, lx2, lx3);
        if (ndim == 2) {
          std::fprintf(fp, "%g %g\n", block_size.x1min, block_size.x2min);
          std::fprintf(fp, "%g %g\n", block_size.x1max, block_size.x2min);
          std::fprintf(fp, "%g %g\n", block_size.x1max, block_size.x2max);
          std::fprintf(fp, "%g %g\n", block_size.x1min, block_size.x2max);
          std::fprintf(fp, "%g %g\n", block_size.x1min, block_size.x2min);
          std::fprintf(fp, "\n\n");
        }
        if (ndim == 3) {
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
  if (f2) std::fclose(fp);
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
//! \fn void Mesh::NewTimeStep()
//! \brief function that loops over all MeshBlocks and find new timestep
//!        this assumes that phydro->NewBlockTimeStep is already called

void Mesh::NewTimeStep() {
  MeshBlock *pmb = my_blocks(0);

  // prevent timestep from growing too fast in between 2x cycles (even if every MeshBlock
  // has new_block_dt > 2.0*dt_old)
  dt = static_cast<Real>(2.0)*dt;
  // consider first MeshBlock on this MPI rank's linked list of blocks:
  dt = std::min(dt, pmb->new_block_dt_);
  dt_hyperbolic = pmb->new_block_dt_hyperbolic_;
  dt_parabolic = pmb->new_block_dt_parabolic_;
  dt_user = pmb->new_block_dt_user_;

  for (int i=0; i<nblocal; ++i) {
    pmb = my_blocks(i);
    dt = std::min(dt, pmb->new_block_dt_);
    dt_hyperbolic  = std::min(dt_hyperbolic, pmb->new_block_dt_hyperbolic_);
    dt_parabolic  = std::min(dt_parabolic, pmb->new_block_dt_parabolic_);
    dt_user  = std::min(dt_user, pmb->new_block_dt_user_);
  }

#ifdef MPI_PARALLEL
  // pack array, MPI allreduce over array, then unpack into Mesh variables
  Real dt_array[4] = {dt, dt_hyperbolic, dt_parabolic, dt_user};
  MPI_Allreduce(MPI_IN_PLACE, dt_array, 4, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
  dt            = dt_array[0];
  dt_hyperbolic = dt_array[1];
  dt_parabolic  = dt_array[2];
  dt_user       = dt_array[3];
#endif

  if (time < tlim && (tlim - time) < dt) // timestep would take us past desired endpoint
    dt = tlim - time;

  if (STS_ENABLED) {
    Real dt_ratio = dt / dt_parabolic;
    if (sts_max_dt_ratio > 0 && dt_ratio > sts_max_dt_ratio) {
      dt = sts_max_dt_ratio * dt_parabolic;
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollUserBoundaryFunction(BoundaryFace dir, BValFunc my_bc)
//! \brief Enroll a user-defined boundary function

void Mesh::EnrollUserBoundaryFunction(BoundaryFace dir, BValFunc my_bc) {
  std::stringstream msg;
  if (dir < 0 || dir > 5) {
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
//! \fn void Mesh::EnrollUserMGGravityBoundaryFunction(BoundaryFace dir,
//!                                                    MGBoundaryFunc my_bc)
//! \brief Enroll a user-defined Multigrid boundary function

void Mesh::EnrollUserMGGravityBoundaryFunction(BoundaryFace dir, MGBoundaryFunc my_bc) {
  std::stringstream msg;
  if (dir < 0 || dir > 5) {
    msg << "### FATAL ERROR in EnrollBoundaryCondition function" << std::endl
        << "dirName = " << dir << " not valid" << std::endl;
    ATHENA_ERROR(msg);
  }
  MGGravityBoundaryFunction_[static_cast<int>(dir)] = my_bc;
  return;
}

//! \deprecated (felker):
//! * provide trivial overloads for old-style BoundaryFace enum argument
void Mesh::EnrollUserBoundaryFunction(int dir, BValFunc my_bc) {
  EnrollUserBoundaryFunction(static_cast<BoundaryFace>(dir), my_bc);
  return;
}

void Mesh::EnrollUserMGGravityBoundaryFunction(int dir, MGBoundaryFunc my_bc) {
  EnrollUserMGGravityBoundaryFunction(static_cast<BoundaryFace>(dir), my_bc);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollUserRefinementCondition(AMRFlagFunc amrflag)
//! \brief Enroll a user-defined function for checking refinement criteria

void Mesh::EnrollUserRefinementCondition(AMRFlagFunc amrflag) {
  if (adaptive)
    AMRFlag_ = amrflag;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollUserMeshGenerator(CoordinateDirection,MeshGenFunc my_mg)
//! \brief Enroll a user-defined function for Mesh generation

void Mesh::EnrollUserMeshGenerator(CoordinateDirection dir, MeshGenFunc my_mg) {
  std::stringstream msg;
  if (dir < 0 || dir >= 3) {
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
  use_uniform_meshgen_fn_[dir] = false;
  MeshGenerator_[dir] = my_mg;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollUserExplicitSourceFunction(SrcTermFunc my_func)
//! \brief Enroll a user-defined source function

void Mesh::EnrollUserExplicitSourceFunction(SrcTermFunc my_func) {
  UserSourceTerm_ = my_func;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollUserTimeStepFunction(TimeStepFunc my_func)
//! \brief Enroll a user-defined time step function

void Mesh::EnrollUserTimeStepFunction(TimeStepFunc my_func) {
  UserTimeStep_ = my_func;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::AllocateUserHistoryOutput(int n)
//! \brief set the number of user-defined history outputs

void Mesh::AllocateUserHistoryOutput(int n) {
  nuser_history_output_ = n;
  user_history_output_names_ = new std::string[n];
  user_history_func_ = new HistoryOutputFunc[n];
  user_history_ops_ = new UserHistoryOperation[n];
  for (int i=0; i<n; i++) user_history_func_[i] = nullptr;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollUserHistoryOutput(int i, HistoryOutputFunc my_func,
//!                                        const char *name, UserHistoryOperation op)
//! \brief Enroll a user-defined history output function and set its name

void Mesh::EnrollUserHistoryOutput(int i, HistoryOutputFunc my_func, const char *name,
                                   UserHistoryOperation op) {
  std::stringstream msg;
  if (i >= nuser_history_output_) {
    msg << "### FATAL ERROR in EnrollUserHistoryOutput function" << std::endl
        << "The number of the user-defined history output (" << i << ") "
        << "exceeds the declared number (" << nuser_history_output_ << ")." << std::endl;
    ATHENA_ERROR(msg);
  }
  user_history_output_names_[i] = name;
  user_history_func_[i] = my_func;
  user_history_ops_[i] = op;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollUserMetric(MetricFunc my_func)
//! \brief Enroll a user-defined metric for arbitrary GR coordinates

void Mesh::EnrollUserMetric(MetricFunc my_func) {
  UserMetric_ = my_func;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollViscosityCoefficient(ViscosityCoeff my_func)
//! \brief Enroll a user-defined magnetic field diffusivity function

void Mesh::EnrollViscosityCoefficient(ViscosityCoeffFunc my_func) {
  ViscosityCoeff_ = my_func;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollConductionCoefficient(ConductionCoeff my_func)
//! \brief Enroll a user-defined thermal conduction function

void Mesh::EnrollConductionCoefficient(ConductionCoeffFunc my_func) {
  ConductionCoeff_ = my_func;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollFieldDiffusivity(FieldDiffusionCoeff my_func)
//! \brief Enroll a user-defined magnetic field diffusivity function

void Mesh::EnrollFieldDiffusivity(FieldDiffusionCoeffFunc my_func) {
  FieldDiffusivity_ = my_func;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollOrbitalVelocity(OrbitalVelocityFunc my_func)
//  \brief Enroll a user-defined orbital velocity function

void Mesh::EnrollOrbitalVelocity(OrbitalVelocityFunc my_func) {
  OrbitalVelocity_ = my_func;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollOrbitalVelocityDerivative(int i, OrbitalVelocityFunc my_func)
//  \brief Enroll Derivative fuctions of user-defined orbital velocity.

void Mesh::EnrollOrbitalVelocityDerivative(int i, OrbitalVelocityFunc my_func) {
  OrbitalVelocityDerivative_[i] = my_func;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::AllocateRealUserMeshDataField(int n)
//! \brief Allocate Real AthenaArrays for user-defned data in Mesh

void Mesh::AllocateRealUserMeshDataField(int n) {
  if (nreal_user_mesh_data_ != 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Mesh::AllocateRealUserMeshDataField"
        << std::endl << "User Mesh data arrays are already allocated" << std::endl;
    ATHENA_ERROR(msg);
  }
  nreal_user_mesh_data_ = n;
  ruser_mesh_data = new AthenaArray<Real>[n];
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::AllocateIntUserMeshDataField(int n)
//! \brief Allocate integer AthenaArrays for user-defned data in Mesh

void Mesh::AllocateIntUserMeshDataField(int n) {
  if (nint_user_mesh_data_ != 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Mesh::AllocateIntUserMeshDataField"
        << std::endl << "User Mesh data arrays are already allocated" << std::endl;
    ATHENA_ERROR(msg);
  }
  nint_user_mesh_data_ = n;
  iuser_mesh_data = new AthenaArray<int>[n];
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Mesh::ApplyUserWorkBeforeOutput(ParameterInput *pin)
//! \brief Apply MeshBlock::UserWorkBeforeOutput

void Mesh::ApplyUserWorkBeforeOutput(ParameterInput *pin) {
  for (int i=0; i<nblocal; ++i)
    my_blocks(i)->UserWorkBeforeOutput(pin);
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::Initialize(int res_flag, ParameterInput *pin)
//! \brief  initialization before the main loop

void Mesh::Initialize(int res_flag, ParameterInput *pin) {
  bool iflag = true;
  int inb = nbtotal;
  int nthreads = GetNumMeshThreads();

  do {
    if (res_flag == 0) {
#pragma omp parallel for num_threads(nthreads)
      for (int i=0; i<nblocal; ++i) {
        MeshBlock *pmb = my_blocks(i);
        pmb->ProblemGenerator(pin);
        pmb->pbval->CheckUserBoundaries();
      }
    }

    // add initial perturbation for decaying or impulsive turbulence
    if (((turb_flag == 1) || (turb_flag == 2)) && (res_flag == 0))
      ptrbd->Driving();

    // Create send/recv MPI_Requests for all BoundaryData objects
#pragma omp parallel for num_threads(nthreads)
    for (int i=0; i<nblocal; ++i) {
      MeshBlock *pmb = my_blocks(i);
      // BoundaryVariable objects evolved in main TimeIntegratorTaskList:
      pmb->pbval->SetupPersistentMPI();
      // other BoundaryVariable objects:
      if (SELF_GRAVITY_ENABLED == 1)
        pmb->pgrav->gbvar.SetupPersistentMPI();
    }

    // solve gravity for the first time
    if (SELF_GRAVITY_ENABLED == 1)
      pfgrd->Solve(1, 0);
    else if (SELF_GRAVITY_ENABLED == 2)
      pmgrd->Solve(1);

#pragma omp parallel num_threads(nthreads)
    {
      MeshBlock *pmb;
      BoundaryValues *pbval;

      // prepare to receive conserved variables
#pragma omp for private(pmb,pbval)
      for (int i=0; i<nblocal; ++i) {
        pmb = my_blocks(i); pbval = pmb->pbval;
        if (shear_periodic) {
          pbval->ComputeShear(time, time);
        }
        pbval->StartReceivingSubset(BoundaryCommSubset::mesh_init,
                                    pbval->bvars_main_int);
      }

      // send conserved variables
#pragma omp for private(pmb,pbval)
      for (int i=0; i<nblocal; ++i) {
        pmb = my_blocks(i); pbval = pmb->pbval;
        pmb->phydro->hbvar.SwapHydroQuantity(pmb->phydro->u,
                                               HydroBoundaryQuantity::cons);
        pmb->phydro->hbvar.SendBoundaryBuffers();
        if (MAGNETIC_FIELDS_ENABLED)
          pmb->pfield->fbvar.SendBoundaryBuffers();
        // and (conserved variable) passive scalar masses:
        if (NSCALARS > 0) {
          pmb->pscalars->sbvar.var_cc = &(pmb->pscalars->s);
          if (pmb->pmy_mesh->multilevel) {
            pmb->pscalars->sbvar.coarse_buf = &(pmb->pscalars->coarse_s_);
          }
          pmb->pscalars->sbvar.SendBoundaryBuffers();
        }
      }

      // wait to receive conserved variables
#pragma omp for private(pmb,pbval)
      for (int i=0; i<nblocal; ++i) {
        pmb = my_blocks(i); pbval = pmb->pbval;
        pmb->phydro->hbvar.ReceiveAndSetBoundariesWithWait();
        if (MAGNETIC_FIELDS_ENABLED)
          pmb->pfield->fbvar.ReceiveAndSetBoundariesWithWait();
        if (NSCALARS > 0)
          pmb->pscalars->sbvar.ReceiveAndSetBoundariesWithWait();
        if (shear_periodic && orbital_advection==0) {
          pmb->phydro->hbvar.AddHydroShearForInit();
        }
        pbval->ClearBoundarySubset(BoundaryCommSubset::mesh_init,
                                   pbval->bvars_main_int);
      }

      // With AMR/SMR GR send primitives to enable cons->prim before prolongation
      if (GENERAL_RELATIVITY && multilevel) {
        // prepare to receive primitives
#pragma omp for private(pmb,pbval)
        for (int i=0; i<nblocal; ++i) {
          pmb = my_blocks(i); pbval = pmb->pbval;
          pbval->StartReceivingSubset(BoundaryCommSubset::gr_amr,
                                      pbval->bvars_main_int);
        }

        // send primitives
#pragma omp for private(pmb,pbval)
        for (int i=0; i<nblocal; ++i) {
          pmb = my_blocks(i); pbval = pmb->pbval;
          pmb->phydro->hbvar.SwapHydroQuantity(pmb->phydro->w,
                                               HydroBoundaryQuantity::prim);
          pmb->phydro->hbvar.SendBoundaryBuffers();
          if (NSCALARS > 0) {
            pmb->pscalars->sbvar.var_cc = &(pmb->pscalars->r);
            if (pmb->pmy_mesh->multilevel) {
              pmb->pscalars->sbvar.coarse_buf = &(pmb->pscalars->coarse_r_);
            }
            pmb->pscalars->sbvar.SendBoundaryBuffers();
          }
        }

        // wait to receive AMR/SMR GR primitives
#pragma omp for private(pmb,pbval)
        for (int i=0; i<nblocal; ++i) {
          pmb = my_blocks(i); pbval = pmb->pbval;
          pmb->phydro->hbvar.ReceiveAndSetBoundariesWithWait();
          if (NSCALARS > 0) {
            pmb->pscalars->sbvar.ReceiveAndSetBoundariesWithWait();
          }
          pbval->ClearBoundarySubset(BoundaryCommSubset::gr_amr,
                                     pbval->bvars_main_int);
          pmb->phydro->hbvar.SwapHydroQuantity(pmb->phydro->u,
                                               HydroBoundaryQuantity::cons);
          if (NSCALARS > 0) {
            pmb->pscalars->sbvar.var_cc = &(pmb->pscalars->s);
            if (pmb->pmy_mesh->multilevel) {
              pmb->pscalars->sbvar.coarse_buf = &(pmb->pscalars->coarse_s_);
            }
          }
        }
      } // multilevel

      // perform fourth-order correction of midpoint initial condition:
      // (correct IC on all MeshBlocks or none; switch cannot be toggled independently)
      bool correct_ic = my_blocks(0)->precon->correct_ic;
      if (correct_ic)
        CorrectMidpointInitialCondition();

      // Initiate orbital advection
#pragma omp for private(pmb)
      for (int i=0; i<nblocal; ++i) {
        pmb = my_blocks(i);
        if (pmb->porb->orbital_advection_defined) {
          pmb->porb->InitializeOrbitalAdvection();
          if (pmb->porb->orbital_advection_active)
            pmb->porb->orb_bc->SetupPersistentMPI();
        }
      }

      // Now do prolongation, compute primitives, apply BCs
      Hydro *ph;
      Field *pf;
      PassiveScalars *ps;
#pragma omp for private(pmb,pbval,ph,pf,ps)
      for (int i=0; i<nblocal; ++i) {
        pmb = my_blocks(i);
        pbval = pmb->pbval, ph = pmb->phydro, pf = pmb->pfield, ps = pmb->pscalars;
        if (multilevel)
          pbval->ProlongateBoundaries(time, 0.0, pbval->bvars_main_int);

        int il = pmb->is, iu = pmb->ie,
            jl = pmb->js, ju = pmb->je,
            kl = pmb->ks, ku = pmb->ke;
        if (pbval->nblevel[1][1][0] != -1) il -= NGHOST;
        if (pbval->nblevel[1][1][2] != -1) iu += NGHOST;
        if (pmb->block_size.nx2 > 1) {
          if (pbval->nblevel[1][0][1] != -1) jl -= NGHOST;
          if (pbval->nblevel[1][2][1] != -1) ju += NGHOST;
        }
        if (pmb->block_size.nx3 > 1) {
          if (pbval->nblevel[0][1][1] != -1) kl -= NGHOST;
          if (pbval->nblevel[2][1][1] != -1) ku += NGHOST;
        }
        pmb->peos->ConservedToPrimitive(ph->u, ph->w1, pf->b,
                                        ph->w, pf->bcc, pmb->pcoord,
                                        il, iu, jl, ju, kl, ku);
        if (NSCALARS > 0) {
          // r1/r_old for GR is currently unused:
          pmb->peos->PassiveScalarConservedToPrimitive(ps->s, ph->u, ps->r, ps->r,
                                                       pmb->pcoord,
                                                       il, iu, jl, ju, kl, ku);
        }
        // --------------------------
        int order = pmb->precon->xorder;
        if (order == 4) {
          // fourth-order EOS:
          // for hydro, shrink buffer by 1 on all sides
          if (pbval->nblevel[1][1][0] != -1) il += 1;
          if (pbval->nblevel[1][1][2] != -1) iu -= 1;
          if (pbval->nblevel[1][0][1] != -1) jl += 1;
          if (pbval->nblevel[1][2][1] != -1) ju -= 1;
          if (pbval->nblevel[0][1][1] != -1) kl += 1;
          if (pbval->nblevel[2][1][1] != -1) ku -= 1;
          // for MHD, shrink buffer by 3
          //! \todo (felker):
          //! * add MHD loop limit calculation for 4th order W(U)
          // Apply physical boundaries prior to 4th order W(U)
          ph->hbvar.SwapHydroQuantity(ph->w, HydroBoundaryQuantity::prim);
          if (NSCALARS > 0) {
            ps->sbvar.var_cc = &(ps->r);
            if (pmb->pmy_mesh->multilevel) {
              ps->sbvar.coarse_buf = &(ps->coarse_r_);
            }
          }
          pbval->ApplyPhysicalBoundaries(time, 0.0, pbval->bvars_main_int);
          // Perform 4th order W(U)
          pmb->peos->ConservedToPrimitiveCellAverage(ph->u, ph->w1, pf->b,
                                                     ph->w, pf->bcc, pmb->pcoord,
                                                     il, iu, jl, ju, kl, ku);
          if (NSCALARS > 0) {
            pmb->peos->PassiveScalarConservedToPrimitiveCellAverage(
                ps->s, ps->r, ps->r, pmb->pcoord, il, iu, jl, ju, kl, ku);
          }
        }
        // --------------------------
        // end fourth-order EOS

        // Swap Hydro and (possibly) passive scalar quantities in BoundaryVariable
        // interface from conserved to primitive formulations:
        ph->hbvar.SwapHydroQuantity(ph->w, HydroBoundaryQuantity::prim);
        if (NSCALARS > 0) {
          ps->sbvar.var_cc = &(ps->r);
          if (pmb->pmy_mesh->multilevel) {
            ps->sbvar.coarse_buf = &(ps->coarse_r_);
          }
        }

        pbval->ApplyPhysicalBoundaries(time, 0.0, pbval->bvars_main_int);
      }

      // Calc initial diffusion coefficients
#pragma omp for private(pmb,ph,pf)
      for (int i=0; i<nblocal; ++i) {
        pmb = my_blocks(i); ph = pmb->phydro, pf = pmb->pfield;
        if (ph->hdif.hydro_diffusion_defined)
          ph->hdif.SetDiffusivity(ph->w, pf->bcc);
        if (MAGNETIC_FIELDS_ENABLED) {
          if (pf->fdif.field_diffusion_defined)
            pf->fdif.SetDiffusivity(ph->w, pf->bcc);
        }
      }

      if (!res_flag && adaptive) {
#pragma omp for
        for (int i=0; i<nblocal; ++i) {
          my_blocks(i)->pmr->CheckRefinementCondition();
        }
      }
    } // omp parallel

    if (!res_flag && adaptive) {
      iflag = false;
      int onb = nbtotal;
      LoadBalancingAndAdaptiveMeshRefinement(pin);
      amr_updated = true;
      if (nbtotal == onb) {
        iflag = true;
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
  } while (!iflag);

  // calculate the first time step
#pragma omp parallel for num_threads(nthreads)
  for (int i=0; i<nblocal; ++i) {
    my_blocks(i)->phydro->NewBlockTimeStep();
  }

  NewTimeStep();
  return;
}

//----------------------------------------------------------------------------------------
//! \fn MeshBlock* Mesh::FindMeshBlock(int tgid)
//! \brief return the MeshBlock whose gid is tgid

MeshBlock* Mesh::FindMeshBlock(int tgid) {
  if (tgid >= gids_ && tgid <= gide_)
    return my_blocks(tgid - gids_);
  return nullptr;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::SetBlockSizeAndBoundaries(LogicalLocation loc,
//!                RegionSize &block_size, BundaryFlag *block_bcs)
//! \brief Set the physical part of a block_size structure and block boundary conditions

void Mesh::SetBlockSizeAndBoundaries(LogicalLocation loc, RegionSize &block_size,
                                     BoundaryFlag *block_bcs) {
  std::int64_t &lx1 = loc.lx1;
  int &ll = loc.level;
  std::int64_t nrbx_ll = nrbx1 << (ll - root_level);

  // calculate physical block size, x1
  if (lx1 == 0) {
    block_size.x1min = mesh_size.x1min;
    block_bcs[BoundaryFace::inner_x1] = mesh_bcs[BoundaryFace::inner_x1];
  } else {
    Real rx = ComputeMeshGeneratorX(lx1, nrbx_ll, use_uniform_meshgen_fn_[X1DIR]);
    block_size.x1min = MeshGenerator_[X1DIR](rx, mesh_size);
    block_bcs[BoundaryFace::inner_x1] = BoundaryFlag::block;
  }
  if (lx1 == nrbx_ll - 1) {
    block_size.x1max = mesh_size.x1max;
    block_bcs[BoundaryFace::outer_x1] = mesh_bcs[BoundaryFace::outer_x1];
  } else {
    Real rx = ComputeMeshGeneratorX(lx1+1, nrbx_ll, use_uniform_meshgen_fn_[X1DIR]);
    block_size.x1max = MeshGenerator_[X1DIR](rx, mesh_size);
    block_bcs[BoundaryFace::outer_x1] = BoundaryFlag::block;
  }

  // calculate physical block size, x2
  if (mesh_size.nx2 == 1) {
    block_size.x2min = mesh_size.x2min;
    block_size.x2max = mesh_size.x2max;
    block_bcs[BoundaryFace::inner_x2] = mesh_bcs[BoundaryFace::inner_x2];
    block_bcs[BoundaryFace::outer_x2] = mesh_bcs[BoundaryFace::outer_x2];
  } else {
    std::int64_t &lx2 = loc.lx2;
    nrbx_ll = nrbx2 << (ll - root_level);
    if (lx2 == 0) {
      block_size.x2min = mesh_size.x2min;
      block_bcs[BoundaryFace::inner_x2] = mesh_bcs[BoundaryFace::inner_x2];
    } else {
      Real rx = ComputeMeshGeneratorX(lx2, nrbx_ll, use_uniform_meshgen_fn_[X2DIR]);
      block_size.x2min = MeshGenerator_[X2DIR](rx, mesh_size);
      block_bcs[BoundaryFace::inner_x2] = BoundaryFlag::block;
    }
    if (lx2 == (nrbx_ll) - 1) {
      block_size.x2max = mesh_size.x2max;
      block_bcs[BoundaryFace::outer_x2] = mesh_bcs[BoundaryFace::outer_x2];
    } else {
      Real rx = ComputeMeshGeneratorX(lx2+1, nrbx_ll, use_uniform_meshgen_fn_[X2DIR]);
      block_size.x2max = MeshGenerator_[X2DIR](rx, mesh_size);
      block_bcs[BoundaryFace::outer_x2] = BoundaryFlag::block;
    }
  }

  // calculate physical block size, x3
  if (mesh_size.nx3 == 1) {
    block_size.x3min = mesh_size.x3min;
    block_size.x3max = mesh_size.x3max;
    block_bcs[BoundaryFace::inner_x3] = mesh_bcs[BoundaryFace::inner_x3];
    block_bcs[BoundaryFace::outer_x3] = mesh_bcs[BoundaryFace::outer_x3];
  } else {
    std::int64_t &lx3 = loc.lx3;
    nrbx_ll = nrbx3 << (ll - root_level);
    if (lx3 == 0) {
      block_size.x3min = mesh_size.x3min;
      block_bcs[BoundaryFace::inner_x3] = mesh_bcs[BoundaryFace::inner_x3];
    } else {
      Real rx = ComputeMeshGeneratorX(lx3, nrbx_ll, use_uniform_meshgen_fn_[X3DIR]);
      block_size.x3min = MeshGenerator_[X3DIR](rx, mesh_size);
      block_bcs[BoundaryFace::inner_x3] = BoundaryFlag::block;
    }
    if (lx3 == (nrbx_ll) - 1) {
      block_size.x3max = mesh_size.x3max;
      block_bcs[BoundaryFace::outer_x3] = mesh_bcs[BoundaryFace::outer_x3];
    } else {
      Real rx = ComputeMeshGeneratorX(lx3+1, nrbx_ll, use_uniform_meshgen_fn_[X3DIR]);
      block_size.x3max = MeshGenerator_[X3DIR](rx, mesh_size);
      block_bcs[BoundaryFace::outer_x3] = BoundaryFlag::block;
    }
  }

  block_size.x1rat = mesh_size.x1rat;
  block_size.x2rat = mesh_size.x2rat;
  block_size.x3rat = mesh_size.x3rat;

  return;
}


void Mesh::CorrectMidpointInitialCondition() {
  MeshBlock *pmb;
  Hydro *ph;
  PassiveScalars *ps;
#pragma omp for private(pmb, ph, ps)
  for (int nb=0; nb<nblocal; ++nb) {
    pmb = my_blocks(nb);
    ph = pmb->phydro;
    ps = pmb->pscalars;

    // Assume cell-centered analytic value is computed at all real cells, and ghost
    // cells with the cell-centered U have been exchanged
    int il = pmb->is, iu = pmb->ie, jl = pmb->js, ju = pmb->je,
        kl = pmb->ks, ku = pmb->ke;

    // Laplacian of cell-averaged conserved variables, scalar concentrations
    AthenaArray<Real> delta_cons_, delta_s_;

    // Allocate memory for 4D Laplacian
    int ncells4 = NHYDRO;
    int nl = 0;
    int nu = ncells4 - 1;
    delta_cons_.NewAthenaArray(ncells4, pmb->ncells3, pmb->ncells2, pmb->ncells1);

    // Compute and store Laplacian of cell-averaged conserved variables
    pmb->pcoord->Laplacian(ph->u, delta_cons_, il, iu, jl, ju, kl, ku, nl, nu);

    //! \todo (felker):
    //! * assuming uniform mesh with dx1f=dx2f=dx3f, so this factors out
    //! * also, this may need to be dx1v, since Laplacian is cell-center
    Real h = pmb->pcoord->dx1f(il);  // pco->dx1f(i); inside loop
    Real C = (h*h)/24.0;

    // Compute fourth-order approximation to cell-centered conserved variables
    for (int n=nl; n<=nu; ++n) {
      for (int k=kl; k<=ku; ++k) {
        for (int j=jl; j<=ju; ++j) {
          for (int i=il; i<=iu; ++i) {
            // We do not actually need to store all cell-centered cons. variables,
            // but the ConservedToPrimitivePointwise() implementation operates on 4D
            ph->u(n,k,j,i) = ph->u(n,k,j,i) + C*delta_cons_(n,k,j,i);
          }
        }
      }
    }

    // If NSCALARS < NHYDRO, could reuse delta_cons_ allocated memory...
    int ncells4_s = NSCALARS;
    int nl_s = 0;
    int nu_s = ncells4_s -1;
    if (NSCALARS > 0) {
      delta_s_.NewAthenaArray(ncells4_s, pmb->ncells3, pmb->ncells2, pmb->ncells1);
      pmb->pcoord->Laplacian(ps->s, delta_s_, il, iu, jl, ju, kl, ku, nl_s, nu_s);
    }

    // Compute fourth-order approximation to cell-centered conserved variables
    for (int n=nl_s; n<=nu_s; ++n) {
      for (int k=kl; k<=ku; ++k) {
        for (int j=jl; j<=ju; ++j) {
          for (int i=il; i<=iu; ++i) {
            ps->s(n,k,j,i) = ps->s(n,k,j,i) + C*delta_s_(n,k,j,i);
          }
        }
      }
    }
  } // end loop over MeshBlocks

  // begin second exchange of ghost cells with corrected cell-averaged <U>
  // -----------------  (mostly copied from above section in Mesh::Initialize())
  BoundaryValues *pbval;
  // prepare to receive conserved variables
#pragma omp for private(pmb,pbval)
  for (int i=0; i<nblocal; ++i) {
    pmb = my_blocks(i); pbval = pmb->pbval;
    if (shear_periodic) {
      pbval->ComputeShear(time, time);
    }
    // no need to re-SetupPersistentMPI() the MPI requests for boundary values
    pbval->StartReceivingSubset(BoundaryCommSubset::mesh_init,
                                pbval->bvars_main_int);
  }

#pragma omp for private(pmb,pbval)
  for (int i=0; i<nblocal; ++i) {
    pmb = my_blocks(i);
    pmb->phydro->hbvar.SwapHydroQuantity(pmb->phydro->u,
                                         HydroBoundaryQuantity::cons);
    pmb->phydro->hbvar.SendBoundaryBuffers();
    if (MAGNETIC_FIELDS_ENABLED)
      pmb->pfield->fbvar.SendBoundaryBuffers();
    // and (conserved variable) passive scalar masses:
    if (NSCALARS > 0) {
      pmb->pscalars->sbvar.var_cc = &(pmb->pscalars->s);
      if (pmb->pmy_mesh->multilevel) {
        pmb->pscalars->sbvar.coarse_buf = &(pmb->pscalars->coarse_s_);
      }
      pmb->pscalars->sbvar.SendBoundaryBuffers();
    }
  }

  // wait to receive conserved variables
#pragma omp for private(pmb,pbval)
  for (int i=0; i<nblocal; ++i) {
    pmb = my_blocks(i); pbval = pmb->pbval;
    pmb->phydro->hbvar.SwapHydroQuantity(pmb->phydro->u,
                                         HydroBoundaryQuantity::cons);
    pmb->phydro->hbvar.ReceiveAndSetBoundariesWithWait();
    if (MAGNETIC_FIELDS_ENABLED)
      pmb->pfield->fbvar.ReceiveAndSetBoundariesWithWait();
    if (NSCALARS > 0) {
      pmb->pscalars->sbvar.var_cc = &(pmb->pscalars->s);
      if (pmb->pmy_mesh->multilevel) {
        pmb->pscalars->sbvar.coarse_buf = &(pmb->pscalars->coarse_s_);
      }
      pmb->pscalars->sbvar.ReceiveAndSetBoundariesWithWait();
    }
    if (shear_periodic && orbital_advection==0) {
      pmb->phydro->hbvar.AddHydroShearForInit();
    }
    pbval->ClearBoundarySubset(BoundaryCommSubset::mesh_init,
                               pbval->bvars_main_int);
  } // end second exchange of ghost cells
  return;
}

//----------------------------------------------------------------------------------------
//! \fn int Mesh::ReserveTagPhysIDs(int num_phys)
//! \brief Public function for advancing next_phys_id_ counter
//!
//! E.g. if chemistry or radiation elects to communicate additional information with MPI
//! outside the framework of the BoundaryVariable classes
//!
//! Store signed, but positive, integer corresponding to the next unused value to be used
//! as unique ID for a BoundaryVariable object's single set of MPI calls (formerly "enum
//! AthenaTagMPI"). 5 bits of unsigned integer representation are currently reserved
//! for this "phys" part of the bitfield tag, making 0, ..., 31 legal values

int Mesh::ReserveTagPhysIDs(int num_phys) {
  //! \todo (felker):
  //! * add safety checks? input, output are positive, obey <= 31= MAX_NUM_PHYS
  int start_id = next_phys_id_;
  next_phys_id_ += num_phys;
  return start_id;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::ReserveMeshBlockPhysIDs()
//! \brief private member fn, called in Mesh() ctor
//!
//! depending on compile- and runtime options, reserve the maximum number of "int physid"
//! that might be necessary for each MeshBlock's BoundaryValues object to perform MPI
//! communication for all BoundaryVariable objects
//!
//! \todo (felker):
//! * deduplicate this logic, which combines conditionals in MeshBlock ctor

void Mesh::ReserveMeshBlockPhysIDs() {
#ifdef MPI_PARALLEL
  // if (FLUID_ENABLED) {
  // Advance Mesh's shared counter (initialized to next_phys_id=1 if MPI)
  // Greedy reservation of phys IDs (only 1 of 2 needed for Hydro if multilevel==false)
  ReserveTagPhysIDs(HydroBoundaryVariable::max_phys_id);
  //  }
  if (MAGNETIC_FIELDS_ENABLED) {
    ReserveTagPhysIDs(FaceCenteredBoundaryVariable::max_phys_id);
  }
  if (SELF_GRAVITY_ENABLED) {
    ReserveTagPhysIDs(CellCenteredBoundaryVariable::max_phys_id);
  }
  if (NSCALARS > 0) {
    ReserveTagPhysIDs(CellCenteredBoundaryVariable::max_phys_id);
  }
#endif
  return;
}

//----------------------------------------------------------------------------------------
//! \fn FluidFormulation GetFluidFormulation(const std::string& input_string)
//! \brief Parses input string to return scoped enumerator flag specifying boundary
//! condition. Typically called in Mesh() ctor initializer list

FluidFormulation GetFluidFormulation(const std::string& input_string) {
  if (input_string == "true") {
    return FluidFormulation::evolve;
  } else if (input_string == "disabled") {
    return FluidFormulation::disabled;
  } else if (input_string == "background") {
    return FluidFormulation::background;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in GetFluidFormulation" << std::endl
        << "Input string=" << input_string << "\n"
        << "is an invalid fluid formulation" << std::endl;
    ATHENA_ERROR(msg);
  }
}

void Mesh::OutputCycleDiagnostics() {
  const int dt_precision = std::numeric_limits<Real>::max_digits10 - 1;
  const int ratio_precision = 3;
  if (ncycle_out != 0) {
    if (ncycle % ncycle_out == 0) {
      std::cout << "cycle=" << ncycle << std::scientific
                << std::setprecision(dt_precision)
                << " time=" << time << " dt=" << dt;
      if (dt_diagnostics != -1) {
        if (STS_ENABLED) {
          if (UserTimeStep_ == nullptr)
            std::cout << "=dt_hyperbolic";
          // remaining dt_parabolic diagnostic output handled in STS StartupTaskList
        } else {
          Real ratio = dt / dt_hyperbolic;
          std::cout << "\ndt_hyperbolic=" << dt_hyperbolic << " ratio="
                    << std::setprecision(ratio_precision) << ratio
                    << std::setprecision(dt_precision);
          ratio = dt / dt_parabolic;
          std::cout << "\ndt_parabolic=" << dt_parabolic << " ratio="
                    << std::setprecision(ratio_precision) << ratio
                    << std::setprecision(dt_precision);
        }
        if (UserTimeStep_ != nullptr) {
          Real ratio = dt / dt_user;
          std::cout << "\ndt_user=" << dt_user << " ratio="
                    << std::setprecision(ratio_precision) << ratio
                    << std::setprecision(dt_precision);
        }
      } // else (empty): dt_diagnostics = -1 -> provide no additional timestep diagnostics
      std::cout << std::endl;
    }
  }
  return;
}
