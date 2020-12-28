//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file orbital_advection.cpp
//! \brief implementation of the OrbitalAdvection class

// C/C++ headers
#include <algorithm>  // max(), min()
#include <cfloat>     // FLT_MAX, FLT_MIN
#include <iostream>   // cout, endl
#include <sstream>    //
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"

#include "../bvals/bvals.hpp"
#include "../bvals/bvals_interfaces.hpp"
#include "../coordinates/coordinates.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../reconstruct/reconstruction.hpp"
#include "../scalars/scalars.hpp"

// this class header
#include "orbital_advection.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif


// constructor of OrbitalAdvection class
OrbitalAdvection::OrbitalAdvection(MeshBlock *pmb, ParameterInput *pin)
    : pmb_(pmb), pm_(pmb->pmy_mesh), ph_(pmb->phydro),
      pf_(pmb->pfield), pco_(pmb->pcoord), pbval_(pmb->pbval), ps_(pmb->pscalars) {
  // read parameters from input file
  orbital_splitting_order = pm_->orbital_advection;
  orbital_advection_defined = (orbital_splitting_order != 0) ? true : false;
  orbital_advection_active = orbital_advection_defined;

  // check xorder for reconstruction
  xorder = pmb_->precon->xorder;
  xgh = (xorder<=2)? 1 : 2;

  // Read parameters
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    Omega0  = pin->GetOrAddReal("orbital_advection","Omega0",0.0);
    qshear  = pin->GetOrAddReal("orbital_advection","qshear",0.0);
    shboxcoord = pin->GetOrAddInteger("orbital_advection","shboxcoord",1);
    onx = pmb_->block_size.nx2;
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    gm  = pin->GetOrAddReal("problem","GM",0.0);
    Omega0  = pin->GetOrAddReal("orbital_advection","Omega0",0.0);
    onx = pmb_->block_size.nx2;
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    gm  = pin->GetOrAddReal("problem","GM",0.0);
    Omega0  = pin->GetOrAddReal("orbital_advection","Omega0",0.0);
    onx = pmb_->block_size.nx3;
  }

  // set default orbital velocity functions
  if (pmb->pmy_mesh->OrbitalVelocity_ == nullptr) {
    // not user-defined functions
    if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
      OrbitalVelocity = CartOrbitalVelocity;
      OrbitalVelocityDerivative[0] = CartOrbitalVelocity_x;
      OrbitalVelocityDerivative[1] = ZeroOrbitalVelocity;
    } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
      if (pmb->block_size.nx3>1) { // 3D
        OrbitalVelocity = CylOrbitalVelocity3D;
        OrbitalVelocityDerivative[0] = CylOrbitalVelocity3D_r;
        OrbitalVelocityDerivative[1] = CylOrbitalVelocity3D_z;
      } else { // 2D
        OrbitalVelocity = CylOrbitalVelocity2D;
        OrbitalVelocityDerivative[0] = CylOrbitalVelocity2D_r;
        OrbitalVelocityDerivative[1] = ZeroOrbitalVelocity;
      }
    } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
      OrbitalVelocity = SphOrbitalVelocity;
      OrbitalVelocityDerivative[0] = SphOrbitalVelocity_r;
      OrbitalVelocityDerivative[1] = SphOrbitalVelocity_t;
    }
  } else { // user-defined orbital velocity functions
    OrbitalVelocity = pmb->pmy_mesh->OrbitalVelocity_;
    // if derivatives are not given, caluclate them automatically
    if (pmb->pmy_mesh->OrbitalVelocityDerivative_[0] != nullptr)
      OrbitalVelocityDerivative[0] = pmb->pmy_mesh->OrbitalVelocityDerivative_[0];
    if (pmb->pmy_mesh->OrbitalVelocityDerivative_[1] != nullptr)
      OrbitalVelocityDerivative[1] = pmb->pmy_mesh->OrbitalVelocityDerivative_[1];
  }

  // set meshblock size & orbital_direction
  // For orbital_direction ==1, x2 is the orbital direction (2D/3D)
  // For orbital_direction ==2, x3 is the orbital direction (3D)
  orbital_direction = 0;
  if(orbital_advection_defined) {
    nc1 = pmb_->block_size.nx1 + 2*NGHOST;
    nc2 = 1, nc3 = 1;
    if ((std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0)
         || (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0)) {
      orbital_direction = 1;
      if (pmb_->block_size.nx2 > 1) { // 2D or 3D
        nc2 = pmb_->block_size.nx2 + 2*(NGHOST);
      } else { // 1D
        orbital_advection_active = false;
      }
      if (pmb_->block_size.nx3 > 1) { // 3D
        nc3 = pmb_->block_size.nx3 + 2*(NGHOST);
      }
    } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
      orbital_direction = 2;
      if (pmb_->block_size.nx3 > 1) { // 3D
        nc2 = pmb_->block_size.nx2 + 2*(NGHOST);
        nc3 = pmb_->block_size.nx3 + 2*(NGHOST);
      } else if (pmb_->block_size.nx2 > 1) { // 2D
        orbital_advection_active = false;
        nc2 = pmb_->block_size.nx2 + 2*(NGHOST);
      } else { // 1D
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalAdvection Class" << std::endl
            << "Orbital advection requires 2D or 3D in spherical polar coordinates."
            << std::endl
            << "Check <orbital_advection> order parameter in the input file."
            << std::endl;
        ATHENA_ERROR(msg);
      }
    } else {
      std::stringstream msg;
      msg << "### FATAL ERROR in OrbitalAdvection Class" << std::endl
          << "Orbital advection works only in cartesian, cylindrical, "
          << "or spherical_polar coordinates." << std::endl
          << "Check <orbital_advection> order parameter in the input file."
          << std::endl;
      ATHENA_ERROR(msg);
    }

    // check parameters about the orbital motion for using default orbital velocity
    if (pmb->pmy_mesh->OrbitalVelocity_ == nullptr) {
      if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
        if ((Omega0 == 0.0) || (qshear == 0.0)) {
          std::stringstream msg;
          msg << "### FATAL ERROR in OrbitalAdvection Class" << std::endl
              << "The default orbital velocity profile requires non-zero "
              << "Omega0 and qshear." << std::endl
              << "Check <orbital_advection> Omega0 and qshear in the input file."
              << std::endl;
          ATHENA_ERROR(msg);
        }
      } else if ((std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0)
                || (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0)) {
        if (gm == 0.0) {
          std::stringstream msg;
          msg << "### FATAL ERROR in OrbitalAdvection Class" << std::endl
              << "The default orbital velocity profile requires non-zero GM." << std::endl
              << "Check <problem> GM and <orbital_advection> qshear in the input file."
              << std::endl;
          ATHENA_ERROR(msg);
        }
      }
    }

    // check orbital_refinement
    orbital_refinement = false;
    if (orbital_advection_active) {
      if (pm_->adaptive==true) { //AMR
        orbital_refinement = true;
      } else if (pm_->multilevel==true) { //SMR
        if (orbital_direction == 1) { // cartesian or cylindrical
          int64_t nbx = pm_->nrbx2 * (1L << (pmb_->loc.level - pm_->root_level));
          LogicalLocation loc;
          loc.level = pmb_->loc.level;
          loc.lx1   = pmb_->loc.lx1;
          loc.lx3   = pmb_->loc.lx3;
          for (int64_t dlx=1; dlx < nbx; dlx++) {
            loc.lx2 = pmb_->loc.lx2+dlx;
            if (loc.lx2>=nbx) loc.lx2-=nbx;
            //check level of meshblocks at same i, k
            MeshBlockTree *mbt = pm_->tree.FindMeshBlock(loc);
            if(mbt == nullptr || mbt->GetGid() == -1) {
              orbital_refinement = true;
              break;
            }
          }
        } else if (orbital_direction == 2) { // spherical_polar
          int64_t nbx = pm_->nrbx3 * (1L << (pmb_->loc.level - pm_->root_level));
          LogicalLocation loc;
          loc.level = pmb_->loc.level;
          loc.lx1   = pmb_->loc.lx1;
          loc.lx2   = pmb_->loc.lx2;
          for (int64_t dlx=0; dlx < nbx; dlx++) {
            loc.lx3   = pmb_->loc.lx3+dlx;
            if (loc.lx3>=nbx) loc.lx3-=nbx;
            //check level of meshblocks at same i, j
            MeshBlockTree *mbt = pm_->tree.FindMeshBlock(loc);
            if(mbt == nullptr || mbt->GetGid() == -1) {
              orbital_refinement = true;
              break;
            }
          }
        }
      } // SMR
    } // orbital_advection_active

    // check boundary conditions in the orbital direction
    if (orbital_direction == 1) { // cartesian or cylindrical
      if(pin->GetOrAddString("mesh", "ix2_bc", "none")!="periodic"
        ||pin->GetOrAddString("mesh", "ox2_bc", "none")!="periodic") {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalAdvection Class" << std::endl
            << "Orbital advection in Cartesian or cylindrical coordinates requires "
            << "periodic boundary conditions in the x2 direction." << std::endl
            << "Check <mesh> ix2_bc and ox2_bc in the input file." << std::endl;
        ATHENA_ERROR(msg);
      }
    } else if (orbital_direction==2) { // spherical_polar
      if(pin->GetOrAddString("mesh", "ix3_bc", "none")!="periodic"
          ||pin->GetOrAddString("mesh", "ox3_bc", "none")!="periodic") {
        std::stringstream msg;
        msg << "### FATAL ERROR in OrbitalAdvection Class" << std::endl
            << "Orbital advection in spherical_polar coordinates requires "
            << "periodic boundary conditions in the x3 direction." << std::endl
            << "Check <mesh> ix3_bc and ox3_bc in the input file." << std::endl;
        ATHENA_ERROR(msg);
      }
    }

    // initialize orbital_system_conversion_done flag
    orbital_system_conversion_done = 3;

    // check orbital_uniform_mesh
    // TODO(tomo-ono): non-uniform mesh grids are not allowed now.
    orbital_uniform_mesh = pm_->use_uniform_meshgen_fn_[orbital_direction];
    if (!orbital_uniform_mesh) {
      std::stringstream msg;
      msg << "### FATAL ERROR in OrbitalAdvection Class" << std::endl
          << "Orbital advection currently does not support non-uniform mesh"
          << "in the orbital direction." << std::endl;
      ATHENA_ERROR(msg);
    }

    // check orbital_splitting_order
    if (orbital_splitting_order<0 || orbital_splitting_order > 2) {
      std::stringstream msg;
      msg << "### FATAL ERROR in OrbitalAdvection Class" << std::endl
          << "The order of the orbital splitting must be from 0 to 2." << std::endl
          << "Check <orbital_advection> order parameter." << std::endl;
      ATHENA_ERROR(msg);
    }

    // check relativity
    if (RELATIVISTIC_DYNAMICS) {
      std::stringstream msg;
      msg << "### FATAL ERROR in OrbitalAdvection Class."<<std::endl
          << "Neigher SR nor GR supports Orbital Advection."<<std::endl;
      ATHENA_ERROR(msg);
    }

    // memory allocation
    if (orbital_direction == 1) { // cartesian or cylindrical
      vKc.NewAthenaArray(nc3, nc1);
      dvKc1.NewAthenaArray(nc3, nc1);
      dvKc2.NewAthenaArray(nc3, nc1);
      vKf[0].NewAthenaArray(nc3,nc1+1);
      vKf[1].NewAthenaArray(nc3+1,nc1);
      if (orbital_advection_active) {
        orbital_cons.NewAthenaArray(NHYDRO, nc3, nc1, nc2+onx+1);
        orc.NewAthenaArray(nc3, nc1);
        ofc.NewAthenaArray(nc3, nc1);
        if (orbital_refinement) {
          vKc_coarse.NewAthenaArray((nc3+2*NGHOST)/2, (nc1+2*NGHOST)/2);
          ofc_coarse.NewAthenaArray((nc3+2*NGHOST)/2, (nc1+2*NGHOST)/2);
        }
        if (MAGNETIC_FIELDS_ENABLED) {
          orbital_b1.NewAthenaArray(nc3, nc1+1, nc2+onx+1);
          orbital_b2.NewAthenaArray(nc3+1, nc1, nc2+onx+1);
          orf[0].NewAthenaArray(nc3, nc1+1);
          orf[1].NewAthenaArray(nc3+1, nc1);
          off[0].NewAthenaArray(nc3, nc1+1);
          off[1].NewAthenaArray(nc3+1, nc1);
          if (orbital_refinement) {
            vKf_coarse[0].NewAthenaArray((nc3+2*NGHOST)/2, (nc1+2*NGHOST)/2+1);
            off_coarse[0].NewAthenaArray((nc3+2*NGHOST)/2, (nc1+2*NGHOST)/2+1);
            vKf_coarse[1].NewAthenaArray((nc3+2*NGHOST)/2+1, (nc1+2*NGHOST)/2);
            off_coarse[1].NewAthenaArray((nc3+2*NGHOST)/2+1, (nc1+2*NGHOST)/2);
          }
        }
        if (NSCALARS>0) {
          orbital_scalar.NewAthenaArray(NSCALARS, nc3, nc1, nc2+onx+1);
        }
      }
    } else if (orbital_direction==2) { // spherical_polar
      vKc.NewAthenaArray(nc2, nc1);
      dvKc1.NewAthenaArray(nc2, nc1);
      dvKc2.NewAthenaArray(nc2, nc1);
      vKf[0].NewAthenaArray(nc2,nc1+1);
      vKf[1].NewAthenaArray(nc2+1,nc1);
      if (orbital_advection_active) {
        orbital_cons.NewAthenaArray(NHYDRO, nc2, nc1, nc3+onx+1);
        orc.NewAthenaArray(nc2, nc1);
        ofc.NewAthenaArray(nc2, nc1);
        if (orbital_refinement) {
          vKc_coarse.NewAthenaArray((nc2+2*NGHOST)/2, (nc1+2*NGHOST)/2);
          ofc_coarse.NewAthenaArray((nc2+2*NGHOST)/2, (nc1+2*NGHOST)/2);
        }
        if (MAGNETIC_FIELDS_ENABLED) {
          orbital_b1.NewAthenaArray(nc2, nc1+1, nc3+onx+1);
          orbital_b2.NewAthenaArray(nc2+1, nc1, nc3+onx+1);
          orf[0].NewAthenaArray(nc2, nc1+1);
          orf[1].NewAthenaArray(nc2+1, nc1);
          off[0].NewAthenaArray(nc2, nc1+1);
          off[1].NewAthenaArray(nc2+1, nc1);
          if (orbital_refinement) {
            vKf_coarse[0].NewAthenaArray((nc2+2*NGHOST)/2, (nc1+2*NGHOST)/2+1);
            off_coarse[0].NewAthenaArray((nc2+2*NGHOST)/2, (nc1+2*NGHOST)/2+1);
            vKf_coarse[1].NewAthenaArray((nc2+2*NGHOST)/2+1, (nc1+2*NGHOST)/2);
            off_coarse[1].NewAthenaArray((nc2+2*NGHOST)/2+1, (nc1+2*NGHOST)/2);
          }
        }
        if (NSCALARS>0) {
          orbital_scalar.NewAthenaArray(NSCALARS, nc2, nc1, nc3+onx+1);
        }
      }
    }
    w_orb.NewAthenaArray(NHYDRO,nc3,nc2,nc1);
    u_orb.NewAthenaArray(NHYDRO,nc3,nc2,nc1);

    if (orbital_advection_active) {
      pflux.NewAthenaArray(onx+2*NGHOST+1);
      if (orbital_refinement) {
        u_coarse_send.NewAthenaArray(NHYDRO, (nc3+2*NGHOST)/2,
                                     (nc2+2*NGHOST)/2, (nc1+2*NGHOST)/2);
        u_coarse_recv.NewAthenaArray(NHYDRO, (nc3+2*NGHOST)/2,
                                     (nc2+2*NGHOST)/2, (nc1+2*NGHOST)/2);
        u_temp.NewAthenaArray(NHYDRO, nc3, nc2, nc1);
        if (NSCALARS>0) {
          s_coarse_send.NewAthenaArray(NSCALARS, (nc3+2*NGHOST)/2,
                                       (nc2+2*NGHOST)/2, (nc1+2*NGHOST)/2);
          s_coarse_recv.NewAthenaArray(NSCALARS, (nc3+2*NGHOST)/2,
                                       (nc2+2*NGHOST)/2, (nc1+2*NGHOST)/2);
          s_temp.NewAthenaArray(NSCALARS, nc3, nc2, nc1);
        }
        if (MAGNETIC_FIELDS_ENABLED) {
          b1_coarse_send.NewAthenaArray((nc3+2*NGHOST)/2,
                                        (nc2+2*NGHOST)/2, (nc1+2*NGHOST)/2+1);
          b_coarse_recv.x1f.NewAthenaArray((nc3+2*NGHOST)/2,
                                          (nc2+2*NGHOST)/2, (nc1+2*NGHOST)/2+1);
          b_coarse_recv.x2f.NewAthenaArray((nc3+2*NGHOST)/2,
                                          (nc2+2*NGHOST)/2+1, (nc1+2*NGHOST)/2);
          b_coarse_recv.x3f.NewAthenaArray((nc3+2*NGHOST)/2+1,
                                          (nc2+2*NGHOST)/2, (nc1+2*NGHOST)/2);
          b_temp.x1f.NewAthenaArray(nc3, nc2, nc1+1);
          b_temp.x2f.NewAthenaArray(nc3, nc2+1, nc1);
          b_temp.x3f.NewAthenaArray(nc3+1, nc2, nc1);
          if (orbital_direction == 1) {
            b2_coarse_send.NewAthenaArray((nc3+2*NGHOST)/2+1,
                                          (nc2+2*NGHOST)/2, (nc1+2*NGHOST)/2);
          } else if (orbital_direction == 2) {
            b2_coarse_send.NewAthenaArray((nc3+2*NGHOST)/2,
                                          (nc2+2*NGHOST)/2+1, (nc1+2*NGHOST)/2);
          }
        }
      }
      // call OrbitalAdvectionBoundaryVariable
      orb_bc = new OrbitalBoundaryCommunication(this);
    }
  }

  // preparation for shear_periodic boundary
  if ((orbital_advection_active || pm_->shear_periodic)) {
    int pnum = onx+2*NGHOST+1;
    if (pm_->shear_periodic && MAGNETIC_FIELDS_ENABLED) pnum++;
    if (xorder>2) {
      for (int n=0; n<13; n++) {
        d_src[n].NewAthenaArray(pnum);
      }
    }
  }
}


// destructor
OrbitalAdvection::~OrbitalAdvection() {
  if (orbital_advection_active) {
    // destroy OrbitalBoundaryCommunication
    delete orb_bc;
  }
}


//----------------------------------------------------------------------------------------
//! \fn void OrbitalAdvection::InitializeOrbitalAdvection()
//! \brief Setup for OrbitalAdvection in void Mesh::Initialize

void OrbitalAdvection::InitializeOrbitalAdvection() {
  // set orbital velocity
  SetVKc();
  SetVKf();
  SetDvKc();
  if (orbital_advection_active) {
    if (orbital_refinement) {
      SetVKcCoarse();
      if (MAGNETIC_FIELDS_ENABLED) {
        SetVKfCoarse();
     }
    }

    if(orbital_uniform_mesh) { // uniform mesh
      // set dx
      // TODO(tomo-ono): For the consistency, not use dx2v or dx3v
      if (orbital_direction == 1) { // cartesian or cylindrical
        dx = (pm_->mesh_size.x2max-pm_->mesh_size.x2min)
             /static_cast<Real>((pm_->nrbx2
               *(1L<<(pmb_->loc.level-pm_->root_level))
              )*pmb_->block_size.nx2);
      } else { // spherical_polar
        dx = (pm_->mesh_size.x3max-pm_->mesh_size.x3min)
             /static_cast<Real>((pm_->nrbx3
               *(1L<<(pmb_->loc.level-pm_->root_level))
              )*pmb_->block_size.nx3);
      }
    }
    // else { // non-uniform mesh
    // }

    //set vK_max, vK_min
    vK_max = -(FLT_MAX);
    vK_min = (FLT_MAX);
    if (orbital_direction == 1) { // cartesian or cylindrical
      // cell center
      for (int k=pmb_->ks; k<= pmb_->ke; k++) {
        for (int i=pmb_->is; i<= pmb_->ie; i++) {
          Real pvk = 1.0/pco_->h2v(i);
          vK_max   = std::max(vK_max, vKc(k,i)*pvk);
          vK_min   = std::min(vK_min, vKc(k,i)*pvk);
        }
      }
      if (MAGNETIC_FIELDS_ENABLED) {
        // x1 surface
        for (int k=pmb_->ks; k<= pmb_->ke; k++) {
          for (int i=pmb_->is; i<= pmb_->ie+1; i++) {
            Real pvk = 1.0/pco_->h2f(i);
            vK_max   = std::max(vK_max, vKf[0](k,i)*pvk);
            vK_min   = std::min(vK_min, vKf[0](k,i)*pvk);
          }
        }
        // x3 surface
        for (int k=pmb_->ks; k<= pmb_->ke+1; k++) {
          for (int i=pmb_->is; i<= pmb_->ie; i++) {
            Real pvk = 1.0/pco_->h2v(i);
            vK_max   = std::max(vK_max, vKf[1](k,i)*pvk);
            vK_min   = std::min(vK_min, vKf[1](k,i)*pvk);
          }
        }
      }
    } else if (orbital_direction == 2) { // spherical_polar
      // cell center
      for (int j=pmb_->js; j<= pmb_->je; ++j) {
        for (int i=pmb_->is; i<= pmb_->ie; ++i) {
          Real pvk = 1.0/(pco_->h2v(i)*pco_->h32v(j));
          vK_max   = std::max(vK_max, vKc(j,i)*pvk);
          vK_min   = std::min(vK_min, vKc(j,i)*pvk);
        }
      }
      if (MAGNETIC_FIELDS_ENABLED) {
        // x1 surface
        for (int j=pmb_->js; j<= pmb_->je; ++j) {
          for (int i=pmb_->is; i<= pmb_->ie+1; ++i) {
            Real pvk = 1.0/(pco_->h2f(i)*pco_->h32v(j));
            vK_max   = std::max(vK_max, vKf[0](j,i)*pvk);
            vK_min   = std::min(vK_min, vKf[0](j,i)*pvk);
          }
        }
        // x2 surface
        for (int j=pmb_->js; j<= pmb_->je+1; ++j) {
          for (int i=pmb_->is; i<= pmb_->ie; ++i) {
            Real pvk = 1.0/(pco_->h2v(i)*pco_->h32f(j));
            vK_max   = std::max(vK_max, vKf[1](j,i)*pvk);
            vK_min   = std::min(vK_min, vKf[1](j,i)*pvk);
          }
        }
      }
    }

    // set min_dt
    int mylevel = pmb_->loc.level;
    int lblevel = mylevel, rblevel = mylevel;
    for(int n=0; n<pbval_->nneighbor; n++) {
      NeighborBlock& nb = pbval_->neighbor[n];
      if (nb.ni.type != NeighborConnect::face) break;
      if (orbital_direction == 1) {
        if (nb.fid == BoundaryFace::inner_x2)
          lblevel = nb.snb.level;
        else if (nb.fid == BoundaryFace::outer_x2)
          rblevel = nb.snb.level;
      } else if (orbital_direction == 2) {
        if (nb.fid == BoundaryFace::inner_x3)
          lblevel = nb.snb.level;
        else if (nb.fid == BoundaryFace::outer_x3)
          rblevel = nb.snb.level;
      }
    }

    min_dt = (FLT_MAX);
    // restrictions from meshblock size
    if(orbital_uniform_mesh) { // uniform mesh
      if(vK_max>0.0) {
        if(lblevel > mylevel)
          min_dt = std::min(min_dt, dx*(onx/2-xgh)/vK_max);
        else
          min_dt = std::min(min_dt, dx*(onx-xgh)/vK_max);
      }
      if(vK_min<0.0) {
        if(rblevel > mylevel)
          min_dt = std::min(min_dt, -dx*(onx/2-xgh)/vK_min);
        else
          min_dt = std::min(min_dt, -dx*(onx-xgh)/vK_min);
      }
    }
    // else {
    //  }

    // restrictions from derivatives of orbital velocity
    Real coef = 1.0;
    if (orbital_splitting_order == 1) coef = 0.5;
    if(orbital_direction == 1) {
      for(int k=pmb_->ks; k<=pmb_->ke; k++) {
        for(int i=pmb_->is; i<=pmb_->ie; i++) {
          Real vorb_c   = vKc(k,i)/pco_->h2v(i);
          Real vorb_min = vorb_c;
          Real vorb_max = vorb_c;
          for (int neighbor=1; neighbor<=xgh; neighbor++) {
            Real vorb_m = vKc(k,i-neighbor)/pco_->h2v(i-neighbor);
            Real vorb_p = vKc(k,i+neighbor)/pco_->h2v(i+neighbor);
            vorb_min = std::min(vorb_min, std::min(vorb_m, vorb_p));
            vorb_max = std::max(vorb_max, std::max(vorb_m, vorb_p));
          }
          if (nc3>1) {
            for (int neighbor=1; neighbor<=xgh; neighbor++) {
              Real vorb_m = vKc(k-neighbor,i)/pco_->h2v(i);
              Real vorb_p = vKc(k+neighbor,i)/pco_->h2v(i);
              vorb_min = std::min(vorb_min, std::min(vorb_m, vorb_p));
              vorb_max = std::max(vorb_max, std::max(vorb_m, vorb_p));
            }
          }
          if (vorb_min == vorb_max) continue;
          if(orbital_uniform_mesh) { // uniform mesh
            min_dt = std::min(min_dt, coef*dx/(vorb_max-vorb_min));
          }
          // else { // non-uniform mesh
          // }
        }
      }
    } else if (orbital_direction == 2) {
      for(int j=pmb_->js; j<=pmb_->je; j++) {
        for(int i=pmb_->is; i<=pmb_->ie; i++) {
          Real vorb_c   = vKc(j,i)/(pco_->h2v(i)*pco_->h32v(j));
          Real vorb_min = vorb_c;
          Real vorb_max = vorb_c;
          for (int neighbor=1; neighbor<=xgh; neighbor++) {
            Real vorb_m = vKc(j,i-neighbor)/(pco_->h2v(i-neighbor)*pco_->h32v(j));
            Real vorb_p = vKc(j,i+neighbor)/(pco_->h2v(i+neighbor)*pco_->h32v(j));
            vorb_min = std::min(vorb_min, std::min(vorb_m, vorb_p));
            vorb_max = std::max(vorb_max, std::max(vorb_m, vorb_p));
          }
          for (int neighbor=1; neighbor<=xgh; neighbor++) {
            Real vorb_m = vKc(j-neighbor,i)/(pco_->h2v(i)*pco_->h32v(j-neighbor));
            Real vorb_p = vKc(j+neighbor,i)/(pco_->h2v(i)*pco_->h32v(j+neighbor));
            vorb_min = std::min(vorb_min, std::min(vorb_m, vorb_p));
            vorb_max = std::max(vorb_max, std::max(vorb_m, vorb_p));
          }
          if (vorb_min == vorb_max) continue;
          if(orbital_uniform_mesh) { // uniform mesh
            min_dt = std::min(min_dt, coef*dx/(vorb_max-vorb_min));
          }
          // else { // non-uniform mesh
          // }
        }
      }
    }
  } // orbital_advection_active
  return;
}


//----------------------------------------------------------------------------------------
//! \fn Real OrbitalAdvection::NewOrbitalAdvectionDt()
//! \brief Calculate time step for OrbitalAdvection

Real OrbitalAdvection::NewOrbitalAdvectionDt() {
  return min_dt;
}


//----------------------------------------------------------------------------------------
//! \fn void OrbitalAdvection::SetVKc()
//! \brief Calculate Orbital Velocity at cell centers

void OrbitalAdvection::SetVKc() {
  if (orbital_direction == 1) {
    int il = pmb_->is-(NGHOST); int kl = pmb_->ks;
    int iu = pmb_->ie+(NGHOST); int ku = pmb_->ke;
    if (nc3>1) {
      kl -= NGHOST; ku += NGHOST;
    }
    for(int k=kl; k<=ku; k++) {
      Real z_ = pco_->x3v(k);
#pragma omp simd
      for(int i=il; i<=iu; i++) {
        Real x_ = pco_->x1v(i);
        vKc(k,i)   = OrbitalVelocity(this, x_, 0.0, z_);
      }
    }
  } else if (orbital_direction == 2) {
    int il = pmb_->is-(NGHOST); int jl = pmb_->js-NGHOST;
    int iu = pmb_->ie+(NGHOST); int ju = pmb_->je+NGHOST;
    for(int j=jl; j<=ju; j++) {
      Real y_ = pco_->x2v(j);
#pragma omp simd
      for(int i=il; i<=iu; i++) {
        Real x_ = pco_->x1v(i);
        vKc(j,i)   = OrbitalVelocity(this, x_, y_, 0.0);
      }
    }
  }
}


//----------------------------------------------------------------------------------------
//! \fn void OrbitalAdvection::SetVKf()
//! \brief Calculate Orbital Velocity at cell faces

void OrbitalAdvection::SetVKf() {
  if (orbital_direction == 1) {
    int il = pmb_->is-(NGHOST); int kl = pmb_->ks;
    int iu = pmb_->ie+(NGHOST); int ku = pmb_->ke;
    if (nc3>1) {
      kl -= NGHOST; ku += NGHOST;
    }
    for(int k=kl; k<=ku; k++) {
      Real z_ = pco_->x3v(k);
#pragma omp simd
      for(int i=il; i<=iu+1; i++) {
        Real x_ = pco_->x1f(i);
        vKf[0](k,i)   = OrbitalVelocity(this, x_, 0.0, z_);
      }
    }
    for(int k=kl; k<=ku+1; k++) {
      Real z_ = pco_->x3f(k);
#pragma omp simd
      for(int i=il; i<=iu; i++) {
        Real x_ = pco_->x1v(i);
        vKf[1](k,i)   = OrbitalVelocity(this, x_, 0.0, z_);
      }
    }
  } else if (orbital_direction == 2) {
    int il = pmb_->is-(NGHOST); int jl = pmb_->js-NGHOST;
    int iu = pmb_->ie+(NGHOST); int ju = pmb_->je+NGHOST;
    for(int j=jl; j<=ju; j++) {
      Real y_ = pco_->x2v(j);
#pragma omp simd
      for(int i=il; i<=iu+1; i++) {
        Real x_ = pco_->x1f(i);
        vKf[0](j,i)   = OrbitalVelocity(this, x_, y_, 0.0);
      }
    }
    for(int j=jl; j<=ju+1; j++) {
      Real y_ = pco_->x2f(j);
#pragma omp simd
      for(int i=il; i<=iu; i++) {
        Real x_ = pco_->x1v(i);
        vKf[1](j,i)   = OrbitalVelocity(this, x_, y_, 0.0);
      }
    }
  }
}


//----------------------------------------------------------------------------------------
//! \fn void OrbitalAdvection::SetDvKc()
//! \brief Calculate Orbital Velocity derivatives at cell centers

void OrbitalAdvection::SetDvKc() {
  if (orbital_direction == 1) { // cartesian or cylindrical
    int il = pmb_->is-(NGHOST); int kl = pmb_->ks;
    int iu = pmb_->ie+(NGHOST); int ku = pmb_->ke;
    if (nc3>1) {
      kl -= NGHOST; ku += NGHOST;
    }
    if (OrbitalVelocityDerivative[0] == nullptr) {
      // calculate dvK using user-defined vKf
      for(int k=kl; k<=ku; k++) {
#pragma omp simd
        for(int i=il; i<=iu; i++) {
          dvKc1(k,i) = (vKf[0](k,i+1)-vKf[0](k,i))/pco_->dx1f(i);
        }
      }
    } else {
      // set dvK from user-defined dvK
      for(int k=kl; k<=ku; k++) {
        Real z_ = pco_->x3v(k);
#pragma omp simd
        for(int i=il; i<=iu; i++) {
          Real x_ = pco_->x1v(i);
          dvKc1(k,i)   = OrbitalVelocityDerivative[0](this, x_, 0.0, z_);
        }
      }
    }
    if (nc3>1) { // 3D
      if (OrbitalVelocityDerivative[1] == nullptr) {
        // calculate dvK using user-defined vKf
        for(int k=kl; k<=ku; k++) {
#pragma omp simd
          for(int i=il; i<=iu; i++) {
            dvKc2(k,i) = (vKf[1](k+1,i)-vKf[1](k,i))/pco_->dx3f(k);
          }
        }
      } else {
        // set dvK from user-defined dvK
        for(int k=kl; k<=ku; k++) {
          Real z_ = pco_->x3v(k);
#pragma omp simd
          for(int i=il; i<=iu; i++) {
            Real x_ = pco_->x1v(i);
            dvKc2(k,i)   = OrbitalVelocityDerivative[1](this, x_, 0.0, z_);
          }
        }
      }
    } else { // 2D
      for(int k=kl; k<=ku; k++) {
#pragma omp simd
        for(int i=il; i<=iu; i++) {
          dvKc2(k,i)   = 0.0;
        }
      }
    }
  } else if (orbital_direction == 2) { // spherical_polar
    int il = pmb_->is-(NGHOST); int jl = pmb_->js-NGHOST;
    int iu = pmb_->ie+(NGHOST); int ju = pmb_->je+NGHOST;
    if (OrbitalVelocityDerivative[0] == nullptr) {
      // calculate dvK using user-defined vKf
      for(int j=jl; j<=ju; j++) {
#pragma omp simd
        for(int i=il; i<=iu; i++) {
          dvKc1(j,i) = (vKf[0](j,i+1)-vKf[0](j,i))/pco_->dx1f(i);
        }
      }
    } else {
      // set dvK from user-defined dvK
      for(int j=jl; j<=ju; j++) {
        Real y_ = pco_->x2v(j);
#pragma omp simd
        for(int i=il; i<=iu; i++) {
          Real x_ = pco_->x1v(i);
          dvKc1(j,i)   = OrbitalVelocityDerivative[0](this, x_, y_, 0.0);
        }
      }
    }
    if (OrbitalVelocityDerivative[1] == nullptr) {
      // calculate dvK using user-defined vKf
      for(int j=jl; j<=ju; j++) {
#pragma omp simd
        for(int i=il; i<=iu; i++) {
          dvKc2(j,i) = (vKf[1](j+1,i)-vKf[1](j,i))/pco_->dx2f(j);
        }
      }
    } else {
      // set dvK from user-defined dvK
      for(int j=jl; j<=ju; j++) {
        Real y_ = pco_->x2v(j);
#pragma omp simd
        for(int i=il; i<=iu; i++) {
          Real x_ = pco_->x1v(i);
          dvKc2(j,i)   = OrbitalVelocityDerivative[1](this, x_, y_, 0.0);
        }
      }
    }
  }
}


//----------------------------------------------------------------------------------------
//! \fn void OrbitalAdvection::SetVKcCoarse()
//! \brief Calculate Orbital Velocity at cell centers for coarse cells

void OrbitalAdvection::SetVKcCoarse() {
  Coordinates *cpco = pmb_->pmr->pcoarsec;
  if (orbital_direction == 1) {
    int il = pmb_->cis-(NGHOST); int kl = pmb_->cks;
    int iu = pmb_->cie+(NGHOST); int ku = pmb_->cke;
    if (nc3>1) {
      kl -= NGHOST;
      ku += NGHOST;
    }
    for(int k=kl; k<=ku; k++) {
      Real z_ = cpco->x3v(k);
#pragma omp simd
      for(int i=il; i<=iu; i++) {
        Real x_ = cpco->x1v(i);
        vKc_coarse(k,i)   = OrbitalVelocity(this, x_, 0.0, z_);
      }
    }
  } else if (orbital_direction == 2) {
    int il = pmb_->cis-(NGHOST); int jl = pmb_->cjs-NGHOST;
    int iu = pmb_->cie+(NGHOST); int ju = pmb_->cje+NGHOST;
    for(int j=jl; j<=ju; j++) {
      Real y_ = cpco->x2v(j);
#pragma omp simd
      for(int i=il; i<=iu; i++) {
        Real x_ = cpco->x1v(i);
        vKc_coarse(j,i)   = OrbitalVelocity(this, x_, y_, 0.0);
      }
    }
  }
}


//----------------------------------------------------------------------------------------
//! \fn void OrbitalAdvection::SetVKfCoarse()
//! \brief Calculate Orbital Velocity at cell faces for coarse cells
void OrbitalAdvection::SetVKfCoarse() {
  Coordinates *cpco = pmb_->pmr->pcoarsec;
  if (orbital_direction == 1) {
    int il = pmb_->cis-(NGHOST); int kl = pmb_->cks;
    int iu = pmb_->cie+(NGHOST); int ku = pmb_->cke;
    if (nc3>1) {
      kl -= NGHOST;
      ku += NGHOST;
    }
    for(int k=kl; k<=ku; k++) {
      Real z_ = cpco->x3v(k);
#pragma omp simd
      for(int i=il; i<=iu+1; i++) {
        Real x_ = cpco->x1f(i);
        vKf_coarse[0](k,i)   = OrbitalVelocity(this, x_, 0.0, z_);
      }
    }
    for(int k=kl; k<=ku+1; k++) {
      Real z_ = cpco->x3f(k);
#pragma omp simd
      for(int i=il; i<=iu; i++) {
        Real x_ = cpco->x1v(i);
        vKf_coarse[1](k,i)   = OrbitalVelocity(this, x_, 0.0, z_);
      }
    }
  } else if (orbital_direction == 2) {
    int il = pmb_->cis-(NGHOST); int jl = pmb_->cjs-NGHOST;
    int iu = pmb_->cie+(NGHOST); int ju = pmb_->cje+NGHOST;
    for(int j=jl; j<=ju; j++) {
      Real y_ = cpco->x2v(j);
#pragma omp simd
      for(int i=il; i<=iu+1; i++) {
        Real x_ = cpco->x1f(i);
        vKf_coarse[0](j,i)   = OrbitalVelocity(this, x_, y_, 0.0);
      }
    }
    for(int j=jl; j<=ju+1; j++) {
      Real y_ = cpco->x2f(j);
#pragma omp simd
      for(int i=il; i<=iu; i++) {
        Real x_ = cpco->x1v(i);
        vKf_coarse[1](j,i)   = OrbitalVelocity(this, x_, y_, 0.0);
      }
    }
  }
}
