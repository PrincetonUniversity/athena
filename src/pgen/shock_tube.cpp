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
//! \file shock_tube.cpp
//  \brief Problem generator for shock tube problems.  
//
// Problem generator for shock tube (1-D Riemann) problems. Initializes plane-parallel
// shock along x1 (in 1D, 2D, 3D), along x2 (in 2D, 3D), and along x3 (in 3D).
//======================================================================================

// C++ headers
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../hydro/eos/eos.hpp"
#include "../coordinates/coordinates.hpp"

//======================================================================================
//! \fn void Mesh::InitUserMeshProperties(ParameterInput *pin)
//  \brief Init the Mesh properties
//======================================================================================

void Mesh::InitUserMeshProperties(ParameterInput *pin)
{
  return;
}


//======================================================================================
//! \fn void Mesh::TerminateUserMeshProperties(void)
//  \brief Clean up the Mesh properties
//======================================================================================

void Mesh::TerminateUserMeshProperties(void)
{
  // nothing to do
  return;
}


//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the shock tube tests
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  std::stringstream msg;

// parse shock direction: {1,2,3} -> {x1,x2,x3}

  int shk_dir = pin->GetInteger("problem","shock_dir"); 

// parse shock location (must be inside grid)

  Real xshock = pin->GetReal("problem","xshock"); 
  if (shk_dir == 1 && (xshock < pmy_mesh->mesh_size.x1min ||
                       xshock > pmy_mesh->mesh_size.x1max)) {
    msg << "### FATAL ERROR in Problem Generator" << std::endl << "xshock="
        << xshock << " lies outside x1 domain for shkdir=" << shk_dir << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (shk_dir == 2 && (xshock < pmy_mesh->mesh_size.x2min ||
                       xshock > pmy_mesh->mesh_size.x2max)) {
    msg << "### FATAL ERROR in Problem Generator" << std::endl << "xshock="
        << xshock << " lies outside x2 domain for shkdir=" << shk_dir << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (shk_dir == 3 && (xshock < pmy_mesh->mesh_size.x3min ||
                       xshock > pmy_mesh->mesh_size.x3max)) {
    msg << "### FATAL ERROR in Problem Generator" << std::endl << "xshock="
        << xshock << " lies outside x3 domain for shkdir=" << shk_dir << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

// Parse left state read from input file: dl,ul,vl,wl,[pl]

  Real wl[NHYDRO+NFIELD];
  wl[IDN] = pin->GetReal("problem","dl");
  wl[IVX] = pin->GetReal("problem","ul");
  wl[IVY] = pin->GetReal("problem","vl");
  wl[IVZ] = pin->GetReal("problem","wl");
  if (NON_BAROTROPIC_EOS) wl[IEN] = pin->GetReal("problem","pl");
  if (MAGNETIC_FIELDS_ENABLED) {
    wl[NHYDRO  ] = pin->GetReal("problem","bxl");
    wl[NHYDRO+1] = pin->GetReal("problem","byl");
    wl[NHYDRO+2] = pin->GetReal("problem","bzl");
  }

// Parse right state read from input file: dr,ur,vr,wr,[pr]

  Real wr[NHYDRO+NFIELD];
  wr[IDN] = pin->GetReal("problem","dr");
  wr[IVX] = pin->GetReal("problem","ur");
  wr[IVY] = pin->GetReal("problem","vr");
  wr[IVZ] = pin->GetReal("problem","wr");
  if (NON_BAROTROPIC_EOS) wr[IEN] = pin->GetReal("problem","pr");
  if (MAGNETIC_FIELDS_ENABLED) {
    wr[NHYDRO  ] = pin->GetReal("problem","bxr");
    wr[NHYDRO+1] = pin->GetReal("problem","byr");
    wr[NHYDRO+2] = pin->GetReal("problem","bzr");
  }

// Initialize the discontinuity in the Hydro variables ---------------------------------

  switch(shk_dir) {

//--- shock in 1-direction
  case 1:
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        if (pcoord->x1v(i) < xshock) {
          phydro->u(IDN,k,j,i) = wl[IDN];
          phydro->u(IM1,k,j,i) = wl[IVX]*wl[IDN];
          phydro->u(IM2,k,j,i) = wl[IVY]*wl[IDN];
          phydro->u(IM3,k,j,i) = wl[IVZ]*wl[IDN];
          if (NON_BAROTROPIC_EOS) phydro->u(IEN,k,j,i) =
            wl[IEN]/(phydro->peos->GetGamma() - 1.0)
            + 0.5*wl[IDN]*(wl[IVX]*wl[IVX] + wl[IVY]*wl[IVY] + wl[IVZ]*wl[IVZ]);
        } else {
          phydro->u(IDN,k,j,i) = wr[IDN];
          phydro->u(IM1,k,j,i) = wr[IVX]*wr[IDN];
          phydro->u(IM2,k,j,i) = wr[IVY]*wr[IDN];
          phydro->u(IM3,k,j,i) = wr[IVZ]*wr[IDN];
          if (NON_BAROTROPIC_EOS) phydro->u(IEN,k,j,i) =
            wr[IEN]/(phydro->peos->GetGamma() - 1.0)
            + 0.5*wr[IDN]*(wr[IVX]*wr[IVX] + wr[IVY]*wr[IVY] + wr[IVZ]*wr[IVZ]);
        }
      }
    }}
    break;

//--- shock in 2-direction
  case 2:
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      if (pcoord->x2v(j) < xshock) {
        for (int i=is; i<=ie; ++i) {
          phydro->u(IDN,k,j,i) = wl[IDN];
          phydro->u(IM2,k,j,i) = wl[IVX]*wl[IDN];
          phydro->u(IM3,k,j,i) = wl[IVY]*wl[IDN];
          phydro->u(IM1,k,j,i) = wl[IVZ]*wl[IDN];
          if (NON_BAROTROPIC_EOS) phydro->u(IEN,k,j,i) =
            wl[IEN]/(phydro->peos->GetGamma() - 1.0)
            + 0.5*wl[IDN]*(wl[IVX]*wl[IVX] + wl[IVY]*wl[IVY] + wl[IVZ]*wl[IVZ]);
        }
      } else {
        for (int i=is; i<=ie; ++i) {
          phydro->u(IDN,k,j,i) = wr[IDN];
          phydro->u(IM2,k,j,i) = wr[IVX]*wr[IDN];
          phydro->u(IM3,k,j,i) = wr[IVY]*wr[IDN];
          phydro->u(IM1,k,j,i) = wr[IVZ]*wr[IDN];
          if (NON_BAROTROPIC_EOS) phydro->u(IEN,k,j,i) =
            wr[IEN]/(phydro->peos->GetGamma() - 1.0)
            + 0.5*wr[IDN]*(wr[IVX]*wr[IVX] + wr[IVY]*wr[IVY] + wr[IVZ]*wr[IVZ]);
        }
      }
    }}
    break;

//--- shock in 3-direction

  case 3:
    for (int k=ks; k<=ke; ++k) {
      if (pcoord->x3v(k) < xshock) {
        for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          phydro->u(IDN,k,j,i) = wl[IDN];
          phydro->u(IM3,k,j,i) = wl[IVX]*wl[IDN];
          phydro->u(IM1,k,j,i) = wl[IVY]*wl[IDN];
          phydro->u(IM2,k,j,i) = wl[IVZ]*wl[IDN];
          if (NON_BAROTROPIC_EOS) phydro->u(IEN,k,j,i) =
            wl[IEN]/(phydro->peos->GetGamma() - 1.0)
            + 0.5*wl[IDN]*(wl[IVX]*wl[IVX] + wl[IVY]*wl[IVY] + wl[IVZ]*wl[IVZ]);
        }}
      } else {
        for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          phydro->u(IDN,k,j,i) = wr[IDN];
          phydro->u(IM3,k,j,i) = wr[IVX]*wr[IDN];
          phydro->u(IM1,k,j,i) = wr[IVY]*wr[IDN];
          phydro->u(IM2,k,j,i) = wr[IVZ]*wr[IDN];
          if (NON_BAROTROPIC_EOS) phydro->u(IEN,k,j,i) =
            wr[IEN]/(phydro->peos->GetGamma() - 1.0)
            + 0.5*wr[IDN]*(wr[IVX]*wr[IVX] + wr[IVY]*wr[IVY] + wr[IVZ]*wr[IVZ]);
        }}
      }
    }
    break;

  default:
    msg << "### FATAL ERROR in Problem Generator" << std::endl
        << "shock_dir=" << shk_dir << " must be either 1,2, or 3" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

// now set face-centered (interface) magnetic fields -----------------------------------

  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        if (shk_dir==1 && pcoord->x1v(i) < xshock) {
          pfield->b.x1f(k,j,i) = wl[NHYDRO  ];
          pfield->b.x2f(k,j,i) = wl[NHYDRO+1];
          pfield->b.x3f(k,j,i) = wl[NHYDRO+2];
        } else if (shk_dir==2 && pcoord->x2v(j) < xshock) {
          pfield->b.x1f(k,j,i) = wl[NHYDRO+2];
          pfield->b.x2f(k,j,i) = wl[NHYDRO  ];
          pfield->b.x3f(k,j,i) = wl[NHYDRO+1];
        } else if (shk_dir==3 && pcoord->x3v(k) < xshock) {
          pfield->b.x1f(k,j,i) = wl[NHYDRO+1];
          pfield->b.x2f(k,j,i) = wl[NHYDRO+2];
          pfield->b.x3f(k,j,i) = wl[NHYDRO];
        }

        if (shk_dir==1 && pcoord->x1v(i) >= xshock) {
          pfield->b.x1f(k,j,i) = wr[NHYDRO  ];
          pfield->b.x2f(k,j,i) = wr[NHYDRO+1];
          pfield->b.x3f(k,j,i) = wr[NHYDRO+2];
        } else if (shk_dir==2 && pcoord->x2v(j) >= xshock) {
          pfield->b.x1f(k,j,i) = wr[NHYDRO+2];
          pfield->b.x2f(k,j,i) = wr[NHYDRO  ];
          pfield->b.x3f(k,j,i) = wr[NHYDRO+1];
        } else if (shk_dir==3 && pcoord->x3v(k) >= xshock)  {
          pfield->b.x1f(k,j,i) = wr[NHYDRO+1];
          pfield->b.x2f(k,j,i) = wr[NHYDRO+2];
          pfield->b.x3f(k,j,i) = wr[NHYDRO];
        }
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) += 0.5*(SQR(pfield->b.x1f(k,j,i)) + SQR(pfield->b.x2f(k,j,i)) +
             SQR(pfield->b.x3f(k,j,i)));
        }
      }
    }}

// end by adding bi.x1 at ie+1, bi.x2 at je+1, and bi.x3 at ke+1

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      pfield->b.x1f(k,j,ie+1) = pfield->b.x1f(k,j,ie);
    }}
    for (int k=ks; k<=ke; ++k) {
    for (int i=is; i<=ie; ++i) {
      pfield->b.x2f(k,je+1,i) = pfield->b.x2f(k,je,i);
    }}
    for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie; ++i) {
      pfield->b.x3f(ke+1,j,i) = pfield->b.x3f(ke,j,i);
    }}
  }

  return;
}


//======================================================================================
//! \fn void MeshBlock::UserWorkInLoop(void)
//  \brief User-defined work function for every time step
//======================================================================================

void MeshBlock::UserWorkInLoop(void)
{
  // nothing to do
  return;
}

