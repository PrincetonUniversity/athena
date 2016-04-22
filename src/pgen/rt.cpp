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
//! \file rt.c
//  \brief Problem generator for RT instabilty.
//
// Note the gravitational acceleration is hardwired to be 0.1. Density difference is
// hardwired to be 2.0 in 2D, and is set by the input parameter <problem>/rhoh in 3D
// (default value is 3.0). This reproduces 2D results of Liska & Wendroff, 3D results of
// Dimonte et al.
// 
// FOR 2D HYDRO:
// Problem domain should be -1/6 < x < 1/6; -0.5 < y < 0.5 with gamma=1.4 to match Liska
// & Wendroff. Interface is at y=0; perturbation added to Vy. Gravity acts in y-dirn.
// Special reflecting boundary conditions added in x2 to improve hydrostatic eqm
// (prevents launching of weak waves) Atwood number A=(d2-d1)/(d2+d1)=1/3. Options:
//     iprob = 1  -- Perturb V2 using single mode
//     iprob != 1 -- Perturb V2 using multiple mode
//
// FOR 3D:
// Problem domain should be -.05 < x < .05; -.05 < y < .05, -.1 < z < .1, gamma=5/3 to
// match Dimonte et al.  Interface is at z=0; perturbation added to Vz. Gravity acts in
// z-dirn. Special reflecting boundary conditions added in x3.  A=1/2.  Options:
//     iprob = 1 -- Perturb V3 using single mode
//     iprob = 2 -- Perturb V3 using multiple mode
//     iprob = 3 -- B rotated by "angle" at interface, multimode perturbation
//
// REFERENCE: R. Liska & B. Wendroff, SIAM J. Sci. Comput., 25, 995 (2003)
//======================================================================================

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../bvals/bvals.hpp"
#include "../hydro/eos/eos.hpp"
#include "../hydro/srcterms/srcterms.hpp"
#include "../coordinates/coordinates.hpp"
#include "../utils/utils.hpp"

void ProjectPressureInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                    FaceField &b, int is, int ie, int js, int je, int ks, int ke);
void ProjectPressureOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                    FaceField &b, int is, int ie, int js, int je, int ks, int ke);
void ProjectPressureInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                    FaceField &b, int is, int ie, int js, int je, int ks, int ke);
void ProjectPressureOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                    FaceField &b, int is, int ie, int js, int je, int ks, int ke);

// made global to share with BC functions
static Real gm1;
static Real grav_acc;

//======================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//======================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  if (mesh_size.nx3 == 1) {  // 2D problem
    // Enroll special BCs
    EnrollUserBoundaryFunction(INNER_X2, ProjectPressureInnerX2);
    EnrollUserBoundaryFunction(OUTER_X2, ProjectPressureOuterX2);
  }
  else { // 3D problem
    // Enroll special BCs
    EnrollUserBoundaryFunction(INNER_X3, ProjectPressureInnerX3);
    EnrollUserBoundaryFunction(OUTER_X3, ProjectPressureOuterX3);
  }
  return;
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the Rayleigh-Taylor instability test
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  long int iseed = -1;
  Real gamma = phydro->peos->GetGamma();
  gm1 = gamma - 1.0;
  
  Real kx = 2.0*(PI)/(pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min);
  Real ky = 2.0*(PI)/(pmy_mesh->mesh_size.x2max - pmy_mesh->mesh_size.x2min);
  Real kz = 2.0*(PI)/(pmy_mesh->mesh_size.x3max - pmy_mesh->mesh_size.x3min);

  // Read perturbation amplitude, problem switch, density ratio
  Real amp = pin->GetReal("problem","amp");
  int iprob = pin->GetInteger("problem","iprob");
  Real drat = pin->GetOrAddReal("problem","drat",3.0);


// 2D PROBLEM ---------------------------------------------------------------

  if (block_size.nx3 == 1) {
    grav_acc = phydro->pf_srcterms->GetG2();
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real den=1.0;
        if (pcoord->x2v(j) > 0.0) den *= drat;

        if (iprob == 1) {
          phydro->u(IM2,k,j,i) = (1.0+cos(kx*pcoord->x1v(i)))*(1.0+cos(ky*pcoord->x2v(j)))/4.0;
        } else {
          phydro->u(IM2,k,j,i) = (ran2(&iseed) - 0.5)*(1.0+cos(ky*pcoord->x2v(j)));
        }

        phydro->u(IDN,k,j,i) = den;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) *= (den*amp);
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = (1.0/gamma + grav_acc*den*(pcoord->x2v(j)))/gm1;
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/den;
        }
      }
    }}

    // initialize interface B, same for all iprob
    if (MAGNETIC_FIELDS_ENABLED) {
      // Read magnetic field strength, angle [in degrees, 0 is along +ve X-axis]
      Real b0 = pin->GetReal("problem","b0");
      for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie+1; i++) {
        pfield->b.x1f(k,j,i) = b0;
      }}}
      for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je+1; j++) {
      for (int i=is; i<=ie; i++) {
        pfield->b.x2f(k,j,i) = 0.0;
      }}}
      for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        pfield->b.x3f(k,j,i) = 0.0;
      }}}
      if (NON_BAROTROPIC_EOS) {
        for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          phydro->u(IEN,k,j,i) += 0.5*b0*b0;
        }}}
      }
    }

// 3D PROBLEM ----------------------------------------------------------------

  } else {
    grav_acc = phydro->pf_srcterms->GetG3();
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real den=1.0;
        if (pcoord->x3v(k) > 0.0) den *= drat;

        if (iprob == 1) {
          phydro->u(IM3,k,j,i) = (1.0+cos(kx*(pcoord->x1v(i))))/8.0
                     *(1.0+cos(ky*pcoord->x2v(j)))*(1.0+cos(kz*pcoord->x3v(k)));
        } else {
          phydro->u(IM3,k,j,i) = amp*(ran2(&iseed) - 0.5)*(1.0+cos(kz*pcoord->x3v(k)));
        }

        phydro->u(IDN,k,j,i) = den;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) *= (den*amp);
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = (1.0/gamma + grav_acc*den*(pcoord->x3v(k)))/gm1;
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/den;
        }
      }
    }}

    // initialize interface B
    if (MAGNETIC_FIELDS_ENABLED) {
      // Read magnetic field strength, angle [in degrees, 0 is along +ve X-axis]
      Real b0 = pin->GetReal("problem","b0");
      Real angle = pin->GetReal("problem","angle");
      angle = (angle/180.)*PI;
      for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie+1; i++) {
        if (pcoord->x3v(k) > 0.0) {
          pfield->b.x1f(k,j,i) = b0;
        } else {
          pfield->b.x1f(k,j,i) = b0*cos(angle);
        }
      }}}
      for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je+1; j++) {
      for (int i=is; i<=ie; i++) {
        if (pcoord->x3v(k) > 0.0) {
          pfield->b.x2f(k,j,i) = 0.0;
        } else {
          pfield->b.x2f(k,j,i) = b0*sin(angle);
        }
      }}}
      for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        pfield->b.x3f(k,j,i) = 0.0;
      }}}
      if (NON_BAROTROPIC_EOS) {
        for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          phydro->u(IEN,k,j,i) += 0.5*b0*b0;
        }}}
      }
    }
  } /* end of 3D initialization */

  return;
}


//--------------------------------------------------------------------------------------
//! \fn void ProjectPressureInnerX2()
//  \brief  Pressure is integated into ghost cells to improve hydrostatic eqm

void ProjectPressureInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                            FaceField &b, int is, int ie, int js, int je, int ks, int ke)
{
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
      if (n==(IVY)) {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          a(IVY,k,js-j,i) = -a(IVY,k,js+j-1,i);  // reflect 2-velocity
        }
      } else if (n==(IEN)) {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          a(IEN,k,js-j,i) = a(IEN,k,js+j-1,i) 
             - a(IDN,k,js+j-1,i)*grav_acc*(2*j-1)*pco->dx2f(j)/gm1;
        }
      } else {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          a(n,k,js-j,i) = a(n,k,js+j-1,i);
        }
      }
    }}
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b2
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
      for (int i=is; i<=ie+1; ++i) {
        b.x1f(k,(js-j),i) =  b.x1f(k,(js+j-1),i);
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
      for (int i=is; i<=ie; ++i) {
        b.x2f(k,(js-j),i) = -b.x2f(k,(js+j  ),i);  // reflect 2-field
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
      for (int i=is; i<=ie; ++i) {
        b.x3f(k,(js-j),i) =  b.x3f(k,(js+j-1),i);
      }
    }}
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void ProjectPressureOuterX2()
//  \brief  Pressure is integated into ghost cells to improve hydrostatic eqm

void ProjectPressureOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                            FaceField &b, int is, int ie, int js, int je, int ks, int ke)
{
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
      if (n==(IVY)) {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          a(IVY,k,je+j,i) = -a(IVY,k,je-j+1,i);  // reflect 2-velocity
        }
      } else if (n==(IEN)) {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          a(IEN,k,je+j,i) = a(IEN,k,je-j+1,i) 
             + a(IDN,k,je-j+1,i)*grav_acc*(2*j-1)*pco->dx2f(j)/gm1;
        }
      } else {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          a(n,k,je+j,i) = a(n,k,je-j+1,i);
        }
      }
    }}
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b2
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
      for (int i=is; i<=ie+1; ++i) {
        b.x1f(k,(je+j  ),i) =  b.x1f(k,(je-j+1),i);
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
      for (int i=is; i<=ie; ++i) {
        b.x2f(k,(je+j+1),i) = -b.x2f(k,(je-j+1),i);  // reflect 2-field
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
      for (int i=is; i<=ie; ++i) {
        b.x3f(k,(je+j  ),i) =  b.x3f(k,(je-j+1),i);
      }
    }}
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void ProjectPressureInnerX3()
//  \brief  Pressure is integated into ghost cells to improve hydrostatic eqm

void ProjectPressureInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                            FaceField &b, int is, int ie, int js, int je, int ks, int ke)
{
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
      if (n==(IVZ)) {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          a(IVZ,ks-k,j,i) = -a(IVZ,ks+k-1,j,i);  // reflect 3-vel
        }
      } else if (n==(IEN)) {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          a(IEN,ks-k,j,i) = a(IEN,ks+k-1,j,i) 
             - a(IDN,ks+k-1,j,i)*grav_acc*(2*k-1)*pco->dx3f(k);
        }
      } else {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          a(n,ks-k,j,i) = a(n,ks+k-1,j,i);
        }
      }
    }}
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b3
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=is; i<=ie+1; ++i) {
        b.x1f((ks-k),j,i) =  b.x1f((ks+k-1),j,i);
      }
    }}

    for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma simd
      for (int i=is; i<=ie; ++i) {
        b.x2f((ks-k),j,i) =  b.x2f((ks+k-1),j,i);
      }
    }}

    for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=is; i<=ie; ++i) {
        b.x3f((ks-k),j,i) = -b.x3f((ks+k  ),j,i);  // reflect 3-field
      }
    }}
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void ProjectPressureOuterX3()
//  \brief  Pressure is integated into ghost cells to improve hydrostatic eqm

void ProjectPressureOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                            FaceField &b, int is, int ie, int js, int je, int ks, int ke)
{
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
      if (n==(IVZ)) {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          a(IVZ,ke+k,j,i) = -a(IVZ,ke-k+1,j,i);  // reflect 3-vel
        }
      } else if (n==(IEN)) {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          a(IEN,ke+k,j,i) = a(IEN,ke-k+1,j,i)
             + a(IDN,ke-k+1,j,i)*grav_acc*(2*k-1)*pco->dx3f(k);
        }
      } else {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          a(n,ke+k,j,i) = a(n,ke-k+1,j,i);
        }
      }
    }}
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b3
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=is; i<=ie+1; ++i) {
        b.x1f((ke+k  ),j,i) =  b.x1f((ke-k+1),j,i);
      }
    }}

    for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma simd
      for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
        b.x2f((ke+k  ),j,i) =  b.x2f((ke-k+1),j,i);
      }
    }}

    for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=is; i<=ie; ++i) {
        b.x3f((ke+k+1),j,i) = -b.x3f((ke-k+1),j,i);  // reflect 3-field
      }
    }}
  }

  return;
}
