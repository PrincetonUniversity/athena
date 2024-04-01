//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 *
 * This program is free software: you can redistribute and/or modify it under the terms
 * of the GNU General Public License (GPL) as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of GNU GPL in the file LICENSE included in the code
 * distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

// C++ headers
#include <algorithm>  // min
#include <cmath>      // sqrt
#include <fstream>
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../crdiffusion/crdiffusion.hpp"
#include "../crdiffusion/mg_crdiffusion.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"


#if !CRDIFFUSION_ENABLED
#error "The implicit CR diffusion solver must be enabled (-crdiff)."
#endif

namespace {
  Real e0;
}

void CRFixedInnerX1(AthenaArray<Real> &dst, Real time, int nvar,
                    int is, int ie, int js, int je, int ks, int ke, int ngh,
                    const MGCoordinates &coord) {
  for (int n=0; n<nvar; n++) {
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=0; i<ngh; i++)
          dst(n,k,j,is-i-1) = 2.0*e0 - dst(n,k,j,is);
      }
    }
  }
  return;
}

void CRFixedOuterX1(AthenaArray<Real> &dst, Real time, int nvar,
                    int is, int ie, int js, int je, int ks, int ke, int ngh,
                    const MGCoordinates &coord) {
  for (int n=0; n<nvar; n++) {
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=0; i<ngh; i++)
          dst(n,k,j,ie+i+1) = 2.0*e0 - dst(n,k,j,ie);
      }
    }
  }
  return;
}

void CRFixedInnerX2(AthenaArray<Real> &dst, Real time, int nvar,
                    int is, int ie, int js, int je, int ks, int ke, int ngh,
                    const MGCoordinates &coord) {
  for (int n=0; n<nvar; n++) {
    for (int k=ks; k<=ke; k++) {
      for (int j=0; j<ngh; j++) {
        for (int i=is; i<=ie; i++)
          dst(n,k,js-j-1,i) = 2.0*e0 - dst(n,k,js,i);
      }
    }
  }
  return;
}

void CRFixedOuterX2(AthenaArray<Real> &dst, Real time, int nvar,
                    int is, int ie, int js, int je, int ks, int ke, int ngh,
                    const MGCoordinates &coord) {
  for (int n=0; n<nvar; n++) {
    for (int k=ks; k<=ke; k++) {
      for (int j=0; j<ngh; j++) {
        for (int i=is; i<=ie; i++)
          dst(n,k,je+j+1,i) = 2.0*e0 - dst(n,k,je,i);
      }
    }
  }
  return;
}

void CRFixedInnerX3(AthenaArray<Real> &dst, Real time, int nvar,
                    int is, int ie, int js, int je, int ks, int ke, int ngh,
                    const MGCoordinates &coord) {
  for (int n=0; n<nvar; n++) {
    for (int k=0; k<ngh; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++)
          dst(n,ks-k-1,j,i) = 2.0*e0 - dst(n,ks,j,i);
      }
    }
  }
  return;
}

void CRFixedOuterX3(AthenaArray<Real> &dst, Real time, int nvar,
                    int is, int ie, int js, int je, int ks, int ke, int ngh,
                    const MGCoordinates &coord) {
  for (int n=0; n<nvar; n++) {
    for (int k=0; k<ngh; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++)
          dst(n,ke+k+1,j,i) = 2.0*e0 - dst(n,ke,j,i);
      }
    }
  }
  return;
}


int AMRCondition(MeshBlock *pmb) {
  if (pmb->block_size.x1min >= 0.25 && pmb->block_size.x1min <=0.251
  &&  pmb->block_size.x2min >= 0.25 && pmb->block_size.x2min <=0.251
  &&  pmb->block_size.x3min >= 0.25 && pmb->block_size.x3min <=0.251) {
    if (pmb->pmy_mesh->ncycle >= pmb->loc.level - 1)
      return 1;
  }
  return 0;
}


//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  e0 = 1.0;
  EnrollUserRefinementCondition(AMRCondition);
  EnrollUserMGCRDiffusionBoundaryFunction(BoundaryFace::inner_x1, CRFixedInnerX1);
  EnrollUserMGCRDiffusionBoundaryFunction(BoundaryFace::outer_x1, CRFixedOuterX1);
  EnrollUserMGCRDiffusionBoundaryFunction(BoundaryFace::inner_x2, CRFixedInnerX2);
  EnrollUserMGCRDiffusionBoundaryFunction(BoundaryFace::outer_x2, CRFixedOuterX2);
  EnrollUserMGCRDiffusionBoundaryFunction(BoundaryFace::inner_x3, CRFixedInnerX3);
  EnrollUserMGCRDiffusionBoundaryFunction(BoundaryFace::outer_x3, CRFixedOuterX3);
}


//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief cosmic ray diffusion test
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real gamma = peos->GetGamma();
  Real r0 = 0.2, rho0 = 1.0;

  for(int k=ks; k<=ke; ++k) {
    Real x3 = pcoord->x3v(k);
    for (int j=js; j<=je; ++j) {
      Real x2 = pcoord->x2v(j);
      for (int i=is; i<=ie; ++i) {
        Real x1 = pcoord->x1v(i);
        Real r2 = SQR(x1)+SQR(x2)+SQR(x3);
        phydro->u(IDN,k,j,i) = 1.0;//rho0/std::min(r2,r0*r0);
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS)
          phydro->u(IEN,k,j,i) = phydro->u(IDN,k,j,i)/(gamma-1.0);
      }
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          pfield->b.x1f(k,j,i) = 1.0;
        }
      }
    }
    if (block_size.nx2 > 1) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je+1; ++j) {
          for (int i=is; i<=ie; ++i) {
            pfield->b.x2f(k,j,i) = 1.0;
          }
        }
      }
    }
    if (block_size.nx3 > 1) {
      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
    }

    pfield->CalculateCellCenteredField(pfield->b,pfield->bcc,pcoord,is,ie,js,je,ks,ke);

    if (NON_BAROTROPIC_EOS) {
      for(int k=ks; k<=ke; ++k) {
        for(int j=js; j<=je; ++j) {
          for(int i=is; i<=ie; ++i) {
            phydro->u(IEN,k,j,i) +=
              0.5*(SQR((pfield->bcc(IB1,k,j,i)))
                 + SQR((pfield->bcc(IB2,k,j,i)))
                 + SQR((pfield->bcc(IB3,k,j,i))));
          }
        }
      }
    }
  }

  for(int k=ks; k<=ke; ++k) {
    for(int j=js; j<=je; ++j) {
      for(int i=is; i<=ie; ++i)
        pcrdiff->ecr(k,j,i) = 1.0;
    }
  }

  return;
}

