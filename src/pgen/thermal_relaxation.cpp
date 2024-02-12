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

// C headers

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
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../mesh/mesh.hpp"
#include "../nr_radiation/integrators/rad_integrators.hpp"
#include "../nr_radiation/radiation.hpp"
#include "../parameter_input.hpp"


int RefinementCondition(MeshBlock *pmb);


//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real tgas, er, sigma;

  er = pin->GetOrAddReal("problem","er",10.0);
  tgas = pin->GetOrAddReal("problem","tgas",1.0);
  sigma = pin->GetOrAddReal("problem","sigma",100.0);
  Real gamma = peos->GetGamma();

  // Initialize hydro variable
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        phydro->u(IDN,k,j,i) = 1.0;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = tgas/(gamma-1.0);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        }
      }
    }
  }
  // Now initialize opacity and specific intensity
  if (NR_RADIATION_ENABLED || IM_RADIATION_ENABLED) {
    int nfreq = pnrrad->nfreq;
    //int nang = pnrrad->nang;
    AthenaArray<Real> ir_cm;
    ir_cm.NewAthenaArray(pnrrad->n_fre_ang);
    Real *ir_lab;
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          ir_lab = &(pnrrad->ir(k,j,i,0));
          for (int n=0; n<pnrrad->n_fre_ang; n++) {
             ir_lab[n] = er;
          }
        }
      }
    }

    for (int k=0; k<ncells3; ++k) {
      for (int j=0; j<ncells2; ++j) {
        for (int i=0; i<ncells1; ++i) {
          for (int ifr=0; ifr < nfreq; ++ifr) {
            pnrrad->sigma_s(k,j,i,ifr) = 0.0;
            pnrrad->sigma_a(k,j,i,ifr) = sigma;
            pnrrad->sigma_pe(k,j,i,ifr) = sigma;
            pnrrad->sigma_p(k,j,i,ifr) = sigma;
          }
        }
      }
    }
  }
  return;
}

void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (adaptive) {
    EnrollUserRefinementCondition(RefinementCondition);
  }
}


int RefinementCondition(MeshBlock *pmb) {
  Coordinates *pco = pmb->pcoord;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je;
  int ks=pmb->ks; // ke=pmb->ke;

  Real xmin = pco->x1f(is);
  Real xmax = pco->x1f(ie+1);
  Real ymin = pco->x2f(js);
  Real ymax = pco->x2f(je+1);
  Real tgas = pmb->phydro->w(IPR,ks,js,is)/pmb->phydro->w(IDN,ks,js,is);
  if (xmin > 0.25 && xmax < 0.75 && ymin > 0.25 && ymax < 0.75 && tgas > 50.0) {
    return 1;
  }
  if (xmin > 0.25 && xmax < 0.75 && ymin > 0.25 && ymax < 0.75 && tgas < 50.0) {
    return -1;
  }
  return 0;
}
