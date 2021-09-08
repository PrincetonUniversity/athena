//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file quirk.cpp
//! \brief Problem generator for the Quirk test (Quirk 1994)
//! Variant is from Hanawa et al (2008), section 3.2. First of five shock tubes
//!
//! Problem generator for a shock tube with odd-even perturbation.
//! Physically such a perturbation should not grow but with high resolution schemes
//! such as Roe/HLLC/HLLD often suffer from unphysical amplification of the perturbation
//! at strong shocks. This problem is well known as the Carbuncle phenomenon.
//! To suppress it, use the new LHLLC/LHLLD solvers (Minoshima et al. 2021).
//========================================================================================

// C headers

// C++ headers
#include <cmath>      // sqrt()
#include <cstdio>     // fopen(), freopen(), fprintf(), fclose()
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"

namespace {
int ishock;
Real gm;
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//! \brief Calculate the difference of the post-shock entropy in odd and even rows
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  if (Globals::my_rank == 0) {
    MeshBlock *pmb = my_blocks(0);
    AthenaArray<Real> &w = pmb->phydro->w;
    Real dodd  = w(IDN, pmb->ks, pmb->js+1, ishock);
    Real podd  = w(IPR, pmb->ks, pmb->js+1, ishock);
    Real deven = w(IDN, pmb->ks, pmb->js,   ishock);
    Real peven = w(IPR, pmb->ks, pmb->js,   ishock);
    Real sodd  = podd  / std::pow(dodd,  gm);
    Real seven = peven / std::pow(deven, gm);
    Real deltas = std::abs(sodd - seven);
    if (deltas > 0.05)
      std::cout << "The scheme suffers from the Carbuncle phenomenon : delta s = "
                << deltas << std::endl;
    else
      std::cout << "The scheme looks stable against the Carbuncle phenomenon : delta s = "
                << deltas << std::endl;
    // open output file and write out errors
    std::string fname;
    fname.assign("carbuncle-diff.dat");
    std::stringstream msg;
    FILE *pfile;

    // The file exists -- reopen the file in append mode
    if ((pfile = std::fopen(fname.c_str(), "r")) != nullptr) {
      if ((pfile = std::freopen(fname.c_str(), "a", pfile)) == nullptr) {
        msg << "### FATAL ERROR in function Mesh::UserWorkAfterLoop"
            << std::endl << "Error output file could not be opened" <<std::endl;
        ATHENA_ERROR(msg);
      }

      // The file does not exist -- open the file in write mode and add headers
    } else {
      if ((pfile = std::fopen(fname.c_str(), "w")) == nullptr) {
        msg << "### FATAL ERROR in function Mesh::UserWorkAfterLoop"
            << std::endl << "Error output file could not be opened" <<std::endl;
        ATHENA_ERROR(msg);
      }
      std::fprintf(pfile, "# |s_odd - s_even|\n");
    }

    // write difference
    std::fprintf(pfile, "%e\n", deltas);
    std::fclose(pfile);
  }
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Problem Generator for the shock tube tests
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  gm = peos->GetGamma();
  Real igm1 = 1.0 / (gm - 1.0);

  Real xshock = 0.4;
  for (ishock = 0; pcoord->x1v(ishock) < xshock; ++ishock) {}
  ishock--;

  Real dl =  3.692;
  Real ul = -0.625;
  Real pl =  26.85;
  Real dr =  1.0;
  Real ur = -5.0;
  Real pr =  0.6;
  Real dd = dl - 0.135;
  Real ud = ul + 0.219;
  Real pd = pl - 1.31;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        if (i <= ishock) {
          phydro->u(IDN,k,j,i) = dl;
          phydro->u(IM1,k,j,i) = dl*ul;
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = 0.0;
          phydro->u(IEN,k,j,i) = pl*igm1 + 0.5*dl*SQR(ul);
        } else {
          phydro->u(IDN,k,j,i) = dr;
          phydro->u(IM1,k,j,i) = dr*ur;
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = 0.0;
          phydro->u(IEN,k,j,i) = pr*igm1 + 0.5*dr*SQR(ur);
        }
        if (i == ishock && j % 2 == 0) {
          phydro->u(IDN,k,j,i) = dd;
          phydro->u(IM1,k,j,i) = dd*ud;
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = 0.0;
          phydro->u(IEN,k,j,i) = pd*igm1 + 0.5*dd*SQR(ud);
        }
      }
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    Real bx = pin->GetReal("problem", "bx");
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie+1; ++i)
          pfield->b.x1f(k,j,i) = bx;
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
        for (int i=is; i<=ie; ++i)
          pfield->b.x2f(k,j,i) = 0.0;
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i)
          pfield->b.x3f(k,j,i) = 0.0;
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i)
          phydro->u(IEN,k,j,i) += 0.5*SQR(bx);
      }
    }
  }
  return;
}
