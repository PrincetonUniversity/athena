//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file puncture_z4c.cpp
//  \brief implementation of functions in the Z4c class for initializing puntures evolution

// C++ standard headers
#include <cmath> // pow

// Athena++ headers
#include "z4c.hpp"
#include "z4c_macro.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"

// twopuncturesc: Stand-alone library ripped from Cactus
#include "TwoPunctures.h"

//----------------------------------------------------------------------------------------
// \!fn void Z4c::ADMTwoPunctures(AthenaArray<Real> & u)
// \brief Initialize ADM vars to two punctures

void Z4c::ADMTwoPunctures(ParameterInput *pin, AthenaArray<Real> & u_adm, ini_data *data)
{
  bool verbose = pin->GetOrAddBoolean("problem", "verbose", 0);
  
  //if(verbose)
  //  Z4c::DebugInfoVars();

  ADM_vars adm;
  SetADMAliases(u_adm, adm);

  MeshBlock * pmb = pmy_block;
  Coordinates * pco = pmb->pcoord;

  // Flat spacetime
  ADMMinkowski(u_adm);

  //--

  // construct initial data set
  if(verbose)
    std::cout << "Generating two puncture data." << std::endl;

  if(verbose)
    std::cout << "Done!" << std::endl;

  // interpolate to ADM variables based on solution
  if(verbose)
    std::cout << "Interpolating current MeshBlock." << std::endl;

  int imin[3] = {0, 0, 0};

  // dimensions of block in each direction
  int n[3] = {(*pmb).block_size.nx1 + 2 * GSIZEI,
              (*pmb).block_size.nx2 + 2 * GSIZEJ,
              (*pmb).block_size.nx3 + 2 * GSIZEK};

  int sz = n[0] * n[1] * n[2];
  // this could be done instead by accessing and casting the Athena vars but
  // then it is coupled to implementation details etc.
  Real *gxx = new Real[sz], *gyy = new Real[sz], *gzz = new Real[sz];
  Real *gxy = new Real[sz], *gxz = new Real[sz], *gyz = new Real[sz];

  Real *Kxx = new Real[sz], *Kyy = new Real[sz], *Kzz = new Real[sz];
  Real *Kxy = new Real[sz], *Kxz = new Real[sz], *Kyz = new Real[sz];

  Real *psi = new Real[sz];
  Real *alp = new Real[sz];

  Real *x = new Real[n[0]];
  Real *y = new Real[n[1]];
  Real *z = new Real[n[2]];

  // need to populate coordinates
  for(int ix_I = 0; ix_I < n[0]; ix_I++){
    x[ix_I] = pco->x1v(ix_I);
  }

  for(int ix_J = 0; ix_J < n[1]; ix_J++){
    y[ix_J] = pco->x2v(ix_J);
  }

  for(int ix_K = 0; ix_K < n[2]; ix_K++){
    z[ix_K] = pco->x3v(ix_K);
  }

  TwoPunctures_Cartesian_interpolation
    (data, // struct containing the previously calculated solution
     imin, // min, max idxs of Cartesian Grid in the three directions
     n,    // <-imax, but this collapses
     n,    // total number of indices in each direction
     x,    // x,         // Cartesian coordinates
     y,    // y,
     z,    // z,
     alp,  // alp,       // lapse
     psi,  // psi,       // conformal factor and derivatives
     NULL, // psix,
     NULL, // psiy,
     NULL, // psiz,
     NULL, // psixx,
     NULL, // psixy,
     NULL, // psixz,
     NULL, // psiyy,
     NULL, // psiyz,
     NULL, // psizz,
     gxx,  // gxx,       // metric components
     gxy,  // gxy,
     gxz,  // gxz,
     gyy,  // gyy,
     gyz,  // gyz,
     gzz,  // gzz,
     Kxx,  // kxx,       // extrinsic curvature components
     Kxy,  // kxy,
     Kxz,  // kxz,
     Kyy,  // kyy,
     Kyz,  // kyz,
     Kzz   // kzz
     );


  int flat_ix;
  double psi4;

  GLOOP3(k,j,i){
    flat_ix = i + n[0]*(j + n[1]*k);

    psi4 = pow(psi[flat_ix], 4);
    adm.psi4(k, j, i) = psi4;

    adm.g_dd(0, 0, k, j, i) = psi4 * gxx[flat_ix];
    adm.g_dd(1, 1, k, j, i) = psi4 * gyy[flat_ix];
    adm.g_dd(2, 2, k, j, i) = psi4 * gzz[flat_ix];
    adm.g_dd(0, 1, k, j, i) = psi4 * gxy[flat_ix];
    adm.g_dd(0, 2, k, j, i) = psi4 * gxz[flat_ix];
    adm.g_dd(1, 2, k, j, i) = psi4 * gyz[flat_ix];

    adm.K_dd(0, 0, k, j, i) = Kxx[flat_ix];
    adm.K_dd(1, 1, k, j, i) = Kyy[flat_ix];
    adm.K_dd(2, 2, k, j, i) = Kzz[flat_ix];
    adm.K_dd(0, 1, k, j, i) = Kxy[flat_ix];
    adm.K_dd(0, 2, k, j, i) = Kxz[flat_ix];
    adm.K_dd(1, 2, k, j, i) = Kyz[flat_ix];

  }


  free(gxx); free(gyy); free(gzz);
  free(gxy); free(gxz); free(gyz);

  free(Kxx); free(Kyy); free(Kzz);
  free(Kxy); free(Kxz); free(Kyz);

  free(psi); free(alp);

  free(x); free(y); free(z);

  if(verbose)
    std::cout << "\n\n<-Z4c::ADMTwoPunctures\n\n";
}

//void Z4c::DebugInfoVars(){
//  // dump some basic info to term
//  printf("\n\n->Z4c::DebugInfoVars\n");
//
//  MeshBlock * pmb = pmy_block;
//  Coordinates * pco = pmb->pcoord;
//
//  printf("\n=Ghost node info:\n");
//  printf("(GSIZEI, GSIZEJ, GSIZEK)=(%d, %d, %d)\n",
//         GSIZEI, GSIZEJ, GSIZEK);
//
//  printf("\n=MeshBlock.block_size [physical nodes]");
//
//  int nxyz[3] = {(*pmb).block_size.nx1,
//                 (*pmb).block_size.nx2,
//                 (*pmb).block_size.nx3};
//
//  printf("(nx1, nx2, nx3)=(%d, %d, %d)\n", nxyz[0], nxyz[1], nxyz[2]);
//  printf("(x1min, x1max)=(%lf, %lf)\n",
//         (*pmb).block_size.x1min,
//         (*pmb).block_size.x1max);
//
//  printf("(x2min, x2max)=(%lf, %lf)\n",
//         (*pmb).block_size.x2min,
//         (*pmb).block_size.x2max);
//
//  printf("(x3min, x3max)=(%lf, %lf)\n",
//         (*pmb).block_size.x3min,
//         (*pmb).block_size.x3max);
//
//
//  printf("\n=Coordinates info [current block with ghosts]\n");
//  printf("(x1min, x1max)=(%lf, %lf)\n",
//         (*pco).x1v(0),
//         (*pco).x1v(nxyz[0] + GSIZEI * 2));
//
//  printf("(x2min, x2max)=(%lf, %lf)\n",
//         (*pco).x2v(0),
//         (*pco).x2v(nxyz[1] + GSIZEJ * 2));
//
//  printf("(x3min, x3max)=(%lf, %lf)\n",
//         (*pco).x3v(0),
//         (*pco).x3v(nxyz[2] + GSIZEK * 2));
//
//
//  printf("\n\n<-Z4c::DebugInfoVars\n");
//
//}
