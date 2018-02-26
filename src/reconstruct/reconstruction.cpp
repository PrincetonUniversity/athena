//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file reconstruction.cpp
//  \brief

// C/C++ headers
#include <sstream>
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "reconstruction.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../eos/eos.hpp"

// constructor

Reconstruction::Reconstruction(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block_ = pmb;
//  Coordinates* = pmb->pcoord;
  int is=pmy_block_->is; int ie=pmy_block_->ie;

  // read and set type of spatial reconstruction
  characteristic_reconstruction = false;
  std::string input_recon = pin->GetOrAddString("time","xorder","2");
  if (input_recon == "1") {
    xorder = 1;
  } else if (input_recon == "2") {
    xorder = 2;
  } else if (input_recon == "2c") {
    xorder = 2;
    characteristic_reconstruction = true;
  } else if (input_recon == "3") {
    xorder = 3;
  } else if (input_recon == "3c") {
    xorder = 3;
    characteristic_reconstruction = true;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in Reconstruction constructor" << std::endl
        << "xorder=" << input_recon << " not valid choice for reconstruction"<< std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // Allocate memory for scratch arrays used in PLM and PPM
  int ncells1 = (ie-is+1) + 2*(NGHOST);
  scr01_i_.NewAthenaArray(ncells1);
  scr02_i_.NewAthenaArray(ncells1);

  scr1_ni_.NewAthenaArray(NWAVE,ncells1);
  scr2_ni_.NewAthenaArray(NWAVE,ncells1);
  scr3_ni_.NewAthenaArray(NWAVE,ncells1);
  scr4_ni_.NewAthenaArray(NWAVE,ncells1);

  if (xorder == 3){
    scr03_i_.NewAthenaArray(ncells1);
    scr04_i_.NewAthenaArray(ncells1);
    scr05_i_.NewAthenaArray(ncells1);
    scr06_i_.NewAthenaArray(ncells1);
    scr07_i_.NewAthenaArray(ncells1);
    scr08_i_.NewAthenaArray(ncells1);
    scr09_i_.NewAthenaArray(ncells1);
    scr10_i_.NewAthenaArray(ncells1);
    scr11_i_.NewAthenaArray(ncells1);
    scr12_i_.NewAthenaArray(ncells1);
    scr13_i_.NewAthenaArray(ncells1);
    scr14_i_.NewAthenaArray(ncells1);

    scr5_ni_.NewAthenaArray(NWAVE,ncells1);
    scr6_ni_.NewAthenaArray(NWAVE,ncells1);
    scr7_ni_.NewAthenaArray(NWAVE,ncells1);
    scr8_ni_.NewAthenaArray(NWAVE,ncells1);
  }

  // Precompute PPM coefficients for each direction
  if (xorder == 3){
    c1i.NewAthenaArray(ncells1);
    c2i.NewAthenaArray(ncells1);
    c3i.NewAthenaArray(ncells1);
    c4i.NewAthenaArray(ncells1);
    c5i.NewAthenaArray(ncells1);
    c6i.NewAthenaArray(ncells1);

    // coeffiencients in x1 for uniform mesh
    if (pmb->block_size.x1rat == 1.0) {
#pragma simd
      for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i){
        c1i(i) = 0.5;
        c2i(i) = 0.5;
        c3i(i) = 0.5;
        c4i(i) = 0.5;
        c5i(i) = 1.0/6.0;
        c6i(i) = -1.0/6.0;
      }

    // coeffcients in x1 for non-uniform mesh
    } else {
#pragma simd
      for (int i=is-(NGHOST)+2; i<=ie+(NGHOST)-1; ++i){
        Real& dx_im2 = pmb->pcoord->dx1f(i-2);
        Real& dx_im1 = pmb->pcoord->dx1f(i-1);
        Real& dx_i   = pmb->pcoord->dx1f(i  );
        Real& dx_ip1 = pmb->pcoord->dx1f(i+1);
        Real qa = dx_im2 + dx_im1 + dx_i + dx_ip1;
        Real qb = dx_im1/(dx_im1 + dx_i);
        Real qc = (dx_im2 + dx_im1)/(2.0*dx_im1 + dx_i);
        Real qd = (dx_ip1 + dx_i)/(2.0*dx_i + dx_im1);
        qb = qb + 2.0*dx_i*qb/qa*(qc-qd);
        Real qe = dx_i/(dx_im1 + dx_i + dx_ip1); // Outermost coefficient in CW eq 1.7
        c1i(i) = qe*(2.0*dx_im1+dx_i)/(dx_ip1 + dx_i); // First term in CW eq 1.7 bracket
        c2i(i) = qe*(2.0*dx_ip1+dx_i)/(dx_im1 + dx_i); // Second term in CW eq 1.7
        c3i(i) = 1.0 - qb;
        c4i(i) = qb;
        c5i(i) = dx_i/qa*qd;
        c6i(i) = -dx_im1/qa*qc;
      }
      // remainder loop: only c1, c2 needed for is-(NGHOST)+1
      Real& dx_im1 = pmb->pcoord->dx1f(is-(NGHOST)  );
      Real& dx_i   = pmb->pcoord->dx1f(is-(NGHOST)+1);
      Real& dx_ip1 = pmb->pcoord->dx1f(is-(NGHOST)+2);
      Real qe = dx_i/(dx_im1 + dx_i + dx_ip1);
      c1i(is-(NGHOST)+1) = qe*(2.0*dx_im1+dx_i)/(dx_ip1 + dx_i);
      c2i(is-(NGHOST)+1) = qe*(2.0*dx_ip1+dx_i)/(dx_im1 + dx_i);
    }

    if (pmb->block_size.x2rat != 1.0) {
      int js = pmb->js;
      int je = pmb->je;
      int ncells2 = pmb->block_size.nx2 + 2*(NGHOST);
      Real qa, qb, qc, qd, qe;
    }
    if (pmb->block_size.x3rat != 1.0) {
      int ks = pmb->ks;
      int ke = pmb->ke;
      int ncells3 = pmb->block_size.nx3 + 2*(NGHOST);
      Real qa, qb, qc, qd, qe;
    }

  }

}

// destructor

Reconstruction::~Reconstruction()
{
  scr01_i_.DeleteAthenaArray();
  scr02_i_.DeleteAthenaArray();

  scr1_ni_.DeleteAthenaArray();
  scr2_ni_.DeleteAthenaArray();
  scr3_ni_.DeleteAthenaArray();
  scr4_ni_.DeleteAthenaArray();

  if (xorder == 3){
    scr03_i_.DeleteAthenaArray();
    scr04_i_.DeleteAthenaArray();
    scr05_i_.DeleteAthenaArray();
    scr06_i_.DeleteAthenaArray();
    scr07_i_.DeleteAthenaArray();
    scr08_i_.DeleteAthenaArray();
    scr09_i_.DeleteAthenaArray();
    scr10_i_.DeleteAthenaArray();
    scr11_i_.DeleteAthenaArray();
    scr12_i_.DeleteAthenaArray();
    scr13_i_.DeleteAthenaArray();
    scr14_i_.DeleteAthenaArray();

    scr5_ni_.DeleteAthenaArray();
    scr6_ni_.DeleteAthenaArray();
    scr7_ni_.DeleteAthenaArray();
    scr8_ni_.DeleteAthenaArray();
  }
}
