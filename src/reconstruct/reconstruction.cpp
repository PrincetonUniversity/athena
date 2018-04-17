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
#include "../coordinates/coordinates.hpp"

// constructor

Reconstruction::Reconstruction(MeshBlock *pmb, ParameterInput *pin) {
  pmy_block_ = pmb;

  // read and set type of spatial reconstruction
  characteristic_reconstruction = false;
  uniform_limiter[0] = true;
  uniform_limiter[1] = true;
  uniform_limiter[2] = true;
  std::string input_recon = pin->GetOrAddString("time","xorder","2");
  if (input_recon == "1") {
    xorder = 1;
  } else if (input_recon == "2") {
    xorder = 2;
  } else if (input_recon == "2c") {
    xorder = 2;
    characteristic_reconstruction = true;
  } else if ((input_recon == "3") || (input_recon == "4")) {
    xorder = 4;
  } else if ((input_recon == "3c") || (input_recon == "4c")) {
    xorder = 4;
    characteristic_reconstruction = true;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in Reconstruction constructor" << std::endl
        << "xorder=" << input_recon << " not valid choice for reconstruction"<< std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // check that there are the necessary number of ghost zones for PPM
  if (xorder == 4) {
    int req_nghost = 3;
    if (MAGNETIC_FIELDS_ENABLED)
      req_nghost += 1;
    if (NGHOST < req_nghost) {
      std::stringstream msg;
      msg << "### FATAL ERROR in Reconstruction constructor" << std::endl
          << "xorder=" << input_recon <<
          " (PPM) reconstruction selected, but nghost=" << NGHOST << std::endl
          << "Reconfigure with --nghost=XXX with XXX > " << req_nghost-1 << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
  }

  // switch to secondary PLM and PPM limiters for nonuniform and/or curvilinear meshes
  if ((COORDINATE_SYSTEM == "cylindrical") || (COORDINATE_SYSTEM == "spherical_polar")) {
    // curvilinear: all directions, regardless of non/uniformity
    uniform_limiter[0]=false;
    uniform_limiter[1]=false;
    uniform_limiter[2]=false;
  }
  // nonuniform geometric spacing or user-defined MeshGenerator, for all coordinate
  // systems, use nonuniform limiter (non-curvilinear will default to Cartesian factors)
  if (pmb->block_size.x1rat != 1.0)
    uniform_limiter[0]=false;
  if (pmb->block_size.x2rat != 1.0)
    uniform_limiter[1]=false;
  if (pmb->block_size.x3rat != 1.0)
    uniform_limiter[2]=false;
  // uniform cartesian,minkowski,sinusoidal,tilted,schwarzschild,kerr-schild,gr_user
  // will use first PLM/PPM limiter without any coordinate terms

  // Allocate memory for scratch arrays used in PLM and PPM
  int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
  scr01_i_.NewAthenaArray(ncells1);
  scr02_i_.NewAthenaArray(ncells1);

  scr1_ni_.NewAthenaArray(NWAVE,ncells1);
  scr2_ni_.NewAthenaArray(NWAVE,ncells1);
  scr3_ni_.NewAthenaArray(NWAVE,ncells1);
  scr4_ni_.NewAthenaArray(NWAVE,ncells1);

  if (xorder == 4) {
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

    // Precompute PPM coefficients in x1-direction ---------------------------------------
    c1i.NewAthenaArray(ncells1);
    c2i.NewAthenaArray(ncells1);
    c3i.NewAthenaArray(ncells1);
    c4i.NewAthenaArray(ncells1);
    c5i.NewAthenaArray(ncells1);
    c6i.NewAthenaArray(ncells1);
    hplus_ratio_i.NewAthenaArray(ncells1);
    hminus_ratio_i.NewAthenaArray(ncells1);

    // coeffiencients in x1 for uniform Cartesian mesh
    if (uniform_limiter[0]) {
#pragma omp simd
      for (int i=(pmb->is)-(NGHOST); i<=(pmb->ie)+(NGHOST); ++i) {
        c1i(i) = 0.5;
        c2i(i) = 0.5;
        c3i(i) = 0.5;
        c4i(i) = 0.5;
        c5i(i) = 1.0/6.0;
        c6i(i) = -1.0/6.0;
      }

    // coeffcients in x1 for non-uniform or cuvilinear mesh
    // (unnecessary work in case of uniform curvilinear mesh)
    } else {
#pragma omp simd
      for (int i=(pmb->is)-(NGHOST)+1; i<=(pmb->ie)+(NGHOST)-1; ++i) {
        Real& dx_im1 = pmb->pcoord->dx1f(i-1);
        Real& dx_i   = pmb->pcoord->dx1f(i  );
        Real& dx_ip1 = pmb->pcoord->dx1f(i+1);
        Real qe = dx_i/(dx_im1 + dx_i + dx_ip1);       // Outermost coeff in CW eq 1.7
        c1i(i) = qe*(2.0*dx_im1+dx_i)/(dx_ip1 + dx_i); // First term in CW eq 1.7
        c2i(i) = qe*(2.0*dx_ip1+dx_i)/(dx_im1 + dx_i); // Second term in CW eq 1.7

        if (i > (pmb->is)-(NGHOST)+1) {  // c3-c6 are not computed in first iteration
          Real& dx_im2 = pmb->pcoord->dx1f(i-2);
          Real qa = dx_im2 + dx_im1 + dx_i + dx_ip1;
          Real qb = dx_im1/(dx_im1 + dx_i);
          Real qc = (dx_im2 + dx_im1)/(2.0*dx_im1 + dx_i);
          Real qd = (dx_ip1 + dx_i)/(2.0*dx_i + dx_im1);
          qb = qb + 2.0*dx_i*qb/qa*(qc-qd);
          c3i(i) = 1.0 - qb;
          c4i(i) = qb;
          c5i(i) = dx_i/qa*qd;
          c6i(i) = -dx_im1/qa*qc;
        }
      }
      // Compute curvilinear geometric factors for limiter (Mignone eq 48)
      for (int i=(pmb->is)-1; i<=(pmb->ie)+1; ++i) {
        if ((COORDINATE_SYSTEM == "cylindrical") ||
            (COORDINATE_SYSTEM == "spherical_polar")) {
          Real h_plus, h_minus;
          Real& dx_i   = pmb->pcoord->dx1f(i);
          Real& xv_i   = pmb->pcoord->x1v(i);
          if (COORDINATE_SYSTEM == "cylindrical") {
            // cylindrical radial coordinate
            h_plus = 3.0 + dx_i/(2.0*xv_i);
            h_minus = 3.0 - dx_i/(2.0*xv_i);
          } else {
            // spherical radial coordinate
            h_plus = 3.0 + (2.0*dx_i*(10.0*xv_i + dx_i))/(20.0*SQR(xv_i) + SQR(dx_i));
            h_minus = 3.0 + (2.0*dx_i*(-10.0*xv_i + dx_i))/(20.0*SQR(xv_i) + SQR(dx_i));
          }
          hplus_ratio_i(i) = (h_plus + 1.0)/(h_minus - 1.0);
          hminus_ratio_i(i) = (h_minus + 1.0)/(h_plus - 1.0);
        } else { // Cartesian, SR, GR
          // h_plus = 3.0;
          // h_minus = 3.0;
          // Ratios are = 2 for Cartesian coords, as in original PPM overshoot limiter
          hplus_ratio_i(i) = 2.0;
          hminus_ratio_i(i) = 2.0;
        }
      }
    }

    // Precompute PPM coefficients in x2-direction ---------------------------------------
    if (pmb->block_size.nx2 > 1) {
      int ncells2 = pmb->block_size.nx2 + 2*(NGHOST);
      c1j.NewAthenaArray(ncells2);
      c2j.NewAthenaArray(ncells2);
      c3j.NewAthenaArray(ncells2);
      c4j.NewAthenaArray(ncells2);
      c5j.NewAthenaArray(ncells2);
      c6j.NewAthenaArray(ncells2);
      hplus_ratio_j.NewAthenaArray(ncells2);
      hminus_ratio_j.NewAthenaArray(ncells2);

      // coeffiencients in x2 for uniform Cartesian mesh
      if (uniform_limiter[1]) {
#pragma omp simd
        for (int j=(pmb->js)-(NGHOST); j<=(pmb->je)+(NGHOST); ++j) {
          c1j(j) = 0.5;
          c2j(j) = 0.5;
          c3j(j) = 0.5;
          c4j(j) = 0.5;
          c5j(j) = 1.0/6.0;
          c6j(j) = -1.0/6.0;
        }

      // coeffcients in x2 for non-uniform or cuvilinear mesh
      // (unnecessary work in case of uniform curvilinear mesh)
      } else {
#pragma omp simd
        for (int j=(pmb->js)-(NGHOST)+2; j<=(pmb->je)+(NGHOST)-1; ++j) {
          Real& dx_jm1 = pmb->pcoord->dx2f(j-1);
          Real& dx_j   = pmb->pcoord->dx2f(j  );
          Real& dx_jp1 = pmb->pcoord->dx2f(j+1);
          Real qe = dx_j/(dx_jm1 + dx_j + dx_jp1);       // Outermost coeff in CW eq 1.7
          c1j(j) = qe*(2.0*dx_jm1+dx_j)/(dx_jp1 + dx_j); // First term in CW eq 1.7
          c2j(j) = qe*(2.0*dx_jp1+dx_j)/(dx_jm1 + dx_j); // Second term in CW eq 1.7

          if (j > (pmb->js)-(NGHOST)+1) {  // c3-c6 are not computed in first iteration
            Real& dx_jm2 = pmb->pcoord->dx2f(j-2);
            Real qa = dx_jm2 + dx_jm1 + dx_j + dx_jp1;
            Real qb = dx_jm1/(dx_jm1 + dx_j);
            Real qc = (dx_jm2 + dx_jm1)/(2.0*dx_jm1 + dx_j);
            Real qd = (dx_jp1 + dx_j)/(2.0*dx_j + dx_jm1);
            qb = qb + 2.0*dx_j*qb/qa*(qc-qd);
            c3j(j) = 1.0 - qb;
            c4j(j) = qb;
            c5j(j) = dx_j/qa*qd;
            c6j(j) = -dx_jm1/qa*qc;
          }
        }
        // Compute curvilinear geometric factors for limiter (Mignone eq 48)
        for (int j=(pmb->js)-1; j<=(pmb->je)+1; ++j) {
          // corrections to PPMx2 only for spherical polar coordinates
          if (COORDINATE_SYSTEM == "spherical_polar") {
            // x2 = theta polar coordinate adjustment
            Real h_plus, h_minus;
            Real& dx_j   = pmb->pcoord->dx2f(j);
            Real& xf_j   = pmb->pcoord->x2f(j);
            Real& xf_jp1   = pmb->pcoord->x2f(j+1);
            Real dmu = cos(xf_j) - cos(xf_jp1);
            Real dmu_tilde = sin(xf_j) - sin(xf_jp1);
            h_plus = (dx_j*(dmu_tilde + dx_j*cos(xf_jp1)))/(
                dx_j*(sin(xf_j) + sin(xf_jp1)) - 2.0*dmu);
            h_minus = -(dx_j*(dmu_tilde + dx_j*cos(xf_j)))/(
                dx_j*(sin(xf_j) + sin(xf_jp1)) - 2.0*dmu);
            hplus_ratio_j(j) = (h_plus + 1.0)/(h_minus - 1.0);
            hminus_ratio_j(j) = (h_minus + 1.0)/(h_plus - 1.0);
          } else {
            // h_plus = 3.0;
            // h_minus = 3.0;
            // Ratios are both = 2, as in orig PPM overshoot limiter
            hplus_ratio_j(j) = 2.0;
            hminus_ratio_j(j) = 2.0;
          }
        }
      }
    }

    // Precompute PPM coefficients in x3-direction
    if (pmb->block_size.nx3 > 1) {
      int ncells3 = pmb->block_size.nx3 + 2*(NGHOST);
      c1k.NewAthenaArray(ncells3);
      c2k.NewAthenaArray(ncells3);
      c3k.NewAthenaArray(ncells3);
      c4k.NewAthenaArray(ncells3);
      c5k.NewAthenaArray(ncells3);
      c6k.NewAthenaArray(ncells3);
      hplus_ratio_k.NewAthenaArray(ncells3);
      hminus_ratio_k.NewAthenaArray(ncells3);

      // coeffiencients in x3 for uniform Cartesian mesh
      if (uniform_limiter[2]) {
#pragma omp simd
        for (int k=(pmb->ks)-(NGHOST); k<=(pmb->ke)+(NGHOST); ++k) {
          c1k(k) = 0.5;
          c2k(k) = 0.5;
          c3k(k) = 0.5;
          c4k(k) = 0.5;
          c5k(k) = 1.0/6.0;
          c6k(k) = -1.0/6.0;
        }

      // coeffcients in x3 for non-uniform or cuvilinear mesh
      // (unnecessary work in case of uniform curvilinear mesh)
      } else {
#pragma omp simd
        for (int k=(pmb->ks)-(NGHOST)+2; k<=(pmb->ke)+(NGHOST)-1; ++k) {
          Real& dx_km1 = pmb->pcoord->dx3f(k-1);
          Real& dx_k   = pmb->pcoord->dx3f(k  );
          Real& dx_kp1 = pmb->pcoord->dx3f(k+1);
          Real qe = dx_k/(dx_km1 + dx_k + dx_kp1);       // Outermost coeff in CW eq 1.7
          c1k(k) = qe*(2.0*dx_km1+dx_k)/(dx_kp1 + dx_k); // First term in CW eq 1.7
          c2k(k) = qe*(2.0*dx_kp1+dx_k)/(dx_km1 + dx_k); // Second term in CW eq 1.7

          if (k > (pmb->ks)-(NGHOST)+1) {  // c3-c6 are not computed in first iteration
            Real& dx_km2 = pmb->pcoord->dx3f(k-2);
            Real qa = dx_km2 + dx_km1 + dx_k + dx_kp1;
            Real qb = dx_km1/(dx_km1 + dx_k);
            Real qc = (dx_km2 + dx_km1)/(2.0*dx_km1 + dx_k);
            Real qd = (dx_kp1 + dx_k)/(2.0*dx_k + dx_km1);
            qb = qb + 2.0*dx_k*qb/qa*(qc-qd);
            c3k(k) = 1.0 - qb;
            c4k(k) = qb;
            c5k(k) = dx_k/qa*qd;
            c6k(k) = -dx_km1/qa*qc;
          }
        }
        // Compute curvilinear geometric factors for limiter (Mignone eq 48)
        // No corrections in x3 for the built-in Newtonian coordinate systems
        for (int k=(pmb->ks)-1; k<=(pmb->ke)+1; ++k) {
          // h_plus = 3.0;
          // h_minus = 3.0;
          // Ratios are both = 2 for Cartesian and all curviliniear coords
          hplus_ratio_k(k) = 2.0;
          hminus_ratio_k(k) = 2.0;
        }
      }
    }
  }

}

// destructor

Reconstruction::~Reconstruction() {
  scr01_i_.DeleteAthenaArray();
  scr02_i_.DeleteAthenaArray();

  scr1_ni_.DeleteAthenaArray();
  scr2_ni_.DeleteAthenaArray();
  scr3_ni_.DeleteAthenaArray();
  scr4_ni_.DeleteAthenaArray();

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

  c1i.DeleteAthenaArray();
  c2i.DeleteAthenaArray();
  c3i.DeleteAthenaArray();
  c4i.DeleteAthenaArray();
  c5i.DeleteAthenaArray();
  c6i.DeleteAthenaArray();
  hplus_ratio_i.DeleteAthenaArray();
  hminus_ratio_i.DeleteAthenaArray();

  c1j.DeleteAthenaArray();
  c2j.DeleteAthenaArray();
  c3j.DeleteAthenaArray();
  c4j.DeleteAthenaArray();
  c5j.DeleteAthenaArray();
  c6j.DeleteAthenaArray();
  hplus_ratio_j.DeleteAthenaArray();
  hminus_ratio_j.DeleteAthenaArray();

  c1k.DeleteAthenaArray();
  c2k.DeleteAthenaArray();
  c3k.DeleteAthenaArray();
  c4k.DeleteAthenaArray();
  c5k.DeleteAthenaArray();
  c6k.DeleteAthenaArray();
  hplus_ratio_k.DeleteAthenaArray();
  hminus_ratio_k.DeleteAthenaArray();
}
