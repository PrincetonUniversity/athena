//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file rotating_column.cpp
//! \brief Steady, differentially rotating column threaded by an axial magnetic field.
//!
//! An exact smooth steady state of ideal MHD in which v is NOT parallel to B, so that
//! the EMF E = -v x B is a genuinely nonlinear product of two spatially varying fields.
//! In cylindrical coordinates (r, phi, z) about the chosen axis,
//!
//!   rho = d0,   v = Omega(r) r phi_hat,   B = B_z(r) z_hat,
//!   Omega(r) = omega (1 - r^2/R^2)^n,   B_z(r) = b0 + db (1 - r^2/R^2)^n   for r < R,
//!
//! and Omega = 0, B_z = b0 for r >= R. The pressure that balances the centrifugal and
//! magnetic pressure forces, p'(r) = rho Omega^2 r - (B_z^2/2)', is
//!
//!   p(r) = pinf - rho omega^2 R^2/(2(2n+1)) (1 - r^2/R^2)^(2n+1) - (B_z^2 - b0^2)/2.
//!
//! Every quantity is a function of r only and v is azimuthal, so v.grad(.) = 0 for all
//! of them; E = -Omega r B_z r_hat is curl-free, and the state is exactly stationary.
//! The profiles have compact support, so the periodic box is irrelevant as long as the
//! disk r < R fits inside it; (1 - r^2/R^2)^n is C^(n-1) at r = R.
//!
//! The problem is 2D in the plane perpendicular to the axis. Run with several cells
//! along the axis (nx3 > 1 for axis = 3) to exercise the 3D CT code path on a problem
//! whose edge-averaged EMFs do not factor into products of edge averages: the x1-edge
//! EMF E1 = -v2 B3 has <v2 B3> != <v2><B3> along the edge, and curl of the difference
//! does not vanish (it is 4th-order accurate for E3 corners in 2D, where the corner
//! states are point values). Since the exact solution is the initial condition, the
//! error is simply u(t) - u(0), e.g. from two vtk dumps.
//!
//! Input parameters (block <problem>):
//!   axis   : 1, 2 or 3 (default 3): the axis of the column
//!   d0     : density (1.0)
//!   omega  : central angular velocity (8.0)
//!   radius : R (0.5)
//!   nexp   : n, the exponent of the profiles (6)
//!   b0     : axial field outside the column (1.0)
//!   db     : axial field perturbation at the axis (-0.5)
//!   pinf   : pressure outside the column (1.5)
//!
//! Face-averaged B is set with 6-point tensor Gauss-Legendre quadrature over each face;
//! the conserved variables are set to point values (use time/correct_ic = true).
//========================================================================================

// C headers

// C++ headers
#include <cmath>      // sqrt(), pow()
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

#if !MAGNETIC_FIELDS_ENABLED
#error "This problem generator requires magnetic fields"
#endif
#if !NON_BAROTROPIC_EOS
#error "This problem generator requires an adiabatic EOS"
#endif

namespace {
Real d0, omega, rad, b0, db, pinf;
int nexp;

inline Real Shape(Real r2) {  // (1 - r^2/R^2)^n for r < R, else 0
  Real q = 1.0 - r2/SQR(rad);
  return (q > 0.0) ? std::pow(q, nexp) : 0.0;
}
inline Real Omega(Real r2) { return omega*Shape(r2); }
inline Real Bz(Real r2) { return b0 + db*Shape(r2); }
inline Real Pres(Real r2) {
  Real q = 1.0 - r2/SQR(rad);
  Real prot = (q > 0.0) ? d0*SQR(omega)*SQR(rad)/(2.0*(2*nexp + 1))
              *std::pow(q, 2*nexp + 1) : 0.0;
  return pinf - prot - 0.5*(SQR(Bz(r2)) - SQR(b0));
}

// 6-point Gauss-Legendre nodes and weights on [-1, 1]
const int NQ = 6;
const Real gl_x[NQ] = {-0.9324695142031521, -0.6612093864662645, -0.2386191860831969,
                       0.2386191860831969, 0.6612093864662645, 0.9324695142031521};
const Real gl_w[NQ] = {0.1713244923791704, 0.3607615730481386, 0.4679139345726910,
                       0.4679139345726910, 0.3607615730481386, 0.1713244923791704};

// average of Bz(r^2) over the rectangle [xa - ha/2, xa + ha/2] x [xb - hb/2, xb + hb/2]
Real FaceAveBz(Real xa, Real ha, Real xb, Real hb) {
  Real sum = 0.0;
  for (int a=0; a<NQ; ++a) {
    Real pa = xa + 0.5*ha*gl_x[a];
    for (int b=0; b<NQ; ++b) {
      Real pb = xb + 0.5*hb*gl_x[b];
      sum += gl_w[a]*gl_w[b]*Bz(SQR(pa) + SQR(pb));
    }
  }
  return 0.25*sum;
}
} // namespace

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief rotating magnetized column
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  int axis = pin->GetOrAddInteger("problem", "axis", 3);
  d0 = pin->GetOrAddReal("problem", "d0", 1.0);
  omega = pin->GetOrAddReal("problem", "omega", 8.0);
  rad = pin->GetOrAddReal("problem", "radius", 0.5);
  nexp = pin->GetOrAddInteger("problem", "nexp", 6);
  b0 = pin->GetOrAddReal("problem", "b0", 1.0);
  db = pin->GetOrAddReal("problem", "db", -0.5);
  pinf = pin->GetOrAddReal("problem", "pinf", 1.5);
  Real gm1 = peos->GetGamma() - 1.0;

  if (axis < 1 || axis > 3) {
    std::stringstream msg;
    msg << "### FATAL ERROR in rotating_column.cpp ProblemGenerator" << std::endl
        << "problem/axis must be 1, 2 or 3" << std::endl;
    ATHENA_ERROR(msg);
  }
  if (Pres(0.0) <= 0.0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in rotating_column.cpp ProblemGenerator" << std::endl
        << "pressure on the axis is " << Pres(0.0) << "; increase problem/pinf"
        << std::endl;
    ATHENA_ERROR(msg);
  }

  // in-plane coordinates (xa, xb) are the cyclic successors of the axis:
  // axis 3: (x1, x2), axis 1: (x2, x3), axis 2: (x3, x1); v = Omega (-xb, xa)
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real x1 = pcoord->x1v(i), x2 = pcoord->x2v(j), x3 = pcoord->x3v(k);
        Real xa, xb;
        if (axis == 3) {
          xa = x1; xb = x2;
        } else if (axis == 1) {
          xa = x2; xb = x3;
        } else {
          xa = x3; xb = x1;
        }
        Real r2 = SQR(xa) + SQR(xb);
        Real om = Omega(r2);
        Real va = -om*xb, vb = om*xa;
        Real m1, m2, m3;
        if (axis == 3) {
          m1 = d0*va; m2 = d0*vb; m3 = 0.0;
        } else if (axis == 1) {
          m2 = d0*va; m3 = d0*vb; m1 = 0.0;
        } else {
          m3 = d0*va; m1 = d0*vb; m2 = 0.0;
        }
        phydro->u(IDN,k,j,i) = d0;
        phydro->u(IM1,k,j,i) = m1;
        phydro->u(IM2,k,j,i) = m2;
        phydro->u(IM3,k,j,i) = m3;
        // point values at the cell center, including the point value of B_z
        // (consistent with time/correct_ic, which converts u to cell averages)
        phydro->u(IEN,k,j,i) = Pres(r2)/gm1 + 0.5*d0*(SQR(va) + SQR(vb))
                               + 0.5*SQR(Bz(r2));
      }
    }
  }

  // face-averaged field: only the axial component is nonzero
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie+1; i++) {
        pfield->b.x1f(k,j,i) = (axis == 1) ?
            FaceAveBz(pcoord->x2v(j), pcoord->dx2f(j), pcoord->x3v(k), pcoord->dx3f(k))
            : 0.0;
      }
    }
  }
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je+1; j++) {
      for (int i=is; i<=ie; i++) {
        pfield->b.x2f(k,j,i) = (axis == 2) ?
            FaceAveBz(pcoord->x3v(k), pcoord->dx3f(k), pcoord->x1v(i), pcoord->dx1f(i))
            : 0.0;
      }
    }
  }
  for (int k=ks; k<=ke+1; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        pfield->b.x3f(k,j,i) = (axis == 3) ?
            FaceAveBz(pcoord->x1v(i), pcoord->dx1f(i), pcoord->x2v(j), pcoord->dx2f(j))
            : 0.0;
      }
    }
  }
  return;
}
