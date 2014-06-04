// Conserved-to-primitive inversion for adiabatic hydrodynamics in special relativity

// TODO: make conserved inputs const

// Primary header
#include "../fluid.hpp"

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // atan2(), cbrt(), cos(), pow(), sqrt()

// Athena headers
#include "../athena.hpp"         // array access, macros
#include "../athena_arrays.hpp"  // AthenaArray
#include "../mesh.hpp"           // Block

// Variable inverter
// Inputs:
//   cons: conserved quantities
// Outputs:
//   prim: primitives
// Notes:
//   solves quartic equation for velocity explicitly
//   follows relativistic hydro routine in old Athena's convert_var.c:
//     1) given quartic |v|^4 + a3 |v|^3 + a2 |v|^2 + a1 |v| + a0
//     2) construct resolvent cubic x^3 + b2 x^2 + b1 x + b0:
//          b2 = -a2
//          b1 = -4 a0 + a1 a3
//          b0 = 4 a0 a2 - a1^2 - a0 a3^2
//     3) eliminate quadratic term from cubic y^3 + (3 c1) y - 2 c2 = 0:
//          c1 = b1/3 - b2^2/9
//          c2 = -(1/2) b0 + (1/6) b1 b2 - (1/27) b2^3
//          c3 = c1^3 + c2^2
//     4) find real root of new cubic:
//          if c3 >= 0, y0 = (c2 + sqrt(c3))^(1/3) + (c2 - sqrt(c3))^(1/3)
//          otherwise use y0 = 2 (c2^2 + c3)^(1/6) cos((1/3) atan2(sqrt(-c2), c3))
//          formulas are equivalent except for implicit assumptions about branch cuts
//     5) find real root of original (resolvent) cubic:
//          x0 = y0 - b2/3
//     6) solve for (correct) root of original quartic:
//          d1 = 1/2 * (a3 + sqrt(4 x0 - 4 a2 + a3^2))
//          d0 = 1/2 * (x0 - sqrt(x0^2 - 4 a0))
//          then |v|^2 + d1 |v| + d0 = 0
//          |v| = 1/2 * (-d1 + sqrt(d1^2 - 4 d0))
void Fluid::ConservedToPrimitive(AthenaArray<Real> &cons, AthenaArray<Real> &prim)
{
  // Parameters
  const Real max_velocity = 1.0 - 1.0e-15;

  // Extract ratio of specific heats
  const Real gamma_adi = GetGamma();
  const Real gamma_adi_minus_1 = gamma_adi - 1.0;

  // Determine array bounds
  Block *pb = pparent_block;
  int is = pb->is;
  int ie = pb->ie;
  int jl = pb->js;
  int ju = pb->je;
  int kl = pb->ks;
  int ku = pb->ke;
  if (pb->block_size.nx2 > 1)
  {
    jl -= (NGHOST);
    ju += (NGHOST);
  }
  if (pb->block_size.nx3 > 1)
  {
    kl -= (NGHOST);
    ku += (NGHOST);
  }

  // Make array copies for performance reasons
  AthenaArray<Real> cons_copy = cons.ShallowCopy();
  AthenaArray<Real> prim_copy = prim.ShallowCopy();

  // Go through cells
  for (int k = kl; k <= ku; k++)
    for (int j = jl; j <= ju; j++)
    {
#pragma simd
      for (int i = is - (NGHOST); i <= ie + (NGHOST); i++)
      {
        // Extract conserved quantities
        Real &d = cons_copy(IDN,k,j,i);
        Real &e = cons_copy(IEN,k,j,i);
        Real &mx = cons_copy(IVX,k,j,i);
        Real &my = cons_copy(IVY,k,j,i);
        Real &mz = cons_copy(IVZ,k,j,i);

        // Extract primitives
        Real &rho = prim_copy(IDN,k,j,i);
        Real &pgas = prim_copy(IEN,k,j,i);
        Real &vx = prim_copy(IVX,k,j,i);
        Real &vy = prim_copy(IVY,k,j,i);
        Real &vz = prim_copy(IVZ,k,j,i);

        // Calculate total momentum
        Real m_sq = mx*mx + my*my + mz*mz;

        // Case out based on whether momentum vanishes
        if (m_sq >= TINY_NUMBER*TINY_NUMBER)  // generic case, nonzero velocity
        {
          // Step 1: Prepare quartic coefficients
          Real m_abs = sqrt(m_sq);
          Real denom_inverse = 1.0 / (gamma_adi_minus_1*gamma_adi_minus_1
              * (d*d + m_sq));
          Real a3 = -2.0 * gamma_adi * gamma_adi_minus_1 * m_abs * e * denom_inverse;
          Real a2 = (gamma_adi*gamma_adi * e*e + 2.0 * gamma_adi_minus_1 * m_sq -
              gamma_adi_minus_1*gamma_adi_minus_1 * d*d) * denom_inverse;
          Real a1 = -2.0 * gamma_adi * m_abs * e * denom_inverse;
          Real a0 = m_sq * denom_inverse;

          // Step 2: Find resolvent cubic coefficients
          Real b2 = -a2;
          Real b1 = -4.0*a0 + a1*a3;
          Real b0 = 4.0*a0*a2 - a1*a1 - a0*a3*a3;

          // Step 3: Eliminate quadratic term from cubic
          Real c1 = b1/3.0 - b2*b2/9.0;
          Real c2 = -b0/2.0 + b1*b2/6.0 - b2*b2*b2/27.0;
          Real c3 = c1*c1*c1 + c2*c2;

          // Step 4: Find real root of new cubic
          Real y0;
          if (c3 >= 0.0)
            y0 = cbrt(c2 + sqrt(c3)) + cbrt(c2 - sqrt(c3));
          else
            y0 = 2.0 * pow(c2*c2 + c3, 1.0/6.0) * cos(atan2(sqrt(-c3), c2) / 3.0);

          // Step 5: Find real root of original (resolvent) cubic:
          Real x0 = y0 - b2/3.0;

          // Step 6: Solve for (correct) root of original quartic
          Real d1 = (a3 + sqrt(4.0*x0 - 4.0*a2 + a3*a3)) / 2.0;
          Real d0 = (x0 - sqrt(x0*x0 - 4.0*a0)) / 2.0;
          Real v_abs = (-d1 + sqrt(d1*d1 - 4.0*d0)) / 2.0;

          // Ensure velocity is physical
          v_abs = std::max(v_abs, 0.0);
          v_abs = std::min(v_abs, 1.0-1.0e-15);

          // Set primitives
          rho = d * sqrt(1.0 - v_abs*v_abs);
          pgas = gamma_adi_minus_1 * (e - (mx * vx + my * vy + mz * vz) - rho);
          vx = mx * v_abs / m_abs;
          vy = my * v_abs / m_abs;
          vz = mz * v_abs / m_abs;
        }
        else  // vanishing velocity
        {
          rho = d;
          pgas = gamma_adi_minus_1 * (e - rho);
          vx = 0.0;
          vy = 0.0;
          vz = 0.0;
        }
      }
    }
  return;
}
