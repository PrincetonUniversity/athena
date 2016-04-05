// Conserved-to-primitive inversion for adiabatic hydrodynamics in special relativity

// TODO: make inputs const?

// Primary header
#include "eos.hpp"

// C++ headers
#include <cmath>   // atan2(), cbrt(), cos(), sqrt()
#include <cfloat>  // FLT_MIN

// Athena headers
#include "../hydro.hpp"               // Hydro
#include "../../athena.hpp"           // enums, macros, Real
#include "../../athena_arrays.hpp"    // AthenaArray
#include "../../mesh.hpp"             // MeshBlock
#include "../../parameter_input.hpp"  // GetReal()
#include "../../field/field.hpp"      // FaceField

// Constructor
// Inputs:
//   pf: pointer to hydro object
//   pin: pointer to runtime inputs
HydroEqnOfState::HydroEqnOfState(Hydro *pf, ParameterInput *pin)
{
  pmy_hydro_ = pf;
  gamma_ = pin->GetReal("hydro", "gamma");
  density_floor_ = pin->GetOrAddReal("hydro", "dfloor", 1024*FLT_MIN);
  pressure_floor_ = pin->GetOrAddReal("hydro", "pfloor", 1024*FLT_MIN);
  gamma_max_ = pin->GetOrAddReal("hydro", "gamma_max", 1000.0);
}

// Destructor
HydroEqnOfState::~HydroEqnOfState() {}

// Variable inverter
// Inputs:
//   cons: conserved quantities
//   prim_old: primitive quantities from previous half timestep (not used)
//   bb: face-centered magnetic field (not used)
//   pco: pointer to Coordinates
//   is,ie,js,je,ks,ke: index bounds of region to be updated
// Outputs:
//   prim: primitives
//   bb_cc: cell-centered magnetic field (not used)
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
void HydroEqnOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
  const AthenaArray<Real> &prim_old, const FaceField &bb, AthenaArray<Real> &prim,
  AthenaArray<Real> &bb_cc, Coordinates *pco, int is, int ie, int js, int je, int ks,
  int ke)
{
  // Parameters
  const Real max_velocity = std::sqrt(1.0 - 1.0/SQR(gamma_max_));

  // Extract ratio of specific heats
  const Real gamma_adi = GetGamma();
  const Real gamma_adi_minus_1 = gamma_adi - 1.0;

  // Make array copies for performance reasons
  AthenaArray<Real> cons_copy,prim_copy;
  cons_copy.InitWithShallowCopy(cons);
  prim_copy.InitWithShallowCopy(prim);

  // Go through cells
  for (int k = ks; k <= ke; k++)
    for (int j = js; j <= je; j++)
    {
      #pragma simd
      for (int i = is; i <= ie; i++)
      {
        // Extract conserved quantities
        Real &d = cons_copy(IDN,k,j,i);
        Real &e = cons_copy(IEN,k,j,i);
        const Real &mx = cons_copy(IVX,k,j,i);
        const Real &my = cons_copy(IVY,k,j,i);
        const Real &mz = cons_copy(IVZ,k,j,i);

        // Extract primitives
        Real &rho = prim_copy(IDN,k,j,i);
        Real &pgas = prim_copy(IEN,k,j,i);
        Real &vx = prim_copy(IVX,k,j,i);
        Real &vy = prim_copy(IVY,k,j,i);
        Real &vz = prim_copy(IVZ,k,j,i);

        // Calculate total momentum
        Real m_sq = SQR(mx) + SQR(my) + SQR(mz);

        // Case out based on whether momentum vanishes
        if (m_sq >= SQR(TINY_NUMBER))  // generic case, nonzero velocity
        {
          // Step 1: Prepare quartic coefficients
          Real m_abs = std::sqrt(m_sq);
          Real denom_inverse = 1.0 / (SQR(gamma_adi_minus_1) * (SQR(d) + m_sq));
          Real a3 = -2.0 * gamma_adi * gamma_adi_minus_1 * m_abs * e * denom_inverse;
          Real a2 = (SQR(gamma_adi) * SQR(e) + 2.0 * gamma_adi_minus_1 * m_sq -
              SQR(gamma_adi_minus_1) * SQR(d)) * denom_inverse;
          Real a1 = -2.0 * gamma_adi * m_abs * e * denom_inverse;
          Real a0 = m_sq * denom_inverse;

          // Step 2: Find resolvent cubic coefficients
          Real b2 = -a2;
          Real b1 = -4.0*a0 + a1*a3;
          Real b0 = 4.0*a0*a2 - SQR(a1) - a0*SQR(a3);

          // Step 3: Eliminate quadratic term from cubic
          Real c1 = b1/3.0 - SQR(b2)/9.0;
          Real c2 = -b0/2.0 + b1*b2/6.0 - SQR(b2)*b2/27.0;
          Real c3 = SQR(c1)*c1 + SQR(c2);

          // Step 4: Find real root of new cubic
          Real y0;
          if (c3 >= 0.0)
            y0 = cbrt(c2 + std::sqrt(c3)) + cbrt(c2 - std::sqrt(c3));
          else
            y0 = 2.0 * cbrt(SQR(c2) + c3)
                * std::cos(std::atan2(std::sqrt(-c3), c2) / 3.0);

          // Step 5: Find real root of original (resolvent) cubic:
          Real x0 = y0 - b2/3.0;

          // Step 6: Solve for (correct) root of original quartic
          Real d1 = (a3 + std::sqrt(4.0*x0 - 4.0*a2 + SQR(a3))) / 2.0;
          Real d0 = (x0 - std::sqrt(SQR(x0) - 4.0*a0)) / 2.0;
          Real v_abs = (-d1 + std::sqrt(SQR(d1) - 4.0*d0)) / 2.0;

          // Ensure velocity is physical
          v_abs = (v_abs > 0.0) ? v_abs : 0.0;                    // sets NaN to 0
          v_abs = (v_abs < max_velocity) ? v_abs : max_velocity;

          // Set density, correcting only conserved density if floor applied
          Real gamma_rel = 1.0 / std::sqrt(1.0 - SQR(v_abs));
          rho = d / gamma_rel;
          if (rho < density_floor_)
          {
            rho = density_floor_;
            d = gamma_rel * rho;
          }

          // Set velocity
          vx = mx * v_abs / m_abs;
          vy = my * v_abs / m_abs;
          vz = mz * v_abs / m_abs;

          // Set pressure, correcting only energy if floor applied
          pgas = gamma_adi_minus_1 * (e - (mx*vx + my*vy + mz*vz) - rho);
          if (pgas < pressure_floor_)
          {
            pgas = pressure_floor_;
            e = pgas/gamma_adi_minus_1 + (mx*vx + my*vy + mz*vz) + rho;
          }
        }
        else  // vanishing velocity
        {
          // Set density, correcting only conserved density if floor applied
          rho = d;
          if (rho < density_floor_)
          {
            rho = density_floor_;
            d = rho;
          }

          // Set pressure, correcting only energy if floor applied
          pgas = gamma_adi_minus_1 * (e - rho);
          if (pgas < pressure_floor_)
          {
            pgas = pressure_floor_;
            e = pgas/gamma_adi_minus_1 + rho;
          }

          // Set velocity
          vx = 0.0;
          vy = 0.0;
          vz = 0.0;
        }
      }
    }
  return;
}

// Function for converting all primitives to conserved variables
// Inputs:
//   prim: primitives
//   bb_cc: cell-centered magnetic field (unused)
//   pco: pointer to Coordinates
//   is,ie,js,je,ks,ke: index bounds of region to be updated
// Outputs:
//   cons: conserved variables
// Notes:
//   single-cell function exists for other purposes; call made to that function rather
//       than having duplicate code
void HydroEqnOfState::PrimitiveToConserved(const AthenaArray<Real> &prim,
     const AthenaArray<Real> &bb_cc, AthenaArray<Real> &cons, Coordinates *pco, int is,
     int ie, int js, int je, int ks, int ke)
{
  // Calculate reduced ratio of specific heats
  Real gamma_adi_red = gamma_/(gamma_-1.0);

  // Go through all cells
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
    {
      #pragma simd
      for (int i = is; i <= ie; ++i)
      {
        // Extract primitives
        const Real &rho = prim(IDN,k,j,i);
        const Real &pgas = prim(IEN,k,j,i);
        const Real &v1 = prim(IVX,k,j,i);
        const Real &v2 = prim(IVY,k,j,i);
        const Real &v3 = prim(IVZ,k,j,i);

        // Calculate 4-velocity
        Real u0 = 1.0 / std::sqrt(1.0 - SQR(v1) - SQR(v2) - SQR(v3));
        Real u1 = u0 * v1;
        Real u2 = u0 * v2;
        Real u3 = u0 * v3;

        // Extract conserved quantities
        Real &d = cons(IDN,k,j,i);
        Real &e = cons(IEN,k,j,i);
        Real &m1 = cons(IM1,k,j,i);
        Real &m2 = cons(IM2,k,j,i);
        Real &m3 = cons(IM3,k,j,i);

        // Set conserved quantities
        Real wgas = rho + gamma_adi_red * pgas;
        d = rho * u0;
        e = wgas * u0 * u0 - pgas;
        m1 = wgas * u0 * u1;
        m2 = wgas * u0 * u2;
        m3 = wgas * u0 * u3;
      }
    }
  return;
}

// Function for calculating relativistic sound speeds
// Inputs:
//   rho_h: enthalpy per unit volume
//   pgas: gas pressure
//   vx: 3-velocity component v^x
//   gamma_lorentz_sq: Lorentz factor \gamma^2
// Outputs:
//   plambda_plus: value set to most positive wavespeed
//   plambda_minus: value set to most negative wavespeed
// Notes:
//   same function as in adiabatic_hydro_gr.cpp
//   references Mignone & Bodo 2005, MNRAS 364 126 (MB)
void HydroEqnOfState::SoundSpeedsSR(
    Real rho_h, Real pgas, Real vx, Real gamma_lorentz_sq,
    Real *plambda_plus, Real *plambda_minus)
{
  const Real gamma_adi = gamma_;
  Real cs_sq = gamma_adi * pgas / rho_h;                                 // (MB 4)
  Real sigma_s = cs_sq / (gamma_lorentz_sq * (1.0-cs_sq));
  Real relative_speed = std::sqrt(sigma_s * (1.0 + sigma_s - SQR(vx)));
  *plambda_plus = 1.0/(1.0+sigma_s) * (vx + relative_speed);             // (MB 23)
  *plambda_minus = 1.0/(1.0+sigma_s) * (vx - relative_speed);            // (MB 23)
  return;
}
