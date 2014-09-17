// Tests for Riemann solvers

// Primary header
#include "test_riemann_schwarzschild_adiabatic_hydro.hpp"

// gtest headers
#include "gtest/gtest.h"

// Athena headers
#include "../../../src/athena.hpp"         // enums, macros
#include "../../../src/athena_arrays.hpp"  // array_access

// No flow across equator
TEST_F(Geodesic2DTest, EquatorialFlux)
{
  // Prepare expected inputs and outputs
  for (int i = is; i <= ie; i++)
  {
    // Left primitives
    prim_left(IDN,i) = 1.0;
    prim_left(IEN,i) = 1.0;
    prim_left(IM1,i) = 0.0;
    prim_left(IM2,i) = 0.0;
    prim_left(IM3,i) = 0.0;

    // Right primitives
    prim_right(IDN,i) = 1.0;
    prim_right(IEN,i) = 1.0;
    prim_right(IM1,i) = 0.0;
    prim_right(IM2,i) = 0.0;
    prim_right(IM3,i) = 0.0;

    // Expected fluxes
    flux_expected(IDN,i) = 0.0;
    flux_expected(IEN,i) = 0.0;
    flux_expected(IM1,i) = 0.0;
    flux_expected(IM2,i) = 0.0;
    flux_expected(IM3,i) = 0.0;
  }

  // Check for correct fluxes
  pfluid_integrator->RiemannSolver(ks, jm, is, ie, IM1, IM2, IM3, &prim_left,
      &prim_right, &flux);
  for (int n = 0; n < NVAR; n++)
    EXPECT_DOUBLE_TOL(flux_expected(n,0), flux(n,0));
}
