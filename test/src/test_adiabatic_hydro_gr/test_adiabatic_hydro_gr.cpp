// Tests for variable inversion in adiabatic hydro

// Primary header
#include "test_adiabatic_hydro_gr.hpp"

// gtest headers
#include "gtest/gtest.h"

// Athena headers
#include "../../../src/athena.hpp"         // enums, macros, Real
#include "../../../src/athena_arrays.hpp"  // array_access

// Test inversion with no motion in all directions
TEST_F(AdiabaticHydroGRTest1, TestAll)
{
  // Prepare expected outputs
  prim_expected[IDN] = 1.0;
  prim_expected[IEN] = 1.2;
  prim_expected[IM1] = 0.1;
  prim_expected[IM2] = 0.2;
  prim_expected[IM3] = 0.3;

  // Prepare inputs
  set_primitives();
  primitive_to_conserved();

  // Check conserved-to-primitive inversion
  pfluid->ConservedToPrimitive(cons, prim);
  for (int n = 0; n < NVAR; n++)
    EXPECT_DOUBLE_TOL(prim_expected[n], prim(n,0,0,0));
}

// Test inversion with no y- or z-motion
TEST_F(AdiabaticHydroGRTest2, TestX)
{
  // Prepare expected outputs
  prim_expected[IDN] = 1.0;
  prim_expected[IEN] = 1000.0;
  prim_expected[IM1] = 0.0;
  prim_expected[IM2] = 0.0;
  prim_expected[IM3] = 0.0;

  // Prepare inputs
  set_primitives();
  primitive_to_conserved();

  // Check conserved-to-primitive inversion
  pfluid->ConservedToPrimitive(cons, prim);
  for (int n = 0; n < NVAR; n++)
    EXPECT_DOUBLE_TOL(prim_expected[n], prim(n,0,0,0));
}

// Test inversion with imperfect guess
TEST_F(AdiabaticHydroGRTest2, TestImperfect)
{
  // Prepare expected outputs
  prim_expected[IDN] = 0.79482790185472818;
  prim_expected[IEN] = 41.324042099645951;
  prim_expected[IM1] = 0.60683490871340828;
  prim_expected[IM2] = 0.0;
  prim_expected[IM3] = 0.0;

  // Prepare inputs
  primitive_to_conserved();
  prim(IDN,0,0,0) = 1.0;
  prim(IEN,0,0,0) = 10.0;//0.01;
  prim(IM1,0,0,0) = 0.0;
  prim(IM2,0,0,0) = 0.0;
  prim(IM3,0,0,0) = 0.0;

  // Check conserved-to-primitive inversion
  pfluid->ConservedToPrimitive(cons, prim);
  for (int n = 0; n < NVAR; n++)
    EXPECT_DOUBLE_TOL(prim_expected[n], prim(n,0,0,0));
}

