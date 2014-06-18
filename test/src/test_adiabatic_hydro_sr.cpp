// Tests for variable inversion in adiabatic hydro

// Primary header
#include "test_adiabatic_hydro_sr.hpp"

// gtest headers
#include "gtest/gtest.h"

// Athena headers
#include "../../src/athena.hpp"         // enums, macros
#include "../../src/athena_arrays.hpp"  // array_access

// Test inversion with no y- or z- motion
TEST_F(AdiabaticHydroSRTest1, TestX)
{
  // Prepare inputs
  cons(IDN,0,0,0) = 1.0000015219045548e0;
  cons(IEN,0,0,0) = 3.1000062399724630e1;
  cons(IM1,0,0,0) = 3.4782707762878009e-5;
  cons(IM2,0,0,0) = 0.0;
  cons(IM3,0,0,0) = 0.0;

  // Prepare expected outputs
  prim_expected[IDN] = 1.0000015219031373e+00;
  prim_expected[IEN] = 1.0000020292487640e+01;
  prim_expected[IM1] = 1.6836932068180557e-06;
  prim_expected[IM2] = 0.0;
  prim_expected[IM3] = 0.0;

  // Check conserved-to-primitive inversion
  pfluid->ConservedToPrimitive(cons, prim);
  for (int n = 0; n < NVAR; n++)
    EXPECT_DOUBLE_TOL(prim_expected[n], prim(n,0,0,0));
}
