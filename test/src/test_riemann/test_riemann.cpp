// Tests for Riemann solvers

// Primary header
#include "test_riemann.hpp"

// gtest headers
#include "gtest/gtest.h"

// Athena headers
#include "../../../src/athena.hpp"         // enums, macros
#include "../../../src/athena_arrays.hpp"  // array_access

// Shock tube 3 from Mignone & Bodo 2005, MNRAS 364 126 - left side
TEST_F(HLLCSRTest, MB_3_L)
{
  // Prepare left inputs
  prim_left(IDN,0,0,0) = 1.0e1;
  prim_left(IEN,0,0,0) = 1.3333333333333329e1;
  prim_left(IM1,0,0,0) = 0.0;
  prim_left(IM2,0,0,0) = 0.0;
  prim_left(IM3,0,0,0) = 0.0;

  // Prepare right inputs
  prim_right(IDN,0,0,0) = 1.0e1;
  prim_right(IEN,0,0,0) = 1.3333333333333329e1;
  prim_right(IM1,0,0,0) = 0.0;
  prim_right(IM2,0,0,0) = 0.0;
  prim_right(IM3,0,0,0) = 0.0;

  // Prepare expected outputs
  flux_expected[IDN] = 0.0;
  flux_expected[IEN] = -2.5441511085978707e-15;
  flux_expected[IM1] = 1.3333333333333329e1;
  flux_expected[IM2] = 0.0;
  flux_expected[IM3] = 0.0;

  // Check conserved-to-primitive inversion
  pfluid_integrator->RiemannSolver(0, 0, 0, 0, IM1, IM2, IM3, prim_left, prim_right,
      flux);
  for (int n = 0; n < NVAR; n++)
    EXPECT_DOUBLE_TOL(flux_expected[n], flux(n,0,0,0));
}

// Shock tube 3 from Mignone & Bodo 2005, MNRAS 364 126 - right side
TEST_F(HLLCSRTest, MB_3_R)
{
  // Prepare left inputs
  prim_left(IDN,0,0,0) = 1.0e0;
  prim_left(IEN,0,0,0) = 6.6666666661515567e-7;
  prim_left(IM1,0,0,0) = 0.0;
  prim_left(IM2,0,0,0) = 0.0;
  prim_left(IM3,0,0,0) = 0.0;

  // Prepare right inputs
  prim_right(IDN,0,0,0) = 1.0e0;
  prim_right(IEN,0,0,0) = 6.6666666661515567e-7;
  prim_right(IM1,0,0,0) = 0.0;
  prim_right(IM2,0,0,0) = 0.0;
  prim_right(IM3,0,0,0) = 0.0;

  // Prepare expected outputs
  flux_expected[IDN] = 0.0;
  flux_expected[IEN] = 0.0;
  flux_expected[IM1] = 6.6666666661515567e-7;
  flux_expected[IM2] = 0.0;
  flux_expected[IM3] = 0.0;

  // Check conserved-to-primitive inversion
  pfluid_integrator->RiemannSolver(0, 0, 0, 0, IM1, IM2, IM3, prim_left, prim_right,
      flux);
  for (int n = 0; n < NVAR; n++)
    EXPECT_DOUBLE_TOL(flux_expected[n], flux(n,0,0,0));
}

// Shock tube 3 from Mignone & Bodo 2005, MNRAS 364 126 - center
TEST_F(HLLCSRTest, MB_3_C)
{
  // Prepare left inputs
  prim_left(IDN,0,0,0) = 1.0e1;
  prim_left(IEN,0,0,0) = 1.3333333333333329e1;
  prim_left(IM1,0,0,0) = 0.0;
  prim_left(IM2,0,0,0) = 0.0;
  prim_left(IM3,0,0,0) = 0.0;

  // Prepare right inputs
  prim_right(IDN,0,0,0) = 1.0e0;
  prim_right(IEN,0,0,0) = 6.6666666661515567e-7;
  prim_right(IM1,0,0,0) = 0.0;
  prim_right(IM2,0,0,0) = 0.0;
  prim_right(IM3,0,0,0) = 0.0;

  // Prepare expected outputs
  flux_expected[IDN] = 3.1882534747310025e0;
  flux_expected[IEN] = 9.7877309660532692e0;
  flux_expected[IM1] = 6.3241936054311374e0;
  flux_expected[IM2] = 0.0;
  flux_expected[IM3] = 0.0;

  // Check conserved-to-primitive inversion
  pfluid_integrator->RiemannSolver(0, 0, 0, 0, IM1, IM2, IM3, prim_left, prim_right,
      flux);
  for (int n = 0; n < NVAR; n++)
    EXPECT_DOUBLE_TOL(flux_expected[n], flux(n,0,0,0));
}

// Generic Gamma=5/3 test 1 (R*)
TEST_F(HLLCSRTest, Test1)
{
  // Prepare left inputs
  prim_left(IDN,0,0,0) = 1.0e-2;
  prim_left(IEN,0,0,0) = 1.0e-2;
  prim_left(IM1,0,0,0) = 0.0;
  prim_left(IM2,0,0,0) = 0.0;
  prim_left(IM3,0,0,0) = 0.0;

  // Prepare right inputs
  prim_right(IDN,0,0,0) = 1.0e0;
  prim_right(IEN,0,0,0) = 1.0e0;
  prim_right(IM1,0,0,0) = 0.0;
  prim_right(IM2,0,0,0) = 0.0;
  prim_right(IM3,0,0,0) = 0.0;

  // Prepare expected outputs
  flux_expected[IDN] = -3.0737873516350456e-1;
  flux_expected[IEN] = -7.7818413849874979e-1;
  flux_expected[IM1] = 4.6300192719551114e-1;
  flux_expected[IM2] = 0.0;
  flux_expected[IM3] = 0.0;

  // Check conserved-to-primitive inversion
  pfluid_integrator->RiemannSolver(0, 0, 0, 0, IM1, IM2, IM3, prim_left, prim_right,
      flux);
  for (int n = 0; n < NVAR; n++)
    EXPECT_DOUBLE_TOL(flux_expected[n], flux(n,0,0,0));
}

// Generic Gamma=5/3 test 2 (L*)
TEST_F(HLLCSRTest, Test2)
{
  // Prepare left inputs
  prim_left(IDN,0,0,0) = 1.0e-2;
  prim_left(IEN,0,0,0) = 1.0e2;
  prim_left(IM1,0,0,0) = 0.0;
  prim_left(IM2,0,0,0) = 0.0;
  prim_left(IM3,0,0,0) = 0.0;

  // Prepare right inputs
  prim_right(IDN,0,0,0) = 1.0e0;
  prim_right(IEN,0,0,0) = 1.0e0;
  prim_right(IM1,0,0,0) = 0.0;
  prim_right(IM2,0,0,0) = 0.0;
  prim_right(IM3,0,0,0) = 0.0;

  // Prepare expected outputs
  flux_expected[IDN] = 3.8675130206656653e-3;
  flux_expected[IEN] = 6.0433493073291352e1;
  flux_expected[IM1] = 5.0657246377342176e1;
  flux_expected[IM2] = 0.0;
  flux_expected[IM3] = 0.0;

  // Check conserved-to-primitive inversion
  pfluid_integrator->RiemannSolver(0, 0, 0, 0, IM1, IM2, IM3, prim_left, prim_right,
      flux);
  for (int n = 0; n < NVAR; n++)
    EXPECT_DOUBLE_TOL(flux_expected[n], flux(n,0,0,0));
}

// Generic Gamma=5/3 test 3 (R*)
TEST_F(HLLCSRTest, Test3)
{
  // Prepare left inputs
  prim_left(IDN,0,0,0) = 1.0e2;
  prim_left(IEN,0,0,0) = 1.0000000000000430e-2;
  prim_left(IM1,0,0,0) = 0.0;
  prim_left(IM2,0,0,0) = 0.0;
  prim_left(IM3,0,0,0) = 0.0;

  // Prepare right inputs
  prim_right(IDN,0,0,0) = 1.0e0;
  prim_right(IEN,0,0,0) = 1.0e0;
  prim_right(IM1,0,0,0) = 0.0;
  prim_right(IM2,0,0,0) = 0.0;
  prim_right(IM3,0,0,0) = 0.0;

  // Prepare expected outputs
  flux_expected[IDN] = -1.3467281792926080e-2;
  flux_expected[IEN] = -4.6692917494036548e-2;
  flux_expected[IM1] = 9.6777882577215113e-1;
  flux_expected[IM2] = 0.0;
  flux_expected[IM3] = 0.0;

  // Check conserved-to-primitive inversion
  pfluid_integrator->RiemannSolver(0, 0, 0, 0, IM1, IM2, IM3, prim_left, prim_right,
      flux);
  for (int n = 0; n < NVAR; n++)
    EXPECT_DOUBLE_TOL(flux_expected[n], flux(n,0,0,0));
}

// Generic Gamma=5/3 test 4 (L*)
TEST_F(HLLCSRTest, Test4)
{
  // Prepare left inputs
  prim_left(IDN,0,0,0) = 1.0e2;
  prim_left(IEN,0,0,0) = 9.9999999999999986e1;
  prim_left(IM1,0,0,0) = 0.0;
  prim_left(IM2,0,0,0) = 0.0;
  prim_left(IM3,0,0,0) = 0.0;

  // Prepare right inputs
  prim_right(IDN,0,0,0) = 1.0e0;
  prim_right(IEN,0,0,0) = 1.0e0;
  prim_right(IM1,0,0,0) = 0.0;
  prim_right(IM2,0,0,0) = 0.0;
  prim_right(IM3,0,0,0) = 0.0;

  // Prepare expected outputs
  flux_expected[IDN] = 3.0737873516350465e1;
  flux_expected[IEN] = 7.7818413849874972e1;
  flux_expected[IM1] = 4.6300192719551092e1;
  flux_expected[IM2] = 0.0;
  flux_expected[IM3] = 0.0;

  // Check conserved-to-primitive inversion
  pfluid_integrator->RiemannSolver(0, 0, 0, 0, IM1, IM2, IM3, prim_left, prim_right,
      flux);
  for (int n = 0; n < NVAR; n++)
    EXPECT_DOUBLE_TOL(flux_expected[n], flux(n,0,0,0));
}

// Generic Gamma=5/3 test 5 (L*)
TEST_F(HLLCSRTest, Test5)
{
  // Prepare left inputs
  prim_left(IDN,0,0,0) = 1.0e0;
  prim_left(IEN,0,0,0) = 1.0e0;
  prim_left(IM1,0,0,0) = 0.0;
  prim_left(IM2,0,0,0) = 0.0;
  prim_left(IM3,0,0,0) = 0.0;

  // Prepare right inputs
  prim_right(IDN,0,0,0) = 1.0e-2;
  prim_right(IEN,0,0,0) = 1.0e-2;
  prim_right(IM1,0,0,0) = 0.0;
  prim_right(IM2,0,0,0) = 0.0;
  prim_right(IM3,0,0,0) = 0.0;

  // Prepare expected outputs
  flux_expected[IDN] = 3.0737873516350456e-1;
  flux_expected[IEN] = 7.7818413849874979e-1;
  flux_expected[IM1] = 4.6300192719551114e-1;
  flux_expected[IM2] = 0.0;
  flux_expected[IM3] = 0.0;

  // Check conserved-to-primitive inversion
  pfluid_integrator->RiemannSolver(0, 0, 0, 0, IM1, IM2, IM3, prim_left, prim_right,
      flux);
  for (int n = 0; n < NVAR; n++)
    EXPECT_DOUBLE_TOL(flux_expected[n], flux(n,0,0,0));
}

// Generic Gamma=5/3 test 6 (R*)
TEST_F(HLLCSRTest, Test6)
{
  // Prepare left inputs
  prim_left(IDN,0,0,0) = 1.0e0;
  prim_left(IEN,0,0,0) = 1.0e0;
  prim_left(IM1,0,0,0) = 0.0;
  prim_left(IM2,0,0,0) = 0.0;
  prim_left(IM3,0,0,0) = 0.0;

  // Prepare right inputs
  prim_right(IDN,0,0,0) = 1.0e-2;
  prim_right(IEN,0,0,0) = 1.0e2;
  prim_right(IM1,0,0,0) = 0.0;
  prim_right(IM2,0,0,0) = 0.0;
  prim_right(IM3,0,0,0) = 0.0;

  // Prepare expected outputs
  flux_expected[IDN] = -3.8675130206656653e-3;
  flux_expected[IEN] = -6.0433493073291352e1;
  flux_expected[IM1] = 5.0657246377342176e1;
  flux_expected[IM2] = 0.0;
  flux_expected[IM3] = 0.0;

  // Check conserved-to-primitive inversion
  pfluid_integrator->RiemannSolver(0, 0, 0, 0, IM1, IM2, IM3, prim_left, prim_right,
      flux);
  for (int n = 0; n < NVAR; n++)
    EXPECT_DOUBLE_TOL(flux_expected[n], flux(n,0,0,0));
}

// Generic Gamma=5/3 test 7 (L*)
TEST_F(HLLCSRTest, Test7)
{
  // Prepare left inputs
  prim_left(IDN,0,0,0) = 1.0e0;
  prim_left(IEN,0,0,0) = 1.0e0;
  prim_left(IM1,0,0,0) = 0.0;
  prim_left(IM2,0,0,0) = 0.0;
  prim_left(IM3,0,0,0) = 0.0;

  // Prepare right inputs
  prim_right(IDN,0,0,0) = 1.0e2;
  prim_right(IEN,0,0,0) = 1.0000000000000430e-2;
  prim_right(IM1,0,0,0) = 0.0;
  prim_right(IM2,0,0,0) = 0.0;
  prim_right(IM3,0,0,0) = 0.0;

  // Prepare expected outputs
  flux_expected[IDN] = 1.3467281792926080e-2;
  flux_expected[IEN] = 4.6692917494036548e-2;
  flux_expected[IM1] = 9.6777882577215113e-1;
  flux_expected[IM2] = 0.0;
  flux_expected[IM3] = 0.0;

  // Check conserved-to-primitive inversion
  pfluid_integrator->RiemannSolver(0, 0, 0, 0, IM1, IM2, IM3, prim_left, prim_right,
      flux);
  for (int n = 0; n < NVAR; n++)
    EXPECT_DOUBLE_TOL(flux_expected[n], flux(n,0,0,0));
}

// Generic Gamma=5/3 test 8 (R*)
TEST_F(HLLCSRTest, Test8)
{
  // Prepare left inputs
  prim_left(IDN,0,0,0) = 1.0e0;
  prim_left(IEN,0,0,0) = 1.0e0;
  prim_left(IM1,0,0,0) = 0.0;
  prim_left(IM2,0,0,0) = 0.0;
  prim_left(IM3,0,0,0) = 0.0;

  // Prepare right inputs
  prim_right(IDN,0,0,0) = 1.0e2;
  prim_right(IEN,0,0,0) = 9.9999999999999986e1;
  prim_right(IM1,0,0,0) = 0.0;
  prim_right(IM2,0,0,0) = 0.0;
  prim_right(IM3,0,0,0) = 0.0;

  // Prepare expected outputs
  flux_expected[IDN] = -3.0737873516350465e1;
  flux_expected[IEN] = -7.7818413849874972e1;
  flux_expected[IM1] = 4.6300192719551092e1;
  flux_expected[IM2] = 0.0;
  flux_expected[IM3] = 0.0;

  // Check conserved-to-primitive inversion
  pfluid_integrator->RiemannSolver(0, 0, 0, 0, IM1, IM2, IM3, prim_left, prim_right,
      flux);
  for (int n = 0; n < NVAR; n++)
    EXPECT_DOUBLE_TOL(flux_expected[n], flux(n,0,0,0));
}

// Generic Gamma=5/3 test 9 (L)
TEST_F(HLLCSRTest, Test9)
{
  // Prepare left inputs
  prim_left(IDN,0,0,0) = 1.0000000000000537e0;
  prim_left(IEN,0,0,0) = 1.0000000000000888e0;
  prim_left(IM1,0,0,0) = 8.9999999999998870e-1;
  prim_left(IM2,0,0,0) = 0.0;
  prim_left(IM3,0,0,0) = 0.0;

  // Prepare right inputs
  prim_right(IDN,0,0,0) = 2.0000000000000289e0;
  prim_right(IEN,0,0,0) = 2.0000000000000480e0;
  prim_right(IM1,0,0,0) = 7.9999999999999361e-1;
  prim_right(IM2,0,0,0) = 0.0;
  prim_right(IM3,0,0,0) = 0.0;

  // Prepare expected outputs
  flux_expected[IDN] = 2.0647416048350302e0;
  flux_expected[IEN] = 1.6578947368420337e1;
  flux_expected[IM1] = 1.5921052631578204e1;
  flux_expected[IM2] = 0.0;
  flux_expected[IM3] = 0.0;

  // Check conserved-to-primitive inversion
  pfluid_integrator->RiemannSolver(0, 0, 0, 0, IM1, IM2, IM3, prim_left, prim_right,
      flux);
  for (int n = 0; n < NVAR; n++)
    EXPECT_DOUBLE_TOL(flux_expected[n], flux(n,0,0,0));
}

// Generic Gamma=5/3 test 10 (R)
TEST_F(HLLCSRTest, Test10)
{
  // Prepare left inputs
  prim_left(IDN,0,0,0) = 2.0000000000000289e0;
  prim_left(IEN,0,0,0) = 2.0000000000000480e0;
  prim_left(IM1,0,0,0) = -7.9999999999999361e-1;
  prim_left(IM2,0,0,0) = 0.0;
  prim_left(IM3,0,0,0) = 0.0;

  // Prepare right inputs
  prim_right(IDN,0,0,0) = 1.0000000000000537e0;
  prim_right(IEN,0,0,0) = 1.0000000000000888e0;
  prim_right(IM1,0,0,0) = -8.9999999999998870e-1;
  prim_right(IM2,0,0,0) = 0.0;
  prim_right(IM3,0,0,0) = 0.0;

  // Prepare expected outputs
  flux_expected[IDN] = -2.0647416048350302e0;
  flux_expected[IEN] = -1.6578947368420337e1;
  flux_expected[IM1] = 1.5921052631578204e1;
  flux_expected[IM2] = 0.0;
  flux_expected[IM3] = 0.0;

  // Check conserved-to-primitive inversion
  pfluid_integrator->RiemannSolver(0, 0, 0, 0, IM1, IM2, IM3, prim_left, prim_right,
      flux);
  for (int n = 0; n < NVAR; n++)
    EXPECT_DOUBLE_TOL(flux_expected[n], flux(n,0,0,0));
}
