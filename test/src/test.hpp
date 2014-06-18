// Test suite for Athena++

#ifndef TEST_HPP
#define TEST_HPP

// Test headers
#include "gtest/gtest.h"

// Macros

namespace {

// General test class
class GeneralTest : public testing::Test
{
 protected:

  // Tolerances
  double absolute_tol;
  double relative_tol;

  // Constructor
  GeneralTest(double absolute_tol, double relative_tol)
    : absolute_tol(absolute_tol),
      relative_tol(relative_tol)
  {};

  // Destructor
  ~GeneralTest() {};

  // Function acting like test macros
  void EXPECT_DOUBLE_TOL(double expected, double actual)
  {
    double abs_error = relative_tol * expected;
    abs_error = (abs_error < absolute_tol) ? absolute_tol : abs_error;
    EXPECT_NEAR(expected, actual, abs_error);
  }
};

}

#endif  // TEST_HPP
