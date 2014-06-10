// Test suite for Athena++

// Primary header
#include "test.hpp"

// gtest headers
#include "gtest/gtest.h"

// Main function
int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
