#include "Quaca.h"
#include "catch.hpp"

TEST_CASE("Ohmic memory kernel constructors work as expected",
          "[OhmicMemoryKernel]") {

  SECTION("Direct constructor") {
    OhmicMemoryKernel memorykernel(7.8);
    REQUIRE(memorykernel.get_gamma() == 7.8);
  };

  SECTION("ini file constructor") {
    OhmicMemoryKernel memorykernel("../data/test_files/MemoryKernel.ini");
    REQUIRE(memorykernel.get_gamma() == 0.69420);
  };
};

TEST_CASE("Ohmic memory kernel obeys crossing relation",
          "[OhmicMemoryKernel]") {
  OhmicMemoryKernel mk(30.0);
  REQUIRE(mk.mu(1.0) == std::conj(mk.mu(-1.0)));
};
