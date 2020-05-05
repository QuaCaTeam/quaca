#include "Quaca.h"
#include "catch.hpp"

TEST_CASE("Ohmic memory kernel constructors work as expected",
          "[OhmicMemoryKernel]") {

  SECTION("Direct constructor") {
    OhmicMemoryKernel memorykernel(7.8);
    REQUIRE(memorykernel.get_gamma() == 7.8);
  }

  SECTION("json file constructor") {
    OhmicMemoryKernel memorykernel("../data/test_files/MemoryKernel.json");
    REQUIRE(memorykernel.get_gamma() == 0.69420);
  }
}

TEST_CASE("Ohmic memory kernel obeys crossing relation",
          "[OhmicMemoryKernel]") {
  OhmicMemoryKernel mk(30.0);
  REQUIRE(mk.calculate(1.0) == std::conj(mk.calculate(-1.0)));
}
