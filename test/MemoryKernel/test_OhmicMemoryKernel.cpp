#include "Quaca.h"
#include "catch.hpp"

TEST_CASE("Ohmic memory kernel can be counstructed in different ways",
          "[OhmicMemoryKernel]") {
  SECTION("Constructor with input file") {
    OhmicMemoryKernel *memorykernel =
        new OhmicMemoryKernel("../data/test_files/MemoryKernel.ini");
    REQUIRE(memorykernel->get_gamma() == 0.69420);
  };

  SECTION("Constructor with direct input") {
    OhmicMemoryKernel *memorykernel = new OhmicMemoryKernel(7.8);
    REQUIRE(memorykernel->get_gamma() == 7.8);
  };
};

TEST_CASE("Ohmic memory kernel obeys crossing relation",
          "[OhmicMemoryKernel]") {
  OhmicMemoryKernel *mk = new OhmicMemoryKernel(30.0);
  REQUIRE(mk->mu(1.0) == std::conj(mk->mu(-1.0)));
};
