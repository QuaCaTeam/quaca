#include "catch.hpp"
#include "Quaca.h"

TEST_CASE("MuTest")
{
  SECTION("Constructor with input file")
  {
    OhmicMemoryKernel *memorykernel = new OhmicMemoryKernel("../data/test_files/MemoryKernel.ini");
    REQUIRE( memorykernel->get_gamma() == 0.69420);
  };

  SECTION("The memory kernel obeys causality (mu(omega)=mu*(-omega*)")
  {
    MemoryKernel *mk = new OhmicMemoryKernel(30.0);
    REQUIRE( mk->mu(1E0) == mk->mu(-1E0));
  };
};
