#include "catch.hpp"
#include "Quaca.h"

TEST_CASE("MuTest")
{
  SECTION("The memory kernel obeys causality (mu(omega)=mu*(-omega*)")
  {
    MemoryKernel *mk = new OhmicMemoryKernel(30.0);
    REQUIRE( mk->mu(1E0) == mk->mu(-1E0));
  };
};
