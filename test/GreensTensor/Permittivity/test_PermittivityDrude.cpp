#include "catch.hpp"
#include "Quaca.h"

TEST_CASE("Drude model can be constructed in different ways", "[PermittivityDrude]")
{
  SECTION("Constructor with input file")
  {
      PermittivityDrude *perm = new PermittivityDrude("../data/test_files/PermittivityDrude.ini");
      REQUIRE(perm->get_gamma() == 3.5E-2);
      REQUIRE(perm->get_omega_p() == 3.2);
  };
  SECTION("Constructor with direct input")
  {
      PermittivityDrude *perm = new PermittivityDrude(3.5E-2, 3.2);
      REQUIRE(perm->get_gamma() == 3.5E-2);
      REQUIRE(perm->get_omega_p() == 3.2);
  };
};

TEST_CASE("Drude model obeys crossing relation", "[PermittivityDrude]")
{
  PermittivityDrude *perm = new PermittivityDrude(3.5E-2, 3.2);
  REQUIRE( perm->epsilon(3.0) == std::conj(perm->epsilon(-3.0)) );
};
