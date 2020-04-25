#include "Quaca.h"
#include "catch.hpp"

TEST_CASE("Drude permittivity constructors work as expected",
          "[PermittivityDrude]") {
  SECTION("Direct constructor") {
    PermittivityDrude perm("../data/test_files/PermittivityDrude.json");
    REQUIRE(perm.get_gamma() == 3.5E-2);
    REQUIRE(perm.get_omega_p() == 3.2);
  }
  SECTION("json file constructor") {
    auto omega_p = GENERATE(take(5, random(0.0, 1e3)));
    auto gamma = GENERATE(take(5, random(0.0, 1e3)));

    PermittivityDrude perm(omega_p, gamma);
    REQUIRE(perm.get_gamma() == gamma);
    REQUIRE(perm.get_omega_p() == omega_p);
  }
}

TEST_CASE("Drude permittivity obeys crossing relation", "[PermittivityDrude]") {
  PermittivityDrude perm(3.5E-2, 3.2);

  auto omega = GENERATE(take(10, random(-150.4, 150.4)));
  REQUIRE(perm.epsilon(omega) == std::conj(perm.epsilon(-omega)));
}
