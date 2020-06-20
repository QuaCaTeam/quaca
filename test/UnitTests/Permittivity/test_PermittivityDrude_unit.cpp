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
    auto omega_p = GENERATE(0.12,3.21,43.12,103.2);
    auto gamma = GENERATE(1.21,32.12,143.21);

    PermittivityDrude perm(omega_p, gamma);
    REQUIRE(perm.get_gamma() == gamma);
    REQUIRE(perm.get_omega_p() == omega_p);
  }
}

TEST_CASE("Drude permittivity obeys crossing relation", "[PermittivityDrude]") {
  PermittivityDrude perm(3.5E-2, 3.2);

  auto omega = GENERATE(-213.21,-65.34,-2.32,0.021,9.87,89.32);
  REQUIRE(perm.calculate(omega) == std::conj(perm.calculate(-omega)));
  REQUIRE(perm.calculate_times_omega(omega) == std::conj(perm.calculate_times_omega(-omega)));
};
