#include "Quaca.h"
#include "catch.hpp"

TEST_CASE("Lorentz permittivity constructors work as expected",
          "[PermittivityLorentz]") {
  SECTION("Direct constructor") {
    PermittivityLorentz perm("../data/test_files/PermittivityLorentz.json");
    REQUIRE(perm.get_eps_inf() == 1.4);
    REQUIRE(perm.get_alpha_zero() == 6e-9);
    REQUIRE(perm.get_omega_0() == 3.4);
    REQUIRE(perm.get_memory_kernel()->mu(3.5) == 0.69420);
  }
  SECTION("json file constructor") {

    auto mu = std::make_shared<OhmicMemoryKernel>(0.69420);
    auto eps_inf = GENERATE(take(3, random(0.0, 1e3)));
    auto alpha_zero = GENERATE(take(3, random(0.0, 1e3)));
    auto omega_0 = GENERATE(take(3, random(0.0, 1e3)));

    PermittivityLorentz perm(eps_inf, alpha_zero, omega_0, mu);
    REQUIRE(perm.get_eps_inf() == eps_inf);
    REQUIRE(perm.get_alpha_zero() == alpha_zero);
    REQUIRE(perm.get_omega_0() == omega_0);
    REQUIRE(perm.get_memory_kernel()->mu(3.5) == mu->mu(3.5));
  }
}

TEST_CASE("Lorentz permittivity obeys crossing relation",
          "[PermittivityLorentz]") {
  auto omega_0 = GENERATE(take(2, random(0.0, 1e3)));
  auto eps_inf = GENERATE(take(2, random(0.0, 1e3)));
  auto alpha_zero = GENERATE(take(2, random(0.0, 1e3)));
  auto mu = std::make_shared<OhmicMemoryKernel>(0.69420);

  PermittivityLorentz perm(eps_inf, alpha_zero, omega_0, mu);

  auto omega = GENERATE(take(10, random(-150.4, 150.4)));
  REQUIRE(perm.epsilon(omega) == std::conj(perm.epsilon(-omega)));
}
