#include "Quaca.h"
#include "catch.hpp"

TEST_CASE("LorentzNoBath permittivity constructors work as expected",
          "[PermittivityLorentzNoBath]") {
  SECTION("Direct constructor") {
    PermittivityLorentzNoBath perm(
        "../data/test_files/PermittivityLorentzNoBath.ini");
    REQUIRE(perm.get_eps_inf() == 1.4);
    REQUIRE(perm.get_alpha_zero() == 6e-9);
    REQUIRE(perm.get_omega_0() == 3.4);
  };
  SECTION("ini file constructor") {
    auto eps_inf = GENERATE(take(3, random(0.0, 1e3)));
    auto alpha_zero = GENERATE(take(3, random(0.0, 1e3)));
    auto omega_0 = GENERATE(take(3, random(0.0, 1e3)));

    PermittivityLorentzNoBath perm(eps_inf, alpha_zero, omega_0);
    REQUIRE(perm.get_eps_inf() == eps_inf);
    REQUIRE(perm.get_alpha_zero() == alpha_zero);
    REQUIRE(perm.get_omega_0() == omega_0);
  };
};

TEST_CASE("LorentzNoBath permittivity obeys crossing relation",
          "[PermittivityLorentzNoBath]") {
  auto omega_0 = GENERATE(take(2, random(0.0, 1e3)));
  auto eps_inf = GENERATE(take(2, random(0.0, 1e3)));
  auto alpha_zero = GENERATE(take(2, random(0.0, 1e3)));

  PermittivityLorentzNoBath perm(eps_inf, alpha_zero, omega_0);

  auto omega = GENERATE(take(10, random(-150.4, 150.4)));
  REQUIRE(perm.epsilon(omega) == std::conj(perm.epsilon(-omega)));
};
