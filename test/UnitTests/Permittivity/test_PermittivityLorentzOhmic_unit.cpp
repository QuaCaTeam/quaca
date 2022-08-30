#include "Quaca.h"
#include "catch.hpp"

TEST_CASE("LorentzOhmic permittivity constructors work as expected",
          "[PermittivityLorentzOhmic]") {
  SECTION("Direct constructor") {
    double eps_inf = 1.4;
    double alpha_0 = 6e-9;
    double omega_0 = 3.4;
    double gamma = 2.8;

    PermittivityLorentzOhmic perm(eps_inf, alpha_0, omega_0, gamma);

    REQUIRE(perm.get_eps_inf() == eps_inf);
    REQUIRE(perm.get_alpha_0() == alpha_0);
    REQUIRE(perm.get_omega_0() == omega_0);
    REQUIRE(perm.get_gamma() == gamma);
  };
};

TEST_CASE("LorentzOhmic permittivity obeys crossing relation",
          "[PermittivityLorentzOhmic]") {
  double eps_inf = 1.4;
  double alpha_0 = 6e-9;
  double omega_0 = 3.4;
  double gamma = 2.8;

  PermittivityLorentzOhmic perm(eps_inf, alpha_0, omega_0, gamma);

  auto omega = GENERATE(-1.2,0.1);
  REQUIRE(perm.calculate(omega) == std::conj(perm.calculate(-omega)));
};
