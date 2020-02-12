#include "Quaca.h"
#include "catch.hpp"
#include <armadillo>
#include <complex>
#include <iomanip> // std::setprecision

TEST_CASE("Construction of Green's tensor plate vacuum works properly",
          "[GreensTensorPlateVacuum]") {

  SECTION("Construction with direct input") {

    double omega_p = 9;
    double gamma = 0.1;
    PermittivityDrude perm(omega_p, gamma);
    ReflectionCoefficientsLocBulk refl(&perm);

    double v = 1E-5;
    double za = 0.1;
    double beta = 1E4;
    double delta_cut = 20;
    vec::fixed<2> rel_err = {1E-8, 1E-6};

    GreensTensorPlateVacuum Greens(v, za, beta, &refl, delta_cut, rel_err);

    REQUIRE(Greens.get_v() == v);
    REQUIRE(Greens.get_beta() == beta);
    REQUIRE(Greens.get_za() == za);
    REQUIRE(Greens.get_delta_cut() == delta_cut);
    REQUIRE(Greens.get_rel_err_0() == rel_err(0));
    REQUIRE(Greens.get_rel_err_1() == rel_err(1));

    REQUIRE(Greens.get_vacuums_greens_tensor()->get_v() == v);
    REQUIRE(Greens.get_vacuums_greens_tensor()->get_beta() == beta);
  };
};
