#include "Quaca.h"
#include "catch.hpp"
#include <armadillo>
#include <complex>

TEST_CASE("Construction of Green's tensor plate vacuum works properly",
          "[GreensTensorPlateVacuum]") {

  SECTION("Construction with direct input") {

    double omega_p = 9;
    double gamma = 0.1;
    auto perm = std::make_shared<PermittivityDrude>(omega_p, gamma);
    auto refl = std::make_shared<ReflectionCoefficientsLocBulk>(perm);

    double v = 1E-5;
    double za = 0.1;
    double beta = 1E4;
    double delta_cut = 20;
    vec::fixed<2> rel_err = {1E-8, 1E-6};

    GreensTensorPlateVacuum Greens(v, beta, za, refl, delta_cut, rel_err);

    REQUIRE(Greens.get_v() == v);
    REQUIRE(Greens.get_beta() == beta);
    REQUIRE(Greens.get_za() == za);
    REQUIRE(Greens.get_delta_cut() == delta_cut);
    REQUIRE(Greens.get_rel_err_0() == rel_err(0));
    REQUIRE(Greens.get_rel_err_1() == rel_err(1));

    REQUIRE(Greens.get_vacuums_greens_tensor()->get_v() == v);
    REQUIRE(Greens.get_vacuums_greens_tensor()->get_beta() == beta);
  }
}

TEST_CASE("GreensTensorPlateVacuum calculate_tensor function returns the sum "
          "of plate and vacuum Green's tensor",
          "[GreensTensorPlateVacuum]") {

  double omega_p = 9;
  double gamma = 0.1;
  auto perm = std::make_shared<PermittivityDrude>(omega_p, gamma);
  auto refl = std::make_shared<ReflectionCoefficientsLocBulk>(perm);

  double v = 1E-5;
  double za = 0.1;
  double beta = 1E4;
  double delta_cut = 20;
  vec::fixed<2> rel_err = {1E-8, 1E-6};

  GreensTensorPlateVacuum GTPlateVacuum(v, beta, za, refl, delta_cut, rel_err);
  GreensTensorPlate GTPlate(v, beta, za, refl, delta_cut, rel_err);
  GreensTensorVacuum GTVacuum(v, beta, rel_err(0));

  auto k_x = GENERATE(take(1, random(-1e3, 1e3)));
  auto k_y = GENERATE(take(1, random(-1e3, 1e3)));
  auto omega = GENERATE(take(1, random(1., 1e3)));

  cx_mat::fixed<3, 3> TensorPlateVacuum(fill::zeros);
  GTPlateVacuum.calculate_tensor(omega, {-k_x, -k_y}, TensorPlateVacuum);

  cx_mat::fixed<3, 3> TensorVacuum(fill::zeros);
  GTVacuum.calculate_tensor(omega, {-k_x, -k_y}, TensorVacuum);

  cx_mat::fixed<3, 3> TensorPlate(fill::zeros);
  GTPlate.calculate_tensor(omega, {-k_x, -k_y}, TensorPlate);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      REQUIRE(TensorPlateVacuum(i, j) ==
              TensorVacuum(i, j) + TensorPlate(i, j));
    }
  }
}

TEST_CASE("GreensTensorPlateVacuum integrate_k function returns the sum "
          "of plate and vacuum Green's tensor",
          "[GreensTensorPlateVacuum]") {

  double omega_p = 9;
  double gamma = 0.1;
  auto perm = std::make_shared<PermittivityDrude>(omega_p, gamma);
  auto refl = std::make_shared<ReflectionCoefficientsLocBulk>(perm);

  double v = 1E-5;
  double za = 0.1;
  double beta = 1E4;
  double delta_cut = 20;
  vec::fixed<2> rel_err = {1E-8, 1E-6};

  GreensTensorPlateVacuum GTPlateVacuum(v, beta, za, refl, delta_cut, rel_err);
  GreensTensorPlate GTPlate(v, beta, za, refl, delta_cut, rel_err);
  GreensTensorVacuum GTVacuum(v, beta, rel_err(0));

  auto omega = GENERATE(take(1, random(1., 1e3)));

  cx_mat::fixed<3, 3> TensorPlateVacuum(fill::zeros);
  GTPlateVacuum.integrate_k(omega, TensorPlateVacuum, IM, UNIT);

  cx_mat::fixed<3, 3> TensorVacuum(fill::zeros);
  GTVacuum.integrate_k(omega, TensorVacuum, IM, UNIT);

  cx_mat::fixed<3, 3> TensorPlate(fill::zeros);
  GTPlate.integrate_k(omega, TensorPlate, IM, UNIT);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      REQUIRE(TensorPlateVacuum(i, j) ==
              TensorVacuum(i, j) + TensorPlate(i, j));
    }
  }
}
