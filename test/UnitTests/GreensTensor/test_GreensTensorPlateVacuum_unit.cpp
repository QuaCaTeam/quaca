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

TEST_CASE("GreensTensorPlateVacuum calculate_tensor function returns the sum "
          "of plate and vacuum Green's tensor",
          "[GreensTensorPlateVacuum]") {

  double omega_p = 9;
  double gamma = 0.1;
  PermittivityDrude perm(omega_p, gamma);
  ReflectionCoefficientsLocBulk refl(&perm);

  double v = 1E-8;
  double za = 0.1;
  double beta = 0.01;
  double delta_cut = 20;
  vec::fixed<2> rel_err = {1E-8, 1E-6};

  GreensTensorPlateVacuum GTPlateVacuum(v, za, beta, &refl, delta_cut, rel_err);
  GreensTensorPlate GTPlate(v, za, beta, &refl, delta_cut, rel_err);
  GreensTensorVacuum GTVacuum(v, beta, rel_err(0));

  auto k_x = GENERATE(-12.42,0.124);
  auto k_y = GENERATE(-6.543,34.123);
  auto omega = GENERATE(-1.1, 5.3);
  //Only for \omega^2 - k^2 >= 0 GreensTensorVacuum returns non-trivial results
  omega *= sqrt(pow(k_x,2) + pow(k_y,2));

  Options_GreensTensor opts;
  opts.class_pt = &GTPlateVacuum;
  opts.fancy_complex = IM;
  opts.omega = omega;
  opts.kvec(0) = k_x;
  opts.kvec(1) = k_y;
  cx_mat::fixed<3, 3> TensorPlateVacuum(fill::zeros);
  GTPlateVacuum.calculate_tensor(TensorPlateVacuum, opts);

  opts.class_pt = &GTVacuum;
  cx_mat::fixed<3, 3> TensorVacuum(fill::zeros);
  GTVacuum.calculate_tensor(TensorVacuum, opts);

  opts.class_pt = &GTPlate;
  cx_mat::fixed<3, 3> TensorPlate(fill::zeros);
  GTPlate.calculate_tensor(TensorPlate, opts);

  //Ensure non-trivial results
  REQUIRE(!TensorPlateVacuum.is_zero());
  REQUIRE(!TensorVacuum.is_zero());
  REQUIRE(!TensorPlate.is_zero());

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      REQUIRE(TensorPlateVacuum(i, j) ==
              TensorVacuum(i, j) + TensorPlate(i, j));
    };
  };
};

TEST_CASE("GreensTensorPlateVacuum integrate_k function returns the sum "
          "of plate and vacuum Green's tensor",
          "[GreensTensorPlateVacuum]") {

  double omega_p = 9;
  double gamma = 0.1;
  PermittivityDrude perm(omega_p, gamma);
  ReflectionCoefficientsLocBulk refl(&perm);

  double v = 1E-5;
  double za = 0.1;
  double beta = 1E4;
  double delta_cut = 20;
  vec::fixed<2> rel_err = {1E-8, 1E-6};

  GreensTensorPlateVacuum GTPlateVacuum(v, za, beta, &refl, delta_cut, rel_err);
  GreensTensorPlate GTPlate(v, za, beta, &refl, delta_cut, rel_err);
  GreensTensorVacuum GTVacuum(v, beta, rel_err(0));

  auto k_x = GENERATE(-10.32,2.43);
  auto k_y = GENERATE(-6.53,32.43);
  auto omega = GENERATE(3.2);
  //Only for \omega^2 - k^2 >= 0 GreensTensorVacuum returns non-trivial results
  omega *= sqrt(pow(k_x,2) + pow(k_y,2));

  Options_GreensTensor opts_pv;
  opts_pv.class_pt = &GTPlateVacuum;
  opts_pv.fancy_complex = IM;
  opts_pv.omega = omega;
  cx_mat::fixed<3, 3> TensorPlateVacuum(fill::zeros);
  GTPlateVacuum.integrate_k(TensorPlateVacuum, opts_pv);

  Options_GreensTensor opts_vac = opts_pv;
  opts_vac.class_pt = &GTVacuum;
  cx_mat::fixed<3, 3> TensorVacuum(fill::zeros);
  GTVacuum.integrate_k(TensorVacuum, opts_vac);

  Options_GreensTensor opts_plate = opts_pv;
  opts_plate.class_pt = &GTPlate;
  cx_mat::fixed<3, 3> TensorPlate(fill::zeros);
  GTPlate.integrate_k(TensorPlate, opts_plate);

  //Ensure non-trivial results
  REQUIRE(!TensorPlateVacuum.is_zero());
  REQUIRE(!TensorVacuum.is_zero());
  REQUIRE(!TensorPlate.is_zero());

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      REQUIRE(TensorPlateVacuum(i, j) ==
              TensorVacuum(i, j) + TensorPlate(i, j));
    };
  };
};
