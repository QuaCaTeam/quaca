#include <armadillo>
#include <complex>

#include "Quaca.h"
#include "catch.hpp"

TEST_CASE("Construction of Green's tensor works properly",
          "[GreensTensorPlate]") {
  SECTION("Construction with direct input") {
    double omega = 1;
    double omega_p = 9;
    double gamma = 0.1;
    double v = 1E-5;
    double za = 0.1;
    double beta = 1E4;
    double delta_cut = 20;

    PermittivityDrude perm(gamma, omega_p);
    GreensTensorPlate Greens(v, za, beta, &perm, delta_cut);
    REQUIRE(Greens.get_za() == za);
    REQUIRE(Greens.get_delta_cut() == delta_cut);
  };

  SECTION("Construction from .ini file") {
    GreensTensorPlate Greens("../data/test_files/GreensTensorPlate.ini");
    REQUIRE(Greens.get_za() == 0.1);
    REQUIRE(Greens.get_delta_cut() == 20);
  };
};
TEST_CASE("The operations calculate_tensor and the integrand_k_2d coincide",
          "[GreensTensorPlate]") {
  // Here we considered also the volume element from the integration.
  std::complex<double> I(0.0, 1.0);
  auto omega = GENERATE(take(5, random(0., 1e2)));
  auto k_x = GENERATE(take(5, random(0., 1e2)));
  auto k_y = GENERATE(take(5, random(0., 1e2)));
  double omega_p = 9;
  double gamma = 0.1;
  double v = 1E-5;
  double za = 0.1;
  PermittivityDrude perm(gamma, omega_p);
  GreensTensorPlate Greens(v, za, NAN, &perm, NAN);
  struct Options_GreensTensor opts;
  opts.class_pt = &Greens;

  // First, the calculate_tensor operation is used to generate the
  // Green's tensor with fancy_I
  cx_mat::fixed<3, 3> Green(fill::zeros);
  cx_mat::fixed<3, 3> Green_dagger(fill::zeros);
  cx_mat::fixed<3, 3> Green_fancy_I_ct(fill::zeros);

  opts.omega = omega + k_x * v;
  opts.kvec(0) = k_x;
  opts.kvec(1) = k_y;
  Greens.calculate_tensor(Green, opts);

  opts.omega = -omega;
  Greens.calculate_tensor(Green_dagger, opts);

  Green_fancy_I_ct = (Green - Green_dagger) / std::complex<double>(0e0, 2e0);

  // Second, the integrand_k_2d operation is used
  cx_mat::fixed<3, 3> Green_fancy_I_ik2d(fill::zeros);
  opts.fancy_I = true;
  opts.kvec(0) = acos(k_x / sqrt(k_x * k_x + k_y * k_y));
  // loop over all indices
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      opts.indices(0) = i;
      opts.indices(1) = j;
      Green_fancy_I_ik2d(i, j) = Greens.integrate_k_2d(kappa_double, &opts);
      if (i != j) {
        // As the prefactor I can not be evaluated in a purely real integration
        // routine, it was dropped in integrate_k_2d and has to be inserted here
        Green_fancy_I_ik2d(i, j) *= I;
      }
    }
  }
};

TEST_CASE("Scattered Green's tensor works properly", "[GreensTensorPlate]") {

  SECTION("Green's tensor obeys reciprocity") {
    auto omega = GENERATE(take(5, random(0., 1e2)));
    auto k_x = GENERATE(take(5, random(0., 1e2)));
    auto k_y = GENERATE(take(5, random(0., 1e2)));
    GreensTensorPlate Greens("../data/test_files/GreensTensorPlate.ini");
    struct Options_GreensTensor opts;
    opts.class_pt = &Greens;
    cx_mat::fixed<3, 3> Greens_lhs(fill::zeros);
    cx_mat::fixed<3, 3> Greens_rhs(fill::zeros);

    opts.omega = omega;
    opts.kvec(0) = k_x;
    opts.kvec(1) = k_y;
    Greens.calculate_tensor(Greens_lhs, opts);

    opts.omega = omega;
    opts.kvec(0) = -k_x;
    opts.kvec(1) = -k_y;
    Greens.calculate_tensor(Greens_rhs, opts);

    REQUIRE(approx_equal(Greens_lhs, strans(Greens_rhs), "reldiff", 10E-15));
  };

  SECTION("Green's tensor obeys reality condition") {
    auto omega = GENERATE(take(5, random(0., 1e2)));
    auto k_x = GENERATE(take(5, random(0., 1e2)));
    auto k_y = GENERATE(take(5, random(0., 1e2)));
    GreensTensorPlate Greens("../data/test_files/GreensTensorPlate.ini");
    struct Options_GreensTensor opts;
    opts.class_pt = &Greens;
    cx_mat::fixed<3, 3> Greens_lhs(fill::zeros);
    cx_mat::fixed<3, 3> Greens_rhs(fill::zeros);

    opts.omega = omega;
    opts.kvec(0) = k_x;
    opts.kvec(1) = k_y;
    Greens.calculate_tensor(Greens_lhs, opts);

    opts.omega = -omega;
    opts.kvec(0) = k_x;
    opts.kvec(1) = k_y;
    Greens.calculate_tensor(Greens_rhs, opts);

    REQUIRE(approx_equal(Greens_lhs, trans(Greens_rhs), "reldiff", 10E-15));
  };
};

//  struct Options_GreensTensor opts;
//  opts.fancy_R=true;
//  opts.omega = 3.0;
//  opts.kvec = { 10.0 , NAN };
//  opts.indices = { 0 , 0};
//  opts.class_pt = &Greens;
//  cx_mat::fixed<3,3> test(fill::zeros);
//  Greens.integrate_k_1d(test, opts);
//  std::cout << test << std::endl;
////  std::cout << Greens.integrand_k_2d( 1., &opts) << std::endl;
//  std::cout << Greens.integrand_k_1d( 1., &opts) << std::endl;
