#include <armadillo>
#include <complex>

#include "Quaca.h"
#include "catch.hpp"
#include <iomanip> // std::setprecision

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
  auto omega = GENERATE(take(5, random(-1e2, 1e2)));
  auto k_x = GENERATE(take(5, random(0., 1e2)));
  auto k_y = GENERATE(take(5, random(0., 1e2)));
  double omega_p = 9;
  double gamma = 0.1;
  double v = 1e-2;
  double za = 0.1;
  double kappa_double;
  double cos_phi, k;
  std::complex<double> kappa, volume_element;
  PermittivityDrude perm(gamma, omega_p);
  GreensTensorPlate Greens(v, za, NAN, &perm, NAN);
  struct Options_GreensTensor opts;
  opts.class_pt = &Greens;

  // First, the calculate_tensor operation is used to generate the
  // Green's tensor with fancy_I
  cx_mat::fixed<3, 3> Green(fill::zeros);
  cx_mat::fixed<3, 3> Green_fancy_I_ct(fill::zeros);

  opts.kvec(0) = k_x;
  opts.kvec(1) = k_y;
  opts.omega = omega + k_x * v;
  k = sqrt(k_x * k_x + k_y * k_y);
  cos_phi = k_x / k;
  kappa = sqrt(std::complex<double>(k * k - opts.omega * opts.omega, 0.));
  kappa = std::complex<double>(std::abs(kappa.real()), -std::abs(kappa.imag()));
  volume_element = kappa * k / (k - cos_phi * v * opts.omega);

  Greens.calculate_tensor(Green, opts);
  if (opts.omega < 0) {
    volume_element = conj(volume_element);
  }
  Green *= volume_element;
  Green_fancy_I_ct = (Green - trans(Green)) / (2. * I);

  // Second, the integrand_k_2d operation is used
  cx_mat::fixed<3, 3> Green_fancy_I_ik2d(fill::zeros);
  if (kappa.real() == 0.) {
    kappa_double = kappa.imag();
  } else {
    kappa_double = kappa.real();
  }
  opts.omega = omega;
  opts.fancy_I = true;
  opts.kvec(0) = acos(k_x / k);
  // loop over all indices
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      opts.indices(0) = i;
      opts.indices(1) = j;
      Green_fancy_I_ik2d(i, j) =
          (2 * M_PI) * Greens.integrand_k_2d(kappa_double, &opts);
      if (i != j) {
        // As the prefactor I can not be evaluated in a purely real integration
        // routine, it was dropped in integrate_k_2d and has to be inserted here
        Green_fancy_I_ik2d(i, j) *= I;
      }
    }
  }
  // std::cout << "omega=" << omega << '\n';
  // std::cout << "kx*v=" << k_x * v << '\n';
  // std::cout << "k=" << k << '\n';
  // Green_fancy_I_ct.raw_print();
  // std::cout << '\n';
  // Green_fancy_I_ik2d.raw_print();
  // std::cout << '\n';

  REQUIRE(approx_equal(Green_fancy_I_ct, Green_fancy_I_ik2d, "reldiff", 10E-4));
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

TEST_CASE("Integrated Green's tensor works properly", "[GreensTensorPlate]") {

  SECTION("Integral over Green_fancy_I obeys the crossing relation") {
    auto omega = GENERATE(take(1, random(0., 1e2)));
    GreensTensorPlate Greens("../data/test_files/GreensTensorPlate.ini");
    struct Options_GreensTensor opts;
    opts.class_pt = &Greens;
    cx_mat::fixed<3, 3> Greens_lhs(fill::zeros);
    cx_mat::fixed<3, 3> Greens_rhs(fill::zeros);
    // Test of fancy_I
    opts.omega = omega;
    opts.fancy_I = true;
    Greens.integrate_k_1d(Greens_lhs, opts);

    opts.omega = -omega;
    Greens.integrate_k_1d(Greens_rhs, opts);

    REQUIRE(approx_equal(Greens_lhs, -strans(Greens_rhs), "reldiff", 10E-4));
  };
  SECTION("Integral over Green_fancy_R obeys the crossing relation") {
    auto omega = GENERATE(take(1, random(0., 1e2)));
    GreensTensorPlate Greens("../data/test_files/GreensTensorPlate.ini");
    struct Options_GreensTensor opts;
    opts.class_pt = &Greens;
    cx_mat::fixed<3, 3> Greens_lhs(fill::zeros);
    cx_mat::fixed<3, 3> Greens_rhs(fill::zeros);
    // Test of fancy_R
    opts.omega = omega;
    opts.fancy_R = true;
    Greens.integrate_k_1d(Greens_lhs, opts);

    opts.omega = -omega;
    Greens.integrate_k_1d(Greens_rhs, opts);

    REQUIRE(approx_equal(Greens_lhs, strans(Greens_rhs), "reldiff", 10E-4));
    opts.fancy_R = false;
  };
  SECTION("Integral over Green_fancy_I_kv obeys the crossing relation") {
    auto omega = GENERATE(take(1, random(0., 1e2)));
    GreensTensorPlate Greens("../data/test_files/GreensTensorPlate.ini");
    struct Options_GreensTensor opts;
    opts.class_pt = &Greens;
    cx_mat::fixed<3, 3> Greens_lhs(fill::zeros);
    cx_mat::fixed<3, 3> Greens_rhs(fill::zeros);
    // Test of fancy_I_kv
    opts.omega = omega;
    opts.fancy_I_kv = true;
    Greens.integrate_k_1d(Greens_lhs, opts);

    opts.omega = -omega;
    Greens.integrate_k_1d(Greens_rhs, opts);
    REQUIRE(approx_equal(Greens_lhs, strans(Greens_rhs), "reldiff", 10E-4));
  };
  SECTION("Integral over Green_fancy_I_temp obeys the crossing relation") {
    auto omega = GENERATE(take(1, random(0., 1e-1)));
    GreensTensorPlate Greens("../data/test_files/GreensTensorPlate.ini");
    struct Options_GreensTensor opts;
    opts.class_pt = &Greens;
    cx_mat::fixed<3, 3> Greens_lhs(fill::zeros);
    cx_mat::fixed<3, 3> Greens_rhs(fill::zeros);
    cx_mat::fixed<3, 3> Greens_rhs1(fill::zeros);
    cx_mat::fixed<3, 3> Greens_rhs2(fill::zeros);
    // Test of fancy_I_kv
    opts.omega = -omega;
    opts.fancy_I_temp = true;
    Greens.integrate_k_1d(Greens_lhs, opts);

    opts.omega = omega;
    opts.fancy_I_temp = false;
    opts.fancy_I = true;
    Greens.integrate_k_1d(Greens_rhs1, opts);
    opts.fancy_I = false;
    opts.fancy_I_temp = true;
    Greens.integrate_k_1d(Greens_rhs2, opts);
    Greens_rhs = -Greens_rhs1 + Greens_rhs2;
    REQUIRE(approx_equal(Greens_lhs, strans(Greens_rhs), "reldiff", 10E-4));
  };
  SECTION("Integral over Green_fancy_I_kv_temp obeys the crossing relation") {
    // auto omega = GENERATE(take(1, random(1e2, 1e2)));
    auto omega = GENERATE(take(1, random(0., 1e-1)));
    GreensTensorPlate Greens("../data/test_files/GreensTensorPlate.ini");
    struct Options_GreensTensor opts;
    opts.class_pt = &Greens;
    cx_mat::fixed<3, 3> Greens_lhs(fill::zeros);
    cx_mat::fixed<3, 3> Greens_rhs(fill::zeros);
    cx_mat::fixed<3, 3> Greens_rhs1(fill::zeros);
    cx_mat::fixed<3, 3> Greens_rhs2(fill::zeros);
    // Test of fancy_I_kv
    opts.omega = -omega;
    opts.fancy_I_kv_temp = true;
    Greens.integrate_k_1d(Greens_lhs, opts);

    opts.omega = omega;
    opts.fancy_I_kv_temp = false;
    opts.fancy_I_kv = true;
    Greens.integrate_k_1d(Greens_rhs1, opts);
    opts.fancy_I_kv = false;
    opts.fancy_I_kv_temp = true;

    Greens.integrate_k_1d(Greens_rhs2, opts);
    Greens_rhs = -Greens_rhs1 + Greens_rhs2;
    REQUIRE(approx_equal(Greens_lhs, -strans(Greens_rhs), "reldiff", 10E-4));
  };
};

TEST_CASE("Integrated Green's tensor matches asymptotes",
          "[GreensTensorPlate]") {
  SECTION("Low-frequency asymptote") {
    // \omega << \omega_p and \omega << v/z_a
    std::complex<double> I(0.0, 1.0);
    double omega_p = 9;
    double gamma = 0.1;
    double v = 1e-5;
    double za = 0.1;
    auto omega = GENERATE(take(1, random(-0.1 * 1e-6, 0.1 * 1e-6)));
    double delta_cut = 30;
    PermittivityDrude perm(gamma, omega_p);
    GreensTensorPlate Greens(v, za, NAN, &perm, delta_cut);
    struct Options_GreensTensor opts;
    opts.class_pt = &Greens;
    cx_mat::fixed<3, 3> GT_Ana(fill::zeros);
    cx_mat::fixed<3, 3> GT_Num(fill::zeros);
    GT_Ana(0, 0) = 2 * omega * gamma / pow(omega_p, 2) / pow(2 * za, 3);
    GT_Ana(1, 1) = GT_Ana(0, 0);
    GT_Ana(2, 2) = 2. * GT_Ana(0, 0);
    GT_Ana(0, 2) = (2 * 3 * v * gamma / (pow(omega_p, 2) * pow(2 * za, 4))) / I;
    GT_Ana(2, 0) = -GT_Ana(0, 2);
    opts.fancy_I = true;
    opts.omega = omega;
    Greens.integrate_k_1d(GT_Num, opts);
    REQUIRE(approx_equal(GT_Ana, GT_Num, "reldiff", 10E-4));
  };
};
