#include "Quaca.h"
#include "catch.hpp"
#include <iostream>

TEST_CASE("PowerSpectrum constructors work as expected", "[PowerSpectrum]") {
  SECTION("Direct constructor") {
    auto v = GENERATE(take(1, random(0., 1.)));
    double beta = GENERATE(take(1, random(1e-5, 1e3)));
    double omega_a = GENERATE(take(1, random(0., 1e3)));
    double alpha_zero = GENERATE(take(1, random(1e-5, 1.)));

    double relerr_k = 1E-9;
    auto greens = std::make_shared<GreensTensorVacuum>(v, beta, relerr_k);
    auto alpha = std::make_shared<Polarizability>(omega_a, alpha_zero, greens);
    PowerSpectrum powerspectrum(greens, alpha);

    REQUIRE(Approx(powerspectrum.get_greens_tensor()->get_v()).epsilon(1e-6) ==
            v);
    REQUIRE(
        Approx(powerspectrum.get_greens_tensor()->get_beta()).epsilon(1e-6) ==
        beta);
    REQUIRE(Approx(powerspectrum.get_polarizability()->get_omega_a())
                .epsilon(1e-6) == omega_a);
    REQUIRE(Approx(powerspectrum.get_polarizability()->get_alpha_zero())
                .epsilon(1e-6) == alpha_zero);
  }

  SECTION("json file constructor") {
    double omega_a = 1.3;
    double alpha_zero = 6e-9;

    double beta = 5.;
    double v = 0.01;

    PowerSpectrum powerspectrum("../data/test_files/PowerSpectrum.json");

    REQUIRE(Approx(powerspectrum.get_greens_tensor()->get_v()).epsilon(1e-6) ==
            v);
    REQUIRE(
        Approx(powerspectrum.get_greens_tensor()->get_beta()).epsilon(1e-6) ==
        beta);
    REQUIRE(Approx(powerspectrum.get_polarizability()->get_omega_a())
                .epsilon(1e-6) == omega_a);
    REQUIRE(Approx(powerspectrum.get_polarizability()->get_alpha_zero())
                .epsilon(1e-6) == alpha_zero);
  }
}

TEST_CASE("Power spectrum without bath is hermitian", "[PowerSpectrum]") {
  // Generate power spectrum
  auto v = GENERATE(1e-4);
  auto beta = GENERATE(10.);
  auto omega_a = GENERATE(1.6);
  auto alpha_zero = GENERATE(1e-9);

  double relerr_k = 1E-9;
  auto greens = std::make_shared<GreensTensorVacuum>(v, beta, relerr_k);
  auto alpha = std::make_shared<Polarizability>(omega_a, alpha_zero, greens);
  PowerSpectrum powerspectrum(greens, alpha);

  // Matrices to store results
  cx_mat::fixed<3, 3> lhs(fill::zeros);
  cx_mat::fixed<3, 3> rhs(fill::zeros);

  auto omega = GENERATE(-21.65, -5.34, 1.32, 65.32);

  powerspectrum.calculate(omega, lhs, FULL);
  powerspectrum.calculate(omega, rhs, FULL);

  // Ensure non-trivial results
  REQUIRE(!lhs.is_zero());
  REQUIRE(!rhs.is_zero());

  REQUIRE(approx_equal(lhs, trans(rhs), "reldiff", 1e-8));
}

TEST_CASE(
    "Power spectrum without bath reduces to polarizability in the static case",
    "[PowerSpectrum]") {
  // Generate power spectrum
  auto beta = GENERATE(10.);
  auto omega_a = GENERATE(1.6);
  auto alpha_zero = GENERATE(1e-9);
  // Ensure equilibrium result
  double v = 1e-11;

  double relerr_k = 1E-9;
  auto greens = std::make_shared<GreensTensorVacuum>(v, beta, relerr_k);
  auto alpha = std::make_shared<Polarizability>(omega_a, alpha_zero, greens);
  PowerSpectrum powerspectrum(greens, alpha);

  // Matrices to store results
  cx_mat::fixed<3, 3> lhs(fill::zeros);
  cx_mat::fixed<3, 3> rhs(fill::zeros);

  auto omega = GENERATE(0.321, 1.65, 56.21);
  // Compute the powerspectrum
  powerspectrum.calculate(omega, lhs, FULL);

  // Set integration options for the polarizability
  alpha->calculate_tensor(omega, rhs, IM);
  rhs *= 1. / (M_PI * (1. - exp(-beta * omega)));

  // Ensure non-trivial result
  REQUIRE(!lhs.is_zero());
  REQUIRE(!rhs.is_zero());

  REQUIRE(approx_equal(lhs, rhs, "reldiff", 10e-8));
}

TEST_CASE("Power spectrum with bath is hermitian", "[PowerSpectrum]") {
  // Generate power spectrum
  auto v = GENERATE(1e-4);
  auto beta = GENERATE(10.);
  auto omega_a = GENERATE(1.6);
  auto alpha_zero = GENERATE(1e-9);
  double gamma = GENERATE(0.21);

  double relerr_k = 1E-9;
  auto greens = std::make_shared<GreensTensorVacuum>(v, beta, relerr_k);
  auto ohmic_kernel = std::make_shared<OhmicMemoryKernel>(gamma);
  auto alpha = std::make_shared<Polarizability>(omega_a, alpha_zero,
                                                ohmic_kernel, greens);
  PowerSpectrum powerspectrum(greens, alpha);

  // Matrices to store results
  cx_mat::fixed<3, 3> lhs(fill::zeros);
  cx_mat::fixed<3, 3> rhs(fill::zeros);

  auto omega = GENERATE(-.97, -1e-3, 0.43, .9);

  powerspectrum.calculate(omega, lhs, FULL);
  powerspectrum.calculate(omega, rhs, FULL);

  // Ensure non-trivial results
  REQUIRE(!lhs.is_zero());
  REQUIRE(!rhs.is_zero());
  REQUIRE(approx_equal(lhs, trans(rhs), "reldiff", 1e-8));
}

TEST_CASE(
    "Power spectrum with bath reduces to polarizability in the static case",
    "[PowerSpectrum]") {
  // Generate power spectrum
  auto beta = GENERATE(10.);
  auto omega_a = GENERATE(1.6);
  auto alpha_zero = GENERATE(1e-9);
  double gamma = GENERATE(0.21);
  // Ensure equilibrium result
  double v = 1e-14;

  double relerr_k = 1E-9;
  auto greens = std::make_shared<GreensTensorVacuum>(v, beta, relerr_k);
  auto ohmic_kernel = std::make_shared<OhmicMemoryKernel>(gamma);
  auto alpha = std::make_shared<Polarizability>(omega_a, alpha_zero,
                                                ohmic_kernel, greens);
  PowerSpectrum powerspectrum(greens, alpha);

  // Matrices to store results
  cx_mat::fixed<3, 3> lhs(fill::zeros);
  cx_mat::fixed<3, 3> rhs(fill::zeros);

  auto omega = GENERATE(take(1, random(0., 1e2)));

  powerspectrum.calculate(omega, lhs, FULL);

  // Set integration options for the polarizability
  alpha->calculate_tensor(omega, rhs, IM);
  rhs *= 1. / (M_PI * (1. - exp(-beta * omega)));

  // Ensure non-trivial result
  REQUIRE(!lhs.is_zero());
  REQUIRE(!rhs.is_zero());

  REQUIRE(approx_equal(lhs, rhs, "reldiff", 10e-8));
}
