#include "Quaca.h"
#include "catch.hpp"
#include <iomanip>
#include <iostream>

TEST_CASE("PowerSpectrumHarmOsc constructors work as expected",
          "[PowerSpectrumHarmOsc]") {
  SECTION("Direct constructor") {
    auto v = GENERATE(take(1, random(0., 1.)));
    double beta = GENERATE(take(1, random(1e-5, 1e3)));
    double omega_a = GENERATE(take(1, random(0., 1e3)));
    double alpha_zero = GENERATE(take(1, random(1e-5, 1.)));

    double relerr_k = 1E-9;
    GreensTensorVacuum greens(v, beta, relerr_k);
    PolarizabilityNoBath alpha(omega_a, alpha_zero, &greens);
    PowerSpectrumHarmOsc powerspectrum(&greens, &alpha);

    REQUIRE(Approx(powerspectrum.get_greens_tensor()->get_v()).epsilon(1e-6) ==
            v);
    REQUIRE(
        Approx(powerspectrum.get_greens_tensor()->get_beta()).epsilon(1e-6) ==
        beta);
    REQUIRE(Approx(powerspectrum.get_polarizability()->get_omega_a())
                .epsilon(1e-6) == omega_a);
    REQUIRE(Approx(powerspectrum.get_polarizability()->get_alpha_zero())
                .epsilon(1e-6) == alpha_zero);
    REQUIRE(powerspectrum.has_bath == false);
  };

  SECTION("json file constructor") {
    double omega_a = 1.3;
    double alpha_zero = 6e-9;

    double beta = 5.;
    double v = 0.01;

    PowerSpectrumHarmOsc powerspectrum(
        "../data/test_files/PowerSpectrumHarmOsc.json");

    REQUIRE(Approx(powerspectrum.get_greens_tensor()->get_v()).epsilon(1e-6) ==
            v);
    REQUIRE(
        Approx(powerspectrum.get_greens_tensor()->get_beta()).epsilon(1e-6) ==
        beta);
    REQUIRE(Approx(powerspectrum.get_polarizability()->get_omega_a())
                .epsilon(1e-6) == omega_a);
    REQUIRE(Approx(powerspectrum.get_polarizability()->get_alpha_zero())
                .epsilon(1e-6) == alpha_zero);
    REQUIRE(powerspectrum.has_bath == true);
  };
}

TEST_CASE("Power spectrum without bath is hermitian",
          "[PowerSpectrumHarmOsc]") {
  // Generate power spectrum
  auto v = GENERATE(1e-4);
  auto beta = GENERATE(10.);
  auto omega_a = GENERATE(1.6);
  auto alpha_zero = GENERATE(1e-9);

  double relerr_k = 1E-9;
  GreensTensorVacuum greens(v, beta, relerr_k);
  PolarizabilityNoBath alpha(omega_a, alpha_zero, &greens);
  PowerSpectrumHarmOsc powerspectrum(&greens, &alpha);

  // Matrices to store results
  cx_mat::fixed<3, 3> lhs(fill::zeros);
  cx_mat::fixed<3, 3> rhs(fill::zeros);

  auto omega = GENERATE(-21.65,-5.34,1.32,65.32);

  Options_PowerSpectrum opts;
  opts.spectrum = FULL;
  opts.omega = omega;

  powerspectrum.calculate(lhs, opts);
  powerspectrum.calculate(rhs, opts);

  //Ensure non-trivial results
  REQUIRE(!lhs.is_zero());
  REQUIRE(!rhs.is_zero());

  REQUIRE(approx_equal(lhs, trans(rhs), "reldiff", 1e-8));
}

TEST_CASE(
    "Power spectrum without bath reduces to polarizability in the static case",
    "[PowerSpectrumHarmOsc]") {
  // Generate power spectrum
  auto beta = GENERATE(10.);
  auto omega_a = GENERATE(1.6);
  auto alpha_zero = GENERATE(1e-9);
  // Ensure equilibrium result
  double v = 1e-11;

  double relerr_k = 1E-9;
  GreensTensorVacuum greens(v, beta, relerr_k);
  PolarizabilityNoBath alpha(omega_a, alpha_zero, &greens);
  PowerSpectrumHarmOsc powerspectrum(&greens, &alpha);

  // Matrices to store results
  cx_mat::fixed<3, 3> lhs(fill::zeros);
  cx_mat::fixed<3, 3> rhs(fill::zeros);

  auto omega = GENERATE(0.321,1.65,56.21);
  // Compute the powerspectrum
  Options_PowerSpectrum opts_S;
  opts_S.spectrum = FULL;
  opts_S.omega = omega;
  powerspectrum.calculate(lhs, opts_S);

  // Set integration options for the polarizability
  Options_Polarizability opts_alpha;
  opts_alpha.omega = omega;
  opts_alpha.fancy_complex = IM;
  alpha.calculate_tensor(rhs, opts_alpha);
  rhs *= 1. / (M_PI * (1. - exp(-beta * omega)));

  //Ensure non-trivial result
  REQUIRE(!lhs.is_zero());
  REQUIRE(!rhs.is_zero());

  REQUIRE(approx_equal(lhs, rhs, "reldiff", 10e-8));
}

TEST_CASE("Power spectrum with bath is hermitian", "[PowerSpectrumHarmOsc]") {
  // Generate power spectrum
  auto v = GENERATE(1e-4);
  auto beta = GENERATE(10.);
  auto omega_a = GENERATE(1.6);
  auto alpha_zero = GENERATE(1e-9);
  double gamma = GENERATE(0.21);

  double relerr_k = 1E-9;
  GreensTensorVacuum greens(v, beta, relerr_k);
  OhmicMemoryKernel ohmic_kernel(gamma);
  PolarizabilityBath alpha(omega_a, alpha_zero, &ohmic_kernel, &greens);
  PowerSpectrumHarmOsc powerspectrum(&greens, &alpha);

  // Matrices to store results
  cx_mat::fixed<3, 3> lhs(fill::zeros);
  cx_mat::fixed<3, 3> rhs(fill::zeros);

  auto omega = GENERATE(-.97,-1e-3,0.43,.9);

  Options_PowerSpectrum opts;
  opts.spectrum = FULL;
  opts.omega = omega;

  powerspectrum.calculate(lhs, opts);
  powerspectrum.calculate(rhs, opts);

  //Ensure non-trivial results
  REQUIRE(!lhs.is_zero());
  REQUIRE(!rhs.is_zero());
  REQUIRE(approx_equal(lhs, trans(rhs), "reldiff", 1e-8));
}

TEST_CASE(
    "Power spectrum with bath reduces to polarizability in the static case",
    "[PowerSpectrumHarmOsc]") {
  // Generate power spectrum
  auto beta = GENERATE(10.);
  auto omega_a = GENERATE(1.6);
  auto alpha_zero = GENERATE(1e-9);
  double gamma = GENERATE(0.21);
  // Ensure equilibrium result
  double v = 1e-14;

  double relerr_k = 1E-9;
  GreensTensorVacuum greens(v, beta, relerr_k);
  OhmicMemoryKernel ohmic_kernel(gamma);
  PolarizabilityBath alpha(omega_a, alpha_zero, &ohmic_kernel, &greens);
  PowerSpectrumHarmOsc powerspectrum(&greens, &alpha);

  // Matrices to store results
  cx_mat::fixed<3, 3> lhs(fill::zeros);
  cx_mat::fixed<3, 3> rhs(fill::zeros);

  auto omega = GENERATE(take(1, random(0., 1e2)));
  // Compute the powerspectrum
  Options_PowerSpectrum opts_S;
  opts_S.spectrum = FULL;
  opts_S.omega = omega;
  powerspectrum.calculate(lhs, opts_S);

  // Set integration options for the polarizability
  Options_Polarizability opts_alpha;
  opts_alpha.omega = omega;
  opts_alpha.fancy_complex = IM;
  alpha.calculate_tensor(rhs, opts_alpha);
  rhs *= 1. / (M_PI * (1. - exp(-beta * omega)));

  //Ensure non-trivial result
  REQUIRE(!lhs.is_zero());
  REQUIRE(!rhs.is_zero());

  REQUIRE(approx_equal(lhs, rhs, "reldiff", 10e-8));
}
