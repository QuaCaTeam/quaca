#include "Quaca.h"
#include "catch.hpp"
#include <iostream>
#include <iomanip>

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
    REQUIRE(powerspectrum.has_bath==false);
  };

  SECTION("ini file constructor") {
    double omega_a = 1.3;
    double alpha_zero = 6e-9;

    double beta = 5.;
    double v = 0.01;

    PowerSpectrumHarmOsc powerspectrum(
        "../data/test_files/PowerSpectrumHarmOsc.ini");

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

TEST_CASE("Power spectrum without bath is hermitian", "[PowerSpectrumHarmOsc]") {
  // Generate randomnized power spectrum
  auto v = GENERATE(take(1, random(0., 1.)));
  double beta = GENERATE(take(1, random(1., 1e2)));
  double omega_a = GENERATE(take(1, random(0.1, 1e1)));
  double alpha_zero = GENERATE(take(1, random(1e-9, 1e-7)));

  double relerr_k = 1E-9;
  GreensTensorVacuum greens(v, beta, relerr_k);
  PolarizabilityNoBath alpha(omega_a, alpha_zero, &greens);
  PowerSpectrumHarmOsc powerspectrum(&greens, &alpha);

  // Matrices to store results
  cx_mat::fixed<3, 3> lhs(fill::zeros);
  cx_mat::fixed<3, 3> rhs(fill::zeros);

  auto omega = GENERATE(take(5, random(-1e1, 1e1)));

  Options_PowerSpectrum opts;
  opts.spectrum = FULL;
  opts.omega = omega;

  powerspectrum.calculate(lhs, opts);
  powerspectrum.calculate(rhs, opts);

  REQUIRE(approx_equal(lhs, trans(rhs), "reldiff", 1e-8));
}

TEST_CASE("Power spectrum without bath reduces to polarizability in the static case",
          "[PowerSpectrumHarmOsc]") {
  // Generate randomnized power spectrum
  double beta = GENERATE(take(3, random(1e-5, 1e5)));
  double omega_a = GENERATE(take(3, random(0., 1e1)));
  double alpha_zero = GENERATE(take(3, random(1e-9, 1e-7)));
  // Ensure equilibrium result
  double v = 1e-11;

  double relerr_k = 1E-9;
  GreensTensorVacuum greens(v, beta, relerr_k);
  PolarizabilityNoBath alpha(omega_a, alpha_zero, &greens);
  PowerSpectrumHarmOsc powerspectrum(&greens, &alpha);

  // Matrices to store results
  cx_mat::fixed<3, 3> lhs(fill::zeros);
  cx_mat::fixed<3, 3> rhs(fill::zeros);

  auto omega = GENERATE(take(5, random(0., 1e2)));
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

  REQUIRE(approx_equal(lhs, rhs, "reldiff", 10e-8));
}

TEST_CASE("Power spectrum with bath is hermitian", "[PowerSpectrumHarmOsc]") {
  // Generate randomnized power spectrum
  auto v = GENERATE(take(1, random(0., 1.)));
  double beta = GENERATE(take(1, random(1., 1e2)));
  double omega_a = GENERATE(take(1, random(0.1, 1e1)));
  double alpha_zero = GENERATE(take(1, random(1e-9, 1e-7)));
  double gamma = GENERATE(take(1,random(0.,1e1)));

  double relerr_k = 1E-9;
  GreensTensorVacuum greens(v, beta, relerr_k);
  OhmicMemoryKernel ohmic_kernel(gamma);
  PolarizabilityBath alpha(omega_a, alpha_zero, &ohmic_kernel, &greens);
  PowerSpectrumHarmOsc powerspectrum(&greens, &alpha);

  // Matrices to store results
  cx_mat::fixed<3, 3> lhs(fill::zeros);
  cx_mat::fixed<3, 3> rhs(fill::zeros);

  auto omega = GENERATE(take(5, random(-1e1, 1e1)));

  Options_PowerSpectrum opts;
  opts.spectrum = FULL;
  opts.omega = omega;

  powerspectrum.calculate(lhs, opts);
  powerspectrum.calculate(rhs, opts);

  REQUIRE(approx_equal(lhs, trans(rhs), "reldiff", 1e-8));
}

TEST_CASE("Power spectrum with bath reduces to polarizability in the static case",
          "[PowerSpectrumHarmOsc]") {
  // Generate randomnized power spectrum
  double beta = GENERATE(take(1, random(1e-5, 1e5)));
  double omega_a = GENERATE(take(1, random(0., 1e1)));
  double alpha_zero = GENERATE(take(1, random(1e-9, 1e-7)));
  double gamma = GENERATE(take(1,random(0.,1.)));
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
  REQUIRE(approx_equal(lhs, rhs, "reldiff", 10e-8));
}