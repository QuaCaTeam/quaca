#include "Quaca.h"
#include "catch.hpp"
#include <iostream>

TEST_CASE("Magnetic power spectrum is hermitian",
          "[PowerSpectrumHarmOscMagnetic]") {
  // Generate randomnized power spectrum
  double omega_p = 9;
  double gamma = .1;
  double v = 1e-2;
  double za = 0.1;
  double beta = 1000.;
  double delta_cut = 100;
  vec::fixed<2> rel_err = {1E-8, 1E-6};
  PermittivityDrude perm(omega_p, gamma);
  ReflectionCoefficientsLocBulk refl(&perm);
  GreensTensorPlateMagnetic greens(v, za, beta, &refl, delta_cut, rel_err);

  double omega_a = GENERATE(take(1, random(0.1, 1e1)));
  double alpha_zero = GENERATE(take(1, random(1e-9, 1e-7)));
  PolarizabilityNoBathMagnetic alpha(omega_a, alpha_zero, &greens);

  PowerSpectrumHarmOscMagnetic powerspectrum(&greens, &alpha);

  // Matrices to store results
  cx_mat::fixed<3, 3> lhs(fill::zeros);

  auto omega = GENERATE(take(5, random(-1e1, 1e1)));

  Options_PowerSpectrum opts;
  opts.spectrum = FULL;
  opts.omega = omega;

  powerspectrum.calculate(lhs, opts);

  REQUIRE(approx_equal(lhs, trans(lhs), "reldiff", 1e-8));
}

TEST_CASE("Magnetic power spectrum has non-vanishing magnetic contributions",
          "[PowerSpectrumHarmOscMagnetic]") {
  // Generate randomnized power spectrum
  double omega_p = 9;
  double gamma = .1;
  double v = 1e-2;
  double za = 0.1;
  double beta = 1000.;
  double delta_cut = 100;
  vec::fixed<2> rel_err = {1E-8, 1E-6};
  PermittivityDrude perm(omega_p, gamma);
  ReflectionCoefficientsLocBulk refl(&perm);
  GreensTensorPlateMagnetic greens_magnetic(v, za, beta, &refl, delta_cut, rel_err);
  GreensTensorPlate greens(v, za, beta, &refl, delta_cut, rel_err);

  double omega_a = GENERATE(take(1, random(0.1, 1e1)));
  double alpha_zero = GENERATE(take(1, random(1e-9, 1e-7)));
  PolarizabilityNoBathMagnetic alpha_magnetic(omega_a, alpha_zero, &greens_magnetic);
  PolarizabilityNoBath alpha(omega_a, alpha_zero, &greens);

  PowerSpectrumHarmOscMagnetic powerspectrum_magnetic(&greens_magnetic, &alpha_magnetic);
  PowerSpectrumHarmOsc powerspectrum(&greens, &alpha);

  // Matrices to store results
  cx_mat::fixed<3, 3> lhs(fill::zeros);
  cx_mat::fixed<3, 3> rhs(fill::zeros);
  cx_mat::fixed<3, 3> EMPTY(fill::zeros);

  //for frequencies around the plasmon frequency and beyond there should be contributions of the
  //magnetic terms
  auto omega = GENERATE(take(5, random(10., 1e2)));

  Options_PowerSpectrum opts;
  opts.spectrum = FULL;
  opts.omega = omega;

  powerspectrum.calculate(rhs, opts);
  powerspectrum_magnetic.calculate(lhs, opts);

  REQUIRE(!approx_equal(lhs-rhs, EMPTY, "reldiff", 1e-8));
}


TEST_CASE(
    "Magnetic power spectrum reduces to polarizability in the static case",
    "[PowerSpectrumHarmOscMangetic]") {
  // Generate randomnized power spectrum
  double omega_p = 9;
  double gamma = .1;
  double v = 0.;
  double beta = 1000;
  double za = 0.1;
  double delta_cut = 100;
  vec::fixed<2> rel_err = {1E-8, 1E-6};
  PermittivityDrude perm(omega_p, gamma);
  ReflectionCoefficientsLocBulk refl(&perm);
  GreensTensorPlateMagnetic greens(v, za, beta , &refl, delta_cut, rel_err);

  double omega_a = GENERATE(take(1, random(0.1, 1e1)));
  double alpha_zero = GENERATE(take(1, random(1e-9, 1e-7)));
  PolarizabilityNoBathMagnetic alpha(omega_a, alpha_zero, &greens);

  PowerSpectrumHarmOscMagnetic powerspectrum(&greens, &alpha);

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

  std::cout << lhs << rhs << std::endl;
  REQUIRE(approx_equal(lhs, rhs, "reldiff", 10e-8));
}

