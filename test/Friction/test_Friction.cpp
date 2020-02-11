#include "Quaca.h"
#include "catch.hpp"
#include <iomanip>
#include <iostream>

TEST_CASE("Quantum friction constructors work", "[Friction]") {
  SECTION("Constructor with initialisation list works") {
    auto v = GENERATE(take(1, random(0., 1.)));
    double beta = GENERATE(take(1, random(1e-5, 1e3)));
    double omega_a = GENERATE(take(1, random(0., 1e3)));
    double alpha_zero = GENERATE(take(1, random(1e-5, 1.)));
    double relerr_omega = 1E-1;

    double relerr_k = 1E-9;
    GreensTensorVacuum greens(v, beta, relerr_k);
    PolarizabilityNoBath alpha(omega_a, alpha_zero, &greens);
    PowerSpectrumHarmOsc powerspectrum(&greens, &alpha);
    Friction quant_fric(&greens, &alpha, &powerspectrum, relerr_omega);

    REQUIRE(Approx(quant_fric.get_greens_tensor()->get_v()).epsilon(1e-6) == v);
    REQUIRE(Approx(quant_fric.get_greens_tensor()->get_beta()).epsilon(1e-6) ==
            beta);
    REQUIRE(
        Approx(quant_fric.get_polarizability()->get_omega_a()).epsilon(1e-6) ==
        omega_a);
    REQUIRE(Approx(quant_fric.get_polarizability()->get_alpha_zero())
                .epsilon(1e-6) == alpha_zero);
    REQUIRE(Approx(quant_fric.get_powerspectrum()->get_greens_tensor()->get_v())
                .epsilon(1e-6) == v);
    REQUIRE(
        Approx(quant_fric.get_powerspectrum()->get_greens_tensor()->get_beta())
            .epsilon(1e-6) == beta);
    REQUIRE(
        Approx(
            quant_fric.get_powerspectrum()->get_polarizability()->get_omega_a())
            .epsilon(1e-6) == omega_a);
    REQUIRE(Approx(quant_fric.get_powerspectrum()
                       ->get_polarizability()
                       ->get_alpha_zero())
                .epsilon(1e-6) == alpha_zero);
  };

  SECTION("Constructor with ini file works") {
    double omega_a = 1.3;
    double alpha_zero = 6e-9;

    double beta = 5.;
    double v = 0.01;

    Friction quant_fric("../data/test_files/FrictionVacuum.ini");

    REQUIRE(Approx(quant_fric.get_greens_tensor()->get_v()).epsilon(1e-6) == v);
    REQUIRE(Approx(quant_fric.get_greens_tensor()->get_beta()).epsilon(1e-6) ==
            beta);
    REQUIRE(
        Approx(quant_fric.get_polarizability()->get_omega_a()).epsilon(1e-6) ==
        omega_a);
    REQUIRE(Approx(quant_fric.get_polarizability()->get_alpha_zero())
                .epsilon(1e-6) == alpha_zero);
    REQUIRE(Approx(quant_fric.get_powerspectrum()->get_greens_tensor()->get_v())
                .epsilon(1e-6) == v);
    REQUIRE(
        Approx(quant_fric.get_powerspectrum()->get_greens_tensor()->get_beta())
            .epsilon(1e-6) == beta);
    REQUIRE(
        Approx(
            quant_fric.get_powerspectrum()->get_polarizability()->get_omega_a())
            .epsilon(1e-6) == omega_a);
    REQUIRE(Approx(quant_fric.get_powerspectrum()
                       ->get_polarizability()
                       ->get_alpha_zero())
                .epsilon(1e-6) == alpha_zero);
  };
};

TEST_CASE("Analytical results with vacuum Green's tensor gets reproduced",
          "[Friction]") {
  // Units: c=1, 4 pi epsilon_0 = 1, hbar = 1
  double omega_a = .3;
  double alpha_zero = 6e-9;

  double beta = 1e-1;
  double v = 1e-5;
  double analytical_result = -2. / 12. * v * alpha_zero * pow(omega_a, 6) *
                             beta / pow(sinh(omega_a * beta / 2.), 2);
  // double analytical_result = (3./2.)*alpha_zero*pow(omega_a,2)/beta;

  double relerr_omega = 1e-6;
  double epsabs = 0;

  double relerr_k = 1E-9;
  GreensTensorVacuum greens(v, beta, relerr_k);
  PolarizabilityNoBath alpha(omega_a, alpha_zero, &greens);
  PowerSpectrumHarmOsc powerspectrum(&greens, &alpha);
  Friction quant_fric(&greens, &alpha, &powerspectrum, relerr_omega);

  Options_Friction opts;
  opts.class_pt = &quant_fric;
  opts.non_LTE = true;
  // opts.full_spectrum = true;
  double num_result = quant_fric.calculate(opts);
  REQUIRE(Approx(num_result).epsilon(1e-4) == analytical_result);
};

TEST_CASE("Analytical results with scattered Green's tensor gets reproduced",
          "[Friction]") {
  // Units: c=1, 4 pi epsilon_0 = 1, hbar = 1
  double omega_a = 1.3;
  double alpha_zero = 6e-9;
  double za = 0.01;
  double omega_p = 9.;
  double gamma = 0.1;
  double rho;
  rho = gamma * M_PI * 4. / pow(omega_p, 2);
  double beta = 1e6;
  double v = 1e-4;
  double delta_cut = 30;
  double analytical_result = -(63. - 45.) * pow(alpha_zero * rho, 2) *
                             pow(v / M_PI, 3) / pow(2 * za, 10);
  analytical_result += -(6. - 3.) * pow(alpha_zero * rho / beta, 2) *
                       (v / M_PI) / pow(2 * za, 8);

  vec::fixed<2> rel_err = {1E-6, 1E-4};
  double relerr_omega = 1e-2;
  double epsabs = 0;

  PermittivityDrude perm(omega_p, gamma);
  ReflectionCoefficientsLocBulk refl(&perm);
  GreensTensorPlate greens(v, za, beta, &refl, delta_cut, rel_err);
  PolarizabilityNoBath alpha(omega_a, alpha_zero, &greens);
  PowerSpectrumHarmOsc powerspectrum(&greens, &alpha);
  Friction quant_fric(&greens, &alpha, &powerspectrum, relerr_omega);

  Options_Friction opts;
  opts.class_pt = &quant_fric;
  opts.non_LTE = true;
  // opts.full_spectrum = true;
  double num_result = quant_fric.calculate(opts);
  std::cout << "ana=" << analytical_result << std::endl;
  std::cout << "num=" << num_result << std::endl;
  std::cout << "num/ana=" << num_result / analytical_result << std::endl;
  REQUIRE(Approx(num_result).epsilon(1e-2) == analytical_result);
};
