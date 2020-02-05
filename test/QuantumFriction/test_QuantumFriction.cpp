#include "Quaca.h"
#include "catch.hpp"
#include <iomanip>
#include <iostream>

TEST_CASE("Quantum friction constructors work", "[QuantumFriction]") {
  SECTION("Constructor with initialisation list works") {
    auto v = GENERATE(take(1, random(0., 1.)));
    double beta = GENERATE(take(1, random(1e-5, 1e3)));
    double omega_a = GENERATE(take(1, random(0., 1e3)));
    double alpha_zero = GENERATE(take(1, random(1e-5, 1.)));

    GreensTensorVacuum green(v, beta);
    PolarizabilityNoBath alpha(omega_a, alpha_zero, &green);
    PowerSpectrumHarmOsc powerspectrum(&green, &alpha);
    QuantumFriction quant_fric(&green, &alpha, &powerspectrum);

    REQUIRE(Approx(quant_fric.greens_tensor->get_v()).epsilon(1e-6) == v);
    REQUIRE(Approx(quant_fric.greens_tensor->get_beta()).epsilon(1e-6) == beta);
    REQUIRE(Approx(quant_fric.polarizability->get_omega_a()).epsilon(1e-6) ==
            omega_a);
    REQUIRE(Approx(quant_fric.polarizability->get_alpha_zero()).epsilon(1e-6) ==
            alpha_zero);
    REQUIRE(Approx(quant_fric.powerspectrum->greens_tensor->get_v())
                .epsilon(1e-6) == v);
    REQUIRE(Approx(quant_fric.powerspectrum->greens_tensor->get_beta())
                .epsilon(1e-6) == beta);
    REQUIRE(Approx(quant_fric.powerspectrum->polarizability->get_omega_a())
                .epsilon(1e-6) == omega_a);
    REQUIRE(Approx(quant_fric.powerspectrum->polarizability->get_alpha_zero())
                .epsilon(1e-6) == alpha_zero);
  }

  SECTION("Constructor with ini file works") {
    double omega_a = 1.3;
    double alpha_zero = 6e-9;

    double beta = 5.;
    double v = 0.01;

    QuantumFriction quant_fric("../data/test_files/QuantumFrictionVacuum.ini");

    REQUIRE(Approx(quant_fric.greens_tensor->get_v()).epsilon(1e-6) == v);
    REQUIRE(Approx(quant_fric.greens_tensor->get_beta()).epsilon(1e-6) == beta);
    REQUIRE(Approx(quant_fric.polarizability->get_omega_a()).epsilon(1e-6) ==
            omega_a);
    REQUIRE(Approx(quant_fric.polarizability->get_alpha_zero()).epsilon(1e-6) ==
            alpha_zero);
    REQUIRE(Approx(quant_fric.powerspectrum->greens_tensor->get_v())
                .epsilon(1e-6) == v);
    REQUIRE(Approx(quant_fric.powerspectrum->greens_tensor->get_beta())
                .epsilon(1e-6) == beta);
    REQUIRE(Approx(quant_fric.powerspectrum->polarizability->get_omega_a())
                .epsilon(1e-6) == omega_a);
    REQUIRE(Approx(quant_fric.powerspectrum->polarizability->get_alpha_zero())
                .epsilon(1e-6) == alpha_zero);
  }
}

TEST_CASE("Analytical results with vacuum Green's tensor gets reproduced",
          "[QuantumFriction]") {
  // Units: c=1, 4 pi epsilon_0 = 1, hbar = 1
  double omega_a = .3;
  double alpha_zero = 6e-9;

  double beta = 1e-1;
  double v = 1e-5;
  double analytical_result = -2. / 12. * v * alpha_zero * pow(omega_a, 6) *
                             beta / pow(sinh(omega_a * beta / 2.), 2);
  // double analytical_result = (3./2.)*alpha_zero*pow(omega_a,2)/beta;

  double relerr = 1e-6;
  double epsabs = 0;

  GreensTensorVacuum green(v, beta);
  PolarizabilityNoBath alpha(omega_a, alpha_zero, &green);
  PowerSpectrumHarmOsc powerspectrum(&green, &alpha);
  QuantumFriction quant_fric(&green, &alpha, &powerspectrum);

  Options_Friction opts;
  opts.class_pt = &quant_fric;
  opts.non_LTE = true;
  // opts.full_spectrum = true;
  double num_result = quant_fric.calculate(opts, relerr, epsabs);
  REQUIRE(Approx(num_result).epsilon(1e-4) == analytical_result);
}
TEST_CASE("Analytical results with scattered Green's tensor gets reproduced",
          "[QuantumFriction]") {
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

  PermittivityDrude perm(gamma, omega_p);
  ReflectionCoefficientsLocBulk refl(&perm); 
  GreensTensorPlate green(v, za, beta, &refl, delta_cut, rel_err);
  PolarizabilityNoBath alpha(omega_a, alpha_zero, &green);
  PowerSpectrumHarmOsc powerspectrum(&green, &alpha);
  QuantumFriction quant_fric(&green, &alpha, &powerspectrum);

  Options_Friction opts;
  opts.class_pt = &quant_fric;
  opts.non_LTE = true;
  // opts.full_spectrum = true;
  double num_result = quant_fric.calculate(opts, relerr_omega, epsabs);
  std::cout << "ana=" << analytical_result << std::endl;
  std::cout << "num=" << num_result << std::endl;
  std::cout << "num/ana=" << num_result / analytical_result << std::endl;
  REQUIRE(Approx(num_result).epsilon(1e-2) == analytical_result);
}
