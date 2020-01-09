#include "catch.hpp"
#include "Quaca.h"
#include <iostream>


TEST_CASE("Power spectrum of harmonic oscillator works properly")
{
  
  SECTION("Constructor with initialisation list works")
  {
    auto v = GENERATE(take(1, random(0.,1.)));
    double beta = GENERATE(take(1, random(1e-5, 1e3)));
    double omega_a = GENERATE(take(1, random(0., 1e3)));
    double alpha_zero = GENERATE(take(1, random(1e-5, 1.)));

    GreensTensorVacuum green(v, beta);
    PolarizabilityNoBath alpha(omega_a, alpha_zero, &green);
    PowerSpectrumHarmOsc powerspectrum(&green, &alpha);

    REQUIRE(Approx(powerspectrum.greens_tensor->get_v()).epsilon(1e-6) == v);
    REQUIRE(Approx(powerspectrum.greens_tensor->get_beta()).epsilon(1e-6) == beta);
    REQUIRE(Approx(powerspectrum.polarizability->get_omega_a()).epsilon(1e-6) == omega_a);
    REQUIRE(Approx(powerspectrum.polarizability->get_alpha_zero()).epsilon(1e-6) == alpha_zero);
  }

  SECTION("Constructor with ini file works")
  {
    double omega_a = 1.3;
    double alpha_zero = 6e-9;

    double beta = 5.;
    double v = 0.01;

    PowerSpectrumHarmOsc powerspectrum("../data/test_files/PowerSpectrumHarmOsc.ini");

    REQUIRE(Approx(powerspectrum.greens_tensor->get_v()).epsilon(1e-6) == v);
    REQUIRE(Approx(powerspectrum.greens_tensor->get_beta()).epsilon(1e-6) == beta);
    REQUIRE(Approx(powerspectrum.polarizability->get_omega_a()).epsilon(1e-6) == omega_a);
    REQUIRE(Approx(powerspectrum.polarizability->get_alpha_zero()).epsilon(1e-6) == alpha_zero);
  }

  

  SECTION("Power spectrum is hermitian")
  {
    auto omega = GENERATE(take(5,random(-1e3,1e3)));

    PowerSpectrumHarmOsc powerspectrum("../data/test_files/PowerSpectrumHarmOsc.ini");

    cx_mat::fixed<3,3> lhs(fill::zeros);
    cx_mat::fixed<3,3> rhs(fill::zeros);
    powerspectrum.calculate(lhs, omega);
    powerspectrum.calculate(rhs, omega);
    cout << "Power spectrum" <<lhs << rhs << endl;

    REQUIRE(approx_equal(lhs, trans(conj(rhs)), "absdiff", 1e-5));
  }
  
  SECTION("Power spectrum reduces to polarizability in the static case")
  {
    PowerSpectrumHarmOsc powerspectrum("../data/test_files/PowerSpectrumHarmOsc.ini");
    powerspectrum.greens_tensor->set_v(0);
    PolarizabilityNoBath alpha("../data/test_files/PowerSpectrumHarmOsc.ini");
    cx_mat::fixed<3,3> lhs(fill::zeros);
    cx_mat::fixed<3,3> rhs(fill::zeros);
    auto omega = GENERATE(take(5, random(-1e3,1e3)));
    powerspectrum.calculate(lhs, omega);
    alpha.calculate(rhs, omega);
    cout << "Power spectrum" <<lhs << rhs << endl;
    REQUIRE(approx_equal(lhs, rhs, "absdiff", 10e-5));
  }
}

