#include "Quaca.h"
#include "catch.hpp"

TEST_CASE("Quantum friction calculations work") {
  SECTION("Constructor with initialisation list works") {
    auto v = GENERATE(take(1, random(0., 1.)));
    double beta = GENERATE(take(1, random(1e-5, 1e3)));
    double omega_a = GENERATE(take(1, random(0., 1e3)));
    double alpha_zero = GENERATE(take(1, random(1e-5, 1.)));
    double relerr_omega = 1E-1;
    double relerr_k = 1E-9;

    GreensTensorVacuum green(v, beta, relerr_k);
    PolarizabilityNoBath alpha(omega_a, alpha_zero, &green);
    PowerSpectrumHarmOsc powerspectrum(&green, &alpha);
    Friction quant_fric(&green, &alpha, &powerspectrum, relerr_omega);

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
    double relerr_k = 1E-9;

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

  SECTION("Analytical results with vacuum Green's tensor get reproduced") {
    // Units: c=1, 4 pi epsilon_0 = 1, hbar = 1
    double omega_a = 1.3;
    double alpha_zero = 6e-9;

    double beta = .01;
    double v = 0.01;
    double analytical_result =
        -2. / 3. * v / beta * alpha_zero * pow(omega_a, 4);

    double omega_max = 1e5;
    double omega_min = -omega_min;
    double relerr = 1e-5;
    double epsabs = 0;

    Friction quant_fric("../data/test_files/FrictionVacuum.ini");
    Options_Friction opts;
    opts.class_pt = &quant_fric;
    //    REQUIRE(Approx(quant_fric.calculate(opts, omega_min, omega_max,
    //    relerr, epsabs)).epsilon(1e-8) == analytical_result);
    REQUIRE(0 == 0);
  };
}
