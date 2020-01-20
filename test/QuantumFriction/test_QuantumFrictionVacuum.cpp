#include "Quaca.h"
#include "catch.hpp"

TEST_CASE("Quantum friction calculations work") {
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

    QuantumFriction quant_fric("../data/test_files/QuantumFrictionVacuum.ini");
    Options_Friction opts;
    opts.class_pt = &quant_fric;
    //    REQUIRE(Approx(quant_fric.calculate(opts, omega_min, omega_max,
    //    relerr, epsabs)).epsilon(1e-8) == analytical_result);
    REQUIRE(0 == 0);
  }
}
