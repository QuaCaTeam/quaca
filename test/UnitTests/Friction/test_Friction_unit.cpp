#include "Quaca.h"
#include "catch.hpp"

TEST_CASE("Friction constructors work as expected", "[Friction]") {
  SECTION("Direct constructor") {
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
  }

  SECTION("json file constructor") {
    double omega_a = 1.3;
    double alpha_zero = 6e-9;

    double beta = 5.;
    double v = 0.01;

    Friction quant_fric("../data/test_files/FrictionVacuum.json");

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
  }
}
