#include <armadillo>
#include <complex>

#include "Quaca.h"
#include "catch.hpp"

TEST_CASE("Vacuum Greens Tensor works properly") {

  SECTION("Constructor with argument list works") {
    auto v = GENERATE(take(2, random(0., 1.)));
    auto beta = GENERATE(take(2, random(1e-3, 1e3)));
    GreensTensorVacuum Greens(v, beta);

    REQUIRE(Approx(Greens.get_v()).epsilon(1E-6) == v);
    REQUIRE(Approx(Greens.get_beta()).epsilon(1E-6) == beta);
  }

  SECTION("Constructor with ini file works") {
    double v = 0.1;
    double beta = 5;

    GreensTensorVacuum Greens("../data/test_files/GreensTensorVacuum.ini");

    REQUIRE(Approx(Greens.get_v()).epsilon(1E-6) == v);
    REQUIRE(Approx(Greens.get_beta()).epsilon(1E-6) == beta);
  }

  /*!
   * Some basic relations any Green's tensor should fulfill which can
   * be found in docs under: Relations_and_tests.pdf
   */

  GreensTensorVacuum Greens(0.01, 1e3);
  struct Options_GreensTensor opts;
  opts.fancy_I = true;
  opts.class_pt = &Greens;
  cx_mat::fixed<3, 3> Greens_lhs(fill::zeros);
  cx_mat::fixed<3, 3> Greens_rhs(fill::zeros);

  SECTION("Crossing relation in frequency domain see eq. [1]") {
    auto omega = GENERATE(take(5, random(0.0, 1e3)));
    auto k_x = GENERATE(take(5, random(0.0, 1e3)));
    auto k_y = GENERATE(take(5, random(0.0, 1e3)));

    opts.omega = -omega;
    opts.kvec(0) = -k_x;
    opts.kvec(1) = -k_y;
    Greens.calculate_tensor(Greens_lhs, opts);

    opts.omega = omega;
    opts.kvec(0) = k_x;
    opts.kvec(1) = k_y;
    Greens.calculate_tensor(Greens_rhs, opts);

    REQUIRE(
        approx_equal(Greens_lhs, trans(conj(Greens_rhs)), "absdiff", 10E-5));
  }

  SECTION("Reciprocity, see eq. [6]") {
    auto omega = GENERATE(take(5, random(0.0, 1e3)));
    auto k_x = GENERATE(take(5, random(0.0, 1e3)));
    auto k_y = GENERATE(take(5, random(0.0, 1e3)));

    opts.omega = omega;
    opts.kvec(0) = -k_x;
    opts.kvec(1) = -k_y;
    Greens.calculate_tensor(Greens_lhs, opts);

    opts.omega = omega;
    opts.kvec(0) = k_x;
    opts.kvec(1) = k_y;
    Greens.calculate_tensor(Greens_rhs, opts);

    REQUIRE(approx_equal(Greens_lhs, trans(Greens_rhs), "absdiff", 10E-5));
  }

  SECTION("Reality, see eq. [7]") {
    /*!
     * Due to the divergence of the real part only the imaginary
     * part is calculated in any case
     */

    auto omega = GENERATE(take(5, random(0.0, 1e3)));

    opts.omega = omega;
    Greens.integrate_1d_k(Greens_lhs, opts);

    opts.omega = -omega;
    Greens.integrate_1d_k(Greens_rhs, opts);

    REQUIRE(approx_equal(Greens_lhs, -Greens_rhs, "absdiff", 10E-5));
  }
  SECTION("Test the integration routine") {

    /*!
     * Test the integration routine of the vacuums Green's tensor, the
     * analytical results can be found in the docs under: VacuumGreen.pdf, see
     * eq. [21]
     */

    cx_mat::fixed<3, 3> test(fill::zeros);

    auto omega = GENERATE(take(5, random(0.0, 1e3)));
    opts.omega = omega;
    Greens.integrate_1d_k(test, opts);

    double analytic_result_prefactor =
        2. / 3. * pow(opts.omega, 3) / pow(1 - pow(Greens.get_v(), 2), 2);

    REQUIRE(Approx(test(0, 0).real()).epsilon(1E-6) ==
            analytic_result_prefactor);
    REQUIRE(Approx(test(1, 1).real()).epsilon(1E-6) ==
            analytic_result_prefactor * (1 + pow(Greens.get_v(), 2)) /
                (1 - pow(Greens.get_v(), 2)));
    REQUIRE(Approx(test(2, 2).real()).epsilon(1E-6) ==
            analytic_result_prefactor * (1 + pow(Greens.get_v(), 2)) /
                (1 - pow(Greens.get_v(), 2)));
  }
};
