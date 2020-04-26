#include "Quaca.h"
#include "catch.hpp"
#include <armadillo>
#include <complex>

TEST_CASE("Vacuum Green's tensor constructors work as expected",
          "[GreensTensorVacuum]") {

  SECTION("Direct constructor") {
    auto v = GENERATE(take(2, random(0., 1.)));
    auto beta = GENERATE(take(2, random(1e-3, 1e3)));
    auto relerr = GENERATE(take(2, random(1e-9, 1e-1)));

    GreensTensorVacuum Greens(v, beta, relerr);

    REQUIRE(Greens.get_v() == v);
    REQUIRE(Greens.get_beta() == beta);
    REQUIRE(Greens.get_relerr() == relerr);
  }

  SECTION("json file constructor") {
    double v = 0.1;
    double beta = 5;
    double relerr = 1E-9;

    GreensTensorVacuum Greens("../data/test_files/GreensTensorVacuum.json");

    REQUIRE(Approx(Greens.get_v()).epsilon(1E-6) == v);
    REQUIRE(Approx(Greens.get_beta()).epsilon(1E-6) == beta);
    REQUIRE(Greens.get_relerr() == relerr);
  }
}

TEST_CASE("Integrand 1d k is correctly implemented", "[GreensTensorVacuum]") {
  // Generate a Green's tensor with random attributes v and beta
  auto v = GENERATE(take(3, random(1e-10, 1.)));
  auto beta = GENERATE(take(3, random(1., 8e1)));
  double relerr = 1E-9;
  GreensTensorVacuum Greens(v, beta, relerr);

  // Create the matries storing the Green's tensors
  cx_mat::fixed<3, 3> Greens_lhs(fill::zeros);
  cx_mat::fixed<3, 3> Greens_rhs(fill::zeros);

  // Create the variables for the num_result, taking care, that
  //\omega^2 - k^2 >= 0 to stay in the non-trivial regime
  auto omega = GENERATE(take(3, random(-1e1, 1e1)));
  auto k_v = GENERATE(take(3, random(-1., 1.)));
  if (k_v < 0)
    k_v *= omega / (1 + v);
  if (k_v >= 0)
    k_v *= omega / (1 - v);

  // Check the integrand for all possible integration options
  double omega_kv = omega + v * k_v;
  double xi = pow(omega_kv, 2) - pow(k_v, 2);
  double result_x = .5 * xi;
  double result_yz = .5 * (pow(omega_kv, 2) - .5 * xi);

  SECTION("Option: IM") {
    double factor = 1.;
    REQUIRE(Approx(Greens.integrand_k(k_v, omega, {0, 0}, IM, UNIT))
                .epsilon(1E-7) == result_x * factor);
    REQUIRE(Approx(Greens.integrand_k(k_v, omega, {1, 1}, IM, UNIT))
                .epsilon(1E-7) == result_yz * factor);
    REQUIRE(Approx(Greens.integrand_k(k_v, omega, {2, 2}, IM, UNIT))
                .epsilon(1E-7) == result_yz * factor);
  }

  SECTION("Option: IM, KV") {
    double factor = k_v;
    REQUIRE(
        Approx(Greens.integrand_k(k_v, omega, {0, 0}, IM, KV)).epsilon(1E-7) ==
        result_x * factor);
    REQUIRE(
        Approx(Greens.integrand_k(k_v, omega, {1, 1}, IM, KV)).epsilon(1E-7) ==
        result_yz * factor);
    REQUIRE(
        Approx(Greens.integrand_k(k_v, omega, {2, 2}, IM, KV)).epsilon(1E-7) ==
        result_yz * factor);
  }

  SECTION("Option: IM, TEMP") {
    double factor = 1. / (1. - exp(-beta * omega_kv));
    REQUIRE(Approx(Greens.integrand_k(k_v, omega, {0, 0}, IM, TEMP))
                .epsilon(1E-7) == result_x * factor);
    REQUIRE(Approx(Greens.integrand_k(k_v, omega, {1, 1}, IM, TEMP))
                .epsilon(1E-7) == result_yz * factor);
    REQUIRE(Approx(Greens.integrand_k(k_v, omega, {2, 2}, IM, TEMP))
                .epsilon(1E-7) == result_yz * factor);
  }

  SECTION("Option: IM, KV_TEMP") {
    double factor = k_v / (1. - exp(-beta * omega_kv));
    REQUIRE(Approx(Greens.integrand_k(k_v, omega, {0, 0}, IM, KV_TEMP))
                .epsilon(1E-7) == result_x * factor);
    REQUIRE(Approx(Greens.integrand_k(k_v, omega, {1, 1}, IM, KV_TEMP))
                .epsilon(1E-7) == result_yz * factor);
    REQUIRE(Approx(Greens.integrand_k(k_v, omega, {2, 2}, IM, KV_TEMP))
                .epsilon(1E-7) == result_yz * factor);
  }

  SECTION("Option: IM, NON_LTE") {
    double factor =
        1. / (1. - exp(-beta * (omega_kv))) - 1. / (1. - exp(-beta * omega));
    REQUIRE(Approx(Greens.integrand_k(k_v, omega, {0, 0}, IM, NON_LTE))
                .epsilon(1E-7) == result_x * factor);
    REQUIRE(Approx(Greens.integrand_k(k_v, omega, {1, 1}, IM, NON_LTE))
                .epsilon(1E-7) == result_yz * factor);
    REQUIRE(Approx(Greens.integrand_k(k_v, omega, {2, 2}, IM, NON_LTE))
                .epsilon(1E-7) == result_yz * factor);
  }
}

/*!
 * Some basic relations any Green's tensor should fulfill which can
 * be found in docs under: Relations_and_num_results.pdf
 */
TEST_CASE("Crossing relation in frequency domain see eq. [1]",
          "[GreensTensorVacuum]") {
  // Generate a Green's tensor with random attributes v and beta
  auto v = GENERATE(take(1, random(0., 1.)));
  auto beta = GENERATE(take(1, random(1e-5, 1e5)));
  double relerr = 1E-9;

  // Create the variables for the num_result, taking care, that
  //\omega^2 - k^2 >= 0 to stay in the non-trivial regime
  auto k_x = GENERATE(take(3, random(-1e3, 1e3)));
  auto k_y = GENERATE(take(3, random(-1e3, 1e3)));
  auto omega = GENERATE(take(3, random(1., 1e3)));
  double k = sqrt(k_x * k_x + k_y * k_y);
  omega *= k;

  GreensTensorVacuum Greens(v, beta, relerr);

  cx_mat::fixed<3, 3> Greens_lhs(fill::zeros);
  cx_mat::fixed<3, 3> Greens_rhs(fill::zeros);

  Greens.calculate_tensor(omega, {k_x, k_y}, Greens_lhs);
  Greens.calculate_tensor(-omega, {-k_x, -k_y}, Greens_rhs);

  REQUIRE(approx_equal(Greens_lhs, trans(conj(Greens_rhs)), "reldiff", 10E-5));
}

TEST_CASE("Reciprocity, see eq. [6]", "[GreensTensorVacuum]") {
  // Generate a Green's tensor with random attributes v and beta
  auto v = GENERATE(take(1, random(0., 1.)));
  auto beta = GENERATE(take(1, random(1e-5, 1e5)));
  double relerr = 1E-9;

  auto omega = GENERATE(take(5, random(1., 1e3)));
  auto k_x = GENERATE(take(5, random(0.0, 1e3)));
  auto k_y = GENERATE(take(5, random(0.0, 1e3)));

  // Take care that we are looking at the non trivial part of the Green's tensor
  // where \omega^2 - k^2 >= 0
  double k = sqrt(k_x * k_x + k_y * k_y);
  omega *= k;

  GreensTensorVacuum Greens(v, beta, relerr);

  cx_mat::fixed<3, 3> Greens_lhs(fill::zeros);
  cx_mat::fixed<3, 3> Greens_rhs(fill::zeros);

  Greens.calculate_tensor(omega, {k_x, k_y}, Greens_lhs);
  Greens.calculate_tensor(omega, {-k_x, -k_y}, Greens_rhs);

  REQUIRE(approx_equal(Greens_lhs, trans(Greens_rhs), "reldiff", 10E-5));
}

TEST_CASE("Reality, see eq. [7]", "[GreensTensorVacuum]") {
  // Generate a Green's tensor with random attributes v and beta
  auto v = GENERATE(take(1, random(0., 1.)));
  auto beta = GENERATE(take(1, random(1e-5, 1e5)));
  auto omega = GENERATE(take(5, random(-1e3, 1e3)));
  double relerr = 1E-9;

  GreensTensorVacuum Greens(v, beta, relerr);

  cx_mat::fixed<3, 3> Greens_lhs(fill::zeros);
  cx_mat::fixed<3, 3> Greens_rhs(fill::zeros);

  Greens.integrate_k(omega, Greens_lhs, IM, UNIT);
  Greens.integrate_k(-omega, Greens_rhs, IM, UNIT);

  REQUIRE(approx_equal(Greens_lhs, -Greens_rhs, "reldiff", 10E-5));
}

TEST_CASE("Test the integration routine", "[GreensTensorVacuum]") {

  SECTION("Option: IM") {
    // Generate a Green's tensor with random attributes v and beta
    auto v = GENERATE(take(1, random(0., 1.)));
    auto beta = GENERATE(take(1, random(1e-5, 1e5)));
    auto omega = GENERATE(take(1, random(-1e2, 1e2)));
    double relerr = 1E-9;
    GreensTensorVacuum Greens(v, beta, relerr);

    // Matrix to store the analytic results
    cx_mat::fixed<3, 3> ana_result(fill::zeros);
    // Computing the analytical result and storing it in analytic
    double ana_pref = 2. / 3. * pow(omega, 3) / pow(1 - pow(v, 2), 2);
    ana_result(0, 0) = ana_pref;
    ana_result(1, 1) = ana_pref * (1 + pow(v, 2)) / (1 - pow(v, 2));
    ana_result(2, 2) = ana_pref * (1 + pow(v, 2)) / (1 - pow(v, 2));

    // Matrix storing the numerical integration
    cx_mat::fixed<3, 3> num_result(fill::zeros);
    Greens.integrate_k(omega, num_result, IM, UNIT);

    REQUIRE(approx_equal(num_result, ana_result, "reldiff", 10E-5));
  }

  SECTION("Option: IM, KV") {

    // Generate a Green's tensor with random attributes v and beta
    auto v = GENERATE(take(1, random(0., 1.)));
    auto beta = GENERATE(take(1, random(1e-5, 1e5)));
    auto omega = GENERATE(take(1, random(-1e2, 1e2)));
    double relerr = 1E-9;
    GreensTensorVacuum Greens(v, beta, relerr);

    // Matrix to store the analytic results
    cx_mat::fixed<3, 3> ana_result(fill::zeros);
    // Computing the analytical result and storing it in analytic
    double ana_pref = 2. / 3. * pow(omega, 4) * v / pow(1 - pow(v, 2), 3);
    ana_result(0.0) = ana_pref;
    ana_result(1, 1) = ana_pref * (2. + pow(v, 2)) / (1. - pow(v, 2));
    ana_result(2, 2) = ana_result(1, 1);

    // Matrix storing the numerical integration
    cx_mat::fixed<3, 3> num_result(fill::zeros);
    Greens.integrate_k(omega, num_result, IM, KV);

    REQUIRE(approx_equal(num_result, ana_result, "reldiff", 10E-5));
  }

  SECTION("Option: IM, TEMP") {
    // Generate a Green's tensor with random attributes v and beta
    auto v = GENERATE(take(1, random(0., 1.)));
    auto omega = GENERATE(take(3, random(-1e2, 1e2)));
    auto beta = GENERATE(take(3, random(10e-10, 10e-12)));
    beta *= fabs(omega);
    double relerr = 1E-9;
    GreensTensorVacuum Greens(v, beta, relerr);

    // Matrix to store the analytic results
    cx_mat::fixed<3, 3> ana_result(fill::zeros);
    // Computing the analytical result and storing it in analytic
    double ana_pref = pow(omega, 2) / (2. * pow(v, 3) * beta);
    ana_result(0.0) = ana_pref * (2. * v / (1. - pow(v, 2)) - 2. * atanh(v));
    ana_result(1, 1) =
        ana_pref * ((3. * pow(v, 3) - v) / pow(1 - pow(v, 2), 2) + atanh(v));
    ana_result(2, 2) = ana_result(1, 1);

    // Matrix storing the numerical integration
    cx_mat::fixed<3, 3> num_result(fill::zeros);
    Greens.integrate_k(omega, num_result, IM, TEMP);

    REQUIRE(approx_equal(num_result, ana_result, "reldiff", 10E-4));
  }

  SECTION("Option: IM, KV_TEMP") {
    // Generate a Green's tensor with random attributes v and beta
    auto v = GENERATE(take(1, random(0., 1.)));
    auto omega = GENERATE(take(3, random(-1e1, 1e1)));
    auto beta = GENERATE(take(3, random(10e-10, 10e-12)));
    //  beta *= fabs(omega);
    double relerr = 1E-9;
    GreensTensorVacuum Greens(v, beta, relerr);

    // Matrix to store the analytic results
    cx_mat::fixed<3, 3> ana_result(fill::zeros);
    // Computing the analytical result and storing it in analytic
    double ana_pref = std::pow(omega, 3) / (6. * std::pow(v, 4) * beta);
    ana_result(0, 0) =
        ana_pref * (2 * (5 * pow(v, 3) - 3. * v) / pow(1. - pow(v, 2), 2) +
                    6. * std::atanh(v));
    ana_result(1, 1) = ana_pref * ((8. * pow(v, 3) - 3. * v - 13. * pow(v, 5)) /
                                       pow(pow(v, 2) - 1, 3) -
                                   3. * std::atanh(v));
    ana_result(2, 2) = ana_result(1, 1);

    // Matrix storing the numerical integration
    cx_mat::fixed<3, 3> num_result(fill::zeros);
    Greens.integrate_k(omega, num_result, IM, KV_TEMP);

    REQUIRE(approx_equal(num_result, ana_result, "reldiff", 10E-5));
  }
}
