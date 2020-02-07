#include "Quaca.h"
#include "catch.hpp"
#include <armadillo>
#include <complex>

TEST_CASE("Constructors work properly", "[GreensTensorVacuum]") {

  SECTION("Constructor with argument list works") {
    auto v = GENERATE(take(2, random(0., 1.)));
    auto beta = GENERATE(take(2, random(1e-3, 1e3)));
    auto relerr = GENERATE(take(2, random(1e-9, 1e-1)));

    GreensTensorVacuum Greens(v, beta, relerr);

    REQUIRE(Greens.get_v() == v);
    REQUIRE(Greens.get_beta() == beta);
    REQUIRE(Greens.get_relerr() == relerr);
  }

  SECTION("Constructor with ini file works") {
    double v = 0.1;
    double beta = 5;
    double relerr = 1E-9;

    GreensTensorVacuum Greens("../data/test_files/GreensTensorVacuum.ini");

    REQUIRE(Approx(Greens.get_v()).epsilon(1E-6) == v);
    REQUIRE(Approx(Greens.get_beta()).epsilon(1E-6) == beta);
    REQUIRE(Greens.get_relerr() == relerr);
  };
};

TEST_CASE("Integrand 1d k is correctly implemented", "[GreensTensorVacuum]") {
  // Generate a Green's tensor with random attributes v and beta
  auto v = GENERATE(take(1, random(0., 1.)));
  auto beta = GENERATE(take(1, random(1., 1e2)));
  double relerr = 1E-9;
  GreensTensorVacuum Greens(v, beta, relerr);

  // Create the matries storing the Green's tensors
  cx_mat::fixed<3, 3> Greens_lhs(fill::zeros);
  cx_mat::fixed<3, 3> Greens_rhs(fill::zeros);

  // Create the variables for the num_result, taking care, that
  //\omega^2 - k^2 >= 0 to stay in the non-trivial regime
  auto omega = GENERATE(take(10, random(-1e1, 1e1)));
  auto k_v = GENERATE(take(1, random(-1., 1.)));
  if (k_v < 0)
    k_v *= omega / (1 + v);
  if (k_v >= 0)
    k_v *= omega / (1 - v);

  // Check the integrand for all possible integration options
  double omega_kv = omega + v * k_v;
  double xi = pow(omega_kv, 2) - pow(k_v, 2);
  double result_x = .5 * xi;
  double result_yz = .5 * (pow(omega_kv, 2) - .5 * xi);
  double factor = 1.;

  // Create a struct with the integration options
  struct Options_GreensTensor opts;
  opts.class_pt = &Greens;
  opts.omega = omega;

  SECTION("Option: fancy_I") {
    opts.fancy_I = true;
    opts.indices = {0, 0};
    REQUIRE(Approx(Greens.integrand_1d_k(k_v, &opts)).epsilon(1E-7) ==
            result_x * factor);
    opts.indices = {1, 1};
    REQUIRE(Approx(Greens.integrand_1d_k(k_v, &opts)).epsilon(1E-7) ==
            result_yz * factor);
    opts.indices = {2, 2};
    REQUIRE(Approx(Greens.integrand_1d_k(k_v, &opts)).epsilon(1E-7) ==
            result_yz * factor);
  };

  SECTION("Option: fancy_I_kv") {
    opts.fancy_I = false;
    opts.fancy_I_kv = true;
    factor = k_v;
    opts.indices = {0, 0};
    REQUIRE(Approx(Greens.integrand_1d_k(k_v, &opts)).epsilon(1E-7) ==
            result_x * factor);
    opts.indices = {1, 1};
    REQUIRE(Approx(Greens.integrand_1d_k(k_v, &opts)).epsilon(1E-7) ==
            result_yz * factor);
    opts.indices = {2, 2};
    REQUIRE(Approx(Greens.integrand_1d_k(k_v, &opts)).epsilon(1E-7) ==
            result_yz * factor);
  };
  SECTION("Option: fancy_I_temp") {
    opts.fancy_I_kv = false;
    opts.fancy_I_temp = true;
    factor = 1. / (1. - exp(-beta * omega_kv));
    opts.indices = {0, 0};
    REQUIRE(Approx(Greens.integrand_1d_k(k_v, &opts)).epsilon(1E-7) ==
            result_x * factor);
    opts.indices = {1, 1};
    REQUIRE(Approx(Greens.integrand_1d_k(k_v, &opts)).epsilon(1E-7) ==
            result_yz * factor);
    opts.indices = {2, 2};
    REQUIRE(Approx(Greens.integrand_1d_k(k_v, &opts)).epsilon(1E-7) ==
            result_yz * factor);
  };

  SECTION("Option: fancy_I_kv_temp") {
    opts.fancy_I_temp = false;
    opts.fancy_I_kv_temp = true;
    factor = k_v / (1. - exp(-beta * omega_kv));
    opts.indices = {0, 0};
    REQUIRE(Approx(Greens.integrand_1d_k(k_v, &opts)).epsilon(1E-7) ==
            result_x * factor);
    opts.indices = {1, 1};
    REQUIRE(Approx(Greens.integrand_1d_k(k_v, &opts)).epsilon(1E-7) ==
            result_yz * factor);
    opts.indices = {2, 2};
    REQUIRE(Approx(Greens.integrand_1d_k(k_v, &opts)).epsilon(1E-7) ==
            result_yz * factor);
  };

  SECTION("Option: fancy_I_non_LTE") {
    opts.fancy_I_kv_temp = false;
    opts.fancy_I_non_LTE = true;
    factor =
        1. / (1. - exp(-beta * (omega_kv))) - 1. / (1. - exp(-beta * omega));
    opts.indices = {0, 0};
    REQUIRE(Approx(Greens.integrand_1d_k(k_v, &opts)).epsilon(1E-7) ==
            result_x * factor);
    opts.indices = {1, 1};
    REQUIRE(Approx(Greens.integrand_1d_k(k_v, &opts)).epsilon(1E-7) ==
            result_yz * factor);
    opts.indices = {2, 2};
    REQUIRE(Approx(Greens.integrand_1d_k(k_v, &opts)).epsilon(1E-7) ==
            result_yz * factor);
  };
};

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
  GreensTensorVacuum Greens(v, beta, relerr);
  // Create a struct with the integration options
  struct Options_GreensTensor opts;
  opts.fancy_I = true;
  opts.class_pt = &Greens;

  // Create the matries storing the Green's tensors
  cx_mat::fixed<3, 3> Greens_lhs(fill::zeros);
  cx_mat::fixed<3, 3> Greens_rhs(fill::zeros);

  // Create the variables for the num_result, taking care, that
  //\omega^2 - k^2 >= 0 to stay in the non-trivial regime
  auto k_x = GENERATE(take(3, random(-1e3, 1e3)));
  auto k_y = GENERATE(take(3, random(-1e3, 1e3)));
  auto omega = GENERATE(take(3, random(1., 1e3)));
  double k = sqrt(k_x * k_x + k_y * k_y);
  omega *= k;

  // Compute the respective Green's tensors
  opts.omega = -omega;
  opts.kvec(0) = -k_x;
  opts.kvec(1) = -k_y;
  Greens.calculate_tensor(Greens_lhs, opts);

  opts.omega = omega;
  opts.kvec(0) = k_x;
  opts.kvec(1) = k_y;
  Greens.calculate_tensor(Greens_rhs, opts);

  REQUIRE(approx_equal(Greens_lhs, trans(conj(Greens_rhs)), "reldiff", 10E-5));
};

TEST_CASE("Reciprocity, see eq. [6]", "[GreensTensorVacuum]") {
  // Generate a Green's tensor with random attributes v and beta
  auto v = GENERATE(take(1, random(0., 1.)));
  auto beta = GENERATE(take(1, random(1e-5, 1e5)));
  double relerr = 1E-9;
  GreensTensorVacuum Greens(v, beta, relerr);
  // Create a struct with the integration options
  struct Options_GreensTensor opts;
  opts.fancy_I = true;
  opts.class_pt = &Greens;

  // Create the matries storing the Green's tensors
  cx_mat::fixed<3, 3> Greens_lhs(fill::zeros);
  cx_mat::fixed<3, 3> Greens_rhs(fill::zeros);
  auto omega = GENERATE(take(5, random(1., 1e3)));
  auto k_x = GENERATE(take(5, random(0.0, 1e3)));
  auto k_y = GENERATE(take(5, random(0.0, 1e3)));

  // Take care that we are looking at the non trivial part of the Green's tensor
  // where \omega^2 - k^2 >= 0
  double k = sqrt(k_x * k_x + k_y * k_y);
  omega *= k;

  opts.omega = omega;
  opts.kvec(0) = -k_x;
  opts.kvec(1) = -k_y;
  Greens.calculate_tensor(Greens_lhs, opts);

  opts.omega = omega;
  opts.kvec(0) = k_x;
  opts.kvec(1) = k_y;
  Greens.calculate_tensor(Greens_rhs, opts);

  REQUIRE(approx_equal(Greens_lhs, trans(Greens_rhs), "reldiff", 10E-5));
};

TEST_CASE("Reality, see eq. [7]", "[GreensTensorVacuum]") {
  // Generate a Green's tensor with random attributes v and beta
  auto v = GENERATE(take(1, random(0., 1.)));
  auto beta = GENERATE(take(1, random(1e-5, 1e5)));
  auto omega = GENERATE(take(5, random(-1e3, 1e3)));
  double relerr = 1E-9;
  GreensTensorVacuum Greens(v, beta, relerr);

  // Create a struct with the integration options
  Options_GreensTensor opts;
  opts.fancy_I = true;
  opts.class_pt = &Greens;

  // Create the matries storing the Green's tensors
  cx_mat::fixed<3, 3> Greens_lhs(fill::zeros);
  cx_mat::fixed<3, 3> Greens_rhs(fill::zeros);

  opts.omega = omega;
  Greens.integrate_1d_k(Greens_lhs, opts);

  opts.omega = -omega;
  Greens.integrate_1d_k(Greens_rhs, opts);

  REQUIRE(approx_equal(Greens_lhs, -Greens_rhs, "reldiff", 10E-5));
};

TEST_CASE("Test the integration routine", "[GreensTensorVacuum]") {

  SECTION("Option: fancy_I") {
    // Generate a Green's tensor with random attributes v and beta
    auto v = GENERATE(take(1, random(0., 1.)));
    auto beta = GENERATE(take(1, random(1e-5, 1e5)));
    auto omega = GENERATE(take(1, random(-1e2, 1e2)));
    double relerr = 1E-9;
    GreensTensorVacuum Greens(v, beta, relerr);

    // Create a struct with the integration options
    Options_GreensTensor opts;
    opts.class_pt = &Greens;
    opts.omega = omega;

    // Matrix storing the numerical integration
    cx_mat::fixed<3, 3> num_result(fill::zeros);

    // Matrix to store the analytic results
    cx_mat::fixed<3, 3> ana_result(fill::zeros);

    double ana_pref = 0;
    // Integration of the vacuum Green's tensor
    opts.fancy_I = true;
    Greens.integrate_1d_k(num_result, opts);

    // Computing the analytical result and storing it in analytic
    ana_pref = 2. / 3. * pow(omega, 3) / pow(1 - pow(v, 2), 2);
    ana_result(0, 0) = ana_pref;
    ana_result(1, 1) = ana_pref * (1 + pow(v, 2)) / (1 - pow(v, 2));
    ana_result(2, 2) = ana_pref * (1 + pow(v, 2)) / (1 - pow(v, 2));

    REQUIRE(approx_equal(num_result, ana_result, "reldiff", 10E-5));
  };

  SECTION("Option: fancy_I_kv") {

    // Generate a Green's tensor with random attributes v and beta
    auto v = GENERATE(take(1, random(0., 1.)));
    auto beta = GENERATE(take(1, random(1e-5, 1e5)));
    auto omega = GENERATE(take(1, random(-1e2, 1e2)));
    double relerr = 1E-9;
    GreensTensorVacuum Greens(v, beta, relerr);

    // Create a struct with the integration options
    Options_GreensTensor opts;
    opts.class_pt = &Greens;
    opts.omega = omega;

    // Matrix storing the numerical integration
    cx_mat::fixed<3, 3> num_result(fill::zeros);

    // Matrix to store the analytic results
    cx_mat::fixed<3, 3> ana_result(fill::zeros);

    double ana_pref = 0;
    // Integration of the vacuum Green's tensor
    opts.fancy_I = false;
    opts.fancy_I_kv = true;
    Greens.integrate_1d_k(num_result, opts);

    // Computing the analytical result and storing it in analytic
    ana_pref = 2. / 3. * pow(omega, 4) * v / pow(1 - pow(v, 2), 3);
    ana_result(0.0) = ana_pref;
    ana_result(1, 1) = ana_pref * (2. + pow(v, 2)) / (1. - pow(v, 2));
    ana_result(2, 2) = ana_result(1, 1);

    REQUIRE(approx_equal(num_result, ana_result, "reldiff", 10E-5));
  }

  SECTION("Option: fancy_I_temp") {
    // Generate a Green's tensor with random attributes v and beta
    auto v = GENERATE(take(1, random(0., 1.)));
    auto omega = GENERATE(take(3, random(-1e2, 1e2)));
    auto beta = GENERATE(take(3, random(10e-10, 10e-12)));
    beta *= fabs(omega);
    double relerr = 1E-9;
    GreensTensorVacuum Greens(v, beta, relerr);

    // Create a struct with the integration options
    Options_GreensTensor opts;
    opts.class_pt = &Greens;
    opts.omega = omega;

    // Matrix storing the numerical integration
    cx_mat::fixed<3, 3> num_result(fill::zeros);

    // Matrix to store the analytic results
    cx_mat::fixed<3, 3> ana_result(fill::zeros);

    double ana_pref = 0;
    // Integration of the vacuum Green's tensor
    opts.fancy_I_temp = true;
    Greens.integrate_1d_k(num_result, opts);

    // Computing the analytical result and storing it in analytic
    ana_pref = pow(omega, 2) / (2. * pow(v, 3) * beta);
    ana_result(0.0) = ana_pref * (2. * v / (1. - pow(v, 2)) - 2. * atanh(v));
    ana_result(1, 1) =
        ana_pref * ((3. * pow(v, 3) - v) / pow(1 - pow(v, 2), 2) + atanh(v));
    ana_result(2, 2) = ana_result(1, 1);

    REQUIRE(approx_equal(num_result, ana_result, "reldiff", 10E-4));
  }

  SECTION("Option: fancy_I_kv_temp") {
    // Generate a Green's tensor with random attributes v and beta
    auto v = GENERATE(take(1, random(0., 1.)));
    auto omega = GENERATE(take(3, random(-1e1, 1e1)));
    auto beta = GENERATE(take(3, random(10e-10, 10e-12)));
    //  beta *= fabs(omega);
    double relerr = 1E-9;
    GreensTensorVacuum Greens(v, beta, relerr);

    // Create a struct with the integration options
    Options_GreensTensor opts;
    opts.class_pt = &Greens;
    opts.omega = omega;

    // Matrix storing the numerical integration
    cx_mat::fixed<3, 3> num_result(fill::zeros);

    // Matrix to store the analytic results
    cx_mat::fixed<3, 3> ana_result(fill::zeros);

    double ana_pref = 0;

    // Integration of the vacuum Green's tensor
    opts.fancy_I_temp = false;
    opts.fancy_I_kv_temp = true;
    Greens.integrate_1d_k(num_result, opts);

    // Computing the analytical result and storing it in analytic
    ana_pref = std::pow(omega, 3) / (6. * std::pow(v, 4) * beta);
    ana_result(0, 0) =
        ana_pref * (2 * (5 * pow(v, 3) - 3. * v) / pow(1. - pow(v, 2), 2) +
                    6. * std::atanh(v));
    ana_result(1, 1) = ana_pref * ((8. * pow(v, 3) - 3. * v - 13. * pow(v, 5)) /
                                       pow(pow(v, 2) - 1, 3) -
                                   3. * std::atanh(v));
    ana_result(2, 2) = ana_result(1, 1);

    REQUIRE(approx_equal(num_result, ana_result, "reldiff", 10E-5));
  };
};
