#include "Quaca.h"
#include "catch.hpp"
#include <armadillo>
#include <complex>

TEST_CASE("Vacuum Green's tensor constructors work as expected",
          "[GreensTensorVacuum]") {

  SECTION("Direct constructor") {
    auto v = GENERATE(3.2e-3,1.2e-1);
    auto beta = GENERATE(3.65,100.34);
    auto relerr = GENERATE(1e-8,1e-7);

    GreensTensorVacuum Greens(v, beta, relerr);

    REQUIRE(Greens.get_v() == v);
    REQUIRE(Greens.get_beta() == beta);
    REQUIRE(Greens.get_relerr() == relerr);
  };

  SECTION("json file constructor") {
    double v = 0.1;
    double beta = 5;
    double relerr = 1E-9;

    GreensTensorVacuum Greens("../data/test_files/GreensTensorVacuum.json");

    REQUIRE(Approx(Greens.get_v()).epsilon(1E-6) == v);
    REQUIRE(Approx(Greens.get_beta()).epsilon(1E-6) == beta);
    REQUIRE(Greens.get_relerr() == relerr);
  };
};

TEST_CASE("Integrand 1d k is correctly implemented", "[GreensTensorVacuum]") {
  // Generate a Green's tensor with random attributes v and beta
  auto v = GENERATE(1e-4,1e-8);
  auto beta = GENERATE(12.21,1.32,5.23);
  double relerr = 1E-9;
  GreensTensorVacuum Greens(v, beta, relerr);

  // Create the variables for the num_result, taking care, that
  //\omega^2 - k^2 >= 0 to stay in the non-trivial regime
  auto omega = GENERATE(-7.43,0.21,1.76);
  auto k_v = GENERATE(-.9,.1,.8);
  if (k_v < 0)
    k_v *= omega / (1 + v);
  if (k_v >= 0)
    k_v *= omega / (1 - v);

  // Check the integrand for all possible integration options
  double omega_kv = omega + v * k_v;
  double xi = pow(omega_kv, 2) - pow(k_v, 2);
  //Stores the prefactors of the integration options
  double factor = 1.;
  cx_mat::fixed<3,3> LHS(fill::zeros);
  cx_mat::fixed<3,3> RHS(fill::zeros);
  

  //Set the values of the analytical result
  LHS(0,0) = .5 * xi;
  LHS(1,1) = .5 * (pow(omega_kv, 2) - .5 * xi);
  LHS(2,2) = LHS(1,1);

  double result_x = LHS(0,0).real();
  double result_yz = LHS(1,1).real();

  //Matrix with only zero entries
  cx_mat::fixed<3,3> zero_mat(fill::zeros);

  // Create a struct with the integration options
  struct Options_GreensTensor opts;
  opts.class_pt = &Greens;
  opts.omega = omega;

  SECTION("Option: IM") {
    opts.fancy_complex = IM;
    for(int i = 0; i < 3; ++i) {
      for(int j = 0; j < 3; ++j) {
	opts.indices(0) = i;
	opts.indices(1) = j;
	RHS(i,j) = Greens.integrand_k(k_v, &opts);
      }
    }
    //Ensure that result is non-trivial
    REQUIRE(!approx_equal(RHS,zero_mat,"reldiff",1e-1));

    REQUIRE(approx_equal(LHS,RHS,"reldiff",1e-12));
  }

  SECTION("Option: IM, KV") {
    opts.fancy_complex = IM;
    opts.weight_function = KV;
    factor = k_v;
    for(int i = 0; i < 3; ++i) {
      for(int j = 0; j < 3; ++j) {
	opts.indices(0) = i;
	opts.indices(1) = j;
	RHS(i,j) = Greens.integrand_k(k_v, &opts);
      }
    }
    //Ensure that result is non-trivial
    REQUIRE(!approx_equal(RHS,zero_mat,"reldiff",1e-1));

    REQUIRE(approx_equal(factor*LHS,RHS,"reldiff",1e-12));
  };

  SECTION("Option: IM, TEMP") {
    opts.fancy_complex = IM;
    opts.weight_function = TEMP;
    factor = 1. / (1. - exp(-beta * omega_kv));
    for(int i = 0; i < 3; ++i) {
      for(int j = 0; j < 3; ++j) {
	opts.indices(0) = i;
	opts.indices(1) = j;
	RHS(i,j) = Greens.integrand_k(k_v, &opts);
      }
    }
    //Ensure that result is non-trivial
    REQUIRE(!approx_equal(RHS,zero_mat,"reldiff",1e-1));

    REQUIRE(approx_equal(factor*LHS,RHS,"reldiff",1e-12));
  };

  SECTION("Option: IM, KV_TEMP") {
    opts.fancy_complex = IM;
    opts.weight_function = KV_TEMP;
    factor = k_v / (1. - exp(-beta * omega_kv));
    for(int i = 0; i < 3; ++i) {
      for(int j = 0; j < 3; ++j) {
	opts.indices(0) = i;
	opts.indices(1) = j;
	RHS(i,j) = Greens.integrand_k(k_v, &opts);
      }
    }
    //Ensure that result is non-trivial
    REQUIRE(!approx_equal(RHS,zero_mat,"reldiff",1e-1));

    REQUIRE(approx_equal(factor*LHS,RHS,"reldiff",1e-12));
  };

  SECTION("Option: IM, NON_LTE") {
    opts.fancy_complex = IM;
    opts.weight_function = NON_LTE;
    factor =
        1. / (1. - exp(-beta * (omega_kv))) - 1. / (1. - exp(-beta * omega));
    for(int i = 0; i < 3; ++i) {
      for(int j = 0; j < 3; ++j) {
	opts.indices(0) = i;
	opts.indices(1) = j;
	RHS(i,j) = Greens.integrand_k(k_v, &opts);
      }
    }
    //Ensure that result is non-trivial
    if(approx_equal(RHS,zero_mat,"reldiff",1e-1))
    {
      std::cout << omega << '\t' << k_v << '\t' << beta << '\t' << v << std::endl;
      std::cout << RHS << std::endl;
    }
    REQUIRE(!approx_equal(RHS,zero_mat,"reldiff",1e-1));

    REQUIRE(approx_equal(factor*LHS,RHS,"reldiff",1e-12));
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
  opts.fancy_complex = IM;
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
  opts.fancy_complex = IM;
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
  opts.fancy_complex = IM;
  opts.class_pt = &Greens;

  // Create the matries storing the Green's tensors
  cx_mat::fixed<3, 3> Greens_lhs(fill::zeros);
  cx_mat::fixed<3, 3> Greens_rhs(fill::zeros);

  opts.omega = omega;
  Greens.integrate_k(Greens_lhs, opts);

  opts.omega = -omega;
  Greens.integrate_k(Greens_rhs, opts);

  REQUIRE(approx_equal(Greens_lhs, -Greens_rhs, "reldiff", 10E-5));
};

TEST_CASE("Test the integration routine", "[GreensTensorVacuum]") {

  SECTION("Option: IM") {
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
    opts.fancy_complex = IM;
    Greens.integrate_k(num_result, opts);

    // Computing the analytical result and storing it in analytic
    ana_pref = 2. / 3. * pow(omega, 3) / pow(1 - pow(v, 2), 2);
    ana_result(0, 0) = ana_pref;
    ana_result(1, 1) = ana_pref * (1 + pow(v, 2)) / (1 - pow(v, 2));
    ana_result(2, 2) = ana_pref * (1 + pow(v, 2)) / (1 - pow(v, 2));

    REQUIRE(approx_equal(num_result, ana_result, "reldiff", 10E-5));
  };

  SECTION("Option: IM, KV") {

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
    opts.fancy_complex = IM;
    opts.weight_function = KV;
    Greens.integrate_k(num_result, opts);

    // Computing the analytical result and storing it in analytic
    ana_pref = 2. / 3. * pow(omega, 4) * v / pow(1 - pow(v, 2), 3);
    ana_result(0.0) = ana_pref;
    ana_result(1, 1) = ana_pref * (2. + pow(v, 2)) / (1. - pow(v, 2));
    ana_result(2, 2) = ana_result(1, 1);

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
    opts.fancy_complex = IM;
    opts.weight_function = TEMP;
    Greens.integrate_k(num_result, opts);

    // Computing the analytical result and storing it in analytic
    ana_pref = pow(omega, 2) / (2. * pow(v, 3) * beta);
    ana_result(0.0) = ana_pref * (2. * v / (1. - pow(v, 2)) - 2. * atanh(v));
    ana_result(1, 1) =
        ana_pref * ((3. * pow(v, 3) - v) / pow(1 - pow(v, 2), 2) + atanh(v));
    ana_result(2, 2) = ana_result(1, 1);

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
    opts.fancy_complex = IM;
    opts.weight_function = KV_TEMP;
    Greens.integrate_k(num_result, opts);

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
