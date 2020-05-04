#include <armadillo>
#include <complex>

#include "Quaca.h"
#include "catch.hpp"

TEST_CASE("PolarizabilityNoBath constructors work as expected",
          "[PolarizabilityNoBath]") {

  SECTION("Direct constructor") {
    PolarizabilityBath pol(1.3, 6E-9, nullptr);
    REQUIRE(pol.get_omega_a() == 1.3);
    REQUIRE(pol.get_alpha_zero() == 6E-9);
    // Since no internal bath is implemented the get_mu should return 0
    // for any omega
    std::complex<double> test = pol.get_mu(3.0);
    REQUIRE(test.real() == 0.);
  }

  SECTION("json file constructor") {
    PolarizabilityBath pol("../data/test_files/PolarizabilityNoBath.json");
    REQUIRE(pol.get_omega_a() == 1.3);
    REQUIRE(pol.get_alpha_zero() == 6E-9);
    // Since no internal bath is implemented the get_mu should return 0
    // for any omega
    std::complex<double> test = pol.get_mu(3.0);
    REQUIRE(test.real() == 0.);
  }
}

TEST_CASE("PolarizabilityNoBath integrand works as expected",
          "[PolarizabilityNoBath]") {
  // define greens tensor
  double v = 0.1;
  double beta = 10;
  double relerr_k = 1E-9;
  auto greens = std::make_shared<GreensTensorVacuum>(v, beta, relerr_k);

  // define polarizability
  double omega_a = 3.0;
  double alpha_zero = 2.4;
  PolarizabilityBath pol(omega_a, alpha_zero, greens);

  // frequency to evaluate
  double omega = 3.0;

  // calculate as a reference the normal way
  cx_mat::fixed<3, 3> alpha;
  pol.calculate_tensor(omega, alpha, IM);

  // calculation with calculate_tensor give the same result as integrand
  // Matrix to store the result
  cx_mat::fixed<3, 3> alpha_int(fill::zeros);

  // loop over all indices
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      alpha_int(i, j) = pol.integrand_omega(omega, {(double)i, (double)j}, IM);
    }
  }

  // Ensure non-trivial results
  REQUIRE(!alpha.is_zero());
  REQUIRE(!alpha_int.is_zero());

  if (!approx_equal(alpha, alpha_int, "reldiff", 1e-6)) {
    std::cout << alpha << alpha_int << alpha - alpha_int << std::endl;
  }
  REQUIRE(approx_equal(alpha, alpha_int, "reldiff", 1e-6));
};
