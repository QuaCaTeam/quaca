#include <armadillo>
#include <complex>

#include "Quaca.h"
#include "catch.hpp"

TEST_CASE("PolarizabilityBath constructors work as expected",
          "[PolarizabilityBath]") {

  SECTION("Direct constructor") {
    auto mu = std::make_shared<OhmicMemoryKernel>(0.69420);
    PolarizabilityBath pol(1.3, 6E-9, mu, nullptr);
    REQUIRE(pol.get_omega_a() == 1.3);
    REQUIRE(pol.get_alpha_zero() == 6E-9);

    std::complex<double> test = pol.get_mu(3.0);
    REQUIRE(test.real() == 0.69420);
  }

  SECTION("json file constructor") {
    PolarizabilityBath pol("../data/test_files/PolarizabilityBath.json");
    REQUIRE(pol.get_omega_a() == 1.3);
    REQUIRE(pol.get_alpha_zero() == 6E-9);

    // test if we read memory kernel correctly
    std::complex<double> test = pol.get_mu(3.0);
    REQUIRE(test.real() == 0.69420);
  }
}

TEST_CASE("PolarizabilityBath integrand works as expected",
          "[PolarizabilityBath]") {
  // define greens tensor
  double v = 0.1;
  double beta = 10;
  double relerr_k = 1E-9;
  auto greens = std::make_shared<GreensTensorVacuum>(v, beta, relerr_k);

  // define polarizability
  double omega_a = 3.0;
  double alpha_zero = 2.4;
  double gamma = 2.0;
  auto mu = std::make_shared<OhmicMemoryKernel>(gamma);
  PolarizabilityBath pol(omega_a, alpha_zero, mu, greens);

  // frequency to evaluate
  double omega = 3.0;

  // calculate as a reference the normal way
  cx_mat::fixed<3, 3> alpha;
  pol.calculate_tensor(omega, alpha, IM);

  // calculation with calculate_tensor give the same result as integrand
  // Tensor to store the result
  cx_mat::fixed<3,3> alpha_int(fill::zeros);

  // loop over all indices
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      opts.indices(0) = i;
      opts.indices(1) = j;
      alpha_int(i,j) = pol.integrand_omega(omega, {(double) i, (double) j}, IM);
    }
  }
  //Ensure non-trivial results
  REQUIRE(!alpha.is_zero());
  REQUIRE(!alpha_int.is_zero());

  REQUIRE(approx_equal(alpha,alpha_int,"reldiff",1e-6));
}
