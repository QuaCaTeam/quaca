#include <armadillo>
#include <complex>

#include "Quaca.h"
#include "catch.hpp"

TEST_CASE("PolarizabilityBath constructors work as expected",
          "[PolarizabilityBath]") {

  SECTION("Direct constructor") {
    OhmicMemoryKernel mu(0.69420);
    PolarizabilityBath pol(1.3, 6E-9, &mu, NULL);
    REQUIRE(pol.get_omega_a() == 1.3);
    REQUIRE(pol.get_alpha_zero() == 6E-9);

    std::complex<double> test = pol.get_mu(3.0);
    REQUIRE(test.real() == 0.69420);
  };

  SECTION("ini file constructor") {
    PolarizabilityBath pol("../data/test_files/PolarizabilityBath.ini");
    REQUIRE(pol.get_omega_a() == 1.3);
    REQUIRE(pol.get_alpha_zero() == 6E-9);

    // test if we read memory kernel correctly
    std::complex<double> test = pol.get_mu(3.0);
    REQUIRE(test.real() == 0.69420);
  };
};

TEST_CASE("PolarizabilityBath integrand works as expected",
          "[PolarizabilityBath]") {
  // define greens tensor
  double v = 0.1;
  double beta = 10;
  double relerr_k = 1E-9;
  GreensTensorVacuum greens(v, beta, relerr_k);

  // define polarizability
  double omega_a = 3.0;
  double alpha_zero = 2.4;
  double gamma = 2.0;
  OhmicMemoryKernel mu(gamma);
  PolarizabilityBath pol(omega_a, alpha_zero, &mu, &greens);

  // frequency to evaluate
  double omega = 3.0;

  // define options struct for integrand
  Options_Polarizability opts;
  opts.class_pt = &pol;
  opts.fancy_complex = IM;
  opts.indices(0) = 0;
  opts.indices(1) = 0;

  // calculate as a reference the normal way
  cx_mat::fixed<3, 3> alpha;
  opts.omega = omega;
  pol.calculate_tensor(alpha, opts);

  // calculation with calculate_tensor give the same result as integrand

  // loop over all indices
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      opts.indices(0) = i;
      opts.indices(1) = j;
      pol.calculate_tensor(alpha, opts);
      REQUIRE(pol.integrand_omega(omega, &opts) ==
              alpha(opts.indices(0), opts.indices(1)));
    }
  }
};
