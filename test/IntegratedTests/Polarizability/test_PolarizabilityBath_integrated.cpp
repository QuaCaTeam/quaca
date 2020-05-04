#include <complex>

#include "Quaca.h"
#include "catch.hpp"

TEST_CASE("Integrated PolarizabilityBath fulfills the omega_cut much smaller "
          "than omega_a asymptote",
          "[PolarizabilityBath]") {
  // define greens tensor
  auto v = GENERATE(1e-4, 1e-8);
  auto beta = GENERATE(1e-4, 1e4);
  double relerr_k = 1e-9;
  auto greens = std::make_shared<GreensTensorVacuum>(v, beta, relerr_k);

  // define polarizability
  auto omega_a = GENERATE(1.6);
  auto alpha_zero = GENERATE(1e-8);
  auto gamma = GENERATE(0.32, 12.43);
  auto mu = std::make_shared<OhmicMemoryKernel>(gamma);
  Polarizability pol(omega_a, alpha_zero, mu, greens);

  double omega_min = 0.0;
  double fact = 1e-3;
  double relerr = 1e-12;
  double abserr = 0;
  double omega_max = fact * omega_a;

  cx_mat::fixed<3, 3> result(fill::zeros);
  cx_mat::fixed<3, 3> asymp(
      fill::zeros); // double for result and asymptotic value

  // loop over indices
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      result(i, j) = pol.integrate_omega({(double)i, (double)j}, IM, omega_min,
                                         omega_max, relerr, abserr);
    }
  }

  asymp(0, 0) = alpha_zero * alpha_zero * pow(omega_max, 4) / 2.0 * 1.0 /
                    (3 * (1.0 - v * v) * (1.0 - v * v)) +
                alpha_zero * gamma * fact * fact / (2.0);
  asymp(1, 1) = alpha_zero * alpha_zero * pow(omega_max, 4) / 2.0 *
                    (1.0 + v * v) / (3 * pow((1.0 - v * v), 3)) +
                alpha_zero * gamma * fact * fact / (2.0);
  asymp(2, 2) = asymp(1, 1);

  // Ensure non-trivial results
  REQUIRE(!result.is_zero());
  REQUIRE(!asymp.is_zero());

  REQUIRE(approx_equal(result, asymp, "reldiff", 1e-3));
  // Ensure that the absolute error is below the error due to the series
  // expansio
  REQUIRE(approx_equal(result, asymp, "absdiff", pow(omega_max, 2)));
};

TEST_CASE("Integrated PolarizabilityBath fulfills the omega_cut much larger "
          "than omega_a asymptote",
          "[PolarizabilityBath]") {
  // define greens tensor
  double v = 1e-1;
  double beta = 1e3;
  double relerr_k = 1E-12;
  auto greens = std::make_shared<GreensTensorVacuum>(v, beta, relerr_k);

  // define polarizability
  double omega_a = 4.0;
  double alpha_zero = 1e-8;
  double gamma = 1.0;
  auto mu = std::make_shared<OhmicMemoryKernel>(gamma);
  Polarizability pol(omega_a, alpha_zero, mu, greens);

  double omega_min = 0.0;
  double omega_max = omega_a * 1e4;
  double relerr = 1e-15;
  double abserr = 0;

  double asymp =
      alpha_zero * omega_a * omega_a / 2.0 *
      (M_PI -
       2.0 * atan((gamma * gamma - 2.0 * omega_a * omega_a) /
                  (gamma * sqrt(4.0 * omega_a * omega_a - gamma * gamma)))) /
      sqrt(4 * omega_a * omega_a - gamma * gamma);
  cx_mat::fixed<3, 3> result(fill::zeros);
  // Create unitary matrix
  cx_mat::fixed<3, 3> asymp_mat(fill::eye);
  asymp_mat *= asymp;

  // loop over indices
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      result(i, j) = pol.integrate_omega({(double)i, (double)j}, IM, omega_min,
                                         omega_a - 1e-3, relerr, abserr);
      result(i, j) +=
          pol.integrate_omega({(double)i, (double)j}, IM, omega_a - 1e-3,
                              omega_a + 1e-3, relerr, abserr);
      result(i, j) +=
          pol.integrate_omega({(double)i, (double)j}, IM, omega_a + 1e-3,
                              omega_max, relerr, abserr);
    }
  }
  // Ensure non-trivial result
  REQUIRE(!result.is_zero());
  REQUIRE(!asymp_mat.is_zero());

  if (!approx_equal(result, asymp_mat, "reldiff", 1e-4)) {
    std::cout << result << result / asymp_mat << std::endl;
  }
  REQUIRE(approx_equal(result, asymp_mat, "reldiff", 1e-4));
  // Ensure that the error is below is error-bound due to the series expansion
  REQUIRE(approx_equal(result, asymp_mat, "absdiff", sqrt(alpha_zero)));
};
