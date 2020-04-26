#include <complex>

#include "Quaca.h"
#include "catch.hpp"

TEST_CASE("Integrated PolarizabilityBath fulfills the omega_cut much smaller "
          "than omega_a asymptote",
          "[PolarizabilityBath]") {
  // define greens tensor
  auto v = GENERATE(take(3, random(0.0, 1.0)));
  auto beta = GENERATE(take(3, random(1e-4, 1e4)));
  double relerr_k = 1e-9;
  auto greens = std::make_shared<GreensTensorVacuum>(v, beta, relerr_k);

  // define polarizability
  auto omega_a = GENERATE(take(3, random(1e-1, 1e1)));
  auto alpha_zero = GENERATE(take(3, random(0.0, 1e-8)));
  auto gamma = GENERATE(take(3, random(0.0, 1e2)));
  auto mu = std::make_shared<OhmicMemoryKernel>(gamma);
  PolarizabilityBath pol(omega_a, alpha_zero, mu, greens);

  double omega_min = 0.0;
  double fact = 1e-3;
  double omega_max = fact * omega_a;
  double relerr = 1e-12;
  double abserr = 0;

  double result, asymp; // double for result and asymptotic value
  double toterr;

  // loop over indices
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      result = pol.integrate_omega({(double)i, (double)j}, IM, omega_min,
                                   omega_max, relerr, abserr);

      /*
       * error tolerance includes absolute error from integration, which is
       * result*relerr and error from series expansion in omega_a. we roughly
       * estimate the series error with (omega_max)^2
       */
      toterr = result * relerr + omega_max * omega_max;

      // check diagonal entries
      if (i == j) {
        if (i == 0) {
          asymp = alpha_zero * alpha_zero * pow(omega_max, 4) / 2.0 * 1.0 /
                      (3 * (1.0 - v * v) * (1.0 - v * v)) +
                  alpha_zero * gamma * fact * fact / (2.0);
          REQUIRE(Approx(result).margin(toterr) == asymp);
        } else {
          asymp = alpha_zero * alpha_zero * pow(omega_max, 4) / 2.0 *
                      (1.0 + v * v) / (3 * pow((1.0 - v * v), 3)) +
                  alpha_zero * gamma * fact * fact / (2.0);
          REQUIRE(Approx(result).margin(toterr) == asymp);
        }
      } else {
        REQUIRE(result == 0); // off-diagonals are zero
      }
    }
  }
}

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
  double alpha_zero = 1e-4;
  double gamma = 1.0;
  auto mu = std::make_shared<OhmicMemoryKernel>(gamma);
  PolarizabilityBath pol(omega_a, alpha_zero, mu, greens);

  double omega_min = 0.0;
  double omega_max = omega_a * 1e4;
  double relerr = 1e-15;
  double abserr = 0;

  double asymp =
      alpha_zero * omega_a * M_PI / 2.0 *
      (M_PI -
       2.0 * atan((gamma * gamma - 2.0 * omega_a * omega_a) /
                  (gamma * sqrt(4.0 * omega_a * omega_a - gamma * gamma)))) /
      sqrt(4 * omega_a * omega_a - gamma * gamma);
  double result;

  // loop over indices
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      result = pol.integrate_omega({(double)i, (double)j}, IM, omega_min,
                                   omega_a - 1e-3, relerr, abserr);
      result += pol.integrate_omega({(double)i, (double)j}, IM, omega_a - 1e-3,
                                    omega_a + 1e-3, relerr, abserr);
      result += pol.integrate_omega({(double)i, (double)j}, IM, omega_a + 1e-3,
                                    omega_max, relerr, abserr);

      /*
       * error tolerance includes absolute error from integration, which is
       * result*relerr we estimate the truncation error to sqrt(alpha_zero)
       */
      double toterr = result * relerr + sqrt(alpha_zero);

      // check diagonal entries
      if (i == j) {
        REQUIRE(Approx(result).margin(toterr) == asymp);
      } else {
        REQUIRE(result == 0); // off-diagonals are zero
      }
    }
  }
}
