#include <armadillo>
#include <complex>

#include "Quaca.h"
#include "catch.hpp"

TEST_CASE("Integrated PolarizabilityNoBath fulfills the omega_cut much smaller "
          "than omega_a asymptote",
          "[PolarizabilityNoBath]") {
  // define greens tensor
  auto v = GENERATE(1e-4,1e-8);
  auto beta = GENERATE(1e-3,1.,1e2);
  double relerr_k = 1E-9;
  GreensTensorVacuum greens(v, beta, relerr_k);

  // define polarizability
  auto omega_a = GENERATE(0.21,1.7);
  auto alpha_zero = GENERATE(1e-8,1e-9);
  PolarizabilityNoBath pol(omega_a, alpha_zero, &greens);

  Options_Polarizability opts;
  opts.fancy_complex = IM;
  opts.class_pt = &pol;

  double omega_min = 0.0;
  double omega_max = 1e-3 * omega_a; // omega_max much smaller than omega_a
  double relerr = 1e-13;
  double abserr = 0;

  cx_mat::fixed<3,3> result(fill::zeros);
  cx_mat::fixed<3,3> asymp(fill::zeros);
  asymp(0,0) = alpha_zero * alpha_zero * pow(omega_max, 4) / 2.0 * 1.0 /
                  (3 * (1.0 - v * v) * (1.0 - v * v));
  asymp(1,1) =  alpha_zero * alpha_zero * pow(omega_max, 4) / 2.0 *
                  (1.0 + v * v) / (3 * pow((1.0 - v * v), 3));
  asymp(2,2) = asymp(1,1);




  // loop over indices
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      opts.indices(0) = i;
      opts.indices(1) = j;
      result(i,j) = pol.integrate_omega(opts, omega_min, omega_max, relerr, abserr);
    }
  }

  //Ensure non-trivial result
  REQUIRE(!result.is_zero());
  REQUIRE(!asymp.is_zero());

  REQUIRE(approx_equal(result,asymp, "reldiff", 1e-4));
  //Ensure that the total error is above the error due to the series expansio
  REQUIRE(approx_equal(result,asymp, "absdiff", pow(omega_max,2)));
}

TEST_CASE("Integrated PolarizabilityNoBath fulfills the omega_cut much larger "
          "than omega_a asymptote",
          "[PolarizabilityNoBath]") {
  // define greens tensor
  auto v = GENERATE(1e-4,1e-8);
  double beta = 1e5;
  double relerr_k = 1E-9;
  GreensTensorVacuum greens(v, beta, relerr_k);

  // define polarizability
  auto omega_a = GENERATE(0.21,1.5);
  double alpha_zero = 1e-10;
  PolarizabilityNoBath pol(omega_a, alpha_zero, &greens);

  Options_Polarizability opts;
  opts.fancy_complex = IM;
  opts.class_pt = &pol;

  double omega_min = 0.0;
  double omega_max = omega_a * 1e5;
  double relerr = 1e-13;
  double abserr = 0.;

  double asymp = alpha_zero * omega_a * M_PI / 2.0;
  cx_mat::fixed<3,3> result(fill::zeros);
  //create unitary matrix
  cx_mat::fixed<3,3> asymp_mat(fill::eye);
  asymp_mat *= asymp;


  // loop over indices
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      opts.indices(0) = i;
      opts.indices(1) = j;
      result(i,j) =
          pol.integrate_omega(opts, omega_min, omega_a - 1e-3, relerr, abserr);
      result(i,j) += pol.integrate_omega(opts, omega_a - 1e-3, omega_a + 1e-3,
                                    relerr, abserr);
      result(i,j) +=
          pol.integrate_omega(opts, omega_a + 1e-3, omega_max, relerr, abserr);
    }
  }

  //Ensure non-trivial results
  REQUIRE(!result.is_zero());
  REQUIRE(!asymp_mat.is_zero());

  REQUIRE(approx_equal(result,asymp_mat,"reldiff", 1e-4));
  //Ensure that the error is above the error due to the series expansin
  REQUIRE(approx_equal(result,asymp_mat,"absdiff",sqrt(alpha_zero)));
}
